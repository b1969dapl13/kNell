# File: export_json.sage.py
"""
elliptic_curves.db から個別曲線JSONとindex.jsonを生成する (v2)。

使い方:
    sage -python export_json.sage.py --db elliptic_curves.db --out data/
    sage -python export_json.sage.py --db elliptic_curves.db --out data/ --target-points 100
"""
import cmath
import math
import numpy as np

import argparse
import json
import os
import re
import sqlite3
from fractions import Fraction
from itertools import product
from math import ceil, gcd
from sage.all import EllipticCurve, QQ, Integer


# ---------------------------------------------------------------
# 逆変換: (y1, y2) -> (a:b:c)
# ---------------------------------------------------------------
def inverse_transform(y1, y2, k, N):
    """最終ワイエルシュトラス形の点 (y1, y2) から元の射影座標 (a:b:c) を復元。"""
    k = Fraction(k)
    N = Fraction(N)

    # Step 9 inverse
    x1_star = 2 * k * (k**3 + 1) * N + 4 * k**2
    x1 = y1 + x1_star
    x2 = y2

    # Step 6 inverse
    denom6 = k**3 + 1
    if denom6 == 0:
        return None
    w1 = -x1 / (2 * denom6)
    w2 = -x2 / denom6

    # Step 5 inverse
    beta2 = w1**2 - k**2 * N**2 + 4 * k
    beta1 = (-2 * (k**3 + 1) * w1 * N - 2 * k * w1
             - 2 * (k**4 + k) * N**2 - 2 * k**2 * N + 4 * k**3 + 4)
    if beta2 == 0:
        return None
    v1 = (w2 - beta1) / (2 * beta2)

    # Step 4 inverse
    L_v1 = (k**3 + 1) - ((k**3 + 1) * N + k) * v1
    v2 = L_v1 + v1**2 * w1

    # Step 2 inverse
    u2 = v1
    alpha2 = -k**2 * u2 * N + (-k**3 * u2 + k**2 - u2)
    alpha1 = ((k * u2**2 - k**3 * u2 + u2) * N
              + (2 * k**2 * u2**2 + k**3 - k * u2 - 1))
    if alpha2 == 0:
        return None
    u1 = (v2 - alpha1) / (2 * alpha2)

    # Step 1 inverse
    a = u2 - k * u1
    b = u1
    c = Fraction(1)

    def lcm(x, y):
        return x * y // gcd(x, y)

    L = lcm(lcm(a.denominator, b.denominator), c.denominator)
    a_int = int(a * L)
    b_int = int(b * L)
    c_int = int(c * L)
    g = gcd(gcd(abs(a_int), abs(b_int)), abs(c_int))
    if g > 0:
        a_int //= g
        b_int //= g
        c_int //= g
    return (a_int, b_int, c_int)


# ---------------------------------------------------------------
# 符号パターン
# ---------------------------------------------------------------
def sign_of(v):
    if v > 0:
        return 1
    if v < 0:
        return -1
    return 0


def sign_pattern(abc):
    """[a,b,c] -> [sign(a), sign(b), sign(c)]"""
    if abc is None:
        return None
    return [sign_of(v) for v in abc]


def is_all_same_sign(abc):
    """全符号一致 (全て正 or 全て負, ゼロは不可)"""
    if abc is None:
        return False
    s = sign_pattern(abc)
    if 0 in s:
        return False
    return s[0] == s[1] == s[2]


# ---------------------------------------------------------------
# 点数目標から M を決定
# ---------------------------------------------------------------
def determine_M(rank, torsion_order, target=100, M_max=10):
    if rank == 0:
        return 0
    t = max(torsion_order, 1)
    needed = max(1, ceil(target / t))
    M = 0
    while (2 * M + 1) ** rank < needed and M < M_max:
        M += 1
    return M


def coord_digits(x):
    f = Fraction(x)
    return max(len(str(abs(f.numerator))), len(str(abs(f.denominator))))


# ---------------------------------------------------------------
# ねじれ部分群の構造化
# ---------------------------------------------------------------
def parse_torsion_structure(struct_str):
    """例: '(2, 6)' -> [2, 6]"""
    if not struct_str:
        return []
    m = re.search(r'\(([0-9,\s]+)\)', struct_str)
    if not m:
        return []
    return [int(s.strip()) for s in m.group(1).split(',') if s.strip()]


def enumerate_torsion_points(E, torsion_gens, orders):
    if not torsion_gens or not orders:
        return []
    if len(torsion_gens) != len(orders):
        result = []
        for T in E.torsion_points():
            if T.is_zero():
                continue
            result.append((None, T))
        return result

    coeff_ranges = [range(o) for o in orders]
    seen = {}  # (x, y) -> (coeffs, point) で重複除去
    for coeffs in product(*coeff_ranges):
        if all(c == 0 for c in coeffs):
            continue
        try:
            T = sum((c * g for c, g in zip(coeffs, torsion_gens)), E(0))
        except Exception:
            continue
        if T.is_zero():
            continue
        key = (T.xy()[0], T.xy()[1])
        if key not in seen:
            seen[key] = (list(coeffs), T)
    return list(seen.values())


def is_torsion_generator(torsion_coeffs):
    """torsion_coeffs が単位ベクトル (生成元単独) か"""
    if torsion_coeffs is None:
        return False
    nonzero = [(i, c) for i, c in enumerate(torsion_coeffs) if c != 0]
    return len(nonzero) == 1 and nonzero[0][1] == 1


# ---------------------------------------------------------------
# 点ロール判定
# ---------------------------------------------------------------
def classify_point(free_coeffs, torsion_coeffs):
    """
    free_coeffs: [m_1, ..., m_r] or [] (rank=0)
    torsion_coeffs: [c_1, ...] or None

    Returns one of:
      "free_generator"     : +P_i 単独 (他成分0)
      "combination"        : 自由部分のみで複数項 or -P_i
      "torsion_generator"  : T_j 単独
      "torsion_other"      : ねじれ部分の合成
      "mixed"              : 自由 + ねじれ
    """
    free_nonzero = [c for c in free_coeffs if c != 0] if free_coeffs else []
    has_torsion = torsion_coeffs is not None and any(c != 0 for c in torsion_coeffs)

    if not free_nonzero and not has_torsion:
        # ありえない（単位元は除外済み）
        return "combination"

    if not free_nonzero and has_torsion:
        if is_torsion_generator(torsion_coeffs):
            return "torsion_generator"
        return "torsion_other"

    if free_nonzero and not has_torsion:
        if len(free_nonzero) == 1 and free_nonzero[0] == 1:
            return "free_generator"
        return "combination"

    # 自由 + ねじれ
    return "mixed"


def point_label(free_coeffs, torsion_coeffs):
    """人間可読ラベル: '+P_1 + 2P_2 + T_1' など"""
    parts = []
    if free_coeffs:
        for i, c in enumerate(free_coeffs):
            if c == 0:
                continue
            sign = "+" if c > 0 else "-"
            mag = abs(c)
            if mag == 1:
                parts.append(f"{sign}P_{i+1}")
            else:
                parts.append(f"{sign}{mag}P_{i+1}")
    if torsion_coeffs:
        for j, c in enumerate(torsion_coeffs):
            if c == 0:
                continue
            if c == 1:
                parts.append(f"+T_{j+1}")
            else:
                parts.append(f"+{c}T_{j+1}")
    if not parts:
        return "O"
    s = " ".join(parts)
    # 先頭が "+" なら除去
    if s.startswith("+"):
        s = s[1:].lstrip()
    return s


# ---------------------------------------------------------------
# y₁範囲の計算（plot_yj.txtの手法を移植）
# ---------------------------------------------------------------
def compute_y1_ranges(k, N, a_samples=400, a_max=1e8):
    """
    各ブランチ j=0,1,2 について、bⱼ(a) > 0 かつ a > 0 となる範囲での y₁ の [min, max] を計算
    """
    from fractions import Fraction
    import cmath
    
    k_f = float(k)
    N_f = float(N)
    
    # カルダノ公式で bⱼ(a) を計算
    def A3(a): return -k_f
    def A2(a): return k_f**2*N_f*a - k_f**2 + k_f*N_f - a
    def A1(a): return k_f**3*N_f*a - k_f**2*a**2 + k_f*N_f*a**2 + k_f**2*N_f - 3*k_f*a + N_f*a - 1
    def A0(a): return k_f**2*N_f*a**2 - k_f*a**3 - k_f**2*a + k_f*N_f*a - a**2 - k_f
    
    omega = cmath.exp(2j*cmath.pi/3.0)
    
    def bj_func(a, j):
        """ブランチ j の bⱼ(a) を返す"""
        a = complex(a)
        A3v = complex(A3(a))
        A2v = complex(A2(a))
        A1v = complex(A1(a))
        A0v = complex(A0(a))
        
        if abs(A3v) < 1e-15:
            return None
        
        p2 = A2v/A3v
        p1 = A1v/A3v
        p0 = A0v/A3v
        shift = p2/3.0
        p = p1 - p2*p2/3.0
        q = 2.0*p2**3/27.0 - p2*p1/3.0 + p0
        s = cmath.sqrt((q/2.0)**2 + (p/3.0)**3)
        C = (-q/2.0 + s)**(1.0/3.0)
        
        if abs(C) < 1e-15:
            return -shift
        
        wj = omega**j
        t = wj*C - p/(3.0*wj*C)
        return t - shift
    
    def y1_ab(a, b):
        """y₁(a, b) を計算"""
        a = complex(a)
        b = complex(b)
        num = ((-4*k_f**5 - 4*k_f**4*N_f - 8*k_f**2 - 4*k_f*N_f)*a**2
               + (-4*k_f**6 - 4*k_f**5*N_f - 8*k_f**3 - 4*k_f**2*N_f + 4)*a*b
               + (-4*k_f**3*N_f - 4*N_f)*a
               + (4*k_f)*b**2
               + (-4*k_f**5 - 4*k_f**4*N_f - 4*k_f**2 - 4*k_f*N_f)*b
               + (4*k_f**3 + 4))
        den = (k_f*b + a)**2
        if abs(den) < 1e-12:
            return None
        return num/den
    
    # a のサンプリング（対数スケール + 線形スケール）
    try:
        import numpy as np
        a_vals = np.unique(np.concatenate([
            np.geomspace(1e-6, a_max, a_samples),
            np.linspace(1e-6, 1000, a_samples),
            np.linspace(1e-6, 10, a_samples // 2)
        ]))
        a_vals = a_vals[a_vals > 0]
        a_vals = list(a_vals)
    except Exception as e_np:
        # numpy が使えない場合は簡易サンプリング
        print(f"    numpy unavailable for k={k}, N={N}: {e_np}, using fallback sampling")
        a_vals = []
        # 対数スケール
        for i in range(a_samples):
            a_vals.append(1e-6 * (10 ** (i * 14.0 / a_samples)))
        # 線形スケール（小さい範囲を密に）
        for i in range(a_samples):
            a_vals.append(1e-6 + i * (1000.0 / a_samples))
        for i in range(a_samples // 2):
            a_vals.append(1e-6 + i * (10.0 / (a_samples // 2)))
        a_vals = sorted(set(a_vals))
    
    ranges = []
    for j in range(3):
        valid_y1 = []
        for a_val in a_vals:
            if a_val <= 0:
                continue
            try:
                b = bj_func(a_val, j)
                if b is None:
                    continue
                # 実数かつ正の条件を両方満たす
                if abs(b.imag) < 1e-6 and b.real > 1e-9:
                    y1 = y1_ab(a_val, b.real)
                    if y1 is not None and abs(y1.imag) < 1e-6:
                        valid_y1.append(y1.real)
            except Exception:
                continue
        
        if valid_y1:
            ranges.append({
                "branch": j,
                "y1_min": float(min(valid_y1)),
                "y1_max": float(max(valid_y1))
            })
        else:
            ranges.append({
                "branch": j,
                "y1_min": None,
                "y1_max": None
            })
    
    if all(r["y1_min"] is None for r in ranges):
        print(f"    Warning: k={k}, N={N} - no valid y1 ranges found for any branch")
        print(f"             Tested {len(a_vals)} a-values, conditions: a>0, b.real>1e-9, b.imag<1e-6, y1.imag<1e-6")
    
    return ranges

# ---------------------------------------------------------------
# y1 帯 (Cardano 閉形式 3 分岐) の計算
# ---------------------------------------------------------------
def _make_coeff_funcs(k, N):
    k = float(k); N = float(N)
    def A3(a): return -k
    def A2(a): return k**2*N*a - k**2 + k*N - a
    def A1(a): return (k**3*N*a - k**2*a**2 + k*N*a**2
                       + k**2*N - 3*k*a + N*a - 1)
    def A0(a): return (k**2*N*a**2 - k*a**3 - k**2*a
                       + k*N*a - a**2 - k)
    return A3, A2, A1, A0


def _make_branch_funcs(k, N):
    A3f, A2f, A1f, A0f = _make_coeff_funcs(k, N)
    omega = cmath.exp(2j * cmath.pi / 3.0)

    def make_bj(j):
        wj = omega ** j

        def bj(a):
            a = complex(a)
            A3 = complex(A3f(a))
            A2 = complex(A2f(a))
            A1 = complex(A1f(a))
            A0 = complex(A0f(a))
            if A3 == 0:
                return complex('nan')
            p2 = A2 / A3
            p1 = A1 / A3
            p0 = A0 / A3
            shift = p2 / 3.0
            p = p1 - p2 * p2 / 3.0
            q = 2.0 * p2**3 / 27.0 - p2 * p1 / 3.0 + p0
            s = cmath.sqrt((q / 2.0)**2 + (p / 3.0)**3)
            C = (-q / 2.0 + s) ** (1.0 / 3.0)
            if C == 0:
                return -shift
            t = wj * C - p / (3.0 * wj * C)
            return t - shift
        return bj

    return [make_bj(0), make_bj(1), make_bj(2)]


def _y1_ab(a, b, k, N):
    a = complex(a); b = complex(b)
    k = float(k); N = float(N)
    num = ((-4*k**5 - 4*k**4*N - 8*k**2 - 4*k*N) * a**2
           + (-4*k**6 - 4*k**5*N - 8*k**3 - 4*k**2*N + 4) * a*b
           + (-4*k**3*N - 4*N) * a
           + (4*k) * b**2
           + (-4*k**5 - 4*k**4*N - 4*k**2 - 4*k*N) * b
           + (4*k**3 + 4))
    den = (k * b + a) ** 2
    if den == 0:
        return complex('nan')
    return num / den


def compute_y1_bands(k, N,
                    tol_imag=1e-7,
                    tol_den=1e-8,
                    tol_y1_imag=1e-6):
    """各分岐 j=0,1,2 について、b_j(a) が実かつ正となる a での
    y1(a, b_j(a)) の [min, max] を返す。空なら None を入れる。"""
    try:
        b_funcs = _make_branch_funcs(k, N)
    except Exception:
        return []

    a_samples = sorted(set(
        list(np.geomspace(1e-4, 1e8, 4000))
        + list(np.linspace(1e-4, 5.0, 1500))
        + list(np.linspace(1e-4, 0.5, 1500))
    ))

    bands = []
    for j in range(3):
        y1_vals = []
        for a in a_samples:
            try:
                bj = b_funcs[j](a)
            except Exception:
                continue
            if not (math.isfinite(bj.real) and math.isfinite(bj.imag)):
                continue
            if abs(bj.imag) >= tol_imag:
                continue
            b_real = bj.real
            if b_real <= 0:
                continue
            den = k * b_real + a
            if abs(den) <= tol_den:
                continue
            try:
                y1v = _y1_ab(a, b_real, k, N)
            except Exception:
                continue
            if not (math.isfinite(y1v.real) and math.isfinite(y1v.imag)):
                continue
            if abs(y1v.imag) >= tol_y1_imag:
                continue
            y1_vals.append(y1v.real)

        if y1_vals:
            bands.append({
                "branch": j,
                "y1_min": float(min(y1_vals)),
                "y1_max": float(max(y1_vals)),
                "count": len(y1_vals),
            })
        else:
            bands.append({
                "branch": j,
                "y1_min": None,
                "y1_max": None,
                "count": 0,
            })
    return bands

# ---------------------------------------------------------------
# 1曲線のエクスポート
# ---------------------------------------------------------------
def export_curve(row, target_points, digit_limit):
    k = row["k"]
    N = row["n"]
    status = row["status"]
    out = {
        "k": int(k),
        "n": int(N),
        "status": status,
    }

    if status != "ok":
        if row["note"]:
            out["note"] = row["note"]
        return out

    ainvs_raw = json.loads(row["a_invariants"]) if row["a_invariants"] else None
    if ainvs_raw is None:
        out["status"] = "error"
        out["note"] = "missing a_invariants"
        return out

    ainvs = [Integer(int(s)) if isinstance(s, str) else Integer(s)
             for s in ainvs_raw]

    out["a_invariants"] = [str(a) for a in ainvs]
    out["discriminant"] = row["discriminant"]
    out["j_invariant"] = row["j_invariant"]
    out["conductor"] = row["conductor"]
    out["rank"] = row["rank"]
    out["rank_bounds"] = (json.loads(row["rank_bounds"])
                          if row["rank_bounds"] else None)
    out["torsion_order"] = row["torsion_order"]
    out["torsion_struct_raw"] = row["torsion_struct"]
    out["torsion_structure"] = parse_torsion_structure(row["torsion_struct"])
    out["proof_note"] = "gens computed with proof=False"

    # α, β
    k_i = Integer(k)
    N_i = Integer(N)
    alpha = (k_i**3 + 1) * N_i + 3 * k_i
    beta = k_i**2 * (k_i**3 + 1) * N_i + (k_i**6 + 3 * k_i**3 + 1)
    out["alpha"] = str(alpha)
    out["beta"] = str(beta)

    try:
        E = EllipticCurve(QQ, [QQ(s) for s in ainvs_raw])
    except Exception as e:
        out["status"] = "error"
        out["note"] = f"EllipticCurve construction failed: {e}"
        return out

    # 自由部分の生成元
    gens_raw = json.loads(row["gens"]) if row["gens"] else []
    try:
        free_gens = [E(QQ(gx), QQ(gy)) for gx, gy in gens_raw]
    except Exception as e:
        out["status"] = "error"
        out["note"] = f"generator construction failed: {e}"
        return out

    rank = len(free_gens)
    out["free_generators"] = [
        {
            "index": i + 1,
            "label": f"P_{i+1}",
            "xy": [str(g.xy()[0]), str(g.xy()[1])],
        }
        for i, g in enumerate(free_gens)
    ]

    # ねじれ部分群の生成元
    try:
        T_subgroup = E.torsion_subgroup()
        torsion_gens_raw = T_subgroup.gens()
        torsion_gens = []
        for tg in torsion_gens_raw:
            try:
                pt = E(tg)
            except Exception:
                pt = tg.element() if hasattr(tg, 'element') else tg
            if pt.is_zero():
                continue
            torsion_gens.append(pt)
    except Exception:
        torsion_gens = []

    # === 修正ここから ===
    # Sage の T.invariants() と T.gens() は順序が一致しない場合があるため、
    # 各生成元の実位数を直接取る
    orders = [int(g.order()) for g in torsion_gens]
    # === 修正ここまで ===

    out["torsion_generators"] = [
        {
            "index": j + 1,
            "label": f"T_{j+1}",
            "order": int(orders[j]),
            "xy": [str(tg.xy()[0]), str(tg.xy()[1])],
        }
        for j, tg in enumerate(torsion_gens)
    ]
    
    # y1 帯 (Cardano 閉形式)
    try:
        out["y1_bands"] = compute_y1_bands(k, N)
    except Exception as e:
        out["y1_bands"] = []
        out["y1_bands_note"] = f"compute failed: {e}"

    # 全ねじれ点 (生成元の整数結合)
    torsion_pts = enumerate_torsion_points(E, torsion_gens, orders)

    # 自由部分の係数を決定
    M = determine_M(rank, max(row["torsion_order"] or 1, 1),
                    target=target_points)
    out["M"] = M

    points = []

    # まず純粋ねじれ点 (自由部分はゼロ)
    for torsion_coeffs, T in torsion_pts:
        x, y = T.xy()
        x_f = Fraction(int(x.numerator()), int(x.denominator()))
        y_f = Fraction(int(y.numerator()), int(y.denominator()))
        abc = inverse_transform(x_f, y_f, k, N)
        free_coeffs = [0] * rank
        role = classify_point(free_coeffs, torsion_coeffs)
        label = point_label(free_coeffs, torsion_coeffs)
        sp = sign_pattern(abc)
        pt = {
            "free_coeffs": free_coeffs,
            "torsion_coeffs": torsion_coeffs,
            "role": role,
            "label": label,
            "xy": [str(x), str(y)],
            "abc": [str(v) for v in abc] if abc else None,
            "sign_pattern": sp,
            "all_same_sign": is_all_same_sign(abc),
            "order": int(T.order()),
        }
        points.append(pt)

    # 自由部分 + (ねじれを足したバリアント)
    if rank > 0:
        coeff_ranges = [range(-M, M + 1)] * rank
        torsion_variants = [(None, None)] + [(tc, T) for tc, T in torsion_pts]

        for coeffs in product(*coeff_ranges):
            if all(c == 0 for c in coeffs):
                continue
            try:
                P = sum((c * g for c, g in zip(coeffs, free_gens)), E(0))
            except Exception:
                continue

            for torsion_coeffs, T in torsion_variants:
                if T is None:
                    Q = P
                    tc_record = None
                else:
                    Q = P + T
                    tc_record = torsion_coeffs

                if Q.is_zero():
                    continue
                try:
                    x, y = Q.xy()
                except Exception:
                    continue

                x_f = Fraction(int(x.numerator()), int(x.denominator()))
                y_f = Fraction(int(y.numerator()), int(y.denominator()))

                if max(coord_digits(x_f), coord_digits(y_f)) > digit_limit:
                    continue

                abc = inverse_transform(x_f, y_f, k, N)
                free_coeffs = list(coeffs)
                role = classify_point(free_coeffs, tc_record)
                label = point_label(free_coeffs, tc_record)
                sp = sign_pattern(abc)
                pt = {
                    "free_coeffs": free_coeffs,
                    "torsion_coeffs": tc_record,
                    "role": role,
                    "label": label,
                    "xy": [str(x), str(y)],
                    "abc": [str(v) for v in abc] if abc else None,
                    "sign_pattern": sp,
                    "all_same_sign": is_all_same_sign(abc),
                }
                points.append(pt)

    out["points"] = points
    out["points_count"] = len(points)
    out["has_same_sign_point"] = any(p["all_same_sign"] for p in points)
    
    # y₁範囲の計算を追加
    try:
        y1_ranges = compute_y1_ranges(k, N)
        out["y1_ranges"] = y1_ranges
    except Exception as e:
        out["y1_ranges"] = []
        # 失敗してもエラーにはしない（オプショナル機能）
        print(f"    Warning: y1_ranges computation failed for k={k}, N={N}: {e}")
    
    return out


# ---------------------------------------------------------------
# メイン
# ---------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default="elliptic_curves.db")
    parser.add_argument("--out", default="data")
    parser.add_argument("--target-points", type=int, default=100)
    parser.add_argument("--digit-limit", type=int, default=200,
                        help="座標の桁数上限（超えたら線形結合を打ち切り）")
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    conn = sqlite3.connect(args.db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    cur.execute("SELECT MIN(k), MAX(k), MIN(n), MAX(n) FROM curves")
    kmin, kmax, nmin, nmax = cur.fetchone()

    index = {
        "k_range": [int(kmin), int(kmax)],
        "n_range": [int(nmin), int(nmax)],
        "cells": {},
    }

    cur.execute("SELECT * FROM curves ORDER BY k, n")
    rows = cur.fetchall()
    total = len(rows)

    for i, row in enumerate(rows, 1):
        k = row["k"]
        N = row["n"]
        key = f"{k},{N}"

        fname = f"curve_k{k:+05d}_N{N:+05d}.json"
        fpath = os.path.join(args.out, fname)

        if args.skip_existing and os.path.exists(fpath):
            index["cells"][key] = {
                "status": row["status"],
                "rank": row["rank"],
                "has_same_sign_point": False,  # 不明のまま
            }
            print(f"[{i}/{total}] skip {key}")
            continue

        try:
            data = export_curve(row, args.target_points, args.digit_limit)
        except Exception as e:
            print(f"[{i}/{total}] ERROR {key}: {e}")
            data = {
                "k": int(k), "n": int(N),
                "status": "error",
                "note": f"export failed: {e}",
            }

        with open(fpath, "w") as f:
            json.dump(data, f, ensure_ascii=False, separators=(",", ":"))

        index["cells"][key] = {
            "status": row["status"],
            "rank": row["rank"] if row["rank"] is not None else None,
            "has_same_sign_point": data.get("has_same_sign_point", False),
        }
        pc = data.get("points_count", 0)
        print(f"[{i}/{total}] {key} status={row['status']} "
              f"points={pc} same_sign={data.get('has_same_sign_point', False)}")

    with open(os.path.join(args.out, "index.json"), "w") as f:
        json.dump(index, f, ensure_ascii=False, indent=2)

    conn.close()
    print(f"\nDone. Output: {args.out}/")


if __name__ == "__main__":
    main()