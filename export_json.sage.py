# File: export_json.sage.py
"""
elliptic_curves.db から個別曲線JSONとindex.jsonを生成する (v3)。

v3 変更点:
- k>0, N>0 のとき、a:b:c が全て正となる有理点を見つけるため adaptive 探索を実装
- 卵成分（E(R) の有界連結成分）に乗る点が生成可能なら、M を拡張して粘る
- 卵成分判定は短形式 Y^2 = X^3 + AX + B に変換後、Q 演算のみで厳密に判定

使い方:
    sage -python export_json.sage.py --db elliptic_curves.db --out data/
    sage -python export_json.sage.py --db elliptic_curves.db --out data/ --target-points 100
    sage -python export_json.sage.py --db elliptic_curves.db --out data/ --m-extended-max 15
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
from sage.all import EllipticCurve, QQ, Integer, PolynomialRing


# ---------------------------------------------------------------
# 逆変換: (y1, y2) -> (a:b:c)
# ---------------------------------------------------------------
def inverse_transform(y1, y2, k, N):
    """最終ワイエルシュトラス形の点 (y1, y2) から元の射影座標 (a:b:c) を復元。"""
    k = Fraction(k)
    N = Fraction(N)

    x1_star = 2 * k * (k**3 + 1) * N + 4 * k**2
    x1 = y1 + x1_star
    x2 = y2

    denom6 = k**3 + 1
    if denom6 == 0:
        return None
    w1 = -x1 / (2 * denom6)
    w2 = -x2 / denom6

    beta2 = w1**2 - k**2 * N**2 + 4 * k
    beta1 = (-2 * (k**3 + 1) * w1 * N - 2 * k * w1
             - 2 * (k**4 + k) * N**2 - 2 * k**2 * N + 4 * k**3 + 4)
    if beta2 == 0:
        return None
    v1 = (w2 - beta1) / (2 * beta2)

    L_v1 = (k**3 + 1) - ((k**3 + 1) * N + k) * v1
    v2 = L_v1 + v1**2 * w1

    u2 = v1
    alpha2 = -k**2 * u2 * N + (-k**3 * u2 + k**2 - u2)
    alpha1 = ((k * u2**2 - k**3 * u2 + u2) * N
              + (2 * k**2 * u2**2 + k**3 - k * u2 - 1))
    if alpha2 == 0:
        return None
    u1 = (v2 - alpha1) / (2 * alpha2)

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
    if abc is None:
        return None
    return [sign_of(v) for v in abc]


def is_all_same_sign(abc):
    if abc is None:
        return False
    s = sign_pattern(abc)
    if 0 in s:
        return False
    return s[0] == s[1] == s[2]


def is_all_positive(abc):
    """a, b, c 全て正"""
    if abc is None:
        return False
    return all(v > 0 for v in abc)


# ---------------------------------------------------------------
# 卵成分判定（厳密、Q演算のみ）
# ---------------------------------------------------------------
def _to_short_weierstrass(E):
    """E を短形式 Y^2 = X^3 + AX + B に変換し、(E_short, phi, A, B, Delta) を返す。
    phi は E -> E_short の射。失敗時 None。"""
    try:
        E_short = E.short_weierstrass_model()
        phi = E.isomorphism_to(E_short)
        A = QQ(E_short.a4())
        B = QQ(E_short.a6())
        Delta = -16 * (4 * A**3 + 27 * B**2)
        return E_short, phi, A, B, Delta
    except Exception:
        return None


def has_egg_component(E):
    """E(R) が卵成分（有界連結成分）を持つか。"""
    info = _to_short_weierstrass(E)
    if info is None:
        return False
    _, _, _, _, Delta = info
    return Delta > 0


def is_on_egg_component(P, short_info):
    """点 P が E(R) の卵成分上にあるか厳密判定。
    short_info: _to_short_weierstrass(E) の返り値。"""
    if short_info is None:
        return False
    E_short, phi, A, B, Delta = short_info
    if Delta <= 0:
        return False
    if P.is_zero():
        return False
    try:
        P_short = phi(P)
    except Exception:
        return False
    if P_short.is_zero():
        return False
    try:
        xP = QQ(P_short.xy()[0])
    except Exception:
        return False
    # x_P < sqrt(-A/3) ⟺ x_P ≤ 0 または 3·x_P² + A ≤ 0
    if xP <= 0:
        return True
    if 3 * xP**2 + A <= 0:
        return True
    return False


def can_reach_egg(E, free_gens, torsion_pts):
    """生成元の整数結合（+ねじれ）で卵成分上の点が得られる可能性があるか。
    - 単独 free generator が卵にある
    - 単独 torsion が卵にある
    - free + torsion で卵にある
    のいずれか。"""
    short_info = _to_short_weierstrass(E)
    if short_info is None:
        return False, None
    _, _, _, _, Delta = short_info
    if Delta <= 0:
        return False, short_info

    for g in free_gens:
        if is_on_egg_component(g, short_info):
            return True, short_info
    for _, T in torsion_pts:
        if is_on_egg_component(T, short_info):
            return True, short_info
    # free + torsion の単純な組み合わせ
    for g in free_gens:
        for _, T in torsion_pts:
            try:
                Q = g + T
                if is_on_egg_component(Q, short_info):
                    return True, short_info
            except Exception:
                continue
    # 卵成分があるなら、倍々していけば理論上は到達可能（連結成分は群の指数2の部分群）
    # ただし生成元が全て無限成分かつ卵2-torsionがない場合は到達不能
    return False, short_info


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
    seen = {}
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
    if torsion_coeffs is None:
        return False
    nonzero = [(i, c) for i, c in enumerate(torsion_coeffs) if c != 0]
    return len(nonzero) == 1 and nonzero[0][1] == 1


# ---------------------------------------------------------------
# 点ロール判定
# ---------------------------------------------------------------
def classify_point(free_coeffs, torsion_coeffs):
    free_nonzero = [c for c in free_coeffs if c != 0] if free_coeffs else []
    has_torsion = torsion_coeffs is not None and any(c != 0 for c in torsion_coeffs)

    if not free_nonzero and not has_torsion:
        return "combination"

    if not free_nonzero and has_torsion:
        if is_torsion_generator(torsion_coeffs):
            return "torsion_generator"
        return "torsion_other"

    if free_nonzero and not has_torsion:
        if len(free_nonzero) == 1 and free_nonzero[0] == 1:
            return "free_generator"
        return "combination"

    return "mixed"


def point_label(free_coeffs, torsion_coeffs):
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
    if s.startswith("+"):
        s = s[1:].lstrip()
    return s


# ---------------------------------------------------------------
# y₁範囲の計算（plot_yj.txtの手法を移植）
# ---------------------------------------------------------------
def compute_y1_ranges(k, N, a_samples=400, a_max=1e8):
    k_f = float(k)
    N_f = float(N)

    def A3(a): return -k_f
    def A2(a): return k_f**2*N_f*a - k_f**2 + k_f*N_f - a
    def A1(a): return k_f**3*N_f*a - k_f**2*a**2 + k_f*N_f*a**2 + k_f**2*N_f - 3*k_f*a + N_f*a - 1
    def A0(a): return k_f**2*N_f*a**2 - k_f*a**3 - k_f**2*a + k_f*N_f*a - a**2 - k_f

    omega = cmath.exp(2j*cmath.pi/3.0)

    def bj_func(a, j):
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

    try:
        a_vals = np.unique(np.concatenate([
            np.geomspace(1e-6, a_max, a_samples),
            np.linspace(1e-6, 1000, a_samples),
            np.linspace(1e-6, 10, a_samples // 2)
        ]))
        a_vals = a_vals[a_vals > 0]
        a_vals = list(a_vals)
    except Exception as e_np:
        print(f"    numpy unavailable for k={k}, N={N}: {e_np}, using fallback sampling")
        a_vals = []
        for i in range(a_samples):
            a_vals.append(1e-6 * (10 ** (i * 14.0 / a_samples)))
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
# 点生成ヘルパー（指定 M シェルでの点列挙）
# ---------------------------------------------------------------
def _build_point_record(P, free_coeffs, torsion_coeffs, k, N, digit_limit,
                       short_info=None, include_order=False):
    """P から1点分のレコードを作る。座標桁が制限を超えたら None。"""
    if P.is_zero():
        return None
    try:
        x, y = P.xy()
    except Exception:
        return None

    x_f = Fraction(int(x.numerator()), int(x.denominator()))
    y_f = Fraction(int(y.numerator()), int(y.denominator()))

    if max(coord_digits(x_f), coord_digits(y_f)) > digit_limit:
        return None

    abc = inverse_transform(x_f, y_f, k, N)
    role = classify_point(free_coeffs, torsion_coeffs)
    label = point_label(free_coeffs, torsion_coeffs)
    sp = sign_pattern(abc)
    rec = {
        "free_coeffs": list(free_coeffs) if free_coeffs else [],
        "torsion_coeffs": torsion_coeffs,
        "role": role,
        "label": label,
        "xy": [str(x), str(y)],
        "abc": [str(v) for v in abc] if abc else None,
        "sign_pattern": sp,
        "all_same_sign": is_all_same_sign(abc),
        "all_positive": is_all_positive(abc),
    }
    if short_info is not None:
        rec["on_egg"] = is_on_egg_component(P, short_info)
    if include_order:
        try:
            rec["order"] = int(P.order())
        except Exception:
            pass
    return rec


def _enumerate_shell(M_outer, M_inner, free_gens, torsion_variants,
                    E, k, N, digit_limit, short_info):
    """係数の max(|c_i|) == M_outer のシェル（M_inner より外側の殻）の点を列挙。
    M_inner は「既に列挙済みの最大 M」（None なら最初から全部）。"""
    rank = len(free_gens)
    if rank == 0:
        return []

    coeff_ranges = [range(-M_outer, M_outer + 1)] * rank
    results = []
    for coeffs in product(*coeff_ranges):
        if all(c == 0 for c in coeffs):
            continue
        # シェル判定: 少なくとも1成分が ±M_outer
        max_abs = max(abs(c) for c in coeffs)
        if M_inner is not None and max_abs <= M_inner:
            continue
        if max_abs != M_outer:
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
                try:
                    Q = P + T
                except Exception:
                    continue
                tc_record = torsion_coeffs

            if Q.is_zero():
                continue

            rec = _build_point_record(Q, list(coeffs), tc_record, k, N,
                                     digit_limit, short_info=short_info)
            if rec is not None:
                results.append(rec)
    return results


# ---------------------------------------------------------------
# 1曲線のエクスポート
# ---------------------------------------------------------------
def export_curve(row, target_points, digit_limit, m_extended_max=15):
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

    # 短形式情報（卵成分判定用）
    short_info = _to_short_weierstrass(E)
    if short_info is not None:
        _, _, A_short, B_short, Delta_short = short_info
        out["short_weierstrass"] = {
            "A": str(A_short),
            "B": str(B_short),
            "discriminant": str(Delta_short),
            "has_egg_component": Delta_short > 0,
        }
    else:
        out["short_weierstrass"] = None

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
            "on_egg": is_on_egg_component(g, short_info),
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

    orders = [int(g.order()) for g in torsion_gens]

    out["torsion_generators"] = [
        {
            "index": j + 1,
            "label": f"T_{j+1}",
            "order": int(orders[j]),
            "xy": [str(tg.xy()[0]), str(tg.xy()[1])],
            "on_egg": is_on_egg_component(tg, short_info),
        }
        for j, tg in enumerate(torsion_gens)
    ]

    # y1 帯
    try:
        out["y1_bands"] = compute_y1_bands(k, N)
    except Exception as e:
        out["y1_bands"] = []
        out["y1_bands_note"] = f"compute failed: {e}"

    # 全ねじれ点
    torsion_pts = enumerate_torsion_points(E, torsion_gens, orders)

    # 卵成分到達可能性
    egg_reachable, _ = can_reach_egg(E, free_gens, torsion_pts)
    out["egg_reachable"] = egg_reachable

    # 自由部分の係数を決定
    M = determine_M(rank, max(row["torsion_order"] or 1, 1),
                    target=target_points)
    out["M"] = M

    points = []

    # まず純粋ねじれ点
    for torsion_coeffs, T in torsion_pts:
        rec = _build_point_record(T, [0] * rank, torsion_coeffs, k, N,
                                 digit_limit, short_info=short_info,
                                 include_order=True)
        if rec is not None:
            points.append(rec)

    # 自由部分 + (ねじれを足したバリアント)
    torsion_variants = [(None, None)] + [(tc, T) for tc, T in torsion_pts]

    if rank > 0:
        # M=1, 2, ..., M まで段階的にシェルを足す
        for M_curr in range(1, M + 1):
            shell_pts = _enumerate_shell(
                M_curr, M_curr - 1 if M_curr > 0 else None,
                free_gens, torsion_variants, E, k, N, digit_limit, short_info
            )
            points.extend(shell_pts)

    # ----------------------------------------------------------
    # adaptive 拡張: k>0, N>0 で正符号点が未発見 かつ 卵到達可能なら粘る
    # ----------------------------------------------------------
    search_status = "initial"
    found_positive_at_M = None

    # 初期 M で正符号点が見つかったか
    if any(p["all_positive"] for p in points):
        search_status = "found_in_initial"
        # どの M で見つかったかを記録
        for p in points:
            if p["all_positive"]:
                fc = p.get("free_coeffs") or []
                m_used = max((abs(c) for c in fc), default=0)
                if found_positive_at_M is None or m_used < found_positive_at_M:
                    found_positive_at_M = m_used

    should_extend = (
        k > 0 and N > 0
        and rank > 0
        and not any(p["all_positive"] for p in points)
        and egg_reachable
        and M < m_extended_max
    )

    if should_extend:
        search_status = "extending"
        M_final = M
        for M_curr in range(M + 1, m_extended_max + 1):
            shell_pts = _enumerate_shell(
                M_curr, M_curr - 1,
                free_gens, torsion_variants, E, k, N, digit_limit, short_info
            )
            points.extend(shell_pts)
            M_final = M_curr
            if any(p["all_positive"] for p in shell_pts):
                found_positive_at_M = M_curr
                search_status = "found_in_extension"
                break
        else:
            search_status = "extension_exhausted"
        out["M_extended"] = M_final
    elif k > 0 and N > 0 and not any(p["all_positive"] for p in points):
        if not egg_reachable:
            search_status = "no_egg_access"
        elif rank == 0:
            search_status = "rank_zero"
        else:
            search_status = "skipped"

    out["search_status"] = search_status
    out["found_positive_at_M"] = found_positive_at_M

    out["points"] = points
    out["points_count"] = len(points)
    out["has_same_sign_point"] = any(p["all_same_sign"] for p in points)
    out["has_positive_point"] = any(p["all_positive"] for p in points)

    # y₁範囲の数値計算
    try:
        y1_ranges = compute_y1_ranges(k, N)
        out["y1_ranges"] = y1_ranges
    except Exception as e:
        out["y1_ranges"] = []
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
    parser.add_argument("--m-extended-max", type=int, default=15,
                        help="adaptive 拡張時の M 上限 (k>0,N>0 かつ卵到達可能時)")
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
                "has_same_sign_point": False,
                "has_positive_point": False,
            }
            print(f"[{i}/{total}] skip {key}")
            continue

        try:
            data = export_curve(row, args.target_points, args.digit_limit,
                              m_extended_max=args.m_extended_max)
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
            "has_positive_point": data.get("has_positive_point", False),
            "egg_reachable": data.get("egg_reachable", False),
            "search_status": data.get("search_status", None),
        }
        pc = data.get("points_count", 0)
        ss = data.get("search_status", "-")
        print(f"[{i}/{total}] {key} status={row['status']} rank={row['rank']} "
              f"points={pc} pos={data.get('has_positive_point', False)} "
              f"egg={data.get('egg_reachable', False)} search={ss}")

    with open(os.path.join(args.out, "index.json"), "w") as f:
        json.dump(index, f, ensure_ascii=False, indent=2)

    conn.close()
    print(f"\nDone. Output: {args.out}/")


if __name__ == "__main__":
    main()