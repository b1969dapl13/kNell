# File: scan_kN_yj_cardano.sage
import os, math, cmath
import numpy as np
import matplotlib.pyplot as plt

# ====== 走査範囲 ======
K_MAX = 12     # k = 1..K_MAX
N_MAX = 30    # N = 1..N_MAX
OUTDIR = "scan_kN"
os.makedirs(OUTDIR, exist_ok=True)

TOL_IMAG_B = 1e-7
TOL_IMAG_Y = 1e-6
TOL_DEN    = 1e-8

# ====== 係数 ======
def make_coeff_funcs(k, N):
    k = float(k); N = float(N)
    def A3(a): return -k
    def A2(a): return k**2*N*a - k**2 + k*N - a
    def A1(a): return k**3*N*a - k**2*a**2 + k*N*a**2 + k**2*N - 3*k*a + N*a - 1
    def A0(a): return k**2*N*a**2 - k*a**3 - k**2*a + k*N*a - a**2 - k
    return A3, A2, A1, A0

# ====== 閉形式 b_j(a) ======
def make_branch_funcs(k, N):
    A3f, A2f, A1f, A0f = make_coeff_funcs(k, N)
    omega = cmath.exp(2j*cmath.pi/3.0)

    def make_bj(j):
        wj = omega**j
        def bj(a):
            a = complex(a)
            A3 = complex(A3f(a)); A2 = complex(A2f(a))
            A1 = complex(A1f(a)); A0 = complex(A0f(a))
            p2 = A2/A3; p1 = A1/A3; p0 = A0/A3
            shift = p2/3.0
            p = p1 - p2*p2/3.0
            q = 2.0*p2**3/27.0 - p2*p1/3.0 + p0
            s = cmath.sqrt((q/2.0)**2 + (p/3.0)**3)
            C = (-q/2.0 + s)**(1.0/3.0)
            # 退化: C ~ 0  →  t^3 + p t ≈ 0  →  t ∈ {0, ±sqrt(-p)}
            if abs(C) < 1e-14:
                if abs(p) < 1e-14:
                    return -shift
                r = cmath.sqrt(-p)
                t_candidates = [0.0+0j, r, -r]
                return t_candidates[j] - shift
            return wj*C - p/(3.0*wj*C) - shift
        return bj

    return [make_bj(0), make_bj(1), make_bj(2)]

# ====== y1(a,b) ======
def y1_ab(a, b, k, N):
    a = complex(a); b = complex(b); k = float(k); N = float(N)
    num = ((-4*k**5 - 4*k**4*N - 8*k**2 - 4*k*N)*a**2
           + (-4*k**6 - 4*k**5*N - 8*k**3 - 4*k**2*N + 4)*a*b
           + (-4*k**3*N - 4*N)*a
           + (4*k)*b**2
           + (-4*k**5 - 4*k**4*N - 4*k**2 - 4*k*N)*b
           + (4*k**3 + 4))
    den = (k*b + a)**2
    if den == 0:
        return complex('nan')
    return num/den

# ====== 漸近係数 c_j (a -> ∞) ======
def make_c_branches(k, N):
    k = float(k); N = float(N)
    A3 = -k; A2 = N*k**2 - 1; A1 = N*k - k**2; A0 = -k
    p2 = A2/A3; p1 = A1/A3; p0 = A0/A3
    shift = p2/3.0
    p = p1 - p2*p2/3.0
    q = 2.0*p2**3/27.0 - p2*p1/3.0 + p0
    s = cmath.sqrt(complex((q/2.0)**2 + (p/3.0)**3))
    C = (complex(-q/2.0) + s)**(1.0/3.0)
    omega = cmath.exp(2j*cmath.pi/3.0)

    # 退化ケース
    if abs(C) < 1e-14:
        if abs(p) < 1e-14:
            ts = [0.0+0j, 0.0+0j, 0.0+0j]
        else:
            r = cmath.sqrt(-complex(p))
            ts = [0.0+0j, r, -r]
        return [t - shift for t in ts]

    return [omega**j * C - p/(3.0*omega**j * C) - shift for j in range(3)]

def y1_limit_safe(c, k, N):
    k = float(k); N = float(N); c = complex(c)
    den = (k*c + 1)**2
    if abs(den) < 1e-12:
        return None
    num = ((-4*k**5 - 4*k**4*N - 8*k**2 - 4*k*N)
           + (-4*k**6 - 4*k**5*N - 8*k**3 - 4*k**2*N + 4)*c
           + 4*k*c**2)
    return num/den

# ====== 端点分類 ======
def classify_extreme(a_star, a_arr_pos, A_MIN, A_MAX):
    a_lo, a_hi = float(a_arr_pos.min()), float(a_arr_pos.max())
    if a_star <= A_MIN*2:        return "a->0+ (infinitesimal)"
    if a_star >= A_MAX/2:        return "a->inf (asymptotic)"
    if abs(math.log10(a_star) - math.log10(a_lo)) < 0.05:
        return f"left edge of {{b>0}} a~{a_lo:.3g}"
    if abs(math.log10(a_star) - math.log10(a_hi)) < 0.05:
        return f"right edge of {{b>0}} a~{a_hi:.3g}"
    return "interior (finite)"

# ====== 1 つの (k,N) を処理 ======
def process_pair(k_val, N_val, a_samples):
    A_MIN, A_MAX = float(a_samples.min()), float(a_samples.max())
    b_funcs = make_branch_funcs(k_val, N_val)
    c_roots = make_c_branches(k_val, N_val)

    # 実数値のみ残す（複素解は捨てる）
    real_data = {j: {'a': [], 'b': []} for j in range(3)}
    for a in a_samples:
        for j in range(3):
            bj = b_funcs[j](a)
            if abs(bj.imag) < TOL_IMAG_B:
                real_data[j]['a'].append(float(a))
                real_data[j]['b'].append(float(bj.real))

    # a>0 で常に負のブランチを検出（y1 は描かない）
    always_neg = []
    for j in range(3):
        b_arr = np.array(real_data[j]['b'])
        if b_arr.size >= max(50, 0.5*len(a_samples)) and np.all(b_arr < 0):
            always_neg.append(j)

    # y1 計算
    valid = {j: {'a': [], 'b': [], 'y1': []} for j in range(3)}
    for j in range(3):
        if j in always_neg:
            continue
        for a, b in zip(real_data[j]['a'], real_data[j]['b']):
            den = k_val*b + a
            if abs(den) < TOL_DEN:
                continue
            y1v = y1_ab(a, b, k_val, N_val)
            if np.isfinite(y1v.real) and abs(y1v.imag) < TOL_IMAG_Y:
                valid[j]['a'].append(a)
                valid[j]['b'].append(b)
                valid[j]['y1'].append(y1v.real)

    # ---- プロット ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = ['tab:blue', 'tab:orange', 'tab:green']

    # (左) b_j(a)
    axL = axes[0]
    for j in range(3):
        if not real_data[j]['a']:
            continue
        tag = " [always<0: y1 skipped]" if j in always_neg else ""
        axL.plot(real_data[j]['a'], real_data[j]['b'],
                 color=colors[j], lw=1.4, label=f'b_{j}(a){tag}')
    axL.set_xscale('log'); axL.set_yscale('symlog', linthresh=1)
    axL.axhline(0, color='gray', lw=0.5)
    axL.set_xlabel('a'); axL.set_ylabel('b_j(a)  (real part only)')
    axL.set_title(f'k={k_val}, N={N_val}: real branches of b_j(a)')
    axL.grid(True, which='both', alpha=0.3)
    axL.legend(fontsize=8)

    # (右) y1
    axR = axes[1]
    band_info = []

    strip_y_positions = {0: 0.965, 1: 0.945, 2: 0.925}
    strip_thickness   = 0.012

    def contiguous_runs(mask, x):
        runs = []
        in_run = False; x0 = None
        for i, m in enumerate(mask):
            if m and not in_run:
                in_run = True; x0 = x[i]
            elif not m and in_run:
                in_run = False
                runs.append((x0, x[i-1]))
        if in_run:
            runs.append((x0, x[-1]))
        return runs

    bneg_legend_added = {0: False, 1: False, 2: False}
    for j in range(3):
        if j in always_neg or not real_data[j]['a']:
            continue
        a_full = np.array(real_data[j]['a'])
        b_full = np.array(real_data[j]['b'])
        order = np.argsort(a_full)
        a_full = a_full[order]; b_full = b_full[order]

        neg_runs = contiguous_runs(b_full < 0, a_full)
        for (xl, xr) in neg_runs:
            if xl == xr:
                continue
            label = None
            if not bneg_legend_added[j]:
                label = f'b_{j}<0 region'
                bneg_legend_added[j] = True
            axR.axvspan(xl, xr, color=colors[j], alpha=0.08,
                        hatch='//', edgecolor=colors[j], linewidth=0.0,
                        label=label)
            axR.fill_between(
                [xl, xr],
                strip_y_positions[j] - strip_thickness/2,
                strip_y_positions[j] + strip_thickness/2,
                transform=axR.get_xaxis_transform(),
                color=colors[j], alpha=0.85, linewidth=0)

    for j in range(3):
        if j in always_neg or not valid[j]['a']:
            continue
        a_arr  = np.array(valid[j]['a'])
        b_arr  = np.array(valid[j]['b'])
        y1_arr = np.array(valid[j]['y1'])

        axR.plot(a_arr, y1_arr, color=colors[j], lw=1.4,
                 label=f'y1 on j={j}')

        lim = y1_limit_safe(c_roots[j], k_val, N_val)
        if lim is not None and abs(lim.imag) < 1e-6:
            axR.axhline(lim.real, color=colors[j], ls='--', lw=0.9,
                        label=f'lim j={j}: {lim.real:.3g}')

        mask_pos = b_arr > 0
        if np.any(mask_pos):
            a_pos  = a_arr[mask_pos]
            y1_pos = y1_arr[mask_pos]
            imin = int(np.argmin(y1_pos)); imax = int(np.argmax(y1_pos))
            ymin, ymax = float(y1_pos[imin]), float(y1_pos[imax])
            a_at_min, a_at_max = float(a_pos[imin]), float(a_pos[imax])
            cls_min = classify_extreme(a_at_min, a_pos, A_MIN, A_MAX)
            cls_max = classify_extreme(a_at_max, a_pos, A_MIN, A_MAX)
            band_info.append(dict(
                j=j, ymin=ymin, ymax=ymax,
                a_min=a_at_min, a_max=a_at_max,
                cls_min=cls_min, cls_max=cls_max,
                lim=(lim.real if lim is not None else None)))

            axR.axhspan(ymin, ymax, color=colors[j], alpha=0.10)
            axR.axvline(a_at_min, color=colors[j], ls=':',  lw=0.9)
            axR.axvline(a_at_max, color=colors[j], ls='-.', lw=0.9)
            axR.plot([a_at_min],[ymin], marker='v', ms=8,
                     color=colors[j], mec='k', mew=0.6,
                     label=f'min={ymin:.3g}@a={a_at_min:.2g}')
            axR.plot([a_at_max],[ymax], marker='^', ms=8,
                     color=colors[j], mec='k', mew=0.6,
                     label=f'max={ymax:.3g}@a={a_at_max:.2g}')

    axR.set_xscale('log'); axR.set_yscale('symlog', linthresh=10)
    axR.axhline(0, color='gray', lw=0.5)
    axR.set_xlabel('a'); axR.set_ylabel('y1  (symlog)')
    axR.set_title(f'k={k_val}, N={N_val}: y1 on branches\n'
                  '(shaded+hatched: a-range where b_j(a)<0; top strips: per-j)')
    axR.grid(True, which='both', alpha=0.3)
    axR.legend(fontsize=6, loc='best')

    plt.suptitle(
        f'k={k_val}, N={N_val}  | always-negative branches j={always_neg}',
        fontsize=11, y=1.02)
    plt.tight_layout()
    fname = os.path.join(OUTDIR, f'kN_k{k_val:02d}_N{N_val:02d}.png')
    plt.savefig(fname, dpi=120, bbox_inches='tight')
    plt.close(fig)

    return dict(k=k_val, N=N_val,
                always_neg=always_neg,
                band_info=band_info,
                fname=fname)

# ====== a サンプリング ======
a_samples = np.unique(np.concatenate([
    np.geomspace(1e-4, 1e8, 2500),
    np.linspace(1e-4, 5, 800),
    np.linspace(1e-4, 0.5, 800),
]))
a_samples.sort()

# ====== 走査 ======
summary_rows = []
failed_pairs = []
print(f"# scan k=1..{K_MAX}, N=1..{N_MAX}")
print(f"# saving plots into ./{OUTDIR}/\n")

for k_val in range(1, K_MAX+1):
    for N_val in range(1, N_MAX+1):
        try:
            res = process_pair(k_val, N_val, a_samples)
        except Exception as e:
            print(f"!! skip k={k_val}, N={N_val}: {type(e).__name__}: {e}")
            failed_pairs.append((k_val, N_val, repr(e)))
            continue

        line = f"k={k_val:2d} N={N_val:2d} | always_neg j={res['always_neg']}"
        for bi in res['band_info']:
            line += (f"\n    j={bi['j']}: y1 in [{bi['ymin']:.4g}, {bi['ymax']:.4g}]"
                     f"  argmin@a={bi['a_min']:.3g} ({bi['cls_min']})"
                     f"  argmax@a={bi['a_max']:.3g} ({bi['cls_max']})"
                     f"  lim={bi['lim']}")
        print(line)
        summary_rows.append(res)

# ====== CSV サマリー ======
csv_path = os.path.join(OUTDIR, 'summary.csv')
with open(csv_path, 'w') as f:
    f.write("k,N,branch_j,always_neg_js,y1_min,a_at_min,cls_min,"
            "y1_max,a_at_max,cls_max,y1_limit\n")
    for res in summary_rows:
        an = "|".join(map(str, res['always_neg']))
        if not res['band_info']:
            f.write(f"{res['k']},{res['N']},,{an},,,,,,,\n")
            continue
        for bi in res['band_info']:
            f.write(f"{res['k']},{res['N']},{bi['j']},{an},"
                    f"{bi['ymin']},{bi['a_min']},{bi['cls_min']},"
                    f"{bi['ymax']},{bi['a_max']},{bi['cls_max']},"
                    f"{bi['lim']}\n")
print(f"\nSaved: {csv_path}")

if failed_pairs:
    print("\n# failed pairs:")
    for kk, NN, msg in failed_pairs:
        print(f"  k={kk}, N={NN}: {msg}")