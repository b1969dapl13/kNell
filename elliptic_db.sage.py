# File: elliptic_db.py
# SageMath script.
#   sage elliptic_db.py
# or in Sage REPL:
#   sage: load("elliptic_db.py")
#   sage: run(k_range=range(-5,6), n_range=range(-10,11), timeout_sec=60)
#   sage: summary()
#   sage: query(3, -4)
#   sage: rerun_undetermined(timeout_sec=300)

import sqlite3
import json
import time
import traceback
import warnings
import multiprocessing as mp

from sage.all import EllipticCurve, ZZ, QQ

# Sage 内部の古い deprecation 警告を抑制
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    module="sage.schemes.elliptic_curves",
)

DB_PATH = "elliptic_curves.db"


# ---------------------------------------------------------------
# DB 初期化
# ---------------------------------------------------------------
def init_db(path=DB_PATH):
    conn = sqlite3.connect(path)
    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE IF NOT EXISTS curves (
        k              INTEGER NOT NULL,
        n              INTEGER NOT NULL,
        status         TEXT    NOT NULL,   -- 'ok' / 'singular' / 'timeout' / 'crashed' / 'error'
        a_invariants   TEXT,               -- JSON list [a1,a2,a3,a4,a6]
        discriminant   TEXT,
        j_invariant    TEXT,
        conductor      TEXT,
        rank           INTEGER,
        rank_bounds    TEXT,               -- JSON "[low, high]" when not determined
        gens           TEXT,               -- JSON list of [x, y]
        torsion_order  INTEGER,
        torsion_struct TEXT,
        torsion_pts    TEXT,               -- JSON list of [x, y]
        timeout_sec    REAL,
        elapsed_sec    REAL,
        note           TEXT,
        PRIMARY KEY (k, n)
    )
    """)
    conn.commit()
    return conn


def already_done(conn, k, N, retry_timeout=False, retry_error=False):
    """既に十分な結果がある？ status に応じて再挑戦の可否を判定。"""
    cur = conn.cursor()
    cur.execute("SELECT status FROM curves WHERE k=? AND n=?", (k, N))
    row = cur.fetchone()
    if row is None:
        return False
    status = row[0]
    if status in ("ok", "singular"):
        return True
    if status == "timeout" and not retry_timeout:
        return True
    if status in ("error", "crashed") and not retry_error:
        return True
    return False


def save_result(conn, rec):
    cur = conn.cursor()
    cur.execute("""
    INSERT OR REPLACE INTO curves
    (k, n, status, a_invariants, discriminant, j_invariant, conductor,
     rank, rank_bounds, gens, torsion_order, torsion_struct, torsion_pts,
     timeout_sec, elapsed_sec, note)
    VALUES (:k, :n, :status, :a_invariants, :discriminant, :j_invariant, :conductor,
            :rank, :rank_bounds, :gens, :torsion_order, :torsion_struct, :torsion_pts,
            :timeout_sec, :elapsed_sec, :note)
    """, rec)
    conn.commit()


# ---------------------------------------------------------------
# 曲線構成
#   y^2 = x^3 + alpha^2 * x^2 + 8 alpha beta * x + 16 beta^2
#   alpha = (k^3+1) N + 3 k
#   beta  = k^2 (k^3+1) N + (k^6 + 3 k^3 + 1)
# ---------------------------------------------------------------
def build_curve(k, N):
    k = ZZ(k); N = ZZ(N)
    alpha = (k**3 + 1) * N + 3 * k
    beta  = k**2 * (k**3 + 1) * N + (k**6 + 3 * k**3 + 1)
    a2 = alpha**2
    a4 = 8 * alpha * beta
    a6 = 16 * beta**2
    return alpha, beta, [0, a2, 0, a4, a6]   # [a1, a2, a3, a4, a6]


def short_weierstrass_discriminant(a2, a4, a6):
    """y^2 = x^3 + a2 x^2 + a4 x + a6 の判別式 (×16 を省いた素の式)."""
    return (-4 * a2**3 * a6
            + a2**2 * a4**2
            + 18 * a2 * a4 * a6
            - 4 * a4**3
            - 27 * a6**2)


# ---------------------------------------------------------------
# サブプロセスで重い計算
#   - rank / gens / conductor
#   - 子プロセスが SIGSEGV しても親プロセスは生存
# ---------------------------------------------------------------
def _heavy_worker(ainvs_int, pipe):
    try:
        import warnings as _w
        _w.filterwarnings("ignore", category=DeprecationWarning)

        from sage.all import EllipticCurve, QQ
        E = EllipticCurve(QQ, ainvs_int)
        out = {}

        try:
            out["conductor"] = str(E.conductor())
        except Exception as e:
            out["conductor_err"] = f"{type(e).__name__}: {e}"

        # 階数 (PARI)
        rank_known = False
        try:
            r = E.rank(algorithm="pari")
            out["rank"] = int(r)
            rank_known = True
        except Exception as e1:
            try:
                lo, hi = E.rank_bounds()
                out["rank_bounds"] = [int(lo), int(hi)]
                if lo == hi:
                    out["rank"] = int(lo)
                    rank_known = True
            except Exception as e2:
                out["rank_err"] = f"rank: {e1}; bounds: {e2}"

        # 生成元
        try:
            if rank_known:
                gens = E.gens(proof=False)
            else:
                gens = E.gens()
            out["gens"] = [[str(P[0]), str(P[1])] for P in gens]
        except Exception as e:
            out["gens_err"] = f"{type(e).__name__}: {e}"

        pipe.send(("ok", out))
    except Exception as e:
        pipe.send(("error",
                   f"{type(e).__name__}: {e}\n{traceback.format_exc()}"))
    finally:
        pipe.close()


def heavy_with_timeout(ainvs_int, timeout_sec):
    """サブプロセスで重い計算。タイムアウト / クラッシュ耐性付き。"""
    parent_conn, child_conn = mp.Pipe(duplex=False)
    proc = mp.Process(target=_heavy_worker, args=(ainvs_int, child_conn))
    proc.start()
    child_conn.close()
    proc.join(timeout_sec)

    if proc.is_alive():
        proc.terminate()
        proc.join(2.0)
        if proc.is_alive():
            proc.kill()
            proc.join()
        return "timeout", None

    if parent_conn.poll():
        try:
            return parent_conn.recv()
        except EOFError:
            return "crashed", f"exitcode={proc.exitcode}"
    return "crashed", f"exitcode={proc.exitcode}"


# ---------------------------------------------------------------
# 1 ケースの解析
# ---------------------------------------------------------------
def analyze(k, N, timeout_sec=60.0, verbose=True):
    t0 = time.time()
    rec = {
        "k": int(k), "n": int(N),
        "status": "error",
        "a_invariants": None, "discriminant": None, "j_invariant": None,
        "conductor": None,
        "rank": None, "rank_bounds": None, "gens": None,
        "torsion_order": None, "torsion_struct": None, "torsion_pts": None,
        "timeout_sec": timeout_sec, "elapsed_sec": None, "note": None,
    }
    try:
        alpha, beta, ainvs = build_curve(k, N)
        a2_, a4_, a6_ = ainvs[1], ainvs[3], ainvs[4]

        # alpha=0 / beta=0 は自明解が退化しているので singular 扱い
        if alpha == 0 or beta == 0:
            rec["status"] = "singular"
            rec["note"] = f"alpha={alpha}, beta={beta}"
            rec["a_invariants"] = json.dumps([int(a) for a in ainvs])
            rec["elapsed_sec"] = time.time() - t0
            return rec

        # EllipticCurve を作る前に判別式を先回りチェック
        disc_short = short_weierstrass_discriminant(a2_, a4_, a6_)
        if disc_short == 0:
            rec["status"] = "singular"
            rec["a_invariants"] = json.dumps([int(a) for a in ainvs])
            rec["discriminant"] = "0"
            rec["note"] = "discriminant=0 (precomputed)"
            rec["elapsed_sec"] = time.time() - t0
            return rec

        try:
            E = EllipticCurve(QQ, ainvs)
        except ArithmeticError as e:
            rec["status"] = "singular"
            rec["a_invariants"] = json.dumps([int(a) for a in ainvs])
            rec["note"] = f"singular: {e}"
            rec["elapsed_sec"] = time.time() - t0
            return rec

        rec["a_invariants"] = json.dumps([int(a) for a in ainvs])
        rec["discriminant"] = str(E.discriminant())
        rec["j_invariant"]  = str(E.j_invariant())

        # ねじれ群（軽い）
        try:
            T = E.torsion_subgroup()
            rec["torsion_order"]  = int(T.order())
            rec["torsion_struct"] = str(T.invariants())
            rec["torsion_pts"]    = json.dumps(
                [[str(P[0]), str(P[1])] for P in T.points() if not P.is_zero()]
            )
        except Exception as e:
            rec["note"] = f"torsion: {type(e).__name__}: {e}"

        # 重い計算は子プロセス
        ainvs_int = [int(a) for a in ainvs]
        tag, data = heavy_with_timeout(ainvs_int, timeout_sec)

        if tag == "ok":
            rec["status"] = "ok"
            rec["conductor"] = data.get("conductor")
            if "rank" in data:
                rec["rank"] = data["rank"]
            if "rank_bounds" in data:
                rec["rank_bounds"] = json.dumps(data["rank_bounds"])
            if "gens" in data:
                rec["gens"] = json.dumps(data["gens"])
            notes = []
            for key in ("conductor_err", "rank_err", "gens_err"):
                if key in data:
                    notes.append(f"{key}={data[key]}")
            if notes:
                prev = rec.get("note")
                joined = "; ".join(notes)
                rec["note"] = joined if not prev else f"{prev}; {joined}"
        elif tag == "timeout":
            rec["status"] = "timeout"
            rec["note"] = f"heavy computation timed out at {timeout_sec}s"
        elif tag == "crashed":
            rec["status"] = "crashed"
            rec["note"] = f"subprocess crashed: {data}"
        else:
            rec["status"] = "error"
            rec["note"] = str(data)

    except Exception as e:
        rec["status"] = "error"
        rec["note"] = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"

    rec["elapsed_sec"] = time.time() - t0

    if verbose:
        ng = "-"
        if rec["gens"]:
            try:
                ng = str(len(json.loads(rec["gens"])))
            except Exception:
                ng = "?"
        msg = (f"  [{rec['status']}] k={k}, N={N}, "
               f"rank={rec['rank']}, bounds={rec['rank_bounds']}, "
               f"tors={rec['torsion_struct']}, ngens={ng}, "
               f"t={rec['elapsed_sec']:.2f}s")
        if rec["note"]:
            msg += f"  // {rec['note'][:120]}"
        print(msg)
    return rec


# ---------------------------------------------------------------
# (k, N) の列挙と優先順位
# ---------------------------------------------------------------
def enumerate_pairs(k_range, n_range, exclude_k=(-1)):  # exclude_k=(-1, 0, 1)
    """
    優先度 = 2|k| + |N| (小さいほど先).
    k = -1, 0, 1 では k^3+1 == 0 や k == 0 で曲線が退化する組み合わせが多い。
    """
    pairs = []
    for k in k_range:
        if k in exclude_k:
            continue
        for N in n_range:
            priority = abs(k) * 2 + abs(N)
            pairs.append((priority, k, N))
    pairs.sort()
    return pairs


# ---------------------------------------------------------------
# メインループ
# ---------------------------------------------------------------
def run(k_range=range(-5, 6), n_range=range(-10, 11),
        timeout_sec=60.0,
        retry_timeout=False, retry_error=False,
        db_path=DB_PATH, max_cases=None,
        exclude_k=(-1, 0, 1)):
    conn = init_db(db_path)
    pairs = enumerate_pairs(k_range, n_range, exclude_k=exclude_k)
    processed = 0
    skipped = 0
    print(f"Total candidates: {len(pairs)}  DB: {db_path}")
    try:
        for prio, k, N in pairs:
            if max_cases is not None and processed >= max_cases:
                break
            if already_done(conn, k, N,
                            retry_timeout=retry_timeout,
                            retry_error=retry_error):
                skipped += 1
                continue
            print(f"[{processed+1}] k={k}, N={N} (priority={prio})")
            rec = analyze(k, N, timeout_sec=timeout_sec)
            save_result(conn, rec)
            processed += 1
    except KeyboardInterrupt:
        print("\n\n*** KeyboardInterrupt: 中断されました ***")
        print(f"processed={processed}, skipped={skipped}")
    finally:
        conn.close()
    print(f"\nDone. processed={processed}, skipped(existing)={skipped}")


# ---------------------------------------------------------------
# 未確定だけ再挑戦
# ---------------------------------------------------------------
def rerun_undetermined(timeout_sec=300, db_path=DB_PATH):
    """status='ok' だが rank が None のものを再計算."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    targets = cur.execute(
        "SELECT k, n FROM curves WHERE status='ok' AND rank IS NULL"
    ).fetchall()
    conn.close()
    print(f"re-running {len(targets)} undetermined cases "
          f"(timeout={timeout_sec}s)...")
    conn = init_db(db_path)
    processed = 0
    try:
        for k, N in targets:
            print(f"  [{processed+1}/{len(targets)}] k={k}, N={N}")
            rec = analyze(k, N, timeout_sec=timeout_sec)
            save_result(conn, rec)
            processed += 1
    except KeyboardInterrupt:
        print(f"\n\n*** KeyboardInterrupt: 中断されました ({processed}/{len(targets)} 完了) ***")
    finally:
        conn.close()


def rerun_failures(timeout_sec=300, db_path=DB_PATH,
                   statuses=("timeout", "crashed", "error")):
    """timeout/crashed/error を再計算."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    qmarks = ",".join("?" * len(statuses))
    # timeout_sec と elapsed_sec も取得して判定に使う
    targets = cur.execute(
        f"SELECT k, n, timeout_sec, elapsed_sec FROM curves WHERE status IN ({qmarks})",
        statuses
    ).fetchall()
    conn.close()
    
    # 前回の試行時間が今回の制限時間以上なら無駄なのでスキップ
    filtered = []
    skipped = 0
    for k, N, prev_timeout, prev_elapsed in targets:
        # 前回のタイムアウト設定が今回より長い、または実際の経過時間が今回の制限以上
        if prev_timeout and prev_timeout >= timeout_sec:
            if prev_elapsed and prev_elapsed >= timeout_sec * 0.95:  # 95%以上使い切っていたら
                skipped += 1
                continue
        filtered.append((k, N))
    
    print(f"re-running {len(filtered)} failed cases (skipped {skipped} with sufficient timeout) "
          f"(timeout={timeout_sec}s) statuses={statuses}...")
    
    conn = init_db(db_path)
    processed = 0
    try:
        for k, N in filtered:
            print(f"  [{processed+1}/{len(filtered)}] k={k}, N={N}")
            rec = analyze(k, N, timeout_sec=timeout_sec)
            save_result(conn, rec)
            processed += 1
    except KeyboardInterrupt:
        print(f"\n\n*** KeyboardInterrupt: 中断されました ({processed}/{len(filtered)} 完了) ***")
    finally:
        conn.close()


# ---------------------------------------------------------------
# 集計・検索ユーティリティ
# ---------------------------------------------------------------
def summary(db_path=DB_PATH):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    print("== status counts ==")
    for row in cur.execute(
        "SELECT status, COUNT(*) FROM curves "
        "GROUP BY status ORDER BY 2 DESC"):
        print(f"  {row[0]:12s}: {row[1]}")

    print("\n== rank distribution (status='ok', rank determined) ==")
    for row in cur.execute(
        "SELECT rank, COUNT(*) FROM curves "
        "WHERE status='ok' AND rank IS NOT NULL "
        "GROUP BY rank ORDER BY rank"):
        print(f"  rank {row[0]}: {row[1]}")

    print("\n== torsion distribution (status='ok') ==")
    for row in cur.execute(
        "SELECT torsion_struct, COUNT(*) FROM curves "
        "WHERE status='ok' "
        "GROUP BY torsion_struct ORDER BY 2 DESC"):
        print(f"  tors={row[0]}: {row[1]}")

    print("\n== rank not determined (status='ok' with only bounds) ==")
    rows = cur.execute(
        "SELECT k, n, rank_bounds FROM curves "
        "WHERE status='ok' AND rank IS NULL"
    ).fetchall()
    if not rows:
        print("  (none)")
    for row in rows:
        print(f"  k={row[0]}, N={row[1]}, bounds={row[2]}")

    print("\n== high-rank examples (rank >= 2) ==")
    for row in cur.execute(
        "SELECT k, n, rank, torsion_struct, gens FROM curves "
        "WHERE status='ok' AND rank >= 2 "
        "ORDER BY rank DESC, ABS(k)+ABS(n) LIMIT 50"):
        print(f"  k={row[0]:>3}, N={row[1]:>3}, rank={row[2]}, "
              f"tors={row[3]}, gens={row[4]}")

    print("\n== singular cases ==")
    for row in cur.execute(
        "SELECT k, n, note FROM curves WHERE status='singular' "
        "ORDER BY ABS(k)+ABS(n)"):
        print(f"  k={row[0]:>3}, N={row[1]:>3}: {row[2]}")

    conn.close()


def query(k, N, db_path=DB_PATH):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute("SELECT * FROM curves WHERE k=? AND n=?", (k, N))
    row = cur.fetchone()
    conn.close()
    if row is None:
        print(f"(k={k}, N={N}) not in DB.")
        return None
    for key in row.keys():
        print(f"  {key:14s}: {row[key]}")
    return dict(row)


def export_csv(out_path="elliptic_curves.csv", db_path=DB_PATH):
    """CSV へ書き出し (Excel/解析用)。"""
    import csv
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    rows = cur.execute("SELECT * FROM curves ORDER BY ABS(k)+ABS(n), k, n")
    headers = [d[0] for d in rows.description] if False else None
    # rows.description は使えないので一度全部取得
    rows = cur.execute("SELECT * FROM curves ORDER BY ABS(k)+ABS(n), k, n").fetchall()
    if not rows:
        print("DB is empty.")
        conn.close()
        return
    headers = rows[0].keys()
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in rows:
            w.writerow([r[h] for h in headers])
    print(f"exported {len(rows)} rows -> {out_path}")
    conn.close()


# ---------------------------------------------------------------
# fork で起動（spawn だと Sage の初期化が重い）
# ---------------------------------------------------------------
if mp.get_start_method(allow_none=True) is None:
    mp.set_start_method("fork")


# ---------------------------------------------------------------
# スクリプト直接起動時のデフォルト動作
# ---------------------------------------------------------------
# File: elliptic_db.sage.py (末尾の __main__ ブロックのみ)
if __name__ == "__main__":
    # ============================================================
    # 使いたいブロックだけ、そのブロックを囲む """ """ を外す。
    # 複数同時に有効化してもよい（上から順に実行される）。
    # ============================================================

    # [3] 退化チェックは analyze 内に任せて全 k を計算
    
    run(k_range=range(0, 21),
        n_range=range(-30, 61),
        timeout_sec=8.0,
        exclude_k=())
    

    # [4] さらに広い範囲で計算（時間がかかる）
    """
    run(k_range=range(-20, 21),
        n_range=range(-30, 31),
        timeout_sec=30.0)
    """

    # [5] timeout のみを長めの制限時間で再試行
    """
    rerun_failures(timeout_sec=6000, statuses=("timeout",))
    """

    # [6] timeout / crashed / error をまとめて再試行
    """
    rerun_failures(timeout_sec=600)
    """

    # [7] 段階的に timeout を伸ばして再試行（推奨）
    """
    rerun_failures(timeout_sec=120,  statuses=("timeout",))
    rerun_failures(timeout_sec=600,  statuses=("timeout",))
    rerun_failures(timeout_sec=3600, statuses=("timeout",))
    """

    # [8] rank が未確定 (status='ok' だが rank=NULL) のものを再計算
    """
    rerun_undetermined(timeout_sec=600)
    """

    # [9] 既存の timeout も対象に含めて通常ループを回す
    """
    run(k_range=range(-12, 13),
        n_range=range(-15, 16),
        timeout_sec=600.0,
        retry_timeout=True)
    """

    # [10] crashed / error のみ再試行
    """
    rerun_failures(timeout_sec=300, statuses=("crashed", "error"))
    """

    # [11] 単独のケースを確認したい
    """
    query(3, -4)
    query(2, 7)
    """

    # [12] 単独のケースを強制的に再計算して DB へ保存
    """
    _rec = analyze(3, -4, timeout_sec=600)
    _conn = init_db()
    save_result(_conn, _rec)
    _conn.close()
    """

    # [13] CSV へエクスポート
    """
    export_csv("elliptic_curves.csv")
    """

    # [14] 集計表示（最後に呼ぶのが便利）
    summary()