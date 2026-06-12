# File: compact.py
import io
import sys
import tokenize


def compact_python(source: str) -> str:
    """
    Pythonコードを動作を維持したまま可能な限り短くする。
    - コメント（# ...）は保持する
    - 文字列リテラル（docstring含む）は中身を一切変更しない
    - 空行・余分な空白を削除し、インデントを1スペース化
    - 同インデントの単純文を `;` で連結
    """
    src = source.replace("\r\n", "\n").replace("\r", "\n")
    if not src.endswith("\n"):
        src += "\n"

    # ソース全体を一括でトークナイズ
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except (tokenize.TokenError, IndentationError):
        # トークナイズ失敗時は原文をそのまま返す
        return source

    # 行ごとに「コードトークン列」「コメント」を集める
    # 複数行にまたがる文字列トークンは「開始行」に属させる
    line_tokens: dict[int, list[tokenize.TokenInfo]] = {}
    line_comments: dict[int, list[str]] = {}
    line_indent: dict[int, int] = {}  # 各論理行のインデント段数

    current_indent = 0  # INDENT/DEDENT で更新

    for tok in tokens:
        ttype = tok.type
        if ttype == tokenize.INDENT:
            current_indent += 1
            continue
        if ttype == tokenize.DEDENT:
            current_indent -= 1
            continue
        if ttype in (tokenize.ENCODING, tokenize.ENDMARKER, tokenize.NL):
            continue
        if ttype == tokenize.NEWLINE:
            continue

        start_line = tok.start[0]

        if ttype == tokenize.COMMENT:
            line_comments.setdefault(start_line, []).append(tok.string)
            line_indent.setdefault(start_line, current_indent)
            continue

        line_tokens.setdefault(start_line, []).append(tok)
        line_indent.setdefault(start_line, current_indent)

    # 行番号順に処理
    all_lines = sorted(set(line_tokens.keys()) | set(line_comments.keys()))

    out_lines: list[tuple[int, str, bool]] = []
    # (インデント段数, 内容, 連結可能フラグ)

    for ln in all_lines:
        indent = line_indent.get(ln, 0)
        toks = line_tokens.get(ln, [])
        cmts = line_comments.get(ln, [])

        if not toks:
            # コメントのみの行
            out_lines.append((indent, " ".join(cmts), False))
            continue

        code = _compact_tokens(toks)

        if cmts:
            code = code + " " + " ".join(cmts)
            mergeable = False
        else:
            mergeable = _is_simple_statement(code)

        out_lines.append((indent, code, mergeable))

    # 同インデントで連続する単純文を `;` で連結
    merged = _merge_simple(out_lines)

    # 文字列化
    result_lines = [(" " * ind) + body for ind, body in merged]
    return "\n".join(result_lines) + "\n"


def _compact_tokens(toks: list) -> str:
    """トークン列の間の余分な空白を最小化して文字列化する。"""
    out = ""
    prev_type = None
    prev_str = ""
    for t in toks:
        s = t.string
        typ = t.type
        if prev_type is None:
            out = s
        else:
            need = _needs_space(prev_type, prev_str, typ, s)
            out += (" " if need else "") + s
        prev_type = typ
        prev_str = s
    return out


def _is_word_like(typ: int, s: str) -> bool:
    if typ == tokenize.NAME:
        return True
    if typ == tokenize.NUMBER:
        return True
    # 文字列でもプレフィックス付き（f, r, b, u, rb等）は識別子と隣接禁止
    if typ == tokenize.STRING and s and s[0] not in ("'", '"'):
        return True
    return False


def _needs_space(pt: int, ps: str, ct: int, cs: str) -> bool:
    # 識別子・数値・キーワード同士の境界はスペース必須
    if _is_word_like(pt, ps) and _is_word_like(ct, cs):
        return True
    # 数値の直後に `.` が続く場合、`1.x` のような曖昧を避けるためスペース
    # 例: `1 .bit_length()` のようなケース
    if pt == tokenize.NUMBER and ct == tokenize.OP and cs == ".":
        # 整数末尾が数字のみなら `.` は小数点と解釈される可能性
        if ps.isdigit():
            return True
    # 文字列の直後に識別子（暗黙連結防止）
    if pt == tokenize.STRING and ct == tokenize.NAME:
        return True
    # 識別子の直後に文字列プレフィックス付き文字列
    if pt == tokenize.NAME and ct == tokenize.STRING:
        if cs and cs[0] not in ("'", '"'):
            return True
    return False


def _is_simple_statement(code: str) -> bool:
    """`;` 連結対象にできる単純文か判定する。"""
    s = code.strip()
    if not s:
        return False
    if s.endswith(":"):
        return False
    if s.startswith("@"):
        return False
    # 制御・宣言系キーワード
    head = s.split(" ", 1)[0].split("(", 1)[0]
    block_kw = {
        "if", "elif", "else", "for", "while", "try", "except",
        "finally", "with", "def", "class", "async",
    }
    if head in block_kw:
        return False
    return True


def _merge_simple(items: list) -> list:
    """同インデントで連続する単純文を `;` で連結する。"""
    result = []
    buf_indent = -1
    buf = []

    def flush():
        nonlocal buf_indent, buf
        if buf:
            if len(buf) == 1:
                result.append((buf_indent, buf[0]))
            else:
                result.append((buf_indent, ";".join(buf)))
            buf = []
            buf_indent = -1

    for indent, code, mergeable in items:
        if mergeable:
            if buf_indent == -1:
                buf_indent = indent
                buf.append(code)
            elif indent == buf_indent:
                buf.append(code)
            else:
                flush()
                buf_indent = indent
                buf.append(code)
        else:
            flush()
            result.append((indent, code))
    flush()
    return result


def main():
    if len(sys.argv) < 2:
        src = sys.stdin.read()
        sys.stdout.write(compact_python(src))
        return

    in_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) >= 3 else None

    with open(in_path, "r", encoding="utf-8") as f:
        src = f.read()

    result = compact_python(src)

    if out_path:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(result)
        before = len(src)
        after = len(result)
        ratio = (after / before * 100) if before else 0
        print(f"{before} -> {after} bytes ({ratio:.1f}%)", file=sys.stderr)
    else:
        sys.stdout.write(result)


if __name__ == "__main__":
    main()