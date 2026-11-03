# -*- coding: utf-8 -*-
"""附錄 A「全書筆誤總表」：把 ref/md/chNN_formulas.md 的 12 張筆誤表合併成 errata.tex。

用法（在 notes/ 裡）：python appendix/make_errata.py
ref/ 不進 git，所以產生的 errata.tex 要 commit；講義編譯時不需要 ref/。
"""
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# 章：(各章筆誤表的列數, 講義裡的章末總表)
CHAPTERS = {
    1: (13, "1.7"), 2: (16, "2.5"), 3: (12, "3.9"), 4: (16, "4.5"),
    5: (10, "5.5"), 6: (12, "6.4"), 7: (13, "7.4"), 8: (7, "8.6"),
    9: (10, "9.8"), 10: (11, "10.7"), 11: (8, "11.8"), 12: (11, "12.8"),
}
HEADING = "## 本章的筆誤一覽（給 2004 年版）"


def split_row(line):
    """以 | 切欄；$...$ 裡的 | 是絕對值，不是欄的分隔。"""
    cells, cur, in_math = [], "", False
    for ch in line.strip().strip("|"):
        if ch == "$":
            in_math = not in_math
        if ch == "|" and not in_math:
            cells.append(cur.strip())
            cur = ""
        else:
            cur += ch
    assert not in_math, line
    cells.append(cur.strip())
    return cells


def to_tex(cell):
    """數學原樣保留；文字部分處理 LaTeX 特殊字元與 markdown 的 *強調*。"""
    parts = re.split(r"(\$[^$]*\$)", cell)
    for i in range(0, len(parts), 2):
        s = re.sub(r"([%&#_])", r"\\\1", parts[i])
        parts[i] = s
    s = "".join(parts)
    return re.sub(r"(?<![\w$])\*([A-Za-z ]+)\*", r"\\emph{\1}", s)


def read_rows(ch):
    lines = (ROOT / f"ref/md/ch{ch:02d}_formulas.md").read_text(encoding="utf-8").split("\n")
    start = lines.index(HEADING)
    rows = []
    for line in lines[start + 1:]:
        if line.startswith("|") and not line.startswith(("| 位置", "|---")):
            rows.append(split_row(line))
        elif rows and not line.strip():
            break
    for r in rows:
        assert len(r) == 3, r
    return rows


def main():
    out = [
        "% 由 appendix/make_errata.py 產生，不要手改；來源是 ref/md/chNN_formulas.md 的筆誤表。",
        r"\clearpage",
        r"\appendix",
        r"\section{全書筆誤總表}\label{app-errata}",
    ]
    tables, total = [], 0
    for ch, (expected, _) in CHAPTERS.items():
        rows = read_rows(ch)
        assert len(rows) == expected, (ch, len(rows), expected)
        total += expected
        tables.append(rows)
    out += [
        f"把第 1--12 章章末的筆誤表合在一起，共 {total} 條，全部針對 2004 年版。"
        r"頁碼與式號都是書上的；標「〔注〕」的條目，在該章本文另有說明。"
        r"只是記號或寫法上的小問題（各章校對稿裡的「其他〔注〕」）不收在這裡。",
        r"\begin{center}",
        r"\small\setlength{\tabcolsep}{5pt}",
        r"\begin{tabular}{@{}rrl@{\qquad}rrl@{}}",
        r"\toprule",
        r"章 & 條數 & 章末總表 & 章 & 條數 & 章末總表\\",
        r"\midrule",
    ]
    for ch in range(1, 7):
        a, b = CHAPTERS[ch], CHAPTERS[ch + 6]
        out.append(f"{ch} & {a[0]} & {a[1]} 節 & {ch + 6} & {b[0]} & {b[1]} 節\\\\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{center}", ""]
    for ch, rows in zip(CHAPTERS, tables):
        out += [
            r"\needspace{7\baselineskip}",
            rf"\subsection{{第 {ch} 章（{len(rows)} 條）}}",
            r"{\footnotesize\renewcommand{\arraystretch}{1.1}\setlength{\tabcolsep}{4pt}",
            r"\begin{longtable}{@{}>{\raggedright\arraybackslash}p{0.21\linewidth}"
            r">{\raggedright\arraybackslash}p{0.34\linewidth}"
            r">{\raggedright\arraybackslash}p{0.41\linewidth}@{}}",
            r"\toprule 位置 & 書上 & 應為\\ \midrule\endfirsthead",
            r"\toprule 位置 & 書上 & 應為\\ \midrule\endhead",
            r"\bottomrule\endfoot",
        ]
        out += [" & ".join(to_tex(c) for c in r) + r"\\" for r in rows]
        out += [r"\end{longtable}}", ""]
    (HERE / "errata.tex").write_bytes("\n".join(out).encode("utf-8"))
    print("errata.tex:", total, "條")


if __name__ == "__main__":
    main()
