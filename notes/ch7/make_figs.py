r"""
make_figs.py —— 隨機漫步導讀第 7 章（Rudnick & Gaspari《Elements of the Random Walk》第 7 章，動畫 W701–W706）的圖與表

在 notes/ 或 notes/ch7/ 裡執行：
    python ch7/make_figs.py
會先檢查 ch7/part7.tex 寫到的每一個數字（check_numbers_ch7；最後再把 part7.tex 裡出現的每一個小數和檢查過的數比對一遍），
再寫出 ch7/figs.tex（part7.tex 第一行 \input{ch7/figs}）。圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算）。
數學一律由 codes/utils/walk_saw.py 算出；例子與動畫 W701–W706 相同。巨集名稱都以 \Seven 開頭（不和前幾章的 \Fig…、\Two…、
\Three…、\Four…、\Five…、\Six… 撞名）；part7.tex 的 \label 用 seven-…。顏色沿用導言區：cBar 精確值、cWalk 路線、cExact 更正、
cGauss 近似、cWrong 筆誤。
"""
from __future__ import annotations

import math
import re
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent            # notes/ch7
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_saw as SW  # noqa: E402

OUT = HERE / "figs.tex"
TEX = HERE / "part7.tex"
F = Fraction
NS_ALL = tuple(range(1, 15))                       # 精確列舉 N = 1..14
OK = set()                                         # 檢查過的小數字串（part7.tex 的小數都要在這裡）


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣），並記進 OK。"""
    s = f"{float(x):.{d}f}"
    OK.add(s.lstrip("-"))
    return s


def _fr(x):
    x = F(x)
    if x.denominator == 1:
        return str(x.numerator)
    sgn = "-" if x < 0 else ""
    return rf"{sgn}\frac{{{abs(x.numerator)}}}{{{x.denominator}}}"


def _close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) < tol


def _blocks(m, n):
    """弧線圖 m（n 個點）在哪些偶數切口切得開；回傳區塊大小（弧的條數）的 tuple。"""
    cuts = [c for c in range(2, n, 2) if all(not (a <= c < b) for a, b in m)]
    edges = [0] + cuts + [n]
    return tuple((edges[i + 1] - edges[i]) // 2 for i in range(len(edges) - 1))


# =============================================================================
# 數字檢查（part7.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch7():
    # ---- W701 路徑和
    assert SW.split_terms(4, 2, 0, 1) == [(-2, 1, 1), (0, 2, 2), (2, 1, 1)] and SW.split_sum(4, 2, 0, 1) == 6
    assert SW.count_brute(4, 0, 1) == 6 == math.comb(4, 2)
    t2 = SW.split_terms(4, 2, (0, 0), 2)
    assert len(t2) == 9
    cen = [a * b for r, a, b in t2 if r == (0, 0)]
    axis = [a * b for r, a, b in t2 if r != (0, 0) and 0 in r]
    diag = [a * b for r, a, b in t2 if 0 not in r]
    assert cen == [16] and axis == [1] * 4 and diag == [4] * 4
    assert all(a == b for r, a, b in t2) and {a for r, a, b in t2} == {4, 1, 2}
    assert 16 + 4 * 1 + 4 * 4 == 36 == SW.count_brute(4, (0, 0), 2) == math.comb(4, 2) ** 2 == SW.count2(4, 0, 0)
    assert SW.EX_SPLIT["eq27_1d"] == (2, 3) and SW.EX_SPLIT["eq27_2d"] == (4, 9) and SW.count1(3, 1) == 3
    assert SW.count2(3, 1, 0) == 9 and 2 * 3 == 6 and 4 * 9 == 36
    assert SW.split_sum(4, 3, 0, 1) == 6 and SW.split_sum(4, 3, (0, 0), 2) == 36
    assert SW.multi_split(6, [2, 2, 2], 0, 1) == 20 == math.comb(6, 3) == SW.count_brute(6, 0, 1)
    assert SW.multi_split(6, [2, 2, 2], (0, 0), 2) == 400 == math.comb(6, 3) ** 2
    # Ex. 7.1：M1 = 3、M2 = 6（三維，a = 1）→ 分量變異數 1、2、3
    assert SW.EX71 == {"M1": 3, "M2": 6, "var": (1, 2, 3)}
    assert tuple(_s(SW.gauss1(M, 0), 4) for M in (3, 6, 9)) == ("0.3989", "0.2821", "0.2303")
    assert _close(SW.gauss1(9, 0), 1 / math.sqrt(2 * math.pi * 3))
    assert all(abs(SW.gauss_conv(3, 6, x) - SW.gauss1(9, x)) < 1e-8 for x in (-4, -1, 0, 0.7, 2.5))
    assert _close(SW.conv_ratio(3, 6), 1, 1e-8)
    br = SW.conv_ratio(3, 6, book=True)
    assert (_s(br, 4), _s(br ** 3, 4)) == ("0.5774", "0.1925") and _close(br, 3 ** -0.5, 1e-8)
    # (7.11) 的分母 3Ma²：一個分量的變異數 M/2（M = 6：3），格子漫步 M/3（2）
    assert SW.lattice_var(6) == 2 and SW.summary()["W701_var711"] == ("2.0000", "3.0000")

    # ---- W702 自迴避漫步
    assert len(SW.saw_list(2)) == 12 and 4 ** 2 - 12 == 4
    s3 = SW.saw_list(3)
    assert len(s3) == 36 and 4 ** 3 == 64 and len([w for w in s3 if w[1] == (1, 0)]) == 9
    c = [SW.saw_c(n) for n in range(15)]
    assert c == SW.SAW_KNOWN_C[:15]
    assert (c[2], c[3], c[4], c[8], c[14]) == (12, 36, 100, 5916, 2374444)
    for n in range(1, 11):
        bc, _ = SW.saw_brute(n)
        assert bc == c[n]
    assert _s(SW.saw_frac(14), 4) == "0.0088" and _s(SW.saw_R2(14), 2) == "42.79"
    for n in NS_ALL:                                       # 表裡的每一格
        _s(SW.saw_frac(n), 4), _s(SW.saw_R2(n), 4)
    assert _s(SW.mu_eff(14), 4) == "2.6936" and _close(F(c[14], c[13]), SW.mu_eff(14))
    OK.add("2.638")                                         # 連接常數（只引用）
    sw, rw = SW.sample_saw(10), SW.sample_rw(10)
    assert SW.SAMPLE_SEEDS == {"saw": 702, "rw": 702}
    assert sw[-1] == (-3, 3) and sum(v * v for v in sw[-1]) == 18 and len(set(sw)) == 11
    assert rw[-1] == (1, -1) and sum(v * v for v in rw[-1]) == 2 and len(rw) - len(set(rw)) == 3

    # ---- W703 1 − f 展開
    T4 = SW.ie_terms(4)
    assert T4 == [256, 228, 96, 28, 4] == SW.ie_terms_subsets(4)
    assert SW.order1_formula(4) == 228 == 3 * 4 ** 2 * 4 + 36 and 3 * 4 ** 2 * 4 == 192
    assert SW.count2(2, 0, 0) == 4 and SW.count2(4, 0, 0) == 36
    assert SW.ie_partial(4) == [256, 28, 124, 96, 100] and 256 - 228 == 28 and 28 + 96 == 124 and 124 - 28 == 96
    assert SW.order2_types(4) == {"sequential": 0, "nested": 16, "crossed": 32, "shared": 48} and 16 + 32 + 48 == 96
    P6, P8 = SW.ie_partial(6), SW.ie_partial(8)
    assert P6[-1] == 780 and len(P6) == 10
    assert len(P8) == 17 and P8[:4] == [65536, -119332, 190940, -222036] and P8[-1] == 5916
    assert tuple(SW.chord_count(k) for k in range(1, 5)) == (1, 3, 15, 105)
    assert len(SW.matchings(2)) == 3 and len(SW.matchings(3)) == 15
    assert SW.SECOND_ORDER == {"a": ("nested", ((1, 4), (2, 3))), "b": ("sequential", ((1, 2), (3, 4))),
                               "c": ("crossed", ((1, 3), (2, 4)))}

    # ---- W704 生成函數、Dyson
    assert SW.EX_CONV == {"H": [1, 2, 4, 8, 16], "conv": [1, 4, 12, 32, 80]}
    assert SW.convolve([1, 2, 4, 8, 16], [1, 2, 4, 8, 16])[:5] == [1, 4, 12, 32, 80] == [(n + 1) * 2 ** n for n in range(5)]
    gs, gc = SW.gf_check(0.1)
    assert _s(gc, 4) == "1.5625" and _close(gs, gc, 1e-9) and _close(1 / (1 - 0.2) ** 2, 1.5625)
    OK.update({"0.1"})
    assert SW.fourier_conv_error() < 1e-12
    g0, sg = SW.DYSON_TOY
    assert (g0, sg) == (2, F(1, 4)) and sg * g0 == F(1, 2)
    assert SW.dyson_partial(g0, sg, 5) == [2, 3, F(7, 2), F(15, 4), F(31, 8)] and SW.dyson_closed(g0, sg) == 4
    assert SW.dyson_book35(g0, sg) == F(4, 3) and _close(SW.dyson_book33(g0, sg), F(10, 3))
    full, irr, D = SW.dyson_counts(4)
    assert full == [1, 1, 3, 15, 105] == D and irr == [0, 1, 2, 10, 74]
    assert (1 + 2, 3 + 2 + 10, 15 + 2 * 3 + 10 + 74) == (3, 15, 105)
    assert SW.irreducible_count(3) == 10 and 15 - 10 == 5 and SW.reversal_classes(2) == 3 and SW.reversal_classes(3) == 11
    blk = {}
    for m in SW.matchings(3):
        b = _blocks(m, 6)
        assert (len(b) == 1) == SW.is_irreducible(m)
        blk[b] = blk.get(b, 0) + 1
    assert blk == {(3,): 10, (1, 1, 1): 1, (1, 2): 2, (2, 1): 2}
    assert sum(SW.is_irreducible(m) for m in SW.matchings(2)) == 2
    # 幾何級數的收斂：一維接觸作用 f = 0.05、q = 0.5
    sgz = [0.05 * SW.sigma1_1d(z) * SW.g0_1d(z, 0.5) for z in (0.4, 0.49)]
    assert (_s(sgz[0], 2), _s(sgz[1], 2)) == ("0.28", "1.80")
    OK.update({"0.05", "0.5", "0.4", "0.49"})

    # ---- W705 冪次、平均場
    pt = SW.power_table()
    assert [str(v) for v in pt[1]] == ["1", "-1/2", "-2", "-7/2", "-5"]
    assert [str(v) for v in pt[2]] == ["1", "0", "-1", "-2", "-3"]
    assert [str(v) for v in pt[3]] == ["1", "1/2", "0", "-1/2", "-1"]
    assert [str(v) for v in pt[4]] == ["1"] * 5 and [str(v) for v in pt[5]] == ["1", "3/2", "2", "5/2", "3"]
    assert all(SW.power_exp(l, d) == SW.power_exp_count(l, d) for d in range(1, 6) for l in range(5))
    assert tuple(_s(SW.loop_slope(k), 2) for k in ("l1d1", "l1d3", "l2d1")) == ("-0.50", "0.50", "-2.00")
    sig = tuple(_s(SW.sigma1_1d(z), 4) for z in (0.4, 0.49, 0.499))
    assert sig == ("1.6667", "5.0252", "15.8193") and abs(SW.sigma1_1d(0.49) - SW.sigma1_1d_quad(0.49)) < 1e-9
    OK.add("0.499")
    r = SW.sigma1_1d(0.499) / SW.sigma1_1d(0.49)
    assert abs(r / math.sqrt(10) - 1) < 0.01                                 # 近 10 倍 → 大 √10 倍
    assert _s(SW.eq747_residual(0.4, 1, book=True), 4) == "0.3333" and abs(SW.eq747_residual(0.4, 1)) < 1e-12
    assert _close(SW.eq755(0.25, 0.2, 0.25), 5) and _close(1 / (0.25 - 0.2), 20)
    assert abs(SW.eq748_eigen(2.0) + 4) < 1e-4 and abs(SW.eq748_eigen(2.0, book=True) - 4) < 1e-4
    mf = [SW.cmf_ratio(*cs) for cs in SW.MF_CASES]
    assert SW.MF_CASES == ((10, 100, 1), (100, 10 ** 4, 1), (1000, 10 ** 6, 1))
    assert tuple(_s(v, 4) for v in mf) == ("0.3385", "0.3643", "0.3675") and _s(math.exp(-1), 4) == "0.3679"
    assert all(_close(SW.cmf_asym(*cs), math.exp(-1)) for cs in SW.MF_CASES)
    assert tuple(f"{100 * (1 - v / math.exp(-1)):.1f}" for v in mf) == ("8.0", "1.0", "0.1")
    OK.update({"8.0", "1.0", "0.1"})
    assert abs(SW.cmf_ratio(10, 100, 0) - 1) < 1e-9 and abs(SW.cmf_book59(10, 100, 0) / 100.0 ** -11 - 1) < 1e-6
    e75 = tuple(_s(SW.ex75(n), 4) for n in (10, 100, 1000, 10000))
    assert e75 == ("0.6899", "0.9341", "0.9922", "0.9992")
    assert (_s(1 - 2 / (1000 * 0.25), 3), _s(1 - 2 / (10000 * 0.25), 4)) == ("0.992", "0.9992")
    eb = tuple(_s(SW.ex75_big(n), 4) for n in (100, 1000, 10000))
    assert eb == ("0.0452", "0.0596", "0.0621") and _s((1.5 - 1) ** 2 / 4, 4) == "0.0625"

    # ---- W706 Flory
    Rs = SW.flory_Rstar(1000, 2)
    assert _s(Rs, 2) == "177.83" and _close(Rs ** 4, 1e9, 1e-3) and abs(SW.flory_Rstar_num(1000, 2) / Rs - 1) < 1e-8
    t1, t2_ = -Rs ** 2 / 1000, -1000 ** 2 / Rs ** 2
    assert _s(t1, 2) == _s(t2_, 2) == "-31.62" and _s(SW.flory_F(Rs, 1000, 2), 2) == "-63.25"
    assert tuple(_s(SW.flory_Rstar(1000, d), 2) for d in (1, 2, 3, 4)) == ("793.70", "177.83", "68.43", "35.50")
    assert [SW.flory_nu(d) for d in (1, 2, 3, 4)] == [1, F(3, 4), F(3, 5), F(1, 2)]
    assert tuple(_s(SW.flory_slope(d), 4) for d in (1, 2, 3, 4)) == ("1.0000", "0.7500", "0.6000", "0.5000")
    OK.add("0.6")
    assert SW.NU_QUOTED == {1: 1.0, 2: 0.75, 3: 0.588, 4: 0.5}
    OK.add("0.588")
    assert tuple(_s(SW.two_nu_eff(n), 4) for n in (6, 8, 10, 12, 14)) == ("1.4306", "1.4372", "1.4442", "1.4505", "1.4556")
    assert tuple(_s(SW.two_nu_eff(n), 4) for n in (7, 9, 11, 13)) == ("1.4459", "1.4528", "1.4562", "1.4593")


def check_tex_decimals():
    """part7.tex（去掉註解）裡的每一個小數都必須是上面檢查過的數。"""
    if not TEX.exists():
        print("（尚無 part7.tex，略過小數比對）")
        return
    txt = re.sub(r"(?<!\\)%.*", "", TEX.read_text(encoding="utf-8"))
    txt = re.sub(r"[\d.]+\\linewidth|[\d.]+pt\b|\\arraystretch\}\{[\d.]+\}|\(2\.7\)", "", txt)   # 版面尺寸、第 2 章的式號
    eqno = re.compile(r"7\.\d{1,2}")                     # 本章的式號、節號、圖號、練習號（7.1–7.64）
    bad = sorted({m for m in re.findall(r"\d+\.\d+", txt) if m not in OK and not eqno.fullmatch(m)})
    assert not bad, f"part7.tex 有未檢查的小數：{bad}"


# =============================================================================
# 圖與表
# =============================================================================
def fig_grid():
    """二維 4 步回原點：9 個中間點上寫 C(2; r)；下面是一維。"""
    t2 = SW.split_terms(4, 2, (0, 0), 2)
    out = [r"\begin{tikzpicture}[x=0.5cm,y=0.5cm,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[gray!30] (-2.5,-2.5) grid (2.5,2.5);",
           r"\draw[gray!60,->] (-2.7,0) -- (2.8,0); \draw[gray!60,->] (0,-2.7) -- (0,2.8);"]
    for (x, y), a, b in t2:
        fill = {4: "cBar!45", 2: "cBar!25", 1: "cBar!10"}[a]
        out.append(rf"\node[circle,draw=cBar,fill={fill},inner sep=1.2pt,minimum size=11pt] at ({x},{y}) {{${a}$}};")
    out.append(r"\node at (0,-3.2) {$\sum_{\vec r}C(2;\vec r)^2=36$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_gauss():
    """Ex. 7.1：變異數 1、2 的高斯摺積成變異數 3（點：數值摺積）。"""
    xs = np.linspace(-6, 6, 73)
    out = [r"\begin{tikzpicture}[x=0.36cm,y=6.2cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (-6.4,0) -- (6.6,0) node[right] {$x$};",
           r"\draw[->] (0,0) -- (0,0.47);",
           r"\foreach \x in {-6,-3,3,6} {\draw (\x,0) -- (\x,-0.012); \node[below] at (\x,-0.01) {$\x$};}",
           r"\foreach \y in {0.1,0.2,0.3,0.4} {\draw (0,\y) -- (-0.2,\y); \node[left] at (-0.2,\y) {$\y$};}"]
    leg = []
    for k, (M, sty, var) in enumerate(((3, "cBar,densely dashed,thick", 1), (6, "cWalk,densely dashed,thick", 2),
                                       (9, "cExact,very thick", 3))):
        pts = " ".join(f"({x:.3f},{SW.gauss1(M, x):.4f})" for x in xs)
        out.append(rf"\draw[{sty}] plot coordinates {{{pts}}};")
        yl = 0.44 - 0.055 * k
        leg.append(rf"\draw[{sty}] (2.0,{yl:.3f}) -- (3.2,{yl:.3f}) node[right,black] {{變異數 {var}：{_s(SW.gauss1(M, 0), 4)}}};")
    out += leg
    for x in range(-5, 6):
        out.append(rf"\fill ({x},{SW.gauss_conv(3, 6, x):.4f}) circle (1.1pt);")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_split():
    e = SW.EX_SPLIT
    rows = [
        r"4 步，$M=2$ & $1+4+1=" + str(SW.split_sum(4, 2, 0, 1)) + r"$ & $16+4+16=" + str(SW.split_sum(4, 2, (0, 0), 2)) + r"$\\",
        rf"4 步，$M=3$ & ${e['eq27_1d'][0]}\cdot{e['eq27_1d'][1]}={SW.split_sum(4, 3, 0, 1)}$ & "
        rf"${e['eq27_2d'][0]}\cdot{e['eq27_2d'][1]}={SW.split_sum(4, 3, (0, 0), 2)}$\\",
        rf"6 步，$2+2+2$ & {SW.multi_split(6, [2, 2, 2], 0, 1)} & {SW.multi_split(6, [2, 2, 2], (0, 0), 2)}\\",
        rf"暴力列舉 & {SW.count_brute(4, 0, 1)}，{SW.count_brute(6, 0, 1)} & "
        rf"{SW.count_brute(4, (0, 0), 2)}，{SW.count_brute(6, (0, 0), 2)}\\",
    ]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"回原點 & 一維 & 二維\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_twostep():
    """二維 16 條 2 步漫步：12 條自迴避（橙）、4 條原地折返（灰、打叉）。"""
    saws = set(SW.saw_list(2))
    out = [r"\begin{tikzpicture}[x=0.17cm,y=0.17cm,font=\scriptsize]"]
    names = ("第一步向東", "向北", "向西", "向南")
    for g, (ax, ay) in enumerate(SW.STEPS2):
        for k, (bx, by) in enumerate(SW.STEPS2):
            ox = (g * 4 + k) * 5.0 + g * 2.6
            pts = ((0, 0), (ax, ay), (ax + bx, ay + by))
            out.append(rf"\draw[gray!25] ({ox - 2.2:.1f},-2.2) rectangle ({ox + 2.2:.1f},2.2);")
            out.append(rf"\fill ({ox:.1f},0) circle (0.45);")
            if pts in saws:
                out.append(rf"\draw[cWalk,thick,-{{Stealth[length=2.5pt]}}] ({ox:.1f},0) -- ({ox + ax:.1f},{ay}) -- "
                           rf"({ox + ax + bx:.1f},{ay + by});")
            else:
                px, py = -ay * 0.45, ax * 0.45
                out.append(rf"\draw[gray,thick] ({ox:.1f},0) -- ({ox + ax:.1f},{ay}) -- ({ox + ax + px:.2f},{ay + py:.2f}) -- "
                           rf"({ox + px:.2f},{py:.2f});")
                out.append(rf"\node[cWrong,font=\small,inner sep=0pt] at ({ox + 1.6 * ax + px:.2f},{1.6 * ay + py:.2f}) {{$\times$}};")
        out.append(rf"\node at ({(g * 4 + 1.5) * 5.0 + g * 2.6:.1f},3.4) {{{names[g]}：3＋1}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_saw():
    def block(ns):
        r0 = r"$N$ & " + " & ".join(str(n) for n in ns) + r"\\"
        r1 = r"$c_N$ & " + " & ".join(str(SW.saw_c(n)) for n in ns) + r"\\"
        r2 = r"$c_N/4^N$ & " + " & ".join(_s(SW.saw_frac(n), 4) for n in ns) + r"\\"
        r3 = r"$\langle R^2\rangle$ & " + " & ".join(_s(SW.saw_R2(n), 4) for n in ns) + r"\\"
        return [r0, r"\midrule", r1, r2, r3]
    return "\n".join([r"\begin{tabular}{@{}l*{7}{r}@{}}", r"\toprule", *block(range(1, 8)), r"\midrule",
                      *block(range(8, 15)), r"\bottomrule", r"\end{tabular}"])


def fig_loglog():
    """log⟨R²⟩ 對 log N：自迴避（點）、普通漫步（虛線）、斜率 3/2 的參考線（點線）。"""
    pts = [(math.log2(n), math.log2(float(SW.saw_R2(n)))) for n in NS_ALL]
    x14, y14 = pts[-1]
    out = [r"\begin{tikzpicture}[x=1.0cm,y=0.6cm,>=Stealth,font=\scriptsize]",
           r"\draw[->] (0,0) -- (4.3,0) node[right] {$N$};",
           r"\draw[->] (0,0) -- (0,6.4) node[above] {$\langle R^2\rangle$};",
           r"\foreach \x/\l in {0/1,1/2,2/4,3/8,4/16} {\draw (\x,0) -- (\x,-0.12); \node[below] at (\x,-0.1) {$\l$};}",
           r"\foreach \y/\l in {0/1,2/4,4/16,6/64} {\draw (0,\y) -- (-0.06,\y); \node[left] at (-0.06,\y) {$\l$};}",
           r"\draw[gray,densely dashed] (0,0) -- (4,4) node[below=1pt] {普通漫步 $N$};"]
    xa, xb = 1.0, 4.0
    out.append(rf"\draw[cWalk,densely dotted,thick] ({xa},{y14 + 1.5 * (xa - x14):.3f}) -- ({xb},{y14 + 1.5 * (xb - x14):.3f}) "
               r"node[above] {斜率 $\frac32$};")
    for x, y in pts:
        out.append(rf"\fill[cBar] ({x:.3f},{y:.3f}) circle (1.4pt);")
    out.append(rf"\node[cBar,left] at ({x14 - 0.1:.3f},{y14 + 0.1:.3f}) {{自迴避 {_s(SW.saw_R2(14), 2)}}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def _chord(ox, pairs, n, sp, col, lab=None, dot=1.2):
    out = [rf"\draw[thick] ({ox:.3f},0) -- ({ox + (n + 1) * sp:.3f},0);"]
    for a, b in pairs:
        out.append(rf"\draw[{col},densely dashed,thick] ({ox + a * sp:.3f},0) arc[start angle=180,end angle=0,radius={(b - a) * sp / 2:.3f}];")
    for i in range(1, n + 1):
        out.append(rf"\fill ({ox + i * sp:.3f},0) circle ({dot}pt);")
    if lab:
        out.append(rf"\node[below] at ({ox + (n + 1) * sp / 2:.3f},-0.05) {{{lab}}};")
    return out


def fig_chords():
    """一階（Fig. 7.2）與二階的三種圖（Fig. 7.3）。"""
    so = SW.SECOND_ORDER
    sp = 0.33
    out = [r"\begin{tikzpicture}[font=\scriptsize]"]
    out += _chord(0, ((1, 2),), 2, sp, "cWrong", "一階 (7.19)")
    ox = 3 * sp + 0.45
    lab = {"a": "(a) 巢狀", "b": "(b) 並排", "c": "(c) 交叉"}
    for key in ("a", "b", "c"):
        out += _chord(ox, so[key][1], 4, sp, "cWrong", lab[key])
        ox += 5 * sp + 0.45
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_third():
    """三階的 15 種弧線圖（練習 7.2）：不可約（橙）10 種、可約（灰）5 種。"""
    sp, dx, dy = 0.18, 1.6, 0.85
    out = [r"\begin{tikzpicture}[font=\tiny]"]
    for i, m in enumerate(SW.matchings(3)):
        ox, oy = (i % 5) * dx, -(i // 5) * dy
        col = "cWalk" if SW.is_irreducible(m) else "gray"
        body = _chord(ox, m, 6, sp, col, dot=0.9)
        out.append(rf"\begin{{scope}}[yshift={oy:.2f}cm]")
        out += body
        out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_ie():
    T4, P4, P6, P8 = SW.ie_terms(4), SW.ie_partial(4), SW.ie_partial(6), SW.ie_partial(8)
    named = ((r"$T_K$（$N=4$）", T4, False), (r"$P_K$（$N=4$）", P4, True), (r"$P_K$（$N=6$）", P6, True),
             (r"$P_K$（$N=8$）", P8, True))
    assert (P4[-1], P6[-1], P8[-1]) == (SW.saw_c(4), SW.saw_c(6), SW.saw_c(8))
    out = [r"\begin{tabular}{@{}l*{9}{r}@{}}", r"\toprule"]
    for b, Ks in enumerate((range(0, 9), range(9, 17))):
        if b:
            out.append(r"\midrule")
        out += [r"$K$ & " + " & ".join(str(K) for K in Ks) + (r" & \\" if len(Ks) < 9 else r"\\"), r"\midrule"]
        for lab, col, last in named:
            cells = []
            for K in Ks:
                if K < len(col):
                    v = str(col[K])
                    cells.append(rf"$\mathbf{{{v}}}$" if last and K == len(col) - 1 else f"${v}$")
                else:
                    cells.append("")
            if any(cells):
                out.append(lab + " & " + " & ".join(cells) + (r" & \\" if len(Ks) < 9 else r"\\"))
    return "\n".join(out + [r"\bottomrule", r"\end{tabular}"])


def fig_dyson():
    """Dyson 方程的「火車」：g = g0 + g0Σg0 + …（長條另見 fig_dyson_bars）。"""
    out = [r"\begin{tikzpicture}[font=\small,baseline=(current bounding box.center)]"]
    x = 0.0
    out.append(r"\node[anchor=east] at (0,0) {$g=$};")

    def line(x0, n):
        seg, rr = 0.55, 0.22
        o = []
        xx = x0
        for i in range(n + 1):
            o.append(rf"\draw[thick] ({xx:.2f},0) -- ({xx + seg:.2f},0);")
            xx += seg
            if i < n:
                o.append(rf"\draw[cWrong,thick,fill=cWrong!8] ({xx + rr:.2f},0) circle ({rr}) node {{\scriptsize$\Sigma$}};")
                xx += 2 * rr
        return o, xx

    for n in range(3):
        o, x2 = line(x + 0.1, n)
        out += o
        x = x2 + 0.1
        out.append(rf"\node at ({x + 0.2:.2f},0) {{$+$}};")
        x += 0.4
    out.append(rf"\node[anchor=west] at ({x - 0.1:.2f},0) {{$\cdots=\dfrac{{g_0}}{{1-\Sigma g_0}}=\dfrac{{1}}{{g_0^{{-1}}-\Sigma}}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_dyson_bars():
    """玩具數字 g0 = 2、Σ = 1/4 的部分和長條，逼近 4。"""
    ps = SW.dyson_partial(*SW.DYSON_TOY, 5)
    lim = SW.dyson_closed(*SW.DYSON_TOY)
    w, ys = 0.36, 0.36
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]", r"\begin{scope}"]
    for i, v in enumerate(ps):
        out.append(rf"\draw[cBar,fill=cBar!25] ({i * (w + 0.12):.2f},0) rectangle ({i * (w + 0.12) + w:.2f},{float(v) * ys:.3f});")
        out.append(rf"\node[below] at ({i * (w + 0.12) + w / 2:.2f},0) {{${_fr(v)}$}};")
    xe = 5 * (w + 0.12)
    out.append(rf"\draw[cExact,densely dashed,thick] (-0.1,{float(lim) * ys:.3f}) -- ({xe:.2f},{float(lim) * ys:.3f}) "
               rf"node[right] {{${lim}$}};")
    out.append(rf"\draw (-0.1,0) -- ({xe:.2f},0);")
    out.append(rf"\node[below] at ({xe / 2:.2f},-0.5) {{$g_0=2,\ \Sigma=\frac14$ 的部分和}};")
    out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_dyson():
    H = SW.EX_CONV["H"]
    cv = SW.convolve(H, H)[:5]
    ps = SW.dyson_partial(*SW.DYSON_TOY, 5)
    full, irr, D = SW.dyson_counts(4)
    gs, gc = SW.gf_check(0.1)
    rows = [
        r"摺積 (7.20)：$N$ & 0 & 1 & 2 & 3 & 4 & \\",
        r"\quad $H_1=H_2=2^K$ & " + " & ".join(map(str, H)) + r" & \\",
        r"\quad $\sum_MH_1(M)H_2(N-M)$ & " + " & ".join(map(str, cv)) + rf" & $z=0.1$：{_s(gc, 4)}\\",
        r"\midrule",
        r"Dyson：項數 $n$ & 1 & 2 & 3 & 4 & 5 & $\infty$\\",
        r"\quad $g_0\sum_{k<n}(\Sigma g_0)^k$ & " + " & ".join(f"${_fr(v)}$" for v in ps)
        + rf" & ${SW.dyson_closed(*SW.DYSON_TOY)}$\\",
    ]
    a = "\n".join([r"\begin{tabular}[c]{@{}lrrrrrl@{}}", r"\toprule", *rows, r"\bottomrule", r"\end{tabular}"])
    rows2 = [f"{k} & {full[k]} & " + ("---" if k == 0 else str(irr[k])) + f" & {D[k]}\\\\" for k in range(5)]
    b = "\n".join([r"\begin{tabular}[c]{@{}crrr@{}}", r"\toprule", r"計數版 $k$ & $(2k-1)!!$ & 不可約 $I_k$ & $\frac1{1-I(x)}$\\",
                   r"\midrule", *rows2, r"\bottomrule", r"\end{tabular}"])
    return a + "\n\\qquad\n" + b


def tab_power():
    pt = SW.power_table()
    rows = [f"{d} & " + " & ".join(f"${_fr(v)}$" for v in pt[d]) + r"\\" for d in range(1, 6)]
    return "\n".join([r"\begin{tabular}{@{}c|rrrrr@{}}", r"\toprule", r"$d$ \textbackslash\ $l$ & 0 & 1 & 2 & 3 & 4\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def tab_mf():
    mf = [SW.cmf_ratio(*cs) for cs in SW.MF_CASES]
    e75 = [SW.ex75(n) for n in (10, 100, 1000, 10000)]
    eb = [SW.ex75_big(n) for n in (100, 1000, 10000)]
    rows = [
        r"練習 7.4：$(m,\mathcal N)$ & $(10,10^2)$ & $(10^2,10^4)$ & $(10^3,10^6)$ & & 極限\\",
        r"\quad $\mathcal Nz_c^mC_{MF}(m)$ & " + " & ".join(_s(v, 4) for v in mf) + rf" & & {_s(math.exp(-1), 4)}\\",
        r"\midrule",
        r"練習 7.5：$\mathcal N$ & $10$ & $10^2$ & $10^3$ & $10^4$ & \\",
        r"\quad $\frac z{z_c}=\frac12$：$\mathcal N(1-\frac z{z_c})G_{MF}$ & " + " & ".join(_s(v, 4) for v in e75) + r" & 1\\",
        r"\quad $\frac z{z_c}=\frac32$：$\frac{\ln G_{MF}}{\mathcal N}$ & --- & " + " & ".join(_s(v, 4) for v in eb)
        + r" & 0.0625\\",
    ]
    return "\n".join([r"\begin{tabular}{@{}lccccc@{}}", r"\toprule", *rows, r"\bottomrule", r"\end{tabular}"])


def fig_flory():
    """d = 2、m = 1000、u = 1：−R²/m、−um²/R² 與和；最大在 R* = 177.83。"""
    m, d = 1000, 2
    Rs = np.arange(50, 601, 5.0)
    Rstar = SW.flory_Rstar(m, d)
    Fstar = SW.flory_F(Rstar, m, d)
    out = [r"\begin{tikzpicture}[x=0.0085cm,y=0.022cm,>=Stealth,font=\scriptsize]",
           r"\draw[->] (0,0) -- (640,0) node[right] {$R$};",
           r"\draw[->] (0,-160) -- (0,12);",
           r"\foreach \x in {100,200,300,400,500,600} {\draw (\x,0) -- (\x,3); \node[above] at (\x,2) {\x};}",
           r"\foreach \y in {-150,-100,-50} {\draw (0,\y) -- (-8,\y); \node[left] at (-8,\y) {$\y$};}",
           r"\begin{scope}\clip (0,-160) rectangle (640,0);"]
    for f, sty in ((lambda R: -R * R / m, "cBar,thick"), (lambda R: -m * m / R ** d, "cWalk,thick"),
                   (lambda R: SW.flory_F(R, m, d), "cExact,very thick")):
        pts = " ".join(f"({R:.0f},{max(f(R), -400):.2f})" for R in Rs)
        out.append(rf"\draw[{sty}] plot coordinates {{{pts}}};")
    out.append(r"\end{scope}")
    out.append(rf"\draw[gray,densely dashed] ({Rstar:.2f},-160) -- ({Rstar:.2f},0);")
    out.append(rf"\fill[cExact] ({Rstar:.2f},{Fstar:.2f}) circle (1.6pt);")
    out.append(rf"\fill ({Rstar:.2f},{-Rstar ** 2 / m:.2f}) circle (1.4pt);")
    out.append(rf"\node[below right] at ({Rstar:.2f},-150) {{$R^*={_s(Rstar, 2)}$}};")
    for k, (col, txt) in enumerate((("cWalk", r"排斥 $-um^2/R^2$"), ("cBar", r"熵 $-R^2/m$"),
                                    ("cExact", rf"和，最大 ${_s(Fstar, 2)}$"),
                                    ("black", rf"交點：兩項都是 ${_s(-Rstar ** 2 / m, 2)}$"))):
        out.append(rf"\node[{col},anchor=west] at (320,{-32 - 17 * k}) {{{txt}}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_flory():
    ds = (1, 2, 3, 4)
    known = {1: "1", 2: r"$\frac34$（精確）", 3: r"0.588（數值）", 4: r"$\frac12$"}
    rows = [f"{d} & ${_fr(SW.flory_nu(d))}$ & {_s(SW.flory_slope(d), 4)} & {_s(SW.flory_Rstar(1000, d), 2)} & {known[d]}\\\\"
            for d in ds]
    a = "\n".join([r"\begin{tabular}{@{}ccccl@{}}", r"\toprule", r"$d$ & $\nu_F=\frac3{d+2}$ & 數值斜率 & $R^*$（$m=10^3$） & 已知值\\",
                   r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])
    ns = (6, 8, 10, 12, 14)
    b = "\n".join([r"\begin{tabular}{@{}l*{5}{c}c@{}}", r"\toprule", r"二維 $N$ & " + " & ".join(map(str, ns)) + r" & $\infty$\\",
                   r"\midrule", r"$2\nu_{\rm eff}$ & " + " & ".join(_s(SW.two_nu_eff(n), 4) for n in ns) + r" & $\frac32$\\",
                   r"\bottomrule", r"\end{tabular}"])
    return a + "\n\\\\[6pt]\n" + b


ERRATA = [
    (r"p.~169--173，(7.6)(7.11)、Ex.~7.1", r"$(2\pi Ma^2)^{-3/2}$", r"$(2\pi Ma^2/3)^{-3/2}$"),
    (r"p.~173，(7.11)", r"$e^{-3|\vec x-\vec r_1|^2/3M_1a^2}$、$\cdots/3M_ka^2$", r"分母 $2M_1a^2$、$2M_ka^2$"),
    (r"p.~176，(7.21) 第二、三行", r"$Z^{N-M_1}$、$Z^{K_2}$", r"$z^{N-M_1}$、$z^{K_2}$"),
    (r"p.~178，(7.27) 前", r"Figure 7.3(a)", r"Figure 7.3(b)（並排）"),
    (r"p.~183，(7.33)", r"$g_0[\Sigma_1+\Sigma_2]^2g_0$", r"$g_0[\Sigma_1+\Sigma_2]g_0[\Sigma_1+\Sigma_2]g_0$"),
    (r"p.~183，(7.35)", r"$\frac{g_0}{1+(\Sigma_1+\Sigma_2)g_0}$", r"$\frac{g_0}{1-(\Sigma_1+\Sigma_2)g_0}$"),
    (r"p.~188，(7.47)", r"$za^2\nabla^2G-(z/z_c-1)G=a^d\delta$", r"$(1-z/z_c)G-za^2\nabla^2G=a^d\delta$"),
    (r"p.~188，(7.48)", r"$\psi_i(\vec r)=e^{-\vec q\cdot\vec r}$", r"$\psi_i(\vec r)=e^{i\vec q\cdot\vec r}/\sqrt V$"),
    (r"p.~188，(7.50) 後", r"(7.43dd)", r"(7.43d)"),
    (r"p.~189，(7.55) 前", r"積分 $=1/(z_0-z)$", r"積分 $=z_c/(z_0-z)$"),
    (r"p.~190，(7.58) 與前一句", r"$e^{-\mathcal Nz\lambda/z_c}$；$z^m\cdots\lambda^n$", r"$e^{+\mathcal Nz\lambda/z_c}$；$(z/z_c)^m\cdots\lambda^m$"),
    (r"p.~190，(7.59)", r"$e^{-\mathcal N(f\lambda^2+\mathcal N\lambda)+m\ln\lambda}$", r"$e^{-\mathcal N(f\lambda^2+\lambda)+m\ln\lambda}$"),
    (r"p.~191，(7.61) 左邊", r"$C_{MF}(N;\vec x,\vec y)$", r"$C_{MF}(m;\vec x,\vec y)$"),
]


def tab_errata():
    assert len(ERRATA) == 13
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("SevenFigGrid", fig_grid),
    ("SevenFigGauss", fig_gauss),
    ("SevenTabSplit", tab_split),
    ("SevenFigTwoStep", fig_twostep),
    ("SevenTabSaw", tab_saw),
    ("SevenFigLogLog", fig_loglog),
    ("SevenFigChords", fig_chords),
    ("SevenFigThird", fig_third),
    ("SevenTabIE", tab_ie),
    ("SevenFigDyson", fig_dyson),
    ("SevenFigDysonBars", fig_dyson_bars),
    ("SevenTabDyson", tab_dyson),
    ("SevenTabPower", tab_power),
    ("SevenTabMF", tab_mf),
    ("SevenFigFlory", fig_flory),
    ("SevenTabFlory", tab_flory),
    ("SevenTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch7()
    parts = ["% figs.tex —— 由 ch7/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch7/make_figs.py（第 7 章，動畫 W701–W706）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    check_tex_decimals()                                   # 圖表裡用到的小數也已記進 OK
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
