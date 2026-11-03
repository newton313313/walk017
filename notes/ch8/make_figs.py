r"""
make_figs.py —— 隨機漫步導讀第 8 章（Rudnick & Gaspari《Elements of the Random Walk》第 8 章，動畫 W801–W804）的圖與表

在 notes/ 或 notes/ch8/ 裡執行：
    python ch8/make_figs.py
會先檢查 ch8/part8.tex 寫到的每一個數字（check_numbers_ch8；最後再把 part8.tex 裡出現的每一個小數和檢查過的數比對一遍），
再寫出 ch8/figs.tex（part8.tex 第一行 \input{ch8/figs}）。圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算）。
數學一律由 codes/utils/walk_scale.py 算出（第 7 章的 Flory 指數用 walk_saw）；例子與動畫 W801–W804 相同。巨集名稱都以 \Eight 開頭
（不和前幾章的 \Fig…、\Two…、…、\Seven… 撞名）；part8.tex 的 \label 用 eight-…。顏色沿用導言區：cBar 精確值、cWalk 路線、
cExact 更正、cGauss 近似、cWrong 筆誤。
"""
from __future__ import annotations

import math
import re
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent            # notes/ch8
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_saw as SW  # noqa: E402
from utils import walk_scale as WS  # noqa: E402   （套件模組：walk_scale 的相對 import 才找得到 walk_gf、walk_saw）
from utils import walk_gf as GF  # noqa: E402

OUT = HERE / "figs.tex"
TEX = HERE / "part8.tex"
F = Fraction
OK = set()                                         # 檢查過的小數字串（part8.tex 的小數都要在這裡）


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


def _all(vals, d, expect):
    got = tuple(_s(v, d) for v in vals)
    assert got == tuple(expect), (got, expect)


# =============================================================================
# 數字檢查（part8.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch8():
    OK.update({"0.1", "0.01", "0.001", "1.5", "0.4", "0.45", "0.49", "0.05", "0.96", "0.99", "0.6", "0.588"})
    # ---- W801 長尾巴
    assert WS.pareto_alpha(1) == 4 and WS.delta_of(1) == 1 and WS.delta_of(1.5) == 0.5
    assert [WS.moment2_trunc(1, R) for R in (10, 100, 1000)] == [9, 99, 999]
    assert _s(WS.B_pareto(1), 4) == "0.7854" and _close(WS.B_pareto(1), math.pi / 4, 1e-12)
    _all([WS.chi_ratio(q, 1) for q in (0.1, 0.01, 0.001)], 4, ("0.7687", "0.7837", "0.7852"))
    assert (_s(WS.chi_pareto(0.01, 1), 5), _s(WS.chi_small(0.01, 1), 5)) == ("0.99216", "0.99215")
    assert (_s(WS.chi_slope(1), 2), _s(WS.chi_slope(1.5), 2)) == ("1.00", "1.49")
    assert tuple(f"{WS.median_width(N, 1):.0f}" for N in WS.WIDTH_N) == ("18", "172", "1822")
    _all([WS.median_width(N, 1) / N for N in WS.WIDTH_N], 2, ("1.77", "1.72", "1.82"))
    _all([WS.median_width(N, None) for N in WS.WIDTH_N], 2, ("2.86", "8.98", "28.06"))
    _all([WS.median_width(N, None) / math.sqrt(N) for N in WS.WIDTH_N], 2, ("0.90", "0.90", "0.89"))
    assert (_s(WS.width_pred(1, 1), 4), _s(WS.width_pred(1, None), 4)) == ("1.7785", "0.8881")
    assert _s(GF.cauchy_radial_median(1.0), 4) == "2.2644" and _s(GF.CHI3_MEDIAN, 4) == "1.5382"
    assert (_s(WS.width_slope(1), 2), _s(WS.width_slope(None), 2)) == ("1.01", "0.50")
    _all([WS.levy_density(r, 1, 1) / WS.levy_tail(r, 1, 1) for r in (10, 30, 100)], 4, ("0.9803", "0.9978", "0.9998"))
    _all([WS.levy_density(r, 1, 1.5) / WS.levy_tail(r, 1, 1.5) for r in (10, 30, 100)], 4, ("1.1755", "1.0316", "1.0051"))
    assert (_s(WS.levy_slope(100, 1, 1), 2), _s(WS.levy_slope(100, 1, 1.5), 2)) == ("-4.00", "-4.51")
    assert WS.tail_exponent(1) == -4 and WS.tail_exponent(F(1, 2)) == F(-9, 2)
    assert _close(WS.levy_tail(1, 1, 1), 1 / math.pi ** 2)                     # δ = 1：N/(π² r⁴)
    assert WS.SUBST_Y == 0.5 and _s(math.gamma(1.5), 4) == "0.8862"
    I1, I2 = (abs(WS.subst_integral(r, 0.5, WS.SUBST_EPS)) for r in (1, 2))
    assert _s(I2 / I1, 4) == "0.3536" == _s(2 ** -1.5, 4)
    # ---- W802 標度形式
    assert WS.zc(1) == F(1, 2) and WS.zc(3) == F(1, 6) and WS.AB(3) == (6, 1)
    assert tuple(f"{WS.g_lattice(0.5 - e, 0.0):.0f}" for e in WS.COLLAPSE_EPS) == ("5", "10", "50")
    assert WS.COLLAPSE_EPS == (0.1, 0.05, 0.01) and WS.COLLAPSE_X == (0, 1, 2, 4)
    exp = {0.1: ("0.5000", "0.4172", "0.2819", "0.1317"), 0.05: ("0.5000", "0.4085", "0.2652", "0.1145"),
           0.01: ("0.5000", "0.4017", "0.2529", "0.1027")}
    for e in WS.COLLAPSE_EPS:
        _all([WS.collapse(e, x) for x in WS.COLLAPSE_X], 4, exp[e])
    _all([WS.collapse_limit(x) for x in WS.COLLAPSE_X], 4, ("0.5000", "0.4000", "0.2500", "0.1000"))
    assert all(_close(WS.collapse_limit(x), 2 / (4 + x * x)) for x in WS.COLLAPSE_X)
    OK.update({"0.25", "0.5"})
    assert (_s(WS.xi_88c(0.96), 4), _s(WS.xi_lattice(0.96), 4)) == ("2.0412", "2.0205")
    assert (_s(WS.xi_88c(0.99), 4), _s(WS.xi_lattice(0.99), 4)) == ("4.0825", "4.0722")
    assert _close(WS.xi_88c(0.99) / WS.xi_88c(0.96), 2)
    assert (_s(WS.G_ratio(0.96, 8), 4), _s(WS.G_ratio(0.99, 8), 4)) == ("1.0231", "1.0178")
    assert _s(WS.F3(0), 4) == "0.4775" and _close(WS.F3(0), 3 / (2 * math.pi))
    # ---- W803 最陡下降
    assert WS.SD_N == (20, 100) and WS.SD_X == {20: (4, 8), 100: (8, 12)}
    exp = {20: (("0.6818", "0.6832", "0.4668"), ("0.2098", "0.2179", "0.0475")),
           100: (("0.7282", "0.7285", "0.5306"), ("0.4894", "0.4902", "0.2403"))}
    for N in WS.SD_N:
        for x, e in zip(WS.SD_X[N], exp[N]):
            _all([WS.ratio_exact(N, x), WS.sd_ratio(N, x, 0.5, 2, 0.5), WS.sd_ratio(N, x, 0.5, 2, 0.5, book=True)], 4, e)
    assert _close(WS.sd_K(0.5, 2, 0.5), 0.5) and _close(WS.sd_K(0.5, 2, 0.5, book=True), 1)
    assert _close(WS.sd_K_dim(3), 1.5) and _close(WS.sd_K_dim(1), 0.5)
    # 算到底：N = 20、x = 8：δ = (21/4)^{−2} = 16/441，指數 −64/21 + 32/21 = −32/21
    dl = WS.sd_delta(20, 8, 0.5, 2, 0.5)
    assert _close(dl, 16 / 441) and _close(-2 * 8 * math.sqrt(dl), -64 / 21) and _close(21 * dl / 0.5, 32 / 21)
    assert _s(math.exp(-32 / 21), 4) == "0.2179" and _s(math.exp(-64 / 21), 4) == "0.0475"
    assert [WS.fisher(nu) for nu in WS.FISHER_NU] == [2, F(5, 2), 4]
    assert SW.NU_QUOTED[3] == 0.588 and _s(1 / (1 - 0.588), 2) == "2.43"
    assert WS.gamma_exp(0, F(1, 2)) == 1 and WS.gamma_exp(WS.ETA_2D, WS.NU_2D) == F(43, 32) == WS.GAMMA_QUOTED
    assert WS.ETA_2D == F(5, 24) and WS.NU_2D == F(3, 4) and _s(43 / 32, 5) == "1.34375"
    x2, p2 = WS.growth(WS.GAMMA_QUOTED, 1 / WS.MU_QUOTED)
    assert _s(x2, 3) == "2.638" and p2 == F(11, 32) and _s(p2, 5) == "0.34375"
    assert WS.growth(1, F(1, 6)) == (6, 0)
    assert _s(WS.p_numeric(WS.GAMMA_QUOTED, 1 / WS.MU_QUOTED), 4) == "0.3436"
    assert WS.RATIO_N == (10, 12, 14, 16)
    _all([WS.ratio_method(N)[1] for N in WS.RATIO_N], 4, ("2.6592", "2.6529", "2.6483", "2.6454"))
    _all([WS.ratio_method(N)[2] for N in WS.RATIO_N], 4, ("1.3552", "1.3519", "1.3496", "1.3483"))
    assert WS.theta(WS.GAMMA_QUOTED, WS.NU_2D) == F(11, 24) and WS.theta(1, F(1, 2)) == 0
    assert WS.ex85_s(2, F(3, 4)) == F(2, 3)
    _all([WS.return_slope(d, 200) for d in (1, 2, 3)], 2, ("-0.50", "-1.00", "-1.50"))
    # ---- W804 重新標度
    assert WS.kin_factor(F(1, 2), F(1, 4)) == 1 and WS.int_factor(F(1, 2), F(1, 4), 3) == 2
    assert WS.int_factor_book86(F(1, 2), F(1, 4), 3) == F(1, 2)
    k, v = WS.rescale_check(0.5, 0.25, 3)
    assert (_s(k, 2), _s(v, 2)) == ("1.00", "2.00")
    assert [WS.sa_amplitude(100, d) for d in range(1, 6)] == [1000, 100, 10, 1, F(1, 10)] and WS.UPPER_CRITICAL_D == 4
    assert all(_close(WS.int_factor(math.sqrt(1 / 100), 1 / 100, d), 100 ** ((4 - d) / 2)) for d in range(1, 6))
    assert (_s(WS.flory_l1(0.01, 3), 4), _s(1 / WS.flory_l1(0.01, 3), 2)) == ("0.0631", "15.85")
    assert _close(1 / WS.flory_l1(0.01, 3), 100 ** 0.6) and _close(WS.flory_l1_num(0.01, 3), WS.flory_l1(0.01, 3), 1e-9)
    assert [WS.nu_rescale(d) for d in range(1, 5)] == [1, F(3, 4), F(3, 5), F(1, 2)] == [SW.flory_nu(d) for d in range(1, 5)]
    assert WS.ex88_relevant(3, 2) and not WS.ex88_relevant(5, 5) and WS.ex88_relevant(5, 3) and WS.ex88_relevant(3, 4)
    assert WS.ex88_nu(2) == F(3, 4) and _close(WS.ex88_amplitude(100, 2), 100)


def check_tex_decimals():
    """part8.tex（去掉註解）裡的每一個小數都必須是上面檢查過的數。"""
    if not TEX.exists():
        print("（尚無 part8.tex，略過小數比對）")
        return
    txt = re.sub(r"(?<!\\)%.*", "", TEX.read_text(encoding="utf-8"))
    txt = re.sub(r"[\d.]+\\linewidth|[\d.]+pt\b|[\d.]+em\b|\\arraystretch\}\{[\d.]+\}|\(2\.(20|31)\)|\(7\.16\)", "", txt)   # 版面尺寸、第 2、7 章的式號
    eqno = re.compile(r"8\.\d{1,2}")                     # 本章的式號、節號、圖號、練習號（8.1–8.25）
    bad = sorted({m for m in re.findall(r"\d+\.\d+", txt) if m not in OK and not eqno.fullmatch(m)})
    assert not bad, f"part8.tex 有未檢查的小數：{bad}"


# =============================================================================
# 圖與表
# =============================================================================
def fig_chi():
    """log₁₀(1 − χ(q)) 對 log₁₀ q：δ = 1（y = 1）、δ = 1/2（y = 1.5）的數值點，與步長 1 的高斯步伐（斜率 2）。"""
    lq = np.linspace(-3, -1, 9)
    out = [r"\begin{tikzpicture}[x=2.1cm,y=0.68cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[gray!25] (-3,-5) grid (-1,0);",
           r"\draw[->] (-3,-5) -- (-0.8,-5) node[right] {$q$};",
           r"\draw[->] (-3,-5) -- (-3,0.35) node[above] {$1-\chi(q)$};",
           r"\foreach \x in {-3,-2,-1} {\node[below] at (\x,-5) {$10^{\x}$};}",
           r"\foreach \y in {-5,-4,-3,-2,-1,0} {\node[left] at (-3,\y) {$10^{\y}$};}"]
    for y, col, mark, lab in ((1, "cBar", "circle (1.5pt)", rf"$\delta=1$：斜率 {_s(WS.chi_slope(1), 2)}"),
                              (1.5, "cWalk", "circle (1.5pt)", rf"$\delta=1/2$：斜率 {_s(WS.chi_slope(1.5), 2)}")):
        pts = [(x, math.log10(1 - WS.chi_pareto(10 ** x, y))) for x in lq]
        out.append(rf"\draw[{col},thin] plot coordinates {{{' '.join(f'({a:.3f},{b:.3f})' for a, b in pts)}}};")
        out += [rf"\fill[{col}] ({a:.3f},{b:.3f}) {mark};" for a, b in pts]
        yl = -0.15 if y == 1 else -0.6
        out.append(rf"\fill[{col}] (-2.88,{yl}) circle (1.5pt); \node[{col},anchor=west] at (-2.84,{yl}) {{{lab}}};")
    g = [(x, math.log10(1 - math.sin(10 ** x) / 10 ** x)) for x in np.linspace(-2.1, -1, 12)]
    out.append(rf"\draw[cExact,densely dashed,thick] plot coordinates {{{' '.join(f'({a:.3f},{b:.3f})' for a, b in g)}}};")
    out.append(r"\draw[cExact,densely dashed,thick] (-2.94,-1.05) -- (-2.82,-1.05) node[anchor=west] {步長 1：斜率 2};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_width():
    Ns = WS.WIDTH_N
    rows = [rf"{N} & {WS.median_width(N, 1):.0f} & {_s(WS.median_width(N, 1) / N, 2)} & {_s(WS.median_width(N, None), 2)} & "
            rf"{_s(WS.median_width(N, None) / math.sqrt(N), 2)}\\" for N in Ns]
    rows.append(rf"理論 & & {_s(WS.width_pred(1, 1), 4)} & & {_s(WS.width_pred(1, None), 4)}\\")
    rows.append(rf"斜率 & {_s(WS.width_slope(1), 2)} & & {_s(WS.width_slope(None), 2)} & \\")
    return "\n".join([r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule",
                      r"$N$ & 長尾巴 & $\div N$ & 步長 1 & $\div\sqrt N$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_tail():
    rs = (10, 30, 100)
    r1 = r"$\delta=1$ & " + " & ".join(_s(WS.levy_density(r, 1, 1) / WS.levy_tail(r, 1, 1), 4) for r in rs) + \
         rf" & ${_s(WS.levy_slope(100, 1, 1), 2)}$\\"
    r2 = r"$\delta=\frac12$ & " + " & ".join(_s(WS.levy_density(r, 1, 1.5) / WS.levy_tail(r, 1, 1.5), 4) for r in rs) + \
         rf" & ${_s(WS.levy_slope(100, 1, 1.5), 2)}$\\"
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$r$ & 10 & 30 & 100 & 斜率\\", r"\midrule", r1, r2, r"\bottomrule", r"\end{tabular}"])


def fig_rawg():
    """一維 g(z;q) = 1/(1 − 2z cos q)，z = 0.4, 0.45, 0.49。"""
    qs = np.linspace(-1, 1, 81)
    out = [r"\begin{tikzpicture}[x=1.6cm,y=0.068cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (-1.1,0) -- (1.2,0) node[right] {$q$};",
           r"\draw[->] (0,0) -- (0,55) node[above] {$g(z;q)$};",
           r"\foreach \x in {-1,1} {\draw (\x,0) -- (\x,-1.2); \node[below] at (\x,-1) {$\x$};}",
           r"\foreach \y in {10,30,50} {\draw (0,\y) -- (-0.04,\y); \node[left] at (-0.04,\y) {$\y$};}"]
    for e, col in zip(WS.COLLAPSE_EPS, ("cBar!40", "cBar!70", "cBar")):
        z = 0.5 - e
        pts = " ".join(f"({q:.3f},{WS.g_lattice(z, q):.3f})" for q in qs)
        out.append(rf"\draw[{col},thick] plot coordinates {{{pts}}};")
    out.append(r"\node[cBar,anchor=west] at (0.06,50) {$z=0.49$};")
    out.append(r"\node[cBar!70,anchor=west] at (0.42,12.5) {$0.45$};")
    out.append(r"\node[cBar!40,anchor=west] at (0.72,7.5) {$0.4$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_collapse():
    """塌縮：(z_c − z) g 對 x = q/√(z_c − z)，三個 z 與極限 2/(4 + x²)。"""
    xs = np.linspace(-4, 4, 81)
    out = [r"\begin{tikzpicture}[x=0.42cm,y=6.8cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (-4.3,0) -- (4.6,0) node[right] {$x$};",
           r"\draw[->] (0,0) -- (0,0.58) node[above] {$(z_c-z)\,g$};",
           r"\foreach \x in {-4,-2,2,4} {\draw (\x,0) -- (\x,-0.012); \node[below] at (\x,-0.01) {$\x$};}",
           r"\foreach \y in {0.1,0.25,0.5} {\draw (0,\y) -- (-0.12,\y); \node[left] at (-0.12,\y) {$\y$};}"]
    for e, col in zip(WS.COLLAPSE_EPS, ("cBar!40", "cBar!70", "cBar")):
        pts = " ".join(f"({x:.3f},{WS.collapse(e, x):.4f})" for x in xs)
        out.append(rf"\draw[{col},thick] plot coordinates {{{pts}}};")
    pts = " ".join(f"({x:.3f},{WS.collapse_limit(x):.4f})" for x in xs)
    out.append(rf"\draw[cExact,densely dashed,thick] plot coordinates {{{pts}}};")
    out.append(r"\node[cExact,anchor=west] at (1.6,0.36) {極限 $\frac{2}{4+x^2}$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_collapse():
    head = r"$z_c-z$ & $z$ & $g(z;0)$ & " + " & ".join(f"$x={x}$" for x in WS.COLLAPSE_X) + r"\\"
    zs = {0.1: "0.4", 0.05: "0.45", 0.01: "0.49"}
    rows = [rf"{e} & {zs[e]} & {WS.g_lattice(0.5 - e, 0.0):.0f} & " +
            " & ".join(_s(WS.collapse(e, x), 4) for x in WS.COLLAPSE_X) + r"\\" for e in WS.COLLAPSE_EPS]
    lim = r"$\to0$ & & & " + " & ".join(_s(WS.collapse_limit(x), 4) for x in WS.COLLAPSE_X) + r"\\"
    return "\n".join([r"\begin{tabular}{@{}lll*{4}{c}@{}}", r"\toprule", head, r"\midrule", *rows, r"\midrule",
                      lim, r"\bottomrule", r"\end{tabular}"])


def tab_xi():
    rows = [rf"{x} & {_s(WS.xi_88c(x), 4)} & {_s(WS.xi_lattice(x), 4)} & {_s(WS.G_ratio(x, 8), 4)}\\" for x in WS.XI_CASES]
    return "\n".join([r"\begin{tabular}{@{}cccc@{}}", r"\toprule",
                      r"$z/z_c$ & $\xi$（(8.8c)） & $\xi$（格子） & $G_{\text{格子}}/G_{\rm Yukawa}$（$r=8$）\\", r"\midrule",
                      *rows, r"\bottomrule", r"\end{tabular}"])


def fig_sd():
    """一維 N = 20：C(20;x)/C(20;0) 的長條、更正的 e^{−x²/42}（綠）、書上的 e^{−x²/21}（紅虛線）。"""
    N = 20
    out = [r"\begin{tikzpicture}[x=0.16cm,y=2.6cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (-21.5,0) -- (22.5,0) node[right] {$x$};",
           r"\draw[->] (-21.5,0) -- (-21.5,1.15);",
           r"\foreach \x in {-20,-10,10,20} {\draw (\x,0) -- (\x,-0.03); \node[below] at (\x,-0.02) {$\x$};}",
           r"\foreach \y in {0.5,1} {\draw (-21.5,\y) -- (-21.9,\y); \node[left] at (-21.9,\y) {$\y$};}"]
    for x in range(-N, N + 1, 2):
        h = WS.ratio_exact(N, x)
        out.append(rf"\fill[cBar!35] ({x - 0.8},0) rectangle ({x + 0.8},{h:.4f});")
    xs = np.linspace(-20, 20, 121)
    for book, sty in ((False, "cExact,very thick"), (True, "cWrong,densely dashed,thick")):
        pts = " ".join(f"({x:.2f},{WS.sd_ratio(N, abs(x), 0.5, 2, 0.5, book=book):.4f})" for x in xs)
        out.append(rf"\draw[{sty}] plot coordinates {{{pts}}};")
    out.append(r"\node[cExact,anchor=west] at (8.5,0.85) {更正 $e^{-x^2/42}$};")
    out.append(r"\node[cWrong,anchor=west] at (8.5,0.7) {書上 $e^{-x^2/21}$};")
    out.append(r"\node[cBar,anchor=west] at (8.5,0.55) {精確（長條）};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_sd():
    rows = []
    for N in WS.SD_N:
        for x in WS.SD_X[N]:
            rows.append(rf"{N} & {x} & {_s(WS.ratio_exact(N, x), 4)} & {_s(WS.sd_ratio(N, x, 0.5, 2, 0.5), 4)} & "
                        rf"{_s(WS.sd_ratio(N, x, 0.5, 2, 0.5, book=True), 4)}\\")
    return "\n".join([r"\begin{tabular}{@{}rrccc@{}}", r"\toprule",
                      r"$N$ & $x$ & 精確 & 更正 $e^{-x^2/2(N+1)}$ & 書上 $e^{-x^2/(N+1)}$\\", r"\midrule",
                      *rows, r"\bottomrule", r"\end{tabular}"])


def tab_fisher():
    nus = list(WS.FISHER_NU)
    cells = [f"${_fr(nu)}$" for nu in nus] + [r"0.588"]
    outs = [f"${_fr(WS.fisher(nu))}$" for nu in nus] + [_s(1 / (1 - SW.NU_QUOTED[3]), 2)]
    what = ["高斯", "三維 Flory", "二維自迴避", "三維自迴避"]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$\nu$ & " + " & ".join(cells) + r"\\",
                      r"$\frac1{1-\nu}$ & " + " & ".join(outs) + r"\\",
                      r"\midrule", r" & " + " & ".join(what) + r"\\", r"\bottomrule", r"\end{tabular}"])


def tab_ratio():
    Ns = WS.RATIO_N
    return "\n".join([r"\begin{tabular}{@{}l*{4}{c}c@{}}", r"\toprule",
                      r"$N$ & " + " & ".join(map(str, Ns)) + r" & 已知值\\", r"\midrule",
                      r"$\hat\mu_N$ & " + " & ".join(_s(WS.ratio_method(N)[1], 4) for N in Ns) + rf" & {_s(WS.MU_QUOTED, 3)}\\",
                      r"$\hat\gamma_N$ & " + " & ".join(_s(WS.ratio_method(N)[2], 4) for N in Ns) + rf" & ${_fr(WS.GAMMA_QUOTED)}$\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_amp():
    ds = range(1, 6)
    amp = []
    for d in ds:
        a = WS.sa_amplitude(100, d)
        amp.append(str(a) if F(a).denominator == 1 else _s(a, 1))
    return "\n".join([r"\begin{tabular}{@{}l*{5}{c}@{}}", r"\toprule",
                      r"$d$ & " + " & ".join(map(str, ds)) + r"\\", r"\midrule",
                      r"$\frac{4-d}{2}$ & " + " & ".join(f"${_fr(F(4 - d, 2))}$" for d in ds) + r"\\",
                      r"$N^{(4-d)/2}$（$N=100$） & " + " & ".join(amp) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_nu():
    known = {1: "1", 2: r"$\frac34$（精確）", 3: r"0.588（數值）", 4: r"$\frac12$"}
    rows = [rf"{d} & ${_fr(WS.nu_rescale(d))}$ & ${_fr(SW.flory_nu(d))}$ & {known[d]}\\" for d in range(1, 5)]
    return "\n".join([r"\begin{tabular}{@{}cccl@{}}", r"\toprule",
                      r"$d$ & 放大鏡 $\frac3{d+2}$ & 第 7 章拔河 $\nu_F$ & 已知值\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    (r"p.~195，(8.6)$\to$(8.7)", r"replace $k$ by $w=k/r$", r"$w=kr$（$k=w/r$）"),
    (r"p.~196，第一段", r"$z_c$ is the coordination number", r"$z_c=1/(\text{配位數})=\frac1{2d}$"),
    (r"p.~198，(8.18)", r"$\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}z=0$", r"$-\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}z=0$"),
    (r"p.~198，(8.20) 第一行", r"$\exp\big(-B|\vec r|\delta^\nu-\frac{N+1}{z_c}\big)$",
     r"$z_c^{-(N+1)}\exp\big(-B|\vec r|\delta^\nu+\frac{(N+1)\delta}{z_c}\big)$"),
    (r"p.~198，(8.20) 第二行", r"$\exp[-((\nu z_c)^\nu B)^{\frac1{1-\nu}}\cdots]$", r"$\exp[-(1-\nu)((\nu z_c)^\nu B)^{\frac1{1-\nu}}\cdots]$"),
    (r"p.~200，(8.22) 前", r"$\vec r\,'(s')=l_1\vec r(l_2s)$", r"$\vec r\,'(s')=l_1\vec r(s'/l_2)$"),
    (r"p.~200，Ex.~8.6", r"rescaled by $l_1^d/l_2$", r"$l_1^d/l_2^2$"),
]


def tab_errata():
    assert len(ERRATA) == 7
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("EightFigChi", fig_chi),
    ("EightTabWidth", tab_width),
    ("EightTabTail", tab_tail),
    ("EightFigRawG", fig_rawg),
    ("EightFigCollapse", fig_collapse),
    ("EightTabCollapse", tab_collapse),
    ("EightTabXi", tab_xi),
    ("EightFigSD", fig_sd),
    ("EightTabSD", tab_sd),
    ("EightTabFisher", tab_fisher),
    ("EightTabRatio", tab_ratio),
    ("EightTabAmp", tab_amp),
    ("EightTabNu", tab_nu),
    ("EightTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch8()
    parts = ["% figs.tex —— 由 ch8/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch8/make_figs.py（第 8 章，動畫 W801–W804）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    check_tex_decimals()                                   # 圖表裡用到的小數也已記進 OK
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
