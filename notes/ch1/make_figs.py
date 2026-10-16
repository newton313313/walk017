r"""
make_figs.py —— 隨機漫步導讀第 1 章（Rudnick & Gaspari《Elements of the Random Walk》第 1 章，動畫 W101–W112）的圖與表

在 notes/ 或 notes/ch1/ 裡執行：
    python ch1/make_figs.py
會先檢查 ch1/part1.tex 寫到的每一個數字（check_numbers_ch1），再寫出 ch1/figs.tex（part1.tex 第一行 \input{ch1/figs}）。
圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots）。數學一律由 codes/utils/walk1d.py 算出；例子與動畫 W101–W112 相同。
顏色（導言區定義）：cHead 正面／向右、cTail 反面／向左、cWalk 主例路線、cBar 精確值長條、cGauss 高斯近似、cExact 精確、cWrong 筆誤。
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent            # notes/ch1
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk1d as W  # noqa: E402

OUT = HERE / "figs.tex"


def f4(x) -> str:
    return f"{float(x):.4f}"


def sci(x: float, digits: int = 4) -> str:
    """科學記號（TeX）：2.4329\times10^{18}。"""
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    return f"{m:.{digits}f}\\times10^{{{e}}}"


def thin(n: int) -> str:
    """整數加千分位細空白：3\,628\,800。"""
    s = f"{n:,}"
    return s.replace(",", r"\,")


# =============================================================================
# 數字檢查（part1.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch1():
    W.check()                                                     # 模組本身的全部檢查（含書上的數字與場景卡）
    # 1.1
    assert W.MAIN_COINS == "HTHHTHTTHH" and W.MAIN_POS == [0, 1, 0, 1, 2, 1, 2, 1, 0, 1, 2]
    assert [W.sequences_with_heads(3, n) for n in (3, 2, 1, 0)] == [["HHH"], ["HHT", "HTH", "THH"],
                                                                     ["HTT", "THT", "TTH"], ["TTT"]]
    assert W.walk_count(5, 1) == 10 and W.walk_prob(5, 1) == Fraction(5, 16)
    assert W.walk_count(6, 0) == 20 and W.walk_prob(6, 0) == Fraction(5, 16)
    assert max(W.distribution(6).items(), key=lambda kv: kv[1]) == (0, 20)
    assert [W.biased_prob(3, n, Fraction(2, 3)) for n in range(4)] == [Fraction(k, 27) for k in (1, 6, 12, 8)]
    assert W.biased_prob(10, 7, Fraction(2, 3)) == Fraction(120 * 2 ** 7, 3 ** 10)
    assert round(float(W.biased_prob(10, 7, Fraction(2, 3))), 4) == 0.2601
    assert max(range(11), key=lambda n: W.biased_prob(10, n, Fraction(2, 3))) == 7
    # 1.2
    assert W.mean(W.finals(W.ensemble(20, 10, 1))) == Fraction(1, 10)
    assert W.moment_exact(10, 6) == 12160 == W.moment_formula(10, 6)
    assert [W.moment_formula(N, 6) for N in (1, 2, 3)] == [1, 32, 183]
    # 1.3
    assert len(str(math.factorial(100))) == 158
    assert math.factorial(20) == 2432902008176640000 and sci(W.stirling(20)) == "2.4228\\times10^{18}"
    assert 1 / 384 == (1 / 2) ** 4 / 24 and 1 / 8 == (1 / 2) ** 2 / 2          # cosh(y/2) 的展開係數
    for N in range(1, 8):                                                      # cosh(y/2)^N 的 y^4 係數 → (1.26)
        c4 = Fraction(N, 384) + Fraction(N * (N - 1), 128)
        assert 24 * c4 * 16 == W.moment_formula(N, 4) and 2 * Fraction(N, 8) * 4 == N
    assert round(W.stirling_ratio(20), 5) == 1.00418
    c100 = math.comb(100, 50)
    approx100 = 2 ** 100 * math.sqrt(2 / (math.pi * 100))
    assert sci(c100) == "1.0089\\times10^{29}" and sci(approx100) == "1.0114\\times10^{29}"
    assert round(approx100 / c100, 4) == 1.0025
    # 1.4
    assert all((y - 1) ** 3 + (y + 1) ** 3 - 2 * y ** 3 == 6 * y for y in range(-5, 6))
    assert all((y - 1) ** 4 + (y + 1) ** 4 - 2 * y ** 4 == 12 * y * y + 2 for y in range(-5, 6))
    # 1.5：回到原點的生成函數 G(z,0) = 1/sqrt(1-4z^2)
    assert [math.comb(2 * k, k) for k in range(6)] == [1, 2, 6, 20, 70, 252]
    z = 0.3
    series = sum(W.walk_count(N, 0) * z ** N for N in range(0, 200))
    assert abs(series - 1 / math.sqrt(1 - 4 * z * z)) < 1e-12
    integral = sum(1 / (1 - 2 * z * math.cos(-math.pi + 2 * math.pi * k / 400)) for k in range(400)) / 400
    assert abs(integral - 1 / math.sqrt(1 - 4 * z * z)) < 1e-12
    assert abs(W.cauchy_coeff(lambda u: 1 / (1 - 4 * u * u) ** 0.5, 10, 0.3).real - 252) < 1e-6
    # 1.6
    assert round(W.sd_approx(10, 3), 2) == 62.2 and abs(W.sd_exact(10, 3) - 66) < 1e-9
    assert round(W.sd_approx(10, 3) / W.sd_exact(10, 3), 4) == 0.9424
    assert round(W.stirling_ratio(3), 4) == 1.0281
    assert round(W.sd_approx(10 ** 6, 2) / W.sd_exact(10 ** 6, 2), 4) == 1.0422
    return True


# =============================================================================
# 圖與表
# =============================================================================
def fig_main_walk():
    pos = W.MAIN_POS
    out = [r"\begin{tikzpicture}[x=0.82cm,y=0.42cm,>=Stealth]",
           r"\draw[gray!25] (0,-3) grid (10,3);",
           r"\draw[->] (0,0) -- (10.7,0) node[right]{\small 步數 $k$};",
           r"\draw[->] (0,-3.3) -- (0,3.9) node[above]{\small 位置 $x$};"]
    for k in range(2, 11, 2):
        out.append(rf"\node[below,font=\scriptsize] at ({k},0) {{${k}$}};")
    for y in (-2, 2):
        out.append(rf"\node[left,font=\scriptsize] at (0,{y}) {{${y}$}};")
    coin = {"H": r"\coinH", "T": r"\coinT"}
    for k, c in enumerate(W.MAIN_COINS, 1):
        out.append(rf"\node at ({k - 0.5},3.45) {{{coin[c]}}};")
    pts = " -- ".join(f"({k},{x})" for k, x in enumerate(pos))
    out.append(rf"\draw[cWalk,very thick] {pts};")
    for k, x in enumerate(pos):
        out.append(rf"\fill[cWalk] ({k},{x}) circle (2.2pt);")
    out.append(rf"\node[right,font=\small] at (10,{pos[-1]}) {{$d={pos[-1]}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def tab_three_coins():
    coin = {"H": r"\coinH", "T": r"\coinT"}
    rows = [r"\begin{tabular}{cccc}", r"\toprule", r"正面數 $n$ & 序列 & 個數 $\binom3n$ & 終點 $d=2n-3$\\", r"\midrule"]
    for n in (3, 2, 1, 0):
        seqs = W.sequences_with_heads(3, n)
        cells = r"\quad ".join("".join(coin[c] for c in s) for s in seqs)
        rows.append(rf"${n}$ & {cells} & ${len(seqs)}$ & ${2 * n - 3}$\\")
    rows += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(rows) + "\n"


def tab_ten_dist():
    ds = list(range(0, 11, 2))
    head = " & ".join(f"${d}$" for d in ds)
    c = " & ".join(f"${W.walk_count(10, d)}$" for d in ds)
    p = " & ".join(f"${f4(W.walk_prob(10, d))}$" for d in ds)
    g = " & ".join(f"${f4(W.gauss_prob(10, d))}$" for d in ds)
    rows = [r"\begin{tabular}{l" + "c" * len(ds) + "}", r"\toprule",
            rf"終點 $d$ & {head}\\", r"\midrule",
            rf"走法數 $C(10,d)$ & {c}\\",
            rf"機率 $P(10,d)=C/1024$ & {p}\\",
            rf"高斯近似 $\frac{{2}}{{\sqrt{{20\pi}}}}e^{{-d^2/20}}$ & {g}\\",
            r"\bottomrule", r"\end{tabular}"]
    return "\n".join(rows) + "\n"


def fig_delta_grid():
    n = 4
    out = [r"\begin{tikzpicture}[x=1.05cm,y=1.05cm]"]
    for j in range(n):
        for k in range(n):
            x, y = k, -j
            if j == k:
                out.append(rf"\filldraw[fill=cHead!35,draw=cHead!80!black] ({x},{y}) rectangle ++(0.95,0.95);")
                out.append(rf"\node[font=\scriptsize] at ({x + 0.475},{y + 0.475}) {{$\Delta_{j + 1}^2=l^2$}};")
            else:
                out.append(rf"\filldraw[fill=gray!12,draw=gray!60] ({x},{y}) rectangle ++(0.95,0.95);")
                out.append(rf"\node[font=\scriptsize] at ({x + 0.475},{y + 0.475}) {{$\Delta_{j + 1}\Delta_{k + 1}$}};")
    for k in range(n):
        out.append(rf"\node[font=\scriptsize,gray] at ({k + 0.475},1.2) {{$k={k + 1}$}};")
        out.append(rf"\node[font=\scriptsize,gray] at (-0.45,{-k + 0.475}) {{$j={k + 1}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def tab_moments():
    Ns = range(1, 7)
    rows = [r"\begin{tabular}{l" + "c" * 6 + "}", r"\toprule",
            "$N$ & " + " & ".join(f"${N}$" for N in Ns) + r"\\", r"\midrule",
            r"$\overline{x^2}=N$ & " + " & ".join(f"${W.moment_exact(N, 2)}$" for N in Ns) + r"\\",
            r"$\overline{x^4}=3N^2-2N$ & " + " & ".join(f"${W.moment_exact(N, 4)}$" for N in Ns) + r"\\",
            r"$\overline{x^6}=15N^3-30N^2+16N$ & " + " & ".join(f"${W.moment_exact(N, 6)}$" for N in Ns) + r"\\",
            r"$\overline{x^4}/(\overline{x^2})^2=3-2/N$ & "
            + " & ".join(f"${_frac(Fraction(W.moment_formula(N, 4), N * N))}$" for N in Ns) + r"\\",
            r"\bottomrule", r"\end{tabular}"]
    return "\n".join(rows) + "\n"


def _frac(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else rf"\tfrac{{{q.numerator}}}{{{q.denominator}}}"


def fig_kurtosis():
    out = [r"\begin{tikzpicture}[x=0.62cm,y=1.1cm,>=Stealth]",
           r"\draw[->] (0,0) -- (10.8,0) node[right]{\small $N$};",
           r"\draw[->] (0,0) -- (0,3.5);",
           r"\draw[cGauss,dashed,thick] (0,3) -- (10.5,3) node[right,font=\small]{$3$（高斯）};"]
    for N in range(1, 11):
        r = 3 - 2 / N
        out.append(rf"\fill[cBar] ({N},{r:.4f}) circle (2.4pt);")
        out.append(rf"\node[below,font=\scriptsize] at ({N},0) {{${N}$}};")
    for y in (1, 2):
        out.append(rf"\draw (-0.08,{y}) -- (0.08,{y}) node[left=2pt,font=\scriptsize] {{${y}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def fig_gamma():
    sx, sy = 0.3, 0.75                     # cm per unit（w；10^5）
    ws = [0.25 * k for k in range(0, 121)]
    f = " ".join(f"({w:.2f},{W.gamma_integrand(w, 10) / 1e5:.4f})" for w in ws)
    g = " ".join(f"({w:.2f},{W.gamma_gauss(w, 10) / 1e5:.4f})" for w in ws)
    peak = W.gamma_integrand(10, 10)
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\draw[->] (0,0) -- (31.5,0) node[right]{\small $w$};",
           r"\draw[->] (0,0) -- (0,5.3);"]
    for w in (5, 10, 15, 20, 25, 30):
        out.append(rf"\draw ({w},0.06) -- ({w},-0.06) node[below,font=\scriptsize]{{${w}$}};")
    for y in (1, 2, 3, 4):
        out.append(rf"\draw (0.3,{y}) -- (-0.3,{y}) node[left,font=\scriptsize]{{${y}\times10^5$}};")
    out.append(rf"\draw[cBar,thick] plot coordinates {{{f}}};")
    out.append(rf"\draw[cGauss,thick,dashed] plot coordinates {{{g}}};")
    out.append(rf"\node[above right,font=\scriptsize] at (10,{peak / 1e5:.3f}) {{峰高 $10^{{10}}e^{{-10}}\approx{thin(round(peak))}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def tab_stirling():
    rows = [r"\begin{tabular}{rccc}", r"\toprule",
            r"$n$ & $n!$ & $\sqrt{2\pi n}\,n^ne^{-n}$ & 比值 $n!/\text{Stirling}$\\", r"\midrule"]
    for n in (1, 2, 5, 10, 20, 100):
        fact = math.factorial(n)
        st = W.stirling(n) if n < 150 else None
        fs = thin(fact) if n <= 10 else sci(float(fact))
        ss = f"{st:.4f}" if n <= 2 else (f"{st:.3f}" if n == 5 else (thin(round(st)) if n == 10 else sci(st)))
        rows.append(rf"${n}$ & ${fs}$ & ${ss}$ & ${W.stirling_ratio(n):.4f}$\\")
    rows += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(rows) + "\n"


def fig_binom_gauss():
    sx, sy = 0.36, 15.0
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\draw[->] (-11,0) -- (11.4,0) node[right]{\small $d$};"]
    for d in range(-10, 11, 2):
        p = float(W.walk_prob(10, d))
        out.append(rf"\filldraw[fill=cBar!45,draw=cBar] ({d - 0.4},0) rectangle ({d + 0.4},{p:.4f});")
        out.append(rf"\node[below,font=\scriptsize] at ({d},0) {{${d}$}};")
    xs = [x / 4 for x in range(-44, 45)]
    pts = " ".join(f"({x:.2f},{W.gauss_prob(10, x):.4f})" for x in xs)
    out.append(rf"\draw[cGauss,thick] plot[smooth] coordinates {{{pts}}};")
    out.append(r"\draw (-11,0.25) -- (-10.8,0.25) node[left,font=\scriptsize] {$0.25$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def fig_pascal():
    rows = W.pascal_rows(8)
    out = [r"\begin{tikzpicture}[x=0.5cm,y=0.62cm,>=Stealth]"]
    for N, row in enumerate(rows):
        for y, c in row.items():
            style = "font=\\small"
            if (N, y) == (6, 0):
                style = r"font=\small\bfseries,text=cWrong"
            elif N == 5 and y in (-1, 1):
                style = r"font=\small,text=cBar"
            out.append(rf"\node[{style}] (n{N}_{y + 10}) at ({y},{-N}) {{${c}$}};")
        out.append(rf"\node[font=\scriptsize,gray] at (-10.5,{-N}) {{$N={N}$}};")
    out.append(r"\draw[->,cBar,thick] (n5_9) -- (n6_10);")
    out.append(r"\draw[->,cBar,thick] (n5_11) -- (n6_10);")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def fig_diffusion():
    N = 30
    sx, sy = 0.2, 25.0
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\draw[->] (-21,0) -- (21.5,0) node[right]{\small $y$};"]
    for y in range(-20, 21, 2):
        p = float(W.walk_prob(N, y))
        out.append(rf"\filldraw[fill=cBar!45,draw=cBar] ({y - 0.75},0) rectangle ({y + 0.75},{p:.4f});")
    for y in range(-20, 21, 10):
        out.append(rf"\node[below,font=\scriptsize] at ({y},0) {{${y}$}};")
    xs = [x / 2 for x in range(-42, 43)]
    pts = " ".join(f"({x:.1f},{2 * W.heat_kernel(x, N):.4f})" for x in xs)
    out.append(rf"\draw[cGauss,thick] plot[smooth] coordinates {{{pts}}};")
    out.append(r"\draw (-21,0.1) -- (-20.6,0.1) node[left,font=\scriptsize] {$0.1$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out) + "\n"


def fig_chi_power():
    sx, sy = 1.25, 2.6
    qs = [-math.pi + 2 * math.pi * k / 240 for k in range(241)]
    pts = " ".join(f"({q:.4f},{math.cos(q) ** 10:.4f})" for q in qs)
    g0 = " ".join(f"({q:.4f},{math.exp(-5 * q * q):.4f})" for q in qs if abs(q) <= 1.3)
    gl = " ".join(f"({q:.4f},{math.exp(-5 * (q + math.pi) ** 2):.4f})" for q in qs if q <= -math.pi + 1.3)
    gr = " ".join(f"({q:.4f},{math.exp(-5 * (q - math.pi) ** 2):.4f})" for q in qs if q >= math.pi - 1.3)
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\draw[->] (-3.4,0) -- (3.6,0) node[right]{\small $q$};",
           r"\draw[->] (0,0) -- (0,1.25);",
           rf"\draw[cBar,thick] plot coordinates {{{pts}}};",
           rf"\draw[cGauss,thick,dashed] plot coordinates {{{g0}}};",
           rf"\draw[cGauss,thick,dashed] plot coordinates {{{gl}}};",
           rf"\draw[cGauss,thick,dashed] plot coordinates {{{gr}}};",
           r"\node[below,font=\scriptsize] at (-3.1416,0) {$-\pi$};",
           r"\node[below,font=\scriptsize] at (3.1416,0) {$\pi$};",
           r"\node[below left,font=\scriptsize] at (0,0) {$0$};",
           r"\draw (-0.06,1) -- (0.06,1) node[right,font=\scriptsize] {$1$};",
           r"\end{tikzpicture}"]
    return "\n".join(out) + "\n"


def fig_contour():
    s = 4.0
    out = [rf"\begin{{tikzpicture}}[x={s}cm,y={s}cm,>=Stealth]",
           r"\draw[->] (-0.45,0) -- (0.75,0) node[right]{\small $\mathrm{Re}\,z$};",
           r"\draw[->] (0,-0.45) -- (0,0.5) node[above]{\small $\mathrm{Im}\,z$};",
           r"\draw[cExact,thick,postaction={decorate},decoration={markings,mark=at position 0.15 with {\arrow{Stealth}}}]"
           r" (0,0) circle (0.25);",
           r"\node[cExact,font=\scriptsize] at (-0.2,0.28) {$|z|=\tfrac14$};",
           r"\draw[cWrong,thick] (0.47,-0.03) -- (0.53,0.03) (0.47,0.03) -- (0.53,-0.03);",
           r"\node[cWrong,above,font=\scriptsize] at (0.5,0.03) {極點 $z=\tfrac12$};",
           r"\fill (0,0) circle (0.6pt);",
           r"\end{tikzpicture}"]
    return "\n".join(out) + "\n"


def fig_saddle():
    N, a = 10, 2
    zs = W.sd_saddle_exact(N, a)
    xs = [0.62 + 0.34 * k / 80 for k in range(81)]
    re_pts = " ".join(f"({x:.4f},{math.exp(W.sd_p(x, N, a)) / 1000:.4f})" for x in xs)
    ys = [-0.15 + 0.3 * k / 80 for k in range(81)]
    im_pts = " ".join(f"({y:.4f},{abs(math.exp(W.sd_p(complex(zs, y), N, a).real)) / 1000:.4f})" for y in ys)
    v0 = math.exp(W.sd_p(zs, N, a))
    out = [r"\begin{tikzpicture}[>=Stealth]",
           r"\begin{scope}[xshift=-8.4cm,x=14cm,y=2.4cm]",
           r"\draw[->] (0.6,0) -- (0.99,0) node[right]{\small $x$};",
           r"\draw[->] (0.6,0) -- (0.6,1.95);",
           rf"\draw[cBar,thick] plot[smooth] coordinates {{{re_pts}}};",
           rf"\draw[dashed,gray] ({zs:.4f},0) -- ({zs:.4f},{v0 / 1000:.4f});",
           rf"\fill[cWrong] ({zs:.4f},{v0 / 1000:.4f}) circle (1.6pt);",
           rf"\node[below,font=\scriptsize] at ({zs:.4f},0) {{$z^*=\tfrac{{11}}{{13}}$}};",
           r"\node[below,font=\scriptsize] at (0.7,0) {$0.7$};",
           r"\node[below,font=\scriptsize] at (0.95,0) {$0.95$};",
           r"\draw (0.6,1) -- (0.607,1) node[left,font=\scriptsize] {$1000$};",
           r"\node[font=\small] at (0.79,2.1) {沿實軸：谷};",
           r"\end{scope}",
           r"\begin{scope}[xshift=9.4cm,x=16cm,y=14cm]",
           r"\draw[->] (-0.16,0) -- (0.17,0) node[right]{\small $y$};",
           r"\draw[->] (0,0) -- (0,0.33);",
           rf"\draw[cGauss,thick] plot[smooth] coordinates {{{im_pts}}};",
           rf"\fill[cWrong] (0,{v0 / 1000:.4f}) circle (1.6pt);",
           r"\node[below,font=\scriptsize] at (-0.15,0) {$-0.15$};",
           r"\node[below,font=\scriptsize] at (0.15,0) {$0.15$};",
           r"\draw (-0.004,0.1) -- (0.004,0.1) node[left=3pt,font=\scriptsize] {$100$};",
           r"\node[font=\small] at (0,0.36) {沿 $z=z^*+iy$：峰};",
           r"\end{scope}",
           r"\end{tikzpicture}"]
    assert round(v0) == 265 and math.exp(W.sd_p(0.62, N, a)) > 1000 > v0
    return "\n".join(out) + "\n"


ERRATA = [
    ("p. 6，(1.13) 前", r"每步位移的平方叫 $l$", r"$l^2$"),
    ("p. 8，(1.20) 前", r"$x=N-2m$", r"$x=l(2m-N)$（奇數次動差為 0，不影響結果）"),
    ("p. 8--9，(1.22)(1.23)", r"右邊沒有 $e^{Ny/2}$", r"乘上 $e^{Ny/2}$（$y=0$ 時等於 1）"),
    ("p. 11，(1.34)(1.36)", r"$\exp[N(\ln N-N)]$", r"$\exp[N\ln N-N]$"),
    ("p. 12，(1.37)", r"$+\tfrac12\ln(2\pi n)$、$+\tfrac12\ln(2\pi(N-n))$", r"兩項都是 $-$"),
    ("p. 12，(1.38)", r"$+\tfrac12(2m/N)^2$", r"$-\tfrac12(2m/N)^2$"),
    ("p. 12，(1.39)(1.40)", r"$n$；$-\tfrac12\ln 2\pi N$；$2^N/\sqrt{2\pi N}$",
     r"$2m=d/l$；$-\tfrac12\ln\tfrac{\pi N}{2}$；$2\cdot2^N/\sqrt{2\pi N}$"),
    ("p. 14，(1.48)", r"$\exp\big((x-y)^2/2l^2N\big)$", r"$\exp\big(-(x-y)^2/2l^2N\big)$，$\alpha=1/\sqrt{2\pi l^2}$"),
    ("p. 16，(1.54) 後", r"上下限換成 $\pm\infty$ 誤差可忽略", r"還要加 $q=\pm\pi/l$ 的鞍點（因子 2 或 0）"),
    ("p. 17，(1.58) 後", r"$N!\,d^N/dz^N$", r"$\tfrac1{N!}\,d^N/dz^N$"),
    ("p. 20，練習 1.7", r"$d^4=3(d_N^2)^2$", r"$\overline{d^4}=3N^2-2N\approx3(\overline{d^2})^2$"),
    ("p. 21，(1:S-9)", r"$\exp(N|p''|\cdots)$", r"$\exp(-\tfrac N2|p''|\cdots)$"),
    ("p. 23，(1:S-15)", r"前因子 $\tfrac1{2\pi}$", r"$\tfrac1{\sqrt{2\pi}}$"),
]


def tab_errata():
    rows = [r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule"]
    rows += [rf"{w} & {a} & {b}\\" for w, a, b in ERRATA]
    rows += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(rows) + "\n"


FIGURES = [
    ("FigMainWalk", fig_main_walk),
    ("TabThreeCoins", tab_three_coins),
    ("TabTenDist", tab_ten_dist),
    ("FigDeltaGrid", fig_delta_grid),
    ("TabMoments", tab_moments),
    ("FigKurtosis", fig_kurtosis),
    ("FigGamma", fig_gamma),
    ("TabStirling", tab_stirling),
    ("FigBinomGauss", fig_binom_gauss),
    ("FigPascal", fig_pascal),
    ("FigDiffusion", fig_diffusion),
    ("FigChiPower", fig_chi_power),
    ("FigContour", fig_contour),
    ("FigSaddle", fig_saddle),
    ("TabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch1()
    parts = ["% figs.tex —— 由 ch1/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch1/make_figs.py（第 1 章，動畫 W101–W112）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
