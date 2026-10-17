r"""
make_figs.py —— 隨機漫步導讀第 2 章（Rudnick & Gaspari《Elements of the Random Walk》第 2 章，動畫 W201–W212）的圖與表

在 notes/ 或 notes/ch2/ 裡執行：
    python ch2/make_figs.py
會先檢查 ch2/part2.tex 寫到的每一個數字（check_numbers_ch2），再寫出 ch2/figs.tex（part2.tex 第一行 \input{ch2/figs}）。
圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots）。數學一律由 codes/utils/walk_gf.py 算出；例子與動畫 W201–W212 相同。
巨集名稱都以 \Two 開頭（不和第 1 章的 \Fig…、\Tab… 撞名）。顏色沿用導言區：cBar 精確值、cGauss 高斯近似、cExact 更正、cWrong 筆誤、cWalk 足跡。
"""
from __future__ import annotations

import cmath
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent            # notes/ch2
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_gf as GF  # noqa: E402

OUT = HERE / "figs.tex"
F = Fraction


def thin(n: int) -> str:
    """整數加千分位細空白：1\,172\,556。"""
    return f"{n:,}".replace(",", r"\,")


def frac_tex(q: Fraction) -> str:
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    s = "-" if q < 0 else ""
    return rf"{s}\tfrac{{{abs(q.numerator)}}}{{{q.denominator}}}"


def poly_tex(p) -> str:
    """係數 list（x^0 起，Fraction）→ TeX，例如 \tfrac{3}{2}x^{2}-\tfrac{1}{2}。"""
    terms = []
    for k in range(len(p) - 1, -1, -1):
        c = Fraction(p[k])
        if c == 0:
            continue
        mono = "" if k == 0 else ("x" if k == 1 else f"x^{{{k}}}")
        a = abs(c)
        coef = "" if (a == 1 and k > 0) else frac_tex(a)
        sign = "-" if c < 0 else "+"
        terms.append((sign, coef + mono))
    s = "".join(f"{sg}{t}" for sg, t in terms)
    return s[1:] if s.startswith("+") else s


# =============================================================================
# 數字檢查（part2.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch2():
    GF.check()                                                     # 模組本身的全部檢查（書上的數字、筆誤、場景卡）
    # 2.1.1
    assert GF.exp_coeffs(4) == [1, 1, F(1, 2), F(1, 6), F(1, 24)]
    assert GF.half_binom_series(3) == [1, F(1, 2), F(3, 8), F(5, 16)]
    L = GF.legendre_gf(5)
    assert poly_tex(L[2]) == r"\tfrac{3}{2}x^{2}-\tfrac{1}{2}" and poly_tex(L[5]).startswith(r"\tfrac{63}{8}x^{5}")
    assert all(GF.poly_eval(p, -1) == (-1) ** n for n, p in enumerate(L))                 # 練習：P_n(−1) = (−1)^n
    assert round(GF.bessel_series(0, 2.5), 4) == -0.0484 and abs(GF.bessel_integral(0, 2.5) - GF.bessel_series(0, 2.5)) < 1e-9
    assert round(math.exp(2), 3) == 7.389
    assert [round(p, 4) for p in GF.grand_probs(2, 1.0, 5)] == [0.1353, 0.2707, 0.2707, 0.1804, 0.0902, 0.0361]
    mw = GF.mortal_walkers(2000, 0.9, seed=2)
    assert round(sum(n for n, _ in mw) / 2000, 2) == 9.34 and sum(1 for _, p in mw if p == (0, 0)) == 296
    assert round(GF.mortal_origin_exact(0.9), 4) == 0.1452 and 296 / 2000 == 0.148
    assert max(n for n, _ in mw) == 78
    assert round(GF.ellipk(0.9), 4) == 2.2805
    # 2.1.2
    g2, g3, g4 = (GF.count_grid(N) for N in (2, 3, 4))
    assert g2[(0, 0)] == 4 and g2[(1, 1)] == 2 and g2[(2, 0)] == 1 and g3[(1, 0)] == 9 == GF.count_sq(3, 1, 0)
    assert g4[(0, 0)] == 36 == 6 * 6 and g4[(2, 0)] == 16 == 4 * 4 and g4[(1, 1)] == 24 == 4 * 6
    assert GF.count_sq(6, 1, 1) == 300 == 15 * 20
    assert abs(GF.chi_1d(math.pi / 3) - 1) < 1e-12
    # 2.1.3
    assert GF.returns_sq(10) == 63504 == 252 ** 2 and GF.count_sq(10, 1, 1) == 52920 == 210 * 252
    assert (round(GF.gauss_sc(10)), round(2 * GF.gauss_sc(10)), GF.returns_sc(10)) == (630848, 1261696, 1172556)
    assert round(GF.returns_sc(10) / GF.gauss_sc(10), 2) == 1.86
    assert [round(GF.returns_sc(N) / (2 * GF.gauss_sc(N)), 3) for N in (10, 20)] == [0.929, 0.964]
    assert (round(GF.gauss_sq(10)), round(GF.asym_returns_sq(10))) == (33377, 66754)
    # 2.1.4
    assert round((3 / (2 * math.pi)) ** 1.5, 4) == 0.3299
    for kind in ("sc", "bcc", "fcc"):
        assert round(scaled_return(kind, 40), 4) == 0.3238
    assert round(math.exp(-12 / (2 * math.sqrt(3) * 10)), 4) == 0.7072 and round(math.exp(-0.2), 4) == 0.8187
    assert round(4 * math.pi, 3) == 12.566
    assert round(GF.chi_cauchy(1.0), 4) == 0.3679
    assert round(1 - GF.cauchy_radial_cdf(10.0), 4) == 0.1265                  # 練習：一步超過 10 r0
    assert round(GF.CHI3_MEDIAN, 3) == 1.538 and round(GF.cauchy_radial_median(), 3) == 2.264
    # 2.1.5
    P20 = GF.prob_grid_sq(20)
    assert round(float(P20[(0, 0)]), 6) == 0.031045 and round(2 / (20 * math.pi), 6) == 0.031831
    assert GF.variance(P20) == 20
    assert (round(GF.g_fourier(0.4, GF.chi_1d(0.3)), 5), round(GF.g_cont(0.4, 0.3), 5), round(GF.g_cont_book(0.4, 0.3), 5)) \
        == (4.24213, 4.23729, 6.09756)
    assert (round(GF.G1_exact(0.45, 0), 4), round(GF.G1_continuum(0.45, 0), 4)) == (2.2942, 2.357)
    # 2.2
    assert [round(GF.gauss_moment(n), 4) for n in (0, 2, 4, 6)] == [1.7725, 0.8862, 1.3293, 3.3234]
    assert round(GF.gauss_gen(1, 2), 3) == 4.818
    assert GF.gauss_moment_frac(4) / GF.gauss_moment_frac(2) ** 2 == 3
    assert GF.gauss_moment_frac(6) / GF.gauss_moment_frac(2) ** 3 == 15
    assert abs(GF.gauss_moment(2, 3.0) - math.sqrt(math.pi / 3) / 6) < 1e-15 and round(GF.gauss_moment(2, 3.0), 4) == 0.1706
    # 2.3
    assert (round(GF.CELL_VOLUME["bcc"], 4), round(GF.CELL_VOLUME["fcc"], 4)) == (0.7698, 0.7071)
    bs = GF.reciprocal(*GF.BCC_BASIS)                                             # bcc 的倒晶格是 fcc
    for b in bs:
        nz = sorted(abs(t) for t in b)
        assert abs(nz[0]) < 1e-12 and abs(nz[1] - nz[2]) < 1e-12 and abs(nz[1] - math.pi * math.sqrt(3)) < 1e-12
    assert all(abs(GF.roots_sum(m, 6)) < 1e-12 for m in range(1, 6)) and abs(GF.roots_sum(0, 6) - 6) < 1e-12
    # 2.4
    assert [round(GF.coeff_two_poles(N, 1, 1, 0.5, 1.0)) for N in range(6)] == [3, 5, 9, 17, 33, 65]
    assert round(GF.coeff_two_poles(10, 1, 0, 0.5, 1.0)) == 2048
    assert math.comb(20, 10) == 184756 and round(4 ** 10 / math.sqrt(10 * math.pi), 2) == 187078.97
    assert round(2 ** 10 * math.sqrt(2 / (10 * math.pi)), 1) == 258.4
    assert [round(GF.returns_sq(N) / GF.asym_returns_sq(N), 4) for N in (10, 100)] == [0.9513, 0.995]
    assert (round(GF.coeff_log_exact(10, 1.0, 0.5), 1), round(GF.coeff_log_asym(10, 1.0, 0.5), 1),
            round(GF.coeff_log_asym(10, 1.0, 0.5, book=True), 1), round(GF.coeff_log_asym(10, 1.0, 0.5, terms=1), 1)) \
        == (7418.1, 7317.4, 4478.3, 4715.7)
    jac = [(2 ** (N + 1) + (-1) ** N) // 3 for N in range(8)]                      # 練習：1/((1−2z)(1+z))
    ser = [1, 1]
    for N in range(2, 8):
        ser.append(ser[-1] + 2 * ser[-2])
    assert jac == ser == [1, 1, 3, 5, 11, 21, 43, 85]
    c = GF.series_coeffs_power(1.5, 1.0, 100)[100]                                # 練習：(1−z)^{−3/2} ~ 2√(N/π)
    assert round(c, 2) == 11.33 and round(2 * math.sqrt(100 / math.pi), 2) == 11.28


NMAX_UNIV = 40
_FCC_CACHE: list[float] = []


def fcc_return_prob(N: int) -> float:
    """fcc 的 P(N;0) = C(N;0)/12^N，N ≤ 40：均分求積（每軸 M = 41 點，對次數 ≤ 40 的三角多項式精確），
    一次算完所有 N（逐步 count_grid 到 N = 40 太慢）。N ≤ 10 與 GF.returns_fcc 的精確整數核對。"""
    if not _FCC_CACHE:
        M = NMAX_UNIV + 1
        cs = [math.cos(-math.pi + 2 * math.pi * j / M) for j in range(M)]
        acc = [0.0] * (NMAX_UNIV + 1)
        for a in cs:
            for b in cs:
                for c in cs:
                    x = (a * b + b * c + c * a) / 3                 # χ/12
                    v = 1.0
                    for N in range(NMAX_UNIV + 1):
                        acc[N] += v
                        v *= x
        _FCC_CACHE.extend(t / M ** 3 for t in acc)
        for n in range(11):
            assert abs(_FCC_CACHE[n] - GF.returns_fcc(n) / 12 ** n) < 1e-12
    return _FCC_CACHE[N]


def scaled_return(kind: str, N: int) -> float:
    if kind == "fcc":
        p = fcc_return_prob(N)
    else:
        f = {"sc": GF.returns_sc, "bcc": GF.returns_bcc}[kind]
        p = f(N) / GF.CHI_SMALL_Q[kind][0] ** N
    g = 2 if GF.BIPARTITE[kind] else 1
    return N ** 1.5 * p / (g * GF.CELL_VOLUME[kind])


# =============================================================================
# 圖與表
# =============================================================================
def tab_legendre():
    L = GF.legendre_gf(5)
    rows = []
    for n, p in enumerate(L):
        i1 = GF.poly_int11(p)
        i2 = GF.poly_int11(GF.poly_mul(p, p))
        rows.append(rf"{n} & ${poly_tex(p)}$ & ${frac_tex(i1)}$ & ${frac_tex(i2)}$ & ${frac_tex(F(2, 2 * n + 1))}$\\")
    return "\n".join([r"\begin{tabular}{@{}clccc@{}}", r"\toprule",
                      r"$n$ & $P_n(x)$ & $\int_{-1}^1P_n\,dx$ & $\int_{-1}^1P_n^2\,dx$ & $\frac{2}{2n+1}$\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def tab_grand():
    Z = [GF.ideal_gas_Z(N, 2) for N in range(6)]
    P = GF.grand_probs(2, 1.0, 5)
    return "\n".join([r"\begin{tabular}{@{}c" + "c" * 6 + "@{}}", r"\toprule",
                      r"$N$ & " + " & ".join(str(N) for N in range(6)) + r"\\", r"\midrule",
                      r"$Z(N)=2^N/N!$ & " + " & ".join(f"${frac_tex(z)}$" for z in Z) + r"\\",
                      r"$Z(N)z^N/Q$（$z=1$）& " + " & ".join(f"{p:.4f}" for p in P) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def fig_mortal():
    mw = GF.mortal_walkers(2000, 0.9, seed=2)
    hist = [sum(1 for n, _ in mw if n == k) for k in range(21)]
    sx, sy = 0.36, 0.012
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\draw[->] (-0.8,0) -- (21.2,0) node[right,font=\scriptsize] {$n$};",
           r"\draw[->] (-0.8,0) -- (-0.8,240) node[above,font=\scriptsize] {人數};"]
    for k, h in enumerate(hist):
        out.append(rf"\fill[cBar!55] ({k - 0.38},0) rectangle ({k + 0.38},{h});")
    pts = " ".join(f"({k},{2000 * 0.1 * 0.9 ** k:.2f})" for k in range(21))
    out.append(rf"\draw[cGauss,thick] plot coordinates {{{pts}}};")
    out.append(rf"\foreach \k in {{0,5,10,15,20}} \node[below,font=\scriptsize] at (\k,0) {{\k}};")
    out.append(r"\foreach \v in {100,200} {\draw (-0.95,\v) -- (-0.65,\v); \node[left,font=\scriptsize] at (-0.9,\v) {\v};}")
    out.append(r"\node[right,font=\scriptsize,cGauss] at (6,200) {理論 $2000\times0.1\times0.9^n$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_count_grid():
    out = [r"\begin{tikzpicture}[>=Stealth]"]
    # (a) Fig. 2.3：鄰居
    s = 0.9
    out.append(r"\begin{scope}[shift={(0,0)}]")
    out.append(rf"\draw[gray!40] (-1.5*{s},-1.5*{s}) grid[step={s}] (1.5*{s},1.5*{s});")
    for k, (x, y, lab) in enumerate([(-1, 0, "w_1"), (0, -1, "w_2"), (1, 0, "w_3"), (0, 1, "w_4")]):
        out.append(rf"\fill[gray] ({x * s},{y * s}) circle (2.6pt);")
        out.append(rf"\draw[->,gray!70!black,thick] ({x * s * 0.85},{y * s * 0.85}) -- ({x * s * 0.2},{y * s * 0.2});")
        anchor = "below left" if x + y < 0 else "above right"
        out.append(rf"\node[{anchor},font=\scriptsize] at ({x * s},{y * s}) {{$\vec {lab}$}};")
    out.append(r"\fill (0,0) circle (3pt); \node[below right,font=\scriptsize] at (0,0) {$\vec y$};")
    out.append(r"\node[font=\small] at (0,-2.0) {(a) 最後一步來自四個鄰居};")
    out.append(r"\end{scope}")
    # (b)(c) 數字格
    for shift, N, title in ((4.2, 2, "(b) $N=2$"), (9.3, 4, "(c) $N=4$")):
        g = GF.count_grid(N)
        c = 0.5
        out.append(rf"\begin{{scope}}[shift={{({shift},0)}}]")
        out.append(rf"\draw[gray!30] ({-N - 0.5}*{c},{-N - 0.5}*{c}) grid[step={c}] ({N + 0.5}*{c},{N + 0.5}*{c});")
        for (x, y), v in sorted(g.items()):
            col = "cWalk" if (x, y) == (0, 0) else "black"
            out.append(rf"\node[font=\scriptsize,text={col}] at ({x * c},{y * c}) {{{v}}};")
        out.append(rf"\node[font=\small] at (0,{-N * c - 0.75}) {{{title}，總和 ${4 ** N}$}};")
        out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_chi():
    rows = [("一維", r"$2\cos q$", 2, r"$2-q^2$", "1", "是"),
            ("正方", r"$2(\cos q_x+\cos q_y)$", 4, r"$4-q^2$", "1", "是"),
            ("簡單立方 sc", r"$2(\cos q_x+\cos q_y+\cos q_z)$", 6, r"$6-q^2$", "1", "是"),
            ("體心立方 bcc", r"$8\cos\frac{q_x}{\sqrt3}\cos\frac{q_y}{\sqrt3}\cos\frac{q_z}{\sqrt3}$", 8, r"$8-\frac43q^2$",
             f"{GF.CELL_VOLUME['bcc']:.4f}", "是"),
            ("面心立方 fcc", r"$4\big[\cos\frac{q_x}{\sqrt2}\cos\frac{q_y}{\sqrt2}+\cdots\big]$", 12, r"$12-2q^2$",
             f"{GF.CELL_VOLUME['fcc']:.4f}", "否")]
    for name, _, Z, _, _, _ in rows:
        key = {"一維": "1d", "正方": "sq", "簡單立方 sc": "sc", "體心立方 bcc": "bcc", "面心立方 fcc": "fcc"}[name]
        assert GF.CHI_SMALL_Q[key][0] == Z
    body = [rf"{a} & {b} & {c} & {d} & {e} & {f}\\" for a, b, c, d, e, f in rows]
    return "\n".join([r"\begin{tabular}{@{}llcccc@{}}", r"\toprule",
                      r"格子 & $\chi(\vec q)$（$a=1$） & $\chi(\vec0)$ & 小 $q$ & 單胞體積 $\Omega$ & 二部圖\\",
                      r"\midrule", *body, r"\bottomrule", r"\end{tabular}"])


def tab_returns():
    Ns = range(0, 11)
    fs = [("一維 $\\binom{N}{N/2}$", lambda N: GF.count1(N, 0)), ("正方", GF.returns_sq), ("sc", GF.returns_sc),
          ("bcc", GF.returns_bcc), ("fcc", GF.returns_fcc)]
    rows = [name + " & " + " & ".join(thin(f(N)) for N in Ns) + r"\\" for name, f in fs]
    return "\n".join([r"{\small\setlength{\tabcolsep}{3.2pt}\begin{tabular}{@{}l" + "r" * 11 + "@{}}", r"\toprule",
                      r"$N$ & " + " & ".join(str(N) for N in Ns) + r"\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}}"])


def _heat_tikz(fn, n: int, size: float) -> list[str]:
    out = []
    h = size / n
    for i in range(n):
        for j in range(n):
            x = -math.pi + (j + 0.5) * 2 * math.pi / n
            y = -math.pi + (i + 0.5) * 2 * math.pi / n
            v = max(-1.0, min(1.0, fn(x, y)))
            col = f"orange!{round(100 * v)}!white" if v >= 0 else f"blue!{round(-70 * v)}!white"
            out.append(rf"\fill[{col}] ({j * h:.3f},{i * h:.3f}) rectangle ({(j + 1) * h:.3f},{(i + 1) * h:.3f});")
    return out


def fig_heat():
    out = [r"\begin{tikzpicture}"]
    size, n = 3.2, 20
    for k, N in enumerate((1, 2, 10)):
        out.append(rf"\begin{{scope}}[shift={{({k * 4.3},0)}}]")
        out += _heat_tikz(lambda x, y, N=N: (GF.chi_sq(x, y) / 4) ** N, n, size)
        out.append(rf"\draw[gray] (0,0) rectangle ({size},{size});")
        out.append(rf"\node[below,font=\scriptsize] at ({size / 2},0) {{$N={N}$}};")
        out.append(rf"\node[left,font=\scriptsize] at (0,0) {{$-\pi$}}; \node[right,font=\scriptsize] at ({size},{size}) {{$\pi$}};")
        out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def sci(x: float, digits: int = 3) -> str:
    e = int(math.floor(math.log10(abs(x))))
    return rf"${x / 10 ** e:.{digits}f}\times10^{{{e}}}$"


def tab_gauss():
    def fmt(v, N):
        return thin(round(v)) if N <= 10 else sci(v)
    sc_rows, sq_rows = [], []
    for N in (2, 4, 6, 8, 10, 20):
        ex, g = GF.returns_sc(N), GF.gauss_sc(N)
        sc_rows.append(rf"{N} & {fmt(ex, N)} & {fmt(g, N)} & {fmt(2 * g, N)} & {ex / (2 * g):.3f}\\")
    for N in (2, 10, 20, 50, 100):
        ex, g, a = GF.returns_sq(N), GF.gauss_sq(N), GF.asym_returns_sq(N)
        sq_rows.append(rf"{N} & {fmt(ex, N)} & {fmt(g, N)} & {fmt(a, N)} & {ex / a:.4f}\\")
    return "\n".join([r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule",
                      r"$N$ & 精確 $C(N;\vec0)$ & 書上的高斯 & 乘 2（奇偶性） & 精確／乘 2\\", r"\midrule",
                      r"\multicolumn{5}{@{}l}{\small 簡單立方：高斯 $6^N(3/2\pi N)^{3/2}$，即 (2.20) 取 $\vec r=0$}\\",
                      *sc_rows, r"\midrule",
                      r"\multicolumn{5}{@{}l}{\small 正方（主例）：高斯 $4^N/\pi N$}\\", *sq_rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_univ():
    rows = []
    for kind, name in (("sc", "sc"), ("bcc", "bcc"), ("fcc", "fcc")):
        f = {"sc": GF.returns_sc, "bcc": GF.returns_bcc, "fcc": GF.returns_fcc}[kind]
        Z = GF.CHI_SMALL_Q[kind][0]
        ps = " & ".join(f"{f(N) / Z ** N:.5f}" for N in (2, 4, 6, 8, 10))
        ss = " & ".join(f"{scaled_return(kind, N):.4f}" for N in (10, 20, 40))
        rows.append(rf"{name} & {ps} & {ss}\\")
    return "\n".join([r"{\small\setlength{\tabcolsep}{4pt}\begin{tabular}{@{}lccccc|ccc@{}}", r"\toprule",
                      r" & \multicolumn{5}{c|}{$P(N;\vec0)=C(N;\vec0)/Z^N$} & \multicolumn{3}{c}{$N^{3/2}P/(g\Omega)$}\\",
                      r"$N$ & 2 & 4 & 6 & 8 & 10 & 10 & 20 & 40\\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}}"])


def fig_univ():
    sx, sy = 0.16, 22.0
    y0 = 0.22
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           rf"\draw[->] (0,{y0}) -- (42,{y0}) node[right,font=\scriptsize] {{$N$}};",
           rf"\draw[->] (0,{y0}) -- (0,0.345) node[above,font=\scriptsize] {{$N^{{3/2}}P(N;\vec0)/(g\Omega)$}};"]
    lim = (3 / (2 * math.pi)) ** 1.5
    out.append(rf"\draw[dashed,gray] (0,{lim:.4f}) -- (42,{lim:.4f}) node[right,font=\scriptsize] {{$(3/2\pi)^{{3/2}}$}};")
    styles = {"sc": ("cBar", "circle"), "bcc": ("cGauss", "square*"), "fcc": ("cExact", "triangle*")}
    for kind, (col, mark) in styles.items():
        pts = " ".join(f"({N},{scaled_return(kind, N):.4f})" for N in range(2, 41, 2))
        out.append(rf"\draw[{col},thick] plot[mark={mark},mark size=1.1pt] coordinates {{{pts}}};")
    out.append(r"\foreach \n in {10,20,30,40} \node[below,font=\scriptsize] at (\n,0.22) {\n};")
    out.append(r"\foreach \v in {0.25,0.30} {\draw (-0.4,\v) -- (0.4,\v); \node[left,font=\scriptsize] at (-0.4,\v) {\v};}")
    out.append(r"\node[cBar,font=\scriptsize,right] at (24,0.255) {sc};")
    out.append(r"\node[cGauss,font=\scriptsize,right] at (24,0.243) {bcc};")
    out.append(r"\node[cExact,font=\scriptsize,right] at (24,0.231) {fcc};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_coin():
    rows = []
    for N, d in ((10, 0), (10, 2), (20, 0), (20, 4), (40, 0)):
        ex, ga, co = GF.coin_quartic(N, d)
        rows.append(rf"{N} & {d} & {ex:.6f} & {ga:.6f} & {co:.6f} & {100 * (ga - ex) / ex:+.2f}\% & {100 * (co - ex) / ex:+.3f}\%\\")
    return "\n".join([r"\begin{tabular}{@{}rrccccc@{}}", r"\toprule",
                      r"$N$ & $d$ & 精確 $P(N,d)$ & 高斯 $\times2$ & 加 $q^4$ 修正 & 高斯的誤差 & 修正後的誤差\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def tab_median():
    rows = []
    for N in (1, 10, 100):
        mg = GF.median_distance("gauss", N, 400, seed=5)
        mc = GF.median_distance("cauchy", N, 400, seed=5)
        rows.append(rf"{N} & {mg:.2f} & {GF.CHI3_MEDIAN * math.sqrt(N):.2f} & {mc:.2f} & {GF.cauchy_radial_median() * N:.2f}\\")
    return "\n".join([r"\begin{tabular}{@{}rcccc@{}}", r"\toprule",
                      r" & \multicolumn{2}{c}{高斯步長（Ex.~2.7）} & \multicolumn{2}{c}{Cauchy 步長 (2.26)}\\",
                      r"$N$ & 模擬 & $1.538\sqrt N$ & 模擬 & $2.264N$\\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def fig_diffusion():
    out = [r"\begin{tikzpicture}"]
    c = 0.11
    for k, N in enumerate((4, 10, 20)):
        P = GF.prob_grid_sq(N)
        vmax = float(max(P.values()))
        out.append(rf"\begin{{scope}}[shift={{({k * 5.0},0)}}]")
        out.append(rf"\fill[gray!8] ({-20.5 * c},{-20.5 * c}) rectangle ({20.5 * c},{20.5 * c});")
        for (x, y), v in sorted(P.items()):
            t = round(100 * float(v) / vmax)
            if t >= 2:
                out.append(rf"\fill[orange!{t}!white] ({(x - 0.5) * c:.3f},{(y - 0.5) * c:.3f}) rectangle ({(x + 0.5) * c:.3f},{(y + 0.5) * c:.3f});")
        out.append(rf"\draw[gray] ({-20.5 * c},{-20.5 * c}) rectangle ({20.5 * c},{20.5 * c});")
        out.append(rf"\node[below,font=\scriptsize] at (0,{-20.5 * c}) {{$N={N}$，$P(N;\vec0)={float(P[(0, 0)]):.4f}$}};")
        out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_moments():
    rows = []
    for n in (0, 2, 4, 6, 8):
        q = GF.gauss_moment_frac(n)
        rows.append(rf"{n} & ${frac_tex(q)}\sqrt\pi$ & {GF.gauss_moment(n):.4f} & {GF.gauss_moment_numeric(n):.4f}\\")
    return "\n".join([r"\begin{tabular}{@{}cccc@{}}", r"\toprule",
                      r"$n$ & $\frac{n!}{(n/2)!\,4^{n/2}}\sqrt\pi$ & 數值 & 數值積分\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_roots():
    out = [r"\begin{tikzpicture}[>=Stealth,scale=0.9]"]
    for k, (m, L, title, shift) in enumerate(((0, 0.55, "$m=0$：和 $=8$", 0.0), (1, 0.55, "$m=1$：和 $=0$", 6.0),
                                               (3, 1.2, "$m=3$：和 $=0$", 9.6))):
        p = complex(shift, 0)
        pts = [p]
        for j in range(8):
            p += L * cmath.exp(2j * math.pi * m * j / 8)
            pts.append(p)
        for a, b in zip(pts, pts[1:]):
            out.append(rf"\draw[->,teal!70!black,thick] ({a.real:.3f},{a.imag:.3f}) -- ({b.real:.3f},{b.imag:.3f});")
        cx = shift + (2.2 if m == 0 else 0.6)
        out.append(rf"\node[font=\scriptsize] at ({cx:.2f},-0.9) {{{title}}};")
        assert abs(pts[-1] - pts[0] - (8 * L if m == 0 else 0)) < 1e-9
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_poles():
    out = [r"\begin{tikzpicture}[>=Stealth,scale=1.15]"]
    # (a) 兩個極點
    out.append(r"\begin{scope}")
    out.append(r"\draw[->,gray] (-1.3,0) -- (1.5,0) node[right,font=\scriptsize] {$\mathrm{Re}\,z$};")
    out.append(r"\draw[->,gray] (0,-1.2) -- (0,1.2) node[above,font=\scriptsize] {$\mathrm{Im}\,z$};")
    out.append(r"\draw[cExact,dashed] (0,0) circle (0.5);")
    out.append(r"\node[cWrong,font=\small] at (0.5,0) {$\times$}; \node[cWrong,above,font=\scriptsize] at (0.5,0.05) {$\frac12$};")
    out.append(r"\node[gray,font=\small] at (1,0) {$\times$}; \node[gray,above,font=\scriptsize] at (1,0.05) {$1$};")
    out.append(r"\node[font=\scriptsize,align=center] at (0.1,-1.55) {(a) $\frac1{1/2-z}+\frac1{1-z}$：近的極點說了算};")
    out.append(r"\end{scope}")
    # (b) 一維回到原點：±1/2 的分支點
    out.append(r"\begin{scope}[shift={(4.3,0)}]")
    out.append(r"\draw[->,gray] (-1.3,0) -- (1.5,0); \draw[->,gray] (0,-1.2) -- (0,1.2);")
    out.append(r"\draw[cExact,dashed] (0,0) circle (0.5);")
    out.append(r"\draw[cWrong,thick,decorate,decoration={zigzag,amplitude=1pt,segment length=3pt}] (0.5,0) -- (1.4,0);")
    out.append(r"\draw[cWrong,thick,decorate,decoration={zigzag,amplitude=1pt,segment length=3pt}] (-0.5,0) -- (-1.3,0);")
    out.append(r"\fill[cWrong] (0.5,0) circle (1.6pt) (-0.5,0) circle (1.6pt);")
    out.append(r"\node[font=\scriptsize,align=center] at (0.1,-1.55) {(b) $(1-4z^2)^{-1/2}$：$z=\pm\frac12$ 一樣近};")
    out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_asym():
    rows = []
    for N in (5, 10, 20, 50, 100):
        ex = math.comb(2 * N, N)
        asy = 4.0 ** N / math.sqrt(math.pi * N)
        ex2 = GF.returns_sq(2 * N)
        asy2 = GF.asym_returns_sq(2 * N)
        rows.append(rf"{N} & {ex / asy:.4f} & {2 * N} & {ex2 / asy2:.4f}\\")
    return "\n".join([r"\begin{tabular}{@{}rc|rc@{}}", r"\toprule",
                      r"$N$ & $\binom{2N}{N}\Big/\frac{4^N}{\sqrt{\pi N}}$ & $N$ & $C_{\text{正方}}(N;\vec0)\Big/\frac{2\cdot4^N}{\pi N}$\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    ("p.~26，(2.1) 後", r"把右邊對 $y$ 展開", r"把左邊對 $z$ 展開"),
    ("p.~29，(2.6)", r"$\chi(k)=2\cos q$", r"$\chi(q)=2\cos q$"),
    ("p.~34，(2.20) 前", "中心以外沒有可觀的貢獻", r"角落也一樣大：乘 2 或 0"),
    ("p.~35，(2.22)", r"$e^{-3r^2/2aN}$", r"$e^{-3r^2/2a^2N}$"),
    ("p.~36，(2.24)", r"左邊 $d\vec q$", r"$d^3r$"),
    ("p.~37，Ex.~2.7", r"$\sqrt{2/\pi}\,r_0^{-3}$（積分 $4\pi$）", r"$(2\pi)^{-3/2}r_0^{-3}$"),
    ("p.~39，(2.37) 後", "(2.20) 的左邊", "(2.20) 的右邊"),
    ("p.~40，(2.38)", r"$za^2\nabla^2G-(2dz-1)G=\delta$", r"$za^2\nabla^2G+(2dz-1)G=-\delta$"),
    ("p.~41，(2:S1-4)", r"$-\frac{b^2}{2a}$", r"$-\frac{b^2}{4a}$"),
    ("p.~42，(2:S1-5)", r"$-a(x-\frac b{2a})+\frac{b^2}{4a}$", r"$-a(x-\frac b{2a})^2+\frac{b^2}{4a}$"),
    ("p.~43", r"$l_2,l_2,l_3$", r"$l_1,l_2,l_3$"),
    ("p.~44", r"布里淵區體積寫成 $\Omega$", r"$\Omega_{BZ}=(2\pi)^3/a^3$"),
    ("p.~45，(2:S2-23)", r"右邊 $f(\vec R)$", r"$F(\vec R)$"),
    ("p.~47，(2:S2-34) 後", r"$m_i$ 不等於 $N_i$", r"$m_i$ 不是 $N_i$ 的倍數"),
    ("p.~50，(2:S3-11)", r"$+\ln z_c$", r"$-\ln z_c$"),
    ("p.~50，(2:S3-13)", r"$y_c$；$N^{y-1}$", r"$z_c$；$N^{-y-1}$"),
]


def tab_errata():
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("TwoTabLegendre", tab_legendre),
    ("TwoTabGrand", tab_grand),
    ("TwoFigMortal", fig_mortal),
    ("TwoFigCountGrid", fig_count_grid),
    ("TwoTabChi", tab_chi),
    ("TwoTabReturns", tab_returns),
    ("TwoFigHeat", fig_heat),
    ("TwoTabGauss", tab_gauss),
    ("TwoTabUniv", tab_univ),
    ("TwoFigUniv", fig_univ),
    ("TwoTabCoin", tab_coin),
    ("TwoTabMedian", tab_median),
    ("TwoFigDiffusion", fig_diffusion),
    ("TwoTabMoments", tab_moments),
    ("TwoFigRoots", fig_roots),
    ("TwoFigPoles", fig_poles),
    ("TwoTabAsym", tab_asym),
    ("TwoTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch2()
    parts = ["% figs.tex —— 由 ch2/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch2/make_figs.py（第 2 章，動畫 W201–W212）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
