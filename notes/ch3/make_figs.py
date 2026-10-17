r"""
make_figs.py —— 隨機漫步導讀第 3 章（Rudnick & Gaspari《Elements of the Random Walk》第 3 章，動畫 W301–W308）的圖與表

在 notes/ 或 notes/ch3/ 裡執行：
    python ch3/make_figs.py
會先檢查 ch3/part3.tex 寫到的每一個數字（check_numbers_ch3），再寫出 ch3/figs.tex（part3.tex 第一行 \input{ch3/figs}）。
圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots）。數學一律由 codes/utils/walk_rec.py（必要時 walk_gf.py）算出；例子與動畫 W301–W308 相同。
巨集名稱都以 \Three 開頭（不和第 1 章的 \Fig…、第 2 章的 \Two… 撞名）；part3.tex 的 \label 用 three-…。
顏色沿用導言區：cBar 精確值、cGauss 近似、cExact 更正、cWrong 筆誤、cWalk 足跡。
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent            # notes/ch3
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_gf as GF  # noqa: E402
from utils import walk_rec as R  # noqa: E402

OUT = HERE / "figs.tex"
F = Fraction
NS = (10, 100, 1000, 10000)
NMAX = 10000


def thin(n: int) -> str:
    """整數加千分位細空白：63\,504。"""
    return f"{n:,}".replace(",", r"\,")


def _r(xs, k):
    return [round(x, k) for x in xs]


def _simpson(f, a, b, n=4000):
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) + 4 * sum(f(a + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(a + 2 * i * h) for i in range(1, n // 2)))


# =============================================================================
# 數字檢查（part3.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch3():
    # ---- 3.1 Pólya 的問題（W301）
    assert 2 ** 4 == 16 and R.never_return_1d(4) == 6 and 16 - 6 == 10
    for k in ("1d", "sq", "sc"):
        fr = R.simulate_first_returns(k, R.SIM_M, R.SIM_N, R.SIM_SEED[k])
        assert sum(t is not None for t in fr) == R.SIM_RETURNED[k]
        assert sum(t is not None and t <= 100 for t in fr) == R.SIM_RETURNED_100[k]
    assert [R.SIM_RETURNED[k] for k in ("1d", "sq", "sc")] == [977, 676, 325]
    assert round(1 - R.escape("sc"), 4) == 0.3405 and round(R.G_C["sc"], 6) == 1.516386
    # ---- 3.2 數到訪次數（W302）
    assert R.visit_counts_1d(4, 1) == {0: 6, 1: 6, 2: 4} == R.visit_counts_gf_1d(4, 1)
    assert R.visit_poly({0: 6, 1: 6, 2: 4}) == [16, -14, 4]
    assert R.counts_from_poly([16, -14, 4]) == {0: 6, 1: 6, 2: 4} and -(-14 + 2 * 4) == 6 and 2 * 4 // 2 == 4
    assert R.virial_terms_1d(4, 1) == [16, 14, 4] and 14 == 1 * 2 ** 3 + 3 * 2 and 4 == 1 * 2 * 2
    assert GF.count1(1, 1) == 1 and GF.count1(3, 1) == 3 and GF.count1(2, 0) == 2
    groups = {0: [], 1: [], 2: []}
    for pos in R.walks_1d(3):
        groups[pos.count(1)].append("".join("+" if b > a else "-" for a, b in zip(pos, pos[1:])))
    assert groups == {0: ["-+-", "--+", "---"], 1: ["+++", "+--", "-++"], 2: ["++-", "+-+"]}
    assert R.visit_counts_1d(4, 0) == {1: 6, 2: 6, 3: 4}
    # ---- 3.3 一維（W303）
    assert [R.never_return_1d(N) for N in range(11)] == [1, 2, 2, 4, 6, 12, 20, 40, 70, 140, 252] == R.never_return_1d_gf(10)
    assert R.never_return_1d(10) == math.comb(10, 5) == 252 and R.escape_1d(10) == F(252, 1024)
    assert round(252 / 1024, 4) == 0.2461 and R.never_return_1d(9) == 2 * math.comb(8, 4) == 140
    G1d = [GF.count1(k, 0) for k in range(11)]
    assert G1d == [1, 0, 2, 0, 6, 0, 20, 0, 70, 0, 252]
    z = 0.3                                                         # (3.18)(3.19)(3.20)：第二行要有 1/2π；分母 α² − α/z + 1
    I = _simpson(lambda q: 1 / (1 - 2 * z * math.cos(q)), 0, 2 * math.pi)
    assert abs(I / (2 * math.pi) - 1 / math.sqrt(1 - 4 * z * z)) < 1e-12 and I > 6
    ap, am = (1 / z + math.sqrt(1 / z ** 2 - 4)) / 2, (1 / z - math.sqrt(1 / z ** 2 - 4)) / 2
    assert abs(ap * am - 1) < 1e-12 and abs(ap + am - 1 / z) < 1e-12 and am < 1 < ap
    assert abs(am ** 2 - am / 2 + 1 / z) > 1                       # 書上的分母在 α_− 不是 0
    assert 1 / 0.7 ** 2 - 4 < 0                                     # 1/2 < z < 1：兩個極點都在單位圓上
    assert abs(1 / (z * (ap - am)) - 1 / math.sqrt(1 - 4 * z * z)) < 1e-12        # (3.21)
    assert abs(GF.coeff_power_asym(10, 0.5, 0.5) / 2 ** 10 - math.sqrt(2 / (10 * math.pi))) < 1e-14
    assert round(math.sqrt(2), 4) == 1.4142
    assert [round(float(R.escape_1d(N)), 4) for N in (10, 100, 1000)] == [0.2461, 0.0796, 0.0252]
    assert [round(R.escape_1d_asym(N), 4) for N in (10, 100, 1000)] == [0.2523, 0.0798, 0.0252]
    assert [round(R.escape_1d_book(N), 4) for N in (10, 100, 1000)] == [0.1784, 0.0564, 0.0178]
    PB = [R.escape_biased_series(0.6, N) for N in (2, 4, 10, 20, 50, 100)]
    assert _r(PB, 4) == [0.52, 0.4048, 0.294, 0.2427, 0.2094, 0.2016]
    assert round(1 / R.escape_biased(0.6), 9) == 5 and round(R.escape_biased(0.6), 12) == 0.2
    assert round(1 - R.escape_biased(0.6), 12) == 0.8 and R.escape_biased(1.0) == 1
    # ---- 3.4 二維（W304）
    assert [R.returns_count("sq", N) for N in range(0, 11, 2)] == [1, 4, 36, 400, 4900, 63504]
    fs = R.first_return_probs("sq", 10)
    assert _r([fs[n] for n in (2, 4, 6, 8, 10)], 4) == [0.25, 0.0781, 0.043, 0.0286, 0.021]
    assert [fs[n] * 4 ** n for n in (2, 4)] == [4, 20] and round(sum(fs), 4) == 0.4208
    assert round(R.no_return_probs("sq", 10)[10], 4) == 0.5792
    assert round(R.no_return_sq_book(10), 3) == 1.364 and round(R.log10_steps_for(0.01)) == 136
    assert abs(2 / math.pi * GF.ellipk(0.999999) - (math.log(8) - math.log(1e-6)) / math.pi) < 1e-5   # G ≈ (ln 8 − ln(1−x))/π
    assert 1000 - R.SIM_RETURNED["sq"] == 324 and round(R.no_return_probs("sq", 1000)[1000], 4) == 0.3218
    zc, N = 0.25, 1000                                              # 鞍點 z* = z_c(N+1)/(N+2)；f''(z*) ≈ 1/δ²
    zs = zc * (N + 1) / (N + 2)
    fpp = 1 / (zc - zs) ** 2 + (N + 1) / zs ** 2
    assert abs(fpp * (zc - zs) ** 2 - 1) < 2e-3 and abs((zc - zs) / (zc / N) - 1) < 2e-3
    # ---- 3.5 維度（W305）
    assert _r([R.K_d(d) for d in (1, 2, 3, 4)], 3) == [2.0, 6.283, 12.566, 19.739]
    assert abs(R.K_d(2) - 2 * math.pi) < 1e-12 and abs(R.K_d(3) - 4 * math.pi) < 1e-12
    T = [[R.tail_integral(d, e) for e in (0.1, 0.01, 0.001)] for d in (1, 2, 3)]
    assert _r(T[0], 6) == [9, 99, 999] and _r(T[1], 3) == [2.303, 4.605, 6.908] and _r(T[2], 6) == [0.9, 0.99, 0.999]
    c1, c2 = 6, 1                                                   # c2 = 1：χ_sc(q,0,0) = 4 + 2cos q ≈ 6 − q²
    assert abs((6 - GF.chi_sc(1e-3, 0, 0)) / 1e-6 - c2) < 1e-6
    B1 = R.K_d(3) / (2 * math.pi) ** 3 / c1 * (c1 ** 2 / c2) ** 1.5 * (math.pi / 2)           # (3.28)(3.32)，a → c1
    assert abs(B1 - 9 / math.pi) < 1e-12 and round(B1, 3) == 2.865
    A = R.sc_correction_coeff()
    assert round(A, 5) == 0.28696 and round(1 / R.G_C["sc"], 6) == 0.659463
    P3 = R.no_return_probs("sc", 1000)
    assert _r([P3[100], P3[1000]], 6) == [0.688062, 0.668534]
    assert _r([1 / R.G_C["sc"] + A / math.sqrt(n) for n in (100, 1000)], 6) == [0.688159, 0.668537]
    # ---- 3.6 Watson（W306）
    assert [round(R.G_C[k], 6) for k in ("sc", "bcc", "fcc")] == [1.516386, 1.393204, 1.344661]
    assert [round(R.escape(k), 6) for k in ("sc", "bcc", "fcc")] == [0.659463, 0.71777, 0.743682]
    assert abs(R.watson_sc_bessel() - R.G_C["sc"]) < 1e-10 and abs(R.watson_bcc_series() - R.G_C["bcc"]) < 1e-10
    assert abs(R.watson_fcc_series() - R.G_C["fcc"]) < 1e-3
    assert math.floor(R.G_C["bcc"] * 1e6) / 1e6 == 1.393203                 # 書上印的是截斷值
    assert [R.Z[k] for k in ("sc", "bcc", "fcc")] == [6, 8, 12]
    assert round(1 - R.escape("sc"), 6) == 0.340537 and round(R.G_C["sc"] - 1, 6) == 0.516386
    assert _r([R.f_m("sc", m) for m in range(5)], 4) == [0.6595, 0.2246, 0.0765, 0.026, 0.0089]
    # ---- 3.7 走過幾個點（W307–W308）
    assert R.new_site_probs_brute("1d", 4) == [1, 1, F(1, 2), F(1, 2), F(3, 8)]
    assert [R.distinct_sites_brute("1d", n) for n in range(5)] == [1, 2, F(5, 2), 3, F(27, 8)]
    assert sum(len(set(p)) for p in R.walks_1d(4)) == 54 == 16 * F(27, 8)
    assert [R.never_return_1d(n) / 2 ** n for n in range(5)] == [1, 1, 0.5, 0.5, 0.375]
    assert R.distinct_sites_1d(4) == 9 * F(6, 16) and R.distinct_sites_1d(3) == 6 * F(2, 4)
    assert round(float(R.distinct_sites_1d(100)), 3) == 15.997 and round(R.footprint_1d_asym(100), 3) == 15.958
    assert [R.distinct_sites_brute("sq", n) for n in range(6)] == [1, 2, F(11, 4), F(7, 2), F(267, 64), F(155, 32)]
    assert R.TRAIL_SITES == 159 and round(R.distinct_sites("sq", 400)[400], 2) == 162.57
    D3 = R.distinct_sites("sc", 1000)
    assert [round(D3[n] / n, 4) for n in (10, 100, 1000)] == [0.9051, 0.7229, 0.6782]
    assert round(2 * A, 4) == 0.5739 and round(1000 / R.G_C["sc"] + 2 * A * math.sqrt(1000), 1) == 677.6
    assert round(D3[1000], 1) == 678.2
    # ---- 練習
    assert round(R.escape("bcc"), 6) == 0.71777 and round(sum(R.escape_biased_series(0.6, N) for N in range(101)), 2) == 24.95
    assert sum(R.f_m("sc", m) for m in range(400)) > 1 - 1e-12
    Gx = [GF.count1(k, 1) for k in range(6)]                       # 一維從 0 到 +1：z + 3z³ + 10z⁵
    assert Gx == [0, 1, 0, 3, 0, 10]


# =============================================================================
# 表與圖
# =============================================================================
def tab_sim():
    rows = []
    for name, k, inf in (("一維（直線）", "1d", "1"), ("二維（正方格子）", "sq", "1"),
                         ("三維（簡單立方）", "sc", f"{1 - R.escape('sc'):.4f}")):
        P = R.no_return_probs(k, 1000)
        rows.append(rf"{name} & {R.SIM_RETURNED_100[k]} & {1000 * (1 - P[100]):.1f} & {R.SIM_RETURNED[k]} & "
                    rf"{1000 * (1 - P[1000]):.1f} & {inf}\\")
    return "\n".join([r"\begin{tabular}{@{}lrrrrc@{}}", r"\toprule",
                      r" & \multicolumn{2}{c}{100 步內回過家} & \multicolumn{2}{c}{1000 步內回過家} & 走無窮多步\\",
                      r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
                      r"格子 & 模擬 & 理論 & 模擬 & 理論 & 回家的機率\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_virial():
    """Fig. 3.1 的樣子：從 x 到 y 的路線經過大圓點 w 一次、兩次、三次。"""
    out = [r"\begin{tikzpicture}[>=Stealth,walk/.style={decorate,decoration={snake,amplitude=0.7pt,segment length=5pt},"
           r"cWalk,thick},dot/.style={circle,fill=black,inner sep=1.1pt}]"]
    for i, (sign, loops) in enumerate((("-t", 0), ("+t^2", 1), ("-t^3", 2))):
        x0 = 4.2 * i
        out.append(rf"\begin{{scope}}[xshift={x0}cm]")
        out.append(r"\node[dot,label=below:{\scriptsize$\vec x$}] (x) at (0,0) {};")
        out.append(r"\node[dot,label=below:{\scriptsize$\vec y$}] (y) at (3,0) {};")
        out.append(r"\coordinate (w) at (1.5,0.9);")
        out.append(r"\draw[walk] (x) .. controls (0.4,0.9) .. (w);")
        out.append(r"\draw[walk] (w) .. controls (2.6,0.9) .. (y);")
        for j in range(loops):
            a = 60 + 60 * j if loops == 2 else 90
            out.append(rf"\draw[walk] (w) .. controls +({a - 30}:1.1) and +({a + 30}:1.1) .. (w);")
        out.append(r"\node[circle,fill=cWrong,inner sep=2.6pt] at (w) {};")
        out.append(rf"\node[font=\small] at (1.5,-0.75) {{${sign}$}};")
        out.append(r"\end{scope}")
    out.append(r"\node[font=\small] at (12.6,0.3) {$+\cdots$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_never():
    Ns = range(11)
    nr = [R.never_return_1d(N) for N in Ns]
    rows = [r"從不回家 & " + " & ".join(str(v) for v in nr) + r"\\",
            r"全部 $2^N$ & " + " & ".join(thin(2 ** N) for N in Ns) + r"\\",
            r"$P_N$ & " + " & ".join(f"{v / 2 ** N:.3f}" for N, v in zip(Ns, nr)) + r"\\",
            r"第 $N$ 步在家 & " + " & ".join(str(R.CT(N)) if N % 2 == 0 else "--" for N in Ns) + r"\\"]
    return "\n".join([r"{\small\setlength{\tabcolsep}{3.6pt}\begin{tabular}{@{}l" + "r" * 11 + "@{}}", r"\toprule",
                      r"$N$ & " + " & ".join(str(N) for N in Ns) + r"\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}}"])


def tab_no_return():
    P1, P2, P3 = (R.no_return_probs(k, NMAX) for k in ("1d", "sq", "sc"))
    A, g = R.sc_correction_coeff(), 1 / R.G_C["sc"]
    rows = []
    for N in NS:
        rows.append(rf"$10^{{{round(math.log10(N))}}}$ & {P1[N]:.4f} & {R.escape_1d_asym(N):.4f} & "
                    rf"{P2[N]:.4f} & {R.no_return_sq_book(N):.4f} & {R.no_return_sq_refined(N):.4f} & "
                    rf"{P3[N]:.5f} & {g + A / math.sqrt(N):.5f}\\")
    return "\n".join([r"{\small\setlength{\tabcolsep}{4pt}\begin{tabular}{@{}rcccccccc@{}}", r"\toprule",
                      r" & \multicolumn{2}{c}{一維} & \multicolumn{3}{c}{二維（正方）} & \multicolumn{2}{c}{三維（sc）}\\",
                      r"\cmidrule(lr){2-3}\cmidrule(lr){4-6}\cmidrule(lr){7-8}",
                      r"$N$ & 精確 & $\sqrt{2/\pi N}$ & 精確 & $\pi/\ln N$ & $\frac{\pi}{\ln N+\ln8+\gamma}$ & 精確 & "
                      r"$0.659463+\frac{0.28696}{\sqrt N}$\\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}}"])


def fig_tail():
    sx, sy = 6.0, 0.54
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\fill[cExact!15] (0,0) rectangle (1,1);",
           r"\draw[->] (0,0) -- (1.1,0) node[right,font=\scriptsize] {$q$};",
           r"\draw[->] (0,0) -- (0,5.6) node[above,font=\scriptsize] {$q^{d-1}/q^{2}$};"]

    def pts(f, q0):
        return " ".join(f"({q0 + (1 - q0) * i / 60:.4f},{f(q0 + (1 - q0) * i / 60):.4f})" for i in range(61))

    out.append(rf"\draw[cBar,thick] plot coordinates {{{pts(lambda q: q ** -2, 5 ** -0.5)}}};")
    out.append(rf"\draw[cGauss,thick] plot coordinates {{{pts(lambda q: 1 / q, 0.2)}}};")
    out.append(r"\draw[cExact,thick] (0,1) -- (1,1);")
    out.append(r"\draw[cBar,->] (0.4472,5) -- (0.4472,5.5); \draw[cGauss,->] (0.2,5) -- (0.2,5.5);")
    out.append(r"\node[cBar,font=\scriptsize,right] at (0.56,5.0) {$d=1$：$1/q^2$，面積 $\to\infty$};")
    out.append(r"\node[cGauss,font=\scriptsize,right] at (0.56,4.3) {$d=2$：$1/q$，面積 $\to\infty$（對數）};")
    out.append(r"\node[cExact,font=\scriptsize,right] at (0.56,3.6) {$d=3$：$1$，面積有限（綠色）};")
    out.append(r"\foreach \x in {0.5,1} \node[below,font=\scriptsize] at (\x,0) {\x};")
    out.append(r"\node[below,font=\scriptsize] at (0,0) {0};")
    out.append(r"\foreach \v in {1,2,3,4,5} {\draw (-0.012,\v) -- (0.012,\v); \node[left,font=\scriptsize] at (-0.012,\v) {\v};}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_gzero():
    xs = (0.9, 0.99, 0.999)
    rows = [r"一維 $1/\sqrt{1-x^2}$ & " + " & ".join(f"{R.G0('1d', x):.3f}" for x in xs) + r" & $\infty$\\",
            r"正方 $\frac2\pi K(x)$ & " + " & ".join(f"{R.G0('sq', x):.3f}" for x in xs) + r" & $\infty$\\",
            r"sc $\int_0^\infty e^{-t}I_0(xt/3)^3dt$ & " + " & ".join(f"{R.G0('sc', x):.3f}" for x in xs)
            + rf" & {R.G_C['sc']:.3f}\\"]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$G(z;0)$，$z=x\,z_c$ & $x=0.9$ & $0.99$ & $0.999$ & $1$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_lattices():
    forms = {"sc": r"$\frac{\sqrt6}{32\pi^3}\Gamma(\frac1{24})\Gamma(\frac5{24})\Gamma(\frac7{24})\Gamma(\frac{11}{24})$",
             "bcc": r"$\frac{\Gamma(1/4)^4}{4\pi^3}$", "fcc": r"$\frac{9\,\Gamma(1/3)^6}{2^{14/3}\pi^4}$"}
    rows = []
    for k in ("sc", "bcc", "fcc"):
        G = R.G_C[k]
        rows.append(rf"{k} & {R.Z[k]} & {forms[k]} & {G:.6f} & {1 / G:.6f} & {1 - 1 / G:.6f} & {G - 1:.6f}\\")
    return "\n".join([r"{\small\begin{tabular}{@{}lclcccc@{}}", r"\toprule",
                      r"格子 & 鄰居 & $G(z_c;0)$ 的閉式 & 數值 & 逃脫 $1/G$ & 回家 $1-1/G$ & 平均回家次數 $G-1$\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}}"])


def tab_fm():
    R0 = 1 - R.escape("sc")
    th = [R.f_m("sc", m) for m in range(4)] + [R0 ** 4]
    fin = R.return_count_dist("sc", 1000, 3)
    fin = fin + [1 - sum(fin)]
    sim = [c / R.SIM_M for c in R.SIM_COUNTS_SC]
    rows = [r"$f_m=R^m(1-R)$（走無窮多步） & " + " & ".join(f"{v:.4f}" for v in th) + r"\\",
            r"1000 步內的精確值 & " + " & ".join(f"{v:.4f}" for v in fin) + r"\\",
            r"1000 人的模擬 & " + " & ".join(f"{v:.3f}" for v in sim) + r"\\"]
    return "\n".join([r"\begin{tabular}{@{}lccccc@{}}", r"\toprule",
                      r"回家 $m$ 次 & $0$ & $1$ & $2$ & $3$ & $\ge4$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_foot():
    sx, sy = 2.6, 3.7
    D = {k: R.distinct_sites(k, NMAX) for k in ("1d", "sq", "sc")}
    grid = sorted({round(10 ** (1 + 3 * i / 90)) for i in range(91)})
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           r"\draw[->] (1,0) -- (4.15,0) node[right,font=\scriptsize] {$n$};",
           r"\draw[->] (1,0) -- (1,1.05) node[above,font=\scriptsize] {$\delta_n/n$};"]
    for k, col in (("1d", "cBar"), ("sq", "cGauss"), ("sc", "cExact")):
        pts = " ".join(f"({math.log10(n):.4f},{D[k][n] / n:.4f})" for n in grid)
        out.append(rf"\draw[{col},very thick] plot coordinates {{{pts}}};")
    a1 = " ".join(f"({math.log10(n):.4f},{R.footprint_1d_asym(n) / n:.4f})" for n in grid)
    out.append(rf"\draw[cBar,dashed] plot coordinates {{{a1}}};")
    a2 = " ".join(f"({math.log10(n):.4f},{R.footprint_2d_book(n) / n:.4f})" for n in grid if n >= 25)
    out.append(rf"\draw[cGauss,dashed] plot coordinates {{{a2}}};")
    a3 = " ".join(f"({math.log10(n):.4f},{R.footprint_2d_mw(n) / n:.4f})" for n in grid)
    out.append(rf"\draw[cGauss,dotted,thick] plot coordinates {{{a3}}};")
    g = 1 / R.G_C["sc"]
    out.append(rf"\draw[cExact,dashed] (1,{g:.4f}) -- (4,{g:.4f});")
    out.append(r"\foreach \e in {1,2,3,4} {\draw (\e,-0.012) -- (\e,0.012); \node[below,font=\scriptsize] at (\e,0) {$10^{\e}$};}")
    out.append(r"\foreach \v in {0.2,0.4,0.6,0.8,1} {\draw (0.97,\v) -- (1.03,\v); \node[left,font=\scriptsize] at (0.97,\v) {\v};}")
    out.append(rf"\node[cExact,font=\scriptsize,below] at (3.3,{g - 0.01:.3f}) {{三維 sc：$\to 1/G=0.6595$（虛線）}};")
    out.append(r"\node[cGauss,font=\scriptsize,right] at (3.0,0.2) {二維：$\pi/\ln n$（虛線）、$\pi/\ln 8n$（點線）};")
    out.append(r"\node[cBar,font=\scriptsize,right] at (2.0,0.27) {一維：$\sqrt{8/\pi n}$（虛線）};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_foot():
    D = {k: R.distinct_sites(k, NMAX) for k in ("1d", "sq", "sc")}
    A, g = R.sc_correction_coeff(), 1 / R.G_C["sc"]
    rows = []
    for n in NS:
        rows.append(rf"$10^{{{round(math.log10(n))}}}$ & {D['1d'][n]:.2f} & {R.footprint_1d_asym(n):.2f} & "
                    rf"{D['sq'][n]:.1f} & {R.footprint_2d_book(n):.1f} & {R.footprint_2d_mw(n):.1f} & "
                    rf"{D['sc'][n]:.1f} & {g * n + 2 * A * math.sqrt(n):.1f}\\")
    return "\n".join([r"{\small\setlength{\tabcolsep}{4pt}\begin{tabular}{@{}rcccccccc@{}}", r"\toprule",
                      r" & \multicolumn{2}{c}{一維} & \multicolumn{3}{c}{二維（正方）} & \multicolumn{2}{c}{三維（sc）}\\",
                      r"\cmidrule(lr){2-3}\cmidrule(lr){4-6}\cmidrule(lr){7-8}",
                      r"$n$ & 精確 $\delta_n$ & $\sqrt{8n/\pi}$ & 精確 $\delta_n$ & $\pi n/\ln n$ & $\pi n/\ln 8n$ & 精確 $\delta_n$ & "
                      r"$0.6595n+0.574\sqrt n$\\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}}"])


ERRATA = [
    (r"p.~52--53，(3.3)(3.5)", r"$\prod_{i=1}^{N}$", r"$\prod_{i=0}^{N}$（起點也算腳印）"),
    (r"p.~56，(3.18) 第二行", r"沒有 $\frac{1}{2\pi}$", r"前面乘 $\frac{1}{2\pi}$"),
    (r"p.~56，(3.19) 第二行", r"$\alpha^2-\alpha/2+1/z$", r"$\alpha^2-\alpha/z+1$"),
    (r"p.~57，(3.20) 後", r"$z<1$", r"$0<z<\frac12$"),
    (r"p.~57，(3.23)", r"$z_c^{-N}N^{-1/2}/\Gamma(\frac12)$", r"$z_c^{-N-1/2}N^{-1/2}/\Gamma(\frac12)$"),
    (r"p.~57，(3.24)", r"$N^{-1/2}/\Gamma(\frac12)=1/\sqrt{\pi N}$", r"$\sqrt{2}N^{-1/2}/\Gamma(\frac12)=\sqrt{2/\pi N}$"),
    (r"p.~59，(3.27) 後", r"sc 的 $c_2=1/2$", r"$c_2=1$"),
    (r"p.~60", r"$f''(z_c)=1/\delta^2$", r"$f''(z_c-\delta)\approx1/\delta^2$"),
    (r"p.~61，(3.32) 與前一句", r"$a$（$q=Qa(z_c-z)^{1/2}/c_2^{1/2}$）", r"$c_1$"),
    (r"p.~66，(3.53) 後", r"$(z_c-z)^{1/2}$ 的係數", r"$(z_c-z)^{-1/2}$ 的係數"),
    (r"p.~66，(3.56)", r"$\sum_{n'=2}^{\infty}$", r"$\sum_{n'=2}^{n}$"),
    (r"p.~67，(3.57)", r"積分下限 1", r"積分下限 2"),
]


def tab_errata():
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("ThreeTabSim", tab_sim),
    ("ThreeFigVirial", fig_virial),
    ("ThreeTabNever", tab_never),
    ("ThreeTabNoReturn", tab_no_return),
    ("ThreeFigTail", fig_tail),
    ("ThreeTabGzero", tab_gzero),
    ("ThreeTabLattices", tab_lattices),
    ("ThreeTabFm", tab_fm),
    ("ThreeFigFoot", fig_foot),
    ("ThreeTabFoot", tab_foot),
    ("ThreeTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch3()
    parts = ["% figs.tex —— 由 ch3/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch3/make_figs.py（第 3 章，動畫 W301–W308）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
