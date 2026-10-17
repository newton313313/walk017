r"""
make_figs.py —— 隨機漫步導讀第 5 章（Rudnick & Gaspari《Elements of the Random Walk》第 5 章，動畫 W501–W507）的圖與表

在 notes/ 或 notes/ch5/ 裡執行：
    python ch5/make_figs.py
會先檢查 ch5/part5.tex 寫到的每一個數字（check_numbers_ch5），再寫出 ch5/figs.tex（part5.tex 第一行 \input{ch5/figs}）。
圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算 Bessel 函數）。數學一律由 codes/utils/walk_var.py 算出；
例子與動畫 W501–W507 相同。巨集名稱都以 \Five 開頭（不和第 1 章的 \Fig…、第 2–4 章的 \Two…、\Three…、\Four… 撞名）；
part5.tex 的 \label 用 five-…。顏色沿用導言區：cBar 精確值、cGauss 近似、cExact 更正、cWrong 筆誤、cWalk 路線。
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import integrate as _int
from scipy import special as _sp

HERE = Path(__file__).resolve().parent            # notes/ch5
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_var as V  # noqa: E402

OUT = HERE / "figs.tex"
F = Fraction
PI = math.pi
TS = (1, 3, 5, 8)                                  # Fig. 5.3、Fig. 5.6 的四個時間
TCOL = {1: "cBar", 3: "cWalk", 5: "cExact", 8: "violet!75!black"}


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣）。"""
    return f"{float(x):.{d}f}"


def _g(x, k=4):
    """k 位有效數字（0.0008387 這類小數也照寫）。"""
    s = f"{float(x):.{k}g}"
    return s if "e" not in s else f"{float(x):.{k + 3}f}".rstrip("0")


def _mc(d):
    vals = V.rt_mc(d, 3.0, 100000, 5)
    return V.mc_fraction(vals, 0, 1.5), float(np.mean(np.abs(vals) > 3.0))


# =============================================================================
# 數字檢查（part5.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch5():
    # ---- 5.1.1 不公正的硬幣（W501）
    p = F(3, 5)
    P10 = V.biased_pmf(10, p)
    assert P10 == V.biased_dp(10, p) and V.moments(P10) == (2, F(48, 5)) and float(F(48, 5)) == 9.6
    assert V.biased_mean(10, p) == 2 and V.biased_var(10, p) == F(48, 5) and V.moments(V.biased_pmf(10, F(1, 2)))[1] == 10
    assert 2 * p - 1 == F(1, 5) and 1 - (2 * p - 1) ** 2 == 4 * p * (1 - p) == F(24, 25)
    assert [_s(P10[x], 4) for x in (0, 2, 4)] == ["0.2007", "0.2508", "0.2150"] and P10[0] != P10[4]
    assert [_s(2 * V.biased_gauss(float(x), 10, 0.6), 4) for x in (0, 2, 4)] == ["0.2091", "0.2575", "0.2091"]
    assert P10[2] == math.comb(10, 6) * p ** 6 * (1 - p) ** 4 and math.comb(10, 6) == 210
    assert _s(0.6 ** 6, 6) == "0.046656" and _s(0.4 ** 4, 4) == "0.0256" and P10[0] == 252 * F(6, 25) ** 5
    for N, mu, var, pk, gk in ((10, 2, "9.6", "0.2508", "0.2575"), (20, 4, "19.2", "0.1797", "0.1821"),
                               (40, 8, "38.4", "0.1279", "0.1288")):
        P = V.biased_pmf(N, p)
        m, v = V.moments(P)
        assert m == mu and str(float(v)) == var and _s(max(P.values()), 4) == pk
        assert max(P, key=P.get) == mu and _s(2 * V.biased_gauss(float(mu), N, 0.6), 4) == gk
    assert V.taylor_residual(0.6) < 1e-5 and V.taylor_residual(0.6, book=True) > 1e3 * V.taylor_residual(0.6)
    # Ex. 5.1
    assert V.diffusion_D(0.6, kind="naive") == 0.5 and abs(V.diffusion_D(0.6, kind="careful") - 0.48) < 1e-15
    assert _s(math.sqrt(10), 3) == "3.162" and _s(math.sqrt(9.6), 3) == "3.098" and abs(2 * 0.48 * 10 - 9.6) < 1e-12
    assert _s(V.width_ratio(0.6), 4) == "0.9798" and abs(V.width_ratio(0.9) - 0.6) < 1e-15
    assert round(100 * (1 - V.width_ratio(0.6))) == 2 and round(100 * (1 - V.width_ratio(0.9))) == 40
    # ---- 5.1.2 先走一步再亂走
    m, v = V.moments(dict(V.shift_pmf(10, 0.2)))
    assert abs(m - 2) < 1e-12 and abs(v - 10) < 1e-12 and V.shift_pmf(10, 0.2) == V.shift_dp(10, 0.2)
    # ---- 5.2.1 兩個狀態（W502）
    Q = F(3, 4)
    h = F(1, 2)
    LR = [V.persistent_dp(N, Q) for N in range(5)]
    assert LR[0] == {0: (h, h)} and LR[1] == {-1: (h, 0), 1: (0, h)}
    assert LR[2] == {-2: (F(3, 8), 0), 0: (F(1, 8), F(1, 8)), 2: (0, F(3, 8))}
    assert LR[3] == {-3: (F(9, 32), 0), -1: (F(1, 8), F(3, 32)), 1: (F(3, 32), F(1, 8)), 3: (0, F(9, 32))}
    assert LR[4] == {-4: (F(27, 128), 0), -2: (F(15, 128), F(9, 128)), 0: (F(13, 128), F(13, 128)),
                     2: (F(9, 128), F(15, 128)), 4: (0, F(27, 128))}
    assert all(V.persistent_pmf(N, Q) == V._persistent_brute(N, Q) for N in range(1, 7))
    assert Q * h == F(3, 8) and (1 - Q) * h == F(1, 8) and 2 * F(3, 8) * 4 == 3
    assert [V.persistent_var(N, Q) for N in range(1, 5)] == [1, 3, F(11, 2), F(33, 4)]
    assert [V.persistent_var(N, h) for N in range(1, 5)] == [1, 2, 3, 4]
    assert [V.persistent_var(N, F(1, 4)) for N in range(1, 5)] == [1, 1, F(3, 2), F(7, 4)]
    assert all(V.persistent_var(N, q) == V.persistent_var_formula(N, q) for N in range(1, 12) for q in (Q, F(1, 4)))
    l1, l2 = V.eigs(0.1, 0.75)
    assert _s(l1.real, 4) == "0.9848" and _s(l2.real, 4) == "0.5077" and l1.imag == l2.imag == 0
    c1, c2 = V.eigs(0.5, 0.75)
    assert _s(c1.real, 4) == "0.6582" and _s(abs(c1.imag), 4) == "0.2584" and abs(c1 - c2.conjugate()) < 1e-15
    assert 0.75 * math.sin(0.5) > 0.25 > 0.75 * math.sin(0.1)
    for N, s in ((5, "0.5699"), (20, "0.0388")):
        a, b, c = V.total_fourier(N, 0.3, 0.75), V.T_total(N, 0.3, 0.75), V.total_eig(N, 0.3, 0.75)
        assert _s(a.real, 4) == s and abs(a - b) < 1e-12 and abs(a - c) < 1e-12
    T = V.T_matrix(0.3, 0.75)
    assert abs(np.trace(T) - 2 * 0.75 * math.cos(0.3)) < 1e-15 and abs(np.linalg.det(T) - (2 * 0.75 - 1)) < 1e-15
    # ---- 5.2.2 持續長度、二階遞迴（W503）
    assert _s(V.lambda1_small(0.1, 0.75), 3) == "0.985"
    assert [_s(V.persistent_var(N, 0.75), 2) for N in (10, 100, 1000)] == ["26.00", "296.00", "2996.00"]
    assert [V.persistent_var_asym(N, 0.75) for N in (10, 100, 1000)] == [30, 300, 3000]
    assert [_s(V.persistent_var(N, 0.25), 3) for N in (10, 100, 1000)] == ["3.777", "33.778", "333.778"]
    assert [_s(V.persistent_var_asym(N, 0.25), 3) for N in (10, 100, 1000)] == ["3.333", "33.333", "333.333"]
    assert abs(V.persistent_var(1000, 0.75) - 3000 + 4) < 1e-9 and abs(V.persistent_var(1000, 0.25) - 1000 / 3 - 4 / 9) < 1e-9
    assert abs(V.persistent_var(1000, 0.75) - float(V.persistent_var_formula(1000, F(3, 4)))) < 1e-8
    assert _s(V.persistent_pmf(100, 0.75)[0], 5) == "0.04587" and _s(2 * V.gauss_persistent(0.0, 100, 0.75), 5) == "0.04607"
    assert [V.persistence_length(q) for q in (Q, h, F(1, 4))] == [3, 1, F(1, 3)] and 1000 // 3 == 333
    assert V.persistent_var(10, F(1)) == 100 and V.persistent_var(10, F(0)) == 0 and V.persistent_var(11, F(0)) == 1
    assert all(abs(V.recursion2_residual(N, q, 0.75)) < 1e-12 for N in (3, 6, 10) for q in (0.1, 0.3))
    assert _s(abs(V.eig_residual(l1, 0.1, 0.75, book=True)), 3) == "1.000" and abs(V.eig_residual(l1, 0.1, 0.75)) < 1e-12
    assert abs(abs(V.eig_residual(l1, 0.1, 0.75, book=True)) - 2 * abs(1 - 2 * 0.75)) < 1e-12
    # ---- 5.2.3 電報方程（W504）
    lat = []
    for dt in (0.01, 0.005):
        xs, P = V.tele_lattice(1.0, dt)
        i = int(np.argmin(np.abs(xs)))
        assert abs(xs[i]) < 1e-12
        lat.append(P[i] / (2 * dt))
    assert [_s(x, 4) for x in lat] == ["0.3392", "0.3380"] and _s(2 * lat[1] - lat[0], 4) == "0.3368"
    S0 = [float(V.tele_cont(0.0, t)) for t in TS]
    assert [_s(x, 4) for x in S0] == ["0.3368", "0.2199", "0.1738", "0.1388"]
    assert _s(math.exp(-1) / 2, 4) == "0.1839" and _s(_sp.i0(1.0), 4) == "1.2661" and _s(_sp.i1(1.0), 4) == "0.5652"
    assert abs(math.exp(-1) / 2 * (_sp.i0(1.0) + _sp.i1(1.0)) - S0[0]) < 1e-14
    fr = [float(V.tele_cont(t * (1 - 1e-9), t)) for t in TS]
    assert [_g(x) for x in fr] == ["0.2759", "0.06223", "0.01179", "0.0008387"]
    assert all(abs(f - math.exp(-t) / 2 * (1 + t / 2)) < 1e-6 for f, t in zip(fr, TS))
    assert [_g(V.tele_spike(t)) for t in TS] == ["0.1839", "0.02489", "0.003369", "0.0001677"]
    assert [_s(1 - math.exp(-t), 4) for t in TS] == ["0.6321", "0.9502", "0.9933", "0.9997"]
    assert all(abs(V.tele_mass(t) - 1) < 1e-9 for t in TS) and _s(V.tele_mass(3, v=2.0, book=True), 4) == "1.9502"
    assert abs(V.tele_fourier(0.7, 2.0) - V.tele_fourier_numeric(0.7, 2.0)) < 1e-9
    assert _s(math.exp(-3), 4) == "0.0498"
    mc1 = _mc(1)
    assert _s(mc1[0], 4) == "0.6060" and _s(mc1[1], 4) == "0.0488" and _s(V.radial_fraction(1, 0, 1.5, 3), 4) == "0.6051"
    assert V.tele_D() == 0.5 and V.tele_lp() == 1.0 and V.tele_D() == 0.5 * 1.0 * V.tele_lp()
    assert [_s(V.tele_gauss(0.0, t), 4) for t in TS] == ["0.3989", "0.2303", "0.1784", "0.1410"]
    assert [_s(V.tele_var(t), 4) for t in TS] == ["0.5677", "2.5012", "4.5000", "7.5000"]
    assert all(abs(V.tele_var(t) - V.tele_var_numeric(t)) < 1e-6 for t in (1, 3))
    # ---- 5.2.4 二維、三維（W505）
    assert _s(V.angle_avg_2d_closed(1, 1).real, 4) == "0.4472" and abs(V.angle_avg_2d(1, 1) - V.angle_avg_2d_closed(1, 1)) < 1e-12
    assert _s(V.angle_avg_3d_closed(1, 1).real, 4) == "0.4636" and abs(V.angle_avg_3d(1, 1) - V.angle_avg_3d_closed(1, 1)) < 1e-9
    assert V.sm_residual() < 1e-12
    g = V.angle_avg_2d_closed(0.6, 0.9)
    assert abs(g / (1 - g) - 1 / (math.sqrt(1.6 ** 2 + 0.81) - 1)) < 1e-12
    assert [_s(V.rt2_radial(r, 3), 4) for r in (1, 2, 2.5, 2.9)] == ["0.2978", "0.4167", "0.3941", "0.4052"]
    rr = np.array([0.5, 1.0, 2.5, 2.9])
    assert np.allclose(V.rt2_pole(rr, 3) + V.rt2_cut(rr, 3), V.rt2_cont(rr, 3), rtol=1e-13)
    assert abs(V.rt2_mass(3) - 1) < 1e-5 and _s(1 - math.exp(-3), 4) == "0.9502"
    assert _s(V.rt2_mass(3, book=True), 4) == "1.2658" and _s(V.rt2_mass(1, book=True), 4) == "1.9340"
    assert all(abs(V.rt3_mass(t) - 1) < 1e-5 for t in TS)
    assert _s(3 * math.exp(-3), 4) == "0.1494" and _s(1 - 4 * math.exp(-3), 4) == "0.8009"
    assert abs(_int.quad(lambda u: u * math.log((1 + u) / (1 - u)), 0, 1)[0] - 1) < 1e-9
    assert _s(V.rt3_pole_book(3.5, 3), 4) == "0.0113" and _s(V.rt3_cut_book(3.5, 3), 4) == "-0.0113"
    assert abs(V.rt3_5106(2.0, 3, book=False) + V.rt3_pole_book(2.0, 3)) < 1e-12
    assert abs(V.rt3_5106(2.0, 3, book=True) + V.rt3_pole_book(2.0, 3)) > 1e-3
    mc2, mc3 = _mc(2), _mc(3)
    assert _s(mc2[0], 4) == "0.3335" and _s(V.radial_fraction(2, 0, 1.5, 3), 4) == "0.3310"
    assert _s(mc3[0], 4) == "0.2880" and _s(V.radial_fraction(3, 0, 1.5, 3), 4) == "0.2881"
    assert _s(mc2[1], 4) == "0.0493" and _s(mc3[1], 4) == "0.0497"
    P3 = {t: [V.rt3_radial(f * t, t) for f in (0.25, 0.5, 0.75, 0.9)] for t in TS}
    assert [_s(x, 4) for x in P3[1]] == ["0.1472", "0.5334", "1.0275", "1.3699"]
    assert [_s(x, 4) for x in P3[3]] == ["0.1729", "0.4558", "0.4828", "0.3594"]
    assert [_s(x, 4) for x in P3[5]] == ["0.1869", "0.3630", "0.2169", "0.0926"]
    assert [_s(x, 4) for x in P3[8]] == ["0.1948", "0.2411", "0.0622", "0.0118"]
    assert [_g(math.exp(-t)) for t in TS] == ["0.3679", "0.04979", "0.006738", "0.0003355"]
    # ---- 5.3.1 等待時間（W506）
    assert _s(V.survival(3), 4) == "0.0498" and _s(V.survival(3, "gamma2"), 4) == "0.1991"
    assert abs(V.survival(3, "gamma2") - 4 * math.exp(-3)) < 1e-15
    assert abs(V.survival_numeric(3, "gamma2") - V.survival(3, "gamma2")) < 1e-8 and V.mean_wait("gamma2") == 2
    A, B = V.AB(0.7)
    assert _s(A, 4) == _s(B, 4) == "0.5882" and abs(B - 1 / 1.7) < 1e-15
    A2, B2 = V.AB(0.7, "gamma2")
    assert _s(A2, 4) == "0.3460" and _s(B2, 4) == "0.9343" and abs(A2 - (1 - 0.7 * B2)) < 1e-15
    assert all(abs(x - y) < 1e-8 for x, y in zip(V.AB_numeric(0.7, "gamma2"), (A2, B2)))
    assert abs(B2 / (1 - A2) - 1 / 0.7) < 1e-12
    pm1 = [V.steps_prob(m, 3) for m in range(7)]
    pm2 = [V.steps_prob(m, 3, "gamma2") for m in range(7)]
    assert [_s(x, 4) for x in pm1] == ["0.0498", "0.1494", "0.2240", "0.2240", "0.1680", "0.1008", "0.0504"]
    assert [_s(x, 4) for x in pm2[:5]] == ["0.1991", "0.4481", "0.2689", "0.0720", "0.0108"]
    assert abs(pm2[1] - 9 * math.exp(-3)) < 1e-15 and 9 / 2 + 27 / 6 == 9
    assert all(abs(V.steps_prob_laplace(m, 3, k) - V.steps_prob(m, 3, k)) < 1e-7 for m in (0, 2) for k in ("poisson", "gamma2"))
    assert V.mean_steps(3) == 3 and _s(V.mean_steps(3, "gamma2"), 4) == "1.2506"
    assert abs(V.mean_steps(3, "gamma2") - (1.5 - 0.25 + math.exp(-6) / 4)) < 1e-15
    # ---- 5.3.2 等待＋步伐（W507）
    assert [_s(V.ctrw_bessel(x, 3), 4) for x in (0, 1, 2)] == ["0.2430", "0.1968", "0.1118"]
    assert all(abs(V.ctrw_bessel(x, 3) - V.ctrw_lattice(x, 3)) < 1e-12 for x in (0, 1, 2))
    assert abs(V.ctrw_transform(0.7, 0.4) - V.ex58(0.7, 0.4)) < 1e-12
    assert _s(V.ctrw_lattice(0, 10), 4) == "0.1278" and _s(V.ctrw_gauss(0, 10), 4) == "0.1262"
    assert _s(V.ctrw_lattice(0, 10, "gamma2"), 4) == "0.1805" and _s(V.ctrw_gauss(0, 10, "gamma2"), 4) == "0.1784"
    assert [_s(V.ctrw_var(t), 2) for t in (10, 40)] == ["10.00", "40.00"]
    assert [_s(math.sqrt(V.ctrw_var(t)), 3) for t in (10, 40)] == ["3.162", "6.325"]
    assert [_s(V.ctrw_var(t, "gamma2"), 2) for t in (10, 40)] == ["4.75", "19.75"]
    assert all(abs(V.ctrw_var(t, k) - V.mean_steps(t, k)) < 1e-8 for t in (10, 40) for k in ("poisson", "gamma2"))
    assert _s(V.mortal_total(10), 4) == "0.3679" and _s(math.log(2) / 0.1, 3) == "6.931"
    assert abs(V.mortal_series(10) - V.mortal_total(10)) < 1e-12 and abs(V.mortal_laplace(10) - V.mortal_total(10)) < 1e-8


# =============================================================================
# 表與圖
# =============================================================================
def tab_biased():
    p = F(3, 5)
    P10 = V.biased_pmf(10, p)
    r1 = [r"精確 (5.1) & " + " & ".join(_s(P10[x], 4) for x in (0, 2, 4)) + r"\\",
          r"$2\times$(5.7) & " + " & ".join(_s(2 * V.biased_gauss(float(x), 10, 0.6), 4) for x in (0, 2, 4)) + r"\\"]
    Ps = {N: V.biased_pmf(N, p) for N in (10, 20, 40)}
    mom = {N: V.moments(P) for N, P in Ps.items()}
    r2 = [r"平均 $(2p-1)N$ & " + " & ".join(str(int(mom[N][0])) for N in Ps) + r"\\",
          r"變異數 $4Np(1-p)$ & " + " & ".join(str(float(mom[N][1])) for N in Ps) + r"\\",
          r"峰值（精確） & " + " & ".join(_s(max(P.values()), 4) for P in Ps.values()) + r"\\",
          r"峰值（$2\times$高斯） & " + " & ".join(_s(2 * V.biased_gauss(float(mom[N][0]), N, 0.6), 4) for N in Ps) + r"\\"]
    a = "\n".join([r"\begin{tabular}[t]{@{}lccc@{}}", r"\toprule", r"$N=10$ & $x=0$ & $x=2$ & $x=4$\\", r"\midrule",
                   *r1, r"\bottomrule", r"\end{tabular}"])
    b = "\n".join([r"\begin{tabular}[t]{@{}lccc@{}}", r"\toprule", r"Fig.~5.1 & $N=10$ & $N=20$ & $N=40$\\", r"\midrule",
                   *r2, r"\bottomrule", r"\end{tabular}"])
    return a + "\n\\qquad\n" + b


def _fr(x):
    x = F(x)
    return "0" if x == 0 else (str(x.numerator) if x.denominator == 1 else rf"\frac{{{x.numerator}}}{{{x.denominator}}}")


def tab_lr():
    Q = F(3, 4)
    rows = []
    for N in range(5):
        d = V.persistent_dp(N, Q)
        cells = [(f"$({_fr(d[x][0])},{_fr(d[x][1])})$" if x in d else "") for x in range(-4, 5)]
        rows.append(f"{N} & " + " & ".join(cells) + r"\\")
    return "\n".join([r"{\small\setlength{\tabcolsep}{3pt}\begin{tabular}{@{}c" + "c" * 9 + "@{}}", r"\toprule",
                      r"$N\backslash x$ & " + " & ".join(f"${x}$" for x in range(-4, 5)) + r"\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}}"])


def tab_pvar():
    Ns = (10, 100, 1000)
    rows = []
    for Q, lab, d in ((0.75, r"\tfrac34", 2), (0.25, r"\tfrac14", 3)):
        rows.append(rf"$Q={lab}$ 精確 & " + " & ".join(_s(V.persistent_var(N, Q), d) for N in Ns) + r"\\")
        rows.append(r"\quad $NQ/(1-Q)$ & " + " & ".join(
            (str(int(V.persistent_var_asym(N, Q))) if Q == 0.75 else _s(V.persistent_var_asym(N, Q), 3)) for N in Ns)
            + r"\\")
    sm = [r"$Q$ & $N=1$ & $2$ & $3$ & $4$ & $n=\frac{Q}{1-Q}$\\", r"\midrule"]
    for Q in (F(3, 4), F(1, 2), F(1, 4)):
        sm.append(f"${_fr(Q)}$ & " + " & ".join(f"${_fr(V.persistent_var(N, Q))}$" for N in range(1, 5))
                  + f" & ${_fr(V.persistence_length(Q))}$" + r"\\")
    a = "\n".join([r"\begin{tabular}[t]{@{}cccccc@{}}", r"\toprule", *sm, r"\bottomrule", r"\end{tabular}"])
    b = "\n".join([r"\begin{tabular}[t]{@{}lccc@{}}", r"\toprule", r"變異數 & $N=10$ & $N=100$ & $N=1000$\\", r"\midrule",
                   *rows, r"\bottomrule", r"\end{tabular}"])
    return a + "\n\\qquad\n" + b


def fig_tele():
    out = [r"\begin{tikzpicture}[x=1.3cm,y=7.5cm,>=Stealth]",
           r"\draw[->] (0,0) -- (9.5,0) node[right,font=\scriptsize] {$x$};",
           r"\draw[->] (0,0) -- (0,0.46) node[above,font=\scriptsize] {$\Sigma(x,t)$};",
           r"\foreach \x in {0,...,9} {\draw (\x,0) -- (\x,-0.008); \node[below,font=\scriptsize] at (\x,-0.006) {$\x$};}",
           r"\foreach \y in {0.1,0.2,0.3,0.4} {\draw (0,\y) -- (-0.08,\y); \node[left,font=\scriptsize] at (-0.08,\y) {$\y$};}"]
    xs = np.linspace(0, 9, 61)
    gs = V.tele_gauss(xs, 8.0)
    out.append(r"\draw[gray,thick,dashed] " + " -- ".join(f"({x:.3f},{y:.5f})" for x, y in zip(xs, gs)) + ";")
    for t in TS:
        xs = np.linspace(0, t * (1 - 1e-9), 70)
        ys = V.tele_cont(xs, float(t))
        pts = " -- ".join(f"({x:.4f},{y:.5f})" for x, y in zip(xs, ys))
        c = TCOL[t]
        out += [rf"\draw[{c},thick] {pts} -- ({t},0);",
                rf"\draw[->,{c},very thick] ({t},0) -- ({t},0.355);",
                rf"\node[above,font=\scriptsize,{c},align=center,inner sep=1pt] at ({t},0.355) "
                rf"{{$t={t}$\\${_g(V.tele_spike(t))}$}};"]
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_tele():
    rows = [(r"$\Sigma(0,t)$（(5.61)）", [_s(V.tele_cont(0.0, float(t)), 4) for t in TS]),
            (r"高斯 $1/\sqrt{4\pi Dt}$（$D=\frac12$）", [_s(V.tele_gauss(0.0, float(t)), 4) for t in TS]),
            (r"前緣後方 $\Sigma(vt^-,t)=\frac{e^{-t}}2(1+\frac t2)$", [_g(V.tele_cont(t * (1 - 1e-9), float(t))) for t in TS]),
            (r"每個尖峰 $\frac12e^{-t}$", [_g(V.tele_spike(t)) for t in TS]),
            (r"連續部分總量 $1-e^{-t}$", [_s(1 - math.exp(-t), 4) for t in TS]),
            (r"變異數 $t-\frac{1-e^{-2t}}2$", [_s(V.tele_var(t), 4) for t in TS])]
    body = [f"{a} & " + " & ".join(b) + r"\\" for a, b in rows]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$v=\gamma=1$ & $t=1$ & $t=3$ & $t=5$ & $t=8$\\", r"\midrule", *body,
                      r"\bottomrule", r"\end{tabular}"])


def fig_rt3():
    ymax = 1.5
    out = [r"\begin{tikzpicture}[x=1.3cm,y=1.9cm,>=Stealth]",
           r"\draw[->] (0,0) -- (8.5,0) node[right,font=\scriptsize] {$r$};",
           rf"\draw[->] (0,0) -- (0,{ymax + 0.25}) node[above,font=\scriptsize] {{$P_3(r,t)=4\pi r^2\rho$}};",
           r"\foreach \x in {0,...,8} {\draw (\x,0) -- (\x,-0.03); \node[below,font=\scriptsize] at (\x,-0.02) {$\x$};}",
           r"\foreach \y in {0.5,1,1.5} {\draw (0,\y) -- (-0.08,\y); \node[left,font=\scriptsize] at (-0.08,\y) {$\y$};}",
           r"\draw[gray!30] (0,0.5) -- (8.2,0.5) (0,1) -- (8.2,1) (0,1.5) -- (8.2,1.5);",
           r"\begin{scope}", rf"\clip (0,0) rectangle (8.3,{ymax});"]
    for t in TS:
        ph = np.linspace(2e-3, math.asin(0.99), 110)
        rs = t * np.sin(ph)
        ys = [V.rt3_radial(float(r), float(t)) for r in rs]
        pts = " -- ".join(f"({r:.4f},{min(y, 3.0):.4f})" for r, y in zip(rs, ys))
        out.append(rf"\draw[{TCOL[t]},thick] {pts};")
    out.append(r"\end{scope}")
    for t in TS:
        c = TCOL[t]
        out += [rf"\draw[->,{c},very thick] ({t},0) -- ({t},1.3);",
                rf"\node[above,font=\scriptsize,{c},align=center,inner sep=1pt,fill=white] at ({t},1.3) "
                rf"{{$t={t}$\\${_g(math.exp(-t))}$}};"]
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_rt3():
    rows = []
    for t in TS:
        rows.append(f"$t={t}$ & " + " & ".join(_s(V.rt3_radial(f * t, float(t)), 4) for f in (0.25, 0.5, 0.75, 0.9))
                    + f" & {_g(math.exp(-t))}" + r"\\")
    return "\n".join([r"\begin{tabular}[t]{@{}lccccc@{}}", r"\toprule",
                      r"$P_3$ & $r=\frac t4$ & $\frac t2$ & $\frac{3t}4$ & $0.9t$ & 球殼\\", r"\midrule",
                      *rows, r"\bottomrule", r"\end{tabular}"])


def tab_dims():
    mc = {d: _mc(d) for d in (1, 2, 3)}
    th = [V.radial_fraction(d, 0, 1.5, 3) for d in (1, 2, 3)]
    rows = [r"沒轉過：$e^{-3}$ & " + " & ".join([_s(math.exp(-3), 4)] * 3) + r"\\",
            r"沒轉過：模擬 & " + " & ".join(_s(mc[d][1], 4) for d in (1, 2, 3)) + r"\\",
            r"$r<1.5$：公式 & " + " & ".join(_s(x, 4) for x in th) + r"\\",
            r"$r<1.5$：模擬 & " + " & ".join(_s(mc[d][0], 4) for d in (1, 2, 3)) + r"\\",
            r"前緣 & " + _g(V.tele_cont(3 * (1 - 1e-9), 3.0)) + r" & $1/\sqrt{\ }$ 發散 & $\ln$ 發散\\"]
    return "\n".join([r"\begin{tabular}[t]{@{}lccc@{}}", r"\toprule",
                      r"$t=3$ & 一維 & 二維 & 三維\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_ctrw():
    ms = range(7)
    r1 = r"固定速率 $p_0=1$ & " + " & ".join(_s(V.steps_prob(m, 3), 4) for m in ms) + f" & {_s(V.mean_steps(3), 0)}" + r"\\"
    r2 = (r"兩段等待 & " + " & ".join(_s(V.steps_prob(m, 3, "gamma2"), 4) for m in ms)
          + f" & {_s(V.mean_steps(3, 'gamma2'), 4)}" + r"\\")
    a = "\n".join([r"{\small\setlength{\tabcolsep}{4pt}\begin{tabular}{@{}l" + "c" * 8 + "@{}}", r"\toprule",
                   r"$t=3$ 走了 $m$ 步 & " + " & ".join(f"${m}$" for m in ms) + r" & 平均\\", r"\midrule", r1, r2,
                   r"\bottomrule", r"\end{tabular}}"])
    K = ("poisson", "gamma2")
    rows = [r"$P(0,10)$：精確 & " + " & ".join(_s(V.ctrw_lattice(0, 10, k), 4) for k in K) + r"\\",
            r"$P(0,10)$：高斯 $\sqrt{B(0)/2\pi t}$ & " + " & ".join(_s(V.ctrw_gauss(0, 10, k), 4) for k in K) + r"\\",
            r"變異數（精確）$t=10\to40$ & " + " & ".join(
                f"{_s(V.ctrw_var(10, k), 2)} $\\to$ {_s(V.ctrw_var(40, k), 2)}" for k in K) + r"\\",
            r"寬度 $t=10\to40$ & " + " & ".join(
                f"{_s(math.sqrt(V.ctrw_var(10, k)), 3)} $\\to$ {_s(math.sqrt(V.ctrw_var(40, k)), 3)}" for k in K) + r"\\"]
    b = "\n".join([r"\begin{tabular}{@{}lcc@{}}", r"\toprule",
                   r"一維格子 $C=\frac12$ & 固定速率 $B(0)=1$ & 兩段等待 $B(0)=2$\\", r"\midrule", *rows,
                   r"\bottomrule", r"\end{tabular}"])
    return a + "\n\\\\[6pt]\n" + b


ERRATA = [
    (r"p.~96，(5.3)", r"$l(2p-1)\frac{\partial P}{\partial y}$", r"$l(1-2p)\frac{\partial P}{\partial y}$"),
    (r"p.~103，(5.37)", r"$\lambda^2-2Q\cos ql\,\lambda+(1-2Q)=0$", r"$\lambda^2-2Q\cos ql\,\lambda-(1-2Q)=0$"),
    (r"p.~107，(5.57) 前", r"when $x<t$", r"$x<vt$"),
    (r"p.~107，(5.60)(5.61)", r"$\frac12\big(2\gamma+\frac{\partial}{\partial t}\big)$",
     r"$\frac{1}{2v}\big(2\gamma+\frac{\partial}{\partial t}\big)$"),
    (r"p.~112，(5.83) 後", r"$b\to i\gamma$", r"$b\to i\gamma/v$"),
    (r"p.~112，(5.84)", r"$2\gamma\sqrt{\pi/2}\,e^{-\gamma t}\cosh(\gamma\sqrt{v^2t^2-r^2})$",
     r"$2\gamma\,e^{-\gamma t}\cosh(\gamma\sqrt{v^2t^2-r^2}/v)$"),
    (r"p.~113，(5.88)", r"$\frac{2\pi}{x}\delta(x-y)$", r"$\frac1x\delta(x-y)$"),
    (r"p.~113，(5.89)", r"兩項前面各有 $2\pi$", r"兩個 $2\pi$ 都去掉"),
    (r"p.~115，(5.96) 前", r"$\frac1{vq}\big(\frac\pi2-i\ln((1+x)(1-x))\big)$",
     r"$\frac1{vq}\big(\frac\pi2-\frac i2\ln\frac{1+x}{1-x}\big)$"),
    (r"p.~117，(5.106)", r"$e^{-\gamma ut\cot u}$", r"$e^{+\gamma ut\cot u}$"),
    (r"p.~124，(5.133) 後", r"width proportional to $t$", r"變異數 $\propto t$（寬度 $\propto\sqrt t$）"),
]


def tab_errata():
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("FiveTabBiased", tab_biased),
    ("FiveTabLR", tab_lr),
    ("FiveTabPVar", tab_pvar),
    ("FiveFigTele", fig_tele),
    ("FiveTabTele", tab_tele),
    ("FiveFigRT", fig_rt3),
    ("FiveTabRT", tab_rt3),
    ("FiveTabDims", tab_dims),
    ("FiveTabCTRW", tab_ctrw),
    ("FiveTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch5()
    parts = ["% figs.tex —— 由 ch5/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch5/make_figs.py（第 5 章，動畫 W501–W507）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
