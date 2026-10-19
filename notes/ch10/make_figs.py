r"""
make_figs.py —— 隨機漫步導讀第 10 章（Rudnick & Gaspari《Elements of the Random Walk》第 10 章，動畫 W1001–W1005）的圖與表

在 notes/ 或 notes/ch10/ 裡執行：
    python ch10/make_figs.py
會先檢查 ch10/part10.tex 寫到的每一個數字（check_numbers_ch10；最後再把 part10.tex 裡出現的每一個小數和檢查過的數比對一遍），
再寫出 ch10/figs.tex（part10.tex 第一行 \input{ch10/figs}）。圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算）。
數學一律由 codes/utils/walk_cond.py 算出（z_c 用 walk_crit.ZC_2D）；例子與動畫 W1001–W1005 相同。巨集名稱都以 \Ten 開頭
（不和前幾章的 \Fig…、\Two…、…、\Nine… 撞名）；part10.tex 的 \label 用 ten-…。顏色沿用導言區：cBar 精確值、cWalk 路線、
cExact 更正、cGauss 近似、cWrong 筆誤。環面的列舉、週期盒子的 FFT 都由 walk_cond 的 lru_cache 快取，同一次執行只算一次（全部約 3 秒）。
"""
from __future__ import annotations

import math
import re
import sys
import warnings
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import special

HERE = Path(__file__).resolve().parent            # notes/ch10
sys.path.insert(0, str(HERE.parent.parent / "codes"))
warnings.filterwarnings("ignore")                  # scipy 的 IntegrationWarning 不影響結果

from utils import walk_cond as WK  # noqa: E402   （套件模組：walk_cond 的相對 import 才找得到 walk_saw、walk_crit、walk_rec、walk_gf）
from utils import walk_crit as WC  # noqa: E402

OUT = HERE / "figs.tex"
TEX = HERE / "part10.tex"
F = Fraction
OK = set()                                         # 檢查過的小數字串（part10.tex 的小數都要在這裡）
MODES = ("beta", "gamma", "alpha", "nu")


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣），並記進 OK。"""
    s = f"{float(x):.{d}f}"
    OK.add(s.lstrip("-"))
    return s


def _e(x, d=1):
    """科學記號 a.b×10^c（walk_cond 的 .1e 寫法），尾數記進 OK。"""
    mant, ex = f"{float(x):.{d}e}".split("e")
    OK.add(mant.lstrip("-"))
    return rf"{mant}\times10^{{{int(ex)}}}"


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


def _matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def _matvec(A, v):
    return [sum(a * b for a, b in zip(row, v)) for row in A]


def _eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


# =============================================================================
# 數字檢查（part10.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch10():
    OK.update({"0.2", "0.6", "1.0", "0.5", "0.1", "0.01", "7.57", "7.5", "9.28", "9.76", "9.77"})
    zc = WC.ZC_2D
    assert _s(zc, 4) == "0.3791"
    # ---- W1001 關在盒子裡：環面上的自迴避漫步
    c3, c4, c5 = (WK.torus_counts(L) for L in WK.TORI)
    assert list(c3) == [1, 4, 12, 32, 80, 160, 256, 296, 168] and sum(c3) == 1009
    assert list(c4) == [1, 4, 12, 36, 96, 264, 624, 1536, 3072, 6240, 9840, 15912, 17472, 19200, 10656, 5712]
    assert sum(c4) == 90677 and len(c5) == 25 and c5[5] == 280 and c5[-1] == 236456 and sum(c5) == 22094601
    assert list(WK.SAW_INF[:16]) == [1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932, 881500,
                                     2374444, 6416596]
    for L, c in zip(WK.TORI, (c3, c4, c5)):
        assert WK.torus_first_diff(L) == (L, 4)
        assert list(c[:L]) == list(WK.SAW_INF[:L]) and WK.SAW_INF[L] - c[L] == 4
    assert list(WK.torus_counts_plain(4)) == list(c4)                       # 第二種算法：不用對稱的 DFS
    assert all(WK.torus_counts_brute(4, N) == c4[N] for N in range(7))     # 第三種：4^N 條全部列舉
    assert WK.GRAY_Q4 == 16 * c4[15] == 91392
    ex = {k: WK.torus_example(k) for k in ("short", "mid", "full")}
    assert [(e["N"], e["steps"], e["wraps"]) for e in ex.values()] == [(3, "WNE", 0), (8, "WSSEENNE", 2),
                                                                       (15, "WSSEEENNNEEESSW", 3)]
    assert all(e["start"] == (1, 1) for e in ex.values()) and ex["mid"]["lift"][-1] == (3, 1)
    assert len(set(ex["full"]["sites"])) == 16
    g02, n02, _ = WK.torus_moments(4, 0.2)
    assert _close(g02, sum(c * 0.2 ** N for N, c in enumerate(c4)), 1e-12)
    assert _close(g02 * n02, sum(N * c * 0.2 ** N for N, c in enumerate(c4)), 1e-12)
    assert (_s(g02, 4), _s(g02 * n02, 4)) == ("2.8782", "4.1445")
    assert [_s(4 * 0.2, 1), _s(12 * 0.04, 2), _s(36 * 0.008, 3)] == ["0.8", "0.48", "0.288"]
    assert [_s(2 * 12 * 0.04, 2), _s(3 * 36 * 0.008, 3)] == ["0.96", "0.864"]
    _all([WK.torus_mean_N(4, z) for z in WK.MEAN_N_Z], 4, ("1.4400", "4.8899", "11.7409"))
    assert all(WK.book101_mean_N(4, z) == 0 for z in WK.MEAN_N_Z)
    # ---- W1002 平均場
    assert (WK.MF_N0, WK.MF_R, WK.MF_U, WK.MF_N0S) == (100, 0.5, 1.0, (100, 1000, 10000))
    _all([WK.gamma_n(n) for n in WK.MF_NS], 4, ("1.7183", "1.7634", "1.8131", "1.8624", "1.8676"))
    assert _close(WK.gamma_n(1), WK.gamma_n1_direct(), 1e-8) and _close(WK.gamma_n(0.1), WK.gamma_n(0.1, method="pcf"), 1e-6)
    gmf = WK.gamma_mf(0.5)
    assert (_s(gmf, 4), _s(WK.gamma_mf(0.5, book=True), 4)) == ("1.8682", "0.9341")
    gs, g0, gbook = WK.gamma_gauss_check()
    assert all(_close(x, 2, 1e-9) for x in gs) and _close(g0, 2, 1e-9) and _close(gbook, 1, 1e-9)
    s1, s2 = WK.s_integral(1.0), WK.s_integral(-1.0)
    assert _s(s1.imag, 4) == "5.6853" and _close(s1.imag, 2 * math.pi * math.exp(-0.1), 1e-6) and abs(s2) < 1e-8
    a, b, c = math.sqrt(math.pi / 100), math.exp(6.25), special.erfc(2.5)       # 閉式 I = √(π/N0u) e^{N0r²/4u} erfc(r√(N0/u)/2)
    assert (_s(a, 4), _s(b, 2), _e(c, 3), _s(a * b * c, 5)) == ("0.1772", "518.01", r"4.070\times10^{-4}", "0.03736")
    assert _close(a * b * c, WK.mf_I(0.5)) and _close(WK.mf_I(0.5), WK.mf_I(0.5, method="quad"), 1e-9)
    assert _close(50 * a * b * c, gmf) and _close(0.25 / 4 * 100, 6.25)
    _all([WK.gamma_mf(0.5, 1, N0) for N0 in WK.MF_N0S], 4, ("1.8682", "1.9844", "1.9984"))
    _all([math.exp(WK.ln_gamma_mf(-0.5, 1, N0) - WK.ln_gamma_mf_asym(-0.5, 1, N0)) for N0 in WK.MF_N0S], 4,
         ("0.9998", "1.0000", "1.0000"))
    _all([WK.ln_gamma_mf(-0.5, 1, N0) / N0 for N0 in WK.MF_N0S], 4, ("0.0912", "0.0665", "0.0630"))
    assert _s(0.25 / 4, 4) == "0.0625"
    gneg = WK.gamma_mf(-0.5)
    assert _s(gneg, 1) == "9179.7" and WK.MF_NS_NEG == (0.01, 0.001, 0.0001, 1e-6)
    _all([WK.gamma_n(n, -0.5) / gneg for n in WK.MF_NS_NEG], 4, ("0.3299", "0.8311", "0.9801", "0.9998"))
    _all(WK.fn3_check(), 4, ("1.8682", "1.8682", "1.8682", "0.9341", "0.9341"))
    # ---- W1003 自旋波
    assert (WK.WALL_R, WK.WALL_U, WK.WALL_X) == (-1.0, 1.0, (1, 2, 4))
    _all([WK.wall_M(x) for x in WK.WALL_X], 4, ("0.6089", "0.8884", "0.9930"))
    assert _s(math.sqrt(2 / abs(WK.WALL_R)), 4) == "1.4142"
    e_bvp, e_res = WK.wall_errors()
    assert _e(e_bvp) == r"2.8\times10^{-6}" and e_res < 1e-7
    assert WK.ex101_grad_check() < 1e-7 and WK.ex102_coeffs() < 1e-11
    l1, l2 = WK.sw_lambdas()
    assert (WK.SW_R, WK.SW_U, WK.SW_Q2) == (-1, 1, 1) and (l1, l2) == (F(1, 2), F(3, 2))
    M = WK.SW_M
    X, Xi = WK.sw_X(), WK.sw_Xinv()
    assert M == (F(2, 3), F(2, 3), F(1, 3)) and sum(m * m for m in M) == 1
    assert [[18 * v for v in row] for row in X] == [[17, 8, 4], [8, 17, 4], [4, 4, 11]]
    assert [[27 * v for v in row] for row in Xi] == [[38, -16, -8], [-16, 38, -8], [-8, -8, 50]]
    I3 = _eye(3)
    assert all(X[i][j] == l1 * I3[i][j] + M[i] * M[j] for i in range(3) for j in range(3))
    assert all(Xi[i][j] == 2 * I3[i][j] - F(4, 3) * M[i] * M[j] for i in range(3) for j in range(3))
    assert _matmul(X, Xi) == I3 == _matmul(Xi, X)
    assert _matvec(X, M) == [l2 * m for m in M]
    for v in ((1, -1, 0), (1, 1, -4)):
        assert sum(a * m for a, m in zip(v, M)) == 0 and _matvec(X, v) == [l1 * a for a in v]
    assert sum(X[i][i] for i in range(3)) == F(45, 18) == 2 * l1 + l2
    errs, ok = WK.sw_check()
    assert ok and max(errs.values()) < 1e-11
    tbl = WK.mode_table()
    assert [t[:4] for t in tbl] == [(3, 2, 1, F(3, 8)), (2, 1, 1, F(3, 4)), (1, 0, 1, F(3, 2)), (0, -1, 1, F(3))]
    _all([t[4] for t in tbl], 4, ("1.6330", "1.1547", "0.8165", "0.5774"))
    assert _close(tbl[-1][4], math.sqrt(l1 / l2)) and _s(math.sqrt(l2 / l1), 4) == "1.7321"
    gnum, gplus, gminus = WK.gauss_n2_check()
    _all((gnum, gminus), 4, ("3.9891", "3.2988"))
    assert _close(gnum, gplus, 1e-8) and WK.GAUSS_M == (0.6, 0.8) and WK.GAUSS_H == (0.5, 0.0)
    M2, h = (F(3, 5), F(4, 5)), (F(1, 2), F(0))
    Xi2 = WK.sw_Xinv(M=M2)
    hXh = sum(h[i] * Xi2[i][j] * h[j] for i in range(2) for j in range(2))
    assert Xi2[0][0] == F(38, 25) and hXh == F(19, 50) and (_s(hXh, 2), _s(hXh / 4, 3)) == ("0.38", "0.095")
    pref = math.pi / math.sqrt(l1 * l2)
    assert _s(pref, 4) == "3.6276" and _close(pref * math.exp(float(hXh) / 4), gplus, 1e-9)
    assert _close(pref * math.exp(-float(hXh) / 4), gminus, 1e-9)
    # ---- W1004 相關洞
    assert _e(WK.case1_check()) == r"1.2\times10^{-11}"
    Dax = WK.D_axis()
    assert (WK.HOLE_L, WK.HOLE_R, WK.HOLE_D, len(Dax)) == (64, -1.0, (0, 1, 2, 4, 8, 16, 32), 33)
    _all([Dax[d] for d in WK.HOLE_D], 4, ("0.1083", "0.0613", "0.0359", "0.0167", "0.0065", "0.0016", "-0.0001"))
    _all([1 - Dax[d] for d in WK.HOLE_D], 4, ("0.8917", "0.9387", "0.9641", "0.9833", "0.9935", "0.9984", "1.0001"))
    assert all(np.diff(Dax) < 0) and round(100 * Dax[0]) == 11 and round(100 * Dax[1]) == 6
    assert all(_close(Dax[d], WK.D_direct(64, -1.0, (d, 0, 0)), 1e-9) for d in (0, 1, 4))     # 第二種算法
    c0, cc = WK.coulomb_fit()
    assert (_s(c0, 4), _s(cc, 4), _s(1 / (4 * math.pi), 4), _s(4 * math.pi * cc, 4)) == ("0.0035", "0.0801", "0.0796", "1.0064")
    d46 = WK.D_book46()
    assert _s(d46[0], 4) == "0.0541" and _close(2 * d46[0], Dax[0])
    am_book, am = WK.book45_argmax()
    assert _close(am_book, 0.5) and _close(am, 1.0)
    _all([WK.D_cont(p) for p in WK.CONT_RHO], 4, ("0.0602", "0.0374", "0.0198", "0.0099"))
    assert all(_close(WK.D_cont(p), WK.D_cont_quad(p), 1e-8) for p in WK.CONT_RHO)
    _all([WK.D_inf((p, 0, 0)) for p in WK.CONT_RHO], 4, ("0.0648", "0.0394", "0.0202", "0.0100"))
    _all([1 / (4 * math.pi * p) for p in WK.CONT_RHO], 4, ("0.0796", "0.0398", "0.0199", "0.0099"))
    assert _s(WK.SC_EWALD, 6) == "2.837297" and _s(WK.SC_EWALD / (4 * math.pi * 64), 4) == "0.0035"
    assert all(abs(Dax[p] - WK.D_inf((p, 0, 0)) - WK.box_shift(64, -1.0, (p, 0, 0))) < 5e-6 for p in WK.CONT_RHO)
    assert (WK.VIS_N, WK.VIS_R) == (100, (1, 2, 4, 8, 16))
    vis = [WK.sc_visits(100, R) for R in WK.VIS_R]
    ef = [WK.erfc_form(R, 100) for R in WK.VIS_R]
    assert [_s(v, 4) for v in vis[:4]] + [_s(vis[4], 6)] == ["0.4509", "0.1931", "0.0610", "0.0103", "0.000183"]
    assert [_s(v, 4) for v in ef[:4]] + [_s(ef[4], 6)] == ["0.4121", "0.1744", "0.0586", "0.0100", "0.000174"]
    _all([v / e for v, e in zip(vis, ef)], 3, ("1.094", "1.108", "1.042", "1.028", "1.052"))
    assert [round(100 * (min(v / e for v, e in zip(vis, ef)) - 1)), round(100 * (max(v / e for v, e in zip(vis, ef)) - 1))] == [3, 11]
    assert _close(WK.sc_visits_dp(20, 2)[0], WK.sc_visits(20, 2), 1e-12)
    assert [_e(WK.erfc_form_book(R, 100)) for R in (8, 16)] == [r"4.3\times10^{-4}", r"1.7\times10^{-28}"]
    gr = [WK.green_sc(R) for R in WK.VIS_R]
    _all(gr, 4, ("0.5164", "0.2573", "0.1217", "0.0599", "0.0299"))
    _all([g * 2 * math.pi * R / 3 for g, R in zip(gr, WK.VIS_R)], 4, ("1.0815", "1.0779", "1.0198", "1.0041", "1.0009"))
    assert _s(WK.green_sc(0), 4) == "1.5164" and _close(gr[0], WK.green_sc(0) - 1, 1e-6)
    assert _close(WK.green_sc_series(4), gr[2], 2e-3)
    assert _s(math.sqrt(2 * 101 / 3), 1) == "8.2"                     # erfc 的引數 R√(3/(2(N+1))) = 1
    e80 = WK.eq1080_check()                                            # (10.80)：級數＝書上左邊 × z_c^{−N}＝右邊
    assert max(e80) - min(e80) < 1e-9 * max(e80)
    assert _close(*WK.eq1083_check(), 1e-9)                            # (10.82)→(10.83) 的 Fourier 積分＝erfc 閉式
    # ---- W1005 Redner–Reynolds
    assert WK.RHO_Z[1] == zc
    rho_exp = {3: ("0.1391", "0.3374", "0.5284", "0.6823"), 4: ("0.0900", "0.3056", "0.5585", "0.7338"),
               5: ("0.0603", "0.2769", "0.5764", "0.7565")}
    peak_exp = {3: ("0.4443", "0.4289"), 4: ("0.4380", "0.5846"), 5: ("0.4272", "0.7242")}
    for L in WK.TORI:
        _all([WK.torus_rho(L, z) for z in WK.RHO_Z], 4, rho_exp[L])
        assert all(_close(WK.torus_rho(L, z), WK.torus_rho_fd(L, z), 1e-6) for z in WK.RHO_Z)
        zp, hp = WK.torus_peak(L)
        assert (_s(zp, 4), _s(hp, 4)) == peak_exp[L] and _close(WK.torus_var_fd(L, zp), hp, 1e-4)
        assert _close(WK.torus_rho(L, 1e4), F(L * L - 1, L * L), 1e-3)
    assert [F(L * L - 1, L * L) for L in WK.TORI] == [F(8, 9), F(15, 16), F(24, 25)]
    assert WK.MF_RHO_N0 == 10 ** 6 and WK.MF_RHO_R == (-0.5, -0.1, -0.01, 0.1, 0.5)
    rho = [WK.mf_rho(r) for r in WK.MF_RHO_R]
    assert [_s(x, 4) for x in rho[:3]] + [_e(x) for x in rho[3:]] == ["0.2500", "0.0500", "0.0050",
                                                                      r"1.0\times10^{-5}", r"2.0\times10^{-6}"]
    assert all(_close(WK.mf_rho(r), WK.mf_rho_fd(r), 1e-6) for r in (-0.5, -0.1)) and _close(WK.mf_rho_quad(-0.5), 0.25, 1e-9)
    assert _s(math.log(WK.mf_rho(-0.01) / WK.mf_rho(-0.1)) / math.log(0.1), 4) == "1.0000"
    chi = [WK.mf_chi(r) for r in (-0.1, 0.1, 0.5)]
    assert [_s(chi[0], 4), _e(chi[1]), _e(chi[2])] == ["0.5000", r"1.0\times10^{-4}", r"4.0\times10^{-6}"]
    assert _close(WK.mf_chi(-0.1), WK.mf_chi_fd(-0.1), 1e-6)
    assert _s(WK.book55_rho(-0.5), 4) == "-0.2500"
    E = {k: WK.rr_exponents(*WK.RR_CASES[k]) for k in ("mf", "d2", "d3")}
    assert [E["mf"][k] for k in MODES] == [1, 0, 0, F(1, 2)] and [E["d2"][k] for k in MODES] == [F(1, 2)] * 3 + [F(3, 4)]
    assert [_s(E["d3"][k], 3) for k in MODES] == ["0.764", "0.236", "0.236", "0.588"]
    for k, (nu, d) in WK.RR_CASES.items():
        t = WK.table101(2 - d * nu, nu)                                  # Table 10.1 的 (α, ν) 寫法 = (10.73) 的 (ν, d)
        assert all(_close(t[q], E[k][q]) for q in MODES)
        s_, two_a, dnu, dcnu = WK.rr_laws(E[k], d)
        assert _close(s_, 2) and _close(two_a, dnu) and _close(dcnu, 4 * E[k]["nu"])
    assert WK.rr_laws(E["d2"], 2)[1:] == (F(3, 2), F(3, 2), 3)
    assert [_s(x, 3) for x in WK.rr_laws(E["d3"], 3)[1:]] == ["1.764", "1.764", "2.352"]
    assert WK.dc_mean_field() == 4


def check_tex_decimals():
    """part10.tex（去掉註解）裡的每一個小數都必須是上面檢查過的數。"""
    if not TEX.exists():
        print("（尚無 part10.tex，略過小數比對）")
        return
    txt = re.sub(r"(?<!\\)%.*", "", TEX.read_text(encoding="utf-8"))
    txt = re.sub(r"[\d.]+\\linewidth|[\d.]+pt\b|[\d.]+em\b|[\d.]+cm\b|\\arraystretch\}\{[\d.]+\}", "", txt)   # 版面尺寸
    eqno = re.compile(r"10\.\d{1,2}")                    # 本章的式號、節號、圖號、練習號（10.1–10.84）
    bad = sorted({m for m in re.findall(r"\d+\.\d+", txt) if m not in OK and not eqno.fullmatch(m)})
    assert not bad, f"part10.tex 有未檢查的小數：{bad}"


# =============================================================================
# 圖與表
# =============================================================================
ARROW1 = r"postaction={decorate},decoration={markings,mark=at position 0.56 with {\arrow{Stealth[length=2.6pt]}}}"
ARROW2 = (r"postaction={decorate},decoration={markings,mark=at position 0.44 with {\arrow{Stealth[length=2.6pt]}},"
          r"mark=at position 0.62 with {\arrow{Stealth[length=2.6pt]}}}")


def _pt(x0, s):
    return lambda x, y: f"({x0 + s * x:.3f},{s * y:.3f})"


def _torus_panel(name, x0, s):
    """L×L 環面（格點畫在格子中心，盒子 [−½, L−½]²）：對邊的同向小箭頭＝黏在一起；盒子裡的步（橙）與越界的半步（藍）。"""
    ex = WK.torus_example(name)
    L = ex["L"]
    P = _pt(x0, s)
    lo, hi = -0.5, L - 0.5
    out = [r"\draw[gray!70," + ARROW1 + "] " + P(lo, lo) + " -- " + P(lo, hi) + ";",
           r"\draw[gray!70," + ARROW1 + "] " + P(hi, lo) + " -- " + P(hi, hi) + ";",
           r"\draw[gray!70," + ARROW2 + "] " + P(lo, lo) + " -- " + P(hi, lo) + ";",
           r"\draw[gray!70," + ARROW2 + "] " + P(lo, hi) + " -- " + P(hi, hi) + ";"]
    out += [rf"\fill[gray!55] {P(x, y)} circle (0.7pt);" for y in range(L) for x in range(L)]
    for a, b, cross in ex["segments"]:
        out.append(rf"\draw[{'cBar' if cross else 'cWalk'},very thick,line cap=round] {P(*a)} -- {P(*b)};")
    ex_, ey = (v % L for v in ex["lift"][-1])
    out.append(rf"\fill[cExact] {P(*ex['start'])} circle (1.7pt);")
    out.append(rf"\draw[cWalk!70!black,fill=white,thick] {P(ex_, ey)} circle (1.5pt);")
    return out


def _lift_panel(x0, s):
    """'mid' 攤開到無限平面（walk_cond 的 lift）：虛線是盒子和它下面的複本（剪到 y ≥ −1.75）。"""
    ex = WK.torus_example("mid")
    L = ex["L"]
    P = _pt(x0, s)
    out = [r"\begin{scope}", rf"\clip {P(-0.75, -1.75)} rectangle {P(L - 0.25, L - 0.25)};",
           rf"\draw[gray!60,densely dashed] {P(-0.5, -0.5)} rectangle {P(L - 0.5, L - 0.5)};",
           rf"\draw[gray!60,densely dashed] {P(-0.5, -L - 0.5)} rectangle {P(L - 0.5, -0.5)};"]
    out += [rf"\fill[gray!40] {P(x, y)} circle (0.7pt);" for y in range(-1, L) for x in range(L)]
    out.append(r"\end{scope}")
    out.append(r"\draw[cWalk,very thick,line join=round,line cap=round] " + " -- ".join(P(*p) for p in ex["lift"]) + ";")
    out.append(rf"\fill[cExact] {P(*ex['lift'][0])} circle (1.7pt);")
    out.append(rf"\draw[cWalk!70!black,fill=white,thick] {P(*ex['lift'][-1])} circle (1.5pt);")
    return out


def fig_torus():
    """4×4 環面上的三條例子（walk_cond.torus_example 的 segments）＋ 'mid' 攤開到平面。"""
    s, gap = 0.38, 2.85
    lab = {"short": r"$N=3$：\texttt{WNE}\\不越界",
           "mid": r"$N=8$：\texttt{WSSEENNE}\\越界 2 次",
           "full": rf"$N=15$：走滿 16 格\\越界 {WK.torus_example('full')['wraps']} 次"}
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]"]
    for k, name in enumerate(("short", "mid", "full")):
        out += _torus_panel(name, gap * k, s)
        out.append(rf"\node[below,align=center,inner sep=1pt] at ({gap * k + 1.5 * s:.3f},{-1.85 * s:.3f}) {{{lab[name]}}};")
    x0 = 3 * gap + 0.15
    out += _lift_panel(x0, s)
    end = WK.torus_example("mid")["lift"][-1]
    out.append(rf"\node[below,align=center,inner sep=1pt] at ({x0 + 1.5 * s:.3f},{-1.85 * s:.3f}) "
               rf"{{$N=8$ 攤開到平面\\終點 $({end[0]},{end[1]})$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_torus():
    c4, inf = WK.torus_counts(4), WK.SAW_INF

    def cell(v, N):
        if N <= 3:
            return rf"{{\color{{cExact}}{v}}}"
        return rf"\textbf{{{v}}}" if N == 4 else str(v)

    out = [r"\begin{tabular}{@{}l*{8}{r}@{}}", r"\toprule"]
    for k, Ns in enumerate((range(8), range(8, 16))):
        if k:
            out.append(r"\midrule")
        out += [r"$N$ & " + " & ".join(map(str, Ns)) + r"\\", r"\midrule",
                r"$4\times4$ 環面 $c_N$ & " + " & ".join(cell(c4[N], N) for N in Ns) + r"\\",
                r"無限平面 $c_N$ & " + " & ".join(cell(inf[N], N) for N in Ns) + r"\\"]
    out += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(out)


def fig_mf():
    """f(λ) = e^{−N0(rλ/2+uλ²/4)}（N0 = 100、u = 1），r = 0.5（藍，從 0 遞減）與 r = −0.5（橙，λ* = |r|/u 有峰），各自正規化到峰值 1。"""
    lam_s = 0.5 / WK.MF_U
    xp = [float(x) for x in np.linspace(0, 0.3, 46)] + [0.45, 0.6, 0.8, 1.0]
    xn = [float(x) for x in np.linspace(0, 1.0, 60)]
    yp = [WK.mf_f(x, 0.5) / WK.mf_f(0.0, 0.5) for x in xp]
    yn = [WK.mf_f(x, -0.5) / WK.mf_f(lam_s, -0.5) for x in xn]
    assert len(xp) <= 60 and len(xn) <= 60 and max(yn) <= 1 + 1e-12 and _close(yp[0], 1)
    out = [r"\begin{tikzpicture}[x=5.0cm,y=1.9cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (0,0) -- (1.12,0) node[right] {$\lambda$};",
           r"\draw[->] (0,0) -- (0,1.2) node[above] {$f(\lambda)/$峰值};",
           r"\foreach \x in {0,0.5,1} {\draw (\x,0) -- (\x,-0.03) node[below] {$\x$};}",
           r"\foreach \y in {0.5,1} {\draw (0,\y) -- (-0.015,\y) node[left] {$\y$};}",
           rf"\draw[gray!60,densely dotted] ({lam_s},0) -- ({lam_s},1.06);"]
    for xs, ys, col in ((xp, yp, "cBar"), (xn, yn, "cWalk")):
        out.append(rf"\draw[{col},thick] plot coordinates {{{' '.join(f'({a:.4f},{b:.4f})' for a, b in zip(xs, ys))}}};")
    out.append(r"\node[cBar,anchor=west,align=left] at (0.035,0.66) {$r=0.5$\\（$\vec M\approx0$）};")
    out.append(r"\node[cWalk,anchor=west,align=left] at (0.66,0.8) {$r=-0.5$：峰在\\$\lambda^*=|r|/u$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_n():
    gn = [_s(WK.gamma_n(n), 4) for n in WK.MF_NS] + [_s(WK.gamma_mf(0.5), 4)]
    gneg = WK.gamma_mf(-0.5)
    rn = [_s(WK.gamma_n(n, -0.5) / gneg, 4) for n in WK.MF_NS_NEG]
    return "\n".join([r"\begin{tabular}{@{}l*{6}{c}@{}}", r"\toprule",
                      r"$r=0.5$，$n=$ & 3 & 2 & 1 & 0.1 & 0.01 & $\to0$\\", r"\midrule",
                      r"$\Gamma_n$ & " + " & ".join(gn) + r"\\", r"\midrule",
                      r"$r=-0.5$，$n=$ & $10^{-2}$ & $10^{-3}$ & $10^{-4}$ & $10^{-6}$ & & \\", r"\midrule",
                      r"$\Gamma_n/\Gamma_{MF}$ & " + " & ".join(rn) + r" & & \\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_asym():
    Ns = WK.MF_N0S
    g = [_s(WK.gamma_mf(0.5, 1, N0), 4) for N0 in Ns]
    ra = [_s(math.exp(WK.ln_gamma_mf(-0.5, 1, N0) - WK.ln_gamma_mf_asym(-0.5, 1, N0)), 4) for N0 in Ns]
    ln = [_s(WK.ln_gamma_mf(-0.5, 1, N0) / N0, 4) for N0 in Ns]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$N_0$ & 100 & 1000 & $10^4$ & 極限\\", r"\midrule",
                      r"$r=0.5$：$\Gamma_{MF}$ & " + " & ".join(g) + r" & $\frac1r=2$\\",
                      r"$r=-0.5$：精確 $\div$ 漸近 & " + " & ".join(ra) + r" & 1\\",
                      r"$r=-0.5$：$\ln\Gamma_{MF}/N_0$ & " + " & ".join(ln) + rf" & $\frac{{r^2}}{{4u}}={_s(0.25 / 4, 4)}$\\",
                      r"\bottomrule", r"\end{tabular}"])


def fig_wall():
    """牆邊的狀態方程解 M = √(|r|/u) tanh(x√(|r|/2))（r = −1、u = 1），60 點；切線在癒合長度 √(2/|r|) 碰到 M0。"""
    xs = [float(x) for x in np.linspace(0, 6, 60)]
    hl = math.sqrt(2 / abs(WK.WALL_R))
    out = [r"\begin{tikzpicture}[x=0.8cm,y=1.9cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\fill[gray!22] (-0.32,0) rectangle (0,1.18);",
           r"\draw[gray!70!black,thick] (0,0) -- (0,1.18);",
           r"\node[gray!60!black,font=\tiny] at (-0.16,0.12) {牆};",
           r"\draw[->] (0,0) -- (6.5,0) node[right] {$x$};",
           r"\foreach \x in {2,4,6} {\draw (\x,0) -- (\x,-0.03) node[below] {$\x$};}",
           r"\foreach \y in {0.5,1} {\draw (-0.32,\y) -- (-0.4,\y) node[left] {$\y$};}",
           r"\node[left] at (-0.32,1.12) {$M$};",
           r"\draw[cExact,densely dashed] (0,1) -- (6.2,1) node[above left] {$M_0=\sqrt{|r|/u}=1$};",
           rf"\draw[gray!60,densely dotted] (0,0) -- ({hl:.4f},1) -- ({hl:.4f},0);",
           rf"\node[anchor=west,gray!60!black] at ({hl:.4f},0.2) {{癒合長度 $\sqrt{{2/|r|}}={_s(hl, 4)}$}};",
           rf"\draw[cBar,thick] plot coordinates {{{' '.join(f'({x:.3f},{WK.wall_M(x):.4f})' for x in xs)}}};"]
    for x, anc in zip(WK.WALL_X, ("north west", "north west", "north")):
        m = WK.wall_M(x)
        out.append(rf"\fill[cBar] ({x},{m:.4f}) circle (1.4pt); \node[cBar,anchor={anc},inner sep=2pt] at ({x},{m:.4f}) {{{_s(m, 4)}}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_modes():
    tbl = WK.mode_table()
    l1, l2 = WK.sw_lambdas()
    col = lambda n: "$\\to0$" if n == 0 else str(n)
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$n$ & " + " & ".join(col(t[0]) for t in tbl) + r"\\", r"\midrule",
                      rf"橫向 $\lambda_1={_fr(l1)}$ 的個數 & " + " & ".join(f"${t[1]}$" for t in tbl) + r"\\",
                      rf"縱向 $\lambda_2={_fr(l2)}$ 的個數 & " + " & ".join(str(t[2]) for t in tbl) + r"\\",
                      r"$\det X=\lambda_1^{n-1}\lambda_2$ & " + " & ".join(f"${_fr(t[3])}$" for t in tbl) + r"\\",
                      r"$(\det X)^{-1/2}$ & " + " & ".join(_s(t[4], 4) for t in tbl) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def fig_hole():
    """−D(Δ)：64³ 週期格子沿軸（r = −1，33 點，藍）與 Coulomb 擬合 c0 − c/Δ（Δ ≥ 2，橙虛線）。"""
    Dax = WK.D_axis()
    c0, cc = WK.coulomb_fit()
    out = [r"\begin{tikzpicture}[x=0.165cm,y=24cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (0,0) -- (35,0) node[right] {$\Delta$};",
           r"\draw[->] (0,-0.122) -- (0,0.014) node[above] {$-D$};",
           r"\foreach \x in {8,16,24,32} {\draw (\x,0) -- (\x,0.002) node[above] {$\x$};}",
           r"\node[above left,inner sep=1pt] at (0,0) {$0$};",
           r"\foreach \y/\l in {-0.05/-0.05,-0.1/-0.10} {\draw (0,\y) -- (-0.45,\y) node[left] {$\l$};}"]
    pts = [(d, -float(v)) for d, v in enumerate(Dax)]
    fit = [(d, c0 - cc / d) for d in range(2, 33)]
    out.append(r"\draw[cWalk,densely dashed,thick] plot coordinates {" + " ".join(f"({a},{b:.5f})" for a, b in fit) + "};")
    out.append(r"\draw[cBar,thick] plot coordinates {" + " ".join(f"({a},{b:.5f})" for a, b in pts) + "};")
    out += [rf"\fill[cBar] ({a},{b:.5f}) circle (1.1pt);" for a, b in pts if a in WK.HOLE_D]
    out.append(rf"\node[cBar,anchor=west] at (1.3,{-Dax[0]:.4f}) {{$\Delta=0$：$1-D={_s(1 - Dax[0], 4)}$（少 11\%）}};")
    out.append(r"\draw[cBar,thick] (11,-0.072) -- (14,-0.072) node[right,black] {$-D(\Delta)$，$64^3$ 週期格子};")
    out.append(rf"\draw[cWalk,densely dashed,thick] (11,-0.088) -- (14,-0.088) node[right,black] "
               rf"{{$c_0-c/\Delta$，$c={_s(cc, 4)}\approx\frac1{{4\pi}}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_hole():
    Dax = WK.D_axis()
    ps = WK.CONT_RHO
    rows = [(r"連續 $\frac{1-e^{-\sqrt2\rho}}{4\pi\rho}$", [WK.D_cont(p) for p in ps]),
            (r"無限大格子", [WK.D_inf((p, 0, 0)) for p in ps]),
            (r"$64^3$ 週期格子", [Dax[p] for p in ps]),
            (r"Coulomb $\frac1{4\pi\rho}$", [1 / (4 * math.pi * p) for p in ps])]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$D$（$r=-1$），$\rho=$ & " + " & ".join(map(str, ps)) + r"\\", r"\midrule",
                      *[lab + " & " + " & ".join(_s(v, 4) for v in vals) + r"\\" for lab, vals in rows],
                      r"\bottomrule", r"\end{tabular}"])


def tab_visits():
    Rs = WK.VIS_R
    vis = [WK.sc_visits(100, R) for R in Rs]
    ef = [WK.erfc_form(R, 100) for R in Rs]
    f4 = lambda v, R: _s(v, 6 if R == 16 else 4)
    book = ["", "", ""] + [rf"{{\color{{cWrong}}${_e(WK.erfc_form_book(R, 100))}$}}" for R in (8, 16)]
    return "\n".join([r"\begin{tabular}{@{}l*{5}{c}@{}}", r"\toprule",
                      r"$R$ & " + " & ".join(map(str, Rs)) + r"\\", r"\midrule",
                      r"精確 & " + " & ".join(f4(v, R) for v, R in zip(vis, Rs)) + r"\\",
                      r"erfc 形 & " + " & ".join(f4(v, R) for v, R in zip(ef, Rs)) + r"\\",
                      r"比 & " + " & ".join(_s(v / e, 3) for v, e in zip(vis, ef)) + r"\\",
                      r"書上 & " + " & ".join(book) + r"\\", r"\midrule",
                      r"$G(R)\cdot\frac{2\pi R}3$ & " + " & ".join(_s(WK.green_sc(R) * 2 * math.pi * R / 3, 4) for R in Rs) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def fig_rho():
    """環面的 ρ(z) = ⟨N⟩/V（實線）與 Var N/V（虛線），L = 3、4、5 由淺到深，各 60 點；z = z_c 的鉛直虛線。"""
    zc = WC.ZC_2D
    shade = {3: "!40", 4: "!70", 5: ""}
    out = [r"\begin{tikzpicture}[x=4.6cm,y=2.5cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (0,0) -- (1.29,0) node[right] {$z$};",
           r"\draw[->] (0,0) -- (0,1.08);",
           r"\foreach \x in {0,0.6,1.2} {\draw (\x,0) -- (\x,-0.025) node[below] {$\x$};}",
           r"\foreach \y in {0.5,1} {\draw (0,\y) -- (-0.012,\y) node[left] {$\y$};}",
           rf"\draw[gray!60,densely dashed] ({zc:.4f},0) -- ({zc:.4f},0.84) node[above,black] {{$z_c$}};"]
    for L in WK.TORI:
        zs, rho = WK.torus_curve(L)
        var = [WK.torus_moments(L, float(z))[2] / (L * L) for z in zs]
        assert len(zs) == 60
        out.append(rf"\draw[cBar{shade[L]},thick] plot coordinates {{{' '.join(f'({a:.4f},{b:.4f})' for a, b in zip(zs, rho))}}};")
        out.append(rf"\draw[cWalk{shade[L]},densely dashed,thick] plot coordinates {{{' '.join(f'({a:.4f},{b:.4f})' for a, b in zip(zs, var))}}};")
        zp, hp = WK.torus_peak(L)
        out.append(rf"\fill[cWalk{shade[L]}] ({zp:.4f},{hp:.4f}) circle (1.3pt);")
    out.append(r"\node[anchor=north west,align=left,inner sep=1pt] at (0.43,1.07) "
               r"{實線 {\color{cBar}$\rho=\langle N\rangle/V$}，虛線 {\color{cWalk}$\mathrm{Var}\,N/V$}\\由淺到深：$L=3,4,5$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_rho():
    rows = []
    for L in WK.TORI:
        zp, hp = WK.torus_peak(L)
        rows.append(rf"${L}\times{L}$ & " + " & ".join(_s(WK.torus_rho(L, z), 4) for z in WK.RHO_Z)
                    + rf" & {_s(zp, 4)} & {_s(hp, 4)}\\")
    return "\n".join([r"\begin{tabular}{@{}l*{6}{c}@{}}", r"\toprule",
                      r"$\rho$，$z=$ & 0.2 & $z_c$ & 0.6 & 1.0 & 峰 $z^*$ & 峰高\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_exp():
    E = {k: WK.rr_exponents(*WK.RR_CASES[k]) for k in ("mf", "d2", "d3")}
    laws = {k: WK.rr_laws(E[k], d) for k, (_, d) in WK.RR_CASES.items()}
    cell = lambda v, k: _s(v, 3) if k == "d3" else f"${_fr(v)}$"
    lab = {"beta": (r"$\beta'=d\nu-1$", r"$1-\alpha$"), "gamma": (r"$\gamma'=2-d\nu$", r"$\alpha$"),
           "alpha": (r"$\alpha'=2-d\nu$", r"$\alpha$"), "nu": (r"$\nu'=\nu$", r"$\nu$")}
    rows = [f"{lab[q][0]} & {lab[q][1]} & " + " & ".join(cell(E[k][q], k) for k in E) + r"\\" for q in MODES]
    ck = r"\,$\checkmark$"
    rows += [r"\midrule",
             r"$\alpha'+2\beta'+\gamma'$ & 2 & " + " & ".join("2" + ck for k in E if _close(laws[k][0], 2)) + r"\\",
             r"$2-\alpha'$ & & " + " & ".join(cell(laws[k][1], k) for k in E) + r"\\",
             r"$d\,\nu'$ & & " + " & ".join(cell(laws[k][2], k) + ck for k in E) + r"\\",
             r"書上 $d_c\nu'$ & & " + " & ".join((cell(laws[k][3], k) if k == "mf" else
                                                    rf"{{\color{{cWrong}}{cell(laws[k][3], k)}}}") for k in E) + r"\\",
             rf"$d_c=(2-\alpha')/\nu'$ & 4 & {WK.dc_mean_field()} & & \\"]
    head = (r"指數 & Table 10.1 & \shortstack{平均場\\$d=4$} & \shortstack{二維\\$\nu=\frac34$} & "
            r"\shortstack{三維（引用）\\$\nu\approx0.588$}\\")
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule", head, r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    (r"p.~235，(10.1)", r"$\sum_NC(N;\vec x,\vec y)\,z^n$", r"$\sum_NC(N;\vec x,\vec y)\,z^N$"),
    (r"p.~238，(10.14)$\to$(10.15)", r"對 $S$ 的積分是 $i\pi$；$\Gamma_{MF}=\frac{N_0}4\int_0^\infty f\,d\lambda$",
     r"$2\pi i\,\theta(\lambda)$；$\Gamma_{MF}=\frac{N_0}2\int_0^\infty f\,d\lambda$"),
    (r"p.~238，註 3", r"$\lambda\to4\lambda$、$4u\to f$、$2r\to1-z/z_c$", r"$\lambda\to2\lambda$、$u\to f$、$r\to1-z/z_c$"),
    (r"p.~242，(10.32)", r"$\exp\big[-\frac14\vec h\cdot X^{-1}\cdot\vec h\big]$", r"$\exp\big[+\frac14\vec h\cdot X^{-1}\cdot\vec h\big]$"),
    (r"p.~242，(10.33)", r"$\sqrt{\lambda_2/\lambda_1}=\sqrt{\frac{r+q^2+u|\vec M|^2}{r+q^2+3u|\vec M|^2}}$",
     r"$\sqrt{\lambda_1/\lambda_2}=\cdots$（右邊對）"),
    (r"p.~245，(10.45)", r"$e^{-N_0(r\lambda+u\lambda^2)}$", r"$e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}$"),
    (r"p.~245，正文", r"$D$ \emph{increases} with increasing separation", r"$-D$（與 $1-D$）隨距離增加；$D$ 遞減"),
    (r"p.~245，(10.46)", r"$\sum_{\vec q}\frac{\psi_{\vec q}(\vec x)\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}$",
     r"$\sum_{\vec q\ne0}\frac{2|r|\,\psi_{\vec q}(\vec x)\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}$"),
    (r"pp.~247--248，(10.53)(10.55)", r"$e^{-\frac V{a^d}\frac{r^2}{4u}}$；$\ln V-\frac V{a^d}\frac{r^2}{4u}$",
     r"$e^{+\frac V{a^d}\frac{r^2}{4u}}$；$\ln V+\frac V{a^d}\frac{r^2}{4u}$"),
    (r"p.~251，(10.75)", r"$2-\alpha'=d_c\nu'$（hold in general）", r"$2-\alpha'=d\,\nu'$"),
    (r"p.~253，(10.83)", r"$\mathrm{erfc}\big(\frac{z_cR^2}{2(N+1)}\big)$", r"$\mathrm{erfc}\big(\sqrt{\frac{z_cR^2}{4(N+1)}}\big)$"),
]


def tab_errata():
    assert len(ERRATA) == 11
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("TenFigTorus", fig_torus),
    ("TenTabTorus", tab_torus),
    ("TenFigMF", fig_mf),
    ("TenTabN", tab_n),
    ("TenTabAsym", tab_asym),
    ("TenFigWall", fig_wall),
    ("TenTabModes", tab_modes),
    ("TenFigHole", fig_hole),
    ("TenTabHole", tab_hole),
    ("TenTabVisits", tab_visits),
    ("TenFigRho", fig_rho),
    ("TenTabRho", tab_rho),
    ("TenTabExp", tab_exp),
    ("TenTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch10()
    parts = ["% figs.tex —— 由 ch10/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch10/make_figs.py（第 10 章，動畫 W1001–W1005）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    check_tex_decimals()                                   # 圖表裡用到的小數也已記進 OK
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
