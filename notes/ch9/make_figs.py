r"""
make_figs.py —— 隨機漫步導讀第 9 章（Rudnick & Gaspari《Elements of the Random Walk》第 9 章，動畫 W901–W906）的圖與表

在 notes/ 或 notes/ch9/ 裡執行：
    python ch9/make_figs.py
會先檢查 ch9/part9.tex 寫到的每一個數字（check_numbers_ch9；最後再把 part9.tex 裡出現的每一個小數和檢查過的數比對一遍），
再寫出 ch9/figs.tex（part9.tex 第一行 \input{ch9/figs}）。圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算）。
數學一律由 codes/utils/walk_crit.py 算出（二維自迴避的局部指數用 walk_saw）；例子與動畫 W901–W906 相同。巨集名稱都以 \Nine 開頭
（不和前幾章的 \Fig…、\Two…、…、\Eight… 撞名）；part9.tex 的 \label 用 nine-…。顏色沿用導言區：cBar 精確值、cWalk 路線、
cExact 更正、cGauss 近似、cWrong 筆誤。Metropolis（flip_tau、mc_trace）由 walk_crit 的 lru_cache 快取，同一次執行只算一次（全部約 4 秒）。
梯子圖另外用鍵集合的暴力列舉（_ladder_sets）核對 walk_crit.ladder_G 的閉式，並由它決定要畫哪些路徑與閉圈。
"""
from __future__ import annotations

import math
import re
import sys
import warnings
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import special

HERE = Path(__file__).resolve().parent            # notes/ch9
sys.path.insert(0, str(HERE.parent.parent / "codes"))
warnings.filterwarnings("ignore")                  # scipy 的 IntegrationWarning 不影響結果

from utils import walk_saw as SW  # noqa: E402
from utils import walk_crit as WC  # noqa: E402   （套件模組：walk_crit 的相對 import 才找得到 walk_saw、walk_scale）

OUT = HERE / "figs.tex"
TEX = HERE / "part9.tex"
F = Fraction
OK = set()                                         # 檢查過的小數字串（part9.tex 的小數都要在這裡）


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣），並記進 OK。"""
    s = f"{float(x):.{d}f}"
    OK.add(s.lstrip("-"))
    return s


def _g(x, n=4):
    """有效數字 n 位（walk_crit 的 .4g 寫法）；科學記號寫成 a\times10^{b}。尾數記進 OK。"""
    s = f"{float(x):.{n}g}"
    if "e" in s:
        mant, ex = s.split("e")
        OK.add(mant.lstrip("-"))
        return rf"{mant}\times10^{{{int(ex)}}}"
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
# 梯子（2×3，6 點 7 鍵）上的鍵集合列舉：Π(1 + K S_i·S_j) 展開後存活的項
# =============================================================================
LAD_XY = {i: (i % 3, i // 3) for i in range(6)}    # 格點 (x, y) → x + 3y


def _components(sub):
    par = {}

    def root(v):
        while par.setdefault(v, v) != v:
            v = par[v]
        return v

    for i, j in sub:
        par[root(i)] = root(j)
    return len({root(v) for e in sub for v in e})


def _ladder_sets(a=None, b=None):
    """存活的鍵集合 [(鍵, 閉圈數)]：給 a、b 時端點度數 1、其他點 0 或 2（分子）；否則每點 0 或 2（分母 Z）。
    （這張圖每點最多 3 條鍵，所以「偶數次」就是 0 或 2。）"""
    Ns, edges = WC.LADDER
    out = []
    for k in range(len(edges) + 1):
        for sub in combinations(edges, k):
            deg = [0] * Ns
            for i, j in sub:
                deg[i] += 1
                deg[j] += 1
            if all((d == 1) if v in (a, b) else d in (0, 2) for v, d in enumerate(deg)):
                loops = _components(sub) - (0 if a is None else 1)
                out.append((sub, loops))
    return sorted(out, key=lambda t: (t[1], len(t[0])))


def _ladder_G_enum(n, K, pair):
    num = sum(K ** len(s) * n ** c for s, c in _ladder_sets(*pair))
    den = sum(K ** len(s) * n ** c for s, c in _ladder_sets())
    return num / den


def _path_part(sub, a):
    """sub 中與端點 a 相連的鍵（路徑）與其餘的鍵（閉圈）。"""
    if a is None:
        return [], list(sub)
    reach, path, todo = {a}, [], list(sub)
    changed = True
    while changed:
        changed = False
        for e in list(todo):
            if e[0] in reach or e[1] in reach:
                reach.update(e)
                path.append(e)
                todo.remove(e)
                changed = True
    return path, todo


# =============================================================================
# 數字檢查（part9.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch9():
    OK.update({"0.5", "0.9", "0.99", "0.1", "0.01", "0.2", "0.15", "0.05", "0.025", "2.0", "3.0", "0.0001", "0.3"})
    # ---- W901 漫步的指數、自發磁化
    assert [WC.flory_p(d) for d in range(1, 5)] == [2, F(3, 2), F(6, 5), 1]
    assert [2 * WC.NU_KNOWN[d] for d in (1, 2, 4)] == [2, F(3, 2), 1]
    assert (_s(WC.NU_KNOWN[3], 3), _s(2 * WC.NU_KNOWN[3], 3), _s(WC.flory_p(3) / 2, 3)) == ("0.588", "1.176", "0.600")
    assert _s(WC.flory_gap_3d(), 4) == "0.0204"
    _all([SW.two_nu_eff(N) for N in WC.SAW_N], 4, ("1.4442", "1.4505", "1.4556"))
    ms = WC.ms_table()
    assert tuple(ms) == WC.MS_X == (0.5, 0.9, 0.99)
    _all([ms[x][0] for x in WC.MS_X], 4, ("0.9575", "0.5254", "0.1725"))
    _all([ms[x][1] for x in WC.MS_X], 4, ("0.9980", "0.8959", "0.6859"))
    assert all(_close(WC.mf_m(x), WC.mf_m_newton(x), 1e-9) for x in WC.MS_X)         # 兩種算法
    assert _s(WC.mf_m_asym(0.99), 4) == "0.1732" and _close(WC.mf_m_asym(0.99), math.sqrt(0.03))
    assert _s(WC.TC_2D, 4) == "2.2692" and _close(math.sinh(2 / WC.TC_2D), 1, 1e-12)
    assert _close(WC.TC_2D, 2 / math.log(1 + math.sqrt(2)), 1e-12)
    assert (_s(WC.beta_slope(WC.mf_m, 1.0), 4), _s(WC.beta_slope(WC.onsager_M, WC.TC_2D), 4)) == ("0.5000", "0.1250")
    assert [_s(WC.BETA_CITED[k], 3) for k in ("3d_ising", "3d_saw", "book_3d_saw")] == ["0.326", "0.302", "0.308"]
    # ---- W902 O(n)、對稱破缺
    E, E2, M, M2 = WC.rotation_check()
    assert _s(E, 4) == _s(E2, 4) == "-6.3851"
    assert _s(np.linalg.norm(M), 4) == _s(np.linalg.norm(M2), 4) == "6.7943"
    assert (_s(WC.angle_deg(M), 1), _s(WC.angle_deg(M2), 1)) == ("-13.1", "106.9")
    assert _close(WC.angle_deg(M2) - WC.angle_deg(M), 120, 1e-9)
    for K, Z, P, am in ((1.0, "121.2329", ("0.4504", "0.0330", "0.0333"), "0.9337"),
                        (0.2, "17.3497", ("0.1283", "0.2306", "0.2823"), "0.4871")):
        ex, cl = WC.ising_exact(WC.RING4, K), WC.ring4_closed(K)
        assert _s(ex["Z"], 4) == _s(cl[0], 4) == Z and _s(ex["absm"], 4) == _s(cl[1], 4) == am
        assert _close(cl[0], 2 * math.exp(4 * K) + 12 + 2 * math.exp(-4 * K))
        _all([ex["PM"][4], ex["PM"][2], ex["PM"][0]], 4, P)
        _all(cl[2:], 4, P)
        assert all(_close(ex["PM"][m], ex["PM"][-m]) for m in (2, 4)) and abs(ex["m"]) < 1e-12
        assert _close(cl[2], math.exp(4 * K) / cl[0]) and _close(cl[3], 4 / cl[0])
        assert _close(cl[4], (4 + 2 * math.exp(-4 * K)) / cl[0])
    pz, pa = WC.pair_closed(1.0)
    ex = WC.ising_exact(WC.PAIR, 1.0)
    assert _s(pz, 4) == _s(ex["Z"], 4) == "6.1723" and _s(pa, 4) == _s(ex["absm"], 4) == "0.8808"
    assert _close(pz, 2 * (math.e + 1 / math.e)) and _close(pa, math.e / (math.e + 1 / math.e))
    assert (WC.MC_SEED, WC.MC_SWEEPS, WC.MC_REPLICAS, WC.MC_M0, WC.MC_T, WC.MC_L) == (902, 4000, 50, 0.5, 2.0, (4, 6, 8))
    assert WC.MC_TRACES == ((8, 3.0), (4, 2.0), (8, 2.0))
    assert tuple(f"{WC.flip_tau(L):.0f}" for L in WC.MC_L) == ("133", "608", "2151")
    assert [WC.trace_flips(L, T) for L, T in WC.MC_TRACES[1:]] == [30, 2]
    assert _s(WC.onsager_M(2.0), 4) == "0.9113"
    (a1, a2), (e1, e2) = WC.mc_check_L4()
    assert (_s(a1, 4), _s(a2, 4), _s(e1, 4), _s(e2, 4)) == ("0.9182", "0.9189", "-1.7528", "-1.7554")
    assert WC.FIELD_H == (0.1, 0.01)
    _all([WC.ising_exact(WC.RING4, 1.0, H)["m"] for H in WC.FIELD_H], 4, ("0.3506", "0.0367"))
    _all([WC.ising_exact(WC.TORUS4, 1.0, H)["m"] for H in WC.FIELD_H], 4, ("0.9210", "0.1584"))
    assert _s(WC.onsager_M(1.0), 4) == "0.9993"
    # ---- W903 臨界指數、標度律
    b = WC.BOOK_3D
    assert (_s(b["alpha"], 2), _s(b["beta"], 3), _s(b["gamma"], 2), _s(b["delta"], 3), _s(b["nu"], 3), _s(b["eta"], 3)) == \
        ("0.12", "0.308", "1.19", "4.865", "0.602", "0.023") and b["d"] == 3
    a16, a18, lhs = WC.book3d_alpha()
    assert (_s(a16, 3), _s(a18, 3), _s(lhs, 3), _s(a16, 2)) == ("0.194", "0.194", "1.926", "0.19")
    lb = WC.laws(b)
    assert (_s(lb["fisher"], 3), _s(lb["hyper_delta"], 3), _s(lb["widom"], 3)) == ("1.190", "4.865", "1.190")
    assert _s(b["nu"] * (b["d"] - 2 + b["eta"]) / 2, 3) == "0.308"
    c = WC.CITED_3D
    assert (_s(c["nu"], 4), _s(c["gamma"], 3), _s(c["eta"], 3), _s(c["alpha"], 3), _s(c["beta"], 3)) == \
        ("0.5876", "1.157", "0.031", "0.237", "0.302")
    K = WC.K1D
    assert K == 1.0 and _s(math.tanh(K), 4) == "0.7616"
    assert WC.ising1d_G(0, K) == 1
    _all([WC.ising1d_G(r, K) for r in range(1, 5)], 4, ("0.7616", "0.5800", "0.4417", "0.3364"))
    assert all(_close(WC.ising1d_G(r, K), WC.ising1d_transfer_G(r, K), 1e-12) for r in range(9))
    assert all(_close(WC.ising1d_G(r, K), math.exp(-r / WC.ising1d_xi(K)), 1e-12) for r in range(9))
    assert _s(WC.ising1d_xi(K), 4) == "3.6719" and _s(WC.ising1d_chi(K), 4) == "7.3891"
    assert _close(WC.ising1d_chi(K), math.exp(2)) and _close(WC.ising1d_chi(K), 1 + 2 * math.tanh(K) / (1 - math.tanh(K)))
    _all([WC.ising1d_chi_partial(K, R) for R in (5, 10, 20)], 4, ("5.7520", "6.9696", "7.3615"))
    assert _s(WC.ising1d_ring_closed(16, K)[1], 4) == "7.2021"
    assert WC.FISHER1D_K == (1.0, 2.0, 3.0)
    _all([WC.ising1d_xi2_moment(k) / WC.ising1d_xi(k) ** 2 for k in WC.FISHER1D_K], 4, ("1.9877", "1.9998", "2.0000"))
    assert all(_close(WC.ising1d_xi2_moment(k), 2 * math.tanh(k) / (1 - math.tanh(k)) ** 2) for k in WC.FISHER1D_K)
    assert WC.ising1d_xi2_moment(3.0) / WC.ising1d_xi(3.0) > WC.ising1d_xi2_moment(1.0) / WC.ising1d_xi(1.0) * 5
    _all([WC.ising1d_chi(k) / WC.ising1d_xi(k) for k in WC.FISHER1D_K], 4, ("2.0123", "2.0002", "2.0000"))
    E3 = WC.EXPONENTS
    expect = {"mf": (0, F(1, 2), 1, 3, F(1, 2), 0, 4), "ising2d": (0, F(1, 8), F(7, 4), 15, 1, F(1, 4), 2),
              "saw2d": (F(1, 2), F(5, 64), F(43, 32), F(91, 5), F(3, 4), F(5, 24), 2)}
    for k, v in expect.items():
        e = E3[k]
        assert tuple(e[x] for x in ("alpha", "beta", "gamma", "delta", "nu", "eta", "d")) == v
        L = WC.laws(e)
        assert L["rushbrooke"] == 2 and L["fisher"] == e["gamma"] and L["hyper_alpha"] == e["alpha"]
        assert L["hyper_delta"] == e["delta"] and L["widom"] == e["gamma"]
        assert WC.from_nu_eta(e["nu"], e["eta"], e["d"]) == (e["alpha"], e["beta"], e["gamma"], e["delta"])
        D = e["beta"] + e["gamma"]                          # Δ = β + γ（練習 9.1–9.2）
        assert WC.widom_exponents(e["alpha"], D) == (e["beta"], e["gamma"], e["delta"]) and D / e["beta"] == e["delta"]
    assert 2 - 3 * F(1, 2) == F(1, 2) and E3["saw2d"]["beta"] + E3["saw2d"]["gamma"] == F(91, 64)
    assert WC.WIDOM_CASE == (0, F(15, 8)) and WC.widom_exponents(*WC.WIDOM_CASE) == (F(1, 8), F(7, 4), 15)
    for meth in ("exact", "fd"):
        _all(WC.widom_slopes(*WC.WIDOM_CASE, method=meth), 4, ("0.1250", "1.7500", "0.0667"))
    assert F(1, 8) * 14 == F(7, 4)
    # ---- W904 GLW、雙井、練習 9.3–9.4
    assert max(WC.fourier_check()) < 1e-12
    assert WC.LANDAU_R == (1.0, 0.0, -1.0) and WC.LANDAU_H == 0.1
    assert _close(WC.landau_min(1.0), 0) and _close(WC.landau_min(-1.0), 1)
    m1, f1 = WC.landau_closed(-1.0)
    assert _close(m1, 1) and _s(f1, 2) == "-0.25" and _close(WC.landau_F(1.0, -1.0), -0.25)
    assert _s(WC.landau_min(-1.0, 1.0, WC.LANDAU_H), 4) == "1.0467"
    m3, f3 = WC.ex94_closed(-1.0)
    assert (_s(m3, 4), _s(f3, 2)) == ("1.7321", "-0.75") and _close(m3, math.sqrt(3))
    assert _close(WC.ex94_m(-1.0, 0.0), m3, 1e-9)
    assert _s(WC.ex94_scaling(1.0, -0.01), 4) == _s(WC.ex94_scaling(1.0, -0.0001), 4) == "2.1038"
    Mx = WC.ex94_scaling(1.0, -0.01)
    assert _close(Mx ** 3 / 3 - Mx, 1.0, 1e-9)
    _all(WC.ex94_slopes(), 4, ("0.5000", "1.5000", "0.0000", "0.3333", "1.0000"))
    hh = 1e-7                                               # γ = γ' = 1，振幅比 2
    assert abs(WC.ex94_m(0.01, hh) / hh / 100 - 1) < 1e-3
    assert abs((WC.ex94_m(-0.01, hh) - WC.ex94_closed(-0.01)[0]) / hh / 50 - 1) < 1e-3
    assert _close(WC.ex94_m(0.0, 0.001), 0.003 ** (1 / 3), 1e-9)
    assert WC.ZMF_V == (10, 100, 1000) and WC.ZMF_H == 0.01
    for meth in ("quad", "grid"):
        _all([WC.zmf_avg(V, h=WC.ZMF_H, method=meth) for V in WC.ZMF_V], 4, ("0.0869", "0.7576", "1.0042"))
    _all([WC.zmf_avg(V, what="absm") for V in WC.ZMF_V], 4, ("0.8888", "0.9922", "0.9992"))
    _all([WC.zmf_avg(V, what="f") for V in WC.ZMF_V], 4, ("-0.2719", "-0.2397", "-0.2478"))
    assert all(abs(WC.zmf_avg(V)) < 1e-9 for V in WC.ZMF_V)
    assert _s(WC.landau_min(-1.0, 1.0, WC.ZMF_H), 4) == "1.0050"
    # ---- W905 由磁鐵讀出 ⟨R²⟩
    assert _s(WC.ZC_2D, 4) == "0.3791" and _s(1 / WC.ZC_2D, 3) == "2.638"
    assert WC.saw_table() == [(1, 4, 4), (2, 12, 32), (3, 36, 164), (4, 100, 704), (5, 284, 2716), (6, 780, 9808)]
    assert 4 * 4 + 8 * 2 == 32
    chi1, num1, x1 = WC.saw_gf(0.1)
    assert (_s(chi1, 4), _s(num1, 4), _s(x1, 4)) == ("1.5699", "0.9964", "0.6347") and _close(num1 / chi1, x1)
    assert WC.SAW_Z == (0.1, 0.15, 0.2)
    _all([WC.saw_gf(z)[2] for z in WC.SAW_Z], 4, ("0.6347", "1.2799", "2.4593"))
    assert WC.S_GAMMA == F(43, 32) and WC.S_GAMMA2NU == F(91, 32) and WC.COEF_N == (10, 100, 1000)
    _all([WC.coef_ratio(N, WC.S_GAMMA) for N in WC.COEF_N], 4, ("1.0227", "1.0023", "1.0002"))
    _all([WC.coef_ratio(N, WC.S_GAMMA2NU) for N in WC.COEF_N], 4, ("1.2760", "1.0264", "1.0026"))
    assert [_g(WC.coef_ratio_book(N, WC.S_GAMMA, F(91, 5)), 2) for N in (10, 100)] == [r"1.4\times10^{-17}", r"1.9\times10^{-34}"]
    _all([WC.r2_ratio(N, WC.S_GAMMA, F(3, 4), WC.ZC_2D) for N in WC.COEF_N], 4, ("1.2476", "1.0240", "1.0024"))
    zp = WC.zc_pow()
    assert (_s(zp[0], 4), _s(zp[1], 4)) == ("4.2850", "3.6823")
    assert _close(zp[0], WC.ZC_2D ** -1.5) and _close(zp[1], WC.ZC_2D ** (-43 / 32))
    assert WC.SAW_N == (10, 12, 14)
    _all([WC.saw_amp(N) for N in WC.SAW_N], 4, ("0.8299", "0.8224", "0.8168"))
    # ---- W906 n → 0、路徑積分、高斯積分
    _all([WC.hs_rhs(s) for s in (1, 2)], 4, ("0.7788", "0.3679"))
    assert all(_close(WC.hs_rhs(s), WC.hs_lhs(s), 1e-9) for s in (1, 2))
    assert _close(WC.hs_lhs(1), math.exp(-0.25)) and _close(WC.hs_lhs(2), math.exp(-1))
    z72 = WC.hs_book72()
    assert (_s(z72.real, 4), _s(z72.imag, 4)) == ("0.3894", "-0.2395")
    assert abs(z72 - 0.5 * math.exp(-0.25) * (1 - 1j * special.erfi(0.5))) < 1e-9
    assert _s(WC.intersection_factor(), 4) == "0.6065" and _close(WC.intersection_factor(), math.exp(-0.5))
    assert all(_close(WC.visit_weight(k), math.exp(-k * k / 4)) for k in (1, 2, 3))
    assert all(WC.pair_sum(k) == k * k == k + 2 * math.comb(k, 2) for k in range(1, 6))
    w = [WC.weak_saw_count(8, u) for u in WC.WEAK_U]
    assert WC.WEAK_N == 8 and WC.WEAK_U == (0.0, 1.0, 4.0, math.inf)
    assert w[0] == 4 ** 8 == 65536 and (_s(w[1], 2), _s(w[2], 2)) == ("23695.11", "7988.15") and w[3] == 5916
    assert _g(WC.weak_saw_count(8, 1.0, sign=+1), 2) == r"6.2\times10^{5}"
    assert (_s(WC.kernel_mass_book(0.1, 1), 4), _s(WC.kernel_mass_book(0.1, 3), 5)) == ("0.3162", "0.00609")
    I = WC.gauss_I()
    assert _s(I, 4) == _s(WC.gauss_I_quad(), 4) == "5.0627" and _close(I, 2 * math.pi / math.sqrt(3) * math.exp(1 / 3))
    for (p, q), v in (((1, 1), F(1, 3)), ((2, 0), F(2, 3)), ((2, 2), F(2, 3)), ((4, 0), F(4, 3))):
        assert _close(WC.gauss_moment(WC.GAUSS_T, p, q), v, 1e-9)
    assert F(4, 9) + F(2, 9) == F(2, 3) and 3 * F(2, 3) ** 2 == F(4, 3)
    C2 = np.linalg.inv(np.array([[2.0, -1.0], [-1.0, 3.0]]))
    assert (_s(WC.wick4(C2, 0, 0, 1, 1), 2), _s(WC.wick4_book(C2, 0, 0, 1, 1), 2)) == ("0.32", "0.44")
    assert _close(WC.gauss_moment(((2.0, -1.0), (-1.0, 3.0)), 2, 2), 0.32, 1e-9)
    assert [WC.wick_count(m) for m in range(1, 5)] == [1, 3, 15, 105] == [WC.wick_count_enum(m) for m in range(1, 5)]
    assert WC.TROTTER_EPS == (0.1, 0.05, 0.025) and WC.TROTTER_M == (10, 20, 40)
    t1 = [WC.trotter_err(e, 1) for e in WC.TROTTER_EPS]
    t2 = [WC.trotter_err(e, 2) for e in WC.TROTTER_EPS]
    assert [_g(x) for x in t1] == ["0.005784", "0.001502", "0.0003827"]
    assert [_g(x) for x in t2] == ["0.0009465", "0.000124", r"1.589\times10^{-5}"]
    _all([t1[0] / t1[1], t1[1] / t1[2]], 2, ("3.85", "3.92"))
    _all([t2[0] / t2[1], t2[1] / t2[2]], 2, ("7.63", "7.81"))
    _all([t1[k] / WC.trotter_lead(e) for k, e in enumerate(WC.TROTTER_EPS)], 4, ("0.9268", "0.9627", "0.9811"))
    assert [_g(WC.trotter_global(M)) for M in WC.TROTTER_M] == ["0.02366", "0.0117", "0.005821"]
    assert max(WC.ex95_check()) < 1e-10
    assert WC.LADDER_K == 0.5 and WC.LADDER_PAIRS == ((0, 3), (0, 5))
    exp = {(0, 3): ("0.6027", "0.5610", "0.65625"), (0, 5): ("0.3562", "0.3171", "0.40625")}
    for p in WC.LADDER_PAIRS:
        g1, g2, g0 = (WC.ladder_G(n, 0.5, p) for n in (1, 2, 0))
        assert (_s(g1, 4), _s(g2, 4), _s(g0, 5)) == exp[p]
        assert _close(g1, WC.ladder_G_ising(0.5, p)) and _close(g2, WC.ladder_G_xy(0.5, p))
        assert _close(g0, WC.ladder_saw_sum(0.5, p))
        for n in (0, 1, 2, 3):
            for k in (0.3, 0.5):
                assert _close(_ladder_G_enum(n, k, p), WC.ladder_G(n, k, p), 1e-12)
        assert sorted(len(s) for s, c in _ladder_sets(*p) if c == 0) == WC.ladder_saws(*p)
    assert [(len(s), c) for s, c in _ladder_sets(0, 3)] == [(1, 0), (3, 0), (5, 0), (5, 1)]
    assert [(len(s), c) for s, c in _ladder_sets(0, 5)] == [(3, 0), (3, 0), (3, 0), (5, 0)]
    assert [(len(s), c) for s, c in _ladder_sets()] == [(0, 0), (4, 1), (4, 1), (6, 1)]
    assert _close(0.5 + 0.5 ** 3 + 0.5 ** 5, 0.65625) and _close(3 * 0.5 ** 3 + 0.5 ** 5, 0.40625)


def check_tex_decimals():
    """part9.tex（去掉註解）裡的每一個小數都必須是上面檢查過的數。"""
    if not TEX.exists():
        print("（尚無 part9.tex，略過小數比對）")
        return
    txt = re.sub(r"(?<!\\)%.*", "", TEX.read_text(encoding="utf-8"))
    txt = re.sub(r"[\d.]+\\linewidth|[\d.]+pt\b|[\d.]+em\b|[\d.]+cm\b|\\arraystretch\}\{[\d.]+\}", "", txt)   # 版面尺寸
    eqno = re.compile(r"9\.\d{1,2}")                     # 本章的式號、節號、圖號、練習號（9.1–9.77）
    bad = sorted({m for m in re.findall(r"\d+\.\d+", txt) if m not in OK and not eqno.fullmatch(m)})
    assert not bad, f"part9.tex 有未檢查的小數：{bad}"


# =============================================================================
# 圖與表
# =============================================================================
def fig_ms():
    """M_s 對 T/T_c：平均場 m = tanh(m/T)（橙）與二維 Ising 的 Onsager 精確式（藍），各 ≤ 60 點；圓點是表上的三個溫度。"""
    xs = np.concatenate([np.linspace(0.05, 0.9, 20), 1 - np.geomspace(0.09, 1e-7, 30)])
    out = [r"\begin{tikzpicture}[x=4.2cm,y=2.6cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (0,0) -- (1.3,0) node[right] {$T/T_c$};",
           r"\draw[->] (0,0) -- (0,1.14) node[above] {$M_s$};",
           r"\node[below] at (0,-0.03) {$0$};",
           r"\foreach \x in {0.5,1} {\draw (\x,0) -- (\x,-0.03) node[below] {$\x$};}",
           r"\foreach \y in {0.5,1} {\draw (0,\y) -- (-0.02,\y) node[left] {$\y$};}",
           r"\draw[gray!50,densely dotted] (1,0) -- (1,1.05);"]
    for fn, col in ((WC.mf_m, "cGauss"), (lambda x: WC.onsager_M(x * WC.TC_2D), "cBar")):
        pts = [(0.0, 1.0)] + [(x, fn(x)) for x in xs] + [(1.0, 0.0), (1.2, 0.0)]
        assert len(pts) <= 60
        out.append(rf"\draw[{col},thick] plot coordinates {{{' '.join(f'({a:.4f},{b:.4f})' for a, b in pts)}}};")
        out += [rf"\fill[{col}] ({x},{fn(x):.4f}) circle (1.4pt);" for x in WC.MS_X]
    out.append(r"\draw[cBar,thick] (0.05,0.34) -- (0.13,0.34) node[right] {二維 Ising，$\beta=\frac18$};")
    out.append(r"\draw[cGauss,thick] (0.05,0.2) -- (0.13,0.2) node[right] {平均場，$\beta=\frac12$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_flory():
    ds = range(1, 5)
    known = [_s(2 * WC.NU_KNOWN[d], 3) if d == 3 else f"${_fr(2 * WC.NU_KNOWN[d])}$" for d in ds]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$d$ & " + " & ".join(map(str, ds)) + r"\\", r"\midrule",
                      r"Flory $p=\frac6{d+2}$ & " + " & ".join(f"${_fr(WC.flory_p(d))}$" for d in ds) + r"\\",
                      r"已知的 $2\nu$ & " + " & ".join(known) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_ms():
    ms = WC.ms_table()
    mf = " & ".join(_s(ms[x][0], 4) for x in WC.MS_X) + rf" & {_s(WC.beta_slope(WC.mf_m, 1.0), 4)}"
    on = " & ".join(_s(ms[x][1], 4) for x in WC.MS_X) + rf" & {_s(WC.beta_slope(WC.onsager_M, WC.TC_2D), 4)}"
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$T/T_c$ & " + " & ".join(map(str, WC.MS_X)) + r" & 斜率\\", r"\midrule",
                      rf"平均場 & {mf}\\", rf"Onsager & {on}\\", r"\bottomrule", r"\end{tabular}"])


def fig_pm():
    """4 個自旋一圈的 P(M)：K = 1（藍）與 K = 0.2（橙）並排的長條。"""
    out = [r"\begin{tikzpicture}[x=0.6cm,y=4.0cm,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw (-5,0) -- (5,0);",
           r"\node[anchor=south west] at (-5,0.5) {$P(M)$};"]
    for K, col, dx in ((1.0, "cBar", -0.3), (0.2, "cWalk", 0.3)):
        PM = WC.ising_exact(WC.RING4, K)["PM"]
        for M in (-4, -2, 0, 2, 4):
            out.append(rf"\fill[{col}!80] ({M + dx - 0.26:.2f},0) rectangle ({M + dx + 0.26:.2f},{PM[M]:.4f});")
    p4, p0 = WC.ising_exact(WC.RING4, 1.0)["PM"][4], WC.ising_exact(WC.RING4, 0.2)["PM"][0]
    out.append(rf"\node[cBar,above] at (3.7,{p4:.4f}) {{{_s(p4, 4)}}};")
    out.append(rf"\node[cWalk,above] at (0.3,{p0:.4f}) {{{_s(p0, 4)}}};")
    out.append(r"\foreach \m in {-4,-2,0,2,4} {\node[below] at (\m,0) {$\m$};}")
    out.append(r"\node[below] at (5,0) {$M$};")
    out.append(r"\fill[cBar!80] (-2.6,0.47) rectangle (-2.3,0.5); \node[anchor=west] at (-2.25,0.485) {$K=1$（$T=1$）};")
    out.append(r"\fill[cWalk!80] (-2.6,0.4) rectangle (-2.3,0.43); \node[anchor=west] at (-2.25,0.415) {$K=0.2$（$T=5$）};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_trace():
    """三條 Metropolis 軌跡 m(sweep)（各取 200 點）：L = 8、T = 3.0；L = 4、T = 2.0；L = 8、T = 2.0。"""
    lab = {(8, 3.0): r"$L=8$，$T=3.0$\\（高溫）",
           (4, 2.0): rf"$L=4$，$T=2.0$\\翻 {WC.trace_flips(4, 2.0)} 次",
           (8, 2.0): rf"$L=8$，$T=2.0$\\翻 {WC.trace_flips(8, 2.0)} 次"}
    out = [r"\begin{tikzpicture}[x=1.2cm,y=0.42cm,font=\scriptsize,baseline=(current bounding box.center)]"]
    for k, (L, T) in enumerate(WC.MC_TRACES):
        y0 = -2.5 * k
        sw, m = WC.mc_trace(L, T)
        sw, m = sw[::2], m[::2]
        assert len(sw) == 200
        out.append(rf"\draw[gray!50] (0,{y0 - 1:.1f}) rectangle (4,{y0 + 1:.1f});")
        out.append(rf"\draw[gray!50,densely dotted] (0,{y0:.1f}) -- (4,{y0:.1f});")
        pts = " ".join(f"({s / 1000:.3f},{y0 + v:.3f})" for s, v in zip(sw, m))
        out.append(rf"\draw[cBar,thin] plot coordinates {{{pts}}};")
        out.append(rf"\node[anchor=west,align=left,inner sep=2pt] at (4.05,{y0:.1f}) {{{lab[(L, T)]}}};")
    out.append(r"\node[left,inner sep=1pt] at (0,1) {$1$}; \node[left,inner sep=1pt] at (0,0) {$m$}; "
               r"\node[left,inner sep=1pt] at (0,-1) {$-1$};")
    yb = -2.5 * 2 - 1
    out.append(rf"\foreach \x/\l in {{0/0,2/2000,4/4000}} {{\node[below] at (\x,{yb:.1f}) {{\l}};}}")
    out.append(rf"\node[below] at (3.0,{yb:.1f}) {{sweep}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_flip():
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r"$L$（$T=2.0$） & " + " & ".join(map(str, WC.MC_L)) + r"\\", r"\midrule",
                      r"平均翻轉等待（sweep） & " + " & ".join(f"{WC.flip_tau(L):.0f}" for L in WC.MC_L) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_field():
    r4 = " & ".join(_s(WC.ising_exact(WC.RING4, 1.0, H)["m"], 4) for H in WC.FIELD_H)
    t4 = " & ".join(_s(WC.ising_exact(WC.TORUS4, 1.0, H)["m"], 4) for H in WC.FIELD_H)
    return "\n".join([r"\begin{tabular}{@{}lcc@{}}", r"\toprule",
                      r"$\langle m\rangle$（$K=1$） & $H=0.1$ & $H=0.01$\\", r"\midrule",
                      rf"4 個自旋一圈 & {r4}\\", rf"$4\times4$ 週期 & {t4}\\",
                      rf"\multicolumn{{3}}{{@{{}}l@{{}}}}{{無限大格子：$M_s={_s(WC.onsager_M(1.0), 4)}$}}\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_nineone():
    """書上 Table 9.1（照抄；兩處更正用紅／綠標出）。"""
    b = WC.BOOK_3D
    a_fix = rf"{{\color{{cWrong}}{_s(b['alpha'], 2)}}}$\to${{\color{{cExact}}$\approx${_s(WC.book3d_alpha()[0], 2)}}}"
    rows = [(r"比熱 $T\gtrless T_c$", r"$|T-T_c|^{-\alpha}$", r"$\alpha,\alpha'$", "0", a_fix),
            (r"磁化 $h=0$", r"$|T_c-T|^{\beta}$", r"$\beta$", r"$\frac12$", _s(b["beta"], 3)),
            (r"磁化率 $T>T_c$", r"$|T-T_c|^{-\gamma}$", r"$\gamma$", "1", _s(b["gamma"], 2)),
            (r"磁化率 $T<T_c$", r"$|T-T_c|^{-\gamma'}$", r"$\gamma'$", "1", r"無定義"),
            (r"磁化 $T=T_c$", r"$|h|^{1/\delta}$", r"$\delta$", "3", _s(b["delta"], 3)),
            (r"關聯長度（{\color{cWrong}$\xi$}$\to${\color{cExact}$\xi^2$}）", r"$|T-T_c|^{-\nu}$", r"$\nu,\nu'$", r"$\frac12$",
             _s(b["nu"], 3)),
            (r"$C(T_c,0;\rho)$", r"$\rho^{-d+2-\eta}$", r"$\eta$", "0", _s(b["eta"], 3))]
    body = [" & ".join(r) + r"\\" for r in rows]
    return "\n".join([r"\begin{tabular}{@{}llccc@{}}", r"\toprule", r"量 & 行為 & 指數 & 平均場 & 三維（書上）\\", r"\midrule",
                      *body, r"\bottomrule", r"\end{tabular}"])


def tab_laws():
    """三組精確指數（分數）與書上的三維欄；四條標度律＋Widom 逐格驗算（$\checkmark$／$\times$）。"""
    cols = [WC.EXPONENTS[k] for k in ("mf", "ising2d", "saw2d")]
    b = WC.BOOK_3D
    bfmt = {"alpha": 2, "beta": 3, "gamma": 2, "delta": 3, "nu": 3, "eta": 3}
    rows = []
    for key, sym in (("alpha", r"\alpha"), ("beta", r"\beta"), ("gamma", r"\gamma"), ("delta", r"\delta"),
                     ("nu", r"\nu"), ("eta", r"\eta")):
        rows.append(rf"${sym}$ & " + " & ".join(f"${_fr(e[key])}$" for e in cols) + rf" & {_s(b[key], bfmt[key])}\\")
    rows.append(r"$d$ & " + " & ".join(str(e["d"]) for e in cols) + rf" & {b['d']}\\")
    yes, no = r"{\color{cExact}$\checkmark$}", r"{\color{cWrong}$\times$}"
    law_rows = []
    lb = WC.laws(b)
    for key, lab, target in (("rushbrooke", r"$\alpha+2\beta+\gamma\overset?=2$", lambda e: 2),
                             ("fisher", r"$(2-\eta)\nu\overset?=\gamma$", lambda e: e["gamma"]),
                             ("hyper_alpha", r"$2-d\nu\overset?=\alpha$", lambda e: e["alpha"]),
                             ("hyper_delta", r"$\frac{d+2-\eta}{d-2+\eta}\overset?=\delta$", lambda e: e["delta"]),
                             ("widom", r"$\beta(\delta-1)\overset?=\gamma$", lambda e: e["gamma"])):
        cells = []
        for e in cols:
            v = WC.laws(e)[key]
            assert v == target(e)
            cells.append(f"${_fr(v)}$ {yes}")
        ok_b = abs(lb[key] - target(b)) < 0.002
        assert ok_b == (key not in ("rushbrooke", "hyper_alpha"))
        cells.append(f"{_s(lb[key], 3)} {yes if ok_b else no}")
        law_rows.append(lab + " & " + " & ".join(cells) + r"\\")
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r" & 平均場 & 二維 Ising & 二維 $n=0$ & 書上三維\\", r"\midrule", *rows, r"\midrule", *law_rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_chain():
    """一維 Ising 鏈（K = 1）的 G(r) = t^{|r|}（長條，r = −8..8）與 e^{−|r|/ξ}（虛線）。"""
    K = WC.K1D
    xi = WC.ising1d_xi(K)
    out = [r"\begin{tikzpicture}[x=0.33cm,y=1.55cm,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw (-8.8,0) -- (8.8,0) node[right] {$r$};"]
    for r in range(-8, 9):
        out.append(rf"\fill[cBar!70] ({r - 0.34:.2f},0) rectangle ({r + 0.34:.2f},{WC.ising1d_G(r, K):.4f});")
    xs = np.linspace(-8.4, 8.4, 57)
    out.append(r"\draw[cExact,densely dashed,thick] plot coordinates {" +
               " ".join(f"({x:.3f},{math.exp(-abs(x) / xi):.4f})" for x in xs) + "};")
    out.append(r"\foreach \x in {-8,-4,0,4,8} {\node[below] at (\x,0) {$\x$};}")
    out.append(rf"\node[cBar,anchor=west] at (1.6,1.0) {{$G(r)=t^{{|r|}}$，$t=\tanh1={_s(math.tanh(K), 4)}$}};")
    out.append(rf"\node[cExact,anchor=west] at (2.6,0.78) {{虛線 $e^{{-|r|/\xi}}$，$\xi={_s(xi, 4)}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_well():
    """Landau 雙井 F = r m²/2 + m⁴/4 − h m（u = 1）：r = 1、0、−1（藍，由淺到深）與 r = −1、h = 0.1（橙虛線）。"""
    ms = np.linspace(-2, 2, 81)
    out = [r"\begin{tikzpicture}[x=1.35cm,y=1.4cm,>=Stealth,font=\scriptsize,baseline=(current bounding box.center)]",
           r"\draw[->] (-2.15,0) -- (2.35,0) node[right] {$m$};",
           r"\draw[->] (0,-0.6) -- (0,1.68) node[above] {$F(m)$};",
           r"\foreach \x in {-2,-1,1,2} {\draw (\x,0) -- (\x,-0.04);}",
           r"\node[below] at (-2,-0.02) {$-2$}; \node[below] at (2,-0.02) {$2$};",
           r"\foreach \y in {1} {\draw (0,\y) -- (-0.05,\y) node[left] {$\y$};}",
           r"\begin{scope}", r"\clip (-2.1,-0.6) rectangle (2.1,1.5);"]
    for r, col in zip(WC.LANDAU_R, ("cBar!35", "cBar!65", "cBar")):
        pts = [(m, WC.landau_F(m, r)) for m in ms if WC.landau_F(m, r) < 2.2]
        out.append(rf"\draw[{col},thick] plot coordinates {{{' '.join(f'({a:.3f},{b:.4f})' for a, b in pts)}}};")
    h = WC.LANDAU_H
    pts = [(m, WC.landau_F(m, -1.0, 1.0, h)) for m in ms if WC.landau_F(m, -1.0, 1.0, h) < 2.2]
    out.append(rf"\draw[cWalk,densely dashed,thick] plot coordinates {{{' '.join(f'({a:.3f},{b:.4f})' for a, b in pts)}}};")
    out.append(r"\end{scope}")
    m1, f1 = WC.landau_closed(-1.0)
    out.append(rf"\fill[cExact] ({-m1:.3f},{f1:.4f}) circle (1.6pt); \fill[cExact] ({m1:.3f},{f1:.4f}) circle (1.6pt);")
    out.append(rf"\node[cExact,below] at ({-m1:.3f},{f1 - 0.03:.4f}) {{$m^*=\pm1$，$F_{{\min}}={_s(f1, 2)}$}};")
    mh = WC.landau_min(-1.0, 1.0, h)
    fh = WC.landau_F(mh, -1.0, 1.0, h)
    out.append(rf"\shade[ball color=cWalk] ({mh:.4f},{fh + 0.07:.4f}) circle (2.6pt);")
    out.append(rf"\node[cWalk,anchor=west] at ({mh + 0.1:.3f},{fh - 0.02:.4f}) {{$h={h}$：$m^*={_s(mh, 4)}$}};")
    out.append(r"\node[cBar!50,anchor=south] at (-1.28,1.5) {$r=1$};")         # 曲線在 F = 1.5 出框的位置
    out.append(r"\node[cBar!75,anchor=south east] at (1.57,1.5) {$r=0$};")
    out.append(r"\node[cBar,anchor=south] at (1.91,1.5) {$r=-1$};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_zmf():
    Vs = WC.ZMF_V
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r"$V$（$r=-1$，$u=1$） & " + " & ".join(map(str, Vs)) + r"\\", r"\midrule",
                      r"$\langle m\rangle$（$h=0.01$） & " + " & ".join(_s(WC.zmf_avg(V, h=WC.ZMF_H), 4) for V in Vs) + r"\\",
                      r"$\langle|m|\rangle$（$h=0$） & " + " & ".join(_s(WC.zmf_avg(V, what="absm"), 4) for V in Vs) + r"\\",
                      r"$-\ln Z/V$（$h=0$） & " + " & ".join(f"${_s(WC.zmf_avg(V, what='f'), 4)}$" for V in Vs) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_ex94():
    names = (r"\beta", r"\Delta", r"\alpha", r"1/\delta", r"\gamma")
    theory = (F(1, 2), F(3, 2), 0, F(1, 3), 1)
    return "\n".join([r"\begin{tabular}{@{}l*{5}{c}@{}}", r"\toprule",
                      r"練習 9.4 & " + " & ".join(f"${n}$" for n in names) + r"\\", r"\midrule",
                      r"數值斜率 & " + " & ".join(_s(v, 4) for v in WC.ex94_slopes()) + r"\\",
                      r"理論 & " + " & ".join(f"${_fr(t)}$" for t in theory) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_saw():
    rows = WC.saw_table()
    return "\n".join([r"\begin{tabular}{@{}l*{6}{r}@{}}", r"\toprule",
                      r"$N$ & " + " & ".join(str(N) for N, _, _ in rows) + r"\\", r"\midrule",
                      r"$c_N$ & " + " & ".join(str(cN) for _, cN, _ in rows) + r"\\",
                      r"$c_N\langle R^2\rangle_N$ & " + " & ".join(str(s) for _, _, s in rows) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_coef():
    Ns = WC.COEF_N
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r"精確 $\div$ 漸近，$N=$ & " + " & ".join(map(str, Ns)) + r"\\", r"\midrule",
                      r"$s=\gamma=\frac{43}{32}$ & " + " & ".join(_s(WC.coef_ratio(N, WC.S_GAMMA), 4) for N in Ns) + r"\\",
                      r"$s=\gamma+2\nu=\frac{91}{32}$ & " + " & ".join(_s(WC.coef_ratio(N, WC.S_GAMMA2NU), 4) for N in Ns) + r"\\",
                      r"係數比 $\div$ (9.41) & " + " & ".join(_s(WC.r2_ratio(N, WC.S_GAMMA, F(3, 4), WC.ZC_2D), 4)
                                                          for N in Ns) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def _mini(sub, a, b, x0, y0, label):
    """一個小梯子（格距 0.42 cm）：灰色底圖、路徑（橙）、閉圈（紅，帶 n）、端點（藍）、下方標權重。"""
    s = 0.42
    xy = {v: (x0 + s * px, y0 + s * py) for v, (px, py) in LAD_XY.items()}
    out = [rf"\draw[gray!35] ({xy[i][0]:.2f},{xy[i][1]:.2f}) -- ({xy[j][0]:.2f},{xy[j][1]:.2f});" for i, j in WC.LADDER[1]]
    path, loops = _path_part(sub, a)
    for es, col in ((path, "cWalk"), (loops, "cWrong")):
        out += [rf"\draw[{col},very thick] ({xy[i][0]:.2f},{xy[i][1]:.2f}) -- ({xy[j][0]:.2f},{xy[j][1]:.2f});" for i, j in es]
    out += [rf"\fill[gray!60] ({px:.2f},{py:.2f}) circle (0.8pt);" for px, py in xy.values()]
    for v in (a, b):
        if v is not None:
            out.append(rf"\fill[cBar] ({xy[v][0]:.2f},{xy[v][1]:.2f}) circle (1.6pt);")
    out.append(rf"\node[below,inner sep=1.5pt] at ({x0 + s:.2f},{y0:.2f}) {{{label}}};")
    return out


def fig_ladder():
    """2×3 梯子：Π(1 + K S_i·S_j) 展開後存活的項。上兩列是 G_03、G_05 的分子（路徑＋不相交的閉圈），下列是分母 Z（閉圈）。"""
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]"]
    rows = (((0, 3), r"$G_{03}$ 的分子"), ((0, 5), r"$G_{05}$ 的分子"), ((None, None), r"分母 $Z$"))
    for k, ((a, b), lab) in enumerate(rows):
        y0 = -1.08 * k
        out.append(rf"\node[anchor=west] at (0,{y0 + 0.21:.2f}) {{{lab}}};")
        sets = _ladder_sets(a, b) if a is not None else _ladder_sets()
        for j, (sub, nl) in enumerate(sets):
            e = len(sub)
            w = ("n" if nl else "") + ("K" if e == 1 else (f"K^{{{e}}}" if e else ""))
            out += _mini(sub, a, b, 1.75 + 1.3 * j, y0, f"${w or '1'}$")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_ladder():
    rows = []
    for p, saw in zip(WC.LADDER_PAIRS, (r"$K+K^3+K^5$", r"$3K^3+K^5$")):
        rows.append(rf"$G_{{{p[0]}{p[1]}}}$ & {_s(WC.ladder_G_ising(0.5, p), 4)} & {_s(WC.ladder_G_xy(0.5, p), 4)} & "
                    rf"{_s(WC.ladder_G(0, 0.5, p), 5)} & {saw}\\")
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"$K=0.5$ & $n=1$ & $n=2$ & $n\to0$ & 自迴避路徑\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_trotter():
    t1 = [WC.trotter_err(e, 1) for e in WC.TROTTER_EPS]
    t2 = [WC.trotter_err(e, 2) for e in WC.TROTTER_EPS]
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r"$\epsilon$ & " + " & ".join(map(str, WC.TROTTER_EPS)) + r"\\", r"\midrule",
                      r"一階誤差 & " + " & ".join(f"${_g(x)}$" for x in t1) + r"\\",
                      r"乘上 $e^{\frac{i\epsilon^2}4[\phi,\nabla^2]}$ & " + " & ".join(f"${_g(x)}$" for x in t2) + r"\\", r"\midrule",
                      r"$t=1$ 切成 $M$ 段 & " + " & ".join(f"$M={M}$" for M in WC.TROTTER_M) + r"\\",
                      r"\quad 總誤差 & " + " & ".join(f"${_g(WC.trotter_global(M))}$" for M in WC.TROTTER_M) + r"\\",
                      r"\bottomrule", r"\end{tabular}"])


def tab_weak():
    w = [WC.weak_saw_count(8, u) for u in WC.WEAK_U]
    cells = [f"{w[0]:.0f}", _s(w[1], 2), _s(w[2], 2), f"{w[3]:.0f}"]
    book = _g(WC.weak_saw_count(8, 1.0, sign=+1), 2)
    return "\n".join([r"\begin{tabular}{@{}l*{4}{r}r@{}}", r"\toprule",
                      r"$u$ & 0 & 1 & 4 & $\infty$ & 書上 $u=1$\\", r"\midrule",
                      r"$\sum_we^{-(u/2)n_w}$ & " + " & ".join(cells) + rf" & {{\color{{cWrong}}${book}$}}\\",
                      r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    (r"p.~210，Table 9.1 三維欄", r"$\alpha=\alpha'=0.12$", r"$\approx0.19$（(9.16)(9.18) 都得 0.194）"),
    (r"p.~210，Table 9.1 第二個關聯長度", r"$\xi=\int C\rho^2d^d\rho/\!\int C\,d^d\rho$", r"$\xi^2=\int C\rho^2d^d\rho/\!\int C\,d^d\rho$"),
    (r"p.~218，(9.37) 之後", r"關聯長度 $\propto K_1/(z_c-z)^\nu$", r"$\xi^2\approx K_1/(z_c-z)^{2\nu}$"),
    (r"p.~218，(9.38)$\to$(9.39)", r"$z_c^{-(N+s)}N^{\delta-1}/\Gamma(s)$", r"$z_c^{-(N+s)}N^{s-1}/\Gamma(s)$"),
    (r"p.~219，(9.40)$\to$(9.41)", r"must be equal to $V\Gamma_N/z_0$", r"$V\Gamma_N/z_0^N$"),
    (r"p.~219，(9.41)", r"$K_1z_c^{-\gamma}\,\Gamma(\gamma)N^{2\nu}/\Gamma(\gamma+2\nu)$",
     r"$K_1z_c^{-2\nu}\,\Gamma(\gamma)N^{2\nu}/\Gamma(\gamma+2\nu)$"),
    (r"p.~225，(9.70)", r"$e^{\nabla^2\epsilon}$ 的核 $(2\pi)^{-d/2}e^{-\frac d2|\Delta\vec R|^2/\epsilon}$",
     r"$e^{\frac12\nabla^2\epsilon}$ 的核 $(2\pi\epsilon)^{-d/2}e^{-|\Delta\vec R|^2/2\epsilon}$"),
    (r"p.~226，(9.72)", r"$\int_0^\infty$", r"$\int_{-\infty}^{\infty}$"),
    (r"p.~226，(9.75)", r"$-\int\frac d2\big|\frac{d\vec R}{ds}\big|^2ds\ \frac u4\iint\delta$",
     r"$-\int\frac d2\big|\frac{d\vec R}{ds}\big|^2ds-\frac u4\iint\delta$"),
    (r"p.~230，(9:S-19)", r"$\partial^3/\partial x_r\partial x_s\partial x_t$；括號裡的 $\langle x_rx_s\rangle$",
     r"$\partial^3/\partial h_r\partial h_s\partial h_t$；$\langle x_tx_u\rangle$"),
]


def tab_errata():
    assert len(ERRATA) == 10
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("NineFigMs", fig_ms),
    ("NineTabFlory", tab_flory),
    ("NineTabMs", tab_ms),
    ("NineFigPM", fig_pm),
    ("NineFigTrace", fig_trace),
    ("NineTabFlip", tab_flip),
    ("NineTabField", tab_field),
    ("NineTabNineOne", tab_nineone),
    ("NineTabLaws", tab_laws),
    ("NineFigChain", fig_chain),
    ("NineFigWell", fig_well),
    ("NineTabZmf", tab_zmf),
    ("NineTabExNineFour", tab_ex94),
    ("NineTabSaw", tab_saw),
    ("NineTabCoef", tab_coef),
    ("NineFigLadder", fig_ladder),
    ("NineTabLadder", tab_ladder),
    ("NineTabTrotter", tab_trotter),
    ("NineTabWeak", tab_weak),
    ("NineTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch9()
    parts = ["% figs.tex —— 由 ch9/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch9/make_figs.py（第 9 章，動畫 W901–W906）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    check_tex_decimals()                                   # 圖表裡用到的小數也已記進 OK
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
