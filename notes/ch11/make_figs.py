r"""
make_figs.py —— 隨機漫步導讀第 11 章（Rudnick & Gaspari《Elements of the Random Walk》第 11 章，動畫 W1101–W1106）的圖與表

在 notes/ 或 notes/ch11/ 裡執行：
    python ch11/make_figs.py
會先檢查 ch11/part11.tex 寫到的每一個數字（check_numbers_ch11；最後再把 part11.tex 裡出現的每一個小數和檢查過的數比對一遍），
再寫出 ch11/figs.tex（part11.tex 第一行 \input{ch11/figs}）。圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算）。
數學一律由 codes/utils/walk_rg.py 算出（W1106 的對照用 walk_crit.EXPONENTS）；例子與動畫 W1101–W1106 相同。巨集名稱都以 \Eleven 開頭
（不和前幾章的 \Fig…、\Two…、…、\Ten… 撞名）；part11.tex 的 \label 用 eleven-…。顏色沿用導言區：cBar 精確值、cWalk 路線、
cExact 更正、cGauss 近似、cWrong 筆誤。全部約 2 秒（盒子計數與 400 條漫步由 walk_rg 快取）。
"""
from __future__ import annotations

import itertools
import math
import re
import sys
import warnings
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent            # notes/ch11
sys.path.insert(0, str(HERE.parent.parent / "codes"))
warnings.filterwarnings("ignore")                  # scipy 的 IntegrationWarning 不影響結果

from utils import walk_rg as RG    # noqa: E402   （套件模組：walk_rg 的相對 import 才找得到 walk_saw、walk_crit）
from utils import walk_crit as WC  # noqa: E402

OUT = HERE / "figs.tex"
TEX = HERE / "part11.tex"
F = Fraction
OK = set()                                         # 檢查過的小數字串（part11.tex 的小數都要在這裡）
EXP = ("alpha", "beta", "gamma", "delta", "nu")


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


def _rg2_enum(N):
    """⟨R_g²⟩：二維 4^N 條普通漫步全部列舉（第三種算法，只給小 N）。"""
    steps = ((1, 0), (-1, 0), (0, 1), (0, -1))
    tot = F(0)
    for w in itertools.product(steps, repeat=N):
        p = np.vstack([[0, 0], np.cumsum(w, axis=0)]).astype(np.int64)
        s1, s2 = p.sum(axis=0), int((p * p).sum())
        tot += F(s2, N + 1) - F(int(s1 @ s1), (N + 1) ** 2)
    return tot / 4 ** N


# =============================================================================
# 數字檢查（part11.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch11():
    OK.update({"0.5", "0.3", "0.7", "0.25", "0.008", "0.125", "0.01", "0.001", "0.05", "0.1", "0.03", "0.003",
               "1.5", "0.588", "0.9", "0.4", "2.25", "4.5", "2.5", "9.28", "0.6", "1.2"})
    # ---- W1101 碎形
    assert _s(RG.LN4_LN3, 4) == "1.2619" and _s(RG.LN2_LN3, 4) == "0.6309"
    _all([RG.koch_counts(k)[1] for k in range(6)], 4, ["1.0000", "1.3333", "1.7778", "2.3704", "3.1605", "4.2140"])
    _all([RG.cantor_counts(k)[1] for k in range(6)], 4, ["1.0000", "0.6667", "0.4444", "0.2963", "0.1975", "0.1317"])
    assert [RG.koch_counts(k)[0] for k in range(6)] == [1, 4, 16, 64, 256, 1024]
    assert len(RG.koch_points(4)) == 257 and len(RG.cantor_intervals(5)) == 32
    assert sum(b - a for a, b in RG.cantor_intervals(5)) == F(2, 3) ** 5
    kp = RG.koch_points(4)
    assert np.all(np.hypot(*kp[65:].T) >= 1 / 3 - 1e-12)            # 〔注〕第 2 段以後都在 R = 1/3 外（30° 錐）
    assert np.all(np.abs(np.degrees(np.arctan2(kp[1:, 1], kp[1:, 0]))) <= 30 + 1e-9)
    _all([n for _, _, n in RG.mass_table("koch")], 4, ["1.0000", "0.2500", "0.0625", "0.0156", "0.0039"])
    _all([n for _, _, n in RG.mass_table("cantor")], 4, ["1.0000", "0.5000", "0.2500", "0.1250", "0.0625"])
    assert (_s(RG.mass_radius_dim("koch"), 4), _s(RG.mass_radius_dim("cantor"), 4)) == ("1.2619", "0.6309")
    assert _s(RG.mass_radius_dim("koch", [2.0 ** -j for j in range(9)]), 4) == "1.2571"
    assert [n for _, n in RG.box_counts("koch", RG.BOX_K, 3)] == [22, 83, 325, 1291]
    assert (_s(RG.box_dim("koch"), 4), _s(RG.box_dim("cantor"), 4)) == ("1.2362", "0.6254")
    _all([RG.replica_dim(*RG.REPLICAS[k]) for k in ("line", "square", "koch", "cantor", "sierpinski")], 4,
         ["1.0000", "2.0000", "1.2619", "0.6309", "1.5850"])
    assert all(RG.rg2_exact(N) == _rg2_enum(N) for N in range(1, 6))
    _all([RG.rw_rg2(N) for N in RG.RW_NS], 2, ["11.23", "44.65", "163.53", "696.53"])
    _all([RG.rg2_exact(N) for N in RG.RW_NS], 2, ["10.83", "42.83", "170.83", "682.83"])
    assert (_s(RG.rw_df(), 4), _s(RG.rw_df(exact=True), 4)) == ("2.0258", "2.0069")
    assert RG.saw_df(2) == F(4, 3) and _s(RG.saw_df(3), 2) == "1.70" and RG.saw_df(4) == 2
    _all([RG.saw_df_eff(N) for N in RG.SAW_NEFF], 4, ["1.3848", "1.3788", "1.3740"])
    # ---- W1102 被動的重整化
    assert _s(RG.landau_m(-0.25, 0.0), 4) == "0.5000"
    _all([RG.landau_m(0.0, h) for h in (0.008, 0.125)], 4, ["0.2000", "0.5000"])
    _all([RG.landau_scaled(x) for x in (1, 2, 10)], 4, ["0.6823", "1.0000", "2.0000"])
    _all([RG.landau_m(r, h) for r, h in RG.COLLAPSE_POINTS], 4, ["0.5000", "0.7071", "1.0000"])
    assert all(_close(h / r ** 1.5, 2) for r, h in RG.COLLAPSE_POINTS)
    f, book, right = RG.book1116_limit(1e-4, 1e-3)
    assert (_s(f, 4), _s(right, 4)) == ("9.9667", "10.0000") and _close(book, 1e-5, 1e-12)
    assert RG.mf_rescale_check() < 1e-12 and all(RG.gauss_rescale_check(d) < 1e-12 for d in (1, 2, 3, 4))
    g3, g4 = RG.gauss_exponents(3), RG.gauss_exponents(4)
    assert [g3[k] for k in ("beta", "gamma", "delta", "alpha")] == [F(1, 4), 1, 5, F(1, 2)]
    assert [g4[k] for k in ("beta", "gamma", "delta", "alpha")] == [F(1, 2), 1, 3, 0]
    assert all(g4[k] == WC.EXPONENTS["mf"][k] for k in EXP)
    assert RG.scaling_combos(4) == (F(3, 2), F(3, 2), F(1, 2), F(1, 2))
    assert RG.scaling_combos(3) == (F(3, 2), F(5, 4), F(1, 2), F(1, 4))
    # ---- W1103 一維 Gaussian 鏈
    num, closed = RG.square_check()
    assert (_s(num, 6), _s(closed, 6)) == ("1.224819", "1.224819")
    rp, ap, _ = RG.decimate_gauss(0.5, 1.0)
    assert (_s(rp, 1), _s(ap, 1), _s(rp / ap, 2), _s(RG.rr_gauss(0.5), 2)) == ("0.9", "0.4", "2.25", "2.25")
    assert RG.ring_decimation_check(8, 0.5, 1.0, 0.0) < 1e-10 and RG.ring_decimation_check() < 1e-10
    assert _s(RG.ring_decimation_check(8, 0.5, 1.0, 0.0, book=True), 4) == "1.8432"
    assert _close(RG.ring_decimation_check(8, 0.5, 1.0, 0.0, book=True), 2 * math.log(2 * math.pi / 2.5))
    assert _s(RG.rr_curve()[-1][0], 4) == "1.7417" and _close(RG.rr_curve()[-1][1], 10)
    _all(RG.gauss_iterate(RG.COBWEB_R0, 4), 4, ["0.0500", "0.2025", "0.8510", "4.1282", "33.5553"])
    assert RG.gauss_fixed_points() == (0.0, -3.0)
    tab = RG.small_r_table()
    _all([t[1] for t in tab], 4, ["1.0000", "1.0001", "1.0005", "1.0021", "1.0086"])
    _all([t[2] for t in tab], 4, ["1.0000"] * 5)
    _all([t[3] for t in tab], 4, ["1.0000", "1.0001", "1.0005", "1.0021", "1.0085"])
    assert _s(RG.hh_gauss(0.0, 1.0), 4) == "2.8284" and _close(RG.hh_gauss(0.0, 1.0), 2 ** 1.5)
    assert all(RG.gauss_m_invariance(r, h) < 1e-12 for r, h in ((0.3, 0.2), (0.01, 0.001)))
    b47 = RG.book1147_m(0.01, 0.001, 1)
    assert (_s(b47[0], 4), _s(b47[1], 4), _s(b47[2], 4)) == ("0.0500", "0.1000", "0.1000")
    _all(RG.chain_lambdas(), 4, ["2.0000", "1.5000", "0.5000"])
    # ---- W1104 一維 Ising
    assert (_s(2 * math.cosh(2.0), 4), _s(math.cosh(2.0), 4)) == ("7.5244", "3.7622")
    assert tuple(_s(x, 4) for x in RG.book1153(1.0)) == ("3.7622", "7.5244", "7.5244")
    assert (_s(RG.ising_A(1.0), 4), _s(RG.ising_Kp(1.0), 4)) == ("3.8793", "0.6625") and RG.ising_pair_check(1.0) < 1e-12
    _all(RG.ising_iterate(1.0, 5), 4, ["1.0000", "0.6625", "0.3501", "0.1137", "0.0128", "0.0002"])
    _all(RG.ex118(1.0, 0.1), 4, ["3.8897", "0.6602", "0.0482", "0.1964"])
    assert RG.ex118_check(1.0, 0.1) < 1e-12
    A, K1, K2, Hp = RG.ex118(1.0, 0.1)
    assert _close(RG.ring_Z_brute(8, 1.0, 0.1), A ** 4 * RG.ring_Z_brute(4, K1, Hp), 1e-6)   # Z_8(K,H) = A^4 Z_4(K1',H')
    Z = {N: (RG.ring_Z_brute(N, 1.0), RG.ring_Z_transfer(N, 1.0), RG.ring_Z_decimation(N, 1.0),
             RG.ring_Z_decimation(N, 1.0, book=True)) for N in RG.RING_NS}
    assert all(_close(z[0], z[1], 1e-6 * z[1]) and _close(z[2], z[1], 1e-6 * z[1]) for z in Z.values())
    assert [_s(Z[N][1], 4) for N in RG.RING_NS] == ["121.2329", "9160.4387", "68584536.5912"]
    assert [_s(Z[N][3], 4) for N in RG.RING_NS] == ["18.5498", "143.2705", "11787.1848"]
    Ks = RG.ising_iterate(1.0, 2)                                    # 4 自旋環的連乘 A(K_1)^2 A(K_2) · 2e^{K_3}
    assert (_s(RG.ising_A(1.0) ** 2, 4), _s(RG.ising_A(Ks[1]), 4)) == ("15.0488", "2.8383")
    assert _close(RG.ising_A(Ks[1]), 2 * math.exp(Ks[2]))
    assert _close(RG.ising_A(1.0) ** 2 * RG.ising_A(Ks[1]) * 2 * math.exp(Ks[2]), Z[4][1], 1e-9)
    c, s = 2 * math.cosh(1.0), 2 * math.sinh(1.0)
    assert (_s(c, 4), _s(s, 4), _s(c ** 4, 4), _s(s ** 4, 4)) == ("3.0862", "2.3504", "90.7140", "30.5189")
    _all(RG.x_iterate(RG.x_of_K(1.0), 5), 4, ["0.1353", "0.2658", "0.4965", "0.7966", "0.9747", "0.9997"])
    assert all(_close(RG.x_map(RG.x_of_K(K)), RG.x_of_K(RG.ising_Kp(K))) for K in (0.3, 1.0, 2.0))
    assert RG.x_fixed_points() == (0.0, 1.0) and (RG.x_slope(0.0), RG.x_slope(1.0)) == (2.0, 0.0)
    _all([RG.xi_exact(K) for K in RG.XI_K], 4, ["3.6719", "27.2960", "201.7140"])
    _all([math.exp(2 * K) / 2 for K in RG.XI_K], 4, ["3.6945", "27.2991", "201.7144"])
    _all([RG.xi_ratio(K) for K in RG.XI_K], 4, ["0.4969", "0.4999", "0.5000"])
    assert all(max(abs(e) for e in RG.tanh_rg_check(K)) < 1e-9 for K in RG.XI_K)
    _all([RG.free_energy_ratio(K) for K in RG.FREE_K], 4, ["1.1269", "1.0091", "1.0000"])
    assert all(_close(RG.free_energy_series(K), math.log(2 * math.cosh(K))) for K in RG.FREE_K)
    assert _s(RG.free_energy_ratio(1.0, 16), 4) == "1.1277"
    assert _close(RG.free_energy_ratio(1.0, 16), math.log(Z[16][0]) / 16)
    assert _s(math.sqrt(2), 4) == "1.4142" and _close(RG.ising_A(10.0) / math.exp(10.0), math.sqrt(2), 1e-8)
    # ---- W1105 動量殼
    assert (RG.shell_fraction(2), RG.shell_fraction(3)) == (F(3, 4), F(7, 8))
    at = RG.alpha_table()
    _all([t[1] for t in at], 4, ["0.0000", "0.7854", "1.5708", "2.3562", "3.1416"])
    _all([t[2] for t in at], 4, ["0.0000", "0.2929", "1.0000", "1.7071", "2.0000"])
    _all([t[3] for t in at], 4, ["0.0000", "0.3084", "1.2337", "2.7758", "4.9348"])
    _all([t[4] + 0.0 for t in at], 4, ["0.0000", "-0.5858", "-2.0000", "-3.4142", "-4.0000"])
    err, real, book = RG.fourier_check()
    assert err < 1e-12 and (_s(real, 4), _s(book, 4)) == ("3.0394", "-4.3508")
    r3, h3 = RG.shell_rg(1.0, 1.0, 3)
    assert r3 == 4 and _s(h3, 4) == "5.6569" and all(_close(a, b) for a, b in zip((r3, h3), RG.shell_rg_steps(1.0, 1.0, 3)))
    fl = [RG.flow(0.01, 0.01, 3, t) for t in RG.FLOW_ELL]
    assert [(_s(r, 4), _s(h, 4)) for r, h in fl] == [("0.0739", "0.1218"), ("0.5460", "1.4841")]
    assert all(_close(a, b, 1e-8) for t in RG.FLOW_ELL for a, b in zip(RG.flow(0.01, 0.01, 3, t), RG.flow_numeric(0.01, 0.01, 3, t)))
    assert all(_close(a, b) for a, b in zip(RG.flow(0.3, 0.2, 3, math.log(2)), RG.shell_rg(0.3, 0.2, 3)))
    assert all(RG.pde_residual(d) < 1e-7 and RG.char_check(d) < 1e-12 for d in (1, 3, 4))
    land = lambda r, h: math.sqrt(r) * RG.landau_scaled(h / r ** 1.5)           # Landau 的組合放進 d = 3 的 (11.87)
    e = 1e-5
    res = (-0.5 * land(0.4, 0.3) + 2 * 0.4 * (land(0.4 + e, 0.3) - land(0.4 - e, 0.3)) / (2 * e)
           + 2.5 * 0.3 * (land(0.4, 0.3 + e) - land(0.4, 0.3 - e)) / (2 * e))
    assert abs(res) > 0.01
    # ---- W1106 不動點與指數
    assert RG.toy_eigen() == [(2, (2, 1)), (-3, (1, -2))]
    assert all(_close(a, b, 1e-9) for c in RG.TOY_STARTS for a, b in zip(RG.toy_flow(*c, 0.7), RG.toy_flow_expm(*c, 0.7)))
    xs = [RG.xi_from_flow(0.01, L) for L in (2, 1)]
    assert [(_s(a, 4), _s(b, 4)) for a, b in xs] == [("2.3026", "10.0000"), ("4.6052", "100.0000")]
    assert all(_close(a, b, 1e-6) for L in (2, 1) for a, b in zip(RG.xi_from_flow(0.01, L), RG.xi_from_flow_numeric(0.01, L)))
    want = {"gauss3": (F(1, 2), F(1, 4), 1, 5, F(1, 2)), "gauss4": (0, F(1, 2), 1, 3, F(1, 2)),
            "ising2d": (0, F(1, 8), F(7, 4), 15, 1), "saw2d": (F(1, 2), F(5, 64), F(43, 32), F(91, 5), F(3, 4))}
    for k, (Lr, Lh, d) in RG.RG_CASES.items():
        E = RG.rg_exponents(Lr, Lh, d)
        assert tuple(E[q] for q in EXP) == want[k] and E == RG.rg_exponents_widom(Lr, Lh, d)
        assert E["alpha"] + 2 * E["beta"] + E["gamma"] == 2
        assert RG.lambdas_from(E["nu"], E["beta"], d) == (Lr, Lh)
    for k, w in (("gauss4", "mf"), ("ising2d", "ising2d"), ("saw2d", "saw2d")):
        assert all(RG.rg_exponents(*RG.RG_CASES[k])[q] == WC.EXPONENTS[w][q] for q in EXP)
    assert all(RG.rg_exponents(*RG.gauss_lambdas(d), d)[q] == RG.gauss_exponents(d)[q] for d in (3, 4) for q in EXP)


def check_tex_decimals():
    """part11.tex（去掉註解）裡的每一個小數都必須是上面檢查過的數。"""
    if not TEX.exists():
        print("（尚無 part11.tex，略過小數比對）")
        return
    txt = re.sub(r"(?<!\\)%.*", "", TEX.read_text(encoding="utf-8"))
    txt = re.sub(r"[\d.]+\\linewidth|[\d.]+pt\b|[\d.]+em\b|[\d.]+cm\b|\\arraystretch\}\{[\d.]+\}", "", txt)   # 版面尺寸
    eqno = re.compile(r"11\.\d{1,3}")                    # 本章的式號、節號、圖號、練習號（11.1–11.114）
    bad = sorted({m for m in re.findall(r"\d+\.\d+", txt) if m not in OK and not eqno.fullmatch(m)})
    assert not bad, f"part11.tex 有未檢查的小數：{bad}"


# =============================================================================
# 圖與表
# =============================================================================
ARROW = r"postaction={decorate},decoration={markings,mark=at position 0.55 with {\arrow{Stealth[length=3pt]}}}"


def _poly(pts, sx=1.0, sy=1.0, x0=0.0, y0=0.0, d=3):
    return " -- ".join(f"({x0 + sx * x:.{d}f},{y0 + sy * y:.{d}f})" for x, y in pts)


def fig_koch():
    """Koch 曲線第 0–4 代（一列）；第 4 代畫出以左端為圓心、半徑 L/3（虛線）與 L（點線）的圓，前四分之一（綠）在 L/3 圓內。"""
    W, gap = 2.95, 0.33
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]"]
    for k in range(5):
        x0 = k * (W + gap)
        pts = RG.koch_points(k)
        if k == 4:
            q = len(pts) // 4 + 1
            out.append(r"\begin{scope}")
            out.append(rf"\clip ({x0 - 0.05:.3f},-0.12) rectangle ({x0 + W + 0.12:.3f},{0.98:.3f});")
            out.append(rf"\draw[cExact,densely dashed] ({x0:.3f},0) circle ({W / 3:.3f});")
            out.append(rf"\draw[gray,densely dotted] ({x0:.3f},0) circle ({W:.3f});")
            out.append(r"\end{scope}")
            out.append(r"\draw[cBar,line join=round] " + _poly(pts[q - 1:], W, W, x0) + ";")
            out.append(r"\draw[cExact,thick,line join=round] " + _poly(pts[:q], W, W, x0) + ";")
            out.append(rf"\fill ({x0:.3f},0) circle (1.2pt);")
            out.append(rf"\node[cExact,anchor=south] at ({x0 + W / 6:.3f},0.42) {{$R=L/3$}};")
        else:
            out.append(r"\draw[cBar,line join=round] " + _poly(pts, W, W, x0) + ";")
        out.append(rf"\node[anchor=north] at ({x0 + W / 2:.3f},-0.1) {{$k={k}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_cantor():
    """Cantor 集第 0–5 代，一列一代。"""
    W, dy = 6.3, 0.34
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]"]
    for k in range(6):
        y = -k * dy
        segs = " ".join(rf"({W * float(a):.3f},{y:.3f}) -- ({W * float(b):.3f},{y:.3f})" for a, b in RG.cantor_intervals(k))
        out.append(rf"\draw[cWalk,line width=2.6pt] {segs};")
        out.append(rf"\node[anchor=east] at (-0.1,{y:.3f}) {{$k={k}$}};")
        out.append(rf"\node[anchor=west] at ({W + 0.1:.3f},{y:.3f}) {{{2 ** k}}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_fractal():
    km = {k: n for k, _, n in RG.mass_table("koch")}
    cm = {k: n for k, _, n in RG.mass_table("cantor")}
    rows = []
    for k in range(6):
        nk, lk, _ = RG.koch_counts(k)
        nc, lc, _ = RG.cantor_counts(k)
        mk = _s(km[k], 4) if k in km else ""
        mc = _s(cm[k], 4) if k in cm else ""
        rows.append(rf"{k} & {nk} & {_s(lk, 4)} & {nc} & {_s(lc, 4)} & {mk} & {mc}\\")
    return "\n".join([r"\begin{tabular}{@{}c*{6}{r}@{}}", r"\toprule",
                      r" & \multicolumn{2}{c}{Koch} & \multicolumn{2}{c}{Cantor} & \multicolumn{2}{c}{$n(R=3^{-k})$}\\",
                      r"$k$ & 段數 & 總長 & 段數 & 總長 & Koch & Cantor\\", r"\midrule", *rows, r"\midrule",
                      rf"斜率 & & $\to\infty$ & & $\to0$ & {_s(RG.mass_radius_dim('koch'), 4)} & {_s(RG.mass_radius_dim('cantor'), 4)}\\",
                      r"\bottomrule", r"\end{tabular}"])


def _axes(xmax, ymax, sx, sy, xt, yt, xl, yl, xmin=0.0, ymin=0.0, fmt="{:g}"):
    out = [rf"\draw[->,gray] ({xmin * sx:.3f},{ymin * sy:.3f}) -- ({xmax * sx + 0.25:.3f},{ymin * sy:.3f}) node[right,black] {{{xl}}};",
           rf"\draw[->,gray] ({max(xmin, 0.0) * sx:.3f},{ymin * sy:.3f}) -- "
           rf"({max(xmin, 0.0) * sx:.3f},{ymax * sy + 0.25:.3f}) node[above,black] {{{yl}}};"]
    yx = max(xmin, 0.0) * sx
    for t in xt:
        out.append(rf"\draw[gray] ({t * sx:.3f},{ymin * sy:.3f}) -- ++(0,-0.06) node[below,black] {{{fmt.format(t)}}};")
    for t in yt:
        out.append(rf"\draw[gray] ({yx:.3f},{t * sy:.3f}) -- ++(-0.06,0) node[left,black] {{{fmt.format(t)}}};")
    return out


def fig_collapse():
    """塌縮：左 m 對 h（r = 0.25, 0.5, 1），右 m/r^{1/2} 對 h/r^{3/2}（x ≤ 4）疊成一條 y + y³ = x。"""
    sty = {0.25: "cBar,thick", 0.5: "cWalk,thick,densely dashed", 1.0: "cExact,thick,densely dotted"}
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    sx, sy = 2.5, 1.7
    out += _axes(2, 1.25, sx, sy, (0, 0.5, 1, 1.5, 2), (0.5, 1), r"$h$", r"$m$")
    for r in RG.COLLAPSE_R:
        out.append(rf"\draw[{sty[r]}] " + _poly(RG.landau_curve(r), sx, sy) + ";")
    for (r, h) in RG.COLLAPSE_POINTS:
        m = RG.landau_m(r, h)
        out.append(rf"\fill[black] ({h * sx:.3f},{m * sy:.3f}) circle (1.3pt);")
    for i, (r, c) in enumerate(((0.25, "cBar"), (0.5, "cWalk"), (1, "cExact"))):
        out.append(rf"\node[{c},anchor=west] at ({1.4 * sx:.3f},{(0.48 - 0.16 * i) * sy:.3f}) {{$r={r}$}};")
    X0 = 6.3
    sx2, sy2 = 1.25, 1.45
    out.append(rf"\begin{{scope}}[xshift={X0}cm]")
    out += _axes(4, 1.4, sx2, sy2, (0, 1, 2, 3, 4), (0, 0.5, 1), r"$x=h/r^{3/2}$", r"$y=m/r^{1/2}$")
    for r in RG.COLLAPSE_R:
        out.append(rf"\draw[{sty[r]}] " + _poly(RG.collapse_curve(r, 40, 4 * r ** 1.5), sx2, sy2) + ";")
    out.append(rf"\draw[gray,densely dotted] ({2 * sx2:.3f},0) -- ({2 * sx2:.3f},{sy2:.3f}) -- (0,{sy2:.3f});")
    out.append(rf"\fill[black] ({2 * sx2:.3f},{sy2:.3f}) circle (1.5pt);")
    out.append(rf"\node[anchor=west,align=left] at ({2.1 * sx2:.3f},{0.62 * sy2:.3f}) {{三個黑點\\都到 $(2,1)$}};")
    out.append(rf"\node[anchor=west] at ({3.0 * sx2:.3f},{1.5 * sy2:.3f}) {{$y+y^3=x$}};")
    out.append(r"\end{scope}")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_chain():
    """Fig. 11.6–11.7（更正的座標 r/a）：r'' = r(r+4a)/a、對角虛線、r_0 = 0.05 的蛛網。"""
    sx, sy = 2.55, 0.29
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(2, 10, sx, sy, (0, 0.5, 1, 1.5, 2), (2, 4, 6, 8, 10), r"$r/a$", r"$r''/a$")
    out.append(rf"\draw[gray,densely dashed] (0,0) -- ({2 * sx:.3f},{2 * sy:.3f}) node[above,black] {{$r''=r$}};")
    out.append(r"\draw[cBar,thick] " + _poly(RG.rr_curve(), sx, sy) + ";")
    xe, ye = RG.rr_curve()[-1]
    out.append(rf"\node[cBar,anchor=east] at ({xe * sx - 0.05:.3f},{ye * sy:.3f}) {{$r''=r(r+4a)/a$}};")
    cw = RG.cobweb_points()
    out.append(r"\draw[cWalk,thick,->] " + _poly(cw, sx, sy) + ";")
    for x, y in cw[1::2]:
        out.append(rf"\fill[cWalk] ({x * sx:.3f},{y * sy:.3f}) circle (1.1pt);")
    rs = RG.gauss_iterate(RG.COBWEB_R0, 3)
    out.append(rf"\node[cWalk,anchor=west] at ({rs[2] * sx + 0.05:.3f},{rs[3] * sy:.3f}) "
               rf"{{$r_3={_s(rs[3], 4)}$}};")
    out.append(rf"\node[cWalk,anchor=west,align=left] at ({0.06 * sx:.3f},{7.2 * sy:.3f}) "
               rf"{{$r_0={_s(rs[0], 2)}\to{_s(rs[1], 4)}$\\$\to{_s(rs[2], 4)}\to{_s(rs[3], 4)}$}};")
    out.append(r"\fill (0,0) circle (1.4pt);")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_small_r():
    rows = [rf"{k} & {_s(a, 4)} & {_s(b, 4)} & {_s(c, 4)}\\" for k, a, b, c in RG.small_r_table()]
    return "\n".join([r"\begin{tabular}{@{}crrr@{}}", r"\toprule",
                      r"$k$ & $\dfrac{r^{(k)}}{4^kr}$ & $\dfrac{h^{(k)}}{2^{3k/2}h}$ & $\dfrac{\text{自旋因子}}{2^{k/2}}$\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def fig_xmap():
    """Ex. 11.9：x' = 2x/(1+x²) 的蛛網（x_0 = e^{−2}）。"""
    s = 3.5
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(1, 1, s, s, (0, 0.5, 1), (0.5, 1), r"$x=e^{-2K}$", r"$x'$")
    out.append(rf"\draw[gray,densely dashed] (0,0) -- ({s:.3f},{s:.3f});")
    xs = np.linspace(0, 1, 40)
    out.append(r"\draw[cBar,thick] " + _poly([(x, RG.x_map(x)) for x in xs], s, s) + ";")
    it = RG.x_iterate(RG.x_of_K(1.0), 5)
    pts = [(it[0], 0.0)]
    for k in range(5):
        pts.append((it[k], it[k + 1]))
        if k < 4:
            pts.append((it[k + 1], it[k + 1]))
    out.append(r"\draw[cWalk,thick,->] " + _poly(pts, s, s) + ";")
    out.append(r"\fill (0,0) circle (1.5pt);")
    out.append(rf"\fill ({s:.3f},{s:.3f}) circle (1.5pt);")
    out.append(r"\node[anchor=north west,align=left] at (0.08,-0.3) {$x=0$：$T=0$，斜率 2（不穩定）};")
    out.append(rf"\node[anchor=west,align=left] at ({s + 0.08:.3f},{s:.3f}) {{$x=1$：$T=\infty$\\斜率 0（穩定）}};")
    out.append(rf"\node[cBar,anchor=west] at ({0.55 * s:.3f},{0.22 * s:.3f}) {{藍：$x'=\frac{{2x}}{{1+x^2}}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_ising():
    rows = []
    for N in RG.RING_NS:
        z = RG.ring_Z_transfer(N, 1.0)
        b = RG.ring_Z_decimation(N, 1.0, book=True)
        rows.append(rf"${N}$ & {_s(z, 4)}\,$\checkmark$ & {{\color{{cWrong}}{_s(b, 4)}}}\\")
    return "\n".join([r"\begin{tabular}{@{}rrr@{}}", r"\toprule",
                      r"$N$ & 暴力＝轉移矩陣＝連乘 & 書上指數\\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def fig_shell():
    """Fig. 11.8 的二維版：布里淵區每邊縮一半，外圈殼層（3/4 的模態）被積掉。"""
    a = 1.25
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]",
           rf"\fill[cWalk!25,even odd rule] ({-a},{-a}) rectangle ({a},{a}) ({-a / 2},{-a / 2}) rectangle ({a / 2},{a / 2});",
           rf"\draw[cWalk,thick] ({-a},{-a}) rectangle ({a},{a});",
           rf"\draw[cBar,thick] ({-a / 2},{-a / 2}) rectangle ({a / 2},{a / 2});",
           rf"\draw[gray,->] ({-a - 0.2},0) -- ({a + 0.3},0) node[right,black] {{$q_x$}};",
           rf"\draw[gray,->] (0,{-a - 0.2}) -- (0,{a + 0.3}) node[above,black] {{$q_y$}};",
           rf"\node[anchor=north west] at ({a},0) {{$\frac\pi c$}};",
           rf"\node[anchor=north west,cBar] at ({a / 2},0) {{$\frac\pi{{2c}}$}};",
           rf"\node[align=center] at (0,{-a - 0.55}) {{殼層（橙）被積掉：\\$1-2^{{-d}}={_fr(RG.shell_fraction(2))}$}};",
           r"\end{tikzpicture}"]
    return "\n".join(out)


def _ell_max(r0, h0, d=3):
    return min(math.log(1 / abs(r0)) / 2, math.log(1.5 / h0) / ((d + 2) / 2))


FLOW_STARTS = ((0.01, 0.01), (-0.01, 0.01), (0.03, 0.003), (-0.03, 0.003))


def fig_flow():
    """(11.82)(11.83) 的流線（d = 3）：起點 (±0.01, 0.01)、(±0.03, 0.003)，走到 |r| = 1 或 h = 1.5。"""
    sx, sy = 2.1, 2.1
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(1, 1.5, sx, sy, (-1, -0.5, 0, 0.5, 1), (0.5, 1, 1.5), r"$r$", r"$h$", xmin=-1)
    for r0, h0 in FLOW_STARTS:
        L = RG.flow_line(r0, h0, 3, _ell_max(r0, h0))
        out.append(rf"\draw[cBar,thick,{ARROW}] " + _poly(L, sx, sy) + ";")
        out.append(rf"\fill[cBar] ({r0 * sx:.3f},{h0 * sy:.3f}) circle (0.9pt);")
    for t in RG.FLOW_ELL:
        r, h = RG.flow(0.01, 0.01, 3, t)
        out.append(rf"\fill[cWalk] ({r * sx:.3f},{h * sy:.3f}) circle (1.4pt) node[right] {{$\ell={int(t)}$}};")
        _s(r, 4), _s(h, 4)
    out.append(rf"\node[anchor=west,align=left] at ({0.62 * sx:.3f},{1.0 * sy:.3f}) {{$h\propto|r|^{{5/4}}$\\都離開原點}};")
    out.append(r"\fill (0,0) circle (1.5pt);")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_alpha():
    rows = [rf"{k} & {_s(q, 4)} & {_s(a, 4)} & {_s(b, 4)} & {{\color{{cWrong}}{_s(c, 4) if k else '0'}}}\\"
            for k, q, a, b, c in RG.alpha_table()]
    rows = [r.replace("{0.0000}", "{0}").replace(" 0.0000 ", " 0 ") for r in rows]
    return "\n".join([r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule",
                      r"$k$ & $q=\frac{2\pi k}8$ & $a(1-\cos q)$ & $\frac12aq^2$ & 書上\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_toy():
    """示意矩陣 A = [[1,2],[2,−2]] 的流向圖：相關方向 (2,1)（Λ = 2）、無關方向 (1,−2)（Λ = −3）。"""
    s = 1.0
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,scale=0.85,baseline=(current bounding box.center)]",
           r"\draw[cExact,densely dashed] (-3.1,-1.55) -- (3.1,1.55) node[above left,align=right] {相關 $(2,1)$\\$\Lambda=2$};",
           r"\draw[cWrong,densely dashed] (-0.85,1.7) -- (1.1,-2.2) node[right,align=left] {無關 $(1,-2)$\\$\Lambda=-3$（臨界線）};"]
    for c1, c2 in RG.TOY_STARTS:
        L = RG.toy_flow_line(c1, c2)
        out.append(rf"\draw[cBar,thick,{ARROW}] " + _poly(L, s, s) + ";")
        x, y = L[0]
        out.append(rf"\fill[cBar] ({x:.3f},{y:.3f}) circle (1pt);")
    out.append(r"\fill (0,0) circle (1.8pt) node[left=2pt] {不動點};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_exp():
    cases = ("gauss3", "gauss4", "ising2d", "saw2d")
    E = {k: RG.rg_exponents(*RG.RG_CASES[k]) for k in cases}
    head = (r" & \shortstack{Gaussian\\$d=3$} & \shortstack{Gaussian\\$d=4$} & \shortstack{二維\\Ising} & "
            r"\shortstack{二維\\自迴避}\\")
    lam = [r"$\Lambda_r$ & " + " & ".join(f"${_fr(RG.RG_CASES[k][0])}$" for k in cases) + r"\\",
           r"$\Lambda_h$ & " + " & ".join(f"${_fr(RG.RG_CASES[k][1])}$" for k in cases) + r"\\"]
    lab = {"alpha": r"$\alpha=2-\frac d{\Lambda_r}$", "beta": r"$\beta=\frac{d-\Lambda_h}{\Lambda_r}$",
           "gamma": r"$\gamma=\frac{2\Lambda_h-d}{\Lambda_r}$", "delta": r"$\delta=\frac{\Lambda_h}{d-\Lambda_h}$",
           "nu": r"$\nu=\frac1{\Lambda_r}$"}
    rows = [lab[q] + " & " + " & ".join(f"${_fr(E[k][q])}$" for k in cases) + r"\\" for q in EXP]
    ck = r"$\checkmark$"
    tail = [r"$\alpha+2\beta+\gamma$ & " + " & ".join("2" for k in cases if E[k]["alpha"] + 2 * E[k]["beta"] + E[k]["gamma"] == 2) + r"\\",
            r"＝第 9 章 & & " + " & ".join(ck for _ in cases[1:]) + r"\\"]
    return "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule", head, r"\midrule", *lam, r"\midrule", *rows,
                      r"\midrule", *tail, r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    (r"p.~263，(11.16) 後", r"$f(h/r^{3/2})\to h^{1/3}r$", r"$f(h/r^{3/2})\to h^{1/3}r^{-1/2}$"),
    (r"p.~265，(11.31)", r"$\big(\frac{2\pi}{r+2a}\big)^{N/2}$", r"$\big(\frac{2\pi}{r+2a}\big)^{N/4}$"),
    (r"pp.~267--268，Fig.~11.6--11.7", r"軸 $r/\sqrt a$、$r''/\sqrt a$", r"$r/a$、$r''/a$"),
    (r"p.~267，正文", r"the reduced temperature, $t''$", r"$r''$"),
    (r"p.~269，(11.46)(11.47)", r"$m=2^{-n/2}f(r2^{2n},h2^{3n/2})$", r"$m=2^{+n/2}f(r2^{2n},h2^{3n/2})$"),
    (r"p.~270，(11.53)", r"$\cosh K(S_2+S_4)=2\cosh2K$（$S_2=S_4$）", r"$2\cosh K(S_2+S_4)=2\cosh2K$"),
    (r"p.~271，(11.63)", r"$\cdots A(K_j)^{N/2^{j+1}}\cdots$", r"$\cdots A(K_j)^{N/2^j}\cdots$"),
    (r"p.~273，(11.66)(11.67)", r"$\alpha(q)=-2a(1-\cos qc)\sim ac^2q^2$", r"$\alpha(q)=a(1-\cos qc)\sim\frac12ac^2q^2$"),
]


def tab_errata():
    assert len(ERRATA) == 8
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("ElevenFigKoch", fig_koch),
    ("ElevenFigCantor", fig_cantor),
    ("ElevenTabFractal", tab_fractal),
    ("ElevenFigCollapse", fig_collapse),
    ("ElevenFigChain", fig_chain),
    ("ElevenTabSmallR", tab_small_r),
    ("ElevenFigXmap", fig_xmap),
    ("ElevenTabIsing", tab_ising),
    ("ElevenFigShell", fig_shell),
    ("ElevenFigFlow", fig_flow),
    ("ElevenTabAlpha", tab_alpha),
    ("ElevenFigToy", fig_toy),
    ("ElevenTabExp", tab_exp),
    ("ElevenTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch11()
    parts = ["% figs.tex —— 由 ch11/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch11/make_figs.py（第 11 章，動畫 W1101–W1106）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    check_tex_decimals()                                   # 圖表裡用到的小數也已記進 OK
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
