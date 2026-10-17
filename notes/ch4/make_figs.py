r"""
make_figs.py —— 隨機漫步導讀第 4 章（Rudnick & Gaspari《Elements of the Random Walk》第 4 章，動畫 W401–W408）的圖與表

在 notes/ 或 notes/ch4/ 裡執行：
    python ch4/make_figs.py
會先檢查 ch4/part4.tex 寫到的每一個數字（check_numbers_ch4），再寫出 ch4/figs.tex（part4.tex 第一行 \input{ch4/figs}）。
圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots）。數學一律由 codes/utils/walk_bc.py（必要時 walk_gf.py、walk_rec.py）算出；
例子與動畫 W401–W408 相同。巨集名稱都以 \Four 開頭（不和第 1 章的 \Fig…、第 2 章的 \Two…、第 3 章的 \Three… 撞名）；
part4.tex 的 \label 用 four-…。顏色沿用導言區：cBar 精確值、cGauss 近似、cExact 更正、cWrong 筆誤與吸收牆、cWalk 路線、cTail 鏡像路線。
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent            # notes/ch4
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_bc as BC  # noqa: E402
from utils import walk_gf as GF  # noqa: E402
from utils import walk_rec as R  # noqa: E402

OUT = HERE / "figs.tex"
F = Fraction
S3 = math.sqrt(3)
PI = math.pi


def _r(xs, k):
    return [round(x, k) for x in xs]


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣）。"""
    return f"{x:.{d}f}"


def _simpson(f, a, b, n=4000):
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) + 4 * sum(f(a + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
                    + 2 * sum(f(a + 2 * i * h) for i in range(1, n // 2)))


# =============================================================================
# 數字檢查（part4.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch4():
    # ---- 4.1.1 鏡像走路人（W401）
    W = BC.walks(4, 1)
    assert len(W) == 2 ** 4 == 16 and sum(BC.touches_wall(w) for w in W) == 10
    surv = BC.survivors_brute(4, 1)
    assert surv == BC.survivors_image(4, 1) == BC.survivors_dp(4, 1) == {1: 2, 3: 3, 5: 1} and sum(surv.values()) == 6
    assert [GF.count1(4, d) for d in (0, 2, 4, 6)] == [6, 4, 1, 0]                       # C(4; x−1)、C(4; x+1)
    assert [GF.count1(4, x - 1) - GF.count1(4, x + 1) for x in (1, 3, 5)] == [2, 3, 1]
    assert sum(w[-1] > 0 for w in W) == 11 == 6 + 4 + 1 and len(BC.reflection_pairs(4, 1)) == 5 and 11 - 5 == 6
    assert sum(GF.count1(4, x + 1) for x in (1, 3, 5)) == 5                               # 從 −1 出發、終點 > 0
    assert BC.reflect_first([1, 0, 1, 2, 3]) == [-1, 0, 1, 2, 3]
    assert BC.reflect_first([1, 2, 1, 0, 1]) == [-1, -2, -1, 0, 1]
    assert BC.survival_prob(4, 1) == F(6, 16)
    assert BC.survival_prob(10, 3) == F(21, 32) == F(210 + 252 + 210, 1024) and 21 / 32 == 0.65625
    assert [math.comb(10, k) for k in (4, 5, 6)] == [210, 252, 210]
    assert _s(float(BC.survival_prob(100, 3)), 4) == "0.2356" and _s(float(BC.survival_prob(1000, 3)), 5) == "0.07557"
    assert [_s(BC.survival_asym(N, 3), d) for N, d in ((10, 4), (100, 4), (1000, 5))] == ["0.7569", "0.2394", "0.07569"]
    assert [_s(BC.survival_3d(N, 3), 4) for N in (100, 1000)] == ["0.3958", "0.1305"]
    assert abs(BC.survival_3d(100, 3) - BC.survival_3d_dp(100, 3)) < 1e-12
    assert [_s(BC.survival_3d_asym(N, 3), 4) for N in (100, 1000)] == ["0.4146", "0.1311"]
    assert abs(BC.survival_3d_asym(300, 3) - BC.survival_asym(100, 3)) < 1e-12            # sc 的 N 步 ≈ 一維的 N/3 步
    assert BC.wall_reduction_1d(10) == F(1067, 512) == BC.wall_reduction_brute(10) and _s(1067 / 512, 3) == "2.084"
    assert _s(float(BC.wall_reduction_1d(100)), 3) == "7.499"
    assert [_s(BC.wall_reduction_asym(N), 3) for N in (10, 100)] == ["2.523", "7.979"]
    N = 30.0                                                                            # (4.4c)(4.6)(4.7)(4.8)
    assert abs(_simpson(lambda x: x * math.exp(-3 * x * x / (2 * N)), 0, 80) - N / 3) < 1e-9
    a = 1.7
    lhs = _simpson(lambda x: math.exp(-3 * (x - a) ** 2 / (2 * N)) - math.exp(-3 * (x + a) ** 2 / (2 * N)), 0, 80)
    assert abs(lhs - _simpson(lambda x: math.exp(-3 * x * x / (2 * N)), -a, a)) < 1e-9
    assert abs(2 * (N / 3) * (2 * PI * N / 3) ** -0.5 - math.sqrt(2 * N / (3 * PI))) < 1e-12
    # ---- 4.1.2 兩面牆（W402）
    seq = [BC.box_dp(N, 6, 3) for N in range(5)]
    assert seq == [{3: 1}, {2: 1, 4: 1}, {1: 1, 3: 2, 5: 1}, {2: 3, 4: 3}, {1: 3, 3: 6, 5: 3}]
    assert [sum(d.values()) for d in seq] == [1, 2, 4, 6, 12]
    d10 = BC.box_dp(10, 6, 3)
    assert [d10.get(x, 0) for x in range(1, 6)] == [81, 0, 162, 0, 81] == [BC.box_images(10, 6, 3, x) for x in range(1, 6)]
    assert [round(BC.box_modes(10, 6, 3, x)) for x in range(1, 6)] == [81, 0, 162, 0, 81]
    assert sum(d10.values()) == 324 and 2 ** 10 == 1024 and _s(324 / 1024, 3) == "0.316"
    assert _r(BC.box_eigs(6), 4) == _r([S3, 1, 0, -1, -S3], 4) and BC.box_eig_residual(6) < 1e-12
    assert round(S3 ** 10) == 243 and F(2, 6) * (243 + 243) == 162 and F(2, 6) * (F(243, 2) + F(243, 2)) == 81
    assert abs(math.sin(PI / 6) - 0.5) < 1e-15 and abs(math.sin(5 * PI / 6) - 0.5) < 1e-15
    assert BC.box_image_points(6, 3) == ([-9, 3, 15], [-15, -3, 9])
    assert GF.count1(10, 0) == 252 and GF.count1(10, 6) == GF.count1(10, -6) == 45 and 252 - 45 - 45 == 162
    for z, m, x2 in ((0.3, 1, 2), (0.2, 2, 1), (0.25, 3, 4)):                           # (4.12)(4.13) 更正：sin q x2
        assert abs(BC.box_amplitude(z, 6, m, x2) - BC.box_amplitude_formula(z, 6, m, x2)) < 1e-12
    # ---- 4.1.3 長短時間（W403）
    tot = [BC.box_total(N, 6) for N in range(11)]
    assert tot == [5, 8, 14, 24, 42, 72, 126, 216, 378, 648, 1134]
    assert all(tot[N + 2] == 3 * tot[N] for N in range(1, 9)) and tot[2] != 3 * tot[0] and 5 * 2 == 10
    assert _s(BC.decay_rate(6), 4) == "0.1438" and abs(math.cos(PI / 6) - S3 / 2) < 1e-15
    for L, m in ((10, 1), (10, 3), (20, 1), (6, 5)):                                    # (4.16) 的格子版
        assert abs(sum(math.sin(m * PI * x / L) for x in range(1, L)) - 1 / math.tan(m * PI / (2 * L))) < 1e-12
    assert abs(sum(math.sin(2 * PI * x / 10) for x in range(1, 10))) < 1e-12
    ps = [sum(1 / (2 * k + 1) ** 2 for k in range(M + 1)) for M in range(4)]
    assert _r(ps, 4) == [1, 1.1111, 1.1511, 1.1715]
    assert _s(BC.odd_inv_squares(1000), 4) == _s(PI ** 2 / 8, 4) == "1.2337" and _s(BC.zeta2(1000), 4) == "1.6449"
    assert abs(PI ** 2 / 6 * 3 / 4 - PI ** 2 / 8) < 1e-15 and abs(8 / PI ** 2 * PI ** 2 / 8 - 1) < 1e-15
    assert abs(BC.box_total_modes(0, 10) - 9) < 1e-9                                    # 短時間：2^N (L − 1)
    rt = {L: [BC.decay_rate(L, k) for k in ("lattice", "cont", "book")] for L in (10, 20)}
    assert [f"{v:.4g}" for v in rt[10][:2] + rt[20][:2]] == ["0.05018", "0.04935", "0.01239", "0.01234"]
    assert [f"{rt[L][2]:.3g}" for L in (10, 20)] == ["0.0987", "0.0247"]
    assert _s(rt[10][0] / rt[20][0], 2) == "4.05" and abs(rt[10][1] / rt[20][1] - 4) < 1e-12
    assert _s(2 * 10 ** 2 / PI ** 2, 2) == "20.26" and _s(1 / rt[10][0], 2) == "19.93"
    sv = [float(BC.box_survival(N, 10)) for N in (50, 100, 200)]
    assert _s(sv[0], 6) == "0.072103" and _s(sv[1], 6) == "0.005865" and f"{sv[2]:.3e}" == "3.881e-05"
    lo = [BC.box_survival_long(100, 10, k) for k in ("lattice", "cont", "book")]
    assert _s(lo[0], 6) == "0.005861" and _s(lo[1], 6) == "0.005830" and f"{lo[2]:.2e}" == "4.19e-05"
    L = 10
    low = 2 / L * (1 / math.tan(PI / (2 * L))) ** 2 * math.cos(PI / L) ** 100 / (L - 1)
    assert abs(lo[0] / low - 1) < 1e-12
    assert _s(8 / PI ** 2, 4) == "0.8106" and _s(PI ** 2 / 2, 3) == "4.935"
    # ---- 4.1.4 球與立方體（W404）
    assert [round(BC.j0_zero(n) / PI, 9) for n in (1, 2, 3)] == [1, 2, 3]
    assert [_s(BC.sphere_weight(n, 1), 3) for n in (1, 2)] == ["2.546", "0.637"] and _s(8 / PI, 3) == "2.546"
    assert all(abs(BC.sphere_weight(n, 1) - BC.sphere_weight_numeric(n, 1)) < 1e-9 for n in (1, 2))
    assert _s(BC.sphere_weight_book(1, 1), 3) == "0.516"
    assert abs(8 / PI * PI ** 2 / 6 - 4 / 3 * PI) < 1e-14 and _s(4 / 3 * PI, 3) == "4.189"
    assert len(BC.ball_sites(8)) == 2103 and _s(BC.ball_rate(8), 4) == "0.0242"
    assert _s(BC.sphere_rate(8), 4) == "0.0257" and _s(BC.sphere_rate(8, book=True), 3) == "0.154"
    assert round(BC.sphere_rate(8, book=True) / BC.sphere_rate(8)) == 6
    assert [BC.cube_total(N, 6) for N in (0, 1, 2)] == [125, 600, 2970] and 125 * 6 - 6 * 25 == 600 and 5 ** 3 == 125
    assert abs(BC.cube_total_modes(2, 6) * 36 - 2970) < 1e-8
    P = BC.PREFACTOR
    assert [_s(P[k], 3) for k in ("interval", "cube", "sphere")] == ["0.811", "0.533", "0.608"]
    assert abs(P["cube"] - P["interval"] ** 3) < 1e-15 and abs(P["sphere"] - (8 / PI) / (4 / 3 * PI)) < 1e-15
    assert _s(BC.same_volume_rate_ratio(), 3) == "1.155"
    assert abs(BC.same_volume_rate_ratio() - 3 * (3 / (4 * PI)) ** (2 / 3)) < 1e-12
    # ---- 4.2.1 穩定態、Fick、帳篷（W405）
    assert _s(BC.laplacian_coeff(), 4) == "0.1667"
    w = BC.tent_lattice(5, 2)
    assert w == [0, F(6, 5), F(12, 5), F(8, 5), F(4, 5), 0] == [BC.tent_formula(5, 2, x) for x in range(6)]
    assert [float(v) for v in w[1:5]] == [1.2, 2.4, 1.6, 0.8]
    assert w[1] - w[0] == F(6, 5) and w[5] - w[4] == F(-4, 5) and F(-4, 5) - F(6, 5) == -2 == -1 / F(1, 2)
    assert w[2] == F(1, 2) * (w[1] + w[3]) + 1 and w[1] == F(1, 2) * (w[0] + w[2])
    assert BC.wall_fluxes(5, 2) == (F(3, 5), F(2, 5)) == (w[1] / 2, w[4] / 2) and F(3, 5) + F(2, 5) == 1
    assert BC.absorb_left(5, 2) == BC.absorb_left_solve(5, 2) == F(3, 5) and abs(BC.absorb_left_dp(5, 2) - 0.6) < 1e-12
    assert BC.absorb_left(1000, 2) == F(998, 1000)
    # ---- 4.2.2 點源（W406）
    assert _s(BC.time_integral(2.0), 5) == "0.03979" and abs(1 / (8 * PI) - BC.time_integral(2.0)) < 1e-9
    assert abs(BC.K_const(2.0) - BC.K_const_numeric(2.0)) < 1e-9 and abs(BC.K_const(2.0) - math.sqrt(2 * PI)) < 1e-15
    G = [BC.sc_green(n, 0, 0) for n in (0, 1, 2, 3, 4, 5, 10)]
    assert _r(G[:6], 4) == [1.5164, 0.5164, 0.2573, 0.1653, 0.1217, 0.0966] and f"{G[6]:.4g}" == "0.04787"
    nG = [n * g for n, g in zip((1, 2, 3, 4, 5, 10), G[1:])]
    assert _r(nG, 4) == [0.5164, 0.5147, 0.4958, 0.4869, 0.4830, 0.4787]
    assert _s(BC.sc_green_asym(1), 4) == "0.4775" == _s(3 / (2 * PI), 4) and _s(BC.sc_green_book(1), 4) == "0.2387"
    assert abs(G[1] - (G[0] - 1)) < 1e-9 and abs(G[0] - R.G_C["sc"]) < 1e-9 and abs(BC.sc_laplace_residual(2, 1, 0)) < 1e-10
    assert _s(R.G_C["sc"], 6) == "1.516386"
    # ---- 4.2.3 吸收牆旁的點源（W406）
    assert abs(BC.plane_flux(3, 1)) < 1e-8 and abs(BC.plane_flux(0.5, 1) + 1) < 1e-8
    assert BC.plane_flux_exact(3, 1) == 0 and BC.plane_flux_exact(0.5, 1) == -1
    assert abs(BC.dipole_flux_integral(1.0)) < 1e-6
    assert abs(BC.conc_wall(50, 3, 0, 1) / BC.conc_dipole(50, 3, 0, 1) - 1) < 1e-2
    # ---- 4.2.4–4.2.5 吸收球、電容（W407）
    assert _s(BC.sphere_flux(1), 3) == "12.566" and all(abs(BC.sphere_flux_numeric(r, 1) - 4 * PI) < 1e-6 for r in (2, 5))
    assert BC.conc_sphere(2, 1) == 0.5
    rmax = BC.cell_rmax(**BC.CELL_EXAMPLE)
    assert BC.CELL_EXAMPLE == {"D": 2e-9, "cinf": 0.25, "q": 1.0} and f"{rmax:.1e}" == "3.9e-05" and round(rmax * 1e6) == 39
    assert abs(rmax - math.sqrt(3 * 2e-9 * 0.25 / 1.0)) < 1e-18 and 2 ** 3 == 8
    assert BC.kelvin_image(4, 1) == (0.25, -0.25) and BC.hit_prob(4, 1) == 0.25 and abs(BC.hit_prob_numeric(4, 1) - 0.25) < 1e-8
    assert abs(BC.sphere_surface_max(4, 1)) < 1e-12 and 1 - 0.25 == 0.75
    assert BC.capacitance_sphere(1) == 1
    # ---- 4.2.6 Berg–Purcell（W408，書外）
    bp = BC.BP_EXAMPLE
    assert bp == {"r0": 1e-6, "a": 1e-9}
    nh = BC.bp_half(**bp)
    assert round(nh) == 3142 and _s(100 * BC.bp_area_fraction(nh, **bp), 3) == "0.079"
    assert _r([BC.berg_purcell(n, **bp) for n in (100, 1000, 3142, 1e4, 1e5)], 3) == [0.031, 0.241, 0.5, 0.761, 0.970]
    assert all(abs(BC.berg_purcell(n, **bp) - BC.berg_purcell_series(n, **bp)) < 1e-12 for n in (100, 1e4))
    assert abs(BC.disk_capacitance(1) - 2 / PI) < 1e-15 and abs(BC.disk_capacitance_numeric(1) - 2 / PI) < 1e-9
    assert _s(2 / PI, 4) == "0.6366"
    # ---- 4.3 吸收邊界（W408）
    cat = [BC.catalan(k) for k in range(6)]
    assert cat == [1, 1, 2, 5, 14, 42] == [BC.half_counts(2 * k, 0, 0) for k in range(6)]
    assert all(BC.half_counts(n, x, y) == BC.half_counts_image(n, x, y) for n in range(9) for x in range(4) for y in range(4))
    assert [math.comb(2 * k, k) - math.comb(2 * k, k + 1) for k in range(6)] == cat
    assert abs(BC.half_gf_residual(0.3, 1)) < 1e-12
    assert [BC.half_visits_critical(x, y) for x, y in ((0, 0), (1, 3), (2, 2))] == [2, 4, 6]
    assert abs(float(BC.tent_lattice(400, 2)[3]) - 4) < 0.05                             # 帳篷 L → ∞ 的極限
    # ---- 練習
    z = 0.3                                                                             # 4.1：半直線的鏡像
    lam = (1 - math.sqrt(1 - 4 * z * z)) / (2 * z)
    for x1, x2 in ((1, 2), (3, 3)):
        ser = sum(BC.survivors_dp(n, x1).get(x2, 0) * z ** n for n in range(120))
        assert abs(ser - (lam ** abs(x1 - x2) - lam ** (x1 + x2)) / math.sqrt(1 - 4 * z * z)) < 1e-12
        assert abs(ser - BC.half_gf(z, x1 - 1, x2 - 1)) < 1e-12
    N = 12.0                                                                            # 4.2：對 y、z 積分
    gy = _simpson(lambda y: math.exp(-3 * y * y / (2 * N)), -60, 60)
    assert abs(gy * gy * (2 * PI * N / 3) ** -1.5 - (2 * PI * N / 3) ** -0.5) < 1e-12
    q, x = 0.7, 2.3                                                                     # 4.4
    assert abs(math.sin(q * (x + 1)) + math.sin(q * (x - 1)) - 2 * math.cos(q) * math.sin(q * x)) < 1e-14
    for m in range(1, 6):
        for mp in range(1, 6):
            s = sum(2 / 6 * math.sin(m * PI * x / 6) * math.sin(mp * PI * x / 6) for x in range(1, 6))
            assert abs(s - (m == mp)) < 1e-12
    assert abs((PI ** 2 / 8) ** 3 * (8 / PI ** 2) ** 3 - 1) < 1e-15 and _s(512 / PI ** 6, 3) == "0.533"   # 4.5


# =============================================================================
# 表與圖
# =============================================================================
def fig_mirror():
    """Fig. 4.1 的樣子：第一次碰牆之前的那段對 0 鏡射（兩個例子）。"""
    pairs = {tuple(p): q for p, q in BC.reflection_pairs(4, 1)}
    out = [r"\begin{tikzpicture}[x=1.1cm,y=0.48cm,>=Stealth,dot/.style={circle,inner sep=1.2pt}]"]
    for i, (path, lab) in enumerate((([1, 0, 1, 2, 3], ("above right", "below right")),
                                      ([1, 2, 1, 0, 1], ("anchor=north west,inner sep=1pt", "anchor=south west,inner sep=1pt")))):
        img = pairs[tuple(path)]
        assert img == BC.reflect_first(path)
        k0 = next(k for k, v in enumerate(path) if v <= 0)
        arrow = r"\to".join(str(v) for v in path) + r"\ \leftrightarrow\ " + r"\to".join(str(v) for v in img)
        out += [rf"\begin{{scope}}[xshift={7.4 * i}cm]",
                r"\draw[gray!25] (0,-2) grid (4,3);",
                r"\draw[cWrong,very thick] (0,0) -- (4.25,0) node[right,font=\scriptsize] {懸崖 $0$};",
                r"\foreach \y in {-2,-1,0,1,2,3} \node[left,font=\scriptsize] at (0,\y) {$\y$};",
                r"\foreach \k in {0,...,4} \node[below,font=\scriptsize] at (\k,-2) {$\k$};",
                r"\node[font=\scriptsize] at (2,-3.3) {步數 $k$};",
                rf"\node[font=\scriptsize] at (2,3.7) {{${arrow}$}};"]
        seg_img = " -- ".join(f"({k},{v})" for k, v in enumerate(img[:k0 + 1]))
        seg_org = " -- ".join(f"({k},{v})" for k, v in enumerate(path[:k0 + 1]))
        rest = " -- ".join(f"({k},{v})" for k, v in enumerate(path) if k >= k0)
        out += [rf"\draw[cTail,very thick,dashed] {seg_img};",
                rf"\draw[cWalk,very thick] {seg_org};",
                rf"\draw[cWalk,very thick] {rest};",
                rf"\draw[cTail,thick,densely dashed] {rest};"]
        for k, v in enumerate(path):
            out.append(rf"\node[dot,fill=cWalk] at ({k},{v}) {{}};")
        for k, v in enumerate(img[:k0]):
            out.append(rf"\node[dot,fill=cTail] at ({k},{v}) {{}};")
        out += [rf"\node[{lab[0]},font=\scriptsize,cWalk] at (0.06,{0.9 if i else 1}) {{起點 $x_0=1$}};",
                rf"\node[{lab[1]},font=\scriptsize,cTail] at (0.06,{-0.9 if i else -1}) {{鏡像點 $-x_0=-1$}};",
                rf"\node[circle,draw=cWrong,thick,inner sep=2.6pt] at ({k0},0) {{}};",
                rf"\node[below right,font=\scriptsize,cWrong] at ({k0},0) {{第一次碰牆}};",
                r"\end{scope}"]
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_survive():
    Ns = (10, 100, 1000)
    sp = [BC.survival_prob(N, 3) for N in Ns]
    row1 = r"一維：精確存活機率 & $21/32=0.65625$ & " + " & ".join(
        _s(float(v), 4 if N < 1000 else 5) for N, v in zip(Ns[1:], sp[1:])) + r"\\"
    row2 = r"一維：$x_0\sqrt{2/\pi N}$ & " + " & ".join(
        _s(BC.survival_asym(N, 3), 4 if N < 1000 else 5) for N in Ns) + r"\\"
    row3 = r"sc 格子、吸收平面：精確 & -- & " + " & ".join(_s(BC.survival_3d(N, 3), 4) for N in Ns[1:]) + r"\\"
    row4 = r"三維：$\sqrt{6x_0^2/\pi N}$（(4.5) 後） & -- & " + " & ".join(
        _s(BC.survival_3d_asym(N, 3), 4) for N in Ns[1:]) + r"\\"
    wr = [BC.wall_reduction_1d(N) for N in (10, 100)]
    row5 = (r"一維 (4.8)：牆吃掉的（$\div2^N$）精確 & $\frac{1067}{512}=" + _s(float(wr[0]), 3) + "$ & "
            + _s(float(wr[1]), 3) + r" & \\")
    row6 = r"一維 (4.8)：$\sqrt{2N/\pi}$ & " + " & ".join(_s(BC.wall_reduction_asym(N), 3) for N in (10, 100)) + r" & \\"
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r"起點 $x_0=3$、牆在 $0$ & $N=10$ & $N=100$ & $N=1000$\\", r"\midrule",
                      row1, row2, row3, row4, r"\midrule", row5, row6, r"\bottomrule", r"\end{tabular}"])


def fig_decay():
    L, sx, sy, ymin = 10, 0.06, 0.54, -6
    out = [rf"\begin{{tikzpicture}}[x={sx}cm,y={sy}cm,>=Stealth]",
           rf"\draw[->] (0,{ymin}) -- (214,{ymin}) node[right,font=\scriptsize] {{$N$}};",
           rf"\draw[->] (0,{ymin}) -- (0,0.7) node[above,font=\scriptsize] {{存活比例（對數刻度）}};",
           r"\foreach \n in {0,50,100,150,200} {\draw (\n,-6.06) -- (\n,-5.94); "
           r"\node[below,font=\scriptsize] at (\n,-6) {$\n$};}",
           r"\foreach \e in {0,-1,-2,-3,-4,-5,-6} {\draw (-1.5,\e) -- (1.5,\e); "
           r"\node[left,font=\scriptsize] at (-1.5,\e) {$10^{\e}$};}",
           r"\draw[gray!25] (0,-1) -- (200,-1) (0,-2) -- (200,-2) (0,-3) -- (200,-3) (0,-4) -- (200,-4) (0,-5) -- (200,-5);"]
    for kind, style in (("lattice", "cBar,thick"), ("cont", "cExact,very thick,dashed"), ("book", "cWrong,thick,dash dot")):
        y0, y1 = (math.log10(BC.box_survival_long(N, L, kind)) for N in (0, 200))
        N1 = 200 if y1 >= ymin else 200 * (ymin - y0) / (y1 - y0)
        out.append(rf"\draw[{style}] (0,{y0:.4f}) -- ({N1:.2f},{y0 + (y1 - y0) * N1 / 200:.4f});")
    for N in range(0, 201, 10):
        out.append(rf"\fill[cBar] ({N},{math.log10(float(BC.box_survival(N, L))):.4f}) circle (1.7pt);")
    out += [r"\node[right,font=\scriptsize,cBar] at (88,-0.15) {點：精確計數 $\div 2^N(L-1)$；藍線：格子最低模態 $\propto\cos^N(\pi/L)$};",
            r"\node[right,font=\scriptsize,cExact] at (88,-0.8) {綠虛線：更正 $\frac{8}{\pi^2}e^{-\pi^2N/2L^2}$（和藍線幾乎重合）};",
            r"\node[right,font=\scriptsize,cWrong] at (88,-1.45) {紅點劃線：書上 (4.24) $\frac{8}{\pi^2}e^{-\pi^2N/L^2}$（兩倍陡）};",
            r"\end{tikzpicture}"]
    return "\n".join(out)


def tab_rates():
    rows = []
    for L in (10, 20):
        a, b, c = (BC.decay_rate(L, k) for k in ("lattice", "cont", "book"))
        rows.append(rf"{L} & {a:.4g} & {b:.4g} & {c:.3g}\\")
    ra = BC.decay_rate(10) / BC.decay_rate(20)
    rows.append(rf"\midrule 牆距加倍，衰減變慢 & {ra:.2f} 倍 & 4 倍 & 4 倍\\")
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r"牆距 $L$ & 格子 $-\ln\cos(\pi/L)$ & 更正 $\pi^2/2L^2$ & {\color{cWrong}書上 $\pi^2/L^2$}\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def tab_shapes():
    P = BC.PREFACTOR
    rows = [rf"線段 $0<x<L$（一維） & $2^NL$ & $2^N\frac{{8L}}{{\pi^2}}e^{{-\pi^2N/2L^2}}$ & $8/\pi^2={P['interval']:.3f}$ & $\pi^2/2L^2$\\",
            rf"立方體（邊長 $L$，練習 4.5） & $6^NL^3$ & $6^N\frac{{512L^3}}{{\pi^6}}e^{{-\pi^2N/2L^2}}$ & $512/\pi^6={P['cube']:.3f}$ & $\pi^2/2L^2$\\",
            rf"球（半徑 $R$） & $6^N\cdot\frac43\pi R^3$ & $6^N\frac{{8R^3}}{{\pi}}e^{{-\pi^2N/6R^2}}$ & $6/\pi^2={P['sphere']:.3f}$ & $\pi^2/6R^2$\\"]
    return "\n".join([r"{\small\setlength{\tabcolsep}{4pt}\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                      r"區域（吸收牆） & 短時間總數 & 長時間總數 & 前因子 $\div$ 體積 & 每步衰減率\\", r"\midrule",
                      *rows, r"\bottomrule", r"\end{tabular}}"])


def fig_tent():
    L, x0 = 5, 2
    w = BC.tent_lattice(L, x0)
    fl, fr = BC.wall_fluxes(L, x0)
    out = [r"\begin{tikzpicture}[x=1.75cm,y=0.8cm,>=Stealth]",
           r"\draw[->] (-0.3,0) -- (5.6,0) node[right,font=\scriptsize] {$x$};",
           r"\foreach \x in {0,...,5} \node[below,font=\scriptsize] at (\x,0) {$\x$};"]
    for x in range(1, L):
        v = float(w[x])
        out.append(rf"\fill[cBar!20,draw=cBar] ({x - 0.2},0) rectangle ({x + 0.2},{v:.4f});")
        out.append(rf"\node[font=\scriptsize,cBar] at ({x},{v / 2:.4f}) {{$\frac{{{w[x].numerator}}}{{{w[x].denominator}}}$}};")
    out += [rf"\draw[cExact,thick] (0,0) -- ({x0},{float(w[x0]):.4f}) -- ({L},0);",
            r"\draw[cWrong,line width=2pt] (0,0) -- (0,2.9) (5,0) -- (5,2.9);",
            r"\node[above,font=\scriptsize,cWrong] at (0,2.9) {吸收牆};",
            r"\node[above,font=\scriptsize,cWrong] at (5,2.9) {吸收牆};",
            rf"\draw[->,cWalk,very thick] ({x0},3.45) -- ({x0},{float(w[x0]) + 0.12:.4f});",
            rf"\node[above,font=\scriptsize,cWalk] at ({x0},3.45) {{每步在 $x_0={x0}$ 放 1 人（$s_0=1$）}};",
            r"\node[font=\scriptsize,cExact] at (0.5,1.45) {斜率 $\frac65$};",
            r"\node[font=\scriptsize,cExact] at (3.75,1.6) {斜率 $-\frac45$};",
            rf"\node[right,font=\scriptsize,cExact] at ({x0 + 0.12},{float(w[x0]) + 0.35:.4f}) {{斜率跳 $-s_0/\beta=-2$}};",
            r"\draw[->,cWrong,thick] (0.85,-0.75) -- (0.05,-0.75);",
            rf"\node[below,font=\scriptsize,cWrong] at (0.45,-0.75) {{流進左牆 $w(1)/2=\frac{{{fl.numerator}}}{{{fl.denominator}}}$}};",
            r"\draw[->,cWrong,thick] (4.15,-0.75) -- (4.95,-0.75);",
            rf"\node[below,font=\scriptsize,cWrong] at (4.55,-0.75) {{流進右牆 $w(4)/2=\frac{{{fr.numerator}}}{{{fr.denominator}}}$}};",
            r"\end{tikzpicture}"]
    return "\n".join(out)


def tab_green():
    ns = (0, 1, 2, 3, 4, 5, 10)
    G = [BC.sc_green(n, 0, 0) for n in ns]
    gs = [f"{g:.4f}" if n <= 5 else f"{g:.5f}" for n, g in zip(ns, G)]
    rows = [r"$w=G(z_c;n\hat x)$ & " + " & ".join(gs) + r" & $\to0$\\",
            r"$n\,G$ & -- & " + " & ".join(f"{n * g:.4f}" for n, g in zip(ns[1:], G[1:]))
            + rf" & $\frac{{3}}{{2\pi}}={BC.sc_green_asym(1):.4f}$\\"]
    return "\n".join([r"{\small\begin{tabular}{@{}l" + "c" * 8 + "@{}}", r"\toprule",
                      r"$n$（沿 $x$ 軸） & 0 & 1 & 2 & 3 & 4 & 5 & 10 & $\infty$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}}"])


def tab_bp():
    bp = BC.BP_EXAMPLE
    ns = (100, 1000, 3142, 10 ** 4, 10 ** 5)
    head = ["100", "1000", "3142", "$10^4$", "$10^5$"]
    rows = [r"$J/J_{\max}=\frac{na}{na+\pi r_0}$ & " + " & ".join(f"{BC.berg_purcell(n, **bp):.3f}" for n in ns) + r"\\",
            r"受器占表面 $\frac{na^2}{4r_0^2}$ & " + " & ".join(
                f"{100 * BC.bp_area_fraction(n, **bp):.2g}" + r"\%" for n in ns) + r"\\"]
    return "\n".join([r"\begin{tabular}{@{}lccccc@{}}", r"\toprule",
                      r"受器數 $n$（$r_0=1\,\mu$m、$a=1$ nm） & " + " & ".join(head) + r"\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    (r"p.~75，(4.12)(4.13)", r"$\sqrt{2/L}\,\sin qx_1$", r"$\sqrt{2/L}\,\sin qx_2$"),
    (r"p.~75，(4.15)〔注〕", r"$\sum_{n=0}^{\infty}$", r"$\sum_{n=1}^{\infty}$（$n=0$ 項為 0）"),
    (r"p.~75 末行", r"$0\le z\le L$", r"$0\le x\le L$"),
    (r"p.~77，(4.22)", r"$\frac14\sum_{k=0}^\infty\frac1{k^2}$、$\frac34\sum_{n=0}^\infty\frac1{n^2}$", r"和從 1 起（$k=0$ 項發散）"),
    (r"p.~77，(4.24)", r"$e^{-\pi^2\Delta^2N/L^2}$", r"$e^{-\pi^2\Delta^2N/2L^2}$"),
    (r"p.~80，(4.31)", r"$\frac{8\pi}{R}\cdots\frac2R\frac{j_0j_0}{1-2z+z\Delta^2k^2}$",
     r"$\frac{8\pi}{R}\cdots k_{n_0}^2\frac{j_0j_0}{1-6z+z\Delta^2k^2}$"),
    (r"p.~81，(4.35)", r"$e^{-\Delta^2\pi^2N/R^2}$", r"$e^{-\Delta^2\pi^2n^2N/6R^2}$"),
    (r"p.~81，(4.37)", r"$(6e^{-\Delta^2\pi^2/R^2})^N$", r"$(6e^{-\Delta^2\pi^2/6R^2})^N$"),
    (r"p.~83，(4.40)", r"$\frac{a^2}{3}\nabla^2w$", r"$\frac{a^2}{6}\nabla^2w$"),
    (r"p.~83，(4.41)", r"$\frac{3}{a^2}S$", r"$\frac{6}{a^2}S$"),
    (r"p.~87，(4.58) 前", r"wall in the $x,z$ plane", r"$y,z$ 平面（$x=0$）"),
    (r"p.~91，(4.64) 前", r"$\rho=-4\pi\nabla^2\phi$", r"$\nabla^2\phi=-4\pi\rho$"),
    (r"p.~92", r"$na\ll r_0$", r"$na^2\ll r_0^2$（$a\ll r_0/\sqrt n$）"),
    (r"p.~93", r"$a\ll r/\sqrt n$", r"$a\ll r_0/\sqrt n$"),
    (r"p.~93，(4:S-1) 後", r"starting at $y$ and ending at $y$", r"starting at $x$ and ending at $y$"),
    (r"p.~93，(4:S-2)", r"$G(z;0,y)=zG(z;x,\Delta)+\delta_{x,\Delta}$", r"$G(z;x,0)=zG(z;x,\Delta)+\delta_{x,0}$"),
]


def tab_errata():
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("FourFigMirror", fig_mirror),
    ("FourTabSurvive", tab_survive),
    ("FourFigDecay", fig_decay),
    ("FourTabRates", tab_rates),
    ("FourTabShapes", tab_shapes),
    ("FourFigTent", fig_tent),
    ("FourTabGreen", tab_green),
    ("FourTabBP", tab_bp),
    ("FourTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch4()
    parts = ["% figs.tex —— 由 ch4/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch4/make_figs.py（第 4 章，動畫 W401–W408）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
