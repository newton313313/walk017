r"""
make_figs.py —— 隨機漫步導讀第 6 章（Rudnick & Gaspari《Elements of the Random Walk》第 6 章，動畫 W601–W607）的圖與表

在 notes/ 或 notes/ch6/ 裡執行：
    python ch6/make_figs.py
會先檢查 ch6/part6.tex 寫到的每一個數字（check_numbers_ch6），再寫出 ch6/figs.tex（part6.tex 第一行 \input{ch6/figs}）。
圖都用 TikZ 畫，座標與直方圖計數由這裡算好（不用 pgfplots、不在 TeX 裡算）。數學一律由 codes/utils/walk_shape.py 算出；
例子與動畫 W601–W607 相同。巨集名稱都以 \Six 開頭（不和前幾章的 \Fig…、\Two…、\Three…、\Four…、\Five… 撞名）；
part6.tex 的 \label 用 six-…。顏色沿用導言區：cBar 精確值／k=1、cWalk 路線／k=2、cExact 更正／k=3、cGauss 近似、cWrong 筆誤。
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent            # notes/ch6
sys.path.insert(0, str(HERE.parent.parent / "codes"))

from utils import walk_shape as SH  # noqa: E402

OUT = HERE / "figs.tex"
F = Fraction
PI = math.pi
KCOL = ("cBar", "cWalk", "cExact")                 # k = 1, 2, 3
DS = (100, 20, 3)                                  # 1/d 展開的三個維度


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣）。"""
    return f"{float(x):.{d}f}"


def _fr(x):
    x = F(x)
    if x.denominator == 1:
        return str(x.numerator)
    sgn = "-" if x < 0 else ""
    return rf"{sgn}\frac{{{abs(x.numerator)}}}{{{x.denominator}}}"


def _close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) < tol


def _angle(vec):
    """二維主軸的角度（度，0–180）。"""
    a = math.degrees(math.atan2(vec[1], vec[0])) % 180.0
    return a


def _hist(vals, lo, hi, nb):
    """密度歸一化（除以全部樣本數與格寬）的直方圖。"""
    cnt, edges = np.histogram(vals, bins=nb, range=(lo, hi))
    w = (hi - lo) / nb
    return cnt / (len(vals) * w), edges


# =============================================================================
# 數字檢查（part6.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch6():
    assert (SH.SEED_CLOUD, SH.SEED_FIG63, SH.SEED_FIG64, SH.SEED_FIG65, SH.SEED_GAUSS) == (601, 603, 643, 605, 606)
    assert SH.N_D == 50 and SH.M_D == {100: 2000, 20: 2000, 3: 20000}

    # ---- W601 小例子「東、東、北、西」
    ex = SH.small_example()
    assert list(SH.SMALL_STEPS) == [(1, 0), (1, 0), (0, 1), (-1, 0)]
    assert ex["pts"] == [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1)] and len(ex["pts"]) == 5
    assert ex["com"] == [F(6, 5), F(2, 5)]
    assert ex["d2"] == [F(8, 5), F(1, 5), F(4, 5), F(1), F(2, 5)] and sum(ex["d2"]) == F(20, 5) == 4
    assert ex["Rg2"] == F(4, 5) == sum(ex["d2"]) / 5
    # 一大群：1000 條 100 步
    assert SH.cloud().shape == (1000, 101, 2)
    cs = SH.cloud_stats()
    Tm = cs["Tmean"]
    assert (_s(Tm[0, 0], 2), _s(Tm[1, 1], 2), _s(Tm[0, 1], 2)) == ("8.35", "8.59", "0.03")
    assert cs["exact_T11"] == F(100 * 102, 12 * 101) == SH.trace_a(100) / 2 and _s(cs["exact_T11"], 2) == "8.42"
    assert tuple(_s(v, 2) for v in cs["eig_of_mean"]) == ("8.60", "8.34")
    assert (_s(cs["R1"], 2), _s(cs["R2"], 2), _s(cs["ratio"], 2)) == ("14.14", "2.80", "5.05")
    f63 = SH.fig63()
    assert len(f63) == 6 and all(len(p) == 1001 for p in f63)
    assert tuple(_s(r, 1) for r in SH.fig63_ratios()) == ("2.1", "3.0", "6.7", "5.5", "14.5", "4.6")

    # ---- W602 迴轉半徑張量
    T = ex["T"]
    assert T == [[F(14, 25), F(3, 25)], [F(3, 25), F(6, 25)]] == SH.tensor_frac(ex["pts"])
    xs = [F(p[0]) - ex["com"][0] for p in ex["pts"]]
    ys = [F(p[1]) - ex["com"][1] for p in ex["pts"]]
    assert xs == [F(-6, 5), F(-1, 5), F(4, 5), F(4, 5), F(-1, 5)] and ys == [F(-2, 5)] * 3 + [F(3, 5)] * 2
    assert [x * x * 25 for x in xs] == [36, 1, 16, 16, 1] and [y * y * 25 for y in ys] == [4, 4, 4, 9, 9]
    assert [x * y * 25 for x, y in zip(xs, ys)] == [12, 2, -8, 12, -3]
    assert F(36 + 1 + 16 + 16 + 1, 5 * 25) == T[0][0] and F(4 + 4 + 4 + 9 + 9, 5 * 25) == T[1][1]
    assert F(12 + 2 - 8 + 12 - 3, 5 * 25) == T[0][1]
    assert math.isqrt(8 ** 2 + 4 * 3 ** 2) ** 2 == 100 and (20 + 10) // 2 == 15 and (20 - 10) // 2 == 5
    assert 14 + 6 == 20 and 14 - 6 == 8 and (14 - 15) * 3 + 3 * 1 == 0          # (3,1) 是 15 的特徵向量
    assert ex["eig"] == (F(15, 25), F(5, 25)) == (F(3, 5), F(1, 5)) and ex["vec"] == (3, 1)
    assert _s(math.degrees(math.atan(1 / 3)), 2) == "18.43" and _s(math.degrees(math.atan(1 / 3)) + 90, 2) == "108.43"
    trT, trT2 = T[0][0] + T[1][1], T[0][0] ** 2 + 2 * T[0][1] ** 2 + T[1][1] ** 2
    detT = T[0][0] * T[1][1] - T[0][1] ** 2
    assert trT == F(4, 5) == ex["Rg2"] and trT2 == F(2, 5) == F(9 + 1, 25) and detT == F(3, 25) == F(3, 5) * F(1, 5)
    assert (trT ** 2 - trT2) / 2 == detT and trT ** 2 == F(16, 25) and trT2 == F(10, 25)
    inv = SH.invariants(np.array(T, dtype=float))
    assert _close(inv[0], 0.8) and _close(inv[1], 0.4) and _close(inv[3], 0.12)
    rng = np.random.default_rng(6)
    A3 = rng.normal(size=(3, 3))
    T3 = A3 @ A3.T
    assert abs(SH.det_from_traces(T3) - np.linalg.det(T3)) < 1e-9
    phis = np.radians(np.linspace(0, 180, 18001))
    ext = SH.extent(np.array(T, dtype=float), phis)
    assert _s(np.degrees(phis[np.argmax(ext)]), 2) == "18.43" and _close(ext.max(), 0.6, 1e-8)
    assert _s(np.degrees(phis[np.argmin(ext)]), 2) == "108.43" and _close(ext.min(), 0.2, 1e-8)
    assert (_s(math.sqrt(0.6), 4), _s(math.sqrt(0.2), 4)) == ("0.7746", "0.4472")
    # 主例 Fig. 6.4（種子 643）
    p64 = SH.fig64()
    assert p64.shape == (1001, 2) and p64[:, 0].min() == -49 and p64[:, 0].max() == 4
    assert p64[:, 1].min() == -24 and p64[:, 1].max() == 7
    st = SH.fig64_stats()
    assert tuple(_s(v, 1) for v in st["com"]) == ("-17.9", "-8.9")
    assert (_s(st["T"][0, 0], 1), _s(st["T"][1, 1], 1), _s(st["T"][0, 1], 1)) == ("207.5", "82.9", "95.0")
    assert tuple(_s(v, 1) for v in st["eig"]) == ("258.8", "31.6") and tuple(_s(v, 2) for v in st["R"]) == ("16.09", "5.62")
    assert _s(st["ratio"], 2) == "8.18" and _s(st["angle"], 1) == "28.4" and _s(st["asph"], 4) == "0.6117"
    assert (_s(st["Trot"][0, 0], 1), _s(st["Trot"][1, 1], 1), _s(st["Trot"][0, 1], 1)) == ("94.1", "196.3", "101.4")
    assert np.allclose(st["eig_rot"], st["eig"]) and np.allclose(st["Tkr"], st["T"])
    vals, vecs = SH.principal(st["Trot"])
    assert _s(_angle(vecs[:, 0]), 1) == "58.4" and _close(SH.asph(vals), st["asph"])
    assert np.allclose(SH.tensor(SH.rotate(p64, PI / 6)), st["Trot"])
    assert np.allclose(SH.tensor_pairs(p64), st["T"])

    # ---- W603 非球度
    A = SH.asph
    assert A([F(1), F(1), F(1)]) == 0 and A([F(1), F(0), F(0)]) == 1 == F(1 + 1 + 0, 2 * 1)
    assert A([F(4), F(1), F(1)]) == F(1, 4) == F(9 + 9 + 0, 2 * 36) and A([F(1), F(1), F(0)]) == F(1, 4) == F(0 + 1 + 1, 2 * 4)
    assert A([F(3, 5), F(1, 5)]) == F(1, 4) == F(2, 5) ** 2 / F(4, 5) ** 2 == ex["asph"]
    for lam in ([F(4), F(1), F(1)], [F(7), F(2), F(1, 3)], [F(3, 5), F(1, 5)]):
        d = len(lam)
        lhs = sum((lam[i] - lam[j]) ** 2 for i in range(d) for j in range(i + 1, d))
        assert lhs == d * sum(v * v for v in lam) - sum(lam) ** 2
    Ad = {d: SH.A_d(d) for d in (2, 3, 4, 10)}
    assert Ad == {2: F(4, 7), 3: F(10, 19), 4: F(1, 2), 10: F(4, 9)}
    assert [_s(Ad[d], 4) for d in (2, 3, 4, 10)] == ["0.5714", "0.5263", "0.5000", "0.4444"]
    assert abs(float(SH.A_d(10 ** 7)) - 0.4) < 1e-6 and F(4 + 2 * 3, 4 + 5 * 3) == F(10, 19)
    s3 = SH.sim3d_stats()
    ev3 = SH.sim3d()
    assert ev3.shape == (20000, 3)
    assert _s(s3["A_ratio"], 4) == "0.5231" and _s(s3["A_mean"], 4) == "0.3939"
    assert _s(SH.A_exact(100, 3, "lattice"), 4) == "0.5235" and _s(F(10, 19), 4) == "0.5263"
    assert tuple(int(c) for c in s3["hist"]) == (119, 616, 1183, 1631, 1890, 1795, 1808, 1621, 1669, 1529, 1505,
                                                   1305, 1065, 1000, 655, 389, 162, 56, 2, 0)
    assert sum(int(c) for c in s3["hist"]) == 20000 and _s(100 * s3["frac_lt_half"], 1) == "69.3"
    a_ind = SH.asph(ev3)
    tr = ev3.sum(axis=1)
    assert _close(np.sum(a_ind * tr ** 2) / np.sum(tr ** 2), s3["A_ratio"], 1e-9)   # A_d＝(Tr T)² 加權的平均
    ap = {d: SH.Aprime(d) for d in (2, 3, 4)}
    assert [_s(ap[d], 6) for d in (2, 3, 4)] == ["0.396400", "0.394274", "0.393829"] and _s(ap[3], 4) == "0.3943"
    assert _close((5 - 3.5 * SH.ZETA3) / 2, ap[2], 1e-6) and _close(4 - 3 * SH.ZETA3, ap[4], 1e-6)
    assert _close(SH.Aprime_laplace(3), ap[3], 1e-6)
    small = [F(3, 5), F(1, 5)]
    assert _close(SH.asph_book16(small, 4), 0.2) and _close(SH.asph_book16([4 * v for v in small], 4), 0.8)
    assert A([4 * v for v in small]) == A(small)

    # ---- W604 d ≫ N：樓梯
    assert _s(SH.orth_fraction(10, 100), 4) == "0.6282" and _s(SH.orth_approx(10, 100), 4) == "0.6065"
    d, N = 100, 10
    assert _close(math.prod(2 * d - 2 * k for k in range(N)) / (2 * d) ** N, SH.orth_fraction(N, d))
    assert _close(math.prod(1 - k / d for k in range(N)), SH.orth_fraction(N, d)) and _s(math.exp(-0.5), 4) == "0.6065"
    sc = SH.staircase(3)
    assert sc.tolist() == [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]
    assert SH.staircase_com(3) == [F(3, 4), F(1, 2), F(1, 4)] == [F(4 - n, 4) for n in (1, 2, 3)]
    Ts = SH.tensor_frac(sc)
    assert [[16 * v for v in row] for row in Ts] == [[3, 2, 1], [2, 4, 2], [1, 2, 3]] and Ts == SH.kramers(3)
    dev = [F(v) - F(3, 4) for v in (0, 1, 1, 1)]
    assert dev == [F(-3, 4), F(1, 4), F(1, 4), F(1, 4)] and sum(v * v for v in dev) / 4 == F(1, 4) * (F(9, 16) + F(3, 16)) == F(3, 16)
    e3 = SH.eigs_exact(3)
    assert tuple(_s(v, 4) for v in e3) == ("0.4268", "0.1250", "0.0732")
    assert np.allclose(e3, [(2 + math.sqrt(2)) / 8, 1 / 8, (2 - math.sqrt(2)) / 8])
    assert np.allclose(e3, sorted(np.linalg.eigvalsh(np.array(Ts, dtype=float)), reverse=True))
    for N in (3, 4, 10):
        ex_l = [1 / (4 * (N + 1) * math.sin(k * PI / (2 * (N + 1))) ** 2) for k in range(1, N + 1)]
        assert np.allclose(SH.eigs_exact(N), ex_l)
    assert tuple(_s(SH.eigs_asym(3, k), 4) for k in (1, 2, 3)) == ("0.4053", "0.1013", "0.0450")
    assert round(100 * (e3[0] - SH.eigs_asym(3, 1)) / e3[0]) == 5
    Tsn = np.array(Ts, dtype=float)
    for k, v in ((1, (1, math.sqrt(2), 1)), (2, (1, 0, -1)), (3, (1, -math.sqrt(2), 1))):
        v = np.array(v) / np.linalg.norm(v)
        assert np.allclose(Tsn @ v, e3[k - 1] * v)
        sm = SH.sine_mode(3, k)
        assert np.allclose(np.abs(sm @ v), 1.0)
        w = np.array([math.sin(k * PI * n / 4) for n in (1, 2, 3)])
        assert np.allclose(np.abs(w / np.linalg.norm(w) @ v), 1.0)
    assert _s(SH.zeta2, 4) == "1.6449" and _s(SH.zeta2 ** 2, 4) == "2.7058" and _s(SH.zeta4, 4) == "1.0823"
    assert _close(SH.zeta2, PI ** 2 / 6, 1e-9) and _close(SH.zeta4, PI ** 4 / 90, 1e-9) and _close(SH.zeta4 / SH.zeta2 ** 2, 0.4)
    sa = {N: SH.staircase_asph(N) for N in (3, 10, 100)}
    assert sa[3] == F(13, 25) and [_s(sa[N], 4) for N in (3, 10, 100)] == ["0.5200", "0.4150", "0.4002"]
    tra3 = sum(Ts[i][i] for i in range(3))
    assert tra3 == F(5, 8) and _close(sum(v * v for v in e3) / float(tra3) ** 2, 13 / 25)
    assert F(1, 1) / F(1, 9) == 9 and F(1, 4) / F(1, 9) == F(9, 4) == 2.25
    assert tuple(_s(v, 2) for v in s3["ratios"]) == ("11.95", "2.69")
    assert tuple(_s(v, 2) for v in s3["ratio_of_ratios"]) == ("13.69", "2.98")
    assert _s(1.01 * 9 / 0.75, 2) == "12.12" and _s(0.91 * 2.25 / 0.75, 2) == "2.73"
    for N in (3, 6):                                   # 練習 6.3(b)(c) 的編號：j = 1..N+1 會得 (N+2−n)/(N+1)
        for n in range(1, N + 1):
            assert F(sum(1 for j in range(0, N + 1) if j >= n), N + 1) == F(N + 1 - n, N + 1)
            assert F(sum(1 for j in range(1, N + 2) if j >= n), N + 1) == F(N + 2 - n, N + 1)

    # ---- W605 Kramers 恆等式、高斯步伐
    K4 = SH.kramers(4)
    M4 = [[4, 3, 2, 1], [3, 6, 4, 2], [2, 4, 6, 3], [1, 2, 3, 4]]
    assert [[25 * v for v in row] for row in K4] == M4 and SH.kramers_count(4) == M4
    assert all(M4[m - 1][n - 1] == m * (5 - n) for m in range(1, 5) for n in range(m, 5))
    assert sum(K4[i][i] for i in range(4)) == F(20, 25) == F(4, 5) == SH.trace_a(4)
    assert SH.tensor_kramers_frac(SH.SMALL_STEPS) == T
    vx, vy = (1, 1, 0, -1), (0, 0, 1, 0)
    Av = [sum(M4[i][j] * vx[j] for j in range(4)) for i in range(4)]
    assert Av == [6, 7, 3, -1] and sum(a * b for a, b in zip(Av, vx)) == 6 + 7 + 0 + 1 == 14
    assert sum(a * b for a, b in zip(Av, vy)) == 3 and M4[2][2] == 6
    nd = np.array(SH.kramers_book(4, "nodiag"), dtype=float)
    assert abs(np.array(vx) @ nd @ np.array(vx)) < 1e-12 and abs(np.array(vy) @ nd @ np.array(vy)) < 1e-12
    assert np.allclose(np.array(SH.kramers_book(4, "noscale"), dtype=float), 25 * np.array(K4, dtype=float))
    for N in (4, 9):
        K = np.array(SH.kramers(N), dtype=float)
        x = np.arange(1, N + 1) / (N + 1)
        assert np.allclose(K, np.minimum.outer(x, x) * (1 - np.maximum.outer(x, x)))     # x_<(1−x_>)
        ev = SH.eigs_exact(N)
        for k in range(1, N + 1):
            sm = SH.sine_mode(N, k)
            assert np.allclose(K @ sm, ev[k - 1] * sm)
        assert np.allclose(SH.kramers_np(N), K)
    e4 = SH.eigs_exact(4)
    assert tuple(_s(v, 4) for v in e4) == ("0.5236", "0.1447", "0.0764", "0.0553")
    assert tuple(_s(SH.eigs_asym(4, k), 4) for k in (1, 2, 3, 4)) == ("0.5066", "0.1267", "0.0563", "0.0317")
    assert _close(sum(e4), 0.8)
    assert len(SH.wick_pairings(2)) == 3 and _close(SH.gauss_moment(2), 3, 1e-8)
    for N in (4, 100):
        for dd in (2, 3):
            ta, ta2 = SH.trace_a(N), SH.trace_a2(N)
            assert SH.moments_exact(N, dd, "gauss") == ((ta ** 2 + 2 * ta2) / dd ** 2, ta ** 2 / dd ** 2, ta2 / dd ** 2)
    assert SH.moments_enum(3, 2) == SH.moments_exact(3, 2, "lattice")
    Nb = 10 ** 6
    assert abs(float(SH.trace_a(Nb)) / (Nb / 6) - 1) < 1e-5 and abs(float(SH.trace_a2(Nb)) / (Nb ** 2 / 90) - 1) < 1e-5
    assert _close(PI ** 4 / 36 + 2 * PI ** 4 / 90, PI ** 4 / 20) and F(1, 36) + F(2, 90) == F(1, 20)
    assert F(1, 20) / F(1, 36) == F(9, 5) and F(1, 20) / F(1, 90) == F(9, 2)
    for dd in range(1, 12):
        assert F(9 - 5 + 2 * dd, 9 + 5 * (dd - 1)) == SH.A_d(dd) == SH.A_from_moments(dd, F(1, 20), F(1, 36), F(1, 90))
    assert F(1, 36) + F(1, 90) == F(7, 180)                                   # 照書上 (6.35b)
    assert _s(SH.A_exact(4, 3), 4) == "0.6013" and _s(SH.A_exact(100, 3), 4) == "0.5265"
    sgm = SH.sim_gauss_moments()[:3]
    assert tuple(_s(v, 2) for v in sgm) == ("56.09", "31.31", "12.56")
    assert tuple(_s(v, 2) for v in SH.moments_exact(100, 3)) == ("56.67", "31.48", "12.60")
    assert _s(SH.A_from_moments(3, *sgm), 4) == "0.5261"
    assert 2 + 16 + 18 == 36 and F(36, 720) == F(1, 20)
    assert SH.A_S21(3) == F(10, 19) and SH.A_S21(3, True) == F(25, 19) == F(5, 2) * F(10, 19)
    assert F(9, 2) - F(5, 2) == 2 and F(90, 36) == F(5, 2)

    # ---- W606 1/d 展開
    assert [SH.wick_count(m) for m in range(1, 6)] == [1, 3, 15, 105, 945]
    assert all(SH.wick_count(m) == math.factorial(2 * m) // (2 ** m * math.factorial(m)) for m in range(1, 6))
    assert all(len(SH.wick_pairings(m)) == SH.wick_count(m) for m in range(1, 5))
    assert all(_close(SH.gauss_moment(m), SH.wick_count(m), 1e-6 * SH.wick_count(m)) for m in range(1, 6))
    assert _s(SH.resolvent_weight(), 4) == "3.1416"
    lead = [1 / (PI ** 2 * k ** 2) for k in (1, 2, 3)]
    assert [_s(v, 4) for v in lead] == ["0.1013", "0.0253", "0.0113"]
    assert all(abs(SH.sum_identity(n) - 0.75) < 1e-6 for n in (1, 2, 3)) and F(1, 2) * (1 + F(1, 2)) == F(3, 4)
    assert [_s(SH.shift_1d(d), 4) for d in DS] == ["0.0075", "0.0375", "0.2500"]
    assert [_s(math.sqrt(2 / d), 4) for d in DS] == ["0.1414", "0.3162", "0.8165"]
    od = {d: SH.onedee_table(d) for d in DS}
    for d in DS:
        e, trd = SH.eig_sim(d)
        assert e.shape == (SH.M_D[d], 3) and trd.shape == (SH.M_D[d],)
        assert np.allclose(od[d]["lead"], lead) and np.allclose(od[d]["pred_mean"], [SH.mean_1d(k, d, 50) / 51 for k in (1, 2, 3)])
        assert np.allclose(od[d]["pred_var"], [SH.var_1d(k, d, 50) / 51 ** 2 for k in (1, 2, 3)])
    want = {100: (("1.01", "1.00", "1.00"), ("0.98", "0.97", "0.94")),
            20: (("1.01", "1.01", "1.00"), ("1.02", "0.87", "0.68")),
            3: (("1.01", "0.90", "0.75"), ("1.01", "0.47", "0.23"))}
    for d in DS:
        assert tuple(_s(v, 2) for v in od[d]["mean_ratio"]) == want[d][0]
        assert tuple(_s(v, 2) for v in od[d]["var_ratio"]) == want[d][1]
    book = {(100, "var"): (0.97, 0.93, 0.95), (3, "var"): (0.98, 0.48, 0.23), (3, "mean"): (1.01, 0.91, 0.75)}
    for (d, kind), vals in book.items():
        ours = od[d]["var_ratio"] if kind == "var" else od[d]["mean_ratio"]
        assert all(abs(float(_s(a, 2)) - b) <= 0.05 + 1e-12 for a, b in zip(ours, vals))

    # ---- W607 主半徑的分布、Fixman、自迴避
    al = 1 / PI ** 2
    assert _s(al, 4) == "0.1013" and _s(al / 3, 4) == "0.0338" and _s(SH.gamma_density(al / 3, al, 3), 2) == "7.16"
    xg = np.linspace(1e-6, 2.0, 200001)
    gd = SH.gamma_density(xg, al, 3)
    assert abs(xg[np.argmax(gd)] - al / 3) < 1e-4
    dx = xg[1] - xg[0]
    assert abs(np.sum(xg * gd) * dx - al) < 1e-4 and abs(np.sum((xg - al) ** 2 * gd) * dx - 2 * al ** 2 / 3) < 1e-4
    e3sim, tr3 = SH.eig_sim(3)
    assert _s(e3sim[:, 0].mean(), 4) == "0.1280" and _s(od[3]["pred_mean"][0], 4) == "0.1267"
    fm = SH.fixman_mean_var(3)
    assert _close(fm[0], 1 / 6) and _close(fm[1], 1 / 135) and F(1, 45 * 3) == F(1, 135)
    assert (_s(fm[0], 4), _s(fm[1], 4), _s(math.sqrt(fm[1]), 4)) == ("0.1667", "0.0074", "0.0861")
    assert _close(SH.zeta2 / PI ** 2, 1 / 6, 1e-9) and _close(2 / 3 * SH.zeta4 / PI ** 4, 1 / 135, 1e-9)
    xf = np.linspace(0.08, 0.15, 71)
    yf = SH.fixman_density(xf, 3)
    assert _s(xf[np.argmax(yf)], 3) == "0.111" and _s(yf.max(), 2) == "6.20"
    assert tuple(_s(v, 3) for v in SH.fixman_density(np.array([0.1, 0.2, 0.3]), 3)) == ("6.078", "3.281", "1.011")
    assert tuple(_s(v, 4) for v in SH.fixman_moments_numeric(3)) == ("1.0000", "0.1667", "0.0074")
    assert (_s(tr3.mean(), 4), _s(tr3.var(), 4)) == ("0.1671", "0.0076")
    assert tuple(_s(v, 4) for v in SH.fixman_mean_var(3, 50)) == ("0.1666", "0.0074")
    assert np.max(np.abs(SH.fixman_density(np.array([0.05, 0.1, 0.2, 0.3]), 3, book=True))) < 1e-2
    for s in (0.3, 2.0, 7.0):
        assert abs(SH.fixman_K(s, 3) - SH.fixman_K_product(s, 3, 4000)) < 1e-3
    assert _close(SH.A_saw(3), float(F(10, 19)) + 0.008) and _s(SH.A_saw(3), 4) == "0.5343"
    # 練習 6.4、6.5、6.1
    for N in (3, 4, 7):
        tri = 2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
        assert np.allclose(tri @ ((N + 1) * np.array(SH.kramers(N), dtype=float)), np.eye(N))
    for N in (3, 5):
        assert SH.closed_cov(N) == [[(1 if i == j else 0) - F(1, N) for j in range(N)] for i in range(N)]
        assert all(sum(row) == 0 for row in SH.closed_cov(N))
    xs0 = [p[0] for p in ex["pts"]]
    assert sum((xs0[i] - xs0[j]) ** 2 for i in range(5) for j in range(i + 1, 5)) == 14 and F(14, 25) == T[0][0]


# =============================================================================
# 表與圖
# =============================================================================
def _cross(x, y, r, opt="thick"):
    return rf"\draw[{opt}] ({x - r:.4f},{y - r:.4f}) -- ({x + r:.4f},{y + r:.4f}) ({x - r:.4f},{y + r:.4f}) -- ({x + r:.4f},{y - r:.4f});"


def _axes_ellipse(cx, cy, ang, R1, R2, lab2=(0.93, "left")):
    """半長 2R_i 的兩條主軸與橢圓；lab2 = 2R_2 標籤的 (pos, anchor)。"""
    u = (math.cos(math.radians(ang)), math.sin(math.radians(ang)))
    v = (-u[1], u[0])
    out = [rf"\draw[cExact,thick,rotate around={{{ang:.3f}:({cx:.4f},{cy:.4f})}}] ({cx:.4f},{cy:.4f}) "
           rf"ellipse [x radius={2 * R1:.4f}, y radius={2 * R2:.4f}];"]
    for (w, R, name, (pos, anc)) in ((u, R1, "1", (0.93, "below right")), (v, R2, "2", lab2)):
        a = (cx - 2 * R * w[0], cy - 2 * R * w[1])
        b = (cx + 2 * R * w[0], cy + 2 * R * w[1])
        out.append(rf"\draw[cBar,thick] ({a[0]:.4f},{a[1]:.4f}) -- ({b[0]:.4f},{b[1]:.4f})"
                   rf" node[pos={pos},{anc},font=\scriptsize,cBar,inner sep=1pt] {{$2R_{name}$}};")
    return out


def fig_small():
    ex = SH.small_example()
    cx, cy = (float(v) for v in ex["com"])
    R1, R2 = (math.sqrt(float(v)) for v in ex["eig"])
    ang = math.degrees(math.atan2(1, 3))
    pts = ex["pts"]
    out = [r"\begin{tikzpicture}[x=1.45cm,y=1.45cm,>=Stealth]",
           r"\draw[gray!25,step=1] (-0.4,-0.6) grid (2.75,1.45);"]
    out += _axes_ellipse(cx, cy, ang, R1, R2, lab2=(0.06, "right"))
    for (x, y) in pts:
        out.append(rf"\draw[gray,densely dashed] ({cx:.3f},{cy:.3f}) -- ({x},{y});")
    for a, b in zip(pts[:-1], pts[1:]):
        out.append(rf"\draw[->,cWalk,very thick,shorten >=2.5pt,shorten <=2.5pt] {a} -- {b};")
    pos = ("below left", "below", "below right", "above right", "left")     # (1,1) 正好在短軸上，標在左邊
    for (x, y), d2, p in zip(pts, ex["d2"], pos):
        out.append(rf"\fill ({x},{y}) circle (1.7pt) node[{p},font=\scriptsize,inner sep=2pt] {{${_fr(d2).replace(chr(92) + 'frac', chr(92) + 'tfrac')}$}};")
    out.append(_cross(cx, cy, 0.07))
    out.append(rf"\node[font=\scriptsize,anchor=north west,inner sep=1pt] at ({cx + 0.04:.3f},{cy - 0.04:.3f}) {{重心}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def _compress(p):
    """去掉同一直線上的中間點（畫起來一樣，檔案小很多）。"""
    keep = [p[0]]
    for i in range(1, len(p) - 1):
        d1 = p[i] - p[i - 1]
        d2 = p[i + 1] - p[i]
        if not np.array_equal(d1, d2):
            keep.append(p[i])
    keep.append(p[-1])
    return keep


def fig_main():
    p = SH.fig64()
    st = SH.fig64_stats()
    cx, cy = (float(v) for v in st["com"])
    R1, R2 = (float(v) for v in st["R"])
    pts = _compress(p)
    out = [r"\begin{tikzpicture}[x=0.105cm,y=0.105cm,>=Stealth]"]
    for i in range(0, len(pts), 60):
        seg = pts[i:i + 61]
        out.append(r"\draw[cWalk!85!black,line width=0.35pt,line join=round] "
                   + " -- ".join(f"({int(a)},{int(b)})" for a, b in seg) + ";")
    out += _axes_ellipse(cx, cy, float(st["angle"]), R1, R2)
    out.append(r"\draw[thick,fill=white] (0,0) circle (1.8pt);")
    out.append(_cross(cx, cy, 1.1, "very thick"))
    out.append(r"\draw[|-|,gray] (-48,-29) -- (-38,-29) node[midway,below,font=\scriptsize] {10};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_main():
    ex = SH.small_example()
    st = SH.fig64_stats()
    vals, vecs = SH.principal(st["Trot"])
    T = ex["T"]
    c = ", ".join(_s(v, 1) for v in st["com"])
    t0 = ", ".join(_s(v, 1) for v in (st["T"][0, 0], st["T"][1, 1], st["T"][0, 1]))
    t1 = ", ".join(_s(v, 1) for v in (st["Trot"][0, 0], st["Trot"][1, 1], st["Trot"][0, 1]))
    eg = ", ".join(_s(v, 1) for v in st["eig"])
    R = ", ".join(_s(v, 2) for v in st["R"])
    rows = [
        r"點數 $N+1$ & 5 & 1001 & 1001\\",
        rf"重心 & $(\frac65,\frac25)$ & $({c})$ & $({c})$\\",
        rf"$T_{{11}},\ T_{{22}},\ T_{{12}}$ & ${_fr(T[0][0])},\ {_fr(T[1][1])},\ {_fr(T[0][1])}$ & {t0} & {t1}\\",
        rf"$R_1^2,\ R_2^2$ & $\frac35,\ \frac15$ & {eg} & {eg}\\",
        rf"$R_1,\ R_2$ & {_s(math.sqrt(0.6), 4)}, {_s(math.sqrt(0.2), 4)} & {R} & {R}\\",
        rf"主軸角 & ${_s(math.degrees(math.atan(1 / 3)), 2)}^\circ$ & ${_s(st['angle'], 1)}^\circ$ & ${_s(_angle(vecs[:, 0]), 1)}^\circ$\\",
        rf"非球度 $\big(\frac{{R_1^2-R_2^2}}{{R_1^2+R_2^2}}\big)^2$ & $\frac14$ & {_s(st['asph'], 4)} & {_s(SH.asph(vals), 4)}\\",
    ]
    return "\n".join([r"\begin{tabular}{@{}lccc@{}}", r"\toprule",
                      r" & 東東北西 & 主例（種子 643） & 主例轉 $30^\circ$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_A():
    ds = (2, 3, 4, 10)
    r1 = r"$A_d=\frac{4+2d}{4+5d}$ & " + " & ".join(f"${_fr(SH.A_d(d))}$" for d in ds) + r" & $\frac25$\\"
    r2 = r"\quad 小數 & " + " & ".join(_s(SH.A_d(d), 4) for d in ds) + r" & 0.4000\\"
    r3 = r"$A'_d$（個別值的平均） & " + " & ".join(_s(SH.Aprime(d), 6) for d in (2, 3, 4)) + r" & --- & $\to\frac25$\\"
    return "\n".join([r"\begin{tabular}{@{}lccccc@{}}", r"\toprule",
                      r"$d$ & 2 & 3 & 4 & 10 & $\infty$\\", r"\midrule", r1, r2, r3,
                      r"\bottomrule", r"\end{tabular}"])


def fig_hist():
    s3 = SH.sim3d_stats()
    h = [int(c) for c in s3["hist"]]
    out = [r"\begin{tikzpicture}[x=9cm,y=0.0016cm,>=Stealth]"]
    for i, c in enumerate(h):
        if c:
            out.append(rf"\draw[cBar,fill=cBar!25] ({0.05 * i:.2f},0) rectangle ({0.05 * (i + 1):.2f},{c});")
    out += [r"\draw[->] (0,0) -- (1.06,0) node[right,font=\scriptsize] {$\hat A$};",
            r"\draw[->] (0,0) -- (0,2150) node[above,font=\scriptsize] {條數};",
            r"\foreach \x in {0,0.2,0.4,0.6,0.8,1} {\draw (\x,0) -- (\x,-40); \node[below,font=\scriptsize] at (\x,-30) {$\x$};}",
            r"\foreach \y in {500,1000,1500} {\draw (0,\y) -- (-0.01,\y); \node[left,font=\scriptsize] at (-0.01,\y) {\y};}"]
    for v, sty, lab, anc in ((s3["A_ratio"], "cExact,very thick,dashed", r"平均的比 " + _s(s3["A_ratio"], 4), "south west"),
                             (s3["A_mean"], "cWalk,very thick,densely dotted", r"個別值的平均 " + _s(s3["A_mean"], 4), "south east")):
        out.append(rf"\draw[{sty}] ({v:.4f},0) -- ({v:.4f},2050) node[anchor={anc},font=\scriptsize,inner sep=1pt] {{{lab}}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_stair():
    e3, e4 = SH.eigs_exact(3), SH.eigs_exact(4)
    rows = [r"樓梯 $N=3$：精確 & " + " & ".join(_s(v, 4) for v in e3) + r" & \\",
            r"\quad 大 $N$ 公式 $\frac{N+1}{\pi^2k^2}$ & " + " & ".join(_s(SH.eigs_asym(3, k), 4) for k in (1, 2, 3)) + r" & \\",
            r"Kramers $N=4$：精確 & " + " & ".join(_s(v, 4) for v in e4) + r"\\",
            r"\quad 大 $N$ 公式 & " + " & ".join(_s(SH.eigs_asym(4, k), 4) for k in (1, 2, 3, 4)) + r"\\"]
    a = "\n".join([r"\begin{tabular}{@{}lcccc@{}}", r"\toprule", r"特徵值 & $k=1$ & $k=2$ & $k=3$ & $k=4$\\", r"\midrule",
                   *rows, r"\bottomrule", r"\end{tabular}"])
    sa = " & ".join(_s(SH.staircase_asph(N), 4) for N in (3, 10, 100))
    b = "\n".join([r"\begin{tabular}[t]{@{}lcccc@{}}", r"\toprule", r"樓梯 & $N=3$ & 10 & 100 & $\infty$\\", r"\midrule",
                   r"$\frac{\mathrm{Tr}\,\mathbf T^2}{(\mathrm{Tr}\,\mathbf T)^2}$ & " + sa + r" & 0.4000\\",
                   r"\bottomrule", r"\end{tabular}"])
    s3 = SH.sim3d_stats()
    r, rr = s3["ratios"], s3["ratio_of_ratios"]
    c = "\n".join([r"\begin{tabular}[t]{@{}ll@{}}", r"\toprule", r"三維的比 & $R_1^2:R_2^2:R_3^2$\\", r"\midrule",
                   r"(6.53)，$d\to\infty$ & $9:2.25:1$\\",
                   rf"模擬 $\langle R_1^2\rangle:\langle R_2^2\rangle:\langle R_3^2\rangle$ & ${_s(r[0], 2)}:{_s(r[1], 2)}:1$\\",
                   rf"模擬 $\langle R_1^2/R_3^2\rangle:\langle R_2^2/R_3^2\rangle$ & ${_s(rr[0], 2)}:{_s(rr[1], 2)}:1$\\",
                   r"書上 (6.54) & {\color{cWrong}$13.76:3.03:1$}\\",
                   r"\bottomrule", r"\end{tabular}"])
    return a + "\n\\\\[6pt]\n" + b + "\n\\qquad\n" + c


def fig_oned():
    W, H, gap = 4.75, 2.3, 0.75
    rng = {100: (0.5, 1.6, (0.6, 0.8, 1.0, 1.2, 1.4)), 20: (0.0, 2.5, (0, 0.5, 1, 1.5, 2, 2.5)),
           3: (0.0, 4.0, (0, 1, 2, 3, 4))}
    lead = np.array([1 / (PI ** 2 * k ** 2) for k in (1, 2, 3)])
    out = [r"\begin{tikzpicture}[>=Stealth]"]
    for j, d in enumerate(DS):
        x0 = j * (W + gap)
        lo, hi, ticks = rng[d]
        e = SH.eig_sim(d)[0] / lead
        hs = [_hist(e[:, k], lo, hi, 25) for k in range(3)]
        ymax = max(hh.max() for hh, _ in hs) * 1.08

        def X(v):
            return x0 + (v - lo) / (hi - lo) * W

        def Y(v):
            return v / ymax * H

        for k, (hh, edges) in enumerate(hs):
            path = [f"({X(edges[0]):.3f},0)"]
            for i, c in enumerate(hh):
                path += [f"({X(edges[i]):.3f},{Y(c):.3f})", f"({X(edges[i + 1]):.3f},{Y(c):.3f})"]
            path.append(f"({X(edges[-1]):.3f},0)")
            out.append(rf"\draw[{KCOL[k]},thick,fill={KCOL[k]},fill opacity=0.12] " + " -- ".join(path) + ";")
        out.append(rf"\draw[->] ({x0:.3f},0) -- ({x0 + W + 0.15:.3f},0);")
        out.append(rf"\draw ({x0:.3f},0) -- ({x0:.3f},{H:.3f});")
        for t in ticks:
            out.append(rf"\draw ({X(t):.3f},0) -- ({X(t):.3f},-0.08) node[below,font=\scriptsize] {{${t}$}};")
        out.append(rf"\draw[gray,thick,dashed] ({X(1):.3f},0) -- ({X(1):.3f},{H:.3f});")
        s1 = 1 + 0.75 / d
        out.append(rf"\draw[black,densely dotted,thick] ({X(s1):.3f},0) -- ({X(s1):.3f},{H:.3f});")
        out.append(rf"\node[anchor=north east,font=\scriptsize,align=right] at ({x0 + W:.3f},{H:.3f}) "
                   rf"{{$d={d}$\\{SH.M_D[d]} 條}};")
    out.append(rf"\node[font=\scriptsize] at ({1.5 * W + gap:.3f},-0.75) {{$\Lambda_k\big/\frac{{N+1}}{{\pi^2k^2}}$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_oned():
    rows = []
    for d in DS:
        t = SH.onedee_table(d)
        f4 = lambda arr: ", ".join(_s(v, 4) for v in arr)   # noqa: E731
        f2 = lambda arr: ", ".join(_s(v, 2) for v in arr)   # noqa: E731
        rows.append(rf"{d}（{SH.M_D[d]}） & {f4(t['mean'])} & {f4(t['pred_mean'])} & {f4(t['std'])} & "
                    rf"{f2(t['mean_ratio'])} & {f2(t['var_ratio'])}\\")
    return "\n".join([r"{\setlength{\tabcolsep}{3.5pt}\begin{tabular}{@{}lccccc@{}}", r"\toprule",
                      r"$d$（條數） & 平均，$k=1,2,3$ & (6.61) & 標準差 & 平均$\div$(6.61) & 變異數$\div$(6.62)\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}}"])


def _panel(x0, W, H, ymax, hist, curve, ccol, vlines, ylab):
    hh, edges = hist
    xs, ys = curve

    def X(v):
        return x0 + v / 0.6 * W

    def Y(v):
        return v / ymax * H

    out = []
    for i, c in enumerate(hh):
        if c > 0:
            out.append(rf"\draw[cBar!70,fill=cBar!20] ({X(edges[i]):.3f},0) rectangle ({X(edges[i + 1]):.3f},{Y(c):.3f});")
    out.append(rf"\draw[{ccol},very thick] " + " -- ".join(f"({X(a):.3f},{Y(b):.3f})" for a, b in zip(xs, ys)) + ";")
    for v, sty, lab in vlines:
        out.append(rf"\draw[{sty}] ({X(v):.3f},0) -- ({X(v):.3f},{H:.3f}) node[above,font=\scriptsize,inner sep=1pt] {{{lab}}};")
    out.append(rf"\draw[->] ({x0:.3f},0) -- ({x0 + W + 0.2:.3f},0);")
    out.append(rf"\draw[->] ({x0:.3f},0) -- ({x0:.3f},{H + 0.25:.3f}) node[above,font=\scriptsize] {{{ylab}}};")
    for t in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
        out.append(rf"\draw ({X(t):.3f},0) -- ({X(t):.3f},-0.07) node[below,font=\scriptsize] {{${t:g}$}};")
    for t in (2, 4, 6):
        out.append(rf"\draw ({x0:.3f},{Y(t):.3f}) -- ({x0 - 0.07:.3f},{Y(t):.3f}) node[left,font=\scriptsize] {{${t}$}};")
    return out


def fig_fix():
    W, H, gap = 6.6, 2.8, 1.3
    e3, tr3 = SH.eig_sim(3)
    al = 1 / PI ** 2
    xs = np.concatenate([np.linspace(0.0, 0.06, 30), np.linspace(0.062, 0.6, 90)])
    g = SH.gamma_density(xs, al, 3)
    xf = np.concatenate([np.linspace(0.003, 0.06, 30), np.linspace(0.062, 0.6, 90)])
    fx = SH.fixman_density(xf, 3)
    assert len(xs) <= 120 and len(xf) <= 120
    out = [r"\begin{tikzpicture}[>=Stealth]"]
    out += _panel(0.0, W, H, 8.0, _hist(e3[:, 0], 0, 0.6, 30), (xs, g), "cGauss",
                  ((float(e3[:, 0].mean()), "cWalk,thick,dashed", "模擬平均 " + _s(e3[:, 0].mean(), 4)),),
                  r"密度")
    out.append(r"\node[font=\scriptsize,anchor=west] at (2.4,1.9) {\color{cGauss}Gamma (6.73)：峰 0.0338、高 7.16};")
    out.append(r"\node[font=\scriptsize,anchor=west] at (2.4,1.45) {\color{cBar!80!black}模擬 $R_1^2/(N+1)$};")
    x1 = W + gap
    out += _panel(x1, W, H, 8.0, _hist(tr3, 0, 0.6, 30), (xf, fx), "cExact",
                  ((1 / 6, "cWalk,thick,dashed", r"平均 $\frac16$"),), r"密度")
    out.append(rf"\node[font=\scriptsize,anchor=west] at ({x1 + 2.6:.3f},1.9) {{\color{{cExact}}Fixman (6.74)：峰 0.111、高 6.20}};")
    out.append(rf"\node[font=\scriptsize,anchor=west] at ({x1 + 2.6:.3f},1.45) {{\color{{cBar!80!black}}模擬 $\mathrm{{Tr}}\,\mathbf T/(N+1)$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


ERRATA = [
    (r"p.~135，(6.16)", r"$(d-1)\sum_{k=1}^{N}R_k^2$", r"$(d-1)\big(\sum_{k=1}^{d}R_k^2\big)^2$"),
    (r"p.~138，練習 6.3(b)(c)", r"(b) $j<n$／$j\ge n$ 配 (c) $\sum_{j=1}^{N+1}$", r"點編號 $j=0..N$：$\sum_{j=0}^{N}$"),
    (r"p.~139，(6.22)", r"$\alpha<\beta$；$\alpha>\beta$", r"$\alpha\le\beta$；$\alpha>\beta$（含對角線）"),
    (r"p.~139，(6.23)", r"$a_{x_1,x_2}=x_<(N+1-x_>)$", r"$a_{x_1,x_2}=\frac{x_<(N+1-x_>)}{(N+1)^2}$"),
    (r"p.~141，(6.33a)", r"$\sum_{\alpha,\gamma}a_{\alpha,\alpha}a_{\delta,\delta}$", r"$\sum_{\alpha,\delta}a_{\alpha,\alpha}a_{\delta,\delta}$"),
    (r"p.~141，(6.35a)", r"$\frac{1}{d^2}\frac{(N+1)^2}{\pi^4d^2}\big[\big(\sum_n\frac1{n^2}\big)+\cdots\big]$",
     r"$\frac{1}{d^2}\frac{(N+1)^2}{\pi^4}\big[\big(\sum_n\frac1{n^2}\big)^2+\cdots\big]$"),
    (r"p.~141，(6.35b)", r"$\frac{\pi^4}{36}+\frac{\pi^4}{90}$", r"$\frac{\pi^4}{36}+\frac{2\pi^4}{90}$"),
    (r"p.~146，(6.52)", r"$\sum_i\delta(\lambda-\lambda_i)$", r"$\pi\sum_i\delta(\lambda-\lambda_i)$"),
    (r"p.~146，(6.54)", r"$\langle R_1^2\rangle:\langle R_2^2\rangle:\langle R_3^2\rangle=13.76:3.03:1$",
     r"$\approx11.9:2.7:1$（書上的數接近比值的平均）"),
    (r"p.~153，(6.76)", r"$K(R)=\prod\big(1-\frac{2is\alpha_n}{d}\big)^{-d/2}$", r"$K(s)=\prod\big(1+\frac{2is\alpha_n}{d}\big)^{-d/2}$"),
    (r"p.~159，(6:S2-1)", r"$\cdots\big/\big(\frac{\langle T_{11}^2\rangle}{\langle T_{11}T_{22}\rangle}+d-1\big)$",
     r"分母再乘 $\frac{\langle T_{11}T_{22}\rangle}{\langle T_{12}^2\rangle}$"),
    (r"p.~164，(6:S2-25) 前", r"Adding (6:S2-16), (6:S2-21), and (6:S2-23)", r"(6:S2-24)"),
]


def tab_errata():
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("SixFigSmall", fig_small),
    ("SixFigMain", fig_main),
    ("SixTabMain", tab_main),
    ("SixTabA", tab_A),
    ("SixFigHist", fig_hist),
    ("SixTabStair", tab_stair),
    ("SixFigOneDee", fig_oned),
    ("SixTabOneDee", tab_oned),
    ("SixFigFixman", fig_fix),
    ("SixTabErrata", tab_errata),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch6()
    parts = ["% figs.tex —— 由 ch6/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch6/make_figs.py（第 6 章，動畫 W601–W607）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
