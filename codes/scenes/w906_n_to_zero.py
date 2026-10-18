"""
W906 n → 0 的魔術：O(n) 模型就是自迴避漫步＋第 9 章回顧
（§9.5–9.6，pp. 219–231，(9.42)–(9.77)，(9:S-1)–(9:S-22)，Ex. 9.5–9.6）

問：「n = 0 個分量的磁鐵」聽起來荒謬——它為什麼會數自迴避漫步？目標 (9.45)。
  1. (9.46) 輔助場：e^{−u|S|⁴/4} = (1/√(πu))∫dφ e^{−φ²/u − iφ|S|²}（0.7788、0.3679）；(⋯)^n → Z_n → 1；
     (9:S-11) 高斯積分 5.0627；⟨x₁x₂⟩ = 1/3、⟨x₁²⟩ = 2/3；Wick 2/3、4/3；配對數 1、3、15、105。
  2. Trotter (9.66)、Ex. 9.6：一階誤差 ∝ ε²、加交換子項 ∝ ε³；切片 → 路徑積分 (9.67)–(9.71)（示意）。
  3. (9.72)–(9.74)：走過 k 次 → e^{−uk²/4}，每相交一次 e^{−u/2} = 0.6065；筆誤框 1（∫_0^∞）；
     N = 8：65536、23695.11、7988.15、5916；筆誤框 2（(9.75) 少了負號 → 6.2e5）。
  4. 2×3 梯子：G_03、G_05（n = 1、2、→ 0）＝自迴避路徑的和；(9.76)(9.77)。
  5. 第 9 章回顧（六列）；收束。
畫面上的數字都由 utils.walk_crit（c_8 由 utils.walk_saw）算出，check_math() 核對。
模組層級的畫圖函式：saw_demo、cross_demo（格子上的路）、pairing_row（Wick 弧線圖）、time_slices（路徑積分示意）、
weak_bars（長條）、ladder_graph／ladder_paths／path_mob（梯子與自迴避路徑）。
"""
import math
import sys
import warnings
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_crit as WC
from utils import walk_saw as SW
from utils import lattice_draw as L
from utils import walk_draw as D

HS_S = (1, 2)                                     # |S|² = 1、2
SAW_DEMO = ((-2, -2), (-1, -2), (-1, -1), (0, -1), (0, 0), (-1, 0), (-1, 1), (0, 1), (1, 1))      # 8 步自迴避
CROSS_DEMO = ((-2, 0), (-1, 0), (0, 0), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))    # (−1, 0) 走過兩次
PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
PAIR_TEX = (r"C_{rs}C_{tu}", r"C_{rt}C_{su}", r"C_{ru}C_{ts}")
SLICE_H = (0.15, 0.55, 0.3, 0.8, 0.45, 0.95, 0.7)    # 路徑積分示意的高度（相對值，無數字）
PATH_COLOR = YELLOW
START_COLOR = "#7CE38B"
END_COLOR = RED


def f2(x):
    return f"{float(x):.2f}"


def f4(x):
    return f"{float(x):.4f}"


def f5(x):
    return f"{float(x):.5f}"


def g4(x):
    """4 位有效數字；科學記號寫成 TeX。"""
    s = f"{float(x):.4g}"
    if "e" in s:
        m, e = s.split("e")
        return r"%s\times10^{%d}" % (m, int(e))
    return s


def sci1(x):
    m, e = f"{float(x):.1e}".split("e")
    return r"%s\times10^{%d}" % (m, int(e))


def fcount(v):
    return f"{v:.0f}" if float(v).is_integer() else f"{v:.2f}"


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        C = np.linalg.inv(np.array(WC.GAUSS_T, float))
        nb = {"hs": [WC.hs_rhs(s) for s in HS_S],
              "gI": WC.gauss_I(), "gIq": WC.gauss_I_quad(),
              "mom": [WC.gauss_moment(WC.GAUSS_T, a, b) for a, b in ((1, 1), (2, 0), (2, 2), (4, 0))],
              "C": C, "wick": [WC.wick4(C, 0, 0, 1, 1), WC.wick4(C, 0, 0, 0, 0)],
              "wc": [WC.wick_count(m) for m in range(1, 5)], "wce": [WC.wick_count_enum(m) for m in range(1, 5)],
              "t1": [WC.trotter_err(e) for e in WC.TROTTER_EPS], "t2": [WC.trotter_err(e, 2) for e in WC.TROTTER_EPS],
              "vw": [WC.visit_weight(k) for k in (1, 2)], "ifac": WC.intersection_factor(), "b72": WC.hs_book72(),
              "weak": [WC.weak_saw_count(WC.WEAK_N, u) for u in WC.WEAK_U],
              "weak_book": WC.weak_saw_count(WC.WEAK_N, 1.0, sign=+1),
              "G": {p: [WC.ladder_G(n, WC.LADDER_K, p) for n in (1, 2, 0)] for p in WC.LADDER_PAIRS},
              "Gi": {p: WC.ladder_G_ising(WC.LADDER_K, p) for p in WC.LADDER_PAIRS},
              "Gx": {p: WC.ladder_G_xy(WC.LADDER_K, p) for p in WC.LADDER_PAIRS},
              "tau": [WC.flip_tau(size) for size in WC.MC_L], "a3d": WC.book3d_alpha()[0], "slopes": WC.ex94_slopes()}
    return nb


def check_math():
    nb = numbers()
    assert [f4(v) for v in nb["hs"]] == ["0.7788", "0.3679"]
    assert abs(nb["hs"][0] - math.exp(-0.25)) < 1e-9 and abs(nb["hs"][1] - math.exp(-1)) < 1e-9
    assert all(abs(WC.hs_lhs(s) - v) < 1e-9 for s, v in zip(HS_S, nb["hs"]))
    assert f4(nb["gI"]) == "5.0627" and abs(nb["gI"] - 2 * math.pi / math.sqrt(3) * math.exp(1 / 3)) < 1e-12
    assert f4(nb["gIq"]) == "5.0627" and WC.GAUSS_T == ((2.0, -1.0), (-1.0, 2.0)) and WC.GAUSS_H == (1.0, 0.0)
    mq = [Fraction(1, 3), Fraction(2, 3), Fraction(2, 3), Fraction(4, 3)]
    assert all(abs(v - float(q)) < 1e-12 for v, q in zip(nb["mom"], mq))
    C = nb["C"]
    assert abs(C[0][1] - 1 / 3) < 1e-12 and abs(C[0][0] - 2 / 3) < 1e-12
    assert abs(nb["wick"][0] - float(mq[2])) < 1e-12 and abs(nb["wick"][1] - float(mq[3])) < 1e-12
    assert Fraction(4, 9) + Fraction(2, 9) == mq[2] and 3 * Fraction(2, 3) ** 2 == mq[3]
    assert nb["wc"] == [1, 3, 15, 105] == nb["wce"] and len(PAIRINGS) == nb["wc"][1]
    assert [g4(v) for v in nb["t1"]] == ["0.005784", "0.001502", "0.0003827"]
    assert [g4(v) for v in nb["t2"]] == ["0.0009465", "0.000124", r"1.589\times10^{-5}"]
    assert [f2(nb["t1"][i] / nb["t1"][i + 1]) for i in (0, 1)] == ["3.85", "3.92"]
    assert [f2(nb["t2"][i] / nb["t2"][i + 1]) for i in (0, 1)] == ["7.63", "7.81"]
    assert WC.TROTTER_EPS == (0.1, 0.05, 0.025) and WC.trotter_mats.__wrapped__.__defaults__ == (6, 906)
    assert abs(nb["vw"][0] - math.exp(-0.25)) < 1e-15 and abs(nb["vw"][1] - math.exp(-1)) < 1e-15
    assert f4(nb["ifac"]) == "0.6065" and abs(nb["ifac"] - math.exp(-0.5)) < 1e-12
    assert all(WC.pair_sum(k) == k + 2 * math.comb(k, 2) for k in range(1, 6))
    b = nb["b72"]
    assert f"{b.real:.4f}" == "0.3894" and f"{b.imag:.4f}" == "-0.2395" and f4(math.exp(-0.25)) == "0.7788"
    assert WC.WEAK_N == 8 and WC.WEAK_U[:3] == (0.0, 1.0, 4.0) and math.isinf(WC.WEAK_U[3])
    assert [fcount(v) for v in nb["weak"]] == ["65536", "23695.11", "7988.15", "5916"]
    assert nb["weak"][0] == 4 ** WC.WEAK_N and nb["weak"][-1] == SW.saw_c(WC.WEAK_N)
    assert sci1(nb["weak_book"]) == r"6.2\times10^{5}" and nb["weak_book"] > nb["weak"][0]
    assert WC.LADDER[0] == 6 and len(WC.LADDER[1]) == 7 and WC.LADDER_K == 0.5 and WC.LADDER_PAIRS == ((0, 3), (0, 5))
    g03, g05 = nb["G"][(0, 3)], nb["G"][(0, 5)]
    assert [f4(g03[0]), f4(g03[1]), f5(g03[2])] == ["0.6027", "0.5610", "0.65625"]
    assert [f4(g05[0]), f4(g05[1]), f5(g05[2])] == ["0.3562", "0.3171", "0.40625"]
    K = WC.LADDER_K
    assert g03[2] == K + K ** 3 + K ** 5 and g05[2] == 3 * K ** 3 + K ** 5
    for p in WC.LADDER_PAIRS:
        g = nb["G"][p]
        assert abs(nb["Gi"][p] - g[0]) < 1e-12 and abs(nb["Gx"][p] - g[1]) < 1e-9
        assert abs(g[2] - WC.ladder_saw_sum(K, p)) < 1e-15
        assert sorted(len(q) - 1 for q in ladder_paths(*p)) == WC.ladder_saws(*p)
    assert WC.ladder_saws(0, 3) == [1, 3, 5] and WC.ladder_saws(0, 5) == [3, 3, 3, 5]
    assert len(set(SAW_DEMO)) == len(SAW_DEMO) and len(set(CROSS_DEMO)) == len(CROSS_DEMO) - 1
    # 第 9 章回顧
    assert [f"{t:.0f}" for t in nb["tau"]] == ["133", "608", "2151"]
    assert f2(WC.BOOK_3D["alpha"]) == "0.12" and f2(nb["a3d"]) == "0.19"
    assert all(WC.flory_p(d) == Fraction(6, d + 2) for d in (1, 2, 3, 4))
    assert WC.EXPONENTS["mf"]["beta"] == Fraction(1, 2) and WC.EXPONENTS["ising2d"]["beta"] == Fraction(1, 8)
    assert abs(WC.ising1d_chi(1.0) - math.exp(2.0)) < 1e-12 and WC.landau_closed(-1.0) == (1.0, -0.25)
    assert [f4(x) for x in nb["slopes"][:3]] == ["0.5000", "1.5000", "0.0000"]
    assert all(abs(WC.intersection_factor(u) - math.exp(-u / 2)) < 1e-12 for u in (0.5, 1.0, 4.0))


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框 [x0,x1]×[y0,y1]（只縮不放），再置中在方框裡。"""
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


# ---------------------------------------------------------------- 畫圖元件
def saw_demo(cell=0.62):
    """−2..2 的格點＋一條 8 步自迴避漫步（起點 i 綠、終點 j 紅）→ VGroup(格點, 足跡, 標籤)。"""
    dots = L.lattice_dots(2, cell, radius=0.04)
    walk = L.trail(SAW_DEMO, cell, color=PATH_COLOR, stroke_width=4, end_radius=0.09)
    labs = VGroup(tx("i", 30, START_COLOR).next_to(L.lat_point(*SAW_DEMO[0], cell), DL, buff=0.08),
                  tx("j", 30, END_COLOR).next_to(L.lat_point(*SAW_DEMO[-1], cell), UR, buff=0.08))
    return VGroup(dots, walk, labs)


def cross_demo(cell=0.62):
    """有一個格點走過兩次的路：格點、足跡、圈出重複的格點＋「k = 2」→ VGroup(格點, 足跡, 圈, 標籤)。"""
    dots = L.lattice_dots(2, cell, radius=0.04)
    walk = L.trail(CROSS_DEMO, cell, color=PATH_COLOR, stroke_width=4, end_radius=0.09)
    p = L.lat_point(*CROSS_DEMO[1], cell)
    ring = Circle(radius=0.2).move_to(p).set_stroke(RED, 3)
    lab = tx("k=2", 28, RED).next_to(ring, UL, buff=0.05)
    return VGroup(dots, walk, ring, lab)


def pairing_row(dx=0.45):
    """⟨x_r x_s x_t x_u⟩ = 三種兩兩配對（弧線圖，同第 7 章）→ VGroup(左邊, 圖1, +, 圖2, +, 圖3)；每張圖 VGroup(點, 弧, 標籤)。"""
    parts = [tx(r"\langle x_rx_sx_tx_u\rangle=", 32)]
    for k, (pairs, lab) in enumerate(zip(PAIRINGS, PAIR_TEX)):
        pts = [i * dx * RIGHT for i in range(4)]
        dots = VGroup(*[Dot(p, radius=0.055).set_color(WHITE) for p in pts])
        arcs = VGroup(*[ArcBetweenPoints(pts[a], pts[b], angle=-0.9 * PI).set_stroke(TEAL, 3) for a, b in pairs])
        name = tx(lab, 24, GREY_A).next_to(dots, DOWN, buff=0.18)
        if k:
            parts.append(tx("+", 32))
        parts.append(VGroup(dots, arcs, name))
    return VGroup(*parts).arrange(RIGHT, buff=0.3)


def time_slices(width=5.0, height=2.4):
    """路徑積分示意：M+1 條時間切片（虛線）＋一條折線依序穿過（固定的示意高度）→ VGroup(切片, 折線, 點, 標籤, 插入的完備集)。"""
    M = len(SLICE_H) - 1
    xs = [width * k / M for k in range(M + 1)]
    slices = VGroup(*[DashedLine(np.array([x, 0, 0]), np.array([x, height, 0]), dash_length=0.08).set_stroke(GREY_C, 1.5)
                      for x in xs])
    pts = [np.array([x, h * height, 0]) for x, h in zip(xs, SLICE_H)]
    path = VMobject().set_points_as_corners(pts).set_stroke(PATH_COLOR, 3)
    dots = VGroup(*[Dot(p, radius=0.06).set_color(PATH_COLOR) for p in pts])
    labs = VGroup(tx(r"\vec R_0", 24).next_to(pts[0], LEFT, buff=0.12), tx(r"\vec R_M", 24).next_to(pts[-1], RIGHT, buff=0.12),
                  tx("0", 22, GREY_B).next_to(slices[0], DOWN, buff=0.1), tx("t", 22, GREY_B).next_to(slices[-1], DOWN, buff=0.1),
                  tx(r"\epsilon=t/M", 22, GREY_B).next_to(slices, DOWN, buff=0.1))
    ins = tx(r"\sum_{\vec R}|\vec R\rangle\langle\vec R|", 26, GREY_A).next_to(slices, UP, buff=0.12)
    return VGroup(slices, path, dots, labs, ins)


def weak_bars(vals, ux=1.35, height=2.4):
    """N = 8 的加權條數（u = 0、1、4、∞）：長條＋數值＋u 標籤 → VGroup(長條, 底線, 數值, u 標籤)。"""
    chart = D.bar_chart({k: v for k, v in enumerate(vals)}, ux=ux, uy=height / max(vals), width=0.6)
    bars = chart[0]
    bars[-1].set_fill(YELLOW, 0.8).set_stroke(YELLOW, 1.5)
    base = Line(-0.5 * ux * RIGHT, (len(vals) - 0.5) * ux * RIGHT).set_stroke(GREY_B, 2)
    nums = VGroup(*[tx(fcount(v), 24).next_to(b, UP, buff=0.08) for v, b in zip(vals, bars)])
    us = VGroup(*[tx("u=%s" % (r"\infty" if math.isinf(u) else f"{u:g}"), 24, GREY_A).next_to(b, DOWN, buff=0.15)
                  for u, b in zip(WC.WEAK_U, bars)])
    return VGroup(bars, base, nums, us)


def ladder_graph(cell=1.25):
    """2×3 梯子（格點 (x, y) → x + 3y）：鍵、格點、編號 → VGroup(鍵, 點, 編號)。"""
    pts = [L.lat_point(i % 3, i // 3, cell) for i in range(WC.LADDER[0])]
    bonds = VGroup(*[Line(pts[i], pts[j]).set_stroke(GREY_B, 3) for i, j in WC.LADDER[1]])
    dots = VGroup(*[Dot(p, radius=0.1).set_color(WHITE) for p in pts])
    labs = VGroup(*[tx(str(i), 24, GREY_A).next_to(p, DOWN if i < 3 else UP, buff=0.14) for i, p in enumerate(pts)])
    return VGroup(bonds, dots, labs)


def _dfs(path, b, nbr, out):
    if path[-1] == b:
        out.append(list(path))
        return
    for w in nbr[path[-1]]:
        if w not in path:
            _dfs(path + [w], b, nbr, out)


def ladder_paths(a, b):
    """梯子上從 a 到 b 的全部自迴避路徑（格點序列），由短到長。"""
    nbr = {i: [] for i in range(WC.LADDER[0])}
    for i, j in WC.LADDER[1]:
        nbr[i].append(j)
        nbr[j].append(i)
    out = []
    _dfs([a], b, nbr, out)
    return sorted(out, key=len)


def path_mob(graph, path, color=PATH_COLOR):
    """沿著梯子格點（畫好之後的位置）的一條折線（一個 VMobject）。"""
    return VMobject().set_points_as_corners([graph[1][v].get_center().copy() for v in path]).set_stroke(color, 8, 0.85)


# ---------------------------------------------------------------- 各段
def part_question(s):
    saw = saw_demo()
    e45 = tx(r"\lim_{n\to0}\frac1n\sum_\alpha\langle S_{i\alpha}S_{j\alpha}\rangle=\sum_NC(N;\vec R_i,\vec R_j)\,z^N", 42)
    tag = mixed("目標 (9.45)：", "$n\\to0$", " 的磁鐵", font_size=26, color=GREY_A)
    lab_l = mixed("$O(n)$", " 的關聯函數", font_size=24, color=BLUE_B)
    lab_r = mixed("$i$", " 到 ", "$j$", " 的 ", "$N$", " 步自迴避漫步條數", font_size=24, color=YELLOW)
    tag.next_to(e45, UP, buff=0.4, aligned_edge=LEFT)
    labs = VGroup(lab_l, lab_r).arrange(RIGHT, buff=0.9).next_to(e45, DOWN, buff=0.35)
    right = VGroup(tag, e45, labs)
    place_in(saw, -6.7, -3.3, -2.5, 2.7)
    place_in(right, -2.9, 6.7, -1.4, 2.2)
    assert_in_frame(saw, right)
    assert_apart(saw, right)
    assert_apart(lab_l, lab_r)
    assert_column(tag, e45, labs)
    s.play(FadeIn(saw[0]), run_time=0.6)
    s.play(ShowCreation(saw[1][0]), FadeIn(saw[1][1:]), FadeIn(saw[2]), run_time=1.2)
    s.say("問：「", "$n=0$", " 個分量的磁鐵」聽起來荒謬——", wait=1.3)
    s.say("它為什麼會數自迴避漫步？", wait=1.3)
    s.play(FadeIn(tag), Write(e45), run_time=1.2)
    s.say("目標 (9.45)：", "$n\\to0$", " 時，關聯函數除以 ", "$n$", wait=1.5)
    s.play(FadeIn(labs), run_time=0.6)
    s.say("等於從 ", "$i$", " 到 ", "$j$", " 的自迴避漫步，每一步乘 ", "$z$", wait=1.6)
    s.clear_stage()


def part_hs(s):
    nb = numbers()
    hs = nb["hs"]
    e46 = mixed("(9.46) ", "$e^{-\\frac u4|\\vec S|^4}=\\frac{1}{\\sqrt{\\pi u}}\\int_{-\\infty}^{\\infty}d\\phi\\,"
                           "e^{-\\phi^2/u-i\\phi|\\vec S|^2}$", font_size=30)
    num = mixed("數值（", "$u=1$", "）：", "$|\\vec S|^2=1$", " → ", "$%s=e^{-1/4}$" % f4(hs[0]), "；",
                "$|\\vec S|^2=2$", " → ", "$%s=e^{-1}$" % f4(hs[1]), font_size=26, color=D.EXACT_COLOR)
    chain = VGroup(mixed("四次項變成 ", "$i\\phi|\\vec S|^2$", "：對 ", "$\\vec S$", " 只剩二次", font_size=28),
                   mixed("每個分量是同一個高斯積分 → ", "$Z_n=(\\cdots)^n$", font_size=28),
                   mixed("$n\\to0$", "：", "$Z_n\\to1$", font_size=30, color=YELLOW)).arrange(DOWN, buff=0.28)
    gI = mixed("(9:S-11) ", "$\\int d^2x\\,e^{-\\frac12\\vec x\\cdot\\hat T\\vec x+\\vec h\\cdot\\vec x}="
                            "\\frac{2\\pi}{\\sqrt{\\det\\hat T}}e^{\\frac12\\vec h\\cdot\\hat T^{-1}\\vec h}$", font_size=28)
    place_in(e46, -6.7, 6.7, 2.2, 3.1)
    place_in(num, -6.7, 6.7, 1.45, 1.95)
    place_in(chain, -6.7, 6.7, -0.85, 1.15)
    place_in(gI, -6.7, 6.7, -2.9, -1.2)
    assert_in_frame(e46, num, chain, gI)
    assert_column(e46, num, chain, gI)
    s.play(FadeIn(e46), run_time=0.8)
    s.say("(9.46)：輔助場 ", "$\\phi$", " 把四次項換掉", wait=1.4)
    s.play(FadeIn(num), run_time=0.6)
    s.say("數值核對：", "$|\\vec S|^2=1,\\ 2$", " → ", "$%s,\\ %s$" % (f4(hs[0]), f4(hs[1])), wait=1.4)
    s.play(FadeIn(chain[0]), run_time=0.6)
    s.say("指數裡只剩 ", "$|\\vec S|^2$", "：對 ", "$\\vec S$", " 是高斯積分", wait=1.4)
    s.play(FadeIn(chain[1]), run_time=0.6)
    s.say("$n$", " 個分量各給一樣的因子，乘起來是 ", "$(\\cdots)^n$", wait=1.4)
    s.play(FadeIn(chain[2]), run_time=0.6)
    s.say("$n\\to0$", " 時 ", "$Z_n\\to1$", "：分母消失，只剩分子", wait=1.5)
    s.play(FadeIn(gI), run_time=0.8)
    s.say("(9:S-11)：二次型的高斯積分有閉式", wait=1.3)
    s.clear_stage()


def part_wick(s):
    nb = numbers()
    gT = mixed("$\\hat T=\\begin{pmatrix}2&-1\\\\-1&2\\end{pmatrix}$", "，", "$\\vec h=(1,0)$", "：",
               "$\\frac{2\\pi}{\\sqrt3}e^{1/3}=%s$" % f4(nb["gI"]), font_size=28)
    quad = mixed("（二重數值積分：", "$%s$" % f4(nb["gIq"]), "）", font_size=22, color=GREY_B)
    top = VGroup(gT, quad).arrange(RIGHT, buff=0.4)
    mom = mixed("(9:S-18) ", "$\\langle x_1x_2\\rangle=(\\hat T^{-1})_{12}=\\frac13$", "，", "$\\langle x_1^2\\rangle=\\frac23$",
                font_size=28)
    wick = mixed("(9:S-20) ", "$\\langle x_1^2x_2^2\\rangle=\\frac49+\\frac29=\\frac23$", "，",
                 "$\\langle x_1^4\\rangle=3\\cdot\\big(\\frac23\\big)^2=\\frac43$", font_size=28)
    arcs = pairing_row()
    cnt = mixed("配對數 ", "$\\frac{(2m)!}{2^mm!}=%s$" % ",\\ ".join(str(c) for c in nb["wc"]), "（", "$m=1..4$",
                "；第 7 章的弧線圖）", font_size=26)
    place_in(top, -6.7, 6.7, 2.05, 3.1)
    place_in(mom, -6.7, 6.7, 1.2, 1.8)
    place_in(wick, -6.7, 6.7, 0.3, 0.95)
    place_in(arcs, -6.7, 6.7, -1.75, 0.0)
    place_in(cnt, -6.7, 6.7, -2.9, -2.1)
    assert_in_frame(top, mom, wick, arcs, cnt)
    assert_column(top, mom, wick, arcs, cnt)
    s.play(FadeIn(gT), run_time=0.8)
    s.say("例：", "$2\\times2$", " 的 ", "$\\hat T$", "、", "$\\vec h=(1,0)$", "：", "$%s$" % f4(nb["gI"]), "（數值積分相同）",
          wait=1.4)
    s.play(FadeIn(quad), run_time=0.5)
    s.play(FadeIn(mom), run_time=0.6)
    s.say("二階矩就是逆矩陣：", "$\\langle x_1x_2\\rangle=\\frac13$", "、", "$\\langle x_1^2\\rangle=\\frac23$", wait=1.5)
    s.play(FadeIn(arcs[0]), run_time=0.5)
    s.play(LaggedStart(*[FadeIn(m) for m in arcs[1:]], lag_ratio=0.3), run_time=1.2)
    s.say("Wick：四階矩＝所有兩兩配對的和", wait=1.4)
    s.play(FadeIn(wick), run_time=0.6)
    s.say("$\\langle x_1^2x_2^2\\rangle=\\frac49+\\frac29$", "；", "$\\langle x_1^4\\rangle=3\\cdot(\\frac23)^2$", wait=1.5)
    s.play(FadeIn(cnt), run_time=0.6)
    s.say("配對數 ", "、".join(str(c) for c in nb["wc"]), "——第 7 章的弧線圖", wait=1.5)
    s.clear_stage()


def part_trotter(s):
    nb = numbers()
    t1, t2 = nb["t1"], nb["t2"]
    e66 = mixed("(9.66) ", "$e^{(\\frac12\\nabla^2-i\\phi)\\epsilon}=e^{-i\\phi\\epsilon}e^{\\frac12\\nabla^2\\epsilon}e^{O(\\epsilon^2)}$",
                font_size=30)
    nxt = mixed("Ex. 9.6 下一項：", "$e^{\\frac{i\\epsilon^2}{4}[\\phi,\\nabla^2]}$", font_size=26, color=GREY_A)
    hd = [tx(r"\epsilon", 28, YELLOW), zh("一階誤差", font_size=22, color=YELLOW), zh("比值", font_size=22, color=YELLOW),
          zh("乘上交換子項", font_size=22, color=YELLOW), zh("比值", font_size=22, color=YELLOW)]
    rows = [hd]
    for k, e in enumerate(WC.TROTTER_EPS):
        r1 = f2(t1[k - 1] / t1[k]) if k else VectorizedPoint()
        r2 = f2(t2[k - 1] / t2[k]) if k else VectorizedPoint()
        rows.append([f"{e:g}", g4(t1[k]), r1, g4(t2[k]), r2])
    rows.append([zh("趨勢", font_size=22, color=GREEN), tx(r"\propto\epsilon^2", 28, GREEN), tx(r"\to4", 28, GREEN),
                 tx(r"\propto\epsilon^3", 28, GREEN), tx(r"\to8", 28, GREEN)])
    tab = D.tex_table(rows, [1.0, 1.9, 1.0, 2.3, 1.0], row_h=0.66, font_size=26)
    sl = time_slices()
    place_in(e66, -6.7, 6.7, 2.4, 3.1)
    place_in(nxt, -6.7, 6.7, 1.75, 2.25)
    place_in(tab, -6.7, 0.5, -2.9, 1.45)
    place_in(sl, 0.9, 6.7, -2.6, 1.45)
    assert_in_frame(e66, nxt, tab, sl)
    assert_column(e66, nxt, tab)
    assert_column(e66, nxt, sl)
    assert_apart(tab, sl)
    assert_table(tab, 5)
    cells = tab[0]
    r1 = "、".join(f2(t1[i] / t1[i + 1]) for i in (0, 1))
    r2 = "、".join(f2(t2[i] / t2[i + 1]) for i in (0, 1))
    s.play(FadeIn(e66), run_time=0.8)
    s.say("(9.66)：", "$\\frac12\\nabla^2$", " 和 ", "$-i\\phi$", " 不對易，不能直接拆開", wait=1.5)
    s.say("Trotter：一小段 ", "$\\epsilon$", " 可以拆，誤差是 ", "$O(\\epsilon^2)$", wait=1.4)
    s.play(FadeIn(nxt), run_time=0.6)
    s.play(FadeIn(cells[:5]), ShowCreation(tab[1]), run_time=0.6)
    s.play(FadeIn(cells[5:20]), run_time=1.0)
    s.say("6 點環、種子 906：", "$\\epsilon$", " 減半，一階誤差除以 %s → 4：" % r1, "$\\propto\\epsilon^2$", wait=1.8)
    s.say("乘上交換子項：除以 %s → 8：" % r2, "$\\propto\\epsilon^3$", wait=1.5)
    s.play(FadeIn(cells[20:]), run_time=0.6)
    s.play(ShowCreation(sl[0]), FadeIn(sl[4]), run_time=0.6)
    s.play(ShowCreation(sl[1]), FadeIn(sl[2]), FadeIn(sl[3]), run_time=1.0)
    s.say("切成 ", "$M$", " 小段，中間插入 ", "$\\sum_{\\vec R}|\\vec R\\rangle\\langle\\vec R|$", wait=1.5)
    s.say("(9.67)–(9.71)：每段走一小步——變成對路徑的積分", wait=1.4)
    s.clear_stage()


def part_cross(s):
    nb = numbers()
    vw, ifac = nb["vw"], nb["ifac"]
    b = nb["b72"]
    e72 = mixed("(9.72)–(9.74) 走過 ", "$k$", " 次：", "$\\frac{1}{\\sqrt{\\pi u}}\\int_{-\\infty}^{\\infty}d\\phi\\,"
                "e^{-\\phi^2/u-ik\\phi}=e^{-uk^2/4}$", font_size=28)
    kl = mixed("$k=1$", "：", "$e^{-u/4}$", "；", "$k=2$", "：", "$e^{-u}$", "；每相交一次 ",
               "$e^{-u}/(e^{-u/4})^2=e^{-u/2}=%s$" % f4(ifac), "（", "$u=1$", "）", font_size=26)
    sm = mixed("$\\exp\\Big(-\\frac u4\\sum_{l,l'}\\delta\\Big)$", "：", "$k^2=k+2\\binom k2$", font_size=26)
    demo = cross_demo()
    box = D.erratum_box("p. 226，(9.72)", r"\int_0^\infty", r"\int_{-\infty}^{\infty}", width=5.6, font_size=40)
    place_in(e72, -6.7, 6.7, 2.45, 3.1)
    place_in(kl, -6.7, 6.7, 1.75, 2.25)
    place_in(sm, -6.7, 6.7, 0.6, 1.5)
    place_in(demo, -6.7, -2.6, -2.9, 0.3)
    place_in(box, -2.2, 6.7, -2.5, 0.05)
    assert_in_frame(e72, kl, sm, demo, box)
    assert_column(e72, kl, sm, demo)
    assert_column(e72, kl, sm, box)
    assert_apart(demo, box)
    s.play(FadeIn(e72), run_time=0.8)
    s.play(FadeIn(demo[0]), ShowCreation(demo[1][0]), FadeIn(demo[1][1:]), run_time=1.2)
    s.play(ShowCreation(demo[2]), FadeIn(demo[3]), run_time=0.6)
    s.say("(9.72)–(9.74)：格點走過 ", "$k$", " 次，對 ", "$\\phi$", " 積分得 ", "$e^{-uk^2/4}$", wait=1.5)
    s.say("$k=1$", " 是 ", "$e^{-u/4}$", "、", "$k=2$", " 是 ", "$e^{-u}$", wait=1.3)
    s.play(FadeIn(kl), run_time=0.6)
    s.say("每相交一次多乘 ", "$e^{-u/2}=%s$" % f4(ifac), "（", "$u=1$", "）", wait=1.4)
    s.play(FadeIn(sm), run_time=0.6)
    s.say("$k^2=k+2\\binom k2$", "：每次經過各自的權重＋每一對相交", wait=1.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (9.72) 從 0 積到 ", "$\\infty$", "；應該從 ", "$-\\infty$", " 積起", wait=1.4)
    s.say("從 0 積得 ", "$%s-%si$" % (f"{b.real:.4f}", f"{abs(b.imag):.4f}"), "，不是 ", "$e^{-1/4}=%s$" % f4(vw[0]),
          wait=1.5)
    s.clear_stage()


def part_weak(s):
    nb = numbers()
    w = nb["weak"]
    title = mixed("$N=%d$" % WC.WEAK_N, "：全部 ", "$4^{%d}=%s$" % (WC.WEAK_N, fcount(w[0])), " 條漫步，每相交一次乘 ",
                  "$e^{-u/2}$", font_size=28)
    bars = weak_bars(w)
    note = VGroup(mixed("$u=0$", "：每一條都算", font_size=26),
                  mixed("$u\\to\\infty$", "：只剩自迴避", font_size=26),
                  mixed("$c_%d=%s$" % (WC.WEAK_N, fcount(w[-1])), font_size=26, color=YELLOW)).arrange(DOWN, aligned_edge=LEFT,
                                                                                                    buff=0.3)
    box = D.erratum_box("p. 226，(9.75)",
                        r"-\int_0^t\frac d2\Big|\frac{d\vec R}{ds}\Big|^2ds\ \frac u4\iint\delta",      # \!\! 與 tipa 衝突
                        r"-\int_0^t\frac d2\Big|\frac{d\vec R}{ds}\Big|^2ds-\frac u4\iint\delta", width=11.0, font_size=32)
    place_in(title, -6.7, 6.7, 2.5, 3.1)
    place_in(bars, -6.7, 0.4, -0.95, 2.2)
    place_in(note, 1.0, 6.7, -0.8, 2.0)
    place_in(box, -6.7, 6.7, -2.9, -1.25)
    assert_in_frame(title, bars, note, box)
    assert_column(title, bars, box)
    assert_column(title, note, box)
    assert_apart(bars, note)
    for a, c in zip(bars[2], bars[2][1:]):
        assert_apart(a, c, 0.03)
    ws = "、".join(fcount(v) for v in w[:3])
    s.play(FadeIn(title), run_time=0.8)
    s.say("$N=%d$" % WC.WEAK_N, "：全部 ", "$4^{%d}$" % WC.WEAK_N, " 條漫步，每相交一次乘 ", "$e^{-u/2}$", wait=1.4)
    s.play(ShowCreation(bars[1]), LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars[0]], lag_ratio=0.25),
           FadeIn(bars[3]), run_time=1.5)
    s.play(FadeIn(bars[2]), run_time=0.6)
    s.say("$u=0,\\ 1,\\ 4$", "：", ws, wait=1.4)
    s.play(FadeIn(note), run_time=0.6)
    s.say("$u\\to\\infty$", "：只剩自迴避的 ", "$c_%d=%s$" % (WC.WEAK_N, fcount(w[-1])), " 條", wait=1.4)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (9.75) 兩項之間少了一個負號", wait=1.3)
    s.say("少了負號：每相交乘 ", "$e^{+u/2}$", "，", "$u=1$", " 得 ", "$%s$" % sci1(nb["weak_book"]), "——越相交越多",
          wait=1.8)
    s.clear_stage()


def part_ladder(s):
    nb = numbers()
    g03, g05 = nb["G"][(0, 3)], nb["G"][(0, 5)]
    graph = ladder_graph()
    wl = mixed("權重 ", "$\\prod(1+K\\vec S_i\\cdot\\vec S_j)$", "，", "$|\\vec S|^2=n$", "：每個閉圈乘 ", "$n$", font_size=26)
    f03 = tx(r"G_{03}=\frac{K+K^3+K^5+nK^5}{1+2nK^4+nK^6}", 38)
    hd = [tx("K=0.5", 26, YELLOW), mixed("$n=1$", "（Ising）", font_size=20, color=YELLOW),
          mixed("$n=2$", "（xy）", font_size=20, color=YELLOW), tx(r"n\to0", 26, YELLOW), zh("自迴避路徑", font_size=20, color=YELLOW)]
    assert f"{WC.LADDER_K:g}" == "0.5"
    rows = [hd, [tx("G_{03}", 28), f4(g03[0]), f4(g03[1]), f5(g03[2]), r"K+K^3+K^5"],
            [tx("G_{05}", 28), f4(g05[0]), f4(g05[1]), f5(g05[2]), r"3K^3+K^5"]]
    tab = D.tex_table(rows, [1.4, 1.75, 1.6, 1.5, 2.2], row_h=0.72, font_size=28)
    for i in (1, 2):
        D.table_cell(tab, i, 3, 5).set_color(YELLOW)
    e76 = mixed("(9.76)(9.77) ", "$G=\\sum_NC(N)\\,e^{-rN/2}$", "，", "$e^{-r/2}=z/z_c$", "：關聯函數＝自迴避漫步的生成函數",
                font_size=26)
    place_in(graph, -6.7, -2.7, 0.75, 3.0)
    place_in(VGroup(wl, f03).arrange(DOWN, buff=0.35), -2.4, 6.7, 1.35, 3.1)
    place_in(tab, -2.4, 6.7, -1.35, 1.05)
    place_in(e76, -6.7, 6.7, -2.9, -1.75)
    term_pt = np.array([graph.get_center()[0], 0.25, 0])
    terms03 = [tx(t, 32, YELLOW).move_to(term_pt) for t in ("K", "K+K^3", "K+K^3+K^5")]
    terms05 = [tx(t, 32, YELLOW).move_to(term_pt) for t in ("K^3", "2K^3", "3K^3", "3K^3+K^5")]
    assert_in_frame(graph, wl, f03, tab, e76, *terms03, *terms05)
    assert_column(wl, f03, tab, e76)
    assert_apart(graph, VGroup(wl, f03, tab))
    assert_apart(graph, terms03[-1])
    assert_apart(VGroup(*terms03, *terms05), tab)
    assert_apart(VGroup(*terms03, *terms05), e76)
    assert_table(tab, 5)
    p03, p05 = ladder_paths(0, 3), ladder_paths(0, 5)
    dots = graph[1]
    s.play(ShowCreation(graph[0]), FadeIn(dots), FadeIn(graph[2]), run_time=1.0)
    s.say("小梯子：", "$2\\times3$", "、6 個點、7 條鍵（格點 ", "$(x,y)\\to x+3y$", "）", wait=1.3)
    s.play(FadeIn(wl), run_time=0.6)
    s.say("權重展開後只剩「路」和「閉圈」；", "$|\\vec S|^2=n$", "：每個閉圈乘 ", "$n$", wait=1.8)
    s.play(dots[0].animate.set_color(START_COLOR), dots[3].animate.set_color(END_COLOR), run_time=0.5)
    s.say("$G_{03}$", "：從 0 到 3 的自迴避路徑，長 1、3、5", wait=0.2)
    prev = None
    for path, term in zip(p03, terms03):
        pm = path_mob(graph, path)
        anims = [ShowCreation(pm), FadeIn(term)] + ([FadeOut(prev[0]), FadeOut(prev[1])] if prev else [])
        s.play(*anims, run_time=0.8)
        s.wait(0.3)
        prev = (pm, term)
    s.play(FadeIn(f03), FadeOut(prev[0]), run_time=0.8)
    s.say("分子多出的 ", "$nK^5$", "：一條路再加一個閉圈", wait=1.4)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), FadeIn(tab[0][5:10]), run_time=0.8)
    s.say("$K=0.5$", "：Ising ", "$%s$" % f4(g03[0]), "、xy ", "$%s$" % f4(g03[1]), "、", "$n\\to0$", " ",
          "$%s$" % f5(g03[2]), wait=1.5)
    s.say("$n\\to0$", "：閉圈全部消失，只剩 ", "$K+K^3+K^5$", wait=1.4)
    s.play(FadeOut(prev[1]), dots[3].animate.set_color(WHITE), dots[5].animate.set_color(END_COLOR), run_time=0.5)
    s.say("對角 ", "$G_{05}$", "：四條路，長 3、3、3、5", wait=0.2)
    prev = None
    for path, term in zip(p05, terms05):
        pm = path_mob(graph, path)
        anims = [ShowCreation(pm), FadeIn(term)] + ([FadeOut(prev[0]), FadeOut(prev[1])] if prev else [])
        s.play(*anims, run_time=0.7)
        s.wait(0.2)
        prev = (pm, term)
    s.play(FadeIn(tab[0][10:]), FadeOut(prev[0]), run_time=0.8)
    s.say("$n\\to0$", "：", "$%s=3K^3+K^5$" % f5(g05[2]), "——又是自迴避路徑的和", wait=1.4)
    s.play(FadeIn(e76), run_time=0.8)
    s.say("(9.76)(9.77)：", "$e^{-r/2}=z/z_c$", "——關聯函數＝自迴避漫步的生成函數", wait=1.5)
    s.say("書上 p. 227 的 ", "$r/2=1-T_c/T$", " 和 (9.21) 差 2 倍（只差常數倍）", wait=1.5)
    s.clear_stage()


def part_recap(s):
    nb = numbers()
    fz = 24
    tau = "、".join(f"{t:.0f}" for t in nb["tau"])
    b_mf, b_2d = WC.EXPONENTS["mf"]["beta"], WC.EXPONENTS["ising2d"]["beta"]
    rows = [[tx("W901", fz), mixed("$p=2\\nu=\\frac{6}{d+2}$", "；", "$\\beta=\\frac{%d}{%d}$" % (b_mf.numerator, b_mf.denominator),
                                   " 對 ", "$\\frac{%d}{%d}$" % (b_2d.numerator, b_2d.denominator), "（平均場對二維 Ising）",
                                   font_size=fz)],
            [tx("W902", fz), mixed("$\\langle M\\rangle=0$", "，但翻轉等待 ", tau, " sweep：對稱破缺", font_size=fz)],
            [tx("W903", fz), mixed("$\\chi=e^{2K}$", "、四條標度律、書上 ", "$\\alpha=%s\\to%s$" % (f2(WC.BOOK_3D["alpha"]),
                                                                                         f2(nb["a3d"])), font_size=fz)],
            [tx("W904", fz), mixed("雙井 ", "$m^*=\\sqrt{-r/u}$", "：", "$\\beta=\\frac12,\\ \\Delta=\\frac32,\\ \\alpha=0$",
                                   font_size=fz)],
            [tx("W905", fz), mixed("係數一除：", "$\\langle R^2\\rangle\\propto N^{2\\nu}$", "，前因子 ", "$z_c^{2\\nu}$",
                                   "（不是 ", "$z_c^{\\gamma}$", "）", font_size=fz)],
            [tx("W906", fz), mixed("每相交一次 ", "$e^{-u/2}$", "；", "$n\\to0$", " 只剩自迴避", font_size=fz)]]
    tab = D.tex_table(rows, [1.6, 10.4], row_h=0.78, font_size=fz)
    title = zh("第 9 章回顧", font_size=30, color=YELLOW)
    place_in(title, -6.7, 6.7, 2.6, 3.1)
    place_in(tab, -6.7, 6.7, -2.9, 2.35)
    assert_in_frame(title, tab)
    assert_apart(title, tab)
    assert_table(tab, 2)
    s.play(FadeIn(title), ShowCreation(tab[1][0]), ShowCreation(tab[1][2]), run_time=0.8)
    s.say("第 9 章回顧", wait=0.3)
    for i in range(6):
        s.play(FadeIn(tab[0][2 * i]), FadeIn(tab[0][2 * i + 1]), run_time=0.6)
        s.wait(0.7)
    s.say("$n\\to0$", " 時每個閉圈都帶一個 ", "$n$", " 而消失", wait=1.6)
    s.say("剩下的就是一條從 ", "$i$", " 到 ", "$j$", "、不碰自己的路", wait=1.6)
    s.say("$O(n)$", " 的關聯函數就是自迴避漫步的生成函數", wait=1.8)


class W906NToZero(ChuScene):
    def construct(self):
        self.show_header("W906　n→0 的魔術：O(n) 模型就是自迴避漫步")
        check_math()
        part_question(self)
        part_hs(self)
        part_wick(self)
        part_trotter(self)
        part_cross(self)
        part_weak(self)
        part_ladder(self)
        part_recap(self)
        self.wait(1.0)
