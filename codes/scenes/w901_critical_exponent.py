"""
W901 漫步的指數就是臨界指數（§9.1–9.2 開頭，pp. 203–205，(9.1)–(9.3)，Fig. 9.1）

問：自迴避漫步越走越大、磁鐵加熱就失去磁性——為什麼說這是同一個數學問題？
  1. 左：20 步二維自迴避漫步（SW.sample_saw(20)，種子 702）一步一步長出來；右：8×8 自旋 T = 1.5（64 上）→ T = 3.5（33 上、31 下）。
  2. (9.1)(9.2)：Flory p = 6/(d+2) 對已知 2ν（d = 1..4）；三維差 2%；二維精確列舉 2ν_eff = 1.4442、1.4505、1.4556 → 3/2。
  3. (9.3)、Fig. 9.1：平均場 m = tanh(m/T) 與 Onsager 的 M_s(T)；表 T/T_c = 0.5、0.9、0.99；斜率 β = 1/2、1/8（書上 1/3 是三維）。
  4. de Gennes：p = 2ν；漫步 ↔ n = 0 的 O(n) 磁鐵（對照表，無數字）。
畫面上的數字都由 utils.walk_crit（與 walk_saw）算出，check_math() 核對。
模組層級的畫圖函式：grid_lines、saw_box（足跡框）、spin_arrows（自旋箭頭）、ms_axes、ms_curve（磁化曲線）。
"""
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

SAW_STEPS = 20
SNAP_TS = WC.SNAP_T                    # (1.5, 3.5)
NEFF_N = (10, 12, 14)
DIMS = (1, 2, 3, 4)
UP_COLOR = "#F2A33A"                   # 自旋朝上：暖色
DOWN_COLOR = "#4A90D9"                 # 自旋朝下：冷色
MF_COLOR = TEAL_B
ONS_COLOR = YELLOW
GOOD = "#7CE38B"
BAD = RED_B


def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else r"\frac{%d}{%d}" % (q.numerator, q.denominator)


def ms_points(fn, tc, n=50):
    """(T/T_c, M_s) 折線的點：0..1 之間用 x = 1 − (1 − k/n)^4（越靠近 T_c 越密），再接 (1, 0)、(1.2, 0)。"""
    pts = [(0.0, 1.0)]
    for k in range(1, n):
        x = 1 - (1 - k / n) ** 4
        pts.append((x, fn(x * tc)))
    return pts + [(1.0, 0.0), (1.2, 0.0)]


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = {"saw": [tuple(int(v) for v in p) for p in SW.sample_saw(SAW_STEPS)],
              "snap": [np.array(WC.mc_snapshot(8, T)) for T in SNAP_TS],
              "flory": [WC.flory_p(d) for d in DIMS],
              "two_nu": [2 * WC.NU_KNOWN[d] for d in DIMS],
              "gap": WC.flory_gap_3d(),
              "neff": [SW.two_nu_eff(N) for N in NEFF_N],
              "ms": [(WC.mf_m(x), WC.onsager_M(x * WC.TC_2D)) for x in WC.MS_X],
              "ms_tab": WC.ms_table(),
              "beta": (WC.beta_slope(WC.mf_m, 1.0), WC.beta_slope(WC.onsager_M, WC.TC_2D)),
              "mf_pts": ms_points(WC.mf_m, 1.0),
              "ons_pts": ms_points(WC.onsager_M, WC.TC_2D)}
    return nb


def spin_counts(sp):
    return int((sp == 1).sum()), int((sp == -1).sum())


def check_math():
    nb = numbers()
    p = nb["saw"]
    assert len(p) == SAW_STEPS + 1 and len(set(p)) == SAW_STEPS + 1 and p[0] == (0, 0)
    assert all(abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1 for a, b in zip(p, p[1:]))
    lo, hi = nb["snap"]
    assert lo.shape == hi.shape == (8, 8)
    assert spin_counts(lo) == (64, 0) and spin_counts(hi) == (33, 31)
    assert f4(WC.TC_2D) == "2.2692" and SNAP_TS[0] < WC.TC_2D < SNAP_TS[1]
    assert nb["flory"] == [Fraction(2), Fraction(3, 2), Fraction(6, 5), Fraction(1)]
    assert [frac_tex(q) for q in nb["flory"]] == ["2", r"\frac{3}{2}", r"\frac{6}{5}", "1"]
    tn = nb["two_nu"]
    assert tn[0] == 2 and tn[1] == Fraction(3, 2) and tn[3] == 1 and f"{tn[2]:.3f}" == "1.176"
    assert f"{float(nb['flory'][2]) / 2:.3f}" == "0.600" and f"{WC.NU_KNOWN[3]:.3f}" == "0.588"
    assert f4(nb["gap"]) == "0.0204" and f"{100 * nb['gap']:.0f}" == "2"
    assert abs(float(nb["flory"][2]) / tn[2] - 1 - nb["gap"]) < 1e-12
    assert all(WC.flory_p(d) == 2 * SW.flory_nu(d) for d in DIMS)
    assert [f4(v) for v in nb["neff"]] == ["1.4442", "1.4505", "1.4556"]
    want = [("0.9575", "0.9980"), ("0.5254", "0.8959"), ("0.1725", "0.6859")]
    assert [(f4(a), f4(b)) for a, b in nb["ms"]] == want
    assert [(f4(a), f4(b)) for a, b in (nb["ms_tab"][x] for x in WC.MS_X)] == want
    assert abs(WC.mf_m(0.9) - WC.mf_m_newton(0.9)) < 1e-12
    assert [f4(v) for v in nb["beta"]] == ["0.5000", "0.1250"]
    assert WC.BETA_CITED["3d_ising"] == 0.326
    for key in ("mf_pts", "ons_pts"):
        pts = nb[key]
        assert len(pts) <= 60
        assert all(0 <= x <= 1.2 and 0 <= y <= 1 for x, y in pts)


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


def row_cells(tab, i, ncols):
    return tab[0][i * ncols:(i + 1) * ncols]


def dense_points(ax, pts, k=6):
    """折線上的密集取樣點（畫面座標），用來檢查文字沒有壓到曲線。"""
    out = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        for t in np.linspace(0, 1, k, endpoint=False):
            out.append(ax.c2p(x0 + t * (x1 - x0), y0 + t * (y1 - y0)))
    return out + [ax.c2p(*pts[-1])]


def assert_clear(mob, pts, gap=0.05):
    """mob 的外框（外擴 gap）裡沒有任何一個點。"""
    x0, x1 = mob.get_left()[0] - gap, mob.get_right()[0] + gap
    y0, y1 = mob.get_bottom()[1] - gap, mob.get_top()[1] + gap
    for p in pts:
        assert not (x0 <= p[0] <= x1 and y0 <= p[1] <= y1), (mob, p)


# ---------------------------------------------------------------- 畫圖元件
def grid_lines(x0, x1, y0, y1, cell, origin, color=GREY_D):
    """穿過格點的淡格線：x = x0..x1、y = y0..y1，各往外延半格。"""
    g = VGroup()
    for x in range(x0, x1 + 1):
        g.add(Line(L.lat_point(x, y0 - 0.5, cell, origin), L.lat_point(x, y1 + 0.5, cell, origin)))
    for y in range(y0, y1 + 1):
        g.add(Line(L.lat_point(x0 - 0.5, y, cell, origin), L.lat_point(x1 + 0.5, y, cell, origin)))
    return g.set_stroke(color, 1)


def saw_box(pts, side, center):
    """足跡框：VGroup(框, 格線, 足跡 VGroup(線, 起點, 終點), 起點到終點的虛線)；足跡的外框置中。"""
    center = np.array(center, dtype=float)
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    span = max(max(xs) - min(xs) + 1, max(ys) - min(ys) + 1)
    cell = (side - 0.4) / span
    mid = np.array([(max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2, 0.0])
    origin = center - mid * cell
    frame = Square(side_length=side).move_to(center).set_stroke(GREY_B, 1.5)
    grid = grid_lines(min(xs), max(xs), min(ys), max(ys), cell, origin)
    path = L.trail(pts, cell, origin, stroke_width=4, end_radius=0.08)
    rvec = DashedLine(path[1].get_center(), path[2].get_center()).set_stroke(WHITE, 2.5)
    return VGroup(frame, grid, path, rvec)


def spin_arrows(spins, cell, center):
    """L×L 的 ±1 自旋 → VGroup(L² 支箭頭)；朝上暖色、朝下冷色。
    朝下的箭頭就是朝上的箭頭轉 180°（點的順序相同），所以 Transform(path_arc=PI) 會把它們原地翻過來。"""
    n = spins.shape[0]
    center = np.array(center, dtype=float)
    g = VGroup()
    for i in range(n):
        for j in range(n):
            c = center + ((j - (n - 1) / 2) * RIGHT + ((n - 1) / 2 - i) * UP) * cell
            a = Arrow(c + 0.36 * cell * DOWN, c + 0.36 * cell * UP, buff=0, thickness=4,
                      max_width_to_length_ratio=0.16)
            if spins[i, j] > 0:
                a.set_fill(UP_COLOR, 1)
            else:
                a.rotate(PI, about_point=c).set_fill(DOWN_COLOR, 1)
            g.add(a)
    return g


def ms_axes(width=5.4, height=3.7):
    ax = D.plot_axes((0, 1.25, 0.25), (0, 1.1, 0.25), width=width, height=height, x_ticks=(0, 0.5, 1.0),
                     x_tick_labels=("0", "0.5", "1"), y_ticks=(0.5, 1.0), y_tick_labels=("0.5", "1"),
                     font_size=22)
    xl = tx(r"T/T_c", 26, GREY_B).next_to(ax.x_axis.get_end(), DOWN, buff=0.15)
    yl = tx(r"M_s", 26, GREY_B).next_to(ax.y_axis.get_end(), RIGHT, buff=0.15)
    return ax, VGroup(xl, yl)


def ms_curve(ax, pts, color):
    """一條磁化曲線＝一個 VMobject 折線（≤ 60 點）。"""
    return VMobject().set_points_as_corners([ax.c2p(x, y) for x, y in pts]).set_stroke(color, 3)


# ---------------------------------------------------------------- 各段
def part_intro(s):
    nb = numbers()
    left = saw_box(nb["saw"], 3.8, [-3.4, 0.2, 0])
    frame_r = Square(side_length=3.8).move_to([3.4, 0.2, 0]).set_stroke(GREY_B, 1.5)
    lo = spin_arrows(nb["snap"][0], 0.42, [3.4, 0.2, 0])
    hi = spin_arrows(nb["snap"][1], 0.42, [3.4, 0.2, 0])
    (u0, d0), (u1, d1) = spin_counts(nb["snap"][0]), spin_counts(nb["snap"][1])
    tl = zh("自迴避漫步：越走越大", font_size=30).next_to(left[0], UP, buff=0.2)
    tr = zh("磁鐵：加熱就失去磁性", font_size=30).next_to(frame_r, UP, buff=0.2)
    bl = mixed(f"$N={SAW_STEPS}$", " 步；虛線 ", "$\\vec R$", "：起點到終點", font_size=26).next_to(left[0], DOWN, buff=0.2)
    br0 = mixed(f"$T={SNAP_TS[0]}<T_c$", f"：朝上 {u0}、朝下 {d0}", font_size=26).next_to(frame_r, DOWN, buff=0.2)
    br1 = mixed(f"$T={SNAP_TS[1]}>T_c$", f"：朝上 {u1}、朝下 {d1}", font_size=26).next_to(frame_r, DOWN, buff=0.2)
    assert_in_frame(left, frame_r, lo, hi, tl, tr, bl, br0, br1)
    assert_apart(left, frame_r)
    for a, b in ((tl, left), (tr, frame_r), (bl, left), (br0, frame_r), (br1, frame_r), (tl, tr), (bl, br0)):
        assert_apart(a, b)
    s.say("問：自迴避漫步越走越大、磁鐵加熱就失去磁性")
    s.play(FadeIn(left[0]), FadeIn(left[1]), FadeIn(tl), FadeIn(frame_r), FadeIn(tr), run_time=0.8)
    s.say("為什麼說這是同一個數學問題？", wait=2.0)
    s.say("左：", f"$N={SAW_STEPS}$", " 步的二維自迴避漫步，不准踩回走過的點")
    s.play(FadeIn(left[2][1]), run_time=0.4)
    s.play(ShowCreation(left[2][0]), run_time=4.5, rate_func=linear)
    s.play(FadeIn(left[2][2]), ShowCreation(left[3]), FadeIn(bl), run_time=0.8)
    s.say("$\\langle R^2\\rangle$", "：", "$R^2$", " 對所有走法的平均；", "$N$", " 越大越大", wait=2.5)
    s.say("右：", "$8\\times8$", " 個自旋，", f"$T={SNAP_TS[0]}$", " 時全部朝上：有磁化")
    s.play(FadeIn(lo, lag_ratio=0.03), FadeIn(br0), run_time=1.5)
    s.wait(2.0)
    s.say("加熱到 ", f"$T={SNAP_TS[1]}$", f"：{u1} 個朝上、{d1} 個朝下，磁化消失")
    s.play(Transform(lo, hi, path_arc=PI), FadeTransform(br0, br1), run_time=2.0)
    s.wait(2.0)
    s.say(f"$T_c={f4(WC.TC_2D)}$", "（二維 Ising）以上沒有自發磁化", wait=3.0)
    s.say("兩邊都有一個「越來越大」或「消失」的冪次：臨界指數", wait=3.0)
    s.clear_stage()


def part_flory(s):
    nb = numbers()
    fl, tn = nb["flory"], nb["two_nu"]
    gap = f"{100 * nb['gap']:.0f}"
    e1 = VGroup(mixed("(9.1)", "$\\langle R^2\\rangle\\propto N^{p}$", font_size=30),
                mixed("(9.2)　Flory：", "$p=6/(d+2)$", font_size=30)).arrange(RIGHT, buff=1.2)
    rows = [["d"] + [str(d) for d in DIMS],
            [mixed("Flory ", "$p$", font_size=28)] + [frac_tex(q) for q in fl],
            [mixed("已知 ", "$2\\nu$", font_size=28)] + [frac_tex(tn[0]), frac_tex(tn[1]), f"{tn[2]:.3f}", frac_tex(tn[3])],
            [zh("比較", font_size=28)] + [zh("精確", font_size=28, color=GOOD), zh("精確", font_size=28, color=GOOD),
                                          zh(f"差 {gap}%", font_size=28, color=BAD), zh("精確", font_size=28, color=GOOD)]]
    tab = D.tex_table(rows, [2.6, 1.6, 1.6, 1.6, 1.6], row_h=0.85, font_size=32)
    for j in range(1, 5):
        D.table_cell(tab, 2, j, 5).set_color(YELLOW)
    lab = zh("二維精確列舉（第 7 章）的局部指數", font_size=26, color=GREY_B)
    ntab = D.tex_table([["N"] + [str(N) for N in NEFF_N] + [r"\to"],
                        [r"2\nu_{\rm eff}"] + [f4(v) for v in nb["neff"]] + [r"\frac32"]],
                       [2.0, 1.6, 1.6, 1.6, 1.1], row_h=0.72, font_size=32)
    D.table_cell(ntab, 1, 4, 5).set_color(YELLOW)
    bottom = VGroup(lab, ntab).arrange(DOWN, buff=0.2)
    place_in(e1, -6.7, 6.7, 2.55, 3.1)
    place_in(tab, -6.7, 6.7, -0.95, 2.45)
    place_in(bottom, -6.7, 6.7, -2.9, -1.1)
    assert_in_frame(e1, tab, bottom)
    assert_column(e1, tab, bottom)
    assert_table(tab, 5)
    assert_table(ntab, 5)
    s.say("(9.1)：自迴避漫步 ", "$\\langle R^2\\rangle\\propto N^{p}$", "，", "$p$", " 只看維度 ", "$d$")
    s.play(FadeIn(e1), run_time=1.0)
    s.say("(9.2)：Flory 的估計 ", "$p=\\frac{6}{d+2}$", "（第 7 章）", wait=2.5)
    s.play(FadeIn(row_cells(tab, 0, 5)), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in row_cells(tab, 1, 5)], lag_ratio=0.15), run_time=1.0)
    s.say("Flory：", "$d=1,2,3,4$", " 得 ", "$p=2,\\ \\frac32,\\ \\frac65,\\ 1$", wait=2.5)
    s.play(LaggedStart(*[FadeIn(c) for c in row_cells(tab, 2, 5)], lag_ratio=0.15), run_time=1.0)
    s.say("已知的 ", "$2\\nu$", "：一、二、四維和 Flory 一模一樣", wait=3.0)
    s.play(LaggedStart(*[FadeIn(c) for c in row_cells(tab, 3, 5)], lag_ratio=0.15), run_time=1.0)
    s.say("三維：Flory ", f"$\\nu={float(fl[2]) / 2:.3f}$", "，已知 ", f"$\\nu={WC.NU_KNOWN[3]:.3f}$", "（只引用）",
          wait=2.5)
    s.say(f"差 {gap}%：書上說「一個百分點左右」", wait=2.5)
    s.play(FadeIn(lab), FadeIn(row_cells(ntab, 0, 5)), ShowCreation(ntab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in row_cells(ntab, 1, 5)], lag_ratio=0.15), run_time=1.0)
    s.say("二維精確列舉：" + "、".join(f4(v) for v in nb["neff"]), wait=2.5)
    s.say("$N$", " 越大越靠近 ", "$\\frac32$", "：二維的 Flory 值是對的", wait=3.0)
    s.clear_stage()


def part_magnet(s):
    nb = numbers()
    ax, axl = ms_axes()
    graph = VGroup(ax, axl)
    place_in(graph, -6.7, -0.4, -2.9, 2.9)
    mf = ms_curve(ax, nb["mf_pts"], MF_COLOR)
    ons = ms_curve(ax, nb["ons_pts"], ONS_COLOR)
    tc_line = DashedLine(ax.c2p(1, 0), ax.c2p(1, 1.08)).set_stroke(GREY_B, 2)
    tc_lab = tx("T_c", 24, GREY_B).next_to(tc_line, UP, buff=0.06)
    dots = VGroup()
    for x, (a, b) in zip(WC.MS_X, nb["ms"]):
        dots.add(Dot(ax.c2p(x, a), radius=0.06).set_fill(MF_COLOR, 1), Dot(ax.c2p(x, b), radius=0.06).set_fill(ONS_COLOR, 1))
    leg = VGroup(mixed("平均場 ", "$m=\\tanh(m/T)$", font_size=24, color=MF_COLOR),
                 zh("二維 Ising（Onsager）", font_size=24, color=ONS_COLOR)).arrange(DOWN, buff=0.15, aligned_edge=LEFT)
    leg.move_to(ax.c2p(0.42, 0.28))
    e93 = mixed("(9.3)　", "$M_s\\propto(T_c-T)^{\\beta}$", font_size=32)
    eon = tx(r"M_s=\Big(1-\sinh^{-4}\frac{2}{T}\Big)^{1/8}", 32, ONS_COLOR)
    rows = [["T/T_c", zh("平均場", font_size=26, color=MF_COLOR), zh("二維 Ising", font_size=26, color=ONS_COLOR)]]
    for x, (a, b) in zip(WC.MS_X, nb["ms"]):
        rows.append([str(x), f4(a), f4(b)])
    rows.append([mixed("斜率 ", "$\\beta$", font_size=26), f4(nb["beta"][0]), f4(nb["beta"][1])])
    tab = D.tex_table(rows, [1.7, 1.7, 1.9], row_h=0.6, font_size=30)
    for j in (1, 2):
        D.table_cell(tab, 4, j, 3).set_color(YELLOW)
    concl = mixed("$\\beta=\\frac12$", "（平均場），", "$\\beta=\\frac18$", "（二維）", font_size=28, color=YELLOW)
    right = VGroup(e93, eon, tab, concl).arrange(DOWN, buff=0.3)
    place_in(right, 0.2, 6.7, -2.9, 3.1)
    assert_in_frame(graph, right, tc_lab, leg)
    assert_apart(graph, right)
    curve_pts = dense_points(ax, nb["mf_pts"]) + dense_points(ax, nb["ons_pts"])
    for m in (leg, tc_lab, axl[1]):
        assert_clear(m, curve_pts)
    assert_column(e93, eon, tab, concl)
    assert_table(tab, 3)
    s.say("磁鐵：", "$M_s$", " 是自發磁化（不加外場時的磁化）")
    s.play(ShowCreation(ax), FadeIn(axl), FadeIn(e93), run_time=1.0)
    s.say("(9.3)：靠近 ", "$T_c$", " 時 ", "$M_s\\propto(T_c-T)^{\\beta}$", wait=2.5)
    s.play(ShowCreation(mf), FadeIn(leg[0]), ShowCreation(tc_line), FadeIn(tc_lab), run_time=1.5)
    s.say("平均場：解 ", "$m=\\tanh(m/T)$", "（", "$T$", " 以 ", "$T_c$", " 為單位）", wait=2.5)
    s.play(ShowCreation(ons), FadeIn(leg[1]), FadeIn(eon), run_time=1.5)
    s.say("二維 Ising：Onsager 的精確解，", f"$T_c={f4(WC.TC_2D)}$", wait=2.5)
    s.play(FadeIn(dots), FadeIn(row_cells(tab, 0, 3)), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][3:12]], lag_ratio=0.1), run_time=1.0)
    s.say("同一個 ", "$T/T_c$", "：二維的磁化留得久，到 ", "$T_c$", " 前才猛掉", wait=2.5)
    s.play(LaggedStart(*[FadeIn(c) for c in row_cells(tab, 4, 3)], lag_ratio=0.2), run_time=0.8)
    s.say("在 ", "$T_c$", " 附近量 ", "$\\ln M_s$", " 對 ", "$\\ln(1-T/T_c)$", " 的斜率", wait=2.5)
    s.play(FadeIn(concl), run_time=0.8)
    s.say("$\\beta=\\frac12$", "（平均場）、", "$\\beta=\\frac18$", "（二維 Ising）", wait=2.5)
    s.say("書上說 ", "$\\beta\\approx\\frac13$", f"：那是三維（Ising 約 {WC.BETA_CITED['3d_ising']}，只引用）", wait=3.0)
    s.clear_stage()


def part_degennes(s):
    head = mixed("de Gennes：", "$p=2\\nu$", "，", "$\\nu$", " 是一個假想磁鐵的指數", font_size=32)
    fs = 28
    rows = [[zh("漫步（自迴避）", font_size=fs, color=YELLOW), VectorizedPoint(),
             mixed("磁鐵（", "$n=0$", " 的 ", "$\\mathrm{O}(n)$", " 模型）", font_size=fs, color=YELLOW)],
            [mixed("步數 ", "$N$", font_size=fs), r"\leftrightarrow",
             mixed("$\\frac{1}{z_c-z}$", "（離臨界點多遠）", font_size=fs)],
            [r"\langle R^2\rangle\propto N^{2\nu}", r"\leftrightarrow", r"\xi^2\propto(z_c-z)^{-2\nu}"],
            [mixed("條數 ", "$c_N\\propto\\mu^N N^{\\gamma-1}$", font_size=fs), r"\leftrightarrow",
             r"\chi\propto(z_c-z)^{-\gamma}"],
            [zh("不准踩到自己", font_size=fs), r"\leftrightarrow", mixed("$n=0$", " 個分量", font_size=fs)]]
    tab = D.tex_table(rows, [4.8, 1.0, 5.6], row_h=0.92, font_size=34)
    place_in(head, -6.7, 6.7, 2.5, 3.1)
    place_in(tab, -6.7, 6.7, -2.8, 2.1)
    assert_in_frame(head, tab)
    assert_apart(head, tab)
    assert_table(tab, 3)
    s.say("de Gennes：", "$p=2\\nu$", "，", "$\\nu$", " 是一個假想磁鐵的指數")
    s.play(FadeIn(head), run_time=0.8)
    s.say("那個磁鐵的自旋有 ", "$n=0$", " 個分量：", "$\\mathrm{O}(n)$", " 模型（W902）", wait=2.5)
    s.play(FadeIn(row_cells(tab, 0, 3)), ShowCreation(tab[1]), run_time=0.8)
    s.play(FadeIn(row_cells(tab, 1, 3)), run_time=0.8)
    s.say("$z$", "：每一步的逸度（權重）；", "$z_c=1/\\mu$", " 是臨界點", wait=2.5)
    s.say("步數越多，就是越靠近臨界點", wait=2.5)
    s.play(FadeIn(row_cells(tab, 2, 3)), run_time=0.8)
    s.say("$\\langle R^2\\rangle$", " 對應關聯長度的平方 ", "$\\xi^2$", wait=2.5)
    s.play(FadeIn(row_cells(tab, 3, 3)), run_time=0.8)
    s.say("$c_N$", "：", "$N$", " 步自迴避漫步的條數（書上 §9.4 叫 ", "$\\Gamma_N$", "）", wait=2.5)
    s.say("條數 ", "$c_N$", " 對應磁化率 ", "$\\chi$", wait=2.5)
    s.play(FadeIn(row_cells(tab, 4, 3)), run_time=0.8)
    s.say("不准踩到自己，對應「", "$n=0$", " 個分量」", wait=2.5)
    s.say("漫步的 ", "$\\nu$", " 和磁鐵的 ", "$\\beta$", " 都是臨界指數", wait=2.5)
    s.say("第 9 章：為什麼「", "$n=0$", " 的磁鐵」會數自迴避漫步", wait=3.5)


class W901CriticalExponent(ChuScene):
    def construct(self):
        self.show_header("W901　漫步的指數就是臨界指數")
        check_math()
        part_intro(self)
        part_flory(self)
        part_magnet(self)
        part_degennes(self)
        self.wait(1.0)
