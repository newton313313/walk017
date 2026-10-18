"""
W701 漫步就是「路徑的和」（§7.1，pp. 167–171，(7.1)–(7.8)，Ex. 7.1）

問：從原點走 4 步回到原點，一維有幾條？二維呢？
  1. (7.3) 切成兩段：一維 N=4、M=2，中間點 r = −2, 0, 2：1·1 + 2·2 + 1·1 = 6（＝暴力列舉）。
  2. 二維：中心 4·4 = 16、四個 (±2,0),(0,±2) 各 1·1、四個 (±1,±1) 各 2·2 → 36 = C(4,2)²；
     M = N−1 退化成 (2.7)：3+3 = 6、4·9 = 36；切 2+2+2（(7.5)）：20、400。
  3. 高斯 (7.6)、Ex. 7.1：M1 = 3、M2 = 6 → 變異數 1+2 = 3；峰高 0.3989、0.2821 → 摺積 0.2303；筆誤框。
  4. 切成很多段 (7.7)(7.8)：折線變成光滑曲線，指數像「速度平方」。
畫面上的數字都由 utils.walk_saw 算出，check_math() 核對。
"""
import itertools
import math
import sys
import warnings
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_saw as SW
from utils import lattice_draw as L
from utils import walk_draw as D

IN_COL = TEAL_B
OUT_COL = ORANGE
MID_COL = YELLOW
G_COLS = (TEAL_B, BLUE_B, ORANGE)


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = dict(
            t1=SW.split_terms(4, 2, 0, 1), brute1=SW.count_brute(4, 0, 1),
            t2=SW.split_terms(4, 2, (0, 0), 2), brute2=SW.count_brute(4, (0, 0), 2),
            e1=SW.split_terms(4, 3, 0, 1), e2=SW.split_terms(4, 3, (0, 0), 2),
            m1=SW.multi_split(6, (2, 2, 2), 0, 1), m2=SW.multi_split(6, (2, 2, 2), (0, 0), 2),
            peaks=tuple(f4(SW.gauss1(M, 0.0)) for M in (3, 6, 9)),
            conv0=f4(SW.gauss_conv(3, 6, 0.0)), conv15=float(SW.gauss_conv(3, 6, 1.5)),
            ratio=f4(SW.conv_ratio(3, 6)), ratio_book=float(SW.conv_ratio(3, 6, book=True)),
            var=tuple(SW.lattice_var(M) for M in (3, 6, 9)),
        )
    return nb


def classify2(t2):
    """二維 N=4、M=2 的中間點分三類：中心、軸上、對角。"""
    cen = [t for t in t2 if t[0] == (0, 0)]
    axis = [t for t in t2 if t[0] != (0, 0) and 0 in t[0]]
    diag = [t for t in t2 if 0 not in t[0]]
    return cen, axis, diag


def check_math():
    nb = numbers()
    t1 = nb["t1"]
    assert t1 == [(-2, 1, 1), (0, 2, 2), (2, 1, 1)] == SW.EX_SPLIT["1d_terms"]
    assert sum(a * b for _, a, b in t1) == 6 == nb["brute1"] == SW.EX_SPLIT["1d_total"]
    assert 2 ** 4 == 16 and len(one_d_walks(4)) == 16
    assert sum(1 for p in one_d_walks(4) if p[-1] == 0) == 6
    for r, a, b in t1:                                            # 畫出來的進／出路線條數
        assert len(in_paths(r)) == a and len(out_paths(r)) == b
    cen, axis, diag = classify2(nb["t2"])
    assert len(cen) == 1 and cen[0][1] * cen[0][2] == 16 == SW.EX_SPLIT["2d_center"]
    assert len(axis) == 4 and all(a * b == 1 for _, a, b in axis) and SW.EX_SPLIT["2d_axis"] == 4
    assert len(diag) == 4 and all(a * b == 4 for _, a, b in diag) and SW.EX_SPLIT["2d_diag"] == 16
    assert cen[0][1] == 4 and all(a == 1 for _, a, _ in axis) and all(a == 2 for _, a, _ in diag)
    assert 16 + 4 * 1 + 4 * 4 == 36 == nb["brute2"] == math.comb(4, 2) ** 2 == SW.EX_SPLIT["2d_total"]
    assert [(a, b) for _, a, b in nb["e1"]] == [(3, 1), (3, 1)] and SW.EX_SPLIT["eq27_1d"] == (2, 3)
    assert len(nb["e2"]) == 4 and all((a, b) == (9, 1) for _, a, b in nb["e2"]) and SW.EX_SPLIT["eq27_2d"] == (4, 9)
    assert 3 + 3 == 6 and 4 * 9 == 36
    assert nb["m1"] == 20 == SW.EX_SPLIT["split222_1d"] and nb["m2"] == 400 == SW.EX_SPLIT["split222_2d"]
    assert SW.EX71["M1"] == 3 and SW.EX71["M2"] == 6 and SW.EX71["var"] == (1, 2, 3)
    assert nb["var"] == (1, 2, 3) and 1 + 2 == 3
    assert nb["peaks"] == ("0.3989", "0.2821", "0.2303") and nb["conv0"] == "0.2303"
    assert abs(nb["conv15"] - float(SW.gauss1(9, 1.5))) < 1e-6          # 摺積還是高斯（變異數 3）
    assert nb["ratio"] == "1.0000"
    assert f4(nb["ratio_book"]) == "0.5774" and f4(nb["ratio_book"] ** 3) == "0.1925" == f4(3 ** -1.5)


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def one_d_walks(n):
    """一維 n 步的全部走法：位置序列（從 0 開始）。"""
    out = []
    for steps in itertools.product((1, -1), repeat=n):
        out.append(tuple(itertools.accumulate((0,) + steps)))
    return out


def in_paths(r):
    return [p for p in one_d_walks(2) if p[-1] == r]


def out_paths(r):
    return [tuple(r + x for x in p) for p in one_d_walks(2) if r + p[-1] == 0]


def poly(points, color, width=4):
    return VMobject().set_points_as_corners(points).set_stroke(color, width)


def mid_diagram(r, center, ux=0.55, uy=0.3):
    """一維 4 步、第 2 步在 r 的全部路線：走法圖（橫軸步數、縱軸位置）。"""
    org = center - 2 * ux * RIGHT
    dots = VGroup(*[Dot(D.path_point(k, x, ux, uy, org), radius=0.028).set_fill(GREY_D, 1)
                    for k in range(5) for x in range(-2, 3)])
    zero = DashedLine(D.path_point(0, 0, ux, uy, org), D.path_point(4, 0, ux, uy, org),
                      dash_length=0.06).set_stroke(GREY_C, 1.2)
    ins = VGroup(*[poly([D.path_point(k, x, ux, uy, org) for k, x in enumerate(p)], IN_COL)
                   for p in in_paths(r)])
    outs = VGroup(*[poly([D.path_point(k + 2, x, ux, uy, org) for k, x in enumerate(p)], OUT_COL)
                    for p in out_paths(r)])
    mid = Dot(D.path_point(2, r, ux, uy, org), radius=0.09).set_fill(MID_COL, 1)
    ends = VGroup(Dot(D.path_point(0, 0, ux, uy, org), radius=0.06).set_fill(WHITE, 1),
                  Dot(D.path_point(4, 0, ux, uy, org), radius=0.06).set_fill(WHITE, 1))
    return VGroup(VGroup(dots, zero), ins, outs, mid, ends)


def smooth_pt(t, x0=-4.6, x1=4.6, y0=-0.3):
    return np.array([x0 + (x1 - x0) * t, y0 + 0.7 * math.sin(2 * PI * 1.2 * t) + 0.25 * math.sin(7.0 * t), 0])


def chord_path(n, color=YELLOW, dots=True):
    pts = [smooth_pt(i / n) for i in range(n + 1)]
    line = VMobject().set_points_as_corners(pts).set_stroke(color, 4 if n < 40 else 3)
    g = VGroup(line)
    if dots:
        g.add(VGroup(*[Dot(p, radius=0.06).set_fill(WHITE, 1) for p in pts]))
    return g


# ---------------------------------------------------------------- 各段
def part_question(s):
    nl = D.number_line(4, unit=0.55)
    nl.shift(np.array([-3.4, 0.6, 0]))
    zero = D.line_point(nl, 0).copy()
    w1 = Dot(zero, radius=0.11).set_fill(YELLOW, 1)
    arcs = VGroup()
    pos = [0, 1, 2, 1, 0]
    for k in range(4):
        a = zero + pos[k] * 0.55 * RIGHT
        b = zero + pos[k + 1] * 0.55 * RIGHT
        ang = -0.6 * PI if b[0] > a[0] else 0.95 * PI          # 去程矮、回程高，都在數線上方
        arcs.add(ArcBetweenPoints(a, b, angle=ang).set_stroke(YELLOW, 3))
    lab1 = at(zh("一維", font_size=30, color=GREY_A), -3.4, 2.0)
    org2 = np.array([3.4, 0.4, 0])
    dots2 = L.lattice_dots(2, cell=0.6, origin=org2, radius=0.05)
    tr = L.trail([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)], cell=0.6, origin=org2, color=YELLOW, stroke_width=4)
    lab2 = at(zh("二維", font_size=30, color=GREY_A), 3.4, 2.0)
    q = at(zh("？", font_size=60, color=YELLOW), 0, 0.5)
    assert_in_frame(nl, arcs, lab1, dots2, tr, lab2, q)
    for a, b in ((nl, dots2), (nl, q), (q, dots2), (lab1, arcs), (lab2, dots2)):
        assert_apart(a, b)
    s.play(FadeIn(nl), FadeIn(lab1), FadeIn(dots2), FadeIn(lab2), FadeIn(w1), run_time=0.8)
    s.say("從原點走 4 步，再回到原點", wait=0.3)
    s.play(LaggedStart(*[ShowCreation(a) for a in arcs], lag_ratio=0.6), run_time=1.6)
    s.play(ShowCreation(tr[0]), FadeIn(tr[1:]), run_time=1.4)
    s.say("一維有幾條？二維呢？先猜一猜", wait=0.3)
    s.play(FadeIn(q, scale=1.3), run_time=0.6)
    s.wait(2.0)
    s.say("這一幕的方法：把路線切成兩段再數", wait=2.4)
    s.clear_stage()


def part_split_1d(s):
    nb = numbers()
    eq = at(tx(r"C(N;\vec x,\vec y)=\sum_{\vec r}C(M;\vec x,\vec r)\,C(N-M;\vec r,\vec y)", 36), -0.6, 2.55)
    tag = tx(r"(7.3)", 26, GREY_B).next_to(eq, RIGHT, buff=0.5)
    nl = D.number_line(4, unit=0.55, label_every=1, font_size=20)
    nl.shift(1.55 * UP)
    mids = VGroup(*[Dot(D.line_point(nl, r), radius=0.09).set_fill(MID_COL, 1)
                    for r, _, _ in nb["t1"]])
    xs = (-4.4, 0.0, 4.4)
    diags, labs = VGroup(), VGroup()
    for (r, a, b), x in zip(nb["t1"], xs):
        diags.add(mid_diagram(r, np.array([x, -0.1, 0])))
        labs.add(at(tx(r"r=%d:\ %d\cdot%d=%d" % (r, a, b, a * b), 30), x, -1.2))
    terms = "+".join(r"%d\cdot%d" % (a, b) for _, a, b in nb["t1"])
    tot = sum(a * b for _, a, b in nb["t1"])
    total = at(mixed("$%s=%d$" % (terms, tot), "＝暴力列舉 16 條中回到原點的 %d 條" % nb["brute1"],
                     font_size=30), 0, -2.25)
    total[1].set_color(GREY_A)
    leg = at(VGroup(Line(ORIGIN, 0.5 * RIGHT).set_stroke(IN_COL, 4), zh("進來", font_size=22, color=IN_COL),
                    Line(ORIGIN, 0.5 * RIGHT).set_stroke(OUT_COL, 4), zh("出去", font_size=22, color=OUT_COL))
             .arrange(RIGHT, buff=0.15), 5.3, 1.55)
    for m in (eq, total):
        fit(m, 12.0)
    assert_in_frame(eq, tag, nl, mids, diags, labs, total, leg)
    for a, b in ((eq, nl), (nl, diags), (diags, labs), (labs, total), (nl, leg), (eq, leg), (tag, leg)):
        assert_apart(a, b)
    for i in range(3):
        for j in range(i + 1, 3):
            assert_apart(diags[i], diags[j], 0.3)
    s.play(Write(eq), FadeIn(tag), run_time=1.4)
    s.say("先走 ", "$M$", " 步到中間點 ", "$\\vec r$", "，再走 ", "$N-M$", " 步", wait=2.6)
    s.say("像轉機：經過哪個城市都可以，全部加起來", wait=2.8)
    s.play(FadeIn(nl), run_time=0.6)
    s.say("一維、", "$N=4$", "、", "$M=2$", "：第 2 步只能在 ", "$r=-2,0,2$", wait=0.3)
    s.play(LaggedStart(*[FadeIn(m, scale=1.5) for m in mids], lag_ratio=0.3), FadeIn(leg), run_time=1.0)
    s.wait(1.8)
    for i in range(3):
        g = diags[i]
        s.play(FadeIn(g[0]), FadeIn(g[4]), FadeIn(g[3]), run_time=0.5)
        s.play(ShowCreation(g[1]), run_time=0.6)
        s.play(ShowCreation(g[2]), run_time=0.6)
        s.play(FadeIn(labs[i]), run_time=0.4)
        if i == 1:
            s.say("到 0 有 2 種走法，從 0 回來也有 2 種：", "$2\\cdot2$", wait=2.4)
    s.say("進來的走法數乘出去的走法數，再對中間點加總", wait=2.6)
    s.play(FadeIn(total), run_time=0.8)
    s.say("一維 4 步回原點：6 條，和暴力列舉一樣", wait=2.8)
    s.clear_stage()


def part_split_2d(s):
    nb = numbers()
    cen, axis, diag = classify2(nb["t2"])
    org = np.array([-4.1, -0.1, 0])
    cell = 0.8
    vals = {r: a for r, a, _ in nb["t2"]}
    grid = L.count_grid(vals, cell=cell, origin=org, font_size=34, zero_dots=True, n=2)
    nums = {m.lat_xy: m for m in grid if isinstance(m, Tex)}
    start = Circle(radius=0.3).set_stroke(GREY_B, 2).move_to(org)
    head = at(mixed("$C(2;\\vec r)$", "：走 2 步到 ", "$\\vec r$", " 的走法數", font_size=26), -4.1, 2.45)
    fit(head, 4.6)
    rows = [
        mixed("中心：", "$%d\\cdot%d=%d$" % (cen[0][1], cen[0][2], cen[0][1] * cen[0][2]), font_size=30),
        mixed("四個 ", "$(\\pm2,0),(0,\\pm2)$", "：各 ", "$1\\cdot1$", font_size=30),
        mixed("四個 ", "$(\\pm1,\\pm1)$", "：各 ", "$2\\cdot2$", font_size=30),
        tx(r"16+4\cdot1+4\cdot4=%d=\binom42^2" % nb["brute2"], 34),
        mixed("$M=N-1$", "（(2.7)）：一維 ", "$3+3=6$", "、二維 ", "$4\\cdot9=36$", font_size=28),
        mixed("切 ", "$2+2+2$", "（(7.5)，", "$N=6$", "）：一維 %d、二維 %d" % (nb["m1"], nb["m2"]), font_size=28),
    ]
    for m, y in zip(rows, (2.45, 1.65, 0.85, -0.05, -1.05, -2.0)):
        at(m, 2.6, y)
        fit(m, 7.0)
    rows[3].set_color(YELLOW)
    allr = VGroup(*rows)
    assert_in_frame(grid, head, allr, start)
    for m in rows:
        assert_apart(grid, m)
        assert_apart(head, m)
    for i in range(5):
        assert_apart(rows[i], rows[i + 1])
    assert_apart(head, grid)
    s.play(FadeIn(head), FadeIn(grid), ShowCreation(start), run_time=1.0)
    s.say("二維：先數走 2 步到每個中間點的走法數", wait=2.6)
    c_num = nums[(0, 0)]
    s.play(c_num.animate.set_color(YELLOW), FadeIn(rows[0]), run_time=0.7)
    s.say("中心：去 4 種、回 4 種，", "$4\\cdot4=16$", wait=2.4)
    ax_nums = [nums[t[0]] for t in axis]
    s.play(*[m.animate.set_color(TEAL_B) for m in ax_nums], FadeIn(rows[1]), run_time=0.7)
    s.say("四個軸上的點：各 ", "$1\\cdot1$", wait=2.2)
    dg_nums = [nums[t[0]] for t in diag]
    s.play(*[m.animate.set_color(ORANGE) for m in dg_nums], FadeIn(rows[2]), run_time=0.7)
    s.say("四個對角的點：各 ", "$2\\cdot2$", wait=2.2)
    s.play(Write(rows[3]), run_time=1.0)
    s.say("加起來 36：正好是 ", "$\\binom42^2$", "，兩枚硬幣各自回原點", wait=3.0)
    s.play(FadeIn(rows[4]), run_time=0.6)
    s.say("切在最後一步（", "$M=N-1$", "）：就是第 2 章的 (2.7)", wait=2.8)
    s.say("一維 2 個鄰居各 3 條、二維 4 個鄰居各 9 條", wait=2.8)
    s.play(FadeIn(rows[5]), run_time=0.6)
    s.say("切成三段也一樣：每個中間點都要加總", wait=2.8)
    s.clear_stage()


def part_gauss(s):
    nb = numbers()
    eq = at(tx(r"C(M;\vec x,\vec r)\approx\Big(\frac{3}{2\pi Ma^2}\Big)^{3/2}"
               r"e^{-3|\vec x-\vec r|^2/2Ma^2}", 34), -0.6, 2.55)
    tag = tx(r"(7.6)", 26, GREY_B).next_to(eq, RIGHT, buff=0.5)
    ax = D.plot_axes((-6, 6, 1), (0, 0.45, 0.1), width=5.6, height=3.5, x_ticks=(-6, -3, 0, 3, 6),
                     y_ticks=(0.1, 0.2, 0.3, 0.4), y_tick_labels=("0.1", "0.2", "0.3", "0.4"))
    ax.move_to(np.array([-3.55, -0.45, 0]))
    graphs = VGroup(*[ax.get_graph(lambda x, M=M: float(SW.gauss1(M, x)), x_range=(-6, 6, 0.05))
                      .set_stroke(c, 4) for M, c in zip((3, 6, 9), G_COLS)])
    pk = nb["peaks"]
    r1 = mixed("Ex. 7.1：三維，", "$M_1=3$", "、", "$M_2=6$", "（", "$a=1$", "）", font_size=28)
    legs = VGroup()
    for c, tex, txt in ((G_COLS[0], "$M_1=3$", "：峰高 " + pk[0]), (G_COLS[1], "$M_2=6$", "：峰高 " + pk[1]),
                        (G_COLS[2], "$M_1+M_2=9$", "：峰高 " + pk[2])):
        legs.add(VGroup(Line(ORIGIN, 0.5 * RIGHT).set_stroke(c, 5), mixed(tex, txt, font_size=28))
                 .arrange(RIGHT, buff=0.2))
    r5 = mixed("一個分量的變異數：", "$1+2=3$", font_size=28, color=YELLOW)
    for m, y in ((r1, 1.75), (legs[0], 1.1), (legs[1], 0.55), (legs[2], 0.0), (r5, -0.65)):
        at(m, 3.55, y)
        fit(m, 5.9)
    for m in legs:
        m.align_to(legs[0], LEFT)
    box = D.erratum_box("p. 169，(7.6)、Ex. 7.1", r"\big(2\pi M a^2\big)^{-3/2}",
                        r"\big(2\pi M a^2/3\big)^{-3/2}", width=5.9)
    at(box, 3.55, -2.05)
    fit(eq, 11.6)
    assert_in_frame(eq, tag, ax, graphs, r1, legs, r5, box)
    for a, b in ((eq, r1), (eq, ax), (ax, r1), (ax, legs), (ax, r5), (ax, box), (r1, legs[0]), (legs[0], legs[1]),
                 (legs[1], legs[2]), (legs[2], r5), (r5, box)):
        assert_apart(a, b)
    s.play(Write(eq), FadeIn(tag), run_time=1.4)
    s.say("每一段很長時，走法數接近高斯（三維、格距 ", "$a$", "）", wait=2.8)
    s.play(FadeIn(ax), FadeIn(r1), run_time=0.8)
    s.say("Ex. 7.1：先走 3 步、再走 6 步，看一個直角分量", wait=0.3)
    s.play(ShowCreation(graphs[0]), FadeIn(legs[0]), run_time=1.0)
    s.play(ShowCreation(graphs[1]), FadeIn(legs[1]), run_time=1.0)
    s.say("一個分量的變異數：", "$M/3$", "，也就是 1 和 2", wait=2.6)
    s.play(ShowCreation(graphs[2]), FadeIn(legs[2]), run_time=1.2)
    s.say("對中間點積分（摺積）：峰高 %s" % pk[2], wait=2.6)
    s.play(FadeIn(r5), run_time=0.6)
    s.say("高斯的摺積還是高斯：變異數相加，", "$1+2=3$", wait=3.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上的前因子少了一個 ", "$3$", wait=2.4)
    s.say("照書上，每摺積一次多出 ", "$3^{-3/2}=0.1925$", "（一個分量 0.5774）", wait=3.2)
    s.clear_stage()


def part_continuum(s):
    e7 = at(tx(r"\exp\Big[-\frac{3}{2a^2}\sum_i\frac{|\Delta\vec r_i|^2}{\Delta N_i}\Big]", 34), -2.6, 2.6)
    t7 = tx(r"(7.7)", 26, GREY_B).next_to(e7, RIGHT, buff=0.4)
    e8 = at(tx(r"C\propto\prod_s\int d\vec r(s)\,\exp\Big[-\frac{3}{2a^2}\int_0^N"
               r"\Big|\frac{d\vec r}{ds}\Big|^2ds\Big]", 34), 1.6, 1.45)
    t8 = tx(r"(7.8)", 26, GREY_B).next_to(e8, RIGHT, buff=0.4)
    paths = [chord_path(3), chord_path(6), chord_path(12), chord_path(60, YELLOW, dots=False)]
    p3 = paths[0][1]
    a, b = p3[0].get_center().copy(), p3[1].get_center().copy()
    brace = Line(a, b).set_stroke(TEAL_B, 6)
    dlab = tx(r"\Delta\vec r_i,\ \Delta N_i", 28, TEAL_B).next_to(brace.get_center(), DOWN + RIGHT, buff=0.15)
    n1 = at(mixed("像「速度平方」：", "$\\Big|\\frac{d\\vec r}{ds}\\Big|^2$", "，動能", font_size=28,
                  color=YELLOW), -3.2, -1.9)
    n2 = at(mixed("注意：", "$\\Delta N_i\\gg1$", " 高斯才成立；", "$d\\vec r/ds$", " 只是記號", font_size=26,
                  color=GREY_A), 2.9, -2.6)
    for m in (e8, n1, n2):
        fit(m, 7.4)
    assert_in_frame(e7, t7, e8, t8, *paths, dlab, n1, n2)
    for p in paths:
        for m in (e7, t7, e8, t8, n1, n2):
            assert_apart(p, m)
    for x, y in ((e7, e8), (t7, e8), (n1, n2), (dlab, n1), (dlab, n2)):
        assert_apart(x, y)
    s.play(Write(e7), FadeIn(t7), run_time=1.2)
    s.say("切成很多段：每段一個高斯，指數加起來", wait=0.3)
    s.play(ShowCreation(paths[0][0]), FadeIn(paths[0][1]), run_time=1.2)
    s.play(ShowCreation(brace), FadeIn(dlab), run_time=0.8)
    s.say("每一段的位移 ", "$\\Delta\\vec r_i$", "、步數 ", "$\\Delta N_i$", wait=2.4)
    s.play(FadeOut(brace), FadeOut(dlab), run_time=0.4)
    cur = paths[0]
    for nxt in paths[1:]:
        s.play(Transform(cur, nxt), run_time=1.2)
        s.wait(0.6)
    s.say("段越切越細：折線變成光滑的曲線 ", "$\\vec r(s)$", wait=2.4)
    s.play(Write(e8), FadeIn(t8), run_time=1.4)
    s.say("對每個中間點積分：所有路徑的和（路徑積分）", wait=3.0)
    s.play(FadeIn(n1), run_time=0.6)
    s.say("指數裡是「速度平方」：像動能", wait=2.6)
    s.play(FadeIn(n2), run_time=0.6)
    s.say("但每段其實要很長，高斯才成立；", "$d\\vec r/ds$", " 只是記號", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("漫步的條數＝所有路徑的和", wait=2.4)
    s.say("切段、求和、再切下去，就是路徑積分", wait=3.0)


class W701PathSum(ChuScene):
    def construct(self):
        self.show_header("W701　漫步就是「路徑的和」")
        check_math()
        part_question(self)
        part_split_1d(self)
        part_split_2d(self)
        part_gauss(self)
        part_continuum(self)
        part_close(self)
        self.wait(1.0)
