"""
W702 自迴避漫步：會記得走過的地方（§7.2–7.2.1，pp. 171–172，(7.9)(7.10)）

問：如果醉漢不准踩到自己走過的點，10 步之後離起點會比較遠還是比較近？
  1. 三種自迴避：嚴格（踩到就停）、溫和（越靠近越可能停）、「真」自迴避（偏向不踩，但一直走）；本章用嚴格的。
     權重 (7.9) Π_{i>j} e^{−v(|w_i − w_j|)}：重合時 v = ∞ → 權重 0 → 那條路徑被刪掉。
  2. 逐條畫：N=2：16 條中 12 條（4 條原地折返）；N=3：64 條中 36 條（第一步向東的 9 條 ×4）；N=4：256 條中 100 條。
  3. 精確列舉 N = 1..14：c_N、c_N/4^N、⟨R²⟩（普通漫步 = N）。
  4. 對比（種子 702）：自迴避的 10 步終點 (−3,3)、R² = 18；普通漫步終點 (1,−1)、R² = 2（中間踩回 3 次）。
畫面上的數字都由 utils.walk_saw 算出，check_math() 核對。
"""
import sys
from fractions import Fraction
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

NS = (1, 2, 3, 4, 6, 8, 10, 12, 14)
SAW_COL = YELLOW
RW_COL = BLUE_B
BAD_COL = RED
OK_COL = "#7CE38B"
BACKTRACK = (((0, 0), (1, 0), (0, 0)), ((0, 0), (0, 1), (0, 0)), ((0, 0), (-1, 0), (0, 0)), ((0, 0), (0, -1), (0, 0)))


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def revisits(pts):
    seen, n = set(), 0
    for p in pts:
        n += p in seen
        seen.add(p)
    return n


@lru_cache(maxsize=None)
def numbers():
    saw = [tuple(p) for p in SW.sample_saw(10)]
    rw = [tuple(p) for p in SW.sample_rw(10)]
    return dict(
        l2=[tuple(map(tuple, w)) for w in SW.saw_list(2)], l3=[tuple(map(tuple, w)) for w in SW.saw_list(3)],
        n4=len(SW.saw_list(4)),
        c=[SW.saw_c(N) for N in NS], frac=[f4(SW.saw_frac(N)) for N in NS], R2=[f4(SW.saw_R2(N)) for N in NS],
        R2_raw=[float(SW.saw_R2(N)) for N in NS],
        saw=saw, rw=rw,
    )


def check_math():
    nb = numbers()
    assert len(nb["l2"]) == 12 and 4 ** 2 == 16 and 16 - len(BACKTRACK) == 12
    assert len(nb["l3"]) == 36 and 4 ** 3 == 64 and 9 * 4 == 36
    assert all(w[1] == (1, 0) for w in nb["l3"][:9]) and sum(w[1] == (1, 0) for w in nb["l3"]) == 9
    assert nb["n4"] == 100 and 4 ** 4 == 256
    for w in nb["l2"] + nb["l3"]:
        assert len(set(w)) == len(w)                               # N+1 個點兩兩不同
    for w in BACKTRACK:
        assert len(set(w)) < len(w)
    assert nb["c"] == [4, 12, 36, 100, 780, 5916, 44100, 324932, 2374444]
    assert nb["frac"] == ["1.0000", "0.7500", "0.5625", "0.3906", "0.1904", "0.0903", "0.0421", "0.0194", "0.0088"]
    assert nb["R2"] == ["1.0000", "2.6667", "4.5556", "7.0400", "12.5744", "19.0128", "26.2425", "34.1870",
                        "42.7864"]
    assert all(SW.saw_frac(N) == Fraction(SW.saw_c(N), 4 ** N) for N in NS)
    assert all(r > N for r, N in zip(nb["R2_raw"][1:], NS[1:]))                # 比普通漫步的 N 大
    saw, rw = nb["saw"], nb["rw"]
    assert len(saw) == 11 and len(set(saw)) == 11 and saw[-1] == (-3, 3) and r2(saw[-1]) == 18
    assert len(rw) == 11 and rw[-1] == (1, -1) and r2(rw[-1]) == 2 and revisits(rw) == 3
    assert SW.SAMPLE_SEEDS == {"saw": 702, "rw": 702}


def r2(p):
    return p[0] ** 2 + p[1] ** 2


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def dots_box(x0, x1, y0, y1, cell, org, radius=0.03, color=GREY_D):
    return VGroup(*[Dot(L.lat_point(x, y, cell, org), radius=radius).set_fill(color, 1)
                    for y in range(y1, y0 - 1, -1) for x in range(x0, x1 + 1)])


def cross(p, size=0.12, color=BAD_COL, width=4):
    return VGroup(Line(p + size * UL, p + size * DR), Line(p + size * UR, p + size * DL)).set_stroke(color, width)


def tile(points, box, cell, center, color=SAW_COL, bad=False, frame=1.15):
    """小路線圖：淡格點＋一條折線＋起點綠點（box = (x0, x1, y0, y1)）。"""
    x0, x1, y0, y1 = box
    org = center - cell * np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0])
    dots = dots_box(x0, x1, y0, y1, cell, org)
    pts = [L.lat_point(x, y, cell, org) for x, y in points]
    if bad:                                                       # 原地折返：去回畫成兩條平行線
        a, b = pts[0], pts[1]
        off = 0.035 * rotate_vector(normalize(b - a), PI / 2)
        line = VGroup(Line(a + off, b + off), Line(b - off, a - off)).set_stroke(GREY_B, 3)
    else:
        line = VMobject().set_points_as_corners(pts).set_stroke(color, 3)
    start = Dot(pts[0], radius=0.06).set_fill(OK_COL, 1)
    box_m = Square(side_length=frame).set_stroke(GREY_D, 1.2).move_to(center)
    g = VGroup(box_m, dots, line, start)
    if bad:
        g.add(cross(center, size=0.4, width=5))
    return g


def icon_strict(center, cell=0.45):
    org = center - cell * np.array([0.5, 1.0, 0])
    pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
    dots = dots_box(-1, 2, -1, 2, cell, org)
    tr = L.trail(pts, cell, org, SAW_COL, 4, ends=False)
    p0, p1 = L.lat_point(0, 1, cell, org), L.lat_point(0, 0, cell, org)
    bad = DashedLine(p0, p1, dash_length=0.07).set_stroke(BAD_COL, 3)
    return VGroup(dots, tr, bad, cross(p1, 0.13), Dot(p1, radius=0.06).set_fill(OK_COL, 1))


def icon_soft(center, cell=0.45):
    org = center - cell * np.array([0.5, 1.0, 0])
    pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
    dots = dots_box(-1, 2, -1, 2, cell, org)
    halos = VGroup(*[Circle(radius=0.32).set_stroke(width=0).set_fill(BAD_COL, 0.22)
                     .move_to(L.lat_point(x, y, cell, org)) for x, y in pts[:-1]])
    tr = L.trail(pts, cell, org, SAW_COL, 4, ends=False)
    return VGroup(dots, halos, tr, Dot(L.lat_point(0, 0, cell, org), radius=0.06).set_fill(OK_COL, 1))


def icon_true(center, cell=0.45):
    org = center - cell * np.array([0.5, 1.0, 0])
    pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
    dots = dots_box(-1, 2, -1, 2, cell, org)
    tr = L.trail(pts, cell, org, SAW_COL, 4, ends=False)
    here = L.lat_point(0, 1, cell, org)
    good = VGroup(*[Arrow(here, L.lat_point(x, y, cell, org), buff=0.06).set_stroke(OK_COL, 5)
                    for x, y in ((-1, 1), (0, 2))])
    weak = Arrow(here, L.lat_point(0, 0, cell, org), buff=0.06).set_stroke(GREY_B, 2)
    return VGroup(dots, tr, good, weak, Dot(L.lat_point(0, 0, cell, org), radius=0.06).set_fill(OK_COL, 1))


# ---------------------------------------------------------------- 各段
def part_question(s):
    cell = 0.6
    org = np.array([-0.3, -0.6, 0])
    dots = dots_box(-2, 3, -1, 3, cell, org, radius=0.045, color=GREY_C)
    pts = [(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)]
    tr = L.trail(pts, cell, org, SAW_COL, 5, ends=False)
    start = Dot(L.lat_point(0, 0, cell, org), radius=0.09).set_fill(OK_COL, 1)
    p0, p1 = L.lat_point(0, 1, cell, org), L.lat_point(0, 0, cell, org)
    bad = DashedLine(p0, p1, dash_length=0.08).set_stroke(BAD_COL, 4)
    x = cross(p1, 0.17, BAD_COL, 6)
    q = at(zh("？", font_size=64, color=YELLOW), 3.8, 0.4)
    assert_in_frame(dots, tr, q)
    assert_apart(dots, q)
    s.play(FadeIn(dots), FadeIn(start), run_time=0.6)
    s.say("醉漢不准踩到自己走過的點", wait=0.3)
    s.play(ShowCreation(tr), run_time=1.6)
    s.play(ShowCreation(bad), run_time=0.5)
    s.play(ShowCreation(x), run_time=0.4)
    s.wait(1.2)
    s.say("10 步之後，離起點會比較遠，還是比較近？", wait=0.3)
    s.play(FadeIn(q, scale=1.3), run_time=0.5)
    s.wait(2.6)
    s.clear_stage()


def part_kinds(s):
    xs = (-4.5, 0.0, 4.5)
    icons = VGroup(icon_strict(np.array([xs[0], 1.7, 0]), 0.6), icon_soft(np.array([xs[1], 1.7, 0]), 0.6),
                   icon_true(np.array([xs[2], 1.7, 0]), 0.6))
    labs = VGroup(*[at(m, x, 0.0) for m, x in zip(
        (mixed("嚴格：踩到就停", font_size=28), mixed("溫和：越靠近越可能停", font_size=28),
         mixed("「真」自迴避：偏向不踩", font_size=28)), xs)])
    sub = at(zh("但一直走", font_size=24, color=GREY_A), xs[2], -0.5)
    tag = at(zh("本章", font_size=24, color=BLACK), xs[0], -0.55)
    tag_bg = SurroundingRectangle(tag, buff=0.08).set_fill(YELLOW, 1).set_stroke(width=0)
    w = at(tx(r"\prod_{i>j}e^{-v(|\vec w_i-\vec w_j|)}", 40), -2.6, -1.75)
    wt = tx(r"(7.9)", 26, GREY_B).next_to(w, LEFT, buff=0.4)
    w2 = at(mixed("重合：", "$v=\\infty\\ \\Rightarrow\\ e^{-v}=0$", font_size=32), 3.0, -1.75)
    w2[1].set_color(BAD_COL)
    for m in labs:
        fit(m, 4.2)
    assert_in_frame(icons, labs, sub, tag_bg, w, wt, w2)
    for i in range(3):
        assert_apart(icons[i], labs[i])
        for j in range(i + 1, 3):
            assert_apart(icons[i], icons[j], 0.3)
            assert_apart(labs[i], labs[j], 0.2)
    for a, b in ((labs, w), (labs, w2), (w, w2), (sub, w2), (tag_bg, w), (wt, tag_bg), (sub, labs[2])):
        assert_apart(a, b)
    s.say("「自迴避」有三種說法", wait=0.3)
    for i in range(3):
        s.play(FadeIn(icons[i]), FadeIn(labs[i]), *([FadeIn(sub)] if i == 2 else []), run_time=0.7)
        s.wait(0.4)
    s.say("嚴格：踩到自己就停（那條路線不算）", wait=2.8)
    s.say("溫和：離走過的點越近，越可能停", wait=2.8)
    s.say("「真」自迴避：偏向不踩，但一直走下去", wait=2.8)
    s.play(FadeIn(tag_bg), FadeIn(tag), run_time=0.5)
    s.say("本章用嚴格的：", "$N+1$", " 個點兩兩不同", wait=2.6)
    s.play(Write(w), FadeIn(wt), run_time=1.2)
    s.say("每一對點乘一個權重 ", "$e^{-v}$", wait=2.4)
    s.play(FadeIn(w2), run_time=0.6)
    s.say("重合時 ", "$v=\\infty$", "：權重 0，那條路徑被刪掉", wait=2.8)
    s.clear_stage()


def part_enumerate(s):
    nb = numbers()
    cell = 0.22
    # N = 2：先看 4 條原地折返
    bad = VGroup(*[tile(w, (-2, 2, -2, 2), cell, np.array([x, 0.9, 0]), bad=True)
                   for w, x in zip(BACKTRACK, (-2.7, -0.9, 0.9, 2.7))])
    bad_lab = at(mixed("原地折返：", "$4$", " 條", font_size=30, color=GREY_A), 0, -0.6)
    assert_in_frame(bad, bad_lab)
    assert_apart(bad, bad_lab)
    s.say("一條條數：", "$N=2$", " 有 ", "$4^2=16$", " 條普通漫步", wait=0.3)
    s.play(LaggedStart(*[FadeIn(t) for t in bad], lag_ratio=0.2), run_time=1.2)
    s.play(FadeIn(bad_lab), run_time=0.5)
    s.say("其中 4 條原地折返：第 2 步踩回起點", wait=2.4)
    s.play(FadeOut(bad), FadeOut(bad_lab), run_time=0.6)
    xs = [-3.5 + 1.4 * i for i in range(6)]
    good = VGroup(*[tile(w, (-2, 2, -2, 2), cell, np.array([xs[i % 6], 1.75 - 1.4 * (i // 6), 0]))
                    for i, w in enumerate(nb["l2"])])
    lab2 = at(mixed("$N=2$", "：", "$4^2=16$", " 條中 ", "$%d$" % len(nb["l2"]), " 條", font_size=32), 0, -1.1)
    assert_in_frame(good, lab2)
    assert_apart(good, lab2)
    s.play(LaggedStart(*[FadeIn(t) for t in good], lag_ratio=0.08), run_time=1.6)
    s.play(FadeIn(lab2), run_time=0.5)
    s.say("剩下 12 條：", "$16-4=12$", wait=3.0)
    s.play(FadeOut(good), FadeOut(lab2), run_time=0.6)
    # N = 3：第一步向東的 9 條
    t3 = VGroup(*[tile(w, (-1, 3, -2, 2), cell, np.array([-4.6 + 1.4 * (i % 3), 1.65 - 1.4 * (i // 3), 0]),
                       frame=1.2) for i, w in enumerate(nb["l3"][:9])])
    r1 = at(mixed("$N=3$", "：第一步向東的 9 條", font_size=30), 2.6, 1.7)
    r2 = at(mixed("$\\times4$", "（四個方向）", "$=%d$" % len(nb["l3"]), font_size=30), 2.6, 0.9)
    r3 = at(mixed("$4^3=64$", " 條中 ", "$%d$" % len(nb["l3"]), " 條", font_size=30, color=YELLOW), 2.6, 0.1)
    r4 = at(mixed("$N=4$", "：", "$4^4=256$", " 條中 ", "$%d$" % nb["n4"], " 條", font_size=30, color=YELLOW),
            2.6, -1.1)
    rows = VGroup(r1, r2, r3, r4)
    assert_in_frame(t3, rows)
    for m in rows:
        assert_apart(t3, m)
    for a, b in ((r1, r2), (r2, r3), (r3, r4)):
        assert_apart(a, b)
    s.say("$N=3$", "：只畫第一步向東的", wait=0.3)
    s.play(LaggedStart(*[FadeIn(t) for t in t3], lag_ratio=0.1), FadeIn(r1), run_time=1.6)
    s.wait(1.2)
    s.play(FadeIn(r2), run_time=0.5)
    s.say("轉 90 度就得到其他三個方向：", "$9\\times4=36$", wait=2.6)
    s.play(FadeIn(r3), run_time=0.5)
    s.say("64 條只剩 36 條", wait=2.6)
    s.play(FadeIn(r4), run_time=0.5)
    s.say("$N=4$", "：256 條只剩 100 條，刪掉的越來越多", wait=2.8)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    hdr = [r"N"] + [str(N) for N in NS]
    rows = [hdr, [r"c_N"] + [str(c) for c in nb["c"]], [r"c_N/4^N"] + nb["frac"],
            [r"\langle R^2\rangle"] + nb["R2"]]
    tab = D.tex_table(rows, [1.55] + [1.25] * 9, row_h=0.48, font_size=24)
    tab.move_to(np.array([0, 2.12, 0]))
    ax = D.plot_axes((0, 15, 1), (0, 45, 5), width=7.6, height=2.9, x_ticks=(0, 2, 4, 6, 8, 10, 12, 14),
                     y_ticks=(10, 20, 30, 40), x_label="N", y_label=r"\langle R^2\rangle")
    ax.move_to(np.array([-1.6, -1.25, 0]))
    pts = VGroup(*[Dot(ax.c2p(N, r), radius=0.07).set_fill(SAW_COL, 1) for N, r in zip(NS, nb["R2_raw"])])
    line = DashedLine(ax.c2p(0, 0), ax.c2p(14.5, 14.5), dash_length=0.1).set_stroke(RW_COL, 3)
    leg1 = VGroup(Dot(radius=0.07).set_fill(SAW_COL, 1), mixed("自迴避的 ", "$\\langle R^2\\rangle$", font_size=26))
    leg1.arrange(RIGHT, buff=0.2)
    leg2 = VGroup(DashedLine(ORIGIN, 0.6 * RIGHT, dash_length=0.1).set_stroke(RW_COL, 3),
                  mixed("普通漫步 ", "$=N$", font_size=26)).arrange(RIGHT, buff=0.2)
    at(leg1, 4.6, -0.6)
    at(leg2, 4.6, -1.3)
    leg2.align_to(leg1, LEFT)
    rects = [SurroundingRectangle(VGroup(*[D.table_cell(tab, i, j, 10) for j in range(10)]), buff=0.04)
             .set_stroke(YELLOW, 2) for i in (1, 2, 3)]
    assert_in_frame(tab, ax, pts, leg1, leg2, *rects)
    for a, b in ((tab, ax), (ax, leg1), (ax, leg2), (leg1, leg2), (tab, leg1)):
        assert_apart(a, b)
    s.say("精確列舉：電腦一條條數（二維正方格子）", wait=0.3)
    s.play(FadeIn(tab[1]), FadeIn(VGroup(*[D.table_cell(tab, i, 0, 10) for i in range(4)])),
           FadeIn(VGroup(*[D.table_cell(tab, 0, j, 10) for j in range(1, 10)])), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(D.table_cell(tab, 1, j, 10)) for j in range(1, 10)], lag_ratio=0.12), run_time=1.4)
    s.play(ShowCreation(rects[0]), run_time=0.5)
    s.say("$c_N$", "：每多一步乘的倍數比 4 小", wait=2.8)
    s.say("$N=14$", "：有 2374444 條", wait=2.8)
    s.play(FadeOut(rects[0]), LaggedStart(*[FadeIn(D.table_cell(tab, 2, j, 10)) for j in range(1, 10)],
                                          lag_ratio=0.12), run_time=1.4)
    s.play(ShowCreation(rects[1]), run_time=0.5)
    s.say("存活的比例 ", "$c_N/4^N$", "：14 步只剩 0.0088", wait=2.8)
    s.play(FadeOut(rects[1]), LaggedStart(*[FadeIn(D.table_cell(tab, 3, j, 10)) for j in range(1, 10)],
                                          lag_ratio=0.12), run_time=1.4)
    s.play(ShowCreation(rects[2]), run_time=0.5)
    s.say("終點的平方距離 ", "$\\langle R^2\\rangle$", "：普通漫步剛好是 ", "$N$", wait=2.8)
    s.play(FadeOut(rects[2]), FadeIn(ax), run_time=0.7)
    s.play(ShowCreation(line), FadeIn(leg2), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(p, scale=1.5) for p in pts], lag_ratio=0.15), FadeIn(leg1), run_time=1.4)
    s.say("自迴避的點全在虛線上方：14 步是 42.7864，不是 14", wait=3.4)
    s.say("而且越走越往上翹：走得比普通漫步遠", wait=2.6)
    s.clear_stage()


def part_compare(s):
    nb = numbers()
    saw, rw = nb["saw"], nb["rw"]
    cell = 0.42

    def panel(pts, center):
        xs_, ys_ = [p[0] for p in pts], [p[1] for p in pts]
        box = (min(xs_) - 1, max(xs_) + 1, min(ys_) - 1, max(ys_) + 1)
        org = center - cell * np.array([(box[0] + box[1]) / 2, (box[2] + box[3]) / 2, 0])
        return org, dots_box(*box, cell, org, radius=0.035, color=GREY_C)

    o1, d1 = panel(saw, np.array([-3.4, 0.45, 0]))
    o2, d2 = panel(rw, np.array([3.4, 0.45, 0]))
    t1 = L.trail(saw, cell, o1, SAW_COL, 4)
    t2 = L.trail(rw, cell, o2, RW_COL, 4)
    seen, again = set(), {}
    for p in rw:
        if p in seen:
            again[p] = again.get(p, 0) + 1
        seen.add(p)
    rings = VGroup()
    for p, k in again.items():
        ring = Circle(radius=0.16).set_stroke(BAD_COL, 3).move_to(L.lat_point(*p, cell, o2))
        num = tx(r"\times%d" % k, 22, BAD_COL).move_to(ring.get_center() + 0.5 * cell * (RIGHT + DOWN))
        rings.add(VGroup(ring, num))                            # 標在格子中央，不壓到格點
    assert sum(again.values()) == revisits(rw)
    h1 = at(zh("自迴避（10 步）", font_size=28, color=SAW_COL), -3.4, 2.75)
    h2 = at(zh("普通漫步（10 步）", font_size=28, color=RW_COL), 3.4, 2.75)
    e1 = saw[-1]
    e2 = rw[-1]
    b1 = at(mixed("終點 ", "$(%d,%d)$" % e1, "：", "$R^2=%d$" % r2(e1), font_size=30), -3.4, -2.3)
    b2 = at(mixed("終點 ", "$(%d,%d)$" % e2, "：", "$R^2=%d$" % r2(e2), "（踩回 %d 次）" % revisits(rw),
                  font_size=30), 3.4, -2.3)
    assert_in_frame(d1, d2, t1, t2, h1, h2, b1, b2, rings)
    for a, b in ((d1, d2), (d1, h1), (d2, h2), (d1, b1), (d2, b2), (b1, b2), (h1, h2)):
        assert_apart(a, b)
    s.say("同一顆種子：自迴避和普通漫步各走 10 步", wait=0.3)
    s.play(FadeIn(d1), FadeIn(d2), FadeIn(h1), FadeIn(h2), run_time=0.7)
    s.play(ShowCreation(t1[0]), FadeIn(t1[1]), run_time=2.4)
    s.play(FadeIn(t1[2]), FadeIn(b1), run_time=0.6)
    s.say("自迴避：終點 ", "$(-3,3)$", "，", "$R^2=18$", wait=2.4)
    s.play(ShowCreation(t2[0]), FadeIn(t2[1]), run_time=2.4)
    s.play(FadeIn(t2[2]), FadeIn(rings), FadeIn(b2), run_time=0.6)
    s.say("普通漫步：中間踩回 3 次（紅圈），終點 ", "$R^2=2$", wait=2.8)
    s.say("同樣 10 步，自迴避走得遠多了", wait=2.6)
    s.say("不能回頭，只好往外走", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("自迴避漫步被自己的足跡推著走", wait=2.6)
    s.say("條數越來越少，卻走得越來越遠", wait=3.0)


class W702SelfAvoiding(ChuScene):
    def construct(self):
        self.show_header("W702　自迴避漫步：會記得走過的地方")
        check_math()
        part_question(self)
        part_kinds(self)
        part_enumerate(self)
        part_table(self)
        part_compare(self)
        part_close(self)
        self.wait(1.0)
