"""
W603 非球度：0 是圓球、1 是一條線（§6.1.4，pp. 133–135，(6.10)–(6.18)，Fig. 6.5，Ex. 6.1）

問：三維漫步平均起來有多扁？在 0（圓球）到 1（一條線）之間猜一個數。
  1. (6.10)(6.13)＋個別非球度；形狀例子：圓球 0、線 1、雪茄與扁盤 1/4、W602 小例子 1/4、Fig. 6.4 主例 0.6117。
  2. 揭曉 (6.12)(6.14)(6.15)：A_d = (4+2d)/(4+5d)；表 d = 2, 3, 4, 10, ∞。
  3. Fig. 6.5：三維格子 2 萬條 100 步的直方圖；平均非球度 0.5231（精確 0.5235，N→∞ 0.5263）；小於 0.5 占 69.3%。
  4. (6.16)–(6.18)：個別值的平均 A'_3 = 0.394274（模擬 0.3939）；筆誤框 (6.16)。
畫面上的數字都由 utils.walk_shape 算出，check_math() 核對。
"""
import math
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
from utils import walk_shape as SH
from utils import walk_draw as D

F = Fraction
SHAPES = [("圓球", (F(1), F(1), F(1))), ("線", (F(1), F(0), F(0))), ("雪茄", (F(4), F(1), F(1))),
          ("扁盤", (F(1), F(1), F(0)))]
DS = (2, 3, 4, 10)
A_COL = ORANGE           # 平均的比 A_d
AP_COL = "#7CE38B"       # 個別值的平均 A'_d
BAR_COL = BLUE_C
LOW_COL = TEAL_B


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return ("-" if q < 0 else "") + r"\frac{%d}{%d}" % (abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def numbers():
    se = SH.small_example()
    fs = SH.fig64_stats()
    st = SH.sim3d_stats()
    v = [float(x) for x in se["eig"]]
    return dict(shape_vals=[SH.asph(list(t)) for _, t in SHAPES], small=se["eig"], small_a=SH.asph(list(se["eig"])),
                fig64_eig=(f"{fs['eig'][0]:.1f}", f"{fs['eig'][1]:.1f}"), fig64_a=f4(fs["asph"]),
                Ad={d: SH.A_d(d) for d in DS}, Ad_dec={d: f4(SH.A_d(d)) for d in DS}, hist=[int(h) for h in st["hist"]],
                A_sim=f4(st["A_ratio"]), A_ex=f4(SH.A_exact(100, 3, "lattice")), A_inf=f4(SH.A_d(3)),
                lt=f"{100 * st['frac_lt_half']:.1f}", Ap={d: f"{SH.Aprime(d):.6f}" for d in (2, 3, 4)},
                Ap_sim=f4(st["A_mean"]), A_sim_raw=st["A_ratio"], Ap_sim_raw=st["A_mean"],
                book=(f"{SH.asph_book16(v, 4):.1f}", f"{SH.asph_book16([4 * x for x in v], 4):.1f}"),
                fixed=(SH.asph(list(se["eig"])), SH.asph([4 * x for x in se["eig"]])),
                walk64=SH.fig64())


def check_math():
    nb = numbers()
    assert nb["shape_vals"] == [0, 1, F(1, 4), F(1, 4)]
    assert nb["small"] == (F(3, 5), F(1, 5)) and nb["small_a"] == F(1, 4)
    assert nb["fig64_eig"] == ("258.8", "31.6") and nb["fig64_a"] == "0.6117"
    for t in (F(1), F(4), F(3, 5)):                                       # 與大小無關
        assert SH.asph([t * x for x in SHAPES[2][1]]) == F(1, 4)
    v = [F(1), F(2), F(5)]                                                 # (6.10)(6.13)：兩種算法
    lhs = sum((v[i] - v[j]) ** 2 for i in range(3) for j in range(i + 1, 3))
    assert lhs == 3 * sum(x * x for x in v) - sum(v) ** 2
    assert nb["Ad"] == {2: F(4, 7), 3: F(10, 19), 4: F(1, 2), 10: F(4, 9)}
    assert nb["Ad_dec"] == {2: "0.5714", 3: "0.5263", 4: "0.5000", 10: "0.4444"}
    assert abs(float(SH.A_d(10 ** 9)) - 0.4) < 1e-8
    assert nb["hist"] == [119, 616, 1183, 1631, 1890, 1795, 1808, 1621, 1669, 1529, 1505, 1305, 1065, 1000, 655, 389,
                          162, 56, 2, 0]
    assert sum(nb["hist"]) == 20000 and max(nb["hist"]) <= 2000 and int(np.argmax(nb["hist"])) == 4
    assert (nb["A_sim"], nb["A_ex"], nb["A_inf"]) == ("0.5231", "0.5235", "0.5263")
    assert nb["lt"] == "69.3" and sum(nb["hist"][:10]) / 20000 == SH.sim3d_stats()["frac_lt_half"]
    assert nb["Ap"] == {2: "0.396400", 3: "0.394274", 4: "0.393829"} and nb["Ap_sim"] == "0.3939"
    assert nb["book"] == ("0.2", "0.8") and nb["fixed"] == (F(1, 4), F(1, 4))
    w = nb["walk64"]
    assert len(w) == 1001


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def eq_label(tex_mob, label, buff=0.35):
    return VGroup(tex_mob, tx(r"\text{" + label + "}", 22, GREY_B)).arrange(RIGHT, buff=buff)


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def walk_icon(P, width, center, color=YELLOW, stroke=1.5):
    P = np.asarray(P, dtype=float)
    lo, hi = P.min(axis=0), P.max(axis=0)
    u = width / max(hi[0] - lo[0], 1e-9)
    c = (lo + hi) / 2
    pts = np.zeros((len(P), 3))
    pts[:, 0] = center[0] + u * (P[:, 0] - c[0])
    pts[:, 1] = center[1] + u * (P[:, 1] - c[1])
    return VMobject().set_points_as_corners(pts).set_stroke(color, stroke)


def shape_icon(k, center):
    c = np.array([center[0], center[1], 0])
    if k == 0:
        body = Circle(radius=0.45).set_stroke(WHITE, 2).set_fill(GREY_B, 0.3)
        return VGroup(body, Ellipse(width=0.9, height=0.25).set_stroke(GREY_B, 1.5)).move_to(c)
    if k == 1:
        return Line(0.75 * LEFT, 0.75 * RIGHT).set_stroke(WHITE, 6).move_to(c)
    if k == 2:
        return Ellipse(width=1.5, height=0.75).set_stroke(WHITE, 2).set_fill(GREY_B, 0.3).move_to(c)
    if k == 3:
        return VGroup(Ellipse(width=1.4, height=0.42).set_stroke(WHITE, 2).set_fill(GREY_B, 0.3),
                      Ellipse(width=1.4, height=0.42).set_stroke(GREY_C, 1.5).shift(0.06 * DOWN)).move_to(c)
    se = SH.small_example()
    if k == 4:
        return walk_icon(se["pts"], 0.9, c, YELLOW, 3)
    return walk_icon(numbers()["walk64"], 1.7, c, YELLOW, 1.0)


# ---------------------------------------------------------------- 各段
def part_question(s):
    y0 = 0.0
    x0, x1 = -4.5, 4.5
    line = Line(np.array([x0, y0, 0]), np.array([x1, y0, 0])).set_stroke(GREY_B, 3)
    ticks = VGroup(*[Line(np.array([x, y0 - 0.1, 0]), np.array([x, y0 + 0.1, 0])).set_stroke(GREY_B, 3)
                     for x in (x0, 0, x1)])
    nums = VGroup(*[tx(t, 30, GREY_A).move_to(np.array([x, y0 - 0.45, 0])) for t, x in (("0", x0), ("0.5", 0), ("1", x1))])
    ball = shape_icon(0, (x0, 1.2))
    rod = shape_icon(1, (x1, 1.2))
    lb = zh("圓球", font_size=26, color=GREY_A).next_to(nums[0], DOWN, buff=0.2)
    lr = zh("一條線", font_size=26, color=GREY_A).next_to(nums[2], DOWN, buff=0.2)
    q = tx(r"?", 60, YELLOW).move_to(np.array([-1.5, y0 + 0.9, 0]))
    arr = Arrow(q.get_bottom() + 0.05 * DOWN, np.array([-1.5, y0 + 0.12, 0]), buff=0, thickness=3).set_color(YELLOW)
    qa = VGroup(q, arr)
    ask = at(zh("三維漫步平均起來有多扁？", font_size=34, color=YELLOW), 0, 2.5)
    assert_in_frame(line, nums, ball, rod, lb, lr, qa, ask)
    for a, b in ((ask, ball), (ask, rod), (ask, qa), (ball, line), (rod, line), (lb, lr), (qa, nums)):
        assert_apart(a, b)
    s.play(FadeIn(ask), run_time=0.6)
    s.say("問：三維漫步平均起來有多扁？", wait=1.5)
    s.play(ShowCreation(line), FadeIn(ticks), FadeIn(nums), run_time=0.8)
    s.play(FadeIn(ball), FadeIn(lb), FadeIn(rod), FadeIn(lr), run_time=0.8)
    s.say("0 是圓球、1 是一條線", wait=2.6)
    s.play(FadeIn(qa), run_time=0.5)
    s.play(qa.animate.shift(3.0 * RIGHT), run_time=1.2)
    s.play(qa.animate.shift(1.5 * LEFT), run_time=0.8)
    s.say("在 0 到 1 之間猜一個數", wait=3.0)
    s.clear_stage()


def part_shapes(s):
    nb = numbers()
    e1 = eq_label(tx(r"\sum_{i<j}\big(R_i^2-R_j^2\big)^2=d\,\mathrm{Tr}\,\mathbf T^2-\big(\mathrm{Tr}\,\mathbf T\big)^2",
                     32), "(6.10)(6.13)")
    e2 = VGroup(zh("個別非球度＝", font_size=28), tx(r"\frac{\sum_{i<j}\big(R_i^2-R_j^2\big)^2}"
                                                       r"{(d-1)\Big(\sum_i R_i^2\Big)^2}", 34)).arrange(RIGHT, buff=0.2)
    at(e1, 0, 2.6)
    at(e2, 0, 1.4)
    xs = (-5.55, -3.35, -1.15, 1.05, 3.25, 5.55)
    names = ["圓球", "線", "雪茄", "扁盤", "W602 小例子", "Fig. 6.4 主例"]
    triples = [r"(1,1,1)", r"(1,0,0)", r"(4,1,1)", r"(1,1,0)", r"\Big(\frac{3}{5},\frac{1}{5}\Big)",
               r"(%s,\,%s)" % nb["fig64_eig"]]
    vals = [frac_tex(v) for v in nb["shape_vals"]] + [frac_tex(nb["small_a"]), nb["fig64_a"]]
    cols = VGroup()
    for k in range(6):
        nm = zh(names[k], font_size=24, color=GREY_A).move_to(np.array([xs[k], 0.45, 0]))
        ic = shape_icon(k, (xs[k], -0.4))
        tr = fit(tx(triples[k], 26, GREY_A).move_to(np.array([xs[k], -1.35, 0])), 2.05)
        vl = tx(vals[k], 34, YELLOW).move_to(np.array([xs[k], -2.15, 0]))
        cols.add(VGroup(nm, ic, tr, vl))
    assert_in_frame(e1, e2, cols)
    assert_apart(e1, e2)
    assert_apart(e2, cols)
    for k in range(6):
        for i in range(3):
            assert_apart(cols[k][i], cols[k][i + 1])
        if k < 5:
            assert_apart(cols[k], cols[k + 1])
    s.play(Write(e1), run_time=1.2)
    s.say("特徵值兩兩相減再平方：用不變量就能算", wait=2.9)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("除以 ", "$(d-1)\\big(\\sum R_i^2\\big)^2$", "：介於 0 和 1，與大小無關", wait=2.9)
    s.say("括號裡是三個主迴轉半徑的平方 ", "$(R_1^2,R_2^2,R_3^2)$", wait=2.4)
    s.play(FadeIn(cols[0]), FadeIn(cols[1]), run_time=0.8)
    s.say("圓球 ", "$(1,1,1)$", "：0；一條線 ", "$(1,0,0)$", "：1", wait=2.9)
    s.play(FadeIn(cols[2]), FadeIn(cols[3]), run_time=0.8)
    s.say("雪茄 ", "$(4,1,1)$", " 和扁盤 ", "$(1,1,0)$", "：都是 ", "$\\tfrac14$", wait=2.9)
    s.play(FadeIn(cols[4]), FadeIn(cols[5]), run_time=0.8)
    s.say("W602 的小例子：", "$\\tfrac14$", "；Fig. 6.4 的主例：%s" % nb["fig64_a"], wait=2.9)
    s.clear_stage()


def part_formula(s):
    nb = numbers()
    e3 = eq_label(tx(r"A_d=\frac{\sum_{i<j}\big\langle\big(R_i^2-R_j^2\big)^2\big\rangle}"
                     r"{(d-1)\Big\langle\Big(\sum_i R_i^2\Big)^2\Big\rangle}=\frac{4+2d}{4+5d}", 36),
                  "(6.12)(6.14)(6.15)")
    at(e3, 0, 1.75)
    rows = [["d"] + [str(d) for d in DS] + [r"\infty"],
            ["A_d"] + [frac_tex(nb["Ad"][d]) for d in DS] + [r"\frac{2}{5}"],
            [VectorizedPoint()] + [nb["Ad_dec"][d] for d in DS] + ["0.4000"]]
    tab = D.tex_table(rows, [1.4, 1.5, 1.5, 1.5, 1.5, 1.5], row_h=0.78, font_size=32)
    at(tab, 0, -1.0)
    col3 = VGroup(*[D.table_cell(tab, i, 2, 6) for i in (1, 2)])
    assert_in_frame(e3, tab)
    assert_apart(e3, tab)
    assert_table(tab, 6)
    s.play(Write(e3), run_time=1.8)
    s.say("平均非球度：分子、分母各自先平均，再相除", wait=2.9)
    s.say("揭曉：", "$A_d=\\tfrac{4+2d}{4+5d}$", "（W605 會推導）", wait=2.9)
    s.play(FadeIn(tab), run_time=1.0)
    s.play(col3.animate.set_color(YELLOW), run_time=0.6)
    s.say("三維：", "$\\tfrac{10}{19}=%s$" % nb["Ad_dec"][3], "，大約一半扁", wait=2.9)
    s.say("二維 ", "$\\tfrac47$", "、四維 ", "$\\tfrac12$", "、十維 ", "$\\tfrac49$", wait=2.4)
    s.say("維度越高越小；", "$d\\to\\infty$", " 時趨近 ", "$\\tfrac25$", wait=2.9)
    s.clear_stage()


def hist_axes():
    ax = D.plot_axes((0, 1, 0.1), (0, 2000, 500), width=8.4, height=3.4, x_ticks=(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     y_ticks=(500, 1000, 1500, 2000), x_tick_labels=["0", "0.2", "0.4", "0.6", "0.8", "1"])
    ax.move_to(np.array([-1.35, -0.45, 0]))
    return ax


def hist_bars(ax, hist):
    bars = VGroup()
    for i, h in enumerate(hist):
        p0, p1 = ax.c2p(i * 0.05, 0), ax.c2p((i + 1) * 0.05, max(h, 1))
        r = Rectangle(width=(p1[0] - p0[0]) * 0.92, height=max(p1[1] - p0[1], 1e-3))
        r.move_to((p0 + p1) / 2).set_fill(BAR_COL, 0.75).set_stroke(BAR_COL, 1)
        bars.add(r)
    return bars


def vline(ax, x, color, top=2000):
    return DashedLine(ax.c2p(x, 0), ax.c2p(x, top), dash_length=0.1).set_stroke(color, 3.5)


def part_hist(s):
    nb = numbers()
    ax = hist_axes()
    bars = hist_bars(ax, nb["hist"])
    ttl = at(mixed("Fig. 6.5：三維格子 2 萬條、各 100 步", font_size=28), -1.35, 2.55)
    xl = zh("個別非球度", font_size=22, color=GREY_B).next_to(ax.x_axis.get_end(), DOWN, buff=0.45)
    la = vline(ax, nb["A_sim_raw"], A_COL)
    r1 = mixed("$A_3$", " 平均的比", font_size=26, color=A_COL)
    r2 = mixed("模擬 %s" % nb["A_sim"], font_size=26, color=A_COL)
    r3 = mixed("$N=100$", " 精確 %s" % nb["A_ex"], font_size=24)
    r4 = mixed("$N\\to\\infty$", "：%s" % nb["A_inf"], font_size=24)
    r5 = mixed("小於 0.5：%s%%" % nb["lt"], font_size=26, color=LOW_COL)
    panel = VGroup(r1, r2, r3, r4, r5).arrange(DOWN, buff=0.32, aligned_edge=LEFT)
    panel.move_to(np.array([5.0, 0.2, 0]))
    assert_in_frame(ax, bars, ttl, xl, panel, la)
    for a, b in ((ttl, ax), (ttl, panel), (ax, panel), (xl, panel), (bars, panel)):
        assert_apart(a, b)
    s.play(FadeIn(ax), FadeIn(ttl), FadeIn(xl), run_time=0.8)
    s.say("Fig. 6.5：三維格子 2 萬條，各走 100 步", wait=0.2)
    s.play(LaggedStart(*[GrowFromEdge(r, DOWN) for r in bars], lag_ratio=0.06), run_time=2.0)
    s.say("每條算一個非球度：從 0 到 1 都有", wait=2.6)
    s.say("最多的在 0.20 到 0.25 之間", wait=2.2)
    s.play(ShowCreation(la), FadeIn(r1), FadeIn(r2), run_time=0.8)
    s.say("平均非球度（平均的比）：模擬 %s" % nb["A_sim"], wait=2.6)
    s.play(FadeIn(r3), FadeIn(r4), run_time=0.6)
    s.say("$N=100$", " 的精確值 %s；" % nb["A_ex"], "$N\\to\\infty$", "：", "$\\tfrac{10}{19}=%s$" % nb["A_inf"],
          wait=2.9)
    s.play(*[bars[i].animate.set_fill(LOW_COL, 0.85).set_stroke(LOW_COL, 1) for i in range(10)], FadeIn(r5),
           run_time=0.8)
    s.say("小於 0.5 的占 %s%%：平均不代表典型" % nb["lt"], wait=2.9)
    return dict(ax=ax, bars=bars, panel=panel, la=la, ttl=ttl, xl=xl)


def part_aprime(s, g):
    nb = numbers()
    ax = g["ax"]
    e4 = eq_label(tx(r"A'_d=\Big\langle\frac{\sum_{i<j}\big(R_i^2-R_j^2\big)^2}{(d-1)\big(\sum_i R_i^2\big)^2}"
                     r"\Big\rangle", 32), "(6.16)--(6.18)")
    at(e4, -1.35, 2.45)
    lp = vline(ax, nb["Ap_sim_raw"], AP_COL)
    p1 = mixed("$A'_3$", " 個別值的平均", font_size=26, color=AP_COL)
    p2 = mixed("$A'_3=%s$" % nb["Ap"][3], font_size=26, color=AP_COL)
    p3 = mixed("模擬 %s" % nb["Ap_sim"], font_size=26, color=AP_COL)
    p4 = mixed("二維 %s" % nb["Ap"][2], font_size=24)
    p5 = mixed("四維 %s" % nb["Ap"][4], font_size=24)
    panel = VGroup(p1, p2, p3, p4, p5).arrange(DOWN, buff=0.32, aligned_edge=LEFT)
    panel.move_to(np.array([5.0, 0.2, 0]))
    assert_in_frame(e4, panel, lp)
    for a, b in ((e4, ax), (e4, panel), (ax, panel)):
        assert_apart(a, b)
    s.play(FadeOut(g["ttl"]), FadeOut(g["panel"]), run_time=0.5)
    s.play(FadeIn(e4), run_time=0.8)
    s.say("另一種平均：先算每條的值，再平均（", "$A'_d$", "）", wait=2.6)
    s.play(ShowCreation(lp), FadeIn(p1), FadeIn(p2), FadeIn(p3), run_time=0.8)
    s.say("$A'_3=%s$" % nb["Ap"][3], "，模擬 %s：比 " % nb["Ap_sim"], "$\\tfrac{10}{19}$", " 小很多", wait=2.9)
    s.say("比值的平均 ", "$\\ne$", " 平均的比", wait=2.4)
    s.play(FadeIn(p4), FadeIn(p5), run_time=0.6)
    s.say("二維 %s、四維 %s：幾乎不隨維度變" % (nb["Ap"][2], nb["Ap"][4]), wait=2.9)
    s.clear_stage()


def part_erratum(s):
    nb = numbers()
    box = D.erratum_box("p. 135，(6.16)", r"(d-1)\sum_{k=1}^{N}R_k^2", r"(d-1)\Big(\sum_{k=1}^{d}R_k^2\Big)^2",
                        width=8.6)
    at(box, 0, 1.75)
    hdr = lambda t: zh(t, font_size=26, color=YELLOW)
    rows = [[VectorizedPoint(), hdr("原本"), hdr("整條放大 2 倍")],
            [zh("照書上", font_size=26, color=D.WRONG_COLOR), nb["book"][0], nb["book"][1]],
            [zh("更正後", font_size=26, color=D.RIGHT_COLOR), frac_tex(nb["fixed"][0]), frac_tex(nb["fixed"][1])]]
    tab = D.tex_table(rows, [2.2, 2.0, 3.0], row_h=0.72, font_size=32)
    at(tab, 0, -1.5)
    note = at(mixed("W602 小例子 ", "$(R_1^2,R_2^2)=\\big(\\frac35,\\frac15\\big)$", "，放大 2 倍 → 乘 4",
                    font_size=24, color=GREY_A), 0, 0.05)
    assert_in_frame(box, tab, note)
    for a, b in ((box, note), (note, tab), (box, tab)):
        assert_apart(a, b)
    assert_table(tab, 3)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (6.16) 的分母寫錯了：少了平方，也加錯範圍", wait=2.9)
    s.say("正確的分母：特徵值總和的平方，只加到 ", "$d$", wait=2.6)
    s.play(FadeIn(note), FadeIn(tab), run_time=0.8)
    s.say("照書上，整條漫步放大 2 倍，值從 %s 變 %s" % nb["book"], wait=2.9)
    s.say("就不是形狀了；更正後放大前後都是 ", "$\\tfrac14$", wait=2.9)
    s.clear_stage()


def part_close(s):
    s.say("三維漫步平均「一半扁」（", "$\\tfrac{10}{19}$", "）", wait=2.6)
    s.say("但個別的形狀，從近乎圓球到近乎直線都有", wait=3.0)


class W603Asphericity(ChuScene):
    def construct(self):
        self.show_header("W603　非球度：0 是圓球、1 是一條線")
        check_math()
        part_question(self)
        part_shapes(self)
        part_formula(self)
        g = part_hist(self)
        part_aprime(self, g)
        part_erratum(self)
        part_close(self)
        self.wait(1.0)
