"""
W703 把自迴避展開：1 − f（§7.2.2，pp. 172–176，(7.11)–(7.19)，Fig. 7.1–7.4，Ex. 7.2）

問：「不准踩到自己」很難算。能不能先算普通漫步，再一點一點扣掉？
  1. (7.17) e^{−v} = 1 − f：f 只在 x ≈ 0 附近不是 0；嚴格自迴避：重合時 f = 1。每一對點乘 (1 − f) 再展開：
     0 階＝普通漫步；一階 (7.19)＝一條線＋一個圈＋一條虛線（Fig. 7.2）。筆誤框 (7.11)。
  2. 弧線圖：二階 3 種＝Fig. 7.3 (a) 巢狀、(b) 並排、(c) 交叉；三階 15 種（Ex. 7.2）；(2k−1)!! = 1, 3, 15, 105。
  3. 二維 N=4 逐階：T_k = 256, 228, 96, 28, 4（T_1 = 3·4²·4 + 36）；部分和 256 → 28 → 124 → 96 → 100（排容原理）；
     二階 96 = 巢狀 16 + 交叉 32 + 共用一點 48。
  4. N=8：65536 → −119332 → 190940 → −222036 → … → 5916（17 項）：低階的部分和連正負都不對。
畫面上的數字都由 utils.walk_saw 算出，check_math() 核對。
"""
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

ARC_COL = ORANGE
POS_COL = "#F2A33A"
NEG_COL = "#4A90D9"
OK_COL = "#7CE38B"
KIND_ZH = {"nested": "巢狀", "sequential": "並排", "crossed": "交叉"}


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        xs = np.linspace(-20, 20, 20001)
        g = SW.gauss1(6, xs, book711=True)
        var_book = float(np.sum(xs ** 2 * g) / np.sum(g))
    return dict(
        T4=list(SW.ie_terms(4)), P4=list(SW.ie_partial(4)), o1=SW.order1_formula(4), o2=dict(SW.order2_types(4)),
        P8=list(SW.ie_partial(8)), c8=SW.saw_c(8), c4=SW.saw_c(4),
        m2=list(SW.matchings(2)), m3=list(SW.matchings(3)), cc=[SW.chord_count(k) for k in (1, 2, 3, 4)],
        var_book=var_book, var_lat=SW.lattice_var(6),
    )


def check_math():
    nb = numbers()
    assert nb["T4"] == [256, 228, 96, 28, 4] and nb["P4"] == [256, 28, 124, 96, 100]
    assert nb["o1"] == 228 == 3 * 4 ** 2 * 4 + 36 == nb["T4"][1]
    assert sum(4 ** (4 - 2) * math.comb(2, 1) ** 2 for _ in range(3)) == 192 and math.comb(4, 2) ** 2 == 36
    assert nb["o2"] == {"sequential": 0, "nested": 16, "crossed": 32, "shared": 48}
    assert 16 + 32 + 48 == 96 == nb["T4"][2]
    assert 256 - 228 + 96 - 28 + 4 == 100 == nb["P4"][-1] == nb["c4"]
    P8 = nb["P8"]
    assert len(P8) == 17 and P8[:4] == [65536, -119332, 190940, -222036] and P8[4] == 248468
    assert P8[-1] == 5916 == nb["c8"]
    assert nb["cc"] == [1, 3, 15, 105] and len(nb["m2"]) == 3 and len(nb["m3"]) == 15
    assert all(SW.chord_count(k) == math.prod(range(1, 2 * k, 2)) for k in (1, 2, 3, 4))
    so = SW.SECOND_ORDER
    assert so["a"] == ("nested", ((1, 4), (2, 3))) and so["b"] == ("sequential", ((1, 2), (3, 4)))
    assert so["c"] == ("crossed", ((1, 3), (2, 4)))
    assert sorted(m for _, m in so.values()) == sorted(nb["m2"])
    assert f"{nb['var_book']:.4f}" == "3.0000" and 6 / 2 == 3                 # 照 (7.11)：M/2
    assert nb["var_lat"] == 2 == 6 // 3                                        # 格子漫步：M/3


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def chord_diagram(m, width, origin, color=ARC_COL, dot_r=0.05, stroke=3, dashed=True):
    """弧線圖：水平實線（漫步）上 2k 個點，每對 (a, b) 一條上方的虛線弧。origin＝實線的中點。"""
    origin = np.array(origin, dtype=float)
    n = 2 * len(m)
    sp = width / n
    xs = [origin + (-width / 2 + sp * (i + 0.5)) * RIGHT for i in range(n)]
    line = Line(origin - width / 2 * RIGHT, origin + width / 2 * RIGHT).set_stroke(WHITE, stroke)
    arcs = VGroup()
    for a, b in m:
        arc = ArcBetweenPoints(xs[a - 1], xs[b - 1], angle=-PI).set_stroke(color, stroke)
        arcs.add(DashedVMobject(arc, num_dashes=max(6, int(5 * (b - a)))) if dashed else arc)
    dots = VGroup(*[Dot(p, radius=dot_r).set_fill(WHITE, 1) for p in xs])
    return VGroup(line, arcs, dots)


def first_order_picture(center, r=0.6, gap=0.45):
    """Fig. 7.2：x → u1 → 繞一圈 → u2 → y，u1、u2 之間一條虛線（f）。"""
    c = np.array(center, dtype=float)
    cc = c + r * UP

    def on(theta):
        return cc + r * np.array([math.cos(theta), math.sin(theta), 0])

    th0 = -PI / 2 - gap
    loop = ParametricCurve(lambda t: on(th0 - t * (2 * PI - 2 * gap)), t_range=(0, 1, 0.01)).set_stroke(WHITE, 3)
    u1, u2 = on(th0), on(-PI / 2 + gap)
    xp, yp = c + 2.0 * LEFT + 0.0 * UP, c + 2.0 * RIGHT
    l1 = Line(xp, u1).set_stroke(WHITE, 3)
    l2 = Line(u2, yp).set_stroke(WHITE, 3)
    dash = DashedLine(u1, u2, dash_length=0.06).set_stroke(ARC_COL, 4)
    dots = VGroup(*[Dot(p, radius=0.06).set_fill(WHITE, 1) for p in (xp, u1, u2, yp)])
    labs = VGroup(tx(r"\vec x", 28).next_to(xp, DOWN, buff=0.15), tx(r"\vec u_1", 28).next_to(u1, DOWN + LEFT, buff=0.08),
                  tx(r"\vec u_2", 28).next_to(u2, DOWN + RIGHT, buff=0.08), tx(r"\vec y", 28).next_to(yp, DOWN, buff=0.15))
    return VGroup(VGroup(l1, loop, l2), dash, dots, labs)


def signed_bars(vals, xs, ux, uy, origin, width=0.7):
    """有正有負的長條（往下畫負的），附數值標籤。"""
    origin = np.array(origin, dtype=float)
    bars, labs = VGroup(), VGroup()
    for x, v in zip(xs, vals):
        h = v * uy
        col = POS_COL if v > 0 else NEG_COL
        r = Rectangle(width=width * ux, height=max(abs(h), 0.02)).set_fill(col, 0.8).set_stroke(col, 1.5)
        r.move_to(origin + x * ux * RIGHT + h / 2 * UP)
        bars.add(r)
        t = tx(str(v), 22, col)
        t.next_to(r, UP if v > 0 else DOWN, buff=0.08)
        labs.add(t)
    return bars, labs


# ---------------------------------------------------------------- 各段
def part_question(s):
    e = at(tx(r"c_N=4^N-(\cdots)+(\cdots)-(\cdots)+\cdots", 44), 0, 0.8)
    x0 = e.get_left()[0]                                     # 用前綴的寬度找出 c_N、4^N 的位置
    xa = x0 + tx(r"c_N", 44).get_width() / 2
    xb = x0 + tx(r"c_N=4^N", 44).get_width() - tx(r"4^N", 44).get_width() / 2
    ya = e.get_bottom()[1] - 0.4
    a = at(mixed("自迴避的條數", font_size=26, color=YELLOW), xa, ya)
    b = at(mixed("普通漫步", font_size=26, color=TEAL_B), xb, ya - 0.55)       # 錯開兩行，不互相重疊
    assert_in_frame(e, a, b)
    for p, q in ((a, b), (e, a), (e, b)):
        assert_apart(p, q)
    s.say("「不准踩到自己」很難算", wait=2.0)
    s.play(Write(e), run_time=1.4)
    s.play(FadeIn(a), FadeIn(b), run_time=0.5)
    s.say("能不能先算普通漫步，再一點一點扣掉？", wait=2.8)
    s.clear_stage()


def part_expand(s):
    nb = numbers()
    e17 = at(tx(r"e^{-v}=1-f", 40), -4.0, 2.55)
    t17 = tx(r"(7.17)", 24, GREY_B).next_to(e17, RIGHT, buff=0.35)
    ax = D.plot_axes((-3, 3, 1), (0, 1.2, 0.5), width=4.4, height=2.1, x_ticks=(-2, 0, 2), y_ticks=(1,),
                     x_label="x")
    ax.move_to(np.array([-4.0, 0.55, 0]))
    curve = ax.get_graph(lambda x: 1 - math.exp(-4 * math.exp(-x * x / 0.25)), x_range=(-3, 3, 0.02)).set_stroke(
        ARC_COL, 4)
    flab = tx(r"f(x)", 28, ARC_COL).next_to(ax.c2p(0.55, 0.85), RIGHT, buff=0.1)
    hard = Dot(ax.c2p(0, 1), radius=0.08).set_fill(RED, 1)
    hlab = at(mixed("嚴格：重合時 ", "$f=1$", font_size=24, color=RED), -4.0, -1.15)
    ex = at(tx(r"\prod_{i>j}(1-f_{ij})=1-\sum f_{ij}+\sum\sum f_{ij}f_{kl}-\cdots", 32), 2.75, 2.5)
    r0 = at(mixed("0 階：", "$1$", "，就是普通漫步", font_size=28), 2.75, 1.45)
    r1 = at(mixed("1 階：扣掉「某一對點重合」", font_size=28), 2.75, 0.75)
    r2 = at(mixed("2 階：加回「兩對同時重合」……", font_size=28), 2.75, 0.05)
    for m in (ex, r0, r1, r2):
        fit(m, 6.6)
    for m in (r1, r2):
        m.align_to(r0, LEFT)
    assert_in_frame(e17, t17, ax, curve, flab, hard, hlab, ex, r0, r1, r2)
    for a, b in ((e17, ax), (ax, hlab), (t17, ex), (ax, ex), (ax, r0), (ax, r1), (ax, r2), (ex, r0), (r0, r1),
                 (r1, r2), (hlab, r2)):
        assert_apart(a, b)
    s.play(Write(e17), FadeIn(t17), run_time=1.0)
    s.say("把每一對點的權重寫成 ", "$e^{-v}=1-f$", wait=2.4)
    s.play(FadeIn(ax), ShowCreation(curve), FadeIn(flab), run_time=1.2)
    s.say("$f$", " 只在兩點靠得很近（", "$x\\approx0$", "）時不是 0", wait=2.6)
    s.play(FadeIn(hard, scale=1.5), FadeIn(hlab), run_time=0.6)
    s.say("嚴格自迴避：重合時 ", "$f=1$", "，否則 0", wait=2.4)
    s.play(Write(ex), run_time=1.4)
    s.say("每一對點乘 ", "$(1-f)$", "，再乘開", wait=2.4)
    s.play(FadeIn(r0), run_time=0.5)
    s.say("0 階：全部取 1，就是普通漫步", wait=2.2)
    s.play(FadeIn(r1), run_time=0.5)
    s.play(FadeIn(r2), run_time=0.5)
    s.say("1 階扣掉重合的，2 階加回扣太多的……", wait=2.8)
    s.clear_stage()
    # (7.19) 一階
    e19 = at(tx(r"-\sum_{\vec u_1,\vec u_2}C(M_1;\vec x,\vec u_1)\,C(M_2;\vec u_1,\vec u_2)\,"
                r"C(M_3;\vec u_2,\vec y)\,f(|\vec u_1-\vec u_2|)", 32), -0.4, 2.5)
    t19 = tx(r"(7.19)", 24, GREY_B).next_to(e19, RIGHT, buff=0.35)
    pic = first_order_picture(np.array([-3.6, 0.1, 0]))
    cd = chord_diagram(((1, 2),), 3.2, np.array([3.4, 0.2, 0]))
    cd_labs = VGroup(tx(r"\vec u_1", 26).next_to(cd[2][0], DOWN, buff=0.15),
                     tx(r"\vec u_2", 26).next_to(cd[2][1], DOWN, buff=0.15))
    arrow = Arrow(np.array([-1.0, 0.5, 0]), np.array([1.4, 0.5, 0]), buff=0).set_stroke(GREY_B, 3)
    cap7 = at(zh("一條線＋一個圈＋一條虛線（Fig. 7.2）", font_size=26, color=GREY_A), -3.6, -1.0)
    box = D.erratum_box("p. 173，(7.11)", r"e^{-3|\vec x-\vec r_1|^2/3M_1a^2}", r"e^{-3|\vec x-\vec r_1|^2/2M_1a^2}",
                        width=6.6, font_size=44)
    at(box, 3.0, -1.75)
    fit(e19, 12.0)
    assert_in_frame(e19, t19, pic, cd, cd_labs, arrow, cap7, box)
    for a, b in ((e19, pic), (e19, cd), (pic, cd), (pic, cap7), (cd, box), (cd_labs, box), (cap7, box), (pic, arrow),
                 (cd, arrow), (t19, cd)):
        assert_apart(a, b)
    s.play(Write(e19), FadeIn(t19), run_time=1.6)
    s.say("一階 (7.19)：走到 ", "$\\vec u_1$", "、繞一圈回到 ", "$\\vec u_2$", "、再走到 ", "$\\vec y$", wait=0.3)
    s.play(ShowCreation(pic[0]), FadeIn(pic[2]), FadeIn(pic[3]), run_time=1.8)
    s.play(ShowCreation(pic[1]), FadeIn(cap7), run_time=0.6)
    s.say("虛線是 ", "$f$", "：", "$\\vec u_1$", " 和 ", "$\\vec u_2$", " 重合才算", wait=2.6)
    s.play(GrowArrow(arrow), FadeIn(cd[0]), FadeIn(cd[2]), FadeIn(cd_labs), run_time=0.8)
    s.play(ShowCreation(cd[1]), run_time=0.8)
    s.say("把漫步拉直：一條線上兩個點，一條虛線弧", wait=2.6)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (7.11) 的指數分母寫成 ", "$3M_1a^2$", wait=2.4)
    s.say("照書上，一個方向的變異數變成 ", "$M/2$", "（", "$M=6$", "：3）", wait=2.8)
    s.say("格子漫步是 ", "$M/3$", "（", "$M=6$", "：2）", wait=2.4)
    s.clear_stage()


def part_chords(s):
    nb = numbers()
    so = SW.SECOND_ORDER
    xs = (-4.5, 0.0, 4.5)
    big = VGroup()
    labs = VGroup()
    for key, x in zip(("a", "b", "c"), xs):
        kind, m = so[key]
        big.add(chord_diagram(m, 3.4, np.array([x, 0.3, 0])))
        labs.add(at(mixed("(%s) %s" % (key, KIND_ZH[kind]), font_size=28), x, -0.5))
    head = at(mixed("二階：4 個點配成 2 對，", "$3$", " 種（Fig. 7.3）", font_size=30), 0, 2.4)
    assert_in_frame(big, labs, head)
    for i in range(3):
        assert_apart(big[i], labs[i])
        assert_apart(head, big[i])
        for j in range(i + 1, 3):
            assert_apart(big[i], big[j], 0.3)
    s.play(FadeIn(head), run_time=0.6)
    s.say("二階：兩對點同時重合，弧線怎麼擺？", wait=0.3)
    for i in range(3):
        s.play(FadeIn(big[i][0]), FadeIn(big[i][2]), run_time=0.4)
        s.play(ShowCreation(big[i][1]), FadeIn(labs[i]), run_time=0.8)
    s.say("巢狀、並排、交叉：正好 3 種", wait=2.6)
    s.clear_stage()
    small = VGroup()
    for i, m in enumerate(nb["m3"]):
        small.add(chord_diagram(m, 2.2, np.array([-5.2 + 2.6 * (i % 5), 1.45 - 1.3 * (i // 5), 0]), dot_r=0.04,
                                stroke=2.5))
    head3 = at(mixed("三階：6 個點配成 3 對，", "$15$", " 種（Ex. 7.2）", font_size=30), 0, 2.75)
    row = at(mixed("$(2k-1)!!=%d,\\ %d,\\ %d,\\ %d$" % tuple(nb["cc"]), "：和第 6 章的 Wick 配對一樣",
                   font_size=30, color=YELLOW), 0, -2.2)
    fit(row, 12.5)
    assert_in_frame(small, head3, row)
    for i in range(15):
        assert_apart(small[i], row)
        assert_apart(small[i], head3)
        for j in range(i + 1, 15):
            assert_apart(small[i], small[j], 0.1)
    s.play(FadeIn(head3), run_time=0.5)
    s.say("三階：15 種（Ex. 7.2）", wait=0.3)
    s.play(LaggedStart(*[FadeIn(m) for m in small], lag_ratio=0.1), run_time=2.4)
    s.wait(1.2)
    s.play(FadeIn(row), run_time=0.6)
    s.say("$k$", " 階有 ", "$(2k-1)!!$", " 種：1、3、15、105", wait=2.6)
    s.say("和第 6 章的 Wick 配對一樣", wait=2.4)
    s.clear_stage()


def part_n4(s):
    nb = numbers()
    T, P = nb["T4"], nb["P4"]
    ux, uy = 0.95, 3.4 / 256
    org = np.array([-5.7, -2.3, 0])
    chart = D.bar_chart({k: P[k] for k in range(5)}, ux=ux, uy=uy, origin=org, width=0.7, labels=True,
                        font_size=24)
    xl = VGroup(*[tx(r"P_%d" % k, 24, GREY_B).next_to(org + k * ux * RIGHT, DOWN, buff=0.15) for k in range(5)])
    base = Line(org - 0.6 * ux * RIGHT, org + 4.6 * ux * RIGHT).set_stroke(GREY_C, 2)
    target = DashedLine(org - 0.6 * ux * RIGHT + P[-1] * uy * UP, org + 4.6 * ux * RIGHT + P[-1] * uy * UP,
                        dash_length=0.1).set_stroke(OK_COL, 2)
    tlab = mixed("$c_4=%d$" % nb["c4"], font_size=24, color=OK_COL).next_to(target, RIGHT, buff=0.12)
    rows = [
        mixed("$N=4$", "：第 ", "$k$", " 階 ", "$T_k=%s$" % ",\\ ".join(map(str, T)), font_size=28),
        tx(r"T_1=3\cdot4^2\cdot4+36=%d" % nb["o1"], 32),
        zh("三對隔 2 步的點、一對隔 4 步", font_size=24, color=GREY_A),
        tx(r"T_2=%d=%d+%d+%d" % (T[2], nb["o2"]["nested"], nb["o2"]["crossed"], nb["o2"]["shared"]), 32),
        zh("巢狀 16＋交叉 32＋共用一點 48", font_size=24, color=GREY_A),
        tx(r"P_K=\sum_{k\le K}(-1)^kT_k", 32),
        tx(r"256-228+96-28+4=100", 32, YELLOW),
    ]
    for m, y in zip(rows, (2.55, 1.75, 1.2, 0.5, -0.05, -0.95, -1.95)):
        at(m, 3.2, y)
        fit(m, 6.2)
    assert_in_frame(chart, xl, base, target, tlab, *rows)
    for m in rows:
        assert_apart(chart, m)
        assert_apart(tlab, m)
    for i in range(6):
        assert_apart(rows[i], rows[i + 1])
    s.play(FadeIn(rows[0]), run_time=0.6)
    s.say("二維 ", "$N=4$", "：逐階算，第 ", "$k$", " 階記作 ", "$T_k$", wait=2.6)
    s.play(FadeIn(rows[1]), FadeIn(rows[2]), run_time=0.7)
    s.say("一階：能重合的點對隔偶數步", wait=2.4)
    s.say("隔 2 步的三對各 ", "$4^2\\cdot4$", "，隔 4 步的一對 36", wait=2.8)
    s.play(FadeIn(rows[3]), FadeIn(rows[4]), run_time=0.7)
    s.say("二階 96：巢狀 16、交叉 32、共用一點 48", wait=2.8)
    s.play(FadeIn(rows[5]), FadeIn(base), FadeIn(xl), run_time=0.7)
    s.say("部分和：一正一負交錯加上去", wait=0.3)
    for k in range(5):
        s.play(GrowFromEdge(chart[0][k], DOWN), FadeIn(chart[1][k]), run_time=0.5)
    s.play(ShowCreation(target), FadeIn(tlab), FadeIn(rows[6]), run_time=0.8)
    s.say("256、28、124、96、100：最後剛好是自迴避的 100 條", wait=3.0)
    s.say("這就是排容原理", wait=2.0)
    s.clear_stage()


def part_n8(s):
    nb = numbers()
    P = nb["P8"]

    def cell(v):
        return tx(str(v), 24, NEG_COL if v < 0 else (OK_COL if v == P[-1] else WHITE))

    r1 = [tx("K", 24, YELLOW)] + [tx(str(k), 24, YELLOW) for k in range(9)]
    r2 = [tx("P_K", 24, YELLOW)] + [cell(v) for v in P[:9]]
    r3 = [tx("K", 24, YELLOW)] + [tx(str(k), 24, YELLOW) for k in range(9, 17)] + [VectorizedPoint()]
    r4 = [tx("P_K", 24, YELLOW)] + [cell(v) for v in P[9:]] + [VectorizedPoint()]
    tab = D.tex_table([r1, r2, r3, r4], [1.1] + [1.3] * 9, row_h=0.48, font_size=24)
    tab.move_to(np.array([0, 2.12, 0]))
    ux, uy = 1.05, 1.35 / 250000
    org = np.array([-5.4, -0.95, 0])
    xs = [0, 1, 2, 3, 4, 6]
    bars, blabs = signed_bars(P[:5] + [P[-1]], xs, ux, uy, org)
    base = Line(org - 0.6 * ux * RIGHT, org + 6.6 * ux * RIGHT).set_stroke(GREY_C, 2)
    dots = tx(r"\cdots", 30, GREY_B).move_to(org + 5 * ux * RIGHT + 0.25 * UP)
    kl = VGroup(*[tx(r"P_{%d}" % k, 22, GREY_B).next_to(org + x * ux * RIGHT, DOWN if v > 0 else UP, buff=0.12)
                  for k, x, v in zip((0, 1, 2, 3, 4, 16), xs, P[:5] + [P[-1]])])
    note1 = at(mixed("真正的答案：", "$c_8=%d$" % nb["c8"], font_size=28, color=OK_COL), 4.4, 0.0)
    note2 = at(zh("低階的部分和連正負都不對", font_size=26, color=GREY_A), 4.4, -0.75)
    note3 = at(zh("（伏筆：W705）", font_size=24, color=GREY_B), 4.4, -1.4)
    assert_in_frame(tab, bars, blabs, base, dots, kl, note1, note2, note3)
    for a, b in ((tab, bars), (tab, blabs), (bars, note1), (blabs, note1), (blabs, note2), (bars, note2),
                 (note1, note2), (note2, note3), (bars, note3), (blabs, note3)):
        assert_apart(a, b)
    s.say("$N=8$", "：一共 17 項部分和", wait=0.3)
    s.play(FadeIn(tab), run_time=1.0)
    s.wait(1.6)
    s.play(FadeIn(base), FadeIn(kl), FadeIn(dots), run_time=0.5)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN if v > 0 else UP) for b, v in zip(bars, P[:5] + [P[-1]])],
                       lag_ratio=0.2), FadeIn(blabs), run_time=2.0)
    s.say("65536、", "$-119332$", "、190940、", "$-222036$", "……忽正忽負", wait=3.0)
    s.play(FadeIn(note1), run_time=0.6)
    s.say("全部 17 項加完，才得到 5916", wait=2.6)
    s.play(FadeIn(note2), FadeIn(note3), run_time=0.6)
    s.say("低階的部分和連正負都不對：要全部加起來才可信", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("$1-f$", " 展開把自迴避變成「一條線加很多弧」的圖", wait=2.8)
    s.say("只取前幾階，數字可能錯得離譜", wait=3.0)


class W703OneMinusF(ChuScene):
    def construct(self):
        self.show_header("W703　把自迴避展開：1 - f")
        check_math()
        part_question(self)
        part_expand(self)
        part_chords(self)
        part_n4(self)
        part_n8(self)
        part_close(self)
        self.wait(1.0)
