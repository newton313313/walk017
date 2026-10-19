"""
W1203 一迴圈與 (r, u) 平面的流（§12.2.4–12.2.5，pp. 292–297，(12.16)–(12.38)，Fig. 12.11–12.12）

問：兩張一迴圈圖各給多少？
  1. Fig. 12.11：(a) 一個頂點、殼層線自己接成圈 → 改 r；(b) 兩個頂點、兩條殼層線 → 改 u；(12.18)(12.19)；筆誤框 1（(12.17) 前）。
  2. Σ_shell 1 = N d δℓ；動量守恆 (12.26)(12.27)；(12.28)；乘上重新標度 ⇒ (12.37)(12.38)；筆誤框 2；本系列把 d 併進 u。
  3. 不動點：Gaussian (0, 0)、Wilson–Fisher；表 3 列 × 5 欄（精確與一階）。
  4. (r, u) 平面的流：6 條流線、兩個不動點；Gaussian 點的特徵值 2、1。
畫面上的數字都由 utils.walk_eps 算出，check_math() 核對。
模組層級的畫圖函式：frag_dot、self_loop、double_line、loop_a、loop_b、flow_axes、flow_curve。
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
from utils import walk_eps as E
from utils import lattice_draw as L
from utils import walk_draw as D

SHELL_COLOR = TEAL
EPS = 1.0
RMAX, UMAX = 0.45, 0.05
LINE_COLORS = (GREY_B, YELLOW, BLUE_B, GREY_B, YELLOW, BLUE_B)


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q):
    q = Fraction(q)
    sign = "-" if q < 0 else ""
    return r"%s\frac{%d}{%d}" % (sign, abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        starts = E.flow_starts()
        return {"shell": (E.shell_fraction(4, 0.01), E.shell_linear(4, 0.01)),
                "fp": {e: E.fixed_point(e) for e in E.NU_EPS},
                "low": {e: E.fixed_point_lowest(e) for e in E.NU_EPS},
                "starts": starts,
                "lines": [np.array(E.flow_line(r0, u0)) for r0, u0 in starts],
                "gauss": E.eigen(EPS, 1, "gauss"),
                "beta": E.beta(0.1, 0.01, EPS)}


def check_math():
    nb = numbers()
    a, b = nb["shell"]
    assert f4(a) == "0.0392" and f4(b) == "0.0400"
    dr, du = nb["beta"]
    assert abs(dr - (0.2 + 12 * 0.01 / 1.1)) < 1e-12 and abs(du - (0.01 - 36e-4 / 1.21)) < 1e-12
    assert E.NU_EPS == (Fraction(1, 10), Fraction(1, 2), Fraction(1))
    fp, low = nb["fp"], nb["low"]
    assert fp[Fraction(1)] == (Fraction(-1, 7), Fraction(1, 49))
    assert fp[Fraction(1, 2)] == (Fraction(-1, 13), Fraction(2, 169))
    assert fp[Fraction(1, 10)] == (Fraction(-1, 61), Fraction(10, 3721))
    for e in E.NU_EPS:
        assert fp[e] == (-e / (6 + e), e / (6 + e) ** 2)
        assert low[e] == (-e / 6, e / 36)
    assert [(f4(r), f4(u)) for r, u in (fp[Fraction(1)], fp[Fraction(1, 2)], fp[Fraction(1, 10)])] == \
        [("-0.1429", "0.0204"), ("-0.0769", "0.0118"), ("-0.0164", "0.0027")]
    assert [(f4(r), f4(u)) for r, u in (low[Fraction(1)], low[Fraction(1, 2)], low[Fraction(1, 10)])] == \
        [("-0.1667", "0.0278"), ("-0.0833", "0.0139"), ("-0.0167", "0.0028")]
    rw, uw = (float(x) for x in fp[Fraction(1)])
    assert len(nb["starts"]) == 6
    for k, pts in enumerate(nb["lines"]):
        assert len(pts) == 40
        assert np.all(np.abs(pts[:, 0]) <= RMAX + 1e-9) and np.all((pts[:, 1] >= 0) & (pts[:, 1] <= UMAX))
        if k % 3 == 1:
            assert abs(pts[-1, 0] - rw) < 0.002 and abs(pts[-1, 1] - uw) < 0.0002          # 臨界線流進 WF 點
        else:
            assert abs(abs(pts[-1, 0]) - RMAX) < 1e-6 and np.sign(pts[-1, 0]) == (k % 3 - 1)  # 跑出 |r| = 0.45
    assert [f4(x) for x in nb["gauss"]] == ["2.0000", "1.0000"]


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


# ---------------------------------------------------------------- 畫圖元件
def frag_dot(p):
    return Dot(p, radius=0.09).set_color(YELLOW)


def leg(p, direction, length=0.6):
    d = np.array(direction, dtype=float)
    return Line(p, p + length * d / np.linalg.norm(d)).set_stroke(WHITE, 4)


def loop_a(origin):
    """Fig. 12.11(a)：一個頂點、兩隻外腳、殼層線自己接成圈。"""
    o = np.array(origin, dtype=float)
    ring = Circle(radius=0.35).move_to(o + 0.35 * UP).set_stroke(SHELL_COLOR, 4)
    return VGroup(leg(o, LEFT), leg(o, RIGHT), ring, frag_dot(o))


def loop_b(origin):
    """Fig. 12.11(b)：兩個頂點、兩條殼層線、四隻外腳。"""
    o = np.array(origin, dtype=float)
    p, q = o + 0.55 * LEFT, o + 0.55 * RIGHT
    arcs = VGroup(ArcBetweenPoints(p, q, angle=0.9), ArcBetweenPoints(p, q, angle=-0.9)).set_stroke(SHELL_COLOR, 4)
    return VGroup(leg(p, UL), leg(p, DL), leg(q, UR), leg(q, DR), arcs, frag_dot(p), frag_dot(q))


def flow_axes():
    ax = D.plot_axes((-RMAX, RMAX + 0.04, 0.1), (0, UMAX + 0.004, 0.01), width=7.0, height=4.6,
                     x_ticks=(-0.4, -0.2, 0.2, 0.4), x_tick_labels=["-0.4", "-0.2", "0.2", "0.4"],
                     y_ticks=(0.01, 0.02, 0.03, 0.04, 0.05), x_label="r", y_label="u")
    return ax


def flow_curve(ax, pts, color):
    """一條流線（≤ 40 點；去掉重合的點）＋末端箭頭。"""
    ps = [ax.c2p(r, u) for r, u in pts]
    keep = [ps[0]]
    for p in ps[1:]:
        if np.linalg.norm(p - keep[-1]) > 1e-3:
            keep.append(p)
    line = VMobject().set_points_smoothly(keep).set_stroke(color, 4)
    end = keep[-1]
    back = next(p for p in reversed(keep) if np.linalg.norm(end - p) > 0.15)
    u = (end - back) / np.linalg.norm(end - back)
    tip = Arrow(end - 0.3 * u, end, buff=0, thickness=4).set_color(color)
    return VGroup(line, tip)


# ---------------------------------------------------------------- 各段
def part_diagrams(s):
    pa, pb = loop_a(ORIGIN), loop_b(ORIGIN)
    la = VGroup(zh("(a) 一個頂點 → 改 r", font_size=24, color=GREY_A),
                tx(r"\frac{6u}{N(r+1)}\sum_{\vec Q}\sum_{\vec q}SS", 30)).arrange(DOWN, buff=0.2)
    lb = VGroup(zh("(b) 兩個頂點 → 改 u", font_size=24, color=GREY_A),
                tx(r"\frac{36u^2}{N^2(r+1)^2}\sum_{\vec Q,\vec Q'}\sum S\cdots S\,\delta", 30)).arrange(DOWN, buff=0.2)
    ga = VGroup(pa, la).arrange(DOWN, buff=0.35)
    gb = VGroup(pb, lb).arrange(DOWN, buff=0.35)
    figs = VGroup(ga, gb).arrange(RIGHT, buff=1.5)
    tag = zh("Fig. 12.11　(12.18)(12.19)", font_size=22, color=GREY_B)
    col = VGroup(tag, figs).arrange(DOWN, buff=0.35)
    place_in(col, -6.6, 6.6, -0.4, 3.1)
    box = D.erratum_box("p. 293，(12.17) 前", r"\text{Figure 12.16(b)}", r"\text{Figure 12.11(b)}", width=6.5)
    box.move_to(np.array([0, -1.85, 0]))
    assert_in_frame(col, box)
    assert_apart(col, box)
    assert_apart(ga, gb)
    s.say("問：兩張一迴圈圖各給多少？", wait=1.8)
    s.play(FadeIn(tag), FadeIn(pa), run_time=1.0)
    s.say("(a) 一個頂點，兩隻殼層手自己接成圈", wait=2.0)
    s.play(Write(la), run_time=1.2)
    s.say("剩兩隻內區的手：像 ", "$SS$", " 的項——改的是 ", "$r$", wait=2.0)
    s.say("6 就是上一幕的 ", "$6u$", " 頂點；圈給 ", "$\\frac1{r+1}$", wait=2.0)
    s.say("$\\sum_{\\vec Q}$", " 跑遍殼層；殼層裡 ", "$r+q^2\\to r+1$", wait=2.2)
    s.play(FadeIn(pb), run_time=1.0)
    s.say("(b) 兩個頂點，兩條殼層線相連", wait=1.8)
    s.play(Write(lb), run_time=1.2)
    s.say("剩四隻內區的手：改的是 ", "$u$", "；", "$36=6\\times6$", wait=2.2)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上引的圖號寫錯了：應該是 Fig. 12.11(b)", wait=2.0)
    s.clear_stage()


def part_recursion(s):
    nb = numbers()
    a, b = nb["shell"]
    l1 = mixed("$\\sum_{\\vec Q\\in\\rm shell}1=N\\,d\\,\\delta\\ell$", "（", "$d=4$", "、", "$\\delta\\ell=0.01$", "：",
               "$%s$" % f4(a), " vs ", "$%s$" % f4(b), "）", font_size=26)
    l2 = mixed("(12.26) ", "$\\vec Q+\\vec Q'+\\vec q_1+\\vec q_2=0$", "，取 ", "$\\vec q_1+\\vec q_2=0$",
               "$\\Rightarrow$", "(12.27)", font_size=26)
    l3 = VGroup(tx(r"r+\frac{12u}{r+1}d\,\delta\ell", 32), tx(r"u-\frac{36u^2}{(r+1)^2}d\,\delta\ell", 32)
                ).arrange(RIGHT, buff=1.0)
    l3tag = zh("(12.28)　積掉殼層之後", font_size=22, color=GREY_B)
    col = VGroup(l1, l2, VGroup(l3tag, l3).arrange(DOWN, buff=0.2)).arrange(DOWN, buff=0.4)
    place_in(col, -6.6, 6.6, 0.0, 3.1)
    assert_in_frame(col)
    assert_column(*col)
    s.say("殼層裡有多少個 ", "$\\vec Q$", "？回到上一幕的表", wait=1.8)
    s.play(FadeIn(l1), run_time=1.0)
    s.say("殼層的模態數 ", "$\\approx N\\,d\\,\\delta\\ell$", "：0.0392 對 0.0400", wait=2.0)
    s.play(FadeIn(l2), run_time=1.0)
    s.say("動量守恆：外腳的動量加起來是 0，殼層的和只剩一個", wait=2.2)
    s.play(FadeIn(col[2]), run_time=1.2)
    s.say("(a) 給 ", "$r$", " 加上 ", "$\\frac{12u}{r+1}d\\,\\delta\\ell$", "；(b) 給 ", "$u$", " 減掉一項", wait=2.4)
    s.say("兩項都正比於 ", "$\\delta\\ell$", "：殼層越薄，改變越小", wait=2.0)
    box = D.erratum_box("p. 296–297，(12.36)(12.38)(12.39)", r"(4-d)u-36u^2", r"(4-d)u-36d\,u^2",
                        width=6.4, font_size=30)
    eqs = VGroup(tx(r"\frac{dr}{d\ell}=2r+\frac{12u}{1+r}", 34),
                 tx(r"\frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}", 34)).arrange(DOWN, buff=0.3)
    for m in eqs:
        m.set_color(YELLOW)
    etag = mixed("本系列（", "$d$", " 併進 ", "$u$", "，", "$\\epsilon=4-d$", "）", font_size=22, color=GREY_B)
    right = VGroup(etag, eqs).arrange(DOWN, buff=0.25)
    bottom = VGroup(box, right).arrange(RIGHT, buff=0.7)
    place_in(bottom, -6.6, 6.6, -2.85, -0.3)
    assert_in_frame(bottom)
    assert_apart(col, bottom)
    assert_apart(box, right)
    s.say("再乘上 W1201 的重新標度乘數 ", "$e^{2\\delta\\ell}$", "、", "$e^{(4-d)\\delta\\ell}$", wait=2.2)
    s.play(FadeIn(right), run_time=1.2)
    s.say("得到 (12.37)(12.38)：", "$r$", " 和 ", "$u$", " 隨 ", "$\\ell$", " 的變化率", wait=2.2)
    s.say("$2r$", " 來自 ", "$e^{2\\delta\\ell}$", "；", "$\\epsilon u$", " 來自 ", "$e^{(4-d)\\delta\\ell}$", wait=2.4)
    s.play(FadeIn(box), run_time=0.8)
    s.say("(12.28)、(12.40)、(12.41) 都有 ", "$d$", "；(12.38) 卻漏了", wait=2.2)
    s.say("書上留著 ", "$d$", "，本系列併進 ", "$u$", "：12 和 36 不帶 ", "$d$", wait=2.2)
    s.clear_stage()


def part_fixed(s):
    nb = numbers()
    fp, low = nb["fp"], nb["low"]
    head = [r"\epsilon", r"r^*", r"u^*", r"-\epsilon/6", r"\epsilon/36"]
    rows = [head]
    for e, lab in ((Fraction(1), "1"), (Fraction(1, 2), r"\frac12"), (Fraction(1, 10), "0.1")):
        r, u = fp[e]
        r1, u1 = low[e]
        rows.append([lab, r"%s=%s" % (frac_tex(r), f4(r)), r"%s=%s" % (frac_tex(u), f4(u)), f4(r1), f4(u1)])
    tab = D.tex_table(rows, [0.9, 2.6, 2.6, 1.6, 1.6], row_h=0.65, font_size=28)
    ex = VGroup(tx(r"\text{Gaussian: }(0,0)", 32),
                tx(r"\text{Wilson--Fisher: }r^*=-\frac{\epsilon}{6+\epsilon},\quad u^*=\frac{\epsilon}{(6+\epsilon)^2}", 32)
                ).arrange(DOWN, buff=0.3)
    ex[1].set_color(YELLOW)
    ttag = zh("精確（左）與一階（右）", font_size=22, color=GREY_B)
    col = VGroup(ex, ttag, tab).arrange(DOWN, buff=0.35)
    place_in(col, -6.6, 6.6, -2.85, 3.1)
    assert_in_frame(col)
    assert_column(ex, ttag, tab)
    assert_table(tab, 5)
    s.say("不動點：兩個變化率同時為 0", wait=1.8)
    s.play(FadeIn(ex[0]), run_time=0.8)
    s.say("一個是 Gaussian 點：", "$r=u=0$", wait=1.8)
    s.play(Write(ex[1]), run_time=1.4)
    s.say("另一個是 Wilson–Fisher 點：", "$u^*>0$", " 只在 ", "$\\epsilon>0$", " 時存在", wait=2.4)
    s.play(FadeIn(ttag), FadeIn(tab), run_time=1.2)
    s.say("$\\epsilon=1$", "：", "$r^*=-\\frac17$", "、", "$u^*=\\frac1{49}$", wait=2.0)
    s.say("$\\epsilon$", " 越小，WF 點越靠近 Gaussian 點", wait=2.0)
    s.say("一階近似", "$-\\frac\\epsilon6$", "、", "$\\frac\\epsilon{36}$", "：", "$\\epsilon$", " 越小越準", wait=2.2)
    s.say("$\\epsilon\\to0$", " 時兩個不動點合而為一（書上 ", "$u^*\\to\\epsilon/36d$", "）", wait=2.4)
    s.clear_stage()


def part_flow(s):
    nb = numbers()
    ax = flow_axes()
    ax.move_to(np.array([-2.7, 0.05, 0]))
    curves = VGroup(*[flow_curve(ax, pts, c) for pts, c in zip(nb["lines"], LINE_COLORS)])
    rw, uw = (float(x) for x in nb["fp"][Fraction(1)])
    gdot = Dot(ax.c2p(0, 0), radius=0.1).set_color(WHITE)
    wdot = Dot(ax.c2p(rw, uw), radius=0.1).set_color(WHITE)
    glab = zh("Gaussian", font_size=20, color=GREY_A).next_to(gdot, DR, buff=0.05)
    wlab = zh("WF", font_size=22, color=WHITE).next_to(wdot, RIGHT, buff=0.12)
    g1, g2 = nb["gauss"]
    tag = mixed("$\\epsilon=1$", "，", "$(r,u)$", " 平面的流", font_size=24, color=GREY_A)
    info = VGroup(mixed("WF 點 ", "$(%s,\\ %s)$" % (f4(rw), f4(uw)), font_size=26),
                  mixed("黃線：臨界線 ", "$r_0=r_c(u_0)$", font_size=26),
                  mixed("Gaussian 點的特徵值", font_size=26),
                  tx(r"\lambda_r=%s,\quad\lambda_u=\epsilon=%s" % (f4(g1), f4(g2)), 30).set_color(YELLOW)
                  ).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
    right = VGroup(tag, info).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    place_in(right, 1.5, 6.7, -1.5, 2.6)
    assert_in_frame(ax, curves, right, glab, wlab)
    assert_apart(VGroup(ax, curves, glab, wlab), right)
    s.say("把流畫在 ", "$(r,u)$", " 平面上：", "$\\epsilon=1$", wait=1.8)
    s.play(ShowCreation(ax), FadeIn(tag), run_time=1.0)
    s.play(FadeIn(gdot), FadeIn(wdot), FadeIn(glab), FadeIn(wlab), run_time=0.8)
    s.say("兩個白點：Gaussian 點與 Wilson–Fisher 點", wait=1.8)
    s.play(*[ShowCreation(c[0]) for c in curves[:3]], run_time=2.0)
    s.play(*[GrowArrow(c[1]) for c in curves[:3]], run_time=0.5)
    s.say("從 ", "$u_0=0.002$", " 出發的三條：先靠近 WF 點", wait=2.0)
    s.play(*[ShowCreation(c[0]) for c in curves[3:]], run_time=2.0)
    s.play(*[GrowArrow(c[1]) for c in curves[3:]], run_time=0.5)
    s.say("從 ", "$u_0=0.04$", " 出發的也一樣", wait=1.6)
    s.play(FadeIn(info[0]), FadeIn(info[1]), run_time=0.8)
    s.say("中間那條（臨界線）流進 WF 點", wait=1.8)
    s.say("旁邊兩條先靠近、再往 ", "$r\\to\\pm$", " 跑掉：高溫或低溫", wait=2.2)
    s.say("起點不同，臨界線上的系統最後都到同一點：普適性", wait=2.4)
    s.say("WF 點只有一個不穩定方向：只要調一個參數（溫度）", wait=2.4)
    s.play(FadeIn(info[2]), FadeIn(info[3]), run_time=0.8)
    s.say("Gaussian 點：兩個特徵值都是正的", wait=1.8)
    s.say("$\\epsilon>0$", " 時 ", "$u$", " 方向也是相關的——Gaussian 點不穩定", wait=2.4)
    s.say("四維以下，臨界的系統都流向 Wilson–Fisher 點", wait=2.2)
    s.say("——臨界指數由它決定", wait=1.8)


class W1203OneLoopFlow(ChuScene):
    def construct(self):
        self.show_header("W1203　一迴圈與 (r, u) 平面的流")
        check_math()
        part_diagrams(self)
        part_recursion(self)
        part_fixed(self)
        part_flow(self)
        self.wait(1.0)
