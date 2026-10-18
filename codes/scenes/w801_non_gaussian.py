"""
W801 長尾巴的步伐：不再是高斯（§8.1–8.1.1，pp. 193–195，(8.1)–(8.7)，Ex. 8.1）

問：如果醉漢偶爾來一個超大步，走 1000 步之後，還是一團高斯的霧嗎？
  1. 對比足跡（WS.trail(1000, None) 與 WS.trail(1000, 1)，種子 801）；(8.1)(8.2) 截到 R 的二階動差 9、99、999。
  2. (8.3)(8.4) χ ≈ A − B q^{2−δ}：B = π/4 = 0.7854；(1−χ)/q = 0.7687、0.7837、0.7852；log–log 斜率 1.00、1.49、參考線 2。
  3. 寬度表：長尾巴 18、172、1822（除以 N → 1.7785）；步長 1：2.86、8.98、28.06（除以 √N → 0.8881）；斜率 1.01、0.50。
  4. (8.5) 尾巴 r^{−5+δ}：比值 0.9803 … 1.0051、局部斜率 −4.00、−4.51；換變數 (8.6)(8.7)：0.3536 = 2^{−3/2}；Ex. 8.1；筆誤框。
畫面上的數字都由 utils.walk_scale 算出，check_math() 核對。模組層級的畫圖函式：trail_box（足跡框＋比例尺）、loglog_axes。
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
from utils import walk_scale as WS
from utils import walk_draw as D

NS = WS.WIDTH_N                       # 10, 100, 1000
RS = (10, 100, 1000)                  # 截斷半徑
QS = (0.1, 0.01, 0.001)
TAIL_R = (10, 30, 100)
SPAN = {None: (40, 10), 1: (1000, 1000)}     # 框代表的邊長、比例尺長度
UNIT_COLOR = BLUE_B
LEVY_COLOR = YELLOW
Y15_COLOR = TEAL_B
LOG_SHIFT_X, LOG_SHIFT_Y = 3.0, 5.5   # ManimGL 的軸畫在 0：平移座標


def f2(x):
    return f"{x:.2f}"


def f4(x):
    return f"{float(x):.4f}"


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = {"m2": [WS.moment2_trunc(1, R) for R in RS], "m2q": [WS.moment2_trunc_quad(1, R) for R in RS],
              "B": WS.B_pareto(1), "Bq": WS.B_pareto_quad(1),
              "ratio": [WS.chi_ratio(q, 1) for q in QS],
              "slope": [WS.chi_slope(1), WS.chi_slope(1.5)],
              "wl": [WS.median_width(N, 1) for N in NS], "wu": [WS.median_width(N, None) for N in NS],
              "pl": WS.width_pred(1, 1), "pu": WS.width_pred(1, None),
              "ws": [WS.width_slope(1), WS.width_slope(None)],
              "t1": [WS.levy_density(r, 1, 1) / WS.levy_tail(r, 1, 1) for r in TAIL_R],
              "t15": [WS.levy_density(r, 1, 1.5) / WS.levy_tail(r, 1, 1.5) for r in TAIL_R],
              "ts": [WS.levy_slope(100, 1, 1), WS.levy_slope(100, 1, 1.5)],
              "sub": abs(WS.subst_integral(2, WS.SUBST_Y, WS.SUBST_EPS)) / abs(WS.subst_integral(1, WS.SUBST_Y, WS.SUBST_EPS)),
              "trail_u": WS.trail(1000, None), "trail_l": WS.trail(1000, 1)}
        xs = np.linspace(-3, -1, 9)
        nb["pts"] = {y: [(x, math.log10(1 - WS.chi_pareto(10 ** x, y))) for x in xs] for y in (1, 1.5)}
    return nb


def check_math():
    nb = numbers()
    assert [f"{v:.0f}" for v in nb["m2"]] == ["9", "99", "999"]
    assert all(abs(a - b) < 1e-6 * max(1, b) for a, b in zip(nb["m2"], nb["m2q"]))
    assert WS.pareto_alpha(1) == 4 and WS.delta_of(1) == 1
    assert f4(nb["B"]) == "0.7854" and abs(nb["B"] - math.pi / 4) < 1e-12 and abs(nb["B"] - nb["Bq"]) < 1e-6
    assert [f4(v) for v in nb["ratio"]] == ["0.7687", "0.7837", "0.7852"]
    assert [f2(v) for v in nb["slope"]] == ["1.00", "1.49"]
    assert [f"{v:.0f}" for v in nb["wl"]] == ["18", "172", "1822"]
    assert [f2(v / N) for v, N in zip(nb["wl"], NS)] == ["1.77", "1.72", "1.82"]
    assert [f2(v) for v in nb["wu"]] == ["2.86", "8.98", "28.06"]
    assert [f2(v / math.sqrt(N)) for v, N in zip(nb["wu"], NS)] == ["0.90", "0.90", "0.89"]
    assert f4(nb["pl"]) == "1.7785" and f4(nb["pu"]) == "0.8881"
    assert [f2(v) for v in nb["ws"]] == ["1.01", "0.50"]
    assert [f4(v) for v in nb["t1"]] == ["0.9803", "0.9978", "0.9998"]
    assert [f4(v) for v in nb["t15"]] == ["1.1755", "1.0316", "1.0051"]
    assert [f2(v) for v in nb["ts"]] == ["-4.00", "-4.51"]
    assert WS.tail_exponent(1) == -4 and WS.tail_exponent(Fraction(1, 2)) == Fraction(-9, 2)
    assert f4(nb["sub"]) == "0.3536" and f4(2 ** -1.5) == "0.3536"
    for key, y in (("trail_u", None), ("trail_l", 1)):
        t = nb[key]
        assert t.shape == (1001, 3)
        span = max(np.ptp(t[:, 0]), np.ptp(t[:, 1]))
        assert span <= SPAN[y][0], span
    for y in (1, 1.5):
        for x, v in nb["pts"][y]:
            assert 0 <= x + LOG_SHIFT_X <= 2 and 0 <= v + LOG_SHIFT_Y <= 5


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


def trail_box(pts, span, bar, side, center, color, bar_label):
    """足跡框：pts (n, 2)；框代表 span×span 的區域（足跡的外框置中），框下畫長 bar 的比例尺。
    回傳 VGroup(框, 足跡, 比例尺, 比例尺數字)。足跡是一個 VMobject。"""
    center = np.array(center, dtype=float)
    k = side / span
    mid = (pts.min(axis=0) + pts.max(axis=0)) / 2
    p3 = [center + np.array([(x - mid[0]) * k, (y - mid[1]) * k, 0.0]) for x, y in pts]
    path = VMobject().set_points_as_corners(p3).set_stroke(color, 1.3)
    frame = Square(side_length=side).move_to(center).set_stroke(GREY_B, 1.5)
    y0 = center[1] - side / 2 - 0.25
    x0 = center[0] - side / 2
    sb = VGroup(Line([x0, y0, 0], [x0 + bar * k, y0, 0]).set_stroke(WHITE, 3),
                Line([x0, y0 - 0.08, 0], [x0, y0 + 0.08, 0]).set_stroke(WHITE, 2),
                Line([x0 + bar * k, y0 - 0.08, 0], [x0 + bar * k, y0 + 0.08, 0]).set_stroke(WHITE, 2))
    lab = tx(bar_label, 26).next_to(sb, DOWN, buff=0.1)
    return VGroup(frame, path, sb, lab)


def loglog_axes(width=5.6, height=4.2):
    ax = D.plot_axes((0, 2.2, 1), (0, 5.3, 1), width=width, height=height, x_ticks=(0, 1, 2),
                     x_tick_labels=("-3", "-2", "-1"), y_ticks=(0.5, 1.5, 2.5, 3.5, 4.5),
                     y_tick_labels=("-5", "-4", "-3", "-2", "-1"), font_size=22)
    xl = tx(r"\log_{10}q", 26, GREY_B).next_to(ax.x_axis.get_end(), DR, buff=0.1)
    yl = tx(r"\log_{10}(1-\chi)", 26, GREY_B).next_to(ax.y_axis.get_end(), RIGHT, buff=0.15)
    return ax, VGroup(xl, yl)


def lp(ax, x, v):
    return ax.c2p(x + LOG_SHIFT_X, v + LOG_SHIFT_Y)


# ---------------------------------------------------------------- 各段
def part_trails(s):
    nb = numbers()
    left = trail_box(nb["trail_u"][:, :2], *SPAN[None], 3.7, [-3.3, 0.15, 0], UNIT_COLOR, "10")
    right = trail_box(nb["trail_l"][:, :2], *SPAN[1], 3.7, [3.3, 0.15, 0], LEVY_COLOR, "1000")
    tl = zh("步長都是 1", font_size=28, color=UNIT_COLOR).next_to(left[0], UP, buff=0.2)
    tr = mixed("長尾巴：", "$P(\\ell>r)=r^{-1}$", font_size=28, color=LEVY_COLOR).next_to(right[0], UP, buff=0.2)
    assert_in_frame(left, right, tl, tr)
    assert_apart(left, right)
    assert_apart(tl, left)
    assert_apart(tr, right)
    s.say("問：醉漢偶爾來一個超大步，走 1000 步之後……", wait=0.5)
    s.play(FadeIn(left[0]), FadeIn(right[0]), FadeIn(tl), FadeIn(tr), run_time=0.8)
    s.say("還是一團高斯的霧嗎？（三維足跡，畫前兩個座標）", wait=0.3)
    s.play(ShowCreation(left[1]), ShowCreation(right[1]), run_time=4.0, rate_func=linear)
    s.play(FadeIn(left[2:]), FadeIn(right[2:]), run_time=0.8)
    s.say("左邊：一團霧，寬度大約 ", "$\\sqrt N$", wait=3.0)
    s.say("右邊：幾團霧，被少數幾個長跳連起來", wait=3.0)
    s.say("比例尺差了 100 倍：右邊的框大得多", wait=3.0)
    s.clear_stage()


def part_moment(s):
    nb = numbers()
    e1 = tx(r"p(\vec r)\propto r^{-\alpha},\qquad \langle r^2\rangle=\int p(\vec r)\,r^2\,d^3r"
            r"\propto\int r^{2}\,r^{-\alpha}\,r^2\,dr=\int r^{4-\alpha}dr", 34)
    e2 = mixed("主例 ", "$\\alpha=4$", "（", "$y=1$", "）：截到 ", "$R$", "，", "$\\langle r^2\\rangle_R=R-1$",
               font_size=30)
    tab = D.tex_table([["R", "10", "100", "1000"],
                       [r"\langle r^2\rangle_R"] + [f"{v:.0f}" for v in nb["m2"]]],
                      [2.0, 1.6, 1.6, 1.6], row_h=0.66, font_size=32)
    for j in range(1, 4):
        D.table_cell(tab, 1, j, 4).set_color(YELLOW)
    concl = zh("R 越大越大：二階動差不存在", font_size=30, color=RED_B)
    place_in(e1, -6.7, 6.7, 2.2, 3.1)
    place_in(e2, -6.7, 6.7, 1.2, 1.9)
    place_in(tab, -6.7, 6.7, -0.9, 0.9)
    place_in(concl, -6.7, 6.7, -2.4, -1.6)
    assert_in_frame(e1, e2, tab, concl)
    assert_column(e1, e2, tab, concl)
    assert_table(tab, 4)
    s.say("(8.1)：步伐的密度是冪次 ", "$r^{-\\alpha}$", wait=0.3)
    s.play(Write(e1), run_time=1.5)
    s.say("二階動差 (8.2)：被積函數 ", "$r^{4-\\alpha}$", wait=2.5)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("$\\alpha=4$", "：", "$\\int dr$", " 線性發散，截到 ", "$R$", " 看看", wait=1.5)
    s.play(FadeIn(tab[0][:4]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][4:]], lag_ratio=0.25), run_time=1.0)
    s.say("$R=10,\\ 100,\\ 1000$", "：", "$\\langle r^2\\rangle_R=9,\\ 99,\\ 999$", wait=2.5)
    s.play(FadeIn(concl), run_time=0.7)
    s.say("沒有有限的變異數：中央極限定理的前提不成立", wait=3.0)
    s.clear_stage()


def part_chi(s):
    nb = numbers()
    e1 = tx(r"\chi(\vec q)\approx A-B|\vec q|^{2-\delta}", 36)
    e2 = tx(r"\delta=1:\ \ B=\frac{\pi}{4}=%s" % f4(nb["B"]), 32)
    tab = D.tex_table([["q"] + [str(q) for q in QS] + [r"\to"],
                       [r"\frac{1-\chi(q)}{q}"] + [f4(v) for v in nb["ratio"]] + [f4(nb["B"])]],
                      [1.6, 1.2, 1.2, 1.2, 1.2], row_h=0.75, font_size=26)
    D.table_cell(tab, 1, 4, 5).set_color(YELLOW)
    left = VGroup(e1, e2, tab).arrange(DOWN, buff=0.4)
    place_in(left, -6.7, -0.4, -1.0, 3.0)
    ax, axl = loglog_axes()
    graph = VGroup(ax, axl)
    place_in(graph, 0.0, 6.7, -2.9, 3.1)
    dots1 = VGroup(*[Dot(lp(ax, x, v), radius=0.06).set_color(YELLOW) for x, v in nb["pts"][1]])
    dots15 = VGroup(*[Dot(lp(ax, x, v), radius=0.06).set_color(Y15_COLOR) for x, v in nb["pts"][1.5]])
    ln1 = Line(dots1[0].get_center(), dots1[-1].get_center()).set_stroke(YELLOW, 2)
    ln15 = Line(dots15[0].get_center(), dots15[-1].get_center()).set_stroke(Y15_COLOR, 2)
    v_end = nb["pts"][1][-1][1]
    ref = DashedLine(lp(ax, -3, v_end - 4), lp(ax, -1, v_end)).set_stroke(GREY_B, 2)
    l1 = mixed("$\\delta=1$", "：斜率 ", "$%s$" % f2(nb["slope"][0]), font_size=24, color=YELLOW)
    l15 = mixed("$\\delta=\\frac12$", "：斜率 ", "$%s$" % f2(nb["slope"][1]), font_size=24, color=Y15_COLOR)
    lr = zh("高斯步伐：斜率 2", font_size=24, color=GREY_B)
    legend = VGroup(l1, l15, lr).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
    place_in(legend, -6.7, -0.4, -2.9, -1.4)
    assert_in_frame(left, graph, legend)
    assert_apart(left, graph)
    assert_apart(left, legend)
    assert_apart(graph, legend)
    assert_table(tab, 5)
    s.say("(8.3) 結構函數：", "$\\chi(\\vec q)\\approx A-B|\\vec q|^{2-\\delta}$", wait=0.3)
    s.play(Write(e1), run_time=1.2)
    s.say("不是 ", "$q^2$", "：小 ", "$q$", " 的展開斷在分數次方", wait=2.5)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("主例 ", "$\\delta=1$", "：(8.4) 的係數 ", "$B=\\pi/4$", wait=1.5)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][5:]], lag_ratio=0.2), run_time=1.0)
    s.say("數值積分：0.7687、0.7837、0.7852 → 0.7854", wait=3.0)
    s.play(ShowCreation(ax), FadeIn(axl), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(d) for d in dots1], lag_ratio=0.1), ShowCreation(ln1), run_time=1.2)
    s.play(LaggedStart(*[FadeIn(d) for d in dots15], lag_ratio=0.1), ShowCreation(ln15), run_time=1.2)
    s.play(ShowCreation(ref), FadeIn(legend), run_time=1.0)
    s.say("log–log 斜率：", "$1.00$", " 和 ", "$1.49$", "，就是 ", "$2-\\delta$", wait=3.0)
    s.say("高斯步伐的斜率是 2：小 ", "$q$", " 時長尾巴的 ", "$1-\\chi$", " 大得多", wait=3.0)
    s.clear_stage()


def part_width(s):
    nb = numbers()
    rows = [[zh("N", font_size=26, color=YELLOW), "10", "100", "1000", zh("理論", font_size=26, color=YELLOW),
             zh("斜率", font_size=26, color=YELLOW)],
            [zh("長尾巴 中位數", font_size=24, color=LEVY_COLOR)] + [f"{v:.0f}" for v in nb["wl"]] + ["", f2(nb["ws"][0])],
            [zh("除以 N", font_size=24, color=LEVY_COLOR)] + [f2(v / N) for v, N in zip(nb["wl"], NS)] + [f4(nb["pl"]), ""],
            [zh("步長 1 中位數", font_size=24, color=UNIT_COLOR)] + [f2(v) for v in nb["wu"]] + ["", f2(nb["ws"][1])],
            [zh("除以 √N", font_size=24, color=UNIT_COLOR)] + [f2(v / math.sqrt(N)) for v, N in zip(nb["wu"], NS)]
            + [f4(nb["pu"]), ""]]
    tab = D.tex_table(rows, [2.6, 1.4, 1.4, 1.4, 1.6, 1.3], row_h=0.66, font_size=30)
    for i, j in ((2, 4), (4, 4), (1, 5), (3, 5)):
        D.table_cell(tab, i, j, 6).set_color(YELLOW)
    head = mixed("各 2000 人走 ", "$N$", " 步（三維），", "$|\\vec R_N|$", " 的中位數", font_size=28, color=GREY_B)
    concl = mixed("寬度 ", "$\\propto N^{1/y}$", "（", "$y=1$", "：", "$\\propto N$", "），不是 ", "$\\sqrt N$",
                  font_size=32, color=YELLOW)
    place_in(head, -6.7, 6.7, 2.5, 3.1)
    place_in(tab, -6.7, 6.7, -1.4, 2.2)
    place_in(concl, -6.7, 6.7, -2.6, -1.8)
    assert_in_frame(head, tab, concl)
    assert_column(head, tab, concl)
    assert_table(tab, 6)
    s.say("寬度怎麼長？各 2000 人，量 ", "$|\\vec R_N|$", " 的中位數", wait=0.3)
    s.play(FadeIn(head), FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:18]], lag_ratio=0.08), run_time=1.2)
    s.say("長尾巴：18、172、1822，每次大約 10 倍", wait=2.5)
    s.say("除以 ", "$N$", "：1.77、1.72、1.82，理論值 1.7785", wait=3.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][18:]], lag_ratio=0.08), run_time=1.2)
    s.say("步長 1：2.86、8.98、28.06，除以 ", "$\\sqrt N$", " 都是 0.90 左右", wait=3.0)
    s.say("log–log 斜率 1.01 對 0.50", wait=2.5)
    s.play(FadeIn(concl), run_time=0.8)
    s.say("寬度長得像 ", "$N^{1/y}$", "：長尾巴走得比 ", "$\\sqrt N$", " 快得多", wait=3.0)
    s.clear_stage()


def part_tail(s):
    nb = numbers()
    e1 = tx(r"C(N;\vec r)\propto|\vec r|^{-5+\delta}", 36)
    tab = D.tex_table([["r", "10", "30", "100", r"\text{slope at }100", r"\to"],
                       [r"\delta=1"] + [f4(v) for v in nb["t1"]] + [f2(nb["ts"][0]), "-4"],
                       [r"\delta=\frac12"] + [f4(v) for v in nb["t15"]] + [f2(nb["ts"][1]), "-4.5"]],
                      [1.4, 1.5, 1.5, 1.5, 2.3, 1.2], row_h=0.72, font_size=30)
    for i in (1, 2):
        D.table_cell(tab, i, 5, 6).set_color(YELLOW)
    lab = zh("數值逆 Fourier 除以 換變數的漸近式（s = 1）", font_size=26, color=GREY_B)
    top = VGroup(e1, lab, tab).arrange(DOWN, buff=0.3)
    place_in(top, -6.7, 6.7, -0.3, 3.1)
    e2 = tx(r"\int_0^\infty e^{ikr}k^y\,dk\ \overset{w=kr}{=}\ r^{-(y+1)}\int_0^\infty e^{iw}w^y\,dw", 32)
    e3 = tx(r"y=\frac12:\ \ \frac{|I(2)|}{|I(1)|}=%s=2^{-3/2}" % f4(nb["sub"]), 30)
    e4 = mixed("Ex. 8.1：", "$\\frac{4\\pi}{r}\\int q\\sin(qr)\\,q^{2-\\delta}dq\\ \\to\\ r^{-1}\\cdot r^{-(4-\\delta)}=r^{-5+\\delta}$",
               font_size=28)
    bot = VGroup(e2, e3, e4).arrange(DOWN, buff=0.25)
    place_in(bot, -6.7, 6.7, -2.9, -0.5)
    box = D.erratum_box("p. 195，(8.6)→(8.7)", r"w=k/r", r"w=kr", width=5.0, font_size=40)
    place_in(box, -6.7, 6.7, 0.6, 2.6)
    assert_in_frame(top, bot, box)
    assert_apart(top, bot)
    assert_apart(box, bot)
    assert_table(tab, 6)
    s.say("尾巴 (8.5)：", "$C(N;\\vec r)\\propto|\\vec r|^{-5+\\delta}$", "，冪次而不是 ", "$e^{-r^2}$", wait=0.3)
    s.play(Write(e1), FadeIn(lab), run_time=1.2)
    s.play(FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:]], lag_ratio=0.08), run_time=1.5)
    s.say("$\\delta=1$", "：0.9803、0.9978、0.9998；", "$\\delta=\\frac12$", " 慢一點也趨近 1", wait=3.0)
    s.say("$r=100$", " 的局部斜率 ", "$-4.00$", "、", "$-4.51$", "：就是 ", "$-5+\\delta$", wait=3.0)
    s.play(FadeIn(bot[0]), run_time=1.0)
    s.say("換變數 (8.6)(8.7)：令 ", "$w=kr$", "，", "$r$", " 整個提到積分外面", wait=3.0)
    s.play(FadeIn(bot[1]), run_time=0.8)
    s.say("$r$", " 加倍，積分變成 ", "$2^{-3/2}=0.3536$", " 倍", wait=2.5)
    s.play(FadeIn(bot[2]), run_time=0.8)
    s.say("Ex. 8.1：三維逆 Fourier 多一個 ", "$1/r$", "，合起來 ", "$r^{-5+\\delta}$", wait=3.0)
    s.play(FadeOut(top), run_time=0.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上寫 ", "$w=k/r$", "；應該是 ", "$w=kr$", wait=2.5)
    s.say("照書上 ", "$k=rw$", "，積分變成 ", "$r^{y+1}\\int e^{ir^2w}w^ydw$", "，", "$r$", " 提不出來", wait=3.0)
    s.say("二階動差一發散，中央極限定理就不管用", wait=2.5)
    s.say("寬度長得像 ", "$N^{1/y}$", "，尾巴是冪次 ", "$r^{-5+\\delta}$", wait=3.0)


class W801NonGaussian(ChuScene):
    def construct(self):
        self.show_header("W801　長尾巴的步伐：不再是高斯")
        check_math()
        part_trails(self)
        part_moment(self)
        part_chi(self)
        part_width(self)
        part_tail(self)
        self.wait(1.0)
