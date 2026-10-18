"""
W704 生成函數與 Dyson 方程（§7.2.3–7.2.4，pp. 176–184，(7.20)–(7.37)，Fig. 7.5–7.10，Ex. 7.3）

問：展開的每一項都是「切段再求和」的摺積，越高階越難算。有沒有辦法把摺積變成乘法？
  1. (7.20)–(7.22) 生成函數把摺積變乘積：H1 = H2 = 2^K → 摺積 1, 4, 12, 32, 80 = (N+1)2^N ↔ 1/(1−2z)²；z = 0.1：1.5625。
     筆誤框 1：(7.21) 第二個 z 寫成大寫 Z。
  2. 費曼規則：實線 g0、虛線 −F、頂點動量守恆；一階 (7.26)；並排兩圈 = g0 Σ1 g0 Σ1 g0（書上 p. 178 寫成 (a)）。
  3. 火車車廂 (7.30)(7.36)：g0 = 2、Σ = 1/4 的部分和 2, 3, 7/2, 15/4, 31/8 → 4；筆誤框 2：(7.35) 分母的號（4/3），(7.33) 照抄得 10/3。
  4. 自能＝不可約（1PI）的圖；計數版的 Dyson：1, 3, 15, 105 = 1/(1 − (x + 2x² + 10x³ + 74x⁴)) 的係數。
畫面上的數字都由 utils.walk_saw 算出，check_math() 核對。模組層級的畫圖函式：chord_diagram（弧線圖／費曼圖）、train（車廂）。
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
from utils import walk_draw as D

ARC_COLOR = YELLOW
CAR_COLOR = TEAL_C
BAR_COLOR = BLUE_C


def f4(x):
    return f"{x:.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    g0, sg = SW.DYSON_TOY
    return {"H": list(SW.EX_CONV["H"]), "conv": SW.convolve(SW.EX_CONV["H"], SW.EX_CONV["H"]),
            "gf": SW.gf_check(0.1), "fourier": SW.fourier_conv_error(),
            "part": SW.dyson_partial(2, Fraction(1, 4), 5), "closed": SW.dyson_closed(g0, sg),
            "b35": SW.dyson_book35(g0, sg), "b33": SW.dyson_book33(g0, sg), "counts": SW.dyson_counts(4)}


def check_math():
    nb = numbers()
    g0, sg = SW.DYSON_TOY
    assert (g0, sg) == (2, Fraction(1, 4)) and sg * g0 == Fraction(1, 2)
    assert nb["H"] == [1, 2, 4, 8, 16] == [2 ** k for k in range(5)]
    assert nb["conv"] == SW.EX_CONV["conv"] == [1, 4, 12, 32, 80] == [(n + 1) * 2 ** n for n in range(5)]
    assert [nb["H"][j] * nb["H"][2 - j] for j in range(3)] == [4, 4, 4] and nb["conv"][2] == 12   # 1·4+2·2+4·1
    s, closed = nb["gf"]
    assert f4(s) == f4(closed) == "1.5625" and abs(closed - 1 / 0.8 ** 2) < 1e-12
    assert nb["fourier"] < 1e-9
    part = nb["part"]
    assert part == [2, 3, Fraction(7, 2), Fraction(15, 4), Fraction(31, 8)]
    assert [part[0]] + [part[n] - part[n - 1] for n in range(1, 5)] == [2 * Fraction(1, 2) ** n for n in range(5)]
    assert nb["closed"] == 4 == 1 / (Fraction(1, 2) - Fraction(1, 4))
    assert nb["b35"] == Fraction(4, 3) and nb["b33"] == Fraction(10, 3)
    full, irr, ser = nb["counts"]
    assert full[1:] == [1, 3, 15, 105] == ser[1:] and irr[1:] == [1, 2, 10, 74]
    so = SW.SECOND_ORDER
    assert [so[k][0] for k in "abc"] == ["nested", "sequential", "crossed"]
    assert [SW.is_irreducible(so[k][1]) for k in "abc"] == [True, False, True]
    assert sorted(SW.matchings(2)) == sorted(so[k][1] for k in "abc")
    assert list(SW.matchings(1)) == [((1, 2),)]


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


def arrow_tip(p, color=WHITE, size=0.16):
    """指向右邊的小三角形（傳播子的箭頭）。"""
    t = Triangle().set_fill(color, 1).set_stroke(width=0).set_height(size)
    return t.rotate(-PI / 2).move_to(p)


def chord_diagram(m, width=3.0, origin=ORIGIN, color=ARC_COLOR, arrows=False, dashed=True, angle=-0.75 * PI):
    """弧線圖：水平實線（漫步）上 2k 個點（編號 1..2k），每對 (a, b) 一條上方的弧（交互作用）。
    arrows=True 時每段實線加箭頭（費曼圖）。回傳 VGroup(線, 點們, 弧們, 箭頭們)；點的座標存在 .pts。"""
    origin = np.array(origin, dtype=float)
    n = 2 * len(m)
    gap = width / (n + 1)
    pts = [origin + (-width / 2 + gap * (i + 1)) * RIGHT for i in range(n)]
    line = Line(origin + width / 2 * LEFT, origin + width / 2 * RIGHT).set_stroke(WHITE, 3)
    dots = VGroup(*[Dot(p, radius=0.06).set_color(WHITE) for p in pts])
    arcs = VGroup()
    for a, b in m:
        arc = ArcBetweenPoints(pts[a - 1], pts[b - 1], angle=angle)
        if dashed:
            span = np.linalg.norm(pts[b - 1] - pts[a - 1])
            arc = DashedVMobject(arc, num_dashes=max(5, int(1.3 * span / 0.13)))
        arcs.add(arc.set_stroke(color, 2.5))
    tips = VGroup()
    if arrows:
        ends = [line.get_start()] + pts + [line.get_end()]
        tips = VGroup(*[arrow_tip((p + q) / 2) for p, q in zip(ends, ends[1:])])
    g = VGroup(line, dots, arcs, tips)
    g.pts = [p.copy() for p in pts]
    g.ends = (line.get_start().copy(), line.get_end().copy())
    return g


def train(ncars, origin, car_w=0.72, link=0.42, color=CAR_COLOR):
    """火車：g0（帶箭頭的連結）、Σ（車廂）交替；origin = 左端。回傳 VGroup(依序的零件們)。"""
    x, y = float(origin[0]), float(origin[1])
    parts = VGroup()

    def add_link(x0):
        ln = Line([x0, y, 0], [x0 + link, y, 0]).set_stroke(WHITE, 3)
        return VGroup(ln, arrow_tip([x0 + link / 2, y, 0], size=0.13))

    parts.add(add_link(x))
    x += link
    for _ in range(ncars):
        body = RoundedRectangle(width=car_w, height=0.4, corner_radius=0.08)
        body.set_fill(color, 0.45).set_stroke(color, 2).move_to([x + car_w / 2, y, 0])
        sig = tx(r"\Sigma", 24).move_to(body)
        wheels = VGroup(*[Circle(radius=0.06).set_fill(GREY_B, 1).set_stroke(GREY_B, 1)
                          .move_to([x + car_w * f, y - 0.26, 0]) for f in (0.25, 0.75)])
        parts.add(VGroup(body, sig, wheels))
        x += car_w
        parts.add(add_link(x))
        x += link
    return parts


# ---------------------------------------------------------------- 各段
def part_question(s):
    eq = tx(r"\sum_{M_1+M_2+M_3=N}\ \sum_{\vec u_1,\vec u_2}C(M_1;\vec x,\vec u_1)\,C(M_2;\vec u_1,\vec u_2)\,"
            r"C(M_3;\vec u_2,\vec y)\,f(|\vec u_1-\vec u_2|)", 36)
    dg = chord_diagram(((1, 2),), width=6.0, origin=(0, 0.1, 0), arrows=True)
    labs = VGroup()
    for k, p in enumerate(dg.pts):
        labs.add(tx(r"\vec u_%d" % (k + 1), 26, GREY_B).next_to(p, DOWN, buff=0.18))
    for t, p in zip((r"\vec x", r"\vec y"), dg.ends):
        labs.add(tx(t, 26, GREY_B).next_to(p, DOWN, buff=0.18))
    seg_x = [dg.ends[0], dg.pts[0], dg.pts[1], dg.ends[1]]
    for k in range(3):
        labs.add(tx("M_%d" % (k + 1), 26, BLUE_B).move_to((seg_x[k] + seg_x[k + 1]) / 2 + 0.38 * DOWN))
    flab = tx("f", 28, ARC_COLOR).next_to(dg[2], UP, buff=0.1)
    pic = VGroup(dg, labs, flab)
    bot = mixed("摺積 ", "$\\sum_{M}H_1(M)\\,H_2(N-M)$", "：能不能變成乘法？", font_size=30, color=YELLOW)
    place_in(eq, -6.7, 6.7, 2.0, 3.1)
    place_in(pic, -4.0, 4.0, -0.7, 1.5)
    place_in(bot, -6.7, 6.7, -2.8, -1.4)
    assert_in_frame(eq, pic, bot)
    assert_column(eq, pic, bot)
    s.say("問：展開的每一項都是「切段再求和」的摺積", wait=0.5)
    s.play(FadeIn(eq), run_time=1.0)
    s.play(ShowCreation(dg[0]), FadeIn(dg[1]), FadeIn(dg[3]), FadeIn(labs), run_time=1.0)
    s.play(ShowCreation(dg[2]), FadeIn(flab), run_time=0.8)
    s.say("一階 (7.19)：三段漫步接起來，中間點還要求和", wait=2.5)
    s.say("越高階，切的段越多、求和越多層", wait=2.0)
    s.play(FadeIn(bot), run_time=0.8)
    s.say("有沒有辦法把摺積變成乘法？（先猜）", wait=3.0)
    s.clear_stage()


def part_gf(s):
    nb = numbers()
    top = VGroup(tx(r"\sum_N z^N\sum_{M}H_1(M)\,H_2(N-M)=h_1(z)\,h_2(z)", 34),
                 tx(r"h(z)=\sum_K z^K H(K)", 30, GREY_B)).arrange(RIGHT, buff=0.9)
    rows = [["N", "0", "1", "2", "3", "4"],
            [r"H_1(N)=H_2(N)=2^N"] + [str(v) for v in nb["H"]],
            [r"\sum_M H_1(M)\,H_2(N-M)"] + [str(v) for v in nb["conv"]]]
    tab = D.tex_table(rows, [4.6, 1.0, 1.0, 1.0, 1.0, 1.0], row_h=0.66, font_size=30)
    note = tx(r"=(N+1)\,2^N", 30, YELLOW)
    gfl = tx(r"h_1(z)=h_2(z)=\frac{1}{1-2z}\ \Rightarrow\ \sum_N (N+1)\,2^N z^N=\frac{1}{(1-2z)^2}", 32)
    zl = mixed("$z=0.1$", "：級數＝", "$\\frac{1}{0.8^2}=%s$" % f4(nb["gf"][1]), font_size=30, color=YELLOW)
    place_in(top, -6.7, 6.7, 2.3, 3.15)
    place_in(tab, -6.7, 3.2, -0.1, 2.0)
    note.next_to(D.table_cell(tab, 2, 5, 6), RIGHT, buff=0.45)
    place_in(gfl, -6.7, 6.7, -1.45, -0.45)
    place_in(zl, -6.7, 6.7, -2.85, -1.75)
    assert_in_frame(top, tab, note, gfl, zl)
    assert_column(top, tab, gfl, zl)
    assert_apart(tab, note)
    assert_table(tab, 6)
    s.say("生成函數：把數列 ", "$H(K)$", " 變成冪級數 ", "$h(z)$", wait=0.5)
    s.play(FadeIn(top[1]), run_time=0.8)
    s.play(Write(top[0]), run_time=1.5)
    s.say("(7.22)：摺積的生成函數＝兩個生成函數相乘", wait=2.5)
    s.play(FadeIn(tab[0][:12]), ShowCreation(tab[1]), run_time=1.0)
    s.say("例：一維 ", "$K$", " 步的走法數 ", "$2^K$", "：1, 2, 4, 8, 16", wait=2.5)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][12:]], lag_ratio=0.2), run_time=1.2)
    box = SurroundingRectangle(D.table_cell(tab, 2, 3, 6), buff=0.08).set_stroke(YELLOW, 2)
    s.play(ShowCreation(box), run_time=0.5)
    s.say("例如 ", "$N=2$", "：", "$1\\cdot4+2\\cdot2+4\\cdot1=12$", wait=2.5)
    s.play(FadeIn(note), FadeOut(box), run_time=0.7)
    s.say("摺積 1, 4, 12, 32, 80 剛好是 ", "$(N+1)\\,2^N$", wait=2.5)
    s.play(FadeIn(gfl), run_time=1.0)
    s.say("另一邊：", "$h(z)=\\frac{1}{1-2z}$", "，兩個相乘", wait=2.0)
    s.play(FadeIn(zl), run_time=0.8)
    s.say("$z=0.1$", "：級數一項項加起來，就是 1.5625", wait=2.5)
    s.play(FadeOut(tab), FadeOut(note), FadeOut(gfl), FadeOut(zl), run_time=0.7)
    fl = mixed("Fourier 也一樣 (7.24)：", "$\\sum_{\\vec r}e^{i\\vec q\\cdot\\vec r}\\sum_{\\vec u}A(\\vec r-\\vec u)\\,"
               "B(\\vec u)=\\tilde A(\\vec q)\\,\\tilde B(\\vec q)$", font_size=30)
    box1 = D.erratum_box("p. 176，(7.21)", r"z^{M_1}H_1(M_1)\,Z^{N-M_1}H_2(N-M_1)",
                         r"z^{M_1}H_1(M_1)\,z^{N-M_1}H_2(N-M_1)", width=13.2)
    place_in(fl, -6.7, 6.7, 0.6, 1.8)
    place_in(box1, -6.7, 6.7, -2.8, 0.1)
    assert_in_frame(fl, box1)
    assert_column(top, fl, box1)
    s.play(FadeIn(fl), run_time=0.8)
    s.say("Fourier 變換也一樣：摺積變乘積", wait=1.5)
    s.play(FadeIn(box1), run_time=0.8)
    s.say("書上 (7.21)：第二個 ", "$z$", " 寫成大寫 ", "$Z$", wait=2.0)
    s.say("大寫 ", "$Z$", " 成了另一個變數，乘積就拆不開", wait=2.5)
    s.clear_stage()


def part_feynman(s):
    leg_prop = VGroup(Line(ORIGIN, 1.1 * RIGHT).set_stroke(WHITE, 3), arrow_tip(0.55 * RIGHT))
    leg_arc = DashedVMobject(ArcBetweenPoints(ORIGIN, 1.1 * RIGHT, angle=-0.75 * PI), num_dashes=8)
    leg_arc.set_stroke(ARC_COLOR, 2.5)
    leg_dot = VGroup(Line(ORIGIN, 1.1 * RIGHT).set_stroke(WHITE, 3), Dot(0.55 * RIGHT, radius=0.08))
    legend = VGroup()
    for icon, lab in ((leg_prop, mixed("傳播子 ", "$g_0=\\frac{1}{1-z\\chi(\\vec q)}$", font_size=28)),
                      (leg_arc, mixed("交互作用 ", "$-F(\\vec Q)$", font_size=28)),
                      (leg_dot, zh("頂點：動量守恆", font_size=28))):
        legend.add(VGroup(icon, lab).arrange(RIGHT, buff=0.35))
    legend.arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    dg = chord_diagram(((1, 2),), width=5.4, origin=ORIGIN, arrows=True)
    seg = [dg.ends[0], dg.pts[0], dg.pts[1], dg.ends[1]]
    qlabs = VGroup(*[tx(t, 26, BLUE_B).move_to((seg[k] + seg[k + 1]) / 2 + 0.38 * DOWN)
                     for k, t in enumerate((r"\vec q", r"\vec Q", r"\vec q"))])
    alab = tx(r"\vec q-\vec Q", 24, ARC_COLOR).next_to(dg[2], UP, buff=0.1)
    first = VGroup(dg, qlabs, alab)
    eq = tx(r"g_0\Big[\sum_{\vec Q}F(\vec q-\vec Q)\,g_0(\vec Q)\Big]g_0=g_0\,\Sigma_1\,g_0", 34)
    dgb = chord_diagram(((1, 2), (3, 4)), width=5.2, origin=ORIGIN, arrows=True)
    bl = zh("(b) 並排", font_size=26, color=GREY_B).next_to(dgb, UP, buff=0.2)
    eqb = tx(r"=g_0\,\Sigma_1\,g_0\,\Sigma_1\,g_0", 34)
    place_in(legend, -6.7, -0.4, 0.75, 3.1)
    place_in(first, 0.6, 6.7, 0.75, 2.9)
    place_in(eq, -6.7, 6.7, -0.65, 0.45)
    low = VGroup(VGroup(dgb, bl), eqb).arrange(RIGHT, buff=0.8)
    place_in(low, -6.7, 6.7, -2.85, -1.05)
    assert_in_frame(legend, first, eq, low)
    assert_apart(legend, first)
    assert_column(legend, eq, low)
    assert_column(first, eq)
    for a, b in zip(legend, legend[1:]):
        assert_apart(a, b)
    s.say("費曼規則：展開的每一項畫成一張圖", wait=0.5)
    s.play(FadeIn(legend[0]), run_time=0.8)
    s.say("實線＝傳播子 ", "$g_0=\\frac{1}{1-z\\chi(\\vec q)}$", wait=2.5)
    s.play(FadeIn(legend[1]), FadeIn(legend[2]), run_time=0.8)
    s.say("虛線＝交互作用 ", "$-F(\\vec Q)$", "；每個頂點動量守恆", wait=2.0)
    s.play(ShowCreation(dg[0]), FadeIn(dg[1]), FadeIn(dg[3]), FadeIn(qlabs), run_time=1.0)
    s.play(ShowCreation(dg[2]), FadeIn(alab), run_time=0.8)
    s.say("一階 (7.26)：進來 ", "$\\vec q$", "，中間那段帶 ", "$\\vec Q$", "，對 ", "$\\vec Q$", " 求和", wait=2.0)
    s.play(Write(eq), run_time=1.5)
    s.say("方括號叫一階自能 ", "$\\Sigma_1$", "：", "$g_0\\,\\Sigma_1\\,g_0$", wait=2.5)
    s.play(ShowCreation(dgb[0]), FadeIn(dgb[1]), FadeIn(dgb[3]), FadeIn(bl), run_time=1.0)
    s.play(ShowCreation(dgb[2]), FadeIn(eqb), run_time=1.0)
    s.say("並排的兩個圈（Fig. 7.3(b)）：每圈各給一個 ", "$\\Sigma_1$", wait=2.5)
    s.say("書上 p. 178 寫成 (a)；(a) 其實是巢狀的", wait=2.0)
    s.clear_stage()


def part_dyson(s):
    nb = numbers()
    part = nb["part"]
    eq = tx(r"g=g_0+g_0\Sigma g_0+g_0\Sigma g_0\Sigma g_0+\cdots=\frac{1}{g_0^{-1}-\Sigma}", 36)
    place_in(eq, -6.7, 6.7, 2.25, 3.15)
    rows, terms = VGroup(), VGroup()
    ys = [1.6, 0.8, 0.0, -0.8, -1.6]
    for n, y in enumerate(ys):
        tr = train(n, np.array([-6.6, y, 0]))
        rows.add(tr)
        terms.add(tx("+" + ftex(part[n] - (part[n - 1] if n else 0)), 28, YELLOW).move_to([-1.15, y, 0]))
    origin = np.array([1.1, -2.2, 0])
    ux, uy = 0.9, 0.85
    bars = D.bar_chart({n + 1: part[n] for n in range(5)}, ux=ux, uy=uy, origin=origin, color=BAR_COLOR)
    base = Line(origin + 0.4 * ux * RIGHT, origin + 5.6 * ux * RIGHT).set_stroke(GREY_B, 2)
    blabs = VGroup(*[tx(ftex(part[n]), 24).next_to(origin + (n + 1) * ux * RIGHT, DOWN, buff=0.12) for n in range(5)])
    lim = DashedLine(origin + 0.4 * ux * RIGHT + 4 * uy * UP, origin + 5.6 * ux * RIGHT + 4 * uy * UP,
                     dash_length=0.1).set_stroke(RED_B, 2.5)
    limlab = tx("4", 30, RED_B).next_to(lim, RIGHT, buff=0.12)
    chart = VGroup(bars, base, blabs, lim, limlab)
    assert_in_frame(eq, rows, terms, chart)
    assert_column(eq, rows)
    assert_column(eq, chart)
    for a, b in zip(rows, rows[1:]):
        assert_apart(a, b)
    for r, t in zip(rows, terms):
        assert_apart(r, t)
    assert_apart(VGroup(rows, terms), chart)
    for a, b in zip(blabs, blabs[1:]):
        assert_apart(a, b)
    s.say("Dyson 方程：同一種車廂一節節接下去", wait=0.5)
    s.play(Write(eq), run_time=1.5)
    s.say("車廂＝自能 ", "$\\Sigma$", "，車廂之間用 ", "$g_0$", " 連起來", wait=1.5)
    s.say("玩具數字：", "$g_0=2$", "、", "$\\Sigma=\\frac14$", "，每多一節乘 ", "$\\Sigma g_0=\\frac12$", wait=1.5)
    s.play(ShowCreation(base), FadeIn(lim), FadeIn(limlab), run_time=0.8)
    for n in range(5):
        s.play(LaggedStart(*[FadeIn(p) for p in rows[n]], lag_ratio=0.15), FadeIn(terms[n]),
               GrowFromEdge(bars[0][n], DOWN), FadeIn(blabs[n]), run_time=1.0)
        if n == 1:
            s.say("部分和：2，再加 1 得 3", wait=1.0)
    s.say("部分和 2, 3, ", "$\\frac72$", ", ", "$\\frac{15}4$", ", ", "$\\frac{31}8$", "，越來越靠近 4", wait=2.5)
    s.say("幾何級數加完：", "$\\frac{1}{1/2-1/4}=4$", "，就是虛線", wait=2.0)
    box2 = D.erratum_box("p. 183，(7.35)", r"\frac{g_0}{1+(\Sigma_1+\Sigma_2)g_0}",
                         r"\frac{g_0}{1-(\Sigma_1+\Sigma_2)g_0}", width=6.0)
    place_in(box2, -6.7, 0.3, -1.6, 1.6)
    assert_in_frame(box2)
    assert_apart(box2, chart)
    s.play(FadeOut(rows), FadeOut(terms), run_time=0.6)
    s.play(FadeIn(box2), run_time=0.8)
    s.say("書上 (7.35)：分母寫成加號", wait=2.0)
    s.say("照書上得 ", "$\\frac43$", "，不是 4", wait=2.5)
    s.say("同頁 (7.33) 少了中間的 ", "$g_0$", "，照抄得 ", "$\\frac{10}3$", wait=2.5)
    s.clear_stage()


def part_irreducible(s):
    nb = numbers()
    full, irr, ser = nb["counts"]
    so = SW.SECOND_ORDER
    tags = {"a": "(a) 巢狀 (7.34)", "b": "(b) 並排", "c": "(c) 交叉 (Fig. 7.9)"}
    panels = VGroup()
    cuts = VGroup()
    for key in "abc":
        dg = chord_diagram(so[key][1], width=3.4, arrows=True)
        title = zh(tags[key], font_size=24, color=GREY_B).next_to(dg, UP, buff=0.18)
        ok = SW.is_irreducible(so[key][1])
        verdict = zh("不可約" if ok else "可約：剪一刀就分開", font_size=24,
                     color=GREEN if ok else RED_B).next_to(dg, DOWN, buff=0.22)
        panels.add(VGroup(dg, title, verdict))
    panels.arrange(RIGHT, buff=0.9)
    rows = [["k", "1", "2", "3", "4"],
            [mixed("全部弧線圖 ", "$(2k-1)!!$", font_size=26)] + [str(v) for v in full[1:]],
            [mixed("不可約的圖 ", "$I_k$", font_size=26)] + [str(v) for v in irr[1:]],
            [mixed("$1/(1-I(x))$", " 的係數", font_size=26)] + [str(v) for v in ser[1:]]]
    tab = D.tex_table(rows, [4.0, 1.1, 1.1, 1.1, 1.1], row_h=0.58, font_size=28)
    for j in range(1, 5):
        D.table_cell(tab, 2, j, 5).set_color(GREEN)
        D.table_cell(tab, 3, j, 5).set_color(YELLOW)
    form = tx(r"\frac{1}{1-(x+2x^2+10x^3+74x^4)}=1+x+3x^2+15x^3+105x^4+\cdots", 32)
    place_in(panels, -6.7, 6.7, 0.95, 3.15)
    for p, key in zip(panels, "abc"):
        if not SW.is_irreducible(so[key][1]):
            c = p[0][0].get_center().copy()                 # 並排：剪在第 2、3 點中間（線的正中央）
            cuts.add(Line(c + 0.3 * DOWN, c + 0.3 * UP).set_stroke(RED_B, 4))
    assert len(cuts) == 1
    place_in(tab, -6.7, 6.7, -1.65, 0.75)
    place_in(form, -6.7, 6.7, -2.85, -1.85)
    assert_in_frame(panels, tab, form)
    assert_column(panels, tab, form)
    assert_table(tab, 5)
    for a, b in zip(panels, panels[1:]):
        assert_apart(a, b)
    for p in panels:
        assert_apart(p[1], p[0])
        assert_apart(p[0], p[2])
    s.say("自能 ", "$\\Sigma$", " 只收「不可約」的圖（1PI）", wait=0.5)
    for p in panels:
        s.play(ShowCreation(p[0][0]), FadeIn(p[0][1]), FadeIn(p[0][3]), FadeIn(p[1]), run_time=0.7)
        s.play(ShowCreation(p[0][2]), run_time=0.6)
    s.say("不可約：剪斷任何一條實線，圖都還連在一起", wait=2.0)
    s.play(FadeIn(panels[0][2]), FadeIn(panels[2][2]), run_time=0.7)
    s.say("二階：巢狀 (7.34) 與交叉 (Fig. 7.9) 都不可約", wait=2.5)
    s.play(ShowCreation(cuts), FadeIn(panels[1][2]), run_time=0.8)
    s.say("並排的剪成兩個一階：它已經在車廂裡了", wait=2.5)
    s.play(FadeIn(tab[0][:10]), ShowCreation(tab[1]), run_time=1.0)
    s.say("計數版的 Dyson：全部的圖＝不可約的圖排成一串", wait=2.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][10:15]], lag_ratio=0.2), run_time=1.0)
    s.say("不可約的個數：1, 2, 10, 74", wait=2.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][15:]], lag_ratio=0.2), FadeIn(form), run_time=1.2)
    s.say("代進 ", "$\\frac{1}{1-I(x)}$", "：展開剛好是 1, 3, 15, 105", wait=2.5)
    s.say("生成函數把摺積變乘積", wait=2.0)
    s.say("幾何級數把無限多張圖收成一個分母：", "$g=\\frac{1}{g_0^{-1}-\\Sigma}$", wait=3.0)


class W704Dyson(ChuScene):
    def construct(self):
        self.show_header("W704　生成函數與 Dyson 方程")
        check_math()
        part_question(self)
        part_gf(self)
        part_feynman(self)
        part_dyson(self)
        part_irreducible(self)
        self.wait(1.0)
