"""
W802 標度形式：所有曲線疊成一條（§8.1.2 前半，pp. 195–197，(8.8)–(8.14)）

問：z 越靠近 z_c，g(z;q) 越尖越高；不同 z 的曲線，能不能換個座標疊成一條？
  1. (8.8)：g = 1/(1 − zχ)、χ = a − bq² → (1/a)/(z_c − z + b q²/a²)；一維 g(z;0) = 5、10、50；筆誤框 z_c = 2d。
  2. 塌縮：g 對 q → (z_c − z)g 對 x = q/√(z_c − z)；表 x = 0, 1, 2, 4 → 極限 2/(4 + x²)；(8.9)(8.10)。
  3. 實空間 (8.11)–(8.13)：三維 ξ = 2.0412（格子 2.0205）、4.0825（4.0722）；G 比值 1.0231、1.0178；F(0) = 0.4775。
  4. F 的兩個假設（小示意曲線）。
畫面上的數字都由 utils.walk_scale 算出，check_math() 核對。模組層級的畫圖函式：g_curves、collapse_curves、small_axes。
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

EPS = WS.COLLAPSE_EPS                 # z_c − z = 0.1, 0.05, 0.01
XS = WS.COLLAPSE_X                    # 0, 1, 2, 4
ZS = tuple(0.5 - e for e in EPS)      # 0.4, 0.45, 0.49
COLORS = (BLUE_B, TEAL_B, YELLOW)
LIMIT_COLOR = RED_B
NPTS = 59


def f4(x):
    return f"{float(x):.4f}"


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = {"g0": [WS.g_lattice(z, 0.0) for z in ZS],
              "col": {e: [WS.collapse(e, x) for x in XS] for e in EPS},
              "lim": [WS.collapse_limit(x) for x in XS],
              "xi": {x: (WS.xi_88c(x), WS.xi_lattice(x)) for x in WS.XI_CASES},
              "Gr": {x: WS.G_ratio(x, 8) for x in WS.XI_CASES},
              "F0": WS.F3(0)}
    return nb


def check_math():
    nb = numbers()
    assert ZS == (0.4, 0.45, 0.49) and EPS == (0.1, 0.05, 0.01) and XS == (0, 1, 2, 4)
    assert WS.zc(1) == Fraction(1, 2) and WS.AB(1) == (2, 1) and WS.zc(3) == Fraction(1, 6) and WS.AB(3) == (6, 1)
    assert [f"{v:.0f}" for v in nb["g0"]] == ["5", "10", "50"]
    assert all(abs(v - 1 / (1 - 2 * z)) < 1e-9 for v, z in zip(nb["g0"], ZS))
    assert [f4(v) for v in nb["col"][0.1]] == ["0.5000", "0.4172", "0.2819", "0.1317"]
    assert [f4(v) for v in nb["col"][0.05]] == ["0.5000", "0.4085", "0.2652", "0.1145"]
    assert [f4(v) for v in nb["col"][0.01]] == ["0.5000", "0.4017", "0.2529", "0.1027"]
    assert [f4(v) for v in nb["lim"]] == ["0.5000", "0.4000", "0.2500", "0.1000"]
    assert all(abs(WS.collapse_limit(x) - 2 / (4 + x * x)) < 1e-12 for x in XS)
    assert [(f4(a), f4(b)) for a, b in nb["xi"].values()] == [("2.0412", "2.0205"), ("4.0825", "4.0722")]
    assert abs(nb["xi"][0.99][0] / nb["xi"][0.96][0] - 2) < 1e-9 and abs((1 - 0.99) / (1 - 0.96) - 0.25) < 1e-9
    assert [f4(v) for v in nb["Gr"].values()] == ["1.0231", "1.0178"]
    assert f4(nb["F0"]) == "0.4775" and abs(nb["F0"] - 3 / (2 * math.pi)) < 1e-12
    for x in WS.XI_CASES:                                    # rG = F(u)：r/ξ = 6u
        r = 8
        u = r * math.sqrt((1 - x) / 6)
        assert abs(r * WS.G_yukawa(x, r) - WS.F3(u)) < 1e-12


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


def axes_with_labels(x_range, y_range, width, height, x_ticks, y_ticks, x_lab, y_lab, y_fmt):
    """plot_axes（縱軸畫在 x = 0）：縱軸刻度數字貼在 x = 0 的縱軸左邊。"""
    ax = D.plot_axes(x_range, y_range, width=width, height=height, x_ticks=x_ticks, font_size=22)
    labs = VGroup()
    for v in y_ticks:
        labs.add(tx(y_fmt(v), 20, GREY_B).next_to(ax.c2p(0, v), LEFT, buff=0.1))
    labs.add(tx(x_lab, 26, GREY_B).next_to(ax.x_axis.get_end(), DOWN, buff=0.2))
    labs.add(tx(y_lab, 26, GREY_B).next_to(ax.y_axis.get_end(), RIGHT, buff=0.15))
    return ax, labs


def smooth(ax, pts, color, dashed=False):
    c = VMobject().set_points_smoothly([ax.c2p(x, y) for x, y in pts]).set_stroke(color, 3)
    if dashed:
        c = DashedVMobject(c, num_dashes=40).set_stroke(color, 3)
    return c


def g_curves(ax):
    qs = np.linspace(-1, 1, NPTS)
    return VGroup(*[smooth(ax, [(q, WS.g_lattice(z, q)) for q in qs], c) for z, c in zip(ZS, COLORS)])


def collapse_curves(ax):
    xs = np.linspace(-4, 4, NPTS)
    cur = VGroup(*[smooth(ax, [(x, WS.collapse(e, x)) for x in xs], c) for e, c in zip(EPS, COLORS)])
    lim = smooth(ax, [(x, WS.collapse_limit(x)) for x in xs], LIMIT_COLOR, dashed=True)
    return cur, lim


def small_axes(width=4.0, height=2.0):
    ax = D.plot_axes((0, 1.0, 0.5), (0, 0.55, 0.25), width=width, height=height)
    lab = VGroup(tx("u", 24, GREY_B).next_to(ax.x_axis.get_end(), DOWN, buff=0.15),
                 tx("F(u)", 24, GREY_B).next_to(ax.y_axis.get_end(), RIGHT, buff=0.15))
    graph = smooth(ax, [(u, WS.F3(u)) for u in np.linspace(0, 1.0, 40)], YELLOW)
    return VGroup(ax, lab, graph)


# ---------------------------------------------------------------- 各段
def part_88(s):
    nb = numbers()
    e1 = tx(r"g(z;\vec q)=\frac{1}{1-z\chi(\vec q)},\qquad \chi(\vec q)\approx a-b|\vec q|^2", 36)
    e2 = tx(r"g\approx\frac{1}{a}\,\frac{1}{z_c-z+\frac{b}{a^2}|\vec q|^2},\qquad z_c=\frac1a=\frac{1}{2d}", 36)
    note = mixed("這裡的 ", "$a$", " 是 ", "$\\chi(\\vec 0)$", "，不是格距（簡單立方 ", "$a=2d$", "、", "$b=1$", "）",
                 font_size=26, color=GREY_B)
    e3 = mixed("一維：", "$g=\\frac{1}{1-2z\\cos q}$", "，", "$z=0.4,\\ 0.45,\\ 0.49$", " → ",
               "$g(z;0)=%s$" % ",\\ ".join(f"{v:.0f}" for v in nb["g0"]), font_size=30)
    col = VGroup(e1, e2, note, e3).arrange(DOWN, buff=0.35)
    place_in(col, -6.7, 6.7, -0.4, 3.1)
    box = D.erratum_box("p. 196，第一段", r"z_c=\text{coordination number}=2d", r"z_c=1/(2d)", width=8.0, font_size=36)
    place_in(box, -6.7, 6.7, -2.9, -0.7)
    assert_in_frame(col, box)
    assert_column(*col, box)
    s.say("問：", "$z$", " 越靠近 ", "$z_c$", "，", "$g(z;q)$", " 越尖越高……", wait=0.5)
    s.say("不同 ", "$z$", " 的曲線，能不能換個座標疊成一條？", wait=1.5)
    s.play(Write(e1), run_time=1.5)
    s.say("(8.8)：小 ", "$q$", " 時 ", "$\\chi\\approx a-bq^2$", wait=3.0)
    s.play(FadeIn(e2), run_time=1.0)
    s.say("分母在 ", "$z=z_c$", "、", "$q=0$", " 變成 0", wait=3.0)
    s.play(FadeIn(note), run_time=0.8)
    s.say("注意：這裡的 ", "$a$", " 是 ", "$\\chi(\\vec 0)$", "，不是格距", wait=3.5)
    s.play(FadeIn(e3), run_time=1.0)
    s.say("一維 ", "$z_c=\\frac12$", "：", "$z=0.4,\\ 0.45,\\ 0.49$", " 時 ", "$g(z;0)=5,\\ 10,\\ 50$", wait=3.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上說 ", "$z_c$", " 等於配位數 ", "$2d$", "；應該是 ", "$1/(2d)$", wait=3.0)
    s.say("照書上一維 ", "$z_c=2$", "，可是 ", "$z>\\frac12$", " 時 ", "$\\sum(2z)^N$", " 早就發散了", wait=3.5)
    s.clear_stage()


def part_collapse(s):
    nb = numbers()
    ax1, lab1 = axes_with_labels((-1, 1, 0.5), (0, 55, 10), 6.0, 4.6, (-1, 0, 1), (10, 20, 30, 40, 50),
                                 "q", "g(z;q)", lambda v: str(v))
    ax2, lab2 = axes_with_labels((-4, 4, 1), (0, 0.6, 0.1), 6.0, 4.6, (-4, -2, 0, 2, 4), (0.1, 0.2, 0.3, 0.4, 0.5),
                                 r"x=q/\sqrt{z_c-z}", r"(z_c-z)\,g", lambda v: f"{v:.1f}")
    for ax, lab in ((ax1, lab1), (ax2, lab2)):
        place_in(VGroup(ax, lab), -6.7, -0.3, -2.9, 3.0)
    dv = ax1.get_origin() - ax2.get_origin()
    ax2.shift(dv)
    lab2.shift(dv)
    gc = g_curves(ax1)
    cc, lim = collapse_curves(ax2)
    legend = VGroup(*[mixed("$z=%s$" % z, font_size=24, color=c) for z, c in zip(ZS, COLORS)]).arrange(RIGHT, buff=0.4)
    rows = [["x"] + [str(x) for x in XS]]
    for e, c in zip(EPS, COLORS):
        rows.append([r"z_c-z=%s" % e] + [f4(v) for v in nb["col"][e]])
    rows.append([r"\frac{2}{4+x^2}"] + [f4(v) for v in nb["lim"]])
    tab = D.tex_table(rows, [1.9, 1.05, 1.05, 1.05, 1.05], row_h=0.62, font_size=26)
    for i, c in enumerate(COLORS + (LIMIT_COLOR,)):
        D.table_cell(tab, i + 1, 0, 5).set_color(c)
    e9 = tx(r"g=|\vec q|^{-2}f\Big(q(z_c-z)^{-1/2}\Big)", 30)
    e10 = tx(r"g=|\vec q|^{-2+\eta}f\Big(q(z_c-z)^{-\nu}\Big)", 30)
    e10n = mixed("普通漫步 ", "$\\nu=\\frac12$", "、", "$\\eta=0$", "；自迴避不同", font_size=24, color=GREY_B)
    right = VGroup(legend, tab, e9, e10, e10n).arrange(DOWN, buff=0.3)
    place_in(right, 0.0, 6.7, -2.9, 3.1)
    assert_in_frame(VGroup(ax1, lab1), VGroup(ax2, lab2), right)
    assert_apart(VGroup(ax1, lab1), right)
    assert_apart(VGroup(ax2, lab2), right)
    assert_column(*right)
    assert_table(tab, 5)
    s.say("一維三條曲線：", "$z=0.4,\\ 0.45,\\ 0.49$", wait=0.3)
    s.play(ShowCreation(ax1), FadeIn(lab1), FadeIn(legend), run_time=1.0)
    s.play(LaggedStart(*[ShowCreation(c) for c in gc], lag_ratio=0.4), run_time=2.0)
    s.say("越靠近 ", "$z_c$", "：峰越高（5、10、50），也越窄", wait=3.0)
    s.say("縱軸乘上 ", "$z_c-z$", "，橫軸除以 ", "$\\sqrt{z_c-z}$", wait=1.0)
    s.play(ReplacementTransform(ax1, ax2), ReplacementTransform(lab1, lab2),
           *[ReplacementTransform(a, b) for a, b in zip(gc, cc)], run_time=2.5)
    s.play(ShowCreation(lim), run_time=1.0)
    s.say("三條曲線幾乎疊在一起，紅色虛線是極限 ", "$\\frac{2}{4+x^2}$", wait=3.0)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][5:]], lag_ratio=0.04), run_time=1.5)
    s.say("$x=2$", "：0.2819、0.2652、0.2529 → 0.2500", wait=3.5)
    s.say("越靠近 ", "$z_c$", "，越貼近同一條曲線", wait=3.0)
    s.play(FadeIn(e9), run_time=0.8)
    s.say("(8.9)：", "$g$", " 是 ", "$|\\vec q|^{-2}$", " 乘上一個只看 ", "$q(z_c-z)^{-1/2}$", " 的函數", wait=3.5)
    s.play(FadeIn(e10), FadeIn(e10n), run_time=0.8)
    s.say("(8.10) 推廣：兩個指數 ", "$\\nu$", "、", "$\\eta$", "；普通漫步 ", "$\\frac12$", "、0", wait=3.5)
    s.clear_stage()


def part_real(s):
    nb = numbers()
    e1 = mixed("令 ", "$\\vec q=\\vec\\kappa/|\\vec r|$", "：", "$G(z;\\vec r)=|\\vec r|^{-d+2-\\eta}\\,F\\Big(|\\vec r|(z_c-z)^{\\nu}\\Big)$",
               font_size=32)
    e2 = mixed("三維：", "$G\\approx\\frac{3}{2\\pi}\\frac{e^{-r/\\xi}}{r}$", "，",
               "$\\xi=\\frac{1}{\\sqrt{6(1-z/z_c)}}\\propto(z_c-z)^{-1/2}$", font_size=30)
    rows = [[r"z/z_c", r"\xi", mixed("格子精確 ", "$\\xi$", font_size=24, color=YELLOW), mixed("$G$", " 比值（", "$r=8$", "）", font_size=24, color=YELLOW)]]
    for x in WS.XI_CASES:
        a, b = nb["xi"][x]
        rows.append([str(x), f4(a), f4(b), f4(nb["Gr"][x])])
    tab = D.tex_table(rows, [1.6, 1.7, 2.2, 2.6], row_h=0.66, font_size=30)
    D.table_cell(tab, 1, 1, 4).set_color(YELLOW)
    D.table_cell(tab, 2, 1, 4).set_color(YELLOW)
    e3 = mixed("$rG=F(u)=\\frac{3}{2\\pi}e^{-6u}$", "，", "$u=r\\sqrt{z_c-z}$", "，", "$F(0)=%s$" % f4(nb["F0"]),
               font_size=30)
    col = VGroup(e1, e2, tab, e3).arrange(DOWN, buff=0.4)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(*col)
    assert_table(tab, 4)
    s.say("回到實空間 (8.11)–(8.13)：令 ", "$\\vec q=\\vec\\kappa/|\\vec r|$", wait=0.3)
    s.play(FadeIn(e1), run_time=1.0)
    s.say("$G$", " 是一個冪次乘上只看 ", "$|\\vec r|(z_c-z)^\\nu$", " 的函數 ", "$F$", wait=3.5)
    s.play(FadeIn(e2), run_time=1.0)
    s.say("三維簡單立方：Yukawa 形式，衰減長度 ", "$\\xi$", wait=3.0)
    s.play(FadeIn(tab[0][:4]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][4:]], lag_ratio=0.1), run_time=1.0)
    s.say("$\\xi=2.0412$", "、", "$4.0825$", "；格子精確 2.0205、4.0722", wait=3.5)
    s.say("$z_c-z$", " 變成 ", "$\\frac14$", "，", "$\\xi$", " 變成 2 倍：", "$\\xi\\propto(z_c-z)^{-1/2}$", wait=3.5)
    s.say("$r=8$", " 時格子精確 ", "$G$", " 除以 Yukawa：1.0231、1.0178", wait=3.0)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("$rG$", " 只看 ", "$u=r\\sqrt{z_c-z}$", "：", "$F(0)=\\frac{3}{2\\pi}=0.4775$", wait=3.5)
    s.clear_stage()


def part_F(s):
    a1 = mixed("$x\\to0$", "：", "$F(x)$", " 趨於常數", font_size=30)
    a2 = mixed("$x\\to\\infty$", "：", "$F(x)$", " 指數衰減", font_size=30)
    head = zh("書上對 F 的兩個假設", font_size=30, color=YELLOW)
    txt = VGroup(head, a1, a2).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    pic = small_axes()
    cap = mixed("三維：", "$\\frac{3}{2\\pi}e^{-6u}$", " 兩條都滿足", font_size=26, color=GREY_B)
    right = VGroup(pic, cap).arrange(DOWN, buff=0.3)
    place_in(txt, -6.7, -0.3, -1.5, 2.0)
    place_in(right, 0.3, 6.7, -2.5, 2.6)
    assert_in_frame(txt, right)
    assert_apart(txt, right)
    s.say("(8.14) 之前，書上對 ", "$F$", " 做了兩個假設", wait=0.3)
    s.play(FadeIn(head), FadeIn(a1), run_time=1.0)
    s.say("$x\\to0$", "：", "$F$", " 趨於常數——回到普通的冪次", wait=3.0)
    s.play(FadeIn(a2), run_time=0.8)
    s.say("$x\\to\\infty$", "：", "$F$", " 指數衰減——超過 ", "$\\xi$", " 就看不到", wait=3.0)
    s.play(ShowCreation(pic[0]), FadeIn(pic[1]), run_time=0.8)
    s.play(ShowCreation(pic[2]), FadeIn(cap), run_time=1.2)
    s.say("三維的 ", "$\\frac{3}{2\\pi}e^{-6u}$", "：兩條都滿足", wait=3.0)
    s.say("靠近 ", "$z_c$", " 只剩一個長度 ", "$\\xi\\propto(z_c-z)^{-\\nu}$", wait=3.0)
    s.say("用 ", "$\\xi$", " 當尺，所有曲線疊成一條——這就是標度形式", wait=3.5)


class W802ScalingForm(ChuScene):
    def construct(self):
        self.show_header("W802　標度形式：所有曲線疊成一條")
        check_math()
        part_88(self)
        part_collapse(self)
        part_real(self)
        part_F(self)
        self.wait(1.0)
