"""
W1205 O(n) 與自迴避漫步（§12.3–12.4，pp. 300–304，(12.61)–(12.72)，Fig. 12.14–12.17，Ex. 12.5–12.8）

問：自旋有 n 個分量時，12 和 36 會變成什麼？
  1. (12.61) 頂點 (S·S)(S·S)（Fig. 12.14：波浪線＝u、兩個點積各一條實線）；Ex. 12.5 (12.63) 的殼層配對 δ_{αβ}/(1+r)。
  2. Fig. 12.15：r 的兩種圖 → 2n + 4（(12.64)）；u 的圖 → 4(n+8)；表 n = 0..3（on_counts 的張量縮併，n = 0 線性外插）。
  3. (12.66)–(12.70)：ν = ½ + (n+2)ε/(4(n+8))；n = 0：9/16 = 0.5625；表 n = 0..3（一階／不展開／已知（引用））；WF 點 n = 0。
  4. Flory (12.71)(12.72)：3/(6 − ε) 的展開 vs 1/16；N = 1000 的三根長條。
  5. ε 展開的四個隱憂（p. 304）；收束。
畫面上的數字都由 utils.walk_eps 算出，check_math() 核對。
模組層級的畫圖函式：wavy、vertex_fig（Fig. 12.14）、r_loop_fig／r_cross_fig（Fig. 12.15）、size_bars。
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
from utils import walk_draw as D

OK_COLOR = "#7CE38B"
WAVE_COLOR = TEAL
LEG_COLOR = WHITE
SHELL_COLOR = ORANGE
NS = (0, 1, 2, 3)
BAR_NU = (Fraction(1, 2), Fraction(9, 16), 0.588)
BAR_COLORS = (GREY_B, YELLOW, OK_COLOR)


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def counts():
    """n = 1..4 的張量縮併（einsum）；n = 0 由 n = 1、2 線性外插。"""
    oc = {n: E.on_counts(n) for n in (1, 2, 3, 4)}
    oc[0] = tuple(2 * a - b for a, b in zip(oc[1], oc[2]))
    return oc


def check_math():
    oc = counts()
    for n in (0, 1, 2, 3, 4):
        assert abs(oc[n][0] - (2 * n + 4)) < 1e-9 and abs(oc[n][1] - 4 * (n + 8)) < 1e-9, n
        assert E.AB(n) == (4 * (n + 2), 4 * (n + 8)) and E.AB(n)[0] == 2 * (2 * n + 4)
    assert [int(round(oc[n][0])) for n in NS] == [4, 6, 8, 10]
    assert [E.AB(n)[0] for n in NS] == [8, 12, 16, 20] and [E.AB(n)[1] for n in NS] == [32, 36, 40, 44]
    assert E.AB(1) == (12, 36) and (round(oc[1][0]), round(oc[1][1])) == (6, 36)
    # §3 ν
    assert E.nu_first(1, 0) == Fraction(9, 16) and f4(Fraction(9, 16)) == "0.5625"
    assert math.floor(Fraction(9, 16) * 1000 + Fraction(1, 2)) == 563              # 書上 0.563 是 0.5625 四捨五入
    assert all(E.nu_first(Fraction(e), 0) == Fraction(1, 2) + Fraction(e, 16) for e in (1, 2, 3))
    assert [f4(E.nu_first(1, n)) for n in NS] == ["0.5625", "0.5833", "0.6000", "0.6136"]
    assert [f4(E.nu_unexpanded(1, n)) for n in NS] == ["0.5714", "0.6000", "0.6250", "0.6471"]
    assert [E.NU_3D_CITED[n] for n in NS] == [0.588, 0.630, 0.672, 0.711]
    assert E.fixed_point(Fraction(1), 0) == (Fraction(-1, 9), Fraction(2, 81))
    # §4 Flory
    assert E.flory_nu(3) == Fraction(3, 5) and f"{float(E.flory_nu(3)):.1f}" == "0.6"
    assert [E.flory_series(k) for k in range(4)] == [Fraction(1, 2), Fraction(1, 12), Fraction(1, 72), Fraction(1, 432)]
    assert E.flory_series(0) + E.flory_series(1) == E.nu_first(1, 1) == Fraction(7, 12) and f4(Fraction(7, 12)) == "0.5833"
    assert abs(sum(float(E.flory_series(k)) * 0.5 ** k for k in range(60)) - 3 / 5.5) < 1e-12
    assert E.nu_unexpanded(Fraction(1, 2)) == Fraction(3, 6) / (1 - Fraction(1, 12)) == Fraction(3, 6 - Fraction(1, 2))
    assert Fraction(1, 12) != Fraction(1, 16)
    assert [f"{E.saw_size(1000, v):.2f}" for v in BAR_NU] == ["31.62", "48.70", "58.08"]


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
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
def wavy(a, b, waves=4, amp=0.08, color=WAVE_COLOR):
    a, b = np.array(a), np.array(b)
    v = b - a
    nrm = np.array([-v[1], v[0], 0]) / np.linalg.norm(v)
    ts = np.linspace(0, 1, 8 * waves + 1)
    pts = [a + t * v + amp * math.sin(2 * PI * waves * t) * nrm for t in ts]
    return VMobject().set_points_smoothly(pts).set_stroke(color, 4)


def vertex_fig(origin=ORIGIN, half=0.9, leg=0.75):
    """Fig. 12.14：波浪線（u）兩端各一個點積，兩條實線＝同一個分量。"""
    o = np.array(origin)
    A, B = o + half * LEFT, o + half * RIGHT
    w = wavy(A, B)
    dots = VGroup(Dot(A, radius=0.08), Dot(B, radius=0.08)).set_color(WHITE)
    dirs = (np.array([-1, 1, 0]), np.array([-1, -1, 0]), np.array([1, 1, 0]), np.array([1, -1, 0]))
    ends = [A + leg * dirs[0] / math.sqrt(2), A + leg * dirs[1] / math.sqrt(2),
            B + leg * dirs[2] / math.sqrt(2), B + leg * dirs[3] / math.sqrt(2)]
    legs = VGroup(*[Line(A if k < 2 else B, e).set_stroke(LEG_COLOR, 4) for k, e in enumerate(ends)])
    labs = VGroup(*[tx(s, 26, GREY_A).next_to(e, LEFT if k < 2 else RIGHT, buff=0.08)
                    for k, (s, e) in enumerate(zip((r"\alpha", r"\alpha", r"\beta", r"\beta"), ends))])
    ulab = tx("u", 26, WAVE_COLOR).next_to(w, UP, buff=0.12)
    return VGroup(w, legs, dots, labs, ulab)


def r_loop_fig(origin=ORIGIN, half=0.7, leg=0.6):
    """Fig. 12.15 左：殼層的兩隻手在同一個點積 → 自己繞一圈（Σ_β δ_ββ = n）。"""
    o = np.array(origin)
    A, B = o + half * LEFT, o + half * RIGHT
    w = wavy(A, B, waves=3)
    loop = Circle(radius=0.38).move_to(A + 0.38 * LEFT).set_stroke(SHELL_COLOR, 4)
    legs = VGroup(Line(B, B + leg * np.array([0.7, 0.7, 0])), Line(B, B + leg * np.array([0.7, -0.7, 0]))
                  ).set_stroke(LEG_COLOR, 4)
    dots = VGroup(Dot(A, radius=0.08), Dot(B, radius=0.08)).set_color(WHITE)
    return VGroup(w, loop, legs, dots)


def r_cross_fig(origin=ORIGIN, half=0.7, leg=0.6):
    """Fig. 12.15 右：殼層的兩隻手在不同點積 → 接起來只給 1。"""
    o = np.array(origin)
    A, B = o + half * LEFT, o + half * RIGHT
    w = wavy(A, B, waves=3)
    arc = ArcBetweenPoints(A, B, angle=-0.75 * PI).set_stroke(SHELL_COLOR, 4)
    legs = VGroup(Line(A, A + leg * np.array([-0.7, -0.7, 0])), Line(B, B + leg * np.array([0.7, -0.7, 0]))
                  ).set_stroke(LEG_COLOR, 4)
    dots = VGroup(Dot(A, radius=0.08), Dot(B, radius=0.08)).set_color(WHITE)
    return VGroup(w, arc, legs, dots)


def size_bars(width=6.0, bar_h=0.42, gap=0.32):
    vals = [E.saw_size(1000, v) for v in BAR_NU]
    k = width / max(vals)
    labs = (r"\sqrt{N}", r"N^{9/16}", r"N^{0.588}")
    out = VGroup()
    for i, (v, lab, c) in enumerate(zip(vals, labs, BAR_COLORS)):
        y = -i * (bar_h + gap)
        bar = Rectangle(width=v * k, height=bar_h).set_fill(c, 0.75).set_stroke(c, 1)
        bar.move_to(np.array([v * k / 2, y, 0]))
        name = tx(lab, 28, c).next_to(bar, LEFT, buff=0.2)
        name.shift((-0.25 - name.get_right()[0]) * RIGHT)
        val = tx(f"{v:.2f}", 26, c).next_to(bar, RIGHT, buff=0.15)
        out.add(VGroup(bar, name, val))
    return out


# ---------------------------------------------------------------- 各段
def part_vertex(s):
    tag = zh("(12.61)　n 分量的四次項", font_size=22, color=GREY_B)
    eq = tx(r"\frac{u}{N}\sum(\vec S\cdot\vec S)(\vec S\cdot\vec S)\,\delta", 36)
    fig = vertex_fig()
    note = zh("Fig. 12.14", font_size=20, color=GREY_B)
    left = VGroup(tag, eq, VGroup(fig, note).arrange(DOWN, buff=0.3)).arrange(DOWN, buff=0.5)
    place_in(left, -6.7, -0.3, -2.85, 3.1)
    tag2 = zh("Ex. 12.5　(12.63)", font_size=22, color=GREY_B)
    eq2 = tx(r"\langle S^\alpha(\vec Q_1)S^\beta(\vec Q_2)\rangle=\frac{\delta_{\alpha\beta}}{1+r}\,\delta(\vec Q_1+\vec Q_2)", 32)
    right = VGroup(tag2, eq2).arrange(DOWN, buff=0.35)
    place_in(right, 0.2, 6.7, -1.0, 1.4)
    assert_in_frame(left, right)
    assert_apart(left, right)
    assert_column(tag, eq, fig)
    s.say("問：自旋有 ", "$n$", " 個分量時，12 和 36 會變成什麼？", wait=2.9)
    s.play(FadeIn(tag), Write(eq), run_time=1.2)
    s.say("(12.61)：四次項是兩個點積相乘", wait=3.0)
    s.play(ShowCreation(fig[0]), FadeIn(fig[2]), FadeIn(fig[4]), FadeIn(note), run_time=1.0)
    s.play(ShowCreation(fig[1]), FadeIn(fig[3]), run_time=1.0)
    s.say("波浪線是 ", "$u$", "；每個點積畫一條實線，分量相同", wait=3.0)
    s.play(FadeIn(tag2), Write(eq2), run_time=1.2)
    s.say("Ex. 12.5：殼層的配對只接同一個分量 ", "$\\alpha=\\beta$", wait=3.0)
    s.clear_stage()


def part_counts(s):
    oc = counts()
    fa = r_loop_fig()
    fb = r_cross_fig()
    la = mixed("同一個點積：", "$\\sum_\\beta\\delta_{\\beta\\beta}=n$", "，2 種", font_size=22, color=SHELL_COLOR)
    lb = mixed("不同點積：", "$1$", "，4 種", font_size=22, color=SHELL_COLOR)
    ga = VGroup(fa, la).arrange(DOWN, buff=0.35)
    gb = VGroup(fb, lb).arrange(DOWN, buff=0.35)
    figs = VGroup(ga, gb).arrange(RIGHT, buff=0.9)
    tot = VGroup(mixed("$r$", "：", "$2n+4$", "（(12.64)，Ex. 12.6）", font_size=26),
                 mixed("$u$", "：", "$4(n+8)$", "（Fig. 12.16–12.17，Ex. 12.7）", font_size=26)).arrange(DOWN, buff=0.3)
    note = zh("Fig. 12.15", font_size=20, color=GREY_B)
    left = VGroup(note, figs, tot).arrange(DOWN, buff=0.4)
    place_in(left, -6.7, 0.3, -2.85, 3.1)
    rows = [["n", "2n+4", "A=4(n+2)", "B=4(n+8)"]] + [[str(n), str(int(round(oc[n][0]))), str(E.AB(n)[0]),
                                                       str(E.AB(n)[1])] for n in NS]
    tab = D.tex_table(rows, [0.7, 1.4, 1.9, 1.9], row_h=0.58, font_size=28)
    for j in range(4):
        tab[0][4 * 2 + j].set_color(YELLOW_A)
    place_in(tab, 0.6, 6.7, -1.6, 2.0)
    assert_in_frame(left, tab)
    assert_apart(left, tab)
    assert_apart(ga, gb)
    assert_column(note, figs, tot)
    assert_table(tab, 4)
    s.play(FadeIn(note), ShowCreation(fa[0]), FadeIn(fa[2:]), run_time=1.0)
    s.play(ShowCreation(fa[1]), FadeIn(la), run_time=1.0)
    s.say("$r$", " 的圖：殼層的兩隻手在同一個點積，繞一圈", wait=3.0)
    s.say("繞一圈要對分量求和：", "$\\sum_\\beta\\delta_{\\beta\\beta}=n$", wait=3.0)
    s.play(ShowCreation(fb[0]), FadeIn(fb[2:]), run_time=0.8)
    s.play(ShowCreation(fb[1]), FadeIn(lb), run_time=1.0)
    s.say("在不同點積：分量已經接死，只給 1，有 4 種", wait=3.0)
    s.play(FadeIn(tot[0]), run_time=0.8)
    s.say("合起來 ", "$2n+4$", "；乘 2 就是 ", "$A=4(n+2)$", wait=2.9)
    s.play(FadeIn(tot[1]), run_time=0.8)
    s.say("$u$", " 的圖數出 ", "$B=4(n+8)$", wait=3.0)
    s.play(FadeIn(tab[1]), FadeIn(tab[0]), run_time=1.0)
    s.say("程式真的做了 ", "$n=1..4$", " 的張量縮併：", "$n=1$", " 回到 12、36", wait=3.0)
    s.clear_stage()


def part_nu(s):
    tag = zh("(12.66)–(12.70)", font_size=22, color=GREY_B)
    eqs = VGroup(tx(r"\frac{dr}{d\ell}=2r+\frac{4(n+2)u}{1+r}", 32),
                 tx(r"\frac{du}{d\ell}=\epsilon u-\frac{4(n+8)u^2}{(1+r)^2}", 32),
                 tx(r"\nu=\frac12+\frac{(n+2)\epsilon}{4(n+8)}", 36).set_color(YELLOW),
                 mixed("$n=0$", "：", "$\\nu=\\frac12+\\frac{\\epsilon}{16}$", font_size=26),
                 tx(r"\epsilon=1:\ \nu=\frac{9}{16}=%s" % f4(Fraction(9, 16)), 32, OK_COLOR)).arrange(DOWN, buff=0.3)
    left = VGroup(tag, eqs).arrange(DOWN, buff=0.3)
    place_in(left, -6.7, -0.6, -2.85, 3.1)
    heads = [Tex("n", font_size=28).set_color(YELLOW), zh("一階", font_size=22, color=YELLOW),
             zh("不展開", font_size=22, color=YELLOW), zh("已知（引用）", font_size=22, color=YELLOW)]
    rows = [heads] + [[str(n), f4(E.nu_first(1, n)), f4(E.nu_unexpanded(1, n)), f"{E.NU_3D_CITED[n]:.3f}"] for n in NS]
    tab = D.tex_table(rows, [0.7, 1.5, 1.5, 1.9], row_h=0.55, font_size=28)
    names = VGroup(*[zh(t, font_size=18, color=GREY_B) for t in ("自迴避", "Ising", "XY", "Heisenberg")])
    for i, m in enumerate(names):
        m.move_to(np.array([5.6 + 0.15 + m.get_width() / 2, -(i + 1.5) * 0.55, 0]))
    cap = mixed("$\\epsilon=1$", "，三維", font_size=22, color=GREY_A).move_to(np.array([2.8, 0.3, 0]))
    wf = mixed("$n=0$", " 的 WF 點：", "$(%s,\\ %s)$" % tuple(ftex(x) for x in E.fixed_point(Fraction(1), 0)),
               font_size=24, color=GREY_A).move_to(np.array([2.8, -3.3, 0]))
    blk = VGroup(cap, tab, names, wf)
    place_in(blk, 0.0, 6.7, -2.85, 3.0)
    assert_in_frame(left, blk)
    assert_apart(left, blk)
    assert_column(tag, *eqs)
    assert_table(tab, 4)
    assert_column(cap, tab, wf)
    for i, m in enumerate(names):
        assert_apart(tab[0][4 * (i + 1) + 3], m)
    s.play(FadeIn(tag), Write(eqs[0]), Write(eqs[1]), run_time=1.4)
    s.say("(12.66)(12.67)：12 換成 ", "$4(n+2)$", "、36 換成 ", "$4(n+8)$", wait=3.0)
    s.play(Write(eqs[2]), run_time=1.0)
    s.say("同一套推導，得到 (12.68)", wait=3.0)
    s.play(FadeIn(eqs[3]), run_time=0.8)
    s.say("$n=0$", " 就是自迴避漫步：", "$\\nu=\\tfrac12+\\tfrac{\\epsilon}{16}$", wait=2.9)
    s.play(FadeIn(eqs[4]), run_time=0.8)
    s.say("(12.70)：", "$\\tfrac{9}{16}=0.5625$", "；書上的 0.563 是四捨五入", wait=3.0)
    s.play(FadeIn(cap), FadeIn(tab[1]), FadeIn(tab[0][:4]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(VGroup(*tab[0][4 * i:4 * i + 4], names[i - 1])) for i in range(1, 5)], lag_ratio=0.3),
           run_time=1.6)
    s.say("$n$", " 越大 ", "$\\nu$", " 越大，和已知值（引用）的趨勢一樣", wait=3.0)
    s.play(FadeIn(wf), run_time=0.8)
    s.say("$n=0$", " 的 Wilson–Fisher 點在 ", "$(-\\tfrac19,\\tfrac{2}{81})$", wait=2.9)
    s.clear_stage()


def part_flory(s):
    tag = zh("(12.71)(12.72)　Flory", font_size=22, color=GREY_B)
    eqs = VGroup(tx(r"\nu_F=\frac{3}{d+2}=\frac35=0.6", 34),
                 tx(r"\frac{3}{6-\epsilon}=\frac12+\frac{\epsilon}{12}+\frac{\epsilon^2}{72}+\frac{\epsilon^3}{432}+\cdots", 34),
                 mixed("一次項 ", "$\\tfrac{1}{12}$", " vs ", "$n=0$", " 的 ", "$\\tfrac{1}{16}$", "：Flory 不是精確的",
                       font_size=26, color=YELLOW)).arrange(DOWN, buff=0.35)
    top = VGroup(tag, eqs).arrange(DOWN, buff=0.3)
    place_in(top, -6.7, 6.7, 0.1, 3.1)
    bars = size_bars()
    btitle = mixed("$N=1000$", " 步的大小 ", "$N^\\nu$", font_size=24, color=GREY_A)
    bot = VGroup(btitle, bars).arrange(DOWN, buff=0.3)
    place_in(bot, -4.5, 4.5, -2.85, -0.2)
    assert_in_frame(top, bot)
    assert_apart(top, bot)
    assert_column(tag, *eqs)
    assert_column(*bars)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.0)
    s.say("Flory：", "$\\nu_F=\\tfrac{3}{d+2}$", "，三維得 ", "$\\tfrac35=0.6$", wait=2.9)
    s.play(Write(eqs[1]), run_time=1.2)
    s.say("Ising 的 ", "$\\tfrac{1}{2-\\epsilon/3}$", " 恰好就是 ", "$\\tfrac{3}{6-\\epsilon}$", wait=3.0)
    s.play(FadeIn(eqs[2]), run_time=0.8)
    s.say("一次項是 ", "$\\tfrac1{12}$", "，不是 ", "$n=0$", " 的 ", "$\\tfrac1{16}$", wait=2.9)
    s.say("但數值接近：", "$0.5833$", " vs ", "$0.5625$", wait=3.0)
    s.play(FadeIn(btitle), *[GrowFromEdge(b[0], LEFT) for b in bars], *[FadeIn(b[1]) for b in bars], run_time=1.2)
    s.play(*[FadeIn(b[2]) for b in bars], run_time=0.6)
    s.say("$N=1000$", "：", "$\\sqrt N=31.62$", "、", "$N^{9/16}=48.70$", wait=3.0)
    s.say("已知（引用）", "$N^{0.588}=58.08$", "：自迴避的高分子撐得更開", wait=3.0)
    s.clear_stage()


def part_caveats(s):
    title = zh("ε 展開的四個隱憂（p. 304）", font_size=30, color=YELLOW)
    items = VGroup(mixed("1. ", "$\\epsilon=1$", " 並不小", font_size=28),
                   zh("2. 級數是漸近的：收斂半徑是 0", font_size=28),
                   mixed("3. 可能有 ", "$e^{-A/\\epsilon}$", " 這種冪級數抓不到的項", font_size=28),
                   zh("4. 重整化群本身沒有被嚴格證明", font_size=28)).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    col = VGroup(title, items).arrange(DOWN, buff=0.55)
    place_in(col, -6.0, 6.0, -2.6, 2.9)
    assert_in_frame(col)
    assert_column(title, *items)
    s.play(FadeIn(title), run_time=0.6)
    notes = [("$\\epsilon=1$", " 是在拿小參數的展開算大參數",),
             ("ε 的級數越算越多項，最後反而發散",),
             ("有些效應根本不是 ", "$\\epsilon$", " 的冪級數"),
             ("整個方法到現在仍缺嚴格的證明",)]
    for m, note in zip(items, notes):
        s.play(FadeIn(m, shift=0.2 * RIGHT), run_time=0.6)
        s.say(*note, wait=2.6)
    s.say("同一套遞迴，", "$n\\to0$", " 就是自迴避漫步", wait=2.9)
    s.say("第 1 章的 ", "$\\sqrt N$", " 在這裡變成 ", "$N^{0.5625}$", "（已知 ", "$N^{0.588}$", "）", wait=3.0)


class W1205OnSaw(ChuScene):
    def construct(self):
        self.show_header("W1205　n 分量與自迴避漫步")
        check_math()
        part_vertex(self)
        part_counts(self)
        part_nu(self)
        part_flory(self)
        part_caveats(self)
        self.wait(1.0)
