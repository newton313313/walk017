"""
W1206 圖解法的重整化群＋全書回顧（§12.5–12.6，pp. 304–317，(12.73)–(12.89)，Fig. 12.18–12.28）

問：不切殼層、直接用費曼圖，也能做重整化嗎？
  1. (12.73) Dyson 級數：質量算子 Σ、雙線＝重整化傳播子（Fig. 12.18–12.20）；動量排序截斷（Fig. 12.21 西瓜圖示意）。
  2. Γ4、u_R = q0^{4−d+2η} g（(12.81)(12.82)）；(12.86) ⇒ dg0/dℓ = εg0 − Bg0² ＝ (12.39)；g0(e^{−5}) = 0.0270 ＝ W1204 的 u(5)。
  3. 冪次數算（Fig. 12.27、(12.89)）：u²r^{d−3} = r(ur^{−(4−d)/2})²；筆誤框 1；數值冪次 −1.0000、−0.5000；反項與 Callan–Symanzik（文字）。
  4. 第 12 章回顧（6 列，無表頭）。
  5. 全書回顧（兩頁各 6 列：E.BOOK_RECAP、E.recap_numbers()、E.scene_range）；收束。
畫面上的數字都由 utils.walk_eps 算出，check_math() 核對。
模組層級的畫圖函式：dyson_row（Fig. 12.18–12.20）、watermelon_fig（Fig. 12.21）、g_axes／g_curve、book_page。
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
G_ELLS = (1, 2, 3, 4, 5)
G0 = 0.005
G_N = 40


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {"g": [E.g_numeric(math.exp(-l), G0, 1.0) for l in G_ELLS],
                "pc": (E.pc_numeric(2), E.pc_numeric(2.5)),
                "recap": E.recap_numbers(),
                "toy": E.toy_coeffs(*E.TOY)}


def check_math():
    nb = numbers()
    # §2 g0(q0)
    assert f4(nb["g"][-1]) == "0.0270" and f4(math.exp(-5)) == "0.0067"
    for l, g in zip(G_ELLS, nb["g"]):
        assert abs(g - E.u_analytic(l, G0, 1.0)) < 1e-9
    assert f4(E.u_analytic(5, 0.005, 1.0)) == "0.0270" and f4(1 / 36) == "0.0278" and E.AB(1)[1] == 36
    # §3 冪次數算
    for d in (2, Fraction(5, 2), 3, Fraction(7, 2)):
        d = Fraction(d)
        assert Fraction(3, 2) * d - d / 2 - 3 == E.pc_exponent(d) == E.pc_rhs(d) == d - 3 == 1 - (4 - d)
    assert E.pc_exponent(2, book=True) == 1 and E.pc_exponent(Fraction(5, 2), book=True) == 2
    assert [f4(v) for v in nb["pc"]] == ["-1.0000", "-0.5000"]
    # §4 第 12 章回顧
    rp, up, _ = nb["toy"]
    assert (f4(rp), f4(up)) == ("0.3000", "0.0075") and E.TOY == (0.2, 0.01)
    assert E.pair_factor(3, 1) == 10 and math.factorial(5) // (math.factorial(3) * math.factorial(2)) == 10
    assert E.ring_aut(3) == 6
    assert E.fixed_point(Fraction(1)) == (Fraction(-1, 7), Fraction(1, 49))
    assert E.nu_first(1) == Fraction(7, 12) and f4(E.nu_first(1)) == "0.5833" and f4(E.nu_first(1, 0)) == "0.5625"
    # §5 全書回顧
    rec = nb["recap"]
    assert rec == E.RECAP_EXPECT
    assert len(E.BOOK_RECAP) == 12 and [c for c, _, _ in E.BOOK_RECAP] == list(range(1, 13))
    for ch, _, tex in E.BOOK_RECAP:
        for v in rec[ch]:
            if ch == 9:
                continue
            if "/" in v:
                a, b = v.split("/")
                forms = (r"\tfrac{%s}{%s}" % (a, b), r"\tfrac%s%s" % (a, b), v)
                assert any(fm in tex for fm in forms), (ch, v)
            else:
                assert v in tex, (ch, v)
    assert f"{(0.600 - 0.588) / 0.588:.4f}" == rec[9][0] == "0.0204" and round(100 * 0.0204) == 2
    assert E.scene_range(1) == "W101–W112" and E.scene_range(12) == "W1201–W1206"
    assert sum(E.BOOK_SCENES.values()) == 87


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def cap_width(m, w):
    if m.get_width() > w:
        m.set_width(w)
    return m


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
def seg(w=0.7):
    return Line(ORIGIN, w * RIGHT).set_stroke(WHITE, 3)


def sigma_blob():
    return VGroup(Circle(radius=0.27).set_fill(GREY_D, 1).set_stroke(ORANGE, 2), tx(r"\Sigma", 24, ORANGE))


def chain(k):
    parts = [seg()]
    for _ in range(k):
        parts += [sigma_blob(), seg()]
    return VGroup(*parts).arrange(RIGHT, buff=0)


def dyson_row():
    """Fig. 12.18–12.20：雙線 ＝ 單線 ＋ 單線 Σ 單線 ＋ …"""
    dbl = VGroup(seg(1.0).shift(0.06 * UP), seg(1.0).shift(0.06 * DOWN)).set_stroke(YELLOW, 3)
    return VGroup(dbl, tx("=", 32), chain(0), tx("+", 32), chain(1), tx("+", 32), chain(2), tx(r"+\cdots", 32)
                  ).arrange(RIGHT, buff=0.22)


def watermelon_fig(half=1.0, ext=0.8):
    """Fig. 12.21 示意：兩個頂點之間三條內線（q1、q2、q3），外線 q。"""
    A, B = half * LEFT, half * RIGHT
    inner = VGroup(ArcBetweenPoints(A, B, angle=-PI / 2), Line(A, B), ArcBetweenPoints(A, B, angle=PI / 2)
                   ).set_stroke(TEAL, 3)
    outer = VGroup(Line(A + ext * LEFT, A), Line(B, B + ext * RIGHT)).set_stroke(WHITE, 3)
    dots = VGroup(Dot(A, radius=0.08), Dot(B, radius=0.08)).set_color(WHITE)
    labs = VGroup(tx("q_1", 22, TEAL).next_to(inner[0], UP, buff=0.08),
                  tx("q_2", 22, TEAL).move_to(0.18 * UP),
                  tx("q_3", 22, TEAL).next_to(inner[2], DOWN, buff=0.08),
                  tx("q", 22).next_to(outer[0], UP, buff=0.1), tx("q", 22).next_to(outer[1], UP, buff=0.1))
    return VGroup(outer, inner, dots, labs)


def g_axes():
    return D.plot_axes((0, 6.4, 1), (0, 0.034, 0.01), width=5.0, height=3.5, x_ticks=(1, 2, 3, 4, 5, 6),
                       y_ticks=(0.01, 0.02, 0.03), y_tick_labels=["0.01", "0.02", "0.03"],
                       x_label=r"\ell=-\ln q_0", y_label="g_0")


def g_curve(ax):
    ps = [ax.c2p(t, E.u_analytic(t, G0, 1.0)) for t in np.linspace(0, 6, G_N)]
    return VMobject().set_points_smoothly(ps).set_stroke(TEAL, 4)


def book_page(chs, fs=24):
    rows = []
    for ch, name, tex in E.BOOK_RECAP:
        if ch not in chs:
            continue
        rows.append([zh(f"{ch}　{name}", font_size=fs, color=YELLOW),
                     zh(E.scene_range(ch), font_size=fs - 4, color=GREY_B),
                     cap_width(tx(tex, 26), 7.6)])
    return D.tex_table(rows, [3.0, 2.3, 7.8], row_h=0.74)


# ---------------------------------------------------------------- 各段
def part_dyson(s):
    tag = zh("(12.73)　質量算子 Σ", font_size=22, color=GREY_B)
    eq = tx(r"\frac{1}{r+q^2}+\frac{1}{r+q^2}\,\Sigma\,\frac{1}{r+q^2}+\cdots=\frac{1}{r+q^2-\Sigma}", 36)
    row = dyson_row()
    note = zh("Fig. 12.18–12.20：雙線＝重整化傳播子", font_size=20, color=GREY_B)
    top = VGroup(tag, eq, VGroup(row, note).arrange(DOWN, buff=0.25)).arrange(DOWN, buff=0.4)
    place_in(top, -6.7, 6.7, 0.0, 3.1)
    wm = watermelon_fig()
    order = tx(r"q<q_0<q_1,\,q_2,\,q_3", 32, YELLOW)
    wnote = mixed("Fig. 12.21 示意：動量依大小排序、截斷在 ", "$q_0$", font_size=20, color=GREY_B)
    bot = VGroup(wm, VGroup(order, wnote).arrange(DOWN, buff=0.25)).arrange(RIGHT, buff=0.8)
    place_in(bot, -6.7, 6.7, -2.85, -0.3)
    assert_in_frame(top, bot)
    assert_apart(top, bot)
    assert_column(tag, eq, row, note)
    s.say("問：不切殼層、直接用費曼圖，也能做重整化嗎？", wait=2.2)
    s.play(FadeIn(tag), Write(eq), run_time=1.4)
    s.say("(12.73)：把所有插入 ", "$\\Sigma$", " 的圖加起來", wait=2.2)
    s.play(FadeIn(row), FadeIn(note), run_time=1.0)
    s.say("雙線＝重整化的傳播子 ", "$1/(r+q^2-\\Sigma)$", wait=2.2)
    s.say("和第 7 章 W704 的 Dyson 方程同一個形狀", wait=2.2)
    s.play(FadeIn(wm), run_time=1.0)
    s.play(Write(order), FadeIn(wnote), run_time=1.0)
    s.say("內線動量依大小排序，截斷在 ", "$q_0$", "：", "$q<q_0<q_1,q_2,q_3$", wait=2.4)
    s.clear_stage()


def part_g(s):
    nb = numbers()
    tag = zh("(12.81)–(12.86)", font_size=22, color=GREY_B)
    eqs = VGroup(mixed("$\\Gamma_4$", " 定出 ", "$u_R=q_0^{4-d+2\\eta}\\,g$", font_size=26),
                 tx(r"q_0\frac{dg_0}{dq_0}=(-4+d-2\eta)g_0+\mathcal{B}g_0^2", 30),
                 tx(r"q_0=e^{-\ell}:\ \frac{dg_0}{d\ell}=\epsilon g_0-\mathcal{B}g_0^2", 30),
                 mixed("＝ (12.39)（", "$\\eta=0$", "，", "$\\mathcal{B}=36$", "）", font_size=24, color=YELLOW)
                 ).arrange(DOWN, buff=0.4)
    left = VGroup(tag, eqs).arrange(DOWN, buff=0.35)
    place_in(left, -6.7, -0.4, -2.85, 3.1)
    ax = g_axes()
    ax.move_to(np.array([3.4, 0.8, 0]))
    cur = g_curve(ax)
    dots = VGroup(*[Dot(ax.c2p(l, g), radius=0.07).set_color(WHITE).set_stroke(TEAL, 2)
                    for l, g in zip(G_ELLS, nb["g"])])
    star = DashedLine(ax.c2p(0, 1 / 36), ax.c2p(6.2, 1 / 36), dash_length=0.08).set_stroke(RED, 2)
    slab = tx(r"g^*=\tfrac{1}{36}=%s" % f4(1 / 36), 22, RED).next_to(ax.c2p(6.2, 1 / 36), UP, buff=0.12).shift(0.9 * LEFT)
    res = VGroup(mixed("$g_0(1)=0.005$", " 積到 ", "$q_0=e^{-5}=%s$" % f4(math.exp(-5)), font_size=22),
                 mixed("$g_0=%s$" % f4(nb["g"][-1]), "＝W1204 的 ", "$u(\\ell=5)$", font_size=22, color=OK_COLOR)
                 ).arrange(DOWN, buff=0.18)
    place_in(res, 0.4, 6.7, -2.85, -1.4)
    assert_in_frame(left, ax, slab, res)
    assert_apart(left, VGroup(ax, res))
    assert_apart(ax, res)
    assert_column(tag, *eqs)
    s.play(FadeIn(tag), FadeIn(eqs[0]), run_time=1.0)
    s.say("四點函數 ", "$\\Gamma_4$", " 定出無因次的耦合 ", "$g$", wait=2.2)
    s.play(Write(eqs[1]), run_time=1.2)
    s.say("(12.86)：", "$g_0$", " 隨截斷 ", "$q_0$", " 怎麼變", wait=2.2)
    s.play(Write(eqs[2]), FadeIn(eqs[3]), run_time=1.2)
    s.say("換變數 ", "$q_0=e^{-\\ell}$", "：正是 (12.39)", wait=2.2)
    s.play(ShowCreation(ax), ShowCreation(star), FadeIn(slab), run_time=1.0)
    s.play(ShowCreation(cur), FadeIn(dots), run_time=1.4)
    s.play(FadeIn(res), run_time=0.8)
    s.say("$g_0(e^{-5})=%s$" % f4(nb["g"][-1]), "：就是 W1204 的 ", "$u(\\ell=5)$", wait=2.4)
    s.say("一樣流向不動點 ", "$g^*=\\tfrac1{36}=%s$" % f4(1 / 36), wait=2.2)
    s.clear_stage()


def part_power(s):
    nb = numbers()
    tag = zh("Fig. 12.27、(12.89)　冪次數算", font_size=22, color=GREY_B)
    lines = VGroup(mixed("第 ", "$n$", " 階含 ", "$(ur^{-(4-d)/2})^n$", "：", "$d<4,\\ r\\to0$", " 時發散", font_size=26),
                   mixed("$u_R\\propto r^{(4-d)/2}$", " 剛好抵消", font_size=26),
                   tx(r"Q\propto r^{1/2}:\ (d^dQ)^3\to r^{3d/2},\quad \delta\to r^{-d/2},\quad G^3\to r^{-3}", 30),
                   tx(r"u^2r^{d-3}=r\,(ur^{-(4-d)/2})^2", 34, YELLOW)).arrange(DOWN, buff=0.32)
    top = VGroup(tag, lines).arrange(DOWN, buff=0.3)
    place_in(top, -6.7, 6.7, 0.1, 3.1)
    box = D.erratum_box("p. 315，(12.89) 後", r"u^2r^{2d-3}=r^{-1}(ur^{-(4-d)/2})^2",
                        r"u^2r^{d-3}=r\,(ur^{-(4-d)/2})^2", width=10.0)
    place_in(box, -5.5, 5.5, -1.95, -0.1)
    cs = mixed("反項（Fig. 12.28）與 Callan–Symanzik 方程：截斷 ", "$q_0$", " 是算法的產物，答案不能依賴它",
               font_size=22, color=GREY_A)
    place_in(cs, -6.7, 6.7, -2.85, -2.2)
    assert_in_frame(top, box, cs)
    assert_column(tag, *lines)
    assert_column(top, box, cs)
    s.play(FadeIn(tag), FadeIn(lines[0]), run_time=1.0)
    s.say("微擾的第 ", "$n$", " 階在 ", "$r\\to0$", " 時發散", wait=2.2)
    s.play(FadeIn(lines[1]), run_time=0.8)
    s.say("換成 ", "$u_R$", " 之後，發散剛好抵消", wait=2.0)
    s.play(Write(lines[2]), run_time=1.2)
    s.say("西瓜圖：三個 ", "$d^dQ$", "、一個 ", "$\\delta$", "、三個 ", "$G$", wait=2.2)
    s.play(Write(lines[3]), run_time=1.0)
    s.say("合起來 ", "$r^{3d/2-d/2-3}=r^{d-3}$", wait=2.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("數值積分量出的冪次：", "$d=2$", "：", "$%s$" % f4(nb["pc"][0]), "、", "$d=2.5$", "：", "$%s$" % f4(nb["pc"][1]),
          wait=2.4)
    s.say("照書上的 ", "$r^{2d-3}$", " 是 ", "$+1$", "、", "$+2$", "：對不上", wait=2.2)
    s.play(FadeIn(cs), run_time=0.8)
    s.say("截斷 ", "$q_0$", " 只是工具：反項把它的影響消掉", wait=2.4)
    s.clear_stage()


def part_ch12(s):
    nb = numbers()
    rp, up, _ = nb["toy"]
    title = zh("第 12 章回顧", font_size=30, color=YELLOW)
    fs = 22
    labels = ["W1201 加上四次項", "W1202 配對遊戲", "W1203 一迴圈的流", "W1204 ε 展開求 ν", "W1205 n 分量與自迴避",
              "W1206 圖解法"]
    items = [mixed("$r'=r+\\frac{12u}{1+r}$", "、", "$u'=u-\\frac{36u^2}{(1+r)^2}$", "（玩具 ", "$%s$" % f4(rp), "、",
                   "$%s$" % f4(up), "）", font_size=fs),
             mixed("只剩連通圖；", "$\\frac{5!}{3!\\,2!}=%d$" % E.pair_factor(3, 1), "、", "$f_3=%d$" % E.ring_aut(3),
                   font_size=fs),
             mixed("Wilson–Fisher 點 ", "$(-\\frac17,\\ \\frac1{49})$", font_size=fs),
             mixed("$\\nu=\\frac12+\\frac{\\epsilon}{12}=%s$" % f4(E.nu_first(1)), font_size=fs),
             mixed("$\\nu=\\frac12+\\frac{(n+2)\\epsilon}{4(n+8)}$", "；", "$n=0$", "：", "$%s$" % f4(E.nu_first(1, 0)),
                   font_size=fs),
             mixed("(12.86)＝(12.39)、", "$u^2r^{d-3}$", font_size=fs)]
    rows = [[zh(t, font_size=fs, color=YELLOW), m] for t, m in zip(labels, items)]
    tab = D.tex_table(rows, [3.6, 9.0], row_h=0.8)
    col = VGroup(title, tab[0]).arrange(DOWN, buff=0.35)
    place_in(col, -6.7, 6.7, -2.85, 3.1)
    assert_in_frame(col)
    assert_table(tab, 2)
    s.play(FadeIn(title), run_time=0.6)
    s.say("第 12 章回顧：六幕的重點", wait=1.4)
    notes = [("積掉殼層、重新標度：只剩 ", "$r$", " 與 ", "$u$", " 在變"),
             ("配對遊戲：指數化之後只要畫連通圖",),
             ("一迴圈的流：臨界系統都流向 WF 點",),
             ("ε 展開：", "$\\nu$", " 第一次不是 ", "$\\tfrac12$"),
             ("$n\\to0$", "：同一套遞迴就是自迴避漫步"),
             ("圖解法：不切殼層，得到同一個方程",)]
    for i, note in enumerate(notes):
        s.play(FadeIn(tab[0][2 * i]), FadeIn(tab[0][2 * i + 1]), run_time=0.5)
        s.say(*note, wait=1.9)
    s.clear_stage()


def part_book(s):
    for page, chs, notes in (
            (1, range(1, 7), [("全書回顧（一）：第 1–6 章，漫步本身",),
                              ("一維 ", "$\\langle x^2\\rangle=N$", "；三維 sc 回家機率 0.3405"),
                              ("牆用鏡像；持續漫步 ", "$Q=\\tfrac12,\\ N=1000$", "：2996")]),
            (2, range(7, 13), [("全書回顧（二）：自迴避、臨界現象、重整化群",),
                               ("Flory 0.600 vs 0.588：只差 2%（0.0204）",),
                               ("最後用 ", "$\\epsilon$", " 展開算出 ", "$\\nu$", "：0.5833、0.5625")])):
        title = zh(f"全書回顧（{'一' if page == 1 else '二'}）", font_size=30, color=YELLOW)
        tab = book_page(set(chs))
        col = VGroup(title, tab[0]).arrange(DOWN, buff=0.35)
        place_in(col, -6.7, 6.7, -2.85, 3.1)
        assert_in_frame(col)
        assert_table(tab, 3)
        s.play(FadeIn(title), run_time=0.6)
        s.say(*notes[0], wait=0.6)
        s.play(LaggedStart(*[FadeIn(tab[0][3 * i:3 * i + 3]) for i in range(6)], lag_ratio=0.4), run_time=3.0)
        s.wait(2.0)
        for note in notes[1:]:
            s.say(*note, wait=3.0)
        if page == 2:
            s.say("全書 12 章、", "$%d$" % sum(E.BOOK_SCENES.values()), " 幕", wait=2.0)
            s.say("從一枚硬幣到一個 ", "$\\epsilon$", "——", wait=2.0)
            s.say("漫步、磁鐵與重整化群說的是同一個故事", wait=3.0)
        else:
            s.clear_stage()


class W1206DiagrammaticRecap(ChuScene):
    def construct(self):
        self.show_header("W1206　圖解法的重整化群＋全書回顧")
        check_math()
        part_dyson(self)
        part_g(self)
        part_power(self)
        part_ch12(self)
        part_book(self)
        self.wait(1.0)
