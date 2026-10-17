"""
W201 生成函數：一個函數裝下整個數列（§2.1.1，pp. 25–27，(2.1)，Ex. 2.1–2.3）

問：無窮多個數，能不能裝進一個式子裡？
  1. e^x 的係數 1, 1, 1/2, 1/6, 1/24：e^x 生成數列 1/n!；第 1 章見過 (1+e^y)^N/2^N（W105）與 1/(1−2z)（W111）。
  2. (2.1) 1/sqrt(z²−2zx+1) = Σ P_n(x) z^n：令 t = 2xz − z²，(1−t)^{−1/2} = 1 + t/2 + 3t²/8 + 5t³/16 + …，
     收集 z 的次方得 P_0..P_3（練習 2.1）；畫圖、先猜 P_n(1)（全是 1）；筆誤框（p. 26：應是把左邊對 z 展開）與 P_5。
  3. 練習 2.2：兩邊從 −1 積到 1：左邊 = 2 → ∫P_0 dx = 2、其餘 0。
  4. 練習 2.3：兩個生成函數相乘再積分 = 2 + (2/3)yz + (2/5)(yz)² + … → ∫P_nP_m dx = 2δ_nm/(2n+1)（正交）。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import boxed, mixed, zh
from utils import walk_gf as GF
from utils import walk_draw as D

P_COLORS = [GREY_A, BLUE_C, "#7CE38B", ORANGE]          # P_0..P_3 的曲線與式子
COEF_COLOR = YELLOW
P = GF.legendre_gf(5)
Z22 = (0.3, 0.6, 0.9)                                    # 練習 2.2 的數值核對


def frac_tex(c, small=False) -> str:
    c = Fraction(c)
    if c.denominator == 1:
        return str(c.numerator)
    sign = "-" if c < 0 else ""
    cmd = r"\tfrac" if small else r"\frac"
    return rf"{sign}{cmd}{{{abs(c.numerator)}}}{{{c.denominator}}}"


def poly_tex(p) -> str:
    """係數 list（x^0 起）→ 例如 \\frac{3x^2-1}{2}。"""
    d = 1
    for c in p:
        den = Fraction(c).denominator
        d = d * den // math.gcd(d, den)
    s = ""
    for k in range(len(p) - 1, -1, -1):
        c = int(Fraction(p[k]) * d)
        if c == 0:
            continue
        mono = "" if k == 0 else ("x" if k == 1 else f"x^{k}")
        body = (str(abs(c)) if abs(c) != 1 or k == 0 else "") + mono
        s += ("-" if c < 0 else ("+" if s else "")) + body
    return s if d == 1 else rf"\frac{{{s}}}{{{d}}}"


def geometric(r, n):
    """1/(1 − rz) 的前 n 個係數（遞推 c_k = r·c_{k−1}）。"""
    c = [1]
    for _ in range(n - 1):
        c.append(c[-1] * r)
    return c


def gram(nmax):
    return [[GF.poly_int11(GF.poly_mul(P[n], P[m])) for m in range(nmax + 1)] for n in range(nmax + 1)]


def check_math():
    F = Fraction
    assert GF.exp_coeffs(4) == [1, 1, F(1, 2), F(1, 6), F(1, 24)]
    assert [frac_tex(c) for c in GF.exp_coeffs(4)] == ["1", "1", r"\frac{1}{2}", r"\frac{1}{6}", r"\frac{1}{24}"]
    for x in (0.5, 1.0):
        assert abs(sum(float(c) * x ** k for k, c in enumerate(GF.exp_coeffs(20))) - math.exp(x)) < 1e-12
    assert geometric(2, 3) == [1, 2, 4]
    h = GF.half_binom_series(3)
    assert h == [1, F(1, 2), F(3, 8), F(5, 16)]
    assert [poly_tex(P[n]) for n in range(4)] == ["1", "x", r"\frac{3x^2-1}{2}", r"\frac{5x^3-3x}{2}"]
    assert poly_tex(P[5]) == r"\frac{63x^5-70x^3+15x}{8}"
    assert P == GF.legendre_bonnet(5)
    assert [h[1] * (-1), 0, h[2] * 4] == P[2]              # z²：½·(−z²) 與 ⅜·(2xz)²
    assert all(GF.poly_eval(P[n], 1) == 1 for n in range(6))
    assert [f"{GF.ex22_integral(z):.3f}" for z in Z22] == ["2.000"] * 3
    assert [GF.poly_int11(P[n]) for n in range(4)] == [2, 0, 0, 0]
    assert GF.ex23_coeffs(3) == [2, F(2, 3), F(2, 5), F(2, 7)]
    g = gram(3)
    assert all(g[n][m] == (F(2, 2 * n + 1) if n == m else 0) for n in range(4) for m in range(4))
    s = sum(float(c) * 0.25 ** k for k, c in enumerate(GF.ex23_coeffs(40)))
    assert abs(GF.ex23_closed(0.5, 0.5) - s) < 1e-9


# ---------------------------------------------------------------- 小工具
def eq(tex, tag=None, font_size=32, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return m
    return VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：無窮多個數，能不能裝進一個式子裡？", wait=3.5)


def part_exp(s):
    cs = GF.exp_coeffs(4)
    seq = Tex(r",\ ".join(frac_tex(c) for c in cs) + r",\ \dots", font_size=44).set_color(COEF_COLOR)
    seq.move_to(np.array([0, 2.5, 0]))
    terms = [Tex("1", font_size=44), Tex(r"\frac{x}{1!}", font_size=44)]
    terms += [Tex(rf"\frac{{x^{k}}}{{{k}!}}", font_size=44) for k in range(2, 5)]
    parts = [Tex("e^x=", font_size=44)]
    for k, t in enumerate(terms):
        if k:
            parts.append(Tex("+", font_size=44))
        parts.append(t)
    parts.append(Tex(r"+\cdots", font_size=44))
    row = VGroup(*parts).arrange(RIGHT, buff=0.18).move_to(np.array([0, 1.2, 0]))
    ylab = row.get_bottom()[1] - 0.6
    labs = VGroup(*[Tex(frac_tex(c), font_size=36).set_color(COEF_COLOR).move_to(np.array([t.get_center()[0], ylab, 0]))
                    for t, c in zip(terms, cs)])
    note = mixed("「", "$e^x$", " 生成數列 ", "$1/n!$", "」", font_size=30)
    note.move_to(np.array([0, -0.85, 0]))
    g = geometric(2, 3)
    w105 = VGroup(Tex(r"\frac{1}{2^N}(1+e^y)^N", font_size=34), zh("W105", font_size=20, color=D.MUTED))
    w111 = VGroup(Tex(rf"\frac{{1}}{{1-2z}}={g[0]}+{g[1]}z+{g[2]}z^2+\cdots", font_size=34),
                  zh("W111", font_size=20, color=D.MUTED))
    for w in (w105, w111):
        w.arrange(DOWN, buff=0.15)
    recap = VGroup(zh("第 1 章已見過兩個：", font_size=26, color=GREY_B), w105, w111).arrange(RIGHT, buff=0.6)
    recap.move_to(np.array([0, -2.2, 0]))
    assert_in_frame(seq, row, labs, note, recap)
    assert_apart(seq, row)
    assert_apart(row, labs)
    assert_apart(labs, note)
    assert_apart(note, recap)

    s.play(FadeIn(seq, shift=0.2 * DOWN), run_time=1.0)
    s.say("例如這一串數：", "$1,1,\\frac12,\\frac16,\\frac1{24},\\dots$", wait=2.5)
    s.play(Write(row), run_time=2.0)
    s.play(LaggedStart(*[FadeIn(m, shift=0.2 * DOWN) for m in labs], lag_ratio=0.2), run_time=1.2)
    s.say("每個數放在 ", "$x^n$", " 前面、全部加起來：就是 ", "$e^x$", wait=3.0)
    s.play(FadeIn(note), run_time=0.8)
    s.say("我們說：", "$e^x$", " 是數列 ", "$1/n!$", " 的生成函數", wait=3.0)
    s.play(FadeIn(recap, shift=0.2 * UP), run_time=1.0)
    s.say("第 1 章已見過兩個：動差生成函數與等比級數", wait=3.5)
    s.clear_stage()


def part_legendre(s):
    h = GF.half_binom_series(3)
    l1 = eq(r"\frac{1}{\sqrt{z^2-2zx+1}}=\sum_{n\ge0}P_n(x)\,z^n", "(2.1)", font_size=36)
    l2 = mixed("令 ", "$t=2xz-z^2$", "：", font_size=26)
    l3 = Tex(r"(1-t)^{-1/2}=1+" + "+".join(frac_tex(h[k]) + "t" + (f"^{k}" if k > 1 else "") for k in range(1, 4))
             + r"+\cdots", font_size=32)
    l4 = mixed("收集 ", "$z^2$", "：", rf"${frac_tex(h[1])}\cdot(-1)+{frac_tex(h[2])}\cdot4x^2={poly_tex(P[2])}$",
               font_size=26)
    pt = [Tex(rf"P_{n}={poly_tex(P[n])}", font_size=32).set_color(P_COLORS[n]) for n in range(4)]
    p01 = VGroup(pt[0], pt[1]).arrange(RIGHT, buff=0.9)
    p23 = VGroup(pt[2], pt[3]).arrange(RIGHT, buff=0.6)
    col = VGroup(l1, l2, l3, l4, p01, p23)
    for m in col:
        fit(m, 6.9)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.32), -6.5, 2.95)

    ax = D.plot_axes((-1.2, 1.2, 0.5), (-1.2, 1.2, 0.5), width=4.2, height=3.4, x_ticks=[-1, 1], y_ticks=[-1, 1],
                     font_size=20, x_label="x")
    ax.move_to(np.array([4.45, 0.9, 0]))
    curves = VGroup(*[ax.get_graph(lambda x, n=n: float(GF.poly_eval(P[n], x)), x_range=[-1, 1, 0.01])
                      .set_stroke(P_COLORS[n], 3) for n in range(4)])
    dot1 = Dot(ax.c2p(1, 1), radius=0.08).set_fill(YELLOW)
    dlab = Tex(r"P_n(1)=1", font_size=28).set_color(YELLOW).next_to(dot1, UP, buff=0.15)
    expl = fit(Tex(r"x=1:\ \frac{1}{1-z}=1+z+z^2+\cdots", font_size=32), 5.0)
    expl.move_to(np.array([4.25, -1.6, 0]))
    err = D.erratum_box("p. 26，(2.1) 後", r"\text{right side},\ y^5", r"\text{left side},\ z^5", width=5.8)
    place_col(err, -6.5, l1.get_bottom()[1] - 0.3)
    p5 = Tex(rf"P_5={poly_tex(P[5])}", font_size=32)
    place_col(p5, -6.3, err.get_bottom()[1] - 0.2)
    assert_in_frame(col, ax, dlab, expl, err, p5)
    assert_apart(col, ax)
    assert_apart(col, expl)
    assert_apart(ax, expl)
    assert_apart(p5, p01)
    assert_apart(err, ax)

    s.say("物理裡常見的例子：Legendre 多項式", wait=0.3)
    s.play(Write(l1[0]), FadeIn(l1[1]), run_time=1.8)
    s.wait(1.5)
    s.say("左邊對 ", "$z$", " 展開，", "$z^n$", " 的係數就叫 ", "$P_n(x)$", wait=2.5)
    s.play(FadeIn(l2), FadeIn(l3), run_time=1.0)
    s.say("令 ", "$t=2xz-z^2$", "，用二項式展開 ", "$(1-t)^{-1/2}$", wait=3.0)
    s.play(FadeIn(l4), run_time=0.8)
    s.say("收集 ", "$z^2$", " 的項：兩個來源加起來", wait=3.5)
    s.play(LaggedStart(*[FadeIn(m) for m in pt], lag_ratio=0.3), run_time=1.5)
    s.say("得到 ", "$P_0$", " 到 ", "$P_3$", "；", "$P_2$", " 就是練習 2.1", wait=3.0)
    s.play(FadeIn(ax), run_time=0.6)
    s.play(LaggedStart(*[ShowCreation(c) for c in curves], lag_ratio=0.4), run_time=2.4)
    s.say("畫在 ", "$-1\\le x\\le1$", " 上：顏色和左邊的式子一樣", wait=2.0)
    s.say("先猜：", "$P_n(1)$", " 等於多少？", wait=3.0)
    s.play(FadeIn(dot1, scale=0.5), FadeIn(dlab), run_time=0.8)
    s.say("全部通過同一點：", "$P_n(1)=1$", wait=1.5)
    s.play(FadeIn(expl, shift=0.2 * UP), run_time=0.8)
    s.say("$x=1$", " 時左邊是 ", "$\\frac{1}{1-z}=1+z+z^2+\\cdots$", wait=3.5)
    s.play(FadeOut(VGroup(l2, l3, l4)), run_time=0.6)
    s.play(FadeIn(err, shift=0.2 * UP), run_time=0.8)
    s.say("書上說：把右邊對 ", "$y$", " 展開到五次，得到 ", "$P_5$", wait=3.0)
    s.say("應該是：把左邊對 ", "$z$", " 展開到五次", wait=2.5)
    s.play(FadeIn(p5, shift=0.2 * UP), run_time=0.8)
    s.say("這樣得到的 ", "$P_5$", "，和 Bonnet 遞迴算的一樣", wait=3.0)
    s.clear_stage()


def part_ex22(s):
    vals = [f"{GF.ex22_integral(z):.3f}" for z in Z22]
    a = Tex(r"\int_{-1}^{1}\frac{dx}{\sqrt{z^2-2zx+1}}=\frac{1+z-|1-z|}{z}=2", font_size=38)
    an = mixed("數值：", "$z=" + r",\ ".join(str(z) for z in Z22) + "$", " 都得到 ", f"${vals[0]}$",
               font_size=24, color=GREY_B)
    b = mixed("右邊：", "$\\int_{-1}^{1}\\sum_{n}P_n(x)z^n\\,dx=\\sum_{n}z^n\\int_{-1}^{1}P_n\\,dx$", font_size=28)
    rows = [["n", "0", "1", "2", "3"], [r"\int_{-1}^{1}P_n\,dx"] + [frac_tex(GF.poly_int11(P[n])) for n in range(4)]]
    tab = D.tex_table(rows, [2.6, 0.9, 0.9, 0.9, 0.9], row_h=0.9, font_size=32)
    D.table_cell(tab, 1, 1, 5).set_color(D.EXACT_COLOR)
    col = VGroup(a, an, b).arrange(DOWN, buff=0.4).move_to(np.array([0, 1.6, 0]))
    tab.shift(np.array([-tab.get_width() / 2, col.get_bottom()[1] - 0.5, 0]))
    assert_in_frame(col, tab)
    assert_apart(col, tab)

    s.say("練習 2.2：把 (2.1) 兩邊對 ", "$x$", " 從 ", "$-1$", " 積到 ", "$1$", wait=0.3)
    s.play(Write(a), run_time=1.8)
    s.wait(1.5)
    s.play(FadeIn(an), run_time=0.6)
    s.say("左邊算出來是 2：每個 ", "$|z|<1$", " 都一樣", wait=3.0)
    s.play(FadeIn(b), run_time=0.8)
    s.say("右邊逐項積分：", "$z^n$", " 的係數是 ", "$\\int P_n\\,dx$", wait=3.0)
    s.play(FadeIn(tab, shift=0.2 * UP), run_time=1.0)
    s.say("比較係數：只有 ", "$\\int P_0\\,dx=2$", "，其餘全是 0", wait=3.5)
    s.clear_stage()


def part_ex23(s):
    c = GF.ex23_coeffs(3)
    g = gram(3)
    a = fit(Tex(r"\int_{-1}^{1}\frac{dx}{\sqrt{y^2-2yx+1}\sqrt{z^2-2zx+1}}"
                r"=\frac{1}{\sqrt{yz}}\ln\frac{1+\sqrt{yz}}{1-\sqrt{yz}}", font_size=34), 12.0)
    a.move_to(np.array([0, 2.55, 0]))
    b = Tex("=" + "+".join([frac_tex(c[0]), frac_tex(c[1]) + "yz"] + [frac_tex(c[k]) + f"(yz)^{k}" for k in (2, 3)])
            + r"+\cdots", font_size=36)
    b.move_to(np.array([0, 1.45, 0]))
    cc = mixed("左邊逐項展開：", "$\\sum_{n,m}y^nz^m\\int_{-1}^{1}P_nP_m\\,dx$", font_size=26)
    nn = mixed("只有 ", "$y^nz^n$", " 的項 ", "$\\Rightarrow$", " 正交：", font_size=26)
    box = boxed(Tex(r"\int_{-1}^{1}P_nP_m\,dx=\frac{2\,\delta_{nm}}{2n+1}", font_size=38), color=D.GAUSS_COLOR, buff=0.2)
    left = VGroup(cc, nn, box).arrange(DOWN, aligned_edge=LEFT, buff=0.45)
    place_col(left, -6.5, 0.6)
    rows = [[r"n\backslash m", "0", "1", "2", "3"]] + [[str(n)] + [frac_tex(g[n][m], small=True) for m in range(4)] for n in range(4)]
    tab = D.tex_table(rows, [1.0, 0.85, 0.85, 0.85, 0.85], row_h=0.6, font_size=28)
    tab.shift(np.array([2.15, 0.25, 0]))
    for n in range(4):
        D.table_cell(tab, n + 1, 0, 5).set_color(YELLOW)
        for m in range(4):
            D.table_cell(tab, n + 1, m + 1, 5).set_color(D.EXACT_COLOR if n == m else GREY_C)
    ttl = Tex(r"\int_{-1}^{1}P_nP_m\,dx", font_size=28).next_to(tab, UP, buff=0.12)
    assert_in_frame(a, b, left, tab, ttl)
    assert_apart(a, b)
    assert_apart(b, ttl)
    assert_apart(left, tab)
    assert_apart(b, left)

    s.say("練習 2.3：兩個生成函數相乘，再積分", wait=0.3)
    s.play(Write(a), run_time=2.0)
    s.wait(2.0)
    s.play(FadeIn(b, shift=0.2 * DOWN), run_time=0.8)
    s.say("展開成 ", "$yz$", " 的冪級數：", "$2,\\frac23,\\frac25,\\frac27,\\dots$", wait=3.0)
    s.say("只出現 ", "$y^nz^n$", "：沒有 ", "$y^2z^3$", " 這種交叉項", wait=3.0)
    s.play(FadeIn(cc), run_time=0.8)
    s.say("左邊逐項展開：", "$y^nz^m$", " 的係數是 ", "$\\int P_nP_m\\,dx$", wait=3.0)
    s.play(FadeIn(nn), FadeIn(box, shift=0.2 * UP), run_time=1.0)
    s.say("所以 ", "$n\\neq m$", " 時積分是 0：", "$P_n$", " 互相正交", wait=3.0)
    s.play(FadeIn(ttl), FadeIn(tab, shift=0.2 * UP), run_time=1.0)
    s.say("對角線 ", "$2,\\frac23,\\frac25,\\frac27$", "，其餘全是 0", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("生成函數把整串數列打包成一個函數", wait=3.0)
    s.say("要哪一項就展開取係數，還能一次證明整串數的性質", wait=4.0)


class W201WhatIsGF(ChuScene):
    def construct(self):
        self.show_header("W201　生成函數：一個函數裝下整個數列")
        check_math()
        part_question(self)
        part_exp(self)
        part_legendre(self)
        part_ex22(self)
        part_ex23(self)
        part_close(self)
        self.wait(1.0)
