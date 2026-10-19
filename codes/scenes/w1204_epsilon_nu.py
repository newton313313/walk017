"""
W1204 ε 展開求 ν（§12.2.5–12.2.6，pp. 296–300，(12.39)–(12.60)，Fig. 12.13，Ex. 12.2–12.4）

問：耦合 u 最後流到哪裡？
  1. (12.39)–(12.41)（r = 0）：u(ℓ) 的閉式＋數值點（U_CASES、ε = 1）→ u* = 1/36；Fig. 12.13 的數線（Ex. 12.3）。
  2. r = t − 6u ⇒ (12.45) ⇒ (12.51) t(ℓ)；t0 = 1e−4、u0 = 0.01、ℓ = 3：0.0203；筆誤框 1（(12.48)）。
  3. (12.52)–(12.60)：ν = 1/(2 − ε/3) ≈ ½ + ε/12 = 7/12；Ex. 12.4 臨界點 r_c；筆誤框 2（Ex. 12.4）。
  4. ν 的三種讀法（表 4 列 × 4 欄）＋ε = 1 的矩陣＋三條 ν(ε) 曲線。
  5. 平均場 0.5 → 一階 0.5833 → 已知 0.630（引用）；筆誤框 3（(12.60) 後）；收束。
畫面上的數字都由 utils.walk_eps 算出，check_math() 核對。
模組層級的畫圖函式：u_axes／u_curve（u(ℓ)）、fig1213（數線）、nu_axes／nu_curve、compare_line。
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
U_COLORS = (TEAL, ORANGE)
NU_COLORS = (YELLOW, ORANGE, TEAL)          # 一階、不展開、矩陣
ELLS = (1, 2, 5)
EPS1 = 1.0
NU_N = 30
U_N = 40


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def f6(x):
    return f"{float(x):.6f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {
            "u": {u0: [E.u_analytic(l, u0, e) for l in ELLS] for u0, e in E.U_CASES},
            "un": {u0: [E.u_numeric(l, u0, e) for l in ELLS] for u0, e in E.U_CASES},
            "ratio": E.u_analytic(2, 1e-6, EPS1) / 1e-6,
            "ratio_n": E.u_numeric(2, 1e-6, EPS1) / 1e-6,
            "t": E.t_closed(3, 1e-4, 0.01, 1.0),
            "tn": E.t_numeric(3, 1e-4, 0.01, 1.0),
            "i_book": E.int_1248(0.001, 0.005, 1.0, 3, book=True),
            "i_ok": E.int_1248(0.001, 0.005, 1.0, 3),
            "i_49": E.int_1249(0.001, 0.005, 1.0, 3),
            "rc": E.r_critical(1e-4, 0.1),
            "rc_lin": E.r_critical_linear(1e-4, 0.1),
            "nu": {e: (E.nu_first(e), E.nu_unexpanded(e), E.nu_eigen(e)) for e in E.NU_EPS},
            "eig": E.eigen(1.0),
            "jac": E.jacobian(*E.fixed_point(1.0), 1.0),
        }


def check_math():
    nb = numbers()
    # §1 u(ℓ)
    assert E.U_CASES == ((0.005, 1.0), (0.06, 1.0))
    assert [f4(v) for v in nb["u"][0.005]] == ["0.0104", "0.0172", "0.0270"]
    assert [f4(v) for v in nb["u"][0.06]] == ["0.0346", "0.0300", "0.0279"]
    for u0 in (0.005, 0.06):
        assert all(abs(a - b) < 1e-8 for a, b in zip(nb["u"][u0], nb["un"][u0]))
    assert f4(1 / 36) == "0.0278"
    assert f4(nb["ratio"]) == "7.3874" and f4(math.e ** 2) == "7.3891" and abs(nb["ratio"] - nb["ratio_n"]) < 1e-5
    # 數線：u = 0 斜率 +ε、u* 斜率 −ε（du/dℓ 的導數）
    us = 1 / 36
    assert abs(E.beta(0, 1e-7, EPS1)[1] / 1e-7 - EPS1) < 1e-4
    h = 1e-7
    slope = (E.beta(0, us + h, EPS1)[1] - E.beta(0, us - h, EPS1)[1]) / (2 * h)
    assert abs(slope + EPS1) < 1e-5
    # §2 t(ℓ)
    assert f4(nb["t"]) == "0.0203" and f4(nb["tn"]) == "0.0203" and abs(nb["t"] - nb["tn"]) < 1e-8
    assert E.AB(1) == (12, 36) and Fraction(12, 36) == Fraction(1, 3)
    assert f4(nb["i_book"]) == "81.7804" and f4(nb["i_ok"]) == "0.2207" and f4(nb["i_49"]) == "0.2207"
    # §3 ν、Ex. 12.4
    assert E.nu_unexpanded(1) == Fraction(3, 5) and E.nu_first(1) == Fraction(7, 12) and f4(Fraction(7, 12)) == "0.5833"
    assert E.nu_first(Fraction(1, 2)) == Fraction(1, 2) + Fraction(1, 24)
    assert f6(nb["rc"]) == "-0.000631" and f6(nb["rc_lin"]) == "-0.000632" and f6(-6 * 1e-4) == "-0.000600"
    assert abs(nb["rc_lin"] - (-12 * 1e-4 / (2 - 0.1))) < 1e-15
    # §4 表
    want = {Fraction(1, 10): ("0.5083", "0.5085", "0.5076"), Fraction(1, 2): ("0.5417", "0.5455", "0.5256"),
            Fraction(1): ("0.5833", "0.6000", "0.5272")}
    assert E.NU_EPS == (Fraction(1, 10), Fraction(1, 2), Fraction(1))
    for e in E.NU_EPS:
        assert tuple(f4(v) for v in nb["nu"][e]) == want[e], e
    J = nb["jac"]
    exact = ((Fraction(5, 3), 14), (Fraction(1, 21), -1))
    assert all(abs(J[i][j] - float(exact[i][j])) < 1e-12 for i in range(2) for j in range(2))
    assert [f4(v) for v in nb["eig"]] == ["1.8968", "-1.2301"]
    assert abs(nb["eig"][0] - (1 + math.sqrt(22)) / 3) < 1e-12
    # 曲線在 ε = 0 相切（都是 ½，斜率 1/12）
    for fn in (E.nu_first, E.nu_unexpanded):
        assert float(fn(0)) == 0.5
    assert abs(E.nu_eigen(0.0) - 0.5) < 1e-12 and abs((E.nu_eigen(1e-4) - 0.5) / 1e-4 - 1 / 12) < 1e-3
    # §5
    assert E.NU_3D_CITED[1] == 0.630


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


def tip_arrow(p_prev, p_end, color, length=0.25):
    u = (p_end - p_prev) / np.linalg.norm(p_end - p_prev)
    return Arrow(p_end - length * u, p_end, buff=0, thickness=4).set_color(color)


# ---------------------------------------------------------------- 畫圖元件
def u_axes():
    return D.plot_axes((0, 6.4, 1), (0, 0.068, 0.01), width=6.0, height=4.4, x_ticks=(1, 2, 3, 4, 5, 6),
                       y_ticks=(0.02, 0.04, 0.06), y_tick_labels=["0.02", "0.04", "0.06"],
                       x_label=r"\ell", y_label="u")


def u_curve(ax, u0, color):
    ps = [ax.c2p(t, E.u_analytic(t, u0, EPS1)) for t in np.linspace(0, 6, U_N)]
    return VMobject().set_points_smoothly(ps).set_stroke(color, 4)


def u_dots(ax, u0, color):
    nb = numbers()
    return VGroup(*[Dot(ax.c2p(t, v), radius=0.07).set_color(WHITE).set_stroke(color, 2)
                    for t, v in zip(ELLS, nb["un"][u0])])


def fig1213(width=5.4):
    """Fig. 12.13：u 軸上的兩個不動點；u = 0 往外推（不穩定）、u* 往內拉（穩定）。"""
    u0, u1 = -0.008, 0.045
    k = width / (u1 - u0)
    X = lambda u: (u - u0) * k * RIGHT
    us = 1 / 36
    axis = Arrow(X(u0), X(u1), buff=0, thickness=2).set_color(GREY_B)
    p0 = Dot(X(0), radius=0.1).set_color(RED)
    ps = Dot(X(us), radius=0.1).set_color(WHITE)
    arrows = VGroup(Arrow(X(-0.001), X(-0.007), buff=0, thickness=4), Arrow(X(0.002), X(0.011), buff=0, thickness=4),
                    Arrow(X(0.016), X(us - 0.002), buff=0, thickness=4),
                    Arrow(X(0.043), X(us + 0.002), buff=0, thickness=4)).set_color(OK_COLOR)
    for a in arrows:
        a.shift(0.28 * UP)
    l0 = tx("0", 24, GREY_A).next_to(p0, DOWN, buff=0.15)
    ls = tx(r"u^*=\tfrac{1}{36}", 24, GREY_A).next_to(ps, DOWN, buff=0.15)
    lu = tx("u", 24, GREY_B).next_to(axis.get_end(), DOWN, buff=0.15)
    s0 = mixed("斜率 ", "$+\\epsilon$", "：不穩定", font_size=20, color=RED).next_to(p0, UP, buff=0.55)
    s1 = mixed("斜率 ", "$-\\epsilon$", "：穩定", font_size=20, color=OK_COLOR).next_to(ps, UP, buff=0.55)
    return VGroup(axis, arrows, p0, ps, l0, ls, lu, s0, s1)


def nu_axes():
    return D.plot_axes((0, 1.1, 0.25), (0, 0.115, 0.025), width=4.6, height=3.8, x_ticks=(0.5, 1),
                       y_ticks=(0, 0.05, 0.1), y_tick_labels=["0.50", "0.55", "0.60"],
                       x_label=r"\epsilon", y_label=r"\nu")


def nu_curve(ax, k, color):
    fns = (lambda e: float(E.nu_first(e)), lambda e: float(1 / (2 - e / 3)), lambda e: E.nu_eigen(e))
    ps = [ax.c2p(e, fns[k](e) - 0.5) for e in np.linspace(0, 1, NU_N)]
    return VMobject().set_points_smoothly(ps).set_stroke(color, 4)


def compare_line(width=10.0):
    v0, v1 = 0.46, 0.67
    k = width / (v1 - v0)
    X = lambda v: (v - v0) * k * RIGHT
    axis = Arrow(X(v0), X(v1), buff=0, thickness=2).set_color(GREY_B)
    ticks = VGroup()
    for v in (0.5, 0.55, 0.6, 0.65):
        ticks.add(Line(X(v) + 0.08 * DOWN, X(v) + 0.08 * UP).set_stroke(GREY_B, 2))
        ticks.add(tx(f"{v:.2f}", 20, GREY_B).next_to(X(v), DOWN, buff=0.2))
    pts = VGroup()
    for v, lab, col in ((0.5, "平均場", GREY_A), (float(E.nu_first(1)), "一階 ε 展開", YELLOW),
                        (E.NU_3D_CITED[1], "三維已知（引用）", OK_COLOR)):
        d = Dot(X(v), radius=0.12).set_color(col)
        name = zh(lab, font_size=22, color=col).next_to(d, UP, buff=0.8)
        val = tx({0.5: "0.5", E.NU_3D_CITED[1]: "0.630"}.get(v, f4(v)), 26, col).next_to(name, UP, buff=0.12)
        pts.add(VGroup(d, name, val))
    arr = VGroup(CurvedArrow(X(0.5) + 0.25 * UP, X(float(E.nu_first(1))) + 0.25 * UP, angle=-PI / 4),
                 CurvedArrow(X(float(E.nu_first(1))) + 0.25 * UP, X(0.630) + 0.25 * UP, angle=-PI / 4)
                 ).set_color(GREY_A)
    return VGroup(axis, ticks, pts, arr)


# ---------------------------------------------------------------- 各段
def part_u(s):
    nb = numbers()
    ax = u_axes()
    ax.move_to(np.array([-3.4, 0.1, 0]))
    curves = VGroup(*[u_curve(ax, u0, c) for (u0, _), c in zip(E.U_CASES, U_COLORS)])
    dots = VGroup(*[u_dots(ax, u0, c) for (u0, _), c in zip(E.U_CASES, U_COLORS)])
    star = DashedLine(ax.c2p(0, 1 / 36), ax.c2p(6.2, 1 / 36), dash_length=0.08).set_stroke(RED, 2)
    slab = tx(r"u^*=\tfrac{1}{36}=%s" % f4(1 / 36), 22, RED).next_to(ax.c2p(6.2, 1 / 36), UP, buff=0.08).shift(0.9 * LEFT + 0.15 * UP)
    clab = VGroup(tx(r"u_0=0.005", 22, U_COLORS[0]).next_to(ax.c2p(0.9, 0.0104), RIGHT, buff=0.2).shift(0.12 * DOWN),
                  tx(r"u_0=0.06", 22, U_COLORS[1]).next_to(ax.c2p(0.3, 0.06), RIGHT, buff=0.15))
    tag = zh("(12.39)–(12.41)　先看 r = 0", font_size=22, color=GREY_B)
    eqs = VGroup(tx(r"\frac{du}{d\ell}=\epsilon u-36u^2", 32),
                 tx(r"u(\ell)=\frac{\epsilon u_0e^{\epsilon\ell}}{\epsilon+36u_0(e^{\epsilon\ell}-1)}", 32)
                 ).arrange(DOWN, buff=0.3)
    rows = [[r"\ell", r"u_0=0.005", r"u_0=0.06"]] + [[str(t), f4(nb["u"][0.005][i]), f4(nb["u"][0.06][i])]
                                                     for i, t in enumerate(ELLS)]
    tab = D.tex_table(rows, [0.8, 1.9, 1.9], row_h=0.42, font_size=26)
    fig = fig1213(5.0)
    right = VGroup(tag, eqs, tab, fig).arrange(DOWN, buff=0.32)
    place_in(right, 0.5, 6.7, -2.85, 3.1)
    assert_in_frame(ax, right, slab, clab)
    assert_apart(VGroup(ax, slab, clab), right)
    assert_column(tag, eqs, tab, fig)
    assert_table(tab, 3)
    assert_apart(clab[0], clab[1])
    s.say("問：耦合 ", "$u$", " 最後流到哪裡？", wait=2.0)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.0)
    s.say("書上留著 ", "$d$", "，本系列併進 ", "$u$", "：", "$36d\\,u^2\\to36u^2$", wait=2.2)
    s.play(Write(eqs[1]), run_time=1.0)
    s.say("Ex. 12.2：分離變數積分，得到 (12.41)", wait=2.3)
    s.play(ShowCreation(ax), ShowCreation(star), FadeIn(slab), run_time=1.0)
    s.play(ShowCreation(curves), FadeIn(clab), run_time=1.6)
    s.say("$\\epsilon=1$", "：小的 ", "$u_0$", " 往上爬、大的 ", "$u_0$", " 往下降", wait=2.3)
    s.play(FadeIn(dots), FadeIn(tab), run_time=1.0)
    s.say("圓點是數值積分，落在閉式上；兩條都趨向 ", "$u^*$", wait=2.4)
    s.play(FadeIn(fig[0]), FadeIn(fig[2:7]), run_time=0.8)
    s.play(*[GrowArrow(a) for a in fig[1]], FadeIn(fig[7]), FadeIn(fig[8]), run_time=1.0)
    s.say("Fig. 12.13：", "$u=0$", " 斜率 ", "$+\\epsilon$", " 不穩定；", "$u^*$", " 斜率 ", "$-\\epsilon$", " 穩定",
          wait=2.4)
    s.say("Ex. 12.3：", "$u_0=10^{-6}$", "，", "$u(2)/u_0=%s\\approx e^2=%s$" % (f4(nb["ratio"]), f4(math.e ** 2)),
          wait=2.4)
    s.clear_stage()


def part_t(s):
    nb = numbers()
    tag = zh("(12.42)–(12.51)", font_size=22, color=GREY_B)
    e1 = mixed("$r=t-6u$", "（書上 ", "$t-6ud$", "）", font_size=26)
    e2 = tx(r"\frac{dt}{d\ell}=2t-12u\,t", 34)
    e3 = tx(r"t(\ell)=t(0)\,e^{2\ell}\Big[\frac{\epsilon}{\epsilon+36u_0(e^{\epsilon\ell}-1)}\Big]^{1/3}", 34)
    e4 = VGroup(tx(r"t_0=10^{-4},\ u_0=0.01,\ \ell=3", 30, OK_COLOR),
                mixed("閉式 ", "$%s$" % f4(nb["t"]), "＝數值 ", "$%s$" % f4(nb["tn"]), font_size=26, color=OK_COLOR)
                ).arrange(DOWN, buff=0.2)
    e5 = mixed("$\\ell$", " 很大：", "$t\\propto e^{(2-\\epsilon/3)\\ell}$", font_size=26, color=YELLOW)
    col = VGroup(tag, e1, e2, e3).arrange(DOWN, buff=0.4)
    place_in(col, -6.7, 0.0, -0.1, 3.1)
    col2 = VGroup(e4, e5).arrange(DOWN, buff=0.45)
    place_in(col2, 0.4, 6.7, 0.4, 2.6)
    box = D.erratum_box("p. 298，(12.48)", r"12d\int_{u_0}^{u}\frac{du'}{\epsilon u'-36du'^2}",
                        r"12d\int_{u_0}^{u}\frac{u'\,du'}{\epsilon u'-36du'^2}", width=9.0)
    place_in(box, -4.8, 4.8, -2.85, -0.4)
    assert_in_frame(col, col2, box)
    assert_column(tag, e1, e2, e3, box)
    assert_column(e4, e5, box)
    assert_apart(col, col2)
    s.say("接著看 ", "$r$", "：先換成 ", "$t=r+6u$", wait=2.0)
    s.play(FadeIn(tag), FadeIn(e1), run_time=0.8)
    s.play(Write(e2), run_time=1.0)
    s.say("(12.45)：", "$t$", " 的流只差一個因子 ", "$1-6u$", wait=2.3)
    s.play(Write(e3), run_time=1.2)
    s.say("代入 ", "$u(\\ell)$", " 積分，得到 (12.51)", wait=2.3)
    s.play(FadeIn(e4), run_time=0.8)
    s.say("$t_0=10^{-4}$", "、", "$u_0=0.01$", "、", "$\\ell=3$", "：閉式＝數值積分", wait=2.3)
    s.play(FadeIn(e5), run_time=0.8)
    s.say("$u\\to u^*=\\tfrac1{36}$", " 之後，", "$t$", " 長得比 ", "$e^{2\\ell}$", " 慢", wait=2.2)
    s.play(FadeIn(box), run_time=0.8)
    s.say("$u_0=0.001\\to0.005$", "（", "$\\epsilon=1,\\ d=3$", "）：照書上 ", "$%s$" % f4(nb["i_book"]), wait=2.4)
    s.say("應為 ", "$%s$" % f4(nb["i_ok"]), "＝(12.49)；(12.47) 左邊也少了 ", "$12d$", wait=2.4)
    s.clear_stage()


def part_nu(s):
    nb = numbers()
    tag = zh("(12.52)–(12.60)", font_size=22, color=GREY_B)
    eqs = VGroup(tx(r"\xi(\ell^*)=1,\qquad r(\ell^*)=\mathrm{Dist}", 32),
                 tx(r"e^{\ell^*}=(T-T_c^*)^{-1/(2-\epsilon/3)}", 32),
                 tx(r"\nu=\frac{1}{2-\epsilon/3}\approx\frac12+\frac{\epsilon}{12}", 34),
                 tx(r"\epsilon=1:\ \nu=\frac{7}{12}=%s" % f4(Fraction(7, 12)), 36, YELLOW)).arrange(DOWN, buff=0.3)
    left = VGroup(tag, eqs).arrange(DOWN, buff=0.3)
    place_in(left, -6.7, -0.5, -1.25, 3.1)
    head = mixed("Ex. 12.4：臨界點移到 ", "$r_c\\approx-6u_0$", font_size=24)
    sub = mixed("（書上 ", "$-6u_0d$", "）", "$u_0=10^{-4},\\ \\epsilon=0.1$", font_size=22, color=GREY_A)
    rows = [[zh("二分法", font_size=22, color=YELLOW), "%s" % f6(nb["rc"])],
            [mixed("特徵向量 ", "$-\\frac{12u_0}{2-\\epsilon}$", font_size=22), f6(nb["rc_lin"])],
            [mixed("一階 ", "$-6u_0$", font_size=22), f6(-6e-4)]]
    tab = D.tex_table(rows, [2.6, 1.9], row_h=0.62, font_size=26)
    box = D.erratum_box("p. 299，Ex. 12.4", r"\text{(12.53) follows from (12.51)}",
                        r"\text{(12.52) follows from (12.51)}", width=11.0)
    place_in(box, -6.0, 6.0, -2.85, -1.45)
    right = VGroup(head, sub, tab[0]).arrange(DOWN, buff=0.32)
    place_in(right, 0.2, 6.7, -1.25, 3.1)
    assert_in_frame(left, right, box)
    assert_apart(left, right)
    assert_column(tag, *eqs)
    assert_column(head, sub, tab[0])
    assert_apart(VGroup(left, right), box)
    assert_table(tab, 2)
    s.say("關聯長度：流到 ", "$\\xi(\\ell^*)=1$", " 為止", wait=2.0)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.0)
    s.say("那時 ", "$r$", " 離臨界面一段固定距離 ", "$\\mathrm{Dist}$", wait=2.3)
    s.play(Write(eqs[1]), run_time=1.0)
    s.say("$t$", " 以 ", "$e^{(2-\\epsilon/3)\\ell}$", " 長大，解出 ", "$e^{\\ell^*}$", wait=2.3)
    s.play(Write(eqs[2]), run_time=1.0)
    s.say("$\\xi=e^{\\ell^*}$", "，所以 ", "$\\nu=1/(2-\\epsilon/3)$", wait=2.3)
    s.play(Write(eqs[3]), run_time=0.8)
    s.say("(12.60)：", "$\\epsilon=1$", " 得 ", "$\\nu=\\tfrac{7}{12}=%s$" % f4(Fraction(7, 12)), wait=2.2)
    s.play(FadeIn(head), FadeIn(sub), run_time=0.8)
    s.play(FadeIn(tab[0]), run_time=0.8)
    s.say("臨界點不在 ", "$r=0$", "：二分法與特徵向量幾乎一樣", wait=2.4)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上寫 (12.53)，但 (12.53) 在練習之後才定義", wait=2.2)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    heads = [r"\epsilon", r"\tfrac12+\tfrac{\epsilon}{12}", r"\tfrac{1}{2-\epsilon/3}", r"1/\lambda_r"]
    labs = VGroup(*[zh(t, font_size=20, color=c) for t, c in zip(("一階", "不展開", "線性化矩陣"), NU_COLORS)])
    rows = [heads] + [[{Fraction(1, 10): "0.1", Fraction(1, 2): "0.5", Fraction(1): "1"}[e]] +
                      [f4(v) for v in nb["nu"][e]] for e in E.NU_EPS]
    widths = [0.9, 1.9, 1.9, 1.7]
    tab = D.tex_table(rows, widths, row_h=0.55, font_size=28)
    for j, c in enumerate(NU_COLORS):
        tab[0][j + 1].set_color(c)
    xs = np.cumsum([0] + widths)
    for j, m in enumerate(labs):
        m.move_to(np.array([(xs[j + 1] + xs[j + 2]) / 2, 0.28, 0]))
    J = nb["jac"]
    mat = tx(r"\epsilon=1:\ \begin{pmatrix}5/3&14\\1/21&-1\end{pmatrix}", 30)
    eig = tx(r"\lambda=%s,\ %s" % (f4(nb["eig"][0]), f4(nb["eig"][1])), 30, TEAL)
    block = VGroup(VGroup(labs, tab), mat, eig).arrange(DOWN, buff=0.35)
    place_in(block, -6.7, -0.4, -2.85, 3.1)
    ax = nu_axes()
    ax.move_to(np.array([3.3, 0.2, 0]))
    curves = VGroup(*[nu_curve(ax, k, c) for k, c in enumerate(NU_COLORS)])
    assert_in_frame(block, ax, curves)
    assert_apart(block, ax)
    assert_table(tab, 4)
    assert_column(VGroup(labs, tab), mat, eig)
    assert J.shape == (2, 2)
    s.play(FadeIn(tab[1]), FadeIn(tab[0][:4]), FadeIn(labs), run_time=1.0)
    s.say("$\\nu$", " 的三種讀法：一階、不展開、線性化矩陣", wait=2.0)
    s.play(LaggedStart(*[FadeIn(tab[0][4 * i:4 * i + 4]) for i in range(1, 4)], lag_ratio=0.3), run_time=1.6)
    s.say("$\\epsilon$", " 越大，三欄分得越開", wait=2.0)
    s.play(FadeIn(mat), FadeIn(eig), run_time=0.8)
    s.say("$\\epsilon=1$", " 的 WF 點：矩陣的特徵值 ", "$%s$" % f4(nb["eig"][0]), "、", "$%s$" % f4(nb["eig"][1]),
          wait=2.2)
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(ShowCreation(curves), run_time=1.6)
    s.say("三條 ", "$\\nu(\\epsilon)$", " 在 ", "$\\epsilon=0$", " 相切：都是 ", "$\\tfrac12+\\tfrac{\\epsilon}{12}$", wait=2.2)
    s.say("三者只在 ", "$\\epsilon$", " 的一次相同——一迴圈只準到 ", "$O(\\epsilon)$", wait=2.4)
    s.clear_stage()


def part_compare(s):
    line = compare_line(10.0)
    place_in(line, -6.0, 6.0, 0.0, 3.0)
    box = D.erratum_box("p. 299，(12.60) 後", r"\nu\sim0.61", r"\nu\approx0.630", width=5.6)
    place_in(box, -3.0, 3.0, -2.85, -0.6)
    assert_in_frame(line, box)
    assert_apart(line, box)
    assert_apart(line[2][0], line[2][1])
    assert_apart(line[2][1], line[2][2])
    s.play(FadeIn(line[0]), FadeIn(line[1]), FadeIn(line[2][0]), run_time=0.8)
    s.say("和已知值比：平均場 ", "$\\nu=\\tfrac12$", wait=2.0)
    s.play(ShowCreation(line[3][0]), FadeIn(line[2][1]), run_time=1.0)
    s.say("一階 ", "$\\epsilon$", " 展開：", "$0.5833$", wait=2.0)
    s.play(ShowCreation(line[3][1]), FadeIn(line[2][2]), run_time=1.0)
    s.say("三維 Ising 已知約 ", "$0.630$", "（引用）：大幅改善", wait=2.2)
    s.play(FadeIn(box), run_time=0.8)
    s.say("$0.630$", " 是文獻值（只引用），不是一迴圈算出來的", wait=2.2)
    s.say("Wilson–Fisher 點把 ", "$\\nu$", " 從平均場的 ", "$\\tfrac12$", " 推到 0.58", wait=2.4)
    s.say("——第一次「算」出不是 ", "$\\tfrac12$", " 的臨界指數", wait=2.6)


class W1204EpsilonNu(ChuScene):
    def construct(self):
        self.show_header("W1204　ε 展開求 ν")
        check_math()
        part_u(self)
        part_t(self)
        part_nu(self)
        part_table(self)
        part_compare(self)
        self.wait(1.0)
