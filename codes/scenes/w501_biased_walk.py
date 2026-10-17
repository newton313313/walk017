"""
W501 有偏向的漫步（§5.1–5.1.1，pp. 95–98，(5.1)–(5.10)，Fig. 5.1，Ex. 5.1）

問：不公正的硬幣（向右的機率 p = 0.6）。醉漢走 10 步，平均到哪裡？分布比公正硬幣寬還是窄？
  1. 精確二項 (5.1)＝逐步遞迴 (5.2)：平均 2、變異數 9.6（公正 10）；三根長條 vs 2×高斯 (5.7)。
  2. Fig. 5.1：N = 10, 20, 40 往右移、變寬；寬度 ∝ √(p(1−p))。
  3. 泰勒展開 (5.3)＋筆誤框；Ex. 5.1：D = 0.5（錯）vs 0.48；寬度比 √(4p(1−p))。
  4. (5.9)(5.10)「先走一步再亂走」：只平移、不變窄（對照表）。
畫面上的數字都由 utils.walk_var 算出，check_math() 核對。
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
from utils import walk_var as V
from utils import walk_draw as D

P_BIAS = 0.6
P_EXACT = Fraction(3, 5)
N0 = 10
NS = (10, 20, 40)
N_COL = {10: BLUE_C, 20: TEAL_C, 40: PURPLE_B}
B_WIND = 0.2
WIND_COLOR = GREEN_C
MARK = (0, 2, 4)


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


@lru_cache(maxsize=None)
def numbers():
    nb = {}
    P = V.biased_pmf(N0, P_EXACT)
    nb["P10"] = P
    nb["mom"] = V.moments(P)
    nb["fair_var"] = V.biased_var(N0, Fraction(1, 2))
    nb["tri"] = [(x, float(P[x]), 2 * float(V.biased_gauss(x, N0, P_BIAS))) for x in MARK]
    fig = {}
    for N in NS:
        Pn = V.biased_pmf(N, P_BIAS)
        xm = max(Pn, key=Pn.get)
        fig[N] = dict(P=Pn, mean=round(V.biased_mean(N, P_BIAS)), mean_f=V.biased_mean(N, P_BIAS), var=V.biased_var(N, P_BIAS), xm=xm, peak=Pn[xm],
                      gpeak=2 * float(V.biased_gauss(xm, N, P_BIAS)))
    nb["fig"] = fig
    nb["Dn"], nb["Dc"] = V.diffusion_D(P_BIAS, kind="naive"), V.diffusion_D(P_BIAS)
    nb["var_n"], nb["var_c"] = 2 * nb["Dn"] * N0, 2 * nb["Dc"] * N0
    nb["wid"] = (math.sqrt(nb["var_n"]), math.sqrt(nb["var_c"]))
    nb["wr"] = (V.width_ratio(0.6), V.width_ratio(0.9))
    sp = V.shift_pmf(N0, B_WIND)
    m = sum(x * w for x, w in sp)
    nb["shift"] = ({x: w for x, w in sp}, m, sum(x * x * w for x, w in sp) - m * m)
    return nb


def check_math():
    nb = numbers()
    P = nb["P10"]
    assert P == V.biased_dp(N0, P_EXACT) and sum(P.values()) == 1                      # (5.1) = (5.2)
    assert all(x % 2 == 0 for x in P) and len(P) == N0 + 1                              # 只有偶數格有人
    assert nb["mom"] == (2, Fraction(48, 5)) and nb["fair_var"] == 10
    assert nb["mom"] == (V.biased_mean(N0, P_EXACT), V.biased_var(N0, P_EXACT))
    assert [f4(a) for _, a, _ in nb["tri"]] == ["0.2007", "0.2508", "0.2150"]
    assert [f4(g) for _, _, g in nb["tri"]] == ["0.2091", "0.2575", "0.2091"]
    fig = nb["fig"]
    assert all(abs(fig[N]["mean_f"] - fig[N]["mean"]) < 1e-12 for N in NS)
    assert [fig[N]["mean"] for N in NS] == [2, 4, 8] and all(fig[N]["xm"] == fig[N]["mean"] for N in NS)
    assert [f"{fig[N]['var']:.1f}" for N in NS] == ["9.6", "19.2", "38.4"]
    assert [f4(fig[N]["peak"]) for N in NS] == ["0.2508", "0.1797", "0.1279"]
    assert [f4(fig[N]["gpeak"]) for N in NS] == ["0.2575", "0.1821", "0.1288"]
    assert all(min(fig[N]["P"]) >= -N for N in NS)
    assert V.taylor_residual(0.6) < 1e-5 < 1e-3 < V.taylor_residual(0.6, book=True)     # (5.3) 的號
    assert nb["Dn"] == 0.5 and f"{nb['Dc']:.2f}" == "0.48"
    assert abs(nb["var_n"] - 10) < 1e-12 and abs(nb["var_c"] - 9.6) < 1e-12
    assert abs(nb["var_c"] - V.biased_var(N0, P_BIAS)) < 1e-12
    assert (f"{nb['wid'][0]:.3f}", f"{nb['wid'][1]:.3f}") == ("3.162", "3.098")
    assert f4(nb["wr"][0]) == "0.9798" and f"{nb['wr'][1]:.1f}" == "0.6"
    assert abs(nb["wr"][0] - nb["wid"][1] / nb["wid"][0]) < 1e-12
    sp, m, var = nb["shift"]
    assert abs(m - 2) < 1e-12 and abs(var - 10) < 1e-9 and abs(sum(sp.values()) - 1) < 1e-12
    assert all(abs(x - round(x)) < 1e-9 and round(x) % 2 == 0 for x in sp)


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


def sample_steps(n=N0, seed=7):
    rng = np.random.default_rng(seed)
    return [1 if u < P_BIAS else -1 for u in rng.random(n)]


def axis_line(ux, origin, x0, x1, every=2, font_size=20):
    base = Line(origin + x0 * ux * RIGHT, origin + x1 * ux * RIGHT).set_stroke(GREY_C, 2)
    nums = VGroup(*[tx(str(x), font_size, GREY_B).move_to(origin + x * ux * RIGHT + 0.28 * DOWN)
                    for x in range(int(math.ceil(x0)), int(x1) + 1) if x % every == 0])
    return VGroup(base, nums)


def smooth_on(ax, fn, x0, x1, color, n=120, stroke=3):
    xs = np.linspace(x0, x1, n)
    return VMobject().set_points_smoothly([ax.c2p(x, fn(x)) for x in xs]).set_stroke(color, stroke)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nl = D.number_line(6, unit=0.8)
    at(nl, 0, -0.6)
    walker = D.walker_dot().move_to(D.line_point(nl, 0) + 0.3 * UP)
    a_r = Arrow(D.line_point(nl, 0) + 1.0 * UP, D.line_point(nl, 2) + 1.0 * UP, buff=0).set_color(D.HEAD_COLOR)
    a_l = Arrow(D.line_point(nl, 0) + 1.0 * UP, D.line_point(nl, -2) + 1.0 * UP, buff=0).set_color(D.TAIL_COLOR)
    l_r = tx(r"p=0.6", 32, D.HEAD_COLOR).next_to(a_r, UP, buff=0.12)
    l_l = tx(r"1-p=0.4", 32, D.TAIL_COLOR).next_to(a_l, UP, buff=0.12)
    coins = VGroup(D.coin("H"), D.coin("T")).arrange(RIGHT, buff=0.6)
    at(coins, 0, 2.4)
    ask = at(tx(r"\langle x\rangle=?\qquad\text{width}=?", 36, YELLOW), 0, -1.9)
    assert_in_frame(nl, a_r, a_l, l_r, l_l, coins, ask)
    assert_apart(l_r, l_l)
    assert_apart(coins, VGroup(l_r, l_l))
    s.play(FadeIn(nl), FadeIn(walker), FadeIn(coins), run_time=0.8)
    s.say("問：一枚不公正的硬幣", wait=1.2)
    s.play(GrowArrow(a_r), GrowArrow(a_l), FadeIn(l_r), FadeIn(l_l), run_time=0.8)
    s.say("正面向右，機率 ", "$p=0.6$", "；反面向左，機率 0.4", wait=2.0)
    x = 0
    for st in sample_steps():
        x += st
        s.play(walker.animate.move_to(D.line_point(nl, x) + 0.3 * UP), run_time=0.28)
    s.say("醉漢走 10 步，平均到哪裡？", wait=1.5)
    s.play(FadeIn(ask), run_time=0.6)
    s.say("分布比公正硬幣寬還是窄？（先猜）", wait=3.0)
    s.clear_stage()


def part_exact(s):
    nb = numbers()
    ux, uy, org = 0.3, 14.0, np.array([-3.5, -2.3, 0])
    ax = axis_line(ux, org, -10.6, 10.6)
    chart = D.bar_chart({x: float(v) for x, v in nb["P10"].items()}, ux=ux, uy=uy, origin=org, width=1.4)
    gauss = D.curve(lambda x: 2 * float(V.biased_gauss(x, N0, P_BIAS)), -10, 10, ux=ux, uy=uy, origin=org)
    e1 = eq_label(tx(r"P(N;x)=\binom{N}{\frac{N+x}{2}}\,p^{\frac{N+x}{2}}(1-p)^{\frac{N-x}{2}}", 30), "(5.1)")
    e2 = eq_label(tx(r"P(N+1;x)=p\,P(N;x-1)+(1-p)\,P(N;x+1)", 28), "(5.2)")
    m1 = mixed("平均 ", "$(2p-1)N=%d$" % nb["mom"][0], font_size=26)
    m2 = mixed("變異數 ", "$4Np(1-p)=%.1f$" % float(nb["mom"][1]), "（公正：", "$%d$" % nb["fair_var"], "）",
               font_size=26)
    for m in (e1, e2):
        fit(m, 6.0)
    at(e1, 3.75, 2.55)
    at(e2, 3.75, 1.65)
    at(m1, 3.75, 0.85)
    at(m2, 3.75, 0.25)
    rows = [["x", "P(10;x)", r"2\times\text{Gauss}"]]
    rows += [[str(x), f4(a), f4(g)] for x, a, g in nb["tri"]]
    tab = D.tex_table(rows, [1.0, 2.0, 2.3], row_h=0.5, font_size=28)
    at(tab, 3.75, -1.45)
    for i in range(1, 4):
        D.table_cell(tab, i, 2, 3).set_color(D.GAUSS_COLOR)
    left = VGroup(ax, chart, gauss)
    assert_in_frame(left, e1, e2, m1, m2, tab)
    for a, b in ((left, e1), (left, e2), (left, m2), (left, tab), (e1, e2), (e2, m1), (m1, m2), (m2, tab)):
        assert_apart(a, b)
    assert_table(tab, 3)
    s.play(FadeIn(ax), Write(e1), run_time=1.2)
    s.say("精確答案：二項分布 (5.1)", wait=2.0)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("也可以一步一步遞迴 (5.2)：兩種算法一模一樣", wait=2.5)
    s.play(LaggedStart(*[GrowFromEdge(r, DOWN) for r in chart[0]], lag_ratio=0.08), run_time=1.5)
    s.say("走 10 步後的位置：只有偶數格有人", wait=2.0)
    s.play(FadeIn(m1), run_time=0.6)
    s.say("平均往右偏了 2 格", wait=2.0)
    s.play(FadeIn(m2), run_time=0.6)
    s.say("變異數 9.6：比公正硬幣的 10 小一點", wait=2.5)
    s.play(*[D.bar_of(chart, x).animate.set_fill(YELLOW, 0.9) for x in MARK], FadeIn(tab[1]),
           *[FadeIn(D.table_cell(tab, i, j, 3)) for i in range(4) for j in range(2)], run_time=0.8)
    s.say("三根長條：", "$P(10;0),\\ P(10;2),\\ P(10;4)$", wait=2.5)
    s.play(ShowCreation(gauss), *[FadeIn(D.table_cell(tab, i, 2, 3)) for i in range(4)], run_time=1.2)
    s.say("疊上高斯 (5.7)：只有偶數格有人，所以 ", "$\\times2$", wait=3.0)
    s.clear_stage()


def part_fig51(s):
    nb = numbers()
    fig = nb["fig"]
    ax = D.plot_axes((-10, 30, 5), (0, 0.3, 0.1), width=7.2, height=4.2, x_ticks=(-10, 0, 10, 20, 30),
                     y_ticks=(0.1, 0.2, 0.3), x_tick_labels=["-10", "0", "10", "20", "30"],
                     y_tick_labels=["0.1", "0.2", "0.3"], x_label="x")
    ax.move_to(np.array([-2.65, -0.35, 0]))
    groups = {}
    for N in NS:
        col = N_COL[N]
        pts = [(x, float(v)) for x, v in fig[N]["P"].items() if -10 <= x <= 30]
        dots = VGroup(*[Dot(ax.c2p(x, y), radius=0.05).set_fill(col, 1) for x, y in pts])
        cur = smooth_on(ax, lambda x, N=N: 2 * float(V.biased_gauss(x, N, P_BIAS)), -10, 30, col)
        groups[N] = (dots, cur)
    hdr = lambda t: zh(t, font_size=22, color=YELLOW)
    rows = [["N", hdr("平均"), hdr("變異數"), hdr("峰值"), r"2\times\text{Gauss}"]]
    for N in NS:
        rows.append([str(N), str(fig[N]["mean"]), f"{fig[N]['var']:.1f}", f4(fig[N]["peak"]), f4(fig[N]["gpeak"])])
    tab = D.tex_table(rows, [0.7, 0.85, 1.05, 1.15, 1.35], row_h=0.55, font_size=24)
    at(tab, 4.2, 1.9)
    for i, N in enumerate(NS, start=1):
        D.table_cell(tab, i, 0, 5).set_color(N_COL[N])
    g57 = eq_label(tx(r"P\approx\frac{1}{\sqrt{8\pi Np(1-p)}}\exp\Big[-\frac{(x-(2p-1)N)^2}{8Np(1-p)}\Big]", 26),
                   "(5.7)", buff=0.2)
    fit(g57, 5.0)
    at(g57, 4.2, -0.55)
    wline = at(mixed("寬度 ", "$\\propto\\sqrt{p(1-p)}$", "　", "$p\\to1$", "：寬度 ", "$\\to0$", font_size=24), 4.2, -1.75)
    fit(wline, 5.0)
    assert_in_frame(ax, tab, g57, wline)
    for a, b in ((ax, tab), (ax, g57), (ax, wline), (tab, g57), (g57, wline)):
        assert_apart(a, b)
    assert_table(tab, 5)
    s.play(FadeIn(ax), FadeIn(tab[1]), *[FadeIn(c) for c in tab[0][:5]], run_time=1.0)
    s.say("Fig. 5.1：同一枚硬幣，走 10、20、40 步", wait=2.0)
    for i, N in enumerate(NS, start=1):
        dots, cur = groups[N]
        s.play(FadeIn(dots, lag_ratio=0.05), ShowCreation(cur), *[FadeIn(c) for c in tab[0][5 * i:5 * i + 5]],
               run_time=1.2)
        s.say("$N=%d$" % N, "：平均 %d、變異數 %.1f" % (fig[N]["mean"], fig[N]["var"]), wait=1.8)
    s.say("峰往右移、越來越寬，也越來越矮", wait=2.5)
    s.play(FadeIn(g57, shift=0.2 * UP), run_time=0.8)
    s.say("(5.7) 的高斯：寬度 ", "$\\sqrt{4Np(1-p)}$", wait=2.5)
    s.play(FadeIn(wline), run_time=0.6)
    s.say("$p\\to1$", "：每一步都向右，沒有隨機性了", wait=3.0)
    s.clear_stage()


def part_taylor(s):
    nb = numbers()
    e1 = eq_label(tx(r"P(N+1,y)=p\,P(N,y-l)+(1-p)\,P(N,y+l)", 30), "(5.2)")
    e2 = eq_label(tx(r"P+\partial_NP\approx P+l(1-2p)\,\partial_yP+\frac{l^2}{2}\,\partial_y^2P", 30), "(5.3)")
    e3 = eq_label(tx(r"\partial_NP=-v\,\partial_yP+\frac{l^2}{2}\,\partial_y^2P,\quad v=(2p-1)\,l", 30), "(5.4)")
    for m, y in ((e1, 2.6), (e2, 1.75), (e3, 0.85)):
        at(m, 0, y)
    box = at(D.erratum_box("p. 96，(5.3)", r"l(2p-1)\,\partial_y P", r"l(1-2p)\,\partial_y P", width=7.0), 0, -1.0)
    dl = at(VGroup(zh("多保留一項：", font_size=26, color=GREY_A),
                   tx(r"D=\frac{l^2}{2}-\frac{v^2}{2}=2p(1-p)", 30, YELLOW),
                   tx(r"\text{(Ex. 5.1)}", 22, GREY_B)).arrange(RIGHT, buff=0.3), 0, -0.05)
    hdr = lambda t: zh(t, font_size=24, color=YELLOW)
    rows = [[hdr("算法"), "D", hdr("變異數 2DN"), hdr("寬度")],
            [zh("直接寫擴散方程", font_size=24), "%.1f" % nb["Dn"], "%d" % round(nb["var_n"]), "%.3f" % nb["wid"][0]],
            [zh("多保留一項", font_size=24), "%.2f" % nb["Dc"], "%.1f" % nb["var_c"], "%.3f" % nb["wid"][1]]]
    tab = D.tex_table(rows, [3.0, 1.2, 2.2, 1.5], row_h=0.5, font_size=28)
    at(tab, 0, -1.3)
    for j in range(1, 4):
        D.table_cell(tab, 1, j, 4).set_color(D.WRONG_COLOR)
        D.table_cell(tab, 2, j, 4).set_color(D.RIGHT_COLOR)
    ratio = at(mixed("寬度比 ", "$\\sqrt{4p(1-p)}$", "：", "$p=0.6$", " 是 %s，" % f4(nb["wr"][0]),
                     "$p=0.9$", " 是 %.1f" % nb["wr"][1], font_size=26), 0, -2.55)
    assert_in_frame(e1, e2, e3, box, dl, tab, ratio)
    for a, b in ((e1, e2), (e2, e3), (e3, box), (e3, dl), (dl, tab), (tab, ratio)):
        assert_apart(a, b)
    assert_table(tab, 4)
    s.play(Write(e1), run_time=1.2)
    s.say("把格子換成連續的 ", "$y$", "：遞迴 (5.2) 的連續版", wait=2.0)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("泰勒展開到二次：一次項的係數是 ", "$l(1-2p)$", wait=2.5)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (5.3) 寫成 ", "$l(2p-1)$", "：號反了", wait=2.5)
    s.say("往右偏、峰往右移：", "$\\partial P/\\partial y$", " 前面要是負號", wait=3.0)
    s.play(FadeIn(e3, shift=0.2 * DOWN), run_time=0.8)
    s.say("(5.4) 已經寫對：峰以速度 ", "$v=(2p-1)l$", " 往右走", wait=2.5)
    s.play(FadeOut(box), FadeIn(tab[1]), *[FadeIn(c) for c in tab[0][:8]], run_time=0.8)
    s.say("Ex. 5.1 的陷阱：直接寫成擴散方程，", "$D=l^2/2=0.5$", wait=2.5)
    s.say("變異數 10、寬度 3.162：和公正硬幣一樣，不對", wait=2.5)
    s.play(FadeIn(dl), *[FadeIn(c) for c in tab[0][8:]], run_time=0.8)
    s.say("多保留一項：", "$D=2p(1-p)=0.48$", "，變異數 9.6", wait=3.0)
    s.play(FadeIn(ratio), run_time=0.6)
    s.say("寬度比 ", "$\\sqrt{4p(1-p)}$", "：", "$p=0.9$", " 時只剩 0.6", wait=2.5)
    s.say("擴散方程只在 ", "$p\\approx\\frac12$", " 時才對", wait=2.5)
    s.clear_stage()


def part_shift(s):
    nb = numbers()
    sp, m, var = nb["shift"]
    e1 = eq_label(tx(r"\chi(\vec q)=e^{-i\vec q\cdot\vec r_B}\chi_0(\vec q)", 32), "(5.9)")
    e2 = eq_label(tx(r"P(N;\vec r)=P_0(N;\vec r-N\vec r_B)", 32), "(5.10)")
    at(e1, -3.4, 2.55)
    at(e2, -3.4, 1.6)
    ux, uy, org = 0.28, 10.0, np.array([-3.75, -2.2, 0])
    ax = axis_line(ux, org, -10.6, 12.6, every=4)
    b_bias = D.bar_chart({x: float(v) for x, v in nb["P10"].items()}, ux=ux, uy=uy, origin=org + 0.45 * ux * RIGHT,
                         width=0.8)
    b_wind = D.bar_chart(sp, ux=ux, uy=uy, origin=org - 0.45 * ux * RIGHT, width=0.8, color=WIND_COLOR)
    hdr = lambda t: zh(t, font_size=24, color=YELLOW)
    rows = [[VectorizedPoint(), hdr("平均"), hdr("變異數")],
            [mixed("風吹 ", "$b=0.2$", "＋公正硬幣", font_size=22, color=WIND_COLOR), "%d" % round(m), "%d" % round(var)],
            [mixed("不公正硬幣 ", "$p=0.6$", font_size=22, color=D.BAR_COLOR), "%d" % nb["mom"][0],
             "%.1f" % float(nb["mom"][1])]]
    tab = D.tex_table(rows, [3.4, 1.1, 1.3], row_h=0.6, font_size=28)
    at(tab, 3.85, 0.6)
    ex1 = at(zh("例：在流動的水裡擴散的染料", font_size=26, color=GREY_A), 3.85, -1.0)
    ex2 = at(zh("整團跟著水流走，寬度照舊", font_size=26, color=GREY_A), 3.85, -1.6)
    left = VGroup(ax, b_bias, b_wind)
    assert_in_frame(e1, e2, left, tab, ex1, ex2)
    for a, b in ((e1, e2), (e2, left), (left, tab), (left, ex1), (tab, ex1), (ex1, ex2), (e1, tab)):
        assert_apart(a, b)
    assert_table(tab, 3)
    s.play(Write(e1), run_time=1.0)
    s.say("另一種偏向：每步先被吹走 ", "$\\vec r_B$", "，再照舊亂走", wait=2.5)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("轉換只多乘一個相位：整團平移 ", "$N\\vec r_B$", "，寬度不變", wait=3.0)
    s.play(FadeIn(ax), FadeIn(b_bias[0]), FadeIn(tab[1]), FadeIn(tab[0][:3]), FadeIn(tab[0][6:]), run_time=0.8)
    s.say("一維對照：不公正硬幣（藍）", wait=1.5)
    s.play(FadeIn(b_wind[0]), FadeIn(tab[0][3:6]), run_time=0.8)
    s.say("每步先被風吹 0.2 格，再擲公正硬幣（綠）", wait=2.5)
    s.say("平均都是 2；變異數卻是 10 對 9.6", wait=2.5)
    s.say("被風吹的只平移、不變窄", wait=2.0)
    s.play(FadeIn(ex1), FadeIn(ex2), run_time=0.8)
    s.say("例：在流動的水裡擴散的染料", wait=2.5)
    s.clear_stage()


def part_close(s):
    s.say("偏向讓高斯往前走，寬度乘上 ", "$\\sqrt{4p(1-p)}$", wait=3.0)
    s.say("「先走一步再亂走」則只平移、不變窄", wait=3.0)


class W501BiasedWalk(ChuScene):
    def construct(self):
        self.show_header("W501　有偏向的漫步：不公正的硬幣")
        check_math()
        part_question(self)
        part_exact(self)
        part_fig51(self)
        part_taylor(self)
        part_shift(self)
        part_close(self)
        self.wait(1.0)
