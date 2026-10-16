"""
W106 Stirling 公式（§1.3–1.3.1，pp. 10–13，式 (1.28)–(1.36)，Fig. 1.3，Exercise 1.4）

N! = Γ(N+1) = ∫ w^N e^{-w} dw：w^10 e^{-w} 的最高點在 w = 10、高 10^10 e^{-10} ≈ 453999、面積 3628800；
在 w = N 附近展開指數 (1.33) → 虛線高斯，面積 √(2πN) N^N e^{-N} = 3598696 (1.34)；筆誤框 p. 11 (1.34)；
ln 10! = 15.1044 vs 15.0961 (1.35)(1.36)；n!/Stirling(n) 表；練習 1.4：ln N! ≈ N ln N − N。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
import math

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

N = 10
F10 = math.factorial(N)
PEAK = W.gamma_integrand(N, N)
ST10 = W.stirling(N)
LN_F, LN_S = W.ln_factorial(N), W.ln_stirling(N)
REL = {n: (W.ln_factorial(n) - W.ln_stirling(n)) / W.ln_factorial(n) for n in (10, 100)}
RATIO_NS = (1, 2, 5, 10, 100)
RATIOS = [W.stirling_ratio(n) for n in RATIO_NS]
SCALE = 1e5                                     # 圖上 y 的單位：10^5


def riemann(f, a, b, h=0.01):
    n = int(round((b - a) / h))
    return h * sum(f(a + (i + 0.5) * h) for i in range(n))


def check_math():
    assert F10 == 3628800 == math.gamma(N + 1)
    assert round(PEAK) == 453999 and abs(PEAK - 1e10 * math.exp(-10)) < 1e-6
    ws = [0.01 * i for i in range(1, 3000)]
    assert abs(max(ws, key=lambda w: W.gamma_integrand(w, N)) - N) < 0.011              # 最高點 w = N
    assert abs(riemann(lambda w: W.gamma_integrand(w, N), 0, 80) / F10 - 1) < 1e-6         # 面積 = 10!
    assert abs(W.gamma_gauss(N, N) - PEAK) < 1e-6
    assert abs(riemann(lambda w: W.gamma_gauss(w, N), -60, 80) / ST10 - 1) < 1e-6          # 高斯面積 (1.34)
    assert round(ST10) == 3598696 and abs(ST10 - math.sqrt(2 * math.pi * N) * N ** N * math.exp(-N)) < 1e-6
    assert (round(LN_F, 4), round(LN_S, 4)) == W.BOOK_LN10 == (15.1044, 15.0961)           # (1.35)(1.36)
    assert f"{REL[10]:.1e}" == "5.5e-04" and f"{REL[100]:.1e}" == "2.3e-06"
    assert [round(r, 4) for r in RATIOS] == [1.0844, 1.0422, 1.0168, 1.0084, 1.0008]
    assert [round(1 + 1 / (12 * n), 4) for n in RATIO_NS] == [1.0833, 1.0417, 1.0167, 1.0083, 1.0008]
    assert all(abs(W.ln_factorial(n) - W.ln_factorial(n - 1) - math.log(n)) < 1e-9 for n in range(2, 60))
    assert abs(W.ln_factorial(10) - sum(math.log(k) for k in range(1, 11))) < 1e-9
    assert [math.factorial(k) for k in range(1, 11)] == [1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800]
    assert math.gamma(N + 1) == N * math.gamma(N) and math.gamma(1) == 1
    assert 0 < 1 - ST10 / F10 < 0.01                                                   # 差不到 1%
    assert round(100 * (RATIOS[0] - 1)) == 8                                           # n = 1 差 8%


def eq_row(tex, num, font_size=36):
    row = VGroup(Tex(tex, font_size=font_size))
    if num:
        row.add(Tex(num, font_size=24).set_color(GREY_B))
    return row.arrange(RIGHT, buff=0.35)


def at_left(m, y, x):
    return m.move_to(np.array([0, y, 0])).align_to(np.array([x, 0, 0]), LEFT)


# ---------------------------------------------------------------- 1. 問題
def ask(self):
    big = Tex(r"10!=1\cdot2\cdot3\cdots10=3628800", font_size=60).move_to(0.9 * UP)
    prods = Tex(r"1,\ 2,\ 6,\ 24,\ 120,\ 720,\ 5040,\ 40320,\ 362880,\ 3628800", font_size=36)
    prods.set_color(GREY_A).move_to(0.6 * DOWN)
    assert_in_frame(big, prods)
    assert_apart(big, prods)
    self.say("$10!=3628800$", "：能不能不乘，就把它估出來？")
    self.play(Write(big), run_time=1.8)
    self.wait(2.5)
    self.say("一路乘上去：", r"$1!,\ 2!,\ \dots,\ 10!$", " 長得飛快")
    self.play(Write(prods), run_time=3.0)
    self.wait(3.0)
    self.play(FadeOut(big), FadeOut(prods))


# ---------------------------------------------------------------- 2.–3. 積分與高斯
def gamma_plot(self):
    t1 = eq_row(r"N!=\Gamma(N+1)=\int_0^\infty w^N e^{-w}\,dw", r"(1.29)\text{--}(1.31)", 38).move_to(2.55 * UP)
    ax = D.plot_axes((0, 30, 5), (0, 5, 1), width=7.0, height=3.5, x_ticks=[5, 10, 15, 20, 25, 30],
                     y_ticks=[2, 4], y_tick_labels=[r"2\times10^5", r"4\times10^5"], x_label="w",
                     y_label=r"w^{10}e^{-w}")
    ax.move_to(np.array([-1.95, -0.75, 0]))
    f = lambda w: W.gamma_integrand(w, N) / SCALE
    g = lambda w: W.gamma_gauss(w, N) / SCALE
    graph = ax.get_graph(f, x_range=[0, 30, 0.1]).set_stroke(D.BAR_COLOR, 4)
    area = Polygon(*[ax.c2p(w, f(w)) for w in np.linspace(0, 30, 151)], ax.c2p(30, 0), ax.c2p(0, 0))
    area.set_fill(D.BAR_COLOR, 0.3).set_stroke(width=0)
    gauss = DashedVMobject(ax.get_graph(g, x_range=[0.01, 30, 0.1]), num_dashes=45).set_stroke(D.GAUSS_COLOR, 3)
    top = Dot(ax.c2p(N, f(N)), radius=0.07).set_fill(YELLOW, 1)
    drop = DashedLine(ax.c2p(N, f(N)), ax.c2p(N, 0), dash_length=0.08).set_stroke(YELLOW, 1.5)
    rx = 2.3
    p1 = at_left(mixed("最高點在 ", "$w=N=10$", font_size=26), 1.5, rx)
    p2 = at_left(Tex(r"10^{10}e^{-10}\approx 453999", font_size=34).set_color(YELLOW), 0.85, rx)
    p3 = at_left(mixed("藍色面積 ", "$=10!=3628800$", font_size=26, color=D.BAR_COLOR), 0.15, rx)
    assert_in_frame(t1, ax, p1, p2, p3)
    assert_apart(t1, ax)
    assert_apart(ax, VGroup(p1, p2, p3))
    self.say("階乘可以寫成一個積分（Gamma 函數）")
    self.play(Write(t1), run_time=1.8)
    self.wait(2.0)
    self.say("分部積分：", r"$\Gamma(N+1)=N\,\Gamma(N)$", "，", r"$\Gamma(1)=1$", "：正好是階乘", wait=3.5)
    self.say("畫出 ", "$w^{10}e^{-w}$", "：一座山（Fig. 1.3）")
    self.play(FadeIn(ax))
    self.play(ShowCreation(graph), run_time=2.0)
    self.say("$w^N$", " 往上衝、", "$e^{-w}$", " 往下掉：相乘成一座山", wait=3.0)
    self.play(FadeIn(top), ShowCreation(drop), FadeIn(p1))
    self.play(Write(p2))
    self.wait(3.0)
    self.say("山底下的面積，就是 ", "$10!$")
    self.play(FadeIn(area), FadeIn(p3))
    self.wait(3.5)
    return ax, graph, area, gauss, t1, VGroup(p1, p2, p3)


def gauss_part(self, ax, gauss, t1, panel):
    t2 = eq_row(r"N\ln w-w\approx N\ln N-N-\frac{(w-N)^2}{2N}", "(1.33)", 38).move_to(2.55 * UP)
    rx = 2.3
    q1 = at_left(mixed("橘色虛線：在 ", "$w=N$", " 附近的高斯", font_size=24, color=D.GAUSS_COLOR), 1.5, rx)
    q2 = at_left(eq_row(r"\sqrt{2\pi N}\,N^Ne^{-N}", "(1.34)", 36), 0.75, rx)
    q3 = at_left(Tex(r"=3598696", font_size=36).set_color(D.GAUSS_COLOR), 0.0, rx + 0.4)
    box = D.erratum_box("p. 11，(1.34)", r"\exp[N(\ln N-N)]", r"\exp[N\ln N-N]", width=4.4, font_size=28)
    box.move_to(np.array([4.45, -1.75, 0]))
    assert_in_frame(t2, q1, q2, q3, box)
    assert_apart(ax, VGroup(q1, q2, q3, box))
    assert_apart(q3, box)
    self.say("把指數在山頂 ", "$w=N$", " 附近展開到二次")
    self.play(FadeOut(t1, shift=0.15 * UP), FadeIn(t2, shift=0.15 * UP), run_time=1.2)
    self.wait(3.0)
    self.say("指數是二次式：山變成一條鐘形曲線")
    self.play(FadeOut(panel), ShowCreation(gauss), run_time=2.0)
    self.play(FadeIn(q1))
    self.wait(1.5)
    self.say("山頂一樣高，寬度也差不多：兩條曲線幾乎疊在一起")
    self.play(Flash(ax.c2p(N, PEAK / SCALE), color=YELLOW, flash_radius=0.35), run_time=1.2)
    self.wait(2.5)
    self.say("鐘形曲線的面積：Stirling 公式")
    self.play(Write(q2), run_time=1.5)
    self.play(Write(q3))
    self.wait(2.0)
    self.say("和 ", "$10!=3628800$", " 只差不到 1%", wait=3.5)
    self.say("書上 (1.34) 指數的括號放錯了")
    self.play(FadeIn(box, shift=0.2 * UP))
    self.wait(4.0)


# ---------------------------------------------------------------- 4. 比較與表
def compare(self):
    lx = -6.3
    c1 = at_left(eq_row(r"\ln 10!=15.1044", "(1.35)", 36), 2.3, lx)
    c2 = at_left(eq_row(r"\ln\big(\sqrt{20\pi}\,10^{10}e^{-10}\big)=15.0961", "(1.36)", 36), 1.3, lx)
    c3 = at_left(mixed("相對誤差 ", "$5.5\\times10^{-4}$", font_size=28), 0.2, lx)
    c4 = at_left(mixed("$N=100$", "：相對誤差 ", "$2.3\\times10^{-6}$", font_size=28), -0.6, lx)
    rows = [["n", r"n!/\text{Stirling}(n)", r"1+\frac{1}{12n}"]]
    for n, r in zip(RATIO_NS, RATIOS):
        rows.append([str(n), f"{r:.4f}", f"{1 + 1 / (12 * n):.4f}"])
    tab = D.tex_table(rows, [0.9, 2.5, 1.7], row_h=0.62, font_size=30)
    tab.move_to(np.array([3.95, 0.6, 0]))
    assert_in_frame(c1, c2, c3, c4, tab)
    assert_apart(VGroup(c1, c2, c3, c4), tab)
    self.say("比一比：取對數，誤差小得看不出來")
    self.play(Write(c1), run_time=1.2)
    self.play(Write(c2), run_time=1.6)
    self.wait(3.0)
    self.play(FadeIn(c3, shift=0.15 * UP))
    self.wait(2.5)
    self.say("$N$", " 越大越準：", "$N=100$", " 時只差百萬分之二")
    self.play(FadeIn(c4, shift=0.15 * UP))
    self.wait(3.5)
    self.say("不取對數的比值 ", "$n!/\\text{Stirling}(n)$", "：大約 ", "$1+1/(12n)$")
    self.play(FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:3])))
    for i in range(1, len(rows)):
        self.play(FadeIn(VGroup(*tab[0][3 * i:3 * i + 3]), shift=0.1 * UP), run_time=0.6)
    self.wait(3.0)
    self.say("連 ", "$n=1$", " 都只差 8%：近似得出奇地好", wait=4.0)


# ---------------------------------------------------------------- 5. 練習 1.4
def exercise(self):
    e1 = Tex(r"\frac{d}{dN}\ln N!\approx\ln N!-\ln(N-1)!=\ln N", font_size=36).move_to(np.array([0, 2.55, 0]))
    e2 = Tex(r"\ln N!\approx\int^N\ln x\,dx=N\ln N-N", font_size=36).move_to(np.array([0, 1.5, 0]))
    note = zh("（(1.28) 的第一項）", font_size=24, color=GREY_B).next_to(e2, RIGHT, buff=0.3)
    ax = D.plot_axes((0, 10.6, 1), (0, 2.6, 1), width=6.0, height=2.9, x_ticks=[1, 5, 10], y_ticks=[1, 2],
                     x_label="x")
    ax.move_to(np.array([-2.3, -1.25, 0]))
    steps = VGroup(*[Polygon(ax.c2p(k - 1, 0), ax.c2p(k, 0), ax.c2p(k, math.log(k)), ax.c2p(k - 1, math.log(k)))
                     .set_fill(D.BAR_COLOR, 0.45).set_stroke(D.BAR_COLOR, 1.5) for k in range(2, 11)])
    lnx = ax.get_graph(lambda x: math.log(x), x_range=[1, 10.5, 0.05]).set_stroke(D.GAUSS_COLOR, 3)
    s1 = mixed("每階：寬 1、高 ", "$\\ln k$", font_size=24, color=D.BAR_COLOR)
    s2 = mixed("面積和 ", "$=\\ln 10!$", font_size=24, color=D.BAR_COLOR)
    s3 = mixed("橘線 ", "$\\ln x$", " 下的面積", font_size=24, color=D.GAUSS_COLOR)
    side = VGroup(s1, s2, s3).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    side.move_to(np.array([0, -1.0, 0])).align_to(np.array([1.9, 0, 0]), LEFT)
    assert_in_frame(e1, e2, note, ax, steps, lnx, side)
    assert_apart(e2, ax)
    assert_apart(VGroup(ax, steps), side)
    self.say("練習 1.4：", "$\\ln N!$", " 每加一步，就多一個 ", "$\\ln N$")
    self.play(Write(e1), run_time=1.6)
    self.wait(2.0)
    self.say("$N$", " 很大時，差一步的差分就像導數", wait=3.0)
    self.say("台階的面積和是 ", "$\\ln N!$", "，差不多就是曲線下的面積")
    self.play(FadeIn(ax))
    self.play(LaggedStart(*[GrowFromEdge(s, DOWN) for s in steps], lag_ratio=0.1), FadeIn(s1), FadeIn(s2),
              run_time=1.8)
    self.play(ShowCreation(lnx), FadeIn(s3), run_time=1.4)
    self.wait(3.0)
    self.say("積分得到 Stirling 的第一項 ", "$N\\ln N-N$")
    self.play(Write(e2), FadeIn(note), run_time=1.6)
    self.wait(4.5)


class W106Stirling(ChuScene):
    def construct(self):
        self.show_header("W106　Stirling 公式")
        check_math()
        ask(self)
        ax, graph, area, gauss, t1, panel = gamma_plot(self)
        gauss_part(self, ax, gauss, t1, panel)
        self.clear_stage()
        compare(self)
        self.clear_stage()
        exercise(self)
        self.say("階乘＝一個積分；積分", "$\\approx$", "一條鐘形曲線的面積", wait=4.0)
        self.wait(1.0)
