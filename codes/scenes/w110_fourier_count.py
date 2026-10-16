"""
W110 傅立葉：用積分數走法（§1.5，pp. 15–16，(1.49)–(1.57)）

問：(2cos q)^10 對 q 積分，怎麼會數出走法？
  1. χ(q) = e^{iql} + e^{−iql} = 2cos ql (1.50)；N = 3：e^{3iq} + 3e^{iq} + 3e^{−iq} + e^{−3iq}，係數 1,3,3,1 = C(3,d)。
  2. 取係數：(1/2π)∫ e^{−iqd}e^{iqk} dq = 1（k = d）或 0（畫 cos 2q：正負面積抵消）→ (1.51)；N = 10：d = 2 → 210，d = 3 → 0。
  3. 畫 (2cos q)^10/1024：q = 0 有高峰，q = ±π 也有。只取 q = 0（(1.52)–(1.54)）得 105.8——精確 210，只拿到一半；
     q = ±π 兩個半峰補上另一半（(−1)^10 e^{∓2πi} = 1）：合計 211.5；d = 3 兩端相位相反 → 0。
     註記框：書上「上下限換成 ±∞ 的誤差可以忽略」→ 對離散的計數是因子 2。
  4. 再對 N 加總：g(z,q) = Σ χ^N z^N = 1/(1 − zχ(q)) (1.55)；小 q：1/(1 − 2z + zq²l²) (1.57)。
畫面上的數字都由 utils.walk1d 算出，check_math() 核對。
"""
import cmath
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import boxed, mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

SADDLE2 = TEAL                                  # q = ±π 的鞍點
TERMS = [(3, r"e^{3iq}"), (1, r"3e^{iq}"), (-1, r"3e^{-iq}"), (-3, r"e^{-3iq}")]


def numbers():
    q0, both = W.saddle_q0(10, 2), W.saddle_both(10, 2)
    return {"c2": f"{W.fourier_count(10, 2):.3f}", "c3": f"{abs(W.fourier_count(10, 3)):.3f}",
            "q0": f"{q0:.1f}", "q0x": f"{q0:.2f}", "both": f"{both:.1f}", "bothx": f"{both:.2f}",
            "both3": f"{abs(W.saddle_both(10, 3)):.1f}"}


def check_math():
    assert W.chi_power_coeffs(3) == {-3: 1, -1: 3, 1: 3, 3: 1} == {d: W.walk_count(3, d) for d in (-3, -1, 1, 3)}
    assert [int(t[1][0]) if t[1][0].isdigit() else 1 for t in TERMS] == [W.walk_count(3, d) for d, _ in TERMS]
    assert numbers() == {"c2": "210.000", "c3": "0.000", "q0": "105.8", "q0x": "105.77", "both": "211.5",
                         "bothx": "211.53", "both3": "0.0"}
    assert W.walk_count(10, 2) == 210 and W.walk_count(10, 3) == 0
    for m in range(-4, 5):                     # 正交：(1/2π)∫ e^{iqm} dq = [m = 0]
        M = 16
        val = sum(cmath.exp(1j * m * (-math.pi + 2 * math.pi * k / M)) for k in range(M)) / M
        assert abs(val - (1 if m == 0 else 0)) < 1e-12
    assert abs((-1) ** 10 * cmath.exp(-2j * math.pi) - 1) < 1e-12        # d = 2：兩端相位 +1
    assert abs((-1) ** 10 * cmath.exp(-3j * math.pi) + 1) < 1e-12        # d = 3：兩端相位 −1
    assert abs(W.chi(0) ** 10 / 1024 - 1) < 1e-12 and abs(W.chi(math.pi) ** 10 / 1024 - 1) < 1e-12
    q = 0.3                                    # (1.52)(1.53)：χ^N ≈ 2^N e^{−Nq²/2}（小 q）
    assert abs(W.chi(q) ** 10 / (1024 * math.exp(-10 * q * q / 2)) - 1) < 0.01


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


def lobe(ax, fn, a, b, color, n=30):
    xs = np.linspace(a, b, n)
    pts = [ax.c2p(a, 0)] + [ax.c2p(x, fn(x)) for x in xs] + [ax.c2p(b, 0)]
    return Polygon(*pts).set_fill(color, 0.45).set_stroke(width=0)


def dashed_graph(ax, fn, a, b, color):
    g = ax.get_graph(fn, x_range=[a, b, 0.02])
    return DashedVMobject(g, num_dashes=24).set_stroke(color, 3)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：", "$(2\\cos q)^{10}$", " 對 ", "$q$", " 積分，怎麼會數出走法？", wait=4.0)


def part_expand(s):
    a1 = eq(r"\chi(q)=e^{iql}+e^{-iql}=2\cos ql", "(1.50)", font_size=40)
    a1n = mixed("每個因子＝一步：向右 ", "$e^{iql}$", "、向左 ", "$e^{-iql}$", font_size=26, color=GREY_B)
    a2 = mixed("$l=1,\\ N=3$", "：", "$\\chi(q)^3=(e^{iq}+e^{-iq})^3$", font_size=30)
    terms = [Tex(t, font_size=44) for _, t in TERMS]
    parts = [Tex("=", font_size=44)]
    for k, m in enumerate(terms):
        if k:
            parts.append(Tex("+", font_size=44))
        parts.append(m)
    row = VGroup(*parts).arrange(RIGHT, buff=0.25)
    dl = VGroup(*[Tex(f"d={d}", font_size=26).set_color(GREY_B).next_to(m, DOWN, buff=0.3)
                  for m, (d, _) in zip(terms, TERMS)])
    a3 = mixed("係數 1, 3, 3, 1 ＝ ", "$C(3,d)$", "：走 3 步到 ", "$d$", " 的走法數", font_size=28)
    col = VGroup(a1, a1n, a2, VGroup(row, dl), a3).arrange(DOWN, buff=0.45)
    col.move_to(np.array([0, 0.1, 0]))
    assert_in_frame(col)
    s.play(Write(a1[0]), FadeIn(a1[1]), run_time=1.5)
    s.say("每一步的「特徵函數」：向右 ", "$e^{iql}$", "、向左 ", "$e^{-iql}$", wait=0.3)
    s.play(FadeIn(a1n), run_time=0.6)
    s.wait(3.0)
    s.play(FadeIn(a2), run_time=0.8)
    s.play(FadeIn(row, shift=0.2 * UP), run_time=1.2)
    s.say("乘開：每一項 ", "$e^{iqd}$", " 對應一個終點 ", "$d$", wait=0.3)
    s.play(FadeIn(dl), run_time=0.8)
    s.wait(2.5)
    s.play(*[m.animate.set_color(YELLOW) for m in terms], FadeIn(a3), run_time=1.0)
    s.say("係數 1, 3, 3, 1 正是 ", "$C(3,d)$", "：數走法＝取係數", wait=4.5)
    s.clear_stage()


def part_extract(s):
    nb = numbers()
    b1 = Tex(r"\frac{1}{2\pi}\int_{-\pi}^{\pi}e^{-iqd}\,e^{iqk}\,dq=\begin{cases}1, & k=d\\ 0, & k\neq d\end{cases}",
             font_size=32)
    place_col(b1, -6.4, 2.9)
    ax = D.plot_axes((-3.5, 3.5, 1), (-1.3, 1.3, 1), width=4.4, height=2.0, x_ticks=[-PI, PI],
                     x_tick_labels=[r"-\pi", r"\pi"], font_size=20, x_label="q")
    ax.move_to(np.array([4.2, 1.85, 0]))
    fn = lambda q: math.cos(2 * q)
    graph = ax.get_graph(fn, x_range=[-PI, PI, 0.02]).set_stroke(D.BAR_COLOR, 3)
    cuts = [-PI, -3 * PI / 4, -PI / 4, PI / 4, 3 * PI / 4, PI]
    lobes = VGroup(*[lobe(ax, fn, a, b, GREEN if k % 2 == 0 else RED) for k, (a, b) in enumerate(zip(cuts, cuts[1:]))])
    glab = mixed("$e^{iq(k-d)}$", " 的實部，", "$k-d=2$", font_size=20)
    glab.next_to(ax, DOWN, buff=0.15)
    b2 = boxed(eq(r"C(N,d)=\frac{l}{2\pi}\int_{-\pi/l}^{\pi/l}e^{-iqd}\,\chi(q)^N\,dq", "(1.51)", font_size=38),
               color=D.GAUSS_COLOR, buff=0.2)
    b2.move_to(np.array([0, -0.85, 0]))
    b3 = mixed("數值積分 ", "$N=10$", "：", f"$d=2\\ \\to\\ {nb['c2']}$", "，", f"$d=3\\ \\to\\ {nb['c3']}$",
               font_size=28)
    b3.move_to(np.array([0, -2.45, 0]))
    assert_in_frame(b1, ax, glab, b2, b3)
    assert_apart(b1, ax)
    assert_apart(glab, b2)
    s.say("怎麼取出 ", "$e^{iqd}$", " 的係數？乘上 ", "$e^{-iqd}$", " 再積分", wait=0.3)
    s.play(Write(b1), run_time=1.8)
    s.wait(2.0)
    s.play(FadeIn(ax), ShowCreation(graph), FadeIn(glab), run_time=1.2)
    s.play(FadeIn(lobes), run_time=0.8)
    s.say("$k\\neq d$", "：繞整數圈，正負面積抵消，積分是 0", wait=4.0)
    s.say("$k=d$", "：被積函數是 1，除以 ", "$2\\pi$", " 正好是 1", wait=3.0)
    s.play(FadeIn(b2, shift=0.2 * UP), run_time=1.0)
    s.say("於是：走法數＝一個積分 (1.51)", wait=3.0)
    s.play(FadeIn(b3), run_time=0.8)
    s.say("數值積分：", "$d=2$", " 得 210，", "$d=3$", " 得 0（到不了）", wait=3.5)
    s.clear_stage()


def part_saddles(s):
    nb = numbers()
    ax = D.plot_axes((-3.4, 3.4, 1), (0, 1.15, 0.5), width=6.2, height=2.9, x_ticks=[-PI, 0, PI],
                     x_tick_labels=[r"-\pi", "0", r"\pi"], y_ticks=[1], font_size=22, x_label="q")
    ax.move_to(np.array([-3.25, 0.95, 0]))
    graph = ax.get_graph(lambda q: math.cos(q) ** 10, x_range=[-PI, PI, 0.01]).set_stroke(D.BAR_COLOR, 3)
    ylab = Tex(r"(2\cos q)^{10}/2^{10}", font_size=26).set_color(D.BAR_COLOR)
    ylab.move_to(ax.c2p(-1.6, 1.05))
    g0 = dashed_graph(ax, lambda q: math.exp(-5 * q * q), -1.3, 1.3, D.GAUSS_COLOR)
    gp = VGroup(dashed_graph(ax, lambda q: math.exp(-5 * (q - PI) ** 2), PI - 1.3, PI, SADDLE2),
                dashed_graph(ax, lambda q: math.exp(-5 * (q + PI) ** 2), -PI, -PI + 1.3, SADDLE2))
    t1 = eq(r"\chi^N\approx2^Ne^{-Nq^2/2}", "(1.52)(1.53)", font_size=30)
    t2 = eq(r"C(N,d)\approx\frac{2^N}{\sqrt{2\pi N}}e^{-d^2/2N}", "(1.54)", font_size=30)
    t3 = mixed("$N=10,\\ d=2$", "：", f"${nb['q0x']}$", font_size=28)
    t4 = mixed("精確是 ", "$210$", "：只拿到一半！", font_size=28)
    t5 = mixed("$q=\\pm\\pi$", "：", "$(-1)^{10}e^{\\mp2\\pi i}=1$", font_size=28, color=SADDLE2)
    t6 = Tex(rf"2\times{nb['q0x']}={nb['bothx']}", font_size=32).set_color(D.EXACT_COLOR)
    t7 = mixed("$d=3$", "：", "$(-1)^{10}e^{\\mp3\\pi i}=-1$", "，抵消成 ", f"${nb['both3']}$", font_size=28)
    col = VGroup(t1, t2, t3, t4, t5, t6, t7).arrange(DOWN, aligned_edge=LEFT, buff=0.28)
    for m in col:
        fit(m, 5.6)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.28), 1.0, 2.95)
    n1 = mixed("註：書上說「上下限換成 ", "$\\pm\\infty$", "，誤差可以忽略」", font_size=22)
    n2 = mixed("對離散的計數，漏掉的是 ", "$q=\\pm\\pi$", " 那一半：因子 2", font_size=22)
    nb_txt = VGroup(n1, n2).arrange(DOWN, aligned_edge=LEFT, buff=0.15)
    note = VGroup(RoundedRectangle(width=nb_txt.get_width() + 0.4, height=nb_txt.get_height() + 0.35,
                                   corner_radius=0.12).set_stroke(YELLOW, 2).set_fill(BLACK, 0.6), nb_txt)
    note[0].move_to(nb_txt)
    place_col(note, -6.5, -1.35)
    assert_in_frame(ax, ylab, col, note)
    assert_apart(ax, col)
    assert_apart(note, col)
    assert_apart(note, ax)

    s.say("畫出被積函數 ", "$(2\\cos q)^{10}/1024$", wait=0.3)
    s.play(FadeIn(ax), ShowCreation(graph), FadeIn(ylab), run_time=1.5)
    s.say("$q=0$", " 有一座高峰……", "$q=\\pm\\pi$", " 也有！", wait=3.0)
    s.say("先只看 ", "$q=0$", "：峰附近 ", "$\\chi^N\\approx2^Ne^{-Nq^2/2}$", wait=0.3)
    s.play(ShowCreation(g0), FadeIn(t1), run_time=1.2)
    s.wait(2.0)
    s.say("代入 (1.51)、上下限放寬到 ", "$\\pm\\infty$", "：高斯積分", wait=0.3)
    s.play(FadeIn(t2), run_time=0.8)
    s.wait(2.5)
    s.say("先猜：", "$N=10,\\ d=2$", " 會得到多少？", wait=2.5)
    s.play(FadeIn(t3), run_time=0.8)
    s.say("得到 ", f"${nb['q0']}$", "……可是精確答案是 210", wait=2.5)
    s.play(FadeIn(t4), run_time=0.8)
    s.say("只拿到一半！另一半在哪裡？", wait=2.5)
    s.play(ShowCreation(gp), FadeIn(t5), run_time=1.2)
    s.say("在 ", "$q=\\pm\\pi$", "：", "$\\cos q=-1$", "，十次方又變回正的", wait=3.5)
    s.play(FadeIn(t6), run_time=0.8)
    s.say("兩個半峰合起來又是一整座：合計 ", f"${nb['both']}$", wait=3.5)
    s.play(FadeIn(t7), run_time=0.8)
    s.say("$d=3$", "：兩端的相位相反，和中央抵消成 0（奇偶性）", wait=3.5)
    s.play(FadeIn(note, shift=0.2 * UP), run_time=0.8)
    s.say("書上說放寬上下限誤差可忽略：對離散計數其實差因子 2", wait=4.5)
    s.clear_stage()


def part_generating(s):
    g1 = eq(r"g(z,q)=\sum_{N\ge0}\chi(q)^Nz^N=\frac{1}{1-z\chi(q)}", "(1.55)", font_size=38)
    g2 = mixed("小 ", "$q$", "：", "$\\chi(q)=2\\cos ql\\approx2-q^2l^2$", font_size=28)
    g3 = boxed(eq(r"g(z,q)\approx\frac{1}{1-2z+zq^2l^2}", "(1.57)", font_size=38), color=D.GAUSS_COLOR, buff=0.2)
    g4 = mixed("走法數＝", "$g$", " 裡 ", "$z^N$", " 的係數、再取 ", "$e^{iqd}$", " 的係數", font_size=26, color=GREY_B)
    g5 = zh("下一幕：這個函數的奇點", font_size=26, color=GREY_B)
    col = VGroup(g1, g2, g3, g4, g5).arrange(DOWN, buff=0.5).move_to(np.array([0, 0.1, 0]))
    assert_in_frame(col)
    s.say("再把所有 ", "$N$", " 加起來：一個生成函數裝下全部", wait=0.3)
    s.play(Write(g1[0]), FadeIn(g1[1]), run_time=1.8)
    s.wait(2.5)
    s.play(FadeIn(g2), run_time=0.8)
    s.say("等比級數；小 ", "$q$", " 時 ", "$\\chi\\approx2-q^2l^2$", wait=3.0)
    s.play(FadeIn(g3, shift=0.2 * UP), run_time=1.0)
    s.say("得到 (1.57)：下一幕從它的奇點讀出走法", wait=3.0)
    s.play(FadeIn(g4), FadeIn(g5), run_time=0.8)
    s.wait(3.0)


def part_close(s):
    s.say("傅立葉把「數走法」變成「取係數」", wait=3.5)
    s.say("而且 ", "$q=0$", " 和 ", "$q=\\pm\\pi$", " 兩個鞍點各給一半", wait=3.5)


class W110FourierCount(ChuScene):
    def construct(self):
        self.show_header("W110　傅立葉：用積分數走法")
        check_math()
        part_question(self)
        part_expand(self)
        part_extract(self)
        part_saddles(self)
        part_generating(self)
        part_close(self)
        self.wait(2.0)
