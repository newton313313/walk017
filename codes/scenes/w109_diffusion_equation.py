"""
W109 從遞迴到擴散方程（§1.4，pp. 13–14，(1.44)–(1.48)，Ex. 1.5）

問：格子變細、步數變多，這個遞迴會變成什麼方程？
  1. 二階差分 (1.44)：f(y−l) + f(y+l) − 2f(y) ≈ l² f''(y)；f(y) = y²（取 y = 2）：1 + 9 − 2·4 = 2 = f''，剛好相等。
     畫面：曲線上三點、弦的中點（比曲線高 ½ × 二階差分）。
  2. (1.45) C(N+1) − 2C(N) ≈ l²∂²C → (1.46) P(N+1) − P(N) ≈ (l²/2)∂²P → (1.47) ∂P/∂N = (l²/2)∂²P/∂y²（D = l²/2）。
  3. 動畫：機率全在 y = 0，diffusion_step 走 30 步；放大縱軸、疊上 2·e^{−y²/2N}/√(2πN)（因子 2 同 W107）；
     N = 30：Σ y²P = 30；P(30,0) = 0.1445 對 0.1457；筆誤框 p. 14 (1.48) 指數少了負號，α = 1/√(2πl²)。
  4. 練習 1.5：P = e^{−y²/2N}/√(2πN) 代入，兩邊都等於 (y²/2N² − 1/2N) P。
畫面上的數字都由 utils.walk1d 算出，check_math() 核對。
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
from utils import walk1d as W
from utils import walk_draw as D

NSTEPS = 30
YMAX = 24                                  # 長條畫 y = −24..24（更遠的機率 < 10⁻⁵）
O5 = np.array([-1.7, -2.4, 0])
UXD = 0.19
UY_A, UY_B = 5.0, 16.0                     # 放大前／放大後的縱軸


def history():
    Ps = [{0: Fraction(1)}]
    for _ in range(NSTEPS):
        Ps.append(W.diffusion_step(Ps[-1]))
    return Ps


def numbers():
    P30 = history()[NSTEPS]
    return {"var": str(sum(y * y * p for y, p in P30.items())),
            "P0": f"{float(P30[0]):.4f}", "g0": f"{2 * W.heat_kernel(0, NSTEPS):.4f}"}


def check_math():
    sq = lambda y: y * y
    assert all(W.second_difference(sq, y) == 2 for y in range(-5, 6))
    assert [sq(1), sq(2), sq(3)] == [1, 4, 9] and (sq(1) + sq(3)) / 2 == 5 and 5 - sq(2) == 1
    Ps = history()
    for N, P in enumerate(Ps):
        assert sum(P.values()) == 1
        assert all(P.get(y, 0) == W.walk_prob(N, y) for y in range(-N - 1, N + 2))
        assert sum(y * y * p for y, p in P.items()) == N              # 變異數 = N（l = 1）
    assert numbers() == {"var": "30", "P0": "0.1445", "g0": "0.1457"}
    assert float(sum(p for y, p in Ps[NSTEPS].items() if abs(y) > YMAX)) < 1e-5
    # 練習 1.5：數值微分核對 ∂P/∂N = ½ ∂²P/∂y² = (y²/2N² − 1/2N) P（N = 30，y = 5）
    N, y, h = 30.0, 5.0, 1e-3
    P = W.heat_kernel(y, N)
    rhs = (y * y / (2 * N * N) - 1 / (2 * N)) * P
    dN = (W.heat_kernel(y, N + h) - W.heat_kernel(y, N - h)) / (2 * h)
    dyy = (W.heat_kernel(y + h, N) + W.heat_kernel(y - h, N) - 2 * P) / (h * h)
    assert abs(dN - rhs) < 1e-8 and abs(dyy / 2 - rhs) < 1e-6


# ---------------------------------------------------------------- 小工具
def eq(tex, tag=None, font_size=30, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return m
    return VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


def dist_bars(P, uy):
    vals = {y: float(P.get(y, 0)) for y in range(-YMAX, YMAX + 1)}
    bars = D.bar_chart(vals, ux=UXD, uy=uy, origin=O5, width=0.9)[0]
    for r in bars:
        if r.get_height() < 0.004:          # 這一步沒人（奇偶性）或太小：透明，保留位置讓 Transform 一對一
            r.set_fill(opacity=0).set_stroke(opacity=0)
    return bars


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：格子變細、步數變多，這個遞迴會變成什麼方程？", wait=3.0)


def part_second_difference(s):
    ax = D.plot_axes((0, 3.6, 1), (0, 10, 2), width=4.2, height=3.8, x_ticks=[1, 2, 3],
                     x_tick_labels=["y-l", "y", "y+l"], font_size=22, y_label="f")
    ax.move_to(np.array([-4.3, 0.3, 0]))
    graph = ax.get_graph(lambda t: t * t, x_range=[0, 3.15, 0.05]).set_stroke(D.BAR_COLOR, 3)
    pts = [ax.c2p(t, t * t) for t in (1, 2, 3)]
    dots = VGroup(*[Dot(p, radius=0.07) for p in pts])
    mid = ax.c2p(2, 5)
    chord = Line(pts[0], pts[2]).set_stroke(TEAL, 3)
    mdot = Dot(mid, radius=0.07).set_fill(TEAL, 1)
    gap = Line(pts[1], mid).set_stroke(YELLOW, 5)
    vals = VGroup(Tex("1", font_size=28).next_to(pts[0], DR, buff=0.08),
                  Tex("4", font_size=28).next_to(pts[1], RIGHT, buff=0.12),
                  Tex("9", font_size=28).next_to(pts[2], LEFT, buff=0.12))
    mlab = Tex("5", font_size=28).set_color(TEAL).next_to(mid, LEFT, buff=0.12)

    r1 = eq(r"f(y-l)+f(y+l)-2f(y)\approx l^2f''(y)", "(1.44)", font_size=32)
    r2 = mixed("弦的中點比曲線高出 ", "$\\tfrac12\\times$", " 二階差分", font_size=24)
    r3 = mixed("例：", "$f(y)=y^2$", "，取 ", "$y=2$", "：", "$f=1,\\ 4,\\ 9$", font_size=26)
    r4 = mixed("中點 ", "$\\tfrac12(1+9)=5$", "，比曲線高 ", "$5-4=1$", font_size=26)
    r5 = Tex(r"1+9-2\cdot4=2=f''(y)", font_size=34).set_color(YELLOW)
    r6 = Tex(r"(y-1)^2+(y+1)^2-2y^2=2", font_size=30)
    r7 = zh("二次函數：剛好相等，不只是近似", font_size=24, color=GREY_B)
    col = place_col(VGroup(r1, r2, r3, r4, r5, r6, r7).arrange(DOWN, aligned_edge=LEFT, buff=0.3), -1.3, 2.9)
    assert_in_frame(ax, graph, col, vals, mlab)
    assert_apart(ax, col)

    s.say("先看一個小工具：二階差分", wait=0.3)
    s.play(FadeIn(ax), ShowCreation(graph), run_time=1.2)
    s.play(Write(r1[0]), FadeIn(r1[1]), run_time=1.5)
    s.say("鄰居兩格相加、減掉自己兩倍 ", "$\\approx l^2f''$", wait=3.0)
    s.play(FadeIn(dots), FadeIn(vals), FadeIn(r3), run_time=1.0)
    s.say("例子 ", "$f(y)=y^2$", "：在 ", "$y-l,\\ y,\\ y+l$", " 取三點", wait=2.5)
    s.play(ShowCreation(chord), FadeIn(mdot), FadeIn(mlab), run_time=1.0)
    s.play(ShowCreation(gap), FadeIn(r2), FadeIn(r4), run_time=1.0)
    s.say("弦的中點比曲線高：差距是二階差分的一半", wait=3.5)
    s.play(FadeIn(r5, shift=0.2 * UP), run_time=0.8)
    s.say("二階差分 ", "$=2=f''$", "：曲線彎得越厲害，差距越大", wait=3.5)
    s.play(FadeIn(r6), FadeIn(r7), run_time=0.8)
    s.wait(2.5)
    s.clear_stage()


def part_derivation(s):
    c1 = eq(r"C(N+1,y)=C(N,y-l)+C(N,y+l)", "(1.42)")
    c2 = eq(r"C(N+1,y)-2C(N,y)=C(N,y-l)+C(N,y+l)-2C(N,y)\approx l^2\,\partial_y^2C", "(1.45)")
    c3 = eq(r"P(N+1,y)-P(N,y)\approx\frac{l^2}{2}\,\partial_y^2P", "(1.46)")
    c3n = mixed("（兩邊除以 ", "$2^{N+1}$", "，", "$C/2^N=P$", "）", font_size=22, color=GREY_B)
    c4 = boxed(eq(r"\frac{\partial P}{\partial N}=\frac{l^2}{2}\,\frac{\partial^2P}{\partial y^2}", "(1.47)",
                  font_size=40), color=D.GAUSS_COLOR, buff=0.2)
    c4n = mixed("擴散方程：每步的擴散係數 ", "$D=l^2/2$", font_size=26, color=D.GAUSS_COLOR)
    col = place_col(VGroup(c1, c2, c3, c3n, c4, c4n).arrange(DOWN, aligned_edge=LEFT, buff=0.38), -6.2, 2.9)
    assert_in_frame(col)
    s.say("把二階差分用在遞迴上", wait=0.3)
    s.play(FadeIn(c1), run_time=0.8)
    s.wait(1.5)
    s.say("兩邊減掉 ", "$2C(N,y)$", "：右邊正好是二階差分", wait=0.3)
    s.play(Write(c2[0]), FadeIn(c2[1]), run_time=2.0)
    s.wait(3.0)
    s.say("除以 ", "$2^{N+1}$", "，走法數變成機率", wait=0.3)
    s.play(Write(c3[0]), FadeIn(c3[1]), FadeIn(c3n), run_time=1.5)
    s.wait(3.0)
    s.say("步數很多時，相鄰兩步的差 ", "$\\approx$", " 對 ", "$N$", " 的導數", wait=0.3)
    s.play(FadeIn(c4, shift=0.2 * UP), run_time=1.0)
    s.play(FadeIn(c4n), run_time=0.6)
    s.say("這就是擴散方程 (1.47)：每步擴散係數 ", "$D=l^2/2$", wait=4.0)
    s.clear_stage()


def part_animation(s):
    Ps = history()
    nb = numbers()
    axis = Line(O5 + (YMAX + 1) * UXD * LEFT, O5 + (YMAX + 1) * UXD * RIGHT).set_stroke(D.MUTED, 2)
    ticks = VGroup(*[Tex(str(t), font_size=22).set_color(D.MUTED).next_to(O5 + t * UXD * RIGHT, DOWN, buff=0.12)
                     for t in range(-20, 21, 10)])
    ticks.add(Tex("y", font_size=30).set_color(D.MUTED).next_to(axis, RIGHT, buff=0.15))
    bars = dist_bars(Ps[0], UY_A)
    num = Integer(0, font_size=40)
    cnt = VGroup(Tex("N=", font_size=40), num).arrange(RIGHT, buff=0.12).move_to(np.array([4.6, 2.6, 0]))
    curve = D.curve(lambda y: 2 * W.heat_kernel(y, NSTEPS), -YMAX, YMAX, ux=UXD, uy=UY_B, origin=O5)
    n1 = Tex(r"\sum_y y^2P(30,y)=" + nb["var"], font_size=30)
    n2 = Tex(r"P(30,0)=" + nb["P0"], font_size=30)
    n3 = Tex(r"2\cdot\tfrac{1}{\sqrt{2\pi\cdot30}}=" + nb["g0"], font_size=30).set_color(D.GAUSS_COLOR)
    nums = place_col(VGroup(n1, n2, n3).arrange(DOWN, aligned_edge=LEFT, buff=0.3), 3.5, 1.9)
    box = D.erratum_box("p. 14，(1.48)", r"\exp\Big(\frac{(x-y)^2}{2l^2N}\Big)",
                        r"\exp\Big(-\frac{(x-y)^2}{2l^2N}\Big)", width=5.6, font_size=26)
    sol = Tex(r"\frac{\alpha}{\sqrt N}\exp\Big(-\frac{(x-y)^2}{2l^2N}\Big),\ \ \alpha=\frac{1}{\sqrt{2\pi l^2}}",
              font_size=26)
    err = place_col(VGroup(box, sol).arrange(DOWN, aligned_edge=LEFT, buff=0.2), -6.5, 3.1)
    assert_in_frame(axis, ticks, bars, cnt, curve, nums, err, dist_bars(Ps[NSTEPS], UY_B))
    assert_apart(err, cnt)
    assert_apart(err, nums)
    assert err.get_bottom()[1] > O5[1] + UY_B * 2 * W.heat_kernel(0, NSTEPS) + 0.1     # 不壓到放大後的山

    s.say("動畫：一開始，機率全部在 ", "$y=0$", wait=0.3)
    s.play(ShowCreation(axis), FadeIn(ticks), FadeIn(cnt), GrowFromEdge(bars, DOWN), run_time=1.2)
    s.wait(1.5)
    s.say("每一步：每格的機率一半往左、一半往右", wait=0.3)
    for N in range(1, NSTEPS + 1):
        num.set_value(N)
        s.play(Transform(bars, dist_bars(Ps[N], UY_A)), run_time=1.0 if N <= 3 else 0.3)
        if N <= 3:
            s.wait(0.6)
    s.say("30 步之後：攤開成一座小山（奇偶交替）", wait=3.0)
    s.play(Transform(bars, dist_bars(Ps[NSTEPS], UY_B)), run_time=1.5)
    s.say("放大縱軸，疊上 ", "$2\\times$", " 高斯密度（因子 2 同 W107）", wait=0.3)
    s.play(ShowCreation(curve), run_time=1.5)
    s.wait(2.0)
    s.play(FadeIn(n1), run_time=0.6)
    s.say("變異數 ", "$\\overline{y^2}=N=30$", "：和 W104 的 ", "$\\overline{x^2}=Nl^2$", " 一致", wait=3.5)
    s.play(FadeIn(n2), FadeIn(n3), run_time=0.8)
    s.say("中央那根：", f"${nb['P0']}$", " 對高斯的 ", f"${nb['g0']}$", wait=3.5)
    s.play(FadeIn(err, shift=0.2 * DOWN), run_time=0.8)
    s.say("書上 (1.48) 指數少了負號：少了它，曲線會往外爆開", wait=4.0)
    s.clear_stage()


def part_exercise(s):
    e0 = mixed("練習 1.5：把 ", "$P=\\frac{1}{\\sqrt{2\\pi N}}e^{-y^2/2N}$", " 代入 (1.47)（", "$l=1$", "）",
               font_size=26)
    e1 = Tex(r"\frac{\partial P}{\partial N}=-\frac{1}{2N}P+\frac{y^2}{2N^2}P"
             r"=\Big(\frac{y^2}{2N^2}-\frac{1}{2N}\Big)P", font_size=32)
    e2 = Tex(r"\frac{\partial P}{\partial y}=-\frac{y}{N}P,\qquad"
             r"\frac{\partial^2P}{\partial y^2}=\Big(\frac{y^2}{N^2}-\frac1N\Big)P", font_size=32)
    e3 = Tex(r"\frac12\,\frac{\partial^2P}{\partial y^2}=\Big(\frac{y^2}{2N^2}-\frac{1}{2N}\Big)P",
             font_size=32)
    e4 = zh("兩邊相等：高斯正是擴散方程的解", font_size=28, color=D.EXACT_COLOR)
    col = place_col(VGroup(e0, e1, e2, e3, e4).arrange(DOWN, aligned_edge=LEFT, buff=0.42), -6.0, 2.9)
    hl1 = SurroundingRectangle(e1, buff=0.1).set_stroke(D.GAUSS_COLOR, 2)
    hl3 = SurroundingRectangle(e3, buff=0.1).set_stroke(D.GAUSS_COLOR, 2)
    assert_in_frame(col, hl1, hl3)
    s.say("練習 1.5：驗證高斯真的滿足擴散方程", wait=0.3)
    s.play(FadeIn(e0), run_time=0.8)
    s.wait(2.0)
    s.say("左邊：對 ", "$N$", " 微分（前因子和指數都有 ", "$N$", "）", wait=0.3)
    s.play(Write(e1), run_time=1.5)
    s.wait(3.0)
    s.say("右邊：對 ", "$y$", " 微分兩次", wait=0.3)
    s.play(Write(e2), run_time=1.5)
    s.wait(2.5)
    s.play(Write(e3), run_time=1.2)
    s.play(ShowCreation(hl1), ShowCreation(hl3), run_time=0.8)
    s.say("兩邊都等於 ", "$\\big(\\tfrac{y^2}{2N^2}-\\tfrac{1}{2N}\\big)P$", wait=3.5)
    s.play(FadeIn(e4, shift=0.2 * UP), run_time=0.6)
    s.wait(2.0)


def part_close(s):
    s.say("隨機漫步放大來看，就是擴散", wait=3.5)


class W109DiffusionEquation(ChuScene):
    def construct(self):
        self.show_header("W109　從遞迴到擴散方程")
        check_math()
        part_question(self)
        part_second_difference(self)
        part_derivation(self)
        part_animation(self)
        part_exercise(self)
        part_close(self)
        self.wait(1.0)
