"""
W211 補充 2：倒晶格與布里淵區（§2.3，pp. 42–47，(2:S2-1)–(2:S2-34)）

問：格子上的函數做傅立葉，為什麼 k 只要在一個方塊（布里淵區）裡取就夠了？
  1. 混疊：n = −4..4 上 cos(πn/4) 與 cos(9πn/4) 完全一樣（兩條連續曲線在 9 個整數點重合）
     → f(k) = Σ_R F(R) e^{ik·R}（(2:S2-21)）是週期函數 f(k + K) = f(k)（(2:S2-22)），因為 K·R = 2π×整數（(2:S2-8)）。
  2. 倒晶格 (2:S2-5)–(2:S2-8)：b_i = (2π/Ω) a_j × a_k，a_i·b_j = 2πδ_ij；二維示意：正方格子（格距 1）、倒晶格（間距 2π）、
     布里淵區 [−π, π]²（L.bz_frame）。單胞體積 sc 1、bcc 0.7698、fcc 0.7071（GF.CELL_VOLUME，GF.reciprocal 核對）；
     布里淵區體積 (2π)³/Ω。書上 p. 43「l2, l2 and l3」→ l1, l2, l3（字幕帶過）。
  3. 反演 (2:S2-23)–(2:S2-25)：F(R) = Ω/(2π)³ ∫_BZ f(k) e^{−ik·R} d³k ——就是 (2.17) 的 (a/2π)³。
  4. 有限格子（週期邊界）：k = Σ (m_i/N_i) b_i，共 N0 = N1N2N3 個；Σ_R e^{ik·R} = N0 δ_{k,0}。
     N = 8 的單位根：m = 0 排成一直線（和 8），m = 1、3 頭接尾圍成封閉多邊形（和 0，GF.roots_sum）。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_gf as GF
from utils import lattice_draw as L
from utils import walk_draw as D

NR = 8
MS = (0, 1, 3)
BASES = {"sc": GF.SC_BASIS, "bcc": GF.BCC_BASIS, "fcc": GF.FCC_BASIS}


def numbers():
    return {"omega": {k: f"{GF.CELL_VOLUME[k]:.4f}" for k in BASES},
            "sums": {m: GF.roots_sum(m, NR) for m in MS}}


def check_math():
    nb = numbers()
    assert nb["omega"] == {"sc": "1.0000", "bcc": "0.7698", "fcc": "0.7071"}
    for k, B in BASES.items():
        b = GF.reciprocal(*B)
        om = abs(GF.triple(*B))
        assert abs(om - GF.CELL_VOLUME[k]) < 1e-12
        assert all(abs(GF.dot(B[i], b[j]) - 2 * math.pi * (i == j)) < 1e-12 for i in range(3) for j in range(3))
        assert abs(abs(GF.triple(*b)) - (2 * math.pi) ** 3 / om) < 1e-9                # 布里淵區體積
        assert all(abs(math.sqrt(GF.dot(v, v)) - 1) < 1e-12 for v in B) or k != "sc"   # sc 的基底長 1
    for n in range(-4, 5):                                                             # 混疊
        assert abs(math.cos(9 * math.pi * n / 4) - math.cos(math.pi * n / 4)) < 1e-12
    assert abs(math.cos(9 * math.pi * 0.5 / 4) - math.cos(math.pi * 0.5 / 4)) > 1         # 非整數點分得開
    assert abs(nb["sums"][0] - NR) < 1e-12 and abs(nb["sums"][1]) < 1e-12 and abs(nb["sums"][3]) < 1e-12
    assert all(abs(GF.roots_sum(m, NR)) < 1e-12 for m in range(1, NR))


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def tagged(tex, tag, font_size=30, color=WHITE):
    return VGroup(tx(tex, font_size, color), tx(r"\text{" + tag + "}", 22, GREY_B)).arrange(RIGHT, buff=0.3)


def dots_grid(n, step, origin, radius=0.05, color=GREY_B):
    o = np.array(origin, dtype=float)
    return VGroup(*[Dot(o + i * step * RIGHT + j * step * UP, radius=radius).set_fill(color, 1)
                    for i in range(-n, n + 1) for j in range(-n, n + 1)])


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：格子上的函數做傅立葉", wait=1.5)
    s.say("為什麼 ", "$\\vec k$", " 只要在一個方塊（布里淵區）裡取就夠了？", wait=3.0)


def part_alias(s):
    ax = D.plot_axes((-4.5, 4.5, 1), (-1.3, 1.3, 1), width=11.0, height=2.9, x_ticks=list(range(-4, 5)),
                     y_ticks=[-1, 1], font_size=20, x_label="n")
    ax.move_to(np.array([0, 1.35, 0]))
    c1 = ax.get_graph(lambda x: math.cos(math.pi * x / 4), x_range=[-4.5, 4.5, 0.02]).set_stroke(D.BAR_COLOR, 4)
    c2 = ax.get_graph(lambda x: math.cos(9 * math.pi * x / 4), x_range=[-4.5, 4.5, 0.005]).set_stroke(D.GAUSS_COLOR, 2.5)
    pts = VGroup(*[Dot(ax.c2p(n, math.cos(math.pi * n / 4)), radius=0.09).set_fill(YELLOW, 1) for n in range(-4, 5)])
    l1 = tx(r"\cos\tfrac{\pi n}{4}", 28, D.BAR_COLOR)
    l2 = tx(r"\cos\tfrac{9\pi n}{4}=\cos\Big(\tfrac{\pi n}{4}+2\pi n\Big)", 28, D.GAUSS_COLOR)
    legs = VGroup(l1, l2).arrange(RIGHT, buff=0.8).next_to(ax, DOWN, buff=0.2)
    f1 = tagged(r"f(\vec k)=\sum_{\vec R}F(\vec R)\,e^{i\vec k\cdot\vec R}", "(2:S2-21)")
    f2 = tagged(r"f(\vec k+\vec K)=f(\vec k)", "(2:S2-22)")
    f3 = mixed("因為 ", "$\\vec K\\cdot\\vec R=2\\pi\\times$", "整數", "（(2:S2-8)）", font_size=28)
    bot = VGroup(f1, VGroup(f2, f3).arrange(RIGHT, buff=0.6)).arrange(DOWN, buff=0.3)
    bot.move_to(np.array([0, -2.0, 0]))
    assert_in_frame(ax, legs, bot)
    assert_apart(legs, bot)
    s.say("一維格點 ", "$n=-4,\\dots,4$", " 上的波", wait=0.3)
    s.play(FadeIn(ax), ShowCreation(c1), FadeIn(l1), run_time=1.5)
    s.play(FadeIn(pts), run_time=0.8)
    s.say("我們只看得到整數點上的值（黃點）", wait=2.5)
    s.play(ShowCreation(c2), FadeIn(l2), run_time=2.0)
    s.say("$k$", " 加上 ", "$2\\pi$", "：曲線晃得快很多，但在整數點上完全一樣", wait=3.5)
    s.play(Indicate(pts, color=WHITE, scale_factor=1.4), run_time=1.0)
    s.say("格子上分不出 ", "$k$", " 和 ", "$k+2\\pi$", "：這叫混疊", wait=3.0)
    s.play(FadeIn(f1), run_time=1.0)
    s.play(FadeIn(f2), FadeIn(f3), run_time=1.0)
    s.say("所以格子函數的傅立葉轉換是週期函數", wait=3.0)
    s.say("週期 ", "$\\vec K$", "：和每個 ", "$\\vec R$", " 的內積都是 ", "$2\\pi$", " 的整數倍", wait=3.5)
    s.clear_stage()


def part_reciprocal(s):
    f1 = tagged(r"\vec R=n_1\vec a_1+n_2\vec a_2+n_3\vec a_3", "(2:S2-5)")
    f2 = tagged(r"\vec b_i=\frac{2\pi}{\Omega}\,\vec a_j\times\vec a_k", "(2:S2-7)")
    top = VGroup(f1, f2).arrange(RIGHT, buff=0.9)
    fit(top, 13.0).move_to(np.array([0, 2.75, 0]))
    ro = np.array([-4.6, -0.55, 0])
    real = dots_grid(3, 0.5, ro)
    a1 = Arrow(ro, ro + 0.5 * RIGHT, buff=0, thickness=3).set_color(YELLOW)
    a2 = Arrow(ro, ro + 0.5 * UP, buff=0, thickness=3).set_color(YELLOW)
    la = VGroup(tx(r"\vec a_1", 24, YELLOW).next_to(a1, DOWN, buff=0.08),
                tx(r"\vec a_2", 24, YELLOW).next_to(a2, LEFT, buff=0.08))
    t1 = zh("正方格子（格距 1）", font_size=24).next_to(real, UP, buff=0.3)
    u = 0.3
    S = 2 * math.pi * u
    ko = np.array([3.7, -0.75, 0])
    recip = dots_grid(1, S, ko, radius=0.08, color=TEAL)
    bz = L.bz_frame(width=S, origin=ko, font_size=22)
    fill = Square(side_length=S).move_to(ko).set_stroke(width=0).set_fill(TEAL, 0.25)
    b1 = Arrow(ko, ko + S * RIGHT, buff=0, thickness=3).set_color(TEAL)
    b2 = Arrow(ko, ko + S * UP, buff=0, thickness=3).set_color(TEAL)
    lb = VGroup(tx(r"\vec b_1", 24, TEAL).next_to(ko + S * RIGHT, DR, buff=0.08),
                tx(r"\vec b_2", 24, TEAL).next_to(ko + S * UP, UL, buff=0.08))
    t2 = mixed("倒晶格（間距 ", "$2\\pi$", "）", font_size=24).next_to(recip, UP, buff=0.25)
    mid = VGroup(tx(r"\vec a_i\cdot\vec b_j=2\pi\,\delta_{ij}", 30), tx(r"\text{(2:S2-8)}", 22, GREY_B),
                 mixed("$\\Omega=\\vec a_1\\cdot(\\vec a_2\\times\\vec a_3)$", font_size=26),
                 zh("單胞體積", font_size=22, color=GREY_B),
                 mixed("布里淵區 ", "$[-\\pi,\\pi]^2$", font_size=26, color=TEAL)).arrange(DOWN, buff=0.2)
    mid.move_to(np.array([-0.9, -0.55, 0]))
    assert_in_frame(top, real, la, t1, recip, bz, lb, t2, mid)
    assert_apart(top, VGroup(t1, t2))
    assert_apart(VGroup(real, la), mid)
    assert_apart(mid, VGroup(recip, bz, lb))
    s.say("三維的格子：", "$\\vec a_1,\\vec a_2,\\vec a_3$", " 是基本向量", wait=0.3)
    s.play(FadeIn(f1), run_time=1.0)
    s.say("書上 p. 43 把整數 ", "$l_1,l_2,l_3$", " 印成 ", "$l_2,l_2,l_3$", wait=2.5)
    s.play(FadeIn(t1), FadeIn(real), GrowArrow(a1), GrowArrow(a2), FadeIn(la), run_time=1.2)
    s.play(FadeIn(f2), run_time=1.0)
    s.say("倒晶格：用外積造出 ", "$\\vec b_1,\\vec b_2,\\vec b_3$", wait=3.0)
    s.play(FadeIn(mid[0]), FadeIn(mid[1]), run_time=0.8)
    s.say("$\\vec b_j$", " 和其他的 ", "$\\vec a_i$", " 都垂直，和自己的內積是 ", "$2\\pi$", wait=3.5)
    s.play(FadeIn(mid[2]), FadeIn(mid[3]), run_time=0.8)
    s.play(FadeIn(t2), FadeIn(recip), GrowArrow(b1), GrowArrow(b2), FadeIn(lb), run_time=1.2)
    s.say("格距 1 的正方格子：倒晶格的間距是 ", "$2\\pi$", wait=3.0)
    s.play(FadeIn(fill), FadeIn(bz), FadeIn(mid[4]), run_time=1.0)
    s.say("圍著原點的那一格：布里淵區", wait=2.5)
    s.say("$\\vec k$", " 只要在這一格裡取，其他的都是它平移 ", "$\\vec K$", wait=3.5)
    s.clear_stage()


def part_volumes(s):
    nb = numbers()
    rows = [[zh("晶格", font_size=24, color=YELLOW), mixed("單胞體積 ", "$\\Omega$", font_size=24, color=YELLOW),
             tx(r"a=1", 26, YELLOW)],
            [zh("簡單立方 sc", font_size=24), "a^3", "1"],
            [zh("體心立方 bcc", font_size=24), r"\tfrac{4}{3\sqrt3}\,a^3", nb["omega"]["bcc"]],
            [zh("面心立方 fcc", font_size=24), r"a^3/\sqrt2", nb["omega"]["fcc"]]]
    tab = D.tex_table(rows, [2.6, 2.4, 1.5], row_h=0.62, font_size=28)
    tab.move_to(np.array([-6.6 + tab.get_width() / 2, 3.1 - tab.get_height() / 2, 0]))
    v = mixed("布里淵區體積 ", "$=\\dfrac{(2\\pi)^3}{\\Omega}$", font_size=30)
    note = zh("（a＝最近鄰距離；(2:S2-9) 至 (2:S2-20)）", font_size=20, color=GREY_B)
    vn = VGroup(v, note).arrange(DOWN, aligned_edge=LEFT, buff=0.25).next_to(tab, RIGHT, buff=0.6)
    g1 = tagged(r"\int_{BZ}e^{i\vec k\cdot(\vec R-\vec R')}\,d^3k=\frac{(2\pi)^3}{\Omega}\,\delta_{\vec R,\vec R'}",
                "(2:S2-23)", 32)
    g2 = tagged(r"F(\vec R)=\frac{\Omega}{(2\pi)^3}\int_{BZ}f(\vec k)\,e^{-i\vec k\cdot\vec R}\,d^3k", "(2:S2-25)", 34,
                YELLOW)
    g3 = mixed("簡單立方 ", "$\\Omega=a^3$", "：就是 (2.17) 的 ", "$\\big(\\tfrac{a}{2\\pi}\\big)^3$", font_size=28)
    right = VGroup(g1, g2, g3).arrange(DOWN, buff=0.3)
    for m in right:
        fit(m, 12.5)
    right.arrange(DOWN, buff=0.3).move_to(np.array([0, -1.45, 0]))
    v, note = vn
    assert_in_frame(tab, v, note, right)
    assert_apart(VGroup(tab, v, note), right)
    s.say("三維的三種立方晶格（最近鄰距離 ", "$a$", "）", wait=0.3)
    s.play(FadeIn(tab), run_time=1.2)
    s.say("同樣的最近鄰距離：bcc、fcc 排得比較密", wait=3.0)
    s.play(FadeIn(v), FadeIn(note), run_time=0.8)
    s.say("單胞越小，布里淵區越大：兩者相乘是 ", "$(2\\pi)^3$", wait=3.0)
    s.play(FadeIn(g1), run_time=1.0)
    s.say("在布里淵區裡積分：不同的格點互相抵消", wait=3.0)
    s.play(FadeIn(g2), run_time=1.0)
    s.say("所以反演只要積一個布里淵區", wait=2.5)
    s.play(FadeIn(g3), run_time=0.8)
    s.say("這正是 (2.17) 前面的 ", "$(a/2\\pi)^3$", wait=3.0)
    s.clear_stage()


def part_finite(s):
    nb = numbers()
    f1 = mixed("週期邊界：", "$\\vec k=\\sum_i\\frac{m_i}{N_i}\\vec b_i$", "，共 ", "$N_0=N_1N_2N_3$", " 個", font_size=28)
    f2 = tagged(r"\sum_{\vec R}e^{i\vec k\cdot\vec R}=N_0\,\delta_{\vec k,\vec 0}", "(2:S2-34)", 30)
    top = VGroup(f1, f2).arrange(RIGHT, buff=0.6)
    fit(top, 13.2).move_to(np.array([0, 2.75, 0]))
    ro = np.array([-5.1, 0.0, 0])
    circ = L.roots_of_unity(NR, radius=1.3, origin=ro, m=1)
    t0 = mixed("$N=8$", "：", "$e^{2\\pi imn/8}$", font_size=26).next_to(circ, DOWN, buff=0.35)
    chains, labs = VGroup(), VGroup()
    for m, x in zip(MS, (-0.95, 3.2, 5.5)):
        ch = L.root_chain(NR, m, length={0: 0.5, 1: 0.55, 3: 1.1}[m], origin=ORIGIN,
                          color=[TEAL, YELLOW, D.GAUSS_COLOR][MS.index(m)])
        ch.move_to(np.array([x, 0.0, 0]))
        chains.add(ch)
        val = nb["sums"][m]
        lab = mixed("$m=" + str(m) + "$", "：和 ", "$" + str(round(abs(val))) + "$", font_size=26)
        labs.add(lab.move_to(np.array([x, -1.5, 0])))
    assert_in_frame(top, circ, t0, chains, labs)
    assert_apart(circ, chains[0])
    for a, b in zip(chains, chains[1:]):
        assert_apart(a, b)
    s.say("有限的格子（週期邊界）：", "$\\vec k$", " 只能取分立的值", wait=0.3)
    s.play(FadeIn(f1), run_time=1.0)
    s.say("布里淵區裡剛好 ", "$N_0$", " 個：和格點一樣多", wait=3.0)
    s.play(FadeIn(f2), run_time=1.0)
    s.say("一維 ", "$N=8$", " 的例子：", "$e^{ikn}$", " 是 8 次單位根", wait=0.3)
    s.play(ShowCreation(circ[0]), LaggedStart(*[GrowArrow(a) for a in circ[1]], lag_ratio=0.1), FadeIn(t0), run_time=1.6)
    s.wait(1.5)
    for k in range(3):
        s.play(LaggedStart(*[GrowArrow(a) for a in chains[k]], lag_ratio=0.15), run_time=2.0)
        s.play(FadeIn(labs[k]), run_time=0.5)
        if k == 0:
            s.say("$m=0$", "：每個都是 1，排成一直線，和是 8", wait=2.5)
        elif k == 1:
            s.say("$m=1$", "：頭接尾繞一圈回到起點，和是 0", wait=2.5)
        else:
            s.say("$m=3$", "：繞三圈，還是回到起點，和是 0", wait=2.5)
    s.say("只有 ", "$\\vec k=\\vec 0$", "（或差一個 ", "$\\vec K$", "）加起來不抵消", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("格子上的波，", "$\\vec k$", " 和 ", "$\\vec k+\\vec K$", " 分不出來", wait=3.0)
    s.say("所以只需要一個布里淵區", wait=2.5)


class W211ReciprocalLattice(ChuScene):
    def construct(self):
        self.show_header("W211　補充 2：倒晶格與布里淵區")
        check_math()
        part_question(self)
        part_alias(self)
        part_reciprocal(self)
        part_volumes(self)
        part_finite(self)
        part_close(self)
        self.wait(1.0)
