"""
W108 遞迴關係與 Pascal 三角形（§1.4，p. 13，(1.42)–(1.43)）

問：不用公式，能不能一步一步把所有走法數出來？
  1. (1.42)：C(N;x,y) = C(N−1;x,y−l) + C(N−1;x,y+l)——最後一步一定來自左鄰或右鄰。
     pascal_triangle(pascal_rows(8)) 一列列長出來，每格有兩個小箭頭從上一列的鄰居指下來；第 6 列 y = 0：20 = 10 + 10。
  2. Markov：兩條不同的過去（HHTT、THTH）走到同一個現在，下一步一樣是各 ½（一階 Markov 過程）。
  3. (1.43)：6!/(3!3!) = 5!/(2!3!) + 5!/(3!2!)，即 20 = 10 + 10；列和 1, 2, 4, …, 256 = 2^N，第 10 列和 1024。
  4. 機率版：P(N+1,y) = ½[P(N,y−1) + P(N,y+1)]；例 P(6,0) = ½(10/32 + 10/32) = 20/64（下一幕的起點）。
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

NMAX = 8
ROWS = W.pascal_rows(NMAX)
DX, DY = 0.36, 0.56
TRI_O = np.array([-3.0, 2.75, 0])          # 第 0 列（y = 0）的位置
RP = 1.3                                   # 右側面板左緣
RP_W = 5.4                                 # 右側面板寬
PATH_A = [0, 1, 2, 1, 0]                   # HHTT
PATH_B = [0, -1, 0, -1, 0]                 # THTH
MK_O = np.array([1.75, -1.4, 0])           # Markov 小圖原點
MK_UX, MK_UY = 0.75, 0.46


def check_math():
    rows = W.pascal_rows(10)
    assert rows[:NMAX + 1] == ROWS
    for N in range(1, 11):
        for y, c in rows[N].items():
            assert c == W.walk_count(N, y) == rows[N - 1].get(y - 1, 0) + rows[N - 1].get(y + 1, 0)
    assert rows[6][0] == 20 and rows[5][-1] == 10 and rows[5][1] == 10
    f = math.factorial
    assert f(6) // (f(3) * f(3)) == 20 and f(5) // (f(2) * f(3)) == 10 and f(5) // (f(3) * f(2)) == 10
    assert [sum(r.values()) for r in rows] == [2 ** N for N in range(11)]
    assert sum(rows[8].values()) == 256 and sum(rows[10].values()) == 1024
    assert W.walk_prob(5, -1) == W.walk_prob(5, 1) == Fraction(10, 32)
    assert W.walk_prob(6, 0) == Fraction(1, 2) * (W.walk_prob(5, -1) + W.walk_prob(5, 1)) == Fraction(20, 64)
    for N in range(0, 10):                     # 機率版對每一格都成立：P(N+1,y) = ½[P(N,y−1) + P(N,y+1)]
        for y in range(-N - 1, N + 2):
            assert W.walk_prob(N + 1, y) == Fraction(1, 2) * (W.walk_prob(N, y - 1) + W.walk_prob(N, y + 1))
    for N in range(1, 11):                     # (1.43) 的一般形：C(N,n) = C(N−1,n−1) + C(N−1,n)
        for n in range(1, N):
            assert math.comb(N, n) == math.comb(N - 1, n - 1) + math.comb(N - 1, n)
    assert W.positions(W.coins_to_steps("HHTT")) == PATH_A
    assert W.positions(W.coins_to_steps("THTH")) == PATH_B
    assert W.diffusion_step({0: 1}) == {-1: Fraction(1, 2), 1: Fraction(1, 2)}


# ---------------------------------------------------------------- 小工具
def fit(m, width=RP_W):
    if m.get_width() > width:
        m.set_width(width)
    return m


def eq_below(tex, tag=None, font_size=30, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return fit(m)
    return fit(VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(DOWN, aligned_edge=RIGHT, buff=0.1))


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


def tri_pos(N, y):
    return TRI_O + y * DX * RIGHT + N * DY * DOWN


def in_arrows(N, y, color=GREY_B):
    """指向 (N, y) 的小箭頭：從上一列的 y−1、y+1（存在的才畫）。"""
    g = VGroup()
    for py in (y - 1, y + 1):
        if abs(py) <= N - 1:
            a, b = tri_pos(N - 1, py), tri_pos(N, y)
            u = (b - a) / np.linalg.norm(b - a)
            g.add(Arrow(a + 0.17 * u, b - 0.17 * u, buff=0, thickness=2.0).set_color(color))
    return g


def labels():
    nl = VGroup(*[Tex(f"N={N}", font_size=20).set_color(D.MUTED).move_to(np.array([-6.35, tri_pos(N, 0)[1], 0]))
                  for N in range(NMAX + 1)])
    yl = VGroup(*[Tex(str(y), font_size=20).set_color(D.MUTED).move_to(np.array([tri_pos(NMAX, y)[0], -2.3, 0]))
                  for y in range(-NMAX, NMAX + 1, 2)])
    yl.add(zh("位置", font_size=20, color=D.MUTED).move_to(np.array([-6.4, -2.3, 0])))
    return nl, yl


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：不用公式，能不能一步一步把所有走法數出來？", wait=5.0)


def part_triangle(s):
    tri = D.pascal_triangle(ROWS, dx=DX, dy=DY, origin=TRI_O, font_size=24)
    arrows = {(N, y): in_arrows(N, y) for N in range(1, NMAX + 1) for y in ROWS[N]}
    row_arrows = [VGroup()] + [VGroup(*[a for y in sorted(ROWS[N]) for a in arrows[(N, y)]])
                               for N in range(1, NMAX + 1)]
    nl, yl = labels()
    f1 = eq_below(r"C(N;x,y)=C(N-1;x,y-l)+C(N-1;x,y+l)", "(1.42)", font_size=28)
    f2 = fit(mixed("最後一步一定來自 ", "$y-l$", " 或 ", "$y+l$", font_size=26))
    f3 = fit(mixed("（", "$C(N;x,y)$", "：從 ", "$x$", " 出發、", "$N$", " 步後在 ", "$y$", " 的走法數）",
                   font_size=22, color=GREY_B))
    top = place_col(VGroup(f1, f2, f3).arrange(DOWN, aligned_edge=LEFT, buff=0.25), RP, 3.0)
    big = Tex("20=10+10", font_size=48).set_color(YELLOW).move_to(np.array([RP + RP_W / 2, 0.3, 0]))
    assert_in_frame(tri, nl, yl, top, big, *row_arrows[1:])
    assert_apart(tri, top)
    assert_apart(VGroup(*row_arrows[1:]), top)

    s.play(Write(f1[0]), FadeIn(f1[1]), run_time=1.5)
    s.play(FadeIn(f3), run_time=0.6)
    s.say("遞迴 (1.42)：", "$N$", " 步後在 ", "$y$", " 的走法數", wait=4.0)
    s.play(FadeIn(f2, shift=0.2 * DOWN), run_time=0.6)
    s.say("最後一步只有兩種可能：從左鄰或右鄰走過來", wait=4.0)
    s.play(FadeIn(nl[0]), FadeIn(tri[0]), FadeIn(yl), run_time=0.8)
    s.say("一列一列長出來：起點只有 1 種", wait=1.5)
    for N in range(1, NMAX + 1):
        slow = N <= 2
        s.play(FadeIn(nl[N]), LaggedStartMap(GrowArrow, row_arrows[N], lag_ratio=0.08),
               run_time=1.0 if slow else 0.7)
        s.play(FadeIn(tri[N], shift=0.1 * DOWN), run_time=0.6 if slow else 0.4)
        if N == 1:
            s.say("第 1 步：", "$-1$", " 和 ", "$+1$", " 各 1 種", wait=3.0)
        elif N == 2:
            s.say("每一格 ", "$=$", " 左上 ", "$+$", " 右上：", "$1+1=2$", wait=3.5)
        else:
            s.wait(0.3)
    s.wait(1.5)
    s.say("看第 6 列的中央 ", "$y=0$", wait=0.3)
    hl = VGroup(D.pascal_entry(tri, 5, -1), D.pascal_entry(tri, 5, 1), D.pascal_entry(tri, 6, 0))
    s.play(*[m.animate.set_color(YELLOW) for m in hl], arrows[(6, 0)].animate.set_color(YELLOW), run_time=0.8)
    s.play(FadeIn(big, scale=1.2), run_time=0.8)
    s.say("到 ", "$y=0$", " 的 20 種＝從 ", "$-1$", " 來的 10 種＋從 ", "$+1$", " 來的 10 種", wait=5.0)
    return tri, arrows, top, big, hl


def part_markov(s, big):
    ax = D.walk_axes(5, 2, ux=MK_UX, uy=MK_UY, origin=MK_O, x_every=1, y_every=1)
    ax[3].remove(*list(ax[3][:5]))          # 步數刻度的數字會壓到往下走的路線：拿掉（只留位置的數字）
    pa = D.walk_path(PATH_A, ux=MK_UX, uy=MK_UY, origin=MK_O, color=D.WALK_COLOR)
    pb = D.walk_path(PATH_B, ux=MK_UX, uy=MK_UY, origin=MK_O, color=TEAL)
    here = D.path_point(4, 0, MK_UX, MK_UY, MK_O)
    nxt = VGroup()
    for y in (1, -1):
        end = D.path_point(5, y, MK_UX, MK_UY, MK_O)
        nxt.add(DashedLine(here, end).set_stroke(WHITE, 2))
        nxt.add(Tex(r"\tfrac12", font_size=26).next_to(end, RIGHT, buff=0.12))
    title = mixed("一階 ", "Markov", " 過程", font_size=28, color=YELLOW)
    note = zh("下一步只看現在的位置，不管怎麼走到這裡", font_size=22)
    head = place_col(VGroup(title, note).arrange(DOWN, aligned_edge=LEFT, buff=0.15), RP, 1.15)
    assert_in_frame(ax, pa, pb, nxt, head)
    s.play(FadeOut(big), FadeIn(ax), run_time=0.8)
    s.play(FadeIn(head), run_time=0.6)
    s.say("為什麼能這樣一列一列算？", wait=3.5)
    s.play(ShowCreation(pa[0]), FadeIn(pa[1]), run_time=1.2)
    s.play(ShowCreation(pb[0]), FadeIn(pb[1]), run_time=1.2)
    s.say("兩條不同的過去，第 4 步都在 ", "$y=0$", wait=3.5)
    s.play(LaggedStartMap(FadeIn, nxt, lag_ratio=0.3), run_time=1.0)
    s.say("下一步一樣是左右各 ", "$\\tfrac12$", "：只看現在，不看過去", wait=4.0)
    s.say("這叫一階 Markov 過程：遞迴 (1.42) 就靠這一點", wait=5.0)
    s.play(FadeOut(VGroup(ax, pa, pb, nxt, head)), run_time=0.6)


def part_identity(s, tri, hl):
    g1 = eq_below(r"\frac{6!}{3!\,3!}=\frac{5!}{2!\,3!}+\frac{5!}{3!\,2!}", "(1.43)", font_size=34)
    g2 = Tex("20=10+10", font_size=36).set_color(YELLOW)
    g3 = Tex(r"\sum_y C(N,y)=2^N", font_size=34).set_color(TEAL)
    g4 = mixed("第 10 列：", "$\\sum_y C(10,y)=2^{10}=1024$", font_size=26)
    col = place_col(VGroup(g1, g2, g3, g4).arrange(DOWN, aligned_edge=LEFT, buff=0.35), RP, 1.2)
    sums = VGroup(*[Tex(str(sum(ROWS[N].values())), font_size=22).set_color(TEAL)
                    .move_to(np.array([0.6, tri_pos(N, 0)[1], 0])) for N in range(NMAX + 1)])
    assert_in_frame(col, sums)
    assert_apart(sums, col)
    assert_apart(sums, tri)
    s.say("用公式寫：Pascal 恆等式 (1.43)", wait=0.3)
    s.play(Write(g1[0]), FadeIn(g1[1]), run_time=1.5)
    s.play(FadeIn(g2), run_time=0.6)
    s.say("6 步走到 0 的走法 ", "$=$", " 5 步走到 ", "$\\pm1$", " 的走法相加", wait=4.0)
    s.play(LaggedStartMap(FadeIn, sums, shift=0.2 * LEFT, lag_ratio=0.15), run_time=1.6)
    s.say("每一列加起來：1, 2, 4, 8, …, 256", wait=3.5)
    s.play(FadeIn(g3), run_time=0.6)
    s.say("每多走一步，總走法數加倍：", "$2^N$", wait=3.5)
    s.play(FadeIn(g4), run_time=0.6)
    s.say("第 10 列加起來就是 W101 的 1024 個人", wait=4.5)
    s.play(FadeOut(col), *[m.animate.set_color(WHITE) for m in hl], run_time=0.6)
    return sums


def part_probability(s, tri, arrows):
    h1 = fit(Tex(r"P(N+1,y)=\tfrac12\big[P(N,y-1)+P(N,y+1)\big]", font_size=32), RP_W - 0.4)
    h2 = fit(Tex(r"P(6,0)=\tfrac12\Big(\frac{10}{32}+\frac{10}{32}\Big)=\frac{20}{64}", font_size=34))
    h3 = mixed("（每一列除以 ", "$2^N$", "：走法數變成機率）", font_size=22, color=GREY_B)
    h4 = zh("下一幕：從這個遞迴走到擴散方程", font_size=24, color=GREY_B)
    col = place_col(VGroup(h1, h3, h2, h4).arrange(DOWN, aligned_edge=LEFT, buff=0.35), RP, 1.2)
    box = SurroundingRectangle(h1, buff=0.15).set_stroke(D.GAUSS_COLOR, 2)
    assert_in_frame(col, box)
    hl = VGroup(D.pascal_entry(tri, 5, -1), D.pascal_entry(tri, 5, 1), D.pascal_entry(tri, 6, 0))
    s.say("換成機率：每一列除以 ", "$2^N$", wait=0.3)
    s.play(Write(h1), FadeIn(h3), run_time=1.5)
    s.play(ShowCreation(box), run_time=0.6)
    s.say("機率版的遞迴：左上、右上各取一半", wait=4.0)
    s.play(*[m.animate.set_color(D.GAUSS_COLOR) for m in hl], arrows[(6, 0)].animate.set_color(D.GAUSS_COLOR),
           FadeIn(h2), run_time=1.0)
    s.say("例：", "$P(6,0)=\\tfrac12\\big(\\tfrac{10}{32}+\\tfrac{10}{32}\\big)=\\tfrac{20}{64}$", wait=4.5)
    s.play(FadeIn(h4), run_time=0.6)
    s.say("這個式子就是下一幕「擴散方程」的起點", wait=4.5)


def part_close(s):
    s.say("每一格 ", "$=$", " 左上 ", "$+$", " 右上：Pascal 三角形就是走法表", wait=5.0)


class W108PascalRecursion(ChuScene):
    def construct(self):
        self.show_header("W108　遞迴關係與 Pascal 三角形")
        check_math()
        part_question(self)
        tri, arrows, top, big, hl = part_triangle(self)
        part_markov(self, big)
        part_identity(self, tri, hl)
        part_probability(self, tri, arrows)
        part_close(self)
        self.wait(2.0)
