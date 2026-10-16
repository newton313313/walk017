"""
W104 均方位移與 √N 定律（§1.2，pp. 6–7，式 (1.10)–(1.14)）

x = ΣΔ_j，x² = ΣΣ Δ_jΔ_k：N = 4 的 4×4 格子，對角 4 格永遠是 l² = 1、非對角 12 格平均 0 → (1.13) x̄² = N l²；
小核對 Σx² = 64（16 條路線，平均 4）；(1.14) √(x̄²) = l√N；ensemble(2000, 100, seed=7) 的扇形與 ±√N 包絡；冪次律。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
import math
from fractions import Fraction

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

WALK_A, WALK_B = "HTHH", "HHTH"
WALKS4 = W.all_walks(4)
DIST4 = W.distribution(4)
SUMSQ4 = sum(sum(w) ** 2 for w in WALKS4)
BIG = W.ensemble(2000, 100, seed=7)
BIG_PATHS = [W.positions(w) for w in BIG[:80]]
MS100 = W.mean_square(W.finals(BIG))

CELL = 0.72
GRID_O = np.array([-5.0, 1.15, 0.0])
RX = -1.9                                        # 右欄左緣
FAN_O = np.array([-6.2, 0.15, 0.0])
FUX, FUY = 0.085, 0.088


def check_math():
    assert W.coins_to_steps(WALK_A) == W.MAIN_STEPS[:4] == [1, -1, 1, 1]
    assert W.coins_to_steps(WALK_B) == [1, 1, -1, 1]
    assert W.pair_grid_counts(4) == (4, 12)
    for j in range(4):
        for k in range(4):
            avg = Fraction(sum(w[j] * w[k] for w in WALKS4), len(WALKS4))
            assert avg == (1 if j == k else 0)                               # (1.12)
    assert len(WALKS4) == 16 and SUMSQ4 == 64 and Fraction(SUMSQ4, 16) == 4 == W.moment_exact(4, 2)   # (1.13)
    assert [DIST4[d] for d in (-4, -2, 0, 2, 4)] == [1, 4, 6, 4, 1]
    assert [DIST4[d] * d * d for d in (-4, -2, 0, 2, 4)] == [16, 16, 0, 16, 16]
    assert math.isqrt(100) ** 2 == 100 and math.isqrt(100) == 10 and math.isqrt(10000) == 100   # (1.14)
    assert W.moment_exact(100, 2) == 100
    assert len(BIG) == 2000 and float(MS100) == 98.922 and round(math.sqrt(float(MS100)), 3) == 9.946
    assert max(abs(x) for p in BIG_PATHS for x in p) <= 30


def eq_row(tex, num, font_size=36):
    row = VGroup(Tex(tex, font_size=font_size))
    if num:
        row.add(Tex(num, font_size=24).set_color(GREY_B))
    return row.arrange(RIGHT, buff=0.35)


def at_left(m, y, x=RX):
    return m.move_to(np.array([0, y, 0])).align_to(np.array([x, 0, 0]), LEFT)


# ---------------------------------------------------------------- 1. 問題
def ask(self):
    guesses = VGroup(zh("50？", font_size=60, color=GREY_B), zh("100？", font_size=60, color=GREY_B),
                     zh("10？", font_size=60, color=GREY_B)).arrange(RIGHT, buff=1.4).move_to(0.4 * UP)
    self.say("走 100 步，大概離家多遠？", wait=1.5)
    self.play(FadeIn(guesses[0], shift=0.2 * UP))
    self.play(FadeIn(guesses[1], shift=0.2 * UP))
    self.say("50？100？先猜一個", wait=2.0)
    self.play(FadeIn(guesses[2], shift=0.2 * UP))
    self.wait(2.5)
    self.play(FadeOut(guesses))


# ---------------------------------------------------------------- 2. Δ 格子
def cell_values(g, steps):
    out = VGroup()
    for j in range(4):
        for k in range(4):
            v = steps[j] * steps[k]
            t = Tex("+1" if v > 0 else "-1", font_size=24)
            t.set_color(BLACK if j == k else (WHITE if v > 0 else RED_B))
            out.add(t.move_to(D.grid_cell(g, j, k)))
    return out


def coin_tops(coins):
    return VGroup(*[D.coin(c, radius=0.22).move_to(GRID_O + k * CELL * RIGHT + 1.35 * UP)
                    for k, c in enumerate(coins)])


def grid_part(self):
    r1 = at_left(eq_row(r"x=\sum_{j=1}^{N}\Delta_j", "(1.11)"), 2.6)
    r2 = at_left(eq_row(r"x^2=\Big(\sum_j\Delta_j\Big)^2=\sum_j\sum_k\Delta_j\Delta_k", "(1.10)"), 1.45)
    assert_apart(r1, r2)
    g = D.delta_grid(4, cell=CELL, origin=GRID_O)
    collab = VGroup(*[Tex(rf"\Delta_{k + 1}", font_size=26).move_to(GRID_O + k * CELL * RIGHT + 0.68 * UP)
                      for k in range(4)])
    rowlab = VGroup(*[Tex(rf"\Delta_{j + 1}", font_size=26).move_to(GRID_O + j * CELL * DOWN + 0.68 * LEFT)
                      for j in range(4)])
    coins = coin_tops(WALK_A)
    vals = cell_values(g, W.coins_to_steps(WALK_A))
    left = VGroup(g, collab, rowlab, coins)
    assert_in_frame(r1, r2, left)
    assert_apart(left, r1)
    assert_apart(left, r2)
    self.say("每一步的位移 ", "$\\Delta_j=\\pm l$", "：終點是它們的和")
    self.play(Write(r1), run_time=1.2)
    self.wait(2.5)
    self.say("平方之後展開：每一對 ", "$(j,k)$", " 都有一項")
    self.play(Write(r2), run_time=1.6)
    self.wait(3.0)
    self.say("$N=4$", "：格子 ", "$(j,k)$", " 放 ", "$\\Delta_j\\Delta_k$", "（例：正反正正）")
    self.play(FadeIn(g, lag_ratio=0.03), FadeIn(collab), FadeIn(rowlab), FadeIn(coins, lag_ratio=0.2), run_time=1.5)
    self.play(LaggedStart(*[FadeIn(v) for v in vals], lag_ratio=0.04), run_time=1.6)
    self.wait(3.0)
    self.say("$x^2$", " 就是 16 格的總和：這條路線 ", "$x=2$", "，", "$x^2=4$", wait=3.0)
    coins_b = coin_tops(WALK_B)
    vals_b = cell_values(g, W.coins_to_steps(WALK_B))
    self.say("換一條路線（正正反正）：對角格不變，非對角格會變號")
    self.play(*[FadeTransform(a, b) for a, b in zip(coins, coins_b)],
              *[Transform(a, b) for a, b in zip(vals, vals_b)], run_time=1.5)
    self.wait(3.5)
    return g, vals, r1, r2


def average_part(self, g, vals):
    zeros = VGroup()
    for j in range(4):
        for k in range(4):
            t = Tex("1" if j == k else "0", font_size=28).set_color(BLACK if j == k else GREY_B)
            zeros.add(t.move_to(D.grid_cell(g, j, k)))
    r3 = at_left(mixed("對角 4 格：", "$\\Delta_j^2=l^2=1$", font_size=28), 0.6)
    r4 = at_left(mixed("非對角 12 格：", "$\\pm1$", "，彼此獨立，平均 0", font_size=28), -0.1)
    n4 = Tex("(1.12)", font_size=24).set_color(GREY_B).next_to(r4, RIGHT, buff=0.3)
    r5 = at_left(eq_row(r"\overline{x^2}=N\,l^2", "(1.13)", font_size=44), -1.05)
    box = SurroundingRectangle(r5[0], buff=0.12).set_stroke(YELLOW, 2)
    chk = at_left(Tex(r"N=4:\ \frac{1\cdot16+4\cdot4+6\cdot0+4\cdot4+1\cdot16}{16}=\frac{64}{16}=4", font_size=32),
                  -2.2)
    assert_in_frame(r3, r4, n4, r5, box, chk)
    for a, b in ((r3, r4), (r4, r5), (box, chk), (g, r3)):
        assert_apart(a, b)
    self.say("對全部路線平均：", "$\\overline{\\Delta_j\\Delta_k}=0$", "（", "$j\\ne k$", "）")
    self.play(*[Transform(a, b) for a, b in zip(vals, zeros)], run_time=1.5)
    self.play(FadeIn(r3, shift=0.2 * LEFT))
    self.wait(2.5)
    self.play(FadeIn(r4, shift=0.2 * LEFT), FadeIn(n4))
    self.wait(3.0)
    self.say("剩下對角線上的 ", "$N$", " 個 1：")
    self.play(Write(r5), ShowCreation(box), run_time=1.5)
    self.wait(3.5)
    self.say("核對：", "$N=4$", " 的 16 條路線，", "$x^2$", " 加起來是 64")
    self.play(Write(chk), run_time=1.8)
    self.wait(4.0)


# ---------------------------------------------------------------- 3. √N 定律
def sqrt_law(self):
    f = eq_row(r"\sqrt{\overline{x^2}}=l\sqrt{N}", "(1.14)", font_size=56).move_to(1.9 * UP)
    rows = [["N", r"\sqrt{\overline{x^2}}"], ["4", "2"], ["100", "10"], ["10000", "100"]]
    tab = D.tex_table(rows, [2.4, 2.4], row_h=0.62, font_size=36)
    tab.move_to(np.array([0, -0.75, 0]))
    assert_in_frame(f, tab)
    assert_apart(f, tab)
    self.say("開根號：典型的距離")
    self.play(Write(f), run_time=1.5)
    self.wait(2.5)
    self.play(FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:2])))
    for i in (1, 2, 3):
        self.play(FadeIn(VGroup(*tab[0][2 * i:2 * i + 2]), shift=0.15 * UP), run_time=0.8)
        self.wait(0.8)
    self.say("揭曉：走 100 步，大約只離家 10 步", wait=1.0)
    self.play(Indicate(D.table_cell(tab, 2, 1, 2), scale_factor=1.5), run_time=1.2)
    self.wait(2.5)
    self.say("走 10000 步，也只離家大約 100 步", wait=3.5)


# ---------------------------------------------------------------- 4. 扇形
def fan_part(self):
    ax = D.walk_axes(100, 30, ux=FUX, uy=FUY, origin=FAN_O, x_every=20, y_every=10)
    ax[4].next_to(ax[0].get_end(), RIGHT, buff=0.12)          # 軸名移到箭頭旁，不壓刻度 100、30
    ax[5].next_to(ax[1].get_end(), RIGHT, buff=0.12)
    for t in ax[3]:
        assert_apart(ax[4], t)
        assert_apart(ax[5], t)
    fan = D.fan(BIG_PATHS, ux=FUX, uy=FUY, origin=FAN_O)
    env = D.sqrt_envelope(100, FUX, FUY, FAN_O)
    elab = VGroup(Tex(r"+\sqrt{N}", font_size=28), Tex(r"-\sqrt{N}", font_size=28)).set_color(D.GAUSS_COLOR)
    elab[0].next_to(env[0].get_end(), RIGHT, buff=0.1)
    elab[1].next_to(env[1].get_end(), RIGHT, buff=0.1)
    cut = DashedLine(D.path_point(100, -30, FUX, FUY, FAN_O), D.path_point(100, 30, FUX, FUY, FAN_O),
                     dash_length=0.08).set_stroke(WHITE, 1.5)
    info = VGroup(zh("2000 人各走 100 步", font_size=26),
                  zh("（畫出其中 80 條）", font_size=22, color=GREY_B),
                  mixed("第 100 步的樣本：", font_size=24),
                  Tex(r"\overline{x^2}=98.922", font_size=36),
                  Tex(r"\sqrt{\overline{x^2}}\approx 9.946", font_size=36),
                  mixed("理論：", "$\\overline{x^2}=100$", font_size=24, color=D.EXACT_COLOR),
                  mixed("$\\sqrt{100}=10$", font_size=24, color=D.EXACT_COLOR))
    info.arrange(DOWN, buff=0.28, aligned_edge=LEFT).move_to(np.array([0, 0.3, 0])).align_to(np.array([3.3, 0, 0]),
                                                                                              LEFT)
    assert_in_frame(ax, fan, env, elab, info, cut)
    assert_apart(VGroup(ax, elab), info)
    self.say("很多人一起走：", "$N=100$", " 步的路線疊成扇形")
    self.play(FadeIn(ax), FadeIn(info[:2]))
    self.play(LaggedStart(*[ShowCreation(p) for p in fan], lag_ratio=0.03), run_time=4.0)
    self.bring_to_front(ax)                      # 刻度數字壓在扇形上面才看得清楚
    self.wait(1.5)
    self.say("橘色：", "$\\pm\\sqrt{N}$", "；扇形的寬度跟著 ", "$\\sqrt{N}$", " 長")
    self.play(ShowCreation(env), FadeIn(elab), run_time=1.8)
    self.wait(3.5)
    self.say("在第 100 步切一刀，算 2000 人的平方平均")
    self.play(ShowCreation(cut), FadeIn(info[2]))
    self.play(Write(info[3]), run_time=1.2)
    self.play(Write(info[4]), run_time=1.2)
    self.wait(2.0)
    self.play(FadeIn(info[5:], shift=0.15 * UP))
    self.say("樣本 98.922，理論 100：差不到 2%", wait=4.0)


# ---------------------------------------------------------------- 5. 冪次律
def power_law(self):
    ax = D.plot_axes((0, 4.6, 1), (0, 2.6, 1), width=5.2, height=3.4, x_ticks=[2, 4],
                     x_tick_labels=["10^2", "10^4"], y_ticks=[1, 2], y_tick_labels=["10", "100"],
                     x_label="N", y_label=r"\sqrt{\overline{x^2}}")
    ax.move_to(np.array([-3.2, -0.1, 0]))
    one = Tex("1", font_size=20).set_color(GREY_B).next_to(ax.c2p(0, 0), DL, buff=0.08)
    pts = VGroup(*[Dot(ax.c2p(t, t / 2), radius=0.08).set_fill(YELLOW, 1) for t in (0, 2, 4)])
    line = Line(ax.c2p(0, 0), ax.c2p(4.4, 2.2)).set_stroke(YELLOW, 3)
    slope = mixed("斜率 ", "$\\tfrac12$", font_size=26, color=YELLOW).next_to(ax.c2p(3.0, 1.5), UL, buff=0.15)
    f = Tex(r"\sqrt{\overline{x^2}}\propto N^{1/2}", font_size=48).move_to(np.array([3.6, 1.3, 0]))
    t1 = zh("分數次方的「冪次律」", font_size=28).move_to(np.array([3.6, 0.1, 0]))
    t2 = zh("後面幾章（標度、臨界現象）", font_size=24, color=GREY_B).move_to(np.array([3.6, -0.7, 0]))
    t3 = zh("會一直遇到這種式子", font_size=24, color=GREY_B).move_to(np.array([3.6, -1.2, 0]))
    assert_in_frame(ax, one, pts, slope, f, t1, t2, t3)
    assert_apart(VGroup(ax, slope), VGroup(f, t1, t2, t3))
    self.say("兩邊取對數：一條斜率 1/2 的直線")
    self.play(FadeIn(ax), FadeIn(one))
    self.play(LaggedStart(*[FadeIn(p, scale=0.5) for p in pts], lag_ratio=0.4), run_time=1.5)
    self.play(ShowCreation(line), FadeIn(slope))
    self.wait(2.5)
    self.say("距離 ", "$\\propto N^{1/2}$", "：分數次方")
    self.play(Write(f), run_time=1.4)
    self.play(FadeIn(t1))
    self.wait(2.0)
    self.play(FadeIn(t2), FadeIn(t3))
    self.wait(3.5)


class W104SqrtNLaw(ChuScene):
    def construct(self):
        self.show_header("W104　均方位移與 √N 定律")
        check_math()
        ask(self)
        g, vals, r1, r2 = grid_part(self)
        average_part(self, g, vals)
        self.clear_stage()
        sqrt_law(self)
        self.clear_stage()
        fan_part(self)
        self.clear_stage()
        power_law(self)
        self.say("距離只長成步數的平方根：", "$\\sqrt{N}$", " 定律", wait=4.0)
        self.wait(1.0)
