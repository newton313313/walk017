"""
W105 用生成函數算所有動差（§1.2 例題，pp. 7–10，式 (1.15)–(1.27)，Exercise 1.2–1.3）

(1.15) 直接加總；g(y) = (1+e^y)^N / 2^N (1.17)，g(0) = 1，g'(0) = N/2 (1.18)，m^n 平均 = g^(n)(0) (1.19)；
g = e^{Ny/2} cosh(y/2)^N (1.21) → (1.24)；cosh 的冪級數給 (1.25)(1.26)；三種算法的表；比值 3 − 2/N (1.27)；
練習 1.3：Δ_iΔ_jΔ_kΔ_l 平均不為 0 只有「四個相同」與「兩兩成對」，N = 2：2 + 6 = 8。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
import itertools
import math
from fractions import Fraction

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

NS = list(range(1, 7))
M2 = [W.moment_exact(N, 2) for N in NS]
M4 = {"sum": [W.moment_exact(N, 4) for N in NS], "formula": [W.moment_formula(N, 4) for N in NS],
      "cosh": [W.moment_via_cosh(N, 4) for N in NS]}
RATIOS = [Fraction(W.moment_exact(N, 4)) / Fraction(W.moment_exact(N, 2)) ** 2 for N in range(1, 11)]
RATIO_TEX = ["1", "2", r"\tfrac73", r"\tfrac52", r"\tfrac{13}5", r"\tfrac83", r"\tfrac{19}7", r"\tfrac{11}4",
             r"\tfrac{25}9", r"\tfrac{14}5"]
TUPLES = list(itertools.product((1, 2), repeat=4))
ALL_EQ = [t for t in TUPLES if len(set(t)) == 1]
PAIRED = [t for t in TUPLES if len(set(t)) == 2 and all(t.count(v) == 2 for v in set(t))]


def cosh_half_series(N, deg=4):
    """cosh(y/2)^N 的冪級數係數（到 y^deg），用 Fraction 相乘。"""
    base = [Fraction(1, math.factorial(k) * 2 ** k) if k % 2 == 0 else Fraction(0) for k in range(deg + 1)]
    out = [Fraction(1)] + [Fraction(0)] * deg
    for _ in range(N):
        out = [sum(out[i] * base[k - i] for i in range(k + 1)) for k in range(deg + 1)]
    return out


def check_math():
    assert M2 == [1, 2, 3, 4, 5, 6] == [W.moment_formula(N, 2) for N in NS] == [W.moment_via_cosh(N, 2) for N in NS]
    assert M4["sum"] == M4["formula"] == M4["cosh"] == [1, 8, 21, 40, 65, 96]
    assert all(W.moment_formula(N, 4) == N * (3 * N - 2) for N in range(1, 11))               # (1.26)
    assert W.moment_exact(10, 2) == 10 and W.moment_exact(10, 4) == 280
    assert W.mean_heads(10) == 5 and W.mean_heads(2) == 1                                     # (1.18)
    c = [math.comb(2, m) for m in range(3)]
    assert c == [1, 2, 1] and Fraction(sum(c), 4) == 1 and Fraction(sum(m * x for m, x in enumerate(c)), 4) == 1
    for N in range(1, 11):
        s = cosh_half_series(N)
        assert s[2] == Fraction(N, 8) and s[4] == Fraction(N * (3 * N - 2), 384)
        assert 4 * math.factorial(2) * s[2] == N == W.moment_exact(N, 2)                     # (1.25)
        assert 16 * math.factorial(4) * s[4] == N * (3 * N - 2) == W.moment_exact(N, 4)      # (1.26)
    assert RATIOS == [Fraction(1), Fraction(2), Fraction(7, 3), Fraction(5, 2), Fraction(13, 5), Fraction(8, 3),
                      Fraction(19, 7), Fraction(11, 4), Fraction(25, 9), Fraction(14, 5)]
    assert all(r == 3 - Fraction(2, N) for N, r in zip(range(1, 11), RATIOS))               # (1.27)
    assert len(TUPLES) == 16 and len(ALL_EQ) == 2 and len(PAIRED) == 6
    assert len(ALL_EQ) + len(PAIRED) == W.fourth_moment_pairings(2) == 8 == W.moment_exact(2, 4)
    assert all(W.fourth_moment_pairings(N) == N + 3 * N * (N - 1) == 3 * N * N - 2 * N for N in range(1, 11))
    walks2 = W.all_walks(2)
    for t in TUPLES:
        avg = Fraction(sum(w[t[0] - 1] * w[t[1] - 1] * w[t[2] - 1] * w[t[3] - 1] for w in walks2), len(walks2))
        assert avg == (1 if (t in ALL_EQ or t in PAIRED) else 0)


def eq_row(tex, num, font_size=36):
    row = VGroup(Tex(tex, font_size=font_size))
    if num:
        row.add(Tex(num, font_size=24).set_color(GREY_B))
    return row.arrange(RIGHT, buff=0.35)


# ---------------------------------------------------------------- 1. 問題
def ask(self):
    big = Tex(r"\overline{x^4}=\overline{(\Delta_1+\Delta_2+\cdots+\Delta_N)^4}", font_size=52).move_to(0.6 * UP)
    note = mixed("展開有 ", "$N^4$", " 項", font_size=30, color=GREY_B).next_to(big, DOWN, buff=0.5)
    self.say("不展開 ", "$(\\sum\\Delta)^4$", "，能不能直接算出 ", "$\\overline{x^4}$", "？")
    self.play(Write(big), run_time=1.5)
    self.play(FadeIn(note))
    self.wait(3.5)
    self.play(FadeOut(big), FadeOut(note))


# ---------------------------------------------------------------- 2. 生成函數 g(y)
def gen_function(self):
    l1 = eq_row(r"\overline{x^n}=\frac{1}{2^N}\sum_{m=0}^{N}\big((2m-N)\,l\big)^n\frac{N!}{m!\,(N-m)!}", "(1.15)",
                34).move_to(2.45 * UP)
    l2 = eq_row(r"g(y)=\frac{1}{2^N}(1+e^y)^N=\frac{1}{2^N}\sum_{m=0}^{N}\binom{N}{m}e^{my}", "(1.17)", 34)
    l2.move_to(1.05 * UP)
    facts = VGroup(Tex(r"g(0)=1", font_size=34),
                   eq_row(r"g'(0)=\bar m=\frac{N}{2}", "(1.18)", 34),
                   eq_row(r"\overline{m^n}=g^{(n)}(0)", "(1.19)", 34)).arrange(RIGHT, buff=0.9).move_to(0.35 * DOWN)
    ex = Tex(r"N=2:\quad g(y)=\frac{1+2e^y+e^{2y}}{4},\quad g(0)=\frac{1+2+1}{4}=1,\quad g'(0)=\frac{2+2}{4}=1",
             font_size=30).move_to(1.75 * DOWN)
    assert_in_frame(l1, l2, facts, ex)
    for a, b in ((l1, l2), (l2, facts), (facts, ex)):
        assert_apart(a, b)
    self.say("$m$", " 是正面數，", "$x=(2m-N)\\,l$", "：直接加總（式 1.15）")
    self.play(Write(l1), run_time=1.8)
    self.wait(3.5)
    self.say("生成函數：把每一項的 ", "$m$", " 放進 ", "$e^{my}$", " 的指數")
    self.play(Write(l2), run_time=1.8)
    self.wait(3.5)
    self.say("對 ", "$y$", " 微分一次，就乘出一個 ", "$m$")
    self.play(FadeIn(facts[0]))
    self.wait(1.5)
    self.play(FadeIn(facts[1]))
    self.wait(2.0)
    self.play(FadeIn(facts[2]))
    self.wait(3.0)
    self.say("小例子 ", "$N=2$", "：平均正面數 ", "$\\bar m=1=N/2$")
    self.play(Write(ex), run_time=1.8)
    self.wait(4.0)


# ---------------------------------------------------------------- 3. cosh 與 (1.24)
def cosh_part(self):
    a1 = Tex(r"1+e^y=e^{y/2}\big(e^{-y/2}+e^{y/2}\big)=2\,e^{y/2}\cosh(y/2)", font_size=34).move_to(2.5 * UP)
    a2 = eq_row(r"g(y)=e^{Ny/2}\cosh(y/2)^N", "(1.21)", 38).move_to(1.4 * UP)
    a3 = eq_row(r"\overline{x^n}=(-2l)^n\frac{d^n}{dy^n}\cosh(y/2)^N\Big|_{y=0}", "(1.24)", 40).move_to(0.05 * UP)
    box = SurroundingRectangle(a3[0], buff=0.14).set_stroke(YELLOW, 2)
    ser = Tex(r"\cosh(y/2)^N=1+\frac{N}{8}\,y^2+\frac{N(3N-2)}{384}\,y^4+\cdots", font_size=34).move_to(1.45 * DOWN)
    assert_in_frame(a1, a2, a3, box, ser)
    for a, b in ((a1, a2), (a2, box), (box, ser)):
        assert_apart(a, b)
    self.say("關鍵：", "$1+e^y$", " 可以寫成 cosh")
    self.play(Write(a1), run_time=1.8)
    self.wait(3.0)
    self.play(Write(a2), run_time=1.4)
    self.wait(2.5)
    self.say("拿掉 ", "$e^{Ny/2}$", "：每微分一次，乘出的是 ", "$m-N/2=x/2l$")
    self.wait(1.0)
    self.play(Write(a3), run_time=1.8)
    self.play(ShowCreation(box))
    self.wait(3.5)
    self.say("cosh 的冪級數一展開，所有動差一次到手")
    self.play(Write(ser), run_time=1.8)
    self.wait(3.5)
    return VGroup(a1, a2, a3, box, ser)


def results(self, old):
    r1 = eq_row(r"\overline{x^2}=(-2)^2\cdot2!\cdot\frac{N}{8}\,l^2=N\,l^2", "(1.25)", 38)
    r2 = eq_row(r"\overline{x^4}=(-2)^4\cdot4!\cdot\frac{N(3N-2)}{384}\,l^4=l^4N(3N-2)", "(1.26)", 38)
    grp = VGroup(r1, r2).arrange(DOWN, buff=0.6).move_to(0.6 * UP)
    assert_in_frame(grp)
    self.say("係數乘上 ", "$(-2l)^n\\,n!$", "：")
    self.play(FadeOut(old[:4]), old[4].animate.move_to(2.5 * UP))
    self.play(Write(r1), run_time=1.6)
    self.wait(2.5)
    self.play(Write(r2), run_time=1.8)
    self.wait(3.5)
    assert_apart(old[4], grp)


# ---------------------------------------------------------------- 4. 三種算法的表
def table_part(self):
    hdr = [Tex("N", font_size=28), Tex(r"\overline{x^2}", font_size=28)]
    for name in ("加總", "公式", "cosh"):
        hdr.append(mixed("$\\overline{x^4}$", name, font_size=20))
    for m in hdr:
        m.set_color(YELLOW)
    rows = [hdr] + [[str(N), str(M2[i]), str(M4["sum"][i]), str(M4["formula"][i]), str(M4["cosh"][i])]
                    for i, N in enumerate(NS)]
    rows.append(["10", str(W.moment_exact(10, 2)), str(W.moment_exact(10, 4)), str(W.moment_formula(10, 4)),
                 str(W.moment_via_cosh(10, 4))])
    tab = D.tex_table(rows, [0.8, 1.1, 1.6, 1.6, 1.6], row_h=0.46, font_size=28)
    tab.move_to(np.array([-2.2, 0.2, 0]))
    last = VGroup(*[D.table_cell(tab, 7, j, 5) for j in range(5)])
    hl = SurroundingRectangle(last, buff=0.06).set_stroke(D.WALK_COLOR, 2)
    side = VGroup(zh("三種算法：", font_size=26),
                  mixed("直接加總 ", "$(1.15)$", font_size=24, color=GREY_A),
                  mixed("公式 ", "$(1.26)$", font_size=24, color=GREY_A),
                  mixed("cosh 的冪級數 ", "$(1.24)$", font_size=24, color=GREY_A),
                  zh("完全一致", font_size=26, color=D.EXACT_COLOR)).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    side.move_to(np.array([0, 0.5, 0])).align_to(np.array([4.0, 0, 0]), LEFT)
    assert_in_frame(tab, hl, side)
    assert_apart(VGroup(tab, hl), side)
    self.say("表（", "$l=1$", "）：", "$N=1$", " 到 ", "$6$", "，再加上主例 ", "$N=10$")
    self.play(FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:5])))
    for i in range(1, 8):
        self.play(FadeIn(VGroup(*tab[0][5 * i:5 * i + 5]), shift=0.1 * UP), run_time=0.45)
    self.play(FadeIn(side, lag_ratio=0.2), run_time=1.5)
    self.wait(3.0)
    self.say("主例 ", "$N=10$", "：", "$\\overline{x^2}=10$", "，", "$\\overline{x^4}=280$")
    self.play(ShowCreation(hl))
    self.wait(3.5)


# ---------------------------------------------------------------- 5. 比值 → 3；練習 1.3
def ratio_part(self):
    ax = D.plot_axes((0, 11, 1), (0, 3.6, 1), width=7.2, height=3.6, x_ticks=list(range(1, 11)), y_ticks=[1, 2, 3],
                     x_label="N")
    ax.move_to(np.array([-2.6, 0.1, 0]))
    dots = VGroup(*[Dot(ax.c2p(N, float(r)), radius=0.07).set_fill(YELLOW, 1) for N, r in zip(range(1, 11), RATIOS)])
    labs = VGroup(*[Tex(t, font_size=26).set_color(YELLOW).next_to(d, DOWN, buff=0.15)
                    for t, d in zip(RATIO_TEX, dots)])
    three = DashedLine(ax.c2p(0, 3), ax.c2p(10.8, 3), dash_length=0.1).set_stroke(D.GAUSS_COLOR, 2)
    f = Tex(r"\frac{\overline{x^4}}{\big(\overline{x^2}\big)^2}=3-\frac{2}{N}", font_size=42)
    f.move_to(np.array([4.4, 1.6, 0]))
    n27 = Tex("(1.27)", font_size=24).set_color(GREY_B).next_to(f, DOWN, buff=0.2)
    g1 = mixed("$N\\to\\infty$", "：趨近 3", font_size=26, color=D.GAUSS_COLOR).move_to(np.array([4.4, -0.2, 0]))
    g2 = zh("3 是高斯分布的值", font_size=26, color=D.GAUSS_COLOR).move_to(np.array([4.4, -0.75, 0]))
    assert_in_frame(ax, dots, labs, three, f, n27, g1, g2)
    assert_apart(VGroup(ax, labs, three), VGroup(f, n27, g1, g2))
    for a, b in zip(labs, labs[1:]):
        assert_apart(a, b, gap=0.02)
    self.say("比一比：", "$\\overline{x^4}$", " 和 ", "$(\\overline{x^2})^2$", " 的比值")
    self.play(FadeIn(ax), Write(f), FadeIn(n27), run_time=1.5)
    self.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in dots], lag_ratio=0.15),
              LaggedStart(*[FadeIn(t) for t in labs], lag_ratio=0.15), run_time=2.5)
    self.wait(2.5)
    self.say("點列一路爬向虛線 3", wait=0.5)
    self.play(ShowCreation(three), run_time=1.2)
    self.play(FadeIn(g1), FadeIn(g2))
    self.wait(4.0)


def pairing_part(self):
    top = Tex(r"\overline{\Delta_i\Delta_j\Delta_k\Delta_l}\neq0", font_size=40).move_to(np.array([-3.8, 2.55, 0]))
    cells = VGroup()
    for idx, t in enumerate(TUPLES):
        col = D.EXACT_COLOR if t in ALL_EQ else (D.HEAD_COLOR if t in PAIRED else GREY_D)
        cells.add(Tex("".join(map(str, t)), font_size=32).set_color(col)
                  .move_to(np.array([-5.4 + 1.15 * (idx % 4), 1.55 - 0.62 * (idx // 4), 0])))
    lines = VGroup(mixed("四個相同：", "$N$", " 組", font_size=28, color=D.EXACT_COLOR),
                   mixed("兩兩成對：", "$3N(N-1)$", " 組", font_size=28, color=D.HEAD_COLOR),
                   Tex(r"N+3N(N-1)=3N^2-2N", font_size=36),
                   mixed("$N=2$", "：", "$2+6=8=\\overline{x^4}$", font_size=28, color=YELLOW))
    lines.arrange(DOWN, buff=0.35, aligned_edge=LEFT).move_to(np.array([0, 0.9, 0])).align_to(np.array([0.4, 0, 0]),
                                                                                               LEFT)
    cap = zh("i, j, k, l 各取 1 或 2：共 16 組", font_size=22, color=GREY_B).next_to(cells, DOWN, buff=0.35)
    assert_in_frame(top, cells, lines, cap)
    assert_apart(cells, lines)
    assert_apart(top, cells)
    self.say("練習 1.3：直接數，哪些 ", "$\\Delta_i\\Delta_j\\Delta_k\\Delta_l$", " 平均不為 0？")
    self.play(Write(top))
    self.play(FadeIn(cells, lag_ratio=0.05), FadeIn(cap), run_time=1.5)
    self.wait(2.0)
    self.say("灰色：有一個指標只出現奇數次，平均是 0", wait=3.0)
    self.play(*[Indicate(cells[TUPLES.index(t)]) for t in ALL_EQ], FadeIn(lines[0]))
    self.wait(2.5)
    self.play(*[Indicate(cells[TUPLES.index(t)]) for t in PAIRED], FadeIn(lines[1]))
    self.wait(2.5)
    self.say("合計 ", "$3N^2-2N$", "：和 (1.26) 一樣")
    self.play(Write(lines[2]))
    self.play(FadeIn(lines[3]))
    self.wait(4.0)


class W105MomentGenerating(ChuScene):
    def construct(self):
        self.show_header("W105　用生成函數算所有動差")
        check_math()
        ask(self)
        gen_function(self)
        self.clear_stage()
        old = cosh_part(self)
        results(self, old)
        self.clear_stage()
        table_part(self)
        self.clear_stage()
        ratio_part(self)
        self.clear_stage()
        pairing_part(self)
        self.say("生成函數把「加總」變成「微分」，動差一次算完", wait=4.0)
        self.wait(1.0)
