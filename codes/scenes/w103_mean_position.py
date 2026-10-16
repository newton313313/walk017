"""
W103 平均位置是 0（§1.2，pp. 5–6，式 (1.8)–(1.9)）

鎖步系綜：ensemble(20, 10, seed=1) 的 20 人同時走 10 步，終點堆在數線上（w(x)：2,2,3,6,3,1,3）；
(1.8) = (1.9)：總和 2、平均 0.1；左右配對抵消得不完美；對全部 1024 種結果的「超級平均」= 0；
平方平均 252/20 = 63/5 = 12.6（精確值 10，下一幕）。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
from fractions import Fraction

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

M, N = 20, 10
WALKS = W.ensemble(M, N, seed=1)
PATHS = [W.positions(w) for w in WALKS]
FINALS = W.finals(WALKS)
HIST = W.histogram(FINALS)
STEP_MEANS = [W.mean([p[k] for p in PATHS]) for k in range(N + 1)]
DIST = W.distribution(N)
SUMSQ = sum(x * x for x in FINALS)
PAIRS = [(-6, 6), (-4, 4), (-2, 2)]

UNIT = 0.45
NL_Y = -1.05
ROW_Y0, ROW_DY = 2.95, 0.17
STACK_DY = 0.27
COUNT_Y = 0.95
LINE_A, LINE_B = 2.6, 1.75
DOT_COLOR = TEAL_B


def check_math():
    assert len(FINALS) == M and FINALS == [2, 0, 6, 4, -4, -2, -2, 0, -4, -6, 0, -2, 2, 6, 2, 0, 0, 0, -6, 6]
    assert HIST == {-6: 2, -4: 2, -2: 3, 0: 6, 2: 3, 4: 1, 6: 3} and sum(HIST.values()) == M
    assert sum(FINALS) == sum(x * w for x, w in HIST.items()) == 2                 # (1.8) = (1.9)
    assert W.mean(FINALS) == Fraction(2, 20) == Fraction(1, 10) == STEP_MEANS[-1]
    assert [a * HIST[a] + b * HIST[b] for a, b in PAIRS] == [6, -4, 0]
    assert sum(d * c for d, c in DIST.items()) == 0 and all(DIST[d] == DIST[-d] for d in DIST)
    assert SUMSQ == sum(x * x * w for x, w in HIST.items()) == 252
    assert [x * x for x in sorted(HIST)] == [36, 16, 4, 0, 4, 16, 36]
    assert [d * DIST[d] + (-d) * DIST[-d] for d in (10, 8, 6, 4, 2)] == [0] * 5
    assert [DIST[d] for d in (10, 8, 6, 4, 2)] == [1, 10, 45, 120, 210]
    assert W.mean_square(FINALS) == Fraction(252, 20) == Fraction(63, 5) and float(Fraction(63, 5)) == 12.6
    assert W.moment_exact(N, 2) == 10
    assert HIST[0] == 6 and M - HIST[0] == 14
    assert max(abs(x) for p in PATHS for x in p) <= N


def eq_row(tex, num, font_size=36):
    row = VGroup(Tex(tex, font_size=font_size))
    if num:
        row.add(Tex(num, font_size=24).set_color(GREY_B))
    return row.arrange(RIGHT, buff=0.35)


def swap(self, old, new, run_time=1.0):
    anims = [FadeIn(new, shift=0.15 * UP)]
    if old is not None:
        anims.insert(0, FadeOut(old, shift=0.15 * UP))
    self.play(*anims, run_time=run_time)
    return new


# ---------------------------------------------------------------- 1.–2. 鎖步系綜
def lockstep(self):
    nl = D.number_line(N, unit=UNIT).shift(NL_Y * UP)
    guides = VGroup(*[Line(D.line_point(nl, -N), D.line_point(nl, N)).set_stroke(GREY_D, 1, 0.6)
                      .set_y(ROW_Y0 - ROW_DY * i) for i in range(M)])
    dots = VGroup(*[Dot(radius=0.065).set_fill(DOT_COLOR, 1).move_to(g.get_center()) for g in guides])
    clock = zh("第 0 步", font_size=24, color=GREY_B).move_to(np.array([6.0, 2.85, 0]))
    mlab = Tex(r"\bar x=", font_size=30).move_to(np.array([5.75, 2.3, 0]))
    mval = DecimalNumber(0, num_decimal_places=2, font_size=30).next_to(mlab, RIGHT, buff=0.1)
    assert_in_frame(nl, guides, dots, clock, mlab, mval)
    assert_apart(guides, nl)
    assert_apart(guides, VGroup(clock, mlab, mval))
    self.say("平均來說，醉漢最後停在哪裡？", wait=1.0)
    self.play(ShowCreation(nl[0]), FadeIn(nl[1]), FadeIn(nl[2]))
    one = D.walker_dot().move_to(D.line_point(nl, 0))
    self.play(FadeIn(one, scale=0.5))
    for x in W.MAIN_POS[1:]:
        self.play(one.animate.move_to(D.line_point(nl, x)), run_time=0.25)
    self.say("一個人（上一幕的主例）停在 2：看不出「平均」", wait=2.5)
    self.say("找 ", "$M=20$", " 個人一起走", wait=0.3)
    self.play(FadeOut(one), FadeIn(guides, lag_ratio=0.03), run_time=1.2)
    self.play(FadeIn(dots, lag_ratio=0.03), FadeIn(clock), FadeIn(mlab), FadeIn(mval))
    self.say("鎖步系綜：每一拍，20 個人同時各走一步", wait=3.0)
    for k in range(1, N + 1):
        new_clock = zh(f"第 {k} 步", font_size=24, color=GREY_B).move_to(clock)
        self.play(*[d.animate.set_x(D.line_point(nl, p[k])[0]) for d, p in zip(dots, PATHS)],
                  Transform(clock, new_clock), run_time=0.9)
        mval.set_value(float(STEP_MEANS[k]))
        self.wait(0.4)
    self.say("右上角是每一拍的平均位置：一直在 0 附近", wait=3.5)
    return nl, guides, dots, VGroup(clock, mlab, mval)


def stack_up(self, nl, guides, dots, readout):
    rank, targets = {}, []
    for x in FINALS:
        r = rank.get(x, 0)
        rank[x] = r + 1
        targets.append(D.line_point(nl, x) + (0.2 + STACK_DY * r) * UP)
    counts = VGroup(*[Tex(str(HIST[x]), font_size=30).set_color(YELLOW).move_to(np.array([D.line_point(nl, x)[0],
                                                                                           COUNT_Y, 0]))
                      for x in sorted(HIST)])
    wlab = Tex("w(x)", font_size=30).set_color(YELLOW).move_to(np.array([-4.3, COUNT_Y, 0]))
    xlab = Tex("x", font_size=30).set_color(GREY_B).move_to(D.line_point(nl, -N) + np.array([-0.55, -0.32, 0]))
    self.say("把 20 個終點堆在數線上")
    self.play(FadeOut(guides), FadeOut(readout), *[d.animate.move_to(t) for d, t in zip(dots, targets)], run_time=1.8)
    stacks_top = max(t[1] for t in targets) + 0.065
    assert stacks_top + 0.05 < counts.get_bottom()[1]
    assert_in_frame(counts, wlab, xlab)
    self.play(FadeIn(counts, lag_ratio=0.1), FadeIn(wlab), FadeIn(xlab))
    self.wait(1.0)
    self.say("$w(x)$", "：停在 ", "$x$", " 的人數（加起來是 ", "$M=20$", "）", wait=3.5)
    return counts


# ---------------------------------------------------------------- 3. (1.8)＝(1.9)
def mean_formulas(self, dots):
    lst = Tex(r"x_i:\ 2,0,6,4,-4,-2,-2,0,-4,-6,0,-2,2,6,2,0,0,0,-6,6", font_size=30).move_to(LINE_A * UP)
    f18 = eq_row(r"\bar x=\frac{1}{M}\sum_{i=1}^{M}x_i=\frac{2}{20}=0.1", "(1.8)").move_to(LINE_B * UP)
    f19 = eq_row(r"\bar x=\frac{1}{M}\sum_x x\,w(x)", "(1.9)").move_to(LINE_A * UP)
    f3 = Tex(r"=\frac{1}{20}\big[(-6)\cdot2+(-4)\cdot2+(-2)\cdot3+0\cdot6+2\cdot3+4\cdot1+6\cdot3\big]=\frac{2}{20}",
             font_size=30).move_to(LINE_B * UP)
    for m in (lst, f18, f19, f3):
        assert_in_frame(m)
        assert m.get_bottom()[1] > COUNT_Y + 0.2
    self.say("第 ", "$i$", " 個人的終點叫 ", "$x_i$", "：一個一個加起來")
    a = swap(self, None, lst)
    self.wait(3.0)
    self.say("總和是 2，除以人數 20：平均 0.1")
    b = swap(self, None, f18, run_time=1.4)
    self.wait(3.5)
    self.say("也可以先按位置分組：每個位置乘上人數 ", "$w(x)$")
    a = swap(self, a, f19)
    self.wait(2.0)
    b = swap(self, b, f3, run_time=1.4)
    self.wait(4.0)
    return a, b


# ---------------------------------------------------------------- 4. 對稱與超級平均
def symmetry(self, nl, dots, counts, a, b):
    pair = Tex(r"\underbrace{(-6)\cdot2+6\cdot3}_{6}+\underbrace{(-4)\cdot2+4\cdot1}_{-4}"
               r"+\underbrace{(-2)\cdot3+2\cdot3}_{0}=2", font_size=32).move_to(LINE_B * UP)
    assert_in_frame(pair)
    assert pair.get_bottom()[1] > COUNT_Y + 0.15
    self.say("對稱：", "$x$", " 與 ", "$-x$", " 的人數「平均來說」一樣多")
    b = swap(self, b, pair, run_time=1.4)
    for lo, hi in PAIRS:
        sel = [d for d, x in zip(dots, FINALS) if x in (lo, hi)]
        self.play(*[Indicate(d, scale_factor=1.6) for d in sel], run_time=0.8)
    self.wait(2.0)
    self.say("20 個人的系綜太小：抵消得不完美，剩下 0.1", wait=3.5)
    self.say("換另外 20 個人重走，剩下的又會不一樣", wait=3.0)
    o = np.array([0.0, NL_Y + 0.065, 0])
    ideal = D.bar_chart({d: M * c / 2 ** N for d, c in DIST.items()}, ux=UNIT, uy=STACK_DY, origin=o, width=1.3,
                        color=D.BAR_COLOR, fill_opacity=0.3)
    assert ideal[0].get_top()[1] < counts.get_bottom()[1]
    sup = eq_row(r"\bar x=\frac{1}{1024}\sum_d d\,C(10,d)", "").move_to(LINE_A * UP)
    sym = Tex(r"=\frac{1}{1024}\big[(-10+10)\cdot1+(-8+8)\cdot10+(-6+6)\cdot45+(-4+4)\cdot120"
              r"+(-2+2)\cdot210\big]=0", font_size=28).move_to(LINE_B * UP)
    assert_in_frame(sup, sym)
    assert sym.get_bottom()[1] > COUNT_Y + 0.15
    self.say("「超級平均」：對全部 1024 種結果平均")
    self.add(ideal[0])
    self.bring_to_back(ideal[0])
    self.play(FadeIn(ideal[0], lag_ratio=0.05), run_time=1.5)
    self.say("藍色長條：", "$20\\cdot C(10,d)/1024$", "，左右完全對稱", wait=3.0)
    self.say("人數越多，", "$w(x)/M$", " 越接近 ", "$C(10,x)/1024$", wait=3.0)
    a = swap(self, a, sup)
    b = swap(self, b, sym, run_time=1.4)
    self.wait(2.0)
    self.say("$x\\,w(x)$", " 與 ", "$-x\\,w(-x)$", " 完全抵消：", "$\\bar x=0$", wait=3.5)
    return ideal, a, b


# ---------------------------------------------------------------- 5. 0 不代表留在原地
def not_home(self, nl, dots, a, b):
    away = [d for d, x in zip(dots, FINALS) if x != 0]
    home = [d for d, x in zip(dots, FINALS) if x == 0]
    self.say("但 0 不代表留在原地：只有 6 人停在 0")
    self.play(*[d.animate.set_fill(ORANGE, 1) for d in away], *[Indicate(d) for d in home], run_time=1.2)
    self.wait(1.0)
    self.say("其他 14 人都不在 0：要看離 0 多遠", wait=2.5)
    sq = VGroup(*[Tex(str(x * x), font_size=26).set_color(ORANGE)
                  .move_to(np.array([D.line_point(nl, x)[0], NL_Y - 0.8, 0])) for x in sorted(HIST)])
    sqlab = Tex("x^2", font_size=28).set_color(ORANGE).move_to(np.array([-4.3, NL_Y - 0.8, 0]))
    assert_in_frame(sq, sqlab)
    assert_apart(sq, nl)
    self.say("每個人離 0 的距離，取平方：", "$x^2$")
    self.play(FadeIn(sq, lag_ratio=0.1), FadeIn(sqlab))
    self.wait(3.0)
    ms_a = eq_row(r"\overline{x^2}=\frac{1}{M}\sum_x x^2\,w(x)", "").move_to(LINE_A * UP)
    ms_b = Tex(r"=\frac{36\cdot2+16\cdot2+4\cdot3+0\cdot6+4\cdot3+16\cdot1+36\cdot3}{20}=\frac{252}{20}"
               r"=\frac{63}{5}=12.6", font_size=30).move_to(LINE_B * UP)
    assert_in_frame(ms_a, ms_b)
    assert ms_b.get_bottom()[1] > COUNT_Y + 0.15
    self.say("平方的平均：負號不會再互相抵消")
    a = swap(self, a, ms_a)
    b = swap(self, b, ms_b, run_time=1.5)
    self.wait(3.5)
    self.say("這 20 人是 12.6；精確值 ", "$\\overline{x^2}=10$", "（下一幕）", wait=4.0)


class W103MeanPosition(ChuScene):
    def construct(self):
        self.show_header("W103　平均位置是 0")
        check_math()
        nl, guides, dots, readout = lockstep(self)
        counts = stack_up(self, nl, guides, dots, readout)
        a, b = mean_formulas(self, dots)
        ideal, a, b = symmetry(self, nl, dots, counts, a, b)
        not_home(self, nl, dots, a, b)
        self.say("平均位置是 0，但大多數人都不在 0", wait=3.5)
        self.say("要看的是：平方的平均", wait=3.5)
        self.wait(1.0)
