"""
W202 生成函數本身就是答案（§2.1.1–2.1.2，pp. 27–29，(2.2)(2.3)，Ex. 2.5，Fig. 2.1）

問：走路的人每走一步有 10% 的機會消失（存活率 z = 0.9）。平均能走幾步？最後停在原點的機率多大？
  1. 統計力學 (2.2)：Q(z) = Σ Z(N) z^N（fugacity z = e^{−α}）：巨配分函數＝正則配分函數的生成函數；
     練習 2.5 理想氣體 Q = exp(z𝒞T^{3/2}V)；小例子 x = 2：Z(N) = 1, 2, 2, 4/3, 2/3, 4/15，z = 1 的粒子數機率（Poisson）。
  2. 會消失的走路人：剛好 n 步的機率 z^n(1−z)，平均 z/(1−z) = 9；模擬 2000 人（seed 2）對照理論，樣本平均 9.34。
  3. (2.3) G(z; x, y) = Σ z^N C(N; x, y)；主例 1 + 4z² + 36z⁴ + …；一個會消失的人的足跡（seed 25：24 步，停在 (3, −3)）。
  4. 揭曉：停在原點的機率 (1−z) G(z/4; 0, 0) = 0.1452；模擬 296/2000 = 0.148。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import boxed, mixed, zh
from utils import lattice_draw as L
from utils import walk_draw as D
from utils import walk_gf as GF

Z = 0.9
X_GAS = 2
M_WALKERS, SEED_WALKERS = 2000, 2
SEED_TRAIL = 25
SIM_COLOR = D.BAR_COLOR
THEO_COLOR = D.GAUSS_COLOR


def frac_tex(c, small=False) -> str:
    c = Fraction(c)
    if c.denominator == 1:
        return str(c.numerator)
    cmd = r"\tfrac" if small else r"\frac"
    return ("-" if c < 0 else "") + rf"{cmd}{{{abs(c.numerator)}}}{{{c.denominator}}}"


def sim():
    w = GF.mortal_walkers(M_WALKERS, Z, seed=SEED_WALKERS)
    cnt = Counter(n for n, _ in w)
    return [cnt[n] for n in range(10)], sum(n for n, _ in w) / len(w), sum(1 for _, e in w if e == (0, 0))


def numbers():
    hist, mean, home = sim()
    return {
        "Z": [GF.ideal_gas_Z(N, X_GAS) for N in range(6)],
        "p": [f"{p:.4f}" for p in GF.grand_probs(X_GAS, 1.0, 5)],
        "pf": GF.grand_probs(X_GAS, 1.0, 5),
        "Q": f"{GF.grand_Q(1.0, X_GAS):.3f}",
        "surv": [f"{GF.survival_prob(n, Z):.3g}" for n in range(3)],
        "mean_th": f"{GF.survival_mean(Z):.0f}",
        "hist": hist, "theo": [round(M_WALKERS * GF.survival_prob(n, Z)) for n in range(10)],
        "mean": f"{mean:.2f}", "home": home, "home_frac": f"{home / M_WALKERS:.3f}",
        "exact": f"{GF.mortal_origin_exact(Z):.4f}",
        "step": f"{Z / 4:.3f}", "die": f"{1 - Z:.1f}",
        "returns": [GF.returns_sq(N) for N in range(0, 11, 2)],
    }


def check_math():
    nb = numbers()
    F = Fraction
    assert nb["Z"] == [1, 2, 2, F(4, 3), F(2, 3), F(4, 15)]
    assert nb["p"] == ["0.1353", "0.2707", "0.2707", "0.1804", "0.0902", "0.0361"]
    assert nb["Q"] == "7.389" and abs(GF.grand_Q(1.0, X_GAS) - math.exp(2)) < 1e-12
    assert abs(sum(N * p for N, p in enumerate(GF.grand_probs(X_GAS, 1.0, 60))) - 1.0 * X_GAS) < 1e-9   # 平均 zx
    assert nb["surv"] == ["0.1", "0.09", "0.081"] and nb["mean_th"] == "9"
    assert nb["hist"] == [214, 167, 156, 138, 119, 123, 98, 94, 88, 71]
    assert nb["theo"] == [200, 180, 162, 146, 131, 118, 106, 96, 86, 77]
    assert nb["mean"] == f"{GF.MORTAL_MEAN:.2f}" == "9.34"
    assert nb["home"] == GF.MORTAL_AT_ORIGIN == 296 and nb["home_frac"] == "0.148"
    assert nb["exact"] == "0.1452" and abs(GF.mortal_origin_exact(Z) - GF.mortal_origin_series(Z)) < 1e-6
    assert nb["step"] == "0.225" and nb["die"] == "0.1"
    assert nb["returns"] == [1, 4, 36, 400, 4900, 63504]
    t = GF.mortal_trail(Z, seed=SEED_TRAIL)
    assert len(t) - 1 == GF.MORTAL_TRAIL_LEN == 24 and t[0] == (0, 0) and t[-1] == (3, -3)
    assert all(abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1 for a, b in zip(t, t[1:]))
    assert max(max(abs(x), abs(y)) for x, y in t) <= 4


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


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    org = np.array([0, 0.75, 0])
    cell = 1.4
    c = Dot(org, radius=0.13).set_fill(YELLOW, 1)
    arrows = VGroup(*[Arrow(org, L.lat_point(x, y, cell, org), buff=0.15, thickness=3).set_color(GREY_B)
                      for x, y in GF.SQ_STEPS])
    labs = VGroup(*[Tex(nb["step"], font_size=28).next_to(L.lat_point(x, y, cell, org), np.array([x, y, 0]), buff=0.12)
                    for x, y in GF.SQ_STEPS])
    gone = VGroup(Tex(r"\times", font_size=48).set_color(RED), zh("消失 " + nb["die"], font_size=26, color=RED))
    gone.arrange(RIGHT, buff=0.15).move_to(org + np.array([2.9, 1.35, 0]))
    cap = mixed("每一步：四個方向各 ", f"$0.9/4={nb['step']}$", "，消失 ", f"$1-z={nb['die']}$", font_size=28)
    cap.move_to(np.array([0, -1.75, 0]))
    assert_in_frame(arrows, labs, gone, cap)
    assert_apart(labs, cap)
    assert_apart(gone, labs[0])

    s.play(FadeIn(c, scale=0.5), *[GrowArrow(a) for a in arrows], run_time=1.0)
    s.say("問：走路的人每走一步，有 10% 的機會消失", wait=0.3)
    s.play(FadeIn(labs), FadeIn(gone), FadeIn(cap), run_time=0.8)
    s.wait(3.0)
    s.say("平均能走幾步？最後停在原點的機率多大？", wait=3.5)
    s.clear_stage()


def part_statmech(s):
    nb = numbers()
    s1 = eq(r"Q(\alpha)=\sum_{N}Z(N)\,e^{-\alpha N}", "(2.2)", font_size=34)
    s2 = mixed("$z=e^{-\\alpha}$", "（fugacity，", "$\\alpha=-\\mu/kT$", "）", font_size=24)
    s3 = Tex(r"Q(z)=\sum_{N}Z(N)\,z^N", font_size=34)
    s4 = zh("巨配分函數＝正則配分函數的生成函數", font_size=24, color=YELLOW)
    s5 = mixed("練習 2.5：", "$Z(N)=\\frac{1}{N!}\\big(\\mathcal{C}T^{3/2}V\\big)^N$", font_size=24)
    s6 = Tex(r"\Rightarrow\ Q(z)=\exp\big(z\,\mathcal{C}T^{3/2}V\big)", font_size=34)
    col = VGroup(s1, s2, s3, s4, s5, s6)
    for m in col:
        fit(m, 6.4)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.35), -6.5, 2.95)

    ttl = mixed("小例子：", "$x=\\mathcal{C}T^{3/2}V=2$", "，", "$z=1$", font_size=24)
    rows = [["N"] + [str(N) for N in range(6)],
            ["Z(N)"] + [frac_tex(v, small=True) for v in nb["Z"]],
            [r"Z(N)z^N/Q"] + nb["p"]]
    widths = [1.55] + [0.78] * 6
    tab = D.tex_table(rows, widths, row_h=0.55, font_size=22)
    left = 0.55
    place_col(ttl, left, 2.95)
    tab.shift(np.array([left, ttl.get_bottom()[1] - 0.2, 0]))
    base_y = -1.2
    org = np.array([left + widths[0] + widths[1] / 2, base_y, 0])
    bars = D.bar_chart({N: p for N, p in enumerate(nb["pf"])}, ux=0.78, uy=5.0, origin=org, width=0.7)
    for b in bars[0]:
        b.set_fill(SIM_COLOR, 0.75)
    nlabs = VGroup(*[Tex(str(N), font_size=20).set_color(GREY_B).next_to(D.bar_of(bars, N), DOWN, buff=0.1)
                     for N in range(6)])
    qline = mixed(f"$Q=e^2\\approx{nb['Q']}$", "，平均 ", "$\\langle N\\rangle=zx=2$", font_size=24)
    qline.move_to(np.array([left + 3.0, -2.15, 0]))
    assert_in_frame(col, ttl, tab, bars, nlabs, qline)
    assert_apart(col, tab)
    assert_apart(col, bars)
    assert_apart(tab, bars)
    assert_apart(nlabs, qline)

    s.say("先看統計力學：巨配分函數 (2.2)", wait=0.3)
    s.play(Write(s1[0]), FadeIn(s1[1]), run_time=1.5)
    s.wait(2.0)
    s.play(FadeIn(s2), FadeIn(s3), run_time=1.0)
    s.say("令 ", "$z=e^{-\\alpha}$", "（fugacity）：", "$Q$", " 變成 ", "$z$", " 的冪級數", wait=3.0)
    s.play(FadeIn(s4), run_time=0.6)
    s.say("巨配分函數就是正則配分函數的生成函數", wait=3.5)
    s.play(FadeIn(s5), FadeIn(s6), run_time=1.0)
    s.say("練習 2.5：理想氣體的 ", "$Q$", " 是一個指數函數", wait=3.0)
    s.play(FadeIn(ttl), FadeIn(tab), run_time=1.0)
    s.say("小例子 ", "$x=2$", "：", "$Z(N)=2^N/N!$", wait=3.0)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars[0]], lag_ratio=0.15), FadeIn(nlabs), run_time=1.5)
    s.say("$z=1$", "：粒子數 ", "$N$", " 的機率是 ", "$Z(N)z^N/Q$", wait=3.0)
    s.play(FadeIn(qline), run_time=0.8)
    s.say("這是 Poisson 分布，平均 ", "$\\langle N\\rangle=zx=2$", wait=3.5)
    s.clear_stage()


def part_mortal(s):
    nb = numbers()
    sv = nb["surv"]
    f1 = mixed("剛好走 ", "$n$", " 步就消失：", f"$z^n(1-z)={sv[0]},\\ {sv[1]},\\ {sv[2]},\\dots$", font_size=26)
    f2 = mixed("平均 ", f"$\\langle n\\rangle=\\frac{{z}}{{1-z}}={nb['mean_th']}$", font_size=26)
    top = VGroup(f1, f2).arrange(RIGHT, buff=0.8).move_to(np.array([0, 2.6, 0]))
    ux, uy, base_y = 0.68, 2.3 / 214, -1.15
    org = np.array([-3.4, base_y, 0])
    bars = D.bar_chart({n: v for n, v in enumerate(nb["hist"])}, ux=ux, uy=uy, origin=org, width=0.72)
    dots = VGroup(*[Dot(org + n * ux * RIGHT + v * uy * UP, radius=0.07).set_fill(THEO_COLOR, 1)
                    for n, v in enumerate(nb["theo"])])

    def num_row(vals, dy, color, size=22):
        return VGroup(*[Tex(str(v), font_size=size).set_color(color).move_to(org + n * ux * RIGHT + dy * DOWN)
                        for n, v in enumerate(vals)])

    nrow = num_row(range(10), 0.25, GREY_B)
    srow = num_row(nb["hist"], 0.6, SIM_COLOR)
    trow = num_row(nb["theo"], 0.95, THEO_COLOR)
    rl = VGroup(Tex("n", font_size=22).set_color(GREY_B).move_to(org + 0.95 * LEFT + 0.25 * DOWN),
                zh("模擬", font_size=18, color=SIM_COLOR).move_to(org + 0.95 * LEFT + 0.6 * DOWN),
                zh("理論", font_size=18, color=THEO_COLOR).move_to(org + 0.95 * LEFT + 0.95 * DOWN))
    leg = VGroup(
        VGroup(Square(0.22).set_fill(SIM_COLOR, 0.75).set_stroke(SIM_COLOR, 1), zh("模擬 2000 人", font_size=22)
               ).arrange(RIGHT, buff=0.15),
        VGroup(Dot(radius=0.07).set_fill(THEO_COLOR, 1), mixed("理論 ", "$2000\\cdot0.1\\cdot0.9^n$", font_size=22)
               ).arrange(RIGHT, buff=0.15)).arrange(DOWN, aligned_edge=LEFT, buff=0.25)
    place_col(fit(leg, 3.0), 3.75, 1.2)
    mean = mixed("樣本平均 ", f"${nb['mean']}$", "（理論 ", f"${nb['mean_th']}$", "）", font_size=26)
    mean.move_to(np.array([0, -2.68, 0]))
    assert_in_frame(top, bars, dots, nrow, srow, trow, rl, leg, mean)
    assert_apart(top, bars)
    assert_apart(leg, bars)
    assert_apart(nrow, srow)
    assert_apart(srow, trow)
    assert_apart(trow, mean)

    s.say("回到會消失的走路人（p. 28）", wait=0.3)
    s.play(FadeIn(f1), run_time=0.8)
    s.say("剛好走 ", "$n$", " 步就消失的機率：", "$z^n(1-z)$", wait=3.0)
    s.play(FadeIn(f2), run_time=0.8)
    s.say("平均步數 ", "$\\frac{z}{1-z}=9$", wait=2.5)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars[0]], lag_ratio=0.1), FadeIn(nrow), FadeIn(rl[0]),
           run_time=1.8)
    s.play(FadeIn(srow), FadeIn(rl[1]), FadeIn(leg[0]), run_time=0.8)
    s.say("模擬 2000 人：每一種步數各有幾個人", wait=3.0)
    s.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in dots], lag_ratio=0.1), FadeIn(trow), FadeIn(rl[2]),
           FadeIn(leg[1]), run_time=1.5)
    s.say("橘點是理論值 ", "$2000\\cdot0.1\\cdot0.9^n$", wait=3.0)
    s.play(FadeIn(mean), run_time=0.8)
    s.say("樣本平均 9.34，理論 9：吻合", wait=3.5)
    s.clear_stage()


def part_lattice(s):
    nb = numbers()
    rt = nb["returns"]
    c1 = mixed("$C(N;\\vec x,\\vec y)$", "：從 ", "$\\vec x$", " 走 ", "$N$", " 步到 ", "$\\vec y$", " 的走法數",
               font_size=24)
    c2 = eq(r"G(z;\vec x,\vec y)=\sum_{N\ge0}z^N\,C(N;\vec x,\vec y)", "(2.3)", font_size=34)
    terms = [f"{v}z^{{{2 * k}}}" if k > 1 else (f"{v}z^2" if k == 1 else str(v)) for k, v in enumerate(rt)]
    c3a = Tex(r"G(z;\vec 0,\vec 0)=" + "+".join(terms[:4]), font_size=32)
    c3b = Tex("+" + "+".join(terms[4:]) + r"+\cdots", font_size=32)
    c4 = mixed("係數＝回到原點的走法數 ", "$C(N;\\vec 0)$", font_size=22, color=GREY_B)
    col = VGroup(c1, c2, c3a, c3b, c4).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    c3b.shift((c3a.get_right()[0] - c3b.get_right()[0]) * RIGHT)
    for m in col:
        fit(m, 6.6)
    place_col(col, -6.5, 2.7)

    cell, org = 0.55, np.array([3.95, 0.25, 0])
    grid = L.lattice_lines(4, cell, org)
    pts = GF.mortal_trail(Z, seed=SEED_TRAIL)
    tr = L.trail(pts, cell, org, stroke_width=4, end_radius=0.09)
    tl = mixed("走了 24 步後消失，停在 ", "$(3,-3)$", font_size=24)
    tl.next_to(grid, DOWN, buff=0.15)
    ring = Circle(radius=0.2).set_stroke(YELLOW, 3).move_to(L.lat_point(0, 0, cell, org))

    nb_e = nb["exact"]
    r0 = zh("停在原點的機率：", font_size=26)
    r1 = Tex(r"\sum_{N\ge0}z^N(1-z)\,\frac{C(N;\vec 0)}{4^N}", font_size=36)
    r2 = Tex(r"=(1-z)\,G\Big(\frac{z}{4};\vec 0,\vec 0\Big)", font_size=36)
    r3 = boxed(Tex(rf"z=0.9:\ \ {nb_e}", font_size=40), color=D.EXACT_COLOR, buff=0.18)
    r4 = mixed("模擬：", f"${nb['home']}/{M_WALKERS}={nb['home_frac']}$", font_size=28)
    rcol = VGroup(r0, r1, r2, r3, r4).arrange(DOWN, aligned_edge=LEFT, buff=0.38)
    place_col(rcol, -6.5, 2.8)
    assert_in_frame(col, grid, tr, tl, rcol)
    assert_apart(col, grid)
    assert_apart(rcol, grid)

    s.say("回到格子：(2.3) 走法數的生成函數", wait=0.3)
    s.play(FadeIn(c1), run_time=0.8)
    s.play(Write(c2[0]), FadeIn(c2[1]), run_time=1.5)
    s.wait(2.5)
    s.play(FadeIn(c3a), FadeIn(c3b), FadeIn(c4), run_time=1.0)
    s.say("主例：回到原點 ", "$1,4,36,400,4900,63504$", wait=3.0)
    s.play(FadeIn(grid), FadeIn(tr[1]), run_time=0.8)
    s.say("一個會消失的走路人（Fig. 2.1 的樣子）", wait=0.3)
    s.play(ShowCreation(tr[0]), run_time=3.5, rate_func=linear)
    s.play(FadeIn(tr[2], scale=0.5), FadeIn(tl), run_time=0.8)
    s.say("走了 24 步才消失，停在 ", "$(3,-3)$", "，不是原點", wait=3.0)
    s.play(FadeOut(col), ShowCreation(ring), run_time=0.8)
    s.say("揭曉：停在原點的機率是多少？", wait=0.3)
    s.play(FadeIn(r0), Write(r1), run_time=1.8)
    s.say("走 ", "$N$", " 步才消失，乘上走 ", "$N$", " 步回到原點的機率", wait=3.0)
    s.play(FadeIn(r2), run_time=1.0)
    s.say("加起來正好是 ", "$G$", "：只是把 ", "$z$", " 換成 ", "$z/4$", wait=3.0)
    s.play(FadeIn(r3, shift=0.2 * UP), run_time=0.8)
    s.say("$z=0.9$", " 代進去：", f"${nb_e}$", wait=3.0)
    s.play(FadeIn(r4), run_time=0.8)
    s.say("模擬 2000 人：296 人停在原點，", f"${nb['home_frac']}$", wait=3.0)
    s.say("把 ", "$z$", " 換成存活率：生成函數本身就是答案", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("$z$", " 是「每多一個（粒子、單體、一步）就打的折扣」", wait=3.5)
    s.say("生成函數＝所有大小加權的總和，本身就有物理意義", wait=4.0)


class W202GFIsAnswer(ChuScene):
    def construct(self):
        self.show_header("W202　生成函數本身就是答案")
        check_math()
        part_question(self)
        part_statmech(self)
        part_mortal(self)
        part_lattice(self)
        part_close(self)
        self.wait(1.0)
