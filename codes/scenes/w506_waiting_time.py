"""
W506 連續時間的漫步：等一下再走（§5.3，pp. 118–122，(5.107)–(5.123)，Ex. 5.8）

問：走路人每落腳一次就要「等公車」才走下一步。3 分鐘內走了幾步？固定速率和「越等越可能來」哪個走得多？
  1. (5.107)(5.108)：a(τ) = exp(−∫q dt)（V.survival = V.survival_numeric）；a(3) = 0.0498、0.1991；(5.109) ∫aq dt = 1。
  2. (5.110)–(5.121)：走 m 步＝卷積 → Laplace 後 A^m B，A = 1 − λB；λ = 0.7（V.AB = V.AB_numeric）。
  3. t = 3 走了 m 步的機率（V.steps_prob = _laplace = _conv）；平均 3 與 1.2506；(5.122) Σ A^m B = 1/λ。
  4. Ex. 5.8：1/(λ + p0(1 − cos q)) → e^{−p0 t} I_|x|(p0 t)（V.ctrw_bessel = V.ctrw_lattice）。
畫面上的數字都由 utils.walk_var 算出，check_math() 核對。
"""
import math
import sys
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

KINDS = ("poisson", "gamma2")
KCOL = {"poisson": BLUE_C, "gamma2": ORANGE}
KNAME = {"poisson": "固定速率", "gamma2": "兩段等待"}
SEEDS = {"poisson": (3, 6), "gamma2": (21, 27)}
SEED_COUNTS = {"poisson": (3, 5), "gamma2": (1, 2)}
LAM = 0.7
MLAB = {"poisson": 7, "gamma2": 5}          # 標數字的 m = 0..MLAB−1（卡上的範圍）
MBAR = {"poisson": 9, "gamma2": 7}          # 畫出的長條 m = 0..MBAR−1
XS = range(-6, 7)


def f4(x):
    return f"{x:.4f}"


# ---------------------------------------------------------------- 數學
def jump_times(seed, kind, t=3.0):
    """示意：一條時間線上的出發時刻（固定速率：等待 ~ Exp(1)；兩段等待：兩個 Exp(1) 相加，密度 te^{−t}）。"""
    rng = np.random.default_rng(seed)
    ts, s = [], 0.0
    while True:
        w = rng.exponential(1.0) if kind == "poisson" else rng.exponential(1.0) + rng.exponential(1.0)
        s += w
        if s > t:
            return ts
        ts.append(s)


def trapz(f, x):
    return float(np.sum((f[1:] + f[:-1]) * np.diff(x)) / 2)


@lru_cache(maxsize=None)
def numbers():
    nb = {"a3": {k: V.survival(3, k) for k in KINDS}, "a3n": {k: V.survival_numeric(3, k) for k in KINDS}}
    ts = np.linspace(0, 60, 60001)
    nb["norm"] = {k: trapz(np.array([V.wait_density(x, k) for x in ts]), ts) for k in KINDS}
    nb["B0"] = {k: V.mean_wait(k) for k in KINDS}
    nb["AB"] = {k: V.AB(LAM, k) for k in KINDS}
    nb["ABn"] = {k: V.AB_numeric(LAM, k) for k in KINDS}
    nb["P"] = {k: [V.steps_prob(m, 3.0, k) for m in range(MBAR[k])] for k in KINDS}
    nb["mean"] = {k: V.mean_steps(3.0, k) for k in KINDS}
    nb["bes"] = {x: V.ctrw_bessel(x, 3.0) for x in XS}
    return nb


def check_math():
    nb = numbers()
    assert (f4(nb["a3"]["poisson"]), f4(nb["a3"]["gamma2"])) == ("0.0498", "0.1991")
    assert all(abs(nb["a3"][k] - nb["a3n"][k]) < 1e-9 for k in KINDS)
    assert abs(nb["a3"]["poisson"] - math.exp(-3)) < 1e-15 and abs(nb["a3"]["gamma2"] - 4 * math.exp(-3)) < 1e-15
    assert all(abs(nb["norm"][k] - 1) < 1e-6 for k in KINDS)                          # (5.109)
    assert (nb["B0"]["poisson"], nb["B0"]["gamma2"]) == (1, 2)
    assert [f4(x) for x in nb["AB"]["poisson"]] == ["0.5882", "0.5882"]
    assert [f4(x) for x in nb["AB"]["gamma2"]] == ["0.3460", "0.9343"]
    for k in KINDS:
        A, B = nb["AB"][k]
        assert abs(A - nb["ABn"][k][0]) < 1e-6 and abs(B - nb["ABn"][k][1]) < 1e-6
        assert abs(A - (1 - LAM * B)) < 1e-14                                          # (5.119)
        assert abs(B / (1 - A) - 1 / LAM) < 1e-12                                      # (5.122)
    assert f4(1 - LAM * nb["AB"]["gamma2"][1]) == "0.3460"
    assert [f4(x) for x in nb["P"]["poisson"][:7]] == ["0.0498", "0.1494", "0.2240", "0.2240", "0.1680", "0.1008", "0.0504"]
    assert [f4(x) for x in nb["P"]["gamma2"][:5]] == ["0.1991", "0.4481", "0.2689", "0.0720", "0.0108"]
    for k in KINDS:
        for m in range(MLAB[k]):
            p = nb["P"][k][m]
            assert abs(p - V.steps_prob_laplace(m, 3.0, k)) < 1e-9 and abs(p - V.steps_prob_conv(m, 3.0, k)) < 2e-6
    assert (f4(nb["mean"]["poisson"]), f4(nb["mean"]["gamma2"])) == ("3.0000", "1.2506")
    assert [f4(nb["bes"][x]) for x in (0, 1, 2)] == ["0.2430", "0.1968", "0.1118"]
    assert all(abs(nb["bes"][x] - nb["bes"][-x]) < 1e-15 and abs(nb["bes"][x] - V.ctrw_lattice(x, 3.0)) < 1e-13 for x in XS)
    assert abs(V.ctrw_transform(LAM, 1.1) - V.ex58(LAM, 1.1)) < 1e-14
    for k in KINDS:
        assert tuple(len(jump_times(sd, k)) for sd in SEEDS[k]) == SEED_COUNTS[k]


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def eq_label(tex_mob, label, buff=0.4):
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


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def timeline(ts, y, x0=-4.0, x1=4.0, t=3.0, color=YELLOW):
    """一條時間線：灰線＋起點＋每次出發一個點＋右邊「k 步」。"""
    sx = (x1 - x0) / t
    line = Line(np.array([x0, y, 0]), np.array([x1, y, 0])).set_stroke(GREY_C, 2)
    start = Dot(np.array([x0, y, 0]), radius=0.07).set_color(WHITE)
    dots = VGroup(*[Dot(np.array([x0 + s * sx, y, 0]), radius=0.1).set_color(color) for s in ts])
    cnt = mixed("$%d$" % len(ts), " 步", font_size=26, color=color).next_to(line, RIGHT, buff=0.35)
    return VGroup(line, start, dots, cnt)


def labelled_bars(values, nlab, ux, uy, origin, color):
    """長條（m = 0..）＋前 nlab 根標上數字＋底下的 m。"""
    chart = D.bar_chart({m: v for m, v in enumerate(values)}, ux=ux, uy=uy, origin=origin, color=color, width=0.75)
    labs = VGroup(*[tx(f4(values[m]), 18).next_to(chart[0][m], UP, buff=0.07) for m in range(nlab)])
    ms = VGroup(*[tx(str(m), 20, GREY_B).move_to(np.array(origin) + m * ux * RIGHT + 0.25 * DOWN)
                  for m in range(len(values))])
    base = Line(np.array(origin) - 0.6 * ux * RIGHT, np.array(origin) + (len(values) - 0.4) * ux * RIGHT)
    base.set_stroke(GREY_B, 2)
    return VGroup(chart[0], labs, ms, base)


# ---------------------------------------------------------------- 各段
def part_question(s):
    rows = []
    ys = {"poisson": (1.9, 1.15), "gamma2": (-0.45, -1.2)}
    for k in KINDS:
        for sd, y in zip(SEEDS[k], ys[k]):
            rows.append(timeline(jump_times(sd, k), y, color=KCOL[k]))
    names = VGroup(*[zh(KNAME[k], font_size=28, color=KCOL[k]).move_to(np.array([-5.5, sum(ys[k]) / 2, 0]))
                     for k in KINDS])
    axis = VGroup(*[tx(str(t), 22, GREY_B).move_to(np.array([-4.0 + t * 8 / 3, -1.85, 0])) for t in range(4)],
                  zh("分鐘", font_size=22, color=GREY_B).move_to(np.array([4.75, -1.85, 0])))
    notes = VGroup(zh("公車隨時可能來（和等了多久無關）", font_size=22, color=GREY_A),
                   zh("要等兩件事：越等越可能來", font_size=22, color=GREY_A))
    notes[0].move_to(np.array([0, 2.55, 0]))
    notes[1].move_to(np.array([0, 0.2, 0]))
    allm = VGroup(*rows, names, axis, notes)
    assert_in_frame(allm)
    for a, b in ((names[0], rows[0]), (names[1], rows[2]), (notes[0], rows[0]), (notes[1], rows[1]), (notes[1], rows[2]),
                 (rows[3], axis)):
        assert_apart(a, b)
    s.say("問：走路人每落腳一次，就要「等公車」才走下一步", wait=1.5)
    s.play(*[ShowCreation(r[0]) for r in rows], *[FadeIn(r[1]) for r in rows], FadeIn(axis), FadeIn(names), run_time=1.0)
    s.play(FadeIn(notes[0]), LaggedStart(*[FadeIn(d, scale=0.5) for r in rows[:2] for d in r[2]], lag_ratio=0.25),
           run_time=2.0)
    s.play(*[FadeIn(r[3]) for r in rows[:2]], run_time=0.5)
    s.say("固定速率：公車隨時可能來，和已經等了多久無關", wait=2.5)
    s.play(FadeIn(notes[1]), LaggedStart(*[FadeIn(d, scale=0.5) for r in rows[2:] for d in r[2]], lag_ratio=0.25),
           run_time=1.5)
    s.play(*[FadeIn(r[3]) for r in rows[2:]], run_time=0.5)
    s.say("兩段等待：要等兩件事都發生，越等越可能來", wait=2.5)
    s.say("3 分鐘內平均走幾步？哪一種走得多？（先猜）", wait=3.5)
    s.clear_stage()


def part_survival(s):
    nb = numbers()
    ax = D.plot_axes((0, 6, 1), (0, 1.25, 0.25), width=6.0, height=3.9, x_ticks=[1, 2, 3, 4, 5, 6],
                     y_ticks=[0.5, 1.0], font_size=20, x_label="t")
    ax.move_to(np.array([-3.4, -0.5, 0]))
    qs = VGroup(*[DashedVMobject(ax.get_graph(lambda t, k=k: V.hazard(t, k), x_range=(0, 6, 0.05)), num_dashes=40)
                  .set_stroke(KCOL[k], 3) for k in KINDS])
    as_ = VGroup(*[ax.get_graph(lambda t, k=k: V.survival(t, k), x_range=(0, 6, 0.05)).set_stroke(KCOL[k], 4)
                   for k in KINDS])
    legend = zh("虛線：q(t)　實線：a(t)", font_size=22, color=GREY_A).move_to(ax.c2p(4.0, 1.18))
    dots = VGroup(*[Dot(ax.c2p(3, nb["a3"][k]), radius=0.07).set_color(KCOL[k]) for k in KINDS])
    e1 = eq_label(tx(r"a(\tau)=\exp\Big(-\int_0^\tau q(t)\,dt\Big)", 32), "(5.107)(5.108)")
    hz = lambda t, c=YELLOW: zh(t, font_size=24, color=c)
    rows = [[hz("等待"), "q(t)", "a(t)", "a(3)"],
            [hz(KNAME["poisson"], KCOL["poisson"]), "p_0=1", "e^{-t}", f4(nb["a3"]["poisson"])],
            [hz(KNAME["gamma2"], KCOL["gamma2"]), r"\frac{t}{1+t}", "(1+t)e^{-t}", f4(nb["a3"]["gamma2"])]]
    tab = D.tex_table(rows, [1.6, 1.4, 1.8, 1.3], row_h=0.75, font_size=28)
    e2 = eq_label(tx(r"\int_0^\infty a\,q\,dt=1-a(\infty)=1", 30), "(5.109)")
    e3 = mixed("平均等待 ", "$B(0)=\\int_0^\\infty a\\,dt$", "：1 與 2", font_size=26)
    right = VGroup(e1, tab, e2, e3).arrange(DOWN, buff=0.45)
    fit(right, 6.3)
    at(right, 3.5, 0.1)
    left = VGroup(ax, qs, as_, legend, dots)
    assert_in_frame(left, right)
    assert_apart(left, right)
    assert_column(*right)
    assert_table(tab, 4)
    s.play(ShowCreation(ax), FadeIn(e1), run_time=1.0)
    s.say("$q(t)\\,dt$", "：等了 ", "$t$", " 還沒走，接下來 ", "$dt$", " 內出發的機率", wait=2.5)
    s.play(*[ShowCreation(q) for q in qs], FadeIn(legend), run_time=1.2)
    s.say("固定速率 ", "$q=p_0=1$", "；兩段等待 ", "$q=\\frac{t}{1+t}$", "：越等越大", wait=2.5)
    s.play(*[ShowCreation(a) for a in as_], run_time=1.2)
    s.say("$a(\\tau)$", "：等了 ", "$\\tau$", " 還沒出發的機率", wait=2.0)
    s.say("兩段等待一開始很少出發，所以「還在等」的比較多", wait=2.5)
    s.play(FadeIn(tab[0][:4]), ShowCreation(tab[1]), FadeIn(tab[0][4:]), FadeIn(dots), run_time=1.0)
    s.say("$t=3$", "：", "$a(3)=e^{-3}=0.0498$", "；", "$4e^{-3}=0.1991$", wait=3.0)
    s.play(FadeIn(e2), run_time=0.6)
    s.say("$aq$", " 是出發時刻的密度：積分是 1，一定會走", wait=2.5)
    s.play(FadeIn(e3), run_time=0.6)
    s.say("平均等待：固定速率 1，兩段等待 2", wait=2.5)
    s.clear_stage()


def part_laplace(s):
    nb = numbers()
    c1 = eq_label(tx(r"P_0(t)=a(t),\qquad P_m(t)=\int_0^t a(\tau)q(\tau)\,P_{m-1}(t-\tau)\,d\tau", 32),
                  "(5.110)(5.111)")
    c2 = eq_label(tx(r"\tilde P_m(\lambda)=A^mB,\quad A=\int_0^\infty e^{-\lambda t}a\,q\,dt,\quad"
                     r" B=\int_0^\infty e^{-\lambda t}a\,dt", 32), "(5.112)--(5.121)")
    c3 = VGroup(eq_label(tx(r"A=1-\lambda B", 36, YELLOW), "(5.119)"),
                mixed("（分部積分：", "$aq=-a'$", "）", font_size=26, color=GREY_A)).arrange(RIGHT, buff=0.5)
    hz = lambda t, c=YELLOW: zh(t, font_size=24, color=c)
    (Ap, Bp), (Ag, Bg) = nb["AB"]["poisson"], nb["AB"]["gamma2"]
    rows = [[hz("等待"), r"A(\lambda)", r"B(\lambda)", "A(0.7)", "B(0.7)"],
            [hz(KNAME["poisson"], KCOL["poisson"]), r"\frac{p_0}{p_0+\lambda}", r"\frac{1}{p_0+\lambda}", f4(Ap), f4(Bp)],
            [hz(KNAME["gamma2"], KCOL["gamma2"]), r"\frac{1}{(1+\lambda)^2}", r"\frac{2+\lambda}{(1+\lambda)^2}",
             f4(Ag), f4(Bg)]]
    tab = D.tex_table(rows, [1.7, 2.2, 2.2, 1.5, 1.5], row_h=0.8, font_size=28)
    chk = tx(r"1-0.7\times%s=%s" % (f4(Bg), f4(1 - LAM * Bg)), 28, GREY_A)
    col = VGroup(c1, c2, c3, VGroup(tab, chk).arrange(RIGHT, buff=0.6)).arrange(DOWN, buff=0.45)
    fit(col, 13.0)
    at(col, 0, 0.1)
    assert_in_frame(col)
    assert_column(*col)
    assert_table(tab, 5)
    s.play(FadeIn(c1), run_time=1.0)
    s.say("走 ", "$m$", " 步：先等、走、再等……每一段接起來＝卷積", wait=3.0)
    s.play(FadeIn(c2), run_time=1.0)
    s.say("Laplace 轉換把卷積變成乘法：", "$\\tilde P_m=A^mB$", wait=2.5)
    s.say("$A$", "：等完一段、出發；", "$B$", "：最後一段還在等", wait=2.5)
    s.play(FadeIn(c3), run_time=0.8)
    s.say("分部積分（", "$aq=-a'$", "）：", "$A=1-\\lambda B$", wait=2.5)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), FadeIn(tab[0][5:]), run_time=1.0)
    s.say("$\\lambda=0.7$", "：固定速率 ", "$A=B=0.5882$", wait=2.5)
    s.play(FadeIn(chk), run_time=0.6)
    s.say("兩段等待 ", "$A=0.3460$", "、", "$B=0.9343$", "（＝數值積分）", wait=2.5)
    s.say("代回核對 ", "$A=1-\\lambda B$", "：兩種等待都成立", wait=2.5)
    s.clear_stage()


def part_steps(s):
    nb = numbers()
    ux, uy, y0 = 0.62, 4.6, -1.95
    org = {"poisson": np.array([-6.1, y0, 0]), "gamma2": np.array([1.3, y0, 0])}
    charts, titles, means = {}, {}, {}
    for k in KINDS:
        charts[k] = labelled_bars(nb["P"][k], MLAB[k], ux, uy, org[k], KCOL[k])
        mtxt = "3" if k == "poisson" else f4(nb["mean"][k])
        titles[k] = mixed(KNAME[k], "：平均 ", "$%s$" % mtxt, " 步", font_size=28, color=KCOL[k])
        titles[k].move_to(org[k] + (MBAR[k] - 1) / 2 * ux * RIGHT + 3.3 * UP)
        xm = org[k] + nb["mean"][k] * ux * RIGHT
        means[k] = Arrow(xm + 0.9 * DOWN, xm + 0.5 * DOWN, buff=0, thickness=3.0).set_color(WHITE)
    e = eq_label(tx(r"\sum_m A^mB=\frac{B}{1-A}=\frac{1}{\lambda}", 32), "(5.122)")
    at(e, 0, 2.55)
    mlab = tx("m", 24, GREY_B).next_to(charts["gamma2"][3], RIGHT, buff=0.2)
    assert_in_frame(e, mlab, *charts.values(), *titles.values(), *means.values())
    assert_apart(VGroup(charts["poisson"], titles["poisson"], means["poisson"]),
                 VGroup(charts["gamma2"], titles["gamma2"], means["gamma2"]))
    for k in KINDS:
        assert_apart(titles[k], charts[k])
        assert_apart(means[k], charts[k][2], 0.02)
        assert_apart(means[k], charts[k][0], 0.02)
        assert_apart(e, titles[k])
        labs = charts[k][1]
        for a, b in zip(labs, labs[1:]):
            assert_apart(a, b, 0.02)
    for k in KINDS:
        c = charts[k]
        s.play(FadeIn(c[3]), FadeIn(c[2]), FadeIn(titles[k][0]), run_time=0.6)
        s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in c[0]], lag_ratio=0.12), run_time=1.5)
        s.play(FadeIn(c[1]), run_time=0.6)
        if k == "poisson":
            s.say("$t=3$", " 時走了 ", "$m$", " 步的機率（三種算法一致）", wait=2.5)
            s.say("固定速率：Poisson 分布，", "$m=2,3$", " 最可能", wait=2.0)
        else:
            s.say("兩段等待：", "$m=1$", " 最可能", wait=2.0)
    s.play(FadeIn(mlab), *[GrowArrow(means[k]) for k in KINDS], *[FadeIn(titles[k][1:]) for k in KINDS], run_time=0.8)
    s.say("揭曉：平均 3 步對 1.2506 步（箭頭），固定速率走得多", wait=3.0)
    s.say("兩段等待平均要等 2，步數大約只有一半", wait=2.5)
    s.play(FadeIn(e), run_time=0.8)
    s.say("(5.122)：所有 ", "$m$", " 加起來 ", "$=\\frac1\\lambda$", "，反轉後機率和是 1", wait=3.0)
    s.clear_stage()


def part_lattice(s):
    nb = numbers()
    g1 = eq_label(tx(r"\sum_m\chi^mA^mB=\frac{B}{1-\chi A}=\frac{1}{\lambda+p_0(1-\cos q)}", 34), "Ex. 5.8")
    g2 = tx(r"P(x,t)=e^{-p_0t}I_{|x|}(p_0t)", 36, YELLOW)
    top = VGroup(g1, g2).arrange(DOWN, buff=0.35)
    at(top, 0, 2.3)
    ux, uy, org = 0.8, 8.0, np.array([0, -2.15, 0])
    chart = D.bar_chart(nb["bes"], ux=ux, uy=uy, origin=org, color=TEAL_C, width=0.75)
    labs = VGroup(*[tx(f4(nb["bes"][x]), 20).next_to(D.bar_of(chart, x), UP, buff=0.07) for x in (-2, -1, 0, 1, 2)])
    xl = VGroup(*[tx(str(x), 20, GREY_B).move_to(org + x * ux * RIGHT + 0.25 * DOWN) for x in XS])
    base = Line(org + (XS[0] - 0.6) * ux * RIGHT, org + (XS[-1] + 0.6) * ux * RIGHT).set_stroke(GREY_B, 2)
    side = VGroup(mixed("$t=3$", "，", "$p_0=1$", font_size=28),
                  zh("每一格都有人", font_size=26, color=TEAL_B)).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    at(side, -4.6, 0.62)
    assert_in_frame(top, chart, labs, xl, side)
    assert_apart(top, VGroup(chart, labs))
    assert_apart(side, chart)
    assert_apart(side, labs)
    assert_apart(side, top)
    for a, b in zip(labs, labs[1:]):
        assert_apart(a, b, 0.02)
    s.play(FadeIn(g1), run_time=1.0)
    s.say("Ex. 5.8：每次出發再往左或右走一格（", "$\\chi=\\cos q$", "）", wait=2.5)
    s.say("固定速率代入 ", "$A,B$", "：得到 ", "$\\frac{1}{\\lambda+p_0(1-\\cos q)}$", wait=2.5)
    s.play(FadeIn(g2), run_time=0.8)
    s.say("反轉：修正 Bessel 函數 ", "$e^{-p_0t}I_{|x|}(p_0t)$", wait=2.5)
    s.play(ShowCreation(base), FadeIn(xl), FadeIn(side[0]), run_time=0.6)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in chart[0]], lag_ratio=0.08), run_time=1.5)
    s.play(FadeIn(labs), run_time=0.6)
    s.say("$x=0,\\pm1,\\pm2$", "：0.2430、0.1968、0.1118", wait=2.5)
    s.play(FadeIn(side[1]), run_time=0.6)
    s.say("對照：固定走 ", "$N$", " 步，只會落在和 ", "$N$", " 同奇偶的格子", wait=2.5)
    s.say("步數不固定，所以每一格都可能有人（沒有奇偶性）", wait=3.0)


def part_close(s):
    s.say("等待時間的卷積，Laplace 轉換後只是乘法", wait=3.0)
    s.say("走 ", "$m$", " 步的轉換就是 ", "$A^mB$", wait=3.0)


class W506WaitingTime(ChuScene):
    def construct(self):
        self.show_header("W506　連續時間的漫步：等一下再走")
        check_math()
        part_question(self)
        part_survival(self)
        part_laplace(self)
        part_steps(self)
        part_lattice(self)
        part_close(self)
        self.wait(1.0)
