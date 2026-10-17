"""
W504 電報方程：兩個尖峰＋拖在後面的分布（§5.2.2，pp. 103–109，(5.39)–(5.64)，Fig. 5.2–5.3，Ex. 5.3–5.6）

問：走路人以速率 v = 1 前進，平均每單位時間反轉 γ = 1 次。t = 3 時有沒有人在 x = 3 之外？一次都沒轉身的有多少？
  1. (5.39)–(5.48)：R、L 隨波前進＋以 γ 互換；Σ = R + L、Δ = R − L；Ex. 5.4 電報方程；Ex. 5.3 格子極限 0.3392、0.3380 → 0.3368。
  2. 解 (5.61)（更正後）＋筆誤框 (5.60)(5.61)：書上 v = 2、t = 3 積分得 1.9502。
  3. Fig. 5.3：t = 1, 3, 5, 8；尖峰 e^{−t}/2；守恆。
  4. 大 t（Ex. 5.5–5.6）：高斯 D = v²/2γ、變異數 → 2Dt − 0.5、l_p = v/γ。
  5. 揭曉：10 萬人模擬（V.rt_mc(1, 3.0, 100000, 5)）：|x| < 1.5 的比例 0.6060（公式 0.6051），沒轉身 0.0488（e^{−3} = 0.0498）。
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

TS = (1, 3, 5, 8)
T_COL = {1: YELLOW, 3: GREEN_C, 5: TEAL_C, 8: PURPLE_B}
SPIKE_H = {1: 0.37, 3: 0.17, 5: 0.12, 8: 0.08}       # 尖峰箭頭的高度（示意；重量寫在旁邊）
DTS = (0.01, 0.005)
T_MC, M_MC, SEED_MC = 3.0, 100000, 5
LANES, LANE_SEED, LANE_DT = 16, 3, 0.25
R_COL = D.HEAD_COLOR
L_COL = D.TAIL_COLOR


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def g4(x):
    return f"{float(x):.4g}"


def lane_walkers(n=LANES, T=T_MC, dt=LANE_DT, seed=LANE_SEED):
    """示意用的 n 個走路人：每個人的位置（每 dt 一格）與「從沒轉身」。"""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        d = 1 if rng.random() < 0.5 else -1
        flips, t = [], rng.exponential(1.0)
        while t < T:
            flips.append(t)
            t += rng.exponential(1.0)

        def pos(tau, d=d, flips=flips):
            x, cur, dd = 0.0, 0.0, d
            for f in flips:
                if f >= tau:
                    break
                x += dd * (f - cur)
                cur, dd = f, -dd
            return x + dd * (tau - cur)
        out.append(([pos(k * dt) for k in range(int(round(T / dt)) + 1)], not flips))
    return out


def tele_residual(x0=0.7, t0=2.0, h=1e-3):
    """Ex. 5.4：Σ_tt + 2γΣ_t − v²Σ_xx（差分，v = γ = 1）。"""
    S = lambda x, t: float(V.tele_cont(x, t))
    return ((S(x0, t0 + h) - 2 * S(x0, t0) + S(x0, t0 - h)) + h * (S(x0, t0 + h) - S(x0, t0 - h))
            - (S(x0 + h, t0) - 2 * S(x0, t0) + S(x0 - h, t0))) / h ** 2


@lru_cache(maxsize=None)
def numbers():
    nb = {}
    lat = []
    for dt in DTS:
        xl, Pl = V.tele_lattice(1.0, dt)
        lat.append(Pl[int(np.argmin(np.abs(xl)))] / (2 * dt))
    nb["lat"] = (lat[0], lat[1], 2 * lat[1] - lat[0], float(V.tele_cont(0, 1)))
    nb["book"] = V.tele_mass(3, v=2.0, book=True)
    nb["fig"] = {t: dict(c0=float(V.tele_cont(0, t)), front=float(V.tele_cont(t * (1 - 1e-12), t)),
                         spike=V.tele_spike(t), cont=1 - math.exp(-t), mass=V.tele_mass(t),
                         gauss=float(V.tele_gauss(0, t)), var=V.tele_var(t)) for t in TS}
    x1 = V.rt_mc(1, T_MC, M_MC, SEED_MC)
    cnt, edges = np.histogram(x1, bins=30, range=(-3, 3))
    nb["mc"] = dict(frac=V.mc_fraction(x1, 0, 1.5), formula=V.radial_fraction(1, 0, 1.5, T_MC),
                    unturned=float(np.mean(np.abs(x1) > T_MC)), e3=math.exp(-T_MC),
                    hist=cnt / (M_MC * (edges[1] - edges[0])), edges=edges, n_in=int(cnt.sum()))
    nb["lanes"] = lane_walkers()
    return nb


def check_math():
    nb = numbers()
    a, b, rich, formula = nb["lat"]
    assert (f4(a), f4(b), f4(rich), f4(formula)) == ("0.3392", "0.3380", "0.3368", "0.3368")
    assert abs(rich - formula) < 1e-4
    assert abs(tele_residual()) < 1e-6                                                     # Ex. 5.4
    assert f4(nb["book"]) == "1.9502" and abs(V.tele_mass(3, v=2.0) - 1) < 1e-9            # 書上少 1/v
    fig = nb["fig"]
    assert [f4(fig[t]["c0"]) for t in TS] == ["0.3368", "0.2199", "0.1738", "0.1388"]
    assert [g4(fig[t]["front"]) for t in TS] == ["0.2759", "0.06223", "0.01179", "0.0008387"]
    assert [g4(fig[t]["spike"]) for t in TS] == ["0.1839", "0.02489", "0.003369", "0.0001677"]
    assert [f4(fig[t]["cont"]) for t in TS] == ["0.6321", "0.9502", "0.9933", "0.9997"]
    for t in TS:
        assert abs(fig[t]["mass"] - 1) < 1e-9                                               # 守恆
        assert abs(fig[t]["mass"] - 2 * fig[t]["spike"] - fig[t]["cont"]) < 1e-9
        assert float(V.tele_cont(t * 1.0001 + 0.01, t)) == 0                                # 因果性
        assert fig[t]["c0"] < 0.4 and SPIKE_H[t] > fig[t]["front"]
    assert [f4(fig[t]["gauss"]) for t in TS] == ["0.3989", "0.2303", "0.1784", "0.1410"]
    assert all(abs(fig[t]["gauss"] - 1 / math.sqrt(4 * math.pi * V.tele_D() * t)) < 1e-12 for t in TS)
    assert [f4(fig[t]["var"]) for t in TS] == ["0.5677", "2.5012", "4.5000", "7.5000"]
    assert all(abs(V.tele_var(t) - V.tele_var_numeric(t)) < 1e-8 for t in TS)
    assert V.tele_D() == 0.5 and V.tele_lp() == 1.0 and V.tele_D() == 0.5 * 1.0 * V.tele_lp()
    assert [f"{2 * V.tele_D() * t - 0.5:.1f}" for t in TS] == ["0.5", "2.5", "4.5", "7.5"]
    mc = nb["mc"]
    assert (f4(mc["frac"]), f4(mc["formula"]), f4(mc["unturned"]), f4(mc["e3"])) == \
           ("0.6060", "0.6051", "0.0488", "0.0498")
    assert len(mc["hist"]) <= 30 and mc["n_in"] == round(M_MC * (1 - mc["unturned"]))
    assert max(mc["hist"]) < 0.25 and round(100 * mc["e3"]) == 5
    lanes = nb["lanes"]
    assert sum(u for _, u in lanes) == 1 and all(max(abs(x) for x in p) <= T_MC + 1e-9 for p, _ in lanes)
    assert all(abs(abs(p[-1]) - T_MC) < 1e-9 for p, u in lanes if u)


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def eq_label(tex_mob, label, buff=0.35):
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


def cont_curve(ax, t, color, n=120, stroke=3):
    xs = np.linspace(0, t * (1 - 1e-9), n)
    ys = V.tele_cont(xs, t)
    c = VMobject().set_points_smoothly([ax.c2p(x, y) for x, y in zip(xs, ys)]).set_stroke(color, stroke)
    drop = Line(ax.c2p(t, ys[-1]), ax.c2p(t, 0)).set_stroke(color, 1.5, 0.6)
    return VGroup(c, drop)


def gauss_curve(ax, t, color, x1=9.0, n=90):
    xs = np.linspace(0, x1, n)
    pts = [ax.c2p(x, float(V.tele_gauss(x, t))) for x in xs]
    c = VMobject().set_points_smoothly(pts).set_stroke(color, 2.5)
    return DashedVMobject(c, num_dashes=30).set_stroke(color, 2.5)


def spike_arrow(ax, x, h, color, label, font_size=20, side=UP):
    ar = Arrow(ax.c2p(x, 0), ax.c2p(x, h), buff=0, thickness=3).set_color(color)
    lab = tx(label, font_size, color).next_to(ar.get_end(), side, buff=0.08)
    return VGroup(ar, lab)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    lanes = nb["lanes"]
    u, ytop, dy = 1.45, 1.7, 0.22
    xp = lambda x, i: np.array([x * u, ytop - i * dy, 0])
    base = Line(xp(-3.6, LANES - 1) + 0.3 * DOWN, xp(3.6, LANES - 1) + 0.3 * DOWN).set_stroke(GREY_C, 2)
    nums = VGroup(*[tx(str(x), 22, GREY_B).move_to(xp(x, LANES - 1) + 0.6 * DOWN) for x in range(-3, 4)])
    walls = VGroup(*[DashedLine(xp(x, 0) + 0.3 * UP, xp(x, LANES - 1) + 0.3 * DOWN).set_stroke(RED, 2)
                     for x in (-3, 3)])
    dots = VGroup(*[Dot(xp(0, i), radius=0.07).set_fill(WHITE, 1) for i in range(LANES)])
    clock = at(tx(r"t=0", 30, GREY_A), -5.6, 2.3)
    a_r = Arrow(np.array([3.6, 2.6, 0]), np.array([4.6, 2.6, 0]), buff=0).set_color(R_COL)
    a_l = Arrow(np.array([-3.6, 2.6, 0]), np.array([-4.6, 2.6, 0]), buff=0).set_color(L_COL)
    lab = at(mixed("速率 ", "$v=1$", "，反轉率 ", "$\\gamma=1$", font_size=28), 0, 2.6)
    assert_in_frame(base, nums, walls, dots, clock, a_r, a_l, lab)
    for m in (a_r, a_l, dots):
        assert_apart(lab, m)
    assert_apart(clock, dots)
    s.play(FadeIn(base), FadeIn(nums), FadeIn(dots), FadeIn(lab), GrowArrow(a_r), GrowArrow(a_l), FadeIn(clock),
           run_time=0.8)
    s.say("問：走路人以速率 ", "$v=1$", " 一直走，不會停", wait=1.8)
    s.say("平均每單位時間反轉 ", "$\\gamma=1$", " 次（時間隨機）", wait=1.8)
    nsteps = len(lanes[0][0]) - 1
    for k in range(1, nsteps + 1):
        s.play(*[d.animate.move_to(xp(lanes[i][0][k], i)) for i, d in enumerate(dots)],
               Transform(clock, at(tx(r"t=%.2f" % (k * LANE_DT), 30, GREY_A), -5.6, 2.3)),
               run_time=0.3, rate_func=linear)
    unt = [d for d, (_, u0) in zip(dots, lanes) if u0]
    s.play(ShowCreation(walls), *[d.animate.set_fill(RED, 1).scale(1.4) for d in unt], run_time=0.8)
    s.say("$t=3$", " 時有沒有人在 ", "$x=\\pm3$", " 之外？", wait=2.0)
    s.say("一次都沒轉身的（紅）有多少？（先猜）", wait=3.0)
    s.clear_stage()


def part_equations(s):
    nb = numbers()
    a, b, rich, formula = nb["lat"]
    rl = VGroup(tx(r"\partial_tR=-v\,\partial_xR-\gamma R+\gamma L", 30),
                tx(r"\partial_tL=v\,\partial_xL-\gamma L+\gamma R", 30)).arrange(DOWN, buff=0.15, aligned_edge=LEFT)
    rl[0][2].set_color(R_COL)
    rl[1][2].set_color(L_COL)
    e1 = eq_label(rl, "(5.39)--(5.44)")
    e2 = eq_label(tx(r"\Sigma=R+L,\ \Delta=R-L:\quad\partial_t\Sigma=-v\,\partial_x\Delta,\quad"
                     r"\partial_t\Delta=-v\,\partial_x\Sigma-2\gamma\Delta", 30), "(5.45)--(5.48)")
    e3 = VGroup(tx(r"\Sigma_{tt}+2\gamma\,\Sigma_t=v^2\,\Sigma_{xx}", 38, YELLOW),
                zh("電報方程", font_size=28, color=YELLOW), tx(r"\text{(Ex. 5.4)}", 22, GREY_B)).arrange(RIGHT, buff=0.4)
    lt = mixed("Ex. 5.3：", "$Q=1-\\gamma\\delta t,\\ l=v\\,\\delta t$", " 的有慣性漫步，密度＝機率除以 ", "$2l$",
               font_size=26)
    rows = [[r"\delta t", "0.01", "0.005", zh("外插", font_size=24, color=YELLOW), zh("公式", font_size=24, color=YELLOW)],
            [r"\Sigma(0,1)", f4(a), f4(b), f4(rich), f4(formula)]]
    tab = D.tex_table(rows, [1.6, 1.4, 1.4, 1.4, 1.4], row_h=0.6, font_size=30)
    for m, y in ((e1, 2.45), (e2, 1.25), (e3, 0.25), (lt, -0.8), (tab, -2.05)):
        fit(m, 12.8)
        at(m, 0, y)
    assert_in_frame(e1, e2, e3, lt, tab)
    for p, q in ((e1, e2), (e2, e3), (e3, lt), (lt, tab)):
        assert_apart(p, q)
    assert_table(tab, 5)
    s.play(Write(e1), run_time=1.5)
    s.say("右行者 ", "$R$", " 往右跑、左行者 ", "$L$", " 往左跑", wait=2.2)
    s.say("每單位時間有 ", "$\\gamma$", " 的比例轉身，換到另一邊", wait=2.5)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("相加、相減：總數 ", "$\\Sigma$", " 與差 ", "$\\Delta$", wait=2.2)
    s.play(FadeIn(e3, shift=0.2 * DOWN), run_time=0.8)
    s.say("消去 ", "$\\Delta$", "：電報方程，比擴散方程多了 ", "$\\Sigma_{tt}$", wait=2.8)
    s.play(FadeIn(lt), run_time=0.8)
    s.say("Ex. 5.3：它就是 W502 有慣性漫步的連續極限", wait=2.5)
    s.play(FadeIn(tab), run_time=0.8)
    s.say("$\\delta t$", " 減半：0.3392、0.3380，外插 0.3368＝公式", wait=3.0)
    s.clear_stage()


def part_solution(s):
    nb = numbers()
    l1 = tx(r"\Sigma(x,t)=\frac{\gamma e^{-\gamma t}}{2v}\Big[I_0(\gamma s)+\frac{t}{s}\,I_1(\gamma s)\Big]"
            r"\Theta(vt-|x|)", 34)
    l2 = tx(r"+\frac{e^{-\gamma t}}{2}\Big[\delta(x-vt)+\delta(x+vt)\Big],\qquad s=\sqrt{t^2-x^2/v^2}", 34)
    sol = VGroup(l1, l2).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    lab = zh("(5.61)，更正後", font_size=22, color=GREY_B)
    at(sol, 0, 1.9)
    lab.next_to(sol, DOWN, buff=0.15).align_to(sol, RIGHT)
    n1 = at(mixed("拖在後面：轉過身的人", font_size=26, color=GREY_A), -3.3, 0.0)
    n2 = at(mixed("兩個尖峰：一次都沒轉身（", "$e^{-\\gamma t}$", "）", font_size=26, color=GREY_A), 3.0, 0.0)
    box = at(D.erratum_box("p. 107，(5.60)(5.61)", r"\frac{1}{2}\Big(2\gamma+\frac{\partial}{\partial t}\Big)",
                           r"\frac{1}{2v}\Big(2\gamma+\frac{\partial}{\partial t}\Big)", width=8.0), 0, -1.75)
    assert_in_frame(sol, lab, n1, n2, box)
    for p, q in ((sol, lab), (lab, n1), (lab, n2), (sol, n1), (sol, n2), (n1, n2), (n1, box), (n2, box)):
        assert_apart(p, q)
    s.play(Write(l1), run_time=1.6)
    s.say("解 (5.61)：連續部分用 Bessel 函數 ", "$I_0$", "、", "$I_1$", wait=2.5)
    s.play(FadeIn(n1), run_time=0.6)
    s.say("$\\Theta(vt-|x|)$", "：只在 ", "$|x|<vt$", " 裡面有人", wait=2.5)
    s.play(FadeIn(l2, shift=0.2 * DOWN), FadeIn(lab), FadeIn(n2), run_time=0.8)
    s.say("外加兩個 ", "$\\delta$", " 尖峰：一直沒轉身、跑在最前面的人", wait=2.8)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (5.60)(5.61) 少了 ", "$1/v$", wait=2.2)
    s.say("$v=2,\\ t=3$", " 照書上積分得 %s，不是 1" % f4(nb["book"]), wait=3.0)
    s.clear_stage()


def fig_axes():
    ax = D.plot_axes((0, 9, 1), (0, 0.4, 0.1), width=6.8, height=4.4, x_ticks=(0, 1, 3, 5, 8),
                     y_ticks=(0.1, 0.2, 0.3, 0.4), y_tick_labels=["0.1", "0.2", "0.3", "0.4"], x_label="x",
                     y_label=r"\Sigma")
    ax.move_to(np.array([-3.15, -0.3, 0]))
    return ax


def part_fig53(s):
    nb = numbers()
    fig = nb["fig"]
    ax = fig_axes()
    curves = {t: cont_curve(ax, t, T_COL[t]) for t in TS}
    spikes = {t: spike_arrow(ax, t, SPIKE_H[t], T_COL[t], g4(fig[t]["spike"])) for t in TS}
    hdr = lambda t: zh(t, font_size=20, color=YELLOW)
    rows = [["t", r"\Sigma(0,t)", hdr("前緣後方"), hdr("每個尖峰"), hdr("連續部分")]]
    for t in TS:
        rows.append([str(t), f4(fig[t]["c0"]), g4(fig[t]["front"]), g4(fig[t]["spike"]), f4(fig[t]["cont"])])
    tab = D.tex_table(rows, [0.5, 1.15, 1.35, 1.35, 1.25], row_h=0.55, font_size=22)
    at(tab, 4.0, 1.75)
    for i, t in enumerate(TS, start=1):
        D.table_cell(tab, i, 0, 5).set_color(T_COL[t])
    cons = at(mixed("連續部分＋", "$2\\times$", "尖峰＝1", font_size=26), 4.0, -0.6)
    caus = at(mixed("$|x|>vt$", "：0（因果性）", font_size=26), 4.0, -1.3)
    plot = VGroup(ax, *curves.values(), *spikes.values())
    assert_in_frame(plot, tab, cons, caus)
    for p, q in ((plot, tab), (plot, cons), (plot, caus), (tab, cons), (cons, caus)):
        assert_apart(p, q)
    for i in range(len(TS)):
        for j in range(i + 1, len(TS)):
            assert_apart(spikes[TS[i]][1], spikes[TS[j]][1])
    assert_table(tab, 5)
    s.play(FadeIn(ax), FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:5])), run_time=0.8)
    s.say("Fig. 5.3：", "$t=1,3,5,8$", "（只畫 ", "$x\\ge0$", "，左右對稱）", wait=2.0)
    for i, t in enumerate(TS, start=1):
        s.play(ShowCreation(curves[t]), GrowArrow(spikes[t][0]), FadeIn(spikes[t][1]),
               FadeIn(VGroup(*tab[0][5 * i:5 * i + 5])), run_time=1.2)
        if t == 1:
            s.say("$t=1$", "：中間 0.3368，前緣後方還有 0.2759", wait=2.2)
            s.say("箭頭是尖峰：每邊 ", "$\\frac12e^{-t}=0.1839$", wait=2.2)
        elif t == 3:
            s.say("$t=3$", "：尖峰只剩 0.02489，前緣跑到 ", "$x=3$", wait=2.2)
        else:
            s.wait(1.2)
    s.play(FadeIn(caus), run_time=0.6)
    s.say("因果性：", "$|x|>vt$", " 的地方一個人都沒有", wait=2.5)
    s.play(FadeIn(cons), run_time=0.6)
    s.say("守恆：連續部分 ", "$1-e^{-t}$", " 加上兩個尖峰＝1", wait=2.5)
    return VGroup(ax, *curves.values(), *spikes.values()), tab, VGroup(cons, caus)


def part_large_t(s, plot, tab_old, notes):
    nb = numbers()
    fig = nb["fig"]
    ax = plot[0]
    gs = VGroup(*[gauss_curve(ax, t, T_COL[t]) for t in TS])
    hdr = lambda t: zh(t, font_size=20, color=YELLOW)
    rows = [["t", r"\Sigma(0,t)", r"\frac{1}{\sqrt{4\pi Dt}}", hdr("變異數"), r"2Dt-\frac12"]]
    for t in TS:
        rows.append([str(t), f4(fig[t]["c0"]), f4(fig[t]["gauss"]), f4(fig[t]["var"]), f"{2 * V.tele_D() * t - 0.5:.1f}"])
    tab = D.tex_table(rows, [0.5, 1.15, 1.35, 1.25, 1.2], row_h=0.6, font_size=22)
    at(tab, 4.0, 1.6)
    for i, t in enumerate(TS, start=1):
        D.table_cell(tab, i, 0, 5).set_color(T_COL[t])
        D.table_cell(tab, i, 2, 5).set_color(D.GAUSS_COLOR)
    e1 = at(tx(r"D=\frac{v^2}{2\gamma}=%.1f,\quad l_p=\frac{v}{\gamma}=%d,\quad D=\tfrac12\,v\,l_p"
               % (V.tele_D(), round(V.tele_lp())), 28), 4.0, -0.5)
    e2 = at(tx(r"\langle x^2\rangle=\frac{v^2}{\gamma}\Big[t-\frac{1-e^{-2\gamma t}}{2\gamma}\Big]", 28), 4.0, -1.7)
    lab = tx(r"\text{(Ex. 5.5--5.6)}", 20, GREY_B).next_to(e2, DOWN, buff=0.15).align_to(e2, RIGHT)
    for m in (e1, e2):
        fit(m, 5.6)
    assert_in_frame(tab, e1, e2, lab, gs)
    for p, q in ((plot, tab), (plot, e1), (plot, e2), (tab, e1), (e1, e2), (e2, lab), (plot, lab)):
        assert_apart(p, q)
    assert_table(tab, 5)
    s.play(FadeOut(tab_old), FadeOut(notes), run_time=0.6)
    s.play(ShowCreation(gs), run_time=1.5)
    s.say("虛線：擴散方程的高斯，", "$D=v^2/2\\gamma=0.5$", wait=2.5)
    s.say("$t=1$", " 的高斯跑到 ", "$x>1$", "：擴散方程不管因果性", wait=2.5)
    s.play(FadeIn(tab), run_time=0.8)
    s.say("原點：0.3368 對 0.3989，到 ", "$t=8$", " 已經 0.1388 對 0.1410", wait=3.0)
    s.play(FadeIn(e1), run_time=0.6)
    s.say("持續長度 ", "$l_p=v/\\gamma=1$", "：", "$D=\\frac12vl_p$", wait=2.5)
    s.play(FadeIn(e2), FadeIn(lab), run_time=0.6)
    s.say("精確變異數：大 ", "$t$", " 時 ", "$\\to2Dt-0.5$", wait=2.8)
    s.clear_stage()


def part_reveal(s):
    nb = numbers()
    mc = nb["mc"]
    ux, uy, org = 0.95, 11.0, np.array([-3.3, -2.3, 0])
    bars = VGroup()
    for k, h in enumerate(mc["hist"]):
        xc = (mc["edges"][k] + mc["edges"][k + 1]) / 2
        r = Rectangle(width=0.2 * ux * 0.92, height=max(h * uy, 1e-3)).set_fill(D.BAR_COLOR, 0.7).set_stroke(D.BAR_COLOR, 1)
        r.move_to(org + xc * ux * RIGHT + h * uy / 2 * UP)
        r.xc = xc
        bars.add(r)
    base = Line(org + 3.6 * ux * LEFT, org + 3.6 * ux * RIGHT).set_stroke(GREY_C, 2)
    nums = VGroup(*[tx(str(x), 20, GREY_B).move_to(org + x * ux * RIGHT + 0.28 * DOWN) for x in range(-3, 4)])
    cur = D.curve(lambda x: float(V.tele_cont(x, T_MC)), -T_MC * 0.999, T_MC * 0.999, ux=ux, uy=uy, origin=org,
                  color=D.GAUSS_COLOR)
    spk = VGroup(*[Arrow(org + x * ux * RIGHT, org + x * ux * RIGHT + 0.25 * uy * UP, buff=0, thickness=3)
                   .set_color(RED) for x in (-T_MC, T_MC)])
    sl = tx(r"\tfrac12e^{-3}", 24, RED).next_to(spk[1].get_end(), UP, buff=0.08)
    hdr = lambda t: zh(t, font_size=24, color=YELLOW)
    rows = [[VectorizedPoint(), hdr("10 萬人模擬"), hdr("公式")],
            [r"|x|<1.5", f4(mc["frac"]), f4(mc["formula"])],
            [zh("沒轉身", font_size=24), f4(mc["unturned"]), r"e^{-3}=%s" % f4(mc["e3"])]]
    tab = D.tex_table(rows, [1.6, 2.1, 1.9], row_h=0.65, font_size=28)
    at(tab, 3.75, 0.9)
    ans = at(zh("x = 3 之外：一個都沒有", font_size=28, color=YELLOW), 3.75, -1.0)
    plot = VGroup(base, nums, bars, cur, spk, sl)
    assert_in_frame(plot, tab, ans)
    for p, q in ((plot, tab), (plot, ans), (tab, ans)):
        assert_apart(p, q)
    assert_table(tab, 3)
    inner = [r for r in bars if abs(r.xc) < 1.5]
    s.play(FadeIn(base), FadeIn(nums), LaggedStart(*[GrowFromEdge(r, DOWN) for r in bars], lag_ratio=0.05),
           run_time=1.5)
    s.say("揭曉：10 萬人模擬 ", "$t=3$", "（種子固定）", wait=2.0)
    s.play(ShowCreation(cur), GrowArrow(spk[0]), GrowArrow(spk[1]), FadeIn(sl), run_time=1.0)
    s.say("長條＝模擬、橘線＝公式 (5.61)、紅箭頭＝尖峰", wait=2.5)
    s.play(FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:6])), *[r.animate.set_fill(YELLOW, 0.8) for r in inner],
           run_time=0.8)
    s.say("$|x|<1.5$", " 的比例：模擬 0.6060，公式 0.6051", wait=2.5)
    s.play(FadeIn(VGroup(*tab[0][6:])), run_time=0.6)
    s.say("一次都沒轉身：0.0488，", "$e^{-3}=0.0498$", "：大約 5%", wait=2.5)
    s.play(FadeIn(ans), run_time=0.6)
    s.say("沒有人跑得比 ", "$v$", " 快：", "$x=\\pm3$", " 之外是空的", wait=2.5)
    s.clear_stage()


def part_close(s):
    s.say("電報方程是有慣性的擴散：前緣以 ", "$v$", " 跑", wait=2.5)
    s.say("尖峰以 ", "$e^{-\\gamma t}$", " 變小，後面的人慢慢變成高斯", wait=3.0)


class W504Telegrapher(ChuScene):
    def construct(self):
        self.show_header("W504　電報方程：兩個尖峰＋拖在後面的分布")
        check_math()
        part_question(self)
        part_equations(self)
        part_solution(self)
        plot, tab, notes = part_fig53(self)
        part_large_t(self, plot, tab, notes)
        part_reveal(self)
        part_close(self)
        self.wait(1.0)
