"""
W402 夾在兩面吸收牆之間：正弦模態（§4.1.2，pp. 73–76，(4.9)–(4.15)，Ex. 4.1–4.4）

問：懸崖在 0 和 6，醉漢從 3 出發。走 10 步都沒掉下去的路線有幾條？
  1. 逐步遞迴（BC.box_dp(N, 6, 3)）：N = 0..4 共 1、2、4、6、12；每兩步乘 3 → N = 10：81, 0, 162, 0, 81（324／1024）。
  2. 轉移矩陣的特徵值 2cos(mπ/6)（BC.box_eigs）、特徵向量 sin(mπx/6)；(4.15) 的格子版（BC.box_modes）：
     (√3)^10 = 243，x = 3：(2/6)(243+243) = 162；x = 1：(2/6)(½·243+½·243) = 81；連續極限 1 − 2z + zq²。
  3. Ex. 4.3：鏡像點 3+12k（正）、−3+12k（負）（BC.box_image_points、BC.box_images）→ 81、162、81。
  4. 筆誤框 p. 75 (4.12)(4.13)：sin q x1 → sin q x2。
畫面上的數字都由 utils.walk_bc 算出，check_math() 核對。
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
from utils import walk_bc as BC
from utils import walk_draw as D

LW, X0, NB = 6, 3, 10
DEAD = D.WRONG_COLOR
POS_C = ORANGE
NEG_C = BLUE_C
MODE_COLORS = [YELLOW, GREEN_C, TEAL, BLUE_C, PURPLE_B]


def C(N, y):
    """N 步、位移 y 的無限制走法數。"""
    if abs(y) > N or (N + y) % 2:
        return 0
    return math.comb(N, (N + y) // 2)


def image_terms(N, x, K=1):
    """Ex. 4.3 的非零項：[(正負號, 位移, C)]；正鏡像 3+12k、負鏡像 −3+12k。"""
    pos, neg = BC.box_image_points(LW, X0, K)
    out = [(+1, x - p, C(N, x - p)) for p in pos] + [(-1, x - q, C(N, x - q)) for q in neg]
    return sorted([t for t in out if t[2]], key=lambda t: abs(t[1]))


@lru_cache(maxsize=None)
def numbers():
    nb = {"dp": {N: BC.box_dp(N, LW, X0) for N in list(range(5)) + [NB]}}
    nb["eigs"] = BC.box_eigs(LW)
    nb["p243"] = round(nb["eigs"][0] ** NB)
    nb["modes"] = {x: BC.box_modes(NB, LW, X0, x) for x in range(1, LW)}
    nb["images"] = {x: BC.box_images(NB, LW, X0, x) for x in range(1, LW)}
    return nb


def check_math():
    nb = numbers()
    dp = nb["dp"]
    assert dp[0] == {3: 1} and dp[1] == {2: 1, 4: 1} and dp[2] == {1: 1, 3: 2, 5: 1}
    assert dp[3] == {2: 3, 4: 3} and dp[4] == {1: 3, 3: 6, 5: 3}
    assert [sum(dp[N].values()) for N in range(5)] == [1, 2, 4, 6, 12]
    full = {x: dp[NB].get(x, 0) for x in range(1, LW)}
    assert full == {1: 81, 2: 0, 3: 162, 4: 0, 5: 81} and sum(full.values()) == 324 and 2 ** NB == 1024
    assert 324 == 4 * 3 ** 4 and 12 == 3 * 4
    assert [f"{e:.4f}" for e in nb["eigs"]] == ["1.7321", "1.0000", "0.0000", "-1.0000", "-1.7321"]
    assert BC.box_eig_residual(LW) < 1e-12
    assert nb["p243"] == 243 and abs(nb["eigs"][4] ** NB - 243) < 1e-9
    assert [round(math.sin(m * math.pi * X0 / LW)) for m in range(1, LW)] == [1, 0, -1, 0, 1]
    assert abs(math.sin(math.pi / 6) - 0.5) < 1e-12 and abs(math.sin(5 * math.pi / 6) - 0.5) < 1e-12
    assert 2 * (243 + 243) // 6 == 162 and 2 * (243 / 2 + 243 / 2) / 6 == 81
    assert all(abs(nb["modes"][x] - full[x]) < 1e-9 for x in full)
    assert nb["images"] == full
    assert BC.box_image_points(LW, X0) == ([-9, 3, 15], [-15, -3, 9])
    t3, t1 = image_terms(NB, 3), image_terms(NB, 1)
    assert [(s, c) for s, _, c in t3] == [(1, 252), (-1, 45), (-1, 45)]
    assert [(s, y, c) for s, y, c in t1] == [(1, -2, 210), (-1, 4, 120), (-1, -8, 10), (1, 10, 1)]
    assert sum(s * c for s, _, c in t3) == 162 and sum(s * c for s, _, c in t1) == 81
    assert all(sum(s * c for s, _, c in image_terms(NB, x, 3)) == full[x] for x in full)
    for z in (0.2, 0.4):
        for m in (1, 2):
            assert abs(BC.box_amplitude(z, LW, m, 2) - BC.box_amplitude_formula(z, LW, m, 2)) < 1e-9


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


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


def corridor(xmax, x0, unit, walls, y=0.0, xc=0.0):
    """數線 0..xmax（D.number_line 取右半）＋紅色懸崖＋走路人 → VGroup(線, 懸崖, 走路人)。"""
    nl = D.number_line(xmax, unit=unit, label_every=1, font_size=24)
    nl.shift((xc - xmax / 2 * unit) * RIGHT + y * UP)
    p = lambda x: D.line_point(nl, x)
    ticks = VGroup(*[nl[1][k] for k in range(xmax, 2 * xmax + 1)])
    nums = VGroup(*[nl[2][k] for k in range(xmax, 2 * xmax + 1)])
    base = Line(p(0), p(xmax)).set_stroke(GREY_B, 2)
    cliffs = VGroup(*[Line(p(w) + 0.55 * UP, p(w) + 0.55 * DOWN).set_stroke(DEAD, 6) for w in walls])
    labs = VGroup(*[zh("懸崖", font_size=24, color=DEAD).next_to(c, UP, buff=0.1) for c in cliffs])
    walker = D.walker_dot().move_to(p(x0) + 0.35 * UP)
    return VGroup(VGroup(base, ticks, nums), VGroup(cliffs, labs), walker)


def sine_mode(m, xc, yc, unit=0.32, amp=0.5):
    """sin(mπx/6) 在 0..6 的曲線＋格點上的值（圓點）＋兩端紅牆。"""
    c = MODE_COLORS[m - 1]
    o = np.array([xc - 3 * unit, yc, 0])
    base = Line(o, o + 6 * unit * RIGHT).set_stroke(GREY_C, 1.5)
    walls = VGroup(*[Line(o + w * unit * RIGHT + 0.6 * amp * DOWN, o + w * unit * RIGHT + 0.6 * amp * UP)
                     .set_stroke(DEAD, 4) for w in (0, LW)])
    fn = lambda x: math.sin(m * math.pi * x / LW)
    curve = D.curve(fn, 0, LW, ux=unit, uy=amp, origin=o, color=c, stroke_width=2.5, samples=120)
    dots = VGroup(*[Dot(o + x * unit * RIGHT + fn(x) * amp * UP, radius=0.05).set_fill(c, 1) for x in range(LW + 1)])
    return VGroup(base, walls, curve, dots)


# ---------------------------------------------------------------- 各段
def part_question(s):
    pic = corridor(LW, X0, 1.3, (0, LW), y=0.6)
    guess = at(zh("走 10 步都沒掉下去的路線有幾條？", font_size=34, color=YELLOW), 0, -1.3)
    assert_in_frame(pic, guess)
    assert_apart(pic, guess)
    s.play(FadeIn(pic[0]), FadeIn(pic[1]), run_time=1.0)
    s.play(FadeIn(pic[2], shift=0.3 * DOWN), run_time=0.6)
    s.say("問：懸崖在 0 和 6，醉漢從正中間 3 出發", wait=3.0)
    s.play(FadeIn(guess), run_time=0.6)
    s.say("走 10 步都沒掉下去的路線有幾條？先猜", wait=3.0)
    s.clear_stage()


def part_recursion(s):
    nb = numbers()
    dp = nb["dp"]
    unit, xc = 1.0, -0.6
    colx = lambda x: xc + (x - 3) * unit
    ys = {0: 2.6, 1: 2.0, 2: 1.4, 3: 0.8, 4: 0.2, NB: -1.05}
    nl = corridor(LW, X0, unit, (), y=-2.05, xc=xc)[0]
    walls = VGroup(*[Line(np.array([colx(w), -2.05, 0]), np.array([colx(w), 2.95, 0])).set_stroke(DEAD, 5)
                     for w in (0, LW)])
    rows, labels, tots = VGroup(), VGroup(), VGroup()
    for N in (0, 1, 2, 3, 4, NB):
        row = VGroup()
        for x in range(1, LW):
            v = dp[N].get(x, 0)
            if v or N == NB:
                row.add(at(tx(str(v), 30, WHITE if v else GREY_D), colx(x), ys[N]))
        rows.add(row)
        labels.add(at(tx(f"N={N}", 28, GREY_A), -5.6, ys[N]))
        tot = sum(dp[N].values())
        tots.add(at(mixed("共 ", f"${tot}$" if N != NB else f"${tot}\\ /\\ {2 ** NB}$", font_size=26,
                          color=YELLOW), 4.6, ys[N]))
    dots_v = at(tx(r"\vdots", 30, GREY_B), colx(3), -0.45)
    a1 = Arrow(rows[1][0].get_bottom(), rows[2][1].get_top(), buff=0.06, thickness=2).set_color(TEAL)
    a2 = Arrow(rows[1][1].get_bottom(), rows[2][1].get_top(), buff=0.06, thickness=2).set_color(TEAL)
    x3 = at(tx(r"\times3", 30, ORANGE), 6.1, (ys[2] + ys[4]) / 2)
    br = Brace(VGroup(tots[2], tots[4]), RIGHT, buff=0.1).set_color(ORANGE)
    x3.next_to(br, RIGHT, buff=0.08)
    assert_in_frame(nl, walls, rows, labels, tots, dots_v, br, x3)
    for r in rows:
        for t in r:
            assert all(abs(t.get_center()[0] - w.get_center()[0]) > 0.35 for w in walls)
    assert_apart(labels, walls[0])
    assert_apart(walls[1], tots)
    assert_apart(rows[5], nl)
    s.play(FadeIn(nl), ShowCreation(walls), run_time=1.0)
    s.play(FadeIn(labels[0]), FadeIn(rows[0]), FadeIn(tots[0]), run_time=0.8)
    s.say("一步一步數：每一格寫「走到這裡、還活著」的條數", wait=3.0)
    s.play(FadeIn(labels[1]), FadeIn(rows[1]), FadeIn(tots[1]), run_time=0.8)
    s.play(FadeIn(labels[2]), FadeIn(rows[2]), FadeIn(tots[2]), GrowArrow(a1), GrowArrow(a2), run_time=1.0)
    s.say("每一格＝上一列左右兩格相加；牆上永遠是 0", wait=3.0)
    s.play(FadeOut(VGroup(a1, a2)), FadeIn(labels[3]), FadeIn(rows[3]), FadeIn(tots[3]), run_time=0.8)
    s.play(FadeIn(labels[4]), FadeIn(rows[4]), FadeIn(tots[4]), run_time=0.8)
    s.say("$N=4$", "：3、6、3，共 12", wait=2.5)
    s.play(GrowFromCenter(br), FadeIn(x3), run_time=0.8)
    s.say("1、2、1 變成 3、6、3：每兩步剛好乘 3", wait=3.0)
    s.play(FadeIn(dots_v), FadeIn(labels[5]), FadeIn(rows[5]), FadeIn(tots[5]), run_time=1.2)
    s.say("$N=10$", "：81、0、162、0、81，共 324 條", wait=3.0)
    s.say("1024 條裡只剩 324 條活著；為什麼是 3？", wait=3.0)
    s.clear_stage()


def part_eigen(s):
    nb = numbers()
    mat = tx(r"T=\begin{pmatrix}0&1&0&0&0\\1&0&1&0&0\\0&1&0&1&0\\0&0&1&0&1\\0&0&0&1&0\end{pmatrix}", 28)
    at(mat, -4.75, 1.75)
    mlab = zh("鄰居相加（牆外是 0）", font_size=22, color=GREY_A).next_to(mat, DOWN, buff=0.15)
    vals = [r"\sqrt3", "1", "0", "-1", r"-\sqrt3"]
    rows = [["m", "1", "2", "3", "4", "5"], [r"2\cos(m\pi/6)"] + vals,
            [r"\approx"] + [f"{e:.4f}" for e in nb["eigs"]]]
    tab = D.tex_table(rows, [2.3] + [1.15] * 5, row_h=0.5, font_size=26)
    tab.move_to(np.array([2.6, 2.85 - tab.get_height() / 2, 0]))
    eig = eq_label(tx(r"\sin\frac{m\pi(x-1)}{6}+\sin\frac{m\pi(x+1)}{6}=2\cos\frac{m\pi}{6}\,\sin\frac{m\pi x}{6}",
                      28), "(4.9)(4.10)")
    eig.next_to(tab, DOWN, buff=0.35)
    modes = VGroup(*[sine_mode(m, -5.2 + 2.6 * (m - 1), -1.75) for m in range(1, LW)])
    mlabs = VGroup(*[tx(f"m={m}", 24, MODE_COLORS[m - 1]).next_to(modes[m - 1], UP, buff=0.12) for m in range(1, LW)])
    elabs = VGroup(*[tx(r"\times" + (f"({vals[m - 1]})" if vals[m - 1][0] == "-" else vals[m - 1]), 24,
                        MODE_COLORS[m - 1]).next_to(modes[m - 1], DOWN, buff=0.1)
                     for m in range(1, LW)])
    assert_in_frame(mat, mlab, tab, eig, modes, mlabs, elabs)
    assert_table(tab, 6)
    assert_apart(mat, tab)
    assert_apart(mlab, eig)
    assert_apart(tab, eig)
    assert_apart(eig, mlabs)
    assert_apart(mlab, mlabs)
    for i in range(4):
        assert_apart(modes[i], modes[i + 1], 0.2)
    s.play(FadeIn(mat), FadeIn(mlab), run_time=1.0)
    s.say("內點 1 到 5：一步＝乘上 ", "$5\\times5$", " 的轉移矩陣 ", "$T$", wait=3.0)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.2)
    s.say("特徵值 ", "$2\\cos(m\\pi/6)$", "：最大的是 ", "$\\sqrt3$", wait=3.0)
    s.play(FadeIn(eig), run_time=1.0)
    s.say("特徵向量是正弦：", "$\\sin(m\\pi x/6)$", wait=2.5)
    s.play(LaggedStart(*[ShowCreation(md[2]) for md in modes], lag_ratio=0.2),
           LaggedStart(*[FadeIn(VGroup(md[0], md[1], md[3])) for md in modes], lag_ratio=0.2),
           FadeIn(mlabs), run_time=2.4)
    s.say("在 ", "$x=0$", " 和 ", "$x=6$", " 自動是 0：牆的條件白送", wait=3.0)
    s.play(FadeIn(elabs), run_time=0.8)
    s.say("每走一步，第 ", "$m$", " 個模態乘上 ", "$2\\cos(m\\pi/6)$", wait=3.0)
    s.say("兩步乘 ", "$(\\sqrt3)^2=3$", "：這就是「每兩步乘 3」", wait=3.0)
    s.clear_stage()


def part_modes(s):
    nb = numbers()
    p = nb["p243"]
    e1 = eq_label(tx(r"C_N(x;x_0)=\sum_{m=1}^{L-1}\frac2L\sin\frac{m\pi x_0}{L}\sin\frac{m\pi x}{L}"
                     r"\Big(2\cos\frac{m\pi}{L}\Big)^N", 32), "(4.15)")
    at(e1, 0, 2.45)
    e2 = at(tx(r"x_0=3:\ \ \sin\frac{m\pi}{2}=1,\ 0,\ -1,\ 0,\ 1", 30), -2.9, 1.2)
    e3 = at(tx(r"(\pm\sqrt3)^{10}=" + str(p) + r",\qquad 0^{10}=0", 30), 3.0, 1.2)
    e4 = at(tx(r"x=3:\ \ \frac26\,(" + f"{p}+{p}" + r")=" + f"{round(nb['modes'][3])}", 32, YELLOW), -3.0, 0.0)
    e5 = at(tx(r"x=1:\ \ \frac26\Big(\frac12\cdot" + str(p) + r"+\frac12\cdot" + str(p) + r"\Big)="
               + f"{round(nb['modes'][1])}", 32, YELLOW), 2.9, 0.0)
    e6 = eq_label(tx(r"2\cos q\approx2-q^2\ \ \Rightarrow\ \ \frac{1}{1-2z\cos q}\approx\frac{1}{1-2z+zq^2}", 32),
                  "(4.14)(4.15)")
    at(e6, 0, -1.5)
    assert_in_frame(e1, e2, e3, e4, e5, e6)
    assert_apart(e1, e2)
    assert_apart(e2, e3)
    assert_apart(e2, e4)
    assert_apart(e4, e5)
    assert_apart(e4, e6)
    assert_apart(e5, e6)
    s.play(FadeIn(e1), run_time=1.2)
    s.say("把起點拆成模態，各自乘 ", "$N$", " 次，再加回來", wait=3.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("從正中間出發：偶數 ", "$m$", " 的係數是 0", wait=2.5)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("$m=3$", " 的特徵值是 0：只剩 ", "$m=1,5$", wait=2.5)
    s.play(FadeIn(e4), run_time=0.8)
    s.play(FadeIn(e5), run_time=0.8)
    s.say("162 和 81：和逐步遞迴數的一樣", wait=3.0)
    s.play(FadeIn(e6, shift=0.2 * UP), run_time=1.0)
    s.say("連續極限：", "$2\\cos q\\approx2-q^2$", "，分母變成 ", "$1-2z+zq^2$", wait=3.5)
    s.clear_stage()


def part_images(s):
    nb = numbers()
    unit = 0.4
    nl = D.number_line(15, unit=unit, label_every=3, font_size=20)
    at(nl, 0, 1.5)
    p = lambda x: D.line_point(nl, x)
    box = Rectangle(width=LW * unit, height=0.3).set_stroke(width=0).set_fill(GREEN_E, 0.35)
    box.move_to((p(0) + p(LW)) / 2 + 0.15 * UP)
    walls = VGroup(*[Line(p(w), p(w) + 0.45 * UP).set_stroke(DEAD, 5) for w in (0, LW)])
    pos, neg = BC.box_image_points(LW, X0)
    pd = VGroup(*[Dot(p(x) + 0.42 * UP, radius=0.1).set_fill(POS_C, 1) for x in pos])
    nd = VGroup(*[Dot(p(x) + 0.42 * UP, radius=0.1).set_fill(NEG_C, 1) for x in neg])
    ps = VGroup(*[tx("+", 28, POS_C).next_to(d, UP, buff=0.08) for d in pd])
    ns = VGroup(*[tx("-", 28, NEG_C).next_to(d, UP, buff=0.12) for d in nd])
    leg = VGroup(mixed("正鏡像 ", "$3+12k$", "：", "$\\dots,-9,\\ 3,\\ 15,\\dots$", font_size=24, color=POS_C),
                 mixed("負鏡像 ", "$-3+12k$", "：", "$\\dots,-15,\\ -3,\\ 9,\\dots$", font_size=24, color=NEG_C))
    leg.arrange(RIGHT, buff=0.8)
    at(leg, 0, 0.45)
    f = at(eq_label(tx(r"C_N(x;3)=\sum_k\Big[C(N;x-3-12k)-C(N;x+3-12k)\Big]", 30), "Ex. 4.3"), 0, -0.45)

    def line(x):
        ts = image_terms(NB, x)
        lhs = "".join(("+" if sg > 0 and i else ("-" if sg < 0 else "")) + f"C(10;{y})" for i, (sg, y, _) in enumerate(ts))
        rhs = "".join(("+" if sg > 0 and i else ("-" if sg < 0 else "")) + str(c) for i, (sg, _, c) in enumerate(ts))
        return tx(f"x={x}:\\ \\ " + lhs + "=" + rhs + "=" + str(nb["images"][x]), 28, YELLOW)

    l3 = at(line(3), 0, -1.35)
    l1 = at(line(1), 0, -2.2)
    assert_in_frame(nl, walls, pd, nd, ps, ns, leg, f, l3, l1)
    assert_apart(nl, leg)
    assert_apart(leg, f)
    assert_apart(f, l3)
    assert_apart(l3, l1)
    s.play(FadeIn(nl), FadeIn(box), ShowCreation(walls), run_time=1.0)
    s.say("Ex. 4.3：第三種算法——鏡像點", wait=2.0)
    s.play(FadeIn(pd[1]), FadeIn(ps[1]), run_time=0.6)
    s.play(FadeIn(nd[1]), FadeIn(ns[1]), FadeIn(nd[2]), FadeIn(ns[2]), run_time=0.8)
    s.say("對 0 鏡射得 ", "$-3$", "，對 6 鏡射得 ", "$9$", "：都要減掉", wait=3.0)
    s.play(FadeIn(VGroup(pd[0], pd[2], ps[0], ps[2], nd[0], ns[0])), FadeIn(leg), run_time=1.0)
    s.say("鏡像的鏡像又要加回來：週期 ", "$2L=12$", wait=3.0)
    s.play(FadeIn(f), run_time=1.0)
    s.play(FadeIn(l3), run_time=0.8)
    s.play(FadeIn(l1), run_time=0.8)
    s.say("$N=10$", " 只有附近幾項不是 0：又得到 162、81", wait=3.5)
    s.clear_stage()


def part_erratum(s):
    g = eq_label(tx(r"G(z;x_1,x_2)=\sum_q\sqrt{\frac2L}\,\sin qx_1\cdot\frac{\sqrt{2/L}\,\sin qx_2}{1-2z\cos q}", 34),
                 "(4.12)(4.13)")
    at(g, 0, 1.6)
    box = at(D.erratum_box("p. 75，(4.12)(4.13)", r"\sqrt{2/L}\,\sin q x_1", r"\sqrt{2/L}\,\sin q x_2", width=7.0), 0, -0.6)
    assert_in_frame(g, box)
    assert_apart(g, box)
    s.play(FadeIn(g), run_time=1.0)
    s.say("生成函數：對 ", "$x_1$", " 遞迴，展開成 ", "$\\sin qx_1$", " 的和", wait=3.0)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("遞迴移動的是 ", "$x_1$", "；係數裡留下的是另一端 ", "$x_2$", wait=3.5)
    s.say("(4.14) 已經寫對", wait=2.0)
    s.clear_stage()


def part_close(s):
    s.say("有牆的時候，正弦波取代平面波", wait=2.5)
    s.say("每個模態每走一步就乘上 ", "$2\\cos(m\\pi/L)$", wait=3.0)


class W402TwoWalls(ChuScene):
    def construct(self):
        self.show_header("W402　夾在兩面吸收牆之間：正弦模態")
        check_math()
        part_question(self)
        part_recursion(self)
        part_eigen(self)
        part_modes(self)
        part_images(self)
        part_erratum(self)
        part_close(self)
        self.wait(1.0)
