"""
W405 穩定態＝靜電學；兩牆之間的帳篷（§4.2–4.2.1，pp. 82–86，(4.38)–(4.49)，Ex. 4.6）

問：走廊 0..5，兩端是懸崖；每一步都在 x0 = 2 放進一個醉漢。等很久以後，每一格平均站著幾個人？
  1. (4.38)(4.39) 穩定態方程；泰勒展開 (4.40) → (4.41) Poisson 方程；數值核對 1/6（BC.laplacian_coeff）；筆誤框 3 → 6。
  2. Fick 定律 (4.42)–(4.45) 與靜電的對照表。
  3. 帳篷 (4.46)–(4.49)：格子精確解 0, 6/5, 12/5, 8/5, 4/5, 0（BC.tent_lattice ＝ BC.tent_formula，s0 = 1、β = 1/2）。
  4. Ex. 4.6：流進兩面牆 3/5、2/5（BC.wall_fluxes）＝先掉進左邊懸崖的機率 (L − x0)/L（BC.absorb_left）。
畫面上的數字都由 utils.walk_bc 算出，check_math() 核對。
"""
import sys
from fractions import Fraction as F
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

LW, X0 = 5, 2                      # 牆在 0 與 LW，源在 X0
CLIFF = RED
SRC_COLOR = GREEN_C
TENT_COLOR = YELLOW
FILL_TS = (1, 2, 4, 8, 16, 200)    # 穩定態的累積過程（只畫長條高度）
LS = (5, 10, 100)


# ---------------------------------------------------------------- 數學
def frac_tex(q):
    q = F(q)
    if q.denominator == 1:
        return str(q.numerator)
    sign = "-" if q < 0 else ""
    return sign + r"\tfrac{%d}{%d}" % (abs(q.numerator), q.denominator)


def big_frac(q):
    return frac_tex(q).replace(r"\tfrac", r"\frac")


def fill_history():
    """n_{t+1}(x) = ½[n_t(x−1) + n_t(x+1)] + δ_{x,X0}（牆上是 0）；回傳 {t: [n_t(0..LW)]}。"""
    n = [0.0] * (LW + 1)
    out = {}
    for t in range(1, max(FILL_TS) + 1):
        n = [0.0] + [(n[x - 1] + n[x + 1]) / 2 + (1.0 if x == X0 else 0.0) for x in range(1, LW)] + [0.0]
        if t in FILL_TS:
            out[t] = list(n)
    return out


def walker_frames(T=10, seed=3):
    """示意：每步在 X0 放 1 人、其餘的人各自左右跳；回傳每步的 {編號: 位置}（位置 0 或 LW 表示掉下去）。"""
    rng = np.random.default_rng(seed)
    alive, frames, nid = {}, [], 0
    for _ in range(T):
        new = {i: x + (1 if rng.random() < 0.5 else -1) for i, x in alive.items()}
        new[nid] = X0
        nid += 1
        frames.append(new)
        alive = {i: x for i, x in new.items() if 0 < x < LW}
    return frames


@lru_cache(maxsize=None)
def numbers():
    nb = {"lap": BC.laplacian_coeff(), "w": list(BC.tent_lattice(LW, X0)), "flux": tuple(BC.wall_fluxes(LW, X0))}
    nb["left"] = {L: BC.absorb_left(L, X0) for L in LS}
    w = nb["w"]
    nb["slopes"] = (w[1] - w[0], w[LW] - w[LW - 1])
    return nb


def check_math():
    nb = numbers()
    assert f"{nb['lap']:.4f}" == "0.1667" and abs(nb["lap"] - 1 / 6) < 1e-5
    w = nb["w"]
    assert w == [0, F(6, 5), F(12, 5), F(8, 5), F(4, 5), 0]
    assert [f"{float(v):.1f}" for v in w[1:LW]] == ["1.2", "2.4", "1.6", "0.8"]
    assert [frac_tex(v) for v in w[1:LW]] == [r"\tfrac{6}{5}", r"\tfrac{12}{5}", r"\tfrac{8}{5}", r"\tfrac{4}{5}"]
    assert all(F(BC.tent_formula(LW, X0, x)) == w[x] for x in range(LW + 1))           # (4.49) 在格子上剛好精確
    for x in range(1, LW):                                                              # 格子方程
        assert w[x] == (w[x - 1] + w[x + 1]) / 2 + (1 if x == X0 else 0)
    s0, beta = 1, F(1, 2)
    left, right = nb["slopes"]
    assert (left, right) == (F(6, 5), F(-4, 5))
    assert all(w[x + 1] - w[x] == left for x in range(X0)) and all(w[x + 1] - w[x] == right for x in range(X0, LW))
    assert right - left == -s0 / beta == -2                                             # (4.48)
    assert w[X0] == s0 / beta * F(X0 * (LW - X0), LW)
    assert nb["flux"] == (F(3, 5), F(2, 5)) and sum(nb["flux"]) == 1
    assert nb["flux"] == (w[1] / 2, w[LW - 1] / 2)
    assert nb["left"] == {5: F(3, 5), 10: F(4, 5), 100: F(49, 50)}
    assert nb["left"][5] == F(LW - X0, LW) == nb["flux"][0]
    assert F(BC.absorb_left_solve(LW, X0)) == F(3, 5) and f"{BC.absorb_left_dp(LW, X0):.1f}" == "0.6"
    hist = fill_history()
    assert all(abs(hist[max(FILL_TS)][x] - float(w[x])) < 1e-9 for x in range(LW + 1))
    assert all(hist[FILL_TS[k]][X0] <= hist[FILL_TS[k + 1]][X0] for k in range(len(FILL_TS) - 1))
    frames = walker_frames()
    assert max(sum(1 for x in fr.values() if x == y) for fr in frames for y in range(1, LW)) <= 4


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


def cpt(unit, origin, x):
    return np.array(origin, dtype=float) + x * unit * RIGHT


def corridor(unit, origin, font_size=26):
    """走廊 0..LW：內點白點、兩端紅色懸崖、底下寫格點編號 → VGroup(底線, 點, 懸崖, 編號)。"""
    base = Line(cpt(unit, origin, 0), cpt(unit, origin, LW)).set_stroke(GREY_C, 2)
    dots = VGroup(*[Dot(cpt(unit, origin, x), radius=0.07).set_color(WHITE) for x in range(1, LW)])
    cliffs = VGroup(*[Line(cpt(unit, origin, x) + 0.3 * DOWN, cpt(unit, origin, x) + 0.5 * UP).set_stroke(CLIFF, 6)
                      for x in (0, LW)])
    nums = VGroup(*[tx(str(x), font_size, GREY_A).move_to(cpt(unit, origin, x) + 0.55 * DOWN) for x in range(LW + 1)])
    return VGroup(base, dots, cliffs, nums)


def stack_point(unit, origin, x, k):
    return cpt(unit, origin, x) + (0.25 + 0.24 * k) * UP


def frame_targets(frame, unit, origin):
    """同一格的人往上疊：{編號: 畫面座標}。"""
    out = {}
    for y in range(LW + 1):
        ids = sorted(i for i, x in frame.items() if x == y)
        for k, i in enumerate(ids):
            out[i] = stack_point(unit, origin, y, k) if 0 < y < LW else cpt(unit, origin, y) + 0.25 * UP
    return out


def tent_bars(values, ux, uy, origin, labels=False):
    vals = {x: values[x] for x in range(1, LW)}
    return D.bar_chart(vals, ux=ux, uy=uy, origin=origin, width=0.7, labels=labels,
                       label_fmt=big_frac, font_size=28)


# ---------------------------------------------------------------- 各段
def part_question(s):
    u, org = 1.7, np.array([-4.25, -0.6, 0])
    cor = corridor(u, org)
    labs = VGroup(*[zh("懸崖", font_size=24, color=CLIFF).move_to(cpt(u, org, x) + 1.0 * DOWN) for x in (0, LW)])
    src = Arrow(cpt(u, org, X0) + 2.45 * UP, cpt(u, org, X0) + 1.3 * UP, buff=0).set_color(SRC_COLOR)
    src_lab = zh("每步放進 1 人", font_size=26, color=SRC_COLOR).next_to(src, UP, buff=0.1)
    ask = at(tx(r"w(x)=?", 40, YELLOW), 3.4, 2.1)
    assert_in_frame(cor, labs, src, src_lab, ask)
    assert_apart(src_lab, ask)
    s.play(ShowCreation(cor[0]), FadeIn(cor[1]), FadeIn(cor[3]), run_time=0.8)
    s.play(ShowCreation(cor[2]), FadeIn(labs), run_time=0.6)
    s.say("問：走廊 0..5，兩端是懸崖（踩到就消失）", wait=2.0)
    s.play(GrowArrow(src), FadeIn(src_lab), run_time=0.7)
    s.say("每一步都在 ", "$x_0=2$", " 放進一個醉漢", wait=0.5)
    dots = {}
    for fr in walker_frames():
        tg = frame_targets(fr, u, org)
        anims = []
        for i, p in tg.items():
            if i in dots:
                anims.append(dots[i].animate.move_to(p))
            else:
                dots[i] = Dot(p, radius=0.09).set_color(TENT_COLOR)
                anims.append(FadeIn(dots[i], scale=0.5))
        s.play(*anims, run_time=0.45)
        fallen = [i for i, x in fr.items() if x in (0, LW)]
        if fallen:
            s.play(*[FadeOut(dots.pop(i), shift=0.5 * DOWN) for i in fallen], run_time=0.3)
    s.play(FadeIn(ask), run_time=0.6)
    s.say("等很久以後，每一格平均站著幾個人？（先猜）", wait=3.0)
    s.clear_stage()


def part_poisson(s):
    nb = numbers()
    e1 = at(eq_label(tx(r"w(\vec r)=\frac16\sum_{i=1}^{6}w(\vec\rho_i)+S(\vec r)", 36), "(4.38)(4.39)"), 0, 2.45)
    e2 = at(eq_label(tx(r"\frac16\sum_{i=1}^{6}w(\vec r+a\hat e_i)\approx w(\vec r)+\frac{a^2}{6}\nabla^2w", 36),
                     "(4.40)"), 0, 1.2)
    e3 = at(VGroup(tx(r"-\nabla^2w=\frac{6}{a^2}S(\vec r)", 40, YELLOW), zh("Poisson 方程", font_size=28, color=YELLOW),
                   tx(r"\text{(4.41)}", 22, GREY_B)).arrange(RIGHT, buff=0.45), 0, -0.05)
    num = at(mixed("數值核對：", "$\\dfrac{\\frac16\\sum_i f(\\vec r+a\\hat e_i)-f(\\vec r)}{a^2\\nabla^2f}"
                   "=%.4f$" % nb["lap"], "（＝1/6）", font_size=26), 0, -1.45)
    box = at(D.erratum_box("p. 83，(4.40)(4.41)", r"\frac{a^2}{3}\nabla^2 w,\ \ -\nabla^2 w=\frac{3}{a^2}S",
                           r"\frac{a^2}{6}\nabla^2 w,\ \ -\nabla^2 w=\frac{6}{a^2}S", width=9.6), 0, -1.7)
    assert_in_frame(e1, e2, e3, num, box)
    for a, b in ((e1, e2), (e2, e3), (e3, num), (e3, box)):
        assert_apart(a, b)
    s.play(Write(e1), run_time=1.4)
    s.say("簡單立方格子，每步在源放人：", "$S(\\vec r)$", wait=2.5)
    s.say("穩定態：每格的人數＝六個鄰居的平均＋新放進來的", wait=3.0)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=1.0)
    s.say("泰勒展開：", "$\\pm a\\hat e_i$", " 的一次項兩兩抵銷", wait=2.5)
    s.say("二次項：", "$\\frac16\\cdot2\\cdot\\frac{a^2}{2}=\\frac{a^2}{6}$", "，乘上 ", "$\\nabla^2w$", wait=3.0)
    s.play(FadeIn(e3, shift=0.2 * DOWN), run_time=1.0)
    s.say("移項：穩定態的人數滿足 Poisson 方程", wait=2.5)
    s.play(FadeIn(num), run_time=0.8)
    s.say("數值核對：", "$a=10^{-3}$", " 時比值是 0.1667＝1/6", wait=3.0)
    s.play(FadeOut(num), FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上寫成 ", "$a^2/3$", " 與 ", "$3/a^2$", "：少了一半", wait=3.0)
    s.clear_stage()


def part_fick(s):
    f1 = eq_label(tx(r"\partial_t c=-\nabla\cdot\vec j+S", 36), "(4.42)(4.43)")
    f2 = VGroup(tx(r"\vec j=-\beta\nabla c", 36), zh("Fick 定律", font_size=26, color=GREY_A),
                tx(r"\text{(4.44)}", 22, GREY_B)).arrange(RIGHT, buff=0.4)
    f3 = eq_label(tx(r"\partial_t c=0\ \Rightarrow\ -\beta\nabla^2c=S", 36, YELLOW), "(4.45)")
    left = VGroup(f1, f2, f3).arrange(DOWN, buff=0.75, aligned_edge=LEFT)
    at(left, -3.6, 0.4)
    hdr = lambda t: zh(t, font_size=26, color=YELLOW)
    cell = lambda *p: mixed(*p, font_size=24)
    arrow = r"\leftrightarrow"
    rows = [[hdr("擴散"), VectorizedPoint(), hdr("靜電")],
            [cell("濃度 ", "$c$"), arrow, cell("電位 ", "$\\phi$")],
            [cell("源 ", "$S$"), arrow, cell("電荷 ", "$\\rho$")],
            [cell("吸收牆"), arrow, cell("接地導體（", "$c=0$", "）")],
            [cell("反射牆"), arrow, cell("$\\hat n\\cdot\\nabla c=0$", "（法向流量 0）")]]
    tab = D.tex_table(rows, [1.9, 0.8, 3.6], row_h=0.62, font_size=30)
    at(tab, 3.45, 0.4)
    assert_in_frame(left, tab)
    assert_table(tab, 3)
    assert_apart(left, tab)
    s.play(Write(f1), run_time=1.0)
    s.say("連續版：濃度 ", "$c$", " 的變化＝流進來的－流出去的＋源", wait=3.0)
    s.play(FadeIn(f2, shift=0.2 * RIGHT), run_time=0.8)
    s.say("Fick 定律：粒子從濃的地方流向稀的地方", wait=2.5)
    s.play(FadeIn(f3, shift=0.2 * RIGHT), run_time=0.8)
    s.say("穩定態不再變化：", "$-\\beta\\nabla^2c=S$", "，又是 Poisson 方程", wait=3.0)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.2)
    s.say("和靜電一模一樣：濃度是電位、源是電荷", wait=3.0)
    s.say("吸收牆就是接地導體；反射牆是法向流量為 0", wait=3.0)
    s.clear_stage()


def part_tent(s):
    nb = numbers()
    w = nb["w"]
    ux, uy, org = 0.95, 1.2, np.array([-5.9, -1.8, 0])
    hist = fill_history()
    base = Line(cpt(ux, org, 0), cpt(ux, org, LW)).set_stroke(GREY_C, 2)
    cliffs = VGroup(*[Line(cpt(ux, org, x), cpt(ux, org, x) + 0.8 * UP).set_stroke(CLIFF, 6) for x in (0, LW)])
    nums = VGroup(*[tx(str(x), 26, GREY_A).move_to(cpt(ux, org, x) + 0.35 * DOWN) for x in range(LW + 1)])
    src = Arrow(cpt(ux, org, X0) + 4.6 * UP, cpt(ux, org, X0) + 3.85 * UP, buff=0).set_color(SRC_COLOR)
    src_lab = zh("源", font_size=26, color=SRC_COLOR).next_to(src, RIGHT, buff=0.12)
    bars = tent_bars(hist[FILL_TS[0]], ux, uy, org)[0]
    tcount = at(tx(r"t=%d" % FILL_TS[0], 30, GREY_A), -1.6, 2.5)
    final = tent_bars(w, ux, uy, org, labels=True)
    tent = VMobject().set_points_as_corners([cpt(ux, org, x) + float(w[x]) * uy * UP for x in range(LW + 1)])
    tent.set_stroke(TENT_COLOR, 4)
    chart = VGroup(base, cliffs, nums, src, src_lab, final, tent)

    r1 = tx(r"w(x)=\tfrac12\big[w(x-1)+w(x+1)\big]+\delta_{x,2}", 30)
    r2 = eq_label(tx(r"w(x)=\frac{s_0}{\beta}\begin{cases}\dfrac{x(L-x_0)}{L}, & x\le x_0\\[2mm]"
                     r"\dfrac{x_0(L-x)}{L}, & x\ge x_0\end{cases}", 30), "(4.49)", buff=0.3)
    r3 = mixed("$s_0=1,\\ \\beta=\\tfrac12,\\ L=5,\\ x_0=2$", font_size=26)
    r4 = mixed("斜率 ", "$%s$" % frac_tex(nb["slopes"][0]), " 與 ", "$%s$" % frac_tex(nb["slopes"][1]), font_size=26)
    r5 = mixed("折角處斜率跳 ", "$%s-%s=-2=-s_0/\\beta$" % (frac_tex(nb["slopes"][1]), frac_tex(nb["slopes"][0])),
               font_size=26)
    r5b = tx(r"\text{(4.48)}", 22, GREY_B).next_to(r5, RIGHT, buff=0.3)
    right = VGroup(r1, r2, r3, r4, VGroup(r5, r5b)).arrange(DOWN, buff=0.42)
    fit(right, 6.3)
    at(right, 3.45, 0.2)
    assert_in_frame(chart, tcount, right)
    assert_apart(chart, right)
    assert_apart(tcount, right)
    assert_apart(tcount, VGroup(src, src_lab))
    for k in range(4):
        assert_apart(right[k], right[k + 1])
    assert_apart(final[1], src)

    s.play(ShowCreation(base), FadeIn(nums), ShowCreation(cliffs), GrowArrow(src), FadeIn(src_lab), run_time=1.0)
    s.play(FadeIn(r1), run_time=0.8)
    s.say("一維格子：每格＝左右鄰居的一半＋源", wait=2.0)
    s.play(FadeIn(bars), FadeIn(tcount), run_time=0.6)
    s.say("一步一步放人：長條是第 ", "$t$", " 步時每格的平均人數", wait=0.3)
    for t in FILL_TS[1:]:
        s.play(Transform(bars, tent_bars(hist[t], ux, uy, org)[0]),
               Transform(tcount, at(tx(r"t=%d" % t if t < 200 else r"t\to\infty", 30, GREY_A), -1.6, 2.5)),
               run_time=0.9)
        s.wait(0.2)
    s.play(FadeOut(bars), FadeIn(final[0]), FadeIn(final[1]), run_time=0.8)
    s.say("穩定態：", "$w=0,\\ \\tfrac65,\\ \\tfrac{12}5,\\ \\tfrac85,\\ \\tfrac45,\\ 0$", wait=3.0)
    s.play(ShowCreation(tent), run_time=1.2)
    s.say("頂點連起來是兩段直線：一頂帳篷", wait=2.5)
    s.play(FadeIn(r2, shift=0.2 * UP), run_time=1.0)
    s.say("連續的 (4.49)：兩邊各自線性、在牆上是 0", wait=2.5)
    s.play(FadeIn(r3), run_time=0.6)
    s.say("一維每步跳一格：", "$\\langle x^2\\rangle=N=2\\beta N$", "，所以 ", "$\\beta=\\tfrac12$", wait=3.0)
    s.say("代進去，格子上剛好精確：和長條一模一樣", wait=2.5)
    s.play(FadeIn(r4), run_time=0.6)
    s.play(FadeIn(VGroup(r5, r5b)), run_time=0.6)
    s.say("點電荷讓電場跳一下；這裡源讓斜率跳 ", "$-s_0/\\beta=-2$", wait=3.5)
    s.clear_stage()


def part_ruin(s):
    nb = numbers()
    w = nb["w"]
    u, org = 1.4, np.array([-3.5, 1.05, 0])
    cor = corridor(u, org)
    wl = VGroup(*[tx(r"w=%s" % frac_tex(w[x]), 28, TENT_COLOR).move_to(cpt(u, org, x) + 0.85 * UP)
                  for x in range(1, LW)])
    src = Arrow(cpt(u, org, X0) + 2.05 * UP, cpt(u, org, X0) + 1.25 * UP, buff=0).set_color(SRC_COLOR)
    src_lab = tx("1", 28, SRC_COLOR).next_to(src, RIGHT, buff=0.1)
    a_l = CurvedArrow(cpt(u, org, 1) + 0.2 * UP, cpt(u, org, 0) + 0.2 * UP, angle=TAU / 4).set_color(CLIFF)
    a_r = CurvedArrow(cpt(u, org, LW - 1) + 0.2 * UP, cpt(u, org, LW) + 0.2 * UP, angle=-TAU / 4).set_color(CLIFF)
    l_l = tx(r"\frac{w(1)}{2}=%s" % frac_tex(nb["flux"][0]), 30, CLIFF).next_to(cpt(u, org, 0), LEFT, buff=0.35)
    l_r = tx(r"\frac{w(4)}{2}=%s" % frac_tex(nb["flux"][1]), 30, CLIFF).next_to(cpt(u, org, LW), RIGHT, buff=0.35)
    l_l.shift(0.45 * UP)
    l_r.shift(0.45 * UP)
    bal = at(mixed("每步：進來 1 人＝出去 ", "$%s+%s$" % (frac_tex(nb["flux"][0]), frac_tex(nb["flux"][1])),
                   font_size=28), 0, -0.25)
    ex = at(mixed("Ex. 4.6（賭徒破產）：從 2 出發，先掉進左邊懸崖的機率 ",
                  "$\\frac{L-x_0}{L}=%s$" % frac_tex(nb["left"][5]), font_size=28), 0, -1.2)
    fit(ex, 13.0)
    rows = [["L"] + [str(L) for L in LS] + [r"\infty"],
            [r"P_{\text{left}}"] + [frac_tex(nb["left"][L]) for L in LS] + ["1"]]
    tab = D.tex_table(rows, [1.6, 1.3, 1.3, 1.3, 1.3], row_h=0.55, font_size=30)
    at(tab, 0, -2.35)
    top = VGroup(cor, wl, src, src_lab, a_l, a_r, l_l, l_r)
    assert_in_frame(top, bal, ex, tab)
    assert_table(tab, 5)
    for a, b in ((top, bal), (bal, ex), (ex, tab), (l_l, cor[3]), (l_r, cor[3])):
        assert_apart(a, b)
    assert_apart(wl, src)
    assert_apart(wl, VGroup(a_l, a_r))
    s.play(FadeIn(cor), FadeIn(wl), GrowArrow(src), FadeIn(src_lab), run_time=1.0)
    s.say("每一步，牆邊那格有一半的人跨出去", wait=2.0)
    s.play(ShowCreation(a_l), ShowCreation(a_r), FadeIn(l_l), FadeIn(l_r), run_time=1.0)
    s.say("流進左牆 ", "$\\tfrac35$", "、右牆 ", "$\\tfrac25$", wait=2.5)
    s.play(FadeIn(bal), run_time=0.7)
    s.say("進來的和出去的剛好平衡：這才叫穩定態", wait=2.5)
    s.play(FadeIn(ex, shift=0.2 * UP), run_time=0.8)
    s.say("每個人最後都掉下去；掉在左邊的比例就是 ", "$\\tfrac35$", wait=3.0)
    s.say("三種算法都得 0.6：公式、解方程、逐步遞迴", wait=3.0)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.0)
    s.say("$L\\to\\infty$", "：機率 ", "$\\to1$", "，一維走不到無窮遠", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("源源不絕＋亂走＝靜電", wait=2.5)
    s.say("穩定態的濃度，就是 Poisson 方程的解", wait=3.0)


class W405SteadyTent(ChuScene):
    def construct(self):
        self.show_header("W405　穩定態＝靜電學：兩牆之間的帳篷")
        check_math()
        part_question(self)
        part_poisson(self)
        part_fick(self)
        part_tent(self)
        part_ruin(self)
        part_close(self)
        self.wait(1.0)
