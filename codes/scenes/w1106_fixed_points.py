"""
W1106 不動點與臨界指數（§11.3，pp. 277–284，(11.88)–(11.114)，Fig. 11.9–11.10）

問：參數一直流——會流到哪裡停下來？
  1. 區塊自旋 (11.88)–(11.92)：κ → R(κ)、dκ/dℓ = W(κ)、dξ/dℓ = −ξ；臨界曲面映到自己；不動點 W(κ*) = 0（Fig. 11.9 示意）。
  2. 線性化 (11.93)–(11.96)：示意矩陣 A = [[1, 2], [2, −2]]：Λ = 2（(2, 1)，相關）、−3（(1, −2)，無關）；流向圖 5 條；註 6。
  3. (11.97)–(11.103)：K r e^{Λ_r ℓ*} = Dist ⇒ ξ = (K r/Dist)^{−1/Λ_r} ⇒ ν = 1/Λ_r（Fig. 11.10 示意）；r = 0.01 的數值。
  4. (11.104)–(11.108)：指數由 Λ_r、Λ_h 決定；6 列 × 5 欄表；後三欄＝W903；α + 2β + γ = 2；W1103 的鏈 Λ = 2、3/2。
  5. 第 11 章回顧（六列表）；收束。
畫面上的數字都由 utils.walk_rg（對照用 utils.walk_crit.EXPONENTS）算出，check_math() 核對。
模組層級的畫圖函式：fig119（臨界曲面示意）、toy_axes／toy_curve（流向圖）、dist_axes（Fig. 11.10）、check_mark。
"""
import math
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rg as RG
from utils import walk_crit as WC
from utils import walk_draw as D

OK_COLOR = "#7CE38B"
REL_COLOR = ORANGE          # 相關方向
IRR_COLOR = TEAL            # 無關方向（臨界線）
TOY_COLORS = (TEAL, YELLOW, BLUE_B, ORANGE, PINK)
TOY_X = (-3.0, 3.0)
TOY_Y = (-2.5, 1.5)
TOY_U = 0.85
CASES = ("gauss3", "gauss4", "ising2d", "saw2d")
WC_OF = {"gauss4": "mf", "ising2d": "ising2d", "saw2d": "saw2d"}
EXPS = ("alpha", "beta", "gamma", "delta", "nu")
EXP_TEX = {"alpha": r"\alpha", "beta": r"\beta", "gamma": r"\gamma", "delta": r"\delta", "nu": r"\nu"}
XI_LR = (2, 1)
DIST_N = 40


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\tfrac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def numbers():
    return {"eig": RG.toy_eigen(),
            "lines": [RG.toy_flow_line(c1, c2) for c1, c2 in RG.TOY_STARTS],
            "xi": {lr: RG.xi_from_flow(RG.XI_R, lr) for lr in XI_LR},
            "xi_num": {lr: RG.xi_from_flow_numeric(RG.XI_R, lr) for lr in XI_LR},
            "exp": {k: RG.rg_exponents(*RG.RG_CASES[k]) for k in CASES},
            "chain": RG.chain_lambdas()}


def check_math():
    nb = numbers()
    # §2 示意矩陣
    assert RG.TOY_A == ((1, 2), (2, -2)) and all(isinstance(x, Fraction) for row in RG.TOY_A for x in row)
    (l1, v1), (l2, v2) = nb["eig"]
    assert (l1, v1, l2, v2) == (2, (2, 1), -3, (1, -2))
    for lam, v in nb["eig"]:                                           # A v = Λ v
        assert all(sum(RG.TOY_A[i][j] * v[j] for j in range(2)) == lam * v[i] for i in range(2))
    assert len(RG.TOY_STARTS) == 5 and RG.TOY_STARTS[0][0] == 0
    for (c1, c2), pts in zip(RG.TOY_STARTS, nb["lines"]):
        assert len(pts) == 30
        assert all(TOY_X[0] <= x <= TOY_X[1] and TOY_Y[0] <= y <= TOY_Y[1] for x, y in pts)
        assert all(abs(a - b) < 1e-9 for a, b in zip(pts[-1], RG.toy_flow_expm(c1, c2, 1.0)))
    first = nb["lines"][0]
    assert all(abs(2 * x + y) < 1e-12 for x, y in first) and math.hypot(*first[-1]) < 0.12 * math.hypot(*first[0])
    for pts in nb["lines"][1:]:                                         # 先靠近、再沿 (2, 1) 被推開
        dist = [math.hypot(x, y) for x, y in pts]
        assert min(dist) < dist[0] or dist[-1] > dist[0]
        x, y = pts[-1]
        assert abs(2 * x + y) > 3 * abs(x - 2 * y)                     # 終點以 (2, 1) 方向為主
    # §3 (11.97)–(11.103)，K = Dist = 1
    assert RG.XI_R == 0.01 and RG.xi_from_flow.__defaults__ == (1.0, 1.0)
    assert [f4(v) for v in nb["xi"][2]] == ["2.3026", "10.0000"] and [f4(v) for v in nb["xi"][1]] == ["4.6052", "100.0000"]
    for lr in XI_LR:
        assert all(abs(a - b) < 1e-6 * b for a, b in zip(nb["xi"][lr], nb["xi_num"][lr]))
        assert abs(nb["xi"][lr][1] - RG.XI_R ** (-1 / lr)) < 1e-9
    # §4 指數表
    F = Fraction
    want = {"gauss3": (F(1, 2), F(1, 4), 1, 5, F(1, 2)), "gauss4": (0, F(1, 2), 1, 3, F(1, 2)),
            "ising2d": (0, F(1, 8), F(7, 4), 15, 1), "saw2d": (F(1, 2), F(5, 64), F(43, 32), F(91, 5), F(3, 4))}
    assert RG.RG_CASES == {"gauss3": (2, F(5, 2), 3), "gauss4": (2, 3, 4), "ising2d": (1, F(15, 8), 2),
                           "saw2d": (F(4, 3), F(91, 48), 2)}
    for k in CASES:
        e = nb["exp"][k]
        assert tuple(e[x] for x in EXPS) == want[k], k
        assert e["alpha"] + 2 * e["beta"] + e["gamma"] == 2
        assert e == RG.rg_exponents_widom(*RG.RG_CASES[k])
        if k in WC_OF:
            assert all(e[x] == WC.EXPONENTS[WC_OF[k]][x] for x in EXPS)
    assert [f4(v) for v in nb["chain"][:2]] == ["2.0000", "1.5000"] and RG.gauss_lambdas(1) == (2, F(3, 2))
    # §5 回顧
    assert f4(RG.LN4_LN3) == "1.2619" and f4(RG.LN2_LN3) == "0.6309" and RG.saw_df(2) == F(4, 3)
    assert abs(RG.rw_df(exact=True) - 2) < 0.01
    sc = RG.scaling_combos(4)
    assert sc[0] == sc[1] and sc[2] == sc[3]                           # d = 4：Landau 與 Gaussian 的組合相同
    assert f4(RG.hh_gauss(0.0, 1.0)) == f4(2 ** 1.5) and abs(RG.rr_gauss(1e-6) / 1e-6 - 4) < 1e-5
    assert RG.ising_iterate(1.0, 5)[-1] < 1e-3 and f4(RG.xi_ratio(3.0)) == "0.5000"
    assert abs(RG.x_map(0.3) - 0.6 / 1.09) < 1e-15


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


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


def framed(mob, color, buff=0.25):
    return VGroup(SurroundingRectangle(mob, buff=buff).set_stroke(color, 2), mob)


def check_mark(color=OK_COLOR, size=0.3):
    return VMobject().set_points_as_corners([size * np.array([-0.5, 0.0, 0]), size * np.array([-0.15, -0.4, 0]),
                                             size * np.array([0.55, 0.45, 0])]).set_stroke(color, 5)


def tip_arrow(p_prev, p_end, color, length=0.3):
    u = (p_end - p_prev) / np.linalg.norm(p_end - p_prev)
    return Arrow(p_end - length * u, p_end, buff=0, thickness=4).set_color(color)


# ---------------------------------------------------------------- 畫圖元件
def fig119():
    """Fig. 11.9 示意：臨界曲面（畫成一條曲線）、不動點、曲面上流向不動點、曲面外先靠近再離開。"""
    surf_pts = [np.array(p) for p in ((-6.2, -2.0, 0), (-5.0, -1.0, 0), (-3.6, 0.0, 0), (-2.3, 1.1, 0), (-1.2, 2.4, 0))]
    surf = VMobject().set_points_smoothly(surf_pts).set_stroke(IRR_COLOR, 5)
    fp = Dot(surf_pts[2], radius=0.13).set_color(WHITE)
    fpl = tx(r"\kappa^*", 32).next_to(fp, RIGHT, buff=0.15)
    on = VGroup(Arrow(surf_pts[0] + 0.15 * UR, surf_pts[1] + 0.05 * UR, buff=0.1, thickness=4),
                Arrow(surf_pts[4] + 0.1 * DL, surf_pts[3] + 0.05 * DL, buff=0.1, thickness=4)).set_color(IRR_COLOR)
    offs = VGroup()
    for a, m, b in (((-5.6, 0.6, 0), (-3.9, 0.55, 0), (-2.6, -1.6, 0)), ((-2.0, -0.6, 0), (-3.25, -0.35, 0),
                                                                        (-5.0, 1.9, 0))):
        pts = [np.array(a), np.array(m), np.array(b)]
        path = VMobject().set_points_smoothly(pts).set_stroke(REL_COLOR, 3)
        offs.add(VGroup(path, tip_arrow(path.point_from_proportion(0.95), path.get_end(), REL_COLOR)))
    lab = mixed("臨界曲面 ", "$\\xi=\\infty$", font_size=22, color=IRR_COLOR).next_to(surf_pts[4], LEFT, buff=0.25).shift(0.25 * UP)
    note = zh("Fig. 11.9 示意", font_size=22, color=GREY_B).move_to(np.array([-5.6, 2.8, 0]))
    return VGroup(surf, on, fp, fpl, offs, lab, note)


def toy_axes():
    w, h = (TOY_X[1] - TOY_X[0]) * TOY_U, (TOY_Y[1] - TOY_Y[0]) * TOY_U
    return D.plot_axes((TOY_X[0], TOY_X[1], 1), (TOY_Y[0], TOY_Y[1], 0.5), width=w, height=h,
                       x_label=r"\delta\kappa_1", y_label=r"\delta\kappa_2")


def toy_curve(ax, pts, color):
    ps = [ax.c2p(x, y) for x, y in pts]
    line = VMobject().set_points_smoothly(ps).set_stroke(color, 4)
    return VGroup(line, tip_arrow(ps[-2], ps[-1], color, 0.28))


def eigen_line(ax, v, t0, t1, color):
    return DashedLine(ax.c2p(t0 * v[0], t0 * v[1]), ax.c2p(t1 * v[0], t1 * v[1]), dash_length=0.1).set_stroke(color, 2)


def dist_axes():
    return D.plot_axes((0, 5.2, 1), (0, 1.25, 0.25), width=5.4, height=3.6, x_ticks=(1, 2, 3, 4, 5),
                       y_ticks=(0.5, 1), y_tick_labels=["0.5", "1"], x_label=r"\ell", y_label=r"\delta")


def dist_curve(ax, lr, color):
    """K r e^{Λ_r ℓ}（K = 1、r = 0.01）畫到碰到 Dist = 1 的 ℓ*。"""
    ell_star = RG.xi_from_flow(RG.XI_R, lr)[0]
    ps = [ax.c2p(t, RG.XI_R * math.exp(lr * t)) for t in np.linspace(0, ell_star, DIST_N)]
    return VMobject().set_points_smoothly(ps).set_stroke(color, 4)


# ---------------------------------------------------------------- 各段
def part_question(s):
    fig = fig119()
    eqs = VGroup(tx(r"\kappa_i\to\kappa_i'=R_i(\kappa)", 34),
                 tx(r"\frac{d\kappa_i}{d\ell}=W_i(\kappa),\quad s=e^{\ell}", 34),
                 tx(r"\frac{d\xi}{d\ell}=-\xi", 34),
                 tx(r"W_i(\kappa^*)=0", 34).set_color(YELLOW)).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    tag = zh("(11.88)–(11.92)　區塊自旋", font_size=22, color=GREY_B)
    right = VGroup(tag, eqs).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
    place_in(right, 0.8, 6.7, -2.8, 3.0)
    assert_in_frame(fig, right)
    assert_apart(fig, right)
    assert_column(tag, *eqs)
    s.say("問：參數一直流——會流到哪裡停下來？", wait=1.8)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.0)
    s.say("區塊自旋：每做一次，參數 ", "$\\kappa_i$", " 變成 ", "$R_i(\\kappa)$", wait=2.0)
    s.play(Write(eqs[1]), run_time=1.0)
    s.say("連續版：", "$s=e^{\\ell}$", "，參數沿 ", "$W_i(\\kappa)$", " 流動", wait=2.0)
    s.play(Write(eqs[2]), run_time=0.8)
    s.say("長度單位變大，關聯長度就變短：", "$d\\xi/d\\ell=-\\xi$", wait=2.2)
    s.play(FadeIn(fig[6]), ShowCreation(fig[0]), FadeIn(fig[5]), run_time=1.2)
    s.say("$\\xi=\\infty$", " 的點還是 ", "$\\infty$", "：臨界曲面被映到自己", wait=2.2)
    s.play(FadeIn(fig[2]), FadeIn(fig[3]), GrowArrow(fig[1][0]), GrowArrow(fig[1][1]), Write(eqs[3]), run_time=1.2)
    s.say("曲面上的點流向不動點 ", "$W_i(\\kappa^*)=0$", wait=2.0)
    s.play(*[ShowCreation(o[0]) for o in fig[4]], run_time=1.4)
    s.play(*[GrowArrow(o[1]) for o in fig[4]], run_time=0.5)
    s.say("曲面外的點：先靠近，再被推開", wait=2.0)
    s.clear_stage()


def part_linear(s):
    nb = numbers()
    (l1, v1), (l2, v2) = nb["eig"]
    ax = toy_axes()
    ax.move_to(np.array([-3.4, -0.1, 0]))
    rel = eigen_line(ax, v1, -1.5, 1.5, REL_COLOR)
    irr = eigen_line(ax, v2, -0.75, 1.2, IRR_COLOR)
    curves = VGroup(*[toy_curve(ax, pts, c) for pts, c in zip(nb["lines"], TOY_COLORS)])
    fp = Dot(ax.c2p(0, 0), radius=0.1).set_color(WHITE)
    eq = tx(r"\frac{d\,\delta\kappa_i}{d\ell}=\sum_jA_{ij}\,\delta\kappa_j\ \Rightarrow\ "
            r"\delta\kappa=\sum_\alpha C_\alpha e^{\Lambda_\alpha\ell}\vec v_\alpha", 30)
    tag = zh("(11.93)–(11.96)　在不動點附近線性化", font_size=22, color=GREY_B)
    A = RG.TOY_A
    mat = VGroup(tx(r"A=\begin{pmatrix}%s&%s\\%s&%s\end{pmatrix}" % tuple(ftex(x) for row in A for x in row), 34),
                 zh("（示意）", font_size=24, color=GREY_B)).arrange(RIGHT, buff=0.2)
    e1 = mixed("$\\Lambda=%s$" % ftex(l1), "，", "$\\vec v=(%s,%s)$" % tuple(ftex(x) for x in v1), "：相關",
               font_size=26, color=REL_COLOR)
    e2 = mixed("$\\Lambda=%s$" % ftex(l2), "，", "$\\vec v=(%s,%s)$" % tuple(ftex(x) for x in v2), "：無關",
               font_size=26, color=IRR_COLOR)
    note = framed(VGroup(zh("註 6", font_size=22, color=GREY_A),
                         mixed("$\\Lambda>0$", " 相關、", "$\\Lambda<0$", " 無關、", "$\\Lambda=0$", " 邊際", font_size=24)
                         ).arrange(DOWN, buff=0.15), GREY_C, 0.2)
    right = VGroup(mat, e1, e2, note).arrange(DOWN, buff=0.35)
    top = VGroup(tag, eq).arrange(DOWN, buff=0.2)
    place_in(top, -6.7, 6.7, 2.0, 3.1)
    place_in(right, 0.9, 6.7, -2.85, 1.7)
    assert_in_frame(top, ax, curves, right, rel, irr)
    assert_apart(top, VGroup(ax, curves))
    assert_apart(ax, right)
    assert_apart(top, right)
    assert_column(mat, e1, e2, note)
    s.play(FadeIn(tag), Write(eq), run_time=1.4)
    s.say("在不動點附近，偏離 ", "$\\delta\\kappa$", " 很小：流變成線性", wait=2.0)
    s.say("解是特徵向量的組合，各自乘 ", "$e^{\\Lambda_\\alpha\\ell}$", wait=2.0)
    s.play(FadeIn(mat), run_time=0.8)
    s.say("拿一個示意用的 ", "$2\\times2$", " 矩陣來看", wait=1.6)
    s.play(FadeIn(e1), FadeIn(e2), run_time=0.8)
    s.say("特徵值 ", "$2$", " 與 ", "$-3$", "：一個放大、一個縮小", wait=2.0)
    s.play(ShowCreation(ax), FadeIn(fp), ShowCreation(rel), ShowCreation(irr), run_time=1.2)
    s.play(ShowCreation(curves[0][0]), run_time=1.2)
    s.play(GrowArrow(curves[0][1]), run_time=0.4)
    s.say("$c_1=0$", "：沿著臨界線流進不動點", wait=2.0)
    s.play(*[ShowCreation(c[0]) for c in curves[1:]], run_time=2.0)
    s.play(*[GrowArrow(c[1]) for c in curves[1:]], run_time=0.5)
    s.say("其餘先靠近，再沿 ", "$(2,1)$", " 被推開", wait=2.0)
    s.play(FadeIn(note), run_time=0.8)
    s.say("$\\Lambda>0$", " 相關、", "$\\Lambda<0$", " 無關、", "$\\Lambda=0$", " 邊際", wait=2.2)
    s.clear_stage()


def part_nu(s):
    nb = numbers()
    ax = dist_axes()
    ax.move_to(np.array([-3.4, -0.1, 0]))
    curves = VGroup(dist_curve(ax, 2, YELLOW), dist_curve(ax, 1, ORANGE))
    dline = DashedLine(ax.c2p(0, 1), ax.c2p(5.1, 1), dash_length=0.08).set_stroke(RED, 2)
    dlab = tx(r"\mathrm{Dist}=1", 24, RED).next_to(ax.c2p(5.1, 1), UP, buff=0.08).shift(0.5 * LEFT)
    vl = VGroup(*[DashedLine(ax.c2p(nb["xi"][lr][0], 0), ax.c2p(nb["xi"][lr][0], 1), dash_length=0.06)
                  .set_stroke(GREY_A, 1.5) for lr in XI_LR])
    clab = VGroup(tx(r"\Lambda_r=2", 24, YELLOW).next_to(ax.c2p(math.log(70) / 2, 0.7), LEFT, buff=0.15),
                  tx(r"\Lambda_r=1", 24, ORANGE).next_to(ax.c2p(math.log(50), 0.5), LEFT, buff=0.15))
    note = zh("Fig. 11.10 示意（K = 1、r = 0.01）", font_size=20, color=GREY_B).next_to(ax, UP, buff=0.1)
    tag = zh("(11.97)–(11.103)", font_size=22, color=GREY_B)
    eqs = VGroup(tx(r"\delta(\ell)=Kr\,e^{\Lambda_r\ell},\qquad \xi^{(\ell)}=\xi e^{-\ell}", 30),
                 tx(r"\xi^{(\ell^*)}=1:\ Kr\,e^{\Lambda_r\ell^*}=\mathrm{Dist}", 30),
                 tx(r"\xi=e^{\ell^*}=(Kr/\mathrm{Dist})^{-1/\Lambda_r}", 30),
                 tx(r"\nu=1/\Lambda_r", 40).set_color(YELLOW)).arrange(DOWN, buff=0.3)
    rows = [[r"\Lambda_r", r"\ell^*", r"\xi"]] + [[str(lr), f4(nb["xi"][lr][0]), f4(nb["xi"][lr][1])] for lr in XI_LR]
    tab = D.tex_table(rows, [1.1, 1.6, 1.8], row_h=0.5, font_size=28)
    right = VGroup(tag, eqs, tab).arrange(DOWN, buff=0.35)
    place_in(right, 0.4, 6.7, -2.85, 3.1)
    assert_in_frame(ax, note, right, dlab, clab)
    assert_apart(VGroup(ax, note, dlab), right)
    assert_column(tag, *eqs, tab)
    assert_table(tab, 3)
    s.say("離臨界面有多遠，決定了 ", "$\\xi$", wait=1.6)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.2)
    s.say("偏離沿相關方向長大，", "$\\xi$", " 每步縮小 ", "$e^{-\\ell}$", wait=2.0)
    s.play(ShowCreation(ax), FadeIn(note), ShowCreation(dline), FadeIn(dlab), run_time=1.0)
    s.play(ShowCreation(curves), FadeIn(clab), run_time=1.6)
    s.play(Write(eqs[1]), ShowCreation(vl), run_time=1.0)
    s.say("流到 ", "$\\xi^{(\\ell^*)}=1$", " 時，離臨界面的距離是 ", "$\\mathrm{Dist}$", wait=2.2)
    s.play(Write(eqs[2]), run_time=1.0)
    s.play(Write(eqs[3]), run_time=0.8)
    s.say("所以 ", "$\\xi=(Kr/\\mathrm{Dist})^{-1/\\Lambda_r}$", "：", "$\\nu=1/\\Lambda_r$", wait=2.2)
    s.play(FadeIn(tab), run_time=0.8)
    s.say("$r=0.01$", "：", "$\\Lambda_r=2$", " 得 ", "$\\xi=%s$" % f4(nb["xi"][2][1]),
          "；", "$\\Lambda_r=1$", " 得 ", "$%s$" % f4(nb["xi"][1][1]), wait=2.4)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    tag = zh("(11.104)–(11.108)", font_size=22, color=GREY_B)
    f0 = tx(r"f(r,h)=e^{-\ell d}f\big(re^{\Lambda_r\ell},\ he^{\Lambda_h\ell}\big)", 36)
    fs = VGroup(tx(r"\beta=\frac{d-\Lambda_h}{\Lambda_r}", 34), tx(r"-\gamma=\frac{d-2\Lambda_h}{\Lambda_r}", 34),
                tx(r"2-\alpha=\frac{d}{\Lambda_r}", 34), tx(r"\delta=\frac{\Lambda_h}{d-\Lambda_h}", 34),
                tx(r"\nu=\frac1{\Lambda_r}", 34)).arrange(RIGHT, buff=0.55)
    for m in fs:
        m.set_color(YELLOW)
    col1 = VGroup(tag, f0, fs).arrange(DOWN, buff=0.45)
    place_in(col1, -6.7, 6.7, -1.0, 2.9)
    assert_in_frame(col1)
    assert_column(tag, f0, fs)
    s.say("自由能也照同樣的流標度", wait=1.6)
    s.play(FadeIn(tag), Write(f0), run_time=1.2)
    s.say("取 ", "$e^{\\Lambda_r\\ell}r=1$", "，每個指數就讀出來了", wait=2.0)
    s.play(LaggedStart(*[FadeIn(m) for m in fs], lag_ratio=0.3), run_time=1.6)
    s.say("全部只靠兩個數：", "$\\Lambda_r$", " 和 ", "$\\Lambda_h$", wait=2.0)
    s.play(FadeOut(col1), run_time=0.6)

    heads = [zh("Gaussian d = 3", font_size=22, color=YELLOW), zh("Gaussian d = 4", font_size=22, color=YELLOW),
             zh("二維 Ising", font_size=22, color=YELLOW), zh("二維自迴避", font_size=22, color=YELLOW)]
    rows = [[VectorizedPoint()] + heads] + [[EXP_TEX[x]] + [ftex(nb["exp"][k][x]) for k in CASES] for x in EXPS]
    widths = [1.2, 2.6, 2.6, 2.6, 2.6]
    tab = D.tex_table(rows, widths, row_h=0.56, font_size=30)
    lam = VGroup(*[tx(r"(%s,\ %s)" % (ftex(RG.RG_CASES[k][0]), ftex(RG.RG_CASES[k][1])), 24, GREY_A) for k in CASES])
    lamlab = tx(r"(\Lambda_r,\Lambda_h)", 22, GREY_A)
    xs = np.cumsum([0] + widths)
    for j, m in enumerate(lam):
        m.move_to(np.array([(xs[j + 1] + xs[j + 2]) / 2, 0.3, 0]))
    lamlab.move_to(np.array([xs[1] / 2, 0.3, 0]))
    ticks = VGroup(*[check_mark().move_to(np.array([(xs[j + 1] + xs[j + 2]) / 2, -6 * 0.56 - 0.3, 0]))
                     for j in range(4)])
    ticklab = tx(r"\alpha+2\beta+\gamma=2", 22, OK_COLOR).move_to(np.array([xs[1] / 2, -6 * 0.56 - 0.3, 0]))
    block = VGroup(lamlab, lam, tab, ticks, ticklab)
    w903 = mixed("後三欄＝第 9 章 W903 的指數表", font_size=24, color=OK_COLOR)
    chain = mixed("W1103 的鏈：", "$\\log_24=%s$" % f4(nb["chain"][0]), "、",
                  "$\\log_22^{3/2}=%s$" % f4(nb["chain"][1]), "＝", "$d=1$", " 的 ", "$\\Lambda_r,\\Lambda_h$",
                  font_size=24)
    low = VGroup(w903, chain).arrange(DOWN, buff=0.2)
    col = VGroup(block, low).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_table(tab, 5)
    assert_column(lam, tab, ticks, low)
    assert_apart(lamlab, lam)
    assert_apart(ticklab, ticks)
    s.play(FadeIn(tab[1]), FadeIn(tab[0][:5]), FadeIn(lam), FadeIn(lamlab), run_time=1.0)
    s.say("四組 ", "$(\\Lambda_r,\\Lambda_h)$", "，代進同一組公式", wait=1.8)
    s.play(LaggedStart(*[FadeIn(tab[0][5 * i:5 * i + 5]) for i in range(1, 6)], lag_ratio=0.3), run_time=2.0)
    s.say("$d=4$", " 的 Gaussian 就是平均場：", "$\\beta=\\tfrac12$", "、", "$\\delta=3$", wait=2.2)
    s.say("二維 Ising：", "$\\beta=\\tfrac18$", "、", "$\\gamma=\\tfrac74$", "、", "$\\delta=15$", wait=2.0)
    s.play(ShowCreation(ticks), FadeIn(ticklab), FadeIn(w903), run_time=1.0)
    s.say("後三欄就是第 9 章的指數表；每欄 ", "$\\alpha+2\\beta+\\gamma=2$", wait=2.2)
    s.play(FadeIn(chain), run_time=0.8)
    s.say("W1103 的鏈：", "$\\Lambda_r=2$", "、", "$\\Lambda_h=\\tfrac32$", "，正是 ", "$d=1$", " 的值", wait=2.2)
    s.clear_stage()


def part_recap(s):
    title = zh("第 11 章回顧", font_size=30, color=YELLOW)
    fs = 22
    labels = ["W1101 碎形", "W1102 被動重整化", "W1103 Gaussian 鏈", "W1104 Ising 鏈", "W1105 動量殼", "W1106 不動點"]
    items = [mixed("Koch ", "$\\frac{\\ln4}{\\ln3}=%s$" % f4(RG.LN4_LN3), "、Cantor ", "$%s$" % f4(RG.LN2_LN3),
                   "、漫步 ", "$2$", "、自迴避 ", "$%s$" % ftex(RG.saw_df(2)).replace("tfrac", "frac"), font_size=fs),
             mixed("$m=r^{1/2}f(h/r^{3/2})$", "；", "$d=4$", " 與 Gaussian 一致", font_size=fs),
             mixed("$r''=r(r+4a)/a\\approx4r$", "、", "$h''\\approx2^{3/2}h$", font_size=fs),
             mixed("$K'=\\tfrac12\\ln\\cosh2K\\to0$", "、", "$x'=\\frac{2x}{1+x^2}$", "、", "$\\xi\\propto e^{2K}$",
                   font_size=fs),
             mixed("$\\frac{dr}{d\\ell}=2r$", "、", "$\\frac{dh}{d\\ell}=\\frac{d+2}{2}h$", font_size=fs),
             mixed("$\\nu=1/\\Lambda_r$", "、", "$\\beta=(d-\\Lambda_h)/\\Lambda_r$", font_size=fs)]
    rows = [[zh(t, font_size=fs, color=YELLOW), m] for t, m in zip(labels, items)]
    tab = D.tex_table(rows, [3.2, 9.4], row_h=0.78)
    col = VGroup(title, tab[0]).arrange(DOWN, buff=0.35)
    place_in(col, -6.7, 6.7, -2.85, 3.1)
    assert_in_frame(col)
    assert_table(tab, 2)
    s.play(FadeIn(title), run_time=0.6)
    s.say("第 11 章回顧：六幕的重點", wait=1.4)
    notes = [("碎形：放大 ", "$\\mathcal M$", " 倍、份量 ", "$\\mathcal R$", " 倍，維度是比的對數"),
             ("被動的重整化：只看參數怎麼變",),
             ("Gaussian 鏈：積掉一半，", "$r$", " 每次乘 4"),
             ("Ising 鏈：", "$K$", " 流向 0，一維沒有相變"),
             ("動量殼：重整化變成微分方程",),
             ("不動點：指數由特徵值決定",)]
    for i, note in enumerate(notes):
        s.play(FadeIn(tab[0][2 * i]), FadeIn(tab[0][2 * i + 1]), run_time=0.6)
        s.say(*note, wait=1.6)
    s.say("臨界指數只由不動點附近兩個相關特徵值決定", wait=2.4)
    s.say("——下一章把這套用到 ", "$\\mathrm{O}(n)$", " 模型與自迴避漫步", wait=2.6)


class W1106FixedPoints(ChuScene):
    def construct(self):
        self.show_header("W1106　不動點與臨界指數")
        check_math()
        part_question(self)
        part_linear(self)
        part_nu(self)
        part_table(self)
        part_recap(self)
        self.wait(1.0)
