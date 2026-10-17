"""
W507 等待＋步伐＝高斯；會消失的走路人；第 5 章回顧（§5.3，pp. 122–125，(5.124)–(5.133)，Ex. 5.9）

問：把等待和步伐一起算，走很久以後分布是什麼？時間變 4 倍，寬度變幾倍？
  1. (5.124)–(5.128)：Σ χ^n A^n B = B/(1 − χA)（V.ctrw_transform = V.ex58）。
  2. (5.129)–(5.133)：χ ≈ 1 − Cq²、A ≈ 1 − λB(0) → 變異數 2Ct/B(0)；P(0, 10) 長條＋高斯（V.ctrw_lattice、V.ctrw_gauss）；
     筆誤框 width ∝ t → variance ∝ t；變異數 10.00、40.00（寬度 3.162 → 6.325）、4.75、19.75（V.ctrw_var）。
  3. Ex. 5.9：p0 = 1、p1 = 0.1 → 總數 e^{−p1 t}；t = 10：0.3679（三種算法），半衰期 6.931。
  4. 第 5 章回顧表。
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
PATH_SEEDS = (11, 12, 16, 17, 19, 20, 22)
PATH_COLORS = [TEAL, BLUE_C, PURPLE_B, GOLD, PINK, GREEN_C, MAROON_B]
XB = range(-10, 11)
P1 = 0.1
T10, T40 = 10.0, 40.0


def f4(x):
    return f"{x:.4f}"


# ---------------------------------------------------------------- 數學
def ctrw_path(seed, t=T10):
    """示意：固定速率等待（Exp(1)）＋每次左右一格 → 階梯折點 [(時間, 位置)]。"""
    rng = np.random.default_rng(seed)
    s, x, pts = 0.0, 0, [(0.0, 0)]
    while True:
        s += rng.exponential(1.0)
        if s > t:
            pts.append((t, x))
            return pts
        pts.append((s, x))
        x += int(rng.choice([-1, 1]))
        pts.append((s, x))


@lru_cache(maxsize=None)
def numbers():
    nb = {"tr": (V.ctrw_transform(0.7, 1.0), V.ex58(0.7, 1.0))}
    nb["P"] = {k: {x: V.ctrw_lattice(x, T10, k) for x in XB} for k in KINDS}
    nb["G0"] = {k: float(V.ctrw_gauss(0, T10, k)) for k in KINDS}
    nb["B0"] = {k: V.mean_wait(k) for k in KINDS}
    nb["var"] = {k: (V.ctrw_var(T10, k), V.ctrw_var(T40, k)) for k in KINDS}
    nb["mort"] = (V.mortal_total(T10, 1.0, P1), V.mortal_series(T10, 1.0, P1), V.mortal_laplace(T10, 1.0, P1))
    nb["half"] = math.log(2) / P1
    return nb


def check_math():
    nb = numbers()
    assert f4(nb["tr"][0]) == f4(nb["tr"][1]) == "0.8623" and abs(nb["tr"][0] - nb["tr"][1]) < 1e-14
    assert abs(V.ctrw_transform(0.3, 2.0, "gamma2") - V.AB(0.3, "gamma2")[1] / (1 - math.cos(2.0) * V.AB(0.3, "gamma2")[0])) < 1e-15
    assert (f4(nb["P"]["poisson"][0]), f4(nb["G0"]["poisson"])) == ("0.1278", "0.1262")
    assert (f4(nb["P"]["gamma2"][0]), f4(nb["G0"]["gamma2"])) == ("0.1805", "0.1784")
    assert all(abs(nb["P"]["poisson"][x] - V.ctrw_bessel(x, T10)) < 1e-12 for x in XB)
    assert all(sum(nb["P"][k].values()) > 0.99 for k in KINDS)                      # |x| ≤ 10 幾乎是全部
    assert (nb["B0"]["poisson"], nb["B0"]["gamma2"]) == (1, 2)
    assert [f"{v:.2f}" for v in nb["var"]["poisson"]] == ["10.00", "40.00"]
    assert [f"{v:.2f}" for v in nb["var"]["gamma2"]] == ["4.75", "19.75"]
    assert abs(V.ctrw_var(T10, "gamma2") - V.ctrw_var_numeric(T10, "gamma2", X=60)) < 1e-8
    assert abs(V.ctrw_var(T10, "poisson") - V.ctrw_var_numeric(T10, "poisson", X=60)) < 1e-8
    for k in KINDS:                                                                   # 2Ct/B(0)，C = ½（大 t）
        assert abs(nb["var"][k][1] - 0.5 * 2 * T40 / nb["B0"][k]) <= 0.25 + 1e-12
    sig = [math.sqrt(v) for v in nb["var"]["poisson"]]
    assert [f"{x:.3f}" for x in sig] == ["3.162", "6.325"] and abs(sig[1] / sig[0] - 2) < 1e-12
    assert all(f4(v) == "0.3679" for v in nb["mort"]) and abs(nb["mort"][0] - nb["mort"][1]) < 1e-13
    assert abs(nb["mort"][0] - nb["mort"][2]) < 1e-9 and f"{nb['half']:.3f}" == "6.931"
    assert all(max(abs(x) for _, x in ctrw_path(sd)) <= 6 for sd in PATH_SEEDS)


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


def width_arrow(ax, sig, h, color):
    """±σ 的雙向箭頭（兩支從中間往外的箭頭）。"""
    c = ax.c2p(0, h)
    return VGroup(Arrow(c, ax.c2p(-sig, h), buff=0, thickness=2.5),
                  Arrow(c, ax.c2p(sig, h), buff=0, thickness=2.5)).set_color(color)


# ---------------------------------------------------------------- 各段
def part_question(s):
    ax = D.plot_axes((0, 10, 1), (-7, 7, 1), width=6.2, height=4.8, x_ticks=[2, 4, 6, 8, 10],
                     y_ticks=[-6, -3, 3, 6], font_size=20, x_label="t", y_label="x")
    ax.move_to(np.array([-3.3, -0.15, 0]))
    paths = VGroup(*[VMobject().set_points_as_corners([ax.c2p(t, x) for t, x in ctrw_path(sd)]).set_stroke(c, 2.5)
                     for sd, c in zip(PATH_SEEDS, PATH_COLORS)])
    rows = VGroup(mixed("等待：", "$A(\\lambda),\\ B(\\lambda)$", "（W506）", font_size=28),
                  mixed("步伐：", "$\\chi(q)$", "（一維格子 ", "$\\cos q$", "）", font_size=28),
                  zh("合起來，走很久以後？", font_size=28, color=YELLOW),
                  mixed("$t\\to4t$", "：寬度 ", "$\\times\\,?$", font_size=30, color=YELLOW)).arrange(DOWN, buff=0.5,
                                                                                              aligned_edge=LEFT)
    at(rows, 3.4, 0.0)
    assert_in_frame(ax, paths, rows)
    assert_apart(VGroup(ax, paths), rows)
    s.play(ShowCreation(ax), run_time=0.8)
    s.say("問：每一步先等（固定速率），再往左或右走一格", wait=0.3)
    s.play(*[ShowCreation(p) for p in paths], run_time=3.0)
    s.play(FadeIn(rows[:2]), run_time=0.8)
    s.say("W506 算了等待；現在把步伐一起放進來", wait=2.5)
    s.play(FadeIn(rows[2:]), run_time=0.8)
    s.say("走很久以後分布是什麼？時間變 4 倍，寬度變幾倍？（先猜）", wait=3.5)
    s.clear_stage()


def part_transform(s):
    nb = numbers()
    h1 = mixed("走了 ", "$n$", " 步：位移的轉換 ", "$\\chi(q)^n$", "，時間的轉換 ", "$A(\\lambda)^nB(\\lambda)$",
               font_size=30)
    h2 = eq_label(tx(r"\tilde P(q,\lambda)=\sum_{n\ge0}\chi(q)^nA(\lambda)^nB(\lambda)=\frac{B(\lambda)}{1-\chi(q)A(\lambda)}",
                     38), "(5.124)--(5.128)")
    h3 = mixed("固定速率＋一維格子：", "$\\frac{1}{\\lambda+p_0(1-\\cos q)}$", "（W506 的 Ex. 5.8）", font_size=28)
    h4 = mixed("$\\lambda=0.7,\\ q=1$", "：兩邊都是 ", "$%s$" % f4(nb["tr"][0]), font_size=28, color=GREY_A)
    col = VGroup(h1, h2, h3, h4).arrange(DOWN, buff=0.6)
    fit(col, 13.0)
    at(col, 0, 0.2)
    assert_in_frame(col)
    assert_column(*col)
    s.play(FadeIn(h1), run_time=0.8)
    s.say("走了 ", "$n$", " 步：位移的轉換是 ", "$\\chi^n$", "（步伐彼此獨立）", wait=2.5)
    s.say("時間的轉換是 W506 的 ", "$A^nB$", wait=2.0)
    s.play(Write(h2), run_time=1.5)
    s.say("對所有 ", "$n$", " 加起來：幾何級數", wait=2.5)
    s.play(FadeIn(h3), run_time=0.8)
    s.say("固定速率＋左右各一格：就是 Ex. 5.8 的式子", wait=2.5)
    s.play(FadeIn(h4), run_time=0.6)
    s.say("數值核對：", "$\\lambda=0.7,\\ q=1$", " 時兩邊都是 0.8623", wait=2.5)
    s.clear_stage()


def part_gauss(s):
    nb = numbers()
    k1 = eq_label(tx(r"\chi\approx1-Cq^2,\qquad A=1-\lambda B\approx1-\lambda B(0)", 32), "(5.129)(5.130)")
    k2 = eq_label(tx(r"\frac{B}{1-\chi A}\approx\frac{B(0)}{\lambda B(0)+Cq^2}=\frac{1}{\lambda+\frac{C}{B(0)}q^2}", 34),
                  "(5.131)(5.132)")
    k3 = VGroup(eq_label(tx(r"\sigma^2=\frac{2Ct}{B(0)}", 36, YELLOW), "(5.133)"),
                mixed("一維格子 ", "$C=\\frac12$", "，", "$B(0)$", "＝平均等待", font_size=26)).arrange(RIGHT, buff=0.6)
    col = VGroup(k1, k2, k3).arrange(DOWN, buff=0.55)
    at(col, 0, 0.6)
    assert_in_frame(col)
    assert_column(*col)
    s.play(FadeIn(k1), run_time=0.8)
    s.say("小 ", "$q$", "、小 ", "$\\lambda$", "（大距離、長時間）：展開到最低階", wait=2.5)
    s.play(FadeIn(k2), run_time=0.8)
    s.say("$\\frac{1}{\\lambda+Dq^2}$", " 正是擴散方程：", "$D=C/B(0)$", wait=2.5)
    s.play(FadeIn(k3), run_time=0.8)
    s.say("高斯的變異數 ", "$2Ct/B(0)$", "：等得越久，散得越慢", wait=2.5)
    s.play(FadeOut(k1), FadeOut(k2), k3.animate.move_to(np.array([0, 2.6, 0])), run_time=0.8)
    # 長條＋高斯
    ux, uy = 0.26, 11.0
    org = {"poisson": np.array([-3.5, -2.25, 0]), "gamma2": np.array([3.5, -2.25, 0])}
    grp = {}
    for k in KINDS:
        chart = D.bar_chart(nb["P"][k], ux=ux, uy=uy, origin=org[k], color=KCOL[k], width=0.8)
        g = D.curve(lambda x, k=k: float(V.ctrw_gauss(x, T10, k)), -10.5, 10.5, ux, uy, org[k], color=WHITE,
                    stroke_width=3)
        xl = VGroup(*[tx(str(x), 18, GREY_B).move_to(org[k] + x * ux * RIGHT + 0.22 * DOWN) for x in (-10, -5, 0, 5, 10)])
        title = mixed(KNAME[k], "（", "$B(0)=%d$" % nb["B0"][k], "）", font_size=28, color=KCOL[k])
        line = mixed("$P(0,10)=%s$" % f4(nb["P"][k][0]), "，高斯 ", "$%s$" % f4(nb["G0"][k]), font_size=26)
        txt = VGroup(title, line).arrange(DOWN, buff=0.22)
        at(txt, org[k][0], 1.25)
        grp[k] = (chart, g, xl, txt)
        assert_in_frame(chart, g, xl, txt)
        assert_apart(txt, VGroup(chart, g))
        assert_apart(txt, k3)
    assert_apart(VGroup(*grp["poisson"]), VGroup(*grp["gamma2"]))
    for k in KINDS:
        chart, g, xl, txt = grp[k]
        s.play(FadeIn(txt[0]), FadeIn(xl), LaggedStart(*[GrowFromEdge(b, DOWN) for b in chart[0]], lag_ratio=0.04),
               run_time=1.5)
    s.say("$t=10$", " 的精確分布（長條）", wait=1.5)
    s.play(*[ShowCreation(grp[k][1]) for k in KINDS], run_time=1.2)
    s.say("白線：高斯 (5.133)，已經很接近", wait=2.0)
    s.play(*[FadeIn(grp[k][3][1]) for k in KINDS], run_time=0.6)
    s.say("原點：0.1278 對 0.1262；0.1805 對 0.1784", wait=3.0)
    s.say("兩段等待平均等 2：散得慢，原點比較高", wait=2.5)
    s.clear_stage()


def part_width(s):
    nb = numbers()
    ax = D.plot_axes((-15, 15, 5), (0, 0.14, 0.02), width=6.2, height=3.6, x_ticks=[-10, -5, 0, 5, 10],
                     font_size=20, x_label="x")
    ax.y_axis.set_opacity(0)                      # 縱軸在 x = 0 會穿過 σ 標籤；高度不重要，藏起來
    ax.move_to(np.array([-3.3, -0.6, 0]))
    sig = [math.sqrt(v) for v in nb["var"]["poisson"]]
    cols = (YELLOW, TEAL_C)
    curves = VGroup(*[ax.get_graph(lambda x, t=t: float(V.ctrw_gauss(x, t, "poisson")), x_range=(-15, 15, 0.1))
                      .set_stroke(c, 3.5) for t, c in zip((T10, T40), cols)])
    hs = [float(V.ctrw_gauss(sg, t, "poisson")) for sg, t in zip(sig, (T10, T40))]
    arrows = VGroup(*[width_arrow(ax, sg, h, c) for sg, h, c in zip(sig, hs, cols)])
    alabs = VGroup(*[tx(r"\sigma=%.3f" % sg, 22, c).next_to(a, UP, buff=0.06) for sg, a, c in zip(sig, arrows, cols)])
    clabs = VGroup(tx(r"t=10", 24, cols[0]).move_to(ax.c2p(7.5, 0.11)), tx(r"t=40", 24, cols[1]).move_to(ax.c2p(11, 0.045)))
    box = D.erratum_box("p. 124，(5.133) 後", r"\text{width}\propto t",
                        r"\text{variance}\propto t,\ \ \text{width}\propto\sqrt{t}", width=6.2)
    hz = lambda t, c=YELLOW: zh(t, font_size=24, color=c)
    (p10, p40), (g10, g40) = nb["var"]["poisson"], nb["var"]["gamma2"]
    rows = [[hz("變異數"), "t=10", "t=40"],
            [hz(KNAME["poisson"], KCOL["poisson"]), f"{p10:.2f}", f"{p40:.2f}"],
            [hz(KNAME["gamma2"], KCOL["gamma2"]), f"{g10:.2f}", f"{g40:.2f}"]]
    tab = D.tex_table(rows, [2.0, 1.6, 1.6], row_h=0.6, font_size=28)
    wl = mixed("寬度 ", "$%.3f\\to%.3f$" % tuple(sig), "：時間 ", "$\\times4$", "，寬度 ", "$\\times2$", font_size=26)
    right = VGroup(box, tab, wl).arrange(DOWN, buff=0.45)
    fit(right, 6.3)
    at(right, 3.45, 0.0)
    left = VGroup(ax, curves, arrows, alabs, clabs)
    assert_in_frame(left, right)
    assert_apart(left, right)
    assert_column(*right)
    assert_table(tab, 3)
    assert_apart(alabs[1], arrows[0])
    s.play(ShowCreation(ax), FadeIn(box), run_time=1.0)
    s.say("書上說「寬度正比於 ", "$t$", "」：應該是變異數正比於 ", "$t$", wait=3.0)
    s.play(FadeIn(tab[0][:3]), ShowCreation(tab[1]), FadeIn(tab[0][3:]), run_time=0.8)
    s.say("固定速率：變異數 10.00、40.00，正好是 ", "$t$", wait=2.5)
    s.say("兩段等待：4.75、19.75，大約是 ", "$t/2$", wait=2.5)
    s.play(ShowCreation(curves[0]), GrowFromCenter(arrows[0]), FadeIn(alabs[0]), FadeIn(clabs[0]), run_time=1.0)
    s.play(ShowCreation(curves[1]), GrowFromCenter(arrows[1]), FadeIn(alabs[1]), FadeIn(clabs[1]), run_time=1.0)
    s.play(FadeIn(wl), run_time=0.6)
    s.say("揭曉：時間 ", "$\\times4$", "，寬度只 ", "$\\times2$", "（", "$\\sqrt t$", " 定律）", wait=3.0)
    s.clear_stage()


def part_mortal(s):
    nb = numbers()
    m1 = mixed("Ex. 5.9：照走 ", "$p_0=1$", "，每單位時間以 ", "$p_1=0.1$", " 消失", font_size=28)
    m2 = tx(r"A=\frac{p_0}{p_0+p_1+\lambda},\quad B=\frac{1}{p_0+p_1+\lambda}", 32)
    m3 = tx(r"\frac{B}{1-\chi A}=\frac{1}{\lambda+p_1+p_0(1-\chi)}", 32)
    m4 = mixed("$q=0$", "（", "$\\chi=1$", "）：", "$\\frac{1}{\\lambda+p_1}\\ \\Rightarrow\\ e^{-p_1t}$", font_size=28)
    m5 = mixed("$t=10$", "：", "$e^{-1}=%s$" % f4(nb["mort"][0]), "（三種算法）", font_size=26)
    col = VGroup(m1, m2, m3, m4, m5).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    fit(col, 6.6)
    at(col, -3.25, 0.1)
    ax = D.plot_axes((0, 30, 5), (0, 1.1, 0.25), width=5.6, height=3.8, x_ticks=[10, 20, 30], y_ticks=[0.5, 1.0],
                     font_size=20, x_label="t")
    ax.move_to(np.array([3.6, -0.3, 0]))
    cur = ax.get_graph(lambda t: V.mortal_total(t, 1.0, P1), x_range=(0, 30, 0.1)).set_stroke(YELLOW, 4)
    dot = Dot(ax.c2p(T10, nb["mort"][0]), radius=0.07).set_color(RED)
    dl = tx(f4(nb["mort"][0]), 22, RED).next_to(dot, UR, buff=0.08)
    th = nb["half"]
    hl = VGroup(DashedLine(ax.c2p(0, 0.5), ax.c2p(th, 0.5)), DashedLine(ax.c2p(th, 0.5), ax.c2p(th, 0))).set_stroke(GREY_A, 2)
    hlab = tx(r"t_{1/2}=%.3f" % th, 22, GREY_A).move_to(ax.c2p(3.6, 0.41))
    ylab = zh("還在的人數", font_size=22, color=GREY_A).next_to(ax.c2p(0, 1.1), RIGHT, buff=0.2)
    right = VGroup(ax, cur, dot, dl, hl, hlab, ylab)
    assert_in_frame(col, right)
    assert_apart(col, right)
    assert_column(*col)
    s.play(FadeIn(m1), run_time=0.8)
    s.say("會消失的走路人：被吃掉、衰變、被吸收……", wait=2.5)
    s.play(FadeIn(m2), run_time=0.8)
    s.say("等待中也可能消失：", "$A$", "、", "$B$", " 的分母多了 ", "$p_1$", wait=2.5)
    s.play(FadeIn(m3), run_time=0.8)
    s.say("代入 ", "$\\frac{B}{1-\\chi A}$", "：", "$\\lambda$", " 換成 ", "$\\lambda+p_1$", wait=2.5)
    s.play(FadeIn(m4), ShowCreation(ax), FadeIn(ylab), run_time=0.8)
    s.play(ShowCreation(cur), run_time=1.2)
    s.say("總數（", "$q=0$", "）：", "$e^{-p_1t}$", "，衰減率就是 ", "$p_1$", wait=2.5)
    s.play(FadeIn(m5), FadeIn(dot), FadeIn(dl), run_time=0.6)
    s.say("$t=10$", "：0.3679＝逐步數相加＝Laplace 反轉", wait=2.5)
    s.play(ShowCreation(hl), FadeIn(hlab), run_time=0.8)
    s.say("半衰期 ", "$\\ln2/p_1=6.931$", "：和怎麼走完全無關", wait=2.5)
    s.clear_stage()


def part_review(s):
    hz = lambda t, c=WHITE, fs=26: zh(t, font_size=fs, color=c)
    mx = lambda *p: mixed(*p, font_size=26)
    sc = lambda t: tx(r"\text{" + t + "}", 26, GREY_A)
    rows = [[hz("改了什麼", YELLOW), hz("結果", YELLOW), hz("幕", YELLOW)],
            [hz("偏向"), mx("移動的高斯，寬度 ", "$\\times\\sqrt{4p(1-p)}$"), sc("W501")],
            [hz("慣性"), mx("持續長度 ", "$l_p=\\tfrac{Q}{1-Q}\\,l$"), sc("W502--W503")],
            [hz("電報"), mx("尖峰 ", "$e^{-\\gamma t}$", "＋", "$D=\\tfrac{v^2}{2\\gamma}$", " 的高斯"), sc("W504")],
            [hz("二三維"), hz("前緣、因果性、守恆"), sc("W505")],
            [hz("等待"), mx("$\\tfrac{B}{1-\\chi A}$", "，變異數 ", "$\\tfrac{2Ct}{B(0)}$"), sc("W506--W507")]]
    tab = D.tex_table(rows, [2.0, 6.8, 2.8], row_h=0.9, font_size=26)
    at(tab, 0, 0.15)
    assert_in_frame(tab)
    assert_table(tab, 3)
    s.play(FadeIn(tab[0][:3]), ShowCreation(tab[1]), run_time=0.8)
    s.say("第 5 章回顧：改了走路的規則", wait=1.5)
    for i in range(1, 6):
        s.play(FadeIn(tab[0][3 * i:3 * i + 3], shift=0.15 * UP), run_time=0.6)
        s.wait(1.6)
    s.say("改了規則，走路人最後還是回到高斯", wait=2.5)
    s.say("只是平均、寬度與時間尺度換了樣子", wait=3.0)


class W507CtrwGauss(ChuScene):
    def construct(self):
        self.show_header("W507　等待＋步伐＝高斯；會消失的走路人")
        check_math()
        part_question(self)
        part_transform(self)
        part_gauss(self)
        part_width(self)
        part_mortal(self)
        part_review(self)
        self.wait(1.0)
