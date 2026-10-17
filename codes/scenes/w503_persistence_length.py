"""
W503 持續長度（§5.2–5.2.1，pp. 100–103，(5.23)–(5.38)，Ex. 5.2）

問：Q = 3/4 的醉漢走 1000 步，寬度相當於普通醉漢走幾步？
  1. (5.23) λ1 ≈ 1 − Q q²/(2(1−Q))（q = 0.1：0.985 vs 0.9848）→ (5.26)(5.27) 高斯，變異數 NQ/(1−Q)；
     N = 100、x = 0：精確 0.04587 vs 2×(5.27) 0.04607；表：精確 vs 漸近（差 → 4 與 −4/9）。
  2. 持續長度 (5.28)–(5.31)：n = Q/(1−Q)（3、1、1/3），l_p = nl；揭曉：1000 步 ≈ 普通醉漢 3000 步。
  3. 兩個極端：Q = 1 直線（變異數 N²）、Q = 0 原地來回（0、1）。
  4. 二階遞迴 (5.34)–(5.38)＋筆誤框 (5.37)；Ex. 5.2：兩個 λ^N 相加回到 W502 的總數。
畫面上的數字都由 utils.walk_var 算出，check_math() 核對。
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
from utils import walk_var as V
from utils import walk_draw as D

QV = 0.75
QS = (Fraction(3, 4), Fraction(1, 2), Fraction(1, 4))
NS = (10, 100, 1000)
N_BAR = 100
BIG_N, BIG_Y, BIG_SEED = 1000, 120, 4
Q_SMALL = 0.1
Q_EX = 0.3
P_COL = YELLOW           # 有慣性
O_COL = BLUE_C           # 普通


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return ("-" if q < 0 else "") + r"\frac{%d}{%d}" % (abs(q.numerator), q.denominator)


def sample_path(n, Qv, seed):
    rng = np.random.default_rng(seed)
    d = 1 if rng.random() < 0.5 else -1
    pos = [0, d]
    for _ in range(n - 1):
        if rng.random() >= Qv:
            d = -d
        pos.append(pos[-1] + d)
    return pos


def ex52(N, q, Q):
    """Ex. 5.2：總數 S(N) = A1 λ1^N + A2 λ2^N，A1、A2 由 S(0) = 1、S(1) = cos q 決定。"""
    l1, l2 = V.eigs(q, Q)
    a1 = (math.cos(q) - l2) / (l1 - l2)
    return a1 * l1 ** N + (1 - a1) * l2 ** N


@lru_cache(maxsize=None)
def numbers():
    nb = {"l1": V.eigs(Q_SMALL, QV)[0], "l1s": V.lambda1_small(Q_SMALL, QV)}
    xs, Lp, Rp = V.persistent_np(N_BAR, QV)
    P = Lp + Rp
    nb["bars"] = {int(x): float(p) for x, p in zip(xs, P) if abs(x) <= 50 and x % 2 == 0}
    i0 = int(np.nonzero(xs == 0)[0][0])
    nb["p0"] = (float(P[i0]), 2 * float(V.gauss_persistent(0, N_BAR, QV)))
    nb["tab"] = {(Q, N): (V.persistent_var(N, float(Q)), V.persistent_var_asym(N, float(Q))) for Q in QS[::2] for N in NS}
    nb["n"] = {Q: V.persistence_length(Q) for Q in QS}
    nb["big"] = (V.persistent_var(BIG_N, QV), V.persistent_var_asym(BIG_N, QV))
    nb["k"] = BIG_N // 3
    nb["ext"] = (V.persistent_var(10, 1.0), V.persistent_var(10, 0.0), V.persistent_var(11, 0.0))
    nb["paths"] = (sample_path(BIG_N, QV, BIG_SEED), sample_path(BIG_N, 0.5, BIG_SEED + 100))
    nb["book"] = V.eig_residual(nb["l1"], Q_SMALL, QV, book=True)
    nb["ex"] = (ex52(5, Q_EX, QV), V.total_fourier(5, Q_EX, QV))
    return nb


def check_math():
    nb = numbers()
    l1 = nb["l1"]
    assert (f4(l1.real), f"{nb['l1s']:.3f}") == ("0.9848", "0.985") and abs(l1 - nb["l1s"]) < 1e-3
    assert 0.99 < sum(nb["bars"].values()) < 1                                            # |x| ≤ 50 幾乎全部
    assert (f"{nb['p0'][0]:.5f}", f"{nb['p0'][1]:.5f}") == ("0.04587", "0.04607")
    assert max(nb["bars"].values()) == nb["bars"][0] and len(nb["bars"]) <= 60
    tab = nb["tab"]
    Q3, Q1 = QS[0], QS[2]
    assert [f"{tab[(Q3, N)][0]:.2f}" for N in NS] == ["26.00", "296.00", "2996.00"]
    assert [round(tab[(Q3, N)][1], 9) for N in NS] == [30, 300, 3000]
    assert [f"{tab[(Q1, N)][0]:.3f}" for N in NS] == ["3.777", "33.778", "333.778"]
    assert [f"{tab[(Q1, N)][1]:.3f}" for N in NS] == ["3.333", "33.333", "333.333"]
    assert abs(tab[(Q3, 1000)][1] - tab[(Q3, 1000)][0] - 4) < 1e-9                        # 差 → 4
    assert abs(tab[(Q1, 1000)][1] - tab[(Q1, 1000)][0] + 4 / 9) < 1e-9                    # 差 → −4/9
    assert [nb["n"][Q] for Q in QS] == [3, 1, Fraction(1, 3)]
    assert nb["n"][Q3] == Q3 / (1 - Q3)
    var, asym = nb["big"]
    assert f"{var:.0f}" == "2996" and round(asym, 9) == 3000 and nb["k"] == 333 and 9 * nb["k"] == 2997
    assert abs(nb["n"][Q3] * BIG_N - asym) < 1e-9 and abs((BIG_N / 3) * 3 ** 2 - 3000) < 1e-9
    assert nb["ext"] == (100, 0, 1)
    a, b = nb["paths"]
    assert max(map(abs, a)) <= BIG_Y and max(map(abs, b)) <= BIG_Y and len(a) == BIG_N + 1
    for N in range(8):                                                                       # (5.35)
        for q in (0.1, 0.3):
            assert V.recursion2_residual(N, q, QV) < 1e-14
    for q in (0.1, 0.5):                                                                     # (5.36) 的特徵值
        M = np.array([[2 * QV * math.cos(q), 1 - 2 * QV], [1, 0]])
        num = np.linalg.eigvals(M)
        assert all(min(abs(num - lam)) < 1e-12 for lam in V.eigs(q, QV))
    assert V.eig_residual(l1, Q_SMALL, QV) < 1e-14 and f"{nb['book']:.3f}" == "1.000"       # (5.37) 書上差 1
    ex, ft = nb["ex"]
    assert abs(ex - ft) < 1e-12 and f4(ex.real) == "0.5699"


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


def rods(n, length, unit, start, color, stroke=5):
    """n 根同方向的箭頭（長 length 格），頭接尾往右排。"""
    out = VGroup()
    for k in range(n):
        p0 = start + k * length * unit * RIGHT
        out.add(Arrow(p0, p0 + length * unit * RIGHT, buff=0.03, thickness=stroke * 0.6).set_color(color))
    return out


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    a, b = nb["paths"]
    ux, uy, org = 0.0108, 0.0125, np.array([-5.6, -0.75, 0])
    axes = D.walk_axes(BIG_N, BIG_Y, ux=ux, uy=uy, origin=org, x_every=250, y_every=60, font_size=18)
    pa = D.walk_path(a, ux=ux, uy=uy, origin=org, color=P_COL, stroke_width=2, dots=False)
    pb = D.walk_path(b, ux=ux, uy=uy, origin=org, color=O_COL, stroke_width=2, dots=False)
    la = mixed("$Q=\\frac34$", font_size=22, color=P_COL).next_to(D.path_point(BIG_N, a[-1], ux, uy, org), RIGHT,
                                                                   buff=0.12)
    lb = mixed("$Q=\\frac12$", font_size=22, color=O_COL).next_to(D.path_point(BIG_N, b[-1], ux, uy, org), RIGHT,
                                                                   buff=0.12)
    ask = at(mixed("1000 步的寬度＝普通醉漢走 ", "$?$", " 步", font_size=30, color=YELLOW), 0, 2.55)
    assert_in_frame(axes, pa, pb, la, lb, ask)
    assert_apart(la, lb)
    assert_apart(ask, VGroup(axes, pa, pb))
    s.play(FadeIn(axes), run_time=0.6)
    s.say("問：", "$Q=\\frac34$", " 的醉漢走 1000 步（黃）", wait=0.3)
    s.play(ShowCreation(pa), run_time=2.5)
    s.play(FadeIn(la), run_time=0.4)
    s.say("普通醉漢也走 1000 步（藍）", wait=0.3)
    s.play(ShowCreation(pb), run_time=2.5)
    s.play(FadeIn(lb), run_time=0.4)
    s.play(FadeIn(ask), run_time=0.6)
    s.say("黃的寬度，相當於普通醉漢走幾步？（先猜）", wait=3.0)
    s.clear_stage()


def part_gauss(s):
    nb = numbers()
    e1 = eq_label(tx(r"\lambda_1\approx1-\frac{Q}{1-Q}\,\frac{q^2}{2}", 32), "(5.23)")
    n1 = mixed("$q=0.1$", "：", "$%.3f$" % nb["l1s"], "（精確 ", "$%s$" % f4(nb["l1"].real), "）", font_size=26)
    e2 = eq_label(tx(r"\lambda_1^N\approx e^{-N\frac{Q}{1-Q}\frac{q^2}{2}}\ \Rightarrow\ "
                     r"P(N,x)\approx\frac{1}{\sqrt{2\pi N\frac{Q}{1-Q}}}\,e^{-\frac{x^2}{2NQ/(1-Q)}}", 30),
                  "(5.26)(5.27)")
    fit(e2, 12.0)
    at(e1, -2.6, 2.6)
    at(n1, 3.3, 2.6)
    at(e2, 0, 1.35)
    ux, uy, org = 0.105, 46.0, np.array([0, -2.5, 0])
    base = Line(org + 52 * ux * LEFT, org + 52 * ux * RIGHT).set_stroke(GREY_C, 2)
    nums = VGroup(*[tx(str(x), 18, GREY_B).move_to(org + x * ux * RIGHT + 0.25 * DOWN) for x in range(-50, 51, 25)])
    chart = D.bar_chart(nb["bars"], ux=ux, uy=uy, origin=org, width=1.5)
    gauss = D.curve(lambda x: 2 * float(V.gauss_persistent(x, N_BAR, QV)), -50, 50, ux=ux, uy=uy, origin=org)
    lab = at(mixed("$N=100,\\ x=0$", "：精確 %.5f；", "$2\\times$", "(5.27) %.5f" % nb["p0"][1], font_size=22), 3.6, 0.15)
    lab[1].become(zh("：精確 %.5f；" % nb["p0"][0], font_size=22).move_to(lab[1]))
    fit(lab, 6.2)
    at(lab, 3.6, 0.15)
    tl = at(tx(r"N=100", 26, GREY_A), -4.6, 0.1)
    plot = VGroup(base, nums, chart, gauss)
    assert_in_frame(e1, n1, e2, plot, lab, tl)
    for a, b in ((e1, n1), (e1, e2), (n1, e2), (e2, plot), (e2, lab), (lab, plot), (tl, plot)):
        assert_apart(a, b)
    s.play(Write(e1), run_time=1.0)
    s.say("W502 的最大特徵值：", "$q$", " 小時展開到 ", "$q^2$", wait=2.9)
    s.play(FadeIn(n1), run_time=0.6)
    s.say("$q=0.1$", "：近似 0.985，精確 0.9848", wait=2.6)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("$\\lambda_1^N$", " 是高斯的傅立葉轉換：變異數 ", "$N\\frac{Q}{1-Q}$", wait=3.0)
    s.play(FadeIn(base), FadeIn(nums), FadeIn(tl), LaggedStart(*[GrowFromEdge(r, DOWN) for r in chart[0]],
                                                                  lag_ratio=0.04), run_time=1.5)
    s.play(ShowCreation(gauss), run_time=1.0)
    s.say("$N=100$", " 的精確分布（長條）與高斯（橘）", wait=2.6)
    s.play(FadeIn(lab), run_time=0.6)
    s.say("原點：0.04587 對 0.04607（奇偶性的因子 2）", wait=2.9)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    tab_v = nb["tab"]
    Q3, Q1 = QS[0], QS[2]
    rows = [[VectorizedPoint()] + [r"N=%d" % N for N in NS],
            [mixed("$Q=\\frac34$", " 精確", font_size=24)] + [f"{tab_v[(Q3, N)][0]:.2f}" for N in NS],
            [mixed("$N\\frac{Q}{1-Q}$", font_size=24)] + ["%d" % round(tab_v[(Q3, N)][1]) for N in NS],
            [mixed("$Q=\\frac14$", " 精確", font_size=24)] + [f"{tab_v[(Q1, N)][0]:.3f}" for N in NS],
            [mixed("$N\\frac{Q}{1-Q}$", font_size=24)] + [f"{tab_v[(Q1, N)][1]:.3f}" for N in NS]]
    tab = D.tex_table(rows, [2.8, 1.9, 1.9, 2.1], row_h=0.72, font_size=30)
    at(tab, 0, 0.75)
    for i in (2, 4):
        for j in range(1, 4):
            D.table_cell(tab, i, j, 4).set_color(D.GAUSS_COLOR)
    note = at(mixed("差趨近常數：", "$3000-2996=4$", "，", "$333.333-333.778=-\\frac49$", font_size=26), 0, -2.3)
    assert_in_frame(tab, note)
    assert_table(tab, 4)
    assert_apart(tab, note)
    s.play(FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:4])), run_time=0.8)
    s.say("精確變異數（W502 的逐步表）對漸近式 ", "$N\\frac{Q}{1-Q}$", wait=2.6)
    s.play(FadeIn(VGroup(*tab[0][4:12])), run_time=0.8)
    s.say("$Q=\\frac34$", "：26、296、2996 對 30、300、3000", wait=2.9)
    s.play(FadeIn(VGroup(*tab[0][12:20])), run_time=0.8)
    s.say("$Q=\\frac14$", "（老是回頭）：比普通醉漢還窄", wait=2.9)
    s.play(FadeIn(note), run_time=0.6)
    s.say("差趨近一個常數：步數越多，相對誤差越小", wait=2.9)
    s.clear_stage()


def part_length(s):
    nb = numbers()
    e1 = eq_label(tx(r"n=\frac{Q}{1-Q},\qquad l_p=n\,l", 34), "(5.28)--(5.31)")
    e2 = tx(r"\langle x^2\rangle\approx N\,\frac{Q}{1-Q}\,l^2=\frac{N}{n}\,(n\,l)^2", 32)
    at(e1, -2.9, 2.55)
    at(e2, 3.2, 2.55)
    fit(e2, 6.0)
    rows = [["Q"] + [frac_tex(Q) for Q in QS], ["n"] + [frac_tex(nb["n"][Q]) for Q in QS]]
    tab = D.tex_table(rows, [1.0, 1.2, 1.2, 1.2], row_h=0.8, font_size=32)
    at(tab, -4.0, 0.2)
    u = 0.5
    st = np.array([-1.3, 0.75, 0])
    small = rods(9, 1, u, st, O_COL, stroke=4)
    big = rods(3, 3, u, st + 1.2 * DOWN, P_COL, stroke=6)
    ls = tx(r"l", 30, O_COL).next_to(small, RIGHT, buff=0.3)
    lb = tx(r"l_p=3\,l", 30, P_COL).next_to(big, RIGHT, buff=0.3)
    var, asym = nb["big"]
    r1 = at(mixed("$Q=\\frac34$", "，1000 步：", "$\\langle x^2\\rangle\\approx1000\\cdot3=%d$" % round(asym),
                  "（精確 %.0f）" % var, font_size=28), 0, -1.45)
    r2 = at(mixed("＝普通醉漢走 3000 步＝步長 3 的醉漢走 %d 步" % nb["k"], "$\\ (3^2\\times%d=%d)$" % (nb["k"], 9 * nb["k"]),
                  font_size=28, color=YELLOW), 0, -2.3)
    for m in (r1, r2):
        fit(m, 13.0)
    rod_g = VGroup(small, big, ls, lb)
    assert_in_frame(e1, e2, tab, rod_g, r1, r2)
    for a, b in ((e1, e2), (e1, tab), (tab, rod_g), (e2, rod_g), (rod_g, r1), (r1, r2), (tab, r1)):
        assert_apart(a, b)
    assert_table(tab, 4)
    s.play(Write(e1), run_time=1.0)
    s.say("持續長度：同方向走一段，等於一根長 ", "$l_p$", " 的棍子", wait=2.9)
    s.play(FadeIn(tab), run_time=0.8)
    s.say("$Q=\\frac34$", " 是 3 格、普通醉漢 1 格、", "$Q=\\frac14$", " 只有 ", "$\\frac13$", " 格", wait=2.9)
    s.play(LaggedStart(*[GrowArrow(a) for a in small], lag_ratio=0.1), FadeIn(ls), run_time=1.2)
    s.play(LaggedStart(*[GrowArrow(a) for a in big], lag_ratio=0.2), FadeIn(lb), run_time=1.2)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("步長放大 ", "$n$", " 倍、步數變成 ", "$N/n$", "：變異數 ", "$Nn\\,l^2$", wait=2.9)
    s.play(FadeIn(r1), run_time=0.8)
    s.say("揭曉：1000 步的變異數約 3000（精確 2996）", wait=2.9)
    s.play(FadeIn(r2), run_time=0.8)
    s.say("寬度相當於普通醉漢走 3000 步", wait=2.9)
    s.clear_stage()


def part_extremes(s):
    nb = numbers()
    v1, v0a, v0b = nb["ext"]
    N = 11
    ux, uy = 0.42, 0.155
    o1, o2 = np.array([-6.0, -0.75, 0]), np.array([0.9, -0.75, 0])
    ax1 = D.walk_axes(N, N, ux=ux, uy=uy, origin=o1, x_every=5, y_every=5, font_size=18)
    ax2 = D.walk_axes(N, N, ux=ux, uy=uy, origin=o2, x_every=5, y_every=5, font_size=18)
    p1 = D.walk_path(list(range(N + 1)), ux=ux, uy=uy, origin=o1, color=P_COL)
    p2 = D.walk_path([k % 2 for k in range(N + 1)], ux=ux, uy=uy, origin=o2, color=O_COL)
    t1 = at(mixed("$Q=1$", "：直線前進", font_size=28, color=P_COL), -3.4, 2.75)
    t2 = at(mixed("$Q=0$", "：原地來回", font_size=28, color=O_COL), 3.5, 2.75)
    r1 = at(mixed("$\\langle x^2\\rangle=N^2$", "，", "$N=10$", "：%d" % round(v1), font_size=24), -3.4, 2.15)
    r2 = at(mixed("$N=10$", "：%d，" % round(v0a), "$N=11$", "：%d" % round(v0b), font_size=24), 3.5, 2.15)
    left, right = VGroup(ax1, p1), VGroup(ax2, p2)
    assert_in_frame(left, right, t1, t2, r1, r2)
    for a, b in ((left, right), (t1, r1), (t2, r2), (r1, left), (r2, right), (t1, t2)):
        assert_apart(a, b)
    s.play(FadeIn(ax1), FadeIn(ax2), FadeIn(t1), FadeIn(t2), run_time=0.8)
    s.say("兩個極端", wait=1.2)
    s.play(ShowCreation(p1), run_time=1.5)
    s.play(FadeIn(r1), run_time=0.6)
    s.say("$Q=1$", "：第一步決定方向，之後永遠直走", wait=2.9)
    s.say("變異數 ", "$N^2$", "：只有兩根尖峰，永遠不是高斯", wait=2.9)
    s.play(ShowCreation(p2), run_time=1.5)
    s.play(FadeIn(r2), run_time=0.6)
    s.say("$Q=0$", "：每一步都回頭，變異數 0、1、0、1…", wait=2.9)
    s.say("（", "$Q\\to1$", " 時 ", "$l_p\\to\\infty$", "；", "$Q\\to0$", " 時 ", "$l_p\\to0$", "）", wait=2.6)
    s.clear_stage()


def part_recursion(s):
    nb = numbers()
    ex, _ = nb["ex"]
    e1 = eq_label(tx(r"R(N+2)=2Q\cos ql\,R(N+1)+(1-2Q)\,R(N)", 30), "(5.34)(5.35)")
    e2 = eq_label(tx(r"\begin{pmatrix}R(N+2)\\R(N+1)\end{pmatrix}=\begin{pmatrix}2Q\cos ql&1-2Q\\1&0\end{pmatrix}"
                     r"\begin{pmatrix}R(N+1)\\R(N)\end{pmatrix}", 28), "(5.36)")
    e3 = eq_label(tx(r"\lambda^2-2Q\cos ql\,\lambda-(1-2Q)=0", 32), "(5.37)")
    box = D.erratum_box("p. 103，(5.37)", r"\lambda^2-2Q\cos ql\,\lambda+(1-2Q)=0",
                        r"\lambda^2-2Q\cos ql\,\lambda-(1-2Q)=0", width=9.0)
    ex2 = mixed("Ex. 5.2：", "$R=A_1\\lambda_1^N+A_2\\lambda_2^N$", "；", "$q=0.3,\\ N=5$", "：%s" % f4(ex.real),
                font_size=26)
    for m, y in ((e1, 2.65), (e2, 1.5), (e3, 0.35), (box, -0.95), (ex2, -2.35)):
        at(m, 0, y)
    assert_in_frame(e1, e2, e3, box, ex2)
    for a, b in ((e1, e2), (e2, e3), (e3, box), (box, ex2)):
        assert_apart(a, b)
    s.play(Write(e1), run_time=1.2)
    s.say("另一個角度：只看 ", "$R$", "，它滿足二階遞迴", wait=2.9)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("寫成 2×2 矩陣：特徵值和 W502 的 ", "$\\mathbf T$", " 一樣", wait=2.9)
    s.play(FadeIn(e3, shift=0.2 * DOWN), run_time=0.8)
    s.say("特徵方程 (5.37)：兩根就是 ", "$\\lambda_1$", "、", "$\\lambda_2$", wait=2.9)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (5.37) 的常數項號錯了", wait=2.6)
    s.say("把 ", "$\\lambda_1=%s$" % f4(nb["l1"].real), " 代入書上的式子，差 %.3f" % nb["book"], wait=2.9)
    s.play(FadeIn(ex2), run_time=0.8)
    s.say("Ex. 5.2：兩個 ", "$A\\lambda^N$", " 相加，回到 W502 的總數", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("慣性只是把步長放大成持續長度 ", "$l_p=\\frac{Q}{1-Q}\\,l$", wait=3.0)
    s.say("走得夠久，又是高斯", wait=2.9)


class W503PersistenceLength(ChuScene):
    def construct(self):
        self.show_header("W503　持續長度")
        check_math()
        part_question(self)
        part_gauss(self)
        part_table(self)
        part_length(self)
        part_extremes(self)
        part_recursion(self)
        part_close(self)
        self.wait(1.0)
