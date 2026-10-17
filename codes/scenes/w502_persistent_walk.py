"""
W502 有慣性的漫步：左右兩個狀態（§5.2，pp. 98–100，(5.11)–(5.22)）

問：醉漢有 3/4 的機率照上一步的方向走（Q = 3/4）。走 4 步之後，變異數比普通醉漢的 4 大還是小？
  1. 兩個狀態 (5.11)(5.12)：L、R＝最後一步向左／向右；逐步表 N = 0..4（V.persistent_dp）＝暴力列舉 2^N 條路線。
  2. 揭曉：變異數 N = 1..4，Q = 3/4、1/2、1/4（V.persistent_var）。
  3. 轉移矩陣 (5.15)(5.16)、特徵值 (5.17)(5.18)（q = 0.1 實數、q = 0.5 共軛複數）；傅立葉三種算法核對。
  4. 生成函數 (5.19)–(5.25)：伴隨矩陣 / D；大 N 由 λ1 主導。
畫面上的數字都由 utils.walk_var 算出，check_math() 核對。
"""
import cmath
import itertools
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

QF = Fraction(3, 4)
QV = 0.75
NMAX = 4
QS = (Fraction(3, 4), Fraction(1, 2), Fraction(1, 4))
QLAB = {Fraction(3, 4): "", Fraction(1, 2): "（普通）", Fraction(1, 4): "（老是回頭）"}
L_COL = D.TAIL_COLOR
R_COL = D.HEAD_COLOR
Q_EIG = (0.1, 0.5)
EIG_COL = {0.1: GREEN_C, 0.5: ORANGE}
Q_FT = 0.3
N_FT = (5, 20)
PATH_N, PATH_Y, PATH_SEED = 40, 16, 5


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q, small=True):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    cmd = r"\tfrac" if small else r"\frac"
    return ("-" if q < 0 else "") + cmd + "{%d}{%d}" % (abs(q.numerator), q.denominator)


def brute_dp(N, Q):
    """暴力列舉 2^N 條路線：第一步左右各半，之後同方向 Q、回頭 1−Q → {x: (L, R)}。"""
    if N == 0:
        return {0: (Fraction(1, 2), Fraction(1, 2))}
    out = {}
    for steps in itertools.product((-1, 1), repeat=N):
        w = Fraction(1, 2)
        for a, b in zip(steps, steps[1:]):
            w *= Q if a == b else 1 - Q
        x = sum(steps)
        lw, rw = out.get(x, (0, 0))
        out[x] = (lw + w, rw) if steps[-1] < 0 else (lw, rw + w)
    return out


def sample_path(n, Qv, seed):
    rng = np.random.default_rng(seed)
    d = 1 if rng.random() < 0.5 else -1
    pos = [0, d]
    for _ in range(n - 1):
        if rng.random() >= Qv:
            d = -d
        pos.append(pos[-1] + d)
    return pos


def c_series(N, q, Q):
    """1/D 的係數：c_k = 2Q cos q c_{k−1} − (2Q−1) c_{k−2}（c_0 = 1）。"""
    c = [1.0, 2 * Q * math.cos(q)]
    while len(c) <= N:
        c.append(2 * Q * math.cos(q) * c[-1] - (2 * Q - 1) * c[-2])
    return c


def eig_locus(Q, n=240):
    qs = np.linspace(0, PI, n)
    b1 = [V.eigs(q, Q)[0] for q in qs]
    b2 = [V.eigs(q, Q)[1] for q in qs]
    return b1, b2


@lru_cache(maxsize=None)
def numbers():
    nb = {"dp": {N: V.persistent_dp(N, QF) for N in range(NMAX + 1)},
          "brute": {N: brute_dp(N, QF) for N in range(NMAX + 1)},
          "var": {Q: [V.persistent_var(N, Q) for N in range(1, NMAX + 1)] for Q in QS},
          "eig": {q: V.eigs(q, QV) for q in Q_EIG}}
    ft = {}
    for N in N_FT:
        l1, l2 = V.eigs(Q_FT, QV)
        c = c_series(N, Q_FT, QV)
        ft[N] = dict(four=V.total_fourier(N, Q_FT, QV), mat=V.T_total(N, Q_FT, QV), eig=V.total_eig(N, Q_FT, QV),
                     ser=c[N] - (2 * QV - 1) * math.cos(Q_FT) * c[N - 1],
                     dom=l1 ** N * (l1 - (2 * QV - 1) * math.cos(Q_FT)) / (l1 - l2))
    nb["ft"] = ft
    nb["path"] = (sample_path(PATH_N, QV, PATH_SEED), sample_path(PATH_N, 0.5, PATH_SEED + 100))
    return nb


def check_math():
    nb = numbers()
    for N in range(NMAX + 1):                                                             # (5.11)(5.12) = 暴力列舉
        assert nb["dp"][N] == nb["brute"][N]
        assert sum(a + b for a, b in nb["dp"][N].values()) == 1
    dp = nb["dp"]
    assert dp[0] == {0: (Fraction(1, 2), Fraction(1, 2))}
    assert dp[1] == {-1: (Fraction(1, 2), 0), 1: (0, Fraction(1, 2))}
    assert dp[2] == {-2: (Fraction(3, 8), 0), 0: (Fraction(1, 8), Fraction(1, 8)), 2: (0, Fraction(3, 8))}
    assert dp[3] == {-3: (Fraction(9, 32), 0), -1: (Fraction(1, 8), Fraction(3, 32)),
                     1: (Fraction(3, 32), Fraction(1, 8)), 3: (0, Fraction(9, 32))}
    assert dp[4] == {-4: (Fraction(27, 128), 0), -2: (Fraction(15, 128), Fraction(9, 128)),
                     0: (Fraction(13, 128), Fraction(13, 128)), 2: (Fraction(9, 128), Fraction(15, 128)),
                     4: (0, Fraction(27, 128))}
    var = nb["var"]
    assert var[Fraction(3, 4)] == [1, 3, Fraction(11, 2), Fraction(33, 4)]
    assert var[Fraction(1, 2)] == [1, 2, 3, 4]
    assert var[Fraction(1, 4)] == [1, 1, Fraction(3, 2), Fraction(7, 4)]
    assert float(var[Fraction(3, 4)][3]) == 8.25
    for Q in QS:
        P4 = V.persistent_pmf(4, Q)
        assert sum(x * x * w for x, w in P4.items()) == var[Q][3]                         # 平均 0，變異數＝〈x²〉
    l1, l2 = nb["eig"][0.1]
    assert (f4(l1.real), f4(l2.real)) == ("0.9848", "0.5077") and l1.imag == 0 == l2.imag
    z1, z2 = nb["eig"][0.5]
    assert (f4(z1.real), f4(z1.imag), z2 == z1.conjugate()) == ("0.6582", "0.2584", True)
    assert QV * abs(math.sin(0.5)) > 1 - QV > QV * abs(math.sin(0.1))                     # 共軛複數的條件
    assert abs(abs(z1) ** 2 - (2 * QV - 1)) < 1e-12 and 2 * QF - 1 == Fraction(1, 2)
    for q in (0.1, 0.5, 1.3):                                                              # 特徵值＝數值特徵值
        num = np.linalg.eigvals(np.array(V.T_matrix(q, QV)))
        assert all(min(abs(num - lam)) < 1e-12 for lam in V.eigs(q, QV))
        assert abs(np.linalg.det(np.array(V.T_matrix(q, QV))) - (2 * QV - 1)) < 1e-12
    for q, z in ((0.3, 0.4), (1.1, 0.7 + 0.2j)):                                           # (5.19)–(5.22)
        T = np.array(V.T_matrix(q, QV))
        inv = np.linalg.inv(np.eye(2) - z * T)
        Dz = 1 - 2 * z * QV * math.cos(q) + z * z * (2 * QV - 1)
        adj = np.array([[1 - z * QV * cmath.exp(1j * q), z * (1 - QV) * cmath.exp(-1j * q)],
                        [z * (1 - QV) * cmath.exp(1j * q), 1 - z * QV * cmath.exp(-1j * q)]])
        assert np.max(np.abs(inv - adj / Dz)) < 1e-12
        a, b = V.eigs(q, QV)
        assert abs(Dz - (1 - z * a) * (1 - z * b)) < 1e-12
    ft = nb["ft"]
    for N in N_FT:
        r = ft[N]
        assert max(abs(r["four"] - r["mat"]), abs(r["four"] - r["eig"]), abs(r["four"] - r["ser"])) < 1e-12
        assert abs(r["four"].imag) < 1e-12
    assert f4(ft[5]["four"].real) == "0.5699" and f4(ft[20]["four"].real) == "0.0388"
    assert f4(ft[5]["dom"].real) == "0.6116"
    assert (f"{ft[20]['four'].real:.5f}", f"{ft[20]['dom'].real:.5f}") == ("0.03884", "0.03886")
    a, b = nb["path"]
    assert max(map(abs, a)) <= PATH_Y and max(map(abs, b)) <= PATH_Y and len(a) == len(b) == PATH_N + 1


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


def lr_entry(lw, rw, font_size=31):
    """(L, R) 一格：L 藍灰、R 金。"""
    parts = [tx("(", font_size), tx(frac_tex(lw), font_size, L_COL), tx(",", font_size),
             tx(frac_tex(rw), font_size, R_COL), tx(")", font_size)]
    g = VGroup(*parts).arrange(RIGHT, buff=0.05)
    parts[2].align_to(parts[1], DOWN).shift(0.04 * DOWN)
    return g


def dot_c(ax, z, color, radius=0.07):
    return Dot(ax.c2p(z.real, z.imag), radius=radius).set_fill(color, 1)


def poly_c(ax, zs, color, stroke=3):
    return VMobject().set_points_as_corners([ax.c2p(z.real, z.imag) for z in zs]).set_stroke(color, stroke)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    nl = D.number_line(4, unit=0.9)
    at(nl, 0, 1.6)
    walker = D.walker_dot().move_to(D.line_point(nl, 0) + 0.3 * UP)
    prev = Arrow(D.line_point(nl, -1) + 0.3 * UP, D.line_point(nl, 0) + 0.3 * UP, buff=0.15).set_color(GREY_B)
    prev_l = zh("上一步", font_size=22, color=GREY_B).next_to(prev, LEFT, buff=0.15)
    fwd = Arrow(D.line_point(nl, 0) + 0.95 * UP, D.line_point(nl, 1) + 0.95 * UP, buff=0).set_color(R_COL)
    back = Arrow(D.line_point(nl, 0) + 0.95 * UP, D.line_point(nl, -1) + 0.95 * UP, buff=0).set_color(L_COL)
    fwd_l = mixed("照原方向 ", "$Q=\\frac34$", font_size=24, color=R_COL).next_to(fwd, RIGHT, buff=0.2)
    back_l = mixed("回頭 ", "$1-Q=\\frac14$", font_size=24, color=L_COL).next_to(back, LEFT, buff=0.2)
    a, b = nb["path"]
    org = np.array([-5.7, -1.15, 0])
    ux, uy = 0.27, 0.095
    axes = D.walk_axes(PATH_N, PATH_Y, ux=ux, uy=uy, origin=org, x_every=10, y_every=8, font_size=18)
    pa = D.walk_path(a, ux=ux, uy=uy, origin=org, color=YELLOW, stroke_width=3, dot_radius=0.025)
    pb = D.walk_path(b, ux=ux, uy=uy, origin=org, color=BLUE_C, stroke_width=3, dot_radius=0.025)
    la = mixed("$Q=\\frac34$", font_size=22, color=YELLOW).next_to(pa[1][-1], RIGHT, buff=0.15)
    lb = mixed("$Q=\\frac12$", font_size=22, color=BLUE_C).next_to(pb[1][-1], RIGHT, buff=0.15)
    top = VGroup(nl, walker, prev, prev_l, fwd, back, fwd_l, back_l)
    bot = VGroup(axes, pa, pb, la, lb)
    assert_in_frame(top, bot)
    assert_apart(top, bot)
    assert_apart(la, lb)
    assert_apart(fwd_l, back_l)
    s.play(FadeIn(nl), FadeIn(walker), GrowArrow(prev), FadeIn(prev_l), run_time=0.8)
    s.say("問：這個醉漢記得上一步往哪走", wait=1.5)
    s.play(GrowArrow(fwd), FadeIn(fwd_l), GrowArrow(back), FadeIn(back_l), run_time=0.8)
    s.say("有 ", "$\\frac34$", " 的機率照上一步的方向走（", "$Q=\\frac34$", "）", wait=2.5)
    s.play(FadeIn(axes), run_time=0.6)
    s.play(ShowCreation(pa), ShowCreation(pb), run_time=2.5)
    s.play(FadeIn(la), FadeIn(lb), run_time=0.5)
    s.say("一次往同一邊走好幾步：路線拉得比較長", wait=2.0)
    s.say("走 4 步後，變異數比普通醉漢的 4 大還是小？（先猜）", wait=3.0)
    s.clear_stage()


def part_states(s):
    nb = numbers()
    dp = nb["dp"]
    e1 = tx(r"L(N+1,x)=Q\,L(N,x+1)+(1-Q)\,R(N,x+1)", 28)
    e2 = tx(r"R(N+1,x)=Q\,R(N,x-1)+(1-Q)\,L(N,x-1)", 28)
    e1[0].set_color(L_COL)
    e2[0].set_color(R_COL)
    eqs = eq_label(VGroup(e1, e2).arrange(DOWN, buff=0.18), "(5.11)(5.12)")
    at(eqs, -0.9, 2.55)
    x0, dx = 0.2, 0.78
    ys = [1.3, 0.5, -0.3, -1.1, -1.9]
    rows = VGroup()
    nlabs = VGroup()
    for N in range(NMAX + 1):
        r = VGroup(*[at(lr_entry(*dp[N][x]), x0 + x * dx, ys[N]) for x in sorted(dp[N])])
        rows.add(r)
        nlabs.add(at(tx(r"N=%d" % N, 26, GREY_A), -5.6, ys[N]))
    xl = VGroup(*[at(tx(r"x=%d" % x, 20, GREY_B), x0 + x * dx, -2.55) for x in range(-NMAX, NMAX + 1, 2)])
    leg = VGroup(mixed("$L$", "＝最後一步向左", font_size=22, color=L_COL),
                 mixed("$R$", "＝最後一步向右", font_size=22, color=R_COL)).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
    at(leg, 5.35, 0.5)
    hl = SurroundingRectangle(rows[4][2], buff=0.08).set_stroke(YELLOW, 2)
    assert_in_frame(eqs, rows, nlabs, xl, leg)
    assert_apart(eqs, VGroup(rows, leg))
    assert_apart(rows, leg)
    assert_apart(rows, xl)
    assert_apart(nlabs, rows)
    for N in range(NMAX + 1):
        for k in range(len(rows[N]) - 1):
            assert_apart(rows[N][k], rows[N][k + 1], 0.1)
        if N < NMAX:
            assert_apart(rows[N], rows[N + 1], 0.1)
    s.play(Write(eqs), run_time=1.5)
    s.say("兩個狀態：", "$L$", "＝最後一步向左、", "$R$", "＝最後一步向右", wait=2.5)
    s.say("走到 ", "$x$", " 而且向左：一定是從 ", "$x+1$", " 走過來的", wait=2.5)
    s.say("上一步也向左的照走（", "$Q$", "），向右的回頭（", "$1-Q$", "）", wait=2.5)
    s.play(FadeIn(leg), FadeIn(nlabs[0]), FadeIn(rows[0]), run_time=0.8)
    s.say("起點：還沒走，左右各半", wait=1.8)
    for N in range(1, NMAX + 1):
        s.play(FadeIn(nlabs[N]), LaggedStart(*[FadeIn(e, shift=0.15 * DOWN) for e in rows[N]], lag_ratio=0.15),
               run_time=1.0)
        if N == 1:
            s.say("第一步左右各半：", "$(\\frac12,0)$", " 與 ", "$(0,\\frac12)$", wait=1.8)
        elif N == 2:
            s.say("連走兩步同方向：", "$\\frac12\\cdot\\frac34=\\frac38$", wait=1.8)
        else:
            s.wait(1.0)
    s.play(FadeIn(xl), ShowCreation(hl), run_time=0.8)
    s.say("$N=4$", " 回到原點：", "$(\\frac{13}{128},\\frac{13}{128})$", wait=2.5)
    s.say("和暴力列舉 ", "$2^N$", " 條路線一模一樣（已核對）", wait=2.5)
    s.clear_stage()


def part_reveal(s):
    nb = numbers()
    var = nb["var"]
    hdr = [VectorizedPoint()] + [r"N=%d" % N for N in range(1, NMAX + 1)]
    rows = [hdr]
    for Q in QS:
        lab = mixed("$Q=%s$" % frac_tex(Q, small=False), QLAB[Q], font_size=24) if QLAB[Q] else \
            mixed("$Q=%s$" % frac_tex(Q, small=False), font_size=24)
        rows.append([lab] + [frac_tex(v, small=False) for v in var[Q]])
    tab = D.tex_table(rows, [3.4, 1.3, 1.3, 1.3, 1.3], row_h=0.8, font_size=32)
    at(tab, 0, 0.75)
    for j in range(1, 5):
        D.table_cell(tab, 1, j, 5).set_color(YELLOW)
        D.table_cell(tab, 3, j, 5).set_color(L_COL)
    big = D.table_cell(tab, 1, 4, 5)
    box = SurroundingRectangle(big, buff=0.1).set_stroke(YELLOW, 3)
    note = at(mixed("$\\frac{33}{4}=%.2f$" % float(var[Fraction(3, 4)][3]), "：比普通醉漢的 4 大一倍多", font_size=28),
              0, -2.2)
    assert_in_frame(tab, note, box)
    assert_table(tab, 5)
    assert_apart(tab, note)
    s.play(FadeIn(tab[1]), FadeIn(tab[0][:5]), run_time=0.8)
    s.say("揭曉：照逐步表算變異數", wait=1.5)
    s.play(FadeIn(VGroup(*tab[0][10:15])), run_time=0.8)
    s.say("普通醉漢：變異數＝步數", wait=2.0)
    s.play(FadeIn(VGroup(*tab[0][5:10])), run_time=0.8)
    s.say("$Q=\\frac34$", "：1、3、", "$\\frac{11}2$", "、", "$\\frac{33}4$", "，越長越快", wait=2.5)
    s.play(ShowCreation(box), FadeIn(note), run_time=0.8)
    s.say("同方向走一段，等於每一步變長了", wait=2.5)
    s.play(FadeIn(VGroup(*tab[0][15:20])), run_time=0.8)
    s.say("$Q=\\frac14$", "（老是回頭）：只有 ", "$\\frac74$", "，比 4 小", wait=3.0)
    s.clear_stage()


def part_matrix(s):
    nb = numbers()
    l1, l2 = nb["eig"][0.1]
    z1, _ = nb["eig"][0.5]
    m1 = tx(r"\begin{pmatrix}L_{N+1}\\R_{N+1}\end{pmatrix}=\mathbf T\begin{pmatrix}L_N\\R_N\end{pmatrix},\quad"
            r"\mathbf T=\begin{pmatrix}Qe^{-iq}&(1-Q)e^{-iq}\\(1-Q)e^{iq}&Qe^{iq}\end{pmatrix}", 26)
    e1 = fit(eq_label(m1, "(5.15)(5.16)", buff=0.25), 8.4)
    e2 = eq_label(tx(r"\lambda_{1,2}=Q\cos q\pm\sqrt{Q^2\cos^2q-(2Q-1)}", 30), "(5.17)(5.18)")
    cond = mixed("共軛複數：", "$Q|\\sin q|>1-Q$", font_size=24)
    at(e1, -2.35, 2.5)
    at(e2, -2.35, 1.3)
    at(cond, -2.35, 0.45)
    ax = D.plot_axes((-1.1, 1.1, 0.5), (-0.9, 0.9, 0.5), width=4.4, height=3.6,
                     x_label=r"\mathrm{Re}\,\lambda")
    ax.move_to(np.array([4.35, 0.75, 0]))
    ax.add(*[tx(t, 20, GREY_B).next_to(ax.c2p(v, 0), DOWN + (LEFT if v < 0 else RIGHT), buff=0.06)
             for v, t in ((-1, "-1"), (1, "1"))])
    b1, b2 = eig_locus(QV)
    loc = VGroup(poly_c(ax, b1, GREY_B), poly_c(ax, b2, GREY_B))
    d1 = VGroup(dot_c(ax, l1, EIG_COL[0.1]), dot_c(ax, l2, EIG_COL[0.1]))
    d2 = VGroup(dot_c(ax, z1, EIG_COL[0.5]), dot_c(ax, z1.conjugate(), EIG_COL[0.5]))
    g1 = at(mixed("$q=0.1$", "：", "$\\lambda_1=%s,\\ \\lambda_2=%s$" % (f4(l1.real), f4(l2.real)), font_size=24,
                  color=EIG_COL[0.1]), 4.2, -1.45)
    g2 = at(mixed("$q=0.5$", "：", "$\\lambda=%s\\pm%si$" % (f4(z1.real), f4(z1.imag)), font_size=24,
                  color=EIG_COL[0.5]), 4.2, -2.1)
    for g in (g1, g2):
        fit(g, 4.9)
    ft = nb["ft"]
    hdr = lambda t: zh(t, font_size=22, color=YELLOW)
    rows = [[r"q=0.3", hdr("傅立葉"), hdr("乘矩陣"), hdr("特徵值")]]
    for N in N_FT:
        r = ft[N]
        rows.append([r"N=%d" % N, f4(r["four"].real), f4(r["mat"].real), f4(r["eig"].real)])
    tab = D.tex_table(rows, [1.4, 1.5, 1.5, 1.5], row_h=0.5, font_size=26)
    at(tab, -2.35, -1.55)
    right = VGroup(ax, loc, g1, g2)
    assert_in_frame(e1, e2, cond, tab, right)
    for a, b in ((e1, e2), (e2, cond), (cond, tab), (e1, right), (e2, right), (cond, right), (tab, right),
                 (ax, g1), (g1, g2)):
        assert_apart(a, b)
    assert_table(tab, 4)
    s.play(Write(e1), run_time=1.5)
    s.say("轉成傅立葉：每走一步，", "$(L,R)$", " 就乘上矩陣 ", "$\\mathbf T$", wait=2.5)
    s.say("$e^{-iq}$", "、", "$e^{iq}$", "：向左、向右各移一格", wait=2.0)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("走 ", "$N$", " 步就是 ", "$\\mathbf T^N$", "：要看它的特徵值", wait=2.0)
    s.play(FadeIn(ax), ShowCreation(loc), run_time=1.5)
    s.say("$q$", " 從 0 到 ", "$\\pi$", "：兩個特徵值在複平面上走的路", wait=2.0)
    s.play(FadeIn(d1, scale=0.5), FadeIn(g1), run_time=0.8)
    s.say("$q=0.1$", "：兩個實數，", "$\\lambda_1$", " 接近 1", wait=2.0)
    s.play(FadeIn(d2, scale=0.5), FadeIn(g2), FadeIn(cond), run_time=0.8)
    s.say("$q=0.5$", "：", "$Q|\\sin q|>1-Q$", "，變成共軛複數", wait=2.5)
    s.say("共軛時 ", "$|\\lambda|^2=\\det\\mathbf T=2Q-1=\\frac12$", "：在圓上", wait=2.5)
    s.play(FadeIn(tab[1]), FadeIn(tab[0]), run_time=0.8)
    s.say("核對總數：三種算法，", "$N=5$", " 都是 0.5699", wait=2.5)
    s.clear_stage()


def part_gf(s):
    nb = numbers()
    ft = nb["ft"]
    g1 = eq_label(tx(r"(1-\mathbf Tz)^{-1}=\frac{1}{D}\begin{pmatrix}1-zQe^{iq}&z(1-Q)e^{-iq}\\"
                     r"z(1-Q)e^{iq}&1-zQe^{-iq}\end{pmatrix}", 30), "(5.19)--(5.21)")
    g2 = eq_label(tx(r"D=1-2zQ\cos q+z^2(2Q-1)=(1-z\lambda_1)(1-z\lambda_2)", 30), "(5.22)")
    g3 = eq_label(tx(r"[z^N]\frac{1}{D}=c_N=\frac{\lambda_1^{N+1}-\lambda_2^{N+1}}{\lambda_1-\lambda_2}", 30),
                  "(5.24)(5.25)")
    g4 = mixed("從 ", "$L=R=\\frac12$", " 出發：總數 ", "$=c_N-(2Q-1)\\cos q\\,c_{N-1}$", font_size=26)
    dom = mixed("只留 ", "$\\lambda_1$", "：", "$N=5$", " 得 %s（精確 %s）；" % (f4(ft[5]["dom"].real), f4(ft[5]["four"].real)),
                "$N=20$", " 得 %.5f（精確 %.5f）" % (ft[20]["dom"].real, ft[20]["four"].real), font_size=24)
    for m, y in ((g1, 2.35), (g2, 1.0), (g3, -0.1), (g4, -1.3), (dom, -2.3)):
        fit(m, 12.6)
        at(m, 0, y)
    assert_in_frame(g1, g2, g3, g4, dom)
    for a, b in ((g1, g2), (g2, g3), (g3, g4), (g4, dom)):
        assert_apart(a, b)
    s.play(Write(g1), run_time=1.5)
    s.say("生成函數：把每一步的 ", "$(\\mathbf Tz)^N$", " 全部加起來", wait=2.5)
    s.say("2×2 的反矩陣＝伴隨矩陣除以行列式 ", "$D$", wait=2.5)
    s.play(FadeIn(g2, shift=0.2 * DOWN), run_time=0.8)
    s.say("$D$", " 拆成兩個因式：根是 ", "$1/\\lambda_1$", " 與 ", "$1/\\lambda_2$", wait=2.5)
    s.play(FadeIn(g3, shift=0.2 * DOWN), run_time=0.8)
    s.say("部分分式展開，讀出 ", "$z^N$", " 的係數", wait=2.5)
    s.play(FadeIn(g4, shift=0.2 * DOWN), run_time=0.8)
    s.say("從左右各半出發，", "$L+R$", " 的總數就是這一行", wait=2.5)
    s.play(FadeIn(dom), run_time=0.8)
    s.say("大 ", "$N$", " 由最大的 ", "$\\lambda_1$", " 主導：", "$N=20$", " 幾乎一樣", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("記得上一步的醉漢，用兩個狀態就又是「乘矩陣」", wait=3.0)
    s.say("走久了，由最大的特徵值說了算", wait=3.0)


class W502PersistentWalk(ChuScene):
    def construct(self):
        self.show_header("W502　有慣性的漫步：左右兩個狀態")
        check_math()
        part_question(self)
        part_states(self)
        part_reveal(self)
        part_matrix(self)
        part_gf(self)
        part_close(self)
        self.wait(1.0)
