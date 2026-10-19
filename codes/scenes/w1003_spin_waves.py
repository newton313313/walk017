"""
W1003 自旋波：在平均場上加漣漪（§10.3，pp. 239–243，(10.16)–(10.36)，Ex. 10.1–10.4）

問：平均場把所有自旋看成一樣：在它上面加一點漣漪（自旋波），會多出什麼？
  1. (10.18) S = M + σ；(10.19)–(10.21) 零階、一階（= 0，狀態方程 (10.17)）、二階；Ex. 10.2 三階、四階（誤差 < 1e−11）。
  2. (10.17)、Ex. 10.1：牆邊 M(x) = √(|r|/u) tanh(x√(|r|/2))（r = −1、u = 1）：M(1), M(2), M(4) = 0.6089、0.8884、0.9930；
     癒合長度 √2 = 1.4142；有限差分＋牛頓法與 tanh 相差 2.8e−06；週期邊界沒有牆，M 是常數。
  3. Fourier 模態 (10.22)–(10.24)（cos(2πnx/8)，n = 0, 1, 2）；(10.27) X_q；例子 M = (2/3, 2/3, 1/3)：X = (1/18)(…)，
     λ1 = 1/2（橫向，n − 1 = 2 重）、λ2 = 3/2（縱向）；q = 0 時 λ1 = 0。
  4. (10.30)(10.31)：X⁻¹ = I/λ1 − u|M⟩⟨M|/(λ1λ2) = 2I − (4/3)|M⟩⟨M| = (1/27)(…)，X X⁻¹ = X⁻¹ X = I；隨機 M、n = 1..6 全對；
     (10.32) 高斯積分；筆誤框 1（p. 242，(10.32) 指數的負號）：n = 2 數值 3.9891，負號得 3.2988。
  5. n → 0：表（橫向重數 2, 1, 0, −1；det X = 3/8, 3/4, 3/2, 3；(det X)^{−1/2} = 1.6330、1.1547、0.8165、0.5774 = √(λ1/λ2)）；
     筆誤框 2（p. 242，(10.33) 左邊顛倒：書上 √3 = 1.7321）。
畫面上的數字都由 utils.walk_cond 算出（分數用 Fraction），check_math() 核對。
模組層級的畫圖函式：spin_row、mode_curves、mode_diagram。
"""
import math
import sys
import warnings
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_cond as WK
from utils import walk_draw as D

F = Fraction
WALL_X = (1, 2, 4)
WALL_M = ("0.6089", "0.8884", "0.9930")
X18 = ((17, 8, 4), (8, 17, 4), (4, 4, 11))
XI27 = ((38, -16, -8), (-16, 38, -8), (-8, -8, 50))
MODES = ((3, 2, 1, F(3, 8), "1.6330"), (2, 1, 1, F(3, 4), "1.1547"), (1, 0, 1, F(3, 2), "0.8165"),
         (0, -1, 1, F(3), "0.5774"))
SPIN_COLOR = YELLOW
TRANS_COLOR = TEAL_B
LONG_COLOR = "#F2A33A"
GOOD = "#7CE38B"


def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q, small=False):
    q = F(q)
    if q.denominator == 1:
        return str(q.numerator)
    cmd = r"\tfrac" if small else r"\frac"
    sign = "-" if q < 0 else ""
    return rf"{sign}{cmd}{{{abs(q.numerator)}}}{{{q.denominator}}}"


def pmatrix(rows):
    return r"\begin{pmatrix}" + r"\\".join("&".join(str(v) for v in row) for row in rows) + r"\end{pmatrix}"


def sci_tex(x):
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def matvec(A, v):
    return [sum(a * x for a, x in zip(row, v)) for row in A]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return dict(coef=WK.ex102_coeffs(), wall=tuple(WK.wall_errors()), swc=WK.sw_check(),
                    gn2=tuple(WK.gauss_n2_check()), modes=tuple(WK.mode_table()))


def check_math():
    nb = numbers()
    assert nb["coef"] < 1e-11
    r, u, q2 = WK.SW_R, WK.SW_U, WK.SW_Q2
    assert (r, u, q2) == (-1, 1, 1) and tuple(WK.WALL_X) == WALL_X
    assert [f4(WK.wall_M(x)) for x in WALL_X] == list(WALL_M) and f4(WK.wall_M(50.0)) == "1.0000"
    assert all(WK.wall_M(x) == WK.wall_M(x, -1.0, 1.0) for x in WALL_X)
    assert f4(math.sqrt(2 / abs(r))) == "1.4142" and f4(math.sqrt(abs(r) / u)) == "1.0000"
    assert sci_tex(nb["wall"][0]) == r"2.8\times10^{-6}"
    M = tuple(WK.SW_M)
    assert M == (F(2, 3), F(2, 3), F(1, 3)) and sum(m * m for m in M) == 1 == F(abs(r), u)
    X = WK.sw_X()
    assert [list(row) for row in X] == [[F(v, 18) for v in row] for row in X18]
    l1, l2 = WK.sw_lambdas()
    assert (l1, l2) == (F(1, 2), F(3, 2))
    assert l1 == F(r, 2) + F(q2, 2) + F(u, 2) and l2 == l1 + u
    assert matvec(X, M) == [l2 * m for m in M]
    for v in ((1, -1, 0), (1, 1, -4)):
        assert sum(a * b for a, b in zip(v, M)) == 0 and matvec(X, v) == [l1 * a for a in v]
    assert WK.sw_lambdas(q2=F(0))[0] == 0
    Xi = WK.sw_Xinv()
    assert [list(row) for row in Xi] == [[F(v, 27) for v in row] for row in XI27]
    assert 1 / l1 == 2 and F(u) / (l1 * l2) == F(4, 3)
    I3 = [[F(int(i == j)) for j in range(3)] for i in range(3)]
    assert [[2 * I3[i][j] - F(4, 3) * M[i] * M[j] for j in range(3)] for i in range(3)] == [list(row) for row in Xi]
    assert matmul(X, Xi) == I3 == matmul(Xi, X)
    errs, ok = nb["swc"]
    assert ok and max(errs.values()) < 1e-11
    g_num, g_formula, g_book = nb["gn2"]
    assert f4(g_num) == f4(g_formula) == "3.9891" and f4(g_book) == "3.2988"
    modes = nb["modes"]
    assert [tuple(m[:4]) for m in modes] == [m[:4] for m in MODES]
    assert [f4(m[4]) for m in modes] == [m[4] for m in MODES]
    for n, t, lo, det, inv in modes:
        assert t == n - 1 and lo == 1 and det == l1 ** (n - 1) * l2 and abs(inv - float(det) ** -0.5) < 1e-12
    assert f4(float(l1 / l2) ** 0.5) == "0.5774" == f4(1 / math.sqrt(3))
    assert f4(float(l2 / l1) ** 0.5) == "1.7321"


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框 [x0,x1]×[y0,y1]（只縮不放），再置中在方框裡。"""
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def place_left(m, x0, y, x1=6.7):
    """左緣對齊 x0、中心高度 y；太寬就縮到 x1 為止。"""
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    return m.move_to(np.array([x0 + m.get_width() / 2, y, 0]))


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def assert_pairs(*pairs):
    for a, b in pairs:
        assert_apart(a, b)


def unit(v):
    v = np.array(v, dtype=float)
    return v / np.linalg.norm(v)


def arrow(a, b, color, thickness=4, ratio=0.16):
    return Arrow(np.array(a, dtype=float), np.array(b, dtype=float), buff=0, thickness=thickness,
                 max_width_to_length_ratio=ratio).set_fill(color, 1)


# ---------------------------------------------------------------- 畫圖元件
def spin_row(x0, x1, y, count=12, length=0.72):
    """一排朝上的自旋（平均場）→ (VGroup(箭頭們), 底端座標們)。"""
    bases = [np.array([x, y, 0.0]) for x in np.linspace(x0, x1, count)]
    return VGroup(*[arrow(b, b + length * UP, SPIN_COLOR) for b in bases]), bases


def mode_curves(x0, x1, ys, h=0.38):
    """一維 L = 8 的 cos(2πnx/8)，n = 0, 1, 2 → VGroup(VGroup(基線, 曲線, 標籤) × 3)。"""
    g = VGroup()
    for n, y in zip((0, 1, 2), ys):
        base = Line([x0, y, 0], [x1, y, 0]).set_stroke(GREY_D, 1.5)
        pts = [np.array([x0 + (x1 - x0) * t / 8, y + h * math.cos(2 * math.pi * n * t / 8), 0])
               for t in np.linspace(0, 8, 60)]
        c = VMobject().set_points_as_corners(pts).set_stroke(TRANS_COLOR, 3)
        lab = tx(f"n={n}", 26, GREY_A).next_to(base, LEFT, buff=0.25)
        g.add(VGroup(base, c, lab))
    return g


def mode_diagram(c, mlen=1.45):
    """平面上的示意：M（黃）、縱向（沿 M，橘）、垂直於 M 的平面（斜看成橢圓）裡兩個橫向（青）。
    → VGroup(平面, M, 縱向, 橫向1, 橫向2, M 標籤, λ2 標籤, λ1 標籤)"""
    c = np.array(c, dtype=float)
    plane = Ellipse(width=2.3, height=0.75).move_to(c).set_stroke(TRANS_COLOR, 2).set_fill(TRANS_COLOR, 0.15)
    m = arrow(c, c + mlen * UP, SPIN_COLOR, thickness=5, ratio=0.12)
    lon = arrow(c + mlen * UP, c + (mlen + 0.6) * UP, LONG_COLOR)
    t1 = arrow(c, c + 1.0 * RIGHT, TRANS_COLOR)
    t2 = arrow(c, c + 0.8 * unit([-0.85, -0.5, 0]), TRANS_COLOR)
    mlab = tx(r"\vec M", 28, SPIN_COLOR).next_to(c + 0.7 * mlen * UP, LEFT, buff=0.12)
    l2lab = tx(r"\lambda_2", 28, LONG_COLOR).next_to(c + (mlen + 0.3) * UP, RIGHT, buff=0.15)
    l1lab = tx(r"\lambda_1", 28, TRANS_COLOR).next_to(c + 1.15 * RIGHT, RIGHT, buff=0.1)
    return VGroup(plane, m, lon, t1, t2, mlab, l2lab, l1lab)


# ---------------------------------------------------------------- 各段
def part_expand(s):
    row, bases = spin_row(-6.2, -1.8, 1.95)
    angles = [0.45 * math.sin(2 * math.pi * k / 8) for k in range(len(bases))]
    s18 = place_in(tx(r"\vec S=\vec M+\vec\sigma", 44), 0.0, 6.6, 2.0, 2.95)
    lab18 = place_in(mixed("平均場 ", r"$\vec M$", "＋漣漪（自旋波）", r"$\vec\sigma$", font_size=24, color=GREY_B),
                     0.0, 6.6, 1.4, 1.85)
    texs = (r"\mathcal H_{MF}",
            r"\int(r\vec M-\nabla^2\vec M+u\vec M|\vec M|^2-\vec h)\cdot\vec\sigma=0",
            r"\int\Big\{\frac12|\nabla\vec\sigma|^2+\frac r2|\vec\sigma|^2+\frac u2|\vec M|^2|\vec\sigma|^2"
            r"+u(\vec M\cdot\vec\sigma)^2\Big\}")
    ys = (1.0, 0.2, -0.7)
    tags = VGroup(*[zh(t, font_size=28).move_to([-6.1, y, 0]) for t, y in zip(("零階", "一階", "二階"), ys)])
    eqs = VGroup(*[place_left(tx(t, 34), -5.4, y, 6.6) for t, y in zip(texs, ys)])
    why = place_left(zh("（狀態方程 (10.17)）", font_size=24, color=GREY_B), eqs[1].get_right()[0] + 0.3, ys[1], 6.7)
    ex = place_left(mixed("Ex. 10.2：三階 ", r"$u\int(\vec M\cdot\vec\sigma)|\vec\sigma|^2$", "，四階 ",
                          r"$\frac u4\int|\vec\sigma|^4$", font_size=28), -6.4, -1.75, 6.6)
    exn = place_left(mixed("數值內插＝公式，誤差 ", r"$<10^{-11}$", font_size=24, color=GREY_B), -6.4, -2.5, 6.6)
    assert_in_frame(row, s18, lab18, tags, eqs, why, ex, exn)
    assert_pairs((row, s18), (row, lab18), (row, eqs[0]), (s18, lab18), (lab18, eqs[0]), (eqs[0], eqs[1]),
                 (eqs[1], eqs[2]), (eqs[1], why), (eqs[2], ex), (ex, exn), *zip(tags, eqs))

    s.play(FadeIn(row, lag_ratio=0.05), run_time=1.0)
    s.say("平均場把所有自旋看成一樣", wait=1.5)
    s.play(*[Rotate(a, ang, about_point=b) for a, ang, b in zip(row, angles, bases)], run_time=1.5)
    s.say("在它上面加一點漣漪（自旋波），會多出什麼？", wait=2.0)
    s.play(Write(s18), FadeIn(lab18), run_time=1.2)
    s.say("(10.18)：", r"$\vec S=\vec M+\vec\sigma$", "，把能量對 ", r"$\vec\sigma$", " 展開", wait=2.0)
    s.play(FadeIn(tags[0]), FadeIn(eqs[0]), FadeIn(tags[1]), FadeIn(eqs[1]), FadeIn(why), run_time=1.2)
    s.say("零階是平均場；一階是 0：", r"$\vec M$", " 滿足狀態方程 (10.17)", wait=2.5)
    s.play(FadeIn(tags[2]), FadeIn(eqs[2]), run_time=0.8)
    s.say("二階：", r"$\vec\sigma$", " 的二次式，下面把它對角化", wait=2.0)
    s.play(FadeIn(ex), FadeIn(exn), run_time=0.8)
    s.say("Ex. 10.2：三階、四階項；數值內插和公式一致", wait=2.5)
    s.clear_stage()


def part_wall(s):
    nb = numbers()
    r, u = int(WK.SW_R), int(WK.SW_U)
    heal = f4(math.sqrt(2 / abs(r)))
    ax = D.plot_axes((-1, 6, 1), (0, 1.2, 0.5), width=5.6, height=3.4, x_ticks=(0, 1, 2, 4, 6),
                     y_ticks=(0.5, 1), y_tick_labels=["0.5", "1"], x_label="x")
    ax.move_to(np.array([-3.4, 0.1, 0]))
    p0, p1 = ax.c2p(-0.7, 0), ax.c2p(0, 1.2)
    ww, wh = p1[0] - p0[0], p1[1] - p0[1]
    wall = VGroup(Rectangle(width=ww, height=wh).move_to((p0 + p1) / 2).set_stroke(GREY_B, 2).set_fill(GREY_D, 0.85))
    for y in np.arange(p0[1] + 0.1, p1[1] - ww, 0.3):
        wall.add(Line([p0[0], y, 0], [p1[0], y + ww, 0]).set_stroke(GREY_B, 1.5))
    wall_lab = zh("牆", font_size=24).next_to(wall[0], UP, buff=0.12, aligned_edge=LEFT)
    mlab = tx("M", 26, GREY_B).next_to(ax.c2p(0, 1.2), RIGHT, buff=0.15)
    curve = VMobject().set_points_as_corners([ax.c2p(x, WK.wall_M(x)) for x in np.linspace(0, 6, 60)])
    curve.set_stroke(SPIN_COLOR, 4)
    asym = DashedLine(ax.c2p(0, 1), ax.c2p(6, 1), dash_length=0.1).set_stroke(GREY_B, 2)
    asym_lab = tx(rf"M=\sqrt{{|r|/u}}={math.sqrt(abs(r) / u):g}", 22, GREY_A).next_to(ax.c2p(4.6, 1), UP, buff=0.08)
    dots = VGroup(*[Dot(ax.c2p(x, WK.wall_M(x)), radius=0.06).set_fill(LONG_COLOR, 1) for x in WALL_X])
    dlabs = VGroup(*[tx(f4(WK.wall_M(x)), 22, LONG_COLOR).next_to(d, DR, buff=0.06) for x, d in zip(WALL_X, dots)])
    e17 = place_left(tx(r"-\nabla^2\vec M+r\vec M+u\vec M|\vec M|^2-\vec h=0", 32), 0.0, 2.55, 6.7)
    e17n = place_left(mixed("Ex. 10.1：一維、", r"$\vec h=0$", "、牆上 ", "$M(0)=0$", font_size=26), 0.0, 1.85, 6.7)
    sol = place_left(tx(r"M(x)=\sqrt{|r|/u}\,\tanh\big(x\sqrt{|r|/2}\big)", 34), 0.0, 1.1, 6.7)
    par = place_left(mixed(f"$r={r}$", "、", f"$u={u}$", "：癒合長度 ", rf"$\sqrt{{2/|r|}}={heal}$", font_size=26),
                     0.0, 0.3, 6.7)
    vals = place_left(tx(r",\ ".join(f"M({x})={f4(WK.wall_M(x))}" for x in WALL_X), 28), 0.0, -0.4, 6.7)
    errtex = sci_tex(nb["wall"][0])
    bvp1 = place_left(zh("有限差分＋牛頓法解邊值問題，", font_size=26), 0.0, -1.2, 6.7)
    bvp2 = place_left(mixed("和 tanh 只差 ", f"${errtex}$", font_size=26), 0.0, -1.8, 6.7)
    plot = VGroup(ax, wall, wall_lab, mlab, curve, asym, asym_lab, dots, dlabs)
    right = (e17, e17n, sol, par, vals, bvp1, bvp2)
    assert_in_frame(plot, *right)
    for m in right:
        assert_apart(plot, m)
    assert_pairs(*zip(right, right[1:]), (asym_lab, dlabs), (mlab, asym_lab))

    s.play(FadeIn(ax), FadeIn(wall), FadeIn(wall_lab), FadeIn(mlab), run_time=1.0)
    s.play(FadeIn(e17), run_time=0.8)
    s.say("先看有牆的情形：(10.17) 決定 ", r"$\vec M(x)$", wait=2.0)
    s.play(FadeIn(e17n), Write(sol), run_time=1.2)
    s.play(ShowCreation(curve), run_time=1.5)
    s.play(ShowCreation(asym), FadeIn(asym_lab), run_time=0.6)
    s.say("Ex. 10.1：牆上 ", "$M=0$", "，離牆越遠越接近 ", r"$\sqrt{|r|/u}$", wait=2.5)
    s.play(FadeIn(par), FadeIn(dots), FadeIn(dlabs), FadeIn(vals), run_time=0.8)
    s.say(f"$r={r}$", "、", f"$u={u}$", "：", "$M(1),M(2),M(4)$", " 很快接近 1", wait=2.5)
    s.say("癒合長度 ", rf"$\sqrt{{2/|r|}}={heal}$", "：幾個癒合長度外就回到平均場", wait=2.5)
    s.play(FadeIn(bvp1), FadeIn(bvp2), run_time=0.8)
    s.say("有限差分＋牛頓法解邊值問題，和 tanh 只差 ", f"${errtex}$", wait=2.5)
    s.say("週期邊界沒有牆，", r"$\vec M$", " 就是常數（本章接下來都用週期邊界）", wait=2.5)
    s.clear_stage()


def part_modes_matrix(s):
    r, u, q2 = int(WK.SW_R), int(WK.SW_U), int(WK.SW_Q2)
    M = tuple(WK.SW_M)
    X = WK.sw_X()
    l1, l2 = WK.sw_lambdas()
    n = len(M)
    curves = mode_curves(-5.3, -2.3, (2.15, 1.0, -0.15))
    mtitle = place_in(mixed("一維 ", "$L=8$", "：", r"$\cos(2\pi nx/8)$", font_size=24), -6.4, -1.8, 2.7, 3.15)
    psi = place_left(tx(r"\psi_{\vec q}(\vec x)=\prod_i\frac{e^{iq_ix_i}}{\sqrt{L_i}}", 38), -1.2, 1.95, 6.6)
    qi = place_left(tx(r"q_i=\frac{2\pi n_i}{L}", 38), -1.2, 0.55, 6.6)
    xq = r"X_{\vec q}=\Big(\frac r2+\frac12q^2+\frac12u|\vec M|^2\Big)I+u|\vec M\rangle\langle\vec M|"
    X27 = place_in(tx(xq, 40), -6.4, 6.6, -2.6, -1.2)
    assert_in_frame(curves, mtitle, psi, qi, X27)
    assert_pairs((mtitle, curves), (curves, psi), (curves, qi), (psi, qi), (curves, X27), (qi, X27))
    xtop = place_in(tx(xq, 32), -6.4, 6.6, 2.25, 2.95)
    mtex = r"\vec M=(" + ",".join(frac_tex(m, small=True) for m in M) + ")"
    exl = place_in(mixed(f"$r={r}$", "、", f"$u={u}$", "、", f"$q^2={q2}$", "、", f"${mtex}$", "：",
                         r"$|\vec M|^2=1=|r|/u$", "（平均場的鞍點）", font_size=26), -6.4, 6.6, 1.4, 2.0)
    x18 = [[int(v * 18) for v in row] for row in X]
    mat = place_in(tx(r"X=\frac1{18}" + pmatrix(x18), 34), -6.5, -2.9, -1.5, 1.1)
    l1t = place_in(tx(r"\lambda_1=\frac r2+\frac12q^2+\frac12u|\vec M|^2=" + frac_tex(l1), 30), -2.5, 3.6, 0.45, 1.15)
    l1n = place_in(mixed("橫向：", f"$n-1={n - 1}$", " 重（垂直於 ", r"$\vec M$", "）", font_size=24,
                         color=TRANS_COLOR), -2.5, 3.6, -0.15, 0.3)
    l2t = place_in(tx(r"\lambda_2=\lambda_1+u|\vec M|^2=" + frac_tex(l2), 30), -2.5, 3.6, -1.0, -0.35)
    l2n = place_in(mixed("縱向：1 重（沿著 ", r"$\vec M$", "）", font_size=24, color=LONG_COLOR), -2.5, 3.6, -1.6, -1.15)
    diag = mode_diagram(np.array([5.0, -1.35, 0]))
    assert_in_frame(xtop, exl, mat, l1t, l1n, l2t, l2n, diag)
    assert_pairs((xtop, exl), (exl, mat), (exl, l1t), (exl, diag), (mat, l1t), (mat, l1n), (mat, l2t), (mat, l2n),
                 (l1t, l1n), (l1n, l2t), (l2t, l2n), (l1t, diag), (l1n, diag), (l2t, diag), (l2n, diag))
    eig = ",".join(frac_tex(v) for v in (l1, l1, l2))

    s.play(FadeIn(mtitle), *[FadeIn(g[0]) for g in curves], *[FadeIn(g[2]) for g in curves], run_time=0.8)
    s.play(*[ShowCreation(g[1]) for g in curves], run_time=1.5)
    s.play(Write(psi), run_time=1.0)
    s.say("週期邊界上的漣漪：Fourier 模態 (10.22)", wait=2.0)
    s.play(FadeIn(qi), run_time=0.6)
    s.say(r"$q_i=2\pi n_i/L$", "：每個方向放整數個波長", wait=2.0)
    s.play(Write(X27), run_time=1.5)
    s.say("(10.27)：二階項對每個 ", r"$\vec q$", " 變成一個矩陣 ", r"$X_{\vec q}$", wait=2.5)
    s.play(*[FadeOut(m) for m in (mtitle, curves, psi, qi)], ReplacementTransform(X27, xtop), run_time=1.0)
    s.play(FadeIn(exl), run_time=0.6)
    s.say("例子：", f"$r={r}$", "、", f"$u={u}$", "、", f"$q^2={q2}$", "、", r"$|\vec M|^2=1$", "（平均場的鞍點）",
          wait=2.5)
    s.play(FadeIn(mat), run_time=0.8)
    s.say(r"$X=\frac1{18}(\cdots)$", "：特徵值 ", f"${eig}$", wait=2.5)
    s.play(FadeIn(l1t), FadeIn(l1n), FadeIn(diag[0]), GrowArrow(diag[1]), FadeIn(diag[5]), run_time=1.0)
    s.play(GrowArrow(diag[3]), GrowArrow(diag[4]), FadeIn(diag[7]), run_time=0.8)
    s.say("橫向 ", rf"$\lambda_1={frac_tex(l1)}$", "：垂直於 ", r"$\vec M$", " 的平面，", f"$n-1={n - 1}$", " 重",
          wait=2.5)
    s.play(FadeIn(l2t), FadeIn(l2n), GrowArrow(diag[2]), FadeIn(diag[6]), run_time=0.8)
    s.say("縱向 ", rf"$\lambda_2=\lambda_1+u|\vec M|^2={frac_tex(l2)}$", "：沿著 ", r"$\vec M$", "，1 重", wait=2.5)
    s.say("$q=0$", " 時 ", rf"$\lambda_1={WK.sw_lambdas(q2=F(0))[0]}$", "：整體轉動 ", r"$\vec M$", " 不花能量",
          wait=2.5)
    s.clear_stage()


def part_inverse(s):
    nb = numbers()
    u = WK.SW_U
    l1, l2 = WK.sw_lambdas()
    Xi = WK.sw_Xinv()
    a, b = 1 / l1, F(u) / (l1 * l2)
    tag = zh("(10.30)(10.31)", font_size=20, color=GREY_B).move_to([-5.95, 2.45, 0])
    inv1 = place_left(tx(r"X^{-1}=\frac{I}{\lambda_1}-\frac{u|\vec M\rangle\langle\vec M|}{\lambda_1\lambda_2}="
                         + frac_tex(a) + "I-" + frac_tex(b) + r"|\vec M\rangle\langle\vec M|", 36), -4.8, 2.45, 6.6)
    xi27 = [[int(v * 27) for v in row] for row in Xi]
    inv2 = place_in(tx(r"=\frac1{27}" + pmatrix(xi27), 34), -5.2, -0.6, 0.15, 1.75)
    chk = place_in(mixed(r"$XX^{-1}=X^{-1}X=I$", "（精確）", font_size=28), 0.0, 6.6, 0.65, 1.25)
    rnd1 = place_left(mixed("隨機 ", r"$\vec M$", "、", r"$n=1,\dots,6$", "：特徵值、重數 ", "$(n-1,1)$", "、左右逆全對",
                            font_size=26), -6.4, -0.35, 6.6)
    rnd2 = place_left(mixed("（誤差 ", r"$<10^{-11}$", "）", font_size=24, color=GREY_B), -6.4, -0.95, 6.6)
    glab = place_left(mixed("(10.32)：對 ", r"$\vec\sigma$", " 的高斯積分", font_size=24, color=GREY_B), -6.4, -1.55, -0.4)
    gauss = place_in(tx(r"\prod_{\vec q}\frac{\pi^{n/2}}{\sqrt{\lambda_2}\sqrt{\lambda_1^{n-1}}}\,"
                        r"e^{+\frac14\vec h\cdot X^{-1}\cdot\vec h}", 34), -6.4, -0.4, -2.9, -1.85)
    err = D.erratum_box("p. 242，(10.32)", r"\exp\Big[-\frac14\,\vec h\cdot X_{\vec q}^{-1}\cdot\vec h\Big]",
                        r"\exp\Big[+\frac14\,\vec h\cdot X_{\vec q}^{-1}\cdot\vec h\Big]", width=6.6)
    place_in(err, 0.0, 6.7, -2.95, -1.3)
    assert_in_frame(tag, inv1, inv2, chk, rnd1, rnd2, glab, gauss, err)
    assert_pairs((tag, inv1), (inv1, inv2), (inv1, chk), (inv2, chk), (inv2, rnd1), (chk, rnd1), (rnd1, rnd2),
                 (rnd2, glab), (glab, gauss), (gauss, err), (glab, err), (rnd1, err), (rnd2, err))
    g, gb = f4(nb["gn2"][0]), f4(nb["gn2"][2])

    s.play(FadeIn(tag), Write(inv1), run_time=1.5)
    s.say("(10.30)：", "$X$", " 的逆矩陣只要兩個特徵值", wait=2.0)
    s.play(FadeIn(inv2), FadeIn(chk), run_time=1.0)
    s.say(rf"$1/\lambda_1={frac_tex(a)}$", "、", rf"$u/(\lambda_1\lambda_2)={frac_tex(b)}$", "：",
          r"$XX^{-1}=X^{-1}X=I$", "（精確）", wait=2.5)
    s.play(FadeIn(rnd1), FadeIn(rnd2), run_time=0.8)
    s.say("Ex. 10.3、10.4：隨機 ", r"$\vec M$", "、", r"$n=1,\dots,6$", " 全都對", wait=2.5)
    s.play(FadeIn(glab), Write(gauss), run_time=1.2)
    s.say("(10.32)：高斯積分，指數是 ", r"$+\frac14\vec h\cdot X^{-1}\cdot\vec h$", wait=2.5)
    s.play(FadeIn(err), run_time=0.9)
    s.say("$n=2$", " 的數值積分 ", f"${g}$", "，書上的負號得 ", f"${gb}$", wait=2.5)
    s.clear_stage()


def part_n_to_zero(s):
    nb = numbers()
    l1, l2 = WK.sw_lambdas()
    head = ["n", zh("橫向重數", font_size=24, color=YELLOW), zh("縱向重數", font_size=24, color=YELLOW),
            r"\det X=\lambda_1^{n-1}\lambda_2", r"(\det X)^{-1/2}"]
    rows = [head] + [[str(n), str(t), str(lo), frac_tex(det, small=True), f4(inv)] for n, t, lo, det, inv in nb["modes"]]
    nc = 5
    tab = D.tex_table(rows, [1.2, 2.1, 2.1, 3.4, 2.6], row_h=0.66, font_size=30).shift(np.array([-5.7, 2.95, 0]))
    cells = tab[0]
    neg = SurroundingRectangle(VGroup(*cells[4 * nc:5 * nc]), buff=0.06).set_stroke(LONG_COLOR, 2.5)
    inv0 = f4(nb["modes"][-1][4])
    lim = place_in(mixed(r"$n\to0$", "：", rf"$(\det X)^{{-1/2}}={inv0}=\sqrt{{\lambda_1/\lambda_2}}=1/\sqrt3$",
                         font_size=28), -6.4, 6.6, -1.1, -0.55)
    book = f4(float(l2 / l1) ** 0.5)
    err = D.erratum_box("p. 242，(10.33)",
                        r"\sqrt{\frac{\lambda_2}{\lambda_1}}=\sqrt{\frac{r+q^2+u|\vec M|^2}{r+q^2+3u|\vec M|^2}}",
                        r"\sqrt{\frac{\lambda_1}{\lambda_2}}=\sqrt{\frac{r+q^2+u|\vec M|^2}{r+q^2+3u|\vec M|^2}}",
                        width=9.6)
    place_in(err, -4.9, 4.9, -2.9, -1.3)
    assert_in_frame(tab, neg, lim, err)
    assert_table(tab, nc)
    assert_pairs((tab, lim), (neg, lim), (lim, err))

    s.play(FadeIn(tab[1]), FadeIn(VGroup(*cells[:nc])), run_time=0.8)
    s.say(r"$n\to0$", "：每個 ", r"$\vec q$", " 的模態數怎麼變？", wait=2.0)
    s.play(FadeIn(VGroup(*cells[nc:4 * nc]), lag_ratio=0.02), run_time=1.0)
    s.say("$n=3$", "：兩個橫向一個縱向；", "$n=1$", "：只有縱向", wait=2.5)
    s.play(FadeIn(VGroup(*cells[4 * nc:])), ShowCreation(neg), run_time=1.0)
    s.say("$n=0$", "：橫向重數 ", "$-1$", "，負一個模態", wait=2.0)
    s.say("數學上允許，物理上奇妙", wait=2.0)
    s.play(FadeIn(lim), run_time=0.8)
    s.say(r"$(\det X)^{-1/2}\to\sqrt{\lambda_1/\lambda_2}=1/\sqrt3$", wait=2.5)
    s.play(FadeIn(err), run_time=0.9)
    s.say("右邊對、左邊顛倒：書上的左邊是 ", rf"$\sqrt3={book}$", wait=2.5)
    s.say("每個 ", r"$\vec q$", " 有 ", "$n-1$", " 個橫向、1 個縱向模態", wait=2.5)
    s.say(r"$n\to0$", " 時橫向變成 ", "$-1$", " 個，留下 ", r"$\sqrt{\lambda_1/\lambda_2}$", wait=2.5)


class W1003SpinWaves(ChuScene):
    def construct(self):
        self.show_header("W1003　自旋波：在平均場上加漣漪")
        check_math()
        part_expand(self)
        part_wall(self)
        part_modes_matrix(self)
        part_inverse(self)
        part_n_to_zero(self)
        self.wait(1.0)
