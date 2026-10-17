"""
W304 二維：也會回家，但慢得驚人（§3.5，pp. 58–60，(3.25)–(3.30)，Fig. 3.2）

問：在棋盤街道上亂走，100 步內還沒回過家的人剩多少？10000 步呢？
  1. 正方格子回到原點的走法 C(N;0) = 1, 4, 36, 400, 4900, 63504（N = 0, 2, …, 10），u_N = C/4^N；
     「回家的路線＝第一次回家 × 之後再回家」：u_N = Σ f_k u_{N−k}（例：f_4 = 9/64 − 1/16 = 5/64）→
     f_2..f_10 = 0.25, 0.0781, 0.0430, 0.0286, 0.0210，和 0.4208 → P_10 = 0.5792。
  2. 生成函數：(3.28) 只留 q 小的部分；d = 2 (3.30) 在 z_c 對數發散 → P_N ∝ 1/ln N；
     Fig. 3.2：z 平面、鞍點 z_c(1 − 1/N)、垂直的最陡下降路徑；積分 ≈ 常數 × z_c^{−N}/ln N。
  3. 表：精確／π/ln N（N = 10 時 1.3644，大於 1！）／π/(ln N + ln 8 + γ)；對照 W301 的模擬：1000 − 676 = 324 人。
  4. 先猜再揭曉：降到 1% 要 N = e^{100π} ≈ 10^136 步，比宇宙裡的原子數（約 10^80）還多。
收束：二維也一定回家，但「還沒回家」只像 1/ln N 那樣消失。畫面上的數字都由 utils.walk_rec 算出，check_math() 核對。
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rec as R
from utils import walk_gf as GF
from utils import lattice_draw as L
from utils import walk_draw as D

NS2 = tuple(range(0, 11, 2))
CN = [R.returns_count("sq", N) for N in NS2]               # C(N;0)
FS = R.first_return_probs("sq", 10)                        # f_0..f_10
TABLE_N = (10, 100, 1000, 10000)
NSTR = {10: "10", 100: "100", 1000: "1000", 10000: "10^4"}
Q_STEPS = 100                                              # 開頭的足跡：鎖定足跡的前 100 步
_P = {}


def P(N):
    """P_N：N 步都沒回過家的機率（精確）。"""
    if "sq" not in _P:
        _P["sq"] = R.no_return_probs("sq", 10000)
    return _P["sq"][N]


def f4(x):
    return f"{x:.4f}"


def u_str(N, c):
    return "1" if N == 0 else f4(c / 4 ** N)


def check_math():
    assert CN == [1, 4, 36, 400, 4900, 63504] == [GF.returns_sq(N) for N in NS2]
    assert CN == [math.comb(N, N // 2) ** 2 for N in NS2]
    u = [Fraction(c, 4 ** N) for c, N in zip(CN, NS2)]
    assert [u_str(N, c) for N, c in zip(NS2, CN)] == ["1", "0.2500", "0.1406", "0.0977", "0.0748", "0.0606"]
    assert [f4(FS[n]) for n in (2, 4, 6, 8, 10)] == ["0.2500", "0.0781", "0.0430", "0.0286", "0.0210"]
    assert all(FS[n] == 0 for n in range(1, 11, 2)) and FS[0] == 0
    assert f4(sum(FS)) == "0.4208" and f4(1 - sum(FS)) == "0.5792" == f4(P(10))
    uf = R.return_probs("sq", 10)                                             # 更新方程 u_N = Σ f_k u_{N−k}
    assert all(abs(uf[n] - sum(FS[k] * uf[n - k] for k in range(1, n + 1))) < 1e-14 for n in range(1, 11))
    f2e = u[1]
    f4e = u[2] - f2e * u[1]                                                   # u_4 = f_2 u_2 + f_4 u_0
    assert f2e == Fraction(1, 4) and u[2] == Fraction(9, 64) and f2e * u[1] == Fraction(1, 16)
    assert f4e == Fraction(5, 64) and abs(float(f4e) - FS[4]) < 1e-15
    assert [f4(P(N)) for N in TABLE_N] == ["0.5792", "0.4168", "0.3218", "0.2613"]
    assert [f4(R.no_return_sq_book(N)) for N in TABLE_N] == ["1.3644", "0.6822", "0.4548", "0.3411"]
    assert [f4(R.no_return_sq_refined(N)) for N in TABLE_N] == ["0.6335", "0.4326", "0.3285", "0.2647"]
    assert R.no_return_sq_book(10) > 1 and round(100 * P(10000)) == 26
    assert 1000 - R.SIM_RETURNED["sq"] == 324 and round(1000 * P(1000)) == 322
    assert f"{R.log10_steps_for(0.01):.1f}" == "136.4"
    assert abs(R.log10_steps_for(0.01) - 100 * math.pi / math.log(10)) < 1e-9
    assert int(R.log10_steps_for(0.01)) == 136 > 80
    tr = R.trail("sq", R.TRAIL_N, R.TRAIL_SEED)
    assert len(tr) == R.TRAIL_N + 1 and tuple(tr[0]) == (0, 0) and tuple(tr[-1]) == R.TRAIL_END


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def left_at(m, x, y):
    m.move_to(np.array([x + m.get_width() / 2, y, 0]))
    return m


def house(p, w=0.34, h=0.22, roof=0.18, color=YELLOW):
    p = np.array(p, dtype=float)
    pts = [p + w / 2 * LEFT, p + w / 2 * RIGHT, p + w / 2 * RIGHT + h * UP,
           p + (h + roof) * UP, p + w / 2 * LEFT + h * UP]
    return Polygon(*pts).set_fill(color, 0.85).set_stroke(WHITE, 1.5)


# ---------------------------------------------------------------- 各段
def part_question(s):
    pts0 = R.trail("sq", R.TRAIL_N, R.TRAIL_SEED)[:Q_STEPS + 1]
    xs, ys = [p[0] for p in pts0], [p[1] for p in pts0]
    cx, cy = (min(xs) + max(xs)) // 2, (min(ys) + max(ys)) // 2         # 方格置中在足跡上
    pts = [(x - cx, y - cy) for x, y in pts0]
    m = max(max(abs(x), abs(y)) for x, y in pts)
    cell = min(0.45, 2.7 / (m + 1.5))
    org = np.array([-3.4, 0.1, 0])
    grid = L.lattice_lines(m + 1, cell, org)
    home = house(L.lat_point(-cx, -cy, cell, org) + 0.08 * UP)
    tr = L.trail(pts, cell, org, stroke_width=3)
    q1 = mixed("100 步內還沒回過家的人：", "$?$", font_size=28)
    q2 = mixed("10000 步呢：", "$?$", font_size=28)
    note = zh("（一個人走 100 步的足跡：綠＝起點，紅＝終點）", font_size=20, color=GREY_B)
    col = VGroup(q1, q2, note).arrange(DOWN, buff=0.6, aligned_edge=LEFT)
    fit(col, 6.0)
    left_at(col, 0.6, 0.3)
    assert_in_frame(grid, home, tr, col)
    assert_apart(VGroup(grid, tr), col, 0.2)

    s.say("問：在棋盤街道上亂走，100 步內還沒回過家的人剩多少？", wait=0.3)
    s.play(FadeIn(grid), FadeIn(home), run_time=0.8)
    s.play(ShowCreation(tr), run_time=4.0)
    s.play(FadeIn(q1), FadeIn(note), run_time=0.8)
    s.wait(1.5)
    s.say("10000 步呢？", wait=0.3)
    s.play(FadeIn(q2), run_time=0.6)
    s.say("先猜兩個比例！", wait=3.5)
    s.clear_stage()


def part_renewal(s):
    rows = [["N", "C(N;0)", "u_N", "f_N"]] + \
           [[str(N), str(c), u_str(N, c), f4(FS[N]) if N else "0"] for N, c in zip(NS2, CN)]
    tab = D.tex_table(rows, [0.8, 1.6, 1.5, 1.5], row_h=0.5, font_size=28)
    tab.move_to(np.array([-6.6 + tab.get_width() / 2, 3.0 - tab.get_height() / 2, 0]))
    cell = lambda i, j: D.table_cell(tab, i, j, 4)
    s1 = Tex(r"f_2+f_4+\cdots+f_{10}=" + f4(sum(FS)), font_size=34)
    s2 = Tex(r"P_{10}=1-" + f4(sum(FS)) + "=" + f4(P(10)), font_size=38).set_color(YELLOW)
    left_at(s1, -6.4, -1.2)
    left_at(s2, -6.4, -2.2)
    y0 = 1.55
    p0, pk, pN = np.array([0.6, y0, 0]), np.array([2.8, y0, 0]), np.array([6.0, y0, 0])
    axis = Line(p0 + 0.2 * LEFT, pN + 0.3 * RIGHT).set_stroke(GREY_B, 2)
    ticks = VGroup(*[Line(p + 0.08 * DOWN, p + 0.08 * UP).set_stroke(GREY_B, 2) for p in (p0, pk, pN)])
    tl = VGroup(*[Tex(t, font_size=30).set_color(GREY_B).next_to(p, DOWN, buff=0.15)
                  for t, p in (("0", p0), ("k", pk), ("N", pN))])
    arc1 = ArcBetweenPoints(p0, pk, angle=-0.75 * PI).set_stroke(ORANGE, 4)
    arc2 = ArcBetweenPoints(pk, pN, angle=-0.6 * PI).set_stroke(BLUE_C, 4)
    lab1 = mixed("第一次回家 ", "$f_k$", font_size=22, color=ORANGE).next_to(arc1, UP, buff=0.12)
    lab2 = mixed("之後第 ", "$N$", " 步也在家 ", "$u_{N-k}$", font_size=22, color=BLUE_C).next_to(arc2, UP, buff=0.12)
    sch = VGroup(axis, ticks, tl, arc1, arc2, lab1, lab2)
    ren = VGroup(Tex(r"u_N=\sum_{k\ge 1}f_k\,u_{N-k}", font_size=40),
                 zh("（更新方程）", font_size=20, color=GREY_B)).arrange(RIGHT, buff=0.25)
    ren.move_to(np.array([3.3, 0.3, 0]))
    e1 = Tex(r"u_4=f_2u_2+f_4u_0", font_size=34)
    e2 = Tex(r"f_4=\frac{9}{64}-\frac{1}{4}\cdot\frac{1}{4}=\frac{5}{64}=" + f4(FS[4]), font_size=34)
    fit(e2, 6.4)
    left_at(e1, 0.3, -0.95)
    left_at(e2, 0.3, -2.0)
    assert_in_frame(tab, s1, s2, sch, ren, e1, e2)
    assert_apart(tab, s1, 0.15)
    assert_apart(VGroup(tab, s1, s2), VGroup(sch, ren, e1, e2), 0.2)
    assert_apart(sch, ren, 0.1)
    assert_apart(ren, e1, 0.1)
    assert_apart(e1, e2, 0.1)
    assert_apart(lab1, lab2, 0.1)

    s.say("正方格子：走 ", "$N$", " 步回到原點的走法數 ", "$C(N;0)$", wait=0.3)
    s.play(FadeIn(tab[1]), *[FadeIn(cell(0, j)) for j in range(4)], run_time=0.6)
    s.play(LaggedStart(*[FadeIn(VGroup(cell(i, 0), cell(i, 1))) for i in range(1, 7)], lag_ratio=0.25), run_time=2.0)
    s.wait(1.5)
    s.say("$u_N=C(N;0)/4^N$", "：第 ", "$N$", " 步剛好在家的機率", wait=0.3)
    s.play(LaggedStart(*[FadeIn(cell(i, 2)) for i in range(1, 7)], lag_ratio=0.25), run_time=1.5)
    s.wait(2.0)
    s.say("回家的路線＝第一次回家 × 之後再回家", wait=0.3)
    s.play(ShowCreation(axis), FadeIn(ticks), FadeIn(tl), run_time=0.8)
    s.play(ShowCreation(arc1), FadeIn(lab1), run_time=1.0)
    s.play(ShowCreation(arc2), FadeIn(lab2), run_time=1.0)
    s.play(FadeIn(ren), run_time=1.0)
    s.wait(2.5)
    s.say("例：", "$u_4=f_2u_2+f_4u_0$", "，解出第一次在第 4 步回家的機率", wait=0.3)
    s.play(FadeIn(e1), run_time=0.8)
    s.play(FadeIn(e2), run_time=1.0)
    s.wait(2.5)
    s.say("一個一個解出 ", "$f_N$", "：第一次回家在第 ", "$N$", " 步", wait=0.3)
    s.play(LaggedStart(*[FadeIn(cell(i, 3)) for i in range(1, 7)], lag_ratio=0.25), run_time=1.5)
    s.wait(1.5)
    s.say("加起來：10 步內回過家的機率", wait=0.3)
    s.play(FadeIn(s1), run_time=0.8)
    s.wait(1.5)
    s.say("所以 10 步內還沒回家：", f"$P_{{10}}={f4(P(10))}$", wait=0.3)
    s.play(FadeIn(s2), run_time=0.8)
    s.wait(2.5)
    s.clear_stage()


def part_saddle(s):
    e1 = mixed("(3.27)：", "$q$", " 小時 ", r"$1-z\chi(\vec q)\approx c_1(z_c-z)+c_2q^2/c_1$", font_size=24)
    e2 = mixed("(3.28)：只留 ", "$q<q_0$", " 的部分（奇異性都在這裡）", font_size=24)
    e3 = VGroup(Tex(r"G(z;0)\propto\frac{c_1^2}{2c_2}\ln\frac{z_c-z+c_2q_0^2/c_1^2}{z_c-z}", font_size=36),
                zh("(3.30)，d = 2", font_size=20, color=GREY_B)).arrange(RIGHT, buff=0.25)
    e4 = mixed("$z\\to z_c$", "：分母 ", "$z_c-z\\to 0$", "，對數發散", font_size=26)
    e5 = Tex(r"P_N\propto\frac{1}{\ln N}", font_size=44).set_color(YELLOW)
    left = VGroup(e1, e2, e3, e4, e5).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    fit(left, 7.0)
    left_at(left, -6.5, 3.1 - left.get_height() / 2)
    ax = D.plot_axes(x_range=(-1.2, 2.0, 1), y_range=(-1.2, 1.2, 1), width=5.4, height=4.2,
                     x_label=r"\mathrm{Re}\,z", y_label=r"\mathrm{Im}\,z")
    ax.shift(np.array([2.9, 0.2, 0]) - ax.c2p(0, 0))
    zc = Dot(ax.c2p(1, 0), radius=0.07).set_fill(WHITE, 1)
    zcl = Tex("z_c", font_size=30).next_to(ax.c2p(1, 0), UP, buff=0.15)
    cut = L.branch_cut(ax, 1.0, 1.55)
    circ = D.contour_circle(ax, 0.8, color=TEAL)
    sad = Dot(ax.c2p(0.8, 0), radius=0.08).set_fill(ORANGE, 1)
    sadl = Tex(r"z_c(1-1/N)", font_size=24).set_color(ORANGE).next_to(ax.c2p(0.8, 0), DR, buff=0.15)
    path = Line(ax.c2p(0.8, -0.45), ax.c2p(0.8, 0.45)).set_stroke(ORANGE, 6)
    pathl = zh("最陡下降", font_size=20, color=ORANGE).next_to(ax.c2p(0.8, 0.45), RIGHT, buff=0.15)
    title = mixed("Fig. 3.2：", "$z$", " 平面", font_size=22, color=GREY_B).move_to(np.array([3.5, 2.9, 0]))
    fig = VGroup(ax, zc, zcl, cut, circ, sad, sadl, path, pathl, title)
    assert_in_frame(left, fig)
    assert_apart(left, fig, 0.2)
    assert_apart(sadl, pathl, 0.05)
    assert_apart(sadl, zcl, 0.05)
    assert_apart(sadl, cut, 0.02)
    assert_apart(zcl, pathl, 0.05)

    s.say("二維：用生成函數看 ", "$N$", " 很大的時候", wait=0.3)
    s.play(FadeIn(e1), run_time=0.8)
    s.wait(1.5)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("只有 ", "$q$", " 小的地方會出事：分母在 ", "$z=z_c$", " 變成 0", wait=2.5)
    s.play(Write(e3[0]), FadeIn(e3[1]), run_time=2.0)
    s.wait(2.0)
    s.play(FadeIn(e4), run_time=0.8)
    s.say("在 ", "$z_c$", " 對數發散：還沒回家的機率像 ", r"$1/\ln N$", wait=0.3)
    s.play(FadeIn(e5), run_time=0.8)
    s.wait(2.0)
    s.say("怎麼取係數？看 ", "$z$", " 平面（Fig. 3.2）", wait=0.3)
    s.play(FadeIn(ax), FadeIn(title), run_time=0.8)
    s.play(FadeIn(zc), FadeIn(zcl), ShowCreation(cut), run_time=0.8)
    s.play(ShowCreation(circ), run_time=1.2)
    s.play(FadeIn(sad), FadeIn(sadl), ShowCreation(path), FadeIn(pathl), run_time=1.0)
    s.say("鞍點在 ", "$z_c(1-1/N)$", "，積分 ", r"$\approx$", " 常數 ", r"$\times z_c^{-N}/\ln N$", wait=4.0)
    s.say("除以全部的 ", "$z_c^{-N}=4^N$", "：還沒回家的機率 ", r"$\propto 1/\ln N$", wait=3.0)
    s.clear_stage()


def part_table(s):
    book10 = VGroup(Tex(f4(R.no_return_sq_book(10)), font_size=30).set_color(RED),
                    zh("大於 1！", font_size=18, color=RED)).arrange(RIGHT, buff=0.1)
    hdr = ["N", zh("精確", font_size=24, color=YELLOW), r"\pi/\ln N", r"\pi/(\ln N+\ln 8+\gamma)"]
    rows = [hdr]
    for N in TABLE_N:
        rows.append([NSTR[N], f4(P(N)), book10 if N == 10 else f4(R.no_return_sq_book(N)),
                     f4(R.no_return_sq_refined(N))])
    tab = D.tex_table(rows, [1.4, 1.8, 2.8, 4.0], row_h=0.66, font_size=30)
    tab.move_to(np.array([0, 3.0 - tab.get_height() / 2, 0]))
    cell = lambda i, j: D.table_cell(tab, i, j, 4)
    box = SurroundingRectangle(book10, buff=0.08).set_stroke(RED, 3)
    c1 = mixed("對照 W301 的模擬：1000 步內還沒回家 ", f"$1000-{R.SIM_RETURNED['sq']}={1000 - R.SIM_RETURNED['sq']}$",
               " 人", font_size=26)
    c2 = mixed("理論：", f"$P_{{1000}}={f4(P(1000))}$", f"，約 {round(1000 * P(1000))} 人", font_size=26)
    low = VGroup(c1, c2).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    fit(low, 12.0)
    low.move_to(np.array([0, -1.6, 0]))
    assert_in_frame(tab, box, low)
    assert_apart(tab, low, 0.3)

    s.say("精確值和兩個近似比一比", wait=0.3)
    s.play(FadeIn(tab[1]), *[FadeIn(cell(0, j)) for j in range(4)], run_time=0.8)
    s.play(LaggedStart(*[FadeIn(VGroup(cell(i, 0), cell(i, 1))) for i in range(1, 5)], lag_ratio=0.3), run_time=1.8)
    s.wait(1.5)
    s.play(LaggedStart(*[FadeIn(cell(i, 2)) for i in range(1, 5)], lag_ratio=0.3), run_time=1.5)
    s.say("書上的首項 ", r"$\pi/\ln N$", "：", "$N=10$", " 時大於 1，不可能是機率", wait=0.3)
    s.play(ShowCreation(box), run_time=0.8)
    s.wait(2.5)
    s.say("補上常數項 ", r"$\ln 8+\gamma$", "，就準多了", wait=0.3)
    s.play(LaggedStart(*[FadeIn(cell(i, 3)) for i in range(1, 5)], lag_ratio=0.3), run_time=1.5)
    s.wait(2.5)
    s.say("一萬步還有 ", f"${round(100 * P(10000))}\\%$", " 的人沒回過家", wait=2.5)
    s.play(FadeIn(c1), run_time=0.8)
    s.play(FadeIn(c2), run_time=0.8)
    s.say("模擬和理論對得上", wait=2.5)
    s.clear_stage()


def part_guess(s):
    X0, X1, Y, EMAX = -6.0, 6.0, -1.0, 140
    xr = lambda e: X0 + (X1 - X0) * e / EMAX
    line = Arrow(np.array([X0 - 0.2, Y, 0]), np.array([X1 + 0.5, Y, 0]), buff=0).set_color(GREY_B)
    ticks = VGroup(*[Line(np.array([xr(e), Y - 0.08, 0]), np.array([xr(e), Y + 0.08, 0])).set_stroke(GREY_B, 2)
                     for e in range(0, EMAX + 1, 20)])
    tlab = VGroup(*[Tex(rf"10^{{{e}}}", font_size=24).set_color(GREY_B).next_to(np.array([xr(e), Y, 0]), DOWN, buff=0.18)
                    for e in range(0, EMAX + 1, 20)])
    lg = R.log10_steps_for(0.01)

    def marker(e, label, color):
        p = np.array([xr(e), Y + 0.12, 0])
        arr = Arrow(p + 0.75 * UP, p, buff=0).set_color(color)
        lab = label.next_to(arr, UP, buff=0.1)
        return arr, lab

    a1, l1 = marker(4, mixed("1 萬步：還有 ", f"${round(100 * P(10000))}\\%$", font_size=22, color=TEAL), TEAL)
    a2, l2 = marker(80, mixed("宇宙的原子數 ", r"$\approx 10^{80}$", font_size=22, color=BLUE_C), BLUE_C)
    a3, l3 = marker(lg, mixed("降到 ", r"$1\%$", font_size=26, color=YELLOW), YELLOW)
    l1.align_to(np.array([X0 - 0.4, 0, 0]), LEFT)
    l3.align_to(np.array([6.6, 0, 0]), RIGHT)
    qm = Tex("?", font_size=80).move_to(np.array([3.0, 0.4, 0]))
    f = Tex(r"\frac{\pi}{\ln N}=0.01\ \Rightarrow\ \ln N=100\pi\ \Rightarrow\ N=e^{100\pi}\approx 10^{136}",
            font_size=40)
    fit(f, 12.6)
    f.move_to(np.array([0, 2.3, 0]))
    assert_in_frame(line, ticks, tlab, a1, l1, a2, l2, a3, l3, qm, f)
    assert_apart(l1, l2, 0.2)
    assert_apart(l2, l3, 0.2)
    assert_apart(f, VGroup(l1, l2, l3, qm), 0.2)

    s.say("先猜：要走幾步，「還沒回家」才會降到 ", r"$1\%$", "？", wait=0.3)
    s.play(GrowArrow(line), FadeIn(ticks), FadeIn(tlab), run_time=1.2)
    s.play(GrowArrow(a1), FadeIn(l1), run_time=0.8)
    s.play(FadeIn(qm), run_time=0.5)
    s.wait(3.5)
    s.say("解 ", r"$\pi/\ln N=0.01$", wait=0.3)
    s.play(Write(f), run_time=2.0)
    s.wait(1.5)
    s.say(r"$\log_{10}N=100\pi/\ln 10$", f"$={lg:.1f}$", "：", "$N\\approx 10^{136}$", " 步", wait=0.3)
    s.play(FadeOut(qm), GrowArrow(a3), FadeIn(l3), run_time=1.0)
    s.wait(2.0)
    s.say("比宇宙裡的原子數（約 ", "$10^{80}$", "）還多！", wait=0.3)
    s.play(GrowArrow(a2), FadeIn(l2), run_time=1.0)
    s.wait(3.0)


def part_close(s):
    s.say("二維也一定回家，但「還沒回家」只像 ", r"$1/\ln N$", " 那樣消失", wait=3.0)
    s.say("慢得驚人", wait=2.0)


class W304TwoDimensions(ChuScene):
    def construct(self):
        self.show_header("W304　二維：也會回家，但慢得驚人")
        check_math()
        part_question(self)
        part_renewal(self)
        part_saddle(self)
        part_table(self)
        part_guess(self)
        part_close(self)
        self.wait(1.0)
