"""
W303 一維：一定回家（§3.4 例題，pp. 55–57，(3.11)–(3.24)，Ex. 3.2–3.4）

問：十步醉漢，十步之內從來沒回過家的走法有幾條？（先猜）
  1. 令 w⃗ = x⃗：從不回家（終點不限）的生成函數 = g(z;0)/G(z;0)（(3.12)(3.13)(3.16)）；g = 1/(1−2z)、
     G(z;0) = 1/√(1−4z²) = 1 + 2z² + 6z⁴ + 20z⁶ + 70z⁸ + 252z¹⁰ + ⋯（(3.21)）；
     相除 = √((1+2z)/(1−2z)) = 1 + 2z + 2z² + 4z³ + ⋯ + 252z¹⁰ + ⋯（表：一列一個 N，旁邊 2^N）；
     驚喜：252 = C(10,5) = 第 10 步剛好在家的條數（N = 2k：從不回家 = C(2k,k)）。
  2. P_10 = 252/1024 = 0.2461；(3.22) g/G → 1/√(z_c − z)，(3.23)(3.24) P_N ≈ √(2/(πN))；表：精確／更正／書上；
     筆誤框：書上少了 z_c^{−1/2} = √2。
  3. 圖：精確 P_N（N = 2, 4, …, 100）、√(2/πN)（綠）、1/√(πN)（紅虛線）→ P_N → 0：一維一定回家。
  4. 偏硬幣（練習 3.2–3.4）：p = 0.6：G(1;0) = 5 → 逃脫機率 0.2 = |p − q|；P_N（N = 1..100）的長條 → 0.2；
     回家機率 0.8；p = 1：逃脫機率 1。
收束：公平的一維醉漢一定回家，但「還沒回家」只像 1/√N 那樣慢慢消失。
畫面上的數字都由 utils.walk_rec 算出，check_math() 核對。
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

NR = R.never_return_1d_gf(10)                              # [z^N] g/G：N 步都不回家的走法數
GCOEF = [R.returns_count("1d", 2 * k) for k in range(6)]   # G(z;0) 的 z^{2k} 係數
NS = (10, 100, 1000)
PB_N = (1, 2, 4, 10, 20, 50, 100)
P_BIAS = 0.6
WALK_A = [1, 1, -1, 1, 1, -1, -1, 1, 1, -1]                # 從沒回家（綠）
WALK_B = [-1, 1, 1, -1, -1, -1, 1, 1, 1, 1]                # 回過家（灰）
UX, UY = 0.55, 0.36


def positions(steps):
    out = [0]
    for st in steps:
        out.append(out[-1] + st)
    return out


def fmt(N, x):
    return f"{x:.5f}" if N >= 1000 else f"{x:.4f}"


def pb_fmt(x):
    return f"{round(x, 2):g}" if abs(x - round(x, 2)) < 1e-12 else f"{x:.4f}"


def series(coeffs, powers):
    terms = []
    for c, p in zip(coeffs, powers):
        z = "" if p == 0 else ("z" if p == 1 else f"z^{{{p}}}")
        terms.append((str(c) if (c != 1 or p == 0) else "") + z)
    return "+".join(terms)


def check_math():
    assert NR == [1, 2, 2, 4, 6, 12, 20, 40, 70, 140, 252] == [R.never_return_1d(N) for N in range(11)]
    assert all(R.never_return_1d_brute(N) == NR[N] for N in range(11))
    assert GCOEF == [1, 2, 6, 20, 70, 252] == [math.comb(2 * k, k) for k in range(6)]
    # 第二種核對：(Σ NR z^N)² = (1+2z)/(1−2z)、(Σ GCOEF z^{2k})² = 1/(1−4z²)
    assert [sum(NR[i] * NR[n - i] for i in range(n + 1)) for n in range(11)] == [1] + [2 ** (n + 1) for n in range(1, 11)]
    assert [sum(GCOEF[i] * GCOEF[k - i] for i in range(k + 1)) for k in range(6)] == [4 ** k for k in range(6)]
    assert NR[10] == math.comb(10, 5) == R.returns_count("1d", 10) == 252
    assert all(R.never_return_1d(2 * k) == math.comb(2 * k, k) for k in range(1, 6))
    assert R.escape_1d(10) == Fraction(252, 1024) and fmt(10, float(R.escape_1d(10))) == "0.2461"
    ex = [float(R.escape_1d(N)) for N in NS]
    co = [R.escape_1d_asym(N) for N in NS]
    bk = [R.escape_1d_book(N) for N in NS]
    assert [fmt(N, x) for N, x in zip(NS, ex)] == ["0.2461", "0.0796", "0.02523"]
    assert [fmt(N, x) for N, x in zip(NS, co)] == ["0.2523", "0.0798", "0.02523"]
    assert [fmt(N, x) for N, x in zip(NS, bk)] == ["0.1784", "0.0564", "0.01784"]
    assert all(abs(c - math.sqrt(2) * b) < 1e-15 for c, b in zip(co, bk))      # 書上少了 z_c^{−1/2} = √2
    assert f"{ex[2] / co[2]:.5f}" == "0.99975"
    assert abs(1 / math.sqrt(1 - 4 * P_BIAS * (1 - P_BIAS)) - 5) < 1e-12
    assert abs(R.escape_biased(P_BIAS) - 0.2) < 1e-12 and abs(R.escape_biased(P_BIAS) - abs(0.6 - 0.4)) < 1e-12
    PB = [R.escape_biased_series(P_BIAS, N) for N in PB_N]
    assert [pb_fmt(x) for x in PB] == ["1", "0.52", "0.4048", "0.2940", "0.2427", "0.2094", "0.2016"]
    assert all(a > b > 0.2 for a, b in zip(PB, PB[1:])) or PB[0] == 1
    assert f"{1 - R.escape_biased(P_BIAS):.1f}" == "0.8" and R.escape_biased(1.0) == 1.0
    pa, pb = positions(WALK_A), positions(WALK_B)
    assert len(WALK_A) == len(WALK_B) == 10 and all(x != 0 for x in pa[1:]) and 0 in pb[1:]
    assert max(map(abs, pa + pb)) <= 5


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tagged(mob, tag):
    return VGroup(mob, zh(tag, font_size=20, color=GREY_B)).arrange(RIGHT, buff=0.3)


def left_at(m, x, y):
    m.move_to(np.array([x + m.get_width() / 2, y, 0]))
    return m


# ---------------------------------------------------------------- 各段
def part_question(s):
    org = np.array([-6.0, 0.55, 0])
    ax = D.walk_axes(10, 5, ux=UX, uy=UY, origin=org)
    pa = D.walk_path(positions(WALK_A), ux=UX, uy=UY, origin=org, color=D.EXACT_COLOR, stroke_width=4)
    pb = D.walk_path(positions(WALK_B), ux=UX, uy=UY, origin=org, color=GREY_B, stroke_width=3)
    home = zh("家", font_size=22, color=YELLOW).next_to(org, LEFT, buff=0.15)
    r1 = mixed("十步醉漢：", "$2^{10}=1024$", " 條走法", font_size=28)
    r2 = mixed("綠：第 1–10 步都不在 ", "$0$", "（從沒回家）", font_size=22, color=D.EXACT_COLOR)
    r3 = zh("灰：回過家（碰到橫軸）", font_size=22, color=GREY_B)
    r4 = mixed("從沒回過家的走法：", "$?$", font_size=30, color=YELLOW)
    col = VGroup(r1, r2, r3, r4).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    fit(col, 5.6)
    left_at(col, 0.9, 0.6)
    assert_in_frame(ax, pa, pb, home, col)
    assert_apart(VGroup(ax, pa, pb), col)

    s.say("問：十步醉漢，十步之內從來沒回過家的走法有幾條？", wait=0.3)
    s.play(FadeIn(ax), FadeIn(home), run_time=0.8)
    s.play(ShowCreation(pb[0]), FadeIn(pb[1]), run_time=1.5)
    s.play(ShowCreation(pa[0]), FadeIn(pa[1]), run_time=1.5)
    s.play(FadeIn(r1), FadeIn(r2), FadeIn(r3), run_time=0.8)
    s.wait(1.5)
    s.play(FadeIn(r4), run_time=0.6)
    s.say("先猜一個數字！", wait=3.5)
    s.clear_stage()


def part_gf(s):
    rows = [["N", zh("從不回家", font_size=22, color=YELLOW), "2^N"]] + \
           [[str(N), str(NR[N]), str(2 ** N)] for N in range(11)]
    tab = D.tex_table(rows, [0.9, 1.7, 1.3], row_h=0.44, font_size=26)
    tab.move_to(np.array([-6.6 + tab.get_width() / 2, 3.0 - tab.get_height() / 2, 0]))
    cell = lambda i, j: D.table_cell(tab, i, j, 3)
    X0 = -2.2
    f1 = tagged(mixed("從不回家（終點不限）：", r"$\dfrac{g(z;0)}{G(z;0)}$", font_size=26), "(3.12)(3.13)(3.16)")
    f2 = Tex(r"g(z;0)=\frac{z_c}{z_c-z}=\frac{1}{1-2z}", font_size=34)
    gtex = Tex(r"G(z;0)=\frac{1}{\sqrt{1-4z^2}}=" + series(GCOEF, range(0, 11, 2)) + r"+\cdots",
               font_size=32, isolate=["252z^{10}"])
    f3 = tagged(gtex, "(3.21)")
    lhs = r"\frac{g}{G}=\sqrt{\frac{1+2z}{1-2z}}="
    f4 = Tex(lhs + series(NR[:6], range(6)), font_size=32)
    f5 = Tex("+" + series(NR[6:], range(6, 11)) + r"+\cdots", font_size=32, isolate=["252z^{10}"])
    f6 = mixed("驚喜：", r"$252=\binom{10}{5}$", "＝第 10 步剛好在家的條數", font_size=24, color=YELLOW)
    for f in (f1, f2, f3, f6):
        fit(f, 8.8)
    left_at(f1, X0, 2.55)
    left_at(f2, X0, 1.5)
    left_at(f3, X0, 0.5)
    left_at(f4, X0, -0.7)
    lhs_w = Tex(lhs, font_size=32).get_width()
    f5.move_to(np.array([f4.get_left()[0] + lhs_w + f5.get_width() / 2, -1.55, 0]))
    left_at(f6, X0, -2.4)
    fs = [f1, f2, f3, f4, f5, f6]
    assert_in_frame(tab, *fs)
    assert_apart(tab, VGroup(*fs), 0.2)
    for a, b in zip(fs, fs[1:]):
        assert_apart(a, b, 0.05)
    p252 = [gtex.select_part("252z^{10}"), f5.select_part("252z^{10}")]

    s.say("令 ", r"$\vec w=\vec x$", "：到訪起點 0 次，就是從不回家", wait=0.3)
    s.play(FadeIn(f1), run_time=1.0)
    s.wait(2.5)
    s.say("一維：", "$g$", " 數終點不限的走法，", "$G$", " 數回到原點的走法", wait=0.3)
    s.play(FadeIn(f2), run_time=0.8)
    s.play(FadeIn(f3), run_time=1.2)
    s.wait(2.5)
    s.say("相除：", r"$\sqrt{(1+2z)/(1-2z)}$", " 的係數就是從不回家的條數", wait=0.3)
    s.play(FadeIn(f4), run_time=1.2)
    s.play(FadeIn(f5), run_time=1.0)
    s.wait(1.0)
    s.say("列成表：一列一個 ", "$N$", "，旁邊是全部的 ", "$2^N$", wait=0.3)
    s.play(FadeIn(tab[1]), *[FadeIn(cell(0, j)) for j in range(3)], run_time=0.6)
    s.play(LaggedStart(*[FadeIn(VGroup(cell(i, 0), cell(i, 1), cell(i, 2))) for i in range(1, 12)],
                       lag_ratio=0.25), run_time=3.0)
    s.wait(1.5)
    box_t = SurroundingRectangle(cell(11, 1), buff=0.08).set_stroke(YELLOW, 3)
    box_r = SurroundingRectangle(p252[1], buff=0.06).set_stroke(YELLOW, 3)
    box_g = SurroundingRectangle(p252[0], buff=0.06).set_stroke(YELLOW, 3)
    s.say("十步之內從沒回過家：", f"${NR[10]}$", " 條（共 1024 條）", wait=0.3)
    s.play(ShowCreation(box_t), ShowCreation(box_r), run_time=0.8)
    s.wait(2.0)
    s.say("驚喜：", "$G(z;0)$", " 裡也有 252：第 10 步剛好在家的條數", wait=0.3)
    s.play(ShowCreation(box_g), FadeIn(f6), run_time=1.0)
    s.wait(2.5)
    s.say("一般 ", "$N=2k$", "：從不回家的條數 ", r"$=\binom{2k}{k}$", "＝第 ", "$2k$", " 步剛好在家", wait=4.0)
    s.clear_stage()


def part_asym(s):
    ex = [float(R.escape_1d(N)) for N in NS]
    co = [R.escape_1d_asym(N) for N in NS]
    bk = [R.escape_1d_book(N) for N in NS]
    p1 = Tex(rf"P_{{10}}=\frac{{{NR[10]}}}{{{2 ** 10}}}={fmt(10, ex[0])}", font_size=40)
    p2 = tagged(mixed("$N$", " 大：", r"$\dfrac{g}{G}\to\dfrac{1}{\sqrt{z_c-z}}$", font_size=26), "(3.22)")
    p3 = tagged(Tex(r"P_N\approx\sqrt{\frac{2}{\pi N}}", font_size=44).set_color(D.EXACT_COLOR), "(3.23)(3.24) 更正後")
    left = VGroup(p1, p2, p3).arrange(DOWN, buff=0.55, aligned_edge=LEFT)
    fit(left, 6.0)
    left_at(left, -6.4, 3.1 - left.get_height() / 2)
    err = D.erratum_box("p. 57，(3.23)(3.24)", r"N^{-1/2}\frac{1}{\Gamma(1/2)}",
                        r"\sqrt{2}\,N^{-1/2}\frac{1}{\Gamma(1/2)}=\sqrt{\frac{2}{\pi N}}", width=6.3, font_size=28)
    hdr = ["N", zh("精確", font_size=22, color=YELLOW),
           Tex(r"\sqrt{2/\pi N}", font_size=26).set_color(D.EXACT_COLOR),
           Tex(r"1/\sqrt{\pi N}", font_size=26).set_color(D.WRONG_COLOR)]
    rows = [hdr] + [[str(N), fmt(N, e), fmt(N, c), fmt(N, b)] for N, e, c, b in zip(NS, ex, co, bk)]
    tab = D.tex_table(rows, [1.2, 1.7, 1.8, 1.8], row_h=0.62, font_size=28)
    tab.move_to(np.array([6.5 - tab.get_width() / 2, 2.9 - tab.get_height() / 2, 0]))
    cell = lambda i, j: D.table_cell(tab, i, j, 4)
    for i in range(1, 4):
        cell(i, 2).set_color(D.EXACT_COLOR)
        cell(i, 3).set_color(D.WRONG_COLOR)
    lg = zh("綠：更正後的 (3.24)　紅：書上印的", font_size=20, color=GREY_B)
    lg.next_to(tab, DOWN, buff=0.25).align_to(tab, LEFT)
    err.next_to(lg, DOWN, buff=0.4)
    err.move_to(np.array([tab.get_center()[0], err.get_center()[1], 0]))
    assert_in_frame(left, err, tab, lg)
    assert_apart(left, VGroup(tab, lg, err), 0.15)
    assert_apart(lg, err, 0.1)

    s.say("十步之內從沒回家的機率", wait=0.3)
    s.play(Write(p1), run_time=1.5)
    s.wait(2.0)
    s.say("$N$", " 很大時，只看 ", r"$z_c=\tfrac12$", " 附近的奇異點", wait=0.3)
    s.play(FadeIn(p2), run_time=1.0)
    s.wait(2.5)
    s.say("取係數：還沒回家的機率像 ", r"$\sqrt{2/(\pi N)}$", " 那樣變小", wait=0.3)
    s.play(FadeIn(p3), run_time=1.0)
    s.wait(1.5)
    s.play(FadeIn(tab[1]), *[FadeIn(cell(0, j)) for j in range(4)], run_time=0.6)
    s.play(LaggedStart(*[FadeIn(VGroup(*[cell(i, j) for j in range(4)])) for i in range(1, 4)], lag_ratio=0.4),
           FadeIn(lg), run_time=2.0)
    s.say("精確值很快就貼上綠色的 ", r"$\sqrt{2/\pi N}$", wait=3.0)
    s.say("$N=1000$", "：精確值除以 ", r"$\sqrt{2/(\pi N)}$", f" 是 {ex[2] / co[2]:.5f}", wait=3.0)
    s.say("書上印的是紅色那一欄：整整小了 ", r"$\sqrt2$", " 倍", wait=0.3)
    s.play(FadeIn(err), run_time=1.0)
    s.wait(2.0)
    s.say("$(z_c-z)^{-1/2}$", " 的係數是 ", r"$z_c^{-N-1/2}N^{-1/2}/\Gamma(\tfrac12)$",
          "，少了 ", r"$z_c^{-1/2}=\sqrt2$", wait=4.5)
    s.clear_stage()


def part_plot(s):
    ax = D.plot_axes(x_range=(0, 100, 10), y_range=(0, 0.7, 0.1), width=8.6, height=4.6,
                     x_ticks=(20, 40, 60, 80, 100), y_ticks=(0.2, 0.4, 0.6), x_label="N", y_label="P_N")
    ax.shift(np.array([-5.6, -2.2, 0]) - ax.c2p(0, 0))
    dots = VGroup(*[Dot(ax.c2p(N, float(R.escape_1d(N))), radius=0.045).set_fill(WHITE, 1)
                    for N in range(2, 101, 2)])
    green = ax.get_graph(R.escape_1d_asym, x_range=(1.4, 100, 2)).set_stroke(D.EXACT_COLOR, 3)
    red = DashedVMobject(ax.get_graph(R.escape_1d_book, x_range=(1.0, 100, 2)), num_dashes=60)
    red.set_stroke(D.WRONG_COLOR, 3)
    lg = VGroup(
        VGroup(Dot(radius=0.06).set_fill(WHITE, 1), mixed("精確 ", "$P_N$", font_size=22)).arrange(RIGHT, buff=0.2),
        VGroup(Line(ORIGIN, 0.5 * RIGHT).set_stroke(D.EXACT_COLOR, 3),
               Tex(r"\sqrt{2/\pi N}", font_size=30)).arrange(RIGHT, buff=0.2),
        VGroup(DashedLine(ORIGIN, 0.5 * RIGHT, dash_length=0.08).set_stroke(D.WRONG_COLOR, 3),
               Tex(r"1/\sqrt{\pi N}", font_size=30)).arrange(RIGHT, buff=0.2),
    ).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
    left_at(lg, 3.9, 1.6)
    assert_in_frame(ax, dots, green, red, lg)
    assert_apart(ax, lg, 0.2)

    s.say("把 ", "$P_N$", " 畫出來（", r"$N=2,4,\dots,100$", "）", wait=0.3)
    s.play(FadeIn(ax), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(d) for d in dots], lag_ratio=0.05), FadeIn(lg[0]), run_time=2.5)
    s.play(ShowCreation(green), FadeIn(lg[1]), run_time=1.5)
    s.play(ShowCreation(red), FadeIn(lg[2]), run_time=1.5)
    s.say("點落在綠線上，不在紅線上", wait=2.5)
    s.say("$P_N\\to 0$", "：一維醉漢一定回家", wait=3.0)
    s.clear_stage()


def part_biased(s):
    PB = [R.escape_biased_series(P_BIAS, N) for N in PB_N]
    esc = R.escape_biased(P_BIAS)
    b1 = mixed("偏硬幣：向右 ", "$p=0.6$", "，向左 ", "$q=0.4$", font_size=26)
    b2 = Tex(r"G(1;0)=\frac{1}{\sqrt{1-4pq}}=\frac{1}{\sqrt{0.04}}=5", font_size=36)
    b3 = mixed("逃脫機率 ", rf"$=\frac{{1}}{{5}}={esc:.1f}=|p-q|$", font_size=26)
    b4 = mixed("回家機率 ", rf"$1-|p-q|={1 - esc:.1f}$", font_size=26)
    b5 = mixed("$p=1$", "：一直往右，逃脫機率 ", f"${R.escape_biased(1.0):g}$", font_size=26)
    tag = zh("練習 3.2–3.4", font_size=20, color=GREY_B)
    left = VGroup(b1, b2, b3, b4, b5, tag).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    fit(left, 6.2)
    left_at(left, -6.4, 0.15)
    org = np.array([0.75, -2.0, 0])
    ux, uy = 0.85, 3.4
    ch = D.bar_chart(dict(enumerate(PB)), ux=ux, uy=uy, origin=org, labels=True, label_fmt=pb_fmt, font_size=22)
    base = Line(org + 0.5 * LEFT, org + (6 * ux + 0.5) * RIGHT).set_stroke(GREY_B, 2)
    nl = VGroup(*[Tex(str(N), font_size=24).set_color(GREY_B).next_to(org + k * ux * RIGHT, DOWN, buff=0.15)
                  for k, N in enumerate(PB_N)])
    nlab = Tex("N", font_size=28).set_color(GREY_B).move_to(np.array([6.5, nl.get_center()[1], 0]))
    title = mixed("$P_N$", "：", "$N$", " 步都沒回家的機率", font_size=24).move_to(np.array([3.3, 2.3, 0]))
    y02 = org[1] + esc * uy
    dash = DashedLine(np.array([0.3, y02, 0]), np.array([6.2, y02, 0]), dash_length=0.1).set_stroke(ORANGE, 2.5)
    d02 = Tex(f"{esc:.1f}", font_size=28).set_color(ORANGE).next_to(dash, RIGHT, buff=0.08)
    assert_in_frame(left, ch, base, nl, nlab, title, dash, d02)
    assert_apart(left, VGroup(ch, base, nl, title, dash), 0.2)
    assert_apart(title, ch, 0.1)

    s.say("如果硬幣不公平：向右的機率 ", "$p=0.6$", wait=0.3)
    s.play(FadeIn(b1), run_time=0.8)
    s.say("回到原點的生成函數，在 ", "$z=1$", " 是有限的 5", wait=0.3)
    s.play(Write(b2), run_time=1.5)
    s.wait(2.0)
    s.say("逃脫機率＝", "$1/G(1;0)=0.2=|p-q|$", wait=0.3)
    s.play(FadeIn(b3), run_time=0.8)
    s.wait(1.5)
    s.say("有限步的 ", "$P_N$", "：一路往下，停在 0.2", wait=0.3)
    s.play(ShowCreation(base), FadeIn(nl), FadeIn(nlab), FadeIn(title), run_time=0.8)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in ch[0]], lag_ratio=0.2), run_time=2.0)
    s.play(FadeIn(ch[1]), run_time=0.6)
    s.play(ShowCreation(dash), FadeIn(d02), run_time=0.8)
    s.wait(2.0)
    s.say("回家機率只有 0.8：有 20% 的人一去不回", wait=0.3)
    s.play(FadeIn(b4), run_time=0.8)
    s.wait(2.0)
    s.say("$p=1$", "：每一步都往右，永遠回不了家", wait=0.3)
    s.play(FadeIn(b5), FadeIn(tag), run_time=0.8)
    s.wait(2.5)


def part_close(s):
    s.say("公平的一維醉漢一定回家", wait=2.0)
    s.say("但「還沒回家」只像 ", r"$1/\sqrt N$", " 那樣慢慢消失", wait=2.5)
    s.say("有偏向，就有一定比例的人一去不回", wait=3.0)


class W303OneDimension(ChuScene):
    def construct(self):
        self.show_header("W303　一維：一定回家")
        check_math()
        part_question(self)
        part_gf(self)
        part_asym(self)
        part_plot(self)
        part_biased(self)
        part_close(self)
        self.wait(1.0)
