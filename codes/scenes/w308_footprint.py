"""
W308 足跡大小＋第 3 章回顧（§3.8，pp. 65–68，(3.50)–(3.58)，Ex. 3.9–3.12）

問：走 n 步，足跡大約多大？和 n 成正比嗎？
  1. 一維：Δ_n ≈ √(2/(πn)) → δ_n ≈ √(8n/π)（(3.54)）；n = 100：精確 15.997（R.distinct_sites_1d）、近似 15.958、
     200 人模擬平均 15.665（R.SIM_DISTINCT，重算 assert）；圖：精確 δ_n 的點＋√(8n/π)。
  2. 二維：Δ_n ≈ π/ln n → δ_n ≈ πn/ln n；400 步足跡（R.trail、L.trail，踩過的格點用 L.heat_cells 依次數上色）：
     159 個點、終點 (−8, 14)、平均 δ_400 = 162.57；表：精確／πn/ln n（書）／πn/ln(8n)（Montroll–Weiss）；筆誤框 (3.56)。
  3. 三維 (3.58)：δ_n ≈ n/G(z_c;0) = 0.6595n；δ_n/n = 0.9051、0.7229、0.6782；練習 3.11 次一項 2A√n = 0.5739√n
     （A = R.sc_correction_coeff()，把 W305 的 A/√n' 加總）；練習 3.9、3.10、3.12 一行。
  4. 第 3 章回顧表＋同一個生成函數 (z_c/(z_c − z))/G(z;0)。
畫面上的數字都由 utils.walk_rec 算出，check_math() 核對。
"""
import math
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rec as R
from utils import lattice_draw as L
from utils import walk_draw as D

DIM_COLORS = {1: BLUE_C, 2: GREEN_C, 3: ORANGE}
NS2 = (10, 100, 1000, 10000)
CELL = 0.22
TRAIL_ORG = np.array([-3.53, -1.21, 0])


def fmt2(v):
    return f"{v:.1f}" if v >= 1000 else f"{v:#.4g}"


@lru_cache(maxsize=None)
def numbers():
    nb = {"d1": R.distinct_sites("1d", 100), "trail": R.trail("sq", R.TRAIL_N, R.TRAIL_SEED)}
    nb["d1_exact"] = f"{float(R.distinct_sites_1d(100)):.3f}"
    nb["d1_asym"] = f"{R.footprint_1d_asym(100):.3f}"
    nb["d1_sim"] = f"{R.SIM_DISTINCT[0][2] / 200:.3f}"
    dsq = R.distinct_sites("sq", 10000)
    nb["d400"] = f"{dsq[400]:.2f}"
    nb["sq_exact"] = [fmt2(dsq[n]) for n in NS2]
    nb["sq_book"] = [fmt2(R.footprint_2d_book(n)) for n in NS2]
    nb["sq_mw"] = [fmt2(R.footprint_2d_mw(n)) for n in NS2]
    dsc = R.distinct_sites("sc", 1000)
    nb["sc_ratio"] = [f"{dsc[n] / n:.4f}" for n in (10, 100, 1000)]
    A2 = 2 * R.sc_correction_coeff()
    lead, corr = 1000 / R.G_C["sc"], A2 * math.sqrt(1000)
    nb["A2"] = f"{A2:.4f}"
    nb["sc1000"] = (f"{lead:.1f}", f"{corr:.1f}", f"{lead + corr:.1f}", f"{dsc[1000]:.1f}")
    nb["inv"] = f"{1 / R.G_C['sc']:.4f}"
    nb["bcc"] = f"{R.escape('bcc'):.4f}"
    nb["biased"] = f"{sum(R.escape_biased_series(0.6, n) for n in range(101)):.2f}"
    nb["Gc3"] = f"{R.G_C['sc']:.3f}"
    return nb


def check_math():
    nb = numbers()
    assert (nb["d1_exact"], nb["d1_asym"], nb["d1_sim"]) == ("15.997", "15.958", "15.665")
    assert abs(math.sqrt(800 / math.pi) - R.footprint_1d_asym(100)) < 1e-12
    assert abs(nb["d1"][100] - float(R.distinct_sites_1d(100))) < 1e-9
    sims = R.simulate_distinct("1d", 200, 100, R.SIM_DISTINCT_SEED)
    assert len(sims) == 200 and sum(sims) == R.SIM_DISTINCT[0][2] == 3133 and R.SIM_DISTINCT[0][:2] == ("1d", 100)
    tr = nb["trail"]
    assert len(tr) == R.TRAIL_N + 1 == 401 and tr[0] == (0, 0) and tr[-1] == R.TRAIL_END == (-8, 14)
    assert len(set(tr)) == R.TRAIL_SITES == 159
    assert all(abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1 for a, b in zip(tr, tr[1:]))
    xs, ys = [p[0] for p in tr], [p[1] for p in tr]
    assert (min(xs), max(xs), min(ys), max(ys)) == (-13, 6, -6, 17)          # 版面（CELL、TRAIL_ORG）依此安排
    assert nb["d400"] == "162.57"
    assert nb["sq_exact"] == ["7.881", "49.50", "363.6", "2874.6"]
    assert nb["sq_book"] == ["13.64", "68.22", "454.8", "3410.9"]
    assert nb["sq_mw"] == ["7.169", "47.00", "349.6", "2782.7"]
    assert nb["sc_ratio"] == ["0.9051", "0.7229", "0.6782"]
    assert nb["A2"] == "0.5739" and nb["sc1000"] == ("659.5", "18.1", "677.6", "678.2")
    assert nb["inv"] == "0.6595" and nb["bcc"] == "0.7178" and nb["biased"] == "24.95" and nb["Gc3"] == "1.516"
    assert abs(R.escape_biased(0.6) - 0.2) < 1e-12 and f"{1 - nb['d1'][100] / 100:.1f}" == "0.8"
    d1 = nb["d1"]
    assert all(abs(d1[n] - R.footprint_1d_asym(n)) < 0.5 for n in range(4, 101, 4))   # 點貼著綠線
    assert max(d1) < 18 and R.footprint_1d_asym(100) < 18                             # 圖的縱軸到 18


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


def trail_cells(tr):
    """踩過的格點：踩 1、2、3、≥4 次 → 色階 2/5、3/5、4/5、1（L.heat_cells 的橘色）。"""
    cnt = Counter(tr)
    return L.heat_cells({p: 1 + min(c, 4) for p, c in cnt.items()}, cell=CELL, origin=TRAIL_ORG, vmax=5)


# ---------------------------------------------------------------- 各段
def part_question(s):
    law = at(tx(r"\delta_n=\sum_{n'=0}^{n}\Delta_{n'}=\sum_{n'=0}^{n}P_{n'}", 44), 0, 1.9)
    cols = VGroup()
    for d, x, name, p in ((1, -4.4, "一維", r"P_n\approx\sqrt{2/\pi n}"), (2, 0, "二維", r"P_n\approx\pi/\ln n"),
                          (3, 4.4, "三維 sc", r"P_n\to" + numbers()["inv"])):
        cols.add(at(VGroup(zh(name, font_size=30, color=DIM_COLORS[d]), tx(p, 34, DIM_COLORS[d]))
                    .arrange(DOWN, buff=0.3), x, -0.4))
    assert_in_frame(law, cols)
    assert_apart(law, cols)
    for a, b in zip(cols, cols[1:]):
        assert_apart(a, b)
    s.say("問：走 ", "$n$", " 步，足跡大約多大？和 ", "$n$", " 成正比嗎？", wait=0.5)
    s.play(Write(law), run_time=1.5)
    s.say("W307：足跡就是把還沒回家的機率 ", "$P_{n'}$", " 加起來", wait=3.0)
    s.play(LaggedStart(*[FadeIn(c, shift=0.2 * UP) for c in cols], lag_ratio=0.3), run_time=1.5)
    s.say("所以只要知道 ", "$n$", " 很大時 ", "$P_n$", " 怎麼變小", wait=3.0)
    s.clear_stage()


def part_1d(s):
    nb = numbers()
    f1 = tx(r"\Delta_n=P_n\approx\sqrt{\frac{2}{\pi n}}", 34)
    f2 = eq_label(tx(r"\delta_n\approx\sqrt{\frac{2}{\pi}}\sum_{n'\le n}n'^{-1/2}\approx\sqrt{\frac{8n}{\pi}}", 34),
                  "(3.54)")
    at(fit(VGroup(f1, f2).arrange(RIGHT, buff=0.9), 13.2), 0, 2.45)
    lines = VGroup(mixed("$n=100$", "：精確 ", "$" + nb["d1_exact"] + "$", font_size=28),
                   mixed("近似 ", "$\\sqrt{800/\\pi}=" + nb["d1_asym"] + "$", font_size=28),
                   mixed("200 人模擬：平均 ", "$3133/200=" + nb["d1_sim"] + "$", font_size=28))
    lines.arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    lines.move_to(np.array([-6.5 + lines.get_width() / 2, -0.6, 0]))
    ax = D.plot_axes((0, 105, 20), (0, 18, 5), width=5.6, height=3.8, x_ticks=[20, 40, 60, 80, 100],
                     y_ticks=[5, 10, 15], font_size=20, x_label="n", y_label=r"\delta_n")
    ax.move_to(np.array([3.7, -0.55, 0]))
    curve = VMobject().set_points_smoothly([ax.c2p(n, R.footprint_1d_asym(n)) for n in np.linspace(0.01, 100, 120)])
    curve.set_stroke(GREEN_C, 4)
    dots = VGroup(*[Dot(ax.c2p(n, nb["d1"][n]), radius=0.05).set_fill(WHITE, 1) for n in range(0, 101, 4)])
    last = Dot(ax.c2p(100, nb["d1"][100]), radius=0.09).set_fill(YELLOW, 1)
    lab = tx(nb["d1_exact"], 24, YELLOW).next_to(last, UP, buff=0.15)
    leg = VGroup(VGroup(Dot(radius=0.06).set_fill(WHITE, 1), mixed("精確 ", "$\\delta_n$", font_size=22))
                 .arrange(RIGHT, buff=0.15),
                 VGroup(Line(ORIGIN, 0.4 * RIGHT).set_stroke(GREEN_C, 4), tx(r"\sqrt{8n/\pi}", 26, GREEN_C))
                 .arrange(RIGHT, buff=0.15)).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
    leg.move_to(ax.c2p(22, 15.2))
    assert_in_frame(f1, f2, lines, ax, lab, leg)
    for a, b in ((f1, f2), (f2, lines), (f2, ax), (lines, ax), (f1, ax), (leg, lab)):
        assert_apart(a, b)
    assert leg.get_bottom()[1] > ax.c2p(40, R.footprint_1d_asym(40))[1] + 0.1     # 圖例在曲線上方
    s.play(FadeIn(f1), run_time=0.8)
    s.say("一維：", "$\\Delta_n=P_n\\approx\\sqrt{2/\\pi n}$", "（W303）", wait=3.0)
    s.play(FadeIn(f2, shift=0.2 * DOWN), run_time=1.0)
    s.say("加起來：", "$\\sum_{n'\\le n}n'^{-1/2}\\approx2\\sqrt n$", "，所以 ", "$\\delta_n\\approx\\sqrt{8n/\\pi}$",
          wait=3.5)
    s.play(FadeIn(ax), FadeIn(leg), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in dots], lag_ratio=0.08), ShowCreation(curve), run_time=2.5)
    s.say("白點：精確的 ", "$\\delta_n$", "；綠線：", "$\\sqrt{8n/\\pi}$", wait=3.0)
    s.play(FadeIn(last, scale=0.5), FadeIn(lab), FadeIn(lines[0]), FadeIn(lines[1]), run_time=1.0)
    s.say("$n=100$", "：精確 " + nb["d1_exact"] + "，近似 " + nb["d1_asym"], wait=3.0)
    s.play(FadeIn(lines[2]), run_time=0.8)
    s.say("200 人各走 100 步（模擬）：平均踩過 " + nb["d1_sim"] + " 個點", wait=3.0)
    s.say("一維的足跡只像 ", "$\\sqrt n$", " 那樣長：一直在走回頭路", wait=3.0)
    s.clear_stage()


def part_2d(s):
    nb = numbers()
    tr = nb["trail"]
    cells = trail_cells(tr)
    path = L.trail(tr, cell=CELL, origin=TRAIL_ORG, color=WHITE, stroke_width=1.3, end_radius=0.09)
    pic = VGroup(cells, path)
    r0 = mixed("400 步、終點 ", "$(-8,14)$", "：踩過 ", "$159$", " 個點", font_size=26)
    r0b = mixed("400 步的平均：", "$\\delta_{400}=" + nb["d400"] + "$", font_size=26)
    r1 = eq_label(tx(r"\delta_n\approx\sum_{n'=2}^{n}\frac{\pi}{\ln n'}\approx\frac{\pi n}{\ln n}", 32),
                  "(3.55)--(3.57)")
    box = D.erratum_box("p. 66，(3.56)", r"\sum_{n'=2}^{\infty}", r"\sum_{n'=2}^{n}", width=4.6)
    hdr = ["n", zh("精確", font_size=24, color=YELLOW), r"\pi n/\ln n", r"\pi n/\ln(8n)"]
    rows = [hdr] + [[t, a, b, c] for t, a, b, c in zip(["10", "100", "1000", "10^4"], nb["sq_exact"],
                                                         nb["sq_book"], nb["sq_mw"])]
    tab = D.tex_table(rows, [1.3, 1.8, 2.1, 2.3], row_h=0.52, font_size=28)
    note = zh("第三欄：書上的首項；第四欄：加上常數項（Montroll 與 Weiss）", font_size=20, color=GREY_A)
    right = fit(VGroup(r0, r0b, r1).arrange(DOWN, aligned_edge=LEFT, buff=0.3), 8.0)
    right.move_to(np.array([6.6 - right.get_width() / 2, 2.95 - right.get_height() / 2, 0]))
    at(box, 2.6, -0.6)
    tab.move_to(np.array([2.6, 0.2 - tab.get_height() / 2, 0]))
    note.next_to(tab, DOWN, buff=0.15)
    assert_in_frame(pic, right, box, tab, note)
    assert_table(tab, 4)
    for a, b in ((pic, right), (pic, box), (pic, tab), (pic, note), (right, box), (right, tab), (tab, note)):
        assert_apart(a, b)
    s.say("二維：", "$\\Delta_n=P_n\\approx\\pi/\\ln n$", "（W304）", wait=0.5)
    s.play(LaggedStart(*[FadeIn(c) for c in cells], lag_ratio=0.01), run_time=2.5)
    s.play(ShowCreation(path[0]), FadeIn(path[1]), run_time=3.0)
    s.play(FadeIn(path[2], scale=0.5), run_time=0.4)
    s.say("一條 400 步的足跡：格子越亮，踩過的次數越多", wait=3.0)
    s.play(FadeIn(r0), run_time=0.8)
    s.say("這一條踩過 159 個不同的點", wait=2.5)
    s.play(FadeIn(r0b), run_time=0.8)
    s.say("平均是 " + nb["d400"] + "：400 步裡大約四成是新點", wait=3.0)
    s.play(FadeIn(r1, shift=0.2 * DOWN), run_time=1.0)
    s.say("加起來：", "$\\delta_n\\approx\\pi n/\\ln n$", "：只比 ", "$n$", " 小一個 ", "$\\ln n$", wait=3.5)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (3.56) 的和寫到 ", "$\\infty$", "，應該只加到 ", "$n$", wait=3.5)
    s.play(FadeOut(box), run_time=0.5)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), FadeIn(note), run_time=1.2)
    s.say("首項 ", "$\\pi n/\\ln n$", " 偏大很多；", "$\\pi n/\\ln(8n)$", " 準得多", wait=3.5)
    s.say("$\\ln n$", " 長得很慢：二維的足跡幾乎和 ", "$n$", " 成正比", wait=3.0)
    s.clear_stage()


def part_3d(s):
    nb = numbers()
    f1 = at(eq_label(tx(r"\delta_n\approx\frac{n}{G(z_c;\vec 0)}=" + nb["inv"] + r"\,n", 38), "(3.58)"), 0, 2.45)
    rows = [["n", "10", "100", "1000"], [r"\delta_n/n\ (\mathrm{sc})"] + nb["sc_ratio"]]
    tab = D.tex_table(rows, [2.8, 1.6, 1.6, 1.6], row_h=0.6, font_size=30)
    at(tab, 0, 1.0)
    f2 = at(mixed("練習 3.11：", "$\\delta_n\\approx" + nb["inv"] + "\\,n+" + nb["A2"] + "\\sqrt n$", font_size=28),
            0, -0.25)
    a, b, c, e = nb["sc1000"]
    f3 = at(mixed("$n=1000$", "：", f"${a}+{b}={c}$", "，精確 ", f"${e}$", font_size=28), 0, -1.1)
    f4 = mixed("練習 3.9 bcc：", "$" + nb["bcc"] + "\\,n$", "；3.10 偏硬幣：", "$|p-q|\\,n$", "（", "$p=0.6$",
               "：", "$0.2\\,n$", "，", "$n=100$", " 精確 ", "$" + nb["biased"] + "$", "）；3.12 sc：",
               "$" + nb["inv"] + "\\,n$", font_size=22)
    at(fit(f4, 13.2), 0, -2.3)
    assert_in_frame(f1, tab, f2, f3, f4)
    assert_table(tab, 4)
    for x, y in ((f1, tab), (tab, f2), (f2, f3), (f3, f4)):
        assert_apart(x, y)
    s.play(FadeIn(f1), run_time=0.8)
    s.say("三維：", "$P_n\\to1/G(z_c;\\vec 0)$", "，是個常數", wait=3.0)
    s.say("所以足跡和 ", "$n$", " 成正比：每一步有 66% 的機會踩到新點", wait=3.5)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.0)
    s.say("$\\delta_n/n$", " 慢慢降到 " + nb["inv"] + "：差的是 ", "$\\sqrt n$", " 的修正", wait=3.5)
    s.play(FadeIn(f2), run_time=0.8)
    s.say("練習 3.11：把 W305 的 ", "$0.28696/\\sqrt{n'}$", " 加起來", wait=3.0)
    s.play(FadeIn(f3), run_time=0.8)
    s.say("$n=1000$", "：", f"${c}$", " 對 ", f"${e}$", "，差不到 1", wait=3.0)
    s.play(FadeIn(f4), run_time=0.8)
    s.say("其他練習：bcc、偏硬幣都是 ", "$n\\times$", "逃脫機率", wait=3.5)
    s.clear_stage()


def part_recap(s):
    nb = numbers()
    title = at(zh("第 3 章回顧", font_size=32, color=YELLOW), 0, 2.75)
    hdr = [zh("維度", font_size=26, color=YELLOW), mixed("還沒回家 ", "$P_N$", font_size=26, color=YELLOW),
           mixed("足跡 ", "$\\delta_n$", font_size=26, color=YELLOW), r"G(z_c;\vec 0)"]
    rows = [hdr, [zh("一維", font_size=28, color=DIM_COLORS[1]), r"\sqrt{2/\pi N}", r"\sqrt{8n/\pi}", r"\infty"],
            [zh("二維", font_size=28, color=DIM_COLORS[2]), r"\pi/\ln N", r"\pi n/\ln n", r"\infty"],
            [zh("三維 sc", font_size=28, color=DIM_COLORS[3]), nb["inv"], nb["inv"] + r"\,n", nb["Gc3"]]]
    tab = D.tex_table(rows, [2.2, 3.2, 3.2, 2.8], row_h=0.8, font_size=34)
    tab.move_to(np.array([0, 2.3 - tab.get_height() / 2, 0]))
    gf = at(mixed("同一個生成函數 ", "$\\frac{z_c}{z_c-z}\\,\\frac{1}{G(z;\\vec 0)}$", " 管回家也管足跡",
                  font_size=30), 0, -1.75)
    assert_in_frame(title, tab, gf)
    assert_table(tab, 4)
    assert_apart(title, tab)
    assert_apart(tab, gf)
    s.play(FadeIn(title), FadeIn(VGroup(*[D.table_cell(tab, 0, j, 4) for j in range(4)])), ShowCreation(tab[1]),
           run_time=1.0)
    s.say("第 3 章回顧", wait=0.3)
    s.play(LaggedStart(*[FadeIn(VGroup(*[D.table_cell(tab, i, j, 4) for j in range(4)])) for i in (1, 2, 3)],
                       lag_ratio=0.5), run_time=2.0)
    s.say("一、二維一定回家；三維有 66% 的人一去不返", wait=3.5)
    s.say("足跡：", "$\\sqrt n$", "、", "$n/\\ln n$", "、", "$n$", "：維度越高，越少走回頭路", wait=3.5)
    s.play(FadeIn(gf, shift=0.2 * UP), run_time=1.0)
    s.say("全看 ", "$G(z;\\vec 0)$", " 在 ", "$z_c$", " 是否發散", wait=3.0)


def part_close(s):
    s.say("回家與足跡是同一枚硬幣的兩面", wait=2.5)
    s.say("一個 ", "$G(z;\\vec 0)$", " 全部說完", wait=3.0)


class W308Footprint(ChuScene):
    def construct(self):
        self.show_header("W308　足跡大小與第 3 章回顧")
        check_math()
        part_question(self)
        part_1d(self)
        part_2d(self)
        part_3d(self)
        part_recap(self)
        part_close(self)
        self.wait(1.0)
