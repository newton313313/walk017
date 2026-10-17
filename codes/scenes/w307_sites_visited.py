"""
W307 走過幾個不同的點？（§3.8，pp. 63–65，(3.42)–(3.49)，Fig. 3.3）

問：醉漢走 n 步，平均踩過幾個不同的地方？
  1. (3.42) δ_n = Σ_{n'≤n} Δ_{n'}：一維 4 步的 16 條路線（R.walks_1d(4)，新點畫綠點、左上角寫足跡大小）
     Δ_0..Δ_4 = 1, 1, 1/2, 1/2, 3/8（R.new_site_probs_brute）；δ_4 = 54/16 = 27/8（R.distinct_sites_brute）。
  2. Δ_n = P_n：W303 的 1, 2, 2, 4, 6 除以 2^n；倒著走的動畫（0→1→0→1→2 倒過來是 2→1→0→1→0）。
  3. 書上的算法（Fig. 3.3 的示意圖畫在本檔 fig33）：G_BA/(1 + G_1) 對 B 求和 = (z_c/(z_c − z))/G(z;0)；Δ_10 = 252/1024。
  4. 二維核對：暴力列舉 4^n 條（n ≤ 5）＝ R.distinct_sites("sq", 5)。
畫面上的數字都由 utils.walk_rec 算出，check_math() 核對。
"""
import math
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import accumulate
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

NEW_COLOR = "#7CE38B"
EXAMPLE = [0, 1, 0, 1, 2]
CELL_W, CELL_H, UX, UY = 1.62, 1.22, 0.3, 0.12


def is_new(pos, k):
    return pos[k] not in pos[:k]


def ftex(fr):
    fr = Fraction(fr)
    return str(fr.numerator) if fr.denominator == 1 else r"\tfrac{%d}{%d}" % (fr.numerator, fr.denominator)


def gfmt(v):
    return f"{v:.4f}".rstrip("0").rstrip(".")


@lru_cache(maxsize=None)
def numbers():
    walks = R.walks_1d(4)
    nb = {"walks": walks, "cnt": [sum(is_new(w, k) for w in walks) for k in range(5)]}
    nb["delta"] = [Fraction(c, 16) for c in nb["cnt"]]
    nb["dl"] = list(accumulate(nb["delta"]))
    nb["nr"] = [R.never_return_1d(n) for n in range(5)]
    nb["sq"] = [R.distinct_sites_brute("sq", n) for n in range(6)]
    nb["gf"] = [gfmt(v) for v in R.distinct_sites("sq", 5)]
    return nb


def check_math():
    nb = numbers()
    walks = nb["walks"]
    assert len(walks) == 16 and EXAMPLE in walks and nb["cnt"] == [16, 16, 8, 8, 6]
    assert nb["delta"] == R.new_site_probs_brute("1d", 4) == [1, 1, Fraction(1, 2), Fraction(1, 2), Fraction(3, 8)]
    assert sum(len(set(w)) for w in walks) == 54 and Fraction(54, 16) == R.distinct_sites_brute("1d", 4)
    assert nb["dl"] == [1, 2, Fraction(5, 2), 3, Fraction(27, 8)]
    assert nb["nr"] == [1, 2, 2, 4, 6] and [Fraction(c, 2 ** n) for n, c in enumerate(nb["nr"])] == nb["delta"]
    for w in walks:                      # 第 k 步是新點 ⇔ 倒著走 k 步都沒回到出發點
        for k in range(1, 5):
            rev = [w[k - j] for j in range(k + 1)]
            assert is_new(w, k) == all(x != rev[0] for x in rev[1:])
    assert is_new(EXAMPLE, 4) and EXAMPLE[::-1] == [2, 1, 0, 1, 0] and len(set(EXAMPLE)) == 3
    d10 = R.new_site_probs_brute("1d", 10)
    assert [d * 2 ** n for n, d in enumerate(d10)] == R.never_return_1d_gf(10)       # 暴力列舉＝g/G 的係數
    assert d10[10] == Fraction(252, 1024) == R.escape_1d(10) and 252 == math.comb(10, 5) == R.never_return_1d(10)
    assert f"{float(d10[10]):.4f}" == "0.2461"
    assert nb["sq"] == [1, 2, Fraction(11, 4), Fraction(7, 2), Fraction(267, 64), Fraction(155, 32)]
    assert all(abs(float(a) - b) < 1e-12 for a, b in zip(nb["sq"], R.distinct_sites("sq", 5)))
    assert nb["gf"] == ["1", "2", "2.75", "3.5", "4.1719", "4.8438"]


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


def grid_origin(i):
    r, c = divmod(i, 4)
    return np.array([-6.45 + c * CELL_W, 1.7 - r * CELL_H, 0])


def mini_walk(pos, org):
    """一條路線的小走法圖：VGroup(框, 基線, 折線, 腳印們, 足跡大小)；新點是綠色大點。"""
    frame = RoundedRectangle(width=1.5, height=1.12, corner_radius=0.08).set_stroke(GREY_D, 1)
    frame.move_to(org + 2 * UX * RIGHT)
    base = DashedLine(org, org + 4 * UX * RIGHT, dash_length=0.04).set_stroke(GREY_D, 1.5)
    path = VMobject().set_points_as_corners([D.path_point(k, x, UX, UY, org) for k, x in enumerate(pos)])
    path.set_stroke(GREY_B, 2.5)
    dots = VGroup(*[Dot(D.path_point(k, x, UX, UY, org), radius=0.06 if is_new(pos, k) else 0.035)
                    .set_fill(NEW_COLOR if is_new(pos, k) else GREY_C, 1) for k, x in enumerate(pos)])
    foot = tx(str(len(set(pos))), 22, YELLOW).move_to(frame.get_corner(UL) + np.array([0.15, -0.15, 0]))
    return VGroup(frame, base, path, dots, foot)


def fig33():
    """Fig. 3.3 的示意：從 A 走到 B（G_BA），再從 B 繞回 B（G_1）。"""
    A, B = np.array([-5.7, -1.3, 0]), np.array([-2.8, 0.7, 0])
    pts = [A, (-5.3, 0.1), (-4.5, -0.9), (-3.8, 0.5), (-3.3, -0.4), B]
    leg = VMobject().set_points_smoothly([np.array([*p[:2], 0.0]) for p in pts]).set_stroke(YELLOW, 4)
    lp = [B, (-2.3, 2.0), (-1.0, 1.8), (-1.2, 0.5), B]
    loop = VMobject().set_points_smoothly([np.array([*p[:2], 0.0]) for p in lp]).set_stroke(ORANGE, 4)
    dA = Dot(A, radius=0.11).set_fill(L.START_COLOR, 1)
    dB = Dot(B, radius=0.11).set_fill(RED, 1)
    lA = tx("A", 32).next_to(dA, DOWN, buff=0.15)
    lB = tx("B", 32).next_to(dB, LEFT, buff=0.2)
    gl = at(tx("G_{BA}", 32, YELLOW), -4.3, -1.65)
    g1 = at(tx(r"G_1(z;\vec 0)", 30, ORANGE), -0.9, 2.45)
    return VGroup(leg, loop, dA, dB, lA, lB, gl, g1)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nl = at(D.number_line(4, unit=0.9, label_every=1, font_size=26), 0, 0.5)
    walker = D.walker_dot().move_to(D.line_point(nl, 0))
    count = at(mixed("踩過的不同點：", "$1$", font_size=30), 0, -1.2)
    marks = VGroup(Square(side_length=0.36).set_fill(NEW_COLOR, 0.45).set_stroke(NEW_COLOR, 2)
                   .move_to(D.line_point(nl, 0) + 0.55 * UP))
    assert_in_frame(nl, count)
    assert_apart(nl, count)
    s.say("問：醉漢走 ", "$n$", " 步，平均踩過幾個不同的地方？", wait=0.5)
    s.play(FadeIn(nl), FadeIn(walker), FadeIn(marks[0]), FadeIn(count), run_time=1.0)
    seen = {0}
    for k in range(1, 5):
        x = EXAMPLE[k]
        anims = [walker.animate.move_to(D.line_point(nl, x))]
        if x not in seen:
            seen.add(x)
            sq = marks[0].copy().move_to(D.line_point(nl, x) + 0.55 * UP)
            marks.add(sq)
            new_num = tx(str(len(seen)), 30 * 1.3).move_to(count[1])
            anims += [FadeIn(sq, scale=0.5), Transform(count[1], new_num)]
        s.play(*anims, run_time=0.7)
        s.wait(0.3)
    s.say("例：", "$0\\to1\\to0\\to1\\to2$", "：走了 4 步，踩過 3 個不同的點", wait=3.0)
    s.say("定義 ", "$\\delta_n$", "＝", "$n$", " 步平均踩過的不同點數（足跡大小）", wait=3.0)
    s.clear_stage()


def part_example(s):
    nb = numbers()
    walks = VGroup(*[mini_walk(w, grid_origin(i)) for i, w in enumerate(nb["walks"])])
    title = at(zh("一維 4 步的 16 條路線（綠點＝踩到新點）", font_size=24), -3.4, 2.7)
    law = at(eq_label(tx(r"\delta_n=\sum_{n'=0}^{n}\Delta_{n'}", 30), "(3.42)"), 3.55, 2.7)
    hdr = ["n", zh("新點條數", font_size=24, color=YELLOW), r"\Delta_n", r"\delta_n", "P_n"]
    rows = [hdr] + [[str(n), str(nb["cnt"][n]), ftex(nb["delta"][n]), ftex(nb["dl"][n]),
                     r"\tfrac{%d}{%d}" % (nb["nr"][n], 2 ** n)] for n in range(5)]
    tab = D.tex_table(rows, [0.8, 1.9, 1.1, 1.1, 1.2], row_h=0.58, font_size=30)
    tab.move_to(np.array([0.5 + tab.get_width() / 2, 2.1 - tab.get_height() / 2, 0]))
    foot = mixed("16 條的足跡共 54 個點：", "$\\delta_4=\\frac{54}{16}=\\frac{27}{8}$", font_size=24)
    at(foot, 3.55, -2.05)
    assert_in_frame(walks, title, law, tab, foot)
    assert_table(tab, 5)
    for a, b in ((title, walks), (law, tab), (tab, foot), (walks, tab), (walks, foot), (title, law)):
        assert_apart(a, b)
    s.play(FadeIn(title), FadeIn(law), run_time=0.8)
    s.say("$\\Delta_n$", "＝第 ", "$n$", " 步踩到新點的機率；", "$\\delta_n$", " 是它們的和", wait=1.0)
    s.play(LaggedStart(*[FadeIn(VGroup(m[0], m[1], m[2])) for m in walks], lag_ratio=0.05), run_time=2.0)
    s.play(FadeIn(VGroup(*[D.table_cell(tab, 0, j, 5) for j in range(4)])), ShowCreation(tab[1]), run_time=0.8)
    says = [("第 0 步：出發點本身算新點，", "$\\Delta_0=1$"),
            ("第 1 步：16 條都踩到新點，", "$\\Delta_1=1$"),
            ("第 2 步：一半走回 0，", "$\\Delta_2=8/16=\\tfrac12$"),
            ("第 3 步：還是 8 條踩到新點，", "$\\Delta_3=\\tfrac12$"),
            ("第 4 步：只剩 6 條，", "$\\Delta_4=6/16=\\tfrac38$")]
    for n in range(5):
        cells = VGroup(*[D.table_cell(tab, n + 1, j, 5) for j in range(4)])
        s.play(LaggedStart(*[FadeIn(m[3][n], scale=0.5) for m in walks], lag_ratio=0.03), FadeIn(cells),
               run_time=1.2)
        s.say(*says[n], wait=2.2)
    s.play(LaggedStart(*[FadeIn(m[4]) for m in walks], lag_ratio=0.04), run_time=1.2)
    s.say("左上角：每條路線踩過幾個不同的點", wait=2.5)
    s.play(FadeIn(foot), run_time=0.8)
    s.say("加起來 54，平均 ", "$27/8$", "：正好是 ", "$\\Delta_0+\\cdots+\\Delta_4$", wait=3.0)
    pcol = VGroup(*[D.table_cell(tab, i, 4, 5) for i in range(6)])
    s.play(FadeIn(pcol), run_time=1.0)
    s.say("最後一欄：W303 的「還沒回家」1、2、2、4、6 除以 ", "$2^n$", wait=3.5)
    b1 = SurroundingRectangle(VGroup(*[D.table_cell(tab, i, 2, 5) for i in range(1, 6)]), buff=0.08)
    b2 = SurroundingRectangle(VGroup(*[D.table_cell(tab, i, 4, 5) for i in range(1, 6)]), buff=0.08)
    VGroup(b1, b2).set_stroke(YELLOW, 2.5)
    s.play(ShowCreation(b1), ShowCreation(b2), run_time=0.8)
    s.say("兩欄一模一樣：", "$\\Delta_n=P_n$", "！為什麼？", wait=3.0)
    s.play(FadeOut(VGroup(walks, title, foot, law)), run_time=0.8)
    return tab, VGroup(b1, b2)


def part_reverse(s, keep):
    ux, uy = 0.85, 0.38
    fo, ro = np.array([-6.0, 1.3, 0]), np.array([-6.0, -1.85, 0])
    fax = D.walk_axes(4, 2, ux, uy, fo, x_every=1, y_every=1, font_size=20)
    rax = D.walk_axes(4, 2, ux, uy, ro, x_every=1, y_every=1, font_size=20)
    fwd = D.walk_path(EXAMPLE, ux, uy, fo, color=YELLOW, stroke_width=4, dot_radius=0.06)
    rev = D.walk_path(EXAMPLE[::-1], ux, uy, ro, color=TEAL, stroke_width=4, dot_radius=0.06)
    fline = DashedLine(D.path_point(0, 2, ux, uy, fo), D.path_point(4, 2, ux, uy, fo)).set_stroke(NEW_COLOR, 2.5)
    rline = DashedLine(D.path_point(0, 2, ux, uy, ro), D.path_point(4, 2, ux, uy, ro)).set_stroke(NEW_COLOR, 2.5)
    fnew = Dot(D.path_point(4, 2, ux, uy, fo), radius=0.11).set_fill(NEW_COLOR, 1)
    rstart = Dot(D.path_point(0, 2, ux, uy, ro), radius=0.11).set_fill(NEW_COLOR, 1)
    ftitle = mixed("原來：", "$0\\to1\\to0\\to1\\to2$", font_size=24).next_to(fax, UP, buff=0.12).align_to(fax, LEFT)
    rtitle = mixed("倒著走：", "$2\\to1\\to0\\to1\\to0$", font_size=24).next_to(rax, UP, buff=0.12)
    rtitle.align_to(rax, LEFT)
    assert_in_frame(fax, rax, ftitle, rtitle)
    for a, b in ((ftitle, fax), (fax, rtitle), (rtitle, rax), (fax, keep), (rax, keep), (ftitle, keep)):
        assert_apart(a, b)
    s.play(FadeIn(fax), FadeIn(ftitle), run_time=0.8)
    s.play(ShowCreation(fwd[0]), FadeIn(fwd[1]), run_time=1.5)
    s.play(ShowCreation(fline), FadeIn(fnew, scale=0.5), run_time=0.8)
    s.say("第 4 步踩到 2：之前從沒碰過綠色虛線", wait=3.0)
    s.play(FadeIn(rax), FadeIn(rtitle), run_time=0.8)
    s.play(TransformFromCopy(fwd[0], rev[0]), run_time=1.5)
    s.play(FadeIn(rev[1]), ShowCreation(rline), FadeIn(rstart, scale=0.5), run_time=0.8)
    s.say("把路線倒著走：從 2 出發，之後再也沒回到 2", wait=3.0)
    s.say("「第 ", "$n$", " 步踩到新點」＝「倒著走 ", "$n$", " 步都沒回到出發點」", wait=3.5)
    s.say("倒著走也是公平的隨機漫步：所以 ", "$\\Delta_n=P_n$", wait=3.0)
    s.clear_stage()


def part_book(s):
    fig = fig33()
    cap = at(zh("Fig. 3.3 的示意", font_size=22, color=GREY_B), -3.5, -2.45)
    f1 = at(fit(tx(r"G_{BA}-G_{BA}G_1+G_{BA}G_1^2-\cdots=\frac{G_{BA}}{1+G_1(z;\vec 0)}", 30), 6.6), 3.3, 2.45)
    f2 = at(fit(tx(r"\sum_B\frac{G_{BA}}{1+G_1}=\frac{g(z;\vec 0)}{1+G_1(z;\vec 0)}"
                   r"=\frac{z_c}{z_c-z}\,\frac{1}{G(z;\vec 0)}", 30), 6.6), 3.3, 1.1)
    f3 = at(mixed("$1+G_1=G$", "（3.47）；和 (3.16) 是同一個函數", font_size=24), 3.3, -0.05)
    f4 = at(mixed("一維：", "$\\Delta_{10}=\\frac{252}{1024}=0.2461$", font_size=28), 3.3, -1.3)
    assert_in_frame(fig, cap, f1, f2, f3, f4)
    for a, b in ((f1, f2), (f2, f3), (f3, f4), (fig, f1), (fig, f2), (fig, f3), (fig, f4), (fig, cap)):
        assert_apart(a, b)
    leg, loop, dA, dB, lA, lB, gl, g1 = fig
    s.play(FadeIn(dA), FadeIn(lA), FadeIn(dB), FadeIn(lB), FadeIn(cap), run_time=0.8)
    s.say("書上的算法：走到 B 的路線，扣掉之前到過 B 的", wait=0.3)
    s.play(ShowCreation(leg), FadeIn(gl), run_time=1.5)
    s.play(ShowCreation(loop), FadeIn(g1), run_time=1.5)
    s.say("$G_{BA}G_1$", "：先到 B、再繞回 B；排容一路加減下去", wait=3.0)
    s.play(Write(f1), run_time=1.5)
    s.say("等比級數：剩下「最後一步才第一次踩到 B」的路線", wait=3.5)
    s.play(FadeIn(f2, shift=0.2 * DOWN), run_time=1.0)
    s.say("對所有終點 B 加總：", "$\\sum_BG_{BA}=g(z;\\vec 0)=\\frac{z_c}{z_c-z}$", wait=3.5)
    s.play(FadeIn(f3), run_time=0.8)
    s.say("用 ", "$1+G_1=G$", "：正是 (3.16) 從不回家的生成函數", wait=3.5)
    s.play(FadeIn(f4), run_time=0.8)
    s.say("一維 ", "$n=10$", "：第 10 步踩到新點的機率 0.2461", wait=3.0)
    s.clear_stage()


def part_square(s):
    nb = numbers()
    rows = [[zh("步數 n", font_size=24, color=YELLOW)] + [str(n) for n in range(6)],
            [zh("路線數", font_size=24)] + [str(4 ** n) for n in range(6)],
            [zh("暴力列舉", font_size=24)] + [ftex(v) for v in nb["sq"]],
            [zh("生成函數", font_size=24)] + nb["gf"]]
    tab = D.tex_table(rows, [2.4] + [1.5] * 6, row_h=0.72, font_size=30)
    at(tab, 0, 0.2)
    title = mixed("正方格子：", "$n$", " 步平均踩過幾個點 ", "$\\delta_n$", font_size=28).next_to(tab, UP, buff=0.4)
    assert_in_frame(tab, title)
    assert_table(tab, 7)
    assert_apart(title, tab)
    s.play(FadeIn(title), FadeIn(VGroup(*[D.table_cell(tab, i, j, 7) for i in range(2) for j in range(7)])),
           ShowCreation(tab[1]), run_time=1.0)
    s.say("二維核對：把 ", "$4^n$", " 條路線全部列出來（", "$n\\le5$", "）", wait=2.5)
    s.play(FadeIn(VGroup(*[D.table_cell(tab, 2, j, 7) for j in range(7)])), run_time=1.0)
    s.say("數每條踩過幾個不同的點、取平均", wait=2.5)
    s.play(FadeIn(VGroup(*[D.table_cell(tab, 3, j, 7) for j in range(7)])), run_time=1.0)
    s.say("和生成函數算的完全一樣：", "$\\delta_5=155/32=4.8438$", wait=3.5)


def part_close(s):
    s.say("新點的機率＝還沒回家的機率", wait=2.5)
    s.say("回家的問題一解決，足跡的問題也跟著解決", wait=3.0)


class W307SitesVisited(ChuScene):
    def construct(self):
        self.show_header("W307　走過幾個不同的點？")
        check_math()
        part_question(self)
        tab, boxes = part_example(self)
        part_reverse(self, VGroup(tab, boxes))
        part_book(self)
        part_square(self)
        part_close(self)
        self.wait(1.0)
