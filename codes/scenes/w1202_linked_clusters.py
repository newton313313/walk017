"""
W1202 配對遊戲：只剩連通的圖（§12.2.4 後半、§12.7，pp. 289–292、317–322，(12.12)–(12.15)、(12:S-1)–(12:S-7)，
Fig. 12.6–12.11、12:S-2–12:S-5，Ex. 12.1、12.9–12.11）

問：殼層的變數積掉之後，四次頂點會變成什麼？
  1. Fig. 12.6：(a) 一隻手在殼層 4u/N、(b) 兩隻手在殼層 6u/N；Gaussian 平均＝殼層的手兩兩牽起來（(12:S-2)）。
  2. Fig. 12:S-2：3 個片段 15 種牽法 = 1 + 6 + 8；片段圖 1、3、1；f₃ = 6、f₄ = 8。
  3. Fig. 12:S-5：5 個片段 945 種牽法裡 20 種 = 10 張片段圖 × 2；10 張小圖；筆誤框 1（p. 320）。
  4. (12:S-4) 表 6 列；(12:S-3)(12:S-5) 指數化；筆誤框 2（(12:S-7)）。
  5. Ex. 12.1 連結叢集的數值核對：0.9129、二階 0.0150 = 0.0100 + 0.0050、κ_k、(1−2x)^{−1/2}、0.8087。
畫面上的數字都由 utils.walk_eps 算出，check_math() 核對。
模組層級的畫圖函式：vertex_icon、frag_dot、self_loop、double_line、mini_diagram、check_mark。
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
from utils import walk_eps as E
from utils import lattice_draw as L
from utils import walk_draw as D

SHELL_COLOR = TEAL
OK_COLOR = "#7CE38B"
PAIR_CASES = ((2, 1), (0, 2), (3, 1), (4, 1), (2, 2), (0, 3))
METHODS = (("quad", "數值積分"), ("series", "所有圖的和"), ("linked", "連通圖取指數"), ("closed", "閉式"))


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def sizes_key(n1, n2):
    return tuple([2] * n2 + [1] * n1)


def double_pair(m):
    """一個 5 片段配對裡，連成雙線的那一對片段（片段 k 的腳是 2k、2k+1）。"""
    cross = [tuple(sorted((a // 2, b // 2))) for a, b in m if a // 2 != b // 2]
    assert len(cross) == 2 and cross[0] == cross[1]
    return cross[0]


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m5 = [m for m in E.matchings(5) if tuple(E.ring_sizes(m, 5)) == (2, 1, 1, 1)]
        pairs = [double_pair(m) for m in m5]
        return {"c3": E.type_census(3), "d3": E.diagram_census(3), "f3": E.ring_aut(3), "f4": E.ring_aut(4),
                "n5": len(E.matchings(5)), "t5": E.type_census(5)[(2, 1, 1, 1)], "d5": E.diagram_census(5)[(2, 1, 1, 1)],
                "pf": {p: E.pair_factor(*p) for p in PAIR_CASES},
                "dc": {p: E.diagram_census(p[0] + 2 * p[1])[sizes_key(*p)] for p in PAIR_CASES},
                "pairs": sorted(set(pairs)), "pair_mult": [pairs.count(p) for p in sorted(set(pairs))],
                "avg": {k: E.lc_average(0.1, 1.0, k) for k, _ in METHODS},
                "ord2": E.lc_order(0.1, 1.0, 2), "kap": E.lc_cumulants_rings(0.1, 1.0, 4),
                "w": [E.ring_weight(k) for k in range(1, 5)], "egf": E.egf_check(0.1), "multi": E.lc_multi()}


def check_math():
    nb = numbers()
    assert nb["c3"] == {(1, 1, 1): 1, (2, 1): 6, (3,): 8} and sum(nb["c3"].values()) == 15
    assert nb["d3"] == {(1, 1, 1): 1, (2, 1): 3, (3,): 1}
    assert nb["f3"] == 6 == math.factorial(3) and nb["f4"] == 8 and math.factorial(4) == 24
    assert nb["n5"] == 945 and nb["t5"] == 20 and nb["d5"] == 10 and E.pair_factor(3, 1) == 10
    assert len(nb["pairs"]) == 10 and nb["pair_mult"] == [2] * 10
    assert math.factorial(5) // (math.factorial(3) * math.factorial(2)) == 10
    assert math.factorial(5) // (math.factorial(2) * math.factorial(1)) == 60
    assert nb["pf"] == {(2, 1): 6, (0, 2): 3, (3, 1): 10, (4, 1): 15, (2, 2): 45, (0, 3): 15}
    assert nb["dc"] == nb["pf"]
    for (n1, n2), v in nb["pf"].items():
        n = n1 + 2 * n2
        assert v == math.factorial(n) // (math.factorial(n1) * math.factorial(n2) * 2 ** n2)
    assert E.LC == (0.1, 1.0)
    assert all(f4(v) == "0.9129" for v in nb["avg"].values()) and f4((1 + 0.2) ** -0.5) == "0.9129"
    assert [f4(x) for x in nb["ord2"]] == ["0.0150", "0.0100", "0.0050"]
    assert [f4(x) for x in nb["kap"]] == ["-0.1000", "0.0200", "-0.0080", "0.0048"]
    assert nb["w"] == [1, 2, 8, 48]
    assert all(f4(x) == "1.1180" for x in nb["egf"]) and f4(nb["multi"]) == "0.8087"


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def check_mark(color=OK_COLOR, size=0.3):
    return VMobject().set_points_as_corners([size * np.array([-0.5, 0.0, 0]), size * np.array([-0.15, -0.4, 0]),
                                             size * np.array([0.55, 0.45, 0])]).set_stroke(color, 5)


# ---------------------------------------------------------------- 畫圖元件
def vertex_icon(k_shell, origin, arm=0.55):
    """四次頂點：點＋四隻手；右邊的 k 隻在殼層（青色）。"""
    o = np.array(origin, dtype=float)
    angs = [PI * 3 / 4, PI * 5 / 4, -PI / 4, PI / 4]
    g = VGroup()
    for j, a in enumerate(angs):
        col = SHELL_COLOR if j >= 4 - k_shell else WHITE
        g.add(Line(o, o + arm * np.array([math.cos(a), math.sin(a), 0])).set_stroke(col, 5))
    g.add(Dot(o, radius=0.09).set_color(YELLOW))
    return g


def frag_dot(p):
    return Dot(p, radius=0.08).set_color(YELLOW)


def self_loop(p, direction, r=0.14):
    d = np.array(direction, dtype=float)
    d = d / np.linalg.norm(d)
    return Circle(radius=r).move_to(np.array(p) + r * d).set_stroke(SHELL_COLOR, 3)


def double_line(p, q, angle=0.7):
    return VGroup(ArcBetweenPoints(p, q, angle=angle), ArcBetweenPoints(p, q, angle=-angle)).set_stroke(SHELL_COLOR, 3)


def ring_points(k, radius, origin):
    o = np.array(origin, dtype=float)
    return [o + radius * np.array([math.cos(PI / 2 + 2 * PI * j / k), math.sin(PI / 2 + 2 * PI * j / k), 0])
            for j in range(k)]


def mini_diagram(pair, origin, radius=0.42):
    """5 個片段排成五邊形：pair 那兩個連雙線，其餘自己接成圈（≤ 12 個物件）。"""
    o = np.array(origin, dtype=float)
    pts = ring_points(5, radius, o)
    g = VGroup()
    i, j = pair
    g.add(double_line(pts[i], pts[j], 0.45))
    for k, p in enumerate(pts):
        if k not in pair:
            g.add(self_loop(p, p - o))
    g.add(*[frag_dot(p) for p in pts])
    return g


def type_panel(kind, origin):
    """3 個片段：'111' 三個自環、'21' 一個自環＋雙線、'3' 三角環。"""
    o = np.array(origin, dtype=float)
    pts = ring_points(3, 0.5, o)
    g = VGroup()
    if kind == "111":
        g.add(*[self_loop(p, p - o) for p in pts])
    elif kind == "21":
        g.add(self_loop(pts[0], pts[0] - o), double_line(pts[1], pts[2], 0.6))
    else:
        g.add(*[Line(pts[a], pts[(a + 1) % 3]).set_stroke(SHELL_COLOR, 3) for a in range(3)])
    g.add(*[frag_dot(p) for p in pts])
    return g


# ---------------------------------------------------------------- 各段
def part_question(s):
    va = vertex_icon(1, ORIGIN)
    vb = vertex_icon(2, ORIGIN)
    la = VGroup(tx(r"\text{(a)}\ \frac{4u}N", 30), zh("一隻手在殼層", font_size=22, color=GREY_A)).arrange(DOWN, buff=0.15)
    lb = VGroup(tx(r"\text{(b)}\ \frac{6u}N", 30), zh("兩隻手在殼層", font_size=22, color=GREY_A)).arrange(DOWN, buff=0.15)
    ga = VGroup(va, la).arrange(DOWN, buff=0.3)
    gb = VGroup(vb, lb).arrange(DOWN, buff=0.3)
    figs = VGroup(ga, gb).arrange(RIGHT, buff=1.6)
    tag = zh("Fig. 12.6　(12.14)(12.15)", font_size=22, color=GREY_B)
    left = VGroup(tag, figs).arrange(DOWN, buff=0.35)
    place_in(left, -6.6, -0.6, -1.0, 3.0)
    prop = tx(r"\langle S(\vec Q_1)S(\vec Q_2)\rangle=\frac{\delta(\vec Q_1+\vec Q_2)}{r+Q_1^2}", 32)
    ptag = zh("(12:S-2)　殼層的 Gaussian 平均", font_size=22, color=GREY_B)
    p1, p2 = np.array([-0.6, 0, 0]), np.array([0.6, 0, 0])
    pair = VGroup(frag_dot(p1), frag_dot(p2), double_line(p1, p2, 0.8),
                  Line(p1, p1 + 0.5 * LEFT + 0.35 * UP).set_stroke(WHITE, 4),
                  Line(p1, p1 + 0.5 * LEFT + 0.35 * DOWN).set_stroke(WHITE, 4),
                  Line(p2, p2 + 0.5 * RIGHT + 0.35 * UP).set_stroke(WHITE, 4),
                  Line(p2, p2 + 0.5 * RIGHT + 0.35 * DOWN).set_stroke(WHITE, 4))
    right = VGroup(ptag, prop, pair).arrange(DOWN, buff=0.4)
    place_in(right, 0.0, 6.7, -1.0, 3.0)
    assert_in_frame(left, right)
    assert_apart(left, right)
    s.say("問：殼層的變數積掉之後，四次頂點會變成什麼？", wait=2.0)
    s.play(FadeIn(tag), FadeIn(ga), run_time=1.0)
    s.say("頂點畫成一個點＋四隻手；青色的手在殼層", wait=2.0)
    s.say("(a) 一隻手在殼層：", "$\\binom41=4$", " 種選法，係數 ", "$\\frac{4u}N$", wait=2.2)
    s.play(FadeIn(gb), run_time=1.0)
    s.say("(b) 兩隻手在殼層：", "$\\binom42=6$", " 種，係數 ", "$\\frac{6u}N$", wait=2.2)
    s.play(FadeIn(ptag), Write(prop), run_time=1.4)
    s.say("對殼層做 Gaussian 平均：殼層的手兩兩牽起來", wait=2.2)
    s.play(ShowCreation(pair), run_time=1.2)
    s.say("和第 6 章的 Wick 配對同一件事", wait=2.0)
    s.clear_stage()


def part_three(s):
    nb = numbers()
    kinds = (("111", (1, 1, 1), "type 1"), ("21", (2, 1), "type 1 + type 2"), ("3", (3,), "type 3"))
    panels = VGroup()
    for kind, key, name in kinds:
        pic = type_panel(kind, ORIGIN)
        lab = VGroup(zh(name, font_size=22, color=GREY_A),
                     mixed("$%d$" % nb["c3"][key], " 種牽法", font_size=24),
                     mixed("片段圖 ", "$%d$" % nb["d3"][key], " 張", font_size=22, color=GREY_B)).arrange(DOWN, buff=0.12)
        panels.add(VGroup(pic, lab).arrange(DOWN, buff=0.3))
    panels.arrange(RIGHT, buff=1.4)
    tag = zh("Fig. 12:S-2　3 個 (b) 片段，每個 2 隻殼層手", font_size=22, color=GREY_B)
    total = tx(r"15=%d+%d+%d" % (nb["c3"][(1, 1, 1)], nb["c3"][(2, 1)], nb["c3"][(3,)]), 34).set_color(YELLOW)
    col = VGroup(tag, panels, total).arrange(DOWN, buff=0.4)
    place_in(col, -6.6, 6.6, -0.6, 3.1)
    aut = VGroup(tx(r"f_3=3!=%d" % nb["f3"], 32), tx(r"f_4=%d\ne4!=24" % nb["f4"], 32)).arrange(RIGHT, buff=1.2)
    aut.move_to(np.array([0, -1.8, 0]))
    assert_in_frame(col, aut)
    assert_apart(col, aut)
    s.say("每個 (b) 片段有兩隻殼層手；全部牽好之後連成一個個環", wait=2.2)
    s.play(FadeIn(tag), FadeIn(panels[0]), run_time=1.0)
    s.say("type 1：自己的兩隻手相牽", wait=1.6)
    s.play(FadeIn(panels[1]), run_time=0.8)
    s.say("type 2：兩個片段以雙線相連", wait=1.6)
    s.play(FadeIn(panels[2]), run_time=0.8)
    s.say("type 3：三個片段圍成三角環", wait=1.6)
    s.play(Write(total), run_time=0.8)
    s.say("6 隻手的牽法共 ", "$5!!=15$", " 種：", "$1+6+8$", wait=2.2)
    s.say("不同的片段圖只有 1、3、1 張", wait=1.8)
    s.play(FadeIn(aut[0]), run_time=0.8)
    s.say("Ex. 12.11：三角形的 ", "$3!$", " 種換頂點都不破壞配對", wait=2.2)
    s.play(FadeIn(aut[1]), run_time=0.8)
    s.say("4 個片段的環只有 8 種換法，不是 ", "$4!=24$", wait=2.0)
    s.clear_stage()


def part_five(s):
    nb = numbers()
    minis = VGroup(*[mini_diagram(p, ORIGIN) for p in nb["pairs"]])
    minis.arrange_in_grid(2, 5, h_buff=0.9, v_buff=0.5)
    tag = mixed("Fig. 12:S-5　5 個片段 → 3 個 type 1 ＋ 1 個 type 2", font_size=22, color=GREY_B)
    top = VGroup(tag, minis).arrange(DOWN, buff=0.3)
    place_in(top, -6.6, 6.6, 0.0, 3.1)
    cnt = mixed("$%d$" % nb["n5"], " 種牽法裡有 ", "$%d$" % nb["t5"], " 種 ＝ ", "$%d$" % nb["d5"],
                " 張片段圖 × 2 種雙線接法", font_size=26)
    form = tx(r"\frac{5!}{3!\,2!}=%d" % E.pair_factor(3, 1), 34).set_color(YELLOW)
    mid = VGroup(cnt, form).arrange(RIGHT, buff=0.6).move_to(np.array([0, -0.55, 0]))
    box = D.erratum_box("p. 320", r"\frac{5!}{2!\,1!}", r"\frac{5!}{3!\,2!}=10", width=7.0, font_size=30)
    box.move_to(np.array([0, -2.1, 0]))
    assert_in_frame(top, mid, box)
    assert_column(top, mid, box)
    s.say("Ex. 12.9：5 個片段，3 個自己成圈、2 個以雙線相連", wait=2.2)
    s.play(FadeIn(tag), LaggedStart(*[FadeIn(m) for m in minis], lag_ratio=0.12), run_time=2.4)
    s.say("選哪兩個片段連雙線：", "$\\binom52=10$", " 張圖", wait=2.0)
    s.play(FadeIn(cnt), run_time=0.8)
    s.say("電腦列舉全部 945 種牽法：恰好 20 種是這個形狀", wait=2.2)
    s.say("每張圖的雙線有 2 種接法：", "$20=10\\times2$", wait=2.0)
    s.play(Write(form), run_time=0.8)
    s.play(FadeIn(box), run_time=0.8)
    s.say("照書上是 ", "$%d$" % (math.factorial(5) // 2), "；同一頁前面寫對了", wait=2.2)
    s.clear_stage()


def part_general(s):
    nb = numbers()
    head = [r"(n_1,n_2)", "n", r"\frac{n!}{n_1!\,n_2!\,(2!)^{n_2}}", zh("列舉", font_size=24, color=YELLOW), ""]
    rows = [head]
    for p in PAIR_CASES:
        rows.append([r"(%d,%d)" % p, str(p[0] + 2 * p[1]), str(nb["pf"][p]), str(nb["dc"][p]), check_mark(size=0.26)])
    tab = D.tex_table(rows, [1.5, 0.7, 2.3, 1.1, 0.6], row_h=0.5, font_size=26)
    ttag = zh("(12:S-4)　Ex. 12.10", font_size=22, color=GREY_B)
    left = VGroup(ttag, tab).arrange(DOWN, buff=0.25)
    place_in(left, -6.7, -0.4, -2.85, 3.1)
    expo = tx(r"\sum_{n_1,n_2}\frac{\mathcal T_1^{n_1}\mathcal T_2^{n_2}}{n_1!\,n_2!}=e^{\mathcal T_1+\mathcal T_2}", 34)
    etag = zh("(12:S-3)(12:S-5)", font_size=22, color=GREY_B)
    box = D.erratum_box("p. 322，(12:S-7)", r"\mathcal T_1^{n_1}\mathcal T_3^{n_3}\mathcal T_3^{n_3}",
                        r"\mathcal T_1^{n_1}\mathcal T_2^{n_2}\mathcal T_3^{n_3}", width=6.2, font_size=30)
    right = VGroup(etag, expo, box).arrange(DOWN, buff=0.45)
    place_in(right, 0.0, 6.7, -2.6, 2.9)
    assert_in_frame(left, right)
    assert_apart(left, right)
    assert_column(*right)
    assert_table(tab, 5)
    s.say("一般的片段圖：", "$n_1$", " 個 type 1、", "$n_2$", " 個 type 2", wait=2.0)
    s.play(FadeIn(ttag), FadeIn(tab[1]), FadeIn(VGroup(*tab[0][:5])), run_time=1.0)
    s.say("$n!$", " 種排法，同型的片段互換、雙線兩端互換都要除掉", wait=2.4)
    s.play(LaggedStart(*[FadeIn(VGroup(*tab[0][5 * i:5 * i + 5])) for i in range(1, 7)], lag_ratio=0.3), run_time=2.4)
    s.say("6 種形狀：公式和電腦列舉全部一致", wait=2.0)
    s.play(FadeIn(etag), Write(expo), run_time=1.4)
    s.say("把所有圖加起來，階乘剛好變成指數函數的係數", wait=2.2)
    s.say("所有圖的和＝", "$\\exp\\big($", "連通圖的和", "$\\big)$", wait=2.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("$n_1$", " 也要從 0 開始加；", "$\\mathcal T_2$", " 被印成了 ", "$\\mathcal T_3$", wait=2.2)
    s.clear_stage()


def part_check(s):
    nb = numbers()
    rows = [[zh("算法", font_size=24, color=YELLOW), r"\langle e^X\rangle"]]
    rows += [[zh(name, font_size=24), f4(nb["avg"][k])] for k, name in METHODS]
    tab = D.tex_table(rows, [2.8, 1.6], row_h=0.5, font_size=28)
    ttag = mixed("Ex. 12.1：", "$X=-as^2$", "，", "$(a,\\sigma^2)=(0.1,1)$", font_size=22, color=GREY_B)
    closed = tx(r"(1+2a\sigma^2)^{-1/2}=%s" % f4(nb["avg"]["closed"]), 30).set_color(YELLOW)
    left = VGroup(ttag, tab, closed).arrange(DOWN, buff=0.3)
    place_in(left, -6.7, -1.0, -0.2, 3.1)
    tot, con, dis = nb["ord2"]
    o2 = mixed("二階：全部 ", "$%s$" % f4(tot), " ＝ 連通 ", "$%s$" % f4(con), " ＋ 不連通 ", "$%s$" % f4(dis),
               font_size=26)
    kap = mixed("連通圖 ", "$\\kappa_k$", "：", "$" + ",\\ ".join(f4(k) for k in nb["kap"]) + "$", font_size=26)
    right = VGroup(o2, kap).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
    place_in(right, -0.6, 6.7, 0.4, 2.6)
    egf = tx(r"\sum_n(2n-1)!!\frac{x^n}{n!}=\exp\Big(\sum_kw_k\frac{x^k}{k!}\Big)=(1-2x)^{-1/2}", 32)
    w = mixed("$w_k=%s$" % ",".join(str(v) for v in nb["w"]), "；", "$x=0.1$", "：", "$%s$" % f4(nb["egf"][0]),
              "；三個殼層變數：", "$%s$" % f4(nb["multi"]), font_size=26)
    bottom = VGroup(egf, w).arrange(DOWN, buff=0.3)
    place_in(bottom, -6.7, 6.7, -2.85, -0.5)
    assert_in_frame(left, right, bottom)
    assert_apart(left, right)
    assert_apart(left, bottom)
    assert_apart(right, bottom)
    assert_table(tab, 2)
    s.say("數值核對：一個殼層變數 ", "$s$", "，", "$X=-as^2$", wait=2.0)
    s.play(FadeIn(ttag), FadeIn(tab), run_time=1.0)
    s.say("數值積分、所有圖、只用連通圖取指數：都是 ", "$%s$" % f4(nb["avg"]["quad"]), wait=2.4)
    s.play(FadeIn(closed), run_time=0.8)
    s.say("閉式也一樣", wait=1.4)
    s.play(FadeIn(o2), run_time=0.8)
    s.say("二階：不連通的那部分正好是一階平方的一半", wait=2.2)
    s.play(FadeIn(kap), run_time=0.8)
    s.say("連通圖就是累積量 ", "$\\kappa_k$", "：每一個都是一個環", wait=2.2)
    s.play(Write(egf), run_time=1.6)
    s.say("一行的版本：所有配對＝連成環的配對取指數", wait=2.2)
    s.play(FadeIn(w), run_time=0.8)
    s.say("$x=0.1$", "：三種寫法都是 ", "$%s$" % f4(nb["egf"][0]), "；三個變數也成立", wait=2.2)
    s.say("指數化自動把不連通的圖算好——只要畫連通的圖", wait=2.2)
    s.say("Fig. 12.11 的兩張一迴圈圖，就是下一幕的主角", wait=2.0)


class W1202LinkedClusters(ChuScene):
    def construct(self):
        self.show_header("W1202　配對遊戲：只剩連通的圖")
        check_math()
        part_question(self)
        part_three(self)
        part_five(self)
        part_general(self)
        part_check(self)
        self.wait(1.0)
