"""
W705 微擾論為什麼失效、平均場（§7.2.5–7.2.6，pp. 184–191，(7.38)–(7.61)，Ex. 7.4–7.5）

問：自能的每一階都帶一個小小的 f。是不是只要算前幾階就夠了？
  1. 數冪次 (7.38)(7.39)：l 階＝2l−1 條傳播子、l 個 d^dq → (z_c−z)^{1−l(4−d)/2}；表 d = 1..5、l = 0..4（SW.power_table）；
     數值積分量出的冪次 −0.50、0.50、−2.00（SW.loop_slope）。
  2. 一維一階自能 1/√(1−4z²)：z = 0.4、0.49、0.499 → 1.6667、5.0252、15.8193；(7.42) 截到 n_h 階的奇點。
  3. 平均場 (7.51)(7.57)(7.59)(7.60)；筆誤框：(7.59) 指數裡的 𝒩λ。
  4. Ex. 7.4（0.3385、0.3643、0.3675 → e^{−1} = 0.3679）、Ex. 7.5（0.6899 … 0.9992 → 1；0.0452 … 0.0621 → 0.0625）。
畫面上的數字都由 utils.walk_saw 算出，check_math() 核對。模組層級的畫圖函式：chord_diagram（弧線圖）。
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
from utils import walk_saw as SW
from utils import walk_draw as D

ARC_COLOR = YELLOW
ZS = (0.4, 0.49, 0.499)
NS75 = (10, 10 ** 2, 10 ** 3, 10 ** 4)
NS75B = (10 ** 2, 10 ** 3, 10 ** 4)
SLOPES = (("l1d1", 1, 1), ("l1d3", 1, 3), ("l2d1", 2, 1))          # (case, l, d)
CARD_POWER = {1: ["1", "-1/2", "-2", "-7/2", "-5"], 2: ["1", "0", "-1", "-2", "-3"], 3: ["1", "1/2", "0", "-1/2", "-1"],
              4: ["1", "1", "1", "1", "1"], 5: ["1", "3/2", "2", "5/2", "3"]}


def f2(x):
    return f"{x:.2f}"


def f4(x):
    return f"{x:.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = {"power": SW.power_table(), "slope": [SW.loop_slope(c) for c, _, _ in SLOPES],
              "s1": [SW.sigma1_1d(z) for z in ZS], "s1q": [SW.sigma1_1d_quad(z) for z in ZS],
              "mf": [SW.cmf_ratio(*c) for c in SW.MF_CASES], "asym": [SW.cmf_asym(*c) for c in SW.MF_CASES],
              "mf0": SW.cmf_ratio(10, 100, 0), "book0": SW.cmf_book59(10, 100, 0),
              "ex75": [SW.ex75(n) for n in NS75], "ex75b": [SW.ex75_big(n) for n in NS75B]}
    s1 = nb["s1"]
    nb["ratio"] = [s1[1] / s1[0], s1[2] / s1[1]]
    return nb


def check_math():
    nb = numbers()
    pt = nb["power"]
    assert {d: [str(v) for v in row] for d, row in pt.items()} == CARD_POWER
    for d in range(1, 6):
        for l in range(5):
            assert SW.power_exp(l, d) == SW.power_exp_count(l, d) == pt[d][l] == 1 - Fraction(l * (4 - d), 2)
    assert [f2(v) for v in nb["slope"]] == ["-0.50", "0.50", "-2.00"]
    assert all(f2(v) == f2(float(pt[d][l])) for v, (_, l, d) in zip(nb["slope"], SLOPES))
    assert [f4(v) for v in nb["s1"]] == ["1.6667", "5.0252", "15.8193"]
    assert all(abs(a - b) < 1e-6 for a, b in zip(nb["s1"], nb["s1q"]))
    assert all(abs(v - 1 / math.sqrt(1 - 4 * z * z)) < 1e-12 for v, z in zip(nb["s1"], ZS))
    assert [f2(v) for v in nb["ratio"]] == ["3.02", "3.15"] and f2(math.sqrt(10)) == "3.16"
    assert all(SW.power_exp(l + 1, 3) - SW.power_exp(l, 3) == Fraction(-1, 2) for l in range(4))
    assert all(SW.power_exp(l, 4) == 1 for l in range(5)) and all(SW.power_exp(l + 1, 5) > SW.power_exp(l, 5) for l in range(4))
    assert [f4(v) for v in nb["mf"]] == ["0.3385", "0.3643", "0.3675"]
    assert all(f4(v) == "0.3679" for v in nb["asym"]) and f4(math.exp(-1)) == "0.3679"
    assert all(f * m * m == Ns for m, Ns, f in SW.MF_CASES)
    assert [(m, Ns) for m, Ns, _ in SW.MF_CASES] == [(10, 100), (100, 10 ** 4), (1000, 10 ** 6)]
    assert abs(nb["mf0"] - 1) < 1e-9 and nb["book0"] < 1e-20                     # f = 0：更正版回到 1，書上不是
    assert [f4(v) for v in nb["ex75"]] == ["0.6899", "0.9341", "0.9922", "0.9992"]
    assert [f4(v) for v in nb["ex75b"]] == ["0.0452", "0.0596", "0.0621"]
    assert (Fraction(3, 2) - 1) ** 2 / 4 == Fraction(1, 16) and f4(1 / 16) == "0.0625"
    assert SW.is_irreducible(((1, 2),)) and SW.is_irreducible(((1, 4), (2, 3))) and SW.is_irreducible(((1, 4), (2, 5), (3, 6)))


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


def chord_diagram(m, width=3.0, origin=ORIGIN, color=ARC_COLOR, angle=-0.75 * PI):
    """弧線圖：水平實線上 2k 個點（編號 1..2k），每對 (a, b) 一條上方的虛線弧。回傳 VGroup(線, 點們, 弧們)。"""
    origin = np.array(origin, dtype=float)
    n = 2 * len(m)
    gap = width / (n + 1)
    pts = [origin + (-width / 2 + gap * (i + 1)) * RIGHT for i in range(n)]
    line = Line(origin + width / 2 * LEFT, origin + width / 2 * RIGHT).set_stroke(WHITE, 3)
    dots = VGroup(*[Dot(p, radius=0.06).set_color(WHITE) for p in pts])
    arcs = VGroup()
    for a, b in m:
        span = np.linalg.norm(pts[b - 1] - pts[a - 1])
        arc = DashedVMobject(ArcBetweenPoints(pts[a - 1], pts[b - 1], angle=angle), num_dashes=max(5, int(1.3 * span / 0.13)))
        arcs.add(arc.set_stroke(color, 2.5))
    return VGroup(line, dots, arcs)


def eq_tag(tex_mob, tag):
    return VGroup(zh(tag, font_size=22, color=GREY_B), tex_mob).arrange(RIGHT, buff=0.35)


# ---------------------------------------------------------------- 各段
def part_question(s):
    sig = tx(r"\Sigma=", 40)
    items = VGroup(sig)
    for k, (m, w) in enumerate(((((1, 2),), 1.6), (((1, 4), (2, 3)), 2.0), (((1, 4), (2, 5), (3, 6)), 2.6))):
        dg = chord_diagram(m, width=w)
        fl = tx("f" if k == 0 else "f^%d" % (k + 1), 32, ARC_COLOR).next_to(dg, DOWN, buff=0.3)
        items.add(VGroup(dg, fl))
        items.add(tx("+", 40))
    items.add(tx(r"\cdots", 40))
    items.arrange(RIGHT, buff=0.35)
    q = mixed("$f$", " 很小：只算前幾階就夠了嗎？", font_size=34, color=YELLOW)
    place_in(items, -6.7, 6.7, 0.0, 2.6)
    place_in(q, -6.7, 6.7, -2.5, -1.2)
    assert_in_frame(items, q)
    assert_apart(items, q)
    for a, b in zip(items, items[1:]):
        assert_apart(a, b)
    s.say("問：自能的每一階都帶一個小小的 ", "$f$", wait=0.5)
    s.play(FadeIn(items[0]), run_time=0.5)
    for k in range(3):
        g = items[2 * k + 1]
        s.play(ShowCreation(g[0][0]), FadeIn(g[0][1]), ShowCreation(g[0][2]), FadeIn(g[1]), FadeIn(items[2 * k + 2]),
               run_time=0.9)
    s.play(FadeIn(items[-1]), run_time=0.4)
    s.say("第 ", "$l$", " 階帶 ", "$f^l$", "：看起來越高階越小", wait=2.5)
    s.play(FadeIn(q), run_time=0.8)
    s.say("是不是只要算前幾階就夠了？（先猜）", wait=3.0)
    s.clear_stage()


def part_power(s):
    nb = numbers()
    pt = nb["power"]
    l1 = mixed("第 ", "$l$", " 階：", "$2l-1$", " 條傳播子 ", "$\\frac{1}{(z_c-z)+q^2}$", "、", "$l$", " 個 ", "$d^dq$",
               font_size=28)
    l2 = mixed("令 ", "$q\\to(z_c-z)^{1/2}q$", "：", "$(z_c-z)^{-(2l-1)}\\,(z_c-z)^{dl/2}=(z_c-z)^{1-l(4-d)/2}$",
               font_size=28)
    head = VGroup(l1, l2).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    rows = [["d"] + ["l=%d" % l for l in range(5)]]
    for d in range(1, 6):
        rows.append(["d=%d" % d] + [ftex(v) for v in pt[d]])
    tab = D.tex_table(rows, [1.3, 1.2, 1.2, 1.2, 1.2, 1.2], row_h=0.66, font_size=28)
    for d in range(1, 6):
        for l in range(5):
            v = pt[d][l]
            c = YELLOW if d == 4 else (RED_B if v < 0 else (GREY_B if v == 0 else GREEN_B))
            D.table_cell(tab, d, l + 1, 6).set_color(c)
    sl = VGroup(zh("數值積分量出來", font_size=26, color=GREY_B))
    for v, (_, l, d) in zip(nb["slope"], SLOPES):
        sl.add(tx(r"l=%d,\ d=%d:\ \ %s" % (l, d, f2(v)), 30, TEAL_B))
    sl.arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    place_in(head, -6.7, 6.7, 1.75, 3.15)
    place_in(tab, -6.7, 1.6, -2.9, 1.5)
    place_in(sl, 2.2, 6.7, -1.6, 0.9)
    boxes = VGroup(*[SurroundingRectangle(D.table_cell(tab, d, l + 1, 6), buff=0.07).set_stroke(TEAL_B, 2)
                     for _, l, d in SLOPES])
    assert_in_frame(head, tab, sl)
    assert_column(head, tab)
    assert_column(head, sl)
    assert_apart(tab, sl)
    assert_table(tab, 6)
    s.say("數冪次 (7.38)：第 ", "$l$", " 階有 ", "$2l-1$", " 條傳播子、", "$l$", " 個積分", wait=0.5)
    s.play(FadeIn(l1), run_time=1.0)
    s.wait(2.0)
    s.play(FadeIn(l2), run_time=1.0)
    s.say("令 ", "$q\\to(z_c-z)^{1/2}q$", "：每條傳播子給 ", "$(z_c-z)^{-1}$", wait=2.5)
    s.say("每個積分給 ", "$(z_c-z)^{d/2}$", "：合起來 ", "$1-l(4-d)/2$", wait=2.5)
    s.play(FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:]], lag_ratio=0.04), run_time=1.5)
    s.say("表 (7.39)：", "$d=4$", " 每一階都是 1", wait=2.5)
    s.say("$d<4$", "：階數越高，冪次越負（紅色）", wait=2.5)
    s.say("負的冪次：在 ", "$z\\to z_c$", " 時爆掉的奇點", wait=2.5)
    s.say("$d>4$", "：階數越高，冪次反而越大", wait=2.5)
    s.play(FadeIn(sl), ShowCreation(boxes), run_time=1.0)
    s.say("數值積分量出來：", "$-0.50$", "、", "$0.50$", "、", "$-2.00$", "，和表一致", wait=3.0)
    s.clear_stage()


def part_sigma1d(s):
    nb = numbers()
    eq = tx(r"\Sigma_1=\int_{-\pi}^{\pi}\frac{dq}{2\pi}\,\frac{1}{1-2z\cos q}=\frac{1}{\sqrt{1-4z^2}}", 38)
    rows = [["z"] + [str(z) for z in ZS],
            ["z_c-z", "0.1", "0.01", "0.001"],
            [r"\Sigma_1"] + [f4(v) for v in nb["s1"]]]
    tab = D.tex_table(rows, [1.8, 1.9, 1.9, 1.9], row_h=0.62, font_size=30)
    for j in range(1, 4):
        D.table_cell(tab, 2, j, 4).set_color(YELLOW)
    rat = VGroup(zh("離 1/2 每近 10 倍，大了", font_size=26, color=GREY_B),
                 tx(r"\times%s,\ \ \times%s\ \ \to\ \sqrt{10}=%s" % (f2(nb["ratio"][0]), f2(nb["ratio"][1]),
                                                                  f2(math.sqrt(10))), 32)).arrange(DOWN, buff=0.2)
    p42 = VGroup(mixed("(7.42) 截到第 ", "$n_h$", " 階：奇點 ", "$(z_c-z)^{-p}$", "，", "$p=1-n_h(4-d)/2$", font_size=28),
                 mixed("$d<4$", "：每多一階 ", "$p$", " 就變小（", "$d=3$", " 每階少 ", "$\\frac12$", "）",
                       font_size=26, color=RED_B),
                 mixed("$d=4$", "：分界", font_size=26, color=YELLOW),
                 mixed("$d>4$", "：高階反而不要緊（腳註 1）", font_size=26, color=GREEN_B))
    p42.arrange(DOWN, buff=0.18, aligned_edge=LEFT)
    place_in(eq, -6.7, 6.7, 2.2, 3.15)
    place_in(tab, -6.7, 1.0, -0.15, 2.0)
    place_in(rat, 1.5, 6.7, 0.35, 1.5)
    place_in(p42, -6.7, 6.7, -2.9, -0.4)
    assert_in_frame(eq, tab, rat, p42)
    assert_column(eq, tab, p42)
    assert_column(eq, rat, p42)
    assert_apart(tab, rat)
    assert_table(tab, 4)
    assert_column(*p42)
    s.say("一維的一階自能：一個積分，有閉式", wait=0.5)
    s.play(Write(eq), run_time=1.5)
    s.wait(1.5)
    s.play(FadeIn(tab[0][:8]), ShowCreation(tab[1]), run_time=0.8)
    s.say("$z_c=\\frac12$", "：分母 ", "$\\sqrt{1-4z^2}$", " 在那裡變成 0", wait=2.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][8:]], lag_ratio=0.25), run_time=1.0)
    s.say("$z=0.4,\\ 0.49,\\ 0.499$", "：1.6667、5.0252、15.8193", wait=2.5)
    s.play(FadeIn(rat), run_time=0.8)
    s.say("離 ", "$z_c$", " 每近 10 倍，就大約大 ", "$\\sqrt{10}$", " 倍", wait=2.5)
    s.play(FadeIn(p42[0]), run_time=0.8)
    s.say("(7.42)：截到第 ", "$n_h$", " 階，奇點 ", "$(z_c-z)^{-p}$", wait=2.5)
    s.play(FadeIn(p42[1]), run_time=0.7)
    s.say("$d<4$", "：每多一階 ", "$p$", " 就變小", wait=2.5)
    s.say("$f$", " 再小也壓不住：越靠近 ", "$z_c$", "，高階項越大", wait=3.0)
    s.play(FadeIn(p42[2]), FadeIn(p42[3]), run_time=0.8)
    s.say("$d=4$", " 是分界；", "$d>4$", " 高階反而不要緊（腳註 1）", wait=3.0)
    s.clear_stage()


def part_meanfield(s):
    e1 = eq_tag(tx(r"G(\vec x,\vec y)\approx z_c\frac{\psi_0(\vec x)\,\psi_0^*(\vec y)}{z_0-z},\qquad"
                   r"\psi_0=\frac{1}{\sqrt{\mathcal N}}", 32), "(7.51)")
    e2 = eq_tag(tx(r"G_{MF}=\int_0^\infty\exp\Big[\mathcal N\Big(-f\lambda^2-\lambda(1-z/z_c)\Big)\Big]d\lambda", 32),
                "(7.57)")
    e3 = eq_tag(tx(r"C_{MF}(m)=\frac{\mathcal N^m}{z_c^m\,m!}\int_0^\infty e^{-\mathcal N(f\lambda^2+\lambda)+m\ln\lambda}"
                   r"d\lambda\approx\frac{1}{\mathcal N}z_c^{-m}e^{-fm^2/\mathcal N}", 32), "(7.59)(7.60)")
    eqs = VGroup(e1, e2, e3).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    box = D.erratum_box("p. 190，(7.59)", r"e^{-\mathcal N(f\lambda^2+\mathcal N\lambda)+m\ln\lambda}",
                        r"e^{-\mathcal N(f\lambda^2+\lambda)+m\ln\lambda}", width=8.0, font_size=46)
    place_in(eqs, -6.7, 6.7, -0.45, 3.15)
    place_in(box, -6.7, 6.7, -2.9, -0.75)
    assert_in_frame(eqs, box)
    assert_column(*eqs, box)
    hl = SurroundingRectangle(e3[1][-14:], buff=0.08).set_stroke(YELLOW, 2)
    s.say("平均場：只留最低的一個模態 (7.51)", wait=0.5)
    s.play(FadeIn(e1), run_time=1.0)
    s.say("均勻的模態 ", "$\\psi_0=1/\\sqrt{\\mathcal N}$", "：每個格點都一樣", wait=2.5)
    s.play(FadeIn(e2), run_time=1.0)
    s.say("(7.57)：整個問題只剩一個 ", "$\\lambda$", " 的積分", wait=2.5)
    s.play(FadeIn(e3), run_time=1.0)
    s.say("(7.59)：取 ", "$z^m$", " 的係數；大 ", "$m$", " 的近似得 (7.60)", wait=2.5)
    s.play(ShowCreation(hl), run_time=0.6)
    s.say("取 ", "$f=0$", "：回到普通漫步 ", "$z_c^{-m}/\\mathcal N$", wait=2.5)
    s.say("$e^{-fm^2/\\mathcal N}$", "：自迴避讓條數少掉的比例", wait=2.5)
    s.say("$m^2/\\mathcal N$", " 越大，自迴避的效果越強", wait=2.0)
    s.play(FadeIn(box), FadeOut(hl), run_time=0.8)
    s.say("書上 (7.59)：指數裡多了一個 ", "$\\mathcal N$", wait=2.0)
    s.say("照書上，", "$f=0$", " 也回不到普通漫步的 ", "$z_c^{-m}/\\mathcal N$", wait=3.0)
    s.clear_stage()


def part_examples(s):
    nb = numbers()
    ta = D.tex_table([[r"(m,\mathcal N)", "(10,100)", "(100,10^4)", "(1000,10^6)", r"e^{-fm^2/\mathcal N}"],
                      [r"\mathcal N z_c^m C_{MF}(m)"] + [f4(v) for v in nb["mf"]] + [f4(nb["asym"][0])]],
                     [3.2, 1.9, 1.9, 1.9, 1.9], row_h=0.6, font_size=28)
    tb = D.tex_table([[r"\mathcal N", "10", "10^2", "10^3", "10^4", r"\to"],
                      [r"\mathcal N(1-z/z_c)\,G_{MF}"] + [f4(v) for v in nb["ex75"]] + ["1"]],
                     [3.6, 1.5, 1.5, 1.5, 1.5, 1.2], row_h=0.6, font_size=28)
    tc = D.tex_table([[r"\mathcal N", "10^2", "10^3", "10^4", r"\to"],
                      [r"\ln G_{MF}/\mathcal N"] + [f4(v) for v in nb["ex75b"]] + [f4(1 / 16)]],
                     [3.6, 1.5, 1.5, 1.5, 1.6], row_h=0.6, font_size=28)
    for t, n in ((ta, 5), (tb, 6), (tc, 5)):
        D.table_cell(t, 1, n - 1, n).set_color(YELLOW)
    ha = mixed("Ex. 7.4：精確積分 vs (7.60)，", "$fm^2/\\mathcal N=1$", font_size=26, color=GREY_B)
    hb = mixed("Ex. 7.5：", "$z/z_c=\\frac12$", font_size=26, color=GREY_B)
    hc = mixed("對照：", "$z/z_c=\\frac32$", "，", "$(z/z_c-1)^2/4f=\\frac1{16}$", font_size=26, color=GREY_B)
    blocks = VGroup(*[VGroup(h, t).arrange(DOWN, buff=0.18) for h, t in ((ha, ta), (hb, tb), (hc, tc))])
    place_in(blocks[0], -6.7, 6.7, 1.35, 3.15)
    place_in(blocks[1], -6.7, 6.7, -0.75, 1.05)
    place_in(blocks[2], -6.7, 6.7, -2.9, -1.1)
    assert_in_frame(blocks)
    assert_column(*blocks)
    assert_table(ta, 5)
    assert_table(tb, 6)
    assert_table(tc, 5)
    for b in blocks:
        assert_apart(b[0], b[1])
    s.say("Ex. 7.4：精確積分和 (7.60) 比一比", wait=0.5)
    s.play(FadeIn(ha), FadeIn(ta[0][:5]), ShowCreation(ta[1]), run_time=1.0)
    s.say("三組的 ", "$fm^2/\\mathcal N$", " 都是 1：(7.60) 預測 ", "$e^{-1}$", wait=2.0)
    s.play(LaggedStart(*[FadeIn(c) for c in ta[0][5:]], lag_ratio=0.25), run_time=1.2)
    s.say("0.3385、0.3643、0.3675：越來越接近 0.3679", wait=2.5)
    s.say("(7.60) 是大 ", "$m$", " 的近似：", "$m$", " 越大越準", wait=2.0)
    s.play(FadeIn(blocks[1]), run_time=1.0)
    s.say("Ex. 7.5：", "$z/z_c=\\frac12$", "，乘上 ", "$\\mathcal N(1-z/z_c)$", wait=2.0)
    s.say("0.6899、0.9341、0.9922、0.9992：趨近 1", wait=2.5)
    s.say("條數就是 ", "$z_c^{-m}/\\mathcal N$", "：自迴避沒有效果", wait=2.5)
    s.play(FadeIn(blocks[2]), run_time=1.0)
    s.say("對照 ", "$z/z_c=\\frac32$", "：", "$\\ln G_{MF}/\\mathcal N\\to0.0625$", wait=3.0)
    s.say("四維以下，低階微擾一定失效", wait=2.5)
    s.say("平均場只留一個模態，得到 ", "$e^{-fm^2/\\mathcal N}$", wait=3.0)


class W705MeanField(ChuScene):
    def construct(self):
        self.show_header("W705　微擾論為什麼失效、平均場")
        check_math()
        part_question(self)
        part_power(self)
        part_sigma1d(self)
        part_meanfield(self)
        part_examples(self)
        self.wait(1.0)
