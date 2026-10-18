"""
W803 最陡下降：從標度形式讀出分布（§8.1.2 後半，pp. 197–199，(8.15)–(8.20)，Ex. 8.2–8.5）

問：知道 G(z;r) 的標度形式，怎麼拿回「N 步走到 r」的條數？
  1. (8.15)–(8.19)：圍道積分 → 極值化 → δ = z_c − z；筆誤框 1：(8.18) 少了負號。
  2. (8.20) C ∝ exp[−c (r/N^ν)^{1/(1−ν)}]；簡單立方係數 d/2（d = 3：3/2）→ 高斯；筆誤框 2：少了 (1 − ν)。
  3. 一維 N = 20 長條圖＋e^{−x²/42}（更正）與 e^{−x²/21}（書上）；表 0.6818｜0.6832｜0.4668 …
  4. Fisher 指數 2、5/2、4（三維 2.43）；γ = (2 − η)ν：1、43/32；μ ≈ 2.638、p = 11/32；比值法表；Ex. 8.5 θ = 11/24。
畫面上的數字都由 utils.walk_scale 算出，check_math() 核對。模組層級的畫圖函式：chart_n20（長條圖＋兩條曲線＋刻度）。
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
from utils import walk_scale as WS
from utils import walk_saw as SW
from utils import walk_draw as D

HALF = Fraction(1, 2)
NU3 = SW.NU_QUOTED[3]                 # 0.588（只引用）


def f2(x):
    return f"{float(x):.2f}"


def f4(x):
    return f"{float(x):.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    nb = {"bars": {x: WS.ratio_exact(20, x) for x in range(-20, 21, 2)},
          "sd": {N: [(x, WS.ratio_exact(N, x), WS.sd_ratio(N, x, 0.5, 2, 0.5), WS.sd_ratio(N, x, 0.5, 2, 0.5, book=True))
                     for x in WS.SD_X[N]] for N in WS.SD_N},
          "fisher": [WS.fisher(nu) for nu in WS.FISHER_NU], "fisher3": 1 / (1 - NU3),
          "g_ord": WS.gamma_exp(0, HALF), "g_2d": WS.gamma_exp(WS.ETA_2D, WS.NU_2D),
          "growth": WS.growth(WS.GAMMA_QUOTED, 1 / WS.MU_QUOTED),
          "ratio": {N: WS.ratio_method(N) for N in WS.RATIO_N},
          "theta": WS.theta(WS.GAMMA_QUOTED, WS.NU_2D)}
    return nb


def check_math():
    nb = numbers()
    card = {20: [("0.6818", "0.6832", "0.4668"), ("0.2098", "0.2179", "0.0475")],
            100: [("0.7282", "0.7285", "0.5306"), ("0.4894", "0.4902", "0.2403")]}
    for N in WS.SD_N:
        assert [(f4(a), f4(b), f4(c)) for _, a, b, c in nb["sd"][N]] == card[N]
        for x, _, b, c in nb["sd"][N]:
            assert abs(b - math.exp(-x * x / (2 * (N + 1)))) < 1e-12 and abs(c - math.exp(-x * x / (N + 1))) < 1e-12
    assert abs(WS.sd_K(0.5, 2, 0.5) - 0.5) < 1e-12 and abs(WS.sd_K(0.5, 2, 0.5, book=True) - 1) < 1e-12
    assert abs(WS.sd_K_dim(3) - 1.5) < 1e-12 and all(abs(WS.sd_K_dim(d) - d / 2) < 1e-12 for d in (1, 2, 3))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert abs(WS.sd_delta(100, 10, 0.5, 2, 0.5) - WS.sd_delta_num(100, 10, 0.5, 2, 0.5)) < 1e-3
    assert nb["fisher"] == [2, Fraction(5, 2), 4] and NU3 == 0.588 and f2(nb["fisher3"]) == "2.43"
    assert nb["g_ord"] == 1 and nb["g_2d"] == Fraction(43, 32) and WS.GAMMA_QUOTED == Fraction(43, 32)
    assert WS.ETA_2D == Fraction(5, 24) and WS.NU_2D == Fraction(3, 4)
    x, p = nb["growth"]
    assert f"{x:.3f}" == "2.638" and p == Fraction(11, 32) and f"{float(p):.5f}" == "0.34375"
    assert [f4(nb["ratio"][N][1]) for N in WS.RATIO_N] == ["2.6592", "2.6529", "2.6483", "2.6454"]
    assert [f4(nb["ratio"][N][2]) for N in WS.RATIO_N] == ["1.3552", "1.3519", "1.3496", "1.3483"]
    assert nb["theta"] == Fraction(11, 24)
    assert len(nb["bars"]) == 21 and abs(nb["bars"][0] - 1) < 1e-12


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


def chart_n20(ux=0.15, uy=3.4, origin=(-3.6, -2.2, 0)):
    """一維 N = 20：長條（精確比值）＋更正 e^{−x²/42}（綠）＋書上 e^{−x²/21}（紅虛線）＋刻度。"""
    nb = numbers()
    origin = np.array(origin, dtype=float)
    bars = D.bar_chart(nb["bars"], ux=ux, uy=uy, origin=origin, width=1.5)
    good = D.curve(lambda x: math.exp(-x * x / 42), -20, 20, ux, uy, origin, color=D.RIGHT_COLOR, samples=80)
    bad = D.curve(lambda x: math.exp(-x * x / 21), -20, 20, ux, uy, origin, color=D.WRONG_COLOR, samples=80,
                  dashed=True)
    base = Line(origin + 21 * ux * LEFT, origin + 21 * ux * RIGHT).set_stroke(GREY_B, 2)
    ticks = VGroup(*[tx(str(x), 20, GREY_B).next_to(origin + x * ux * RIGHT, DOWN, buff=0.12) for x in (-20, -10, 0, 10, 20)])
    one = VGroup(DashedLine(origin + 21 * ux * LEFT + uy * UP, origin + 21 * ux * RIGHT + uy * UP).set_stroke(GREY_C, 1),
                 tx("1", 20, GREY_B).next_to(origin + 21 * ux * LEFT + uy * UP, LEFT, buff=0.1))
    xl = tx("x", 24, GREY_B).next_to(origin + 21 * ux * RIGHT, RIGHT, buff=0.1)
    return VGroup(bars, good, bad, VGroup(base, ticks, one, xl))


# ---------------------------------------------------------------- 各段
def part_contour(s):
    e15 = tx(r"C(N;\vec r)=\frac{1}{2\pi i}\oint\frac{G(z;\vec r)}{z^{N+1}}\,dz", 34)
    e17 = mixed("極值化 ", "$A\\Big(|\\vec r|(z_c-z)^{\\nu}\\Big)-(N+1)\\ln z$", "，", "$A=\\ln F\\approx-Bx$", font_size=30)
    e18 = mixed("(8.18) ", "$-\\nu|\\vec r|(z_c-z)^{\\nu-1}A'-\\frac{N+1}{z}=0$", font_size=30)
    e19 = mixed("(8.19) ", "$\\delta=z_c-z\\approx\\Big(\\frac{N+1}{\\nu Bz_c|\\vec r|}\\Big)^{1/(\\nu-1)}$", font_size=30)
    col = VGroup(e15, e17, e18, e19).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -0.75, 3.1)
    box = D.erratum_box("p. 198，(8.18)", r"\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}{z}=0",
                        r"-\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}{z}=0", width=11.0, font_size=34)
    place_in(box, -6.7, 6.7, -2.9, -1.0)
    assert_in_frame(col, box)
    assert_column(*col, box)
    s.say("問：知道 ", "$G(z;\\vec r)$", " 的標度形式……", wait=0.5)
    s.say("怎麼拿回「", "$N$", " 步走到 ", "$\\vec r$", "」的條數？", wait=1.0)
    s.play(Write(e15), run_time=1.5)
    s.say("(8.15)：條數是 ", "$z^N$", " 的係數——繞原點的圍道積分", wait=2.5)
    s.play(FadeIn(e17), run_time=1.0)
    s.say("$N$", " 大：被積函數在鞍點最大，只要把指數極值化", wait=2.5)
    s.say("大 ", "$x$", " 時 ", "$F$", " 指數衰減：", "$A=\\ln F\\approx-Bx$", wait=2.5)
    s.play(FadeIn(e18), run_time=0.8)
    s.say("對 ", "$z$", " 微分 (8.18)；", "$A'=-B$", " 時第一項是正的", wait=2.5)
    s.play(FadeIn(e19), run_time=0.8)
    s.say("解出 (8.19)：鞍點離 ", "$z_c$", " 的距離 ", "$\\delta$", wait=2.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (8.18) 少了開頭的負號", wait=2.0)
    s.say("照書上 ", "$A'=-B$", " 時兩項都是負的，沒有解", wait=3.0)
    s.clear_stage()


def part_820(s):
    e20 = tx(r"C(N;\vec r)\propto\exp\Big[-(1-\nu)\big((\nu z_c)^\nu B\big)^{\frac{1}{1-\nu}}"
             r"(N+1)^{-\frac{\nu}{1-\nu}}|\vec r|^{\frac{1}{1-\nu}}\Big]", 34)
    short = tx(r"=\exp\Big[-c\,\big(r/N^{\nu}\big)^{1/(1-\nu)}\Big]", 34)
    chk = mixed("簡單立方 ", "$\\nu=\\frac12,\\ B=2d,\\ z_c=\\frac{1}{2d}$", " → 係數 ", "$\\frac d2$",
                "（", "$d=3$", "：", "$%s$" % ftex(Fraction(WS.sd_K_dim(3)).limit_denominator(10)), "）→ ",
                "$e^{-dr^2/2N}$", font_size=28)
    note = zh("＝第 2 章的高斯", font_size=26, color=GREEN_B)
    col = VGroup(e20, short, chk, note).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -0.6, 3.1)
    box = D.erratum_box("p. 198，(8.20)", r"e^{-((\nu z_c)^\nu B)^{1/(1-\nu)}\cdots}",
                        r"e^{-(1-\nu)((\nu z_c)^\nu B)^{1/(1-\nu)}\cdots}", width=10.0, font_size=36)
    place_in(box, -6.7, 6.7, -2.9, -0.9)
    assert_in_frame(col, box)
    assert_column(*col, box)
    s.say("把鞍點代回去，得到 (8.20)", wait=0.3)
    s.play(Write(e20), run_time=2.0)
    s.say("指數裡只有 ", "$r/N^\\nu$", " 這個組合", wait=2.0)
    s.play(FadeIn(short), run_time=0.8)
    s.say("距離用 ", "$N^\\nu$", " 當尺，指數是 ", "$\\frac{1}{1-\\nu}$", " 次方", wait=2.5)
    s.play(FadeIn(chk), run_time=1.0)
    s.say("檢驗：簡單立方 ", "$\\nu=\\frac12$", "，係數是 ", "$\\frac d2$", wait=2.5)
    s.play(FadeIn(note), run_time=0.6)
    s.say("$e^{-dr^2/2N}$", "：正是第 2 章的高斯", wait=2.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (8.20) 少了 ", "$(1-\\nu)$", " 這個因子", wait=2.5)
    s.say("照書上一維得 ", "$e^{-x^2/N}$", "，變異數少一半", wait=2.5)
    s.say("推導第一行的 ", "$-(N+1)/z_c$", " 也應為 ", "$+(N+1)\\delta/z_c$", wait=3.0)
    s.clear_stage()


def part_1d(s):
    nb = numbers()
    ch = chart_n20()
    rows = [["N", "x", zh("精確", font_size=24, color=YELLOW), zh("更正", font_size=24, color=GREEN_B),
             zh("書上", font_size=24, color=RED_B)]]
    for N in WS.SD_N:
        for x, a, b, c in nb["sd"][N]:
            rows.append([str(N), str(x), f4(a), f4(b), f4(c)])
    tab = D.tex_table(rows, [0.9, 0.7, 1.3, 1.3, 1.3], row_h=0.6, font_size=28)
    for i in range(1, 5):
        D.table_cell(tab, i, 3, 5).set_color(GREEN_B)
        D.table_cell(tab, i, 4, 5).set_color(RED_B)
    leg = VGroup(mixed("綠：", "$e^{-x^2/42}$", "（更正）", font_size=26, color=GREEN_B),
                 mixed("紅：", "$e^{-x^2/21}$", "（書上）", font_size=26, color=RED_B)).arrange(DOWN, buff=0.15, aligned_edge=LEFT)
    head = mixed("一維 ", "$N=20$", "：", "$C(N;x)/C(N;0)$", font_size=28)
    place_in(head, -6.7, -0.3, 2.4, 3.1)
    place_in(ch, -6.7, -0.3, -2.9, 2.2)
    right = VGroup(leg, tab).arrange(DOWN, buff=0.4)
    place_in(right, 0.2, 6.7, -2.6, 2.9)
    assert_in_frame(head, ch, right)
    assert_apart(head, ch)
    assert_apart(ch, right)
    assert_apart(head, right)
    assert_table(tab, 5)
    s.say("一維、", "$N=20$", "：精確的條數比值畫成長條", wait=0.3)
    s.play(FadeIn(head), FadeIn(ch[3]), run_time=0.8)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in ch[0][0]], lag_ratio=0.05), run_time=1.5)
    s.play(ShowCreation(ch[1]), FadeIn(leg[0]), run_time=1.2)
    s.say("更正的 (8.20)：", "$e^{-x^2/2(N+1)}=e^{-x^2/42}$", "，貼著長條", wait=2.5)
    s.play(ShowCreation(ch[2]), FadeIn(leg[1]), run_time=1.2)
    s.say("書上的版本 ", "$e^{-x^2/21}$", "：太窄了", wait=2.5)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][5:]], lag_ratio=0.04), run_time=1.5)
    s.say("$N=20,\\ x=8$", "：精確 0.2098、更正 0.2179、書上 0.0475", wait=3.0)
    s.say("$N=100$", "：更正版和精確值差不到千分之一", wait=2.5)
    s.say("一維在 ", "$z=\\pm\\frac12$", " 各有一個鞍點，完整條數要乘 2（取比值時消掉）", wait=3.0)
    s.clear_stage()


def part_fisher(s):
    nb = numbers()
    frows = [[r"\nu"] + [ftex(nu) for nu in WS.FISHER_NU] + [str(NU3)],
             [r"\frac{1}{1-\nu}"] + [ftex(v) for v in nb["fisher"]] + [f2(nb["fisher3"])]]
    ftab = D.tex_table(frows, [1.6, 1.5, 1.5, 2.0, 2.0], row_h=0.75, font_size=30)
    flab = VGroup(zh("高斯", font_size=22, color=GREY_B), zh("二維自迴避", font_size=22, color=GREY_B),
                  zh("三維（引用）", font_size=22, color=GREY_B))
    for lab, j in zip(flab, (1, 3, 4)):
        lab.next_to(D.table_cell(ftab, 1, j, 5), DOWN, buff=0.35)
    fblock = VGroup(mixed("Fisher：", "$C\\sim e^{-(r/N^\\nu)^{1/(1-\\nu)}}$", font_size=28, color=YELLOW), VGroup(ftab, flab))
    fblock.arrange(DOWN, buff=0.3)
    g1 = mixed("Ex. 8.3：總條數 ", "$\\propto(z_c-z)^{-\\gamma}$", "，", "$\\gamma=(2-\\eta)\\nu$", font_size=28)
    g2 = mixed("普通：", "$(2-0)\\cdot\\frac12=%s$" % ftex(nb["g_ord"]), "；二維自迴避 ",
               "$\\eta=\\frac{5}{24},\\ \\nu=\\frac34\\Rightarrow\\gamma=%s$" % ftex(nb["g_2d"]), font_size=28)
    x, p = nb["growth"]
    g3 = mixed("Ex. 8.4：條數 ", "$\\propto x^NN^p$", "，", "$x=1/z_c=\\mu\\approx%.3f$" % x, "，",
               "$p=\\gamma-1=%s=%.5f$" % (ftex(p), float(p)), font_size=28)
    gblock = VGroup(g1, g2, g3).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    place_in(fblock, -6.7, 6.7, 0.4, 3.1)
    place_in(gblock, -6.7, 6.7, -2.9, 0.0)
    assert_in_frame(fblock, gblock)
    assert_apart(fblock, gblock)
    assert_table(ftab, 5)
    for lab, j in zip(flab, (1, 3, 4)):
        assert_apart(lab, D.table_cell(ftab, 1, j, 5))
    assert_apart(flab[0], flab[1], 0.2)
    assert_apart(flab[1], flab[2], 0.2)
    assert_column(*gblock)
    s.say("指數 ", "$\\frac{1}{1-\\nu}$", "：Fisher 的結果", wait=0.3)
    s.play(FadeIn(fblock[0]), FadeIn(ftab[0][:5]), ShowCreation(ftab[1]), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(c) for c in ftab[0][5:]], lag_ratio=0.15), FadeIn(flab), run_time=1.2)
    s.say("$\\nu=\\frac12$", " 得 2：高斯；二維自迴避 ", "$\\nu=\\frac34$", " 得 4", wait=3.0)
    s.say("三維 ", "$\\nu\\approx0.588$", "（只引用）：2.43", wait=2.5)
    s.play(FadeIn(g1), run_time=0.8)
    s.say("Ex. 8.3：", "$g$", " 在 ", "$q=0$", " 的發散 ", "$\\gamma=(2-\\eta)\\nu$", wait=2.5)
    s.play(FadeIn(g2), run_time=0.8)
    s.say("普通漫步 ", "$\\gamma=1$", "；二維自迴避 ", "$\\frac{43}{32}$", "（已知值，只引用）", wait=3.0)
    s.play(FadeIn(g3), run_time=0.8)
    s.say("Ex. 8.4：條數 ", "$\\propto\\mu^NN^{\\gamma-1}$", "，二維 ", "$p=\\frac{11}{32}$", wait=3.0)
    s.clear_stage()


def part_ratio(s):
    nb = numbers()
    rows = [["N"] + [str(N) for N in WS.RATIO_N] + [r"\to"],
            [r"\hat\mu_N"] + [f4(nb["ratio"][N][1]) for N in WS.RATIO_N] + ["2.638"],
            [r"\hat\gamma_N"] + [f4(nb["ratio"][N][2]) for N in WS.RATIO_N] + [f"{float(WS.GAMMA_QUOTED):.5f}"]]
    tab = D.tex_table(rows, [1.4, 1.5, 1.5, 1.5, 1.5, 1.7], row_h=0.75, font_size=32)
    for i in (1, 2):
        D.table_cell(tab, i, 5, 6).set_color(YELLOW)
    head = mixed("比值法（第 7 章的 ", "$c_N$", "；", "$c_{15},c_{16}$", " 是文獻值）", font_size=28, color=GREY_B)
    ex85 = mixed("Ex. 8.5：小 ", "$r$", " 的條數 ", "$\\propto r^\\theta$", "，", "$\\theta=\\frac{\\gamma-1}{\\nu}$",
                 "（二維自迴避 ", "$%s$" % ftex(nb["theta"]), "）", font_size=28)
    col = VGroup(head, tab, ex85).arrange(DOWN, buff=0.5)
    place_in(col, -6.7, 6.7, -2.6, 2.9)
    assert_in_frame(col)
    assert_column(*col)
    assert_table(tab, 6)
    s.say("用第 7 章數出來的 ", "$c_N$", " 檢查：比值法", wait=0.3)
    s.play(FadeIn(head), FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:]], lag_ratio=0.08), run_time=1.5)
    s.say("$\\hat\\mu$", "：2.6592、2.6529、2.6483、2.6454 → 2.638", wait=3.0)
    s.say("$\\hat\\gamma$", "：1.3552 … 1.3483 → ", "$\\frac{43}{32}=1.34375$", wait=3.0)
    s.play(FadeIn(ex85), run_time=0.8)
    s.say("Ex. 8.5：小 ", "$r$", " 的條數 ", "$\\propto r^\\theta$", "，二維自迴避 ", "$\\theta=\\frac{11}{24}$", wait=3.0)
    s.say("$\\nu=\\frac12$", " 才是高斯；", "$\\nu\\ne\\frac12$", " 時 ", "$C\\sim e^{-(r/N^\\nu)^{1/(1-\\nu)}}$", wait=3.0)
    s.say("自迴避漫步的分布不是高斯", wait=2.5)


class W803SteepestDescent(ChuScene):
    def construct(self):
        self.show_header("W803　最陡下降：從標度形式讀出分布")
        check_math()
        part_contour(self)
        part_820(self)
        part_1d(self)
        part_fisher(self)
        part_ratio(self)
        self.wait(1.0)
