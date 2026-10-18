"""
W903 臨界指數與標度律（§9.2.2 後半–9.2.3，pp. 209–215，(9.13)–(9.27)，Table 9.1，Ex. 9.1–9.2）

問：一個磁鐵有六個臨界指數 α, β, γ, δ, ν, η——它們彼此獨立嗎？
  1. Table 9.1（平均場｜書上三維 n = 0）；筆誤框 1：α = 0.12 → 0.19（2 − 2β − γ = 2 − 3ν = 0.194）。
  2. 一維 Ising 鏈（K = 1、t = 0.7616）：G(r) = t^{|r|} 長條圖、ξ = 3.6719；χ 的部分和 5.7520、6.9696、7.3615 → e² = 7.3891；
     16 個自旋的環 7.2021；Σr²G/ΣG ÷ ξ² = 1.9877、1.9998、2.0000；筆誤框 2（ξ → ξ²）；χ/ξ = 2.0123、2.0002、2.0000 → η = 1。
  3. 標度律 (9.16)–(9.19)：平均場｜二維 Ising｜二維自迴避（n = 0）每格打勾；平均場的 hyperscaling 只在 d = 4。
  4. Widom 形式 (9.20)–(9.23)、Ex. 9.1–9.2：β = 2 − α − Δ、γ = Δ − β、δ = Δ/β；二維 Ising 的斜率 0.1250、1.7500、0.0667。
畫面上的數字都由 utils.walk_crit 算出，check_math() 核對。模組層級的畫圖函式：g_chart（G(r) 長條圖＋虛線 e^{−|r|/ξ}）。
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
from utils import walk_crit as WC
from utils import walk_draw as D

K1 = WC.K1D                            # 1.0
RMAX = 8                               # 長條圖 r = −8..8
CHI_R = (5, 10, 20)
FK = WC.FISHER1D_K                     # (1, 2, 3)
COLS = ("mf", "ising2d", "saw2d")
KEYS = ("alpha", "beta", "gamma", "delta", "nu", "eta")
GREEK = {"alpha": r"\alpha", "beta": r"\beta", "gamma": r"\gamma", "delta": r"\delta", "nu": r"\nu", "eta": r"\eta"}
GOOD = "#7CE38B"
BAR_C = BLUE_C


def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q, small=False):
    q = Fraction(q)
    cmd = r"\tfrac" if small else r"\frac"
    return str(q.numerator) if q.denominator == 1 else cmd + "{%d}{%d}" % (q.numerator, q.denominator)


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = {"b3": WC.book3d_alpha(), "t": math.tanh(K1),
              "G": {r: WC.ising1d_G(r, K1) for r in range(-RMAX, RMAX + 1)},
              "Gt": {r: WC.ising1d_transfer_G(r, K1) for r in range(-RMAX, RMAX + 1)},
              "xi": WC.ising1d_xi(K1), "chi": WC.ising1d_chi(K1),
              "chi_part": [WC.ising1d_chi_partial(K1, R) for R in CHI_R],
              "ring16": WC.ising1d_ring_closed(16, K1)[1],
              "ring16_ex": sum(WC.ising_exact((16, WC.ring_edges(16)), K1)["corr0"]),
              "xi2": [WC.ising1d_xi2_moment(K) / WC.ising1d_xi(K) ** 2 for K in FK],
              "xi2_sum": [WC.ising1d_xi2_moment_sum(K) / WC.ising1d_xi(K) ** 2 for K in FK],
              "chi_xi": [WC.ising1d_chi(K) / WC.ising1d_xi(K) for K in FK],
              "laws": {c: WC.laws(WC.EXPONENTS[c]) for c in COLS},
              "widom": WC.widom_exponents(*WC.WIDOM_CASE),
              "wslopes": WC.widom_slopes(*WC.WIDOM_CASE),
              "wslopes_fd": WC.widom_slopes(*WC.WIDOM_CASE, method="fd")}
    return nb


def law_ok(c):
    e, lw = WC.EXPONENTS[c], numbers()["laws"][c]
    return (lw["rushbrooke"] == 2, lw["fisher"] == e["gamma"], lw["hyper_alpha"] == e["alpha"],
            lw["hyper_delta"] == e["delta"])


def book(k):
    return f"{WC.BOOK_3D[k]:g}"


def check_math():
    nb = numbers()
    assert [book(k) for k in KEYS] == ["0.12", "0.308", "1.19", "4.865", "0.602", "0.023"]
    assert [f"{v:.3f}" for v in nb["b3"]] == ["0.194", "0.194", "1.926"]
    assert abs(nb["b3"][0] - (2 - 2 * WC.BOOK_3D["beta"] - WC.BOOK_3D["gamma"])) < 1e-12
    mf = WC.EXPONENTS["mf"]
    assert [frac_tex(mf[k]) for k in KEYS] == ["0", r"\frac{1}{2}", "1", "3", r"\frac{1}{2}", "0"]
    assert f4(nb["t"]) == "0.7616" and nb["G"][0] == 1
    assert [f4(nb["G"][r]) for r in range(1, 5)] == ["0.7616", "0.5800", "0.4417", "0.3364"]
    assert all(abs(nb["G"][r] - nb["Gt"][r]) < 1e-12 and nb["G"][r] == nb["G"][-r] for r in nb["G"])
    assert f4(nb["xi"]) == "3.6719" and abs(math.exp(-1 / nb["xi"]) - nb["t"]) < 1e-12
    assert f4(nb["chi"]) == "7.3891" and abs(nb["chi"] - math.e ** 2) < 1e-12
    assert [f4(v) for v in nb["chi_part"]] == ["5.7520", "6.9696", "7.3615"]
    assert f4(nb["ring16"]) == f4(nb["ring16_ex"]) == "7.2021"
    assert [f4(v) for v in nb["xi2"]] == ["1.9877", "1.9998", "2.0000"]
    assert all(abs(a - b) < 1e-4 for a, b in zip(nb["xi2"], nb["xi2_sum"]))
    assert [f4(v) for v in nb["chi_xi"]] == ["2.0123", "2.0002", "2.0000"]
    want = {"mf": ("0", "1/2", "1", "3", "1/2", "0"), "ising2d": ("0", "1/8", "7/4", "15", "1", "1/4"),
            "saw2d": ("1/2", "5/64", "43/32", "91/5", "3/4", "5/24")}
    for c in COLS:
        e = WC.EXPONENTS[c]
        assert tuple(str(e[k]) for k in KEYS) == want[c]
        assert all(law_ok(c)) and nb["laws"][c]["widom"] == e["gamma"]
        assert WC.from_nu_eta(e["nu"], e["eta"], e["d"]) == tuple(e[k] for k in ("alpha", "beta", "gamma", "delta"))
    assert [WC.EXPONENTS[c]["d"] for c in COLS] == [4, 2, 2]
    assert 2 - 3 * mf["nu"] == Fraction(1, 2) != mf["alpha"]
    b, g, dl = nb["widom"]
    assert (b, g, dl) == (Fraction(1, 8), Fraction(7, 4), Fraction(15)) and b * (dl - 1) == g
    a0, Dl = WC.WIDOM_CASE
    assert a0 == 0 and Dl == Fraction(15, 8) and b == 2 - a0 - Dl and g == Dl - b and dl == Dl / b and b + g == Dl
    assert [f4(v) for v in nb["wslopes"]] == ["0.1250", "1.7500", "0.0667"] and f4(1 / dl) == "0.0667"
    assert all(abs(x - y) < 1e-3 for x, y in zip(nb["wslopes"], nb["wslopes_fd"]))
    assert Fraction(1, 8) * 14 == Fraction(7, 4)


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


def row_cells(tab, i, ncols):
    return tab[0][i * ncols:(i + 1) * ncols]


def show_table(s, tab, ncols, extra=(), run_time=1.2):
    """先表頭與橫線，再逐列淡入。"""
    s.play(FadeIn(row_cells(tab, 0, ncols)), ShowCreation(tab[1]), *[FadeIn(m) for m in extra], run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][ncols:]], lag_ratio=0.06), run_time=run_time)


# ---------------------------------------------------------------- 畫圖元件
def g_chart(G, xi, center_x, base_y, ux=0.4, uy=2.3):
    """一維 Ising 的 G(r)（r = −8..8）長條圖＋虛線 e^{−|r|/ξ}（兩段，r = 0 有尖角）＋ r 刻度
    → VGroup(長條 VGroup(長條們, 空標籤), 虛線兩段, 刻度, 底線)。"""
    o = np.array([center_x, base_y, 0.0])
    ch = D.bar_chart(G, ux=ux, uy=uy, origin=o, color=BAR_C, width=0.7)
    fn = lambda r: math.exp(-abs(r) / xi)
    dash = VGroup(D.curve(fn, -RMAX - 0.3, 0, ux, uy, o, color=ORANGE, stroke_width=3, samples=60, dashed=True,
                          num_dashes=24),
                  D.curve(fn, 0, RMAX + 0.3, ux, uy, o, color=ORANGE, stroke_width=3, samples=60, dashed=True,
                          num_dashes=24))
    ticks = VGroup(*[tx(str(r), 22, GREY_B).move_to(o + r * ux * RIGHT + 0.25 * DOWN) for r in (-8, -4, 0, 4, 8)])
    base = Line(o + (RMAX + 0.6) * ux * LEFT, o + (RMAX + 0.6) * ux * RIGHT).set_stroke(GREY_B, 1.5)
    return VGroup(ch, dash, ticks, base)


def recolored(bars, R):
    """長條們的複本：|r| ≤ R 的塗黃（部分和 Σ_{|r|≤R} G 用到的長條）。拿來 Transform 整組，不要個別動子物件。"""
    tgt = bars.copy()
    for b, b0 in zip(tgt, bars):
        if abs(b0.bar_x) <= R:
            b.set_fill(YELLOW, 0.85).set_stroke(YELLOW, 1.5)
    return tgt


# ---------------------------------------------------------------- 各段
def part_table91(s):
    nb = numbers()
    mf = WC.EXPONENTS["mf"]
    fs = 28
    hd = [zh(w, font_size=26, color=YELLOW) for w in ("量", "行為", "指數", "平均場")]
    hd.append(mixed("三維 ", "$n=0$", "（書上）", font_size=26, color=YELLOW))
    quant = [zh("比熱", font_size=fs), "M", r"\chi", "M(T_c,h)", r"\xi", r"C(T_c,\rho)"]
    beh = [r"\propto|T-T_c|^{-\alpha}", r"\propto|T_c-T|^{\beta}", r"\propto|T-T_c|^{-\gamma}",
           r"\propto|h|^{1/\delta}", r"\propto|T-T_c|^{-\nu}", r"\propto\rho^{-d+2-\eta}"]
    rows = [hd] + [[q, b, GREEK[k], frac_tex(mf[k]), book(k)] for q, b, k in zip(quant, beh, KEYS)]
    tab = D.tex_table(rows, [2.2, 3.4, 1.2, 1.6, 2.8], row_h=0.58, font_size=fs)
    place_in(tab, -6.7, 6.7, -1.0, 3.1)
    a_cell = D.table_cell(tab, 1, 4, 5)
    mark = SurroundingRectangle(a_cell, buff=0.08).set_stroke(RED, 2.5)
    box = D.erratum_box("p. 210，Table 9.1", r"\alpha=\alpha'=0.12", r"\alpha\approx0.19", width=5.2, font_size=32)
    b3 = nb["b3"]
    c1 = mixed("(9.16)　", r"$2-2\beta-\gamma=2-2(%s)-%s=%s$" % (book("beta"), book("gamma"), f"{b3[0]:.3f}"),
               font_size=26)
    c2 = mixed("(9.18)　", r"$2-3\nu=2-3(%s)=%s$" % (book("nu"), f"{b3[1]:.3f}"), font_size=26)
    calc = VGroup(c1, c2).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    place_in(box, -6.7, -1.4, -2.9, -1.2)
    place_in(calc, -1.1, 6.7, -2.9, -1.2)
    assert_in_frame(tab, mark, box, calc)
    assert_apart(tab, box)
    assert_apart(tab, calc)
    assert_apart(box, calc)
    assert_table(tab, 5)
    s.say("問：一個磁鐵有六個臨界指數 ", "$\\alpha,\\beta,\\gamma,\\delta,\\nu,\\eta$")
    s.say("它們彼此獨立嗎？", wait=1.5)
    show_table(s, tab, 5, run_time=2.5)
    s.say("Table 9.1：靠近 ", "$T_c$", " 時六個量的冪次行為", wait=2.5)
    s.say("平均場：", "$\\alpha=0$", "、", "$\\beta=\\frac12$", "、", "$\\gamma=1$", "、", "$\\delta=3$", "、",
          "$\\nu=\\frac12$", "、", "$\\eta=0$", wait=2.5)
    s.say("右欄：書上給的三維 ", "$n=0$", "（自迴避漫步）的數值", wait=1.5)
    s.play(ShowCreation(mark), run_time=0.6)
    s.say("同一欄的 ", "$\\beta,\\gamma$", " 代 (9.16)：", "$2-2\\beta-\\gamma=0.194$")
    s.play(FadeIn(c1), run_time=0.8)
    s.wait(1.5)
    s.say("$\\nu$", " 代 (9.18)：", "$2-3\\nu$", " 也是 0.194", wait=0.3)
    s.play(FadeIn(c2), run_time=0.8)
    s.wait(1.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("所以書上的 ", "$\\alpha=0.12$", " 應該是 ", "$\\alpha\\approx0.19$", wait=3.0)
    s.clear_stage()


def part_chain(s):
    nb = numbers()
    G, xi, t = nb["G"], nb["xi"], nb["t"]
    title = mixed("一維 Ising 鏈：", "$K=1$", "，", f"$t=\\tanh K={f4(t)}$", font_size=28)
    chart = g_chart(G, xi, -3.15, -0.75)
    gtab = D.tex_table([["r", "0", "1", "2", "3", "4", r"\cdots"],
                        ["G(r)", "1"] + [f4(G[r]) for r in range(1, 5)] + [r"\cdots"]],
                       [1.2, 0.7, 1.1, 1.1, 1.1, 1.1, 0.6], row_h=0.6, font_size=26)
    place_in(title, -6.7, 0.4, 2.55, 3.1)
    place_in(gtab, -6.7, 0.4, -2.9, -1.45)
    r1 = tx(r"G(r)=\langle s_0s_r\rangle=t^{|r|}", 30)
    r2 = tx(r"\xi=-\frac{1}{\ln t}=%s" % f4(xi), 30)
    r3 = VGroup(DashedLine(ORIGIN, 0.7 * RIGHT).set_stroke(ORANGE, 3), tx(r"e^{-|r|/\xi}", 28, ORANGE)
                ).arrange(RIGHT, buff=0.15)
    r4 = mixed("(9.24)–(9.27)　", "$\\chi=\\sum_rG(r)$", font_size=28)
    ptab = D.tex_table([["R", "5", "10", "20", r"\infty"],
                        [r"\sum_{|r|\le R}G"] + [f4(v) for v in nb["chi_part"]] + [f4(nb["chi"])]],
                       [1.6, 1.05, 1.05, 1.05, 1.05], row_h=0.7, font_size=26)
    D.table_cell(ptab, 1, 4, 5).set_color(YELLOW)
    r6 = tx(r"\frac{1+t}{1-t}=e^{2K}=%s" % f4(nb["chi"]), 30, YELLOW)
    r7 = mixed("16 個自旋的環：", f"$\\chi_{{16}}={f4(nb['ring16'])}$", font_size=26)
    right = VGroup(r1, r2, r3, r4, ptab, r6, r7).arrange(DOWN, buff=0.2)
    place_in(right, 0.8, 6.7, -2.9, 3.1)
    assert_in_frame(title, chart, gtab, right)
    assert_column(title, chart, gtab)
    assert_column(r1, r2, r3, r4, ptab, r6, r7)
    assert_apart(chart, right)
    assert_apart(gtab, right)
    assert_table(gtab, 7)
    assert_table(ptab, 5)
    bars = chart[0][0]
    s.say("一維 Ising 鏈可以完全算出來：", "$K=1$")
    s.play(FadeIn(title), ShowCreation(chart[3]), FadeIn(chart[2]), run_time=0.8)
    s.play(GrowFromEdge(bars, DOWN), FadeIn(r1), run_time=1.2)
    s.say("關聯函數 ", "$G(r)=\\langle s_0s_r\\rangle=t^{|r|}$", wait=2.0)
    show_table(s, gtab, 7, run_time=1.0)
    s.say("$G(r)$", "：1、" + "、".join(f4(G[r]) for r in range(1, 4)) + "……", wait=2.5)
    s.play(ShowCreation(chart[1][0]), ShowCreation(chart[1][1]), FadeIn(r2), FadeIn(r3), run_time=1.2)
    s.say("指數衰減，關聯長度 ", f"$\\xi=-1/\\ln t={f4(xi)}$", wait=2.5)
    s.play(FadeIn(r4), run_time=0.8)
    s.say("(9.24)–(9.27)：磁化率 ", "$\\chi=\\sum_rG(r)$", "，一根一根加", wait=1.0)
    s.play(Transform(bars, recolored(bars, CHI_R[0])), FadeIn(row_cells(ptab, 0, 5)), ShowCreation(ptab[1]),
           FadeIn(ptab[0][5:7]), run_time=1.0)
    s.play(Transform(bars, recolored(bars, RMAX)), FadeIn(ptab[0][7:9]), run_time=1.0)
    s.say("部分和：" + "、".join(f4(v) for v in nb["chi_part"]), wait=2.5)
    s.play(FadeIn(ptab[0][9]), FadeIn(r6), run_time=0.8)
    s.say("全部加起來：", "$\\frac{1+t}{1-t}=e^{2K}$", f"＝{f4(nb['chi'])}", wait=2.0)
    s.play(FadeIn(r7), run_time=0.8)
    s.say(f"16 個自旋的環逐一列舉：{f4(nb['ring16'])}（有限大小）", wait=2.0)
    s.clear_stage()

    head = zh("分子多一個 r²：這個比值是什麼？", font_size=28, color=GREY_B)
    tab = D.tex_table([["K"] + [f"{K:g}" for K in FK],
                       [r"\frac{\sum_r r^2G(r)}{\sum_r G(r)}\Big/\xi^2"] + [f4(v) for v in nb["xi2"]],
                       [r"\chi/\xi"] + [f4(v) for v in nb["chi_xi"]]],
                      [3.6, 1.5, 1.5, 1.5], row_h=1.0, font_size=30)
    for j in (1, 2, 3):
        D.table_cell(tab, 1, j, 4).set_color(YELLOW)
    top = VGroup(head, tab).arrange(DOWN, buff=0.25)
    box = D.erratum_box("p. 210，Table 9.1 第二個關聯長度", r"\xi=\frac{\int C\rho^2d^d\rho}{\int C\,d^d\rho}",
                        r"\xi^2=\frac{\int C\rho^2d^d\rho}{\int C\,d^d\rho}", width=7.0, font_size=32)
    concl = mixed("$\\chi/\\xi\\to2$", "：", "$\\chi\\propto\\xi^{2-\\eta}$", "，一維 ", "$\\eta=1$", font_size=30,
                  color=YELLOW)
    place_in(top, -6.7, 6.7, 0.15, 3.1)
    place_in(box, -6.7, 6.7, -1.95, -0.05)
    place_in(concl, -6.7, 6.7, -2.9, -2.2)
    assert_in_frame(top, box, concl)
    assert_column(top, box, concl)
    assert_table(tab, 4)
    s.say("$\\sum r^2G/\\sum G$", " 除以 ", "$\\xi^2$", "：", "$K=1,2,3$")
    show_table(s, tab, 4, extra=(head,), run_time=1.0)
    s.say("、".join(f4(v) for v in nb["xi2"]) + "：是常數", wait=2.0)
    s.say("所以這個比值是長度的平方", wait=1.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上第二個關聯長度寫成 ", "$\\xi=$", "；應該是 ", "$\\xi^2=$", wait=2.5)
    s.say("比值是長度的平方（Table 9.1 第一列寫對了）", wait=2.5)
    s.play(FadeIn(concl), run_time=0.8)
    s.say("$\\chi/\\xi$", "：" + "、".join(f4(v) for v in nb["chi_xi"]), wait=2.0)
    s.say("$\\chi\\propto\\xi^{2-\\eta}$", "：一維 ", "$\\eta=1$", wait=2.5)
    s.clear_stage()


def part_laws(s):
    nb = numbers()
    hd = [zh("指數", font_size=26, color=YELLOW), zh("平均場", font_size=26, color=YELLOW),
          zh("二維 Ising", font_size=26, color=YELLOW), mixed("二維 ", "$n=0$", font_size=26, color=YELLOW)]
    rows = [hd] + [[GREEK[k]] + [frac_tex(WC.EXPONENTS[c][k], small=True) for c in COLS] for k in KEYS]
    etab = D.tex_table(rows, [1.0, 1.4, 1.8, 1.8], row_h=0.64, font_size=28)
    laws = [r"(9.16)\ \ \alpha+2\beta+\gamma=2", r"(9.17)\ \ (2-\eta)\nu=\gamma", r"(9.18)\ \ 2-d\nu=\alpha",
            r"(9.19)\ \ \frac{d+2-\eta}{d-2+\eta}=\delta"]
    marks = {c: law_ok(c) for c in COLS}
    lrows = [[zh("標度律", font_size=26, color=YELLOW), zh("平均場", font_size=24, color=YELLOW),
              zh("Ising", font_size=24, color=YELLOW), mixed("$n=0$", font_size=24, color=YELLOW)]]
    for i, lw in enumerate(laws):
        lrows.append([lw] + [r"\checkmark" if marks[c][i] else r"\times" for c in COLS])
    ltab = D.tex_table(lrows, [3.7, 1.0, 0.9, 0.9], row_h=0.86, font_size=28)
    for i in range(1, 5):
        for j in (1, 2, 3):
            D.table_cell(ltab, i, j, 4).set_color(GOOD)
    dnote = mixed("(9.18)(9.19) 用 ", "$d=4,\\ 2,\\ 2$", font_size=24, color=GREY_B)
    right = VGroup(ltab, dnote).arrange(DOWN, buff=0.25)
    place_in(etab, -6.7, -0.55, -1.7, 3.1)
    place_in(right, -0.25, 6.7, -1.7, 3.1)
    mfa = mixed("平均場在 ", "$d=3$", "：", "$2-3\\cdot\\frac12=\\frac12\\neq0=\\alpha$", font_size=26)
    b3n = mixed("書上三維欄：", f"$\\alpha+2\\beta+\\gamma={nb['b3'][2]:.3f}$", font_size=26)
    low = VGroup(mfa, b3n).arrange(RIGHT, buff=0.8)
    place_in(low, -6.7, 6.7, -2.9, -2.0)
    assert_in_frame(etab, right, low)
    assert_apart(etab, right)
    assert_apart(etab, low)
    assert_apart(right, low)
    assert_table(etab, 4)
    assert_table(ltab, 4)
    s.say("三個模型的指數：平均場、二維 Ising（精確）、二維自迴避")
    show_table(s, etab, 4, run_time=2.0)
    s.say("二維自迴避就是 ", "$n=0$", "：", "$\\nu=\\frac34$", "、", "$\\gamma=\\frac{43}{32}$", wait=2.5)
    s.play(FadeIn(row_cells(ltab, 0, 4)), ShowCreation(ltab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in ltab[0][4:12]], lag_ratio=0.1), run_time=1.2)
    s.say("(9.16)、(9.17)：三欄都成立", wait=2.5)
    s.play(LaggedStart(*[FadeIn(c) for c in ltab[0][12:]], lag_ratio=0.1), FadeIn(dnote), run_time=1.2)
    s.say("(9.18)(9.19) 含維度 ", "$d$", "：平均場要用 ", "$d=4$", wait=2.5)
    s.play(FadeIn(mfa), run_time=0.8)
    s.say("平均場的 hyperscaling 只在 ", "$d=4$", " 成立", wait=2.0)
    s.play(FadeIn(b3n), run_time=0.8)
    s.say("書上三維欄 ", f"$\\alpha+2\\beta+\\gamma={nb['b3'][2]:.3f}$", "：框 1 的由來", wait=2.5)
    s.clear_stage()


def part_widom(s):
    nb = numbers()
    b, g, dl = nb["widom"]
    e20 = mixed("(9.20)–(9.23)　", "$G(r,h)=|r|^{2-\\alpha}f\\Big(\\frac{h}{|r|^{\\Delta}},\\frac{r}{|r|}\\Big)$",
                font_size=32)
    note = zh("（這裡 r 是約化溫度、G 是自由能）", font_size=26, color=GREY_B)
    d1 = mixed("Ex. 9.1：", "$M=-\\frac{\\partial G}{\\partial h}\\propto|r|^{2-\\alpha-\\Delta}$", "　",
               "$\\Rightarrow\\ \\beta=2-\\alpha-\\Delta$", font_size=30)
    d2 = mixed("$\\chi=\\frac{\\partial M}{\\partial h}\\propto|r|^{2-\\alpha-2\\Delta}$", "　",
               "$\\Rightarrow\\ \\gamma=\\Delta-\\beta$", font_size=30)
    d3 = mixed("Ex. 9.2：", "$M(0,h)\\propto h^{\\beta/\\Delta}$", "　", "$\\Rightarrow\\ \\delta=\\Delta/\\beta$",
               font_size=30)
    d4 = mixed("合起來：", "$\\beta(\\delta-1)=\\Delta-\\beta=\\gamma$", "（Widom）", font_size=30, color=YELLOW)
    derv = VGroup(d1, d2, d3, d4).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
    top = VGroup(VGroup(e20, note).arrange(DOWN, buff=0.15), derv).arrange(DOWN, buff=0.5)
    place_in(top, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(top)
    assert_column(e20, note, d1, d2, d3, d4)
    s.say("Widom：自由能是齊次函數 (9.20)–(9.23)")
    s.play(FadeIn(e20), FadeIn(note), run_time=1.0)
    s.say("這裡 ", "$r$", " 是約化溫度（不是距離）、", "$G$", " 是自由能", wait=1.5)
    s.play(FadeIn(d1), run_time=0.8)
    s.say("Ex. 9.1：對 ", "$h$", " 微分一次得 ", "$M$", "，", "$\\beta=2-\\alpha-\\Delta$", wait=2.5)
    s.play(FadeIn(d2), run_time=0.8)
    s.say("再微分一次得 ", "$\\chi$", "：", "$\\gamma=\\Delta-\\beta$", wait=2.0)
    s.play(FadeIn(d3), run_time=0.8)
    s.say("Ex. 9.2：", "$M\\propto h^{1/\\delta}$", "，", "$\\delta=\\Delta/\\beta$", wait=2.0)
    s.play(FadeIn(d4), run_time=0.8)
    s.say("合起來：", "$\\beta(\\delta-1)=\\gamma$", "（Widom 標度律）", wait=2.5)
    s.clear_stage()

    head = mixed("二維 Ising：", "$\\alpha=0$", "、", "$\\Delta=\\frac{15}{8}$", font_size=32, color=YELLOW)
    wt = D.tex_table([[VectorizedPoint(), r"\beta", r"\gamma", r"1/\delta"],
                      [zh("公式", font_size=28), frac_tex(b), frac_tex(g), frac_tex(1 / dl)],
                      [zh("一個具體的 f 量出的斜率", font_size=28)] + [f4(v) for v in nb["wslopes"]]],
                     [4.4, 1.7, 1.7, 1.7], row_h=0.9, font_size=34)
    for j in (1, 2, 3):
        D.table_cell(wt, 2, j, 4).set_color(YELLOW)
    chk = VGroup(tx(r"\beta(\delta-1)=\gamma", 38), tx(r"\frac18\times14=\frac74", 38, GOOD)).arrange(RIGHT, buff=1.2)
    place_in(head, -6.7, 6.7, 2.4, 3.1)
    place_in(wt, -6.7, 6.7, -0.9, 2.1)
    place_in(chk, -6.7, 6.7, -2.9, -1.3)
    assert_in_frame(head, wt, chk)
    assert_column(head, wt, chk)
    assert_table(wt, 4)
    s.play(FadeIn(head), run_time=0.8)
    s.say("二維 Ising：", "$\\alpha=0$", "、", "$\\Delta=\\beta+\\gamma=\\frac{15}{8}$", wait=1.0)
    show_table(s, wt, 4, run_time=1.2)
    s.say("二維 Ising：", "$\\beta=\\frac18$", "、", "$\\gamma=\\frac74$", "、", "$\\delta=15$", wait=2.0)
    s.say("數值斜率：" + "、".join(f4(v) for v in nb["wslopes"]), wait=2.5)
    s.play(FadeIn(chk), run_time=0.8)
    s.say("$\\frac18\\times14=\\frac74$", "：Widom 的關係成立", wait=2.0)
    s.clear_stage()

    head = mixed("知道 ", "$\\nu$", "、", "$\\eta$", " 和維度 ", "$d$", "：", font_size=34)
    fm = VGroup(tx(r"\alpha=2-d\nu", 36), tx(r"\gamma=(2-\eta)\nu", 36),
                tx(r"\beta=\frac{\nu(d-2+\eta)}{2}", 36), tx(r"\delta=\frac{d+2-\eta}{d-2+\eta}", 36))
    fm.arrange_in_grid(2, 2, h_buff=1.2, v_buff=0.5)
    final = VGroup(head, fm).arrange(DOWN, buff=0.5)
    place_in(final, -6.7, 6.7, -2.6, 2.8)
    assert_in_frame(final)
    assert_apart(head, fm)
    s.play(FadeIn(head), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(m) for m in fm], lag_ratio=0.3), run_time=1.5)
    s.say("六個指數只有兩個獨立", wait=2.0)
    s.say("知道 ", "$\\nu$", " 和 ", "$\\eta$", "，其餘都由標度律決定", wait=3.0)


class W903ScalingLaws(ChuScene):
    def construct(self):
        self.show_header("W903　臨界指數與標度律")
        check_math()
        part_table91(self)
        part_chain(self)
        part_laws(self)
        part_widom(self)
        self.wait(1.0)
