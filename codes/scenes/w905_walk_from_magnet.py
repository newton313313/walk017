"""
W905 從磁鐵讀出漫步：⟨R²⟩ ∝ N^{2ν}（§9.4，pp. 218–219，(9.36)–(9.41)）

問：只知道磁鐵的 ξ ∝ (z_c − z)^{−ν}、χ ∝ (z_c − z)^{−γ}，能不能推出漫步的 ⟨R²⟩ ∝ N^{2ν}？
  1. 問題：磁鐵（n = 0）的兩個冪次 → 漫步的 ⟨R²⟩_N。
  2. (9.36) ξ² = Σr²G/ΣG；(9.37) 分子 Σ c_N⟨R²⟩_N z^N、分母 χ = Σ c_N z^N（z 每步的逸度，z_c = 1/μ = 0.3791）；
     c_N、c_N⟨R²⟩_N（N = 1..6）；部分和 ξ² = 0.6347、1.2799、2.4593（z = 0.1、0.15、0.2）。
  3. [z^N](z_c − z)^{−s} ≈ N^{s−1} z_c^{−N−s}/Γ(s)：精確除以漸近 → 1（s = 43/32、91/32）；筆誤框 1（N^{δ−1}）。
  4. (9.38)–(9.41)：係數相除 → ⟨R²⟩_N = K₁Γ(γ)N^{2ν}/(z_c^{2ν}Γ(γ+2ν))；比值 1.2476、1.0240、1.0024；筆誤框 2（z_c^γ）。
  5. 第 7 章的數據：2ν_eff → 3/2、⟨R²⟩_N/N^{3/2} → 常數（N = 1..14 的點圖）；收束。
畫面上的數字都由 utils.walk_crit（2ν_eff 由 utils.walk_saw）算出，check_math() 核對。
模組層級的畫圖函式：xi_plot／xi_curve（ξ² 的部分和）、ratio_axes／ratio_curve（係數比值對 log N）、amp_axes／amp_dots。
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
from utils import walk_saw as SW
from utils import walk_draw as D

NU_2D = Fraction(3, 4)
MAG_COLOR = BLUE_B
WALK_COLOR = YELLOW
S_COLORS = (TEAL, ORANGE)             # s = γ、γ + 2ν
AMP_N = tuple(range(1, 15))           # 點圖：N = 1..14（walk_saw 精確列舉到 14）


def f4(x):
    return f"{float(x):.4f}"


def sci(x):
    """1.4e-17 → 1.4\\times10^{-17}"""
    m, e = f"{float(x):.1e}".split("e")
    return r"%s\times10^{%d}" % (m, int(e))


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
        nb = {"table": WC.saw_table(),
              "xi2": [WC.saw_gf(z)[2] for z in WC.SAW_Z],
              "coef": [[WC.coef_ratio(N, s) for N in WC.COEF_N] for s in (WC.S_GAMMA, WC.S_GAMMA2NU)],
              "book": WC.coef_ratio_book(10, WC.S_GAMMA, Fraction(91, 5)),
              "r2": [WC.r2_ratio(N, WC.S_GAMMA, NU_2D, WC.ZC_2D) for N in WC.COEF_N],
              "zcp": WC.zc_pow(),
              "nu": [SW.two_nu_eff(N) for N in WC.SAW_N],
              "amp": [WC.saw_amp(N) for N in WC.SAW_N],
              "amp_all": [WC.saw_amp(N) for N in AMP_N]}
    return nb


def check_math():
    nb = numbers()
    assert f4(WC.ZC_2D) == "0.3791" and WC.SAW_Z == (0.1, 0.15, 0.2) and WC.COEF_N == (10, 100, 1000)
    assert WC.SAW_N == (10, 12, 14) and WC.saw_gf.__defaults__ == (14,)                  # 部分和 N ≤ 14
    assert [c for _, c, _ in nb["table"]] == [4, 12, 36, 100, 284, 780]
    assert [r for _, _, r in nb["table"]] == [4, 32, 164, 704, 2716, 9808]
    assert [n for n, _, _ in nb["table"]] == [1, 2, 3, 4, 5, 6]
    assert [f4(x) for x in nb["xi2"]] == ["0.6347", "1.2799", "2.4593"]
    assert nb["xi2"][0] < nb["xi2"][1] < nb["xi2"][2]
    assert WC.S_GAMMA == Fraction(43, 32) and WC.S_GAMMA2NU == WC.S_GAMMA + 2 * NU_2D == Fraction(91, 32)
    assert [f4(x) for x in nb["coef"][0]] == ["1.0227", "1.0023", "1.0002"]
    assert [f4(x) for x in nb["coef"][1]] == ["1.2760", "1.0264", "1.0026"]
    assert sci(nb["book"]) == r"1.4\times10^{-17}"
    assert [f4(x) for x in nb["r2"]] == ["1.2476", "1.0240", "1.0024"]
    zp, zb = nb["zcp"]
    assert f4(zp) == "4.2850" and f4(zb) == "3.6823"
    assert abs(zp - WC.ZC_2D ** (-2 * float(NU_2D))) < 1e-12 and abs(zb - WC.ZC_2D ** (-float(WC.S_GAMMA))) < 1e-12
    assert [f4(x) for x in nb["nu"]] == ["1.4442", "1.4505", "1.4556"]
    assert [f4(x) for x in nb["amp"]] == ["0.8299", "0.8224", "0.8168"]
    assert nb["amp_all"][9::2] == nb["amp"] and nb["amp_all"][0] == 1                    # N = 10、12、14；N = 1：4/4
    assert all(0.8 < a <= 1.0 for a in nb["amp_all"])                                    # 點圖的縱軸 0.8..1
    assert AMP_N[0] == 1 and AMP_N[-1] == 14 == max(WC.SAW_N)                           # 字幕「N = 1..14」
    assert 2 * NU_2D == Fraction(3, 2)
    # 比值都 → 1（N 越大越接近）
    assert all(abs(v[0] - 1) > abs(v[1] - 1) > abs(v[2] - 1) for v in nb["coef"] + [nb["r2"]])


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


def framed(mob, color, buff=0.3):
    return VGroup(SurroundingRectangle(mob, buff=buff).set_stroke(color, 2), mob)


# ---------------------------------------------------------------- 畫圖元件
def xi_plot():
    """ξ²（N ≤ 14 的部分和）對 z 的座標：z ∈ [0, 0.4]、ξ² ∈ [0, 3.3]。"""
    return D.plot_axes((0, 0.4, 0.1), (0, 3.3, 1), width=4.4, height=2.3, x_ticks=(0, 0.1, 0.2, 0.3),
                       y_ticks=(1, 2, 3), x_label="z", y_label=r"\xi^2")


def xi_curve(ax, z1=0.2, n=31):
    return VMobject().set_points_smoothly([ax.c2p(z, WC.saw_gf(z)[2]) for z in np.linspace(0, z1, n)])


def ratio_axes():
    """橫軸 log₁₀N − 1 ∈ [0, 2.2]（N = 10..1000）、縱軸 比值 − 1 ∈ [0, 0.32]：ManimGL 的軸畫在 0，所以平移。"""
    return D.plot_axes((0, 2.2, 1), (0, 0.32, 0.1), width=4.4, height=2.4, x_ticks=(0, 1, 2), y_ticks=(0, 0.1, 0.2, 0.3),
                       x_tick_labels=["10", "100", "1000"], y_tick_labels=["1", "1.1", "1.2", "1.3"], x_label="N")


def ratio_curve(ax, s, color):
    """精確除以漸近 coef_ratio(N, s)，N 取 10..1000 的對數等距整數（40 點）＋ N = 10、100、1000 三個點。"""
    Ns = sorted(set(int(round(v)) for v in np.logspace(1, 3, 40)))
    pts = [ax.c2p(math.log10(N) - 1, WC.coef_ratio(N, s) - 1) for N in Ns]
    curve = VMobject().set_points_smoothly(pts).set_stroke(color, 3)
    dots = VGroup(*[Dot(ax.c2p(math.log10(N) - 1, WC.coef_ratio(N, s) - 1), radius=0.06).set_color(color)
                    for N in WC.COEF_N])
    return VGroup(curve, dots)


def amp_axes():
    """⟨R²⟩_N/N^{3/2} 對 N：N ∈ [0, 15]、縱軸 0.8..1.02 平移成 0..0.22（ManimGL 的軸畫在 0）。"""
    return D.plot_axes((0, 15, 1), (0, 0.22, 0.1), width=4.6, height=2.4, x_ticks=(0, 5, 10), y_ticks=(0, 0.1, 0.2),
                       y_tick_labels=["0.8", "0.9", "1"], x_label="N")


def amp_dots(ax, vals):
    """N = 1..14 的點；N = 10、12、14（卡上的三個）黃色，其餘灰色。"""
    return VGroup(*[Dot(ax.c2p(N, a - 0.8), radius=0.06).set_color(YELLOW if N in WC.SAW_N else GREY_B)
                    for N, a in zip(AMP_N, vals)])


# ---------------------------------------------------------------- 各段
def part_question(s):
    mag = VGroup(mixed("磁鐵（", "$n=0$", "）", font_size=28, color=MAG_COLOR),
                 tx(r"\xi\propto(z_c-z)^{-\nu}", 40), tx(r"\chi\propto(z_c-z)^{-\gamma}", 40)).arrange(DOWN, buff=0.3)
    walk = VGroup(zh("自迴避漫步", font_size=28, color=WALK_COLOR), tx(r"\langle R^2\rangle_N\propto N^{2\nu}", 40))
    walk.arrange(DOWN, buff=0.3)
    bm, bw = framed(mag, MAG_COLOR), framed(walk, WALK_COLOR)
    pair = VGroup(bm, bw).arrange(RIGHT, buff=2.6)
    arrow = Arrow(bm.get_right(), bw.get_left(), buff=0.25)
    qm = zh("？", font_size=48, color=WHITE).next_to(arrow, UP, buff=0.1)
    dic = mixed("$N$", " 步", "$\\leftrightarrow$", "$\\frac{1}{z_c-z}$", "：離臨界點多遠", font_size=30, color=GREY_A)
    top = VGroup(VGroup(pair, arrow, qm), dic).arrange(DOWN, buff=0.7)
    place_in(top, -6.7, 6.7, -2.6, 3.0)
    assert_in_frame(top)
    assert_apart(bm, bw)
    assert_apart(pair, dic)
    s.play(FadeIn(bm), run_time=0.8)
    s.say("問：只知道磁鐵的 ", "$\\xi\\propto(z_c-z)^{-\\nu}$", "、", "$\\chi\\propto(z_c-z)^{-\\gamma}$", wait=2.0)
    s.play(GrowArrow(arrow), FadeIn(qm), FadeIn(bw), run_time=0.8)
    s.say("能不能推出漫步的 ", "$\\langle R^2\\rangle\\propto N^{2\\nu}$", "？", wait=2.0)
    s.play(FadeIn(dic), run_time=0.8)
    s.say("W901 的對照：步數 ", "$N$", " 對應 ", "$1/(z_c-z)$", wait=2.0)
    s.clear_stage()


def part_gf(s):
    nb = numbers()
    zc = WC.ZC_2D
    e36 = mixed("(9.36) ", "$\\xi^2=\\frac{\\sum_rr^2G(r)}{\\sum_rG(r)}$", font_size=28)
    e37 = mixed("(9.37) ", "$\\xi^2=\\frac{\\sum_Nc_N\\langle R^2\\rangle_Nz^N}{\\sum_Nc_Nz^N}$", "，",
                "$\\chi=\\sum_Nc_Nz^N$", font_size=28)
    top = VGroup(e36, e37).arrange(RIGHT, buff=0.8)
    rows = [["N"] + [str(n) for n, _, _ in nb["table"]],
            ["c_N"] + [str(c) for _, c, _ in nb["table"]],
            [r"c_N\langle R^2\rangle_N"] + [str(r) for _, _, r in nb["table"]]]
    tab = D.tex_table(rows, [2.4] + [1.25] * 6, row_h=0.6, font_size=30)
    zline = mixed("$z$", "：每一步的逸度，", "$z_c=1/\\mu=%s$" % f4(zc), font_size=26)
    sums = mixed("$N\\le14$", " 的部分和", font_size=24, color=GREY_A)
    ax = xi_plot()
    plot = VGroup(ax, sums).arrange(DOWN, buff=0.15)
    place_in(top, -6.7, 6.7, 2.1, 3.1)
    place_in(tab, -6.7, 6.7, 0.05, 1.85)
    place_in(zline, -6.7, -0.4, -1.8, -0.6)
    place_in(plot, -0.1, 6.7, -2.9, -0.2)
    zc_line = DashedLine(ax.c2p(zc, 0), ax.c2p(zc, 3), dash_length=0.08).set_stroke(RED, 2)
    zc_lab = tx("z_c", 24, RED).next_to(ax.c2p(zc, 3), UP, buff=0.06)
    curve = xi_curve(ax).set_stroke(MAG_COLOR, 3)
    pts = VGroup(*[Dot(ax.c2p(z, v), radius=0.06).set_color(YELLOW) for z, v in zip(WC.SAW_Z, nb["xi2"])])
    labs = VGroup(*[tx(f4(v), 20, YELLOW).next_to(ax.c2p(z, v), UL, buff=0.1) for z, v in zip(WC.SAW_Z, nb["xi2"])])
    assert_in_frame(top, tab, zline, plot, zc_lab, labs)
    assert_column(top, tab, plot)
    assert_apart(tab, zc_lab)
    assert_apart(tab, zline)
    assert_apart(zline, plot)
    assert_apart(e36, e37)
    assert_table(tab, 7)
    for a, b in zip(labs, labs[1:]):
        assert_apart(a, b, 0.02)
    cs = ",\\ ".join(str(c) for _, c, _ in nb["table"])
    rs = "、".join(str(r) for _, _, r in nb["table"])
    s.play(FadeIn(e36), run_time=0.8)
    s.say("(9.36)：關聯長度的平方＝", "$r^2$", " 對 ", "$G(r)$", " 的平均", wait=2.0)
    s.play(FadeIn(e37), run_time=0.8)
    s.say("(9.37)：把 ", "$G(r)$", " 寫成所有自迴避漫步的和", wait=2.0)
    s.say("$c_N$", "：", "$N$", " 步自迴避漫步的條數（書上 §9.4 叫 ", "$\\Gamma_N$", "）", wait=2.2)
    s.play(FadeIn(tab[0][:7]), ShowCreation(tab[1]), FadeIn(tab[0][7:14]), run_time=0.8)
    s.say("$c_N=%s$" % cs, "：分母就是 ", "$\\chi$", wait=2.0)
    s.play(FadeIn(tab[0][14:]), run_time=0.6)
    s.say("分子多乘一個 ", "$\\langle R^2\\rangle_N$", "：", rs, wait=2.2)
    s.play(FadeIn(zline), run_time=0.6)
    s.say("$z$", " 是每一步的逸度（權重）：", "$z_c=1/\\mu=%s$" % f4(zc), wait=2.0)
    s.play(FadeIn(ax), FadeIn(sums), ShowCreation(zc_line), FadeIn(zc_lab), run_time=0.8)
    s.play(ShowCreation(curve), run_time=1.0)
    s.play(FadeIn(pts), FadeIn(labs), run_time=0.6)
    s.say("$z=0.1,\\ 0.15,\\ 0.2$", "：", "$\\xi^2=%s$" % ",\\ ".join(f4(v) for v in nb["xi2"]), wait=2.0)
    s.say("$z$", " 越接近 ", "$z_c$", "，", "$\\xi$", " 越大——走向臨界點", wait=1.8)
    s.say("書上寫 ", "$(z/z_0)^N$", "，所以說 ", "$z_c=z_0/\\mu$", wait=1.8)
    s.say("用每一步的逸度，就是 ", "$z_c=1/\\mu$", wait=1.6)
    s.clear_stage()


def part_coef(s):
    nb = numbers()
    asym = mixed("係數的漸近（第 2 章補充 3）：", "$[z^N](z_c-z)^{-s}\\approx\\frac{N^{s-1}}{\\Gamma(s)}z_c^{-N-s}$",
                 font_size=28)
    rows = [[zh("精確除以漸近", font_size=22, color=YELLOW)] + ["N=%d" % N for N in WC.COEF_N],
            [r"s=\gamma=%s" % ftex(WC.S_GAMMA)] + [f4(v) for v in nb["coef"][0]],
            [r"s=\gamma+2\nu=%s" % ftex(WC.S_GAMMA2NU)] + [f4(v) for v in nb["coef"][1]]]
    tab = D.tex_table(rows, [3.3, 1.5, 1.5, 1.5], row_h=0.75, font_size=28)
    for i in (1, 2):
        D.table_cell(tab, i, 0, 4).set_color(S_COLORS[i - 1])
    ax = ratio_axes()
    cap = zh("精確除以漸近", font_size=22, color=GREY_A)
    plot = VGroup(cap, ax).arrange(DOWN, buff=0.15)
    box = D.erratum_box("p. 218，(9.38)→(9.39)", r"z_c^{-(N+s)}N^{\delta-1}/\Gamma(s)", r"z_c^{-(N+s)}N^{s-1}/\Gamma(s)",
                        width=8.0, font_size=34)
    place_in(asym, -6.7, 6.7, 2.35, 3.1)
    place_in(tab, -6.7, -0.2, 0.0, 2.1)
    place_in(plot, 0.2, 6.7, -0.55, 2.15)
    place_in(box, -6.7, 6.7, -2.9, -0.85)
    curves = [ratio_curve(ax, sv, c) for sv, c in zip((WC.S_GAMMA, WC.S_GAMMA2NU), S_COLORS)]
    assert_in_frame(asym, tab, plot, box)
    assert_column(asym, tab, box)
    assert_column(asym, plot, box)
    assert_apart(tab, plot)
    assert_table(tab, 4)
    c0 = "、".join(f4(v) for v in nb["coef"][0])
    c1 = "、".join(f4(v) for v in nb["coef"][1])
    s.play(FadeIn(asym), run_time=0.8)
    s.say("係數的漸近（第 2 章補充 3）：冪次 ", "$-s$", " 變成 ", "$N^{s-1}$", wait=2.2)
    s.play(FadeIn(tab[0][:4]), ShowCreation(tab[1]), FadeIn(ax), FadeIn(cap), run_time=0.8)
    s.play(FadeIn(tab[0][4:8]), ShowCreation(curves[0]), run_time=1.0)
    s.say("$s=\\gamma=%s$" % ftex(WC.S_GAMMA), "：精確除以漸近＝", c0, wait=2.2)
    s.play(FadeIn(tab[0][8:]), ShowCreation(curves[1]), run_time=1.0)
    s.say("$s=\\gamma+2\\nu=%s$" % ftex(WC.S_GAMMA2NU), "：", c1, wait=2.2)
    s.say("$N$", " 越大越接近 1：漸近式是對的", wait=1.8)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (9.39) 寫 ", "$N^{\\delta-1}$", "；應該是 ", "$N^{s-1}$", wait=2.0)
    s.say("照書上（", "$\\delta=\\frac{91}{5}$", "），", "$N=10$", " 的比值是 ", "$%s$" % sci(nb["book"]), wait=2.2)
    s.say("完全不對：指數必須是 ", "$s-1$", wait=1.6)
    s.clear_stage()


def part_divide(s):
    nb = numbers()
    zp, zb = nb["zcp"]
    e38 = mixed("(9.38) ", "$\\xi^2\\chi=\\sum_Nc_N\\langle R^2\\rangle_Nz^N=K_1K_2(z_c-z)^{-(\\gamma+2\\nu)}$", font_size=28)
    e39 = mixed("$\\chi=\\sum_Nc_Nz^N=K_2(z_c-z)^{-\\gamma}$", font_size=28)
    e41 = mixed("取 ", "$z^N$", " 的係數相除：",
                "$\\langle R^2\\rangle_N=\\frac{K_1}{z_c^{2\\nu}}\\frac{\\Gamma(\\gamma)}{\\Gamma(\\gamma+2\\nu)}N^{2\\nu}$",
                font_size=28, color=YELLOW)
    eqs = VGroup(e38, e39, e41).arrange(DOWN, buff=0.3)
    rows = [[zh("精確比值除以 (9.41)", font_size=22, color=YELLOW)] + ["N=%d" % N for N in WC.COEF_N],
            [r"\nu=%s" % ftex(NU_2D)] + [f4(v) for v in nb["r2"]]]
    tab = D.tex_table(rows, [3.3, 1.5, 1.5, 1.5], row_h=0.75, font_size=28)
    pw = mixed("$\\mu^{3/2}=%s$" % f4(zp), "；書上 ", "$\\mu^{43/32}=%s$" % f4(zb), font_size=26)
    left = VGroup(tab, pw).arrange(DOWN, buff=0.35)
    box = D.erratum_box("p. 219，(9.41)", r"\frac{K_1}{z_c^{\gamma}}\frac{\Gamma(\gamma)}{\Gamma(\gamma+2\nu)}N^{2\nu}",
                        r"\frac{K_1}{z_c^{2\nu}}\frac{\Gamma(\gamma)}{\Gamma(\gamma+2\nu)}N^{2\nu}", width=6.0, font_size=32)
    place_in(eqs, -6.7, 6.7, 0.35, 3.1)
    place_in(left, -6.7, -0.4, -2.9, 0.1)
    place_in(box, -0.1, 6.7, -2.9, 0.1)
    assert_in_frame(eqs, left, box)
    assert_column(*eqs)
    assert_apart(eqs, left)
    assert_apart(eqs, box)
    assert_apart(left, box)
    assert_apart(tab, pw)
    assert_table(tab, 4)
    r2s = "、".join(f4(v) for v in nb["r2"])
    s.play(FadeIn(e38), run_time=0.8)
    s.say("(9.38)：分子 ", "$\\xi^2\\chi$", " 的奇點是 ", "$(z_c-z)^{-(\\gamma+2\\nu)}$", wait=2.0)
    s.play(FadeIn(e39), run_time=0.6)
    s.say("分母 ", "$\\chi$", " 的奇點是 ", "$(z_c-z)^{-\\gamma}$", wait=1.8)
    s.play(FadeIn(e41), run_time=0.8)
    s.say("兩邊都取 ", "$z^N$", " 的係數，再相除", wait=1.8)
    s.say("$N^{\\gamma+2\\nu-1}$", " 除以 ", "$N^{\\gamma-1}$", "＝", "$N^{2\\nu}$", "：這就是 (9.41)", wait=2.0)
    s.play(FadeIn(tab[0][:4]), ShowCreation(tab[1]), run_time=0.6)
    s.play(FadeIn(tab[0][4:]), run_time=0.6)
    s.say("兩個精確係數之比除以 (9.41)：", r2s, " → 1", wait=2.2)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (9.41) 的分母寫 ", "$z_c^{\\gamma}$", "；應該是 ", "$z_c^{2\\nu}$", wait=2.0)
    s.say("$z_c^{-(N+\\gamma+2\\nu)}$", " 除以 ", "$z_c^{-(N+\\gamma)}$", "$=z_c^{-2\\nu}$", wait=2.0)
    s.play(FadeIn(pw), run_time=0.6)
    s.say("二維：", "$\\mu^{3/2}=%s$" % f4(zp), "；書上 ", "$\\mu^{43/32}=%s$" % f4(zb), wait=2.0)
    s.say("書上說關聯長度 ", "$\\propto K_1/(z_c-z)^{\\nu}$", "；要得到 ", "$K_1K_2$", wait=2.0)
    s.say("應是 ", "$\\xi^2\\approx K_1/(z_c-z)^{2\\nu}$", "（只差寫法，不影響結論）", wait=2.0)
    s.clear_stage()


def part_data(s):
    nb = numbers()
    rows = [["N"] + [str(N) for N in WC.SAW_N],
            [r"2\nu_{\rm eff}"] + [f4(v) for v in nb["nu"]],
            [r"\langle R^2\rangle_N/N^{3/2}"] + [f4(v) for v in nb["amp"]]]
    tab = D.tex_table(rows, [3.0, 1.45, 1.45, 1.45], row_h=0.8, font_size=32)
    arrows = VGroup(tx(r"\to%s" % ftex(2 * NU_2D), 32, D.EXACT_COLOR), tx(r"\to\text{const}", 32, D.EXACT_COLOR))
    for k, a in enumerate(arrows):
        a.next_to(tab[1], RIGHT, buff=0.25)                   # 表格橫線的右端之外
        a.set_y(D.table_cell(tab, k + 1, 3, 4).get_center()[1])
    tabg = VGroup(tab, arrows)
    ax = amp_axes()
    cap = mixed("$\\langle R^2\\rangle_N/N^{3/2}$", "，", "$N=1..14$", font_size=22, color=GREY_A)
    plot = VGroup(cap, ax).arrange(DOWN, buff=0.15)
    note = VGroup(mixed("前因子與 ", "$N$", " 無關", font_size=32, color=YELLOW),
                  mixed("$\\Rightarrow\\langle R^2\\rangle\\propto N^{2\\nu}$", "（de Gennes）", font_size=32, color=YELLOW))
    note.arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    src = zh("第 7 章：二維自迴避漫步的精確列舉", font_size=26, color=GREY_A)
    place_in(src, -6.7, 6.7, 2.5, 3.1)
    place_in(tabg, -6.7, 6.7, 0.4, 2.3)
    place_in(plot, -6.7, -0.4, -2.9, 0.15)
    place_in(note, 0.0, 6.7, -2.4, -0.4)
    dots = amp_dots(ax, nb["amp_all"])
    assert_in_frame(src, tabg, plot, note)
    assert_column(src, tabg, plot)
    assert_apart(tabg, note)
    assert_apart(plot, note)
    assert_apart(tab, arrows)
    assert_table(tab, 4)
    nus = ",\\ ".join(f4(v) for v in nb["nu"])
    amps = ",\\ ".join(f4(v) for v in nb["amp"])
    s.play(FadeIn(src), run_time=0.6)
    s.say("回到第 7 章：二維自迴避漫步的精確列舉", wait=1.8)
    s.play(FadeIn(tab[0][:4]), ShowCreation(tab[1]), FadeIn(tab[0][4:8]), run_time=0.8)
    s.say("局部指數 ", "$2\\nu_{\\rm eff}=%s$" % nus, " → ", "$%s$" % ftex(2 * NU_2D), wait=2.0)
    s.play(FadeIn(arrows[0]), FadeIn(tab[0][8:]), run_time=0.8)
    s.say("$\\langle R^2\\rangle_N/N^{3/2}=%s$" % amps, wait=2.0)
    s.play(FadeIn(plot), run_time=0.6)
    s.play(LaggedStart(*[FadeIn(d) for d in dots], lag_ratio=0.12), run_time=1.5)
    s.play(FadeIn(arrows[1]), run_time=0.5)
    s.say("$N=1..14$", " 一路往下，越來越平：慢慢趨於常數", wait=2.0)
    s.play(FadeIn(note), run_time=0.8)
    s.say("前因子與 ", "$N$", " 無關：", "$\\langle R^2\\rangle\\propto N^{2\\nu}$", "——de Gennes", wait=2.0)
    s.say("磁鐵的 ", "$\\xi$", " 和 ", "$\\chi$", " 只要是冪次，取 ", "$z^N$", " 的係數一除", wait=2.2)
    s.say("就得到漫步的 ", "$\\langle R^2\\rangle\\propto N^{2\\nu}$", wait=2.2)


class W905WalkFromMagnet(ChuScene):
    def construct(self):
        self.show_header("W905　從磁鐵讀出漫步：均方距離正比於 N 的 2ν 次方")
        check_math()
        part_question(self)
        part_gf(self)
        part_coef(self)
        part_divide(self)
        part_data(self)
        self.wait(1.0)
