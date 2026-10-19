"""
W1005 Redner–Reynolds：漫步的密度是序參量（§10.3.1，pp. 246–251，(10.49)–(10.75)，Table 10.1）

問：把一條漫步的步數當成粒子數——它的「密度」能不能像磁化一樣當序參量？
  1. (10.49)–(10.52)：Γ = Σ Γ_N e^{−rN}（巨正則，r 是化學勢）、Ω = (1/V) ln Γ、⟨N⟩ = −∂lnΓ/∂r、ρ = ⟨N⟩/V = −∂Ω/∂r。
  2. 環面精確列舉 L = 3, 4, 5：ρ(z) 三條曲線、z_c 虛線；表 z = 0.2, z_c, 0.6, 1.0, ∞；(10.59) Var N/V 的峰（z*、峰高）。
  3. 平均場 (10.53)–(10.58)：ρ_MF（V = 10⁶、u = 1）、斜率 1 ⇒ β' = 1；χ' ⇒ γ' = α' = 0；ν' = 1/2；ρ_MF(r) 的圖；
     筆誤框 1（(10.53)(10.55) 的負號）。
  4. 標度律 (10.65)–(10.75)、Table 10.1（平均場、二維、三維）；驗算 α'+2β'+γ' = 2、2 − α' = dν'、d_c = 4；
     筆誤框 2（(10.75) 的 d_c）。
  5. 第 10 章回顧（五列表）；收束。
畫面上的數字都由 utils.walk_cond（z_c 用 utils.walk_crit.ZC_2D）算出，check_math() 核對。
模組層級的畫圖函式：rho_axes／rho_curves／zc_marker／torus_legend（ρ 對 z）、mf_axes／mf_curve（ρ_MF 對 r）。
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
from utils import walk_crit as WC
from utils import walk_draw as D

TORUS_COLORS = {3: TEAL, 4: YELLOW, 5: ORANGE}
MF_COLOR = BLUE_B
OK_COLOR = "#7CE38B"
CASES = ("mf", "d2", "d3")
MF_CURVE_R = (-0.6, 0.6, 60)         # ρ_MF(r) 的圖：r ∈ [−0.6, 0.6]、60 點（不含 r = 0）
BIG_Z = 1e4                          # z → ∞ 的數值核對


def f4(x):
    return f"{float(x):.4f}"


def sci(x):
    """1.0e-05 → 1.0\\times10^{-5}"""
    m, e = f"{float(x):.1e}".split("e")
    return r"%s\times10^{%d}" % (m, int(e))


def ftex(q, inline=False):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    if inline:
        return "%s%d/%d" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


def rr_str(v, inline=False):
    """Table 10.1 的數：分數照分數寫，三維（引用 ν = 0.588）寫三位小數。"""
    return ftex(v, inline) if isinstance(v, Fraction) else f"{v:.3f}"


def tf(tex):
    """表格裡的分數用 \\tfrac（上下兩列都是分數時才不會碰在一起）。"""
    return tex.replace(r"\frac", r"\tfrac")


def fit_width(m, w):
    """只縮不放（文字）。"""
    if m.get_width() > w:
        m.set_width(w)
    return m


def small_str(x):
    """ρ_MF、χ'：≥ 0.001 用四位小數，更小的用科學記號 1.0e-05。"""
    return f4(x) if abs(x) >= 1e-3 else f"{x:.1e}"


def small_tex(x):
    return f4(x) if abs(x) >= 1e-3 else sci(x)


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rr = {k: WK.rr_exponents(*WK.RR_CASES[k]) for k in CASES}
        nb = {"curve": {L: WK.torus_curve(L) for L in WK.TORI},
              "rho": {L: [float(WK.torus_rho(L, z)) for z in WK.RHO_Z] for L in WK.TORI},
              "full": {L: Fraction(len(WK.torus_counts(L)) - 1, L * L) for L in WK.TORI},
              "rho_big": {L: float(WK.torus_rho(L, BIG_Z)) for L in WK.TORI},
              "peak": {L: WK.torus_peak(L) for L in WK.TORI},
              "mf_rho": [WK.mf_rho(r) for r in WK.MF_RHO_R],
              "mf_chi": [WK.mf_chi(r) for r in WK.MF_RHO_R],
              "mf_curve": [(r, WK.mf_rho(r)) for r in np.linspace(*MF_CURVE_R)],
              "slope": math.log(WK.mf_rho(-0.1) / WK.mf_rho(-0.01)) / math.log(10),
              "book55": WK.book55_rho(-0.5),
              "rr": rr,
              "laws": {k: WK.rr_laws(rr[k], WK.RR_CASES[k][1]) for k in CASES},
              "dc": WK.dc_mean_field(),
              "c4": WK.torus_counts(4),
              "first4": WK.torus_first_diff(4),
              "hole0": float(WK.D_axis()[0]),
              "modes": WK.mode_table()}
    return nb


def check_math():
    nb = numbers()
    # §2 環面：L = 3, 4, 5；z = 0.2, z_c, 0.6, 1.0
    assert WK.TORI == (3, 4, 5) and WK.RHO_Z == (0.2, WC.ZC_2D, 0.6, 1.0) and f4(WC.ZC_2D) == "0.3791"
    assert WK.torus_curve.__defaults__ == (60, 1.2)
    for L in WK.TORI:
        zs, rs = nb["curve"][L]
        assert len(zs) == 60 and zs[0] == 0 and abs(zs[-1] - 1.2) < 1e-12 and 0 <= min(rs) and max(rs) < 1
    assert [f4(v) for v in nb["rho"][3]] == ["0.1391", "0.3374", "0.5284", "0.6823"]
    assert [f4(v) for v in nb["rho"][4]] == ["0.0900", "0.3056", "0.5585", "0.7338"]
    assert [f4(v) for v in nb["rho"][5]] == ["0.0603", "0.2769", "0.5764", "0.7565"]
    assert [str(nb["full"][L]) for L in WK.TORI] == ["8/9", "15/16", "24/25"]
    assert all(abs(nb["rho_big"][L] - nb["full"][L]) < 1e-3 for L in WK.TORI)        # z → ∞：(V − 1)/V
    assert [f4(nb["peak"][L][0]) for L in WK.TORI] == ["0.4443", "0.4380", "0.4272"]
    assert [f4(nb["peak"][L][1]) for L in WK.TORI] == ["0.4289", "0.5846", "0.7242"]
    zst = [nb["peak"][L][0] for L in WK.TORI]
    hts = [nb["peak"][L][1] for L in WK.TORI]
    assert zst[0] > zst[1] > zst[2] > WC.ZC_2D and hts[0] < hts[1] < hts[2]       # 峰往 z_c 移、越來越高
    # §3 平均場：V = 10⁶、u = 1
    assert WK.MF_RHO_N0 == 10 ** 6 and WK.MF_U == 1.0 and WK.mf_rho.__defaults__ == (1.0, 10 ** 6)
    assert WK.mf_chi.__defaults__ == (1.0, 10 ** 6) and WK.MF_RHO_R == (-0.5, -0.1, -0.01, 0.1, 0.5)
    assert [small_str(v) for v in nb["mf_rho"]] == ["0.2500", "0.0500", "0.0050", "1.0e-05", "2.0e-06"]
    assert [small_str(v) for v in nb["mf_chi"]] == ["0.5000", "0.5000", "0.5000", "1.0e-04", "4.0e-06"]
    assert [small_tex(v) for v in nb["mf_rho"][3:]] == [r"1.0\times10^{-5}", r"2.0\times10^{-6}"]
    assert f4(nb["slope"]) == "1.0000" and f4(nb["book55"]) == "-0.2500" and nb["book55"] < 0
    assert all(0 <= p < 0.33 for _, p in nb["mf_curve"]) and 0.0 not in [r for r, _ in nb["mf_curve"]]
    # §4 Table 10.1：(ν, d) = (1/2, 4)、(3/4, 2)、(0.588, 3)
    assert WK.RR_CASES == {"mf": (Fraction(1, 2), 4), "d2": (Fraction(3, 4), 2), "d3": (0.588, 3)}
    tab = {k: [rr_str(nb["rr"][k][e]) for e in ("beta", "gamma", "alpha", "nu")] for k in CASES}
    assert tab == {"mf": ["1", "0", "0", r"\frac{1}{2}"],
                   "d2": [r"\frac{1}{2}", r"\frac{1}{2}", r"\frac{1}{2}", r"\frac{3}{4}"],
                   "d3": ["0.764", "0.236", "0.236", "0.588"]}, tab
    for k in CASES:
        s4, two_a, dnu, dcnu = nb["laws"][k]
        assert abs(s4 - 2) < 1e-12 and abs(two_a - dnu) < 1e-12                          # 兩條標度律都成立
    assert [rr_str(nb["laws"][k][1]) for k in CASES] == ["2", r"\frac{3}{2}", "1.764"]
    assert [rr_str(nb["laws"][k][1], True) for k in CASES] == ["2", "3/2", "1.764"]
    assert nb["laws"]["d2"][3] == 3 and nb["laws"]["d2"][1] == Fraction(3, 2)               # 書上 d_c ν' = 3 ≠ 3/2
    assert nb["dc"] == 4 and nb["rr"]["mf"]["alpha"] == 0 and nb["rr"]["mf"]["nu"] == Fraction(1, 2)
    # §5 回顧
    assert nb["c4"][4] == 96 and WK.SAW_INF[4] == 100 and nb["c4"][15] == 5712 and nb["first4"] == (4, 4)
    assert f4(1 - nb["hole0"]) == "0.8917" and f"{100 * nb['hole0']:.0f}" == "11"
    assert [(n, t, l) for n, t, l, _, _ in nb["modes"]] == [(3, 2, 1), (2, 1, 1), (1, 0, 1), (0, -1, 1)]
    assert nb["rr"]["mf"]["beta"] == 1 and nb["rr"]["mf"]["gamma"] == 0


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框 [x0,x1]×[y0,y1]（只縮不放），再置中在方框裡。只用在文字與表格。"""
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


def framed(mob, color, buff=0.25):
    return VGroup(SurroundingRectangle(mob, buff=buff).set_stroke(color, 2), mob)


# ---------------------------------------------------------------- 畫圖元件
def rho_axes():
    """ρ 對 z：z ∈ [0, 1.2]、ρ ∈ [0, 1]。"""
    return D.plot_axes((0, 1.3, 0.2), (0, 1.05, 0.25), width=5.4, height=3.7, x_ticks=(0, 0.4, 0.8, 1.2),
                       y_ticks=(0.5, 1), y_tick_labels=["0.5", "1"], x_label="z", y_label=r"\rho")


def rho_curves(ax, curves):
    """WK.torus_curve(L) 的 60 點，一條一個 VMobject。"""
    return VGroup(*[VMobject().set_points_smoothly([ax.c2p(z, r) for z, r in zip(*curves[L])])
                    .set_stroke(TORUS_COLORS[L], 4) for L in WK.TORI])


def zc_marker(ax):
    line = DashedLine(ax.c2p(WC.ZC_2D, 0), ax.c2p(WC.ZC_2D, 1.0), dash_length=0.08).set_stroke(GREY_A, 2)
    return VGroup(line, tx("z_c", 26, GREY_A).next_to(line, UP, buff=0.08))


def torus_legend():
    rows = VGroup(*[VGroup(Line(ORIGIN, 0.4 * RIGHT).set_stroke(TORUS_COLORS[L], 5),
                           tx(r"%d\times%d" % (L, L), 24)).arrange(RIGHT, buff=0.12) for L in WK.TORI])
    return rows.arrange(DOWN, buff=0.12, aligned_edge=LEFT)


def mf_axes():
    """ρ_MF 對 r：r ∈ [−0.6, 0.6]、ρ ∈ [0, 0.32]；縱軸在 r = 0（相變點），縱軸的刻度數字放在軸右邊。"""
    ax = D.plot_axes((-0.6, 0.65, 0.2), (0, 0.33, 0.1), width=4.8, height=2.3, x_ticks=(-0.4, 0.4),
                     x_label="r", y_label=r"\rho_{MF}")
    ax.add(VGroup(*[tx(f"{v:.1f}", 20, GREY_B).next_to(ax.c2p(0, v), RIGHT, buff=0.12) for v in (0.1, 0.2, 0.3)]))
    return ax


def mf_curve(ax, pts):
    """ρ_MF(r)，60 點折線（r = 0 有轉折）。"""
    return VMobject().set_points_as_corners([ax.c2p(r, p) for r, p in pts]).set_stroke(MF_COLOR, 4)


# ---------------------------------------------------------------- 各段
def part_question(s):
    eqs = VGroup(tx(r"\Gamma=\sum_N\Gamma_N\,e^{-rN}", 40),
                 tx(r"\Omega=\frac{1}{V}\ln\Gamma", 40),
                 tx(r"\langle N\rangle=-\frac{\partial}{\partial r}\ln\Gamma", 40),
                 tx(r"\rho=\frac{\langle N\rangle}{V}=-\frac{\partial\Omega}{\partial r}", 40)).arrange(DOWN, buff=0.4,
                                                                                                  aligned_edge=LEFT)
    tag = zh("(10.49)–(10.52)", font_size=22, color=GREY_B)
    left = VGroup(tag, eqs).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    dic = VGroup(zh("漫步 vs 氣體、磁鐵", font_size=28, color=YELLOW),
                 mixed("步數 ", "$N$", "$\\leftrightarrow$", "粒子數", font_size=28),
                 mixed("$r$", "$\\leftrightarrow$", "化學勢（", "$e^{-r}=z/z_c$", "）", font_size=28),
                 mixed("密度 ", "$\\rho$", "$\\leftrightarrow$", "磁化？序參量？", font_size=28)).arrange(
        DOWN, buff=0.4, aligned_edge=LEFT)
    box = framed(dic, YELLOW, 0.35)
    row = VGroup(left, box).arrange(RIGHT, buff=1.0)
    place_in(row, -6.7, 6.7, -2.7, 3.1)
    assert_in_frame(row)
    assert_apart(left, box)
    assert_column(tag, *eqs)
    s.say("問：把一條漫步的步數當成粒子數——", wait=1.6)
    s.say("它的「密度」能不能像磁化一樣當序參量？", wait=2.0)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.0)
    s.say("巨正則：", "$\\Gamma=\\sum_N\\Gamma_Ne^{-rN}$", "，", "$r$", " 是化學勢", wait=2.0)
    s.say("$\\Gamma_N$", "：走 ", "$N$", " 步、終點任意的漫步條數", wait=2.0)
    s.play(FadeIn(box[0]), FadeIn(dic[0]), FadeIn(dic[1]), FadeIn(dic[2]), run_time=0.8)
    s.say("這章 ", "$e^{-r}=z/z_c$", "（第 9 章是 ", "$e^{-rN/2}$", "）", wait=2.0)
    s.play(Write(eqs[1]), Write(eqs[2]), Write(eqs[3]), run_time=1.6)
    s.play(FadeIn(dic[3]), run_time=0.5)
    s.say("密度 ", "$\\rho=-\\partial\\Omega/\\partial r$", "：每個格點平均分到幾步", wait=2.2)
    s.clear_stage()


def part_torus(s):
    nb = numbers()
    ax = rho_axes()
    ax.move_to(np.array([-3.5, 0.05, 0]))
    curves = rho_curves(ax, nb["curve"])
    zc = zc_marker(ax)
    leg = torus_legend().next_to(ax.c2p(1.25, 0.04), UL, buff=0)
    zlab = ["0.2", "z_c", "0.6", "1.0"]
    rows = [["z"] + [r"%d\times%d" % (L, L) for L in WK.TORI]]
    rows += [[zlab[k]] + [f4(nb["rho"][L][k]) for L in WK.TORI] for k in range(len(WK.RHO_Z))]
    rows += [[r"\infty"] + [str(nb["full"][L]) for L in WK.TORI]]
    tabA = D.tex_table(rows, [1.0, 1.3, 1.3, 1.3], row_h=0.47, font_size=26)
    titleB = mixed("(10.59) 漲落 ", "$\\mathrm{Var}\\,N/V$", " 的峰", font_size=24, color=GREY_A)
    rowsB = [[zh("漲落峰", font_size=22, color=YELLOW)] + [r"%d\times%d" % (L, L) for L in WK.TORI],
             ["z^*"] + [f4(nb["peak"][L][0]) for L in WK.TORI],
             [zh("峰高", font_size=22)] + [f4(nb["peak"][L][1]) for L in WK.TORI]]
    tabB = D.tex_table(rowsB, [1.0, 1.3, 1.3, 1.3], row_h=0.47, font_size=26)
    right = VGroup(tabA, titleB, tabB).arrange(DOWN, buff=0.3)
    place_in(right, 0.9, 6.7, -2.6, 3.1)
    assert_in_frame(ax, zc, leg, right, curves)
    assert_apart(ax, right)
    assert_table(tabA, 4)
    assert_table(tabB, 4)
    assert_column(tabA, titleB, tabB)
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(*[ShowCreation(c) for c in curves], FadeIn(leg), run_time=2.0)
    s.say("環面上精確列舉：", "$3\\times3$", "、", "$4\\times4$", "、", "$5\\times5$", " 三種盒子", wait=2.0)
    s.say("密度 ", "$\\rho=\\langle N\\rangle/V$", " 隨逸度 ", "$z$", " 一路增加", wait=2.0)
    s.play(ShowCreation(zc[0]), FadeIn(zc[1]), FadeIn(tabA), run_time=1.0)
    s.say("在 ", "$z_c=%s$" % f4(WC.ZC_2D), " 附近，", "$\\rho$", " 從很小變成有限值", wait=2.2)
    s.say("$z\\to\\infty$", "：走滿全部格點，", "$\\rho\\to(V-1)/V$", wait=2.0)
    s.say("$z$", " 小：漫步很短；", "$z$", " 大：漫步擠滿盒子", wait=2.0)
    s.play(FadeIn(titleB), FadeIn(tabB), run_time=0.9)
    s.say("(10.59) 漲落 ", "$\\mathrm{Var}\\,N/V$", " 在 ", "$z^*$", " 有一個峰", wait=2.0)
    s.say("格子越大，峰越高越陡、", "$z^*$", " 往 ", "$z_c$", " 移", wait=2.2)
    s.say("這就是有限大小的相變", wait=1.6)
    s.clear_stage()


def part_mf(s):
    nb = numbers()
    lab = VGroup(mixed("平均場（", "$V=10^6$", "、", "$u=1$", "）", font_size=26),
                 zh("(10.53)–(10.58)", font_size=22, color=GREY_B)).arrange(DOWN, buff=0.15)
    eq = tx(r"\rho_{MF}=\begin{cases}\dfrac{1}{V}\,\dfrac{1}{r}\to0, & r>0\\[2mm] \dfrac{|r|}{2u}, & r<0\end{cases}", 34)
    top = VGroup(lab, eq).arrange(RIGHT, buff=0.8)
    rows = [["r"] + [f"{r:g}" for r in WK.MF_RHO_R],
            [r"\rho_{MF}"] + [small_tex(v) for v in nb["mf_rho"]],
            [r"\chi'=\partial^2\Omega/\partial r^2"] + [small_tex(v) for v in nb["mf_chi"]]]
    tab = D.tex_table(rows, [2.9] + [1.6] * 5, row_h=0.5, font_size=26)
    upper = VGroup(top, tab).arrange(DOWN, buff=0.3)
    place_in(upper, -6.7, 6.7, 0.15, 3.1)
    ax = mf_axes()
    ax.move_to(np.array([-3.85, -1.5, 0]))
    curve = mf_curve(ax, nb["mf_curve"])
    t1 = mixed("斜率 ", "$\\frac{d\\ln\\rho}{d\\ln|r|}=%s$" % f4(nb["slope"]), "$\\Rightarrow\\beta'=1$", font_size=24)
    t2 = mixed("（磁鐵 ", "$\\beta=\\frac{1}{2}$", "：", "$\\rho\\propto|\\vec M|^2$", "）", font_size=24, color=GREY_A)
    t3 = mixed("$\\chi'$", " 不發散 ", "$\\Rightarrow\\gamma'=\\alpha'=0$", font_size=24)
    t4 = mixed("(10.62)–(10.64)：", "$\\nu'=\\frac{1}{2}$", font_size=24)
    txt = VGroup(t1, t2, t3, t4).arrange(DOWN, buff=0.28, aligned_edge=LEFT)
    place_in(txt, -0.6, 6.7, -2.75, -0.15)
    err = D.erratum_box("pp. 247–248，(10.53)(10.55)", r"e^{-\frac{V}{a^d}\frac{r^2}{4u}}",
                        r"e^{+\frac{V}{a^d}\frac{r^2}{4u}}", width=7.0, font_size=40)
    err.move_to(txt)
    assert_in_frame(upper, ax, curve, txt, err)
    assert_table(tab, 6)
    assert_apart(lab, eq)
    assert_column(upper, ax)
    assert_column(upper, txt)
    assert_apart(ax, txt)
    assert_apart(ax, err)
    assert_column(t1, t2, t3, t4)
    s.play(FadeIn(lab), Write(eq), run_time=1.2)
    s.say("平均場：", "$r>0$", " 時 ", "$\\rho\\to0$", "；", "$r<0$", " 時 ", "$\\rho=|r|/2u$", wait=2.2)
    s.play(FadeIn(tab), run_time=0.9)
    s.say("$V=10^6$", "：", "$r=-0.5,-0.1,-0.01$", " 的密度正比於 ", "$|r|$", wait=2.2)
    s.play(ShowCreation(ax), run_time=0.7)
    s.play(ShowCreation(curve), run_time=1.2)
    s.say("$\\rho$", " 就像磁化：在 ", "$r=0$", " 兩邊從 0 變成有限", wait=2.0)
    s.play(FadeIn(t1), run_time=0.5)
    s.say("斜率 %s：" % f4(nb["slope"]), "$\\rho\\propto|r|^1$", "，所以 ", "$\\beta'=1$", wait=2.0)
    s.play(FadeIn(t2), run_time=0.5)
    s.say("磁鐵的 ", "$\\beta=1/2$", "：因為 ", "$\\rho\\propto|\\vec M|^2$", wait=2.0)
    s.play(FadeIn(t3), FadeIn(t4), run_time=0.6)
    s.say("$\\chi'$", " 在 ", "$r=0$", " 不發散：", "$\\gamma'=\\alpha'=0$", "；", "$\\nu'=1/2$", wait=2.4)
    s.play(FadeOut(txt), FadeIn(err, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (10.53)(10.55) 的指數多了一個負號", wait=1.8)
    s.say("照書上的負號：", "$\\rho=%s<0$" % f4(nb["book55"]), "（", "$r=-0.5$", "）", wait=2.2)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    rr, laws = nb["rr"], nb["laws"]
    head = [zh("指數", font_size=24, color=YELLOW),
            mixed("平均場 ", "$d=%d$" % WK.RR_CASES["mf"][1], font_size=22, color=YELLOW),
            mixed("二維 ", "$\\nu=%s$" % tf(rr_str(WK.RR_CASES["d2"][0])), font_size=22, color=YELLOW),
            mixed("三維 ", "$\\nu\\approx%s$" % rr_str(WK.RR_CASES["d3"][0]), font_size=22, color=YELLOW)]
    rows = [head,
            [r"\beta'=d\nu-1=1-\alpha"] + [tf(rr_str(rr[k]["beta"])) for k in CASES],
            [r"\gamma'=\alpha'=2-d\nu=\alpha"] + [tf(rr_str(rr[k]["gamma"])) for k in CASES],
            [r"\nu'=\nu"] + [tf(rr_str(rr[k]["nu"])) for k in CASES]]
    widths = [4.4, 2.5, 2.5, 2.5]
    tab = D.tex_table(rows, widths, row_h=0.64, font_size=28)
    chk = [[tx(r"\alpha'+2\beta'+\gamma'", 28)] + [tx(r"2\ \checkmark", 28, OK_COLOR) for _ in CASES],
           [r"2-\alpha'"] + [tf(rr_str(laws[k][1])) for k in CASES],
           [r"d\,\nu'"] + [tf(rr_str(laws[k][2])) for k in CASES]]
    tab2 = D.tex_table(chk, widths, row_h=0.58, font_size=28)
    dc = mixed("(10.66) 平均場：", "$2-0=d_c\\cdot\\frac{1}{2}$", "$\\Rightarrow d_c=%d$" % nb["dc"], font_size=26)
    err = D.erratum_box("p. 251，(10.75)", r"2-\alpha'=d_c\,\nu'", r"2-\alpha'=d\,\nu'", width=5.6)
    VGroup(dc, err).arrange(RIGHT, buff=0.6)
    col = VGroup(tab, tab2, VGroup(dc, err)).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -2.85, 3.1)
    assert_in_frame(col)
    assert_table(tab, 4)
    assert_table(tab2, 4)
    assert_column(tab, tab2, err)
    assert_column(tab2, dc)
    assert_apart(dc, err)
    s.play(FadeIn(tab[1]), FadeIn(tab[0][:4]), run_time=0.8)
    s.say("標度律 (10.65)–(10.75)：漫步的指數都由 ", "$\\nu$", " 決定", wait=2.0)
    s.play(LaggedStart(*[FadeIn(tab[0][4 * i:4 * i + 4]) for i in range(1, 4)], lag_ratio=0.5), run_time=1.5)
    s.say("Table 10.1：", "$\\beta'=d\\nu-1$", "、", "$\\gamma'=\\alpha'=2-d\\nu$", "、", "$\\nu'=\\nu$", wait=2.4)
    s.say("三維的 ", "$\\nu\\approx0.588$", " 是引用的數值", wait=1.8)
    s.play(FadeIn(tab2), run_time=0.9)
    s.say("驗算：", "$\\alpha'+2\\beta'+\\gamma'=2$", " 三欄都成立", wait=2.0)
    s.say("$2-\\alpha'=d\\,\\nu'$", "：", "、".join(rr_str(laws[k][1], True) for k in CASES),
          " 也都相等", wait=2.0)
    s.play(FadeIn(dc), run_time=0.6)
    s.say("平均場：", "$2-0=d_c\\cdot1/2$", "，所以 ", "$d_c=%d$" % nb["dc"], wait=2.0)
    s.play(FadeIn(err, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (10.75) 寫成 ", "$d_c\\,\\nu'$", "，還說「一般都成立」", wait=2.2)
    s.say("但 ", "$d=2$", " 時 ", "$d_c\\nu'=%s\\neq%s$" % (rr_str(laws["d2"][3]), rr_str(laws["d2"][1], True)), wait=2.2)
    s.clear_stage()


def part_recap(s):
    nb = numbers()
    c4, modes = nb["c4"], nb["modes"]
    n0 = modes[-1]
    title = zh("第 10 章回顧", font_size=34, color=YELLOW)
    fs = 24
    rows = [[zh("幕", font_size=24, color=YELLOW), zh("重點", font_size=24, color=YELLOW)],
            [zh("W1001 關在盒子裡", font_size=fs),
             mixed("$4\\times4$", " 環面 ", "$c_N$", " 從 ", "$N=%d$" % nb["first4"][0],
                   " 起變少（%d vs %d）、走滿 %d 條" % (c4[4], WK.SAW_INF[4], c4[15]), font_size=fs)],
            [zh("W1002 平均場", font_size=fs),
             mixed("$\\Gamma_{MF}=\\frac{N_0}{2}\\int f$", "（書上 ", "$\\frac{N_0}{4}$", "）、先取 ", "$n\\to0$",
                   font_size=fs)],
            [zh("W1003 自旋波", font_size=fs),
             mixed("$\\lambda_1$", "（", "$n-1$", " 重）、", "$\\lambda_2$", "（%d 重）、" % n0[2], "$n\\to0$", " 時 ",
                   "$%d$" % n0[1], " 個模態、", "$\\sqrt{\\lambda_1/\\lambda_2}$", font_size=fs)],
            [zh("W1004 相關洞", font_size=fs),
             mixed("洞 ", "$1-D$", "：", "$%s\\to1$" % f4(1 - nb["hole0"]), "、", "$D\\approx\\frac{1}{4\\pi\\Delta}$",
                   "、erfc 形", font_size=fs)],
            [zh("W1005 密度", font_size=fs),
             mixed("$\\beta'=%s$" % rr_str(nb["rr"]["mf"]["beta"]), "、",
                   "$\\gamma'=\\alpha'=%s$" % rr_str(nb["rr"]["mf"]["gamma"]), "、", "$d_c=%d$" % nb["dc"], font_size=fs)]]
    for row in rows[1:]:
        fit_width(row[1], 9.6)
    tab = D.tex_table(rows, [2.9, 9.9], row_h=0.8, font_size=fs)
    col = VGroup(title, tab).arrange(DOWN, buff=0.35)
    place_in(col, -6.7, 6.7, -2.8, 3.1)
    assert_in_frame(col)
    assert_table(tab, 2)
    s.play(FadeIn(title), FadeIn(tab[1]), FadeIn(tab[0][:2]), run_time=0.8)
    s.say("第 10 章回顧：五幕的重點", wait=1.4)
    notes = [("盒子小、漫步長：從 ", "$N=%d$" % nb["first4"][0], " 起比無限平面少"),
             ("平均場變成一維積分：", "$N_0/2$", "，不是書上的 ", "$N_0/4$"),
             ("自旋波：", "$n\\to0$", " 時橫向模態變成 ", "$%d$" % n0[1], " 個"),
             ("相關洞：頭尾靠近時少 %s%%，遠處是 Coulomb 形" % f"{100 * nb['hole0']:.0f}",),
             ("密度是序參量：", "$\\beta'=%s$" % rr_str(nb["rr"]["mf"]["beta"]), "、",
              "$\\gamma'=\\alpha'=%s$" % rr_str(nb["rr"]["mf"]["gamma"]), "、", "$d_c=%d$" % nb["dc"])]
    for i, note in enumerate(notes, start=1):
        s.play(FadeIn(tab[0][2 * i]), FadeIn(tab[0][2 * i + 1]), run_time=0.6)
        s.say(*note, wait=1.6)
    s.say("漫步的密度在 ", "$z_c$", " 兩邊從 0 變成有限值", wait=2.2)
    s.say("——擠滿盒子的漫步就是 ", "$\\mathrm{O}(n)$", " 磁鐵的有序相", wait=2.6)


class W1005RednerReynolds(ChuScene):
    def construct(self):
        self.show_header("W1005　Redner–Reynolds：漫步的密度是序參量")
        check_math()
        part_question(self)
        part_torus(self)
        part_mf(self)
        part_table(self)
        part_recap(self)
        self.wait(1.0)
