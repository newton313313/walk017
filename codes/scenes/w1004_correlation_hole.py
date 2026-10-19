"""
W1004 相關洞：頭尾互相排斥（§10.3 後半、§10.4，pp. 243–246、251–253，(10.37)–(10.48)、(10.76)–(10.84)，Fig. 10.2–10.3，Ex. 10.5）

問：擠滿盒子的漫步，頭和尾喜不喜歡靠在一起？
  1. (10.37) G = lim (1/n) ∂²lnZ/∂h(x)∂h(x')；(10.38)–(10.43) 只取結構：平均場 × 自旋波的行列式 × 兩種 λ 積分
     → (h·h'/4){G_0 − (∫f dλ) D}。
  2. r ≫ 0 (10.44)：盒子沒影響（32³ 週期 vs 無限大格子，r = 1，相差 1.2e−11）；r ≪ 0 (10.45)：平均場 × (1 − D)；
     字幕（不畫框）：書上 (10.45) 的指數與 (10.15) 不一致。
  3. (10.36)→(10.46)：鞍點 u|M|² = |r|；64³ 週期格子沿軸的 D、−D(Δ) 的圖（Fig. 10.2）與 Coulomb 形 c0 − c/Δ；
     1 − D：0.8917 → 1.0001；筆誤框 1（(10.46)）；字幕（不畫框）：書上說 D「隨距離增加」。
  4. Ex. 10.5、(10.47)(10.48)：連續三維 D(ρ) 與無限大格子；週期盒子只平移 Ewald 常數 0.0035。
  5. §10.4、Fig. 10.3：sc 普通漫步離起點 R 的步數 vs erfc 形、G(R)·2πR/3 → 1；筆誤框 2（(10.83)）；
     兩條 1/R 曲線並排；收束。
畫面上的數字都由 utils.walk_cond 算出，check_math() 核對。相關洞的 D 在程式裡叫 Dv／hole（D 是 walk_draw 模組）。
模組層級的畫圖函式：hole_axes／hole_curve／hole_dots／coulomb_curve（−D 對 Δ）、cont_axes／cont_curves（連續 D 與 1/(4πρ)）、
visit_axes／visit_curves（步數與 erfc 形）。
"""
import math
import sys
import warnings
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_cond as WK
from utils import walk_draw as D

MF_COLOR = BLUE_B          # 平均場
SW_COLOR = GREEN           # 自旋波
LAM_COLOR = ORANGE         # λ 積分
HOLE_COLOR = TEAL          # 相關洞 −D
FIT_COLOR = ORANGE         # Coulomb 形
DOT_COLOR = YELLOW         # 卡上的點
VIS_COLOR = ORANGE         # erfc 形
BOOK83_R = (8, 16)         # 照書上 (10.83) 的兩個 R


def f4(x):
    return f"{float(x):.4f}"


def sci(x):
    """1.2e-11 → 1.2\\times10^{-11}"""
    m, e = f"{float(x):.1e}".split("e")
    return r"%s\times10^{%d}" % (m, int(e))


def vis_str(x):
    """步數表：≥ 0.001 用四位小數，更小的用 0.000183 的寫法（三位有效數字）。"""
    return f4(x) if x >= 1e-3 else f"{float(x):.6f}"


def erfc_smooth(R, N):
    """(10.83) 更正版的 sc 形 (3/(2πR)) erfc(R√(3/(2(N+1))))，R 可以是小數（WK.erfc_form 會把 R 取整數，只用來畫曲線）。"""
    return 3 / (2 * math.pi * R) * math.erfc(R * math.sqrt(3 / (2 * (N + 1))))


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        axis = [float(v) for v in WK.D_axis()]
        c0, c = WK.coulomb_fit()
        nb = {"case1": WK.case1_check(),
              "argmax": WK.book45_argmax(),
              "axis": axis,
              "hole": [axis[d] for d in WK.HOLE_D],
              "c0": c0, "c": c,
              "book46": float(WK.D_book46()[0]),
              "cont": [WK.D_cont(p) for p in WK.CONT_RHO],
              "inf": [WK.D_inf((p, 0, 0)) for p in WK.CONT_RHO],
              "ewald": WK.SC_EWALD / (4 * math.pi * WK.HOLE_L),
              "visits": [WK.sc_visits(WK.VIS_N, R) for R in WK.VIS_R],
              "erfc": [WK.erfc_form(R, WK.VIS_N) for R in WK.VIS_R],
              "green": [WK.green_sc(R) * 2 * math.pi * R / 3 for R in WK.VIS_R],
              "book83": [WK.erfc_form_book(R, WK.VIS_N) for R in BOOK83_R]}
    return nb


def check_math():
    nb = numbers()
    # (10.44)：32³ 週期 vs 無限大格子（r = 1）
    assert WK.case1_check.__defaults__ == (32, 1.0, 4) and f"{nb['case1']:.1e}" == "1.2e-11"
    assert sci(nb["case1"]) == r"1.2\times10^{-11}"
    # 書上 (10.45) 的極大 |r|/2u = 0.5 vs (10.15) 的 |r|/u = 1（r = −1、u = 1）
    a45, a15 = nb["argmax"]
    assert WK.book45_argmax.__defaults__ == (WK.HOLE_R, 1.0)
    assert f"{a45:.4g}" == "0.5" and f"{a15:.4g}" == "1"
    # 64³ 週期格子、r = −1，沿軸 Δ = 0..32
    assert WK.HOLE_L == 64 and WK.HOLE_R == -1.0 and WK.D_axis.__defaults__ == (64, -1.0)
    assert WK.HOLE_D == (0, 1, 2, 4, 8, 16, 32) and len(nb["axis"]) == 33
    assert [f4(v) for v in nb["hole"]] == ["0.1083", "0.0613", "0.0359", "0.0167", "0.0065", "0.0016", "-0.0001"]
    assert all(a > b for a, b in zip(nb["axis"], nb["axis"][1:]))          # D 嚴格遞減（增加的是 −D）
    assert f4(1 - nb["hole"][0]) == "0.8917" and f4(1 - nb["hole"][-1]) == "1.0001"
    assert f"{100 * nb['hole'][0]:.0f}" == "11"                                # 頭尾重合時少 11%
    assert WK.coulomb_fit.__defaults__ == (64, -1.0, 4, 16)
    assert f4(nb["c0"]) == "0.0035" and f4(nb["c"]) == "0.0801" and f4(1 / (4 * math.pi)) == "0.0796"
    hole_top = max(-v for v in nb["axis"])
    assert -0.12 < -nb["axis"][0] and hole_top < 0.015 and nb["c0"] - nb["c"] / 32 < 0.015   # 都在圖的縱軸範圍裡
    # 筆誤框 1：書上 (10.46) 少了分子的 2|r|（= 2）
    assert f4(nb["book46"]) == "0.0541" and abs(2 * abs(WK.HOLE_R) * nb["book46"] - nb["hole"][0]) < 1e-12
    # Ex. 10.5：連續與無限大格子
    assert WK.CONT_RHO == (1, 2, 4, 8) and WK.D_cont.__defaults__ == (-1.0,) and WK.D_inf.__defaults__ == (-1.0, 1.0)
    assert [f4(v) for v in nb["cont"]] == ["0.0602", "0.0374", "0.0198", "0.0099"]
    assert [f4(v) for v in nb["inf"]] == ["0.0648", "0.0394", "0.0202", "0.0100"]
    assert f"{WK.SC_EWALD}" == "2.837297" and f4(nb["ewald"]) == "0.0035" == f4(nb["c0"])
    # §10.4：sc 格子 N = 100
    assert WK.VIS_N == 100 and WK.VIS_R == (1, 2, 4, 8, 16)
    assert [vis_str(v) for v in nb["visits"]] == ["0.4509", "0.1931", "0.0610", "0.0103", "0.000183"]
    assert [vis_str(v) for v in nb["erfc"]] == ["0.4121", "0.1744", "0.0586", "0.0100", "0.000174"]
    assert [f4(v) for v in nb["green"]] == ["1.0815", "1.0779", "1.0198", "1.0041", "1.0009"]
    assert [f"{v:.1e}" for v in nb["book83"]] == ["4.3e-04", "1.7e-28"]
    assert [sci(v) for v in nb["book83"]] == [r"4.3\times10^{-4}", r"1.7\times10^{-28}"]
    assert vis_str(nb["visits"][WK.VIS_R.index(8)]) == "0.0103"
    assert all(abs(erfc_smooth(R, WK.VIS_N) - v) < 1e-15 for R, v in zip(WK.VIS_R, nb["erfc"]))   # 畫曲線的函式 = 模組
    assert 3 / (2 * math.pi) < 0.5 and erfc_smooth(1.0, WK.VIS_N) < 0.5 and 1 / (4 * math.pi * 0.7) < 0.13


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框 [x0,x1]×[y0,y1]（只縮不放），再置中在方框裡。只用在文字。"""
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


def struct_box(title, body, color):
    """(10.38)–(10.43) 的一塊：標題＋式子，外框。"""
    return framed(VGroup(title, body).arrange(DOWN, buff=0.25), color, buff=0.22)


def dashed(points, color, width=2.5, num_dashes=40):
    c = VMobject().set_points_smoothly(points).set_stroke(color, width)
    return DashedVMobject(c, num_dashes=num_dashes).set_stroke(color, width)


# ---------------------------------------------------------------- 畫圖元件
def hole_axes():
    """−D 對 Δ：Δ ∈ [0, 32]、−D ∈ [−0.12, 0.015]。ManimGL 的橫軸畫在 −D = 0（洞在軸下面），橫軸的刻度數字放在軸上方。"""
    ax = D.plot_axes((0, 34, 8), (-0.12, 0.015, 0.04), width=5.8, height=3.3,
                     y_ticks=(-0.12, -0.08, -0.04, 0), y_tick_labels=["-0.12", "-0.08", "-0.04", "0"], y_label="-D")
    labs = VGroup(*[tx(str(v), 20, GREY_B).next_to(ax.c2p(v, 0), UP, buff=0.12) for v in (8, 16, 24, 32)])
    labs.add(tx(r"\Delta", 24, GREY_B).next_to(ax.x_axis.get_end(), UP, buff=0.12))
    ax.add(labs)
    return ax


def hole_curve(ax, axis_vals):
    """格子上的 −D(Δ)，Δ = 0..32（33 點，折線）。"""
    return VMobject().set_points_as_corners([ax.c2p(d, -v) for d, v in enumerate(axis_vals)]).set_stroke(HOLE_COLOR, 4)


def hole_dots(ax, hole):
    return VGroup(*[Dot(ax.c2p(d, -v), radius=0.06).set_color(DOT_COLOR) for d, v in zip(WK.HOLE_D, hole)])


def coulomb_curve(ax, c0, c):
    """Coulomb 形 c0 − c/Δ，Δ ≥ 2 才畫（60 點，虛線）。"""
    return dashed([ax.c2p(d, c0 - c / d) for d in np.linspace(2, 32, 60)], FIT_COLOR, 3)


def cont_axes():
    """連續 D(ρ)：ρ ∈ [0, 16]、D ∈ [0, 0.12]。"""
    return D.plot_axes((0, 18, 4), (0, 0.13, 0.04), width=5.2, height=2.9, x_ticks=(0, 4, 8, 12, 16),
                       y_ticks=(0.04, 0.08, 0.12), x_label=r"\rho", y_label="D")


def cont_curves(ax, inf_vals):
    """連續 D(ρ)（實線，60 點）、1/(4πρ)（虛線，ρ ≥ 0.7）、無限大格子 ρ = 1, 2, 4, 8（點）。"""
    rs = np.linspace(0, 16, 60)
    rs[0] = 1e-6
    cont = VMobject().set_points_smoothly([ax.c2p(p, WK.D_cont(p)) for p in rs]).set_stroke(HOLE_COLOR, 4)
    coul = dashed([ax.c2p(p, 1 / (4 * math.pi * p)) for p in np.linspace(0.7, 16, 60)], GREY_A)
    dots = VGroup(*[Dot(ax.c2p(p, v), radius=0.06).set_color(DOT_COLOR) for p, v in zip(WK.CONT_RHO, inf_vals)])
    lab = tx(r"\frac{1}{4\pi\rho}", 26, GREY_A).next_to(ax.c2p(3.0, 1 / (12 * math.pi)), UR, buff=0.08)
    return VGroup(cont, coul, dots, lab)


def visit_axes():
    """離起點 R 的步數：R ∈ [0, 16]、∈ [0, 0.5]。"""
    return D.plot_axes((0, 18, 4), (0, 0.52, 0.1), width=5.2, height=2.9, x_ticks=(0, 4, 8, 12, 16),
                       y_ticks=(0.2, 0.4), x_label="R")


def visit_curves(ax, visits):
    """erfc 形（實線，R ∈ [1, 16]，60 點）、3/(2πR)（虛線，N → ∞）、精確 R = 1, 2, 4, 8, 16（點）。"""
    Rs = np.linspace(1, 16, 60)
    erf = VMobject().set_points_smoothly([ax.c2p(R, erfc_smooth(R, WK.VIS_N)) for R in Rs]).set_stroke(VIS_COLOR, 4)
    coul = dashed([ax.c2p(R, 3 / (2 * math.pi * R)) for R in Rs], GREY_A)
    dots = VGroup(*[Dot(ax.c2p(R, v), radius=0.06).set_color(DOT_COLOR) for R, v in zip(WK.VIS_R, visits)])
    lab = tx(r"\frac{3}{2\pi R}", 26, GREY_A).next_to(ax.c2p(3.0, 1 / (2 * math.pi)), UR, buff=0.08)
    return VGroup(erf, coul, dots, lab)


# ---------------------------------------------------------------- 各段
def part_question(s):
    eq = tx(r"G(z;\vec x,\vec x\,')=\lim_{n\to0}\frac{1}{n}\,"
            r"\frac{\partial^2\ln Z}{\partial h(\vec x)\,\partial h(\vec x\,')}", 42)
    top = VGroup(eq, zh("(10.37)", font_size=22, color=GREY_B)).arrange(RIGHT, buff=0.5)
    b1 = struct_box(zh("平均場", font_size=26, color=MF_COLOR),
                    tx(r"f(\lambda)=e^{-N_0(\frac{r}{2}\lambda+\frac{u}{4}\lambda^2)}", 32), MF_COLOR)
    b2 = struct_box(zh("自旋波的行列式", font_size=26, color=SW_COLOR),
                    tx(r"\prod_{\vec q}(\det X_{\vec q})^{-1/2}", 32), SW_COLOR)
    b3 = struct_box(mixed("兩種 ", "$\\lambda$", " 積分", font_size=26, color=LAM_COLOR),
                    tx(r"\int_0^\infty d\lambda\,(\cdots)", 32), LAM_COLOR)
    times = VGroup(tx(r"\times", 40), tx(r"\times", 40))
    row = VGroup(b1, times[0], b2, times[1], b3).arrange(RIGHT, buff=0.3)
    res = tx(r"\frac{\vec h\cdot\vec h\,'}{4}\Big\{G_0-\Big(\int_0^\infty f\,d\lambda\Big)D\Big\}", 44)
    res_box = framed(res, YELLOW, buff=0.25)
    top.move_to(2.45 * UP)
    row.move_to(0.45 * UP)
    res_box.move_to(1.95 * DOWN)
    arrow = Arrow(row.get_bottom(), res_box.get_top(), buff=0.12)
    tag = zh("(10.38)–(10.43)", font_size=22, color=GREY_B).next_to(arrow, RIGHT, buff=0.25)
    assert_in_frame(top, row, res_box, tag)
    assert_column(top, row, arrow, res_box)
    assert_apart(b1, b2)
    assert_apart(b2, b3)
    s.say("問：擠滿盒子的漫步，頭和尾喜不喜歡靠在一起？", wait=2.0)
    s.play(Write(top), run_time=1.2)
    s.say("頭在 ", "$\\vec x$", "、尾在 ", "$\\vec x\\,'$", "：對兩點的磁場各微分一次", wait=2.0)
    s.play(LaggedStart(FadeIn(b1, shift=0.2 * UP), FadeIn(times[0]), FadeIn(b2, shift=0.2 * UP),
                       FadeIn(times[1]), FadeIn(b3, shift=0.2 * UP), lag_ratio=0.3), run_time=1.6)
    s.say("(10.38)–(10.43) 只看結構：三塊相乘", wait=2.0)
    s.play(GrowArrow(arrow), FadeIn(tag), run_time=0.6)
    s.play(Write(res), ShowCreation(res_box[0]), run_time=1.2)
    s.say("$G_0$", " 是沒有洞的部分，再減掉平均場乘上 ", "$D$", wait=2.2)
    s.say("$D$", " 就是相關洞：頭尾靠近時少掉的比例", wait=2.0)
    s.clear_stage()


def part_cases(s):
    nb = numbers()
    left = VGroup(mixed("$r\\gg0$", "：漫步比盒子小", font_size=28, color=MF_COLOR),
                  VGroup(tx(r"G_0=\sum_{\vec q}\frac{\psi_{\vec q}(\vec x)\,\psi_{\vec q}(\vec x\,')}{r+q^2}", 36),
                         zh("(10.44)", font_size=20, color=GREY_B)).arrange(DOWN, buff=0.15),
                  mixed("$32^3$", " 週期格子 vs 無限大格子（", "$r=1$", "）", font_size=24),
                  mixed("相差 ", "$" + sci(nb["case1"]) + "$", "：盒子沒影響", font_size=24)).arrange(DOWN, buff=0.35)
    right = VGroup(mixed("$r\\ll0$", "：漫步擠滿盒子", font_size=28, color=LAM_COLOR),
                   VGroup(tx(r"G\propto\int_0^\infty e^{-N_0(\frac{r}{2}\lambda+\frac{u}{4}\lambda^2)}d\lambda\,(1-D)", 36),
                          zh("(10.45)", font_size=20, color=GREY_B)).arrange(DOWN, buff=0.15),
                   mixed("平均場 ", "$\\times\\,(1-D)$", font_size=24)).arrange(DOWN, buff=0.35)
    pl, pr = framed(left, MF_COLOR, 0.3), framed(right, LAM_COLOR, 0.3)
    pair = VGroup(pl, pr).arrange(RIGHT, buff=0.5, aligned_edge=UP)
    a45, a15 = nb["argmax"]
    book = mixed("書上 (10.45) 印成 ", "$e^{-N_0(r\\lambda+u\\lambda^2)}$", "：極大在 ", "$|r|/2u=%s$" % f"{a45:.4g}",
                 "；(10.15) 的 ", "$f$", " 在 ", "$|r|/u=%s$" % f"{a15:.4g}", font_size=24, color=GREY_A)
    place_in(VGroup(pair, book).arrange(DOWN, buff=0.55), -6.7, 6.7, -2.7, 3.1)
    assert_in_frame(pair, book)
    assert_apart(pl, pr)
    assert_apart(pair, book)
    s.play(FadeIn(pl[0]), FadeIn(left[0]), FadeIn(left[1]), run_time=0.8)
    s.say("$r\\gg0$", "：漫步比盒子小，盒子沒影響", wait=2.0)
    s.play(FadeIn(left[2]), FadeIn(left[3]), run_time=0.6)
    s.say("$32^3$", " 週期格子與無限大格子只差 ", "$" + sci(nb["case1"]) + "$", wait=2.0)
    s.play(FadeIn(pr[0]), FadeIn(right[0]), FadeIn(right[1]), run_time=0.8)
    s.play(FadeIn(right[2]), run_time=0.5)
    s.say("$r\\ll0$", "：漫步擠滿盒子，平均場乘上 ", "$(1-D)$", wait=2.0)
    s.play(FadeIn(book), run_time=0.6)
    s.say("書上 (10.45) 的指數印成 ", "$-N_0(r\\lambda+u\\lambda^2)$", "，和 (10.15) 不一致", wait=2.4)
    s.say("書上的極大在 ", "$|r|/2u=%s$" % f"{a45:.4g}", "，正文是 ", "$|r|/u=%s$" % f"{a15:.4g}", wait=2.4)
    s.clear_stage()


def part_saddle(s):
    nb = numbers()
    lab = mixed("鞍點 ", "$u|\\vec M|^2=|r|$", "（", "$r<0$", "）", font_size=28)
    l1 = tx(r"2\lambda_1=r+q^2+u|\vec M|^2\ \to\ q^2", 34)
    l2 = tx(r"2\lambda_2=r+q^2+3u|\vec M|^2\ \to\ q^2+2|r|", 34)
    l3 = tx(r"\frac{1}{q^2}-\frac{1}{q^2+2|r|}=\frac{2|r|}{q^2(q^2+2|r|)}", 40)
    col = VGroup(lab, l1, l2, l3).arrange(DOWN, buff=0.32, aligned_edge=LEFT)
    eq46 = tx(r"D=\sum_{\vec q\neq0}\frac{2|r|\,\psi_{\vec q}(\vec x)\,\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}", 40)
    box46 = VGroup(framed(eq46, YELLOW, 0.28), zh("(10.36) → (10.46)", font_size=22, color=GREY_B))
    box46[1].next_to(box46[0], UP, buff=0.15)
    VGroup(col, box46).arrange(RIGHT, buff=0.7)
    place_in(VGroup(col, box46), -6.7, 6.7, 0.0, 3.1)
    err = D.erratum_box("p. 245，(10.46)",
                        r"D=\sum_{\vec q}\frac{\psi_{\vec q}(\vec x)\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}",
                        r"D=\sum_{\vec q\neq0}\frac{2|r|\,\psi_{\vec q}(\vec x)\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}",
                        width=11.5)
    err.move_to(1.75 * DOWN)
    assert_in_frame(col, box46, err)
    assert_apart(col, box46)
    assert_apart(VGroup(col, box46), err)
    s.play(FadeIn(lab), run_time=0.6)
    s.play(Write(l1), Write(l2), run_time=1.2)
    s.say("鞍點上：橫向 ", "$2\\lambda_1\\to q^2$", "、縱向 ", "$2\\lambda_2\\to q^2+2|r|$", wait=2.2)
    s.play(Write(l3), run_time=1.0)
    s.say("兩者相減，", "$\\vec q=0$", " 不算：這就是 (10.46) 的 ", "$D$", wait=2.0)
    s.play(FadeIn(box46), run_time=0.8)
    s.wait(0.6)
    s.play(FadeIn(err, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (10.46) 漏了分子的 ", "$2|r|$", "，也沒去掉 ", "$\\vec q=0$", wait=2.2)
    s.say("照書上 ", "$D(0)=%s$" % f4(nb["book46"]), "，只有正確值的一半", wait=2.0)
    s.clear_stage()


def part_hole(s):
    nb = numbers()
    hole = nb["hole"]
    lab = mixed("$64^3$", " 週期格子、", "$r=-1$", "，沿 ", "$x$", " 軸", font_size=24, color=GREY_A)
    rows = [[r"\Delta"] + [str(d) for d in WK.HOLE_D], ["D"] + [f4(v) for v in hole]]
    tab = D.tex_table(rows, [0.9] + [1.15] * 7, row_h=0.5, font_size=26)
    tab.shift(np.array([-tab.get_width() / 2, 2.62, 0]))
    lab.next_to(tab, UP, buff=0.12)
    ax = hole_axes()
    ax.move_to(np.array([-3.35, -0.75, 0]))
    curve = hole_curve(ax, nb["axis"])
    dots = hole_dots(ax, hole)
    fit = coulomb_curve(ax, nb["c0"], nb["c"])
    t1 = mixed("Coulomb 形：", "$-D\\approx c_0-\\frac{c}{\\Delta}$", "（", "$\\Delta\\ge2$", "）", font_size=24)
    t2 = mixed("$c_0=%s$" % f4(nb["c0"]), "、", "$c=%s\\approx\\frac{1}{4\\pi}=%s$" % (f4(nb["c"]), f4(1 / (4 * math.pi))),
               font_size=24)
    t3 = mixed("$1-D$", "：", "$%s$" % f4(1 - hole[0]), "（頭尾重合）", "$\\to%s$" % f4(1 - hole[-1]), font_size=24)
    t4 = zh("頭尾靠近時條數少 %s%%：相關洞" % f"{100 * hole[0]:.0f}", font_size=26, color=HOLE_COLOR)
    t1[0].set_color(FIT_COLOR)
    txt = VGroup(t1, t2, t3, t4).arrange(DOWN, buff=0.38, aligned_edge=LEFT)
    place_in(txt, 0.5, 6.7, -2.4, 0.9)
    txt.align_to(np.array([0.6, 0, 0]), LEFT)
    assert_in_frame(lab, tab, ax, txt)
    assert_table(tab, 8)
    assert_column(lab, tab, ax)
    assert_apart(ax, txt)
    assert_column(t1, t2, t3, t4)
    s.play(FadeIn(lab), FadeIn(tab), run_time=0.9)
    s.say("$64^3$", " 週期格子上直接求和：", "$D$", " 沿軸一路遞減", wait=2.2)
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(ShowCreation(curve), FadeIn(dots), run_time=1.4)
    s.say("畫出 ", "$-D(\\Delta)$", "：Fig. 10.2 的樣子，頭尾附近一個洞", wait=2.2)
    s.play(ShowCreation(fit), FadeIn(t1), FadeIn(t2), run_time=1.2)
    s.say("遠處是 Coulomb 形 ", "$c_0-c/\\Delta$", "，", "$c\\approx1/4\\pi$", wait=2.2)
    s.play(FadeIn(t3), run_time=0.6)
    s.say("$1-D$", "：頭尾重合時 %s，遠處 %s" % (f4(1 - hole[0]), f4(1 - hole[-1])), wait=2.0)
    s.play(FadeIn(t4), run_time=0.6)
    s.say("頭尾靠近時條數少 %s%%：這就是相關洞" % f"{100 * hole[0]:.0f}", wait=2.0)
    s.say("書上說 ", "$D$", "「隨距離增加」——其實 ", "$D$", " 遞減", wait=2.2)
    s.say("增加的是 ", "$-D$", "：Fig. 10.2 畫的就是 ", "$-D$", wait=2.0)
    s.clear_stage()


def part_cont(s):
    nb = numbers()
    head = mixed("連續三維（Ex. 10.5、(10.47)(10.48)）", font_size=26, color=GREY_A)
    eq = tx(r"D(\rho)=\frac{1}{4\pi\rho}\Big(1-e^{-\sqrt{2|r|}\,\rho}\Big)\ \to\ \frac{1}{4\pi\rho}", 44)
    note = mixed("（遠處與 ", "$r$", " 無關）", font_size=24, color=GREY_A)
    rows = [[r"\rho"] + [str(p) for p in WK.CONT_RHO],
            [zh("連續", font_size=24, color=HOLE_COLOR)] + [f4(v) for v in nb["cont"]],
            [zh("無限大格子", font_size=24, color=DOT_COLOR)] + [f4(v) for v in nb["inf"]]]
    tab = D.tex_table(rows, [2.4] + [1.35] * 4, row_h=0.55, font_size=28)
    l1 = mixed("遠處就是 Coulomb：", "$D\\propto|\\vec x-\\vec x\\,'|^{2-d}$", "（", "$d=3$", "：", "$1/\\rho$", "）",
               font_size=26)
    l2 = mixed("週期盒子只平移一個常數 ", "$\\frac{%s}{4\\pi L}=%s$" % (f"{WK.SC_EWALD}", f4(nb["ewald"])),
               "（", "$L=%d$" % WK.HOLE_L, "，Ewald 常數）", font_size=26)
    col = VGroup(head, VGroup(eq, note).arrange(RIGHT, buff=0.35), tab, l1, l2).arrange(DOWN, buff=0.38)
    place_in(col, -6.7, 6.7, -2.8, 3.1)
    assert_in_frame(col)
    assert_table(tab, 5)
    assert_column(*col)
    s.play(FadeIn(head), Write(eq), run_time=1.2)
    s.say("Ex. 10.5：連續空間裡的 ", "$D(\\rho)$", " 有閉式", wait=2.0)
    s.play(FadeIn(note), run_time=0.5)
    s.say("遠處 ", "$D\\to1/4\\pi\\rho$", "，和 ", "$r$", " 無關", wait=2.0)
    s.play(FadeIn(tab), run_time=0.9)
    s.play(FadeIn(l1), run_time=0.6)
    s.say("無限大格子也一樣：", "$D\\propto|\\vec x-\\vec x\\,'|^{2-d}$", wait=2.2)
    s.play(FadeIn(l2), run_time=0.6)
    s.say("週期盒子只讓 ", "$D$", " 整體平移 %s：就是上面的 " % f4(nb["ewald"]), "$c_0$", wait=2.4)
    s.clear_stage()


def part_visits(s):
    nb = numbers()
    head = mixed("§10.4、Fig. 10.3：普通漫步「離起點 ", "$R$", " 的步數」（sc 格子、", "$N=%d$" % WK.VIS_N, "）",
                 font_size=26, color=GREY_A)
    eq = tx(r"\sum_{n\le N}P_n(R)\ \approx\ \frac{3}{2\pi R}\,\mathrm{erfc}\Big(R\sqrt{\frac{3}{2(N+1)}}\Big)", 40)
    rows = [["R"] + [str(R) for R in WK.VIS_R],
            [zh("精確", font_size=24, color=DOT_COLOR)] + [vis_str(v) for v in nb["visits"]],
            [zh("erfc 形", font_size=24, color=VIS_COLOR)] + [vis_str(v) for v in nb["erfc"]],
            [tx(r"N\to\infty:\ G(R)\times\frac{2\pi R}{3}", 24, GREY_A)] + [f4(v) for v in nb["green"]]]
    tab = D.tex_table(rows, [3.3] + [1.45] * 5, row_h=0.56, font_size=28)
    err = D.erratum_box("p. 253，(10.83)",
                        r"\frac{1}{4\pi z_cR}\,\mathrm{erfc}\Big(\frac{z_cR^2}{2(N+1)}\Big)",
                        r"\frac{1}{4\pi z_cR}\,\mathrm{erfc}\Big(\sqrt{\frac{z_cR^2}{4(N+1)}}\Big)", width=9.5)
    col = VGroup(head, eq, tab, err).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -2.85, 3.1)
    assert_in_frame(col)
    assert_table(tab, 6)
    assert_column(*col)
    s.play(FadeIn(head), run_time=0.6)
    s.say("普通漫步走 ", "$N$", " 步：停在離起點 ", "$R$", " 處的步數", wait=2.0)
    s.play(Write(eq), run_time=1.2)
    s.play(FadeIn(tab), run_time=0.9)
    s.say("sc 格子、", "$N=%d$" % WK.VIS_N, "：精確值與 erfc 形", wait=2.0)
    s.say("$N\\to\\infty$", "：", "$G(R)\\times2\\pi R/3\\to1$", "，也就是 ", "$\\propto1/R$", wait=2.4)
    s.play(FadeIn(err, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (10.83) 的 erfc 引數少了根號、", "$2$", " 應為 ", "$4$", wait=2.2)
    i8 = WK.VIS_R.index(8)
    s.say("照書上 ", "$R=8$", " 得 ", "$%s$" % sci(nb["book83"][0]), "（精確 %s）" % vis_str(nb["visits"][i8]), wait=2.2)
    s.say("$R=16$", " 得 ", "$%s$" % sci(nb["book83"][1]), "：完全不對", wait=2.0)
    s.clear_stage()


def part_compare(s):
    nb = numbers()
    axl, axr = cont_axes(), visit_axes()
    VGroup(axl, axr).arrange(RIGHT, buff=1.4).move_to(0.35 * DOWN)
    cl, cr = cont_curves(axl, nb["inf"]), visit_curves(axr, nb["visits"])
    tl = mixed("相關洞 ", "$D(\\rho)$", "（自迴避）", font_size=26, color=HOLE_COLOR).next_to(axl, UP, buff=0.3)
    tr = mixed("離起點 ", "$R$", " 的步數（普通漫步）", font_size=26, color=VIS_COLOR).next_to(axr, UP, buff=0.3)
    both = mixed("兩條都是 ", "$\\propto1/R$", "：同一個 Coulomb 形", font_size=28, color=YELLOW)
    both.next_to(VGroup(axl, axr), DOWN, buff=0.35)
    assert_in_frame(axl, axr, tl, tr, both, cl, cr)
    assert_apart(axl, axr)
    assert_apart(tl, tr)
    assert_column(tl, axl, both)
    s.play(ShowCreation(axl), ShowCreation(axr), FadeIn(tl), FadeIn(tr), run_time=1.0)
    s.play(ShowCreation(cl[0]), ShowCreation(cr[0]), run_time=1.4)
    s.play(FadeIn(cl[1]), FadeIn(cr[1]), FadeIn(cl[2]), FadeIn(cr[2]), FadeIn(cl[3]), FadeIn(cr[3]), run_time=0.9)
    s.say("左：相關洞 ", "$D$", "；右：普通漫步離起點 ", "$R$", " 的步數", wait=2.2)
    s.play(FadeIn(both, shift=0.2 * UP), run_time=0.6)
    s.say("兩條都是 ", "$\\propto1/R$", "：同一個 Coulomb 形", wait=2.0)
    s.say("步數的密度 ", "$\\propto1/R$", "，自迴避的懲罰 ", "$-u/R$", wait=2.2)
    s.say("就讓頭尾互相排斥——一個 Coulomb 形的洞", wait=2.4)


class W1004CorrelationHole(ChuScene):
    def construct(self):
        self.show_header("W1004　相關洞：頭尾互相排斥")
        check_math()
        part_question(self)
        part_cases(self)
        part_saddle(self)
        part_hole(self)
        part_cont(self)
        part_visits(self)
        part_compare(self)
        self.wait(1.0)
