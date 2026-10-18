"""
W904 Ginzburg-Landau-Wilson：讓自旋可以伸縮（§9.3，pp. 215–217，(9.28)–(9.35)，Ex. 9.3–9.4）

問：如果自旋不是固定長度的箭頭，而是可以伸縮的「場」，磁鐵的自由能長什麼樣子？
  1. (9.28) GLW 的有效 Hamiltonian（r、梯度、u、外場四項）、(9.29) Z（沒有 k_BT）；固定長度的箭頭 vs 可伸縮的場。
  2. (9.30)(9.31) Fourier；Ex. 9.3 的離散版（8 點環、種子 904）三個等式誤差 < 1e-12；(9.32)(9.33) 只留 q = 0 → 平均場。
  3. 雙井 r = 1、0、−1（u = 1）：m* = ±1、F_min = −0.25；小球滾到一邊；h = 0.1：m* = 1.0467；β = 1/2。
  4. Ex. 9.4（u = 1/3）：m* = √3 = 1.7321、F_min = −0.75；塌縮 ℳ(1) = 2.1038；數值斜率 β、Δ、α、1/δ、γ。
  5. (9.34) 有限體積的平均場積分：⟨m⟩、⟨|m|⟩、−ln Z/V（V = 10、100、1000）；收束。
畫面上的數字都由 utils.walk_crit 算出，check_math() 核對。
模組層級的畫圖函式：spin_row、field_ring（自旋箭頭）、well_curve（雙井）、weight_curve（(9.34) 的被積函數）。
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

WARM = D.HEAD_COLOR                   # 朝上（正）
COLD = D.TAIL_COLOR                   # 朝下（負）
R_COLORS = (BLUE_B, GREY_B, YELLOW)   # r = 1、0、−1
TILT_COLOR = ORANGE
V_COLORS = (TEAL, ORANGE)             # V = 10、100
FIELD_SEED = 904                      # 與 WC.fourier_check() 的預設種子相同
SLOPE_TH = (Fraction(1, 2), Fraction(3, 2), Fraction(0), Fraction(1, 3), Fraction(1))   # β、Δ、α、1/δ、γ


def f2(x):
    return f"{float(x):.2f}"


def f4(x):
    return f"{float(x):.4f}"


def fg(x):
    return f"{float(x):g}"


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
        nb = {"fourier": WC.fourier_check(),
              "lmin": [WC.landau_min(r) for r in WC.LANDAU_R],
              "lc": WC.landau_closed(-1.0),
              "mh": WC.landau_min(-1.0, 1.0, WC.LANDAU_H),
              "ex94": WC.ex94_closed(-1.0),
              "ex94m": WC.ex94_m(-1.0, 0.0),
              "M1": [WC.ex94_scaling(1.0, r) for r in (-0.01, -0.0001)],
              "slopes": WC.ex94_slopes(),
              "zm_h": [WC.zmf_avg(V, h=WC.ZMF_H) for V in WC.ZMF_V],
              "zm_0": [WC.zmf_avg(V) for V in WC.ZMF_V],
              "zabs": [WC.zmf_avg(V, what="absm") for V in WC.ZMF_V],
              "zf": [WC.zmf_avg(V, what="f") for V in WC.ZMF_V],
              "mstar_h": WC.landau_min(-1.0, 1.0, WC.ZMF_H)}
    return nb


def check_math():
    nb = numbers()
    assert WC.LANDAU_R == (1.0, 0.0, -1.0) and WC.LANDAU_H == 0.1 and WC.ZMF_V == (10, 100, 1000) and WC.ZMF_H == 0.01
    assert len(field_values()) == 8 and all(e < 1e-12 for e in nb["fourier"])          # Ex. 9.3：8 點環
    assert nb["lmin"] == [0.0, 0.0, 1.0]                                               # r = 1、0：m* = 0；r = −1：1
    m1, fmin = nb["lc"]
    assert fg(m1) == "1" and f2(fmin) == "-0.25" and abs(WC.landau_F(m1, -1.0) - fmin) < 1e-15
    assert all(abs(WC.landau_closed(r)[0] - math.sqrt(-r)) < 1e-15 for r in (-0.5, -2.0))  # m* = √(−r/u)
    assert f4(nb["mh"]) == "1.0467" and fg(WC.LANDAU_H) == "0.1"
    m3, f3 = nb["ex94"]
    assert f4(m3) == "1.7321" and abs(m3 - math.sqrt(3)) < 1e-12 and f2(f3) == "-0.75"
    assert abs(nb["ex94m"] - m3) < 1e-9 and abs(WC.EX94_U - 1 / 3) < 1e-15
    assert [f4(x) for x in nb["M1"]] == ["2.1038", "2.1038"] and abs(nb["M1"][0] - nb["M1"][1]) < 1e-9
    assert [f4(x) for x in nb["slopes"]] == ["0.5000", "1.5000", "0.0000", "0.3333", "1.0000"]
    assert all(abs(a - float(b)) < 1e-6 for a, b in zip(nb["slopes"], SLOPE_TH))
    b, dl, a, di, g = SLOPE_TH
    assert dl / b == 3 and 1 / di == 3 and b * (dl / b - 1) == g and a == 0              # δ = Δ/β = 3、β(δ−1) = γ
    assert [f4(x) for x in nb["zm_h"]] == ["0.0869", "0.7576", "1.0042"] and f4(nb["mstar_h"]) == "1.0050"
    assert all(abs(x) < 1e-12 for x in nb["zm_0"])
    assert [f4(x) for x in nb["zabs"]] == ["0.8888", "0.9922", "0.9992"] and nb["zabs"][-1] < m1
    assert [f4(x) for x in nb["zf"]] == ["-0.2719", "-0.2397", "-0.2478"]


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


# ---------------------------------------------------------------- 畫圖元件
def field_values():
    """畫圖用的場：與 WC.fourier_check() 同一組 8 個常態亂數（種子 904）。"""
    return np.random.default_rng(FIELD_SEED).normal(size=8)


def spin_arrow(p, vec, color):
    """從 p 出發的短箭頭；太短（< 0.12）就只畫一小段線。"""
    if np.linalg.norm(vec) < 0.12:
        return Line(p, p + vec).set_stroke(color, 3)
    return Arrow(p, p + vec, buff=0, thickness=2.5).set_color(color)


def spin_row(vals, dx=0.6, unit=0.6, fixed=False, x0=0.0):
    """一排格點（灰點，y = 0）＋箭頭：fixed=True 長度固定 unit（只看正負），否則長度 |S_i|·unit；朝上暖色、朝下冷色。"""
    dots, arrows = VGroup(), VGroup()
    for i, v in enumerate(vals):
        p = (x0 + i * dx) * RIGHT
        dots.add(Dot(p, radius=0.05).set_color(GREY_B))
        ln = unit * (math.copysign(1.0, v) if fixed else v)
        arrows.add(spin_arrow(p, ln * UP, WARM if v > 0 else COLD))
    return VGroup(dots, arrows)


def spin_pair(vals, gap=1.8):
    """兩排並列、格點在同一條水平線上：左＝固定長度（O(n) 自旋）、右＝可伸縮的場；標籤對齊在兩排最低點之下。
    回傳 VGroup(左排, 右排, 左標籤, 右標籤)。"""
    dx = 0.6
    fixed = spin_row(vals, dx=dx, fixed=True)
    field = spin_row(vals, dx=dx, x0=(len(vals) - 1) * dx + gap)
    ybot = min(fixed.get_bottom()[1], field.get_bottom()[1]) - 0.3
    lab1 = mixed("$O(n)$", " 自旋：長度固定", font_size=28, color=GREY_A)
    lab2 = mixed("GLW：可以伸縮的場 ", "$\\vec S_i$", font_size=28, color=GREY_A)
    lab1.move_to(np.array([fixed[0].get_center()[0], ybot, 0]), aligned_edge=UP)
    lab2.move_to(np.array([field[0].get_center()[0], ybot, 0]), aligned_edge=UP)
    return VGroup(fixed, field, lab1, lab2)


def field_ring(vals, radius=1.1, unit=0.42):
    """Ex. 9.3 的 8 點環：格點排成一圈，場畫成徑向箭頭（正＝向外、暖色；負＝向內、冷色），長度 |S_i|·unit。"""
    circ = Circle(radius=radius).set_stroke(GREY_D, 1.5)
    dots, arrows = VGroup(), VGroup()
    for i, v in enumerate(vals):
        a = PI / 2 - TAU * i / len(vals)
        d = np.array([math.cos(a), math.sin(a), 0.0])
        dots.add(Dot(radius * d, radius=0.05).set_color(GREY_A))
        arrows.add(spin_arrow(radius * d, v * unit * d, WARM if v > 0 else COLD))
    return VGroup(circ, dots, arrows)


def well_axes():
    """雙井的座標：m ∈ [−2, 2]、F ∈ [−0.5, 1.5]；ManimGL 的軸畫在 0，所以平移成 (m + 2, F + 0.5)。"""
    return D.plot_axes((0, 4.4, 1), (0, 2.2, 0.5), width=6.4, height=4.9, x_ticks=(0, 1, 2, 3, 4),
                       y_ticks=(0, 0.5, 1, 1.5, 2), x_tick_labels=["-2", "-1", "0", "1", "2"],
                       y_tick_labels=["-0.5", "0", "0.5", "1", "1.5"], x_label="m", y_label="F")


def wp(ax, m, fv):
    return ax.c2p(m + 2, fv + 0.5)


def well_curve(ax, r, h=0.0, color=WHITE, fmax=1.5, n=81):
    """F(m) = r m²/2 + m⁴/4 − h m（u = 1）在畫框內（F ≤ fmax）的那一段。"""
    ms = np.linspace(-2, 2, 801)
    ok = ms[WC.landau_F(ms, r, 1.0, h) <= fmax]
    pts = [wp(ax, m, WC.landau_F(m, r, 1.0, h)) for m in np.linspace(ok.min(), ok.max(), n)]
    return VMobject().set_points_smoothly(pts).set_stroke(color, 3)


def weight_curve(ax, V, h, color, n=401):
    """(9.34) 的被積函數 e^{−V(F(m) − F_min)}（r = −1、u = 1；峰高 1）；橫軸同樣平移 m + 2。"""
    fmin = WC.landau_F(WC.landau_min(-1.0, 1.0, h), -1.0, 1.0, h)
    ms = np.linspace(-2, 2, n)
    w = np.exp(-V * (WC.landau_F(ms, -1.0, 1.0, h) - fmin))
    return VMobject().set_points_as_corners([ax.c2p(m + 2, y) for m, y in zip(ms, w)]).set_stroke(color, 2.5)


# ---------------------------------------------------------------- 各段
def part_glw(s):
    spins = spin_pair(field_values())
    fixed, field, lab1, lab2 = spins
    e28 = tx(r"\mathcal H_{GL}=\sum_i\Big\{\frac r2|\vec S_i|^2+\frac12|\nabla\vec S_i|^2"
             r"+\frac u4|\vec S_i|^4-\vec h\cdot\vec S_i\Big\}", 42)
    t1 = mixed("$\\frac r2|\\vec S_i|^2$", "：", "$r$", " 是約化溫度", font_size=26)
    t2 = mixed("$\\frac12|\\nabla\\vec S_i|^2$", "：相鄰的場要對齊", font_size=26)
    t3 = mixed("$\\frac u4|\\vec S_i|^4$", "：不讓場變得太長", font_size=26)
    t4 = mixed("$-\\vec h\\cdot\\vec S_i$", "：外場", font_size=26)
    terms = VGroup(VGroup(t1, t3).arrange(DOWN, aligned_edge=LEFT, buff=0.3),
                   VGroup(t2, t4).arrange(DOWN, aligned_edge=LEFT, buff=0.3)).arrange(RIGHT, buff=1.0)
    e29 = mixed("(9.29) ", "$Z=\\prod_i\\int d^nS_i\\,e^{-\\mathcal H_{GL}}$", font_size=30)
    place_in(spins, -6.7, 6.7, 0.95, 3.1)
    place_in(e28, -6.7, 6.7, 0.0, 0.85)
    place_in(terms, -6.7, 6.7, -1.55, -0.15)
    place_in(e29, -6.7, 6.7, -2.9, -1.75)
    assert_in_frame(spins, e28, terms, e29)
    assert_column(spins, e28, terms, e29)
    assert_apart(fixed, field)
    assert_apart(lab1, lab2)
    assert_apart(VGroup(fixed, field), VGroup(lab1, lab2))
    assert_apart(terms[0], terms[1])
    s.say("問：如果自旋不是固定長度的箭頭", wait=0.2)
    s.play(FadeIn(fixed[0]), LaggedStart(*[GrowArrow(a) if isinstance(a, Arrow) else ShowCreation(a) for a in fixed[1]],
                                         lag_ratio=0.1), FadeIn(lab1), run_time=1.5)
    s.wait(0.8)
    s.say("而是可以伸縮的「場」——磁鐵的自由能長什麼樣子？", wait=0.2)
    s.play(FadeIn(field[0]), LaggedStart(*[GrowArrow(a) if isinstance(a, Arrow) else ShowCreation(a) for a in field[1]],
                                         lag_ratio=0.1), FadeIn(lab2), run_time=1.5)
    s.wait(1.5)
    s.play(Write(e28), run_time=2.0)
    s.say("(9.28)：Ginzburg-Landau-Wilson 的有效 Hamiltonian", wait=1.6)
    s.play(FadeIn(t1), run_time=0.6)
    s.say("$r$", " 是約化溫度：", "$r>0$", " 高溫、", "$r<0$", " 低溫", wait=1.6)
    s.play(FadeIn(t2), FadeIn(t3), run_time=0.6)
    s.say("梯度項讓相鄰的場對齊；四次項不讓場變得太長", wait=2.0)
    s.play(FadeIn(t4), run_time=0.6)
    s.say("$-\\vec h\\cdot\\vec S_i$", "：外場把場拉向 ", "$\\vec h$", " 的方向", wait=1.5)
    s.play(FadeIn(e29), run_time=0.8)
    s.say("(9.29)：每個格點對 ", "$n$", " 個分量積分", wait=1.5)
    s.say("沒有 ", "$k_BT$", "：溫度藏在 ", "$r$", " 裡——所以叫「有效」", wait=2.0)
    s.clear_stage()


def part_fourier(s):
    ring = field_ring(field_values())
    ring_g = VGroup(ring, zh("Ex. 9.3：8 點的環", font_size=24, color=GREY_A)).arrange(DOWN, buff=0.3)
    e30 = mixed("(9.30)(9.31) ", "$\\vec S(\\vec q)=\\frac1N\\sum_i\\vec S_i\\,e^{i\\vec q\\cdot\\vec r_i}$", font_size=30)
    ids = VGroup(tx(r"\sum_iS_i^2=N\sum_q|S(q)|^2", 32),
                 tx(r"\sum_i(S_{i+1}-S_i)^2=N\sum_q4\sin^2\frac q2\,|S(q)|^2", 32),
                 tx(r"\sum_iS_i^4=N\sum_{q_1+q_2+q_3+q_4=0}S(q_1)S(q_2)S(q_3)S(q_4)", 32)).arrange(DOWN, aligned_edge=LEFT,
                                                                                             buff=0.32)
    errs = VGroup(*[mixed("誤差 ", "$<10^{-12}$", font_size=22, color=D.EXACT_COLOR) for _ in ids])
    for e, row in zip(errs, ids):
        e.next_to(ids, RIGHT, buff=0.4)
        e.set_y(row.get_center()[1])
    block = VGroup(ids, errs)
    e33 = mixed("(9.32)(9.33) 只留 ", "$\\vec q=0$", "：", "$\\mathcal H_{MF}=V\\Big[\\frac r2m^2+\\frac u4m^4-hm\\Big]$",
                font_size=30)
    place_in(ring_g, -6.7, -3.7, -0.4, 3.1)
    place_in(e30, -3.4, 6.7, 2.35, 3.1)
    place_in(block, -3.4, 6.7, -0.4, 2.1)
    place_in(e33, -6.7, 6.7, -2.8, -0.9)
    assert_in_frame(ring_g, e30, block, e33)
    assert_column(e30, block, e33)
    assert_apart(ring_g, e30)
    assert_apart(ring_g, block)
    assert_apart(ring_g, e33)
    assert_column(*ids)
    assert_apart(ids, errs)
    s.play(FadeIn(ring_g), run_time=1.0)
    s.play(FadeIn(e30), run_time=1.0)
    s.say("(9.30)(9.31)：把場寫成 Fourier 分量 ", "$\\vec S(\\vec q)$", wait=1.6)
    s.say("Ex. 9.3 的離散版：8 點的環、隨機的場（種子 904）", wait=1.6)
    s.play(FadeIn(ids[0]), run_time=0.6)
    s.say("平方項：", "$\\sum_iS_i^2=N\\sum_q|S(q)|^2$", wait=1.5)
    s.play(FadeIn(ids[1]), run_time=0.6)
    s.say("梯度項多了 ", "$4\\sin^2\\frac q2$", "——小 ", "$q$", " 時就是 ", "$q^2$", wait=1.6)
    s.play(FadeIn(ids[2]), run_time=0.6)
    s.say("四次項：四個 ", "$q$", " 加起來是 0 的和", wait=1.5)
    s.play(LaggedStart(*[FadeIn(e) for e in errs], lag_ratio=0.3), run_time=1.0)
    s.say("三個等式的相對誤差都小於 ", "$10^{-12}$", wait=1.5)
    s.play(FadeIn(e33), run_time=0.8)
    s.say("只留 ", "$\\vec q=0$", "：每個 ", "$\\vec S_i$", " 都換成同一個 ", "$m$", wait=1.6)
    s.say("得到平均場：", "$V$", " 乘上 ", "$m$", " 的多項式", wait=1.5)
    s.say("(9.32) 的 ", "$V$", "、", "$(2\\pi)^d$", " 因子書上沒對好，不影響 (9.33)", wait=2.2)
    s.clear_stage()


def part_well(s):
    nb = numbers()
    m1, fmin = nb["lc"]
    ax = well_axes()
    place_in(ax, -6.7, -0.6, -2.9, 3.0)
    curves = [well_curve(ax, r, color=c) for r, c in zip(WC.LANDAU_R, R_COLORS)]
    zero = DashedLine(wp(ax, -2, 0), wp(ax, 2, 0), dash_length=0.08).set_stroke(GREY_D, 1.5)
    dots = VGroup(*[Dot(wp(ax, sg * m1, fmin), radius=0.07).set_color(YELLOW) for sg in (-1, 1)])
    ball = Dot(wp(ax, 0, 0) + 0.12 * UP, radius=0.12).set_color(WHITE)
    path = VMobject().set_points_smoothly([wp(ax, m, WC.landau_F(m, -1.0)) + 0.12 * UP for m in np.linspace(0, m1, 25)])
    tilt = well_curve(ax, -1.0, h=WC.LANDAU_H, color=TILT_COLOR)
    mh = nb["mh"]
    dot_h = Dot(wp(ax, mh, WC.landau_F(mh, -1.0, 1.0, WC.LANDAU_H)), radius=0.07).set_color(TILT_COLOR)
    lab_f = mixed("(9.33) ", "$F(m)=\\frac r2m^2+\\frac u4m^4-hm$", font_size=26)
    legend = VGroup(*[mixed("$r=%s$" % fg(r), font_size=26, color=c) for r, c in zip(WC.LANDAU_R, R_COLORS)])
    legend.arrange(RIGHT, buff=0.5)
    lab_min = mixed("$r=-1$", "：", "$m^*=\\pm%s$" % fg(m1), "、", "$F_{\\min}=%s$" % f2(fmin), font_size=26, color=YELLOW)
    lab_h = mixed("$h=%s$" % fg(WC.LANDAU_H), "：", "$m^*=%s$" % f4(mh), font_size=26, color=TILT_COLOR)
    lab_b = mixed("$m^*=\\sqrt{-r/u}\\propto(-r)^{1/2}$", " → ", "$\\beta=\\frac12$", font_size=26)
    col = VGroup(lab_f, legend, lab_min, lab_h, lab_b).arrange(DOWN, aligned_edge=LEFT, buff=0.42)
    place_in(col, -0.3, 6.7, -2.6, 3.0)
    assert_in_frame(ax, col)
    assert_apart(ax, col)
    assert_column(*col)
    s.play(FadeIn(ax), FadeIn(lab_f), run_time=0.8)
    s.say("(9.33) 的自由能，先取 ", "$h=0$", "、", "$u=1$", wait=1.4)
    s.play(ShowCreation(curves[0]), FadeIn(legend[0]), run_time=1.0)
    s.say("$r=1$", "：只有一個谷，在 ", "$m=0$", wait=1.4)
    s.play(ShowCreation(curves[1]), FadeIn(legend[1]), run_time=1.0)
    s.say("$r=0$", "：谷底變平——臨界點", wait=1.4)
    s.play(ShowCreation(curves[2]), FadeIn(legend[2]), ShowCreation(zero), run_time=1.0)
    s.say("$r=-1$", "：中間鼓起，兩個谷 ", "$m^*=\\pm1$", wait=1.4)
    s.play(FadeIn(dots), FadeIn(lab_min), run_time=0.8)
    s.say("谷底 ", "$F_{\\min}=-0.25$", "（", "$=-r^2/4u$", "）", wait=1.4)
    s.play(FadeIn(ball), run_time=0.5)
    s.say("小球放在中間的山頂上", wait=0.2)
    s.play(MoveAlongPath(ball, path), run_time=1.5)
    s.say("它只能滾到其中一邊：對稱破缺", wait=1.4)
    s.play(curves[0].animate.set_stroke(opacity=0.3), curves[1].animate.set_stroke(opacity=0.3),
           ShowCreation(tilt), run_time=1.0)
    s.play(FadeIn(dot_h), FadeIn(lab_h), run_time=0.6)
    s.say("外場 ", "$h=0.1$", " 把雙井歪向右：", "$m^*=1.0467$", wait=1.6)
    s.play(FadeIn(lab_b), run_time=0.6)
    s.say("谷底 ", "$m^*=\\sqrt{-r/u}\\propto(-r)^{1/2}$", " → ", "$\\beta=\\frac12$", wait=1.8)
    s.clear_stage()


def part_ex94(s):
    nb = numbers()
    m3, f3 = nb["ex94"]
    e35 = mixed("(9.35) ", "$F=\\frac r2m^2+\\frac{m^4}{12}-hm$", "（", "$u=\\frac13$", "）", font_size=28)
    dF = tx(r"\frac{\partial F}{\partial m}=rm+\frac{m^3}{3}-h=0", 36)
    sol = mixed("$r=-1$", "：", "$m^*=\\sqrt3=%s$" % f4(m3), "、", "$F_{\\min}=%s$" % f2(f3), font_size=28)
    coll = tx(r"m=|r|^{1/2}\,\mathcal M\big(h/|r|^{3/2}\big)", 36)
    left = VGroup(e35, dF, sol, coll).arrange(DOWN, aligned_edge=LEFT, buff=0.32)
    ax = D.plot_axes((0, 3.3, 1), (0, 3.3, 1), width=4.6, height=2.7, x_ticks=(0, 1, 2, 3), y_ticks=(1, 2, 3),
                     x_label="x", y_label=r"\mathcal M")
    note = mixed("$x=h/|r|^{3/2}$", "，", "$\\mathcal M=m/|r|^{1/2}$", font_size=20, color=GREY_B)
    leg = mixed("線 ", "$r=-0.01$", "　點 ", "$r=-0.0001$", font_size=20, color=GREY_A)
    right = VGroup(ax, note, leg).arrange(DOWN, buff=0.18)
    place_in(left, -6.7, -0.5, -0.35, 3.1)
    place_in(right, -0.2, 6.7, -0.45, 3.1)
    c1 = VMobject().set_points_smoothly([ax.c2p(x, WC.ex94_scaling(x, -0.01)) for x in np.linspace(0, 3, 31)])
    c1.set_stroke(TEAL, 3)
    dts = VGroup(*[Dot(ax.c2p(x, WC.ex94_scaling(x, -0.0001)), radius=0.05).set_color(YELLOW)
                   for x in np.linspace(0, 3, 13)])
    M1 = nb["M1"][0]
    guide = VGroup(DashedLine(ax.c2p(1, 0), ax.c2p(1, M1), dash_length=0.06),
                   DashedLine(ax.c2p(0, M1), ax.c2p(1, M1), dash_length=0.06)).set_stroke(GREY_B, 1.5)
    mlab = tx(r"\mathcal M(1)=%s" % f4(M1), 24, YELLOW).move_to(ax.c2p(1.12, 1.45), aligned_edge=LEFT)
    th = [ftex(q) for q in SLOPE_TH]
    rows = [[zh("指數", font_size=24, color=YELLOW), r"\beta", r"\Delta", r"\alpha", r"1/\delta", r"\gamma"],
            [zh("數值斜率", font_size=24)] + [f4(x) for x in nb["slopes"]],
            [zh("理論", font_size=24)] + th]
    tab = D.tex_table(rows, [2.2] + [1.5] * 5, row_h=0.62, font_size=30)
    place_in(tab, -6.7, 6.7, -2.9, -0.8)
    assert_in_frame(left, right, tab, mlab)
    assert_apart(left, right)
    assert_column(*left)
    assert_apart(right, tab)
    assert_apart(left, tab)
    assert_table(tab, 6)
    s.play(FadeIn(e35), run_time=0.8)
    s.say("Ex. 9.4：書上 (9.35) 寫 ", "$\\frac{m^4}{12}$", "，就是 ", "$u=\\frac13$", wait=1.5)
    s.play(FadeIn(dF), run_time=0.8)
    s.say("求極小：", "$\\partial F/\\partial m=rm+m^3/3-h=0$", wait=1.5)
    s.play(FadeIn(sol), run_time=0.8)
    s.say("$r=-1$", "：", "$m^*=\\sqrt3=1.7321$", "、", "$F_{\\min}=-0.75$", wait=1.5)
    s.play(FadeIn(coll), run_time=0.8)
    s.say("塌縮：", "$m=|r|^{1/2}\\mathcal M(h/|r|^{3/2})$", wait=1.5)
    s.play(FadeIn(ax), FadeIn(note), FadeIn(leg), run_time=0.6)
    s.play(ShowCreation(c1), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(d) for d in dts], lag_ratio=0.1), run_time=1.0)
    s.say("$r=-0.01$", " 與 ", "$r=-0.0001$", " 疊成同一條曲線", wait=1.5)
    s.play(ShowCreation(guide), FadeIn(mlab), run_time=0.8)
    s.say("$x=1$", " 時 ", "$\\mathcal M=2.1038$", "：兩個 ", "$r$", " 都一樣", wait=1.5)
    s.say("所以 ", "$\\beta=\\frac12$", "、", "$\\Delta=\\frac32$", "；", "$F_{\\min}\\propto|r|^2$", " → ", "$\\alpha=0$",
          wait=1.8)
    s.say("$\\delta=\\Delta/\\beta=3$", "、", "$\\gamma=1$", "：平均場的指數", wait=1.5)
    s.play(FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:]], lag_ratio=0.1), run_time=1.0)
    s.say("數值斜率和理論一樣", wait=1.5)
    s.say("書上 (b) 寫 ", "$r^{2-\\alpha}$", "，應為 ", "$|r|^{2-\\alpha}$", "（", "$r<0$", "）", wait=2.0)
    s.clear_stage()


def part_zmf(s):
    nb = numbers()
    e34 = mixed("(9.34) ", "$Z_{MF}=\\int dm\\,e^{-VF(m)}$", "（", "$n=1,\\ r=-1,\\ u=1$", "）", font_size=28)
    ax = D.plot_axes((0, 4.4, 1), (0, 1.25, 0.5), width=4.6, height=3.0, x_ticks=(0, 1, 2, 3, 4),
                     x_tick_labels=["-2", "-1", "0", "1", "2"], x_label="m", y_label=r"e^{-VF}")
    leg = VGroup(*[mixed("$V=%d$" % V, "，", "$h=%s$" % fg(WC.ZMF_H), font_size=22, color=c)
                   for V, c in zip(WC.ZMF_V[:2], V_COLORS)]).arrange(DOWN, aligned_edge=LEFT, buff=0.15)
    plot = VGroup(ax, leg).arrange(DOWN, buff=0.25)
    hd = [r"V", r"\langle m\rangle_{h=0.01}", r"\langle m\rangle_{h=0}", r"\langle|m|\rangle_{h=0}", r"(-\ln Z/V)_{h=0}"]
    assert fg(WC.ZMF_H) == "0.01"
    rows = [hd] + [[str(V), f4(a), "0", f4(b), f4(c)] for V, a, b, c in zip(WC.ZMF_V, nb["zm_h"], nb["zabs"], nb["zf"])]
    rows.append([zh("谷底", font_size=24, color=GREEN), f4(nb["mstar_h"]), "0", fg(nb["lc"][0]), f2(nb["lc"][1])])
    tab = D.tex_table(rows, [1.4, 2.3, 1.8, 2.0, 2.4], row_h=0.66, font_size=28)
    for j in range(1, 5):
        D.table_cell(tab, 4, j, 5).set_color(GREEN)
    place_in(e34, -6.7, 6.7, 2.45, 3.1)
    place_in(plot, -6.7, -2.3, -2.9, 2.2)
    place_in(tab, -2.0, 6.7, -2.6, 2.0)
    w = [weight_curve(ax, V, WC.ZMF_H, c) for V, c in zip(WC.ZMF_V[:2], V_COLORS)]
    assert_in_frame(e34, plot, tab)
    assert_column(e34, plot)
    assert_column(e34, tab)
    assert_apart(plot, tab)
    assert_table(tab, 5)
    cells = tab[0]
    s.play(FadeIn(e34), run_time=0.8)
    s.say("(9.34)：有限體積 ", "$V$", " 的平均場——對 ", "$m$", " 積分", wait=1.5)
    s.play(FadeIn(ax), FadeIn(leg[0]), run_time=0.6)
    s.play(ShowCreation(w[0]), run_time=1.0)
    s.say("加一點外場 ", "$h=0.01$", "；", "$V=10$", " 時兩個峰差不多高", wait=1.5)
    s.play(FadeIn(cells[:5]), ShowCreation(tab[1]), FadeIn(cells[5:10]), run_time=0.8)
    s.say("$V=10$", "：", "$\\langle m\\rangle=0.0869$", "，幾乎是 0", wait=1.4)
    s.play(FadeIn(leg[1]), ShowCreation(w[1]), FadeIn(cells[10:15]), run_time=1.0)
    s.say("$V=100$", "：左邊的峰矮了，", "$\\langle m\\rangle=0.7576$", wait=1.4)
    s.play(FadeIn(cells[15:20]), run_time=0.6)
    s.say("$V=1000$", "：", "$\\langle m\\rangle=1.0042$", "，接近谷底 ", "$m^*=1.0050$", wait=1.5)
    s.play(FadeIn(cells[20:25]), run_time=0.6)
    s.say("$h=0$", "：任何 ", "$V$", " 都是 ", "$\\langle m\\rangle=0$", "——兩邊一樣多", wait=1.5)
    s.say("但 ", "$\\langle|m|\\rangle\\to1$", "，", "$-\\ln Z/V\\to F_{\\min}=-0.25$", wait=1.6)
    s.say("先 ", "$V\\to\\infty$", "、再 ", "$h\\to0$", "，才得到 ", "$m^*\\ne0$", wait=1.6)
    s.say("和 W902 的翻轉等待是同一件事", wait=1.5)
    s.say("GLW 把自旋變成可伸縮的場；只留 ", "$\\vec q=0$", " 就是平均場", wait=2.0)
    s.say("雙井給出 ", "$\\beta=\\frac12$", "、", "$\\delta=3$", "、", "$\\alpha=0$", wait=2.0)


class W904LandauWilson(ChuScene):
    def construct(self):
        self.show_header("W904　Ginzburg-Landau-Wilson：讓自旋可以伸縮")
        check_math()
        part_glw(self)
        part_fourier(self)
        part_well(self)
        part_ex94(self)
        part_zmf(self)
        self.wait(1.0)
