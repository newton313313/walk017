"""
W607 主半徑的分布、Fixman、自迴避＋第 6 章回顧（§6.2.6–6.3，pp. 152–154，(6.72)–(6.77)，Fig. 6.15，Ex. 6.8–6.9）

問：三維漫步最大的 R_1² 分布長什麼樣子？三個加起來（迴轉半徑的平方）呢？
  1. (6.72)(6.73)：R_i² 近似 gamma 分布（SH.gamma_density；Ex. 6.8 變異數 2α²/d）；Fig. 6.15：SH.eig_sim(3)[0][:, 0] 的直方圖
     （30 格、寬 0.02）vs α_1 = 1/π² 的曲線（峰 0.0338、高 7.16）；實際平均 0.1280，(6.61) 的 0.1267。
  2. Fixman (6.74)–(6.76)：SH.fixman_density（120 點）vs SH.eig_sim(3)[1] 的直方圖；峰 0.111（6.20）；x = 0.1, 0.2, 0.3：6.078、3.281、1.011；
     平均 1/6、變異數 1/135（SH.fixman_moments_numeric、SH.fixman_mean_var）；筆誤框 (6.76) 的號（book=True 全是 0）。
  3. 自迴避 (6.77)（SH.A_saw）：0.5263 + 0.008 = 0.5343。
  4. 第 6 章回顧表。
畫面上的數字都由 utils.walk_shape 算出，check_math() 核對。
"""
import math
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_shape as SH
from utils import walk_draw as D

ALPHA1 = 1 / math.pi**2
BINS, BW = 30, 0.02
XS_FIX = np.linspace(0.005, 0.6, 120)
XS_AT = np.array([0.1, 0.2, 0.3])
HIST_COLOR = BLUE_C
CURVE_COLOR = YELLOW


def f4(x):
    return f"{x:.4f}"


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    top, tr = SH.eig_sim(3)
    nb = {"r1": top[:, 0], "tr": tr}
    nb["mode"] = 1 / (3 * math.pi**2)
    nb["peak"] = float(SH.gamma_density(nb["mode"], ALPHA1, 3))
    nb["mean1"] = float(top[:, 0].mean())
    nb["pred1"] = float(SH.onedee_table(3)["pred_mean"][0])
    nb["fix"] = SH.fixman_density(XS_FIX, 3)
    nb["fix_book"] = SH.fixman_density(XS_FIX, 3, book=True)
    nb["fix_at"] = SH.fixman_density(XS_AT, 3)
    xg = np.round(np.arange(0.100, 0.1205, 0.001), 3)                     # 峰附近（間隔 0.001）
    pg = SH.fixman_density(xg, 3)
    nb["fix_peak"] = (float(xg[np.argmax(pg)]), float(pg.max()))
    nb["fix_mv"] = SH.fixman_mean_var(3)
    nb["fix_num"] = SH.fixman_moments_numeric(3)
    nb["fix_N"] = SH.fixman_mean_var(3, SH.N_D)
    nb["fix_sim"] = (float(tr.mean()), float(tr.var()))
    nb["saw"] = SH.A_saw(3)
    nb["A3"] = SH.A_d(3)
    nb["Ap3"] = SH.Aprime(3)
    return nb


def trap(y, x):
    return float(np.sum((y[1:] + y[:-1]) * np.diff(x)) / 2)


def density(values):
    cnt, edges = np.histogram(values, bins=BINS, range=(0, BINS * BW))
    return edges, cnt / (len(values) * BW)


def check_math():
    nb = numbers()
    assert f4(nb["mode"]) == "0.0338" and f"{nb['peak']:.2f}" == "7.16" and f4(ALPHA1) == "0.1013"
    x = np.linspace(1e-6, 3, 300001)                                       # Ex. 6.8：平均 α、變異數 2α²/d
    g = SH.gamma_density(x, ALPHA1, 3)
    m0, m1 = trap(g, x), trap(x * g, x)
    assert abs(m0 - 1) < 1e-3 and abs(m1 - ALPHA1) < 1e-3
    assert abs(trap((x - m1) ** 2 * g, x) - 2 * ALPHA1**2 / 3) < 1e-5
    assert f4(nb["mean1"]) == "0.1280" and f4(nb["pred1"]) == "0.1267"
    assert abs(nb["pred1"] - nb["mean1"]) < abs(ALPHA1 - nb["mean1"])
    e1, d1 = density(nb["r1"])
    assert len(d1) == 30 and abs(e1[1] - e1[0] - 0.02) < 1e-12 and d1.max() < 7.6
    assert nb["fix_peak"][0] == 0.111 and f"{nb['fix_peak'][1]:.2f}" == "6.20"
    assert 0.1 < nb["fix_peak"][0] < 0.12 and nb["fix_peak"][1] >= nb["fix"].max() - 1e-9
    assert [f"{v:.3f}" for v in nb["fix_at"]] == ["6.078", "3.281", "1.011"]
    assert abs(nb["fix_mv"][0] - 1 / 6) < 1e-15 and abs(nb["fix_mv"][1] - 1 / 135) < 1e-15
    assert (f4(1 / 6), f4(1 / 135), f4(math.sqrt(1 / 135))) == ("0.1667", "0.0074", "0.0861")
    assert [f4(v) for v in nb["fix_num"]] == ["1.0000", "0.1667", "0.0074"]
    assert [f4(v) for v in nb["fix_sim"]] == ["0.1671", "0.0076"]
    assert [f4(v) for v in nb["fix_N"]] == ["0.1666", "0.0074"]
    assert np.max(np.abs(nb["fix_book"])) < 1e-5                          # 書上的號：R² > 0 全是 0
    e2, d2 = density(nb["tr"])
    assert np.max(np.abs(d2[4:15] - SH.fixman_density(e2[4:15] + 0.01, 3))) < 0.4
    assert len(XS_FIX) <= 120 and nb["fix"].max() < 7
    assert f4(nb["saw"]) == "0.5343" and f4(float(nb["A3"])) == "0.5263" and abs(nb["saw"] - float(nb["A3"]) - 0.008) < 1e-12
    assert f4(nb["Ap3"]) == "0.3943" and nb["A3"] == Fraction(10, 19)
    assert abs(float(SH.staircase_asph(10**6)) - 0.4) < 1e-6 and SH.A_d(10**9) - Fraction(2, 5) < 1e-8


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


def eq_label(tex_mob, label, buff=0.35):
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


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def hist_bars(ax, values, color=HIST_COLOR):
    edges, dens = density(values)
    bars = VGroup()
    for i in range(BINS):
        p0, p1 = ax.c2p(edges[i], 0), ax.c2p(edges[i + 1], dens[i])
        r = Rectangle(width=p1[0] - p0[0], height=max(p1[1] - p0[1], 1e-3))
        r.set_fill(color, 0.55).set_stroke(color, 1).move_to((p0 + p1) / 2)
        bars.add(r)
    return bars


def xy_curve(ax, xs, ys, color, width=4):
    return VMobject().set_points_as_corners([ax.c2p(x, y) for x, y in zip(xs, ys)]).set_stroke(color, width)


def dist_axes(ymax, x_label):
    return D.plot_axes((0, 0.6, 0.1), (0, ymax, 1), width=6.2, height=4.4, x_ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                       y_ticks=list(range(2, int(ymax) + 1, 2)), font_size=20, x_label=x_label)


def vline(ax, x, ytop, color):
    return DashedLine(ax.c2p(x, 0), ax.c2p(x, ytop), dash_length=0.08).set_stroke(color, 2.5)


# ---------------------------------------------------------------- 各段
def part_question(s):
    q = VGroup(mixed("三維、高斯步伐、", "$N=50$", "，2 萬條（W606 的模擬）", font_size=30),
               mixed("最大的 ", "$R_1^2$", "：分布長什麼樣子？", font_size=34, color=YELLOW),
               mixed("三個加起來 ", "$R_g^2=R_1^2+R_2^2+R_3^2$", "：呢？", font_size=34, color=YELLOW)
               ).arrange(DOWN, buff=0.8)
    place_in(q, -6.5, 6.5, -2.5, 2.8)
    assert_in_frame(q)
    assert_column(*q)
    s.play(FadeIn(q[0]), run_time=0.8)
    s.say("W606 看了平均與變異數；這一幕看整個分布", wait=2.5)
    s.play(FadeIn(q[1]), run_time=0.8)
    s.say("問：三維漫步最大的 ", "$R_1^2$", " 分布長什麼樣子？", wait=2.5)
    s.play(FadeIn(q[2]), run_time=0.8)
    s.say("三個加起來（迴轉半徑的平方）呢？（先猜）", wait=3.5)
    s.clear_stage()


def part_gamma(s):
    nb = numbers()
    ax = dist_axes(8, r"R_1^2/(N+1)")
    bars = hist_bars(ax, nb["r1"])
    cur = ax.get_graph(lambda x: float(SH.gamma_density(x, ALPHA1, 3)), x_range=(0.001, 0.6, 0.002))
    cur.set_stroke(CURVE_COLOR, 4)
    dot = Dot(ax.c2p(nb["mode"], nb["peak"]), radius=0.07).set_color(RED)
    l_a = vline(ax, ALPHA1, 7.6, CURVE_COLOR)
    l_m = vline(ax, nb["mean1"], 7.6, TEAL_A)
    plot = VGroup(ax, bars, cur, dot, l_a, l_m)
    f1 = eq_label(tx(r"P(R_i^2)\propto(R_i^2)^{d/2-1}\,e^{-dR_i^2/2\alpha_i}", 34), "(6.72)(6.73)")
    f2 = mixed("Ex. 6.8：平均 ", "$\\alpha_i$", "，變異數 ", "$\\frac{2\\alpha_i^2}{d}$", font_size=26)
    f3 = mixed("曲線：", "$\\alpha_1=\\frac{1}{\\pi^2}=%s$" % f4(ALPHA1), font_size=26, color=CURVE_COLOR)
    f4l = mixed("峰在 ", "$%s$" % f4(nb["mode"]), "，高 ", "$%.2f$" % nb["peak"], "（紅點）", font_size=26)
    f5 = mixed("實際平均 ", "$%s$" % f4(nb["mean1"]), "（青線）", font_size=26, color=TEAL_A)
    f6 = mixed("(6.61) 的 ", "$%s$" % f4(nb["pred1"]), " 比較準", font_size=26, color=GREEN_B)
    right = VGroup(f1, f2, f3, f4l, f5, f6).arrange(DOWN, buff=0.38, aligned_edge=LEFT)
    place_in(plot, -6.75, 0.2, -2.9, 3.1)
    place_in(right, 0.6, 6.75, -2.9, 3.1)
    assert_in_frame(plot, right)
    assert_apart(plot, right)
    assert_column(*right)
    s.play(FadeIn(f1), run_time=0.8)
    s.say("近似：每個 ", "$R_i^2$", " 各自是一個 gamma 分布", wait=2.5)
    s.play(FadeIn(f2), run_time=0.6)
    s.say("平均 ", "$\\alpha_i$", "、變異數 ", "$2\\alpha_i^2/d$", "：維度越高越窄", wait=2.5)
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.05), run_time=1.5)
    s.say("長條：模擬的 ", "$R_1^2/(N+1)$", "（30 格，寬 0.02）", wait=2.0)
    s.play(ShowCreation(cur), FadeIn(f3), run_time=1.2)
    s.play(FadeIn(dot), FadeIn(f4l), run_time=0.6)
    s.say("曲線：", "$\\alpha_1=1/\\pi^2$", "，峰在 0.0338、高 7.16", wait=2.5)
    s.play(ShowCreation(l_a), ShowCreation(l_m), FadeIn(f5), run_time=0.8)
    s.say("形狀大致對，但曲線偏左：實際平均 0.1280", wait=3.0)
    s.say("零階的 ", "$\\alpha_1=0.1013$", " 少算了 1/d 的位移", wait=2.0)
    s.play(FadeIn(f6), run_time=0.6)
    s.say("用 1/d 修正過的 (6.61)：0.1267，準多了", wait=3.0)
    s.clear_stage()


def part_fixman(s):
    nb = numbers()
    ax = dist_axes(7, r"R_g^2/(N+1)")
    bars = hist_bars(ax, nb["tr"])
    cur = xy_curve(ax, XS_FIX, nb["fix"], CURVE_COLOR)
    pk = Dot(ax.c2p(*nb["fix_peak"]), radius=0.07).set_color(RED)
    dots = VGroup(*[Dot(ax.c2p(x, y), radius=0.06).set_color(WHITE) for x, y in zip(XS_AT, nb["fix_at"])])
    book = xy_curve(ax, XS_FIX, nb["fix_book"], RED, width=6)
    plot = VGroup(ax, bars, cur, pk, dots, book)
    k1 = mixed("$R_g^2=\\sum_nR_n^2$", "：獨立的量相加", font_size=26)
    k2 = eq_label(tx(r"K(s)=\langle e^{-isR_g^2}\rangle=\prod_n\Big(1+\frac{2is\alpha_n}{d}\Big)^{-d/2}", 30),
                  "(6.76)")
    k3 = eq_label(tx(r"P(R_g^2)=\frac{1}{2\pi}\int ds\,e^{isR_g^2}K(s)", 30), "(6.74)(6.75)")
    k4 = mixed("Ex. 6.9：", "$K$", " 是各自的乘積；", "$\\alpha_n=\\frac{1}{\\pi^2n^2}$", font_size=24, color=GREY_A)
    k5 = mixed("峰在 ", "$%.3f$" % nb["fix_peak"][0], "，高 ", "$%.2f$" % nb["fix_peak"][1], "（紅點）",
               font_size=24, color=CURVE_COLOR)
    k6 = mixed("$x=0.1,\\,0.2,\\,0.3$", "（白點）：", "$%s$" % r",\ ".join(f"{v:.3f}" for v in nb["fix_at"]),
               font_size=24, color=CURVE_COLOR)
    right1 = VGroup(k1, k2, k3, k4, k5, k6).arrange(DOWN, buff=0.38, aligned_edge=LEFT)
    place_in(plot, -6.75, 0.0, -2.9, 3.1)
    place_in(right1, 0.4, 6.75, -2.9, 3.1)
    m1 = mixed("平均 ", "$\\frac16=%s$" % f4(nb["fix_mv"][0]), "，變異數 ", "$\\frac{1}{135}=%s$" % f4(nb["fix_mv"][1]),
               font_size=26, color=CURVE_COLOR)
    m2 = mixed("（標準差 ", "$%s$" % f4(math.sqrt(nb["fix_mv"][1])), "）", font_size=24, color=GREY_A)
    m3 = mixed("數值積分：", "$%s$" % r",\ ".join(f4(v) for v in nb["fix_num"]), font_size=26)
    m4 = mixed("模擬：", "$%s$" % r",\ ".join(f4(v) for v in nb["fix_sim"]), font_size=26, color=HIST_COLOR)
    m5 = mixed("$N=50$", " 精確：", "$%s$" % r",\ ".join(f4(v) for v in nb["fix_N"]), font_size=26)
    box = D.erratum_box("p. 153，(6.76)", r"K(R)=\prod_{n}\Big(1-\frac{2is\alpha_n}{d}\Big)^{-d/2}",
                        r"K(s)=\prod_{n}\Big(1+\frac{2is\alpha_n}{d}\Big)^{-d/2}", width=6.5)
    right2 = VGroup(m1, m2, m3, m4, m5).arrange(DOWN, buff=0.32, aligned_edge=LEFT)
    place_in(right2, 0.4, 6.75, -0.6, 3.1)
    place_in(box, 0.2, 6.75, -2.95, -0.8)
    assert_in_frame(plot, right1, right2, box)
    assert_apart(plot, right1)
    assert_apart(plot, right2)
    assert_apart(plot, box)
    assert_apart(right2, box)
    assert_column(*right1)
    assert_column(*right2)
    s.play(FadeIn(k1), run_time=0.6)
    s.say("Fixman：把 ", "$R_g^2$", " 看成許多獨立的 ", "$R_n^2$", " 相加", wait=2.5)
    s.play(FadeIn(k2), FadeIn(k4), run_time=0.8)
    s.say("獨立的量相加：特徵函數 ", "$K$", " 是各自的乘積（Ex. 6.9）", wait=2.5)
    s.play(FadeIn(k3), run_time=0.6)
    s.say("再做一次 Fourier 反轉，就得到分布", wait=2.0)
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.05), run_time=1.5)
    s.say("長條：模擬的 ", "$R_g^2/(N+1)$", "；曲線：Fixman 的精確分布", wait=1.0)
    s.play(ShowCreation(cur), run_time=1.5)
    s.play(FadeIn(pk), FadeIn(dots), FadeIn(k5), FadeIn(k6), run_time=0.6)
    s.say("峰在 0.111、高 6.20；", "$x=0.1$", " 時 6.078", wait=2.5)
    s.play(FadeOut(right1), run_time=0.5)
    s.play(FadeIn(right2[:2]), run_time=0.6)
    s.say("平均 ", "$\\frac16$", "、變異數 ", "$\\frac{1}{135}$", "：兩個無窮和", wait=2.5)
    s.play(FadeIn(right2[2:]), run_time=0.8)
    s.say("數值積分、模擬、", "$N=50$", " 的精確值：都對得上", wait=2.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (6.76)：", "$K(R)$", " 應為 ", "$K(s)$", "，而且號錯了", wait=2.0)
    s.play(ShowCreation(book), run_time=0.8)
    s.say("照書上的號，", "$R^2>0$", " 的地方算出來全是 0（紅線）", wait=3.5)
    s.say("分布全跑到負的 ", "$R^2$", "：號一錯，就不是機率了", wait=2.5)
    s.clear_stage()


def part_saw(s):
    nb = numbers()
    a1 = eq_label(tx(r"A_d=\frac{2d+4}{5d+4}+0.008\,(4-d)", 40), "(6.77)")
    a2 = mixed("三維：", "$%s+0.008=%s$" % (f4(float(nb["A3"])), f4(nb["saw"])), font_size=32, color=YELLOW)
    scale = 9.0
    rows = VGroup()
    for i, (name, val, c) in enumerate((("一般漫步", float(nb["A3"]), BLUE_C), ("自迴避", nb["saw"], ORANGE))):
        bar = Rectangle(width=val * scale, height=0.4).set_fill(c, 0.7).set_stroke(c, 1.5)
        bar.move_to(np.array([val * scale / 2, -0.75 * i, 0]))             # 左端都在 x = 0
        lab = zh(name, font_size=26, color=c).next_to(bar.get_left(), LEFT, buff=0.3)
        num = tx(f4(val), 26, c).next_to(bar, RIGHT, buff=0.3)
        rows.add(VGroup(lab, bar, num))
    col = VGroup(a1, a2, rows).arrange(DOWN, buff=0.6)
    place_in(col, -6.5, 6.5, -2.7, 2.9)
    assert_in_frame(col)
    assert_column(*col)
    assert_apart(rows[0], rows[1])
    s.play(FadeIn(a1), run_time=0.8)
    s.say("自迴避漫步：不准走回走過的地方（只引用）", wait=2.5)
    s.play(FadeIn(a2), run_time=0.6)
    s.say("三維：", "$0.5263+0.008=0.5343$", wait=2.0)
    s.play(*[FadeIn(r) for r in rows], run_time=0.8)
    s.say("自迴避讓漫步稍微更扁，但影響不大", wait=3.0)
    s.say("一般漫步的結果已經抓住了大部分", wait=2.5)
    s.clear_stage()


def part_review(s):
    nb = numbers()
    hz = lambda t: zh(t, font_size=26, color=YELLOW)
    mx = lambda *p: mixed(*p, font_size=26)
    sc = lambda t: tx(r"\text{" + t + "}", 26, GREY_A)
    rows = [[hz("形狀"), mx("$\\mathbf{T}$", " 的特徵值 ", "$R_i^2$", "（主迴轉半徑的平方）"), sc("W601--W602")],
            [hz("非球度"), mx("$\\frac{4+2d}{4+5d}$", "；三維 ", "$\\frac{10}{19}$", "，個別的平均 ",
                            "$%s$" % f4(nb["Ap3"])), sc("W603")],
            [mixed("$d\\to\\infty$", font_size=26, color=YELLOW),
             mx("樓梯：", "$R_n^2=\\frac{N+1}{\\pi^2n^2}$", "，非球度 ", "$\\frac25$"), sc("W604")],
            [hz("步伐"), mx("Kramers 恆等式＋Wick 配對，得到 ", "$\\frac{10}{19}$"), sc("W605")],
            [mixed("$1/d$", font_size=26, color=YELLOW), mx("位移 ", "$\\frac{3}{4d}$", "，寬度 ", "$\\sqrt{2/d}$"),
             sc("W606")],
            [hz("分布"), mx("gamma 近似；Fixman 的精確分布（平均 ", "$\\frac16$", "）"), sc("W607")]]
    tab = D.tex_table(rows, [2.0, 8.6, 2.6], row_h=0.9, font_size=26)
    place_in(tab, -6.75, 6.75, -2.8, 2.9)
    assert_in_frame(tab)
    assert_table(tab, 3)
    s.play(ShowCreation(tab[1][0]), ShowCreation(tab[1][2]), run_time=0.6)
    s.say("第 6 章回顧：漫步的形狀", wait=1.0)
    for i in range(6):
        s.play(FadeIn(tab[0][3 * i:3 * i + 3], shift=0.15 * UP), run_time=0.6)
        s.wait(2.0)
    s.say("漫步沒有「標準形狀」，但形狀有統計", wait=3.0)
    s.say("三維平均一半扁，維度越高越確定", wait=3.5)


class W607FixmanSaw(ChuScene):
    def construct(self):
        self.show_header("W607　主半徑的分布、Fixman、自迴避")
        check_math()
        part_question(self)
        part_gamma(self)
        part_fixman(self)
        part_saw(self)
        part_review(self)
        self.wait(1.0)
