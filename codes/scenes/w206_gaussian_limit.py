"""
W206 高斯極限與奇偶性的因子 2（§2.1.3，pp. 34–35，(2.18)–(2.20)）

問：N = 10，用高斯近似估立方格子回到原點的走法數（精確 1172556），會差多少？
  1. (2.18)(2.19)：(cos q_x a + cos q_y a + cos q_z a)^N = exp[N ln(3 − a²q²/2 + …)] = 3^N e^{−Na²q²/6 + …}；
     沿 q_x 軸畫 ((cos q + 2)/3)^10 與 e^{−10q²/6}：幾乎重合。
  2. (2.20) C(N;r) ≈ 6^N (3/2πN)^{3/2} e^{−3r²/2Na²}；N = 10、r = 0：630848；精確 1172556 是它的 1.86 倍——少了一半！
  3. 另一半在布里淵區的角落 q = (π,π,π)/a：χ = −6，(−6)^10 = 6^10；沿體對角線 cos^10 q 兩端都是 1；
     角落的高斯相位 (−1)^{N+x+y+z}：同奇偶 ×2、不同 → 0。注記框（p. 34，(2.20) 前）；更正後 1261696，精確／近似 0.929（N = 20：0.964）。
  4. 正方格子（主例）2·4^N/(πN)：N = 10 → 66754（精確 63504；沒乘 2 的 33377）；比值 0.9513, 0.9753, 0.9901, 0.995 → 1。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import boxed, mixed, zh
from utils import walk_draw as D
from utils import walk_gf as GF

SQ_NS = [10, 20, 50, 100]
CENTER_COLOR = D.GAUSS_COLOR
CORNER_COLOR = TEAL
Q1 = (f"{((math.cos(1) + 2) / 3) ** 10:.4f}", f"{math.exp(-10 / 6):.4f}")      # q = 1 時的兩條曲線


def numbers():
    g10, g20 = GF.gauss_sc(10), GF.gauss_sc(20)
    return {"exact": GF.returns_sc(10), "gauss": round(g10), "ratio": f"{GF.returns_sc(10) / g10:.2f}",
            "p6": 6 ** 10, "coef": f"{(3 / (20 * math.pi)) ** 1.5:.6f}",
            "two": round(2 * g10), "r10": f"{GF.returns_sc(10) / (2 * g10):.3f}",
            "r20": f"{GF.returns_sc(20) / (2 * g20):.3f}",
            "sq_asym": round(GF.asym_returns_sq(10)), "sq_half": round(GF.gauss_sq(10)), "sq_exact": GF.returns_sq(10),
            "sq_ratios": [round(GF.returns_sq(N) / GF.asym_returns_sq(N), 4) for N in SQ_NS]}


def check_math():
    nb = numbers()
    assert nb["exact"] == 1172556 and nb["gauss"] == 630848 and nb["ratio"] == "1.86"
    assert nb["p6"] == 60466176 and nb["coef"] == "0.010433" and Q1 == ("0.1895", "0.1889")
    assert abs(nb["p6"] * (3 / (20 * math.pi)) ** 1.5 - GF.gauss_sc(10)) < 1e-6
    assert nb["two"] == 1261696 == 2 * nb["gauss"] and nb["r10"] == "0.929" and nb["r20"] == "0.964"
    assert abs(GF.gauss_lattice(10, "sc", 0.0, 0) - 2 * GF.gauss_sc(10)) < 1e-6
    assert GF.gauss_lattice(10, "sc", 1.0, 1) == 0 and abs(GF.fourier_count_sc(10, (1, 0, 0))) < 1e-6   # 奇偶不同 → 0
    assert nb["sq_asym"] == 66754 and nb["sq_half"] == 33377 and nb["sq_exact"] == 63504
    assert nb["sq_ratios"] == GF.SQ_RATIOS == [0.9513, 0.9753, 0.9901, 0.995]
    assert abs(GF.chi_sc(math.pi, math.pi, math.pi) + 6) < 1e-12 and (-6) ** 10 == 6 ** 10
    qs = [-math.pi + k * 2 * math.pi / 400 for k in range(401)]
    assert max(abs(((math.cos(q) + 2) / 3) ** 10 - math.exp(-10 * q * q / 6)) for q in qs) < 0.03       # 幾乎重合
    for q in (0.0, 0.7, 2.2, math.pi):                                       # 沿體對角線：(χ/6)^10 = cos^10 q
        assert abs((GF.chi_sc(q, q, q) / 6) ** 10 - math.cos(q) ** 10) < 1e-12
    assert abs(math.cos(math.pi) ** 10 - 1) < 1e-12


# ---------------------------------------------------------------- 小工具
def eq(tex, tag=None, font_size=32, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return m
    return VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


def dashed_graph(ax, fn, a, b, color, n=36):
    return DashedVMobject(ax.get_graph(fn, x_range=[a, b, 0.02]), num_dashes=n).set_stroke(color, 3)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    e1 = Tex(rf"C(10;\vec 0)={nb['exact']}", font_size=48)
    e2 = mixed("簡單立方格子，精確值（W205）", font_size=26, color=GREY_B)
    e3 = mixed("高斯近似 ", "$\\approx\\ ?$", font_size=38, color=CENTER_COLOR)
    col = VGroup(e1, e2, e3).arrange(DOWN, buff=0.45).move_to(np.array([0, 0.4, 0]))
    assert_in_frame(col)
    s.play(FadeIn(e1), FadeIn(e2), run_time=1.0)
    s.say("問：", "$N=10$", "，用高斯近似估立方格子回到原點的走法數", wait=0.3)
    s.play(FadeIn(e3, shift=0.2 * UP), run_time=0.8)
    s.wait(2.0)
    s.say("精確是 1172556：高斯近似會差多少？", wait=3.5)
    s.clear_stage()


def part_expand(s):
    a1 = eq(r"(\cos q_xa+\cos q_ya+\cos q_za)^N=\exp\Big[N\ln\Big(3-\frac12a^2q^2+O(q^4)\Big)\Big]", "(2.18)",
            font_size=34)
    a3 = mixed("每個 ", "$\\cos qa\\approx1-\\frac12q^2a^2$", "，三個加起來 ", "$3-\\frac12a^2q^2$", "；再用 ",
               "$\\ln(3-\\epsilon)\\approx\\ln3-\\frac{\\epsilon}{3}$", font_size=22, color=GREY_B)
    a2 = eq(r"=3^N\,e^{-\frac{N}{6}a^2|\vec q|^2+\cdots}", "(2.19)", font_size=38)
    top = VGroup(fit(a1, 12.6), fit(a3, 12.6), a2).arrange(DOWN, buff=0.22).move_to(np.array([0, 2.0, 0]))
    ax = D.plot_axes((-3.5, 3.5, 1), (0, 1.15, 0.5), width=7.0, height=2.6, x_ticks=[-PI, 0, PI],
                     x_tick_labels=[r"-\pi", "0", r"\pi"], y_ticks=[1], font_size=22, x_label="q_x")
    ax.move_to(np.array([-1.9, -1.15, 0]))
    g1 = ax.get_graph(lambda q: ((math.cos(q) + 2) / 3) ** 10, x_range=[-PI, PI, 0.01]).set_stroke(D.BAR_COLOR, 4)
    g2 = dashed_graph(ax, lambda q: math.exp(-10 * q * q / 6), -PI, PI, CENTER_COLOR, n=40)
    leg = VGroup(
        VGroup(Line(ORIGIN, 0.5 * RIGHT).set_stroke(D.BAR_COLOR, 4),
               Tex(r"\Big(\frac{\cos q+2}{3}\Big)^{10}", font_size=30)).arrange(RIGHT, buff=0.2),
        VGroup(DashedLine(ORIGIN, 0.5 * RIGHT, dash_length=0.08).set_stroke(CENTER_COLOR, 4),
               Tex(r"e^{-10q^2/6}", font_size=30)).arrange(RIGHT, buff=0.2),
        mixed("（沿 ", "$q_x$", " 軸，", "$a=1$", "）", font_size=20, color=GREY_B),
        mixed("$q=1$", "：", f"${Q1[0]}$", " 對 ", f"${Q1[1]}$", font_size=22)).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    place_col(leg, 2.4, -0.15)
    assert_in_frame(top, ax, leg)
    assert_apart(top, ax)
    assert_apart(ax, leg)
    assert_apart(top, leg)

    s.say("立方格子：小 ", "$q$", " 時把 ", "$\\chi^N$", " 寫成指數 (2.18)", wait=0.3)
    s.play(Write(a1[0]), FadeIn(a1[1]), run_time=2.0)
    s.wait(1.5)
    s.play(FadeIn(a3), run_time=0.8)
    s.say("每個 ", "$\\cos$", " 展開到 ", "$q^2$", "：三個加起來 ", "$3-\\frac12a^2q^2$", wait=3.0)
    s.play(FadeIn(a2), run_time=0.8)
    s.say("$N\\ln3-\\frac{N}{6}a^2q^2$", "：", "$3^N$", " 乘上一個高斯 (2.19)", wait=3.0)
    s.play(FadeIn(ax), ShowCreation(g1), FadeIn(leg[0]), run_time=1.5)
    s.play(ShowCreation(g2), FadeIn(leg[1]), FadeIn(leg[2]), run_time=1.2)
    s.say("$N=10$", "：沿 ", "$q_x$", " 軸，兩條曲線幾乎重合", wait=2.5)
    s.play(FadeIn(leg[3]), run_time=0.6)
    s.say("例如 ", "$q=1$", "：", f"${Q1[0]}$", " 對 ", f"${Q1[1]}$", wait=3.0)
    s.clear_stage()


def part_estimate(s):
    nb = numbers()
    b1 = boxed(eq(r"C(N;\vec r)\approx6^N\Big(\frac{3}{2\pi N}\Big)^{3/2}e^{-3|\vec r|^2/2Na^2}", "(2.20)", font_size=40),
               color=CENTER_COLOR, buff=0.2)
    b2 = Tex(r"N=10,\ \vec r=\vec 0:\ \ 6^{10}\Big(\frac{3}{20\pi}\Big)^{3/2}", font_size=40)
    b2b = Tex(rf"={nb['p6']}\times{nb['coef']}\approx{nb['gauss']}", font_size=40)
    b3 = mixed("精確：", f"${nb['exact']}$", font_size=30)
    b4 = VGroup(Tex(rf"\frac{{{nb['exact']}}}{{{nb['gauss']}}}\approx{nb['ratio']}", font_size=44),
                zh("少了一半！", font_size=36, color=RED)).arrange(RIGHT, buff=0.5)
    col = VGroup(b1, b2, b2b, b3, b4).arrange(DOWN, buff=0.4).move_to(np.array([0, 0.3, 0]))
    assert_in_frame(col)
    s.say("對 ", "$\\vec q$", " 積分（高斯積分）得到 (2.20)", wait=0.3)
    s.play(FadeIn(b1, shift=0.2 * UP), run_time=1.2)
    s.wait(2.5)
    s.play(FadeIn(b2), run_time=1.0)
    s.say("$N=10$", "、", "$\\vec r=\\vec 0$", "：代進去", wait=2.0)
    s.play(FadeIn(b2b), run_time=0.8)
    s.say("高斯近似：回到原點約 ", f"${nb['gauss']}$", " 種", wait=2.5)
    s.play(FadeIn(b3), run_time=0.8)
    s.say("先猜：精確值 1172556 是它的幾倍？", wait=3.0)
    s.play(FadeIn(b4[0]), run_time=0.8)
    s.play(FadeIn(b4[1], scale=1.3), run_time=0.6)
    s.say(f"{nb['ratio']} 倍：高斯近似只拿到一半！", wait=3.0)
    s.clear_stage()


def part_corner(s):
    nb = numbers()
    ax = D.plot_axes((0, 3.45, 1), (0, 1.15, 0.5), width=4.8, height=2.3, x_ticks=[0, PI], x_tick_labels=["0", r"\pi"],
                     y_ticks=[1], font_size=22, x_label="q")
    ax.move_to(np.array([-3.85, 1.55, 0]))
    gr = ax.get_graph(lambda q: math.cos(q) ** 10, x_range=[0, PI, 0.01]).set_stroke(D.BAR_COLOR, 4)
    ends = VGroup(Dot(ax.c2p(0, 1), radius=0.09).set_fill(CENTER_COLOR, 1),
                  Dot(ax.c2p(PI, 1), radius=0.09).set_fill(CORNER_COLOR, 1))
    elab = VGroup(zh("中心", font_size=20, color=CENTER_COLOR).next_to(ends[0], UR, buff=0.06),
                  zh("角落", font_size=20, color=CORNER_COLOR).next_to(ends[1], UP, buff=0.1))
    glab = Tex(r"\cos^{10}q", font_size=28).set_color(D.BAR_COLOR).move_to(ax.c2p(1.7, 0.45))
    under = mixed("沿體對角線 ", "$\\vec q=q(1,1,1)$", "：", "$\\chi/6=\\cos q$", font_size=22, color=GREY_B)
    under.next_to(ax, DOWN, buff=0.12)
    c1 = mixed("角落 ", "$\\vec q=(\\pi,\\pi,\\pi)/a$", "：", "$\\chi=-6$", font_size=26)
    c2 = Tex(r"(-6)^{10}=6^{10}", font_size=38).set_color(CORNER_COLOR)
    c3 = mixed("角落的高斯也貢獻一份，相位 ", "$(-1)^{N+x+y+z}$", font_size=24)
    c4 = mixed("同奇偶 ", "$\\Rightarrow\\times2$", "；不同 ", "$\\Rightarrow0$", "（棋盤格）", font_size=24)
    rcol = VGroup(c1, c2, c3, c4)
    for m in rcol:
        fit(m, 6.4)
    place_col(rcol.arrange(DOWN, aligned_edge=LEFT, buff=0.32), 0.15, 2.85)
    err = D.erratum_box("p. 34，(2.20) 前", r"\text{only the zone center}", r"\text{center}+\text{corner}:\ \times2",
                        width=6.2)
    place_col(err, -6.5, -1.2)
    d1 = Tex(rf"2\times{nb['gauss']}={nb['two']}", font_size=38).set_color(D.EXACT_COLOR)
    d2 = mixed("精確／近似 ", f"$={nb['r10']}$", font_size=26)
    d3 = mixed("$N=20$", "：", f"${nb['r20']}$", "（越來越接近 1）", font_size=26)
    dcol = VGroup(d1, d2, d3).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    place_col(dcol, 0.6, -1.0)
    assert_in_frame(ax, under, rcol, err, dcol, elab)
    assert_apart(VGroup(ax, under), rcol)
    assert_apart(under, err)
    assert_apart(err, dcol)
    assert_apart(rcol, dcol)

    s.say("另一半在哪裡？布里淵區的角落！", wait=0.3)
    s.play(FadeIn(ax), ShowCreation(gr), FadeIn(glab), FadeIn(under), run_time=1.5)
    s.play(FadeIn(ends), FadeIn(elab), run_time=0.8)
    s.say("沿體對角線畫 ", "$\\cos^{10}q$", "：兩端都是 1", wait=3.0)
    s.play(FadeIn(c1), FadeIn(c2), run_time=1.0)
    s.say("角落 ", "$\\chi=-6$", "：偶數次方和中心一樣大", wait=3.0)
    s.play(FadeIn(c3), FadeIn(c4), run_time=1.0)
    s.say("角落也有一座高斯，帶相位 ", "$(-1)^{N+x+y+z}$", wait=3.0)
    s.say("同奇偶：兩座相加 ", "$\\times2$", "；不同奇偶：抵消成 0", wait=3.0)
    s.play(FadeIn(err, shift=0.2 * UP), run_time=0.8)
    s.say("書上的 (2.20) 只算了中心：要乘 2", wait=3.0)
    s.say("同第 1 章 W110：一維的另一半在 ", "$q=\\pm\\pi$", wait=3.0)
    s.play(FadeIn(d1), run_time=0.8)
    s.say("更正後：", f"${nb['two']}$", wait=2.5)
    s.play(FadeIn(d2), FadeIn(d3), run_time=0.8)
    s.say("精確／近似 ", f"${nb['r10']}$", "；", "$N=20$", " 時 ", f"${nb['r20']}$", wait=3.0)
    s.clear_stage()


def part_square(s):
    nb = numbers()
    s1 = zh("正方格子（主例）同樣有中心＋角落：", font_size=26)
    s2 = Tex(r"C(N;\vec 0)\approx\frac{2\cdot4^N}{\pi N}", font_size=44).set_color(D.EXACT_COLOR)
    s3 = mixed("$N=10$", "：", f"${nb['sq_asym']}$", "（精確 ", f"${nb['sq_exact']}$", "）", font_size=28)
    s2b = mixed("中心一份 ", "$\\frac{4^N}{\\pi N}$", "，角落再一份 ", "$\\Rightarrow\\times2$", font_size=24)
    s4 = mixed("沒乘 2 的 ", "$\\frac{4^N}{\\pi N}$", "：", f"${nb['sq_half']}$", font_size=26, color=GREY_B)
    rows = [["N"] + [str(N) for N in SQ_NS], [zh("精確／近似", font_size=22)] + [str(r) for r in nb["sq_ratios"]]]
    tab = D.tex_table(rows, [1.9, 1.0, 1.0, 1.0, 1.0], row_h=0.6, font_size=28)
    col = VGroup(s1, s2, s2b, s3, s4, tab).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    place_col(col, -6.5, 2.85)

    x0, y0, W, H = 1.2, -2.0, 5.0, 3.8                      # 比值圖：N ∈ [0, 105]、比值 ∈ [0.94, 1.0]
    px = lambda N: x0 + N / 105 * W
    py = lambda r: y0 + (r - 0.94) / 0.06 * H
    axes = VGroup(Arrow(np.array([x0, y0, 0]), np.array([x0 + W + 0.3, y0, 0]), buff=0, thickness=2).set_color(GREY_B),
                  Arrow(np.array([x0, y0, 0]), np.array([x0, y0 + H + 0.3, 0]), buff=0, thickness=2).set_color(GREY_B))
    ticks = VGroup(*[Tex(str(N), font_size=22).set_color(GREY_B).move_to(np.array([px(N), y0 - 0.25, 0])) for N in SQ_NS],
                   *[Tex(t, font_size=22).set_color(GREY_B).move_to(np.array([x0 - 0.4, py(v), 0]))
                     for v, t in [(0.95, "0.95"), (1.0, "1")]],
                   Tex("N", font_size=26).set_color(GREY_B).move_to(np.array([x0 + W + 0.45, y0 - 0.25, 0])))
    one = DashedLine(np.array([x0, py(1.0), 0]), np.array([x0 + W, py(1.0), 0]), dash_length=0.1).set_stroke(GREY_B, 2)
    dots = VGroup(*[Dot(np.array([px(N), py(r), 0]), radius=0.08).set_fill(D.EXACT_COLOR, 1)
                    for N, r in zip(SQ_NS, nb["sq_ratios"])])
    plot = VGroup(axes, ticks, one, dots)
    assert_in_frame(col, plot)
    assert_apart(col, plot)

    s.say("回到主例：正方格子", wait=0.3)
    s.play(FadeIn(s1), FadeIn(s2), run_time=1.0)
    s.wait(1.5)
    s.play(FadeIn(s2b), run_time=0.8)
    s.say("同樣有中心與角落兩座峰：近似要帶因子 2", wait=3.0)
    s.play(FadeIn(s3), FadeIn(s4), run_time=1.0)
    s.say("$N=10$", "：近似 ", f"${nb['sq_asym']}$", "，精確 ", f"${nb['sq_exact']}$", wait=3.0)
    s.say("沒乘 2 的 ", f"${nb['sq_half']}$", " 只有一半", wait=2.5)
    s.play(FadeIn(tab), FadeIn(axes), FadeIn(ticks), ShowCreation(one), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in dots], lag_ratio=0.3), run_time=1.5)
    s.say("比值：0.9513、0.9753、0.9901、0.995", wait=3.0)
    s.say("$N$", " 越大越接近 1：高斯極限加上因子 2 是對的", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("$N$", " 很大時走法數是高斯：寬度 ", "$\\sqrt N$", wait=3.5)
    s.say("格子的細節只剩奇偶性的因子 2", wait=3.5)


class W206GaussianLimit(ChuScene):
    def construct(self):
        self.show_header("W206　高斯極限與奇偶性的因子 2")
        check_math()
        part_question(self)
        part_expand(self)
        part_estimate(self)
        part_corner(self)
        part_square(self)
        part_close(self)
        self.wait(1.0)
