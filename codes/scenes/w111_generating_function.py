"""
W111 生成函數的奇點與動差（§1.5.1–1.5.2，pp. 16–20，(1.58)–(1.70)，Ex. 1.6–1.7）

問：1/(1−2z) 為什麼在 z = ½ 爆掉？這和走法有什麼關係？
  1. g(z,0) = Σ Γ(N) z^N = 1 + 2z + 4z² + 8z³ + ⋯ = 1/(1−2z)（(1.64)(1.65)，Γ(N) = 2^N 是總走法數）；
     畫 g(z,0)，z ∈ [0, 0.48] 往上衝（灰線是部分和）；奇點 z_c = ½ = 1/每步的選擇數；g ~ (z_c − z)^{−γ}，γ = 1 (1.66)。
  2. Cauchy (1.58)：C_N = (1/2πi)∮ f(z)/z^{N+1} dz；半徑 0.25 的圓周、極點 ½ 在圈外；數值積分 N = 10 → 1024。
     一句話 (1.59)–(1.62)：對 1/(1 − 2z + zq²) 在 z ≈ 1/(2 − q²) 取鞍點，又得到 2^N e^{−Nq²/2}。
  3. 動差：−∂²g/∂q²|₀ = 2z/(1−2z)² = 2z + 8z² + 24z³ + 64z⁴ + ⋯，z^N 的係數 N·2^N = Σ d² C(N,d)（(1.69)(1.70)）
     → 平均 N（練習 1.6）；四階係數 2^N(3N² − 2N)，N = 10：286720/1024 = 280（練習 1.7：大 N 時 ≈ 3(d²平均)²）。
  4. 總結：係數（走法數）←→ 函數的奇點（位置 = 成長率 2^N，形狀 = 指數 γ）。
畫面上的數字都由 utils.walk1d 算出，check_math() 核對。
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
from utils import walk1d as W
from utils import walk_draw as D

TABLE_NS = (1, 2, 3, 4, 10)
G_COLOR = YELLOW


def table_rows():
    d2, d4 = W.d2_gen_coeffs(10), W.d4_gen_coeffs(10)
    return [(str(N), str(d2[N]), str(d2[N] // 2 ** N), str(d4[N]), str(d4[N] // 2 ** N)) for N in TABLE_NS]


def numbers():
    c = W.cauchy_coeff(lambda z: 1 / (1 - 2 * z), 10, 0.25)
    return {"cauchy": f"{c.real:.3f}", "cauchy_im": f"{abs(c.imag):.3f}"}


def check_math():
    assert W.gz0_coeffs(3) == [1, 2, 4, 8]
    assert W.gz0_coeffs(10) == [sum(W.distribution(N).values()) for N in range(11)]       # Γ(N) = 2^N
    assert abs(sum(c * 0.25 ** k for k, c in enumerate(W.gz0_coeffs(80))) - 1 / (1 - 2 * 0.25)) < 1e-12
    assert abs(1 / (1 - 2 * 0.48) - 25) < 1e-9                                        # 圖的上緣
    assert numbers() == {"cauchy": "1024.000", "cauchy_im": "0.000"}
    assert W.d2_gen_coeffs(4) == [0, 2, 8, 24, 64]
    for N in range(11):
        dist = W.distribution(N)
        assert W.d2_gen_coeffs(10)[N] == sum(d * d * c for d, c in dist.items())
        assert W.d4_gen_coeffs(10)[N] == sum(d ** 4 * c for d, c in dist.items())
        assert W.d2_gen_coeffs(10)[N] % 2 ** N == 0 and W.d4_gen_coeffs(10)[N] % 2 ** N == 0
    assert table_rows() == [("1", "2", "1", "2", "1"), ("2", "8", "2", "32", "8"), ("3", "24", "3", "168", "21"),
                            ("4", "64", "4", "640", "40"), ("10", "10240", "10", "286720", "280")]
    assert 3 * 10 ** 2 == 300 and W.moment_exact(10, 4) == 280 and W.moment_exact(10, 2) == 10
    c = W.cauchy_coeff(lambda z: 2 * z / (1 - 2 * z) ** 2, 10, 0.25)                # 2z/(1−2z)² 的 z^10 係數
    assert abs(c - 10240) < 1e-6
    q = 0.2                                    # (1.59)–(1.62)：1/(1−2z+zq²) 的 z^N 係數 (2−q²)^N ≈ 2^N e^{−Nq²/2}
    c = W.cauchy_coeff(lambda z: 1 / (1 - 2 * z + z * q * q), 10, 0.25)
    assert abs(c.real - (2 - q * q) ** 10) < 1e-6
    assert abs((2 - q * q) ** 10 / (1024 * math.exp(-10 * q * q / 2)) - 1) < 0.01
    assert abs(W.walk_saddle_z(10, q)[1] - 1 / (2 - q * q)) < 1e-12


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


def partial_sum(K):
    return lambda z: sum((2 * z) ** n for n in range(K + 1))


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：", "$\\tfrac{1}{1-2z}$", " 為什麼在 ", "$z=\\tfrac12$", " 爆掉？和走法有什麼關係？", wait=3.0)


def part_singularity(s):
    s1 = fit(eq(r"g(z,0)=\sum_{N\ge0}\Gamma(N)z^N=1+2z+4z^2+8z^3+\cdots=\frac{1}{1-2z}", "(1.64)(1.65)"), 12.6)
    s2 = mixed("$\\Gamma(N)=2^N$", "：", "$N$", " 步的總走法數（每步 2 種選擇）", font_size=24, color=GREY_B)
    top = VGroup(s1, s2).arrange(DOWN, buff=0.2).move_to(np.array([0, 2.35, 0]))
    ax = D.plot_axes((0, 0.55, 0.1), (0, 26, 5), width=5.2, height=3.2, x_ticks=[0.1, 0.2, 0.3, 0.4, 0.5],
                     y_ticks=[5, 10, 15, 20, 25], font_size=20, x_label="z")
    ax.move_to(np.array([-3.5, -0.85, 0]))
    graph = ax.get_graph(lambda z: 1 / (1 - 2 * z), x_range=[0, 0.48, 0.004]).set_stroke(G_COLOR, 4)
    sums = VGroup(*[ax.get_graph(partial_sum(K), x_range=[0, 0.5, 0.01]).set_stroke(GREY_B, 2, 0.8)
                    for K in (1, 3, 7)])
    asym = DashedLine(ax.c2p(0.5, 0), ax.c2p(0.5, 26)).set_stroke(RED, 2)
    zc = Tex(r"z_c=\tfrac12", font_size=28).set_color(RED).next_to(ax.c2p(0.5, 24), RIGHT, buff=0.12)
    glab = Tex(r"\frac{1}{1-2z}", font_size=28).set_color(G_COLOR).move_to(ax.c2p(0.36, 15))
    r1 = mixed("奇點 ", "$z_c=\\tfrac12$", "：每步有 2 種選擇", font_size=26)
    r2 = mixed("係數 ", "$2^N=z_c^{-N}$", "：位置決定成長率", font_size=26)
    r3 = eq(r"g(z,0)\sim(z_c-z)^{-\gamma},\ \ \gamma=1", "(1.66)", font_size=32)
    r4 = mixed("$\\Rightarrow$", " 臨界現象的類比（第 8、9 章）", font_size=24, color=GREY_B)
    col = place_col(VGroup(r1, r2, r3, r4).arrange(DOWN, aligned_edge=LEFT, buff=0.35), 0.4, 1.05)
    assert_in_frame(top, ax, zc, glab, col)
    assert_apart(top, ax)
    assert_apart(VGroup(ax, zc), col)

    s.play(Write(s1[0]), FadeIn(s1[1]), run_time=2.0)
    s.play(FadeIn(s2), run_time=0.6)
    s.say("不看位置、只數總數：", "$q=0$", " 的生成函數", wait=3.0)
    s.play(FadeIn(ax), run_time=0.8)
    s.play(LaggedStartMap(ShowCreation, sums, lag_ratio=0.4), run_time=1.5)
    s.say("灰線是部分和：項數越多越接近", wait=1.0)
    s.play(ShowCreation(graph), FadeIn(glab), run_time=1.5)
    s.wait(1.5)
    s.play(ShowCreation(asym), FadeIn(zc), run_time=0.8)
    s.say("越靠近 ", "$z=\\tfrac12$", " 衝得越高：級數在那裡發散", wait=3.0)
    s.play(FadeIn(r1), run_time=0.6)
    s.say("奇點 ", "$z_c=\\tfrac12$", "＝1 ÷ 每步的選擇數", wait=3.0)
    s.play(FadeIn(r2), run_time=0.6)
    s.say("係數長成 ", "$z_c^{-N}=2^N$", "：奇點離原點越近，走法長得越快", wait=3.5)
    s.play(FadeIn(r3), run_time=0.8)
    s.say("一般的奇點 ", "$(z_c-z)^{-\\gamma}$", "：這裡指數 ", "$\\gamma=1$", wait=3.0)
    s.play(FadeIn(r4), run_time=0.6)
    s.say("同樣的形式就是第 8、9 章臨界現象的類比", wait=3.0)
    s.clear_stage()


def part_cauchy(s):
    nb = numbers()
    cp = D.plot_axes((-0.75, 0.75, 0.25), (-0.75, 0.75, 0.25), width=3.8, height=3.8,
                     x_label=r"\mathrm{Re}\,z", y_label=r"\mathrm{Im}\,z")
    cp.move_to(np.array([-4.3, 0.15, 0]))
    cc = D.contour_circle(cp, 0.25)
    pm = D.pole_marker(cp, 0.5 + 0j, r"z=\tfrac12")
    o = Dot(cp.c2p(0, 0), radius=0.06).set_fill(WHITE, 1)
    walker = Dot(cp.c2p(0.25, 0), radius=0.08).set_fill(D.WALK_COLOR, 1)
    c1 = eq(r"C_N=\frac{1}{2\pi i}\oint\frac{f(z)}{z^{N+1}}\,dz", "(1.58)", font_size=36)
    c2 = mixed("$f(z)=\\frac{1}{1-2z}$", "：圓周半徑 ", "$0.25$", "，極點 ", "$\\tfrac12$", " 在圈外", font_size=26)
    c3 = mixed("數值積分（256 點）：", f"$N=10\\ \\to\\ {nb['cauchy']}=2^{{10}}$", font_size=26)
    c4 = mixed("一句話 (1.59)–(1.62)：對 ", "$\\frac{1}{1-2z+zq^2}$", font_size=26)
    c5 = mixed("在 ", "$z\\approx\\frac{1}{2-q^2}$", " 取鞍點 ", "$\\Rightarrow$", " 又得到 ", "$2^Ne^{-Nq^2/2}$",
               font_size=26)
    col = VGroup(c1, c2, c3, c4, c5).arrange(DOWN, aligned_edge=LEFT, buff=0.4)
    for m in col:
        fit(m, 7.4)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.4), -1.3, 2.8)
    hl = SurroundingRectangle(c3, buff=0.1).set_stroke(D.EXACT_COLOR, 2)
    assert_in_frame(cp, pm, col, hl)
    assert_apart(cp, col)
    s.say("反過來：知道函數，怎麼把係數拿回來？", wait=1.0)
    s.play(Write(c1[0]), FadeIn(c1[1]), run_time=1.5)
    s.wait(2.5)
    s.play(FadeIn(cp), FadeIn(o), run_time=0.8)
    s.play(ShowCreation(cc[0]), FadeIn(cc[1]), FadeIn(pm), FadeIn(c2), run_time=1.2)
    s.say("Cauchy 公式：沿著繞原點的小圓走一圈", wait=1.0)
    s.add(walker)
    s.play(MoveAlongPath(walker, cc[0]), run_time=3.0, rate_func=linear)
    s.remove(walker)
    s.say("圓要小：極點 ", "$z=\\tfrac12$", " 必須在圈外", wait=3.5)
    s.play(FadeIn(c3), ShowCreation(hl), run_time=0.8)
    s.say("數值積分 ", "$N=10$", "：得到 1024，正是 ", "$2^{10}$", wait=3.5)
    s.play(FadeIn(c4), FadeIn(c5), run_time=1.0)
    s.say("加上 ", "$q$", "：在 ", "$z\\approx\\tfrac{1}{2-q^2}$", " 取鞍點，又回到 W110 的高斯", wait=5.0)
    s.clear_stage()


def part_moments(s):
    m1 = fit(Tex(r"-\frac{\partial^2g}{\partial q^2}\Big|_{q=0}=\frac{2z}{(1-2z)^2}=2z+8z^2+24z^3+64z^4+\cdots",
                 font_size=30), 12.6)
    m2 = VGroup(mixed("$z^N$", " 的係數 ", "$=N\\,2^N=\\textstyle\\sum_d d^2\\,C(N,d)$", font_size=24),
                Tex("(1.69)(1.70)", font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)
    top = VGroup(m1, m2).arrange(DOWN, buff=0.25)
    top.move_to(np.array([0, 3.1 - top.get_height() / 2, 0]))
    head = ["N", r"\textstyle\sum_d d^2C(N,d)", r"\overline{d^2}", r"\textstyle\sum_d d^4C(N,d)", r"\overline{d^4}"]
    tab = D.tex_table([head] + [list(r) for r in table_rows()], [0.7, 2.0, 0.8, 2.0, 0.9], row_h=0.46,
                      font_size=24)
    y0 = top.get_bottom()[1] - 0.35
    place_col(tab, -6.5, y0)
    k1 = mixed("練習 1.6：", "$\\overline{d^2}=N\\,2^N/2^N=N$", font_size=24)
    k2 = mixed("$N=10$", "：", "$10240/1024=10$", font_size=24)
    k3 = mixed("四階：", "$\\partial_q^4g\\big|_{q=0}$", " 的係數 ", "$2^N(3N^2-2N)$", font_size=24)
    k4 = mixed("$N=10$", "：", "$286720/1024=280$", font_size=24)
    k5 = mixed("練習 1.7：", "$\\overline{d^4}\\approx3\\big(\\overline{d^2}\\big)^2$", "（", "$280$", " 對 ",
               "$3\\times10^2=300$", "）", font_size=24)
    col = VGroup(k1, k2, k3, k4, k5).arrange(DOWN, aligned_edge=LEFT, buff=0.24)
    for m in col:
        fit(m, 6.2)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.24), 0.45, y0)
    assert_in_frame(top, tab, col)
    assert_apart(top, tab)
    assert_apart(tab, col)
    row = lambda i: VGroup(*[D.table_cell(tab, i, j, 5) for j in range(5)])
    s.say("生成函數也能算動差：對 ", "$q$", " 微分兩次", wait=0.3)
    s.play(Write(m1), run_time=2.0)
    s.wait(2.5)
    s.play(FadeIn(m2), run_time=0.8)
    s.say("$z^N$", " 的係數 ", "$N\\,2^N$", " 正是所有走法的 ", "$d^2$", " 加總", wait=4.0)
    s.play(ShowCreation(tab[1]), FadeIn(row(0)), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(row(i)) for i in range(1, 6)], lag_ratio=0.4), run_time=2.0)
    s.say("一列一列核對：除以 ", "$2^N$", " 就是平均", wait=3.0)
    s.play(FadeIn(k1), FadeIn(k2), run_time=0.8)
    s.say("練習 1.6：", "$\\overline{d^2}=N$", "；", "$N=10$", " 得 10", wait=3.0)
    s.play(FadeIn(k3), FadeIn(k4), run_time=0.8)
    s.say("四階一樣做：", "$N=10$", " 得 280（和 W105 一致）", wait=3.0)
    s.play(FadeIn(k5), run_time=0.8)
    s.say("練習 1.7：大 ", "$N$", " 時 ", "$\\overline{d^4}\\approx3(\\overline{d^2})^2$", "（高斯的值）", wait=3.5)
    s.clear_stage()


def part_summary(s):
    left = VGroup(zh("係數：走法數", font_size=30, color=YELLOW),
                  mixed("成長率 ", "$2^N$", font_size=28),
                  mixed("前因子（下一幕）", font_size=28)).arrange(DOWN, buff=0.55)
    right = VGroup(zh("函數的奇點", font_size=30, color=RED),
                   mixed("位置 ", "$z_c=\\tfrac12$", font_size=28),
                   mixed("形狀 ", "$(z_c-z)^{-\\gamma}$", font_size=28)).arrange(DOWN, buff=0.55)
    left.move_to(np.array([-3.3, 0.6, 0]))
    right.move_to(np.array([3.3, 0.6, 0]))
    for k in (1, 2):
        right[k].match_y(left[k])
    boxes = VGroup(*[SurroundingRectangle(g, buff=0.35).set_stroke(c, 2) for g, c in ((left, YELLOW), (right, RED))])
    links = VGroup(*[Tex(r"\longleftrightarrow", font_size=48).move_to(np.array([0, left[k].get_center()[1], 0]))
                     for k in (1, 2)])
    note = mixed("$2^N=z_c^{-N}$", "：奇點的位置決定成長率；形狀決定前因子", font_size=26, color=GREY_B)
    note.move_to(np.array([0, -1.9, 0]))
    assert_in_frame(boxes, links, note)
    assert_apart(boxes, note)
    s.say("一張圖總結", wait=1.0)
    s.play(FadeIn(left[0]), FadeIn(right[0]), ShowCreation(boxes), run_time=1.0)
    s.play(FadeIn(left[1]), FadeIn(right[1]), GrowFromCenter(links[0]), run_time=1.0)
    s.say("奇點的位置 ", "$\\longleftrightarrow$", " 走法的成長率", wait=3.0)
    s.play(FadeIn(left[2]), FadeIn(right[2]), GrowFromCenter(links[1]), run_time=1.0)
    s.say("奇點的形狀（指數 ", "$\\gamma$", "）", "$\\longleftrightarrow$", " 前因子：下一幕", wait=3.5)
    s.play(FadeIn(note), run_time=0.6)
    s.wait(3.5)


def part_close(s):
    s.say("生成函數的奇點，告訴你走法怎麼長大", wait=3.5)


class W111GeneratingFunction(ChuScene):
    def construct(self):
        self.show_header("W111　生成函數的奇點與動差")
        check_math()
        part_question(self)
        part_singularity(self)
        part_cauchy(self)
        part_moments(self)
        part_summary(self)
        part_close(self)
        self.wait(2.0)
