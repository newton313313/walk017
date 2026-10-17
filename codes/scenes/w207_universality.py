"""
W207 普適性：三種晶格、同一條高斯（§2.1.4，pp. 35–36，(2.21)–(2.22)）

問：簡單立方、體心立方、面心立方——走很多步之後，回到原點的機率會不會不一樣？
  1. 三種晶格（L.lattice3d；鄰居 6、8、12 個，最近鄰距離 a = 1）：χ ≈ Z − c2 q²：6 − q²、8 − (4/3)q²、12 − 2q²
     → χ/Z ≈ 1 − q²/6，三個一樣（GF.CHI_SMALL_Q）。
  2. P(N; 0) = C(N; 0)/Z^N 的表（GF.UNIV_TABLE）看起來不一樣；除以單胞體積 Ω 與奇偶因子 g、乘 N^{3/2}：
     N = 10 → 0.307、0.306、0.306；N = 40 → 0.3238；點列圖收斂到 (3/2π)^{3/2} = 0.3299。
  3. (2.21)–(2.22)：q⁴ 項的修正 ∝ 1/N；一維擲幣的實例（GF.coin_quartic）；筆誤框 e^{−3r²/2aN} → e^{−3r²/2a²N}。
  4. 普適性的兩個來源；反例 χ = c0 + c'|q|^{1.3}（下一幕）。
fcc 的 N = 40 用本檔的 fcc_returns（三項式展開取常數項）；GF.returns_fcc 走 count_grid，N = 40 太慢，N ≤ 10 互相核對。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_gf as GF
from utils import lattice_draw as L
from utils import walk_draw as D

KINDS = ("sc", "bcc", "fcc")
COLORS = {"sc": TEAL, "bcc": "#F2A33A", "fcc": "#C792EA"}
NS = (2, 4, 6, 8, 10)
LIM = (3 / (2 * math.pi)) ** 1.5
Y0 = 0.2                                          # 點列圖的縱軸從 0.2 開始（畫 value − 0.2）


def fcc_returns(N):
    """fcc 回到原點的走法數：χ = XY + YZ + ZX（X = x + 1/x，…）的 N 次方取常數項；
    三項式展開 (XY)^i (YZ)^j (ZX)^k = X^{N−j} Y^{i+j} Z^{N−i}，每個變數各自取常數項 C(n, n/2) = GF.count1(n, 0)。"""
    f = math.factorial
    c = lambda n: GF.count1(n, 0)
    return sum(f(N) // (f(i) * f(j) * f(N - i - j)) * c(N - j) * c(i + j) * c(N - i)
               for i in range(N + 1) for j in range(N + 1 - i))


def returns(kind, N):
    return {"sc": GF.returns_sc, "bcc": GF.returns_bcc, "fcc": fcc_returns}[kind](N)


def p_return(kind, N):
    return returns(kind, N) / GF.CHI_SMALL_Q[kind][0] ** N


def parity_g(kind):
    return 2 if GF.BIPARTITE[kind] else 1


def scaled(kind, N):
    return N ** 1.5 * p_return(kind, N) / (parity_g(kind) * GF.CELL_VOLUME[kind])


def numbers():
    nb = {k: {N: str(round(p_return(k, N), 5)) for N in NS} for k in KINDS}
    nb["omega"] = {k: f"{GF.CELL_VOLUME[k]:.4f}".rstrip("0").rstrip(".") for k in KINDS}
    nb["s10"] = {k: f"{scaled(k, 10):.3f}" for k in KINDS}
    nb["s40"] = {k: f"{scaled(k, 40):.4f}" for k in KINDS}
    nb["lim"] = f"{LIM:.4f}"
    nb["coin"] = {(N, d): [f"{v:.6f}" for v in GF.coin_quartic(N, d)] for N, d in ((10, 0), (10, 2), (20, 0))}
    nb["book"] = f"{math.exp(-3 * 4 / (2 * math.sqrt(3) * 10)):.4f}"
    nb["right"] = f"{math.exp(-3 * 4 / (2 * 3 * 10)):.4f}"
    return nb


def check_math():
    nb = numbers()
    for k, (Z, c2) in (("sc", (6, 1)), ("bcc", (8, Fraction(4, 3))), ("fcc", (12, 2))):
        assert GF.CHI_SMALL_Q[k] == (Z, c2) and Fraction(c2) / Z == Fraction(1, 6)
    q = (0.011, -0.007, 0.013)
    q2 = sum(t * t for t in q)
    for k, fn in (("sc", GF.chi_sc), ("bcc", GF.chi_bcc), ("fcc", GF.chi_fcc)):
        Z, c2 = GF.CHI_SMALL_Q[k]
        assert abs(fn(0, 0, 0) - Z) < 1e-12 and abs((Z - fn(*q)) / q2 - float(c2)) < 1e-3
    assert [fcc_returns(N) for N in range(11)] == [GF.returns_fcc(N) for N in range(11)]
    for k in KINDS:
        assert {N: float(nb[k][N]) for N in NS} == GF.UNIV_TABLE[k]
        assert float(nb["s10"][k]) == GF.UNIV_SCALED10[k] and nb["s40"][k] == "0.3238"
    assert nb["omega"] == {"sc": "1", "bcc": "0.7698", "fcc": "0.7071"} and nb["lim"] == "0.3299"
    assert [parity_g(k) for k in KINDS] == [2, 2, 1] and returns("fcc", 3) == 48 and returns("sc", 3) == 0
    assert nb["coin"] == {(10, 0): ["0.246094", "0.252313", "0.246005"], (10, 2): ["0.205078", "0.206577", "0.205268"],
                          (20, 0): ["0.176197", "0.178412", "0.176182"]}
    assert (nb["book"], nb["right"]) == ("0.7072", "0.8187")
    assert abs(27 * (1 / 12) / (10 * 3 ** 2) - 0.025) < 1e-15            # 擲幣：27b/(Na⁴) = 1/(4N)
    vals = [scaled(k, N) - Y0 for k in KINDS for N in range(2, 41, 2)]
    assert 0 < min(vals) and max(vals) < 0.15 and LIM - Y0 < 0.15


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：簡單立方、體心立方、面心立方三種晶格", wait=2.5)
    s.say("走很多步之後，回到原點的機率會不會不一樣？", wait=3.0)


def part_lattices(s):
    xs = (-4.5, 0.0, 4.5)
    titles = [zh("簡單立方 sc：6 個鄰居（面心）", font_size=24), zh("體心立方 bcc：8 個鄰居（角）", font_size=24),
              zh("面心立方 fcc：12 個鄰居（稜中點）", font_size=24)]
    chis = [tx(r"\chi=2(\cos q_x+\cos q_y+\cos q_z)", 24),
            tx(r"\chi=8\cos\tfrac{q_x}{\sqrt3}\cos\tfrac{q_y}{\sqrt3}\cos\tfrac{q_z}{\sqrt3}", 24),
            tx(r"\chi=4\Big[\cos\tfrac{q_x}{\sqrt2}\cos\tfrac{q_y}{\sqrt2}+\cdots\Big]", 24)]
    smalls = [tx(r"\chi\approx 6-q^2", 32), tx(r"\chi\approx 8-\tfrac43\,q^2", 32), tx(r"\chi\approx 12-2q^2", 32)]
    ratios = [tx(r"\frac{\chi}{Z}\approx 1-\frac{q^2}{6}", 32, YELLOW) for _ in KINDS]
    lats = VGroup()
    for k, x, t, c, sm, r in zip(KINDS, xs, titles, chis, smalls, ratios):
        lats.add(L.lattice3d(k, unit=0.8, origin=np.array([x, 1.3, 0]), color=COLORS[k]))
        fit(t, 4.2).move_to(np.array([x, 2.8, 0]))
        fit(c, 4.2).move_to(np.array([x, -0.2, 0]))
        sm.move_to(np.array([x, -1.0, 0]))
        r.move_to(np.array([x, -2.1, 0]))
    box = SurroundingRectangle(VGroup(*ratios), buff=0.15).set_stroke(YELLOW, 2)
    assert_in_frame(*titles, lats, *chis, *smalls, box)
    for a, b in zip(chis, chis[1:]):
        assert_apart(a, b)
    for a, b in zip(titles, titles[1:]):
        assert_apart(a, b)
    s.play(*[FadeIn(t) for t in titles], *[ShowCreation(m[0]) for m in lats], run_time=1.2)
    s.play(*[FadeIn(m[3]) for m in lats], run_time=0.5)
    s.say("三種立方晶格：最近鄰距離都取 ", "$a=1$", wait=0.3)
    for m in lats:
        s.play(ShowCreation(m[1]), FadeIn(m[2]), run_time=0.9)
    s.say("鄰居分別在立方體的面心、角、稜中點：6、8、12 個", wait=3.0)
    s.play(*[FadeIn(c) for c in chis], run_time=1.0)
    s.say("各自的結構函數 ", "$\\chi(\\vec q)$", "：在 W204 見過", wait=3.0)
    s.play(*[FadeIn(sm, shift=0.2 * DOWN) for sm in smalls], run_time=1.0)
    s.say("在 ", "$\\vec q=\\vec 0$", " 附近展開：", "$\\chi\\approx Z-c_2q^2$", wait=3.5)
    s.say("$Z$", "＝鄰居數；", "$q^2$", " 的係數各不相同", wait=3.0)
    s.play(*[FadeIn(r) for r in ratios], run_time=1.0)
    s.play(ShowCreation(box), run_time=0.8)
    s.say("但除以鄰居數 ", "$Z$", "：三個一模一樣！", wait=3.5)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    head = ["N"] + [str(N) for N in NS] + [r"\Omega", "g"]
    rows = [head] + [[tx(r"\mathrm{" + k + "}", 28, COLORS[k])] + [nb[k][N] for N in NS]
                     + [nb["omega"][k], str(parity_g(k))] for k in KINDS]
    tab = D.tex_table(rows, [1.2] + [1.4] * 5 + [1.3, 0.8], row_h=0.55, font_size=28)
    tab.move_to(np.array([0, 2.35 - tab.get_height() / 2, 0]))
    late = VGroup(*[D.table_cell(tab, i, j, 8) for i in range(4) for j in (6, 7)])
    early = VGroup(*[c for c in tab[0] if c not in late.submobjects])
    title = mixed("回到原點的機率 ", "$P(N;\\vec 0)=C(N;\\vec 0)/Z^N$", font_size=28).move_to(np.array([0, 2.85, 0]))
    form = mixed("$\\dfrac{N^{3/2}\\,P(N;\\vec 0)}{g\\,\\Omega}$", "：", "$N=10$", " 時 ",
                 "$" + nb["s10"]["sc"] + "$", "、", "$" + nb["s10"]["bcc"] + "$", "、",
                 "$" + nb["s10"]["fcc"] + "$", font_size=30)
    form.move_to(np.array([0, -0.55, 0]))
    form40 = mixed("$N=40$", "：三個都是 ", "$" + nb["s40"]["sc"] + "$", "；極限 ",
                   "$(3/2\\pi)^{3/2}=" + nb["lim"] + "$", font_size=30).move_to(np.array([0, -1.75, 0]))
    assert_in_frame(title, tab, form, form40)
    assert_apart(title, tab)
    assert_apart(tab, form)
    assert_apart(form, form40)
    s.play(FadeIn(title), run_time=0.6)
    s.play(FadeIn(early), ShowCreation(tab[1]), run_time=1.2)
    s.say("$N=2,4,\\dots,10$", " 步回到原點的機率（", "$Z^N$", "＝所有走法）", wait=3.5)
    s.say("看起來不一樣：鄰居越多，越難剛好回到家", wait=3.0)
    s.play(LaggedStart(*[FadeIn(c) for c in late], lag_ratio=0.1), run_time=1.2)
    s.say("$\\Omega$", "＝單胞體積：每個格點佔多少空間", wait=3.0)
    s.say("$g=2$", "：sc、bcc 只有偶數步回得來（W206 的奇偶性）", wait=3.5)
    s.say("fcc 有三角形迴路，奇數步也回得來：", "$g=1$", wait=3.0)
    s.play(FadeIn(form, shift=0.2 * UP), run_time=1.0)
    s.say("除掉 ", "$g\\,\\Omega$", "、乘上 ", "$N^{3/2}$", "：", "$N=10$", " 已經幾乎一樣", wait=3.5)
    s.play(FadeIn(form40, shift=0.2 * UP), run_time=1.0)
    s.say("$N=40$", "：到小數第四位都一樣", wait=3.0)
    s.clear_stage()


def part_plot(s):
    nb = numbers()
    ax = D.plot_axes((0, 42, 10), (0, 0.15, 0.05), width=9.0, height=4.6, x_ticks=[10, 20, 30, 40],
                     y_ticks=[0.05, 0.10], y_tick_labels=["0.25", "0.30"], font_size=22, x_label="N")
    ax.move_to(np.array([-1.3, -0.35, 0]))
    ylab = mixed("$N^{3/2}P(N;\\vec 0)/(g\\,\\Omega)$", font_size=26).next_to(ax, UP, buff=0.1)
    ylab.align_to(ax, LEFT)
    line = DashedLine(ax.c2p(0, LIM - Y0), ax.c2p(42, LIM - Y0)).set_stroke(GREY_A, 2)
    llab = tx(r"(3/2\pi)^{3/2}=" + nb["lim"], 26, GREY_A).next_to(ax.c2p(31, LIM - Y0), UP, buff=0.12)
    dots = {}
    for k, dx in zip(KINDS, (-0.35, 0.0, 0.35)):
        dots[k] = VGroup(*[Dot(ax.c2p(N + dx, scaled(k, N) - Y0), radius=0.06).set_fill(COLORS[k], 1)
                           for N in range(2, 41, 2)])
    legend = VGroup(*[VGroup(Dot(radius=0.08).set_fill(COLORS[k], 1), tx(r"\mathrm{" + k + "}", 28, COLORS[k]))
                      .arrange(RIGHT, buff=0.15) for k in KINDS]).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    legend.move_to(np.array([4.9, 0.6, 0]))
    end = tx(r"N=40:\ " + nb["s40"]["sc"], 26).next_to(legend, DOWN, buff=0.5).align_to(legend, LEFT)
    assert_in_frame(ax, ylab, llab, legend, end)
    assert_apart(ax, legend)
    assert_apart(ylab, llab)
    s.play(FadeIn(ax), FadeIn(ylab), run_time=1.0)
    s.play(ShowCreation(line), FadeIn(llab), run_time=0.8)
    s.say("縱軸：調整過的回家機率；虛線是高斯給的極限", wait=1.0)
    s.play(FadeIn(legend), *[LaggedStart(*[FadeIn(d, scale=0.5) for d in dots[k]], lag_ratio=0.15) for k in KINDS],
           run_time=3.0)
    s.say("三種顏色的點：一開始分得開，很快就疊在一起", wait=3.0)
    s.play(FadeIn(end), run_time=0.6)
    s.say("一起往同一條虛線爬：", "$N\\to\\infty$", " 時都是 ", "$(3/2\\pi)^{3/2}$", wait=3.5)
    s.clear_stage()


def part_quartic(s):
    nb = numbers()
    e1 = tx(r"\int e^{-iqr-\frac{a^2N}{6}q^2-bNq^4}\,dq\approx\int e^{-iqr-\frac{a^2N}{6}q^2}\big(1-bNq^4\big)\,dq", 30)
    e1t = VGroup(e1, tx(r"\text{(2.21)}", 22, GREY_B)).arrange(RIGHT, buff=0.3)
    e2a = tx(r"=\sqrt{\frac{6\pi}{Na^2}}\,e^{-3r^2/2a^2N}", 30)
    e2b = tx(r"\Big[1-\frac{27b}{Na^4}\Big(1-6\frac{r^2}{Na^2}+3\Big(\frac{r^2}{Na^2}\Big)^2\Big)\Big]", 30)
    e2t = VGroup(VGroup(e2a, e2b).arrange(RIGHT, buff=0.1), tx(r"\text{(2.22)}", 22, GREY_B)).arrange(RIGHT, buff=0.3)
    fit(e1t, 12.5).move_to(np.array([0, 2.75, 0]))
    fit(e2t, 12.5).move_to(np.array([0, 1.6, 0]))
    corr = SurroundingRectangle(e2b, buff=0.08).set_stroke(YELLOW, 2)
    coin = mixed("一維擲幣：", "$\\cos^Nq=e^{-Nq^2/2-Nq^4/12-\\cdots}$", "，即 ", "$a^2=3,\\ b=\\tfrac1{12}$",
                 font_size=26)
    hdr = [tx("(N,d)", 26, YELLOW), zh("精確", font_size=24, color=YELLOW), zh("高斯", font_size=24, color=YELLOW),
           mixed("加 ", "$q^4$", " 修正", font_size=24, color=YELLOW)]
    rows = [hdr] + [[f"({N},{d})"] + nb["coin"][(N, d)] for N, d in ((10, 0), (10, 2), (20, 0))]
    tab = D.tex_table(rows, [1.5, 1.7, 1.7, 1.9], row_h=0.5, font_size=26)
    left = VGroup(coin, tab).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    left.move_to(np.array([-6.6 + left.get_width() / 2, -1.25, 0]))
    box = D.erratum_box("p. 35，(2.22)", r"e^{-3r^2/2aN}", r"e^{-3r^2/2a^2N}", width=4.4)
    note = mixed("$N=10,\\ r=2,\\ a=\\sqrt3$", "：", font_size=22)
    note2 = mixed("照印 ", "$" + nb["book"] + "$", "，正確 ", "$" + nb["right"] + "$", font_size=22)
    err = VGroup(box, note, note2).arrange(DOWN, buff=0.15)
    err.move_to(np.array([6.7 - err.get_width() / 2, -1.25, 0]))
    assert_in_frame(e1t, e2t, left, err)
    assert_apart(e2t, left)
    assert_apart(left, err)
    assert_apart(e2t, err)
    s.say("為什麼只看 ", "$q^2$", "？把下一項 ", "$q^4$", " 也留下來看看", wait=0.3)
    s.play(Write(e1t), run_time=1.8)
    s.say("$e^{-bNq^4}$", " 展開到一次，再做高斯積分", wait=2.5)
    s.play(FadeIn(e2t), run_time=1.2)
    s.play(ShowCreation(corr), run_time=0.8)
    s.say("修正項 ", "$27b/Na^4$", "：正比於 ", "$1/N$", "，步數越多越小", wait=3.5)
    s.play(FadeIn(coin), run_time=0.8)
    s.say("實例：擲硬幣的 ", "$\\cos^Nq$", "，", "$a^2=3$", "、", "$b=\\tfrac1{12}$", wait=3.0)
    s.play(FadeIn(tab), run_time=1.0)
    s.say("加上修正：0.252 變成 0.246，幾乎就是精確值", wait=3.5)
    s.say("$N=20$", "：高斯本身就更準了（修正 ", "$\\propto 1/N$", "）", wait=3.0)
    s.play(FadeIn(err, shift=0.2 * LEFT), run_time=0.8)
    s.say("書上 (2.22) 的指數少了一個 ", "$a$", "：應是 ", "$2a^2N$", wait=3.5)
    s.clear_stage()


def part_sources(s):
    t = zh("普適性的兩個來源", font_size=32, color=YELLOW)
    r1 = mixed("① ", "$\\chi(\\vec q)=c_0+c_2q^2+\\cdots$", "：最低的修正是 ", "$q^2$", font_size=30)
    r2 = mixed("② ", "$N$", " 很大時只有最低兩項重要：", "$q^4$", " 的修正 ", "$\\propto 1/N$", font_size=30)
    r3 = mixed("反例：", "$\\chi=c_0+c'|q|^{1.3}+\\cdots$", "$\\Rightarrow$", "不是高斯（下一幕）", font_size=30,
               color=D.WRONG_COLOR)
    col = VGroup(t, r1, r2, r3).arrange(DOWN, aligned_edge=LEFT, buff=0.55)
    col.move_to(np.array([0, 0.4, 0]))
    assert_in_frame(col)
    s.play(FadeIn(t), run_time=0.6)
    s.play(FadeIn(r1, shift=0.2 * RIGHT), run_time=0.8)
    s.say("晶格的細節只藏在 ", "$c_0$", "、", "$c_2$", " 兩個數裡", wait=3.0)
    s.play(FadeIn(r2, shift=0.2 * RIGHT), run_time=0.8)
    s.say("高次項的貢獻隨 ", "$N$", " 變小：只剩一條高斯", wait=3.0)
    s.play(FadeIn(r3, shift=0.2 * RIGHT), run_time=0.8)
    s.say("如果 ", "$\\chi$", " 一開始就不是 ", "$q^2$", "：普適性就失效", wait=3.0)


def part_close(s):
    s.say("只要 ", "$\\chi(\\vec q)$", " 從 ", "$q^2$", " 開始變，走很多步後一律是同一條高斯", wait=3.5)
    s.say("這就是普適性", wait=2.0)


class W207Universality(ChuScene):
    def construct(self):
        self.show_header("W207　普適性：三種晶格、同一條高斯")
        check_math()
        part_question(self)
        part_lattices(self)
        part_table(self)
        part_plot(self)
        part_quartic(self)
        part_sources(self)
        part_close(self)
        self.wait(1.0)
