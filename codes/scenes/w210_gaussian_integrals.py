"""
W210 補充 1：高斯積分與它的生成函數（§2.2，pp. 41–42，(2:S1-1)–(2:S1-6)）

問：∫ x⁴ e^{−x²} dx = ? 能不能一次算出所有 ∫ xⁿ e^{−ax²} dx？
  1. (2:S1-1) ∫ e^{−ax²} dx = √(π/a)（a = 1：√π = 1.7725；平方後換極座標 I² = π）。
  2. 生成函數 (2:S1-3) ∫ e^{−ax²+bx} dx：配方 (2:S1-4) ax² − bx = a(x − b/2a)² − b²/4a；筆誤框 −b²/2a → −b²/4a
     （a = 1、b = 2：x² − 2x = (x−1)² − 1，照印是 (x−1)² − 2）。動畫：e^{−x²} 右移 1 再乘 e → (2:S1-5) e^{b²/4a}√(π/a) = e√π = 4.818。
  3. 展開比較係數 (2:S1-6)：∫ xⁿ e^{−ax²} dx = √(π/a) n!/((n/2)!(4a)^{n/2})；a = 1、n = 0, 2, 4, 6：1.7725、0.8862、1.3293、3.3234
     （GF.gauss_moment，數值積分 GF.gauss_moment_numeric 核對）；先猜再揭曉 3√π/4。
  4. 呼應：四階／二階² = 3（W105 的峰度）；(2.22) 的修正 27b/Na⁴ 用的就是 ∫q⁴e^{−αq²}（W207）。
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
from combkit.draw.styles import boxed, mixed, zh
from utils import walk_gf as GF
from utils import walk_draw as D

NS = (0, 2, 4, 6)
FRAC_TEX = {Fraction(1): "1", Fraction(1, 2): r"\tfrac12", Fraction(3, 4): r"\tfrac34", Fraction(15, 8): r"\tfrac{15}{8}"}


def numbers():
    return {"sqrtpi": f"{math.sqrt(math.pi):.4f}", "gen": f"{GF.gauss_gen(1, 2):.3f}",
            "mom": [f"{GF.gauss_moment(n):.4f}" for n in NS], "num": [f"{GF.gauss_moment_numeric(n):.4f}" for n in NS],
            "frac": [FRAC_TEX[GF.gauss_moment_frac(n)] for n in NS]}


def check_math():
    nb = numbers()
    assert nb["sqrtpi"] == "1.7725" and nb["gen"] == "4.818" and abs(GF.gauss_gen(1, 2) - math.e * math.sqrt(math.pi)) < 1e-12
    assert nb["mom"] == ["1.7725", "0.8862", "1.3293", "3.3234"] and nb["num"] == nb["mom"]
    assert [GF.gauss_moment_frac(n) for n in NS] == [1, Fraction(1, 2), Fraction(3, 4), Fraction(15, 8)]
    assert all(abs(GF.gauss_moment_numeric(n)) < 1e-12 for n in (1, 3, 5))        # 奇數 n：0
    for x in (-1.3, 0.0, 0.7, 2.5):                                                 # (2:S1-4) 配方，a = 1、b = 2
        assert abs(GF.square_complete(1, 2, x) - (x * x - 2 * x)) < 1e-12
        assert abs(GF.square_complete(1, 2, x, book=True) - (x * x - 2 * x) + 1) < 1e-12   # 照印的差 1
    for a, b in ((1, 2), (2, 1.5), (0.5, -1)):                                      # (2:S1-5) 與數值積分
        num = GF._simpson(lambda x: math.exp(-a * x * x + b * x), -30, 30, 6000)
        assert abs(num - GF.gauss_gen(a, b)) < 1e-9
    assert GF.gauss_moment_frac(4) / GF.gauss_moment_frac(2) ** 2 == 3              # 峰度 3
    al, N, b = 3 * 10 / 6, 10, 1 / 12                                              # W207：a² = 3、N = 10 → α = 5
    ratio = GF.gauss_moment(4, al) / GF.gauss_moment(0, al)
    assert abs(ratio - 3 / (4 * al * al)) < 1e-15 and abs(b * N * ratio - 27 * b / (N * 9)) < 1e-15
    for N, r in ((10, 0.0), (10, 2.0)):                                            # (2.20)：a = N/6、b = −ir
        num = GF._simpson(lambda q: math.cos(q * r) * math.exp(-N * q * q / 6), -12, 12, 4000)
        assert abs(num - math.sqrt(6 * math.pi / N) * math.exp(-3 * r * r / (2 * N))) < 1e-12


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def tagged(tex, tag, font_size=32, color=WHITE):
    return VGroup(tx(tex, font_size, color), tx(r"\text{" + tag + "}", 22, GREY_B)).arrange(RIGHT, buff=0.3)


def area(ax, f, x0, x1, color, opacity=0.45):
    xs = np.linspace(x0, x1, 120)
    pts = [ax.c2p(x0, 0)] + [ax.c2p(x, f(x)) for x in xs] + [ax.c2p(x1, 0)]
    return Polygon(*pts).set_stroke(width=0).set_fill(color, opacity)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：", "$\\int_{-\\infty}^{\\infty}x^4e^{-x^2}dx=\\ ?$", wait=2.5)
    s.say("能不能一次算出所有的 ", "$\\int x^ne^{-ax^2}dx$", "？", wait=3.0)


def part_basic(s):
    nb = numbers()
    t1 = tagged(r"I=\int_{-\infty}^{\infty}e^{-ax^2}\,dx=\sqrt{\frac{\pi}{a}}", "(2:S1-1)")
    t2 = mixed("平方，換極座標（", "$a=1$", "）：", font_size=26)
    t3 = tx(r"I^2=\iint e^{-(x^2+y^2)}\,dx\,dy=\int_0^{2\pi}\int_0^\infty e^{-r^2}r\,dr\,d\theta=\pi", 30)
    t4 = mixed("$a=1$", "：", "$I=\\sqrt\\pi=" + nb["sqrtpi"] + "$", font_size=30, color=YELLOW)
    for m in (t1, t3):
        fit(m, 7.3)
    col = VGroup(t1, t2, t3, t4).arrange(DOWN, aligned_edge=LEFT, buff=0.4)
    col.move_to(np.array([-6.6 + col.get_width() / 2, 0.3, 0]))
    ax = D.plot_axes((-3, 3, 1), (0, 1.2, 0.5), width=5.2, height=3.0, x_ticks=[-2, -1, 1, 2], y_ticks=[1],
                     font_size=20, x_label="x")
    ax.move_to(np.array([3.95, 0.2, 0]))
    g = ax.get_graph(lambda x: math.exp(-x * x), x_range=[-3, 3, 0.02]).set_stroke(D.BAR_COLOR, 3)
    ar = area(ax, lambda x: math.exp(-x * x), -3, 3, D.BAR_COLOR)
    lab = tx(r"\sqrt\pi\approx" + nb["sqrtpi"], 28, YELLOW).move_to(ax.c2p(0, 0.35))
    gl = tx(r"e^{-x^2}", 28, D.BAR_COLOR).next_to(ax.c2p(1.1, math.exp(-1.21)), UR, buff=0.1)
    assert_in_frame(col, ax, gl)
    assert_apart(col, ax)
    s.play(Write(t1), run_time=1.5)
    s.say("最基本的高斯積分", wait=2.5)
    s.play(FadeIn(ax), ShowCreation(g), FadeIn(gl), run_time=1.2)
    s.play(FadeIn(ar), FadeIn(lab), run_time=0.8)
    s.say("曲線下的面積：", "$\\sqrt{\\pi}\\approx" + nb["sqrtpi"] + "$", wait=2.5)
    s.play(FadeIn(t2), FadeIn(t3), run_time=1.2)
    s.say("一句話的證明：平方變成二維，換極座標", wait=3.0)
    s.say("多出來的 ", "$r$", " 讓積分變得能直接算：", "$I^2=\\pi$", wait=3.0)
    s.play(FadeIn(t4), run_time=0.8)
    s.say("$a\\neq1$", "：把 ", "$x$", " 換成 ", "$x/\\sqrt a$", "，多出一個 ", "$1/\\sqrt a$", wait=3.0)
    s.clear_stage()


def part_generating(s):
    nb = numbers()
    g1 = tagged(r"\int_{-\infty}^{\infty}e^{-ax^2+bx}\,dx", "(2:S1-3)")
    g2 = tagged(r"ax^2-bx=a\Big(x-\frac{b}{2a}\Big)^2-\frac{b^2}{4a}", "(2:S1-4)")
    top = VGroup(g1, g2).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    top.move_to(np.array([-6.6 + top.get_width() / 2, 3.15 - top.get_height() / 2, 0]))
    box = D.erratum_box("p. 41，(2:S1-4)", r"-\frac{b^2}{2a}", r"-\frac{b^2}{4a}", width=4.6)
    ex1 = mixed("$a=1,\\ b=2$", "：", "$x^2-2x=(x-1)^2-1$", font_size=24)
    ex2 = mixed("照印：", "$(x-1)^2-2$", "（", "$x=0$", " 就不對）", font_size=24, color=D.WRONG_COLOR)
    err = VGroup(box, ex1, ex2).arrange(DOWN, buff=0.15)
    err.move_to(np.array([6.7 - err.get_width() / 2, 3.15 - err.get_height() / 2, 0]))
    ax = D.plot_axes((-2.5, 4, 1), (0, 3, 1), width=7.0, height=2.9, x_ticks=[-2, -1, 1, 2, 3], y_ticks=[1, 2],
                     font_size=20, x_label="x")
    ax.move_to(np.array([-2.6, -1.25, 0]))
    f0 = lambda x: math.exp(-x * x)
    f1 = lambda x: math.exp(-(x - 1) ** 2)
    f2 = lambda x: math.exp(-x * x + 2 * x)
    c0 = ax.get_graph(f0, x_range=[-2.5, 4, 0.02]).set_stroke(D.BAR_COLOR, 3)
    c1 = ax.get_graph(f1, x_range=[-2.5, 4, 0.02]).set_stroke(GREY_B, 3)
    c2 = ax.get_graph(f2, x_range=[-1.3, 3.3, 0.02]).set_stroke(D.GAUSS_COLOR, 4)
    l2 = tx(r"e^{-x^2+2x}=e\cdot e^{-(x-1)^2}", 26, D.GAUSS_COLOR).next_to(ax.c2p(1, math.e), RIGHT, buff=0.3)
    res = tagged(r"\int e^{-ax^2+bx}\,dx=e^{b^2/4a}\sqrt{\frac{\pi}{a}}", "(2:S1-5)", 32, YELLOW)
    res2 = mixed("$a=1,\\ b=2$", "：", "$e\\sqrt\\pi=" + nb["gen"] + "$", font_size=28)
    right = VGroup(res, res2).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    right.move_to(np.array([6.7 - right.get_width() / 2, -1.75, 0]))
    assert_in_frame(top, err, ax, l2, right)
    assert_apart(top, err)
    assert_apart(VGroup(top, err), ax)
    assert_apart(ax, right)
    assert_apart(l2, right)
    s.say("在指數裡多加一項 ", "$bx$", "：高斯積分的生成函數", wait=0.3)
    s.play(FadeIn(g1), run_time=1.0)
    s.play(FadeIn(g2), run_time=1.0)
    s.say("配方：把 ", "$bx$", " 吸進平方裡", wait=3.0)
    s.play(FadeIn(err, shift=0.2 * LEFT), run_time=0.8)
    s.say("書上的常數項印成 ", "$\\frac{b^2}{2a}$", "：應是 ", "$\\frac{b^2}{4a}$", wait=3.5)
    s.say("代 ", "$x=0$", " 驗算：左邊是 0，照印的右邊卻是 ", "$-1$", wait=3.5)
    s.play(FadeIn(ax), ShowCreation(c0), run_time=1.2)
    s.say("看一個例子：", "$a=1,\\ b=2$", "，從 ", "$e^{-x^2}$", " 開始", wait=2.5)
    s.play(TransformFromCopy(c0, c1), run_time=1.5)
    s.say("往右移 1：", "$e^{-(x-1)^2}$", "，面積不變", wait=2.5)
    s.play(TransformFromCopy(c1, c2), FadeIn(l2), run_time=1.5)
    s.say("再乘上 ", "$e^{1}$", "：面積也乘 ", "$e$", wait=2.5)
    s.play(FadeIn(right, shift=0.2 * UP), run_time=1.0)
    s.say("所以 ", "$\\int e^{-x^2+2x}dx=e\\sqrt\\pi=" + nb["gen"] + "$", wait=3.5)
    s.say("對任何 ", "$a>0$", " 和 ", "$b$", " 都成立：它就是動差的生成函數", wait=3.0)
    s.clear_stage()


def part_moments(s):
    nb = numbers()
    e1 = tx(r"\int e^{-ax^2}e^{bx}\,dx=\sum_{n}\frac{b^n}{n!}\int x^ne^{-ax^2}\,dx", 30)
    e2 = tx(r"e^{b^2/4a}\sqrt{\frac{\pi}{a}}=\sqrt{\frac{\pi}{a}}\sum_{k}\frac{b^{2k}}{k!\,(4a)^k}", 30)
    e3 = VGroup(boxed(tx(r"\int x^ne^{-ax^2}\,dx=\sqrt{\frac{\pi}{a}}\,\frac{n!}{(n/2)!\,(4a)^{n/2}}", 34), color=YELLOW,
                      buff=0.15), mixed("（", "$n$", " 偶數；奇數是 0）", font_size=24), tx(r"\text{(2:S1-6)}", 22, GREY_B))
    e3.arrange(RIGHT, buff=0.25)
    col = VGroup(VGroup(e1, e2).arrange(RIGHT, buff=0.7), e3).arrange(DOWN, buff=0.4)
    fit(col, 13.0).move_to(np.array([0, 3.15 - col.get_height() / 2, 0]))
    rows = [[tx("n", 28, YELLOW)] + [str(n) for n in NS],
            [tx(r"\frac{n!}{(n/2)!\,4^{n/2}}", 26)] + nb["frac"],
            [mixed("公式（", "$a=1$", "）", font_size=22)] + nb["mom"],
            [zh("數值積分", font_size=22)] + nb["num"]]
    tab = D.tex_table(rows, [2.3, 1.5, 1.5, 1.5, 1.5], row_h=0.62, font_size=28)
    tab.move_to(np.array([-6.6 + tab.get_width() / 2, -1.45, 0]))
    hide = VGroup(D.table_cell(tab, 2, 3, 5), D.table_cell(tab, 3, 3, 5))
    q = tx("?", 34, YELLOW).move_to(hide)
    ans = VGroup(tx(r"\int_{-\infty}^{\infty}x^4e^{-x^2}dx", 30, YELLOW),
                 tx(r"=\tfrac34\sqrt\pi\approx" + nb["mom"][2], 30, YELLOW)).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
    ans.next_to(tab, RIGHT, buff=0.5)
    others = VGroup(*[c for c in tab[0] if c not in hide.submobjects])
    assert_in_frame(col, tab, ans)
    assert_apart(col, tab)
    assert_apart(tab, ans)
    s.say("兩邊都展開成 ", "$b$", " 的冪級數", wait=0.3)
    s.play(FadeIn(e1), run_time=1.2)
    s.say("左邊：", "$e^{bx}$", " 展開，每一項帶一個動差", wait=3.0)
    s.play(FadeIn(e2), run_time=1.2)
    s.say("右邊：", "$e^{b^2/4a}$", " 只有 ", "$b$", " 的偶次方", wait=3.0)
    s.play(FadeIn(e3, shift=0.2 * UP), run_time=1.2)
    s.say("比較 ", "$b^n$", " 的係數：所有動差一次到手", wait=3.0)
    s.say("$n$", " 是奇數：被積函數是奇函數，積分是 0", wait=2.5)
    s.play(FadeIn(others), ShowCreation(tab[1]), FadeIn(q), run_time=1.2)
    s.say("$a=1$", "：", "$n=0,2,6$", " 都和數值積分一致", wait=3.0)
    s.say("猜猜看：", "$n=4$", " 是多少？", wait=2.5)
    s.play(FadeOut(q), FadeIn(hide, scale=1.3), run_time=0.8)
    s.play(FadeIn(ans), run_time=0.8)
    s.say("$\\tfrac34\\sqrt{\\pi}\\approx" + nb["mom"][2] + "$", "：開頭的問題答案", wait=3.0)
    s.clear_stage()


def part_echo(s):
    r1 = tx(r"\frac{\langle x^4\rangle}{\langle x^2\rangle^2}=\frac{3/4}{(1/2)^2}=3", 36)
    l1 = zh("W105：高斯的峰度是 3", font_size=26, color=GREY_B)
    r2 = tx(r"\frac{\int q^4e^{-\alpha q^2}dq}{\int e^{-\alpha q^2}dq}=\frac{3}{4\alpha^2},\quad \alpha=\frac{a^2N}{6}"
            r"\ \Rightarrow\ bN\cdot\frac{27}{a^4N^2}=\frac{27b}{Na^4}", 32)
    l2 = mixed("W207：(2.22) 的 ", "$q^4$", " 修正", font_size=26, color=GREY_B)
    r3 = tx(r"\int e^{-\frac N6q^2-iqr}\,dq=\sqrt{\frac{6\pi}{N}}\,e^{-3r^2/2N}", 32)
    l3 = mixed("W206：(2.20) 每個方向一個，", "$a=\\frac N6,\\ b=-ir$", font_size=26, color=GREY_B)
    col = VGroup(r1, l1, r2, l2, r3, l3).arrange(DOWN, buff=0.28)
    for m in col:
        fit(m, 13.0)
    col.arrange(DOWN, buff=0.28).move_to(np.array([0, 0.2, 0]))
    assert_in_frame(col)
    assert_apart(l1, r2)
    assert_apart(l2, r3)
    s.say("這些積分在前面都出現過", wait=0.3)
    s.play(FadeIn(r1), FadeIn(l1), run_time=1.2)
    s.say("四階動差除以二階動差的平方：3", wait=3.0)
    s.play(FadeIn(r2), FadeIn(l2), run_time=1.2)
    s.say("W207 的 ", "$q^4$", " 修正，用的就是四階動差", wait=3.5)
    s.play(FadeIn(r3), FadeIn(l3), run_time=1.2)
    s.say("(2.20) 的高斯：", "$b$", " 取成虛數，", "$e^{b^2/4a}$", " 就是 ", "$e^{-3r^2/2N}$", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("一個生成函數 ", "$e^{b^2/4a}\\sqrt{\\pi/a}$", wait=2.5)
    s.say("一次給出高斯的所有動差", wait=3.0)


class W210GaussianIntegrals(ChuScene):
    def construct(self):
        self.show_header("W210　補充 1：高斯積分與它的生成函數")
        check_math()
        part_question(self)
        part_basic(self)
        part_generating(self)
        part_moments(self)
        part_echo(self)
        part_close(self)
        self.wait(1.0)
