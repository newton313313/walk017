"""
W209 擴散方程與格林函數（§2.1.5，pp. 38–40，(2.32)–(2.38)，Ex. 2.8–2.9）

問：一滴墨水滴在格子中央，機率怎麼散開？
  1. (2.32) P(N; x, y) = (1/2d) Σ_nns P(N−1; x, y − w_i)；偏向的走法 (2.33) 把 1/2d 換成 p_i。
     正方格子 L.heat_cells(GF.prob_grid_sq(N))，N = 0..20：棋盤格一路擴散；Σ r²P = N；P(20; 0) = 0.031045 vs 2/(πN) = 0.031831。
  2. (2.34)(2.35) 加減 P(N−1; y) → 二階差分 → (2.36)(2.37) ∂P/∂N = (a²/6)∇²P（d 維 a²/2d）；
     書上說解是「(2.20) 的左邊去掉 6^N」→ 應為右邊（字幕帶過）。
  3. Ex. 2.8：P = (3/(2πNa²))^{3/2} e^{−3r²/2Na²} 滿足擴散方程（GF.ex28_sides；N = 10、r = 2 兩邊都是 −5.153×10⁻⁴），∫P d³x = 1。
  4. (2.38)：(1 − 2dz)G − za²∇²G = δ；筆誤框（∇² 項的正負號）；一維傅立葉核對 z = 0.4、q = 0.3：4.24213、4.23729、6.09756。
熱圖的顏色每一張按自己的最大值（看形狀）。畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
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
from utils import lattice_draw as L
from utils import walk_draw as D

NMAX = 20
CELL = 0.135
GRID_O = np.array([-3.55, 0.2, 0])
ZQ = (0.4, 0.3)


def p_ex28(N, r):
    return (3 / (2 * math.pi * N)) ** 1.5 * math.exp(-3 * r * r / (2 * N))


def sci(v):
    m, e = f"{v:.3e}".split("e")
    return m + r"\times10^{" + str(int(e)) + "}"


def numbers():
    P = GF.prob_grid_sq(NMAX)
    dN, lap = GF.ex28_sides(10, 2.0)
    z, q = ZQ
    return {"p0": f"{float(P[(0, 0)]):.6f}", "asym": f"{2 / (math.pi * NMAX):.6f}", "var": str(GF.variance(P)),
            "ex28": (sci(dN), sci(lap)),
            "g": [f"{GF.g_fourier(z, GF.chi_1d(q)):.5f}", f"{GF.g_cont(z, q):.5f}", f"{GF.g_cont_book(z, q):.5f}"]}


def check_math():
    nb = numbers()
    assert nb["p0"] == "0.031045" and nb["asym"] == "0.031831" and nb["var"] == "20"
    P = {(0, 0): Fraction(1)}
    for N in range(NMAX + 1):
        G = GF.prob_grid_sq(N)
        assert G == P and GF.variance(G) == N and sum(G.values()) == 1          # (2.32) 一步步算 = 計數 / 4^N
        assert all((x + y - N) % 2 == 0 for (x, y) in G) and max(abs(x) + abs(y) for (x, y) in G) == N
        P = GF.diffuse_sq(P)
    dN, lap = GF.ex28_sides(10, 2.0)
    assert abs(dN - lap) < 1e-15 and nb["ex28"] == (r"-5.153\times10^{-4}", r"-5.153\times10^{-4}")
    for N, r in ((5, 0.0), (10, 2.0), (30, 3.5)):                         # ∂P/∂N 的解析式與差分核對
        h = 1e-5
        num = (p_ex28(N + h, r) - p_ex28(N - h, r)) / (2 * h)
        assert abs(num - GF.ex28_sides(N, r)[0]) < 1e-9
        assert abs(GF.ex28_sides(N, r)[0] - (-1.5 / N + 1.5 * r * r / N ** 2) * p_ex28(N, r)) < 1e-15
    assert abs(GF.radial_integral(lambda r, r0: p_ex28(10, r)) - 1) < 1e-9
    assert nb["g"] == ["4.24213", "4.23729", "6.09756"]


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def tagged(tex, tag, font_size=30, color=WHITE):
    return VGroup(tx(tex, font_size, color), tx(r"\text{" + tag + "}", 22, GREY_B)).arrange(RIGHT, buff=0.3)


def cells(N):
    P = GF.prob_grid_sq(N)
    return L.heat_cells(P, cell=CELL, origin=GRID_O, vmax=float(max(P.values())))


def counter(N):
    t = tx(r"N=" + str(N) + r"\qquad\sum_{\vec r}r^2P(N;\vec r)=" + str(GF.variance(GF.prob_grid_sq(N))), 30)
    return t.move_to(np.array([3.35, 0.45, 0]))


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：一滴墨水滴在格子中央，機率怎麼散開？", wait=3.0)


def part_spread(s):
    nb = numbers()
    e1 = tagged(r"P(N;\vec x,\vec y)=\frac{1}{2d}\sum_{nns}P(N-1;\vec x,\vec y-\vec w_i)", "(2.32)")
    fit(e1, 6.6).move_to(np.array([3.35, 2.55, 0]))
    e2 = mixed("偏向的走法 (2.33)：", "$\\frac{1}{2d}\\ \\to\\ p_i$", "（往各方向的機率）", font_size=24)
    fit(e2, 6.6).move_to(np.array([3.35, 1.55, 0]))
    frame = Square(side_length=(2 * NMAX + 1) * CELL + 0.1).move_to(GRID_O).set_stroke(GREY_D, 1.5)
    note = zh("顏色：每一張按自己的最大值", font_size=20, color=GREY_B).next_to(frame, DOWN, buff=0.08)
    end1 = mixed("$P(20;\\vec 0)=" + nb["p0"] + "$", font_size=30)
    end2 = mixed("高斯（含奇偶的因子 2）", "$\\dfrac{2}{\\pi N}=" + nb["asym"] + "$", font_size=28)
    ends = VGroup(end1, end2).arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to(np.array([3.35, -1.15, 0]))
    cur = cells(0)
    cnt = counter(0)
    assert_in_frame(e1, e2, frame, note, ends, cnt, counter(NMAX))
    assert_apart(frame, VGroup(e1, e2, ends, counter(NMAX)))
    assert_apart(e2, cnt)
    assert_apart(cnt, ends)
    s.play(FadeIn(e1), run_time=1.2)
    s.say("$P$", "＝走到那裡的機率：從四個鄰居各分到 ", "$\\frac14$", wait=3.0)
    s.say("$d$", " 維有 ", "$2d$", " 個鄰居；配位數 ", "$2d$", wait=2.5)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("如果往各方向的機率不同：", "$\\frac{1}{2d}$", " 換成 ", "$p_i$", wait=3.0)
    s.play(ShowCreation(frame), FadeIn(cur), FadeIn(cnt), FadeIn(note), run_time=1.0)
    s.say("$N=0$", "：全部的機率都在原點", wait=2.0)
    s.say("一步一步：每一格平分給四個鄰居", wait=0.2)
    for N in range(1, NMAX + 1):
        new, ncnt = cells(N), counter(N)
        s.play(FadeOut(cur), FadeIn(new), FadeOut(cnt), FadeIn(ncnt), run_time=0.45)
        cur, cnt = new, ncnt
    s.say("棋盤格：", "$x+y$", " 與 ", "$N$", " 不同奇偶的格子永遠是空的", wait=3.0)
    s.say("$\\sum r^2P=N$", "：散開的面積和步數成正比", wait=3.0)
    s.say("散得越開，原點分到的機率就越少", wait=2.5)
    s.play(FadeIn(ends), run_time=1.0)
    s.say("原點的機率：", "$" + nb["p0"] + "$", "，高斯近似 ", "$" + nb["asym"] + "$", wait=3.5)
    s.clear_stage()


def part_equation(s):
    e1 = tagged(r"P(N;\vec y)=\frac16\sum_{i=1}^{6}P(N-1;\vec y-\vec w_i)", "(2.32)")
    e1.move_to(np.array([0, 2.7, 0]))
    lhs = tx(r"P(N;\vec y)-P(N-1;\vec y)", 30)
    mid = tx(r"=\frac16\sum_{j=1}^{3}", 30)
    br = tx(r"\Big[P(N-1;\vec y-\hat e_ja)+P(N-1;\vec y+\hat e_ja)-2P(N-1;\vec y)\Big]", 30)
    e2 = VGroup(lhs, mid, br).arrange(RIGHT, buff=0.12)
    fit(e2, 13.2).move_to(np.array([0, 1.45, 0]))
    b1 = Brace(lhs, DOWN, buff=0.1)
    l1 = VGroup(tx(r"\approx\frac{\partial P}{\partial N}", 30), tx(r"\text{(2.36)}", 22, GREY_B))
    l1.arrange(RIGHT, buff=0.25).next_to(b1, DOWN, buff=0.1)
    b2 = Brace(br, DOWN, buff=0.1)
    l2 = VGroup(tx(r"\approx a^2\frac{\partial^2P}{\partial y_j^2}", 30), zh("（二階差分）", font_size=24))
    l2.arrange(RIGHT, buff=0.15).next_to(b2, DOWN, buff=0.1)
    res = boxed(tx(r"\frac{\partial P}{\partial N}=\frac{a^2}{6}\nabla^2P", 44), color=YELLOW, buff=0.2)
    rest = VGroup(res, tx(r"\text{(2.37)}", 22, GREY_B)).arrange(RIGHT, buff=0.3).move_to(np.array([-1.2, -1.85, 0]))
    dd = mixed("$d$", " 維：", "$\\dfrac{a^2}{2d}\\nabla^2P$", font_size=28).next_to(rest, RIGHT, buff=0.8)
    assert_in_frame(e1, e2, l1, l2, rest, dd)
    assert_apart(e1, e2)
    assert_apart(l1, l2)
    assert_apart(VGroup(l1, l2), rest)
    s.say("從遞迴到微分方程：以立方格子為例（", "$d=3$", "）", wait=0.3)
    s.play(FadeIn(e1), run_time=1.0)
    s.say("兩邊減去 ", "$P(N-1;\\vec y)$", "：把六個鄰居兩兩配對", wait=2.5)
    s.play(FadeIn(e2), run_time=1.2)
    s.say("右邊每一對：左鄰＋右鄰－2×自己", wait=3.0)
    s.play(GrowFromCenter(b2), FadeIn(l2), run_time=1.0)
    s.say("這是二階差分：格距 ", "$a$", " 很小時 ", "$\\approx a^2$", " 乘二階導數", wait=3.5)
    s.play(GrowFromCenter(b1), FadeIn(l1), run_time=1.0)
    s.say("左邊：多走一步的變化 ", "$\\approx\\partial P/\\partial N$", wait=3.0)
    s.play(FadeIn(rest, shift=0.2 * UP), run_time=1.0)
    s.play(FadeIn(dd), run_time=0.6)
    s.say("擴散方程！", "$N$", " 扮演時間的角色", wait=3.0)
    s.say("書上說解是 (2.20) 的「左邊」去掉 ", "$6^N$", "：應是右邊", wait=3.5)
    s.clear_stage()


def part_ex28(s):
    nb = numbers()
    f1 = mixed("練習 2.8：", "$P=\\Big(\\dfrac{3}{2\\pi Na^2}\\Big)^{3/2}e^{-3r^2/2Na^2}$", font_size=32)
    f2 = tx(r"\frac{\partial P}{\partial N}=\Big(-\frac{3}{2N}+\frac{3r^2}{2N^2a^2}\Big)P", 34)
    f3 = tx(r"\frac{a^2}{6}\nabla^2P=\frac{a^2}{6}\Big(P''+\frac{2}{r}P'\Big)=\Big(-\frac{3}{2N}+\frac{3r^2}{2N^2a^2}\Big)P", 34)
    f4 = mixed("$N=10,\\ r=2,\\ a=1$", "：兩邊都是 ", "$" + nb["ex28"][0] + "$", font_size=28)
    f5 = mixed("而且 ", "$\\int P\\,d^3x=1$", "：機率守恆", font_size=28)
    col = VGroup(f1, f2, f3, f4, f5).arrange(DOWN, aligned_edge=LEFT, buff=0.42)
    for m in col:
        fit(m, 12.6)
    col.arrange(DOWN, aligned_edge=LEFT, buff=0.42).move_to(np.array([0, 0.25, 0]))
    assert_in_frame(col)
    s.say("驗證：(2.20) 的高斯滿足擴散方程", wait=0.3)
    s.play(FadeIn(f1), run_time=1.0)
    s.say("把 ", "$6^N$", " 拿掉、變成機率：這就是 (2.20) 的右邊", wait=3.0)
    s.play(FadeIn(f2), run_time=1.0)
    s.say("對 ", "$N$", " 微分：前因子與指數各給一項", wait=3.0)
    s.play(FadeIn(f3), run_time=1.0)
    s.say("三維的拉普拉斯（球對稱）：結果一模一樣", wait=3.0)
    s.play(FadeIn(f4), FadeIn(f5), run_time=1.0)
    s.say("代數字核對：", "$N=10,\\ r=2$", " 兩邊相等", wait=3.0)
    s.say("總機率守恆：墨水只是散開，不會消失", wait=3.0)
    s.clear_stage()


def part_green(s):
    nb = numbers()
    g1 = tagged(r"G(z;\vec y)-z\sum_iG(z;\vec y-\vec w_i)=\delta_{\vec y,\vec 0}", "(2.9)")
    g2 = tx(r"\sum_iG(\vec y-\vec w_i)\approx 2d\,G+a^2\nabla^2G", 30)
    g3 = tagged(r"(1-2dz)\,G-za^2\nabla^2G=\delta", "(2.38)", 34, YELLOW)
    col = VGroup(g1, g2, g3).arrange(DOWN, buff=0.3).move_to(np.array([0, 2.05, 0]))
    rows = [[zh("一維、傅立葉", font_size=22, color=YELLOW), tx(r"g(z;q)", 26, YELLOW), tx(r"z=0.4,\ q=0.3", 24, YELLOW)],
            [zh("精確", font_size=22), r"\frac{1}{1-2z\cos q}", nb["g"][0]],
            [zh("連續極限（正確）", font_size=22, color=D.RIGHT_COLOR), r"\frac{1}{1-2z+zq^2}", nb["g"][1]],
            [zh("照書上的號", font_size=22, color=D.WRONG_COLOR), r"\frac{1}{1-2z-zq^2}", nb["g"][2]]]
    tab = D.tex_table(rows, [3.0, 2.8, 2.4], row_h=0.76, font_size=28)
    tab.move_to(np.array([0, -1.35, 0]))
    D.table_cell(tab, 3, 2, 3).set_color(D.WRONG_COLOR)
    box = D.erratum_box("p. 40，(2.38)", r"za^2\nabla^2G-(2dz-1)G=\delta", r"za^2\nabla^2G+(2dz-1)G=-\delta",
                        width=10.5, font_size=34)
    box.move_to(np.array([0, -1.0, 0]))
    assert_in_frame(col, tab, box)
    assert_apart(col, tab)
    assert_apart(col, box)
    s.say("生成函數也一樣：(2.9) 取連續極限", wait=0.3)
    s.play(FadeIn(g1), run_time=1.0)
    s.play(FadeIn(g2), run_time=1.0)
    s.say("鄰居的和 ", "$\\approx 2d\\,G$", " 加上二階差分", wait=3.0)
    s.play(FadeIn(g3, shift=0.2 * UP), run_time=1.0)
    s.say("$G$", " 是擴散算子的格林函數：右邊是點源 ", "$\\delta$", wait=3.5)
    s.play(FadeIn(box, shift=0.2 * LEFT), run_time=0.8)
    s.say("書上 (2.38) 的 ", "$\\nabla^2$", " 項正負號弄反了", wait=3.0)
    s.play(FadeOut(box), FadeIn(tab), run_time=1.0)
    s.say("傅立葉核對：", "$\\nabla^2\\to-q^2$", "；正確版和精確值只差 0.1%", wait=3.5)
    s.say("照書上的號：差了四成多", wait=2.5)
    s.clear_stage()


def part_close(s):
    s.say("隨機漫步放大來看就是擴散", wait=2.5)
    s.say("生成函數 ", "$G$", " 是擴散算子的格林函數", wait=3.0)


class W209DiffusionGreen(ChuScene):
    def construct(self):
        self.show_header("W209　擴散方程與格林函數")
        check_math()
        part_question(self)
        part_spread(self)
        part_equation(self)
        part_ex28(self)
        part_green(self)
        part_close(self)
        self.wait(1.0)
