"""
W204 平移對稱與傅立葉：一行代數解出 G（§2.1.2，pp. 31–33，(2.10)–(2.13)，Fig. 2.4，Ex. 2.6）

問：(2.9) 每一格都牽扯到鄰居，怎麼一次解出全部？
  1. 平面波（一維格點）f(x) = cos(qx)，q = π/3：x = −3..3 的值 −1, −1/2, 1/2, 1, 1/2, −1/2, −1；
     左右鄰居的和 f(x−1) + f(x+1) = 2cos q·f(x) = f(x)（χ(π/3) = 1）：「和鄰居的關係」變成「乘上 χ(q)」。
  2. 平移對稱 → (2.10) g(z;q) = Σ_x G(z;x−y) e^{iq·(x−y)}；把 (2.9) 乘上相位、對 x 加總 (2.11)：g = zgχ + 1
     → (2.13) g = 1/(1 − zχ(q))，結構函數 χ(q) = Σ_i e^{iq·(w_i−y)}。
  3. 各種格子的 χ（Fig. 2.4）：一維、正方、sc (2.12)、bcc、fcc（練習 2.6）；χ(0)＝鄰居數 2, 4, 6, 8, 12。
  4. 正方格子 χ 的熱圖：中心 4、四個角 −4（伏筆）；q = 0：g = 1/(1−4z) = Σ 4^N z^N。
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
from utils import lattice_draw as L
from utils import walk_draw as D
from utils import walk_gf as GF

QW = math.pi / 3                                         # 平面波的 q
XS = list(range(-3, 4))
R2, R3 = math.sqrt(2), math.sqrt(3)
LATS = [("1d", "一維", r"2\cos q"),
        ("sq", "正方", r"2(\cos q_x+\cos q_y)"),
        ("sc", "簡單立方 (2.12)", r"2(\cos q_xa+\cos q_ya+\cos q_za)"),
        ("bcc", "體心立方 bcc", r"8\cos\tfrac{q_xa}{\sqrt3}\cos\tfrac{q_ya}{\sqrt3}\cos\tfrac{q_za}{\sqrt3}"),
        ("fcc", "面心立方 fcc", r"4\Big[\cos\tfrac{q_xa}{\sqrt2}\cos\tfrac{q_ya}{\sqrt2}+\cos\tfrac{q_ya}{\sqrt2}"
                                r"\cos\tfrac{q_za}{\sqrt2}+\cos\tfrac{q_za}{\sqrt2}\cos\tfrac{q_xa}{\sqrt2}\Big]")]
CHI0 = {"1d": lambda: GF.chi_1d(0), "sq": lambda: GF.chi_sq(0, 0), "sc": lambda: GF.chi_sc(0, 0, 0),
        "bcc": lambda: GF.chi_bcc(0, 0, 0), "fcc": lambda: GF.chi_fcc(0, 0, 0)}


def frac_tex(c) -> str:
    c = Fraction(c)
    if c.denominator == 1:
        return str(c.numerator)
    return ("-" if c < 0 else "") + rf"\tfrac{{{abs(c.numerator)}}}{{{c.denominator}}}"


def wave_vals():
    return [Fraction(round(2 * math.cos(QW * x)), 2) for x in XS]


def chi0(kind) -> str:
    return str(int(round(CHI0[kind]())))


def total_walks(nmax):
    return [sum(GF.count_grid(N).values()) for N in range(nmax + 1)]


def check_math():
    F = Fraction
    v = wave_vals()
    assert v == [-1, F(-1, 2), F(1, 2), 1, F(1, 2), F(-1, 2), -1]
    assert all(abs(math.cos(QW * x) - float(f)) < 1e-12 for x, f in zip(XS, v))
    assert GF.plane_wave_check(QW) < 1e-12 and abs(GF.chi_1d(QW) - 1) < 1e-12
    f = dict(zip(XS, v))
    assert f[0] + f[2] == f[1] == F(1, 2)                                    # 例 x = 1
    assert [chi0(k) for k, _, _ in LATS] == ["2", "4", "6", "8", "12"]
    q3 = (0.3, 0.7, -1.1)
    assert abs(GF.chi_sc(*q3) - GF.chi_from_steps(q3, GF.SC_STEPS)) < 1e-12
    assert abs(GF.chi_bcc(*q3) - GF.chi_from_steps(q3, GF.BCC_STEPS, 1 / R3)) < 1e-12
    assert abs(GF.chi_fcc(*q3) - GF.chi_from_steps(q3, GF.FCC_STEPS, 1 / R2)) < 1e-12
    assert abs(GF.chi_sq(0.4, -0.9) - GF.chi_from_steps((0.4, -0.9), GF.SQ_STEPS)) < 1e-12
    assert abs(GF.chi_sq(math.pi, math.pi) + 4) < 1e-12 and GF.chi_sq(0, 0) == 4
    z, big = 0.05, GF.count_grids(20)                                         # (2.13)：Σ_x G e^{iq·x} = 1/(1−zχ)
    qx, qy = 0.4, -0.9
    g = sum(z ** N * c * math.cos(qx * x + qy * y) for N in range(21) for (x, y), c in big[N].items())
    assert abs(g - GF.g_fourier(z, GF.chi_sq(qx, qy))) < 1e-10
    assert total_walks(4) == [1, 4, 16, 64, 256] == [4 ** N for N in range(5)]


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


# ---------------------------------------------------------------- 各段
def part_question(s):
    e = boxed(eq(r"G(z;\vec x,\vec y)=z\sum_{i=1}^{4}G(z;\vec x,\vec w_i)+\delta_{\vec x,\vec y}", "(2.9)", font_size=42),
              color=D.GAUSS_COLOR, buff=0.25).move_to(np.array([0, 0.4, 0]))
    assert_in_frame(e)
    s.play(FadeIn(e, shift=0.2 * UP), run_time=1.0)
    s.say("問：(2.9) 每一格都牽扯到鄰居，怎麼一次解出全部？", wait=4.0)
    s.clear_stage()


def part_plane_wave(s):
    v = wave_vals()
    ax = D.plot_axes((-3.6, 3.6, 1), (-1.6, 1.6, 0.5), width=6.2, height=3.4, x_ticks=XS, font_size=20, x_label="x")
    ax.move_to(np.array([-3.1, 0.75, 0]))
    curve = DashedVMobject(ax.get_graph(lambda x: math.cos(QW * x), x_range=[-3.5, 3.5, 0.02]), num_dashes=40)
    curve.set_stroke(GREY_B, 2)
    stems = VGroup(*[Line(ax.c2p(x, 0), ax.c2p(x, float(f))).set_stroke(D.BAR_COLOR, 3) for x, f in zip(XS, v)])
    dots = VGroup(*[Dot(ax.c2p(x, float(f)), radius=0.08).set_fill(D.BAR_COLOR, 1) for x, f in zip(XS, v)])
    labs = VGroup(*[Tex(frac_tex(f), font_size=30).next_to(d, UP if f > 0 else DOWN, buff=0.1)
                    for d, f in zip(dots, v)])
    i0, i1, i2 = XS.index(0), XS.index(1), XS.index(2)
    arrs = VGroup(*[Arrow(dots[i].get_center(), dots[i1].get_center(), buff=0.12, thickness=2.5).set_color(YELLOW)
                    for i in (i0, i2)])
    p1 = Tex(r"f(x)=\cos qx,\quad q=\frac{\pi}{3}", font_size=34)
    p2 = Tex(r"f(x-1)+f(x+1)=2\cos q\cdot f(x)", font_size=34)
    p3 = Tex(rf"\chi(q)=2\cos q,\quad \chi\big(\tfrac{{\pi}}{{3}}\big)={round(GF.chi_1d(QW))}", font_size=34)
    p4 = Tex(rf"x=1:\ \ f(0)+f(2)={frac_tex(v[i0])}{frac_tex(v[i2])}={frac_tex(v[i1])}=f(1)", font_size=32)
    p4.set_color(YELLOW)
    p5 = mixed("平移一格只是乘上一個相位：", font_size=24, color=GREY_B)
    p6 = mixed("「和鄰居的關係」變成「乘上 ", "$\\chi(q)$", "」", font_size=26)
    col = VGroup(p1, p2, p3, p4, p5, p6)
    for m in col:
        fit(m, 5.7)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.38), 0.9, 2.7)
    assert_in_frame(ax, labs, col)
    assert_apart(ax, col)
    assert_apart(labs, col)

    s.say("先看一維格點上的平面波 ", "$\\cos qx$", "，取 ", "$q=\\pi/3$", wait=0.3)
    s.play(FadeIn(ax), ShowCreation(curve), FadeIn(p1), run_time=1.5)
    s.play(LaggedStart(*[ShowCreation(st) for st in stems], lag_ratio=0.1),
           LaggedStart(*[FadeIn(d) for d in dots], lag_ratio=0.1), run_time=1.5)
    s.play(FadeIn(labs), run_time=0.8)
    s.say("只看整數點：", "$-1,-\\frac12,\\frac12,1,\\frac12,-\\frac12,-1$", wait=3.0)
    s.play(FadeIn(p2), run_time=0.8)
    s.say("左右鄰居加起來，正好是自己的 ", "$2\\cos q$", " 倍", wait=3.0)
    s.play(FadeIn(p3), run_time=0.8)
    s.say("$q=\\pi/3$", " 時 ", "$\\chi=1$", "：鄰居的和就是自己", wait=2.5)
    s.play(dots[i0].animate.set_fill(YELLOW), dots[i2].animate.set_fill(YELLOW), *[GrowArrow(a) for a in arrs],
           FadeIn(p4), run_time=1.0)
    s.say("例：", "$x=1$", " 的兩個鄰居 ", "$1$", " 和 ", "$-\\frac12$", "，和是 ", "$\\frac12$", wait=3.0)
    s.play(FadeIn(p5), FadeIn(p6), run_time=0.8)
    s.say("每一點都一樣：平移只乘上相位", wait=3.0)
    s.say("「和鄰居的關係」變成「乘上 ", "$\\chi(q)$", "」", wait=3.0)
    s.clear_stage()


def part_derivation(s):
    d0 = mixed("平移對稱：", "$G(z;\\vec x,\\vec y)$", " 只依賴 ", "$\\vec x-\\vec y$", font_size=24)
    d1 = eq(r"g(z;\vec q)=\sum_{\vec x}G(z;\vec x-\vec y)\,e^{i\vec q\cdot(\vec x-\vec y)}", "(2.10)", font_size=36)
    d2 = mixed("把 (2.9) 乘上 ", "$e^{i\\vec q\\cdot(\\vec x-\\vec y)}$", "、對 ", "$\\vec x$", " 加總，相位拆成兩段：",
               font_size=24)
    d3 = Tex(r"e^{i\vec q\cdot(\vec x-\vec y)}=\underbrace{e^{i\vec q\cdot(\vec x-\vec w_i)}}_{\to\ g}"
             r"\ \underbrace{e^{i\vec q\cdot(\vec w_i-\vec y)}}_{\to\ \chi(\vec q)}", font_size=36)
    d4 = eq(r"g=z\,g\,\chi(\vec q)+1", "(2.11)", font_size=40)
    d5 = boxed(eq(r"g(z;\vec q)=\frac{1}{1-z\chi(\vec q)}", "(2.13)", font_size=44), color=D.GAUSS_COLOR, buff=0.2)
    d6 = mixed("結構函數 ", "$\\chi(\\vec q)=\\sum_{i}e^{i\\vec q\\cdot(\\vec w_i-\\vec y)}$", font_size=24)
    col = VGroup(d0, d1, d2, d3, d4, d5, d6).arrange(DOWN, buff=0.16).move_to(np.array([0, 0.1, 0]))
    assert_in_frame(col)

    s.say("格子處處一樣：", "$G$", " 只依賴兩點的差", wait=0.3)
    s.play(FadeIn(d0), run_time=0.8)
    s.play(Write(d1[0]), FadeIn(d1[1]), run_time=1.8)
    s.say("(2.10)：對所有差 ", "$\\vec x-\\vec y$", " 加上相位，就是傅立葉轉換", wait=3.0)
    s.play(FadeIn(d2), run_time=0.8)
    s.play(Write(d3), run_time=1.6)
    s.say("第一段：", "$\\vec x$", " 跑遍格點時 ", "$\\vec x-\\vec w_i$", " 也跑遍：又得到 ", "$g$", wait=3.5)
    s.say("第二段只看鄰居的位置：加起來是 ", "$\\chi(\\vec q)$", wait=3.0)
    s.play(FadeIn(d4, shift=0.2 * DOWN), run_time=1.0)
    s.say("右邊的 ", "$\\delta$", " 變成 1：", "$g=zg\\chi+1$", wait=3.5)
    s.play(FadeIn(d5, shift=0.2 * UP), run_time=1.0)
    s.say("一行代數：", "$g=1/(1-z\\chi)$", wait=3.5)
    s.play(FadeIn(d6), run_time=0.8)
    s.say("$\\chi(\\vec q)$", " 叫結構函數：只由鄰居的位置決定", wait=3.0)
    s.clear_stage()


def part_lattices(s):
    pics = VGroup()
    for k, (kind, name) in enumerate([("sc", "簡單立方 sc"), ("bcc", "體心立方 bcc"), ("fcc", "面心立方 fcc")]):
        org = np.array([-4.3 + 4.3 * k, 2.25, 0])
        pic = L.lattice3d(kind, unit=0.6, origin=org)
        lab = zh(f"{name}：{chi0(kind)} 個鄰居", font_size=20, color=GREY_B).next_to(pic, DOWN, buff=0.12)
        pics.add(VGroup(pic, lab))
    hdr = [zh("格子", font_size=24, color=YELLOW), Tex(r"\chi(\vec q)", font_size=30).set_color(YELLOW),
           Tex(r"\chi(\vec 0)", font_size=30).set_color(YELLOW)]
    rows = [hdr] + [[zh(name, font_size=24), Tex(tex, font_size=32), Tex(chi0(kind), font_size=30).set_color(YELLOW)]
                    for kind, name, tex in LATS]
    tab = D.tex_table(rows, [2.6, 9.0, 1.3], row_h=0.58)
    tab.shift(np.array([-tab.get_width() / 2, pics.get_bottom()[1] - 0.2, 0]))
    rowm = [VGroup(*[D.table_cell(tab, i, j, 3) for j in range(3)]) for i in range(6)]
    assert_in_frame(pics, tab)
    assert_apart(pics, tab)

    s.say("各種格子都一樣做：只差在 ", "$\\chi(\\vec q)$", "（Fig. 2.4）", wait=0.3)
    s.play(FadeIn(tab[1]), FadeIn(rowm[0]), FadeIn(rowm[1]), FadeIn(rowm[2]), run_time=1.2)
    s.say("一維 ", "$2\\cos q$", "；正方格子四個鄰居：", "$2(\\cos q_x+\\cos q_y)$", wait=3.0)
    s.play(LaggedStart(*[FadeIn(p) for p in pics], lag_ratio=0.3), run_time=1.8)
    s.wait(1.0)
    s.play(FadeIn(rowm[3]), run_time=0.8)
    s.say("簡單立方：鄰居在六個面心 (2.12)", wait=3.0)
    s.play(FadeIn(rowm[4]), FadeIn(rowm[5]), run_time=1.0)
    s.say("練習 2.6：bcc（8 個角）與 fcc（12 條稜中點）", wait=3.5)
    s.say("$a$", " 是最近鄰距離；", "$\\chi(\\vec 0)$", "＝鄰居數：2, 4, 6, 8, 12", wait=3.0)
    s.clear_stage()


def part_heatmap(s):
    w, org = 4.6, np.array([-3.4, 0.05, 0])
    hm = L.heatmap(lambda x, y: GF.chi_sq(x, y) / 4, n=32, width=w, origin=org, vmax=1.0)
    frame = L.bz_frame(width=w, origin=org)
    c_lab = Tex(chi0("sq"), font_size=34).move_to(org + np.array([0.3, 0.3, 0]))
    k_val = int(round(GF.chi_sq(math.pi, math.pi)))
    k_labs = VGroup(*[Tex(str(k_val), font_size=30).move_to(org + np.array([sx * (w / 2 - 0.32), sy * (w / 2 - 0.25), 0]))
                      for sx in (1, -1) for sy in (1, -1)])
    tw = total_walks(4)
    h1 = Tex(r"\chi(\vec q)=2(\cos q_x+\cos q_y)", font_size=34)
    h2 = mixed("中心 ", "$\\chi(\\vec 0)=4$", " 最大；四個角 ", "$\\chi=-4$", font_size=24)
    h3 = Tex(r"\vec q=\vec 0:\ \ g(z;\vec 0)=\frac{1}{1-4z}", font_size=36)
    h4 = Tex("=" + "+".join([str(tw[0]), f"{tw[1]}z"] + [f"{tw[k]}z^{k}" for k in range(2, 5)]) + r"+\cdots",
             font_size=34)
    h5 = mixed("係數 ", "$4^N$", "：所有 ", "$N$", " 步走法的總數", font_size=24)
    h6 = Tex(r"g(z;\vec 0)=\sum_{\vec x}G(z;\vec x)=\sum_N z^N\sum_{\vec x}C(N;\vec x)", font_size=28).set_color(GREY_B)
    col = VGroup(h1, h2, h3, h4, h5, h6)
    for m in col:
        fit(m, 5.6)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.36), 0.95, 2.7)
    assert_in_frame(hm, frame, col)
    assert_apart(frame, col)

    s.say("正方格子的 ", "$\\chi$", "：在布里淵區 ", "$[-\\pi,\\pi]^2$", " 上畫出來", wait=0.3)
    s.play(FadeIn(hm), FadeIn(frame), FadeIn(h1), run_time=1.5)
    s.play(FadeIn(c_lab), FadeIn(k_labs), FadeIn(h2), run_time=0.8)
    s.say("橘色為正、藍色為負：中心 4 最大", wait=3.0)
    s.say("四個角 ", "$(\\pm\\pi,\\pm\\pi)$", " 是 ", "$-4$", "（伏筆：W206 會用到）", wait=3.0)
    s.play(FadeIn(h3), run_time=1.0)
    s.say("$\\vec q=\\vec 0$", "：", "$g=1/(1-4z)$", wait=3.0)
    s.play(FadeIn(h4), FadeIn(h5), run_time=1.0)
    s.say("展開：係數 1, 4, 16, 64, 256 ＝ ", "$4^N$", wait=3.0)
    s.play(FadeIn(h6), run_time=0.8)
    s.say("$\\vec q=\\vec 0$", " 把所有終點加起來：", "$N$", " 步走法的總數", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("傅立葉把「和鄰居的關係」變成「乘上 ", "$\\chi(\\vec q)$", "」", wait=3.5)
    s.say("$G$", " 一行代數就解出來：", "$g=1/(1-z\\chi)$", wait=3.5)


class W204StructureFunction(ChuScene):
    def construct(self):
        self.show_header("W204　平移對稱與傅立葉：一行代數解出 G")
        check_math()
        part_question(self)
        part_plane_wave(self)
        part_derivation(self)
        part_lattices(self)
        part_heatmap(self)
        part_close(self)
        self.wait(1.0)
