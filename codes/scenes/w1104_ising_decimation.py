"""
W1104 一維 Ising 的抽格（§11.2 例題，pp. 269–272，(11.50)–(11.63)，Ex. 11.7–11.9）

問：Ising 鏈的奇數自旋加總掉，剩下的自旋之間是什麼耦合？
  1. 9 個 ±1 自旋，奇數格點淡出；(11.51)(11.52) Σ_{S3} e^{K S3(S2+S4)} = 2cosh K(S2+S4)：兩種情況（K = 1）；筆誤框 1（(11.53)）。
  2. Ex. 11.7 (11.54)(11.55)：A = 2√cosh 2K、K' = ½ ln cosh 2K；(11.62) 迭代表 K_n → 0；Ex. 11.8 一行（外場）。
  3. (11.57)–(11.63)：4、8 個自旋的環：暴力列舉 ＝ 轉移矩陣 ＝ 抽格連乘；筆誤框 2（(11.63) 的指數）。
  4. Ex. 11.9(a)：x = e^{−2K}、x' = 2x/(1+x²) 的蛛網；不動點 x = 0（不穩定）、x = 1（穩定）。
  5. Ex. 11.9(b)(c)：ξ ∝ e^{2K}（表）、ln Z/(NK) → 1；收束。
畫面上的數字都由 utils.walk_rg 算出，check_math() 核對。
模組層級的畫圖函式：ising_chain（自旋鏈）、x_axes／x_curve／cobweb_pts（蛛網）。
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
from utils import walk_rg as RG
from utils import walk_draw as D

UP_COLOR = ORANGE
DOWN_COLOR = BLUE_C
NEW_BOND = YELLOW
OK_COLOR = "#7CE38B"
SPINS = (1, 1, -1, 1, 1, -1, -1, 1, -1)       # S_1..S_9（示意）
K0 = 1.0
H0 = 0.1
X_CURVE_N = 40


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


@lru_cache(maxsize=None)
def numbers():
    x0 = RG.x_of_K(K0)
    return {"b53": RG.book1153(K0),
            "anti": RG.ising_A(K0) * math.exp(-RG.ising_Kp(K0)),          # S2 = −S4：A e^{−K'} = 2
            "A": RG.ising_A(K0), "Kp": RG.ising_Kp(K0),
            "Kn": RG.ising_iterate(K0, 5),
            "ex118": RG.ex118(K0, H0),
            "brute": {N: RG.ring_Z_brute(N, K0) for N in (4, 8)},
            "transfer": {N: RG.ring_Z_transfer(N, K0) for N in (4, 8)},
            "decim": {N: RG.ring_Z_decimation(N, K0) for N in (4, 8)},
            "book": {N: RG.ring_Z_decimation(N, K0, book=True) for N in (4, 8)},
            "xs": RG.x_iterate(x0, 5),
            "xi": [RG.xi_exact(K) for K in RG.XI_K],
            "half": [math.exp(2 * K) / 2 for K in RG.XI_K],
            "ratio": [RG.xi_ratio(K) for K in RG.XI_K],
            "free": [RG.free_energy_ratio(K) for K in RG.FREE_K]}


def check_math():
    nb = numbers()
    # §1 (11.51)–(11.53)
    assert [f4(v) for v in nb["b53"]] == ["3.7622", "7.5244", "7.5244"]
    assert f4(nb["anti"]) == "2.0000" and abs(2 * math.cosh(2 * K0) - nb["b53"][2]) < 1e-12
    # §2 Ex. 11.7、(11.62)、Ex. 11.8
    assert f4(nb["A"]) == "3.8793" and f4(nb["Kp"]) == "0.6625" and RG.ising_pair_check(K0) < 1e-12
    assert [f4(v) for v in nb["Kn"]] == ["1.0000", "0.6625", "0.3501", "0.1137", "0.0128", "0.0002"]
    assert all(a > b > 0 for a, b in zip(nb["Kn"], nb["Kn"][1:]))
    assert [f4(v) for v in nb["ex118"]] == ["3.8897", "0.6602", "0.0482", "0.1964"]
    assert RG.ex118_check(K0, H0) < 1e-12 and abs(nb["ex118"][3] - (H0 + 2 * nb["ex118"][2])) < 1e-15
    # §3 (11.57)–(11.63)
    for N, want in ((4, "121.2329"), (8, "9160.4387")):
        assert f4(nb["brute"][N]) == f4(nb["transfer"][N]) == f4(nb["decim"][N]) == want
    assert [f4(nb["book"][N]) for N in (4, 8)] == ["18.5498", "143.2705"]
    # §4 Ex. 11.9(a)
    assert [f4(v) for v in nb["xs"]] == ["0.1353", "0.2658", "0.4965", "0.7966", "0.9747", "0.9997"]
    assert RG.x_fixed_points() == (0.0, 1.0) and f4(RG.x_slope(0.0)) == "2.0000" and f4(RG.x_slope(1.0)) == "0.0000"
    assert all(abs(RG.x_map(x) - math.exp(-2 * RG.ising_Kp(-math.log(x) / 2))) < 1e-12 for x in nb["xs"][:4])
    # §5 Ex. 11.9(b)(c)
    assert RG.XI_K == (1.0, 2.0, 3.0) and RG.FREE_K == (1.0, 2.0, 5.0)
    assert [f4(v) for v in nb["xi"]] == ["3.6719", "27.2960", "201.7140"]
    assert [f4(v) for v in nb["half"]] == ["3.6945", "27.2991", "201.7144"]
    assert [f4(v) for v in nb["ratio"]] == ["0.4969", "0.4999", "0.5000"]
    assert all(max(abs(e) for e in RG.tanh_rg_check(K)) < 1e-9 for K in RG.XI_K)
    assert [f4(v) for v in nb["free"]] == ["1.1269", "1.0091", "1.0000"]


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框（只縮不放），再置中。"""
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
def site_x(i, dx=1.3):
    return (i - 4) * dx


def ising_chain(spins, y=2.0, dx=1.3):
    """9 個 ±1 自旋：VGroup(鍵們, 格點們, 箭頭們, 標籤們)；第 i 個（0 起算）是 S_{i+1}。"""
    pts = [np.array([site_x(i, dx), y, 0]) for i in range(len(spins))]
    bonds = VGroup(*[Line(a, b).set_stroke(GREY_C, 3) for a, b in zip(pts, pts[1:])])
    dots = VGroup(*[Dot(p, radius=0.07).set_color(WHITE) for p in pts])
    arrows = VGroup(*[Arrow(p - 0.42 * sg * UP, p + 0.42 * sg * UP, buff=0, thickness=4)
                      .set_color(UP_COLOR if sg > 0 else DOWN_COLOR) for p, sg in zip(pts, spins)])
    labels = VGroup(*[tx("S_{%d}" % (i + 1), 26, GREY_A).move_to(p + 0.72 * DOWN) for i, p in enumerate(pts)])
    return VGroup(bonds, dots, arrows, labels)


def x_axes():
    return D.plot_axes((0, 1.08, 0.25), (0, 1.08, 0.25), width=4.6, height=4.6, x_ticks=(0, 0.5, 1),
                       y_ticks=(0.5, 1), x_tick_labels=["0", "0.5", "1"], y_tick_labels=["0.5", "1"],
                       x_label="x", y_label="x'")


def x_curve(ax):
    pts = [ax.c2p(x, RG.x_map(x)) for x in np.linspace(0, 1, X_CURVE_N)]
    return VMobject().set_points_smoothly(pts).set_stroke(TEAL, 4)


def cobweb_pts(xs):
    """(x0, 0) → (x0, x1) → (x1, x1) → (x1, x2) → …"""
    out = [(xs[0], 0.0)]
    for a, b in zip(xs, xs[1:]):
        out += [(a, b), (b, b)]
    return out


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    chain = ising_chain(SPINS)
    bonds, dots, arrows, labels = chain
    zeq = tx(r"Z=\sum_{\{S_i=\pm1\}}e^{K\sum_iS_iS_{i+1}}", 36).move_to(0.45 * UP)
    assert_in_frame(chain, zeq)
    assert_apart(chain, zeq)
    s.say("問：Ising 鏈的奇數自旋加總掉，", wait=1.4)
    s.say("剩下的自旋之間是什麼耦合？", wait=1.4)
    s.play(ShowCreation(bonds), FadeIn(dots), run_time=0.8)
    s.play(*[GrowArrow(a) for a in arrows], FadeIn(labels), run_time=1.0)
    s.play(Write(zeq), run_time=1.0)
    s.say("$S_i=\\pm1$", "，", "$K=J/k_BT$", "：只有最近鄰的耦合", wait=2.0)
    odd = [0, 2, 4, 6, 8]
    s.play(*[FadeOut(arrows[i]) for i in odd], *[FadeOut(dots[i]) for i in odd],
           *[labels[i].animate.set_opacity(0.25) for i in odd], *[FadeOut(b) for b in bonds], run_time=1.2)
    new = VGroup(*[Line(dots[i].get_center(), dots[i + 2].get_center()).set_stroke(NEW_BOND, 5)
                   for i in (1, 3, 5)])
    kp = VGroup(*[tx("K'", 26, NEW_BOND).move_to(np.array([site_x(i + 1), 2.0 + 0.35, 0])) for i in (1, 3, 5)])
    s.play(ShowCreation(new), FadeIn(kp), run_time=1.0)
    s.say("$S_3$", " 只連著 ", "$S_2$", " 和 ", "$S_4$", "：先對它求和", wait=1.8)
    eq = tx(r"\sum_{S_3=\pm1}e^{KS_3(S_2+S_4)}=2\cosh K(S_2+S_4)", 36).move_to(0.45 * UP)
    tag = zh("(11.51)(11.52)", font_size=22, color=GREY_B).next_to(eq, LEFT, buff=0.3)
    cases = VGroup(mixed("$S_2=S_4$", "：", "$2\\cosh 2K=%s$" % f4(nb["b53"][2]), font_size=26),
                   mixed("$S_2=-S_4$", "：", "$2\\cosh 0=%s$" % f4(nb["anti"]).rstrip("0").rstrip("."),
                         font_size=26)).arrange(RIGHT, buff=0.8)
    cases = framed(VGroup(zh("K = 1：", font_size=24, color=GREY_A), cases).arrange(RIGHT, buff=0.3), GREY_C, 0.2)
    box = D.erratum_box("p. 270，(11.53)", r"\cosh K(S_2+S_4)=2\cosh 2K\ (S_2=S_4)",
                        r"2\cosh K(S_2+S_4)=2\cosh 2K\ (S_2=S_4)", width=13.2, font_size=28)
    low = VGroup(cases, box).arrange(DOWN, buff=0.4)
    place_in(low, -6.7, 6.7, -2.9, -0.3)
    assert_in_frame(eq, tag, low)
    assert_apart(eq, low)
    assert_apart(cases, box)
    assert_apart(tag, eq)
    s.play(FadeOut(zeq), FadeIn(eq), FadeIn(tag), run_time=1.0)
    s.play(FadeIn(cases), run_time=0.8)
    s.say("兩邊同號：", "$2\\cosh2=%s$" % f4(nb["b53"][2]), "；反號：", "$2$", wait=2.2)
    s.say("只看 ", "$S_2S_4$", " 同號或反號——又是一個 Ising 耦合", wait=2.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (11.53) 漏了 2：", "$K=1$", " 時 ", "$\\cosh2=%s$" % f4(nb["b53"][0]),
          "，不是 ", "$%s$" % f4(nb["b53"][1]), wait=2.6)
    s.clear_stage()


def part_iterate(s):
    nb = numbers()
    A, Kp = nb["A"], nb["Kp"]
    eqs = VGroup(tx(r"Ae^{K'}=2\cosh 2K", 34), tx(r"Ae^{-K'}=2", 34)).arrange(RIGHT, buff=1.0)
    tag = zh("Ex. 11.7　(11.54)(11.55)", font_size=22, color=GREY_B)
    res = VGroup(tx(r"A=2\sqrt{\cosh 2K}=%s" % f4(A), 34), tx(r"K'=\tfrac12\ln\cosh 2K=%s" % f4(Kp), 34)
                 ).arrange(RIGHT, buff=0.9)
    for m in res:
        m.set_color(YELLOW)
    ttl = mixed("(11.62) 一直抽格，", "$K=1$", " 出發：", font_size=24, color=GREY_A)
    rows = [["n"] + [str(n) for n in range(1, 7)], ["K_n"] + [f4(v) for v in nb["Kn"]]]
    tab = D.tex_table(rows, [1.0] + [1.45] * 6, row_h=0.55, font_size=28)
    e118 = nb["ex118"]
    ex = framed(VGroup(zh("Ex. 11.8　加外場 H = 0.1（K = 1）：", font_size=24, color=GREY_A),
                       tx(r"A=%s,\quad K_1'=%s,\quad K_2'=%s,\quad H'=H+2K_2'=%s"
                          % tuple(f4(v) for v in e118), 30)).arrange(DOWN, buff=0.2), GREY_C, 0.2)
    col = VGroup(tag, eqs, res, ttl, tab, ex).arrange(DOWN, buff=0.32)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(tag, eqs, res, ttl, tab, ex)
    assert_table(tab, 7)
    s.play(FadeIn(tag), Write(eqs), run_time=1.2)
    s.say("要求：", "$\\sum_{S_3}e^{KS_3(S_2+S_4)}=Ae^{K'S_2S_4}$", wait=2.2)
    s.say("同號、反號兩種情況，正好解出兩個未知數", wait=2.0)
    s.play(Write(res), run_time=1.2)
    s.say("$K=1$", "：", "$A=%s$" % f4(A), "、", "$K'=%s$" % f4(Kp), "——耦合變弱了", wait=2.2)
    s.play(FadeIn(ttl), FadeIn(tab[1]), FadeIn(tab[0][:2]), FadeIn(tab[0][7:9]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(VGroup(tab[0][j], tab[0][7 + j])) for j in range(2, 7)], lag_ratio=0.3),
           run_time=1.6)
    s.say("$K_n$", " 一路掉到 0：越看越像高溫", wait=2.2)
    s.say("每抽一次，長度單位變兩倍，耦合卻越來越弱", wait=2.2)
    s.play(FadeIn(ex), run_time=0.8)
    s.say("有外場時：鍵 ", "$K_1'$", "，外場多了兩邊各 ", "$K_2'$", wait=2.2)
    s.say("$H'=%s+2\\times%s=%s$" % (H0, f4(e118[2]), f4(e118[3])), "：外場會長大", wait=2.2)
    s.clear_stage()


def part_ring(s):
    nb = numbers()
    f1 = tx(r"Z_N=\prod_{j=1}^{k}A(K_j)^{N/2^{j}}\cdot 2e^{K_{k+1}}", 34)
    f1n = mixed("(11.63)　", "$N=2^k$", " 一直抽到剩 1 個自旋", font_size=22, color=GREY_B)
    f2 = tx(r"Z_N=(2\cosh K)^N+(2\sinh K)^N", 34)
    f2n = zh("轉移矩陣（第 9 章）", font_size=22, color=GREY_B)
    top = VGroup(VGroup(f1, f1n).arrange(DOWN, buff=0.15), VGroup(f2, f2n).arrange(DOWN, buff=0.15)
                 ).arrange(RIGHT, buff=0.8)
    head = ["N", mixed("暴力 ", "$2^N$", " 項", font_size=24, color=YELLOW), zh("轉移矩陣", font_size=24, color=YELLOW),
            zh("抽格連乘", font_size=24, color=YELLOW)]
    rows = [head] + [[str(N), f4(nb["brute"][N]), f4(nb["transfer"][N]), f4(nb["decim"][N])] for N in (4, 8)]
    tab = D.tex_table(rows, [1.0, 2.6, 2.6, 2.6], row_h=0.6, font_size=30)
    box = D.erratum_box("p. 271，(11.63)", r"\cdots A(K_j)^{N/2^{j+1}}\cdots", r"\cdots A(K_j)^{N/2^{j}}\cdots",
                        width=6.0, font_size=28)
    col = VGroup(top, tab, box).arrange(DOWN, buff=0.4)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(top, tab, box)
    assert_apart(top[0], top[1])
    assert_table(tab, 4)
    s.say("配分函數變成連乘：每抽一次，乘上 ", "$A(K_j)$", " 的幾次方", wait=2.0)
    s.play(FadeIn(top[0]), run_time=1.0)
    s.say("第 ", "$j$", " 次抽掉 ", "$N/2^j$", " 個自旋，每個給一個 ", "$A$", wait=2.0)
    s.play(FadeIn(top[1]), run_time=0.8)
    s.say("$2^N$", " 項的和，變成 ", "$k=\\log_2N$", " 個數的連乘", wait=2.0)
    s.play(FadeIn(tab[1]), FadeIn(tab[0][:4]), run_time=0.8)
    s.play(FadeIn(tab[0][4:8]), run_time=0.8)
    s.say("$K=1$", "、4 個自旋的環：三種算法都是 ", "$%s$" % f4(nb["decim"][4]), wait=2.2)
    s.say("和第 9 章 W902 的 4 自旋環同一個數", wait=1.8)
    s.play(FadeIn(tab[0][8:12]), run_time=0.8)
    s.say("8 個自旋：", "$%s$" % f4(nb["decim"][8]), wait=1.8)
    s.play(FadeIn(box), run_time=0.8)
    s.say("照書上的指數：", "$%s$" % f4(nb["book"][4]), "、", "$%s$" % f4(nb["book"][8]), "——少了一半", wait=2.4)
    s.clear_stage()


def part_cobweb(s):
    nb = numbers()
    xs = nb["xs"]
    ax = x_axes()
    ax.move_to(np.array([-3.2, 0.1, 0]))
    curve = x_curve(ax)
    diag = DashedLine(ax.c2p(0, 0), ax.c2p(1, 1), dash_length=0.08).set_stroke(GREY_A, 2)
    web = cobweb_pts(xs)
    segs = VGroup(*[Line(ax.c2p(*a), ax.c2p(*b)).set_stroke(YELLOW, 3) for a, b in zip(web, web[1:])])
    form = VGroup(tx(r"x=e^{-2K}", 34), tx(r"x'=\frac{2x}{1+x^2}", 38)).arrange(RIGHT, buff=0.6)
    vals = VGroup(*[tx(r"x_{%d}=%s" % (n, f4(v)), 30) for n, v in enumerate(xs)])
    vals.arrange_in_grid(n_rows=3, n_cols=2, h_buff=0.7, v_buff=0.22, fill_rows_first=False)
    fp0 = mixed("$x=0$", "：", "$T=0$", "，斜率 2，不穩定", font_size=26)
    fp1 = mixed("$x=1$", "：", "$T=\\infty$", "，斜率 0，穩定", font_size=26)
    fps = VGroup(fp0, fp1).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    right = VGroup(form, vals, fps).arrange(DOWN, buff=0.5)
    place_in(right, 0.2, 6.7, -2.8, 3.0)
    d0 = Dot(ax.c2p(0, 0), radius=0.09).set_color(RED)
    d1 = Dot(ax.c2p(1, 1), radius=0.09).set_color(OK_COLOR)
    assert_in_frame(ax, right, curve)
    assert_apart(ax, right)
    assert_column(form, vals, fps)
    s.say("Ex. 11.9：換成 ", "$x=e^{-2K}$", "，抽格變成 ", "$x'=2x/(1+x^2)$", wait=2.2)
    s.play(ShowCreation(ax), FadeIn(form), run_time=1.0)
    s.play(ShowCreation(curve), ShowCreation(diag), run_time=1.2)
    s.say("$x'=x$", " 的地方是不動點：", "$x=0$", " 和 ", "$x=1$", wait=2.0)
    s.say("$K=1$", " 出發：", "$x_0=e^{-2}=%s$" % f4(xs[0]), wait=1.6)
    s.play(FadeIn(vals[0]), run_time=0.4)
    for k in range(len(segs) // 2):
        s.play(ShowCreation(segs[2 * k]), run_time=0.45)
        s.play(ShowCreation(segs[2 * k + 1]), FadeIn(vals[k + 1]), run_time=0.45)
    s.say("一步步走向 ", "$x=1$", "：", "$%s$" % f4(xs[-1]), wait=2.0)
    s.play(FadeIn(d0), FadeIn(d1), FadeIn(fps), run_time=0.8)
    s.say("$x=0$", " 是 ", "$T=0$", "：斜率 2，一碰就被推開", wait=2.2)
    s.say("$x=1$", " 是 ", "$T=\\infty$", "：斜率 0，什麼都流進來", wait=2.2)
    s.say("只要 ", "$T>0$", "，最後都流到無限高溫", wait=2.2)
    s.clear_stage()


def part_xi(s):
    nb = numbers()
    top = VGroup(mixed("$x$", " 小時 ", "$x'\\approx2x$", "、", "$\\xi'=\\xi/2$", font_size=28),
                 mixed("$\\Rightarrow$", "$\\xi\\propto1/x=e^{2J/k_BT}$", font_size=28)).arrange(RIGHT, buff=0.4)
    rows = [["K", r"\xi", r"e^{2K}/2", r"\xi/e^{2K}"]]
    rows += [["%d" % K, f4(a), f4(b), f4(c)] for K, a, b, c in zip(RG.XI_K, nb["xi"], nb["half"], nb["ratio"])]
    tab = D.tex_table(rows, [1.0, 2.2, 2.2, 2.0], row_h=0.52, font_size=28)
    ttl = mixed("自由能：", "$\\ln Z/(NK)$", font_size=26, color=GREY_A)
    rows2 = [["K"] + ["%d" % K for K in RG.FREE_K], [r"\ln Z/(NK)"] + [f4(v) for v in nb["free"]]]
    tab2 = D.tex_table(rows2, [2.2, 1.6, 1.6, 1.6], row_h=0.52, font_size=28)
    low = VGroup(ttl, tab2).arrange(DOWN, buff=0.2)
    col = VGroup(top, tab, low).arrange(DOWN, buff=0.4)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(top, tab, low)
    assert_table(tab, 4)
    assert_table(tab2, 4)
    s.say("關聯長度呢？", "$x$", " 小時 ", "$x'\\approx2x$", "，長度單位變成兩倍", wait=2.2)
    s.play(FadeIn(top), run_time=1.0)
    s.say("所以 ", "$\\xi'=\\xi/2$", "，", "$\\xi\\propto1/x=e^{2K}$", wait=2.0)
    s.play(FadeIn(tab), run_time=1.0)
    s.say("精確的 ", "$\\xi=-1/\\ln\\tanh K$", "：比值趨向 ", "$\\tfrac12$", wait=2.2)
    s.say("其實 ", "$\\tanh K'=\\tanh^2K$", "，", "$\\xi$", " 每次剛好減半", wait=2.2)
    s.play(FadeIn(low), run_time=1.0)
    s.say("低溫時 ", "$\\ln Z/(NK)\\to1$", "：", "$F\\to-NJ$", "，全部排整齊", wait=2.2)
    s.say("一維 Ising 的流只有兩個不動點：", "$T=0$", " 不穩定", wait=2.2)
    s.say("$T=\\infty$", " 穩定——一維沒有有限溫度的相變", wait=2.6)


class W1104IsingDecimation(ChuScene):
    def construct(self):
        self.show_header("W1104　一維 Ising 的抽格")
        check_math()
        part_question(self)
        part_iterate(self)
        part_ring(self)
        part_cobweb(self)
        part_xi(self)
        self.wait(1.0)
