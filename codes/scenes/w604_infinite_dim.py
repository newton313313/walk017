"""
W604 無限多維：每條漫步都是同一座樓梯（§6.2 開頭，pp. 135–138，(6.19)–(6.20)，(6.53)–(6.54)，Ex. 6.2–6.4）

問：如果空間有 100 維，醉漢走 10 步，會不會走回已經走過的方向？
  1. Ex. 6.2：每一步都換新方向的比例 0.6282（近似 e^{−N²/2d} = 0.6065）→ d ≫ N² 時每條都是同一座樓梯。
  2. N = 3 的樓梯（斜投影）：重心 3/4, 1/2, 1/4；T = (1/16)[[3,2,1],[2,4,2],[1,2,3]]；
     特徵值 0.4268、0.1250、0.0732 vs (N+1)/(π²k²) 0.4053、0.1013、0.0450；特徵向量是正弦。
  3. (6.19)(6.20)：R_n² = (N+1)/(π²n²)，A_∞ = (π⁴/90)/(π²/6)² = 2/5；有限 N：0.5200、0.4150、0.4002。
  4. 和三維比：(6.53) 9 : 2.25 : 1 vs 模擬 11.95 : 2.69 : 1；筆誤框 (6.54)。
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
from utils import lattice_draw as L
from utils import walk_draw as D

F = Fraction
NQ, DQ = 10, 100                 # Ex. 6.2：100 維、10 步
NS = 3                           # 樓梯的步數
NA = (3, 10, 100)
STAIR_COL = TEAL_B
SIM_COL = BLUE_C
POS_COL = "#F2A33A"
NEG_COL = "#4A90D9"
COM_COL = ORANGE
VEC_TEX = (r"(1,\sqrt2,1)", r"(1,0,-1)", r"(1,-\sqrt2,1)")


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def frac_tex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return ("-" if q < 0 else "") + r"\frac{%d}{%d}" % (abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def numbers():
    pts = SH.staircase(NS)
    Tf = SH.tensor_frac([tuple(int(v) for v in p) for p in pts])
    Tn = np.array(Tf, dtype=float)
    st = SH.sim3d_stats()
    ex = SH.eigs_exact(NS)
    return dict(bars={k: float(SH.orth_fraction(k, DQ)) for k in range(1, NQ + 1)},
                orth=f4(SH.orth_fraction(NQ, DQ)), approx=f4(SH.orth_approx(NQ, DQ)),
                pts=pts, com=SH.staircase_com(NS), Tf=Tf, Tn=Tn,
                eig=tuple(f4(v) for v in ex), eig_raw=ex, asym=tuple(f4(SH.eigs_asym(NS, k)) for k in (1, 2, 3)),
                modes=[SH.sine_mode(NS, k) for k in (1, 2, 3)],
                z2=f4(SH.zeta2), z4=f4(SH.zeta4), z2sq=f4(SH.zeta2 ** 2), Ainf=f4(SH.zeta4 / SH.zeta2 ** 2),
                stair={N: SH.staircase_asph(N) for N in NA}, stair_s={N: f4(SH.staircase_asph(N)) for N in NA},
                r653=(SH.eigs_asym(1000, 1) / SH.eigs_asym(1000, 3), SH.eigs_asym(1000, 2) / SH.eigs_asym(1000, 3)),
                r3d=tuple(f"{r:.2f}" for r in st["ratios"]), r3d_raw=st["ratios"],
                ror=tuple(f"{r:.2f}" for r in st["ratio_of_ratios"]))


def check_math():
    nb = numbers()
    prod = F(1)
    for k in range(1, NQ):
        prod *= 1 - F(k, DQ)
    assert prod == SH.orth_fraction(NQ, DQ) and nb["orth"] == "0.6282" and nb["approx"] == "0.6065"
    assert abs(SH.orth_approx(NQ, DQ) - math.exp(-0.5)) < 1e-15
    assert nb["bars"][1] == 1.0 and len(nb["bars"]) == NQ
    assert all(nb["bars"][k + 1] < nb["bars"][k] for k in range(1, NQ))
    assert f4(nb["bars"][NQ]) == nb["orth"]
    assert [tuple(map(int, p)) for p in nb["pts"]] == [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)]
    assert nb["com"] == [F(3, 4), F(1, 2), F(1, 4)]
    assert SH.com_frac([tuple(int(v) for v in p) for p in nb["pts"]]) == nb["com"]
    assert nb["Tf"] == [[F(3, 16), F(2, 16), F(1, 16)], [F(2, 16), F(4, 16), F(2, 16)], [F(1, 16), F(2, 16), F(3, 16)]]
    assert nb["eig"] == ("0.4268", "0.1250", "0.0732") and nb["asym"] == ("0.4053", "0.1013", "0.0450")
    ex = nb["eig_raw"]
    assert np.allclose(ex, [(2 + math.sqrt(2)) / 8, 1 / 8, (2 - math.sqrt(2)) / 8])
    assert np.allclose(SH.principal(nb["Tn"])[0], ex)
    want = [np.array([1, math.sqrt(2), 1]) / 2, np.array([1, 0, -1]) / math.sqrt(2), np.array([1, -math.sqrt(2), 1]) / 2]
    for k in range(3):
        v = nb["modes"][k]
        assert np.allclose(v, want[k]) and np.allclose(nb["Tn"] @ v, ex[k] * v)        # 正弦＝特徵向量
    assert (nb["z2"], nb["z4"], nb["z2sq"], nb["Ainf"]) == ("1.6449", "1.0823", "2.7058", "0.4000")
    assert abs(SH.zeta2 - math.pi ** 2 / 6) < 1e-12 and abs(SH.zeta4 - math.pi ** 4 / 90) < 1e-12
    assert nb["stair"][3] == F(13, 25)
    assert nb["stair_s"] == {3: "0.5200", 10: "0.4150", 100: "0.4002"}
    assert abs(float(SH.staircase_asph(10 ** 6)) - 0.4) < 1e-9
    assert abs(nb["r653"][0] - 9) < 1e-12 and abs(nb["r653"][1] - 2.25) < 1e-12
    assert nb["r3d"] == ("11.95", "2.69") and nb["ror"] == ("13.69", "2.98")


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


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


def cross_mark(p, size=0.13, color=COM_COL, width=4):
    return VGroup(Line(p + size * UL, p + size * DR).set_stroke(color, width),
                  Line(p + size * UR, p + size * DL).set_stroke(color, width))


def signed_bars(values, ux, uy, origin, width=0.6):
    """D.bar_chart 畫絕對值，負的往下翻（bar_chart 只畫正的長條）。"""
    chart = D.bar_chart({k + 1: abs(float(v)) for k, v in enumerate(values)}, ux=ux, uy=uy, origin=origin,
                        width=width, color=POS_COL)
    for k, v in enumerate(values):
        r = D.bar_of(chart, k + 1)
        if v < -1e-12:
            r.shift(abs(float(v)) * uy * DOWN)
            r.set_fill(NEG_COL, 0.75).set_stroke(NEG_COL, 1.5)
    return chart


def sine_group(k, center, ux=0.5, uy=1.3):
    nb = numbers()
    v = nb["modes"][k - 1]
    org = np.array([center[0] - 2 * ux, center[1], 0])
    chart = signed_bars(v, ux, uy, org)
    norm = float(np.linalg.norm(np.sin(k * np.pi * np.arange(1, NS + 1) / (NS + 1))))
    cur = D.curve(lambda a: math.sin(k * math.pi * a / (NS + 1)) / norm, 0, NS + 1, ux=ux, uy=uy, origin=org,
                  color=GREY_A, stroke_width=2.5, dashed=True, num_dashes=24)
    base = Line(org, org + (NS + 1) * ux * RIGHT).set_stroke(GREY_C, 2)
    lab = tx(r"k=%d:\ %s" % (k, VEC_TEX[k - 1]), 26).move_to(np.array([center[0], center[1] + 1.25, 0]))
    return VGroup(base, chart, cur, lab)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    ask = at(zh("100 維空間，醉漢走 10 步：會不會走回已經走過的方向？", font_size=30, color=YELLOW), 0, 2.6)
    prod = eq_label(tx(r"\Big(1-\frac{1}{100}\Big)\Big(1-\frac{2}{100}\Big)\cdots\Big(1-\frac{9}{100}\Big)=%s"
                       % nb["orth"], 32), "Ex. 6.2")
    at(prod, 0, 1.6)
    ux, uy, org = 0.55, 2.4, np.array([-6.1, -2.2, 0])
    chart = D.bar_chart(nb["bars"], ux=ux, uy=uy, origin=org, width=0.7, color=STAIR_COL)
    base = Line(org + 0.4 * ux * RIGHT, org + (NQ + 0.6) * ux * RIGHT).set_stroke(GREY_C, 2)
    one = DashedLine(org + 0.4 * ux * RIGHT + uy * UP, org + (NQ + 0.6) * ux * RIGHT + uy * UP,
                     dash_length=0.08).set_stroke(GREY_B, 1.5)
    one_lab = tx("1", 22, GREY_B).next_to(one, LEFT, buff=0.1)
    ks = VGroup(*[tx(str(k), 20, GREY_B).move_to(org + k * ux * RIGHT + 0.22 * DOWN) for k in range(1, NQ + 1)])
    last = tx(nb["orth"], 24, STAIR_COL).next_to(D.bar_of(chart, NQ), UP, buff=0.1)
    gauss = D.curve(lambda k: math.exp(-k * k / (2 * DQ)), 0.5, NQ + 0.4, ux=ux, uy=uy, origin=org, dashed=True,
                    num_dashes=30)
    kl = mixed("前 ", "$N$", " 步都換新方向的比例", font_size=22, color=GREY_B).next_to(ks, DOWN, buff=0.12)
    plot = VGroup(chart, base, one, one_lab, ks, last, gauss, kl)
    r1 = at(mixed("$\\approx e^{-N^2/2d}=e^{-1/2}=%s$" % nb["approx"], font_size=30, color=D.GAUSS_COLOR), 3.4, 0.35)
    r2 = at(mixed("$d\\gg N^2$", "：幾乎每一步都垂直於", font_size=28), 3.4, -0.6)
    r3 = at(zh("之前所有的步", font_size=28), 3.4, -1.2)
    r4 = at(zh("轉一轉、照鏡子：每條都是同一座樓梯", font_size=26, color=YELLOW), 3.4, -2.1)
    for m in (r1, r2, r3, r4):
        fit(m, 6.2)
    assert_in_frame(ask, prod, plot, r1, r2, r3, r4)
    for a, b in ((ask, prod), (prod, plot), (plot, r1), (plot, r2), (plot, r3), (plot, r4), (r1, r2), (r2, r3),
                 (r3, r4), (prod, r1)):
        assert_apart(a, b)
    s.play(FadeIn(ask), run_time=0.6)
    s.say("問：100 維空間，醉漢走 10 步", wait=1.5)
    s.say("會不會走回已經走過的方向（同一根軸）？", wait=2.6)
    s.play(FadeIn(base), FadeIn(one), FadeIn(one_lab), FadeIn(ks), FadeIn(kl), run_time=0.6)
    s.play(LaggedStart(*[GrowFromEdge(r, DOWN) for r in chart[0]], lag_ratio=0.12), run_time=2.0)
    s.say("第 ", "$k$", " 步碰到舊軸的機率約 ", "$\\frac{k-1}{d}$", "；一步一步乘起來", wait=2.6)
    s.play(Write(prod), FadeIn(last), run_time=1.2)
    s.say("Ex. 6.2：10 步全都換新方向的比例 %s" % nb["orth"], wait=2.9)
    s.play(ShowCreation(gauss), FadeIn(r1), run_time=1.0)
    s.say("近似：", "$e^{-N^2/2d}=%s$" % nb["approx"], "（橘色虛線）", wait=2.6)
    s.play(FadeIn(r2), FadeIn(r3), run_time=0.6)
    s.say("$d\\gg N^2$", "：幾乎每一步都垂直於之前所有的步", wait=2.9)
    s.play(FadeIn(r4), run_time=0.6)
    s.say("轉一轉、照鏡子，每條漫步都是同一個形狀：樓梯", wait=3.0)
    s.clear_stage()


def stair_drawing(unit=2.0, origin=(-5.1, -1.5)):
    nb = numbers()
    org = np.array([origin[0], origin[1], 0.0])
    P = [L.iso_point(*map(float, p), unit=unit, origin=org) for p in nb["pts"]]
    axes = VGroup()
    labs = VGroup()
    for i, (e, d) in enumerate((((1, 0, 0), RIGHT), ((0, 1, 0), UR), ((0, 0, 1), UP))):
        end = L.iso_point(*((1.75 if i == 0 else 1.35) * np.array(e, dtype=float)), unit=unit, origin=org)
        axes.add(DashedLine(org, end, dash_length=0.08).set_stroke(GREY_C, 1.5))
        labs.add(tx(r"x_%d" % (i + 1), 24, GREY_B).next_to(end, d, buff=0.08))
    steps = VGroup(*[Arrow(P[k], P[k + 1], buff=0.04, thickness=4).set_color(YELLOW) for k in range(NS)])
    slabs = VGroup(*[tx(r"\eta_%d" % (k + 1), 24, YELLOW).next_to(steps[k], d, buff=0.1)
                     for k, d in zip(range(NS), (DOWN, RIGHT, RIGHT))])
    dots = VGroup(*[Dot(p, radius=0.06).set_fill(WHITE, 1) for p in P])
    c = L.iso_point(*[float(x) for x in nb["com"]], unit=unit, origin=org)
    for m in slabs:
        assert_apart(m, cross_mark(c))
        for a in labs:
            assert_apart(m, a)
    return dict(axes=axes, labs=labs, steps=steps, slabs=slabs, dots=dots, cross=cross_mark(c),
                all=VGroup(axes, labs, steps, slabs, dots, cross_mark(c)))


def part_stair(s):
    nb = numbers()
    sd = stair_drawing()
    c1 = eq_label(tx(r"\langle x_n\rangle=\Big(\frac{3}{4},\ \frac{1}{2},\ \frac{1}{4}\Big)", 32), "Ex. 6.3(a)")
    at(c1, 2.7, 2.5)
    c1[0].set_color(COM_COL)
    Tt = at(tx(r"\mathbf T=\frac{1}{16}\begin{pmatrix}3&2&1\\2&4&2\\1&2&3\end{pmatrix}", 34), 2.7, 1.15)
    rows = [[VectorizedPoint(), "k=1", "k=2", "k=3"],
            [zh("精確", font_size=26)] + list(nb["eig"]),
            [r"\frac{N+1}{\pi^2k^2}"] + list(nb["asym"])]
    tab = D.tex_table(rows, [1.7, 1.5, 1.5, 1.5], row_h=0.66, font_size=30)
    at(tab, 2.7, -0.95)
    for j in range(1, 4):
        D.table_cell(tab, 2, j, 4).set_color(D.GAUSS_COLOR)
    assert_in_frame(sd["all"], c1, Tt, tab)
    for a, b in ((sd["all"], c1), (sd["all"], Tt), (sd["all"], tab), (c1, Tt), (Tt, tab)):
        assert_apart(a, b)
    assert_table(tab, 4)
    s.play(FadeIn(sd["axes"]), FadeIn(sd["labs"]), FadeIn(sd["dots"][0]), run_time=0.6)
    s.say("$N=3$", " 的樓梯：第 ", "$n$", " 步沿第 ", "$n$", " 根軸走 1", wait=0.3)
    for k in range(NS):
        s.play(GrowArrow(sd["steps"][k]), FadeIn(sd["slabs"][k]), FadeIn(sd["dots"][k + 1]), run_time=0.6)
    s.say("三根互相垂直的箭頭（斜投影）", wait=2.0)
    s.play(ShowCreation(sd["all"][5]), FadeIn(c1), run_time=0.8)
    s.say("重心：", "$\\tfrac34,\\ \\tfrac12,\\ \\tfrac14$", "（Ex. 6.3(a)）", wait=2.6)
    s.play(FadeIn(Tt, shift=0.2 * DOWN), run_time=0.8)
    s.say("迴轉半徑張量：", "$3\\times3$", "，乘 ", "$\\tfrac{1}{16}$", wait=2.6)
    s.play(FadeIn(tab), run_time=0.8)
    s.say("特徵值 %s、%s、%s" % nb["eig"], wait=2.4)
    s.say("大 ", "$N$", " 公式 ", "$\\tfrac{N+1}{\\pi^2k^2}$", "：%s、%s、%s" % nb["asym"], wait=2.9)
    s.say("才 3 步，最大的已經差不多", wait=2.2)
    s.play(FadeOut(Tt), FadeOut(tab), FadeOut(c1), run_time=0.5)
    groups = VGroup(*[sine_group(k, (x, -0.4)) for k, x in zip((1, 2, 3), (0.7, 3.1, 5.5))])
    head = at(mixed("特徵向量：", "$\\sin\\frac{k\\pi\\alpha}{N+1}$", "（", "$\\alpha=1,2,3$", "）", font_size=28), 3.1,
              2.45)
    assert_in_frame(groups, head)
    for a, b in ((sd["all"], groups), (head, groups), (head, sd["all"])):
        assert_apart(a, b)
    for i in range(2):
        assert_apart(groups[i], groups[i + 1])
    s.play(FadeIn(head), run_time=0.6)
    s.play(LaggedStart(*[FadeIn(g) for g in groups], lag_ratio=0.3), run_time=1.5)
    s.say("特徵向量是正弦的形狀：", "$k$", " 越大，振動越多", wait=2.9)
    s.say("和 W605 的 Kramers 矩陣同一組正弦", wait=2.4)
    s.clear_stage()


def part_limit(s):
    nb = numbers()
    e1 = eq_label(tx(r"R_n^2=\frac{N+1}{\pi^2 n^2}", 36), "(6.19)")
    at(e1, -3.2, 2.35)
    e0 = at(tx(r"A_\infty=\frac{\mathrm{Tr}\,\mathbf T^2}{(\mathrm{Tr}\,\mathbf T)^2}", 34), 3.2, 2.35)
    e2 = eq_label(tx(r"A_\infty=\frac{\sum_n n^{-4}}{\big(\sum_n n^{-2}\big)^2}=\frac{\pi^4/90}{(\pi^2/6)^2}"
                     r"=\frac{%s}{%s}=\frac{2}{5}" % (nb["z4"], nb["z2sq"]), 36), "(6.20)")
    at(e2, 0, 0.95)
    e3 = at(tx(r"\sum_n n^{-2}=\frac{\pi^2}{6}=%s" % nb["z2"], 28, GREY_A), 0, -0.25)
    rows = [["N"] + [str(N) for N in NA] + [r"\infty"],
            [r"A"] + [nb["stair_s"][N] for N in NA] + [r"0.4"]]
    tab = D.tex_table(rows, [1.2, 1.6, 1.6, 1.6, 1.2], row_h=0.7, font_size=32)
    at(tab, -1.1, -1.75)
    note = at(mixed("$N=3$", "：", "$\\frac{13}{25}$", font_size=28, color=STAIR_COL), 4.2, -1.75)
    assert_in_frame(e1, e0, e2, e3, tab, note)
    for a, b in ((e1, e0), (e1, e2), (e0, e2), (e2, e3), (e3, tab), (tab, note)):
        assert_apart(a, b)
    assert_table(tab, 5)
    s.play(Write(e1), run_time=1.0)
    s.say("一般 ", "$N$", "：樓梯的 ", "$\\mathbf T$", " 的特徵值，大 ", "$N$", " 時是 ", "$\\tfrac{N+1}{\\pi^2n^2}$", wait=2.9)
    s.say("第 ", "$n$", " 個主迴轉半徑的平方 ", "$\\propto 1/n^2$", wait=2.4)
    s.play(FadeIn(e0), run_time=0.6)
    s.say("$d\\to\\infty$", "：非球度只剩 ", "$\\mathrm{Tr}\\,\\mathbf T^2/(\\mathrm{Tr}\\,\\mathbf T)^2$", wait=2.9)
    s.play(Write(e2), run_time=1.5)
    s.say("$\\sum n^{-4}=%s$" % nb["z4"], "、", "$\\big(\\sum n^{-2}\\big)^2=%s$" % nb["z2sq"], wait=2.6)
    s.play(FadeIn(e3), run_time=0.6)
    s.say("比值正好 ", "$\\tfrac25$", "：和 W603 表上 ", "$d=\\infty$", " 那一格一樣", wait=2.9)
    s.play(FadeIn(tab), FadeIn(note), run_time=0.8)
    s.say("有限步的樓梯：", "$N=3$", " 是 ", "$\\tfrac{13}{25}=%s$" % nb["stair_s"][3], wait=2.6)
    s.say("$N=10$", "：%s；" % nb["stair_s"][10], "$N=100$", "：%s → 0.4" % nb["stair_s"][100], wait=2.9)
    s.clear_stage()


def part_compare(s):
    nb = numbers()
    r1, r2 = nb["r653"]
    s1, s2 = nb["r3d_raw"]
    ux, uy, org = 1.0, 0.24, np.array([-5.6, -0.85, 0])
    stair_vals = {1: r1, 3: r2, 5: 1.0}
    sim_vals = {1.75: s1, 3.75: s2, 5.75: 1.0}
    c1 = D.bar_chart(stair_vals, ux=ux, uy=uy, origin=org, width=0.65, color=STAIR_COL)
    c2 = D.bar_chart(sim_vals, ux=ux, uy=uy, origin=org, width=0.65, color=SIM_COL)
    lab1 = VGroup(*[tx(t, 24, STAIR_COL).next_to(D.bar_of(c1, x), UP, buff=0.08)
                    for t, x in (("9", 1), ("2.25", 3), ("1", 5))])
    lab2 = VGroup(*[tx(t, 24, SIM_COL).next_to(D.bar_of(c2, x), UP, buff=0.08)
                    for t, x in ((nb["r3d"][0], 1.75), (nb["r3d"][1], 3.75), ("1", 5.75))])
    base = Line(org + 0.4 * RIGHT, org + 6.4 * RIGHT).set_stroke(GREY_C, 2)
    gl = VGroup(*[tx(t, 26, GREY_A).move_to(org + (x + 0.375) * ux * RIGHT + 0.3 * DOWN)
                  for t, x in ((r"R_1^2", 1), (r"R_2^2", 3), (r"R_3^2", 5))])
    chart = VGroup(c1, c2, lab1, lab2, base, gl)
    l1 = mixed("(6.53) ", "$d\\to\\infty$", " 樓梯：", "$9:2.25:1$", font_size=26, color=STAIR_COL)
    l2 = mixed("三維格子模擬：", "$%s:%s:1$" % nb["r3d"], font_size=26, color=SIM_COL)
    leg = VGroup(l1, l2).arrange(DOWN, buff=0.35, aligned_edge=LEFT).move_to(np.array([3.6, 1.2, 0]))
    box = D.erratum_box("p. 146，(6.54)",
                        r"\langle R_1^2\rangle:\langle R_2^2\rangle:\langle R_3^2\rangle=13.76:3.03:1",
                        r"\langle R_1^2\rangle:\langle R_2^2\rangle:\langle R_3^2\rangle\approx11.9:2.7:1",
                        width=12.4, font_size=28)
    at(box, 0, -2.2)
    assert_in_frame(chart, leg, box)
    for a, b in ((chart, leg), (chart, box), (leg, box), (l1, l2)):
        assert_apart(a, b)
    for a in lab1:
        for b in lab2:
            assert_apart(a, b)
    s.play(FadeIn(base), FadeIn(gl), run_time=0.5)
    s.play(LaggedStart(*[GrowFromEdge(r, DOWN) for r in c1[0]], lag_ratio=0.2), FadeIn(lab1), FadeIn(l1),
           run_time=1.2)
    s.say("和三維比：(6.53) 樓梯的比 ", "$9:2.25:1$", wait=2.6)
    s.play(LaggedStart(*[GrowFromEdge(r, DOWN) for r in c2[0]], lag_ratio=0.2), FadeIn(lab2), FadeIn(l2),
           run_time=1.2)
    s.say("三維格子模擬：", "$%s:%s:1$" % nb["r3d"], "，差不多", wait=2.9)
    s.say("三維也是「一長、一中、一短」的樓梯形", wait=2.4)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (6.54) 的 13.76 : 3.03 : 1 不是平均的比", wait=2.9)
    s.say("書上的數字接近「比值的平均」", "$\\langle R_1^2/R_3^2\\rangle:\\langle R_2^2/R_3^2\\rangle$", wait=2.9)
    s.say("模擬的比值平均：", "$%s:%s$" % nb["ror"], wait=2.6)
    s.clear_stage()


def part_close(s):
    s.say("維度夠高時，每條漫步都是同一座樓梯，非球度正好 ", "$\\tfrac25$", wait=3.4)


class W604InfiniteDim(ChuScene):
    def construct(self):
        self.show_header("W604　無限多維：每條漫步都是同一座樓梯")
        check_math()
        part_question(self)
        part_stair(self)
        part_limit(self)
        part_compare(self)
        part_close(self)
        self.wait(1.0)
