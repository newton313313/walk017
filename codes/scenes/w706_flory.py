"""
W706 Flory 的論證＋第 7 章回顧（§7.2.7，pp. 191–192，(7.61)–(7.64)）

問：m 個單體的高分子（自迴避漫步），大小 R 怎麼隨 m 變大？
  1. (7.61) 體積 𝒩 = αR^d、u = f/α：C_MF = (2d)^m/𝒩 · e^{−um²/R^d}（排斥）；(7.62) 熵 e^{−R²/m}。拔河圖。
  2. (7.63)(7.64)：d = 2、m = 1000、u = 1 的兩項與和（R ∈ [50, 600]）；R* = 177.83，兩項都是 −31.62、和 −63.25；R ∝ m^{3/(d+2)}。
  3. 表：Flory 1, 3/4, 3/5, 1/2；數值最大化的斜率 1.0000 … 0.5000；已知值 1、3/4、0.588、1/2；
     二維精確列舉的局部斜率 1.4306 … 1.4556 → 3/2。
  4. 第 7 章回顧（六列）。
畫面上的數字都由 utils.walk_saw 算出，check_math() 核對。模組層級的畫圖函式：tug_of_war（拔河圖）、flory_plot、nu_plot。
"""
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
from utils import walk_saw as SW
from utils import lattice_draw as L
from utils import walk_draw as D

M_EX, D_EX, U_EX = 1000, 2, 1.0
ENT_COLOR = BLUE_C
REP_COLOR = ORANGE
SUM_COLOR = YELLOW
YMIN = -250
NS_NU = (6, 8, 10, 12, 14)


def f2(x):
    return f"{x:.2f}"


def f4(x):
    return f"{x:.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


def ent(R):
    return -R * R / M_EX


def rep(R):
    return -U_EX * M_EX * M_EX / R ** D_EX


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    rs = SW.flory_Rstar(M_EX, D_EX)
    return {"rs": rs, "rs_num": SW.flory_Rstar_num(M_EX, D_EX), "F": SW.flory_F(rs, M_EX, D_EX),
            "nu": [SW.flory_nu(d) for d in (1, 2, 3, 4)], "slope": [SW.flory_slope(d) for d in (1, 2, 3, 4)],
            "quoted": [SW.NU_QUOTED[d] for d in (1, 2, 3, 4)], "two_nu": [SW.two_nu_eff(n) for n in NS_NU],
            "saw": SW.sample_saw(10), "split": SW.EX_SPLIT, "c14": SW.saw_c(14), "frac14": SW.saw_frac(14),
            "r2_14": SW.saw_R2(14), "ie4": SW.ie_partial(4), "chords": [SW.chord_count(k) for k in range(1, 5)]}


def check_math():
    nb = numbers()
    rs = nb["rs"]
    assert f2(rs) == "177.83" and abs(nb["rs_num"] - rs) < 1e-3
    assert f2(ent(rs)) == f2(rep(rs)) == "-31.62" and f2(nb["F"]) == f2(ent(rs) + rep(rs)) == "-63.25"
    assert abs(nb["F"] - SW.flory_F(rs, M_EX, D_EX, u=U_EX)) < 1e-12
    assert abs(rs ** (D_EX + 2) - D_EX * U_EX / 2 * M_EX ** 3) < 1e-3 * M_EX ** 3            # R^{d+2} = (du/2) m³
    assert all(SW.flory_F(rs, M_EX, D_EX) >= SW.flory_F(r, M_EX, D_EX) for r in np.linspace(50, 600, 1101))
    assert nb["nu"] == [1, Fraction(3, 4), Fraction(3, 5), Fraction(1, 2)] == [Fraction(3, d + 2) for d in (1, 2, 3, 4)]
    assert [f4(v) for v in nb["slope"]] == ["1.0000", "0.7500", "0.6000", "0.5000"]
    assert nb["quoted"] == [1.0, 0.75, 0.588, 0.5]
    assert [f4(v) for v in nb["two_nu"]] == ["1.4306", "1.4372", "1.4442", "1.4505", "1.4556"]
    assert all(a < b < 1.5 for a, b in zip(nb["two_nu"], nb["two_nu"][1:]))
    pts = nb["saw"]
    assert len(pts) == 11 and len(set(map(tuple, pts))) == 11                               # 自迴避：11 個點兩兩不同
    sp = nb["split"]
    assert sp["1d_total"] == 6 and sp["2d_total"] == 36
    assert nb["c14"] == 2374444 and f4(float(nb["frac14"])) == "0.0088" and f2(float(nb["r2_14"])) == "42.79"
    assert nb["ie4"] == [256, 28, 124, 96, 100] and nb["ie4"][-1] == SW.saw_c(4)
    assert nb["chords"] == [1, 3, 15, 105]
    assert all(SW.power_exp(l, d) == 1 - Fraction(l * (4 - d), 2) for l in range(5) for d in range(1, 6))


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


def assert_curve_clear(curve, mob, gap=0.05):
    """曲線上的每一個點都不在 mob 的外框裡（外框只看曲線經過的點，不看曲線的外框）。"""
    x0, y0 = mob.get_left()[0] - gap, mob.get_bottom()[1] - gap
    x1, y1 = mob.get_right()[0] + gap, mob.get_top()[1] + gap
    for p in curve.get_points():
        assert not (x0 <= p[0] <= x1 and y0 <= p[1] <= y1), (p.round(2), mob)


def fit_cells(tab, ncols, widths, pad=0.25):
    """太寬的格子縮到欄寬以內（只縮不放）。"""
    for k, c in enumerate(tab[0]):
        w = widths[k % ncols] - pad
        if c.get_width() > w:
            c.set_width(w)


def polymer():
    """抽樣的 10 步自迴避漫步（預設種子）＋每個單體一個點＋外接圓與半徑 R。"""
    pts = numbers()["saw"]
    cell = 0.5
    tr = L.trail(pts, cell=cell, origin=ORIGIN, stroke_width=4)
    mono = VGroup(*[Dot(L.lat_point(x, y, cell, ORIGIN), radius=0.08).set_color(WHITE) for x, y in pts])
    c = np.mean([L.lat_point(x, y, cell, ORIGIN) for x, y in pts], axis=0)
    r = max(np.linalg.norm(L.lat_point(x, y, cell, ORIGIN) - c) for x, y in pts) + 0.35
    circ = DashedVMobject(Circle(radius=r).move_to(c), num_dashes=40).set_stroke(GREY_B, 2)
    rad = Arrow(c, c + r * np.array([np.cos(-0.5), np.sin(-0.5), 0]), buff=0, thickness=3).set_color(TEAL_B)
    rl = tx("R", 34, TEAL_B).next_to(rad.get_end(), DR, buff=0.08)
    return VGroup(circ, tr, mono, rad, rl)


def tug_of_war():
    """拔河圖：左隊＝熵（R 越小越好），右隊＝排斥（R 越大越好），中間的繩結＝R。"""
    def team(title, formula, color):
        box = RoundedRectangle(width=3.4, height=1.5, corner_radius=0.15).set_stroke(color, 3).set_fill(color, 0.15)
        body = VGroup(zh(title, font_size=26, color=color), tx(formula, 36, color)).arrange(DOWN, buff=0.18)
        return VGroup(box, body.move_to(box))
    left = team("熵（漫步的亂度）", r"e^{-R^2/m}", ENT_COLOR).move_to(4.6 * LEFT)
    right = team("排斥（自迴避）", r"e^{-um^2/R^d}", REP_COLOR).move_to(4.6 * RIGHT)
    rope = Line(left.get_right(), right.get_left()).set_stroke(GREY_A, 5)
    knot = Square(side_length=0.3).set_fill(RED_B, 1).set_stroke(WHITE, 1.5)
    kl = tx("R", 32).next_to(knot, UP, buff=0.15)
    al = Arrow(0.4 * LEFT + 0.75 * DOWN, 2.4 * LEFT + 0.75 * DOWN, buff=0, thickness=4).set_color(ENT_COLOR)
    ar = Arrow(0.4 * RIGHT + 0.75 * DOWN, 2.4 * RIGHT + 0.75 * DOWN, buff=0, thickness=4).set_color(REP_COLOR)
    tl = mixed("$R$", " 越小越好", font_size=24, color=ENT_COLOR).next_to(al, DOWN, buff=0.15)
    tr = mixed("$R$", " 越大越好", font_size=24, color=REP_COLOR).next_to(ar, DOWN, buff=0.15)
    return VGroup(left, right, rope, VGroup(knot, kl), VGroup(al, tl), VGroup(ar, tr))


def domain(fn, r0=50.0, r1=600.0):
    rs = np.linspace(r0, r1, 2201)
    ok = rs[[fn(r) >= YMIN + 5 for r in rs]]
    return float(ok.min()), float(ok.max())


def flory_plot():
    nb = numbers()
    rs = nb["rs"]
    ax = D.plot_axes((0, 600, 100), (YMIN, 0, 50), width=7.0, height=4.4, y_ticks=[-200, -100], font_size=20,
                     x_label="R")
    ticks = VGroup(*[tx(str(v), 20, GREY_B).move_to(ax.c2p(v, 0) + 0.25 * UP) for v in range(100, 601, 100)])
    curves = VGroup()
    for fn, col in ((ent, ENT_COLOR), (rep, REP_COLOR), (lambda r: ent(r) + rep(r), SUM_COLOR)):
        a, b = domain(fn)
        curves.add(ax.get_graph(fn, x_range=(a, b, 0.5)).set_stroke(col, 4))
    vl = DashedLine(ax.c2p(rs, 0), ax.c2p(rs, nb["F"]), dash_length=0.08).set_stroke(WHITE, 2)
    dots = VGroup(Dot(ax.c2p(rs, ent(rs)), radius=0.08).set_color(WHITE),
                  Dot(ax.c2p(rs, nb["F"]), radius=0.09).set_color(SUM_COLOR))
    rlab = tx(r"R^*=%s" % f2(rs), 24, SUM_COLOR).move_to(ax.c2p(rs, 0) + 0.62 * UP)
    leg = VGroup(mixed("$-R^2/m$", "（熵）", font_size=22, color=ENT_COLOR),
                 mixed("$-um^2/R^d$", "（排斥）", font_size=22, color=REP_COLOR),
                 mixed("兩項的和", font_size=22, color=SUM_COLOR)).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
    leg.move_to(ax.c2p(215, -212))
    return VGroup(ax, ticks, curves, vl, dots, rlab, leg)


def nu_plot():
    nb = numbers()
    # Axes 把軸放在數值 0 的位置，所以座標平移：橫軸 N − 4、縱軸 斜率 − 1.40（刻度數字照原值寫）
    ax = D.plot_axes((0, 12, 2), (0, 0.12, 0.02), width=5.6, height=2.4, x_ticks=[n - 4 for n in NS_NU],
                     x_tick_labels=[str(n) for n in NS_NU], y_ticks=[0, 0.10], y_tick_labels=["1.40", "1.50"],
                     font_size=18, x_label="N")

    def p(n, v):
        return ax.c2p(n - 4, v - 1.40)
    half = DashedLine(p(4, 1.5), p(15.6, 1.5), dash_length=0.08).set_stroke(RED_B, 2)
    hl = tx(r"\frac32", 24, RED_B).next_to(half, RIGHT, buff=0.1)
    dots = VGroup(*[Dot(p(n, v), radius=0.07).set_color(YELLOW) for n, v in zip(NS_NU, nb["two_nu"])])
    labs = VGroup(*[tx(f4(v), 18, YELLOW).next_to(d, UP, buff=0.12) for v, d in zip(nb["two_nu"], dots)])
    return VGroup(ax, half, hl, dots, labs)


# ---------------------------------------------------------------- 各段
def part_question(s):
    poly = polymer()
    q = VGroup(mixed("$m$", " 個單體＝", "$m$", " 步的自迴避漫步", font_size=30),
               mixed("單體不能疊在一起", font_size=30),
               mixed("普通漫步：", "$R\\propto m^{1/2}$", font_size=30),
               mixed("自迴避：", "$R\\propto m^{\\nu}$", "，", "$\\nu=?$", font_size=30, color=YELLOW))
    q.arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    place_in(poly, -6.5, -0.5, -2.8, 3.0)
    place_in(q, 0.3, 6.7, -2.2, 2.4)
    assert_in_frame(poly, q)
    assert_apart(poly, q)
    assert_column(*q)
    s.say("問：", "$m$", " 個單體的高分子，大小 ", "$R$", " 怎麼隨 ", "$m$", " 變大？", wait=0.5)
    s.play(ShowCreation(poly[1][0]), FadeIn(poly[1][1:]), LaggedStart(*[FadeIn(d) for d in poly[2]], lag_ratio=0.1),
           run_time=2.0)
    s.play(ShowCreation(poly[0]), GrowArrow(poly[3]), FadeIn(poly[4]), run_time=1.0)
    s.play(FadeIn(q[0]), FadeIn(q[1]), run_time=0.8)
    s.say("高分子＝自迴避漫步：單體不能疊在一起", wait=2.5)
    s.play(FadeIn(q[2]), run_time=0.6)
    s.say("普通漫步 ", "$\\langle R^2\\rangle=N$", "：", "$R\\propto m^{1/2}$", wait=2.5)
    s.play(FadeIn(q[3]), run_time=0.6)
    s.say("自迴避會把高分子撐得更開：指數 ", "$\\nu$", " 是多少？", wait=3.0)
    s.clear_stage()


def part_tug(s):
    e1 = VGroup(zh("(7.61)", font_size=22, color=GREY_B),
                tx(r"\mathcal N=\alpha R^d,\quad u=f/\alpha:\qquad C_{MF}=\frac{(2d)^m}{\mathcal N}\,e^{-um^2/R^d}", 34))
    e2 = VGroup(zh("(7.62)", font_size=22, color=GREY_B), mixed("熵：", "$e^{-R^2/m}$", font_size=34))
    for e in (e1, e2):
        e.arrange(RIGHT, buff=0.35)
    eqs = VGroup(e1, e2).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    tug = tug_of_war()
    place_in(eqs, -6.7, 6.7, 1.5, 3.15)
    place_in(tug, -6.7, 6.7, -2.9, 0.9)
    assert_in_frame(eqs, tug)
    assert_apart(eqs, tug)
    assert_apart(tug[0], tug[4])
    assert_apart(tug[1], tug[5])
    assert_apart(tug[4], tug[5])
    s.say("(7.61)：把平均場用在半徑 ", "$R$", " 的球裡：", "$\\mathcal N=\\alpha R^d$", wait=0.5)
    s.play(FadeIn(e1), run_time=1.0)
    s.say("$e^{-um^2/R^d}$", "：排斥，", "$R$", " 越大越好", wait=2.5)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("(7.62)：普通漫步走到 ", "$R$", " 的機率 ", "$e^{-R^2/m}$", "：", "$R$", " 越小越好", wait=2.5)
    s.play(FadeIn(tug[2]), FadeIn(tug[3]), run_time=0.6)
    s.play(FadeIn(tug[0]), FadeIn(tug[4]), run_time=0.8)
    s.say("左隊：熵把高分子往回拉", wait=2.0)
    s.play(FadeIn(tug[1]), FadeIn(tug[5]), run_time=0.8)
    s.say("右隊：排斥把高分子撐開", wait=2.0)
    s.play(tug[3].animate.shift(0.4 * RIGHT), run_time=0.6)
    s.play(tug[3].animate.shift(0.6 * LEFT), run_time=0.6)
    s.play(tug[3].animate.shift(0.2 * RIGHT), run_time=0.5)
    s.say("兩邊拔河：繩結停在哪裡？", wait=2.0)
    s.clear_stage()


def part_plot(s):
    nb = numbers()
    pl = flory_plot()
    side = VGroup(VGroup(zh("(7.63)", font_size=20, color=GREY_B), tx(r"F(R)=-\frac{R^2}{m}-\frac{um^2}{R^d}", 32))
                  .arrange(RIGHT, buff=0.25),
                  tx(r"d=2,\ m=1000,\ u=1", 28, GREY_B),
                  tx(r"R^*=%s" % f2(nb["rs"]), 32, SUM_COLOR),
                  mixed("兩項都是 ", "$%s$" % f2(ent(nb["rs"])), font_size=28),
                  mixed("和 ", "$F(R^*)=%s$" % f2(nb["F"]), font_size=28, color=SUM_COLOR),
                  VGroup(zh("(7.64)", font_size=20, color=GREY_B), tx(r"R^{d+2}=\frac{du}{2}\,m^3", 32))
                  .arrange(RIGHT, buff=0.25),
                  tx(r"R\propto m^{3/(d+2)}", 38, YELLOW))
    side.arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    place_in(pl, -6.7, 0.9, -2.9, 3.15)
    place_in(side, 1.4, 6.7, -2.9, 3.0)
    assert_in_frame(pl, side)
    assert_apart(pl, side)
    assert_column(*side)
    for c in pl[2]:
        assert_curve_clear(c, pl[6])
    s.say("(7.63)：兩個因子相乘，取指數再求最大", wait=0.5)
    s.play(FadeIn(side[0]), run_time=0.8)
    s.play(ShowCreation(pl[0]), FadeIn(pl[1]), FadeIn(side[1]), run_time=1.0)
    s.say("例：", "$d=2$", "、", "$m=1000$", "、", "$u=1$", "，", "$R$", " 從 50 到 600", wait=2.0)
    s.play(ShowCreation(pl[2][0]), FadeIn(pl[6][0]), run_time=1.2)
    s.say("藍：熵 ", "$-R^2/m$", "，", "$R$", " 越大扣越多", wait=2.0)
    s.play(ShowCreation(pl[2][1]), FadeIn(pl[6][1]), run_time=1.2)
    s.say("橘：排斥 ", "$-um^2/R^d$", "，", "$R$", " 越小扣越多", wait=2.0)
    s.play(ShowCreation(pl[2][2]), FadeIn(pl[6][2]), run_time=1.2)
    s.say("黃：兩項的和，中間有一個最高點", wait=2.0)
    s.play(ShowCreation(pl[3]), FadeIn(pl[4]), FadeIn(pl[5]), FadeIn(side[2]), run_time=1.0)
    s.say("最大在 ", "$R^*=177.83$", wait=2.0)
    s.play(FadeIn(side[3]), FadeIn(side[4]), run_time=0.8)
    s.say("此時兩項都是 ", "$-31.62$", "，和是 ", "$-63.25$", wait=2.5)
    s.play(FadeIn(side[5]), run_time=0.8)
    s.say("對 ", "$R$", " 微分等於 0：", "$R^{d+2}=\\frac{du}{2}m^3$", wait=2.5)
    s.play(FadeIn(side[6]), run_time=0.8)
    s.say("所以 ", "$R\\propto m^{3/(d+2)}$", "：", "$\\nu=\\frac{3}{d+2}$", wait=2.5)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    q = nb["quoted"]
    rows = [["d", "1", "2", "3", "4"],
            [mixed("Flory ", "$3/(d+2)$", font_size=24)] + [ftex(v) for v in nb["nu"]],
            [zh("數值最大化的斜率", font_size=24)] + [f4(v) for v in nb["slope"]],
            [zh("已知值", font_size=24), tx(ftex(Fraction(q[0]).limit_denominator(10)), 26),
             mixed("$\\frac34$", "（精確）", font_size=22), mixed("$0.588$", "（數值）", font_size=22),
             tx(ftex(Fraction(q[3]).limit_denominator(10)), 26)]]
    widths = [3.6, 1.7, 1.9, 1.9, 1.7]
    tab = D.tex_table(rows, widths, row_h=0.6, font_size=26)
    fit_cells(tab, 5, widths)
    for j in range(1, 5):
        D.table_cell(tab, 1, j, 5).set_color(YELLOW)
        D.table_cell(tab, 3, j, 5).set_color(GREEN_B)
    nu = nu_plot()
    side = VGroup(mixed("二維精確列舉（W702）", font_size=24),
                  mixed("$\\log\\langle R^2\\rangle$", " 對 ", "$\\log N$", " 的局部斜率", font_size=24),
                  mixed("＝", "$2\\nu\\to\\frac32$", font_size=26, color=RED_B)).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
    place_in(tab, -6.7, 6.7, 0.55, 3.15)
    place_in(nu, -6.7, 1.2, -2.9, 0.2)
    place_in(side, 1.7, 6.7, -2.2, -0.2)
    assert_in_frame(tab, nu, side)
    assert_column(tab, nu)
    assert_column(tab, side)
    assert_apart(nu, side)
    assert_table(tab, 5)
    for a, b in zip(nu[4], nu[4][1:]):
        assert_apart(a, b)
    s.say("Flory：", "$\\nu=3/(d+2)$", wait=0.5)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][5:10]], lag_ratio=0.2), run_time=1.0)
    s.say("一維 1、二維 ", "$\\frac34$", "、三維 ", "$\\frac35$", "、四維 ", "$\\frac12$", wait=2.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][10:15]], lag_ratio=0.2), run_time=1.0)
    s.say("數值最大化 ", "$R^*$", "，量 ", "$\\log R^*$", " 對 ", "$\\log m$", " 的斜率：一樣", wait=2.5)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][15:]], lag_ratio=0.2), run_time=1.0)
    s.say("已知值：一維 1、二維 ", "$\\frac34$", "（精確）", wait=2.0)
    s.say("三維 0.588（數值，只引用）、四維 ", "$\\frac12$", wait=2.5)
    s.play(ShowCreation(nu[0]), ShowCreation(nu[1]), FadeIn(nu[2]), FadeIn(side), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(VGroup(d, t)) for d, t in zip(nu[3], nu[4])], lag_ratio=0.25), run_time=1.5)
    s.say("二維精確列舉（W702）：局部斜率＝", "$2\\nu$", wait=2.5)
    s.say("1.4306、1.4372 … 1.4556：慢慢爬向 ", "$\\frac32$", wait=2.5)
    s.say("所以二維 ", "$\\nu=\\frac34$", "：Flory 剛好猜對", wait=2.5)
    s.clear_stage()


def part_recap(s):
    nb = numbers()
    sp = nb["split"]
    ie = nb["ie4"]
    lab = lambda *p: mixed(*p, font_size=28)
    rows = [[r"\text{W701}", lab("路徑和＝切段求和；高斯變異數相加"),
             r"C(4;0)=%d\ (1\text{D}),\ \ %d\ (2\text{D})" % (sp["1d_total"], sp["2d_total"])],
            [r"\text{W702}", lab("自迴避漫步的精確計數"),
             r"c_{14}=%d,\ \ c_{14}/4^{14}=%s,\ \ \langle R^2\rangle=%s\ \text{vs}\ 14" % (nb["c14"], f4(float(nb["frac14"])),
                                                                          f2(float(nb["r2_14"])))],
            [r"\text{W703}", lab("$1-f$", " 展開＝弧線圖"),
             r"(2k-1)!!;\ \ %s" % r",\ ".join(str(v) for v in ie)],
            [r"\text{W704}", lab("生成函數、Dyson 方程"), r"g=\frac{1}{g_0^{-1}-\Sigma}"],
            [r"\text{W705}", lab("冪次數算、平均場"), r"(z_c-z)^{1-l(4-d)/2};\ \ e^{-fm^2/\mathcal N}"],
            [r"\text{W706}", lab("Flory 的拔河"), r"\nu=\frac{3}{d+2}"]]
    widths = [1.4, 5.4, 6.6]
    tab = D.tex_table(rows, widths, row_h=0.8, font_size=30)
    fit_cells(tab, 3, widths)
    for i in range(6):
        D.table_cell(tab, i, 0, 3).set_color(GREY_B)
        D.table_cell(tab, i, 2, 3).set_color(YELLOW)
    rules = VGroup(tab[1][0], tab[1][2])                     # 沒有表頭：不畫第一列下面那條線
    title = zh("第 7 章回顧", font_size=34, color=TEAL_B)
    place_in(title, -6.7, 6.7, 2.65, 3.15)
    body = VGroup(tab[0], rules)
    place_in(body, -6.7, 6.7, -2.8, 2.5)
    assert_in_frame(title, body)
    assert_apart(title, body)
    assert_table(tab, 3)
    caps = [("W701：漫步的條數＝所有路徑的和",),
            ("W702：自迴避的條數越來越少，卻走得越來越遠",),
            ("W703：", "$1-f$", " 展開，部分和最後剛好＝100",),
            ("W704：Dyson 方程把無限多張圖收成一個分母",),
            ("W705：四維以下低階微擾失效；平均場",),
            ("W706：Flory 的拔河，", "$\\nu=\\frac{3}{d+2}$")]
    s.play(FadeIn(title), ShowCreation(rules), run_time=0.8)
    for i in range(6):
        s.play(LaggedStart(*[FadeIn(D.table_cell(tab, i, j, 3)) for j in range(3)], lag_ratio=0.25), run_time=1.0)
        s.say(*caps[i], wait=2.2)
    s.say("熵和排斥拔河，拔出 ", "$R\\propto m^{3/(d+2)}$", wait=2.5)
    s.say("二維剛好是精確的 ", "$\\frac34$", "，三維的 0.6 也很接近 0.588", wait=3.0)


class W706Flory(ChuScene):
    def construct(self):
        self.show_header("W706　Flory 的論證＋第 7 章回顧")
        check_math()
        part_question(self)
        part_tug(self)
        part_plot(self)
        part_table(self)
        part_recap(self)
        self.wait(1.0)
