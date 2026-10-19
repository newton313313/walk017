"""
W1101 碎形：放大之後還一樣（§11.1–11.1.1，pp. 255–261，(11.1)–(11.6)，Fig. 11.1–11.5，Ex. 11.1–11.3）

問：什麼東西放大之後還長得一樣？
  1. 點、直線（平凡）；圓弧越放大越直（有自己的尺度＝半徑）；跳蚤身上還有更小的跳蚤（不引用原詩）。
  2. Koch 曲線第 0–5 代（段數 4^k、總長 (4/3)^k → ∞）；第 5 代前四分之一放大 3 倍＝整條。Cantor 集第 0–5 代（2^k 段、(2/3)^k → 0）。
  3. 質量－半徑（Fig. 11.5）：n(R/3) = n(R)/4 ⇒ p = ln4/ln3 = 1.2619；Cantor 0.6309；盒子計數 22、83、325、1291 → 1.2362（Cantor 0.6254）。
  4. 複本公式 (11.6) p = ln ℛ / ln ℳ：線段、正方形、Koch、Cantor、Sierpinski。
  5. 普通漫步 d_f = 2（400 條：2.0258、公式 2.0069）；自迴避 d_f = 1/ν：4/3、1.70、2。
畫面上的數字都由 utils.walk_rg 算出，check_math() 核對。
"""
import sys
import warnings
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rg as RG
from utils import walk_draw as D

KOCH_COLOR = TEAL_B
ZOOM_COLOR = YELLOW
CANTOR_COLOR = "#7FB8FF"
SEG_TEX = ("1", "4", "16", "64", "256", "1024")
KLEN = ("1", "1.3333", "1.7778", "2.3704", "3.1605", "4.2140")
CLEN = ("1", "0.6667", "0.4444", "0.2963", "0.1975", "0.1317")
KOCH_N = ("1", "0.2500", "0.0625", "0.0156", "0.0039")
CANTOR_N = ("1", "0.5000", "0.2500", "0.1250", "0.0625")
BOXES = (22, 83, 325, 1291)
REP_ROWS = (("line", "線段", "1.0000"), ("square", "正方形", "2.0000"), ("koch", "Koch 曲線", "1.2619"),
            ("cantor", "Cantor 集", "0.6309"), ("sierpinski", "Sierpinski 三角形", "1.5850"))
SAW_EFF = ("1.3848", "1.3788", "1.3740")


def f4(x):
    return f"{float(x):.4f}"


def num(x):
    """1 寫成 1，其他寫四位小數。"""
    return "1" if abs(float(x) - 1) < 1e-12 else f4(x)


def walk_fit(exact=False):
    """ln N 對 ln R_g 的點與最小平方斜率、截距（和 RG.rw_df 同一組 N）。"""
    xs, ys = [], []
    for N in RG.RW_NS:
        rg2 = float(RG.rg2_exact(N)) if exact else float(RG.rw_rg2(N))
        xs.append(0.5 * np.log(rg2))
        ys.append(np.log(N))
    slope, icpt = np.polyfit(xs, ys, 1)
    return xs, ys, slope, icpt


def box_fit():
    bc = RG.box_counts("koch", RG.BOX_K, 3)
    xs = [np.log(1 / e) for e, _ in bc]
    ys = [np.log(n) for _, n in bc]
    slope, icpt = np.polyfit(xs, ys, 1)
    return xs, ys, slope, icpt


def check_math():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(6):
            n, length, eps = RG.koch_counts(k)
            assert str(n) == SEG_TEX[k] == str(4 ** k) and num(length) == KLEN[k] and eps == Fraction(1, 3 ** k)
            cn, clen, _ = RG.cantor_counts(k)
            assert cn == 2 ** k and num(clen) == CLEN[k] and len(RG.cantor_intervals(k)) == 2 ** k
            assert len(RG.koch_points(k)) == 4 ** k + 1
        pts = np.array(RG.koch_points(5), dtype=float)
        assert len(pts) == 1025 and np.allclose(3 * pts[:257], pts[::4])       # 前四分之一放大 3 倍 = 整條
        assert pts[:, 1].min() >= -1e-12 and pts[:, 1].max() > 0.28             # 突起朝上
        mk, mc = RG.mass_table("koch"), RG.mass_table("cantor")
        assert [num(n) for _, _, n in mk] == list(KOCH_N) and [num(n) for _, _, n in mc] == list(CANTOR_N)
        assert all(abs(R - 3.0 ** -k) < 1e-12 for k, R, _ in mk)
        assert f4(RG.LN4_LN3) == "1.2619" == f4(np.log(4) / np.log(3)) and f4(RG.LN2_LN3) == "0.6309"
        assert f4(RG.mass_radius_dim("koch")) == "1.2619" and f4(RG.mass_radius_dim("cantor")) == "0.6309"
        bc = RG.box_counts("koch", RG.BOX_K, 3)
        assert tuple(n for _, n in bc) == BOXES
        assert all(abs(e - 3.0 ** -k) < 1e-12 for (e, _), k in zip(bc, RG.BOX_K))
        assert f4(box_fit()[2]) == f4(RG.box_dim("koch")) == "1.2362" and f4(RG.box_dim("cantor")) == "0.6254"
        for key, _, p in REP_ROWS:
            assert f4(RG.replica_dim(*RG.REPLICAS[key])) == p
        assert RG.REPLICAS == {"line": (3, 3), "square": (9, 3), "koch": (4, 3), "cantor": (2, 3), "sierpinski": (3, 2)}
        trail = RG.rw_trail()
        assert len(trail) == 1501 and 1 / RG.saw_df(2) == Fraction(3, 4)
        assert RG.rg2_exact(4) == Fraction(4 * 6, 6 * 5)
        assert f4(walk_fit()[2]) == f4(RG.rw_df()) == "2.0258"
        assert f4(walk_fit(True)[2]) == f4(RG.rw_df(exact=True)) == "2.0069"
        assert RG.saw_df(2) == Fraction(4, 3) and f"{RG.saw_df(3):.2f}" == "1.70" and RG.saw_df(4) == 2
        assert tuple(f4(RG.saw_df_eff(N)) for N in RG.SAW_NEFF) == SAW_EFF and RG.SAW_NEFF == (10, 12, 14)


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_left(m, x0, y, x1=6.7):
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    return m.move_to(np.array([x0 + m.get_width() / 2, y, 0]))


def fit_width(m, w):
    if m.get_width() > w:
        m.set_width(w)
    return m


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def poly_mob(pts, color=KOCH_COLOR, width=2.5):
    P = [np.array([x, y, 0.0]) for x, y in pts]
    return VMobject().set_points_as_corners(P).set_stroke(color, width)


def koch_mob(gen, x0, x1, y, color=KOCH_COLOR, width=2.5, count=None):
    pts = np.array(RG.koch_points(gen, p0=(x0, y), p1=(x1, y)), dtype=float)
    if count is not None:
        pts = pts[:count]
    return poly_mob(pts, color, width)


def cantor_row(gen, x0, x1, y, color=CANTOR_COLOR, width=8):
    size = x1 - x0
    return VGroup(*[Line([x0 + float(a) * size, y, 0], [x0 + float(b) * size, y, 0]).set_stroke(color, width)
                    for a, b in RG.cantor_intervals(gen)])


def sierpinski(gen, x0, y0, side, color=ORANGE):
    tris = [(np.array([x0, y0]), side)]
    for _ in range(gen):
        new = []
        for p, s in tris:
            h = s / 2
            new += [(p, h), (p + np.array([h, 0]), h), (p + np.array([h / 2, h * np.sqrt(3) / 2]), h)]
        tris = new
    g = VGroup()
    for p, s in tris:
        a, b, c = p, p + np.array([s, 0]), p + np.array([s / 2, s * np.sqrt(3) / 2])
        g.add(Polygon(*[np.array([v[0], v[1], 0]) for v in (a, b, c)]).set_stroke(color, 1).set_fill(color, 0.8))
    return g


def rep_icon(key, center):
    c = np.array(center, dtype=float)
    if key == "line":
        g = VGroup(Line(c + 0.6 * LEFT, c + 0.6 * RIGHT).set_stroke(WHITE, 3))
        for t in (-0.6, -0.2, 0.2, 0.6):
            g.add(Line(c + t * RIGHT + 0.08 * UP, c + t * RIGHT + 0.08 * DOWN).set_stroke(WHITE, 2))
        return g
    if key == "square":
        g = VGroup()
        for i in range(3):
            for j in range(3):
                g.add(Square(0.18).set_stroke(WHITE, 1.5).move_to(c + (i - 1) * 0.18 * RIGHT + (j - 1) * 0.18 * UP))
        return g
    if key == "koch":
        return koch_mob(2, c[0] - 0.6, c[0] + 0.6, c[1] - 0.15, width=2)
    if key == "cantor":
        return VGroup(cantor_row(1, c[0] - 0.6, c[0] + 0.6, c[1] + 0.1, width=5),
                      cantor_row(2, c[0] - 0.6, c[0] + 0.6, c[1] - 0.12, width=5))
    return sierpinski(2, c[0] - 0.3, c[1] - 0.26, 0.6)


# ---------------------------------------------------------------- 1. 放大之後
def part_intro(s):
    xs, y = (-4.5, 0.0, 4.5), 0.9
    dot = Dot([xs[0], y, 0], radius=0.09)
    line = Line([xs[1] - 1.5, y, 0], [xs[1] + 1.5, y, 0]).set_stroke(WHITE, 4)

    def arc_of(R):
        ang = 3.0 / R
        a = Arc(start_angle=PI / 2 - ang / 2, angle=ang, radius=R).set_stroke(ORANGE, 4)
        return a.shift(np.array([xs[2], y + 0.45, 0]) - a.point_from_proportion(0.5))

    arc = arc_of(1.0)
    labs = VGroup(*[zh(t, font_size=28).move_to([x, -1.0, 0]) for t, x in zip(("點", "直線", "圓弧"), xs)])
    zoom = zh("放大 3 倍", font_size=28, color=YELLOW).move_to([0, 2.6, 0])
    zoom9 = zh("再放大 3 倍", font_size=28, color=YELLOW).move_to([0, 2.6, 0])
    arc3, arc9 = arc_of(3.0), arc_of(9.0)
    assert_in_frame(dot, line, arc, labs, zoom)
    for a, b in ((dot, line), (line, arc), (arc, labs[2]), (line, zoom)):
        assert_apart(a, b)

    s.say("什麼東西放大之後還長得一樣？", wait=1.5)
    s.play(FadeIn(dot), ShowCreation(line), ShowCreation(arc), FadeIn(labs), run_time=1.2)
    s.say("點放大還是點，直線放大還是直線", wait=2.0)
    s.play(FadeIn(zoom), Transform(arc, arc3), run_time=1.5)
    s.say("圓弧放大 3 倍：變直了一些", wait=1.5)
    s.play(FadeTransform(zoom, zoom9), Transform(arc, arc9), run_time=1.5)
    s.say("圓弧越放大越直：它有自己的尺度——半徑", wait=2.5)
    s.say("點和直線沒有尺度，但也太平凡了：沒有細節", wait=2.5)
    s.say("跳蚤身上還有更小的跳蚤——一層又一層", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 2. Koch 與 Cantor
def koch_table():
    rows = [["k"] + [str(k) for k in range(6)],
            [zh("段數", font_size=26)] + list(SEG_TEX),
            [mixed("總長 ", r"$(4/3)^k$", font_size=24)] + list(KLEN)]
    return D.tex_table(rows, [2.4] + [1.3] * 6, row_h=0.55, font_size=28).shift(np.array([-5.4, -0.45, 0]))


def part_koch(s):
    x0, x1, y = -3.5, 3.5, 0.3
    curves = [koch_mob(g, x0, x1, y) for g in range(6)]
    tab = koch_table()
    inf = tx(r"\to\infty", 30, YELLOW).next_to(D.table_cell(tab, 2, 6, 7), RIGHT, buff=0.3)
    assert_in_frame(tab, inf, *curves)
    assert_table(tab, 7)
    assert_apart(curves[5], tab)

    s.play(ShowCreation(curves[0]), FadeIn(tab[1]), *[FadeIn(D.table_cell(tab, i, j, 7)) for i in range(3) for j in (0, 1)],
           run_time=1.0)
    s.say("Koch 曲線：從一條線段開始", wait=1.5)
    notes = ("每段換成 4 段，每段長 1/3，中間突起", "再做一次：16 段", "64 段", "256 段", "第 5 代：1024 段，1025 個點")
    for g in range(1, 6):
        s.play(ReplacementTransform(curves[g - 1], curves[g]),
               *[FadeIn(D.table_cell(tab, i, g + 1, 7)) for i in range(3)], run_time=1.2)
        s.say(notes[g - 1], wait=0.7)
    s.play(FadeIn(inf), run_time=0.5)
    s.say("段數 ", "$4^k$", "，總長 ", "$(4/3)^k$", "：越來越長，趨向無窮大", wait=2.5)
    full = curves[5]
    piece = koch_mob(5, x0, x1, y, ZOOM_COLOR, 4, count=257)
    qp = np.array(RG.koch_points(5, p0=(x0, y), p1=(x1, y)), dtype=float)[:257]
    big = poly_mob(qp[0] + 3 * (qp - qp[0]), ZOOM_COLOR, 4)
    s.play(FadeOut(tab), FadeOut(inf), run_time=0.5)
    s.play(ShowCreation(piece), run_time=1.0)
    s.say("取第 5 代的前四分之一", wait=1.0)
    s.play(Transform(piece, big), run_time=2.0)
    s.say("放大 3 倍，和整條疊在一起：一模一樣", wait=2.5)
    s.play(FadeOut(piece), FadeOut(full), run_time=0.6)

    cx0, cx1 = -6.2, -0.4
    ys = [2.45 - 0.88 * k for k in range(6)]
    rows = VGroup(*[cantor_row(k, cx0, cx1, ys[k]) for k in range(6)])
    head = VGroup(zh("段數", font_size=26, color=YELLOW).move_to([1.6, 2.98, 0]),
                  mixed("總長 ", r"$(2/3)^k$", font_size=24, color=YELLOW).move_to([4.2, 2.98, 0]))
    labs = VGroup()
    for k in range(6):
        labs.add(VGroup(tx(f"k={k}", 26, GREY_B).move_to([0.2, ys[k], 0]),
                        tx(str(2 ** k), 28).move_to([1.6, ys[k], 0]),
                        tx(CLEN[k], 28).move_to([4.2, ys[k], 0])))
    zero = tx(r"\to0", 30, YELLOW).next_to(labs[5][2], RIGHT, buff=0.3)
    assert_in_frame(rows, head, labs, zero)
    assert_apart(rows, labs)
    for k in range(5):
        assert_apart(labs[k], labs[k + 1])
    assert_apart(head, labs[0])

    s.say("Cantor 集：去掉中間的三分之一，一直做下去", wait=1.0)
    s.play(FadeIn(head), ShowCreation(rows[0]), FadeIn(labs[0]), run_time=0.8)
    for k in range(1, 6):
        s.play(TransformFromCopy(rows[k - 1], rows[k]), FadeIn(labs[k]), run_time=0.7)
    s.play(FadeIn(zero), run_time=0.5)
    s.say("$2^k$", " 段，總長 ", "$(2/3)^k$", "：趨向 0，卻有無窮多個點", wait=2.5)
    s.say("兩個都是：放大 3 倍，就得到自己的複本", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 3. 質量－半徑
def mass_table():
    mk, mc = RG.mass_table("koch"), RG.mass_table("cantor")
    rtex = ["1", r"\frac13", r"\frac19", r"\frac1{27}", r"\frac1{81}"]
    rows = [["R"] + rtex,
            [zh("Koch", font_size=24)] + [num(n) for _, _, n in mk],
            [zh("Cantor", font_size=24)] + [num(n) for _, _, n in mc]]
    return D.tex_table(rows, [1.5] + [1.15] * 5, row_h=0.62, font_size=26).shift(np.array([-0.55, 3.0, 0]))


def part_mass(s):
    O = np.array([-6.3, -2.0, 0])
    Lk = 5.0
    curve = koch_mob(5, O[0], O[0] + Lk, O[1])
    inner = koch_mob(5, O[0], O[0] + Lk, O[1], ZOOM_COLOR, 3.5, count=257)
    a0, a1 = -8 * DEGREES, 42 * DEGREES
    big = Arc(start_angle=a0, angle=a1 - a0, radius=Lk, arc_center=O).set_stroke(GREY_B, 2.5)
    small = Arc(start_angle=a0, angle=a1 - a0, radius=Lk / 3, arc_center=O).set_stroke(ZOOM_COLOR, 2.5)
    lb = tx("R", 30, GREY_B).next_to(big.get_end(), UR, buff=0.08)
    ls = tx("R/3", 28, ZOOM_COLOR).next_to(small.get_end(), UP, buff=0.1)
    cen = Dot(O, radius=0.06)
    fig = VGroup(curve, inner, big, small, lb, ls, cen)
    tab = mass_table()
    x0 = -0.5
    e1 = place_left(mixed("(11.1)  ", r"$n(R/3)=\frac14\,n(R)$", font_size=28), x0, 0.55)
    e2 = place_left(mixed(r"$n=KR^p$", r"$\Rightarrow$", r"$3^{-p}=\frac14$", font_size=28), x0, -0.3)
    e3 = place_left(tx(r"p=\frac{\ln4}{\ln3}=" + f4(RG.LN4_LN3), 34, YELLOW), x0, -1.2)
    e4 = place_left(mixed("Cantor：", r"$n(R/3)=\frac12\,n(R)$", r"$\Rightarrow$", r"$p=\frac{\ln2}{\ln3}=" + f4(RG.LN2_LN3) + "$",
                          font_size=26), x0, -2.2)
    assert_in_frame(fig, tab, e1, e2, e3, e4)
    assert_table(tab, 6)
    for a, b in ((fig, tab), (fig, e1), (fig, e2), (fig, e3), (fig, e4), (tab, e1), (e1, e2), (e2, e3), (e3, e4)):
        assert_apart(a, b)

    s.say("維度怎麼算？數一數：半徑 ", "$R$", " 裡面有多少份量", wait=1.5)
    s.play(ShowCreation(curve), FadeIn(cen), run_time=1.2)
    s.play(ShowCreation(big), FadeIn(lb), run_time=0.8)
    s.play(ShowCreation(small), FadeIn(ls), ShowCreation(inner), run_time=1.0)
    s.say("半徑 ", "$R/3$", " 裡面剛好是前四分之一", wait=1.5)
    s.play(FadeIn(tab[1]), *[FadeIn(D.table_cell(tab, i, j, 6)) for i in (0, 1) for j in range(6)], run_time=1.0)
    s.play(FadeIn(e1), run_time=0.6)
    s.say("半徑縮成三分之一，份量變成四分之一", wait=1.5)
    s.play(FadeIn(e2), run_time=0.6)
    s.say("假設 ", "$n=KR^p$", "，代進去 ", "$K$", " 消掉", wait=1.5)
    s.play(FadeIn(e3), run_time=0.6)
    s.say("Koch 曲線的碎形維度 ", f"$p={f4(RG.LN4_LN3)}$", wait=1.5)
    s.play(*[FadeIn(D.table_cell(tab, 2, j, 6)) for j in range(6)], FadeIn(e4), run_time=0.8)
    s.say("Cantor 集：份量每次減半，", f"$p={f4(RG.LN2_LN3)}$", wait=2.5)
    s.clear_stage()

    bx, by, slope, icpt = box_fit()
    ax = D.plot_axes((0, 6, 1), (0, 8, 2), width=5.4, height=4.4, x_ticks=(2, 4, 6), y_ticks=(2, 4, 6, 8),
                     x_label=r"\ln(1/\varepsilon)", y_label=r"\ln N")
    ax.move_to([-2.9, 0.15, 0])
    dots = VGroup(*[Dot(ax.c2p(x, y), radius=0.07).set_color(YELLOW) for x, y in zip(bx, by)])
    fit = Line(ax.c2p(1.5, slope * 1.5 + icpt), ax.c2p(6, slope * 6 + icpt)).set_stroke(TEAL_B, 2.5)
    eps = ["\\varepsilon", r"\frac19", r"\frac1{27}", r"\frac1{81}", r"\frac1{243}"]
    rows = [eps, [zh("盒子數", font_size=22)] + [str(n) for _, n in RG.box_counts("koch", RG.BOX_K, 3)]]
    btab = D.tex_table(rows, [1.1, 0.9, 0.9, 0.9, 0.95], row_h=0.62, font_size=26).shift(np.array([1.6, 2.5, 0]))
    t1 = place_left(mixed("斜率 ", f"${f4(RG.box_dim('koch'))}$", "（", r"$\frac{\ln4}{\ln3}=" + f4(RG.LN4_LN3) + "$", "）",
                          font_size=28), 1.6, 0.4)
    t2 = place_left(mixed("Cantor：", f"${f4(RG.box_dim('cantor'))}$", "（", f"${f4(RG.LN2_LN3)}$", "）", font_size=28),
                    1.6, -0.5)
    t3 = place_left(zh("格子不會剛好對齊，有限大小的誤差", font_size=24, color=GREY_B), 1.6, -1.4)
    assert_in_frame(ax, btab, t1, t2, t3)
    assert_table(btab, 5)
    for a, b in ((ax, btab), (ax, t1), (ax, t2), (btab, t1), (t1, t2), (t2, t3)):
        assert_apart(a, b)

    s.say("第二種算法：盒子計數", wait=1.0)
    s.play(FadeIn(ax), FadeIn(btab), run_time=1.0)
    s.say("用邊長 ", r"$\varepsilon$", " 的方格蓋住曲線，數有幾格", wait=1.5)
    s.play(FadeIn(dots, lag_ratio=0.3), run_time=1.0)
    s.play(ShowCreation(fit), FadeIn(t1), run_time=1.0)
    s.say("斜率 ", f"${f4(RG.box_dim('koch'))}$", "：接近 ", f"${f4(RG.LN4_LN3)}$", wait=1.5)
    s.play(FadeIn(t2), FadeIn(t3), run_time=0.6)
    s.wait(1.5)
    s.clear_stage()


# ---------------------------------------------------------------- 4. 複本公式
def part_replicas(s):
    head = [zh("圖形", font_size=24, color=YELLOW), r"\mathcal R", r"\mathcal M", "p"]
    rows = [head]
    for key, name, _ in REP_ROWS:
        nrep, mag = RG.REPLICAS[key]
        rows.append([fit_width(zh(name, font_size=24), 2.7), str(nrep), str(mag), f4(RG.replica_dim(nrep, mag))])
    tab = D.tex_table(rows, [2.9, 1.1, 1.1, 1.6], row_h=0.72, font_size=28).shift(np.array([-5.3, 2.95, 0]))
    icons = VGroup(*[rep_icon(key, [-6.05, 2.95 - 0.72 * (i + 1.5), 0]) for i, (key, _, _) in enumerate(REP_ROWS)])
    f = tx(r"p=\frac{\ln\mathcal R}{\ln\mathcal M}", 44, YELLOW)
    f.move_to([4.3, 1.9, 0])
    tag = tx("(11.6)", 26, GREY_B).next_to(f, UP, buff=0.3)
    l1 = mixed("放大 ", r"$\mathcal M$", " 倍", font_size=28)
    l2 = mixed("得到 ", r"$\mathcal R$", " 份自己", font_size=28)
    l3 = tx(r"\mathcal M^{\,p}=\mathcal R", 34)
    for m, yy in ((l1, 0.5), (l2, -0.2), (l3, -1.0)):
        m.move_to([4.3, yy, 0])
    l4 = mixed("Koch：", r"$1<p<2$", font_size=26).move_to([4.3, -2.0, 0])
    assert_in_frame(tab, icons, f, tag, l1, l2, l3, l4)
    assert_table(tab, 4)
    for a, b in ((tab, f), (tab, l1), (tab, l3), (tab, l4), (icons, tab), (f, l1), (l1, l2), (l2, l3), (l3, l4), (tag, f)):
        assert_apart(a, b)

    s.say("更快的算法：數複本", wait=1.0)
    s.play(FadeIn(f), FadeIn(tag), run_time=0.8)
    s.play(FadeIn(l1), FadeIn(l2), FadeIn(l3), run_time=0.8)
    s.say("放大 ", r"$\mathcal M$", " 倍，得到 ", r"$\mathcal R$", " 份原來的自己", wait=2.0)
    s.play(FadeIn(tab[1]), *[FadeIn(D.table_cell(tab, 0, j, 4)) for j in range(4)], run_time=0.6)
    says = (("線段：放大 3 倍是 3 段，", "$p=1$"), ("正方形：放大 3 倍是 9 塊，", "$p=2$"),
            ("Koch：放大 3 倍是 4 份，", r"$p=\frac{\ln4}{\ln3}$"), ("Cantor：放大 3 倍是 2 份，", r"$p=\frac{\ln2}{\ln3}$"),
            ("Sierpinski：放大 2 倍是 3 份，", r"$p=\frac{\ln3}{\ln2}$"))
    for i in range(5):
        s.play(FadeIn(icons[i]), *[FadeIn(D.table_cell(tab, i + 1, j, 4)) for j in range(4)], run_time=0.7)
        s.say(*says[i], wait=1.1)
    s.play(FadeIn(l4), run_time=0.5)
    s.say("維度介於 1 和 2 之間——介於線和面之間", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 5. 漫步也是碎形
def part_walk(s):
    pts = np.array(RG.rw_trail(), dtype=float)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    u = 5.6 / max(hi - lo)
    ctr = np.array([-5.1, 0.3])
    P = [np.array([*(ctr + u * (p - (lo + hi) / 2)), 0.0]) for p in pts]
    walk = VMobject().set_points_as_corners(P).set_stroke(YELLOW, 1.4)
    ends = VGroup(Dot(P[0], radius=0.07).set_color(GREEN), Dot(P[-1], radius=0.07).set_color(RED))
    wlab = zh("1500 步", font_size=22, color=GREY_B).next_to(walk, DOWN, buff=0.12)

    xs, ys, slope, icpt = walk_fit()
    ax = D.plot_axes((0, 4, 1), (0, 9, 3), width=3.4, height=3.0, x_ticks=(1, 2, 3), y_ticks=(3, 6, 9),
                     x_label=r"\ln R_g", y_label=r"\ln N")
    ax.move_to([-1.5, 1.35, 0])
    dots = VGroup(*[Dot(ax.c2p(x, y), radius=0.06).set_color(YELLOW) for x, y in zip(xs, ys)])
    fit = Line(ax.c2p(0.8, slope * 0.8 + icpt), ax.c2p(3.6, slope * 3.6 + icpt)).set_stroke(TEAL_B, 2.5)
    x0 = 1.0
    r1 = place_left(mixed("$N$", " 步走 ", r"$R\sim N^{1/2}$", r"$\Rightarrow$", r"$N\propto R^2$", font_size=28), x0, 2.6)
    r2 = place_left(tx(r"\langle R_g^2\rangle=\frac{N(N+2)}{6(N+1)}", 32), x0, 1.7)
    r3 = place_left(mixed("400 條：斜率 ", f"${f4(RG.rw_df())}$", font_size=28), x0, 0.75)
    r4 = place_left(mixed("用公式：", f"${f4(RG.rw_df(exact=True))}$", r"$\ \to\ 2$", font_size=28), x0, 0.05)
    s1 = place_left(mixed("自迴避漫步：", r"$d_f=1/\nu$", font_size=28), -3.3, -0.85)
    s2 = place_left(mixed("二維 ", r"$d_f=\frac43$", "（精確列舉 ", "$N=10,12,14$", "：",
                          f"${SAW_EFF[0]},\\ {SAW_EFF[1]},\\ {SAW_EFF[2]}$", "）", font_size=26), -3.3, -1.65)
    s3 = place_left(mixed("三維 ", f"${RG.saw_df(3):.2f}$", "（", r"$\nu=0.588$", "）；四維 ", f"${RG.saw_df(4)}$",
                          "（普通漫步）", font_size=26), -3.3, -2.45)
    assert_in_frame(walk, ends, wlab, ax, r1, r2, r3, r4, s1, s2, s3)
    for a, b in ((walk, ax), (walk, s1), (walk, s2), (walk, s3), (wlab, s3), (ax, r1), (ax, r3), (ax, s1),
                 (r1, r2), (r2, r3), (r3, r4), (r4, s1), (s1, s2), (s2, s3)):
        assert_apart(a, b)

    s.say("漫步也是碎形嗎？", wait=1.0)
    s.play(ShowCreation(walk), FadeIn(ends[0]), run_time=2.5)
    s.play(FadeIn(ends[1]), FadeIn(wlab), run_time=0.4)
    s.say("一條普通漫步：放大一段，看起來還是一條漫步", wait=1.6)
    s.play(FadeIn(r1), run_time=0.6)
    s.say("$N$", " 步走 ", r"$R\sim N^{1/2}$", "：份量 ", r"$N\propto R^2$", wait=1.6)
    s.play(FadeIn(r2), run_time=0.6)
    s.say("迴轉半徑的平方平均有精確公式", wait=1.8)
    s.play(FadeIn(ax), FadeIn(dots, lag_ratio=0.3), run_time=1.0)
    s.play(ShowCreation(fit), FadeIn(r3), FadeIn(r4), run_time=1.0)
    s.say("斜率 ", f"${f4(RG.rw_df())}$", "、公式 ", f"${f4(RG.rw_df(exact=True))}$", "：碎形維度 2", wait=2.5)
    s.play(FadeIn(s1), run_time=0.6)
    s.say("自迴避漫步：", r"$R\sim N^\nu$", "，所以 ", r"$d_f=1/\nu$", wait=1.6)
    s.play(FadeIn(s2), run_time=0.6)
    s.say("二維 ", r"$\nu=\frac34$", "：", r"$d_f=\frac43$", "，精確列舉的有效值慢慢靠近", wait=2.5)
    s.play(FadeIn(s3), run_time=0.6)
    s.say("三維 ", f"${RG.saw_df(3):.2f}$", "；四維回到 ", "$2$", wait=1.6)
    s.say("碎形維度＝放大 ", r"$\mathcal M$", " 倍、份量變 ", r"$\mathcal R$", " 倍的對數比", wait=2.5)
    s.say("重整化群：把這個想法用在漲落的系統上", wait=2.5)


class W1101Fractals(ChuScene):
    def construct(self):
        self.show_header("W1101　碎形：放大之後還一樣")
        check_math()
        part_intro(self)
        part_koch(self)
        part_mass(self)
        part_replicas(self)
        part_walk(self)
        self.wait(1.0)
