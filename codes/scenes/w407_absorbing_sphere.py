"""
W407 吸收球：細胞長不大；電容（§4.2.4–4.2.5，pp. 89–91，(4.61)–(4.66)，Fig. 4.6–4.7）

問：細胞靠擴散吃營養。細胞半徑加倍，吃進來的營養變幾倍？需要的呢？
  1. (4.61)–(4.63)：c = c∞(1 − r0/r)、j = −β∇c；流量 4πβc∞r0，任何半徑的球面都一樣（BC.sphere_flux_numeric = 4π）。
  2. 吃進 ∝ r0、需要 ∝ r0³ → r_max = √(3Dc∞/q)；假設值例子 3.9×10⁻⁵ m（BC.cell_rmax(**BC.CELL_EXAMPLE)）。
  3. 球外的點源（Fig. 4.6）：鏡像在 r0²/R、強度 −r0/R（BC.kelvin_image）；被吃掉 r0/R = 0.25、逃走 0.75。
  4. 電容 (4.64)–(4.66)：−∮j·dS/c∞ = 4πβC，球 C = r0；筆誤框 ρ = −4π∇²φ → ∇²φ = −4πρ。
畫面上的數字都由 utils.walk_bc 算出，check_math() 核對。
"""
import math
import sys
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_bc as BC
from utils import walk_draw as D

CELL_FILL, CELL_STROKE = GREEN_E, GREEN_C
IN_COLOR, NEED_COLOR = TEAL_B, RED_B
FLUX_RS = (1.5, 3.0, 6.0)
R0, RSRC = 1.0, 4.0                         # 球半徑、點源距離（Fig. 4.6）
KC, KS = np.array([-0.5, -0.1, 0]), 0.7     # Fig. 4.6 的球心（畫面）與比例；畫的時候左右翻轉（源在球的左邊）
KBOX = (-1.6, 8.8, -3.8, 4.2)               # 示意路線離開這個框（資料座標；x = 8.8 是畫面左緣）就算逃走


# ---------------------------------------------------------------- 數學
def brownian_trail(seed, step=0.15, nmax=4000):
    """從 (RSRC, 0) 出發的平面示意路線：碰到球就停（hit），離開框就停（out）。"""
    rng = np.random.default_rng(seed)
    p = np.array([RSRC, 0.0])
    pts = [p.copy()]
    for ang in rng.uniform(0, TAU, nmax):
        p = p + step * np.array([math.cos(ang), math.sin(ang)])
        pts.append(p.copy())
        if math.hypot(p[0], p[1]) <= R0:
            return pts, "hit"
        if not (KBOX[0] < p[0] < KBOX[1] and KBOX[2] < p[1] < KBOX[3]):
            pts.pop()
            return pts, "out"
    return pts, "long"


@lru_cache(maxsize=None)
def pick_trails(n_each=2):
    hits, outs = [], []
    for sd in range(200):
        pts, kind = brownian_trail(sd)
        if kind == "hit" and len(hits) < n_each:
            hits.append(pts)
        elif kind == "out" and len(outs) < n_each and len(pts) > 150:
            outs.append(pts)
        if len(hits) == n_each and len(outs) == n_each:
            break
    return tuple(hits), tuple(outs)


def ex_strings(ex):
    """CELL_EXAMPLE → 畫面上的 TeX（數值由字典產生，check_math 核對）。"""
    mant, expo = f"{ex['D']:.0e}".split("e")
    return (r"D=%s\times10^{%d}\,\mathrm{m^2/s}" % (mant, int(expo)),
            r"c_\infty=%g\,\mathrm{mol/m^3}" % ex["cinf"],
            r"q=%g\,\mathrm{mol/(m^3\,s)}" % ex["q"])


@lru_cache(maxsize=None)
def numbers():
    nb = {"flux": {r: BC.sphere_flux_numeric(r, R0) for r in FLUX_RS}, "four_pi": BC.sphere_flux(R0, 1.0, 1.0)}
    nb["rmax"] = BC.cell_rmax(**BC.CELL_EXAMPLE)
    nb["in_ratio"] = BC.sphere_flux(2.0, 1.0, 1.0) / BC.sphere_flux(1.0, 1.0, 1.0)
    nb["need_ratio"] = (2.0 ** 3) / (1.0 ** 3)
    nb["img"] = BC.kelvin_image(RSRC, R0)
    nb["hit"], nb["hit_num"] = BC.hit_prob(RSRC, R0), BC.hit_prob_numeric(RSRC, R0)
    nb["esc"] = (RSRC - R0) / RSRC
    return nb


def check_math():
    nb = numbers()
    assert all(f"{v:.3f}" == "12.566" for v in nb["flux"].values()) and f"{nb['four_pi']:.3f}" == "12.566"
    assert abs(nb["four_pi"] - 4 * math.pi) < 1e-12
    assert all(abs(v - 4 * math.pi) < 1e-6 for v in nb["flux"].values())
    assert BC.conc_sphere(R0, R0, 1.0) == 0 and abs(BC.conc_sphere(2.0, R0, 1.0) - 0.5) < 1e-12
    assert (round(nb["in_ratio"]), round(nb["need_ratio"])) == (2, 8) and abs(nb["in_ratio"] - 2) < 1e-12
    assert BC.CELL_EXAMPLE == {"D": 2e-9, "cinf": 0.25, "q": 1.0}
    assert ex_strings(BC.CELL_EXAMPLE) == (r"D=2\times10^{-9}\,\mathrm{m^2/s}", r"c_\infty=0.25\,\mathrm{mol/m^3}",
                                           r"q=1\,\mathrm{mol/(m^3\,s)}")
    assert f"{nb['rmax']:.1e}" == "3.9e-05" and round(nb["rmax"] * 1e6) == 39
    assert abs(nb["rmax"] - math.sqrt(3 * 2e-9 * 0.25 / 1.0)) < 1e-15
    assert nb["img"] == (0.25, -0.25) and nb["img"] == (R0 ** 2 / RSRC, -R0 / RSRC)
    assert BC.sphere_surface_max(RSRC, R0) < 1e-9
    assert nb["hit"] == 0.25 and f"{nb['hit_num']:.2f}" == "0.25" and abs(nb["hit_num"] - 0.25) < 1e-6
    assert f"{nb['esc']:.2f}" == "0.75" and abs(nb["esc"] + nb["hit"] - 1) < 1e-12
    assert BC.capacitance_sphere(R0) == R0
    assert abs(math.sqrt(3) ** 3 / 3 - math.sqrt(3)) < 1e-12 and 2.12 ** 3 / 3 < 3.2         # 吃進＝需要於 √3
    hits, outs = pick_trails()
    assert len(hits) == 2 and len(outs) == 2


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def eq_label(tex_mob, label, buff=0.4):
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


def cell_circle(center, radius):
    return Circle(radius=radius).move_to(center).set_fill(CELL_FILL, 0.6).set_stroke(CELL_STROKE, 3)


def inward_arrows(center, radius, n=8, length=0.55, gap=0.12, color=TEAL):
    out = VGroup()
    for k in range(n):
        d = np.array([math.cos(TAU * (k + 0.5) / n), math.sin(TAU * (k + 0.5) / n), 0])
        out.add(Arrow(center + (radius + gap + length) * d, center + (radius + gap) * d, buff=0,
                      stroke_width=3).set_color(color))
    return out


def kpt(x, y):
    return KC + KS * (-x * RIGHT + y * UP)


def trail_mob(pts, color, every=3):
    sub = pts[::every] + [pts[-1]]
    m = VMobject().set_points_as_corners([kpt(p[0], p[1]) for p in sub]).set_stroke(color, 1.6, opacity=0.85)
    return VGroup(m, Dot(kpt(*pts[-1]), radius=0.05).set_color(color))


# ---------------------------------------------------------------- 各段
def part_question(s):
    c1, c2 = np.array([-3.4, 0.2, 0]), np.array([2.9, 0.2, 0])
    s1, s2 = cell_circle(c1, 0.6), cell_circle(c2, 1.2)
    rad1 = Line(c1, c1 + 0.6 * RIGHT).set_stroke(WHITE, 2)
    rad2 = Line(c2, c2 + 1.2 * RIGHT).set_stroke(WHITE, 2)
    l1 = tx("r_0", 28).next_to(rad1, UP, buff=0.06)
    l2 = tx("2r_0", 28).next_to(rad2, UP, buff=0.06)
    a1, a2 = inward_arrows(c1, 0.6), inward_arrows(c2, 1.2)
    rng = np.random.default_rng(7)
    food = VGroup()
    while len(food) < 55:
        p = np.array([rng.uniform(-6.4, 6.4), rng.uniform(-2.6, 2.4), 0])
        if np.linalg.norm(p - c1) > 1.5 and np.linalg.norm(p - c2) > 2.1:
            food.add(Dot(p, radius=0.04).set_color(GOLD))
    tag = zh("營養（擴散而來）", font_size=24, color=GOLD).move_to(np.array([0, 2.75, 0]))
    assert_in_frame(s1, s2, a1, a2, food, tag)
    assert_apart(VGroup(s1, a1), VGroup(s2, a2))
    s.play(FadeIn(food), FadeIn(tag), run_time=0.8)
    s.play(FadeIn(s1), ShowCreation(rad1), FadeIn(l1), run_time=0.8)
    s.say("問：細胞靠擴散吃營養（只吃碰到表面的）", wait=0.3)
    s.play(LaggedStart(*[GrowArrow(a) for a in a1], lag_ratio=0.1), run_time=1.2)
    s.wait(1.0)
    s.play(FadeIn(s2), ShowCreation(rad2), FadeIn(l2), LaggedStart(*[GrowArrow(a) for a in a2], lag_ratio=0.1),
           run_time=1.4)
    s.say("細胞半徑加倍，吃進來的營養變幾倍？", wait=2.5)
    s.say("需要的營養又變幾倍？（先猜）", wait=3.0)
    s.clear_stage()


def part_profile(s):
    nb = numbers()
    ax = D.plot_axes((0, 8, 1), (0, 1.2, 0.25), width=5.5, height=3.6, x_ticks=[1, 2, 4, 6, 8], y_ticks=[0.5, 1],
                     y_tick_labels=[r"\tfrac12c_\infty", r"c_\infty"], font_size=20, x_label="r", y_label="c")
    ax.move_to(np.array([-3.5, -0.35, 0]))
    cell = Polygon(ax.c2p(0, 0), ax.c2p(1, 0), ax.c2p(1, 1.2), ax.c2p(0, 1.2)).set_stroke(width=0)
    cell.set_fill(CELL_FILL, 0.45)
    cell_lab = zh("細胞", font_size=20, color=GREEN_B).move_to(ax.c2p(0.5, 0.6))
    cinf = DashedLine(ax.c2p(0, 1), ax.c2p(8, 1)).set_stroke(GREY_B, 2)
    curve = ax.get_graph(lambda r: BC.conc_sphere(r, R0, 1.0), x_range=(1, 8, 0.05)).set_stroke(YELLOW, 4)
    e1 = eq_label(tx(r"c(r)=c_\infty\Big(1-\frac{r_0}{r}\Big)", 32), "(4.61)(4.62)")
    e2 = eq_label(tx(r"\vec j=-\beta\nabla c=-\hat r\,\frac{\beta c_\infty r_0}{r^2}", 32), "(4.63)")
    e3 = tx(r"-\oint\vec j\cdot d\vec S=4\pi r^2\cdot\frac{\beta c_\infty r_0}{r^2}=4\pi\beta c_\infty r_0", 28, YELLOW)
    rows = [["r"] + [f"{r:g}" for r in FLUX_RS],
            [zh("流量", font_size=24)] + [f"{nb['flux'][r]:.3f}" for r in FLUX_RS]]
    tab = D.tex_table(rows, [1.5, 1.4, 1.4, 1.4], row_h=0.55, font_size=28)
    note = mixed("$\\beta=c_\\infty=r_0=1$", "：數值積分都是 ", "$4\\pi=%.3f$" % nb["four_pi"], font_size=22)
    right = VGroup(e1, e2, e3, tab, note).arrange(DOWN, buff=0.38)
    fit(right, 6.2)
    at(right, 3.55, 0.1)
    assert_in_frame(ax, right)
    assert_apart(ax, right)
    assert_table(tab, 4)
    for k in range(4):
        assert_apart(right[k], right[k + 1])
    s.play(ShowCreation(ax), FadeIn(cell), FadeIn(cell_lab), run_time=1.0)
    s.say("把細胞當成吸收球：碰到表面的營養就被吃掉", wait=2.5)
    s.play(ShowCreation(cinf), FadeIn(e1), run_time=0.8)
    s.play(ShowCreation(curve), run_time=1.2)
    s.say("穩定態：球面上 ", "$c=0$", "，遠處 ", "$c=c_\\infty$", "；中間像 ", "$1/r$", wait=3.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("Fick 定律：流量密度指向球心，大小 ", "$\\propto1/r^2$", wait=3.0)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("乘上球面積 ", "$4\\pi r^2$", "：", "$r$", " 消掉了", wait=3.0)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), FadeIn(note), run_time=1.0)
    s.say("任何半徑的球面，流進來的量都一樣：", "$4\\pi\\beta c_\\infty r_0$", wait=3.0)
    s.clear_stage()


def part_cell(s):
    nb = numbers()
    ax = D.plot_axes((0, 2.6, 0.5), (0, 3.2, 1), width=5.4, height=3.6, font_size=20, x_label="r_0")
    ax.move_to(np.array([-3.6, -0.3, 0]))
    r_cross = math.sqrt(3)
    intake = ax.get_graph(lambda u: u, x_range=(0, 2.6, 0.05)).set_stroke(IN_COLOR, 4)
    need = ax.get_graph(lambda u: u ** 3 / 3, x_range=(0, 2.12, 0.02)).set_stroke(NEED_COLOR, 4)
    il = zh("吃進", font_size=24, color=IN_COLOR).move_to(ax.c2p(0.7, 1.8))
    nl = zh("需要", font_size=24, color=NEED_COLOR).move_to(ax.c2p(2.25, 0.7))
    vline = DashedLine(ax.c2p(r_cross, 0), ax.c2p(r_cross, r_cross)).set_stroke(YELLOW, 2)
    xdot = Dot(ax.c2p(r_cross, r_cross), radius=0.07).set_color(YELLOW)
    xl = tx(r"r_{\max}", 26, YELLOW).next_to(ax.c2p(r_cross, 0), DOWN, buff=0.12)
    t1 = mixed("吃進：", "$4\\pi Dc_\\infty r_0$", font_size=26)
    t2 = mixed("需要：", "$\\tfrac43\\pi r_0^3\\,q$", "（", "$q$", "：每單位體積的消耗率）", font_size=26)
    t3 = mixed("半徑加倍：吃進 ", "$\\times%d$" % round(nb["in_ratio"]), "，需要 ",
               "$\\times%d$" % round(nb["need_ratio"]), font_size=28)
    t4 = tx(r"r_{\max}=\sqrt{3Dc_\infty/q}", 36, YELLOW)
    ds, cs, qs = ex_strings(BC.CELL_EXAMPLE)
    exb = VGroup(zh("例（假設值，只看數量級）", font_size=22, color=GREY_A),
                 tx(ds + r",\ \ " + cs, 26), tx(qs, 26),
                 tx(r"r_{\max}\approx%s\times10^{%d}\,\mathrm{m}\approx%d\,\mu\mathrm{m}"
                    % (f"{nb['rmax']:.1e}".split("e")[0], int(f"{nb['rmax']:.1e}".split("e")[1]),
                       round(nb["rmax"] * 1e6)), 30, YELLOW)).arrange(DOWN, buff=0.18)
    frame = SurroundingRectangle(exb, buff=0.18).set_stroke(GREY_B, 1.5)
    ex = VGroup(frame, exb)
    right = VGroup(t1, t2, t3, t4, ex).arrange(DOWN, buff=0.36)
    fit(right, 6.3)
    at(right, 3.5, 0.15)
    assert_in_frame(ax, xl, right, il, nl)
    assert_apart(VGroup(ax, xl), right)
    for k in range(4):
        assert_apart(right[k], right[k + 1])
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(ShowCreation(intake), FadeIn(il), FadeIn(t1), run_time=1.0)
    s.say("吃進來的：", "$4\\pi Dc_\\infty r_0$", "，只跟半徑成正比", wait=3.0)
    s.play(ShowCreation(need), FadeIn(nl), FadeIn(t2), run_time=1.0)
    s.say("需要的：跟體積成正比，", "$\\propto r_0^3$", wait=2.5)
    s.play(FadeIn(t3), run_time=0.7)
    s.say("所以半徑加倍：吃進 ×2，需要 ×8", wait=3.0)
    s.play(ShowCreation(vline), FadeIn(xdot), FadeIn(xl), FadeIn(t4), run_time=1.0)
    s.say("兩條線交叉以後就餓肚子：有個最大半徑", wait=3.0)
    s.play(FadeIn(ex, shift=0.2 * UP), run_time=0.8)
    s.say("假設值代進去：約 ", "$39\\,\\mu\\mathrm{m}$", "（只看數量級）", wait=3.0)
    s.say("大一點的生物就要靠血液、對流來送營養", wait=3.0)
    s.clear_stage()


def part_kelvin(s):
    nb = numbers()
    hits, outs = pick_trails()
    trails = VGroup(*[trail_mob(p, RED_B) for p in hits], *[trail_mob(p, GREY_B) for p in outs])
    ball = cell_circle(KC, R0 * KS)
    src = Dot(kpt(RSRC, 0), radius=0.1).set_color(ORANGE)
    img = Dot(kpt(nb["img"][0], 0), radius=0.07).set_color(BLUE_B)
    src_lab = zh("源", font_size=24, color=ORANGE).next_to(src, UP, buff=0.12)
    img_lab = zh("鏡像源", font_size=22, color=BLUE_B).move_to(KC + 1.15 * DOWN)
    img_ln = Line(img_lab.get_top(), img.get_center()).set_stroke(BLUE_B, 1.5)
    rline = DashedLine(KC, kpt(RSRC, 0)).set_stroke(GREY_B, 1.5)
    r_lab = tx(r"R=%g" % RSRC, 24, GREY_A).next_to(kpt(2.5, 0), DOWN, buff=0.1)
    k1 = mixed("鏡像源：位置 ", "$\\rho=r_0^2/R$", "，強度 ", "$-r_0/R$", font_size=26)
    k2 = tx(r"r_0=%g,\ R=%g:\ \ \rho=%g,\ \ -r_0/R=%g" % (R0, RSRC, nb["img"][0], nb["img"][1]), 30)
    k3 = mixed("球面上處處 ", "$c=0$", font_size=26)
    k4 = mixed("被吃掉＝球面流量 ", "$=%.2f=r_0/R$" % nb["hit_num"], font_size=28, color=RED_B)
    k5 = mixed("逃到無窮遠 ", "$\\frac{R-r_0}{R}=%.2f$" % nb["esc"], font_size=28, color=GREY_A)
    k6 = mixed("高斯定律：總強度 ", "$1-r_0/R$", font_size=24, color=GREY_B)
    right = VGroup(k1, k2, k3, k4, k5, k6).arrange(DOWN, buff=0.4)
    fit(right, 5.4)
    at(right, 3.85, 0.1)
    left = VGroup(trails, ball, src, img, src_lab, img_lab, r_lab)
    assert_in_frame(left, right)
    assert_apart(left, right)
    for k in range(5):
        assert_apart(right[k], right[k + 1])
    assert_apart(img_lab, r_lab)
    s.play(FadeIn(ball), FadeIn(src), FadeIn(src_lab), ShowCreation(rline), FadeIn(r_lab), run_time=1.0)
    s.say("Fig. 4.6：營養從球外一點放出", wait=2.0)
    s.play(*[ShowCreation(t[0]) for t in trails], run_time=3.0)
    s.play(*[FadeIn(t[1]) for t in trails], run_time=0.4)
    s.add(rline, src, src_lab, r_lab)          # 標籤蓋在路線上面
    s.say("有的碰到球被吃掉（紅），有的逃走（灰）", wait=2.5)
    s.play(FadeIn(img), FadeIn(img_lab), ShowCreation(img_ln), FadeIn(k1), run_time=1.0)
    s.say("球面上要 ", "$c=0$", "：在球裡放一個負的鏡像源", wait=2.5)
    s.play(FadeIn(k2), FadeIn(k3), run_time=0.8)
    s.say("$r_0=1,\\ R=4$", "：鏡像在 0.25、強度 ", "$-0.25$", wait=3.0)
    s.play(FadeIn(k4), run_time=0.7)
    s.say("被吃掉的比例＝球面流量＝", "$r_0/R$", wait=2.5)
    s.play(FadeIn(k5), FadeIn(k6), run_time=0.8)
    s.say("其餘 0.75 逃到無窮遠：高斯定律數總強度", wait=3.0)
    s.clear_stage()


def part_capacitor(s):
    c = np.array([-3.6, -0.1, 0])
    inner = Circle(radius=0.8).move_to(c).set_fill(GREY_D, 1).set_stroke(GREY_A, 3)
    outer = DashedVMobject(Circle(radius=2.6).move_to(c), num_dashes=60).set_stroke(GREY_B, 2)
    arrows = inward_arrows(c, 0.8, n=10, length=1.2, gap=0.2, color=TEAL)
    il = mixed("$c=0$", font_size=22).move_to(c)
    ol = mixed("遠處球殼 ", "$c=c_\\infty$", font_size=22, color=GREY_A).next_to(outer, UP, buff=0.1)
    m1 = eq_label(tx(r"\frac{-\oint\vec j\cdot d\vec S}{c_\infty}=4\pi\beta C", 34, YELLOW), "(4.64)--(4.66)")
    m2 = mixed("球：", "$C=r_0\\ \\Rightarrow\\ -\\oint\\vec j\\cdot d\\vec S=4\\pi\\beta r_0c_\\infty$", font_size=24)
    m3 = mixed("吸收體", "$\\leftrightarrow$", "內導體，", "$c_\\infty-c\\leftrightarrow$", "電位", font_size=24,
               color=GREY_A)
    box = D.erratum_box("p. 91，(4.64) 前", r"\rho=-4\pi\nabla^2\phi", r"\nabla^2\phi=-4\pi\rho", width=5.4)
    right = VGroup(m1, m2, m3, box).arrange(DOWN, buff=0.5)
    fit(right, 6.2)
    at(right, 3.5, 0.0)
    left = VGroup(inner, outer, arrows, il, ol)
    assert_in_frame(left, right)
    assert_apart(left, right)
    for k in range(3):
        assert_apart(right[k], right[k + 1])
    s.play(FadeIn(inner), FadeIn(il), ShowCreation(outer), FadeIn(ol), run_time=1.0)
    s.say("Fig. 4.7：吸收體在裡面，遠處維持 ", "$c_\\infty$", wait=2.5)
    s.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.08), run_time=1.2)
    s.play(FadeIn(m1), run_time=0.8)
    s.say("就是一個電容器：內導體＋遠處的球殼", wait=2.0)
    s.say("總流量除以濃度差＝", "$4\\pi\\beta C$", "：電容 ", "$C$", " 只看形狀", wait=3.0)
    s.play(FadeIn(m2), FadeIn(m3), run_time=0.8)
    s.say("球的電容 ", "$C=r_0$", "：回到 ", "$4\\pi\\beta r_0c_\\infty$", wait=3.0)
    s.say("不是球也一樣：吃進的速率只看吸收體的電容", wait=3.0)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("(4.64) 前的 Poisson 方程寫反了；(4.64) 本身是對的", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("擴散送來的營養跟半徑成正比", wait=2.0)
    s.say("需要的跟體積成正比：只靠擴散的細胞長不大", wait=3.0)


class W407AbsorbingSphere(ChuScene):
    def construct(self):
        self.show_header("W407　吸收球：細胞長不大；電容")
        check_math()
        part_question(self)
        part_profile(self)
        part_cell(self)
        part_kelvin(self)
        part_capacitor(self)
        part_close(self)
        self.wait(1.0)
