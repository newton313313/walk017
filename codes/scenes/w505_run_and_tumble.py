"""
W505 二維、三維的「跑一跑、轉個彎」（§5.2.3–5.2.5，pp. 109–118，(5.65)–(5.106)，Fig. 5.4–5.6，Ex. 5.7）

問：大腸桿菌以速率 v 游、平均每單位時間翻滾 γ 次換一個隨機方向（v = γ = 1）。t = 3 時有多少還在第一條直線上？前緣長什麼樣子？
  1. (5.65)(5.69) Γ = γ − γ|0⟩⟨0|；(5.73)–(5.75) 秩一修正的反算子（Ex. 5.7：V.sm_residual）；角度平均 0.4472、0.4636。
  2. 二維 P2 = 2πrρ（V.rt2_radial）＋圓環 e^{−γt}；筆誤框 (5.84)、(5.88)(5.89)（照書上積分 1.2658）。
  3. 三維 P3 = 4πr²ρ（V.rt3_radial），重畫 Fig. 5.6；t = 3 三類人；因果性 (5.101) + (5.103) = 0。
  4. 10 萬隻的模擬（V.rt_mc(d, 3.0, 100000, 5)）對公式；前緣：一維有限、二維 1/√ 發散、三維 ln 發散，都守恆。
畫面上的數字都由 utils.walk_var 算出，check_math() 核對。尖峰（圓環／球殼）的重量寫在表裡與字幕裡。
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
from utils import walk_var as V
from utils import walk_draw as D

T3 = 3.0
R2_PTS = (1, 2, 2.5, 2.9)
FIG_TS = (1, 3, 5, 8)
FIG_COLORS = [YELLOW, ORANGE, TEAL_C, BLUE_C]
U4 = (0.25, 0.5, 0.75, 0.9)
PATH_SEEDS = (1, 2, 3, 4, 6, 7, 8, 9)
PATH_COLORS = [TEAL, BLUE_C, PURPLE_B, GOLD, PINK, GREEN_C, BLUE_B, MAROON_B]
MC_COLOR, F_COLOR = BLUE_C, YELLOW


def f4(x):
    return f"{x:.4f}"


# ---------------------------------------------------------------- 數學
def rt_path2d(seed, t=T3):
    """示意路線：二維 run-and-tumble（v = γ = 1），回傳折點 list。"""
    rng = np.random.default_rng(seed)
    p, pts, tl = np.zeros(2), [np.zeros(2)], t
    while tl > 0:
        ang = rng.uniform(0, 2 * PI)
        step = min(rng.exponential(1.0), tl)
        p = p + step * np.array([math.cos(ang), math.sin(ang)])
        pts.append(p.copy())
        tl -= step
    return pts


def hist_density(vals, nbins=30, hi=T3):
    """[0, hi) 的直方圖（機率密度：個數 ÷ (總數 × 格寬)）；尖峰（放在 vt 外）不算。"""
    h, _ = np.histogram(vals[vals < hi], bins=nbins, range=(0, hi))
    return h / (len(vals) * hi / nbins)


@lru_cache(maxsize=None)
def numbers():
    nb = {"sm": V.sm_residual()}
    nb["g"] = (V.angle_avg_2d_closed(1, 1).real, V.angle_avg_3d_closed(1, 1).real)
    nb["gn"] = (V.angle_avg_2d(1, 1).real, V.angle_avg_3d(1, 1).real)
    nb["p2"] = [float(V.rt2_radial(r, T3)) for r in R2_PTS]
    nb["ring"] = math.exp(-T3)
    nb["m"] = {1: V.tele_mass(T3), 2: V.rt2_mass(T3), 3: V.rt3_mass(T3)}
    nb["m2book"] = V.rt2_mass(T3, book=True)
    nb["fig"] = {t: [float(V.rt3_radial(u * t, t)) for u in U4] for t in FIG_TS}
    nb["shell"] = {t: math.exp(-t) for t in FIG_TS}
    nb["once"], nb["rest"] = 3 * math.exp(-3), 1 - 4 * math.exp(-3)
    nb["pole"], nb["cut"] = V.rt3_pole_book(3.5, 3), V.rt3_cut_book(3.5, 3)
    nb["front1"] = float(V.tele_cont(3 * (1 - 1e-12), 3))
    nb["mc"] = {}
    for d in (2, 3):
        vals = V.rt_mc(d, 3.0, 100000, 5)
        nb["mc"][d] = dict(hist=hist_density(vals), f=V.mc_fraction(vals, 0, 1.5), F=V.radial_fraction(d, 0, 1.5, 3),
                           sh=float(np.mean(vals > 3)))
    return nb


def check_math():
    nb = numbers()
    assert nb["sm"] < 1e-12
    assert [f4(x) for x in nb["g"]] == ["0.4472", "0.4636"] == [f4(x) for x in nb["gn"]]
    assert abs(nb["g"][0] - 1 / math.sqrt(5)) < 1e-12 and abs(nb["g"][1] - math.atan(0.5)) < 1e-12
    assert [f4(x) for x in nb["p2"]] == ["0.2978", "0.4167", "0.3941", "0.4052"]
    assert f4(nb["ring"]) == "0.0498" and f4(nb["m"][2] - nb["ring"]) == "0.9502" == f4(1 - math.exp(-3))
    assert all(f4(nb["m"][d]) == "1.0000" for d in (1, 2, 3))
    assert f4(nb["m2book"]) == "1.2658"
    assert [f4(x) for x in nb["fig"][1]] == ["0.1472", "0.5334", "1.0275", "1.3699"]
    assert [f4(x) for x in nb["fig"][3]] == ["0.1729", "0.4558", "0.4828", "0.3594"]
    assert [f4(x) for x in nb["fig"][5]] == ["0.1869", "0.3630", "0.2169", "0.0926"]
    assert [f4(x) for x in nb["fig"][8]] == ["0.1948", "0.2411", "0.0622", "0.0118"]
    assert [f4(nb["shell"][t]) for t in FIG_TS] == ["0.3679", "0.0498", "0.0067", "0.0003"]
    assert (f4(nb["once"]), f4(nb["rest"])) == ("0.1494", "0.8009")
    assert f4(nb["ring"] + nb["once"] + nb["rest"]) == "1.0000"
    assert (f4(nb["pole"]), f4(nb["cut"])) == ("0.0113", "-0.0113") and abs(nb["pole"] + nb["cut"]) < 1e-5
    assert f"{nb['front1']:.4g}" == "0.06223"
    mc = nb["mc"]
    assert (f4(mc[2]["f"]), f4(mc[2]["F"]), f4(mc[3]["f"]), f4(mc[3]["F"])) == ("0.3335", "0.3310", "0.2880", "0.2881")
    assert (f4(mc[2]["sh"]), f4(mc[3]["sh"])) == ("0.0493", "0.0497")
    assert max(mc[2]["hist"].max(), mc[3]["hist"].max()) < 1.0
    assert all(np.hypot(*np.array(rt_path2d(sd)).T).max() <= T3 + 1e-9 for sd in PATH_SEEDS)


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


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def trunc_graph(ax, fn, xs, ymax, color, sw=3.0, dashed=False):
    """折線圖；超過 ymax 的地方截斷在 ymax（前緣發散）。"""
    pts = [ax.c2p(float(x), min(float(fn(float(x))), ymax)) for x in xs]
    c = VMobject().set_points_as_corners(pts).set_stroke(color, sw)
    if dashed:
        c = DashedVMobject(c, num_dashes=45).set_stroke(color, sw)
    return c


def front_xs(t, n1=100, n2=20, eps0=0.03):
    """[0, t) 的取樣點，在前緣 r → t 附近加密（每條 ≤ 120 點）。"""
    return np.concatenate([np.linspace(0, (1 - eps0) * t, n1), t * (1 - np.geomspace(eps0, 1e-4, n2))[1:]])


def spike_arrow(ax, x, h, color):
    return Arrow(ax.c2p(x, 0), ax.c2p(x, h), buff=0, thickness=3.0).set_color(color)


def hist_bars(ax, dens, hi=T3, color=MC_COLOR):
    w = hi / len(dens)
    bars = VGroup()
    for k, hgt in enumerate(dens):
        p0, p1 = ax.c2p(k * w, 0), ax.c2p((k + 1) * w, float(hgt))
        r = Rectangle(width=p1[0] - p0[0], height=max(p1[1] - p0[1], 1e-3))
        r.set_fill(color, 0.55).set_stroke(color, 1).move_to((p0 + p1) / 2)
        bars.add(r)
    return bars


def path_mob(pts, org, unit, color):
    pp = [org + unit * np.array([x, y, 0]) for x, y in pts]
    line = VMobject().set_points_as_corners(pp).set_stroke(color, 2.5)
    return VGroup(line, Dot(pp[-1], radius=0.05).set_color(color))


# ---------------------------------------------------------------- 各段
def part_question(s):
    org, unit = np.array([-3.4, 0.0, 0]), 0.82
    circle = DashedVMobject(Circle(radius=T3 * unit).move_to(org), num_dashes=60).set_stroke(GREY_B, 2)
    clab = tx(r"r=vt=3", 24, GREY_B).move_to(org + np.array([2.35, 2.35, 0]))
    paths = VGroup(*[path_mob(rt_path2d(sd), org, unit, c) for sd, c in zip(PATH_SEEDS, PATH_COLORS)])
    ang = 0.7
    straight = Line(org, org + T3 * unit * np.array([math.cos(ang), math.sin(ang), 0])).set_stroke(RED, 4)
    src = Dot(org, radius=0.08).set_color(WHITE)
    rows = VGroup(mixed("跑：以速率 ", "$v$", " 直線前進", font_size=28),
                  mixed("翻滾：平均每單位時間 ", "$\\gamma$", " 次", font_size=28),
                  zh("翻滾後：方向完全隨機", font_size=28),
                  mixed("$v=\\gamma=1$", "，看 ", "$t=3$", font_size=28),
                  zh("紅線：一次都沒翻滾", font_size=26, color=RED)).arrange(DOWN, buff=0.42, aligned_edge=LEFT)
    at(rows, 3.3, 0.2)
    assert_in_frame(circle, clab, paths, rows)
    assert_apart(VGroup(circle, clab), rows)
    s.play(FadeIn(src), ShowCreation(circle), FadeIn(clab), run_time=0.8)
    s.say("問：大腸桿菌以速率 ", "$v$", " 游，翻滾後換一個隨機方向", wait=0.3)
    s.play(*[ShowCreation(p[0]) for p in paths], run_time=2.5)
    s.play(*[FadeIn(p[1]) for p in paths], FadeIn(rows[:4]), run_time=0.8)
    s.say("平均每單位時間翻滾 ", "$\\gamma$", " 次（", "$v=\\gamma=1$", "）", wait=2.5)
    s.play(ShowCreation(straight), FadeIn(rows[4]), run_time=1.0)
    s.say("$t=3$", " 時，有多少還在第一條直線上？", wait=2.5)
    s.say("分布在前緣 ", "$r=vt$", " 長什麼樣子？（先猜）", wait=3.0)
    s.clear_stage()


def part_operator(s):
    nb = numbers()
    e1 = eq_label(tx(r"\Gamma=\gamma-\gamma|0\rangle\langle0|", 36), "(5.65)(5.69)")
    e2 = eq_label(tx(r"\big(A-\gamma|0\rangle\langle0|\big)^{-1}=A^{-1}+\frac{\gamma A^{-1}|0\rangle\langle0|A^{-1}}"
                     r"{1-\gamma\langle0|A^{-1}|0\rangle}", 32), "(5.73)--(5.75)")
    e2b = tx(r"A=\lambda+\gamma+i\vec q\cdot\vec v", 28, GREY_A)
    e3 = mixed("Ex. 5.7：角度離散成 64 格，數值核對誤差 ", "$<10^{-12}$", font_size=26)
    top = VGroup(e1, e2, VGroup(e2b, e3).arrange(RIGHT, buff=0.8)).arrange(DOWN, buff=0.3)
    fit(top, 12.5)
    at(top, 0, 1.85)
    hdr = lambda t: zh(t, font_size=24, color=YELLOW)
    rows = [[hdr("維度"), r"\langle0|A^{-1}|0\rangle", r"\lambda=q=1"],
            [hdr("二維"), r"\frac{1}{\sqrt{(\lambda+\gamma)^2+q^2v^2}}", f4(nb["g"][0])],
            [hdr("三維"), r"\frac{1}{vq}\arctan\frac{vq}{\lambda+\gamma}", f4(nb["g"][1])]]
    tab = D.tex_table(rows, [1.2, 4.2, 1.8], row_h=0.78, font_size=28)
    at(tab, -2.7, -1.5)
    e4 = VGroup(tx(r"\langle0|G|0\rangle=\frac{g}{1-\gamma g}", 34, YELLOW),
                tx(r"g=\langle0|A^{-1}|0\rangle", 28, GREY_A)).arrange(DOWN, buff=0.3)
    at(e4, 4.4, -1.5)
    assert_in_frame(top, tab, e4)
    assert_column(*top)
    assert_apart(top, tab)
    assert_apart(top, e4)
    assert_apart(tab, e4)
    assert_table(tab, 3)
    s.play(Write(e1), run_time=1.2)
    s.say("方向是一個「狀態」：翻滾把任何方向變成完全隨機", wait=2.5)
    s.say("$|0\\rangle$", " 是「各方向平均」：", "$\\Gamma$", " 先平均、再以速率 ", "$\\gamma$", " 放回去", wait=2.5)
    s.play(FadeIn(e2), FadeIn(top[2][0]), run_time=1.0)
    s.say("傅立葉＋Laplace 之後要反一個矩陣：秩一修正有公式", wait=2.5)
    s.play(FadeIn(top[2][1]), run_time=0.6)
    s.say("Ex. 5.7：角度切成 64 格當矩陣，誤差小於 ", "$10^{-12}$", wait=2.5)
    s.play(FadeIn(tab[0][:3]), ShowCreation(tab[1]), FadeIn(tab[0][3:]), run_time=1.0)
    s.say("只需要角度平均 ", "$\\langle0|A^{-1}|0\\rangle$", "：二維、三維各有閉式", wait=3.0)
    s.say("$\\lambda=q=1$", "：0.4472、0.4636，和數值角度平均一致", wait=2.5)
    s.play(FadeIn(e4, shift=0.2 * LEFT), run_time=0.8)
    s.say("總數的轉換是 ", "$g/(1-\\gamma g)$", "；反轉的圍道細節放在講義", wait=3.0)
    s.clear_stage()


def part_2d(s):
    nb = numbers()
    ax = D.plot_axes((0, 3.5, 0.5), (0, 1.0, 0.25), width=5.8, height=4.0, x_ticks=[1, 2, 3],
                     y_ticks=[0.25, 0.5, 0.75, 1.0], font_size=20, x_label="r", y_label="P_2")
    ax.move_to(np.array([-3.5, -0.45, 0]))
    xs = np.concatenate([np.linspace(0, 2.9, 60), 3 - np.geomspace(0.1, 1e-5, 40)])
    cur = trunc_graph(ax, lambda r: V.rt2_radial(r, T3), xs, 1.0, YELLOW, 3.5)
    dots = VGroup(*[Dot(ax.c2p(r, p), radius=0.06).set_color(RED) for r, p in zip(R2_PTS, nb["p2"])])
    arr = spike_arrow(ax, T3, 0.75, RED)
    alab = tx(r"e^{-3}", 24, RED).next_to(arr.get_end(), RIGHT, buff=0.12)
    trunc = zh("發散（截斷）", font_size=20, color=GREY_A).next_to(ax.c2p(2.93, 0.95), LEFT, buff=0.1)
    f1 = eq_label(tx(r"P_2=2\pi r\rho=\frac{\gamma r}{v}\,\frac{e^{-\gamma t+\gamma s/v}}{s}", 34), "(5.84)(5.89)", 0.3)
    f2 = mixed("$s=\\sqrt{v^2t^2-r^2}$", "，外加 ", "$r=vt$", " 的圓環 ", "$e^{-\\gamma t}$", font_size=26)
    rows = [["r"] + [str(r) for r in R2_PTS], ["P_2(r,3)"] + [f4(p) for p in nb["p2"]]]
    tab = D.tex_table(rows, [1.5] + [1.15] * 4, row_h=0.6, font_size=26)
    mass = mixed("圓環 ", "$%s$" % f4(nb["ring"]), "＋其餘 ", "$%s$" % f4(nb["m"][2] - nb["ring"]), "＝", "$1$",
                 font_size=28)
    right = VGroup(f1, f2, tab, mass).arrange(DOWN, buff=0.42)
    fit(right, 6.3)
    at(right, 3.45, 1.0)
    b1 = D.erratum_box("p. 112，(5.84)", r"2\gamma\sqrt{\pi/2}\,e^{-\gamma t}\cosh(\gamma\sqrt{v^2t^2-r^2})",
                       r"2\gamma\,e^{-\gamma t}\cosh(\gamma\sqrt{v^2t^2-r^2}/v)", width=10.5)
    b2 = D.erratum_box("p. 113，(5.88)(5.89)", r"\frac{2\pi}{x}\delta(x-y)", r"\frac{1}{x}\delta(x-y)", width=5.4)
    book = mixed("照書上：", "$t=3$", " 的總數 ", "$%s$" % f4(nb["m2book"]), font_size=28, color=GREY_A)
    errs = VGroup(b1, VGroup(b2, book).arrange(RIGHT, buff=0.7)).arrange(DOWN, buff=0.3)
    at(errs, 0, -0.95)
    left = VGroup(ax, cur, arr, alab, trunc)
    assert_in_frame(left, right, errs, dots)
    assert_apart(left, right)
    assert_apart(VGroup(f1, f2), errs)
    assert_column(*right)
    assert_column(*errs)
    assert_table(tab, 5)
    s.play(ShowCreation(ax), FadeIn(f1), run_time=1.0)
    s.say("二維（更正後，極點＋切割合起來）", wait=2.0)
    s.say("$P_2=2\\pi r\\rho$", "：在距離 ", "$r$", " 附近的機率密度", wait=2.0)
    s.play(ShowCreation(cur), FadeIn(f2), run_time=1.8)
    s.play(FadeIn(dots), FadeIn(tab), run_time=0.8)
    s.say("$t=3$", "：", "$r=1,2,2.5,2.9$", " 的值見表", wait=2.5)
    s.play(FadeIn(trunc), run_time=0.5)
    s.say("前緣 ", "$1/\\sqrt{\\ }$", " 發散，但可積（圖上截斷）", wait=2.5)
    s.play(GrowArrow(arr), FadeIn(alab), run_time=0.8)
    s.say("圓環 ", "$r=vt$", "：一次都沒翻滾的人，", "$e^{-3}=0.0498$", wait=2.5)
    s.play(FadeIn(mass), run_time=0.6)
    s.say("連續部分積分是 0.9502，加起來剛好是 1", wait=2.5)
    s.play(FadeOut(left), FadeOut(dots), FadeOut(tab), FadeOut(mass), run_time=0.6)
    s.play(FadeIn(b1, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (5.84) 多了 ", "$\\sqrt{\\pi/2}$", "，cosh 裡少了 ", "$/v$", wait=2.5)
    s.play(FadeIn(errs[1][0], shift=0.2 * UP), run_time=0.8)
    s.say("(5.88)(5.89) 的 ", "$\\delta$", " 多了 ", "$2\\pi$", wait=2.5)
    s.play(FadeIn(errs[1][1]), run_time=0.6)
    s.say("照書上的係數積分，", "$t=3$", " 得 1.2658，不是 1", wait=2.5)
    s.clear_stage()


def part_3d(s):
    nb = numbers()
    ax = D.plot_axes((0, 8.5, 1), (0, 1.6, 0.4), width=6.0, height=4.4, x_ticks=list(range(1, 9)),
                     y_ticks=[0.4, 0.8, 1.2, 1.6], font_size=20, x_label="r", y_label="P_3")
    ax.move_to(np.array([-3.35, -0.3, 0]))
    curves = VGroup(*[trunc_graph(ax, lambda r, t=t: V.rt3_radial(r, t), front_xs(t), 1.6, c, 3.2)
                      for t, c in zip(FIG_TS, FIG_COLORS)])
    arrows = VGroup(*[spike_arrow(ax, t, 0.3, c) for t, c in zip(FIG_TS, FIG_COLORS)])
    legend = VGroup(*[tx(r"t=%d" % t, 24, c) for t, c in zip(FIG_TS, FIG_COLORS)]).arrange(DOWN, buff=0.15,
                                                                                       aligned_edge=LEFT)
    legend.move_to(ax.c2p(7.3, 1.25))
    once = trunc_graph(ax, lambda r: V.rt3_once_radial(r, 3), front_xs(3), 1.6, ORANGE, 2.5, dashed=True)
    head = mixed("$P_3=4\\pi r^2\\rho$", "：(5.98)＋(5.101)＋(5.103)", font_size=26)
    hdr = lambda t: zh(t, font_size=22, color=YELLOW)
    rows = [["t", r"\frac t4", r"\frac t2", r"\frac{3t}4", "0.9t", hdr("球殼")]]
    rows += [[str(t)] + [f4(v) for v in nb["fig"][t]] + [f4(nb["shell"][t])] for t in FIG_TS]
    tab = D.tex_table(rows, [0.6] + [1.05] * 5, row_h=0.5, font_size=24)
    for k, c in enumerate(FIG_COLORS):
        D.table_cell(tab, k + 1, 0, 6).set_color(c)
    hz = lambda t, c=WHITE: zh(t, font_size=22, color=c)
    rows2 = [[hz("t=3", YELLOW), hz("比例", YELLOW), hz("在哪裡", YELLOW)],
             [hz("沒翻滾"), r"e^{-3}=%s" % f4(nb["ring"]), hz("球殼 r=vt")],
             [hz("翻一次"), r"3e^{-3}=%s" % f4(nb["once"]), hz("虛線；前緣 ln 發散")],
             [hz("其餘"), r"1-4e^{-3}=%s" % f4(nb["rest"]), hz("連續的主體")]]
    tab2 = D.tex_table(rows2, [1.2, 2.6, 2.4], row_h=0.48, font_size=24)
    caus = mixed("$r=3.5>vt$", "：(5.101) ", "$=%s$" % f4(nb["pole"]), "，(5.103) ", "$=%s$" % f4(nb["cut"]),
                 font_size=24)
    right = VGroup(head, tab, tab2, caus).arrange(DOWN, buff=0.3)
    fit(right, 6.4)
    at(right, 3.5, 0.1)
    left = VGroup(ax, curves, arrows, legend)
    assert_in_frame(left, right)
    assert_apart(left, right)
    assert_column(*right)
    assert_table(tab, 6)
    assert_table(tab2, 3)
    s.play(ShowCreation(ax), FadeIn(head), run_time=1.0)
    s.say("三維：徑向密度 ", "$P_3=4\\pi r^2\\rho$", "（重畫 Fig. 5.6）", wait=2.0)
    for k in range(4):
        s.play(ShowCreation(curves[k]), GrowArrow(arrows[k]), FadeIn(legend[k]), run_time=1.0)
    s.play(FadeIn(tab[0][:6]), ShowCreation(tab[1]), FadeIn(tab[0][6:]), run_time=1.0)
    s.say("$r=\\frac t4,\\frac t2,\\frac{3t}4,0.9t$", " 的值見表；箭頭是球殼", wait=3.0)
    s.say("走得越久，越像原點附近的高斯；球殼越來越輕", wait=2.5)
    s.play(ShowCreation(once), FadeIn(tab2), run_time=1.2)
    s.say("$t=3$", " 的三類人：沒翻滾、翻一次、其餘", wait=2.5)
    s.say("翻一次＝(5.103) 最後一項：前緣 ", "$\\ln$", " 發散，也可積", wait=2.5)
    s.say("0.0498＋0.1494＋0.8009＝1", wait=2.0)
    s.play(FadeIn(caus), run_time=0.6)
    s.say("因果性：", "$r=3.5>vt$", " 時兩項 ", "$\\pm0.0113$", " 相消為 0", wait=2.5)
    s.clear_stage()


def part_mc(s):
    nb = numbers()
    mc = nb["mc"]
    groups = []
    for d, x0, fn in ((2, -3.4, V.rt2_radial), (3, 3.4, V.rt3_radial)):
        ax = D.plot_axes((0, 3.5, 0.5), (0, 1.0, 0.25), width=5.4, height=2.8, x_ticks=[1, 2, 3],
                         y_ticks=[0.5, 1.0], font_size=20, x_label="r", y_label="P_%d" % d)
        ax.move_to(np.array([x0, -1.0, 0]))
        bars = hist_bars(ax, mc[d]["hist"])
        xs = front_xs(3) if d == 3 else np.concatenate([np.linspace(0, 2.9, 60), 3 - np.geomspace(0.1, 1e-5, 40)])
        cur = trunc_graph(ax, lambda r, fn=fn: fn(r, T3), xs, 1.0, F_COLOR, 3)
        title = zh("二維" if d == 2 else "三維", font_size=28, color=YELLOW)
        l1 = mixed("$r<1.5$", "：模擬 %s，公式 %s" % (f4(mc[d]["f"]), f4(mc[d]["F"])), font_size=24)
        l2 = mixed("圓環" if d == 2 else "球殼", "：模擬 %s，公式 %s" % (f4(mc[d]["sh"]), f4(nb["ring"])), font_size=24)
        txt = VGroup(title, l1, l2).arrange(DOWN, buff=0.22)
        at(txt, x0, 1.85)
        assert_in_frame(ax, bars, cur, txt)
        assert_apart(txt, VGroup(ax, bars))
        assert_column(*txt)
        groups.append((ax, bars, cur, txt))
    assert_apart(VGroup(*groups[0]), VGroup(*groups[1]))
    s.play(*[ShowCreation(g[0]) for g in groups], *[FadeIn(g[3][0]) for g in groups], run_time=1.0)
    s.say("模擬 10 萬隻（種子 5）：", "$t=3$", " 時離原點的距離", wait=1.0)
    s.play(*[LaggedStart(*[FadeIn(b) for b in g[1]], lag_ratio=0.03) for g in groups], run_time=2.0)
    s.play(*[ShowCreation(g[2]) for g in groups], run_time=1.5)
    s.say("藍：直方圖；黃：公式。前緣的尖角也對上了", wait=2.5)
    s.play(*[FadeIn(g[3][1]) for g in groups], run_time=0.6)
    s.say("$r<1.5$", " 的比例：二維模擬 0.3335，公式 0.3310", wait=2.5)
    s.say("三維模擬 0.2880，公式 0.2881", wait=2.0)
    s.play(*[FadeIn(g[3][2]) for g in groups], run_time=0.6)
    s.say("剛好在 ", "$r=vt$", " 上的：0.0493、0.0497，公式 0.0498", wait=2.5)
    s.clear_stage()
    hz = lambda t, c=WHITE: zh(t, font_size=26, color=c)
    rows = [[hz("維度", YELLOW), hz("前緣", YELLOW), hz("總數", YELLOW)],
            [hz("一維"), mixed("有限：", "$\\Sigma(vt^-)=%s$" % f"{nb['front1']:.4g}", font_size=26), f4(nb["m"][1])],
            [hz("二維"), mixed("$1/\\sqrt{\\ }$", " 發散", font_size=26), f4(nb["m"][2])],
            [hz("三維"), mixed("$\\ln$", " 發散", font_size=26), f4(nb["m"][3])]]
    tab = D.tex_table(rows, [1.6, 4.4, 1.8], row_h=0.75, font_size=30)
    at(tab, 0, 0.3)
    assert_in_frame(tab)
    assert_table(tab, 3)
    s.play(FadeIn(tab[0][:3]), ShowCreation(tab[1]), run_time=0.8)
    for k in range(1, 4):
        s.play(FadeIn(tab[0][3 * k:3 * k + 3]), run_time=0.6)
    s.say("前緣：一維有限、二維 ", "$1/\\sqrt{\\ }$", "、三維 ", "$\\ln$", " 發散", wait=2.5)
    s.say("都可積：總數（連續部分＋尖峰）都是 1", wait=2.5)


def part_close(s):
    s.say("維度越高，前緣越尖", wait=2.0)
    s.say("但沒有人跑得比 ", "$v$", " 快，也一個都不少", wait=2.5)


class W505RunAndTumble(ChuScene):
    def construct(self):
        self.show_header("W505　二維、三維的「跑一跑、轉個彎」")
        check_math()
        part_question(self)
        part_operator(self)
        part_2d(self)
        part_3d(self)
        part_mc(self)
        part_close(self)
        self.wait(1.0)
