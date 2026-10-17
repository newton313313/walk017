"""
W208 離開格子：長尾步長讓 √N 失效（§2.1.4，pp. 36–38，(2.23)–(2.31)，Ex. 2.7）

問：如果偶爾會跳很遠很遠，走 N 步後離家多遠？還是 √N 嗎？
  1. (2.23) χ(q) = ∫ e^{iq·r} p(r) d³r；(2.24) 動差展開（書上左邊的 dq 應為 d³r，字幕帶過）；(2.25) χ(0) = 1。
     Ex. 2.7 高斯步長 → χ = e^{−q²r0²/2} → N 步仍是高斯，每軸變異數 N r0²；筆誤框 √(2/π) r0^{−3} → (2π)^{−3/2} r0^{−3}
     （照印的 ∫p d³r = 4π ≈ 12.566，GF.radial_integral）。
  2. 長尾 (2.26) p = r0/(π²(r²+r0²)²)：先做角度積分 (2.28)，再在上半平面繞 r = i r0 的二階極點 → (2.30) χ = e^{−q r0}
     （q = r0 = 1：0.3679）；e^{−q r0} ≈ 1 − r0|q|：從 |q| 開始變（和高斯步長的 χ 並排畫）。
  3. (2.31) P(N; r) = N r0/(π²(r²+N²r0²)²)：寬度乘 N。中位距離（GF.MEDIANS）：高斯 1.45、4.92、15.79；Cauchy 2.5、23.28、236.06。
  4. 兩條 200 步足跡的 (x, y) 投影（種子 25；卡片上的 7 有一次 1365 的大跳，其餘縮成一點，所以換掉）。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import random
import sys
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

SEED = 25
GAUSS_C = D.GAUSS_COLOR
CAUCHY_C = TEAL
BOX = 4.5


def trails():
    tg = GF.trail_3d(GF.steps_gauss_3d(200, 1.0, random.Random(SEED)))
    tc = GF.trail_3d(GF.steps_cauchy_3d(200, 1.0, random.Random(SEED)))
    return tg, tc


def fit_box(pts, center):
    """(x, y) 投影縮放到 BOX×BOX 的方框裡 → (cell, origin, 方框邊長代表的距離)。"""
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    side = max(max(xs) - min(xs), max(ys) - min(ys))
    cell = BOX / side
    cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    return cell, np.array(center) - cell * (cx * RIGHT + cy * UP), side


def numbers():
    med = {(k, N): f"{v:.2f}".rstrip("0").rstrip(".") for k, N, v in
           [(k, N, GF.median_distance(k, N, 400, seed=5)) for k, N, _ in GF.MEDIANS]}
    tg, tc = trails()
    return {"book_int": f"{GF.radial_integral(GF.p_gauss_book):.3f}", "chi1": f"{GF.chi_cauchy(1.0):.4f}",
            "med": med, "thg": {N: f"{GF.CHI3_MEDIAN * math.sqrt(N):.2f}" for N in (1, 10, 100)},
            "thc": {N: f"{GF.cauchy_radial_median() * N:.2f}" for N in (1, 10, 100)},
            "cg": f"{GF.CHI3_MEDIAN:.3f}", "cc": f"{GF.cauchy_radial_median():.3f}",
            "side_g": f"{fit_box(tg, ORIGIN)[2]:.0f}", "side_c": f"{fit_box(tc, ORIGIN)[2]:.0f}"}


def check_math():
    nb = numbers()
    assert nb["book_int"] == "12.566" and abs(GF.radial_integral(GF.p_gauss_book) - 4 * math.pi) < 1e-9
    assert abs(GF.radial_integral(GF.p_gauss) - 1) < 1e-9
    assert nb["chi1"] == "0.3679" and abs(GF.chi_cauchy_numeric(1.0) - GF.chi_cauchy(1.0)) < 2e-5
    assert nb["med"] == {("gauss", 1): "1.45", ("gauss", 10): "4.92", ("gauss", 100): "15.79",
                         ("cauchy", 1): "2.5", ("cauchy", 10): "23.28", ("cauchy", 100): "236.06"}
    assert {k: float(v) for k, v in nb["med"].items()} == {(k, N): v for k, N, v in GF.MEDIANS}
    assert nb["thg"] == {1: "1.54", 10: "4.86", 100: "15.38"} and nb["thc"] == {1: "2.26", 10: "22.64", 100: "226.44"}
    assert (nb["cg"], nb["cc"]) == ("1.538", "2.264")
    for N, r in ((1, 0.5), (10, 7.0)):                                    # (2.31)：形狀和一步一樣，寬度乘 N
        assert abs(GF.P_cauchy_N(r, N) - GF.p_cauchy(r / N) / N ** 3) < 1e-15
    for q in (0.3, 1.0):                                                  # 高斯步長：χ = e^{−q²r0²/2}
        num = GF.radial_integral(lambda r, r0: GF.p_gauss(r, r0) * math.sin(q * r) / (q * r) if r else GF.p_gauss(0, r0))
        assert abs(num - math.exp(-q * q / 2)) < 1e-6
    for d in (1e-3, 1e-4):                                                # e^{−q r0} ≈ 1 − r0|q|：斜率 −1
        assert abs((1 - GF.chi_cauchy(d)) / d - 1) < 1e-3
    tg, tc = trails()
    assert len(tg) == len(tc) == 201 and int(nb["side_c"]) > 4 * int(nb["side_g"])


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def tagged(tex, tag, font_size=30):
    return VGroup(tx(tex, font_size), tx(r"\text{" + tag + "}", 22, GREY_B)).arrange(RIGHT, buff=0.3)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：如果偶爾會跳很遠很遠……", wait=2.0)
    s.say("走 ", "$N$", " 步之後離家多遠？還是 ", "$\\sqrt N$", " 嗎？", wait=3.0)


def part_moments(s):
    nb = numbers()
    e1 = tagged(r"\chi(\vec q)=\int e^{i\vec q\cdot\vec r}\,p(r)\,d^3r", "(2.23)")
    e2 = tagged(r"=\int\Big[1+i\vec q\cdot\vec r-\tfrac12(\vec q\cdot\vec r)^2+\cdots\Big]p(r)\,d^3r", "(2.24)")
    e3 = mixed("$\\int p\\,d^3r=1$", "（(2.25)）", "$\\Rightarrow\\chi(\\vec 0)=1$", "；各向同性：",
               "$\\chi\\approx1-\\tfrac16q^2\\langle r^2\\rangle$", font_size=28)
    col = VGroup(e1, e2, e3).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    fit(col, 12.8).move_to(np.array([0, 3.15 - col.get_height() / 2, 0]))
    g1 = mixed("練習 2.7：", "$p(r)=(2\\pi r_0^2)^{-3/2}e^{-r^2/2r_0^2}$", font_size=28)
    g2 = tx(r"\Rightarrow\ \chi=e^{-q^2r_0^2/2}\ \Rightarrow\ \chi^N=e^{-Nq^2r_0^2/2}", 30)
    g3 = mixed("$N$", " 步仍是高斯：每軸變異數 ", "$Nr_0^2$", font_size=28, color=GAUSS_C)
    ex = VGroup(g1, g2, g3).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    ex.move_to(np.array([-6.6 + ex.get_width() / 2, -1.45, 0]))
    box = D.erratum_box("p. 37，Ex. 2.7", r"\sqrt{2/\pi}\,r_0^{-3}", r"(2\pi)^{-3/2}r_0^{-3}", width=4.6)
    note = mixed("照印的 ", "$\\int p\\,d^3r=4\\pi\\approx" + nb["book_int"] + "$", font_size=24)
    err = VGroup(box, note).arrange(DOWN, buff=0.2)
    err.move_to(np.array([6.7 - err.get_width() / 2, -1.45, 0]))
    assert_in_frame(col, ex, err)
    assert_apart(col, ex)
    assert_apart(col, err)
    assert_apart(ex, err)
    s.say("離開格子：每一步的位移 ", "$\\vec r$", " 是連續的，機率密度 ", "$p(r)$", wait=0.3)
    s.play(Write(e1), run_time=1.5)
    s.say("結構函數變成積分：", "$p(r)$", " 的傅立葉轉換", wait=3.0)
    s.play(FadeIn(e2), run_time=1.0)
    s.say("把 ", "$e^{i\\vec q\\cdot\\vec r}$", " 展開：係數就是步長的各階動差", wait=3.0)
    s.say("（書上 (2.24) 左邊的 ", "$d\\vec q$", " 應為 ", "$d^3r$", "）", wait=2.5)
    s.play(FadeIn(e3), run_time=1.0)
    s.say("一次項平均掉；", "$q^2$", " 項的係數＝二階動差 ", "$\\langle r^2\\rangle$", wait=3.5)
    s.play(FadeIn(g1), run_time=0.8)
    s.say("例：高斯步長", wait=1.5)
    s.play(FadeIn(g2), run_time=1.0)
    s.play(FadeIn(g3), run_time=0.8)
    s.say("$\\chi^N$", " 還是高斯：走 ", "$N$", " 步的距離 ", "$\\propto\\sqrt N$", wait=3.0)
    s.play(FadeIn(err, shift=0.2 * LEFT), run_time=0.8)
    s.say("書上的歸一化常數印錯：照印的總機率是 ", "$4\\pi$", wait=3.5)
    s.clear_stage()


def part_cauchy(s):
    nb = numbers()
    c1 = tagged(r"p(r)=\frac{r_0}{\pi^2(r^2+r_0^2)^2}", "(2.26)", 32)
    c2 = mixed("$\\int p\\,d^3r=1$", "，但 ", "$\\int r^2p\\,d^3r=\\infty$", font_size=28)
    c3 = tagged(r"\chi(q)=\frac{4\pi}{q}\int_0^\infty r\sin(qr)\,p(r)\,dr", "(2.28)", 30)
    c4 = mixed("被積函數是偶函數：延伸到 ", "$-\\infty<r<\\infty$", "，", font_size=26)
    c5 = mixed("在上半平面繞 ", "$r=ir_0$", "（二階極點）", font_size=26)
    c6 = VGroup(tx(r"\chi(q)=e^{-qr_0}", 40).set_color(YELLOW), tx(r"\text{(2.30)}", 22, GREY_B)).arrange(RIGHT, buff=0.3)
    c7 = mixed("$q=r_0=1$", "：", "$e^{-1}=" + nb["chi1"] + "$", font_size=26)
    col = VGroup(c1, c2, c3, c4, c5, c6, c7).arrange(DOWN, aligned_edge=LEFT, buff=0.27)
    fit(col, 7.6).move_to(np.array([-6.6 + col.get_width() / 2, 0.1, 0]))
    ax = D.plot_axes((-2, 2, 1), (-1.5, 2.5, 1), width=4.0, height=4.0, x_label=r"\mathrm{Re}\,r",
                     y_label=r"\mathrm{Im}\,r", font_size=20)
    ax.move_to(np.array([4.3, 0.2, 0]))
    R = 1.8
    o = ax.c2p(0, 0)
    rad = np.linalg.norm(ax.c2p(R, 0) - o)
    path = VGroup(Line(ax.c2p(-R, 0), ax.c2p(R, 0)), Arc(start_angle=0, angle=PI, radius=rad).shift(o))
    path.set_stroke(CAUCHY_C, 3)
    up = D.pole_marker(ax, 1j, r"ir_0", color=YELLOW)
    dn = D.pole_marker(ax, -1j, "", color=GREY_B)
    dnl = tx(r"-ir_0", 24, GREY_B).next_to(ax.c2p(0, -1), RIGHT, buff=0.15)
    assert_in_frame(col, ax, dnl)
    assert_apart(col, ax)
    s.say("長尾的步長：", "$p(r)$", " 只按 ", "$r^{-4}$", " 衰減", wait=0.3)
    s.play(FadeIn(c1), run_time=1.0)
    s.play(FadeIn(c2), run_time=0.8)
    s.say("總機率是 1，二階動差卻是無窮大", wait=3.0)
    s.play(FadeIn(c3), run_time=1.0)
    s.say("先把角度積掉，剩下對 ", "$r$", " 的一維積分", wait=3.0)
    s.play(FadeIn(c4), FadeIn(ax), FadeIn(up), FadeIn(dn), FadeIn(dnl), run_time=1.0)
    s.say("$r^2+r_0^2=0$", "：兩個二階極點 ", "$r=\\pm ir_0$", wait=3.0)
    s.play(FadeIn(c5), ShowCreation(path), run_time=1.5)
    s.play(Indicate(up[0], color=YELLOW, scale_factor=1.6), run_time=1.0)
    s.say("圍道只圍住上面的 ", "$ir_0$", "：留數給出答案", wait=3.0)
    s.play(FadeIn(c6, shift=0.2 * UP), run_time=1.0)
    s.play(FadeIn(c7), run_time=0.6)
    s.say("$\\chi(q)=e^{-qr_0}$", "：", "$q=r_0=1$", " 時是 ", "$" + nb["chi1"] + "$", wait=3.0)
    s.clear_stage()


def part_cusp(s):
    ax = D.plot_axes((-3, 3, 1), (0, 1.2, 0.5), width=7.6, height=3.9, x_ticks=[-2, -1, 1, 2], y_ticks=[0.5, 1],
                     y_tick_labels=["0.5", "1"], font_size=22, x_label="q")
    ax.move_to(np.array([-2.4, 0.15, 0]))
    gg = ax.get_graph(lambda q: math.exp(-q * q / 2), x_range=[-3, 3, 0.02]).set_stroke(GAUSS_C, 4)
    gc = ax.get_graph(lambda q: GF.chi_cauchy(abs(q)), x_range=[-3, 3, 0.01]).set_stroke(CAUCHY_C, 4)
    tg = ax.get_graph(lambda q: 1 - q * q / 2, x_range=[-1.2, 1.2, 0.02])
    tc = ax.get_graph(lambda q: 1 - abs(q), x_range=[-0.85, 0.85, 0.01])
    tg = DashedVMobject(tg.set_stroke(GAUSS_C, 2), num_dashes=24)
    tc = DashedVMobject(tc.set_stroke(CAUCHY_C, 2), num_dashes=24)
    lg = mixed("高斯步長 ", "$e^{-q^2r_0^2/2}\\approx1-\\tfrac12r_0^2q^2$", font_size=26, color=GAUSS_C)
    lc = mixed("Cauchy 步長 ", "$e^{-qr_0}\\approx1-r_0|q|$", font_size=26, color=CAUCHY_C)
    leg = VGroup(lg, lc).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    leg.move_to(np.array([6.7 - leg.get_width() / 2, 1.3, 0]))
    key = mixed("從 ", "$|q|$", " 開始變，不是 ", "$q^2$", font_size=30, color=YELLOW)
    key.next_to(leg, DOWN, buff=0.6).align_to(leg, LEFT)
    assert_in_frame(ax, leg, key)
    assert_apart(ax, leg)
    s.play(FadeIn(ax), run_time=0.8)
    s.play(ShowCreation(gg), FadeIn(lg), run_time=1.5)
    s.say("高斯步長的 ", "$\\chi$", "：頂端是平的，從 ", "$q^2$", " 開始往下", wait=3.0)
    s.play(ShowCreation(gc), FadeIn(lc), run_time=1.5)
    s.say("Cauchy 步長：頂端是尖的！", wait=2.5)
    s.play(ShowCreation(tg), ShowCreation(tc), run_time=1.2)
    s.play(FadeIn(key), run_time=0.8)
    s.say("尖點：", "$\\chi\\approx1-r_0|q|$", "，一開始就是一次方", wait=3.5)
    s.say("W207 的普適性要求 ", "$q^2$", "：這裡它失效了", wait=3.0)
    s.clear_stage()


def part_medians(s):
    nb = numbers()
    f1 = tagged(r"\chi^N=e^{-Nqr_0}\ \Rightarrow\ P(N;\vec r)=\frac{Nr_0}{\pi^2(r^2+N^2r_0^2)^2}", "(2.31)", 32)
    fit(f1, 13.0).move_to(np.array([0, 2.6, 0]))
    f2 = mixed("和一步的 ", "$p(r)$", " 同一個形狀：", "$r_0$", " 換成 ", "$Nr_0$", "，寬度乘 ", "$N$", font_size=28)
    f2.next_to(f1, DOWN, buff=0.3)
    hdr = [tx("N", 28, YELLOW), zh("高斯步長", font_size=24, color=GAUSS_C), tx(nb["cg"] + r"\sqrt N", 28, GAUSS_C),
           zh("Cauchy 步長", font_size=24, color=CAUCHY_C), tx(nb["cc"] + r"\,N", 28, CAUCHY_C)]
    rows = [hdr] + [[str(N), nb["med"][("gauss", N)], nb["thg"][N], nb["med"][("cauchy", N)], nb["thc"][N]]
                    for N in (1, 10, 100)]
    tab = D.tex_table(rows, [1.4, 2.0, 1.9, 2.2, 1.9], row_h=0.55, font_size=28)
    tab.move_to(np.array([0, -0.75, 0]))
    sub = mixed("400 個走路人、離原點距離的中位數（", "$r_0=1$", "）", font_size=22, color=GREY_B)
    sub.next_to(tab, DOWN, buff=0.2)
    hide = D.table_cell(tab, 3, 3, 5)
    q = tx("?", 34, YELLOW).move_to(hide)
    assert_in_frame(f1, f2, tab, sub)
    assert_apart(f2, tab)
    s.say("$N$", " 步之後：", "$\\chi^N$", " 只是把 ", "$r_0$", " 換成 ", "$Nr_0$", wait=0.3)
    s.play(FadeIn(f1), run_time=1.2)
    s.play(FadeIn(f2), run_time=0.8)
    s.say("距離正比於 ", "$N$", "，不是 ", "$\\sqrt N$", wait=3.0)
    others = VGroup(*[c for c in tab[0] if c is not hide])
    s.play(FadeIn(others), ShowCreation(tab[1]), FadeIn(sub), FadeIn(q), run_time=1.2)
    s.say("模擬：每種步長各 400 個人，看走多遠", wait=3.0)
    s.say("高斯步長：", "$N$", " 乘 100，距離只乘 10", wait=3.0)
    s.say("猜猜看：Cauchy 步長走 100 步，中位距離是多少？", wait=3.0)
    s.play(FadeOut(q), FadeIn(hide, scale=1.3), run_time=0.8)
    s.say("$" + nb["med"][("cauchy", 100)] + "$", "，不是 15：", "$N$", " 乘 10，距離也乘 10", wait=3.5)
    s.clear_stage()


def part_trails(s):
    nb = numbers()
    tg, tc = trails()
    centers = (np.array([-3.4, 0.05, 0]), np.array([3.4, 0.05, 0]))
    mobs = VGroup()
    for pts, ctr, col in ((tg, centers[0], GAUSS_C), (tc, centers[1], CAUCHY_C)):
        cell, org, _ = fit_box(pts, ctr)
        frame = Square(side_length=BOX + 0.3).move_to(ctr).set_stroke(GREY_C, 1.5)
        mobs.add(VGroup(frame, L.trail([(p[0], p[1]) for p in pts], cell=cell, origin=org, color=col, stroke_width=2.5)))
    t1 = mixed("高斯步長：方框邊長 ", "$\\approx" + nb["side_g"] + "\\,r_0$", font_size=26, color=GAUSS_C)
    t2 = mixed("Cauchy 步長：方框邊長 ", "$\\approx" + nb["side_c"] + "\\,r_0$", font_size=26, color=CAUCHY_C)
    t1.next_to(mobs[0], UP, buff=0.15)
    t2.next_to(mobs[1], UP, buff=0.15)
    assert_in_frame(mobs, t1, t2)
    assert_apart(mobs[0], mobs[1])
    s.say("各走 200 步，看足跡（只畫 ", "$x$", "、", "$y$", " 兩個方向）", wait=0.3)
    s.play(FadeIn(mobs[0][0]), FadeIn(t1), run_time=0.6)
    s.play(ShowCreation(mobs[0][1][0]), run_time=3.0)
    s.play(FadeIn(mobs[0][1][1:]), run_time=0.4)
    s.say("高斯步長：亂走成一團", wait=2.0)
    s.play(FadeIn(mobs[1][0]), FadeIn(t2), run_time=0.6)
    s.play(ShowCreation(mobs[1][1][0]), run_time=3.0)
    s.play(FadeIn(mobs[1][1][1:]), run_time=0.4)
    s.say("Cauchy 步長：一團一團，中間幾次大跳躍", wait=3.0)
    s.say("這叫 Lévy 飛行：方框大了好幾倍，全靠那幾次跳躍", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("步長的二階動差發散：", "$\\chi(q)$", " 從 ", "$|q|$", " 開始變", wait=3.0)
    s.say("普適性失效：距離 ", "$\\propto N$", "，而不是 ", "$\\sqrt N$", wait=3.0)


class W208OffLattice(ChuScene):
    def construct(self):
        self.show_header("W208　離開格子：長尾步長讓 √N 失效")
        check_math()
        part_question(self)
        part_moments(self)
        part_cauchy(self)
        part_cusp(self)
        part_medians(self)
        part_trails(self)
        part_close(self)
        self.wait(1.0)
