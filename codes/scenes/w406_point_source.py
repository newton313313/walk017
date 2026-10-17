"""
W406 點源的 1/r 與吸收牆的鏡像源（§4.2.2–4.2.3，pp. 86–89，(4.50)–(4.60)，Fig. 4.4–4.5）

問：在原點不斷放出走路人，穩定以後，離源 r 處的濃度怎麼隨 r 變？像高斯那樣急速下降嗎？
  1. (4.53)–(4.56)：高斯對時間積分 → 1/(4πDr)；D = 1、r = 2：0.03979 = 1/(8π)（BC.time_integral）。
  2. 格子（sc，每步在原點放 1 人）：w = G(z_c; r)（BC.sc_green）；nG → 3/(2π) = 0.4775（更正後的 (4.41)）；書上 3/(4π)。
  3. 吸收牆旁的點源 (4.57)(4.58)：鏡像源 → 偶極 2lx/r³；熱圖 BC.conc_wall；筆誤框 x-z plane → y-z plane。
  4. (4.59)(4.60)：穿過遠處平面的總流量 0（BC.dipole_flux_integral、BC.plane_flux）；牆與源之間 −1。
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
from utils import lattice_draw as L
from utils import walk_draw as D

GAUSS_TS = (0.25, 0.5, 1.0, 2.0, 4.0)
CUM_TS = (0.5, 1.0, 2.0, 4.0, 16.0)
GAUSS_COLORS = [BLUE_B, BLUE_C, TEAL_C, GREEN_C, GREY_B]
NS = (0, 1, 2, 3, 4, 5, 10)
LSRC = 1.0                                   # 源在 (l, 0, 0)
HM_X, HM_Y, HM_N = (-2.0, 4.0), (-3.0, 3.0), 30
HM_W, HM_C = 5.2, np.array([-3.7, 0.0, 0])   # 熱圖寬（正方形）與中心
HM_S = HM_W / (HM_X[1] - HM_X[0])
SRC_COLOR, IMG_COLOR = ORANGE, BLUE_C
TRAIL_SEEDS = (10, 11, 13, 15, 16, 17)   # 示意路線（都在 ±19 格內）


# ---------------------------------------------------------------- 數學
def cum_conc(r, T, D_=1.0):
    """∫_0^T (4πDt)^{−3/2} e^{−r²/4Dt} dt = erfc(r/√(4DT))/(4πDr)（只到時間 T 為止放出的人）。"""
    return math.erfc(r / math.sqrt(4 * D_ * T)) / (4 * math.pi * D_ * r)


def cum_conc_numeric(r, T, D_=1.0, n=20000):
    h = T / n
    return sum(BC.heat_kernel(r, (k + 0.5) * h, D_) for k in range(n)) * h


def walk2d(steps, seed, turn=0):
    """平面格子上的示意路線（整條轉 turn 個直角，讓幾條路線往不同方向散開）。"""
    rng = np.random.default_rng(seed)
    moves = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    p, pts = (0, 0), [(0, 0)]
    for k in rng.integers(0, 4, steps):
        dx, dy = moves[(k + turn) % 4]
        p = (p[0] + dx, p[1] + dy)
        pts.append(p)
    return pts


@lru_cache(maxsize=None)
def numbers():
    nb = {"ti": BC.time_integral(2.0), "G": {n: BC.sc_green(n, 0, 0) for n in NS}}
    nb["nG"] = {n: n * nb["G"][n] for n in NS if n > 0}
    nb["asym"], nb["book"] = BC.sc_green_asym(1), BC.sc_green_book(1)
    nb["flux"] = {X: BC.plane_flux(X, LSRC) for X in (3.0, 0.5)}
    nb["dip"] = [BC.dipole_flux_integral(x) for x in (0.5, 1.0, 2.0, 3.0)]
    return nb


def check_math():
    nb = numbers()
    assert f"{nb['ti']:.5f}" == "0.03979" == f"{1 / (8 * math.pi):.5f}"
    assert abs(nb["ti"] - 1 / (8 * math.pi)) < 1e-8
    assert abs(BC.K_const_numeric(4.0) - BC.K_const(4.0)) < 1e-8 and abs(BC.K_const(4.0) - math.sqrt(4 * math.pi)) < 1e-12
    assert abs(cum_conc(2.0, 1e9) - nb["ti"]) < 1e-5 and abs(cum_conc(2.0, 4.0) - cum_conc_numeric(2.0, 4.0)) < 1e-7
    assert all(cum_conc(1.0, CUM_TS[k]) < cum_conc(1.0, CUM_TS[k + 1]) for k in range(len(CUM_TS) - 1))
    assert max(BC.heat_kernel(0.0, t, 1.0) for t in GAUSS_TS) < 0.2 and 1 / (4 * math.pi * 0.4) < 0.2
    G = nb["G"]
    assert [f"{G[n]:.4f}" for n in NS[:6]] == ["1.5164", "0.5164", "0.2573", "0.1653", "0.1217", "0.0966"]
    assert f"{G[10]:.4f}" == "0.0479"
    assert [f"{nb['nG'][n]:.4f}" for n in NS[1:]] == ["0.5164", "0.5147", "0.4958", "0.4869", "0.4830", "0.4787"]
    assert abs(G[1] - (G[0] - 1)) < 1e-9                                                 # 原點的方程
    assert f"{nb['asym']:.4f}" == "0.4775" and abs(nb["asym"] - 3 / (2 * math.pi)) < 1e-12
    assert f"{nb['book']:.4f}" == "0.2387" and abs(nb["asym"] / nb["book"] - 2) < 1e-12
    assert abs(nb["nG"][10] - nb["asym"]) < 0.002 and nb["nG"][10] > nb["asym"]
    assert all(abs(v) < 1e-6 for v in nb["dip"])                                         # (4.60)
    assert abs(nb["flux"][3.0]) < 1e-6 and abs(nb["flux"][0.5] + 1) < 1e-6
    assert BC.plane_flux_exact(3.0, LSRC) == 0 and BC.plane_flux_exact(0.5, LSRC) == -1
    r = math.hypot(10.0, 3.0)
    assert abs(BC.conc_wall(10.0, 3.0, 0, LSRC) / BC.conc_dipole(10.0, 3.0, 0, LSRC) - 1) < 0.01   # 遠場＝偶極
    assert abs(BC.conc_dipole(10.0, 3.0, 0, LSRC) - 2 * LSRC * 10.0 / r ** 3 / (4 * math.pi)) < 1e-12
    assert all(np.abs(np.array(walk2d(160, sd, k))).max() <= 19 for k, sd in enumerate(TRAIL_SEEDS))
    assert abs(BC.conc_wall(0.0, 1.3, 0, LSRC)) < 1e-15 and BC.conc_wall(-0.5, 0.2, 0, LSRC) < 0 < BC.conc_wall(0.5, 0.2, 0, LSRC)


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


def hm_pt(x, y):
    xc = (HM_X[0] + HM_X[1]) / 2
    return HM_C + (x - xc) * HM_S * RIGHT + y * HM_S * UP


def charge(point, color, sign):
    d = Dot(point, radius=0.15).set_fill(color, 1).set_stroke(WHITE, 2)
    return VGroup(d, tx(sign, 28, BLACK).move_to(point))


def wall_scene():
    """熱圖（z = 0 剖面）＋牆＋源＋鏡像源＋遮住牆後的半透明板 → dict。"""
    hm = L.heatmap(lambda x, y: BC.conc_wall(x, y, 0.0, LSRC), x_range=HM_X, y_range=HM_Y, n=HM_N,
                   width=HM_W, origin=HM_C, vmax=0.1)
    frame = Rectangle(width=HM_W, height=HM_W).move_to(HM_C).set_stroke(GREY_C, 1.5)
    wall = Line(hm_pt(0, HM_Y[0]), hm_pt(0, HM_Y[1])).set_stroke(WHITE, 5)
    wall_lab = mixed("牆 ", "$x=0$", font_size=22).next_to(wall, UP, buff=0.08)
    src = charge(hm_pt(LSRC, 0), SRC_COLOR, "+")
    img = charge(hm_pt(-LSRC, 0), IMG_COLOR, "-")
    src_lab = zh("源", font_size=22).next_to(src, DOWN, buff=0.1)
    img_lab = zh("鏡像源", font_size=22).next_to(img, DOWN, buff=0.1)
    cw = (0 - HM_X[0]) * HM_S
    cover = Rectangle(width=cw, height=HM_W).set_stroke(width=0).set_fill(BLACK, 0.7)
    cover.move_to(hm_pt(HM_X[0] / 2, 0))
    cover_lab = zh("牆後（虛構）", font_size=20, color=GREY_A).move_to(hm_pt(HM_X[0] / 2, 2.4))
    return dict(hm=hm, frame=frame, wall=wall, wall_lab=wall_lab, src=src, img=img, src_lab=src_lab,
                img_lab=img_lab, cover=cover, cover_lab=cover_lab)


# ---------------------------------------------------------------- 各段
def part_question(s):
    org = np.array([-3.4, -0.1, 0])
    colors = [TEAL, BLUE_C, PURPLE_B, GOLD, PINK, GREEN_C]
    trails = VGroup(*[L.trail(walk2d(160, TRAIL_SEEDS[k], k), cell=0.11, origin=org, color=c, stroke_width=2, ends=False)
                      for k, c in enumerate(colors)])
    src = charge(org, SRC_COLOR, "+")
    opt1 = mixed("高斯那樣？ ", "$e^{-r^2/4Dt}$", font_size=30)
    opt2 = mixed("還是冪次？ ", "$r^{-?}$", font_size=30)
    opts = VGroup(opt1, opt2).arrange(DOWN, buff=0.7, aligned_edge=LEFT)
    at(opts, 3.4, 0.0)
    assert_in_frame(trails, src, opts)
    assert_apart(trails, opts)
    s.play(FadeIn(src, scale=0.5), run_time=0.6)
    s.say("問：在原點不斷放出走路人（三維；圖是示意）", wait=0.3)
    s.play(*[ShowCreation(t[0]) for t in trails], run_time=3.0)
    s.say("穩定以後，離源 ", "$r$", " 處的濃度怎麼隨 ", "$r$", " 變？", wait=1.5)
    s.play(FadeIn(opt1, shift=0.2 * LEFT), run_time=0.6)
    s.play(FadeIn(opt2, shift=0.2 * LEFT), run_time=0.6)
    s.say("像高斯那樣急速下降嗎？（先猜）", wait=3.0)
    s.clear_stage()


def part_time_integral(s):
    nb = numbers()
    top = at(eq_label(tx(r"c(r)=\int_0^\infty\frac{e^{-r^2/4Dt}}{(4\pi Dt)^{3/2}}\,dt=\frac{1}{4\pi Dr}", 34),
                      "(4.53)--(4.56)"), 0, 2.45)
    ax = D.plot_axes((0, 5, 1), (0, 0.2, 0.05), width=6.0, height=3.5, x_ticks=[1, 2, 3, 4, 5],
                     y_ticks=[0.05, 0.1, 0.15, 0.2], font_size=20, x_label="r", y_label="c")
    ax.move_to(np.array([-3.4, -0.75, 0]))
    gauss = VGroup(*[ax.get_graph(lambda r, t=t: BC.heat_kernel(r, t, 1.0), x_range=(0, 5, 0.05))
                     .set_stroke(GAUSS_COLORS[k], 2.5) for k, t in enumerate(GAUSS_TS)])
    glab = zh("不同時刻 t 的高斯", font_size=22, color=BLUE_B).move_to(ax.c2p(2.9, 0.17))
    cum = [ax.get_graph(lambda r, T=T: cum_conc(r, T), x_range=(0.4, 5, 0.05)).set_stroke(YELLOW, 4) for T in CUM_TS]
    lim = ax.get_graph(lambda r: 1 / (4 * math.pi * r), x_range=(0.4, 5, 0.05)).set_stroke(YELLOW, 4)
    tlab = [at(tx(r"\int_0^{%g}dt" % T, 26, YELLOW), 0, 0) for T in CUM_TS] + [at(tx(r"\int_0^{\infty}dt", 26, YELLOW), 0, 0)]
    for m in tlab:
        m.move_to(ax.c2p(3.8, 0.11))
    dot = Dot(ax.c2p(2.0, nb["ti"]), radius=0.07).set_color(RED)
    dlab = tx(r"%.5f" % nb["ti"], 24, RED).next_to(dot, UP, buff=0.12)
    r1 = tx(r"t=Tr^2\ \Rightarrow\ c(r)=\frac{K}{(4\pi D)^{3/2}\,r}", 30)
    r2 = tx(r"K=\int_0^\infty T^{-3/2}e^{-1/aT}dT=\sqrt{\pi a},\ \ a=4D", 30)
    r3 = mixed("$D=1,\\ r=2$", "：數值積分", font_size=26)
    r4 = tx(r"=%.5f=\frac{1}{8\pi}" % nb["ti"], 32, RED)
    right = VGroup(r1, r2, r3, r4).arrange(DOWN, buff=0.45)
    fit(right, 6.0)
    at(right, 3.55, -0.45)
    assert_in_frame(top, ax, right, glab, dlab)
    assert_apart(top, ax)
    assert_apart(top, right)
    assert_apart(ax, right)
    for k in range(3):
        assert_apart(right[k], right[k + 1])
    s.play(Write(top), run_time=1.5)
    s.say("一次放出的人：", "$t$", " 時刻是三維高斯", wait=2.0)
    s.say("一直放：把所有時刻的高斯加起來＝對 ", "$t$", " 積分", wait=2.5)
    s.play(ShowCreation(ax), run_time=1.0)
    s.play(LaggedStart(*[ShowCreation(g) for g in gauss], lag_ratio=0.3), FadeIn(glab), run_time=2.0)
    s.say("早的高斯又高又窄，晚的又矮又寬", wait=2.5)
    s.play(ShowCreation(cum[0]), FadeIn(tlab[0]), run_time=1.0)
    s.say("黃線：到時刻 ", "$T$", " 為止放出的人，在 ", "$r$", " 處的總和", wait=0.5)
    cur, curl = cum[0], tlab[0]
    for k in range(1, len(CUM_TS)):
        s.play(Transform(cur, cum[k]), Transform(curl, tlab[k]), run_time=0.9)
    s.play(Transform(cur, lim), Transform(curl, tlab[-1]), run_time=1.0)
    s.say("遠處的濃度來自很晚才到的人：所以降得很慢", wait=3.0)
    s.play(FadeIn(r1), run_time=0.8)
    s.say("代換 ", "$t=Tr^2$", "：", "$r$", " 整個提出來，剩下 ", "$1/r$", wait=3.0)
    s.play(FadeIn(r2), run_time=0.8)
    s.say("剩下的積分 ", "$K$", " 是個常數：", "$\\sqrt{\\pi a}$", wait=2.5)
    s.play(FadeIn(r3), FadeIn(r4), FadeIn(dot), FadeIn(dlab), run_time=0.8)
    s.say("數值核對：", "$r=2$", " 時 0.03979，正是 ", "$\\frac{1}{8\\pi}$", wait=3.0)
    s.clear_stage()


def part_lattice(s):
    nb = numbers()
    G, nG = nb["G"], nb["nG"]
    head = at(mixed("簡單立方格子，每步在原點放 1 人：", "$w(\\vec r)=G(z_c;\\vec r)$", "（第 3 章）", font_size=28), 0, 2.6)
    rows = [["n"] + [str(n) for n in NS],
            ["G(n,0,0)"] + [f"{G[n]:.4f}" for n in NS],
            ["n\\,G"] + ["-"] + [f"{nG[n]:.4f}" for n in NS[1:]]]
    tab = D.tex_table(rows, [2.0] + [1.3] * len(NS), row_h=0.6, font_size=28)
    at(tab, 0, 1.0)
    orig = at(mixed("原點：", "$G(\\vec 0)=\\tfrac16\\cdot6\\,G(1)+1\\ \\Rightarrow\\ G(1)=G(0)-1$", font_size=26), 0, -0.55)
    lim = at(VGroup(tx(r"nG\to\frac{6}{4\pi}=\frac{3}{2\pi}=%.4f" % nb["asym"], 34, YELLOW),
                    zh("更正後的 (4.41)", font_size=24, color=GREY_A)).arrange(RIGHT, buff=0.4), 0, -1.55)
    book = at(mixed("照書上的 ", "$3/a^2$", "：", "$\\frac{3}{4\\pi}=%.4f$" % nb["book"], "（差 2 倍）",
                    font_size=26, color=GREY_B), 0, -2.5)
    assert_in_frame(head, tab, orig, lim, book)
    assert_table(tab, len(NS) + 1)
    for a, b in ((head, tab), (tab, orig), (orig, lim), (lim, book)):
        assert_apart(a, b)
    s.play(FadeIn(head), run_time=0.8)
    s.say("格子上：穩定態就是第 3 章的 Green 函數", wait=2.5)
    s.play(FadeIn(tab[0][:len(NS) + 1]), ShowCreation(tab[1]), run_time=0.8)
    s.play(FadeIn(tab[0][len(NS) + 1:2 * (len(NS) + 1)]), run_time=0.8)
    s.say("$G(0)=1.5164$", "：原點平均站著 1.5 個人", wait=2.5)
    s.play(FadeIn(orig), run_time=0.8)
    s.say("原點的方程：鄰居平均＋1，所以 ", "$G(1)=G(0)-1$", wait=3.0)
    s.play(FadeIn(tab[0][2 * (len(NS) + 1):]), run_time=0.8)
    s.say("乘上距離 ", "$n$", "：", "$nG$", " 幾乎不變，所以 ", "$G\\propto1/r$", wait=3.0)
    s.play(FadeIn(lim, shift=0.2 * UP), run_time=0.8)
    s.say("Poisson 方程 ", "$-\\nabla^2w=6\\delta$", " 的解：", "$\\frac{6}{4\\pi r}$", wait=3.0)
    s.play(FadeIn(book), run_time=0.8)
    s.say("照書上的 ", "$3/a^2$", " 會少一半（呼應 W405 的筆誤）", wait=3.0)
    s.clear_stage()


def part_wall(s):
    w = wall_scene()
    q1 = eq_label(tx(r"c(\vec r)\propto\frac{1}{|\vec r-\vec r_s|}-\frac{1}{|\vec r-\vec r_i|}", 32), "(4.57)")
    q2 = tx(r"\vec r_s=(l,0,0),\ \ \vec r_i=(-l,0,0)", 28, GREY_A)
    q3 = VGroup(tx(r"\xrightarrow{\ r\gg l\ }\ \frac{2lx}{r^3}", 34, YELLOW), mixed("偶極，比 ", "$1/r$", " 快", font_size=24),
                tx(r"\text{(4.58)}", 22, GREY_B)).arrange(RIGHT, buff=0.3)
    box = D.erratum_box("p. 87，(4.58) 前", r"x\text{-}z\ \text{plane}", r"y\text{-}z\ \text{plane}\ (x=0)", width=6.0)
    right = VGroup(q1, q2, q3, box).arrange(DOWN, buff=0.5)
    fit(right, 6.4)
    at(right, 3.25, 0.1)
    left = VGroup(w["hm"], w["frame"], w["wall"], w["wall_lab"], w["src"], w["img"], w["src_lab"], w["img_lab"],
                  w["cover_lab"])
    assert_in_frame(left, right)
    assert_apart(left, right)
    for k in range(3):
        assert_apart(right[k], right[k + 1])
    assert_apart(w["wall_lab"], w["cover_lab"])
    s.play(FadeIn(w["frame"]), ShowCreation(w["wall"]), FadeIn(w["wall_lab"]), FadeIn(w["src"]), FadeIn(w["src_lab"]),
           run_time=1.0)
    s.say("源在 ", "$(l,0,0)$", "，吸收牆是 ", "$x=0$", " 平面：牆上濃度 0", wait=3.0)
    s.play(FadeIn(w["img"]), FadeIn(w["img_lab"]), run_time=0.8)
    s.say("像接地導體：在牆後放一個負的鏡像源", wait=2.5)
    s.play(FadeIn(q1), FadeIn(q2), run_time=1.0)
    s.bring_to_back(w["hm"])             # 熱圖放在牆與電荷的後面
    s.play(FadeIn(w["hm"]), run_time=1.5)
    s.say("熱圖（", "$z=0$", " 剖面）：橘是正、藍是負；牆上剛好抵銷", wait=3.0)
    s.play(FadeIn(w["cover"]), FadeIn(w["cover_lab"]), run_time=0.8)
    s.say("牆後是虛構的：真正的濃度只在 ", "$x>0$", wait=2.5)
    s.play(FadeIn(q3, shift=0.2 * DOWN), run_time=0.8)
    s.say("遠處一正一負抵銷成偶極：", "$2lx/r^3$", "，比 ", "$1/r$", " 快", wait=3.0)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上說牆是 ", "$x$", "-", "$z$", " 平面；源在 ", "$x$", " 軸上，牆要垂直 ", "$x$", " 軸", wait=3.5)
    s.play(FadeOut(right), FadeOut(w["wall_lab"]), run_time=0.7)
    return w


def part_flux(s, w):
    nb = numbers()
    planes = VGroup()
    plabs = VGroup()
    for X, c in ((0.5, TEAL), (3.0, YELLOW)):
        ln = DashedLine(hm_pt(X, HM_Y[0]), hm_pt(X, HM_Y[1])).set_stroke(c, 3)
        planes.add(ln)
        plabs.add(tx(r"X=%g" % X, 24, c).next_to(ln, UP, buff=0.08))
    p1 = eq_label(tx(r"\int_0^\infty\Big[\frac{3x^2}{(x^2+R^2)^{5/2}}-\frac{1}{(x^2+R^2)^{3/2}}\Big]R\,dR=0", 30),
                  "(4.59)(4.60)", buff=0.25)
    p2 = mixed("每個點源穿過無限大平面的流量：", "$\\pm\\frac12$", font_size=24)
    hdr = lambda t: zh(t, font_size=24, color=YELLOW)
    rows = [["X", hdr("源"), hdr("鏡像源"), hdr("合計")],
            ["3", r"+\tfrac12", r"-\tfrac12", "%d" % round(nb["flux"][3.0])],
            ["0.5", r"-\tfrac12", r"-\tfrac12", "%d" % round(nb["flux"][0.5])]]
    tab = D.tex_table(rows, [1.3, 1.5, 1.7, 1.5], row_h=0.62, font_size=30)
    p3 = mixed("數值積分：", "$X=3$", " 得 ", "$0$", "；", "$X=0.5$", " 得 ", "$-1$", font_size=24)
    right = VGroup(p1, p2, tab, p3).arrange(DOWN, buff=0.45)
    fit(right, 6.4)
    at(right, 3.25, 0.1)
    assert_in_frame(planes, plabs, right)
    assert_apart(VGroup(w["hm"], planes, plabs), right)
    assert_apart(plabs[0], plabs[1])
    assert_table(tab, 4)
    for k in range(3):
        assert_apart(right[k], right[k + 1])
    s.play(FadeIn(p1), run_time=1.0)
    s.say("遠處平面上偶極的流量：正負兩部分積分剛好是 0", wait=3.0)
    s.play(ShowCreation(planes[1]), FadeIn(plabs[1]), FadeIn(p2), run_time=0.8)
    s.say("精確說：每個點源有一半穿過任何無限大平面", wait=2.5)
    s.play(FadeIn(tab[0][:8]), ShowCreation(tab[1]), run_time=0.8)
    s.say("$X=3$", "：源 ", "$+\\frac12$", "、鏡像 ", "$-\\frac12$", "，合計 0", wait=3.0)
    s.play(ShowCreation(planes[0]), FadeIn(plabs[0]), FadeIn(tab[0][8:]), run_time=0.8)
    s.say("牆和源之間 ", "$X=0.5$", "：合計 ", "$-1$", "，全部流進牆", wait=3.0)
    s.play(FadeIn(p3), run_time=0.6)
    s.say("對照 W401：", "$x$", " 座標本身是一維漫步，終究碰到牆", wait=3.0)


def part_close(s):
    s.say("持續的點源給 ", "$1/r$", "，就像電荷", wait=2.5)
    s.say("一面吸收牆把它變成偶極：最後全被牆吃掉", wait=3.0)


class W406PointSource(ChuScene):
    def construct(self):
        self.show_header("W406　點源的 1/r 與吸收牆的鏡像源")
        check_math()
        part_question(self)
        part_time_integral(self)
        part_lattice(self)
        w = part_wall(self)
        part_flux(self, w)
        part_close(self)
        self.wait(1.0)
