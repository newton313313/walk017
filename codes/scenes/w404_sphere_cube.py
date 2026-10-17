"""
W404 關在球裡與立方體裡（§4.1.4，pp. 79–82，(4.28)–(4.37)，Ex. 4.5）

問：同樣的體積，關在球裡和關在立方體裡，哪一個活得久？
  1. 球的模態 (4.28)–(4.33)：j0(x) = sin x / x，牆上 j0(kR) = 0 → k_nR = π, 2π, 3π（BC.j0_zero）；角度積分只留 l = 0。
  2. 權重 (4.34) 8R³/(πn²)：R = 1 時 2.546、0.637（BC.sphere_weight = 數值積分）；短時間 (4.36) 加總 = (4/3)πR³ = 4.189。
     筆誤框 p. 80 (4.31)：照書上算 n = 1 的權重是 0.516。
  3. 長時間 (4.37)：6^N (8R³/π) e^{−π²N/6R²}；筆誤框 p. 81 (4.35)(4.37)；格子球 R = 8（2103 點）每步 0.0242
     vs 更正 π²/384 = 0.0257 vs 書上 π²/64 = 0.154。
  4. 立方體（Ex. 4.5）：短時間 6^N L³、長時間 6^N (512L³/π⁶) e^{−π²N/2L²}；格子 L = 6：125、2970 條（遞迴＝模態和）。
  5. 同體積：立方體衰減率除以球的 = 3(3/4π)^{2/3} = 1.155 → 球活得久（書外：Faber–Krahn）；前因子除以體積的表。
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

DEAD = D.WRONG_COLOR
SPH = BLUE_C
CUBE = ORANGE
RB = 8                                     # 格子球半徑
LC = 6                                     # 立方體邊長
NW = 8                                     # 權重長條數


@lru_cache(maxsize=None)
def numbers():
    nb = {"zeros": [BC.j0_zero(n) for n in (1, 2, 3)]}
    nb["w"] = [BC.sphere_weight(n, 1) for n in range(1, NW + 1)]
    nb["wnum"] = [BC.sphere_weight_numeric(n, 1) for n in (1, 2)]
    nb["wbook"] = BC.sphere_weight_book(1, 1)
    nb["vol"] = 4 / 3 * math.pi
    nb["ball"] = BC.ball_rate(RB)
    nb["sites"] = len(BC.ball_sites(RB))
    nb["rc"] = BC.sphere_rate(RB)
    nb["rb"] = BC.sphere_rate(RB, book=True)
    nb["cube"] = {N: BC.cube_total(N, LC) for N in (0, 2)}
    nb["ratio"] = BC.same_volume_rate_ratio()
    nb["pref"] = dict(BC.PREFACTOR)
    return nb


def check_math():
    nb = numbers()
    assert all(abs(z - n * math.pi) < 1e-9 for n, z in zip((1, 2, 3), nb["zeros"]))
    assert all(abs(BC.j0(z)) < 1e-9 for z in nb["zeros"]) and BC.j0(0.0) == 1
    w = nb["w"]
    assert [f"{v:.3f}" for v in w[:2]] == ["2.546", "0.637"]
    assert all(abs(a - b) < 1e-6 for a, b in zip(w[:2], nb["wnum"]))
    assert all(abs(w[n - 1] - 8 / (math.pi * n * n)) < 1e-12 for n in range(1, NW + 1))
    assert f"{nb['wbook']:.3f}" == "0.516"
    assert f"{nb['vol']:.3f}" == "4.189" and abs(8 / math.pi * math.pi ** 2 / 6 - nb["vol"]) < 1e-12
    assert abs(math.fsum(BC.sphere_weight(n, 1) for n in range(1, 100001)) - nb["vol"]) < 1e-4
    assert f"{nb['ball']:.4f}" == "0.0242" and nb["sites"] == 2103
    assert f"{nb['rc']:.4f}" == "0.0257" and abs(nb["rc"] - math.pi ** 2 / (6 * 64)) < 1e-15
    assert f"{nb['rb']:.3f}" == "0.154" and abs(nb["rb"] / nb["rc"] - 6) < 1e-12
    assert nb["cube"] == {0: 125, 2: 2970} and (LC - 1) ** 3 == 125
    assert all(abs(BC.cube_total_modes(N, LC) * 6 ** N - nb["cube"][N]) < 1e-6 for N in (0, 2))
    assert f"{nb['ratio']:.3f}" == "1.155" and abs(nb["ratio"] - 3 * (3 / (4 * math.pi)) ** (2 / 3)) < 1e-12
    p = nb["pref"]
    assert [f"{p[k]:.3f}" for k in ("interval", "cube", "sphere")] == ["0.811", "0.533", "0.608"]
    assert abs(p["cube"] - 512 / math.pi ** 6) < 1e-15 and abs(p["cube"] - p["interval"] ** 3) < 1e-15
    assert abs(p["sphere"] - 6 / math.pi ** 2) < 1e-15
    assert abs(p["sphere"] - BC.sphere_weight(1, 1) / nb["vol"]) < 1e-12
    assert abs(BC.cube_long(0, LC) / LC ** 3 - p["cube"]) < 1e-12


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


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


def ball_section(radius, center, rings=14):
    """球的剖面：最低模態 j0(πr/R) 的濃淡（同心圓）＋紅色吸收牆。"""
    c = np.array(center)
    shade = VGroup()
    for k in range(rings, 0, -1):
        r = radius * k / rings
        val = BC.j0(math.pi * (k - 0.5) / rings)
        shade.add(Circle(radius=r).move_to(c).set_stroke(width=0).set_fill(SPH, 0.85 * val))
    wall = Circle(radius=radius).move_to(c).set_stroke(DEAD, 5)
    return VGroup(shade, wall)


def same_volume_pair(R, y, xs=(-3.0, 3.0)):
    """同體積的球（剖面）與立方體（正面）。"""
    side = R * (4 * math.pi / 3) ** (1 / 3)
    ball = Circle(radius=R).move_to(np.array([xs[0], y, 0])).set_stroke(DEAD, 5).set_fill(SPH, 0.35)
    cube = Square(side_length=side).move_to(np.array([xs[1], y, 0])).set_stroke(DEAD, 5).set_fill(CUBE, 0.3)
    lb = zh("球", font_size=30, color=SPH).next_to(ball, DOWN, buff=0.2)
    lc = zh("立方體", font_size=30, color=CUBE).next_to(cube, DOWN, buff=0.2)
    dots = VGroup(*[D.walker_dot(radius=0.08).move_to(m.get_center()) for m in (ball, cube)])
    return VGroup(ball, cube, lb, lc, dots)


# ---------------------------------------------------------------- 各段
def part_question(s):
    pair = same_volume_pair(1.25, 0.7)
    q = at(mixed("同樣的體積：", "$\\frac43\\pi R^3=a^3$", font_size=30, color=GREY_A), 0, -1.7)
    assert_in_frame(pair, q)
    assert_apart(pair, q)
    s.play(FadeIn(pair[:4]), run_time=1.0)
    s.play(FadeIn(pair[4], scale=0.5), FadeIn(q), run_time=0.8)
    s.say("牆都是吸收牆：碰到就消失", wait=2.5)
    s.say("問：同樣的體積，關在球裡和立方體裡", wait=2.5)
    s.say("哪一個活得久？先猜", wait=3.0)
    s.clear_stage()


def part_modes(s):
    nb = numbers()
    sec = ball_section(1.6, (-4.6, 0.4, 0))
    rarrow = Arrow(sec[1].get_center(), sec[1].get_center() + 1.6 * RIGHT, buff=0, thickness=2).set_color(WHITE)
    rlab = tx("R", 28).next_to(rarrow, UP, buff=0.08)
    wlab = mixed("牆上 ", "$c=0$", font_size=24, color=DEAD).next_to(sec, DOWN, buff=0.2)
    ax = D.plot_axes((0, 10.5, 1), (-0.4, 1.15, 0.5), width=6.6, height=3.4, y_ticks=[1], x_label="x")
    ax.move_to(np.array([2.9, 0.05, 0]))
    pts = [ax.c2p(x, BC.j0(x)) for x in np.linspace(1e-6, 10.2, 200)]
    crv = VMobject().set_points_smoothly(pts).set_stroke(YELLOW, 4)
    zs = VGroup(*[Dot(ax.c2p(z, 0), radius=0.08).set_fill(DEAD, 1) for z in nb["zeros"]])
    zl = VGroup(*[tx(t, 26, DEAD).next_to(d, dr, buff=0.04)
                  for t, d, dr in zip((r"\pi", r"2\pi", r"3\pi"), zs, (UR, UL, UR))])
    f1 = tx(r"j_0(x)=\frac{\sin x}{x}", 34, YELLOW)
    f1.next_to(ax, UP, buff=0.25)
    f2 = eq_label(tx(r"j_0(k_nR)=0\ \Rightarrow\ k_nR=n\pi", 32), "(4.28)--(4.33)")
    f2.next_to(ax, DOWN, buff=0.3)
    assert_in_frame(sec, rarrow, rlab, wlab, ax, crv, zs, zl, f1, f2)
    assert_apart(sec, ax)
    assert_apart(wlab, f2)
    assert_apart(f1, ax)
    for lab in zl:
        assert lab.get_left()[0] > ax.get_left()[0]
    s.play(FadeIn(sec), GrowArrow(rarrow), FadeIn(rlab), run_time=1.2)
    s.play(FadeIn(wlab), run_time=0.6)
    s.say("球的剖面：顏色深淺是最低模態的濃度", wait=3.0)
    s.play(FadeIn(ax), FadeIn(f1), run_time=1.0)
    s.play(ShowCreation(crv), run_time=1.6)
    s.say("從中心放人：角度積分只留 ", "$l=0$", "，只看半徑 (4.30)", wait=3.5)
    s.play(FadeIn(zs), FadeIn(zl), FadeIn(f2), run_time=1.0)
    s.say("牆上要是 0：", "$k_nR=\\pi,\\ 2\\pi,\\ 3\\pi,\\dots$", wait=3.0)
    s.say("和線段的 ", "$\\sin(n\\pi x/L)$", " 一模一樣的條件", wait=2.5)
    s.clear_stage()


def part_weights(s):
    nb = numbers()
    w = nb["w"]
    e1 = eq_label(mixed("第 ", "$n$", " 個模態的權重 ", "$\\dfrac{8R^3}{\\pi n^2}$", font_size=30), "(4.34)")
    at(e1, 0, 2.5)
    ux, uy = 0.62, 1.15
    org = np.array([-6.1, -2.35, 0])
    bars = D.bar_chart({n: w[n - 1] for n in range(1, NW + 1)}, ux=ux, uy=uy, origin=org, color=SPH, width=0.8)
    labs = VGroup(*[tx(f"{w[k]:.3f}", 22).next_to(bars[0][k], UP, buff=0.08) for k in range(2)])
    nlab = VGroup(*[tx(str(n), 20, GREY_B).next_to(org + n * ux * RIGHT, DOWN, buff=0.1) for n in range(1, NW + 1)])
    cap = mixed("$R=1$", "，", "$n=1,2,\\dots$", font_size=22, color=GREY_A).next_to(nlab, DOWN, buff=0.08)
    sm = eq_label(tx(r"\sum_{n\ge1}\frac{8R^3}{\pi n^2}=\frac{8R^3}{\pi}\cdot\frac{\pi^2}{6}=\frac43\pi R^3", 32),
                  "(4.36)")
    at(sm, 3.15, 1.15)
    vol = mixed("$R=1$", "：", f"${nb['vol']:.3f}$", "＝球的體積", font_size=28, color=D.EXACT_COLOR)
    vol.next_to(sm, DOWN, buff=0.25)
    box = D.erratum_box("p. 80，(4.31)", r"\frac{2}{R}\,\frac{j_0\,j_0}{1-2z+z\Delta^2k^2}",
                        r"k_{n_0}^2\,\frac{j_0\,j_0}{1-6z+z\Delta^2k^2}", width=6.2)
    at(box, 3.15, -1.75)
    assert_in_frame(e1, bars, labs, nlab, cap, sm, vol, box)
    assert_apart(e1, labs)
    assert_apart(e1, sm)
    assert_apart(bars, sm)
    assert_apart(bars, box)
    assert_apart(labs, sm)
    assert_apart(sm, vol)
    assert_apart(vol, box)
    s.play(FadeIn(e1), run_time=1.0)
    s.say("從中心出發、終點不限：每個模態分到的權重", wait=3.0)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars[0]], lag_ratio=0.1), FadeIn(nlab), FadeIn(cap),
           run_time=1.6)
    s.play(FadeIn(labs), run_time=0.6)
    s.say("$R=1$", "：2.546、0.637……像 ", "$1/n^2$", " 那樣變小", wait=3.0)
    s.play(FadeIn(sm), run_time=1.0)
    s.play(FadeIn(vol), run_time=0.6)
    s.say("短時間全部加起來：", "$\\sum1/n^2=\\pi^2/6$", wait=3.0)
    s.say("得到 ", "$\\frac43\\pi R^3$", "：就是球的體積，像沒有牆", wait=3.0)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上的 (4.31)：分母該是 ", "$1-6z$", "，前因子該是 ", "$k_{n_0}^2$", wait=3.5)
    s.say("照書上算，", "$n=1$", " 的權重得 0.516，不是 2.546", wait=3.0)
    s.clear_stage()


def part_long(s):
    nb = numbers()
    e1 = eq_label(tx(r"(6-\Delta^2k^2)^N\approx6^N\,e^{-N\Delta^2k^2/6}", 32), "(4.37)")
    e2 = tx(r"\Rightarrow\ \ 6^N\,\frac{8R^3}{\pi}\,e^{-\pi^2N/6R^2}", 34, YELLOW)
    top = at(VGroup(e1, e2).arrange(RIGHT, buff=0.7), 0, 2.45)
    box = D.erratum_box("p. 81，(4.35)(4.37)", r"e^{-\Delta^2\pi^2N/R^2}", r"e^{-\Delta^2\pi^2n^2N/6R^2}", width=5.4)
    at(box, -3.6, 0.6)
    ttl = mixed("每步衰減率（", "$R=8$", "）", font_size=26, color=GREY_A)
    k = 17.0
    rows = [(mixed("格子球（", f"${nb['sites']}$", " 個格點）", font_size=22, color=SPH), nb["ball"], f"{nb['ball']:.4f}", SPH),
            (mixed("更正 ", "$\\pi^2/6R^2$", font_size=22, color=D.EXACT_COLOR), nb["rc"], f"{nb['rc']:.4f}", D.EXACT_COLOR),
            (mixed("書上 ", "$\\pi^2/R^2$", font_size=22, color=DEAD), nb["rb"], f"{nb['rb']:.3f}", DEAD)]
    chart = VGroup()
    for i, (lab, v, txt, c) in enumerate(rows):
        y = -0.55 - 0.75 * i
        lab.move_to(np.array([0, y, 0])).align_to(np.array([3.0, 0, 0]), RIGHT)
        bar = Rectangle(width=v * k, height=0.36).set_stroke(c, 1.5).set_fill(c, 0.7)
        bar.move_to(np.array([3.2 + v * k / 2, y, 0]))
        val = tx(txt, 24, c).next_to(bar, RIGHT, buff=0.1)
        chart.add(VGroup(lab, bar, val))
    ttl.next_to(chart, UP, buff=0.3).align_to(chart, LEFT)
    six = mixed("差 ", "$6$", " 倍", font_size=28, color=DEAD).next_to(chart, DOWN, buff=0.3).align_to(chart[2][1], RIGHT)
    assert_in_frame(top, box, ttl, chart, six)
    assert_apart(top, box)
    assert_apart(top, ttl)
    assert_apart(box, chart)
    assert_apart(box, ttl)
    for i in range(2):
        assert_apart(chart[i], chart[i + 1])
    s.play(FadeIn(e1), run_time=1.0)
    s.say("長時間：每一步，最低模態乘上 ", "$(6-\\Delta^2k^2)/6$", wait=3.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("$k=\\pi/R$", "：衰減率 ", "$\\pi^2/6R^2$", "——分母的 6 來自六個方向", wait=3.5)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (4.35)(4.37) 的指數漏了 ", "$n^2$", " 和 ", "$1/6$", wait=3.0)
    s.play(FadeIn(ttl), LaggedStart(*[FadeIn(c) for c in chart], lag_ratio=0.3), run_time=1.6)
    s.say("格子核對：半徑 8 的格子球，2103 個格點", wait=3.0)
    s.play(FadeIn(six), run_time=0.6)
    s.say("0.0242 接近更正的 0.0257；書上的大了 6 倍", wait=3.0)
    s.clear_stage()


def part_cube(s):
    nb = numbers()
    pic = L.lattice3d("sc", unit=1.1, origin=np.array([-4.5, 0.2, 0]), color=CUBE)
    plab = zh("Ex. 4.5：立方體（斜投影示意）", font_size=22, color=GREY_A).next_to(pic, DOWN, buff=0.35)
    l1 = mixed("三個方向各取奇數模態：", "$\\sin\\frac{n_x\\pi x}{L}\\sin\\frac{n_y\\pi y}{L}\\sin\\frac{n_z\\pi z}{L}$",
               font_size=26)
    l2 = mixed("短時間：", "$6^NL^3$", font_size=28)
    l3 = mixed("長時間：", "$6^N\\,\\frac{512L^3}{\\pi^6}\\,e^{-\\pi^2N/2L^2}$", font_size=28, color=YELLOW)
    l4 = mixed("衰減率 ", "$\\pi^2/2L^2$", "：和一維線段一樣", font_size=26, color=GREY_A)
    rows = [["N", zh("逐步遞迴", font_size=22, color=YELLOW), zh("模態和", font_size=22, color=YELLOW)]]
    rows += [[str(N), str(nb["cube"][N]), str(round(BC.cube_total_modes(N, LC) * 6 ** N))] for N in (0, 2)]
    tab = D.tex_table(rows, [1.0, 1.9, 1.9], row_h=0.48, font_size=26)
    thd = mixed("格子 ", f"$L={LC}$", "：活著的走法總數", font_size=24, color=GREY_A)
    col = VGroup(l1, l2, l3, l4, thd, tab).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    col.move_to(np.array([2.2, 0.15, 0]))
    assert_in_frame(pic, plab, col)
    assert_table(tab, 3)
    assert_apart(pic, col)
    assert_apart(plab, col)
    s.play(ShowCreation(pic[0]), FadeIn(pic[1:]), FadeIn(plab), run_time=1.2)
    s.say("立方體：三個方向各自是一條有兩面牆的線段", wait=3.0)
    s.play(FadeIn(l1), run_time=0.8)
    s.play(FadeIn(l2), FadeIn(l3), run_time=1.0)
    s.say("$\\frac{512}{\\pi^6}=\\big(\\frac{8}{\\pi^2}\\big)^3$", "：三條線段的前因子相乘", wait=3.0)
    s.play(FadeIn(l4), run_time=0.6)
    s.say("每個方向只在 1/3 的步數裡動：衰減率和線段相同", wait=3.0)
    s.play(FadeIn(thd), FadeIn(tab), run_time=1.0)
    s.say("$L=6$", "：125 個內點；走 2 步還活著 2970 條，兩種算法一致", wait=3.5)
    s.clear_stage()


def part_reveal(s):
    nb = numbers()
    pair = same_volume_pair(0.85, 1.75, xs=(-4.6, -1.9))
    pair.remove(pair[4])
    r1 = mixed("同體積：立方體的衰減率除以球的", font_size=26)
    r2 = tx(r"=3\Big(\frac{3}{4\pi}\Big)^{2/3}=" + f"{nb['ratio']:.3f}", 34, ORANGE)
    r3 = mixed("$\\Rightarrow$", " 球活得久", font_size=30, color=SPH)
    rr = VGroup(r1, r2, r3).arrange(DOWN, aligned_edge=LEFT, buff=0.25)
    rr.move_to(np.array([3.3, 1.85, 0]))
    ext = zh("書外：同體積裡，球的最低模態最慢（Faber-Krahn 不等式，1923、1925）", font_size=22, color=GREY_B)
    at(ext, 0, 0.0)
    p = nb["pref"]
    rows = [[zh("形狀", font_size=24, color=YELLOW), zh("長時間前因子除以體積", font_size=24, color=YELLOW)],
            [zh("線段", font_size=24), rf"8/\pi^2={p['interval']:.3f}"],
            [zh("立方體", font_size=24, color=CUBE), rf"512/\pi^6={p['cube']:.3f}"],
            [zh("球", font_size=24, color=SPH), rf"6/\pi^2={p['sphere']:.3f}"]]
    tab = D.tex_table(rows, [1.8, 3.8], row_h=0.5, font_size=28)
    tab.move_to(np.array([0, -0.55 - tab.get_height() / 2, 0]))
    assert_in_frame(pair, rr, ext, tab)
    assert_table(tab, 2)
    assert_apart(pair, rr)
    assert_apart(pair, ext)
    assert_apart(rr, ext)
    assert_apart(ext, tab)
    s.play(FadeIn(pair), run_time=0.8)
    s.say("揭曉：同樣的體積比衰減率", wait=2.0)
    s.play(FadeIn(r1), FadeIn(r2), run_time=1.0)
    s.say("立方體的衰減率是球的 1.155 倍", wait=2.5)
    s.play(FadeIn(r3), Indicate(pair[0], color=SPH), run_time=1.0)
    s.say("球活得久：同體積裡球最「圓」，離牆最遠", wait=3.0)
    s.play(FadeIn(ext), run_time=0.6)
    s.say("書外：這是 Faber-Krahn 不等式", wait=2.5)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.0)
    s.say("前因子也跟形狀有關：線段、立方體、球各不相同", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("形狀只改兩件事：衰減率（最低模態）與前因子", wait=3.0)
    s.say("同樣的體積，球最耐關", wait=2.5)


class W404SphereCube(ChuScene):
    def construct(self):
        self.show_header("W404　關在球裡與立方體裡")
        check_math()
        part_question(self)
        part_modes(self)
        part_weights(self)
        part_long(self)
        part_cube(self)
        part_reveal(self)
        part_close(self)
        self.wait(1.0)
