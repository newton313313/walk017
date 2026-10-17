"""
W401 鏡像走路人（§4.1.1，pp. 69–73，(4.1)–(4.8)，Fig. 4.1）

問：醉漢從 x0 = 1 出發，0 是懸崖。走 4 步還活著的路線有幾條？
  1. 16 條小走法圖（BC.walks(4, 1)）：碰過 0 的 10 條標紅、沒碰過的 6 條標綠；終點 > 0 的 11 條裡 5 條碰過牆。
  2. 鏡射（Fig. 4.1）：1→0→1→2→3 第一次碰牆之前那段對 0 鏡射 → −1→0→1→2→3（BC.reflect_first）；
     5 對一一對應（BC.reflection_pairs）；(4.1) C'(N;x,x0) = C(N;x−x0) − C(N;x+x0)：2、3、1，共 6 = 11 − 5。
  3. 主例 x0 = 3：存活機率 vs x0√(2/πN)；三維 sc 吸收平面 vs √(6x0²/πN)（(4.5)）。
  4. (4.8)：牆吃掉的走法數 μ^N √(2N/3π)；一維被吃掉的比例 1 − P_N(x0)（N = 100）只在 x0 ≲ √N 時明顯。
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

DEAD = D.WRONG_COLOR
ALIVE = D.EXACT_COLOR
MIRROR = BLUE_C
MAIN = [1, 0, 1, 2, 3]
NS_1D = (10, 100, 1000)
NS_3D = (100, 1000)
N_BARS = 100


def C(N, y):
    """N 步、位移 y 的無限制走法數。"""
    if abs(y) > N or (N + y) % 2:
        return 0
    return math.comb(N, (N + y) // 2)


@lru_cache(maxsize=None)
def numbers():
    nb = {"walks": BC.walks(4, 1)}
    nb["surv"] = BC.survivors_brute(4, 1)
    nb["pairs"] = BC.reflection_pairs(4, 1)
    nb["p1"] = {N: BC.survival_prob(N, 3) for N in NS_1D}
    nb["a1"] = {N: BC.survival_asym(N, 3) for N in NS_1D}
    nb["p3"] = {N: BC.survival_3d(N, 3) for N in NS_3D}
    nb["a3"] = {N: BC.survival_3d_asym(N, 3) for N in NS_3D}
    nb["eaten"] = {x0: 1 - float(BC.survival_prob(N_BARS, x0)) for x0 in range(1, 31)}
    return nb


def check_math():
    nb = numbers()
    ws = nb["walks"]
    assert len(ws) == 16 and all(w[0] == 1 and len(w) == 5 for w in ws)
    touched = [w for w in ws if BC.touches_wall(w)]
    assert len(touched) == 10 and len(ws) - len(touched) == 6
    assert nb["surv"] == {1: 2, 3: 3, 5: 1} == BC.survivors_image(4, 1) == BC.survivors_dp(4, 1)
    assert sum(1 for w in ws if w[-1] > 0) == 11
    assert sum(1 for w in touched if w[-1] > 0) == 5 == len(nb["pairs"])
    assert BC.reflect_first(MAIN) == [-1, 0, 1, 2, 3]
    assert (MAIN, [-1, 0, 1, 2, 3]) in [(list(a), list(b)) for a, b in nb["pairs"]]
    for a, b in nb["pairs"]:
        assert b[0] == -1 and b[-1] == a[-1] > 0 and BC.reflect_first(list(a)) == list(b)
    rows = {x: (C(4, x - 1), C(4, x + 1)) for x in (1, 3, 5)}
    assert rows == {1: (6, 4), 3: (4, 1), 5: (1, 0)}
    assert {x: a - b for x, (a, b) in rows.items()} == nb["surv"] and sum(nb["surv"].values()) == 6 == 11 - 5
    assert str(nb["p1"][10]) == "21/32" and f"{float(nb['p1'][10]):.5f}" == "0.65625"
    assert [f"{nb['a1'][N]:.4f}" for N in (10, 100)] == ["0.7569", "0.2394"]
    assert f"{float(nb['p1'][100]):.4f}" == "0.2356"
    assert f"{float(nb['p1'][1000]):.5f}" == "0.07557" and f"{nb['a1'][1000]:.5f}" == "0.07569"
    assert abs(nb["a1"][100] - 3 * math.sqrt(2 / (math.pi * 100))) < 1e-12
    assert [f"{nb['p3'][N]:.4f}" for N in NS_3D] == ["0.3958", "0.1305"]
    assert [f"{nb['a3'][N]:.4f}" for N in NS_3D] == ["0.4146", "0.1311"]
    assert abs(nb["a3"][100] - math.sqrt(6 * 9 / (math.pi * 100))) < 1e-12
    assert abs(BC.wall_reduction_asym(100, 3) - math.sqrt(200 / (3 * math.pi))) < 1e-12
    e = nb["eaten"]
    assert e[1] > 0.9 and e[30] < 0.01 and all(e[k] > e[k + 1] for k in range(1, 30))
    assert math.isqrt(N_BARS) ** 2 == N_BARS


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


def mini_walk(pos, origin, color=D.WALK_COLOR, ux=0.4, uy=0.1):
    """小走法圖：座標軸（位置 0 那條橫軸＝紅色懸崖）＋路線 → VGroup(axes, path)。"""
    ax = D.walk_axes(4, 5, ux=ux, uy=uy, origin=origin, font_size=12)
    axes = VGroup(ax[0], ax[1])
    axes[0].set_color(DEAD)
    path = D.walk_path(pos, ux=ux, uy=uy, origin=origin, color=color, stroke_width=3, dot_radius=0.03)
    return VGroup(axes, path)


def corridor(xmax, x0, unit, walls, y=0.0):
    """數線 0..xmax（D.number_line 取右半）＋紅色懸崖＋走路人 → VGroup(線, 懸崖, 走路人)。"""
    nl = D.number_line(xmax, unit=unit, label_every=1, font_size=24)
    nl.shift((-xmax / 2 * unit) * RIGHT + y * UP)
    p = lambda x: D.line_point(nl, x)
    ticks = VGroup(*[nl[1][k] for k in range(xmax, 2 * xmax + 1)])
    nums = VGroup(*[nl[2][k] for k in range(xmax, 2 * xmax + 1)])
    base = Line(p(0), p(xmax) + (0.4 * unit if len(walls) == 1 else 0) * RIGHT).set_stroke(GREY_B, 2)
    cliffs = VGroup(*[Line(p(w) + 0.55 * UP, p(w) + 0.55 * DOWN).set_stroke(DEAD, 6) for w in walls])
    labs = VGroup(*[zh("懸崖", font_size=24, color=DEAD).next_to(c, UP, buff=0.1) for c in cliffs])
    walker = D.walker_dot().move_to(p(x0) + 0.35 * UP)
    return VGroup(VGroup(base, ticks, nums), VGroup(cliffs, labs), walker)


# ---------------------------------------------------------------- 各段
def part_question(s):
    pic = corridor(6, 1, 1.1, (0,), y=0.6)
    guess = at(zh("走 4 步還活著的路線有幾條？", font_size=34, color=YELLOW), 0, -1.4)
    assert_in_frame(pic, guess)
    assert_apart(pic, guess)
    s.play(FadeIn(pic[0]), FadeIn(pic[1]), run_time=1.0)
    s.play(FadeIn(pic[2], shift=0.3 * DOWN), run_time=0.6)
    s.say("問：醉漢從 ", "$x_0=1$", " 出發，0 是懸崖", wait=2.5)
    s.say("踩到 0 就掉下去；每步左右各一半", wait=2.5)
    s.play(FadeIn(guess), run_time=0.6)
    s.say("走 4 步還活著的路線有幾條？先猜", wait=3.0)
    s.clear_stage()


def part_sixteen(s):
    nb = numbers()
    ws = nb["walks"]
    minis = VGroup()
    for i, w in enumerate(ws):
        r, c = divmod(i, 4)
        minis.add(mini_walk(w, np.array([-6.55 + 2.3 * c, 2.4 - 1.45 * r, 0])))
    touched = [BC.touches_wall(w) for w in ws]
    pos_end = [w[-1] > 0 for w in ws]
    surv = nb["surv"]
    t1 = mixed("活著（綠）：", "$6$", " 條", font_size=30, color=ALIVE)
    t2 = tx(rf"x=1:\ {surv[1]}\quad x=3:\ {surv[3]}\quad x=5:\ {surv[5]}", 24, ALIVE)
    t3 = mixed("碰過 0（紅）：", "$10$", " 條", font_size=30, color=DEAD)
    t4 = mixed("終點 ", "$x>0$", "：", "$11$", " 條", font_size=30)
    t5 = mixed("$11=6+5$", "（5 條碰過牆）", font_size=26, color=ORANGE)
    panel = VGroup(t1, t2, t3, t4, t5).arrange(DOWN, aligned_edge=LEFT, buff=0.42)
    panel.move_to(np.array([4.65, 0.3, 0]))
    marks = VGroup(*[SurroundingRectangle(minis[i], buff=0.04).set_stroke(ORANGE, 2)
                     for i in range(16) if touched[i] and pos_end[i]])
    assert len(marks) == 5
    assert_in_frame(minis, panel, marks)
    for i in range(16):
        assert_apart(minis[i], panel)
    for i in (0, 4, 8):
        assert_apart(minis[i], minis[i + 1], 0.1)
        assert_apart(minis[i], minis[i + 4], 0.1)
    s.play(LaggedStart(*[FadeIn(m[0]) for m in minis], lag_ratio=0.05), run_time=1.2)
    s.play(LaggedStart(*[ShowCreation(m[1][0]) for m in minis], lag_ratio=0.08),
           LaggedStart(*[FadeIn(m[1][1]) for m in minis], lag_ratio=0.08), run_time=2.4)
    s.say("從 1 出發走 4 步：全部 ", "$2^4=16$", " 條", wait=2.5)
    s.say("每張小圖的紅色橫線就是懸崖 0", wait=2.0)
    s.play(*[minis[i][1].animate.set_color(DEAD if touched[i] else ALIVE) for i in range(16)], run_time=1.2)
    s.play(FadeIn(t1), FadeIn(t3), run_time=0.8)
    s.say("碰過 0 的 10 條標紅；沒碰過的 6 條標綠", wait=3.0)
    s.play(FadeIn(t2), run_time=0.6)
    s.say("活著的 6 條：停在 1、3、5 的各有 2、3、1 條", wait=3.0)
    s.play(FadeIn(t4), run_time=0.6)
    s.say("換個數法：終點在 0 的右邊的有 11 條", wait=2.5)
    s.play(LaggedStart(*[ShowCreation(m) for m in marks], lag_ratio=0.15), FadeIn(t5), run_time=1.4)
    s.say("其中 5 條碰過牆（橘框）：", "$6=11-5$", wait=3.0)
    s.say("這 5 條要怎麼「直接」數出來？", wait=2.5)
    s.clear_stage()


def part_reflect(s):
    ux, uy = 1.3, 0.55
    org = np.array([-5.4, -0.3, 0])
    ax = D.walk_axes(4, 3, ux=ux, uy=uy, origin=org, font_size=22, x_every=1, y_every=1)
    ax[0].set_color(DEAD)
    orig = D.walk_path(MAIN, ux=ux, uy=uy, origin=org, color=D.WALK_COLOR, stroke_width=5, dot_radius=0.07)
    refl_pos = BC.reflect_first(MAIN)
    seg = D.walk_path(MAIN[:2], ux=ux, uy=uy, origin=org, color=D.WALK_COLOR, stroke_width=5, dot_radius=0.07)
    seg_r = D.walk_path(refl_pos[:2], ux=ux, uy=uy, origin=org, color=MIRROR, stroke_width=5, dot_radius=0.07)
    full_r = D.walk_path(refl_pos, ux=ux, uy=uy, origin=org, color=MIRROR, stroke_width=2.5, dot_radius=0.0)
    lab_o = tx(r"1\to0\to1\to2\to3", 30, D.WALK_COLOR)
    lab_r = tx(r"-1\to0\to1\to2\to3", 30, MIRROR)
    rule = zh("第一次碰牆之前那段對 0 鏡射", font_size=26, color=GREY_A)
    after = zh("之後兩條重合", font_size=26, color=GREY_A)
    panel = VGroup(lab_o, lab_r, rule, after).arrange(DOWN, aligned_edge=LEFT, buff=0.4)
    panel.move_to(np.array([4.4, 0.6, 0]))
    corr = mixed("碰過牆、終點 ", "$>0$", "（從 1 出發）", font_size=24, color=DEAD)
    corr2 = mixed("$\\Updownarrow$", "　一一對應", font_size=24)
    corr3 = mixed("從 ", "$-1$", " 出發、終點 ", "$>0$", font_size=24, color=MIRROR)
    cgrp = VGroup(corr, corr2, corr3).arrange(DOWN, aligned_edge=LEFT, buff=0.22)
    cgrp.next_to(panel, DOWN, buff=0.5).align_to(panel, LEFT)
    assert_in_frame(ax, orig, seg_r, panel, cgrp)
    assert_apart(ax, panel)
    assert_apart(ax, cgrp)
    assert_apart(panel, cgrp)
    s.play(FadeIn(ax), run_time=0.8)
    s.play(ShowCreation(orig[0]), FadeIn(orig[1]), FadeIn(lab_o), run_time=1.6)
    s.say("Fig. 4.1：這條路線第 1 步就踩到 0", wait=2.5)
    s.add(seg)
    s.play(FadeIn(rule), run_time=0.6)
    s.play(Transform(seg, seg_r), run_time=1.6)
    s.play(FadeIn(lab_r), ShowCreation(full_r[0]), run_time=1.2)
    s.say("碰牆之前的那段翻到下面：變成從 ", "$-1$", " 出發", wait=3.0)
    s.play(FadeIn(after), Indicate(orig[0], color=MIRROR), run_time=1.2)
    s.say("碰到 0 之後兩條路線完全一樣", wait=2.5)
    s.play(FadeIn(cgrp, shift=0.2 * UP), run_time=0.8)
    s.say("從 ", "$-1$", " 出發要走到右邊，一定會經過 0", wait=3.0)
    s.say("所以兩邊可以一條對一條地配對", wait=2.5)
    s.clear_stage()


def part_pairs(s):
    nb = numbers()
    pairs = nb["pairs"]
    top, bot, arrs = VGroup(), VGroup(), VGroup()
    for i, (a, b) in enumerate(pairs):
        x = -6.3 + 2.35 * i
        top.add(mini_walk(list(a), np.array([x, 2.45, 0]), color=DEAD))
        bot.add(mini_walk(list(b), np.array([x, 0.75, 0]), color=MIRROR))
        arrs.add(at(tx(r"\updownarrow", 30, GREY_A), x + 1.0, 1.6))
    lt = mixed("從 ", "$1$", " 出發", font_size=22, color=DEAD).next_to(top, RIGHT, buff=0.25)
    lb = mixed("從 ", "$-1$", " 出發", font_size=22, color=MIRROR).next_to(bot, RIGHT, buff=0.25)
    eq = at(eq_label(tx(r"C'(N;x,x_0)=C(N;x-x_0)-C(N;x+x_0)", 34), "(4.1)"), 0, -0.45)
    rows = [["x", "C(4;x-1)", "C(4;x+1)", "C'(4;x,1)"]]
    for x in (1, 3, 5):
        a, b = C(4, x - 1), C(4, x + 1)
        rows.append([str(x), str(a), str(b), f"{a}-{b}={a - b}"])
    tab = D.tex_table(rows, [1.2, 2.2, 2.2, 2.6], row_h=0.45, font_size=26)
    tab.move_to(np.array([-2.4, -1.0 - tab.get_height() / 2, 0]))
    tot = mixed("共 ", "$2+3+1=6=11-5$", font_size=28, color=ALIVE).next_to(tab, RIGHT, buff=0.5)
    assert_in_frame(top, bot, arrs, lt, lb, eq, tab, tot)
    assert_table(tab, 4)
    assert_apart(top, bot)
    assert_apart(bot, eq)
    assert_apart(eq, tab)
    assert_apart(tab, tot)
    for i in range(4):
        assert_apart(top[i], top[i + 1], 0.1)
    s.play(LaggedStart(*[FadeIn(m) for m in top], lag_ratio=0.15), FadeIn(lt), run_time=1.5)
    s.say("上排：碰過牆、終點在右邊的 5 條", wait=2.0)
    s.play(LaggedStart(*[TransformFromCopy(top[i], bot[i]) for i in range(5)], lag_ratio=0.15),
           FadeIn(arrs), FadeIn(lb), run_time=2.2)
    s.say("下排：鏡射後從 ", "$-1$", " 出發的 5 條", wait=2.5)
    s.say("從 ", "$-1$", " 走到 ", "$x>0$", " 的每一條都對得回來", wait=2.5)
    s.play(FadeIn(eq, shift=0.2 * UP), run_time=1.0)
    s.say("所以：活著的＝全部減去從鏡像點 ", "$-x_0$", " 出發的", wait=3.0)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.2)
    s.say("$C(N;y)$", "：無限制、位移 ", "$y$", " 的走法數", wait=2.5)
    s.play(FadeIn(tot), run_time=0.6)
    s.say("$6-4=2$", "、", "$4-1=3$", "、", "$1-0=1$", "：和暴力數的一樣", wait=3.0)
    s.clear_stage()


def part_survival(s):
    nb = numbers()
    f1 = tx(r"P_N(x_0)\approx x_0\sqrt{\frac{2}{\pi N}}", 32)
    f3 = eq_label(tx(r"P_N(x_0)\approx\sqrt{\frac{6x_0^2}{\pi N}}", 32), "(4.5)")
    h1 = mixed("一維，", "$x_0=3$", font_size=26, color=YELLOW)
    h3 = mixed("三維 sc，吸收平面 ", "$x=0$", "，", "$x_0=3$", font_size=26, color=YELLOW)
    p1, a1 = nb["p1"], nb["a1"]
    r1 = [["N", zh("精確", font_size=24, color=YELLOW), zh("近似", font_size=24, color=YELLOW)],
          ["10", f"21/32={float(p1[10]):.5f}", f"{a1[10]:.4f}"],
          ["100", f"{float(p1[100]):.4f}", f"{a1[100]:.4f}"],
          ["1000", f"{float(p1[1000]):.5f}", f"{a1[1000]:.5f}"]]
    r3 = [["N", zh("精確", font_size=24, color=YELLOW), zh("近似", font_size=24, color=YELLOW)]]
    r3 += [[str(N), f"{nb['p3'][N]:.4f}", f"{nb['a3'][N]:.4f}"] for N in NS_3D]
    t1 = D.tex_table(r1, [1.2, 3.0, 1.6], row_h=0.5, font_size=26)
    t3 = D.tex_table(r3, [1.2, 1.6, 1.6], row_h=0.5, font_size=26)
    left = VGroup(h1, f1, t1).arrange(DOWN, buff=0.3)
    right = VGroup(h3, f3, t3).arrange(DOWN, buff=0.3)
    at(left, -3.3, 0.9)
    right.move_to(np.array([3.4, 0, 0])).align_to(left, UP)
    note = at(zh("x 座標只在 1/3 的步數裡動：三維碰平面還是一維的問題", font_size=24, color=GREY_A), 0, -1.6)
    zero = at(mixed("$N\\to\\infty$", "：", "$P_N\\to0$", "，一維一定掉下去", font_size=26, color=DEAD), 0, -2.4)
    assert_in_frame(left, right, note, zero)
    assert_table(t1, 3)
    assert_table(t3, 3)
    assert_apart(left, right)
    assert_apart(left, note)
    assert_apart(right, note)
    assert_apart(note, zero)
    s.say("主例：從 ", "$x_0=3$", " 出發，走 ", "$N$", " 步還活著的機率", wait=0.5)
    s.play(FadeIn(h1), FadeIn(f1), run_time=1.0)
    s.play(FadeIn(t1[0]), ShowCreation(t1[1]), run_time=1.2)
    s.say("精確值用 (4.1) 加總；", "$N$", " 大時近似越來越準", wait=3.0)
    s.play(FadeIn(zero), run_time=0.6)
    s.say("存活機率像 ", "$1/\\sqrt N$", " 那樣趨近 0", wait=2.5)
    s.play(FadeIn(right), run_time=1.2)
    s.say("三維簡單立方、吸收平面：書上的 (4.5)", wait=2.5)
    s.play(FadeIn(note), run_time=0.6)
    s.say("前因子 ", "$\\sqrt3$", " 倍：只有 1/3 的步數往 ", "$x$", " 方向走", wait=3.0)
    s.clear_stage()


def part_eaten(s):
    nb = numbers()
    e = nb["eaten"]
    ux, uy = 0.3, 3.0
    org = np.array([-5.0, -2.2, 0])
    bars = D.bar_chart(e, ux=ux, uy=uy, origin=org, color=DEAD, width=0.75)
    base = Line(org + 0.4 * ux * RIGHT, org + 31 * ux * RIGHT).set_stroke(GREY_B, 2)
    xl = mixed("起點 ", "$x_0$", font_size=24, color=GREY_B).next_to(base, RIGHT, buff=0.15)
    ticks = VGroup(*[tx(str(k), 22, GREY_B).next_to(org + k * ux * RIGHT, DOWN, buff=0.12) for k in (1, 10, 20, 30)])
    xs = org + math.isqrt(N_BARS) * ux * RIGHT
    mark = DashedLine(xs, xs + 3.2 * UP).set_stroke(YELLOW, 2)
    ml = tx(r"x_0=\sqrt N=10", 24, YELLOW).next_to(mark, UP, buff=0.1)
    ttl = mixed("被牆吃掉的比例 ", "$1-P_{100}(x_0)$", "（一維）", font_size=26)
    ttl.move_to(np.array([0, 2.75, 0])).align_to(np.array([6.5, 0, 0]), RIGHT)
    f48 = eq_label(tx(r"\sum_{x_0\ge1}\big(1-P_N(x_0)\big)\,\mu^N\approx\mu^N\sqrt{\frac{2N}{3\pi}}", 30), "(4.8)")
    f48.move_to(np.array([0, 1.75, 0])).align_to(np.array([6.5, 0, 0]), RIGHT)
    lab3 = zh("（三維：所有起點加總）", font_size=22, color=GREY_A).next_to(f48, DOWN, buff=0.15).align_to(f48, RIGHT)
    f48.add(lab3)
    assert_in_frame(bars, base, xl, ticks, mark, ml, ttl, f48)
    assert_apart(ml, ttl)
    assert_apart(ml, f48)
    assert_apart(ttl, f48)
    assert_apart(f48, bars, 0.1)
    s.play(FadeIn(base), FadeIn(ticks), FadeIn(xl), run_time=0.8)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars[0]], lag_ratio=0.05), FadeIn(ttl), run_time=2.0)
    s.say("起點離牆多遠，才不怕懸崖？（", "$N=100$", "）", wait=2.5)
    s.play(ShowCreation(mark), FadeIn(ml), run_time=0.8)
    s.say("離牆 ", "$\\sqrt N$", " 之內出發的，大多活不久", wait=3.0)
    s.play(FadeIn(f48, shift=0.2 * DOWN), run_time=1.0)
    s.say("把所有起點加起來：牆吃掉的走法數", wait=2.5)
    s.say("三維 (4.8)：只有 ", "$\\mu^N\\sqrt{2N/3\\pi}$", "，相對很少", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("碰過牆的路線，正好被從鏡像點出發的路線抵銷", wait=3.0)
    s.say("數活著的路線，只要做一次減法", wait=2.5)


class W401ImageWalker(ChuScene):
    def construct(self):
        self.show_header("W401　鏡像走路人")
        check_math()
        part_question(self)
        part_sixteen(self)
        part_reflect(self)
        part_pairs(self)
        part_survival(self)
        part_eaten(self)
        part_close(self)
        self.wait(1.0)
