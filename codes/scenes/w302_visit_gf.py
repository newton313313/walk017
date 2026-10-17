"""
W302 新的生成函數：用 (1−t)^M 數到訪次數（§3.2–3.3，pp. 51–54，(3.1)–(3.10)，Fig. 3.1，Ex. 3.1）

問：4 步的 16 條路線裡，從沒踩到 +1 的有幾條？踩到一次？兩次？（+1 畫紅圈；腳印算第 0..4 步）
  1. 分三堆（visit_counts_1d(4, 1)）：M = 0：6 條、M = 1：6 條、M = 2：4 條；權重 ∏(1 − tδ_{s_i,S}) = (1−t)^M（(3.3)–(3.5)）；
     (3.1) C'(4,t) = 6 + 6(1−t) + 4(1−t)² = 16 − 14t + 4t²；t 從 0 推到 1，踩到 +1 的路線淡出，只剩 6。
  2. (3.2) 反過來取回三堆：C'(4,1) = 6、−C'_t(1) = −(−14+8) = 6、½C'_tt(1) = 4。
  3. 排容（Fig. 3.1，virial_terms_1d(4, 1) = [16, 14, 4]）：14 = 1·2³ + 3·2（= 6·1 + 4·2）、4 = 1·2·2；t³ 項 0。
  4. 乘 z^N 求和 (3.9)：G_xy − tG_xwG_wy + t²G_xwG_1G_wy − ⋯ = G_xy − tG_xwG_wy/(1 + tG_1)（練習 3.1）；
     (3.10)：微分 M 次取回「到訪恰好 M 次」。
收束：展開是排容，加起來是幾何級數。畫面上的數字都由 utils.walk_rec 算出，check_math() 核對。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rec as R
from utils import walk_gf as GF
from utils import lattice_draw as L
from utils import walk_draw as D

S = 1                                                 # 特別的點 S（Fig. 3.1 的 w⃗）
WALKS4 = R.walks_1d(4)
COUNTS = R.visit_counts_1d(4, S)                      # {M: 條數}
POLY = R.visit_poly(COUNTS)                           # C'(4,t) 的係數（t^0 起）
VIR = R.virial_terms_1d(4, S)                         # Fig. 3.1 各項（不帶號）
DEMO = [0, 1, 0, 1, 2]                                # 開頭的例子：踩到 +1 兩次
UXm, UYm = 0.22, 0.09                                 # 小走法圖
PILE_X = (-4.3, 0.0, 4.3)
BG = "#1C1C1E"


def visits(pos):
    return sum(x == S for x in pos)


def cprime(t):
    return sum(c * t ** k for k, c in enumerate(POLY))


def check_math():
    assert COUNTS == {0: 6, 1: 6, 2: 4} and sum(COUNTS.values()) == 16 == len(WALKS4)
    assert all(sum(visits(p) == M for p in WALKS4) == c for M, c in COUNTS.items())
    assert POLY == [16, -14, 4] and R.counts_from_poly(POLY) == COUNTS
    assert R.visit_counts_gf_1d(4, S) == COUNTS                           # 第二種算法：(3.9)–(3.10)
    for t in (0.0, 0.25, 0.5, 1.0):                                       # 權重 (1−t)^M 加總 = C'(4,t)
        assert abs(sum((1 - t) ** visits(p) for p in WALKS4) - cprime(t)) < 1e-12
        assert abs(sum(c * (1 - t) ** M for M, c in COUNTS.items()) - cprime(t)) < 1e-12
    assert cprime(0) == 16 and cprime(1) == 6 == COUNTS[0]
    assert POLY[1] + 2 * POLY[2] == -14 + 8 and -(POLY[1] + 2 * POLY[2]) == 6 == COUNTS[1]
    assert (2 * POLY[2]) // 2 == 4 == COUNTS[2]                           # ½·d²/dt² = ½·8
    assert VIR == [16, 14, 4] == [abs(a) for a in POLY]
    assert GF.count1(1, 1) * 2 ** 3 == 8 and GF.count1(3, 1) * 2 == 6 and 8 + 6 == VIR[1]
    assert COUNTS[1] * 1 + COUNTS[2] * 2 == VIR[1]                        # 每次到訪都算一次
    assert GF.count1(1, 1) * GF.count1(2, 0) * 2 == 4 == VIR[2] == COUNTS[2]
    assert len(VIR) == 3 and max(map(visits, WALKS4)) == 2                # t³ 項是 0
    assert VIR[0] - VIR[1] + VIR[2] == 6
    assert DEMO in WALKS4 and visits(DEMO) == 2


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tagged(tex, tag, font_size=36):
    return VGroup(Tex(tex, font_size=font_size), zh(tag, font_size=20, color=GREY_B)).arrange(RIGHT, buff=0.3)


def mini(pos, center):
    org = np.array(center) + 2.45 * UXm * LEFT
    ax = D.walk_axes(4, 4, ux=UXm, uy=UYm, origin=org, x_every=4, y_every=4, font_size=12)
    lvl = DashedLine(D.path_point(0, S, UXm, UYm, org), D.path_point(4.5, S, UXm, UYm, org),
                     dash_length=0.06).set_stroke(RED, 1.5)
    path = D.walk_path(pos, ux=UXm, uy=UYm, origin=org, color=WHITE, stroke_width=2.5, dot_radius=0.03)
    hits = VGroup(*[Dot(D.path_point(k, S, UXm, UYm, org), radius=0.055).set_fill(RED, 1)
                    for k in range(5) if pos[k] == S])
    return VGroup(VGroup(ax[0], ax[1]), lvl, path, hits)


def piles():
    out = []
    for M in sorted(COUNTS):
        g = VGroup()
        for k, pos in enumerate([p for p in WALKS4 if visits(p) == M]):
            i, j = divmod(k, 2)
            g.add(mini(pos, np.array([PILE_X[M] + (j - 0.5) * 1.3, 2.0 - 1.0 * i, 0])))
        out.append(g)
    return out


def fig31(k, center, color=TEAL):
    """Fig. 3.1 的示意圖：從 x⃗ 出發，經過大圓點 w⃗ k 次，再到 y⃗。"""
    c = np.array(center, dtype=float)
    P = lambda x, y: c + np.array([x, y, 0.0])
    xs, ys = P(-1.6, -0.2), P(1.6, -0.2)
    pts = [xs, P(-0.85, 0.12), c]
    for a in {1: [], 2: [90], 3: [55, 125]}[k]:
        th = np.radians(a)
        pts += [P(0.5 * np.cos(th - 0.5), 0.5 * np.sin(th - 0.5)), P(0.8 * np.cos(th), 0.8 * np.sin(th)),
                P(0.5 * np.cos(th + 0.5), 0.5 * np.sin(th + 0.5)), c]
    pts += [P(0.85, 0.12), ys]
    curve = VMobject().set_points_smoothly(pts).set_stroke(color, 3)
    dots = VGroup(Dot(xs, radius=0.06), Dot(ys, radius=0.06), Dot(c, radius=0.12).set_fill(RED, 1))
    labs = VGroup(Tex(r"\vec x", font_size=28).next_to(xs, DOWN, buff=0.12),
                  Tex(r"\vec y", font_size=28).next_to(ys, DOWN, buff=0.12),
                  Tex(r"\vec w", font_size=28).set_color(RED).next_to(c, DOWN, buff=0.18))
    sign = Tex({1: "-t", 2: "+t^2", 3: "-t^3"}[k], font_size=38).set_color(YELLOW).move_to(P(0, 1.2))
    return VGroup(curve, dots, labs, sign)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nl = D.number_line(4, unit=0.85, label_every=2).move_to(np.array([0, 1.3, 0]))
    ring = Circle(radius=0.22).move_to(D.line_point(nl, S)).set_stroke(RED, 4)
    slab = mixed("$S=+1$", font_size=24, color=RED).next_to(D.line_point(nl, S), DOWN, buff=0.3)
    note = zh("腳印算第 0、1、2、3、4 步（共 5 個）", font_size=24, color=GREY_B).move_to(np.array([0, 2.75, 0]))
    man = D.walker_dot().move_to(D.line_point(nl, 0))
    prints, seen = VGroup(), {}
    for x in DEMO:
        h = seen.get(x, 0)
        seen[x] = h + 1
        prints.add(Dot(D.line_point(nl, x) + (0.42 + 0.24 * h) * UP, radius=0.08)
                   .set_fill(RED if x == S else YELLOW, 1))
    clab = mixed("踩到 ", "$+1$", " 的次數：", font_size=30)
    cnt = Tex("M=0", font_size=50)
    VGroup(clab, cnt).arrange(RIGHT, buff=0.3).move_to(np.array([0, -1.0, 0]))
    assert_in_frame(nl, ring, slab, note, prints, clab, cnt)
    assert_apart(note, prints)
    assert_apart(VGroup(nl, slab), VGroup(clab, cnt))

    s.say("問：4 步的 16 條路線裡，從沒踩到 ", "$+1$", " 的有幾條？", wait=0.3)
    s.play(FadeIn(nl), run_time=0.6)
    s.play(ShowCreation(ring), FadeIn(slab), run_time=0.8)
    s.wait(1.5)
    s.say("踩到一次的有幾條？踩到兩次的呢？", wait=2.0)
    s.say("先看一條：", r"$0\to1\to0\to1\to2$", wait=0.3)
    s.play(FadeIn(note), FadeIn(clab), FadeIn(cnt), FadeIn(man), FadeIn(prints[0]), run_time=0.8)
    m = 0
    for k in range(1, 5):
        st = DEMO[k] - DEMO[k - 1]
        man.generate_target()
        man.target.move_to(D.line_point(nl, DEMO[k]))
        s.play(MoveToTarget(man, path_arc=-st * 0.7 * PI), run_time=0.5)
        anims = [FadeIn(prints[k], shift=0.1 * DOWN)]
        if DEMO[k] == S:
            m += 1
            anims += [Transform(cnt, Tex(f"M={m}", font_size=50).move_to(cnt)), Indicate(ring, color=RED)]
        s.play(*anims, run_time=0.5)
    s.say("這條踩到 ", "$+1$", " 兩次：", "$M=2$", "（紅腳印疊了兩個）", wait=3.0)
    s.say("把 16 條全部按 ", "$M$", " 分成三堆", wait=1.5)
    s.clear_stage()


def part_piles(s):
    gs = piles()
    labs = [mixed(f"$M={M}$", "：", f"${COUNTS[M]}$", " 條", font_size=26).move_to(np.array([PILE_X[M], 2.85, 0]))
            for M in sorted(COUNTS)]
    wts = [Tex(w, font_size=34).set_color(YELLOW).move_to(np.array([PILE_X[M], -0.7, 0]))
           for M, w in enumerate([r"\times 1", r"\times(1-t)", r"\times(1-t)^2"])]
    f1 = tagged(r"\prod_{i=0}^{4}\Big(1-t\,\delta_{s_i,S}\Big)=(1-t)^{M}", "(3.3)–(3.5)")
    assert POLY[1] < 0
    f2 = tagged(rf"C'(4,t)={COUNTS[0]}+{COUNTS[1]}(1-t)+{COUNTS[2]}(1-t)^2={POLY[0]}-{-POLY[1]}t+{POLY[2]}t^2", "(3.1)")
    for f in (f1, f2):
        fit(f, 12.8)
    f1.move_to(np.array([-6.3 + f1.get_width() / 2, -1.55, 0]))
    f2.move_to(np.array([-6.3 + f2.get_width() / 2, -2.5, 0]))
    tt = ValueTracker(0.0)
    tnum = DecimalNumber(0, num_decimal_places=2, font_size=36)
    cnum = DecimalNumber(16, num_decimal_places=2, font_size=36).set_color(YELLOW)
    rd = VGroup(Tex("t=", font_size=36), tnum, Tex("C'(4,t)=", font_size=36), cnum).arrange(RIGHT, buff=0.15)
    rd[2].shift(0.35 * RIGHT)
    rd[3].shift(0.35 * RIGHT)
    rd.move_to(np.array([6.4 - rd.get_width() / 2, -1.55, 0]))
    assert_in_frame(*gs, *labs, *wts, f1, f2, rd)
    assert_apart(f1, rd, 0.2)
    for a, b in zip(gs, gs[1:]):
        assert_apart(a, b, 0.3)
    assert_apart(VGroup(*gs), VGroup(*wts))
    assert_apart(VGroup(*wts), f1)

    s.say("16 條按踩到 ", "$+1$", " 的次數分三堆（紅點＝踩到 ", "$+1$", "）", wait=0.3)
    for M in sorted(COUNTS):
        s.play(LaggedStart(*[FadeIn(m) for m in gs[M]], lag_ratio=0.15), FadeIn(labs[M]), run_time=1.4)
    s.wait(1.5)
    s.say("數字 6、6、4 就是答案；但能不能「算」出來？", wait=2.5)
    s.say("每條路線給一個權重：每踩到 ", "$+1$", " 一次，就乘一個 ", "$(1-t)$", wait=0.3)
    s.play(Write(f1[0]), FadeIn(f1[1]), run_time=1.5)
    s.play(*[FadeIn(w) for w in wts], run_time=0.8)
    s.wait(2.0)
    s.say("全部加起來，就是新的生成函數 (3.1)", wait=0.3)
    s.play(Write(f2[0]), FadeIn(f2[1]), run_time=1.8)
    s.wait(2.5)
    ovs = []
    for M in (1, 2):
        g = gs[M]
        ov = Rectangle(width=g.get_width() + 0.2, height=g.get_height() + 0.2).move_to(g)
        ov.set_stroke(width=0).set_fill(BG, 0)
        ov.add_updater(lambda o, M=M: o.set_fill(BG, 0.88 * (1 - (1 - tt.get_value()) ** M)))
        ovs.append(ov)
    tnum.add_updater(lambda d: d.set_value(tt.get_value()))
    cnum.add_updater(lambda d: d.set_value(cprime(tt.get_value())))
    s.add(*ovs)
    s.say("讓 ", "$t$", " 從 0 推到 1：", wait=0.3)
    s.play(FadeIn(rd), run_time=0.6)
    s.play(tt.animate.set_value(1.0), run_time=3.5)
    s.say("$t=1$", "：踩到 ", "$+1$", " 的路線全部變成 0，只剩 6 條", wait=3.0)
    s.play(tt.animate.set_value(0.0), run_time=1.5)
    s.say("$t=0$", "：每條都算 1，加起來是 16 條", wait=2.0)
    for m in [*ovs, tnum, cnum]:
        m.clear_updaters()
    s.remove(*ovs)
    part_recover(s, gs, labs, wts, f1, f2, rd)


def part_recover(s, gs, labs, wts, f1, f2, rd):
    d0 = Tex(rf"C'(4,1)={POLY[0]}{POLY[1]:+d}{POLY[2]:+d}={cprime(1)}", font_size=30)
    d1 = Tex(rf"-\frac{{d}}{{dt}}C'\Big|_{{t=1}}=-({POLY[1]}+{2 * POLY[2]})={-(POLY[1] + 2 * POLY[2])}", font_size=30)
    d2 = Tex(rf"\frac{{1}}{{2}}\frac{{d^2}}{{dt^2}}C'\Big|_{{t=1}}=\frac{{{2 * POLY[2]}}}{{2}}={POLY[2]}", font_size=30)
    ds = [d0, d1, d2]
    for M, d in enumerate(ds):
        fit(d, 4.0)
        d.move_to(np.array([PILE_X[M], -2.4, 0]))
    dy = -1.45 - f2.get_center()[1]
    assert_in_frame(*ds)
    for a, b in zip(ds, ds[1:]):
        assert_apart(a, b, 0.15)
    assert f2.get_bottom()[1] + dy > ds[1].get_top()[1] + 0.05

    s.say("反過來：從多項式取回三堆的條數 (3.2)", wait=0.3)
    s.play(FadeOut(f1), FadeOut(rd), f2.animate.shift(dy * UP), run_time=1.0)
    for M, d in enumerate(ds):
        s.play(Write(d), run_time=1.2)
        s.play(Indicate(labs[M]), run_time=0.8)
        s.wait(1.0)
    s.say("一般：", r"$C(M)=\tfrac{(-1)^M}{M!}\,\tfrac{d^M}{dt^M}C'\big|_{t=1}$", "（微分 ", "$M$", " 次）", wait=3.5)
    s.clear_stage()


def part_inclusion(s):
    figs = [fig31(k, np.array([PILE_X[k - 1], 1.65, 0])) for k in (1, 2, 3)]
    c1 = VGroup(mixed("在第 ", "$N_1$", " 步到 ", "$+1$", "，之後隨便走", font_size=20),
                Tex(rf"N_1=1:\ {GF.count1(1, 1)}\cdot 2^3={GF.count1(1, 1) * 8}", font_size=30),
                Tex(rf"N_1=3:\ {GF.count1(3, 1)}\cdot 2={GF.count1(3, 1) * 2}", font_size=30),
                Tex(rf"8+6={VIR[1]}", font_size=36).set_color(YELLOW)).arrange(DOWN, buff=0.28)
    c2 = VGroup(mixed("到 ", "$+1$", "、走 2 步回到 ", "$+1$", "、再走 1 步", font_size=20),
                Tex(rf"{GF.count1(1, 1)}\cdot{GF.count1(2, 0)}\cdot 2={VIR[2]}", font_size=36).set_color(YELLOW)
                ).arrange(DOWN, buff=0.28)
    c3 = VGroup(zh("4 步內到訪不了 3 次", font_size=20), Tex("0", font_size=36).set_color(YELLOW)).arrange(DOWN, buff=0.28)
    cols = [c1, c2, c3]
    for k, c in enumerate(cols):
        fit(c, 4.0)
        c.move_to(np.array([PILE_X[k], 0.75 - c.get_height() / 2, 0]))
    bl = Tex(rf"C'(4,t)={VIR[0]}-{VIR[1]}t+{VIR[2]}t^2", font_size=44,
             isolate=[str(VIR[0]), f"{VIR[1]}t", f"{VIR[2]}t^2"])
    bl.move_to(np.array([0, -2.3, 0]))
    assert_in_frame(*figs, *cols, bl)
    for a, b in zip(figs, figs[1:]):
        assert_apart(a, b, 0.2)
    for f, c in zip(figs, cols):
        assert_apart(f, c, 0.15)
    assert_apart(VGroup(*cols), bl, 0.15)

    s.say("把 ", "$(1-t)^M$", " 展開：每一項都有意思（Fig. 3.1）", wait=0.3)
    s.play(FadeIn(bl), run_time=0.8)
    s.wait(1.5)
    s.say("示意圖：從 ", r"$\vec x$", " 出發，經過 ", r"$\vec w$", " 一次、兩次、三次，再到 ", r"$\vec y$", wait=0.3)
    s.play(LaggedStart(*[FadeIn(f) for f in figs], lag_ratio=0.4), run_time=2.4)
    s.wait(1.5)
    s.say("$-t$", " 項的 14：在第 ", "$N_1$", " 步到 ", "$+1$", "，之後隨便走", wait=0.3)
    s.play(FadeIn(c1[0]), run_time=0.6)
    s.play(FadeIn(c1[1]), run_time=0.8)
    s.play(FadeIn(c1[2]), run_time=0.8)
    s.play(FadeIn(c1[3]), Indicate(bl.select_part(f"{VIR[1]}t")), run_time=1.0)
    s.wait(1.0)
    s.say("每次到訪都算一次：", r"$6\cdot1+4\cdot2=14$", wait=3.5)
    s.say("$+t^2$", " 項的 4：到 ", "$+1$", "、走 2 步再回到 ", "$+1$", "、再隨便走 1 步", wait=0.3)
    s.play(FadeIn(c2[0]), run_time=0.6)
    s.play(FadeIn(c2[1]), Indicate(bl.select_part(f"{VIR[2]}t^2")), run_time=1.0)
    s.wait(2.5)
    s.say("$-t^3$", " 項：4 步內到訪不了 3 次，所以是 0", wait=0.3)
    s.play(FadeIn(c3), run_time=0.8)
    s.wait(1.5)
    s.say("排容：", r"$16-14+4=6$", "，正是從沒踩到 ", "$+1$", " 的條數", wait=3.5)
    part_series(s, figs, cols, bl)


def part_series(s, figs, cols, bl):
    terms = [Tex(t, font_size=34).set_color(TEAL_B) for t in
             (r"-t\,G_{xw}G_{wy}", r"+t^2G_{xw}G_1G_{wy}", r"-t^3G_{xw}G_1^2G_{wy}")]
    for k, t in enumerate(terms):
        t.move_to(np.array([PILE_X[k], 0.45, 0]))
    e1 = Tex(r"G(z,t;\vec x,\vec y,\vec w)=G_{xy}-tG_{xw}G_{wy}+t^2G_{xw}G_1G_{wy}-\cdots", font_size=38)
    tag1 = zh("(3.9)", font_size=20, color=GREY_B)
    e2 = Tex(r"=G_{xy}-t\,\frac{G_{xw}G_{wy}}{1+tG_1(z;\vec w,\vec w)}", font_size=38)
    lhs_w = Tex(r"G(z,t;\vec x,\vec y,\vec w)", font_size=38).get_width()
    e1.move_to(np.array([-0.4, -0.65, 0]))
    tag1.next_to(e1, RIGHT, buff=0.3)
    e2.next_to(e1, DOWN, buff=0.45).align_to(e1, LEFT).shift((lhs_w + 0.05) * RIGHT)
    e2.set_color(YELLOW)
    e3 = Tex(r"\frac{(-1)^M}{M!}\frac{\partial^M G}{\partial t^M}\Big|_{t=1}"
             r"=\frac{G_{xw}\,G_1^{M-1}\,G_{wy}}{(1+G_1)^{M+1}}", font_size=40)
    tag3 = mixed("(3.10)，", "$M\\ge 1$", font_size=20, color=GREY_B)
    e4 = Tex(r"M=0:\ G_{xy}-\frac{G_{xw}G_{wy}}{1+G_1}", font_size=36)
    up = 2.7
    assert_in_frame(*terms, e1, tag1, e2)
    for f, t in zip(figs, terms):
        assert_apart(f, t, 0.05)
    assert_apart(VGroup(*terms), e1, 0.15)
    e3.move_to(np.array([-0.6, -0.75, 0]))
    tag3.next_to(e3, RIGHT, buff=0.3)
    e4.next_to(e3, DOWN, buff=0.5).align_to(e3, LEFT)
    assert_in_frame(e1.copy().shift(up * UP), e2.copy().shift(up * UP), e3, tag3, e4)
    assert e2.get_bottom()[1] + up > e3.get_top()[1] + 0.2

    s.say("乘上 ", "$z^N$", " 對 ", "$N$", " 求和：分段走的步數相加，卷積變乘積", wait=0.3)
    s.play(FadeOut(VGroup(*cols)), FadeOut(bl), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(t, shift=0.2 * UP) for t in terms], lag_ratio=0.4), run_time=2.0)
    s.wait(2.0)
    s.say("$G_1=G-1$", "：從 ", r"$\vec w$", " 回到 ", r"$\vec w$", "，至少走一步", wait=2.5)
    s.play(Write(e1), FadeIn(tag1), run_time=2.0)
    s.wait(1.5)
    s.say("公比 ", "$-tG_1$", " 的幾何級數，加起來是（練習 3.1）", wait=0.3)
    s.play(Write(e2), run_time=1.8)
    s.wait(3.0)
    s.play(FadeOut(VGroup(*figs)), FadeOut(VGroup(*terms)), run_time=0.8)
    s.play(e1.animate.shift(up * UP), tag1.animate.shift(up * UP), e2.animate.shift(up * UP), run_time=1.0)
    s.say("(3.10)：對 ", "$t$", " 微分 ", "$M$", " 次，取回「到訪恰好 ", "$M$", " 次」的生成函數", wait=0.3)
    s.play(Write(e3), FadeIn(tag3), run_time=2.0)
    s.wait(2.0)
    s.say("$M=0$", "：就是 ", "$t=1$", "，從沒到訪 ", r"$\vec w$", wait=0.3)
    s.play(FadeIn(e4), run_time=1.0)
    s.wait(2.5)


def part_close(s):
    s.say("$(1-t)^M$", " 把「到訪幾次」變成 ", "$t$", " 的多項式", wait=2.5)
    s.say("展開是排容，加起來是幾何級數", wait=3.0)


class W302VisitGF(ChuScene):
    def construct(self):
        self.show_header("W302　新的生成函數：用 (1−t) 的 M 次方數到訪次數")
        check_math()
        part_question(self)
        part_piles(self)
        part_inclusion(self)
        part_close(self)
        self.wait(1.0)
