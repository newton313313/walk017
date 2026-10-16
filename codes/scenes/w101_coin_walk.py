"""
W101 擲硬幣走路（§1.1，pp. 1–2，式 (1.1)）

十次擲幣 HTHHTHTTHH：數線上的人跟著走、下方同步畫走法圖，6 正 4 反 → d = 2；三個應用；
1024 人各擲一套（全部 2^10 種結果各一人）的終點分布 C(10,d)；每一條指定路線的機率都是 1/1024。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
from fractions import Fraction

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

N = 10
COINS = W.MAIN_COINS
POS = W.MAIN_POS
HEADS = COINS.count("H")
D_END = POS[-1]
DIST = W.distribution(N)
TOTAL = 2 ** N
PCT0 = round(100 * DIST[0] / TOTAL, 1)

COIN_Y, NL_Y = 2.55, 1.05
AX_O = np.array([-6.1, -1.45, 0.0])
UX, UY = 0.55, 0.4
BAR_O = np.array([0.0, -1.5, 0.0])
BAR_UY = 3.1 / 252


def check_math():
    assert COINS == "HTHHTHTTHH" and POS == [0, 1, 0, 1, 2, 1, 2, 1, 0, 1, 2]
    assert W.positions(W.coins_to_steps(COINS)) == POS
    assert HEADS == 6 and N - HEADS == 4
    assert D_END == 2 * HEADS - N == 2                                  # (1.1)，l = 1
    assert [DIST[d] for d in range(-N, N + 1, 2)] == [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]
    assert sum(DIST.values()) == TOTAL == 1024
    assert max(DIST, key=DIST.get) == 0 and DIST[0] == 252 and DIST[2] == 210 and DIST[N] == 1
    assert PCT0 == 24.6
    assert all(W.walk_count(N, d) == 0 for d in range(-N + 1, N, 2))     # 奇數位置沒有人
    assert all(d == 2 * ((d + N) // 2) - N for d in DIST)                 # d = 2n − 10
    assert Fraction(1, 2) ** N == Fraction(1, TOTAL)
    assert W.walk_prob(N, 0) == Fraction(252, 1024)


# ---------------------------------------------------------------- 1. 問題
def ask(self):
    self.say("擲 10 次硬幣：正面向右一步，反面向左一步", wait=3.0)
    self.say("最後最可能停在哪裡？", wait=3.0)


# ---------------------------------------------------------------- 2. 走十步
def placeholder(c):
    circ = Circle(radius=c[0].get_width() / 2).set_stroke(GREY_B, 2).move_to(c)
    return VGroup(circ, Tex("?", font_size=26).set_color(GREY_B).move_to(c))


def legend():
    rows = VGroup()
    for face, tex, direc in (("H", "+1", RIGHT), ("T", "-1", LEFT)):
        col = D.HEAD_COLOR if face == "H" else D.TAIL_COLOR
        c = D.coin(face, radius=0.2)
        t = Tex(tex, font_size=30).set_color(col)
        a = Arrow(ORIGIN, 0.6 * direc, buff=0, thickness=2.5).set_color(col)
        rows.add(VGroup(c, t, a).arrange(RIGHT, buff=0.18))
    rows.arrange(DOWN, buff=0.18, aligned_edge=LEFT)
    return rows.move_to(np.array([5.3, 2.6, 0]))


def coin_walk(self):
    row = D.coin_row(COINS, radius=0.27).move_to(COIN_Y * UP)
    holders = [placeholder(c) for c in row]
    nl = D.number_line(N, unit=0.5).shift(NL_Y * UP)
    man = D.walker_dot().move_to(D.line_point(nl, 0))
    ax = D.walk_axes(N, 3, ux=UX, uy=UY, origin=AX_O, x_every=2, y_every=1)
    start = Dot(D.path_point(0, 0, UX, UY, AX_O), radius=0.06).set_fill(D.WALK_COLOR, 1)
    leg = legend()
    assert_in_frame(row, nl, ax, leg)
    assert_apart(row, leg)
    assert_apart(row, nl)
    assert_apart(nl, ax)
    self.play(LaggedStart(*[FadeIn(h) for h in holders], lag_ratio=0.08), FadeIn(leg), run_time=1.4)
    self.play(ShowCreation(nl[0]), FadeIn(nl[1]), FadeIn(nl[2]), FadeIn(man, scale=0.5))
    ask(self)
    self.play(FadeIn(ax), FadeIn(start))
    self.say("上：數線上的人　下：走法圖（橫軸步數、縱軸位置）", wait=3.0)
    notes = {0: "第 1 枚：正面，向右走一步", 1: "第 2 枚：反面，向左走一步", 2: "以後都一樣：一枚硬幣，一步"}
    for k in range(N):
        s = W.MAIN_STEPS[k]
        col = D.HEAD_COLOR if s > 0 else D.TAIL_COLOR
        lab = Tex("+1" if s > 0 else "-1", font_size=24).set_color(col).next_to(row[k], DOWN, buff=0.1)
        p0 = D.path_point(k, POS[k], UX, UY, AX_O)
        p1 = D.path_point(k + 1, POS[k + 1], UX, UY, AX_O)
        seg = Line(p0, p1).set_stroke(D.WALK_COLOR, 4)
        dot = Dot(p1, radius=0.06).set_fill(D.WALK_COLOR, 1)
        if k in notes:
            self.say(notes[k])
        self.play(FadeTransform(holders[k], row[k]), FadeIn(lab, shift=0.1 * DOWN),
                  man.animate.move_to(D.line_point(nl, POS[k + 1])), ShowCreation(seg),
                  run_time=1.4 if k < 2 else 0.9)
        self.add(dot)
        self.wait(1.0 if k < 2 else 0.25)
    end_lab = Tex("d=2", font_size=30).set_color(D.WALK_COLOR)
    end_lab.next_to(D.path_point(N, D_END, UX, UY, AX_O), RIGHT, buff=0.15)
    self.play(FadeIn(end_lab), Indicate(man, scale_factor=1.6))
    self.say("位置依序：", "$0,1,0,1,2,1,2,1,0,1,2$", wait=3.5)
    return row


def tally(self, row):
    hidx = [k for k in range(N) if COINS[k] == "H"]
    tidx = [k for k in range(N) if COINS[k] == "T"]
    hrow = D.coin_row("H" * len(hidx), radius=0.27, buff=0.08)
    trow = D.coin_row("T" * len(tidx), radius=0.27, buff=0.08)
    ht = Tex(r"n=6", font_size=36).set_color(D.HEAD_COLOR)
    tt = Tex(r"N-n=4", font_size=36).set_color(D.TAIL_COLOR)
    block = VGroup(VGroup(hrow, ht).arrange(RIGHT, buff=0.3),
                   VGroup(trow, tt).arrange(RIGHT, buff=0.3)).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    block.move_to(np.array([3.6, -0.5, 0]))
    form = Tex(r"d=l(2n-N)=2\cdot 6-10=2", font_size=40)
    form.next_to(block, DOWN, buff=0.45).align_to(block, LEFT)
    eqn = Tex("(1.1)", font_size=26).set_color(GREY_B).next_to(form, DOWN, buff=0.15).align_to(form, RIGHT)
    box = SurroundingRectangle(form, buff=0.12).set_stroke(D.WALK_COLOR, 2)
    assert_in_frame(block, form, eqn, box)
    self.say("數一數：", "$6$", " 個正面、", "$4$", " 個反面")
    self.play(*[TransformFromCopy(row[k], hrow[i]) for i, k in enumerate(hidx)], run_time=1.3)
    self.play(Write(ht))
    self.play(*[TransformFromCopy(row[k], trow[i]) for i, k in enumerate(tidx)], run_time=1.3)
    self.play(Write(tt))
    self.wait(1.5)
    self.say("向右 ", "$n$", " 步、向左 ", "$N-n$", " 步：終點 ", "$d=l(2n-N)$")
    self.play(Write(form), FadeIn(eqn), run_time=1.6)
    self.play(ShowCreation(box))
    self.wait(1.5)
    self.say("步長 ", "$l$", " 在畫面上取 ", "$1$", "：這一次停在 ", "$d=2$", wait=3.5)


# ---------------------------------------------------------------- 3. 三個應用
def tube_demo():
    tube = Rectangle(width=3.4, height=0.5).set_stroke(GREY_B, 2)
    mol = Dot(radius=0.1).set_fill(TEAL, 1)
    return VGroup(tube, mol)


def chain_demo():
    pts = [np.array([0.7 * x, 0.45 - 0.09 * k, 0]) for k, x in enumerate(POS)]
    line = VMobject().set_points_as_corners(pts).set_stroke(PURPLE_B, 4)
    beads = VGroup(*[Dot(p, radius=0.06).set_fill(WHITE, 1) for p in pts])
    return VGroup(line, beads)


def bet_demo():
    win = VGroup(D.coin("H", 0.25), zh("贏 1 元", font_size=24, color=D.HEAD_COLOR)).arrange(RIGHT, buff=0.15)
    lose = VGroup(D.coin("T", 0.25), zh("輸 1 元", font_size=24, color=D.TAIL_COLOR)).arrange(RIGHT, buff=0.15)
    return VGroup(win, lose).arrange(RIGHT, buff=0.5)


def applications(self):
    self.say("同一個模型，到處都看得到", wait=0.5)
    pics = [tube_demo(), chain_demo(), bet_demo()]
    texts = [zh("分子在一條細管子裡擴散", font_size=28),
             zh("高分子鏈：每一節只有兩種取向", font_size=28),
             zh("擲硬幣的賭局：de Moivre 的賭徒破產問題", font_size=28)]
    for pic, txt, y in zip(pics, texts, (2.1, 0.4, -1.4)):
        pic.move_to(np.array([-4.2, y, 0]))
        txt.move_to(np.array([0, y, 0])).align_to(np.array([-1.3, 0, 0]), LEFT)
        assert_in_frame(pic, txt)
        assert_apart(pic, txt)
    tube, mol = pics[0]
    self.play(FadeIn(tube), FadeIn(mol), FadeIn(texts[0], shift=0.2 * LEFT), run_time=0.9)
    for x in POS[1:]:
        self.play(mol.animate.move_to(tube.get_center() + 0.3 * x * RIGHT), run_time=0.2)
    self.wait(2.0)
    line, beads = pics[1]
    self.play(FadeIn(beads), FadeIn(texts[1], shift=0.2 * LEFT), run_time=0.9)
    self.play(ShowCreation(line), run_time=1.4)
    self.wait(2.0)
    win, lose = pics[2]
    self.play(FadeIn(win), FadeIn(lose), FadeIn(texts[2], shift=0.2 * LEFT), run_time=0.9)
    self.play(Indicate(win), run_time=0.8)
    self.play(Indicate(lose), run_time=0.8)
    self.wait(2.0)
    self.say("共同點：每一步只有兩種選擇，而且彼此獨立", wait=4.0)


# ---------------------------------------------------------------- 4. 1024 人的終點分布
def apart_from_chart(m, ch):
    for part in [*ch[0], *ch[1]]:
        assert_apart(m, part)


def distribution(self):
    top = mixed("1024 個人各擲一套：全部 ", "$2^{10}$", " 種結果，每種一人", font_size=28).move_to(2.7 * UP)
    nl = D.number_line(N, unit=0.5).shift(BAR_O)
    ch = D.bar_chart(DIST, ux=0.5, uy=BAR_UY, origin=BAR_O, width=1.4, labels=True, font_size=22)
    nrow = VGroup(*[Tex(str((d + N) // 2), font_size=20).set_color(D.HEAD_COLOR)
                    .move_to(D.line_point(nl, d) + 0.7 * DOWN) for d in range(-N, N + 1, 2)])
    dlab = Tex("d", font_size=26).set_color(GREY_B).move_to(D.line_point(nl, -N) + np.array([-0.8, -0.33, 0]))
    nlab = Tex("n", font_size=26).set_color(D.HEAD_COLOR).move_to(D.line_point(nl, -N) + np.array([-0.8, -0.7, 0]))
    assert_in_frame(top, nl, ch, nrow, dlab, nlab)
    assert_apart(top, ch)
    self.say("如果不是一個人，而是 1024 個人呢？")
    self.play(FadeIn(top))
    self.play(ShowCreation(nl[0]), FadeIn(nl[1]), FadeIn(nl[2]), FadeIn(dlab))
    self.say("每個人走完十步，站到自己的終點上", wait=0.5)
    self.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in ch[0]], lag_ratio=0.08), run_time=2.2)
    self.play(LaggedStart(*[FadeIn(t) for t in ch[1]], lag_ratio=0.05), run_time=1.0)
    self.say("長條上的數字：停在那個位置的人數（合計 1024）", wait=3.5)
    self.say("正面數 ", "$n$", " 決定終點：", "$d=2n-10$", wait=0.3)
    self.play(FadeIn(nrow, lag_ratio=0.05), FadeIn(nlab))
    self.wait(2.0)
    # 揭曉
    a0 = VGroup(mixed("最可能：", "$d=0$", "，252 人", font_size=26, color=D.EXACT_COLOR),
                mixed("$252/1024\\approx 24.6\\%$", font_size=24, color=D.EXACT_COLOR)).arrange(DOWN, buff=0.12,
                                                                                         aligned_edge=LEFT)
    a2 = mixed("主例 ", "$d=2$", "：210 人", font_size=26, color=D.WALK_COLOR)
    a0.move_to(np.array([0, 1.95, 0])).align_to(np.array([1.8, 0, 0]), LEFT)
    a2.next_to(a0, DOWN, buff=0.3, aligned_edge=LEFT)
    assert_in_frame(a0, a2)
    apart_from_chart(a0, ch)
    apart_from_chart(a2, ch)
    self.say("揭曉：最可能停在原點 ", "$d=0$")
    self.play(D.bar_of(ch, 0).animate.set_fill(D.EXACT_COLOR, 0.9).set_stroke(D.EXACT_COLOR), FadeIn(a0))
    self.wait(2.5)
    self.say("我們剛才的 ", "$d=2$", "：210 人，只比 252 少一點")
    self.play(D.bar_of(ch, 2).animate.set_fill(D.WALK_COLOR, 0.85).set_stroke(D.WALK_COLOR), FadeIn(a2))
    self.wait(2.5)
    zeros = VGroup(*[Tex("0", font_size=20).set_color(RED).move_to(D.line_point(nl, x) + 0.2 * UP)
                     for x in range(-N + 1, N, 2)])
    self.say("奇數位置：一個人都沒有")
    self.play(LaggedStart(*[FadeIn(z, scale=1.5) for z in zeros], lag_ratio=0.08), run_time=1.5)
    self.wait(1.5)
    self.say("走 10 步，", "$d=2n-10$", " 一定是偶數：奇偶性", wait=3.0)
    self.play(FadeOut(top), FadeOut(a0), FadeOut(a2))
    return ch


# ---------------------------------------------------------------- 5. 每條路線 1/1024
def each_path(self, ch):
    rows = VGroup()
    for coins in (COINS, "H" * N):
        r = D.coin_row(coins, radius=0.15, buff=0.05)
        t = Tex(r"=\big(\tfrac12\big)^{10}=\tfrac{1}{1024}", font_size=30)
        rows.add(VGroup(r, t).arrange(RIGHT, buff=0.2))
    rows.arrange(DOWN, buff=0.22, aligned_edge=LEFT)
    rows.move_to(np.array([0, 2.5, 0])).align_to(np.array([-6.5, 0, 0]), LEFT)
    assert_in_frame(rows)
    apart_from_chart(rows, ch)
    self.say("每一條指定的路線，機率都是 ", "$(1/2)^{10}=1/1024$")
    self.play(FadeIn(rows[0], shift=0.2 * DOWN))
    self.wait(1.5)
    self.play(FadeIn(rows[1], shift=0.2 * DOWN))
    self.wait(2.0)
    self.say("全正只有 1 條路線；", "$d=2$", " 有 210 條；", "$d=0$", " 有 252 條")
    self.play(Indicate(ch[1][-1], scale_factor=2.0), Indicate(rows[1][0]), run_time=1.4)
    self.play(Indicate(D.bar_of(ch, 2)), Indicate(rows[0][0]), run_time=1.4)
    self.play(Indicate(D.bar_of(ch, 0)), run_time=1.2)
    self.wait(2.0)
    self.say("停在 ", "$d$", " 的機率＝通往 ", "$d$", " 的路線數 ", "$\\times\\,1/1024$", wait=4.0)
    self.say("每條路線一樣可能", wait=3.0)
    self.say("不一樣的是：有幾條路線通往同一個終點", wait=3.5)


class W101CoinWalk(ChuScene):
    def construct(self):
        self.show_header("W101　擲硬幣走路")
        check_math()
        row = coin_walk(self)
        tally(self, row)
        self.clear_stage()
        applications(self)
        self.clear_stage()
        ch = distribution(self)
        each_path(self, ch)
        self.wait(1.0)
