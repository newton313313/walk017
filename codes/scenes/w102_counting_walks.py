"""
W102 數路線：二項式係數（§1.1，pp. 2–5，式 (1.2)–(1.7)，Fig. 1.1–1.2，Exercise 1.1）

N = 3 的 8 個序列按正面數分組（1, 3, 3, 1）；洗牌論證 3!/2! = 3 → (1.4)；
換成位置 (1.2)(1.6)(1.7)，N = 4 的 1, 4, 6, 4, 1；揭曉 C(10,0) = 252；偏硬幣 p = 2/3 的機率 1, 6, 12, 8 (/27)。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
import itertools
import math
from fractions import Fraction

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk1d as W
from utils import walk_draw as D

SEQ3 = {n: W.sequences_with_heads(3, n) for n in range(4)}
HOME10 = [W.positions(w) for w in W.all_walks(10) if sum(w) == 0]
PERM_GROUPS = {}
for _p in itertools.permutations(["H1", "T", "H2"]):
    PERM_GROUPS.setdefault("".join(t[0] for t in _p), []).append(_p)
DIST4 = W.distribution(4)
C10 = W.walk_count(10, 0)
P10 = W.walk_prob(10, 0)
P_BIASED = [W.biased_prob(3, n, Fraction(2, 3)) for n in range(4)]
P_FAIR = [W.walk_prob(3, 2 * n - 3) for n in range(4)]
PATH_COLORS = {"HTH": D.WALK_COLOR, "HHT": TEAL_B, "THH": PINK}


def check_math():
    assert SEQ3 == {0: ["TTT"], 1: ["HTT", "THT", "TTH"], 2: ["HHT", "HTH", "THH"], 3: ["HHH"]}
    assert [len(SEQ3[n]) for n in range(4)] == [1, 3, 3, 1] == [W.walk_count(3, 2 * n - 3) for n in range(4)]
    assert sum(len(v) for v in SEQ3.values()) == 8 == 2 ** 3
    assert W.positions(W.coins_to_steps("HTH")) == [0, 1, 0, 1]
    assert len(HOME10) == C10 == 252
    assert {k: len(v) for k, v in PERM_GROUPS.items()} == {"HTH": 2, "HHT": 2, "THH": 2}
    assert math.factorial(3) == 6 and math.factorial(3) // math.factorial(2) == 3 == math.comb(3, 2)
    assert [DIST4[d] for d in range(-4, 5, 2)] == [1, 4, 6, 4, 1] and sum(DIST4.values()) == 16 == 2 ** 4
    for N, d in ((4, -2), (4, 0), (10, 0), (10, 2)):
        n = (N + d) // 2                                                   # (1.2)，l = 1
        assert W.walk_count(N, d) == math.factorial(N) // (math.factorial(n) * math.factorial(N - n))   # (1.6)
    assert C10 == math.factorial(10) // math.factorial(5) ** 2 == 252
    assert P10 == Fraction(252, 1024) == Fraction(63, 256) and round(float(P10), 3) == 0.246   # (1.7)
    assert P_BIASED == [Fraction(k, 27) for k in (1, 6, 12, 8)] and sum(P_BIASED) == 1
    assert P_FAIR == [Fraction(k, 8) for k in (1, 3, 3, 1)]


# ---------------------------------------------------------------- 1. 問題
def ask(self):
    o = np.array([-5.6, 0.3, 0])
    ax = D.walk_axes(10, 5, ux=0.5, uy=0.33, origin=o, x_every=2, y_every=5)
    paths = [D.walk_path(HOME10[i], ux=0.5, uy=0.33, origin=o, color=c, stroke_width=3, dots=False)
             for i, c in ((0, TEAL_B), (100, D.WALK_COLOR), (-1, PINK))]
    q = Tex("?", font_size=150).set_color(D.WALK_COLOR).move_to(np.array([3.6, 0.4, 0]))
    assert_in_frame(ax, q, *paths)
    assert_apart(ax, q)
    self.play(FadeIn(ax))
    for p in paths:
        self.play(ShowCreation(p[0]), run_time=1.2)
    self.say("走 10 步又回到原點的路線，一共有幾條？", wait=3.0)
    self.play(Write(q))
    self.say("先猜一個數字", wait=3.5)


# ---------------------------------------------------------------- 2. N = 3：八個序列
def seq_column(n, x):
    d = 2 * n - 3
    head = Tex(f"d={d}", font_size=32).move_to(np.array([x, 2.7, 0]))
    sub = zh(f"{n} 正 {3 - n} 反", font_size=22, color=GREY_B).move_to(np.array([x, 2.27, 0]))
    rows = VGroup(*[D.coin_row(s, radius=0.18, buff=0.06).move_to(np.array([x, 1.65 - 0.55 * i, 0]))
                    for i, s in enumerate(SEQ3[n])])
    cnt = mixed(f"${len(SEQ3[n])}$", " 種", font_size=28, color=YELLOW).move_to(np.array([x, -0.35, 0]))
    return head, sub, rows, cnt


def three_coins(self):
    o = np.array([-6.1, -0.1, 0])
    ux, uy = 0.85, 0.42
    ax = D.walk_axes(3, 3, ux=ux, uy=uy, origin=o, x_every=1, y_every=1)
    row = D.coin_row("HTH", radius=0.25).move_to(np.array([-4.6, 2.45, 0]))
    paths = {s: D.walk_path(W.positions(W.coins_to_steps(s)), ux=ux, uy=uy, origin=o, color=PATH_COLORS[s],
                            stroke_width=4) for s in ("HTH", "HHT", "THH")}
    cols = [seq_column(n, x) for n, x in zip(range(4), (-1.2, 1.2, 3.6, 6.0))]
    total = Tex(r"1+3+3+1=8=2^3", font_size=36).move_to(np.array([2.4, -1.25, 0]))
    assert_in_frame(ax, row, total, *[m for c in cols for m in c])
    assert_apart(ax, row)
    assert_apart(ax, VGroup(*[c[2] for c in cols]))
    self.say("先看 ", "$N=3$", "：例如「正反正」（Fig. 1.1）")
    self.play(FadeIn(row, lag_ratio=0.3), FadeIn(ax))
    self.play(ShowCreation(paths["HTH"][0]), FadeIn(paths["HTH"][1]), run_time=1.5)
    self.wait(2.5)
    self.say("三枚硬幣一共 ", "$2^3=8$", " 種結果：按正面數分組", wait=0.5)
    for head, sub, rows, cnt in cols:
        self.play(FadeIn(head), FadeIn(sub), run_time=0.5)
        self.play(LaggedStart(*[FadeIn(r, shift=0.15 * DOWN) for r in rows], lag_ratio=0.3), run_time=1.0)
    self.wait(1.5)
    self.play(*[FadeIn(c[3]) for c in cols])
    self.play(Write(total))
    self.say("每一組的個數：", "$1,\\ 3,\\ 3,\\ 1$", wait=3.5)
    rows_d1 = cols[2][2]
    boxes = VGroup(*[SurroundingRectangle(r, buff=0.06).set_stroke(PATH_COLORS[s], 2)
                     for r, s in zip(rows_d1, SEQ3[2])])
    self.say("停在 ", "$d=1$", " 的 3 條路線，畫在左邊的走法圖上")
    self.play(ShowCreation(boxes[1]))
    for i, s in ((0, "HHT"), (2, "THH")):
        self.play(ShowCreation(boxes[i]), ShowCreation(paths[s][0]), FadeIn(paths[s][1]), run_time=1.3)
    self.wait(3.5)


# ---------------------------------------------------------------- 3. 洗牌論證
def tagged_row(tags, radius=0.3, step=0.9):
    coins, subs = VGroup(), VGroup()
    for i, tag in enumerate(tags):
        c = D.coin(tag[0], radius).move_to(i * step * RIGHT)
        coins.add(c)
        if len(tag) > 1:
            s = Tex(tag[1], font_size=28)
            s.next_to(c.get_corner(DR), RIGHT, buff=0.02).shift(0.06 * UP)
            subs.add(s)
    return VGroup(coins, subs)


def shuffle(self):
    cols = []
    for x, (pat, perms) in zip((-4.6, -1.5, 1.6), PERM_GROUPS.items()):
        rows = VGroup(*[tagged_row(p).move_to(np.array([x, y, 0])) for p, y in zip(perms, (1.55, 0.45))])
        box = SurroundingRectangle(rows, buff=0.15).set_stroke(PATH_COLORS[pat], 2)
        lab = zh("同一種", font_size=24, color=PATH_COLORS[pat]).next_to(box, DOWN, buff=0.12)
        cols.append((rows, box, lab))
    f1 = Tex(r"3!=3\cdot2\cdot1=6", font_size=38).move_to(np.array([4.9, 1.55, 0]))
    f2 = Tex(r"\frac{3!}{2!}=\frac{6}{2}=3", font_size=38).move_to(np.array([4.9, 0.3, 0]))
    gen = Tex(r"\binom{N}{n}=\frac{N!}{n!\,(N-n)!}", font_size=44).move_to(np.array([0, -1.75, 0]))
    gbox = SurroundingRectangle(gen, buff=0.15).set_stroke(YELLOW, 2)
    eqn = Tex("(1.4)", font_size=26).set_color(GREY_B).next_to(gbox, RIGHT, buff=0.25)
    allm = [m for c in cols for m in c]
    assert_in_frame(f1, f2, gen, gbox, eqn, *allm)
    for c in cols:
        assert_apart(c[1], f1)
        assert_apart(c[2], gbox)
    assert_apart(cols[0][1], cols[1][1])
    assert_apart(cols[1][1], cols[2][1])
    self.say("為什麼是 3 種？把兩個正面編上號碼 1、2（Fig. 1.2）", wait=1.5)
    self.play(LaggedStart(*[FadeIn(r) for c in cols for r in c[0]], lag_ratio=0.2), run_time=2.4)
    self.say("三個「不同」的硬幣排成一列：", "$3!=6$", " 種")
    self.play(Write(f1))
    self.wait(3.0)
    self.say("擦掉號碼：每一欄的兩列變得一模一樣")
    self.play(*[FadeOut(r[1]) for c in cols for r in c[0]])
    self.play(*[ShowCreation(c[1]) for c in cols], *[FadeIn(c[2]) for c in cols])
    self.wait(2.5)
    self.say("兩個正面互換看起來一樣：除以 ", "$2!$")
    self.play(Write(f2))
    self.wait(3.0)
    self.say("一般：", "$N$", " 枚排列，再除掉正面之間、反面之間的互換")
    self.play(Write(gen), ShowCreation(gbox), FadeIn(eqn), run_time=1.5)
    self.wait(3.5)


# ---------------------------------------------------------------- 4. 換成位置；揭曉 C(10,0)
def positions_part(self):
    rows = VGroup()
    for tex, num in ((r"n=\tfrac12\big(N+d/l\big)", "(1.2)"),
                     (r"C(N,d)=\frac{N!}{\big(\frac{N+d/l}{2}\big)!\,\big(\frac{N-d/l}{2}\big)!}", "(1.6)"),
                     (r"P(N,d)=\frac{C(N,d)}{2^N}", "(1.7)")):
        rows.add(VGroup(Tex(tex, font_size=36), Tex(num, font_size=24).set_color(GREY_B)).arrange(RIGHT, buff=0.35))
    rows.arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    rows.move_to(np.array([0, 1.35, 0])).align_to(np.array([-6.4, 0, 0]), LEFT)
    o = np.array([3.9, 0.3, 0])
    ch = D.bar_chart(DIST4, ux=0.45, uy=0.28, origin=o, width=1.6, labels=True, font_size=28)
    dl = VGroup(*[Tex(str(d), font_size=24).set_color(GREY_B).move_to(o + d * 0.45 * RIGHT + 0.25 * DOWN)
                  for d in range(-4, 5, 2)])
    t4 = Tex("N=4", font_size=34).move_to(o + 2.55 * UP)
    s4 = Tex(r"1+4+6+4+1=16=2^4", font_size=32).move_to(o + 0.85 * DOWN)
    r1 = Tex(r"C(10,0)=\frac{10!}{5!\,5!}=252", font_size=40).move_to(np.array([-3.4, -1.95, 0]))
    r2 = Tex(r"P(10,0)=\frac{252}{1024}=\frac{63}{256}\approx0.246", font_size=40).move_to(np.array([3.3, -1.95, 0]))
    assert_in_frame(rows, ch, dl, t4, s4, r1, r2)
    for a, b in ((rows, ch), (rows, t4), (rows, r1), (s4, r2), (r1, r2), (ch, t4)):
        assert_apart(a, b)
    self.say("換成位置：向右 ", "$n$", " 步、向左 ", "$N-n$", " 步，終點 ", "$d$")
    self.play(Write(rows[0]), run_time=1.3)
    self.wait(3.0)
    self.say("走法數 ", "$C(N,d)$", "：把 ", "$n$", " 代進二項式係數")
    self.play(Write(rows[1]), run_time=1.6)
    self.wait(3.5)
    self.say("機率：走法數除以全部 ", "$2^N$", " 種")
    self.play(Write(rows[2]), run_time=1.3)
    self.wait(3.0)
    self.say("$N=4$", "：", "$C(4,d)$", " 依序是 1, 4, 6, 4, 1")
    self.play(FadeIn(t4), LaggedStart(*[GrowFromEdge(b, DOWN) for b in ch[0]], lag_ratio=0.15), FadeIn(dl),
              run_time=1.6)
    self.play(LaggedStart(*[FadeIn(t) for t in ch[1]], lag_ratio=0.1))
    self.play(Write(s4))
    self.wait(3.0)
    self.say("揭曉：走 10 步回到原點的路線有")
    self.play(Write(r1), run_time=1.6)
    self.wait(2.5)
    self.say("機率 ", "$252/1024=63/256$", "，大約四分之一")
    self.play(Write(r2), run_time=1.6)
    self.play(ShowCreation(SurroundingRectangle(r1[-3:], buff=0.08).set_stroke(YELLOW, 2)))
    self.wait(3.0)


# ---------------------------------------------------------------- 5. 練習 1.1：偏的硬幣
def biased(self):
    top = Tex(r"P(n)=\binom{N}{n}p^n(1-p)^{N-n}", font_size=42).move_to(np.array([-1.4, 2.45, 0]))
    tag = zh("練習 1.1", font_size=26, color=GREY_B).next_to(top, RIGHT, buff=0.5)
    o = np.array([-2.6, -1.85, 0])
    ux, uy = 1.25, 4.4
    fair = D.bar_chart({n: float(P_FAIR[n]) for n in range(4)}, ux=ux, uy=uy, origin=o, width=0.66,
                       color=GREY_B, fill_opacity=0.0)
    bias = D.bar_chart({n: float(P_BIASED[n]) for n in range(4)}, ux=ux, uy=uy, origin=o, width=0.46,
                       color=D.HEAD_COLOR, fill_opacity=0.8)
    labs = VGroup(*[Tex(rf"\tfrac{{{k}}}{{27}}", font_size=34).set_color(D.HEAD_COLOR)
                    .next_to(max(bias[0][n], fair[0][n], key=lambda b: b.get_top()[1]), UP, buff=0.1)
                    for n, k in enumerate((1, 6, 12, 8))])
    nl = VGroup(*[Tex(f"n={n}", font_size=26).set_color(GREY_B).move_to(o + n * ux * RIGHT + 0.3 * DOWN)
                  for n in range(4)])
    leg = VGroup(
        VGroup(Square(0.25).set_fill(D.HEAD_COLOR, 0.8).set_stroke(D.HEAD_COLOR, 1),
               Tex(r"p=\tfrac23", font_size=30)).arrange(RIGHT, buff=0.2),
        VGroup(Square(0.25).set_stroke(GREY_B, 2),
               Tex(r"p=\tfrac12", font_size=30)).arrange(RIGHT, buff=0.2)).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    leg.move_to(np.array([4.6, 0.9, 0]))
    assert_apart(leg[0], leg[1])
    tot = Tex(r"\frac{1+6+12+8}{27}=1", font_size=38).move_to(np.array([4.6, -0.7, 0]))
    assert_in_frame(top, tag, fair, bias, labs, nl, leg, tot)
    for a, b in ((top, labs), (leg, labs), (tot, labs), (leg, tot), (leg, bias), (tot, bias)):
        assert_apart(a, b)
    self.say("練習 1.1：硬幣不公平，正面的機率是 ", "$p$")
    self.play(Write(top), FadeIn(tag), run_time=1.5)
    self.wait(3.0)
    self.say("$N=3$", "，", "$p=2/3$", "：灰框是公平硬幣")
    self.play(FadeIn(nl), LaggedStart(*[ShowCreation(b) for b in fair[0]], lag_ratio=0.15), FadeIn(leg[1]))
    self.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bias[0]], lag_ratio=0.15), FadeIn(leg[0]), run_time=1.5)
    self.play(LaggedStart(*[FadeIn(t) for t in labs], lag_ratio=0.15))
    self.wait(2.0)
    self.say("機率依序 ", "$1/27,\\ 6/27,\\ 12/27,\\ 8/27$", "：往正面那邊偏")
    self.wait(3.0)
    self.play(Write(tot))
    self.wait(3.0)


class W102CountingWalks(ChuScene):
    def construct(self):
        self.show_header("W102　數路線：二項式係數")
        check_math()
        ask(self)
        self.clear_stage()
        three_coins(self)
        self.clear_stage()
        shuffle(self)
        self.clear_stage()
        positions_part(self)
        self.clear_stage()
        biased(self)
        self.say("走法數＝二項式係數；機率＝走法數除以 ", "$2^N$", wait=4.0)
        self.wait(1.0)
