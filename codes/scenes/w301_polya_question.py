"""
W301 醉漢回得了家，醉鳥不一定（§3.1、§3.4 開頭，pp. 51, 55）

問：醉漢在街上亂走，一定走得回家嗎？醉鳥在天上亂飛呢？（字幕引角谷靜夫的話）
  1. 小例子：一維 4 步的 16 條路線；第 1–4 步碰過 0 的 10 條（黃），從沒回家的 6 條（灰）；16 − 6 = 10。
  2. 先猜再揭曉：一、二、三維各 1000 人走 1000 步，回過家的人數：100 步 930、579、300 → 1000 步 977、676、325；
     旁邊的理論值 1 − P_N（N = 100：0.9204、0.5832、0.3119；N = 1000：0.9748、0.6782、0.3315）。
  3. Pólya 定理：走無窮多步，回家的機率 一維 1、二維 1、三維（簡單立方）1 − 1/1.516386 = 0.3405。
  4. 本章地圖：W302 → W303 → W304 → W305 → W306 → W307–W308。
收束：維度決定命運。畫面上的數字都由 utils.walk_rec 算出（三個模擬重算一次），check_math() 核對。
"""
import sys
from fractions import Fraction
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

KINDS = ("1d", "sq", "sc")
DIMS = ("一維", "二維", "三維")
LATS = ("直線", "棋盤", "立方格子")
DIM_COLORS = (BLUE_C, TEAL_C, ORANGE)
WALKS4 = R.walks_1d(4)
DRUNK = [1, 1, -1, 1, -1, -1, -1, 1]                       # 開頭醉漢的 8 步：第 6、8 步回到家
BIRD = [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1), (2, 1, 1)]   # 醉鳥越飛越遠
MAP = [("W302", "新的生成函數"), ("W303", "一維"), ("W304", "二維"),
       ("W305", "維度的分水嶺"), ("W306", "三種立方格子"), ("W307–W308", "走過幾個點")]
UX4, UY4 = 0.25, 0.12                                       # 16 張小走法圖
BAR_ORG, BAR_UX, BAR_UY = np.array([-5.0, -1.85, 0]), 1.9, 0.0045
_P = {}


def returned(pos):
    """第 1..N 步碰過 0（回過家）。"""
    return any(x == 0 for x in pos[1:])


def first_return(steps):
    x = 0
    for k, st in enumerate(steps, 1):
        x += st
        if x == 0:
            return k
    return None


def theory(kind, N):
    """1 − P_N：N 步內回過家的精確機率。"""
    if kind not in _P:
        _P[kind] = R.no_return_probs(kind, 1000)
    return 1 - _P[kind][N]


def check_math():
    assert len(WALKS4) == 16 == 2 ** 4 and all(len(p) == 5 for p in WALKS4)
    assert sum(map(returned, WALKS4)) == 10 == 16 - R.never_return_1d(4) and R.never_return_1d(4) == 6
    assert Fraction(10, 16) == Fraction(5, 8) and f"{10 / 16:.3f}" == "0.625"
    assert first_return(DRUNK) == 6 and sum(DRUNK) == 0
    for kind in KINDS:                                               # 模擬重算一次（固定種子）
        fr = R.simulate_first_returns(kind, R.SIM_M, R.SIM_N, R.SIM_SEED[kind])
        assert len(fr) == R.SIM_M == 1000 and R.SIM_N == 1000
        assert sum(t is not None for t in fr) == R.SIM_RETURNED[kind]
        assert sum(t is not None and t <= 100 for t in fr) == R.SIM_RETURNED_100[kind]
    assert [R.SIM_RETURNED_100[k] for k in KINDS] == [930, 579, 300]
    assert [R.SIM_RETURNED[k] for k in KINDS] == [977, 676, 325]
    assert [f"{theory(k, 100):.4f}" for k in KINDS] == ["0.9204", "0.5832", "0.3119"]
    assert [f"{theory(k, 1000):.4f}" for k in KINDS] == ["0.9748", "0.6782", "0.3315"]
    assert round(1000 * theory("1d", 100)) == 920
    assert f"{R.G_C['sc']:.6f}" == "1.516386" and f"{1 - R.escape('sc'):.4f}" == "0.3405"
    assert abs((1 - R.escape("sc")) - (1 - 1 / 1.516386)) < 1e-6
    assert round(100 * R.escape("sc")) == 66


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def house(p, w=0.36, h=0.24, roof=0.2, color=YELLOW):
    p = np.array(p, dtype=float)
    pts = [p + w / 2 * LEFT, p + w / 2 * RIGHT, p + w / 2 * RIGHT + h * UP,
           p + (h + roof) * UP, p + w / 2 * LEFT + h * UP]
    return Polygon(*pts).set_fill(color, 0.85).set_stroke(WHITE, 1.5)


def mini(pos, center):
    """一張小走法圖：VGroup(座標軸, 路線)；另外回傳「回到 0」的黃點。"""
    org = np.array(center) + 2.45 * UX4 * LEFT
    ax = D.walk_axes(4, 4, ux=UX4, uy=UY4, origin=org, x_every=4, y_every=4, font_size=14)
    path = D.walk_path(pos, ux=UX4, uy=UY4, origin=org, color=WHITE, stroke_width=3, dot_radius=0.035)
    hits = VGroup(*[Dot(D.path_point(k, 0, UX4, UY4, org), radius=0.065).set_fill(YELLOW, 1)
                    for k in range(1, 5) if pos[k] == 0])
    return VGroup(VGroup(ax[0], ax[1]), path), hits


def bars(values):
    ch = D.bar_chart(dict(enumerate(values)), ux=BAR_UX, uy=BAR_UY, origin=BAR_ORG, labels=True, font_size=30)
    for k, b in enumerate(ch[0]):
        b.set_fill(DIM_COLORS[k], 0.75).set_stroke(DIM_COLORS[k], 1.5)
    return ch


def map_box(code, title, color=BLUE_C):
    rect = RoundedRectangle(width=3.4, height=1.2, corner_radius=0.15).set_stroke(color, 2.5).set_fill(color, 0.12)
    txt = VGroup(zh(code, font_size=20, color=GREY_B), zh(title, font_size=28)).arrange(DOWN, buff=0.12)
    txt.move_to(rect)
    return VGroup(rect, txt)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nl = D.number_line(5, unit=0.55).move_to(np.array([-3.4, 0.5, 0]))
    home = house(D.line_point(nl, 0) + 0.3 * UP)
    man = D.walker_dot().move_to(D.line_point(nl, 0))
    l1 = zh("醉漢：在一條街上（一維）", font_size=26)
    org3 = np.array([3.4, 0.7, 0])
    cube = L.lattice3d("sc", unit=0.95, origin=org3)
    pts = [L.iso_point(*q, unit=0.95, origin=org3) for q in BIRD]
    flight = VMobject().set_points_as_corners(pts).set_stroke(ORANGE, 4)
    bird = Dot(pts[0], radius=0.1).set_fill(ORANGE, 1)
    qm = Tex("?", font_size=60).set_color(ORANGE).next_to(pts[-1], UR, buff=0.1)
    l2 = zh("醉鳥：在立體空間裡（三維）", font_size=26)
    l1.move_to(np.array([nl.get_center()[0], -1.6, 0]))
    l2.move_to(np.array([cube.get_center()[0], -1.6, 0]))
    assert_in_frame(nl, home, l1, cube, flight, qm, l2)
    assert_apart(VGroup(nl, home), VGroup(cube, flight, qm))
    assert_apart(l1, l2)
    assert_apart(VGroup(cube, flight), l2)

    s.say("問：醉漢在街上亂走，一定走得回家嗎？", wait=0.3)
    s.play(FadeIn(nl), FadeIn(home), FadeIn(l1), run_time=0.8)
    s.play(FadeIn(man), run_time=0.4)
    x = 0
    for st in DRUNK:
        x += st
        man.generate_target()
        man.target.move_to(D.line_point(nl, x))
        s.play(MoveToTarget(man, path_arc=-st * 0.7 * PI), run_time=0.4)
        if x == 0:
            s.play(Indicate(home), run_time=0.6)
    s.say("走著走著，他回到家了，而且不只一次", wait=2.0)
    s.say("醉鳥在天上亂飛呢？", wait=0.3)
    s.play(FadeIn(cube), FadeIn(l2), run_time=1.0)
    s.play(FadeIn(bird), run_time=0.3)
    s.play(ShowCreation(flight), MoveAlongPath(bird, flight.copy()), run_time=2.5)
    s.play(FadeIn(qm), run_time=0.5)
    s.wait(1.0)
    s.say("角谷靜夫：「醉漢終究會走回家，醉鳥卻可能永遠迷路。」", wait=4.5)
    s.say("真的嗎？先從最簡單的一維看起", wait=2.0)
    s.clear_stage()


def part_example(s):
    ret = [returned(p) for p in WALKS4]
    minis, hits = VGroup(), VGroup()
    for k, pos in enumerate(WALKS4):
        i, j = divmod(k, 4)
        m, h = mini(pos, np.array([-5.9 + 1.55 * j, 2.3 - 1.5 * i, 0]))
        minis.add(m)
        hits.add(h)
    n_ret, n_never = sum(ret), R.never_return_1d(4)
    t1 = mixed("一維走 ", "$4$", " 步：共 ", "$2^4=16$", " 條路線", font_size=28)
    t2 = mixed("回過家（第 1–4 步碰過 ", "$0$", "）：", f"${n_ret}$", " 條", font_size=26, color=YELLOW)
    t3 = mixed("從沒回家：", f"${n_never}$", " 條", font_size=26, color=GREY_B)
    t4 = Tex(rf"16-{n_never}={16 - n_never}", font_size=40)
    t5 = mixed("4 步內回過家的比例：", rf"$\frac{{{n_ret}}}{{16}}={n_ret / 16:.3f}$", font_size=24)
    col = VGroup(t1, t2, t3, t4, t5).arrange(DOWN, buff=0.5, aligned_edge=LEFT)
    fit(col, 6.3)
    col.move_to(np.array([0.35 + col.get_width() / 2, 0.2, 0]))
    assert_in_frame(minis, hits, col)
    assert_apart(minis, col)

    s.say("小例子：一維走 4 步，一共 16 條路線", wait=0.3)
    s.play(LaggedStart(*[FadeIn(m) for m in minis], lag_ratio=0.1), run_time=2.5)
    s.play(FadeIn(t1), run_time=0.6)
    s.wait(1.5)
    s.say("橫軸就是「家」：第 1 到 4 步碰過 0 的，標黃", wait=0.3)
    s.play(*[m[1].animate.set_color(YELLOW) for m, r in zip(minis, ret) if r],
           *[FadeIn(h) for h, r in zip(hits, ret) if r], run_time=1.2)
    s.play(FadeIn(t2), run_time=0.6)
    s.wait(2.5)
    s.say("從沒碰過 0 的，標灰", wait=0.3)
    s.play(*[m[1].animate.set_color(GREY_C) for m, r in zip(minis, ret) if not r], run_time=1.0)
    s.play(FadeIn(t3), run_time=0.6)
    s.wait(1.5)
    s.play(FadeIn(t4), run_time=0.6)
    s.say("16 條裡有 10 條回過家", wait=1.0)
    s.play(FadeIn(t5), run_time=0.6)
    s.wait(2.0)
    s.say("走得越久，回家的機會越大嗎？換成二維、三維呢？", wait=3.0)
    s.clear_stage()


def part_guess(s):
    base = Line(BAR_ORG + 1.0 * LEFT, BAR_ORG + (2 * BAR_UX + 1.0) * RIGHT).set_stroke(GREY_B, 2)
    dlab = VGroup()
    for k in range(3):
        g = VGroup(zh(DIMS[k], font_size=26, color=DIM_COLORS[k]), zh(LATS[k], font_size=20, color=GREY_B))
        g.arrange(DOWN, buff=0.08).next_to(BAR_ORG + k * BAR_UX * RIGHT, DOWN, buff=0.18)
        dlab.add(g)
    qs = VGroup(*[Tex("?", font_size=80).move_to(BAR_ORG + k * BAR_UX * RIGHT + 1.2 * UP) for k in range(3)])
    ch100 = bars([R.SIM_RETURNED_100[k] for k in KINDS])
    ch1000 = bars([R.SIM_RETURNED[k] for k in KINDS])
    hdr = [zh(t, font_size=22, color=YELLOW) for t in ("維度", "100 步", "理論", "1000 步", "理論")]
    rows = [hdr] + [[zh(DIMS[k], font_size=24, color=DIM_COLORS[k]), str(R.SIM_RETURNED_100[kd]),
                     f"{theory(kd, 100):.4f}", str(R.SIM_RETURNED[kd]), f"{theory(kd, 1000):.4f}"]
                    for k, kd in enumerate(KINDS)]
    tab = D.tex_table(rows, [1.0, 1.15, 1.35, 1.3, 1.35], row_h=0.62, font_size=28)
    tab.move_to(np.array([0.35 + tab.get_width() / 2, 1.65, 0]))
    cell = lambda i, j: D.table_cell(tab, i, j, 5)
    n1 = zh("模擬：每一維 1000 人，各走 1000 步（固定亂數種子）", font_size=20, color=GREY_B)
    n2 = zh("理論：精確的機率（乘上 1000 就是預期的人數）", font_size=20, color=GREY_B)
    notes = VGroup(n1, n2).arrange(DOWN, buff=0.15, aligned_edge=LEFT)
    notes.next_to(tab, DOWN, buff=0.3).align_to(tab, LEFT)
    n3 = zh("三維：走了 1000 步，也只有約三分之一回過家", font_size=24, color=ORANGE)
    n3.next_to(notes, DOWN, buff=0.55).align_to(tab, LEFT)
    assert_in_frame(base, dlab, qs, ch100, ch1000, tab, notes, n3)
    assert_apart(VGroup(base, dlab, ch1000), VGroup(tab, notes, n3))

    s.say("先猜：每一維 1000 人各走 1000 步，回過家的有幾人？", wait=0.3)
    s.play(ShowCreation(base), FadeIn(dlab), FadeIn(qs), run_time=1.0)
    s.wait(4.0)
    s.say("先看走到第 100 步的時候", wait=0.3)
    s.play(FadeOut(qs), *[GrowFromEdge(b, DOWN) for b in ch100[0]], run_time=1.5)
    s.play(FadeIn(ch100[1]), run_time=0.6)
    s.play(FadeIn(tab[1]), *[FadeIn(cell(i, j)) for i in range(4) for j in range(3)], run_time=1.0)
    s.play(FadeIn(notes), run_time=0.6)
    s.say("理論值是精確的機率：一維 ", "$0.9204$", "，約 920 人", wait=3.5)
    s.say("再走到第 1000 步：", wait=0.3)
    s.play(*[Transform(a, b) for a, b in zip(ch100[0], ch1000[0])],
           *[Transform(a, b) for a, b in zip(ch100[1], ch1000[1])], run_time=2.0)
    s.play(*[FadeIn(cell(i, j)) for i in range(4) for j in (3, 4)], run_time=1.0)
    s.wait(2.0)
    s.say("一維幾乎全部回家；二維約三分之二；三維只有三分之一", wait=3.5)
    s.play(Indicate(VGroup(*[cell(3, j) for j in range(5)]), color=ORANGE), FadeIn(n3), run_time=1.2)
    s.say("三維的長條還會長到 1000 嗎？", wait=3.0)
    s.clear_stage()


def part_polya(s):
    title = VGroup(Tex(r"\text{P\'{o}lya}", font_size=60), zh("定理", font_size=42)).arrange(RIGHT, buff=0.18)
    title.set_color(YELLOW).move_to(np.array([0, 2.55, 0]))
    sub = zh("走無窮多步，回家的機率", font_size=26, color=GREY_B).next_to(title, DOWN, buff=0.22)
    esc = R.escape("sc")
    rows = [[zh("維度", font_size=26, color=YELLOW), zh("格子", font_size=26, color=YELLOW),
             zh("回家的機率", font_size=26, color=YELLOW)],
            [zh("一維", font_size=28), zh("直線", font_size=28), "1"],
            [zh("二維", font_size=28), zh("正方格子", font_size=28), "1"],
            [zh("三維", font_size=28, color=ORANGE), zh("簡單立方", font_size=28, color=ORANGE),
             rf"1-1/{R.G_C['sc']:.6f}={1 - esc:.4f}"]]
    tab = D.tex_table(rows, [1.6, 2.4, 5.4], row_h=0.8, font_size=38)
    tab.move_to(np.array([0, -0.75, 0]))
    cell = lambda i, j: D.table_cell(tab, i, j, 3)
    cell(3, 2).set_color(ORANGE)
    box = SurroundingRectangle(cell(3, 2), buff=0.12).set_stroke(ORANGE, 3)
    assert_in_frame(title, sub, tab, box)
    assert_apart(sub, tab)

    s.say(r"$\text{P\'{o}lya}$", " 定理：走無窮多步，回家的機率是多少？", wait=0.3)
    s.play(FadeIn(title, shift=0.2 * DOWN), FadeIn(sub), run_time=1.0)
    s.play(FadeIn(tab[1]), *[FadeIn(cell(0, j)) for j in range(3)], run_time=0.8)
    s.say("一維：機率 1，一定回家", wait=0.3)
    s.play(*[FadeIn(cell(1, j)) for j in range(3)], run_time=0.8)
    s.wait(1.5)
    s.say("二維：也是 1，一定回家", wait=0.3)
    s.play(*[FadeIn(cell(2, j)) for j in range(3)], run_time=0.8)
    s.wait(1.5)
    s.say("三維：只有 ", f"${1 - esc:.4f}$", "，大約三分之二的人永遠回不了家", wait=0.3)
    s.play(*[FadeIn(cell(3, j)) for j in range(3)], run_time=0.8)
    s.play(ShowCreation(box), run_time=0.8)
    s.wait(3.0)
    s.say("其中 ", f"${R.G_C['sc']:.6f}$", " 是三維的 ", "$G(z_c;0)$", "：W305、W306 會算出來", wait=3.5)
    s.clear_stage()


def part_map(s):
    head = zh("第 3 章地圖", font_size=32, color=YELLOW).move_to(np.array([0, 2.75, 0]))
    pos = [(-4.4, 1.2), (0, 1.2), (4.4, 1.2), (4.4, -1.4), (0, -1.4), (-4.4, -1.4)]
    boxes = [map_box(c, t).move_to(np.array([x, y, 0])) for (c, t), (x, y) in zip(MAP, pos)]
    arrows = [Arrow(boxes[0].get_right(), boxes[1].get_left(), buff=0.1),
              Arrow(boxes[1].get_right(), boxes[2].get_left(), buff=0.1),
              Arrow(boxes[2].get_bottom(), boxes[3].get_top(), buff=0.1),
              Arrow(boxes[3].get_left(), boxes[4].get_right(), buff=0.1),
              Arrow(boxes[4].get_left(), boxes[5].get_right(), buff=0.1)]
    for a in arrows:
        a.set_color(GREY_B)
    assert_in_frame(head, *boxes)
    for a, b in zip(boxes, boxes[1:]):
        assert_apart(a, b, 0.3)

    s.say("這一章的路線圖：", wait=0.3)
    s.play(FadeIn(head), run_time=0.6)
    for k, b in enumerate(boxes):
        anims = [FadeIn(b)]
        if k:
            anims.append(GrowArrow(arrows[k - 1]))
        s.play(*anims, run_time=0.7)
    s.say("先用新的生成函數數「到訪次數」，再看一維、二維", wait=3.5)
    s.say("然後找出維度的分水嶺、比較三種立方格子", wait=3.0)
    s.say("最後：走了 ", "$n$", " 步，平均踩過幾個不同的點？", wait=3.0)


def part_close(s):
    s.say("維度決定命運：一、二維遲早回家", wait=2.5)
    s.say(f"三維有 {round(100 * R.escape('sc'))}% 的人永遠回不了家", wait=3.5)


class W301PolyaQuestion(ChuScene):
    def construct(self):
        self.show_header("W301　醉漢回得了家，醉鳥不一定")
        check_math()
        part_question(self)
        part_example(self)
        part_guess(self)
        part_polya(self)
        part_map(self)
        part_close(self)
        self.wait(1.0)
