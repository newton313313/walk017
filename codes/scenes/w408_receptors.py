"""
W408 布滿受器的球（Berg–Purcell）＋吸收邊界＋第 4 章回顧（§4.2.6、§4.3，pp. 91–94，(4:S-1)–(4:S-3)）

問：細胞表面只有一小部分是受器（吸收點）。吃營養的速率會打幾折？
  1. 書上的論證（p. 92）：電位 ≈ Q/(na)（自己）＋ Q/r0（其他受器）；na ≫ r0 時和整顆球一樣好。筆誤框 na ≪ r0 → na² ≪ r0²。
  2. 書外（Berg & Purcell 1977）：J/Jmax = na/(na + πr0)（電容串聯；BC.berg_purcell ＝ BC.berg_purcell_series）；
     r0 = 1 µm、a = 1 nm：半飽和 n = 3142，只占表面 0.079%。
  3. §4.3 吸收邊界：虛擬格點 G(z; x, −Δ) = 0；鏡像點 −x − 2；Catalan 數 1, 1, 2, 5, 14, 42（BC.half_counts ＝ half_counts_image ＝ catalan）；
     筆誤框 (4:S-2)。
  4. 第 4 章回顧表。
畫面上的數字都由 utils.walk_bc 算出，check_math() 核對。
"""
import itertools
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
REC_COLOR = RED_B
BP_NS = (100, 1000, 3142, 10 ** 4, 10 ** 5)
BP_TEX = ("100", "1000", "3142", "10^4", "10^5")
KS = tuple(range(6))
EXT_COLOR = ORANGE                       # 書外補充的標籤色


# ---------------------------------------------------------------- 數學
def dyck_paths(k):
    """從 0 出發、2k 步回到 0、從不踩到 −1 的路線（位置 list）。"""
    out = []
    for steps in itertools.product((1, -1), repeat=2 * k):
        pos = [0]
        for st in steps:
            pos.append(pos[-1] + st)
        if pos[-1] == 0 and min(pos) >= 0:
            out.append(pos)
    return out


@lru_cache(maxsize=None)
def numbers():
    a, r0 = BC.BP_EXAMPLE["a"], BC.BP_EXAMPLE["r0"]
    nb = {"bp": [BC.berg_purcell(n, a, r0) for n in BP_NS], "half": BC.bp_half(a, r0)}
    nb["area_pct"] = 100 * BC.bp_area_fraction(nb["half"], a, r0)
    nb["free"] = [math.comb(2 * k, k) for k in KS]
    nb["img"] = [math.comb(2 * k, k + 1) for k in KS]
    nb["cat"] = [BC.catalan(k) for k in KS]
    return nb


def check_math():
    nb = numbers()
    a, r0 = BC.BP_EXAMPLE["a"], BC.BP_EXAMPLE["r0"]
    assert (a, r0) == (1e-9, 1e-6)
    assert [f"{v:.3f}" for v in nb["bp"]] == ["0.031", "0.241", "0.500", "0.761", "0.970"]
    assert all(abs(BC.berg_purcell(n, a, r0) - BC.berg_purcell_series(n, a, r0)) < 1e-12 for n in BP_NS)
    assert abs(nb["half"] - math.pi * r0 / a) < 1e-9 and str(round(nb["half"])) == "3142" == BP_TEX[2]
    assert abs(BC.berg_purcell(nb["half"], a, r0) - 0.5) < 1e-12
    assert f"{nb['area_pct']:.3f}" == "0.079"
    assert [str(n) for n in BP_NS[:3]] == list(BP_TEX[:3]) and BP_NS[3:] == (10 ** 4, 10 ** 5)
    assert nb["free"] == [1, 2, 6, 20, 70, 252] and nb["img"] == [0, 1, 4, 15, 56, 210]
    assert nb["cat"] == [1, 1, 2, 5, 14, 42] == [f - i for f, i in zip(nb["free"], nb["img"])]
    assert nb["cat"] == [BC.half_counts(2 * k, 0, 0) for k in KS] == [BC.half_counts_image(2 * k, 0, 0) for k in KS]
    assert len(dyck_paths(3)) == BC.catalan(3) == 5
    assert all(BC.half_counts(N, x, 0) == BC.half_counts_image(N, x, 0) for N in range(9) for x in range(4))


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


def receptor_cell(center, radius, n=20):
    """細胞剖面＋圓周上 n 個受器 → VGroup(cell, receptors)。"""
    cell = Circle(radius=radius).move_to(center).set_fill(CELL_FILL, 0.55).set_stroke(CELL_STROKE, 3)
    recs = VGroup(*[Dot(center + radius * np.array([math.cos(TAU * k / n + 0.1), math.sin(TAU * k / n + 0.1), 0]),
                        radius=0.065).set_color(REC_COLOR) for k in range(n)])
    return VGroup(cell, recs)


def site_line(xs, unit, origin):
    return {x: np.array(origin, dtype=float) + x * unit * RIGHT for x in xs}


def mini_path(pos, ux, uy, origin, color=YELLOW):
    org = np.array(origin, dtype=float)
    pts = [org + k * ux * RIGHT + p * uy * UP for k, p in enumerate(pos)]
    base = Line(org, org + (len(pos) - 1) * ux * RIGHT).set_stroke(GREY_C, 1.5)
    cliff = Line(org + uy * DOWN, org + (len(pos) - 1) * ux * RIGHT + uy * DOWN).set_stroke(RED, 3)
    path = VMobject().set_points_as_corners(pts).set_stroke(color, 3)
    return VGroup(cliff, base, path)


# ---------------------------------------------------------------- 各段
def part_question(s):
    c = np.array([-3.4, -0.1, 0])
    rc = receptor_cell(c, 1.7)
    rad = Line(c, c + 1.7 * np.array([math.cos(-0.6), math.sin(-0.6), 0])).set_stroke(WHITE, 2)
    rl = tx("r_0", 28).next_to(rad.get_center(), UP, buff=0.12)
    lab = zh("受器（吸收點）", font_size=24, color=REC_COLOR).next_to(rc[1][4], UP + RIGHT, buff=0.12)
    rng = np.random.default_rng(5)
    food = VGroup()
    while len(food) < 45:
        p = np.array([rng.uniform(-6.4, 0.2), rng.uniform(-2.7, 2.6), 0])
        if np.linalg.norm(p - c) > 2.0:
            food.add(Dot(p, radius=0.04).set_color(GOLD))
    q = VGroup(mixed("吃營養的速率：", font_size=30), mixed("$\\times\\ ?$", font_size=40, color=YELLOW))
    q.arrange(DOWN, buff=0.4)
    at(q, 3.6, 0.3)
    assert_in_frame(rc, food, lab, q)
    assert_apart(VGroup(rc, food, lab), q)
    s.play(FadeIn(rc[0]), FadeIn(food), ShowCreation(rad), FadeIn(rl), run_time=1.0)
    s.say("W407 把整顆細胞當成吸收球", wait=2.0)
    s.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in rc[1]], lag_ratio=0.05), FadeIn(lab), run_time=1.5)
    s.say("問：其實只有表面的一小部分是受器（吸收點）", wait=2.5)
    s.play(FadeIn(q), run_time=0.6)
    s.say("吃營養的速率會打幾折？（先猜）", wait=3.0)
    s.clear_stage()


def part_book(s):
    c = np.array([-4.1, 0.1, 0])
    rc = receptor_cell(c, 1.6, n=16)
    qlab = mixed("每個受器帶 ", "$Q/n$", font_size=24, color=REC_COLOR).next_to(rc, DOWN, buff=0.3)
    a_lab = mixed("$n$", " 個、半徑 ", "$a$", font_size=24, color=REC_COLOR).next_to(qlab, DOWN, buff=0.15)
    phi = tx(r"\phi\approx", 36)
    t1, t2 = tx(r"\frac{Q}{na}", 36), tx(r"\frac{Q}{r_0}", 36)
    l1 = zh("自己", font_size=22, color=GREY_A).next_to(t1, DOWN, buff=0.15)
    l2 = zh("其他受器（像均勻帶電球）", font_size=22, color=GREY_A).next_to(t2, DOWN, buff=0.15)
    plus = tx("+", 36)
    row = VGroup(phi, VGroup(t1, l1), plus, VGroup(t2, l2)).arrange(RIGHT, buff=0.3, aligned_edge=UP)
    phi.match_y(t1)
    plus.match_y(t1)
    b0 = zh("書上的論證（p. 92）：把受器當成小導體", font_size=24, color=GREY_A)
    b3 = mixed("$na\\gg r_0$", "：第二項主導 ", "$\\Rightarrow$", " 和整顆球一樣好", font_size=26)
    box = D.erratum_box("p. 92–93", r"na\ll r_0,\ \ a\ll r/\sqrt{n}", r"na^2\ll r_0^2,\ \ a\ll r_0/\sqrt{n}", width=6.6)
    right = VGroup(b0, row, b3, box).arrange(DOWN, buff=0.5)
    fit(right, 7.0)
    at(right, 2.95, 0.05)
    left = VGroup(rc, qlab, a_lab)
    assert_in_frame(left, right)
    assert_apart(left, right)
    for k in range(3):
        assert_apart(right[k], right[k + 1])
    assert_apart(l1, l2)
    s.play(FadeIn(rc), FadeIn(b0), run_time=0.8)
    s.say("用靜電想：受器是 ", "$n$", " 個半徑 ", "$a$", " 的小導體，總電荷 ", "$Q$", wait=0.3)
    s.play(FadeIn(qlab), FadeIn(a_lab), run_time=0.7)
    s.wait(1.5)
    s.play(FadeIn(phi), FadeIn(t1), FadeIn(l1), run_time=0.8)
    s.say("某個受器上的電位：自己那份 ", "$Q/(na)$", wait=2.5)
    s.play(FadeIn(plus), FadeIn(t2), FadeIn(l2), run_time=0.8)
    s.say("其他受器加起來像均勻帶電的球：", "$Q/r_0$", wait=2.5)
    s.play(FadeIn(b3), run_time=0.8)
    s.say("$na\\gg r_0$", " 時第一項很小：電容和整顆球差不多", wait=3.0)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("「受器只占小面積」是 ", "$na^2\\ll r_0^2$", "，不是 ", "$na\\ll r_0$", wait=3.0)
    s.say("所以 ", "$na\\gg r_0$", " 和小面積並不矛盾", wait=2.5)
    s.clear_stage()


def part_berg_purcell(s):
    nb = numbers()
    a, r0 = BC.BP_EXAMPLE["a"], BC.BP_EXAMPLE["r0"]
    # 橫軸畫 u = log10(n) − 1.5（Axes 的縱軸在 u = 0，所以從 0 開始）
    ax = D.plot_axes((0, 4, 1), (0, 1.0, 0.25), width=4.8, height=2.6, x_ticks=[0.5, 1.5, 2.5, 3.5],
                     x_tick_labels=["10^2", "10^3", "10^4", "10^5"], y_ticks=[0.5, 1], font_size=20,
                     x_label="n", y_label=r"J/J_{\max}")
    ax.move_to(np.array([-3.5, 0.85, 0]))
    curve = ax.get_graph(lambda u: BC.berg_purcell(10 ** (u + 1.5), a, r0), x_range=(0, 4, 0.05)).set_stroke(YELLOW, 4)
    dots = VGroup(*[Dot(ax.c2p(math.log10(n) - 1.5, v), radius=0.07).set_color(ORANGE)
                    for n, v in zip(BP_NS, nb["bp"])])
    xh = math.log10(nb["half"]) - 1.5
    half = VGroup(DashedLine(ax.c2p(0, 0.5), ax.c2p(xh, 0.5)), DashedLine(ax.c2p(xh, 0.5), ax.c2p(xh, 0)))
    half.set_stroke(GREY_B, 2)
    tag = zh("書外補充：Berg & Purcell 1977", font_size=24, color=EXT_COLOR).next_to(ax, UP, buff=0.15)
    p2 = tx(r"\frac{J}{J_{\max}}=\frac{na}{na+\pi r_0}", 38, YELLOW)
    p3 = VGroup(tx(r"\frac{1}{C}=\frac{1}{r_0}+\frac{\pi}{na}", 32), zh("電容串聯", font_size=24, color=GREY_A)).arrange(RIGHT, buff=0.3)
    p4 = mixed("半飽和：", "$n=\\pi r_0/a$", font_size=26)
    right = VGroup(p2, p3, p4).arrange(DOWN, buff=0.45)
    at(right, 3.4, 1.15)
    p5 = at(mixed("例：", "$r_0=1\\,\\mu\\mathrm{m},\\ a=1\\,\\mathrm{nm}$", "，", "$n=%d$" % round(nb["half"]),
                  " 個受器只占表面 ", "$%.3f\\%%$" % nb["area_pct"], font_size=26), 0, -1.15)
    rows = [["n"] + list(BP_TEX), [r"J/J_{\max}"] + [f"{v:.3f}" for v in nb["bp"]]]
    tab = D.tex_table(rows, [1.8] + [1.25] * 5, row_h=0.55, font_size=28)
    at(tab, 0, -2.3)
    assert_in_frame(ax, tag, right, p5, tab)
    assert_apart(VGroup(ax, tag), right)
    for a_, b_ in ((ax, p5), (right, p5), (p5, tab)):
        assert_apart(a_, b_)
    assert_table(tab, 6)
    s.play(FadeIn(tag), ShowCreation(ax), run_time=1.0)
    s.say("書外補充：Berg 與 Purcell（1977）算得更細", wait=2.5)
    s.play(FadeIn(p2), run_time=0.8)
    s.play(ShowCreation(curve), run_time=1.5)
    s.say("吃營養的速率是整顆球的 ", "$na/(na+\\pi r_0)$", wait=3.0)
    s.play(FadeIn(p3), run_time=0.8)
    s.say("像兩個電容串聯：整顆球 ", "$r_0$", "，受器那層 ", "$na/\\pi$", wait=3.0)
    s.play(FadeIn(p4), ShowCreation(half), run_time=0.8)
    s.play(FadeIn(p5), run_time=0.8)
    s.say("1 微米的細胞、1 奈米的受器：3142 個就達一半", wait=3.0)
    s.say("這些受器只占表面的 0.079%", wait=2.5)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), FadeIn(dots), run_time=1.0)
    s.say("一萬個受器就有 76%：一小撮受器就夠用", wait=3.0)
    s.clear_stage()


def part_absorbing(s):
    nb = numbers()
    xs = list(range(-2, 7))
    unit, org = 0.85, np.array([-3.6, 2.35, 0])
    sp = site_line(xs, unit, org)
    base = Line(sp[0], sp[6]).set_stroke(GREY_C, 2)
    dots = VGroup(*[Dot(sp[x], radius=0.07).set_color(WHITE) for x in range(0, 7)])
    virt = Dot(sp[-1], radius=0.11).set_fill(BLACK, 1).set_stroke(RED, 3)
    imgd = Dot(sp[-2], radius=0.09).set_color(BLUE_C)
    nums = VGroup(*[tx(str(x), 22, GREY_A).next_to(sp[x], DOWN, buff=0.18) for x in xs])
    vlab = mixed("懸崖／虛擬格點 ", "$G=0$", font_size=20, color=RED).next_to(sp[-1], UP, buff=0.18)
    ilab = mixed("鏡像點 ", "$-x-2$", font_size=20, color=BLUE_C).next_to(nums[0], DOWN, buff=0.12)
    line = VGroup(base, dots, nums)
    s1 = eq_label(tx(r"G(z;x,y)=z\big[G(z;x,y-\Delta)+G(z;x,y+\Delta)\big]+\delta_{x,y}", 28), "(4:S-1)")
    s3 = VGroup(tx(r"G(z;x,-\Delta)=0", 30, RED), zh("（虛擬格點）", font_size=22, color=GREY_A),
                tx(r"\text{(4:S-3)}", 22, GREY_B)).arrange(RIGHT, buff=0.3)
    box = D.erratum_box("p. 93，(4:S-2)", r"G(z;0,y)=zG(z;x,\Delta)+\delta_{x,\Delta}",
                        r"G(z;x,0)=zG(z;x,\Delta)+\delta_{x,0}", width=10.5)
    lower = VGroup(s1, s3, box).arrange(DOWN, buff=0.45)
    fit(lower, 12.6)
    at(lower, 0, -0.85)
    assert_in_frame(line, virt, imgd, vlab, ilab, lower)
    assert_apart(VGroup(line, ilab, vlab), lower)
    for k in range(2):
        assert_apart(lower[k], lower[k + 1])
    s.play(ShowCreation(base), FadeIn(dots), FadeIn(nums[2:]), run_time=0.8)
    s.say("§4.3 補充：格點 ", "$n\\ge0$", "，懸崖在 ", "$-1$", wait=2.5)
    s.play(FadeIn(s1), run_time=0.8)
    s.say("一般的格點：左右鄰居各貢獻 ", "$z$", " 倍", wait=2.5)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("在 ", "$y=0$", " 少了左鄰居：書上這一式的下標寫錯了", wait=3.0)
    s.play(FadeIn(virt), FadeIn(nums[1]), FadeIn(vlab), FadeIn(s3), run_time=0.8)
    s.say("加一個虛擬格點 ", "$G(z;x,-\\Delta)=0$", "：整條遞迴就統一了", wait=3.0)
    s.play(FadeIn(imgd), FadeIn(nums[0]), FadeIn(ilab), run_time=0.8)
    s.say("對 ", "$-1$", " 鏡射：鏡像點在 ", "$-x-2$", wait=2.5)
    s.play(FadeOut(lower), run_time=0.6)

    rows = [["k"] + [str(k) for k in KS],
            ["C(2k;0)"] + [str(v) for v in nb["free"]],
            [r"C(2k;-2\to0)"] + [str(v) for v in nb["img"]],
            [zh("差", font_size=26, color=YELLOW)] + [str(v) for v in nb["cat"]]]
    tab = D.tex_table(rows, [3.0] + [1.0] * len(KS), row_h=0.52, font_size=28)
    at(tab, 0, 0.25)
    paths = VGroup(*[mini_path(p, 0.2, 0.2, ORIGIN) for p in dyck_paths(3)]).arrange(RIGHT, buff=0.45, aligned_edge=DOWN)
    at(paths, 0, -1.75)
    plab = mixed("$k=3$", "：", "$%d$" % BC.catalan(3), " 條從 0 回到 0、不碰懸崖的路線", font_size=24).next_to(paths, DOWN, buff=0.2)
    assert_in_frame(tab, paths, plab)
    assert_table(tab, len(KS) + 1)
    assert_apart(VGroup(line, ilab), tab)
    assert_apart(tab, paths)
    s.play(FadeIn(tab[0][:2 * (len(KS) + 1)]), ShowCreation(tab[1]), run_time=0.9)
    s.say("從 0 出發、", "$2k$", " 步回到 0：先不管懸崖", wait=2.0)
    s.play(FadeIn(tab[0][2 * (len(KS) + 1):3 * (len(KS) + 1)]), run_time=0.8)
    s.say("減掉從鏡像點 ", "$-2$", " 出發到 0 的走法", wait=2.5)
    s.play(FadeIn(tab[0][3 * (len(KS) + 1):]), run_time=0.8)
    s.say("剩下的正是 Catalan 數：1, 1, 2, 5, 14, 42", wait=3.0)
    s.play(LaggedStart(*[FadeIn(p) for p in paths], lag_ratio=0.2), FadeIn(plab), run_time=1.2)
    s.say("$k=3$", " 的 5 條：都待在懸崖上面", wait=2.5)
    s.clear_stage()


def part_recap(s):
    h = lambda t: zh(t, font_size=26, color=YELLOW)
    c = lambda *p: mixed(*p, font_size=26)
    rows = [[h("情況"), h("工具與結果"), h("場景")],
            [c("吸收牆"), c("鏡像"), tx(r"\text{W401, W406}", 28)],
            [c("有限區域"), c("正弦／", "$j_0$", " 模態、", "$e^{-\\pi^2N/2L^2}$"), tx(r"\text{W402--W404}", 28)],
            [c("穩定態"), c("Poisson：帳篷、", "$1/r$", "、", "$4\\pi Dr_0c_\\infty$"), tx(r"\text{W405--W407}", 28)]]
    tab = D.tex_table(rows, [2.6, 6.8, 2.8], row_h=0.85, font_size=28)
    at(tab, 0, 0.75)
    bottom = at(mixed("吸收牆＝接地導體", font_size=36, color=YELLOW), 0, -1.6)
    head = at(zh("第 4 章回顧", font_size=30, color=GREY_A), 0, 2.85)
    assert_in_frame(head, tab, bottom)
    assert_table(tab, 3)
    assert_apart(head, tab)
    assert_apart(tab, bottom)
    s.play(FadeIn(head), FadeIn(tab[0][:3]), ShowCreation(tab[1]), run_time=1.0)
    for i in range(1, 4):
        s.play(FadeIn(tab[0][3 * i:3 * i + 3], shift=0.2 * RIGHT), run_time=0.7)
        s.wait(0.8)
    s.say("第 4 章：牆、有限區域、穩定態", wait=2.5)
    s.play(FadeIn(bottom, shift=0.2 * UP), run_time=0.8)
    s.say("三件事背後是同一個類比", wait=2.0)


def part_close(s):
    s.say("吸收牆就是接地導體：鏡像、模態、電容", wait=3.0)
    s.say("靜電學的工具，全都能拿來數走法", wait=3.0)


class W408Receptors(ChuScene):
    def construct(self):
        self.show_header("W408　布滿受器的球＋吸收邊界＋第 4 章回顧")
        check_math()
        part_question(self)
        part_book(self)
        part_berg_purcell(self)
        part_absorbing(self)
        part_recap(self)
        part_close(self)
        self.wait(1.0)
