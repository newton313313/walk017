"""
W1001 關在盒子裡的高分子（§10.1，pp. 233–236，(10.1)，Fig. 10.1）

問：一條高分子關在小盒子裡、和一池單體平衡：盒子的牆會怎樣改變它能走的路？
  1. Fig. 10.1 示意（方框、長鏈、30 個單體；一個單體接上鏈端、一個蒸發回池裡）；(10.1) G(z; x, y) = Σ_N C(N; x, y) z^N；
     筆誤框 1（p. 235，z^n → z^N）：照書上 n → 0 時和 z 無關，4×4 環面永遠是 90677 條、⟨N⟩ = 0。
  2. 4×4 環面（週期邊界，對邊同記號）：'short' WNE、'mid' WSSEENNE（越界 2 次，攤開到無限平面，終點 (3, 1)）、
     'full' WSSEEENNNEEESSW（越界 3 次，走滿 16 格）。
  3. 表：環面 c_N 與無限平面 c_N（N = 0..15）：N ≤ 3 相同、N = 4 少 4 條（EEEE 繞一圈）、走滿 5712 條、共 90677 條；
     3×3、5×5 也在 N = L 少 4 條；16 × 5712 = 91392（四維 Gray 碼）。
  4. 三種 z：⟨N⟩ = 1.4400、4.8899、11.7409（最多 15）。
畫面上的數字都由 utils.walk_cond 算出（z_c 用 walk_crit.ZC_2D），check_math() 核對。
模組層級的畫圖函式：polymer_box、chevron、torus_marks、torus_box、walk_lines、lift_figure、count_table。
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_cond as WK
from utils import walk_crit as WC
from utils import lattice_draw as L
from utils import walk_draw as D

SIDE = 4                                   # 4×4 環面
EXAMPLES = {"short": (3, "WNE", 0), "mid": (8, "WSSEENNE", 2), "full": (15, "WSSEEENNNEEESSW", 3)}
C4 = (1, 4, 12, 36, 96, 264, 624, 1536, 3072, 6240, 9840, 15912, 17472, 19200, 10656, 5712)
CINF = (1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932, 881500, 2374444, 6416596)
Z_CASES = (0.2, WC.ZC_2D, 1.0)
MEAN_N = ("1.4400", "4.8899", "11.7409")
STEP_COLOR = YELLOW
WRAP_COLOR = "#F2A33A"                     # 越界的半段
MARK1_COLOR = TEAL_B                       # 左右邊的記號
MARK2_COLOR = "#C792EA"                    # 上下邊的記號
MONO_COLOR = TEAL_B
SAME_COLOR = "#7CE38B"
LOOP_COLOR = RED
LINK = 0.26                                # 一個單體的長度
CHAIN = ((-1.25, -1.45), (-0.55, -1.05), (-0.85, -0.25), (-1.45, 0.35), (-0.95, 1.05), (-0.15, 0.75),
         (0.25, 0.05), (0.95, -0.35), (1.35, 0.45), (0.85, 1.25), (0.05, 1.45))
LOOP_SEGS = (((0, 1), (1, 1), False), ((1, 1), (2, 1), False), ((2, 1), (3, 1), False),
             ((3, 1), (3.5, 1), True), ((-0.5, 1), (0, 1), True))          # EEEE：繞一圈回到起點


def f4(x):
    return f"{float(x):.4f}"


def check_math():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c4 = tuple(WK.torus_counts(SIDE))
        c3, c5 = WK.torus_counts(3), WK.torus_counts(5)
        means = [WK.torus_mean_N(SIDE, z) for z in Z_CASES]
        book = WK.book101_mean_N(SIDE, 0.3)
    assert c4 == C4 and tuple(WK.SAW_INF[:16]) == CINF and len(c4) == SIDE * SIDE
    assert sum(c4) == 90677 and book == 0
    assert c4[:4] == CINF[:4] and (c4[4], CINF[4]) == (96, 100) and CINF[4] - c4[4] == 4
    assert c4[-1] == 5712 and 16 * c4[-1] == WK.GRAY_Q4 == 91392
    assert max(range(SIDE * SIDE), key=lambda N: c4[N]) == 13          # 環面的 c_N 先多後少
    assert all(CINF[N] < CINF[N + 1] for N in range(15))
    assert (c3[3], WK.SAW_INF[3]) == (32, 36) and (c5[5], WK.SAW_INF[5]) == (280, 284)
    assert all(WK.torus_first_diff(n) == (n, 4) for n in (3, 4, 5))
    for name, (N, steps, wraps) in EXAMPLES.items():
        ex = WK.torus_example(name)
        assert (ex["N"], ex["steps"], ex["wraps"], tuple(ex["start"])) == (N, steps, wraps, (1, 1))
        segs = ex["segments"]
        assert len(segs) == N + wraps <= 18 and sum(1 for sg in segs if sg[2]) == 2 * wraps
        assert tuple(segs[0][0]) == (1, 1) and len(ex["lift"]) == N + 1
    assert len(set(WK.torus_example("full")["sites"])) == SIDE * SIDE
    assert tuple(WK.torus_example("mid")["lift"][-1]) == (3, 1)
    assert [f4(m) for m in means] == list(MEAN_N) and f4(WC.ZC_2D) == "0.3791"
    assert means[0] < means[1] < means[2] < SIDE * SIDE - 1 == 15
    assert f"{Z_CASES[0]:.1f}" == "0.2" and f"{Z_CASES[2]:.1f}" == "1.0"


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框 [x0,x1]×[y0,y1]（只縮不放），再置中在方框裡。"""
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def place_left(m, x0, y, x1=6.7):
    """左緣對齊 x0、中心高度 y；太寬就縮到 x1 為止。"""
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    return m.move_to(np.array([x0 + m.get_width() / 2, y, 0]))


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def unit(v):
    v = np.array(v, dtype=float)
    return v / np.linalg.norm(v)


def r2_points(count):
    """低差異的二維點列（R2 序列），每個座標在 [0, 1)。"""
    g = 1.32471795724474602596
    return [((0.5 + k / g) % 1.0, (0.5 + k / g ** 2) % 1.0) for k in range(1, count + 1)]


def seg_dist(p, a, b):
    p, a, b = (np.array(v, dtype=float) for v in (p, a, b))
    ab = b - a
    t = min(1.0, max(0.0, float(np.dot(p - a, ab) / np.dot(ab, ab))))
    return float(np.linalg.norm(p - a - t * ab))


# ---------------------------------------------------------------- 畫圖元件
def polymer_box(center, half=2.0, n_mono=30):
    """Fig. 10.1 示意：方框、長鏈（尾端多一節會蒸發）、散落的單體 → dict（無數字）。"""
    c = np.array(center, dtype=float)
    P = [np.array(p, dtype=float) for p in CHAIN]
    tail = P[0] + LINK * unit(P[0] - P[1])
    head = P[-1] + LINK * unit(P[-1] - P[-2])
    poly = [tail] + P + [head]

    def to3(q):
        return c + np.array([q[0], q[1], 0.0])

    monos, centers, free = VGroup(), [], None
    for k, (u, v) in enumerate(r2_points(2000)):
        q = np.array([(2 * u - 1) * (half - 0.25), (2 * v - 1) * (half - 0.25)])
        if min(seg_dist(q, a, b) for a, b in zip(poly[:-1], poly[1:])) < 0.34:
            continue
        if any(np.linalg.norm(q - p) < 0.38 for p in centers):
            continue
        if len(monos) == n_mono:
            free = q
            break
        ang = PI * ((k * 0.6180339887) % 1.0)
        d = 0.5 * LINK * np.array([np.cos(ang), np.sin(ang), 0.0])
        monos.add(Line(to3(q) - d, to3(q) + d).set_stroke(MONO_COLOR, 4))
        centers.append(q)
    assert len(monos) == n_mono and free is not None
    box = Square(side_length=2 * half).move_to(c).set_stroke(GREY_B, 3)
    chain = VMobject().set_points_smoothly([to3(p) for p in P]).set_stroke(STEP_COLOR, 5)
    link = Line(to3(tail), to3(P[0])).set_stroke(STEP_COLOR, 5)
    ends = VGroup(Dot(to3(tail), radius=0.08).set_fill(L.START_COLOR, 1),
                  Dot(to3(P[-1]), radius=0.08).set_fill(L.END_COLOR, 1))
    near = int(np.argmin([np.linalg.norm(q - head) for q in centers]))
    return dict(box=box, chain=chain, link=link, ends=ends, monos=monos, near=near,
                p0=to3(P[0]), p_end=to3(P[-1]), head=to3(head), free=to3(free))


def chevron(p, d, size, color):
    """朝 d 方向的小箭頭記號「>」，中心在 p。"""
    d = unit(d)
    perp = np.array([-d[1], d[0], 0.0])
    p = np.array(p, dtype=float)
    pts = [p - 0.5 * size * d + size * perp, p + 0.5 * size * d, p - 0.5 * size * d - size * perp]
    return VMobject().set_points_as_corners(pts).set_stroke(color, 3)


def torus_marks(cell, origin, size=0.12, n=SIDE):
    """對邊黏在一起的記號：左右邊各一個朝上的單箭頭、上下邊各一個朝右的雙箭頭。"""
    lo, hi, mid = -0.5, n - 0.5, (n - 1) / 2
    g = VGroup()
    for x in (lo, hi):
        g.add(chevron(L.lat_point(x, mid, cell, origin), UP, size, MARK1_COLOR))
    for y in (lo, hi):
        p = L.lat_point(mid, y, cell, origin)
        for k in (-1, 1):
            g.add(chevron(p + k * 0.9 * size * RIGHT, RIGHT, size, MARK2_COLOR))
    return g


def torus_box(cell, origin, dot_r=0.06, mark=0.12, n=SIDE):
    """n×n 環面：格點畫在格子中心，盒子是 [−½, n−½]² → VGroup(盒子, 格點, 記號)。"""
    mid = (n - 1) / 2
    box = Square(side_length=n * cell).move_to(L.lat_point(mid, mid, cell, origin)).set_stroke(GREY_B, 2.5)
    dots = VGroup(*[Dot(L.lat_point(x, y, cell, origin), radius=dot_r).set_fill(GREY_A, 1)
                    for y in range(n) for x in range(n)])
    return VGroup(box, dots, torus_marks(cell, origin, mark, n))


def walk_lines(segments, start, cell, origin, width=6, end_r=0.09, color=STEP_COLOR):
    """ex["segments"] 每段一條 Line（越界的半段換色）→ VGroup(線段們, 起點, 終點)。"""
    lines = VGroup(*[Line(L.lat_point(a[0], a[1], cell, origin), L.lat_point(b[0], b[1], cell, origin))
                     .set_stroke(WRAP_COLOR if w else color, width) for a, b, w in segments])
    s0 = Dot(L.lat_point(start[0], start[1], cell, origin), radius=end_r).set_fill(L.START_COLOR, 1)
    e = segments[-1][1]
    s1 = Dot(L.lat_point(e[0], e[1], cell, origin), radius=end_r).set_fill(L.END_COLOR, 1)
    return VGroup(lines, s0, s1)


def lift_figure(ex, cell, origin, n=SIDE):
    """把環面上的漫步攤開到無限平面：背景兩個盒子的複本（上下）→ VGroup(複本, 足跡, 終點座標)。"""
    mid = (n - 1) / 2
    copies = VGroup()
    for oy in (0, -n):
        copies.add(Square(side_length=n * cell).move_to(L.lat_point(mid, mid + oy, cell, origin))
                   .set_stroke(GREY_D, 2))
        copies.add(*[Dot(L.lat_point(x, y + oy, cell, origin), radius=0.035).set_fill(GREY_C, 0.8)
                     for y in range(n) for x in range(n)])
    path = L.trail(ex["lift"], cell, origin, color=STEP_COLOR, stroke_width=4)
    xe, ye = ex["lift"][-1]
    lab = tx(f"({xe},{ye})", 24, GREY_A).next_to(L.lat_point(xe, ye, cell, origin), RIGHT, buff=0.14)
    return VGroup(copies, path, lab)


def count_table(lo, c4):
    """N = lo..lo+7 的 c_N：3 列 × 9 欄（左上角在 ORIGIN）。"""
    cols = range(lo, lo + 8)
    rows = [["N"] + [str(N) for N in cols],
            [mixed("環面 ", "$c_N$", font_size=24)] + [str(c4[N]) for N in cols],
            [mixed("無限平面 ", "$c_N$", font_size=24)] + [str(WK.SAW_INF[N]) for N in cols]]
    return D.tex_table(rows, [2.3] + [1.3] * 8, row_h=0.52, font_size=26)


def tcell(t, i, j):
    return D.table_cell(t, i, j, 9)


# ---------------------------------------------------------------- 各段
def part_polymer(s):
    pb = polymer_box(np.array([-3.75, 0.1, 0]))
    box, chain, link, ends, monos = pb["box"], pb["chain"], pb["link"], pb["ends"], pb["monos"]
    lx = tx(r"\vec x", 30, L.START_COLOR).next_to(ends[0], LEFT, buff=0.08)
    ly = tx(r"\vec y", 30, L.END_COLOR).next_to(ends[1], UP, buff=0.08)
    fig = VGroup(box, chain, link, ends, monos, lx, ly)
    eq = place_in(tx(r"G(z;\vec x,\vec y)=\sum_N C(N;\vec x,\vec y)\,z^N", 40), 0.2, 6.7, 1.75, 2.95)
    zl = mixed("$z$", "：單體的逸度（fugacity）", font_size=28)
    m1 = zh("單體可以接上鏈的兩端，", font_size=28)
    m2 = zh("也可以蒸發回池裡", font_size=28)
    col = VGroup(zl, m1, m2).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    place_in(col, 0.4, 6.7, -0.45, 1.45)
    err = D.erratum_box("p. 235，(10.1)", r"\sum_N C(N;\vec x,\vec y)\,z^n", r"\sum_N C(N;\vec x,\vec y)\,z^N",
                        width=6.3)
    place_in(err, 0.3, 6.7, -2.9, -0.75)
    assert_in_frame(fig, eq, col, err)
    for a, b in ((fig, eq), (fig, col), (fig, err), (eq, col), (col, err)):
        assert_apart(a, b)
    total = sum(WK.torus_counts(SIDE))
    mono = monos[pb["near"]]

    s.say("一條高分子關在小盒子裡，和一池單體平衡", wait=1.5)
    s.play(ShowCreation(box), run_time=0.8)
    s.play(ShowCreation(chain), ShowCreation(link), FadeIn(ends), FadeIn(lx), FadeIn(ly), run_time=1.6)
    s.play(FadeIn(monos, lag_ratio=0.05), run_time=1.2)
    s.say("盒子的牆會怎樣改變它能走的路？", wait=2.5)
    s.play(Write(eq), run_time=1.5)
    s.say("(10.1)：從 ", r"$\vec x$", " 走到 ", r"$\vec y$", "、長 ", "$N$", " 的鏈有 ",
          r"$C(N;\vec x,\vec y)$", " 種", wait=2.5)
    s.play(FadeIn(zl), run_time=0.7)
    s.say("$z$", " 是單體的逸度：每接上一個單體乘一個 ", "$z$", wait=2.0)
    s.play(FadeIn(m1), run_time=0.6)
    s.play(mono.animate.put_start_and_end_on(pb["p_end"], pb["head"]).set_stroke(STEP_COLOR),
           ends[1].animate.move_to(pb["head"]), ly.animate.next_to(pb["head"], UP, buff=0.08), run_time=1.4)
    s.play(FadeIn(m2), run_time=0.6)
    s.play(link.animate.move_to(pb["free"]).set_stroke(MONO_COLOR), ends[0].animate.move_to(pb["p0"]),
           lx.animate.next_to(pb["p0"], LEFT, buff=0.08), run_time=1.4)
    s.say("單體可以接上鏈的兩端，也可以蒸發回池裡", wait=2.0)
    s.play(FadeIn(err), run_time=0.9)
    s.say("書上把 ", "$z$", " 的次方印成 ", "$n$", wait=2.0)
    s.say("$n$", " 是自旋的分量數；照書上 ", r"$n\to0$", " 時和 ", "$z$", " 無關", wait=2.5)
    s.say(r"$4\times4$", " 環面永遠是 ", f"${total}$", " 條，", r"$\langle N\rangle=0$", wait=2.5)
    s.clear_stage()


def part_torus(s):
    cell = 1.05
    origin = np.array([-3.9 - 1.5 * cell, -0.1 - 1.5 * cell, 0])
    tor = torus_box(cell, origin)
    legend = VGroup(
        VGroup(Line(ORIGIN, 0.5 * RIGHT).set_stroke(STEP_COLOR, 6), zh("一般的一步", font_size=24))
        .arrange(RIGHT, buff=0.2),
        VGroup(Line(ORIGIN, 0.5 * RIGHT).set_stroke(WRAP_COLOR, 6), zh("越界的半段（從對邊回來）", font_size=24))
        .arrange(RIGHT, buff=0.2),
    ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
    place_left(legend, -1.1, -2.4, 3.9)
    notes = {"short": "不越界：像在無限平面上", "mid": None, "full": None}
    walks, descs = {}, {}
    for name in ("short", "mid", "full"):
        ex = WK.torus_example(name)
        if name == "short":
            note = notes[name]
        elif name == "mid":
            note = f"越界 {ex['wraps']} 次"
        else:
            note = f"越界 {ex['wraps']} 次，{len(set(ex['sites']))} 格全部走過"
        l1 = mixed(f"$N={ex['N']}$", rf"$\mathtt{{{ex['steps']}}}$", font_size=30, buff=0.4)
        l2 = zh(note, font_size=28)
        dsc = VGroup(l1, l2).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
        descs[name] = place_left(dsc, -1.1, 2.3, 6.7)
        walks[name] = walk_lines(ex["segments"], ex["start"], cell, origin)
    mid_ex = WK.torus_example("mid")
    lcell = 0.36
    lorigin = np.array([4.9 - 1.5 * lcell, -0.5 + 0.5 * lcell, 0])
    lift = lift_figure(mid_ex, lcell, lorigin)
    lift_lab = zh("攤開到無限平面", font_size=24).move_to([4.9, 1.3, 0])
    assert_in_frame(tor, legend, lift, lift_lab, *descs.values(), *walks.values())
    for d in descs.values():
        assert_apart(tor, d)
        assert_apart(d, lift_lab)
    for m in (legend, lift, lift_lab):
        assert_apart(tor, m)
    assert_apart(legend, lift)
    assert_apart(lift, lift_lab)
    xe, ye = mid_ex["lift"][-1]

    s.play(ShowCreation(tor[0]), FadeIn(tor[1]), run_time=1.2)
    s.say("週期邊界：", r"$4\times4$", " 環面，16 個格點畫在格子中心", wait=2.0)
    s.play(FadeIn(tor[2], lag_ratio=0.2), run_time=1.0)
    s.say("同樣記號的兩邊黏在一起：走出右邊就從左邊回來", wait=2.5)
    w = walks["short"]
    s.play(FadeIn(descs["short"]), FadeIn(legend), run_time=0.6)
    s.play(ShowCreation(w[0], lag_ratio=1.0), FadeIn(w[1]), run_time=1.5)
    s.play(FadeIn(w[2]), run_time=0.3)
    s.say("$N=3$", "：不越界，和無限平面上一樣", wait=2.0)
    s.play(FadeOut(w), FadeOut(descs["short"]), run_time=0.5)
    w = walks["mid"]
    s.play(FadeIn(descs["mid"]), run_time=0.5)
    s.play(ShowCreation(w[0], lag_ratio=1.0), FadeIn(w[1]), run_time=2.5)
    s.play(FadeIn(w[2]), run_time=0.3)
    s.say("$N=8$", f"：越界 {mid_ex['wraps']} 次，越界的半段畫成橘色", wait=2.5)
    s.play(FadeIn(lift[0]), FadeIn(lift_lab), run_time=0.8)
    s.play(ShowCreation(lift[1][0]), FadeIn(lift[1][1]), run_time=2.0)
    s.play(FadeIn(lift[1][2]), FadeIn(lift[2]), run_time=0.5)
    s.say("攤開到無限平面：終點 ", f"$({xe},{ye})$", "，背景是盒子的複本", wait=2.5)
    s.play(FadeOut(w), FadeOut(descs["mid"]), FadeOut(lift), FadeOut(lift_lab), run_time=0.6)
    w = walks["full"]
    s.play(FadeIn(descs["full"]), run_time=0.5)
    s.play(ShowCreation(w[0], lag_ratio=1.0), FadeIn(w[1]), run_time=3.5)
    s.play(FadeIn(w[2]), run_time=0.3)
    s.say("$N=15$", "：越界 3 次，16 格全部走過", wait=2.0)
    s.say("小盒子裡的長漫步：一盤義大利麵", wait=2.0)
    s.clear_stage()


def part_table(s):
    c4 = tuple(WK.torus_counts(SIDE))
    c3, c5 = WK.torus_counts(3), WK.torus_counts(5)
    t1 = count_table(0, c4).shift(np.array([-6.35, 2.95, 0]))
    t2 = count_table(8, c4).shift(np.array([-6.35, 1.25, 0]))
    same = SurroundingRectangle(VGroup(*[tcell(t1, i, j) for i in (1, 2) for j in (1, 2, 3, 4)]), buff=0.07)
    same.set_stroke(SAME_COLOR, 2.5)
    diff = SurroundingRectangle(VGroup(tcell(t1, 1, 5), tcell(t1, 2, 5)), buff=0.07).set_stroke(WRAP_COLOR, 2.5)
    full = SurroundingRectangle(tcell(t2, 1, 8), buff=0.07).set_stroke(STEP_COLOR, 2.5)
    nd = WK.torus_first_diff(SIDE)[1]
    lb = mixed("$N=4$", "：", f"${c4[4]}$", " vs ", f"${CINF[4]}$", f"，少了 {CINF[4] - c4[4]} 條直線（",
               r"$\mathtt{EEEE}$", " 等繞一圈回到起點）", font_size=26)
    lc = mixed(f"$N={SIDE * SIDE - 1}$", "（走滿）", f"${c4[-1]}$", " 條；全部 ", f"${sum(c4)}$", " 條", font_size=26)
    ld = mixed("$N=L$", f" 時剛好少 {nd} 條：", r"$3\times3$", " 在 ", "$N=3$", f"（{c3[3]} vs {WK.SAW_INF[3]}）、",
               r"$5\times5$", " 在 ", "$N=5$", f"（{c5[5]} vs {WK.SAW_INF[5]}）", font_size=26)
    place_left(lb, -6.3, -0.95, 4.6)
    place_left(lc, -6.3, -1.65, 4.6)
    place_left(ld, -6.3, -2.4, 6.7)
    lcell = 0.28
    lorigin = np.array([5.7 - 1.5 * lcell, -1.3 - 1.5 * lcell, 0])
    loop = VGroup(torus_box(lcell, lorigin, dot_r=0.03, mark=0.06),
                  walk_lines(LOOP_SEGS, (0, 1), lcell, lorigin, width=4, end_r=0.06, color=LOOP_COLOR))
    assert_in_frame(t1, t2, lb, lc, ld, loop, same, diff, full)
    assert_table(t1, 9)
    assert_table(t2, 9)
    assert_apart(same, diff)
    for a, b in ((t1, t2), (t2, lb), (lb, lc), (lc, ld), (lb, loop), (lc, loop), (ld, loop), (t2, loop)):
        assert_apart(a, b)

    s.play(FadeIn(t1), run_time=1.0)
    s.say("從起點出發的自迴避漫步條數 ", "$c_N$", "：環面 vs 無限平面", wait=2.5)
    s.play(ShowCreation(same), run_time=0.8)
    s.say(r"$N\le3$", " 完全一樣：還碰不到盒子的邊", wait=2.0)
    s.play(ShowCreation(diff), FadeIn(lb), run_time=0.9)
    s.say("$N=4$", "：", f"${c4[4]}$", " vs ", f"${CINF[4]}$", f"，少了 {CINF[4] - c4[4]} 條直線", wait=2.0)
    s.play(FadeIn(loop[0]), run_time=0.5)
    s.play(ShowCreation(loop[1][0], lag_ratio=1.0), FadeIn(loop[1][1]), run_time=1.2)
    s.play(FadeIn(loop[1][2]), run_time=0.3)
    s.say(r"$\mathtt{EEEE}$", " 繞一圈回到起點，不再自迴避", wait=2.0)
    s.play(FadeIn(t2), run_time=1.0)
    s.say("無限平面的 ", "$c_N$", " 一直變多；環面先多後少", wait=2.5)
    s.play(ShowCreation(full), FadeIn(lc), run_time=0.9)
    s.say(f"$N={SIDE * SIDE - 1}$", " 走滿 16 格：", f"${c4[-1]}$", " 條；全部 ", f"${sum(c4)}$", " 條", wait=2.5)
    s.play(FadeIn(ld), run_time=0.8)
    s.say("其他大小也一樣：", "$N=L$", f" 時剛好少 {nd} 條", wait=2.5)
    s.say(r"$4\times4$", " 環面就是四維超立方體", wait=2.0)
    s.say(rf"$16\times{c4[-1]}={WK.GRAY_Q4}$", "：正是四維 Gray 碼的條數", wait=2.5)
    s.clear_stage()


def part_regimes(s):
    cell, yc = 0.5, 0.55
    xs = (-4.4, 0.0, 4.4)
    names = ("short", "mid", "full")
    ztex = (f"z={Z_CASES[0]:.1f}", f"z=z_c={f4(WC.ZC_2D)}", f"z={Z_CASES[2]:.1f}")
    notes = ("漫步比盒子小", "和盒子差不多", f"擠滿（最多 {SIDE * SIDE - 1} 步）")
    groups = []
    for k in range(3):
        origin = np.array([xs[k] - 1.5 * cell, yc - 1.5 * cell, 0])
        tor = torus_box(cell, origin, dot_r=0.045, mark=0.08)
        ex = WK.torus_example(names[k])
        w = walk_lines(ex["segments"], ex["start"], cell, origin, width=4, end_r=0.07)
        lab = tx(ztex[k], 32).move_to([xs[k], 2.15, 0])
        mean = tx(rf"\langle N\rangle={f4(WK.torus_mean_N(SIDE, Z_CASES[k]))}", 32).move_to([xs[k], -0.95, 0])
        note = zh(notes[k], font_size=26).move_to([xs[k], -1.6, 0])
        groups.append((tor, w, lab, mean, note))
    blocks = [VGroup(*g) for g in groups]
    assert_in_frame(*blocks)
    for a, b in zip(blocks, blocks[1:]):
        assert_apart(a, b)
    for tor, w, lab, mean, note in groups:
        assert_apart(lab, tor)
        assert_apart(tor, mean)
        assert_apart(mean, note)

    s.play(*[FadeIn(g[0]) for g in groups], *[FadeIn(g[2]) for g in groups], run_time=1.0)
    s.say("三種 ", "$z$", "：平均步數 ", r"$\langle N\rangle$", " 和盒子比一比", wait=2.0)
    for k, (tor, w, lab, mean, note) in enumerate(groups):
        s.play(ShowCreation(w[0], lag_ratio=1.0), FadeIn(w[1]), run_time=1.2)
        s.play(FadeIn(w[2]), FadeIn(mean), FadeIn(note), run_time=0.6)
        s.say(f"${ztex[k]}$", "：", rf"$\langle N\rangle={MEAN_N[k]}$", "，" + notes[k], wait=2.2)
    s.say("先取 ", r"$n\to0$", "，再讓盒子變大", wait=2.0)
    s.say("順序對調會出錯（W1002 會看到為什麼）", wait=2.5)
    s.say("盒子小、漫步長：漫步就被擠成一盤義大利麵", wait=2.5)
    s.say("本章用 ", r"$\mathrm{O}(n)$", " 模型的平均場與自旋波", wait=2.0)
    s.say("來描述這個「凝聚相」", wait=2.5)


class W1001ConfinedPolymer(ChuScene):
    def construct(self):
        self.show_header("W1001　關在盒子裡的高分子")
        check_math()
        part_polymer(self)
        part_torus(self)
        part_table(self)
        part_regimes(self)
        self.wait(1.0)
