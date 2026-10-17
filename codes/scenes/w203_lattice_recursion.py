"""
W203 遞迴：最後一步來自四個鄰居（§2.1.2，pp. 29–31，(2.4)–(2.9)，Fig. 2.2–2.3）

問：正方格子上走 4 步回到原點，有幾種走法？（先猜；答案 36）
  1. 回顧一維（W110、W111）：(2.4) C(N;x,y)、(2.5) G(z;x,y) 都是積分，(2.6) χ(q) = 2cos q（書上寫成 χ(k)）；一維 C(4;0,0) = 6。
  2. 遞迴 (2.7)（Fig. 2.3）：最後一步一定來自四個鄰居 → C(N;x,y) = Σ_i C(N−1;x,w_i)，C(0;x,y) = δ。
  3. 數字格 N = 0..4 一張張長出來：每格＝上一張四個鄰居的和；原點 1, 0, 4, 0, 36，總和 4^N；
     棋盤格（奇偶性）；彩蛋：兩枚硬幣 u = x+y、v = x−y：36 = 6×6、16 = 4×4。
  4. 乘上 z^N 加總 (2.8)，補上 δ → (2.9) G = zΣ_i G(z;x,w_i) + δ：非局部 → 下一幕用傅立葉。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import boxed, mixed, zh
from utils import lattice_draw as L
from utils import walk_draw as D
from utils import walk_gf as GF

EXAMPLES = ["ENWS", "EWNS", "SSNN"]                      # 開頭的三個 4 步回家的例子
GRIDS = [GF.count_grid(N) for N in range(5)]
CELL, ORG = 0.62, np.array([-3.3, 0.1, 0])


def walk(word):
    pts = [(0, 0)]
    for ch in word:
        dx, dy = GF.SQ_STEPS[GF.SQ_NAMES.index(ch)]
        pts.append((pts[-1][0] + dx, pts[-1][1] + dy))
    return pts


def nbr_sum(N):
    """上一張（N−1）四個鄰居的值與它們的和（＝第 N 張的中心）。"""
    vals = [GRIDS[N - 1].get(w, 0) for w in GF.SQ_STEPS]
    return vals, sum(vals)


def check_math():
    assert all(walk(w)[-1] == (0, 0) and len(w) == 4 for w in EXAMPLES)
    assert GF.count1(4, 0) == 6
    M = 16                                                           # (2.4)：(1/2π)∫(2cos q)^4 dq（均分求積，精確）
    assert abs(sum(GF.chi_1d(-math.pi + 2 * math.pi * k / M) ** 4 for k in range(M)) / M - 6) < 1e-9
    assert GF.chi_1d(0) == 2
    for N in range(1, 5):                                            # (2.7)：每格＝上一張四個鄰居的和
        g, h = GRIDS[N], GRIDS[N - 1]
        cells = {(x, y) for x in range(-5, 6) for y in range(-5, 6)}
        assert all(g.get(c, 0) == sum(h.get((c[0] - dx, c[1] - dy), 0) for dx, dy in GF.SQ_STEPS) for c in cells)
    assert [g.get((0, 0), 0) for g in GRIDS] == [1, 0, 4, 0, 36]
    assert [sum(g.values()) for g in GRIDS] == [1, 4, 16, 64, 256] == [4 ** N for N in range(5)]
    g2 = GRIDS[2]
    assert g2[(0, 0)] == 4 and all(g2[(a, b)] == 2 for a in (1, -1) for b in (1, -1))
    assert all(g2[c] == 1 for c in [(2, 0), (-2, 0), (0, 2), (0, -2)])
    g4 = GRIDS[4]
    want = {(0, 0): 36, (1, 1): 24, (2, 0): 16, (2, 2): 6, (3, 1): 4, (4, 0): 1}
    for (a, b), v in want.items():                                    # 用對稱性展開成全部格點
        for c in {(sa * a, sb * b) for sa in (1, -1) for sb in (1, -1)} | {(sa * b, sb * a) for sa in (1, -1) for sb in (1, -1)}:
            assert g4[c] == v, (c, g4[c], v)
    assert len(g4) == 25
    for N, g in enumerate(GRIDS):                                     # 棋盤格：奇偶不同的格子都是 0
        assert all((x + y - N) % 2 == 0 for (x, y), v in g.items() if v)
    assert nbr_sum(2) == ([1, 1, 1, 1], 4) and nbr_sum(4) == ([9, 9, 9, 9], 36)
    assert GF.count_sq(4, 0, 0) == 36 == GF.count1(4, 0) ** 2 == 6 * 6
    assert GF.count_sq(4, 2, 0) == 16 == GF.count1(4, 2) ** 2 == 4 * 4
    assert all(GF.count_sq(4, x, y) == g4.get((x, y), 0) for x in range(-4, 5) for y in range(-4, 5))
    z, big = 0.05, GF.count_grids(20)                                 # (2.9) 的數值核對（截到 N = 20）
    G = lambda r: sum(z ** N * big[N].get(r, 0) for N in range(21))
    assert abs(G((0, 0)) - (z * sum(G(w) for w in GF.SQ_STEPS) + 1)) < 1e-10


# ---------------------------------------------------------------- 小工具
def eq(tex, tag=None, font_size=32, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return m
    return VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


# ---------------------------------------------------------------- 各段
def part_question(s):
    pics = VGroup()
    for k, w in enumerate(EXAMPLES):
        org = np.array([-3.6 + 3.6 * k, 0.75, 0])
        pics.add(VGroup(L.lattice_lines(2, 0.45, org), L.trail(walk(w), 0.45, org, stroke_width=5, end_radius=0.08),
                        Tex(rf"\text{{{w}}}", font_size=30).move_to(org + 1.55 * DOWN)))
    note = zh("（E、N、W、S：東北西南各一步）", font_size=24, color=GREY_B).move_to(np.array([0, -1.75, 0]))
    assert_in_frame(pics, note)
    assert_apart(pics, note)
    s.say("問：正方格子上走 4 步回到原點，有幾種走法？", wait=0.3)
    s.play(LaggedStart(*[FadeIn(p) for p in pics], lag_ratio=0.3), run_time=1.8)
    s.play(FadeIn(note), run_time=0.6)
    s.wait(3.0)
    s.say("這裡有三種……全部有幾種？先猜一個數字！", wait=3.5)
    s.clear_stage()


def part_recap(s):
    r1 = eq(r"C(N;x,y)=\frac{1}{2\pi}\int_{-\pi}^{\pi}e^{-iq(x-y)}\chi(q)^N\,dq", "(2.4)", font_size=36)
    r2 = eq(r"G(z;x,y)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{e^{-iq(x-y)}}{1-\chi(q)z}\,dq", "(2.5)", font_size=36)
    r3 = eq(r"\chi(q)=2\cos q", "(2.6)", font_size=36)
    r4 = mixed("例：", "$N=4,\\ x=y$", "：", rf"$C(4;0,0)={GF.count1(4, 0)}$", "（記住這個 6）", font_size=26)
    col = VGroup(r1, r2, r3, r4).arrange(DOWN, buff=0.45).move_to(np.array([0, 0.3, 0]))
    assert_in_frame(col)
    s.say("回顧一維（W110、W111）：走法數與生成函數都是積分", wait=0.3)
    s.play(Write(r1[0]), FadeIn(r1[1]), run_time=1.6)
    s.play(Write(r2[0]), FadeIn(r2[1]), run_time=1.6)
    s.wait(2.0)
    s.play(FadeIn(r3), run_time=0.8)
    s.say("$\\chi(q)=2\\cos q$", "：書上 (2.6) 寫成 ", "$\\chi(k)$", "，同一個函數", wait=3.0)
    s.play(FadeIn(r4), run_time=0.8)
    s.say("一維走 4 步回到原點：6 種", wait=2.5)
    s.say("怎麼推廣到二維、三維？先從遞迴下手", wait=3.0)
    s.clear_stage()


def part_recursion(s):
    star = L.neighbor_star(cell=1.25, origin=np.array([-4.5, 0.2, 0]), font_size=30)
    t1 = zh("最後一步一定來自四個鄰居之一", font_size=28, color=YELLOW)
    t2 = eq(r"C(N;\vec x,\vec y)=\sum_{i=1}^{4}C(N-1;\vec x,\vec w_i)", "(2.7)", font_size=38)
    t3 = Tex(r"C(0;\vec x,\vec y)=\delta_{\vec x,\vec y}", font_size=36)
    t4 = mixed("$\\vec x$", "：起點，", "$\\vec y$", "：終點，", "$\\vec w_i$", "：", "$\\vec y$", " 的四個鄰居",
               font_size=24, color=GREY_B)
    col = VGroup(t1, t2, t3, t4).arrange(DOWN, aligned_edge=LEFT, buff=0.5)
    fit(col, 8.2)
    place_col(col, -1.9, 1.75)
    assert_in_frame(star, col)
    assert_apart(star, col)
    s.say("怎麼走到 ", "$\\vec y$", "？看最後一步（Fig. 2.3）", wait=0.3)
    s.play(FadeIn(star[0], scale=0.5), FadeIn(star[1]), FadeIn(star[3]), run_time=1.0)
    s.play(LaggedStart(*[GrowArrow(a) for a in star[2]], lag_ratio=0.2), FadeIn(t1), run_time=1.4)
    s.wait(2.0)
    s.play(Write(t2[0]), FadeIn(t2[1]), run_time=1.6)
    s.say("走到 ", "$\\vec y$", " 的走法＝走到四個鄰居的走法加起來", wait=3.0)
    s.play(FadeIn(t3), FadeIn(t4), run_time=1.0)
    s.say("起點：0 步時只有 ", "$\\vec x$", " 自己是 1", wait=3.0)
    s.clear_stage()


def part_grids(s):
    lines = L.lattice_lines(4, CELL, ORG)
    mobs = [L.count_grid(g, CELL, ORG, font_size=30, zero_dots=True, n=4) for g in GRIDS]
    nlab = [place_col(Tex(f"N={N}", font_size=44).set_color(YELLOW), 0.6, 2.95) for N in range(5)]
    orig = [g.get((0, 0), 0) for g in GRIDS]
    tot = [sum(g.values()) for g in GRIDS]
    rows = [["N", r"C(N;\vec 0)", r"\sum_{\vec r}C(N;\vec r)"]] + [[str(N), str(orig[N]), str(tot[N])] for N in range(5)]
    tab = D.tex_table(rows, [0.9, 2.0, 2.5], row_h=0.5, font_size=28)
    tab.shift(np.array([0.6, nlab[0].get_bottom()[1] - 0.2, 0]))
    hdr = VGroup(*[D.table_cell(tab, 0, j, 3) for j in range(3)])
    rowm = [VGroup(*[D.table_cell(tab, i + 1, j, 3) for j in range(3)]) for i in range(5)]
    D.table_cell(tab, 5, 1, 3).set_color(YELLOW)
    sums = {}
    for N in (2, 4):
        vals, tot_n = nbr_sum(N)
        sums[N] = place_col(Tex("+".join(map(str, vals)) + f"={tot_n}", font_size=36).set_color(YELLOW),
                            0.6, tab.get_bottom()[1] - 0.35)
    ring = Circle(radius=0.27).set_stroke(YELLOW, 3).move_to(L.lat_point(0, 0, CELL, ORG))
    board = VGroup(*[Square(side_length=CELL * 0.96).set_stroke(width=0).set_fill(BLUE_E, 0.5)
                     .move_to(L.lat_point(x, y, CELL, ORG)) for x in range(-4, 5) for y in range(-4, 5) if (x + y) % 2 == 0])
    board.set_z_index(-1)
    c40, c42 = GF.count1(4, 0), GF.count1(4, 2)
    coins = VGroup(mixed("兩枚硬幣：", "$u=x+y,\\ v=x-y$", " 各擲一次", font_size=24),
                   Tex(rf"C(4;\vec 0)={GF.count_sq(4, 0, 0)}={c40}\times{c40}", font_size=32),
                   Tex(rf"C(4;(2,0))={GF.count_sq(4, 2, 0)}={c42}\times{c42}", font_size=32))
    place_col(coins.arrange(DOWN, aligned_edge=LEFT, buff=0.22), 0.6, tab.get_bottom()[1] - 0.3)
    assert_in_frame(lines, tab, sums[2], sums[4], coins, *nlab)
    assert_apart(lines, tab)
    assert_apart(lines, coins)
    assert_apart(tab, coins)

    def nbr_box(g):
        return VGroup(*[SurroundingRectangle(L.grid_number(g, x, y), buff=0.06).set_stroke(YELLOW, 2.5)
                        for x, y in GF.SQ_STEPS])

    s.say("把每一格寫上「走 ", "$N$", " 步到這裡的走法數」", wait=0.3)
    s.play(FadeIn(lines), FadeIn(mobs[0]), FadeIn(nlab[0]), FadeIn(hdr), FadeIn(tab[1]), FadeIn(rowm[0]), run_time=1.2)
    s.wait(1.5)
    s.play(FadeOut(mobs[0]), FadeIn(mobs[1]), FadeTransform(nlab[0], nlab[1]), FadeIn(rowm[1]), run_time=1.0)
    s.say("$N=1$", "：四個鄰居各 1", wait=2.5)
    boxes = nbr_box(mobs[1])
    s.play(ShowCreation(boxes), run_time=0.8)
    s.play(FadeIn(sums[2]), run_time=0.6)
    s.say("每一格＝上一張四個鄰居的和：中心 ", "$1+1+1+1=4$", wait=2.5)
    s.play(FadeOut(boxes), FadeOut(mobs[1]), FadeIn(mobs[2]), FadeTransform(nlab[1], nlab[2]), FadeIn(rowm[2]),
           run_time=1.0)
    s.say("$N=2$", "：中心 4，斜角 2，再往外 1", wait=3.0)
    s.play(FadeOut(sums[2]), FadeOut(mobs[2]), FadeIn(mobs[3]), FadeTransform(nlab[2], nlab[3]), FadeIn(rowm[3]),
           run_time=1.0)
    s.say("$N=3$", "：原點是 0，奇數步回不了家", wait=2.5)
    boxes = nbr_box(mobs[3])
    s.play(ShowCreation(boxes), FadeIn(sums[4]), run_time=0.8)
    s.say("揭曉：中心＝", "$9+9+9+9$", wait=2.5)
    s.play(FadeOut(boxes), FadeOut(mobs[3]), FadeIn(mobs[4]), FadeTransform(nlab[3], nlab[4]), FadeIn(rowm[4]),
           ShowCreation(ring), run_time=1.2)
    s.say("走 4 步回到原點：36 種走法！", wait=3.0)
    s.say("總和 1、4、16、64、256：正好是 ", "$4^N$", wait=3.0)
    s.say("往外依序 24、16、6、4、1，全部加起來 256", wait=3.0)
    s.play(FadeIn(board), run_time=1.0)
    s.say("棋盤格：", "$x+y$", " 與 ", "$N$", " 奇偶不同的格子都是 0", wait=3.5)
    s.play(FadeOut(sums[4]), FadeIn(coins), run_time=1.0)
    s.say("彩蛋：", "$u=x+y$", "、", "$v=x-y$", " 各是一枚硬幣", wait=3.0)
    s.say("$36=6\\times6$", "、", "$16=4\\times4$", "：一維的 6 又出現了", wait=3.5)
    s.clear_stage()


def part_generating(s):
    g0 = mixed("把 (2.7) 乘上 ", "$z^N$", "，對 ", "$N\\ge1$", " 加總：", font_size=26)
    g1 = eq(r"\sum_{N\ge1}z^NC(N;\vec x,\vec y)=z\sum_{i}\sum_{N\ge1}z^{N-1}C(N-1;\vec x,\vec w_i)", "(2.8)",
            font_size=32)
    g2 = mixed("左邊少了 ", "$N=0$", " 的項 ", "$C(0;\\vec x,\\vec y)=\\delta_{\\vec x,\\vec y}$", "：補上", font_size=24)
    g3 = boxed(eq(r"G(z;\vec x,\vec y)=z\sum_{i=1}^{4}G(z;\vec x,\vec w_i)+\delta_{\vec x,\vec y}", "(2.9)",
                  font_size=40), color=D.GAUSS_COLOR, buff=0.22)
    g4 = mixed("只剩 ", "$z$", " 一個變數，但右邊要用到鄰居（非局部）", font_size=24)
    g5 = zh("下一幕：用傅立葉一次解出全部", font_size=24, color=GREY_B)
    col = VGroup(g0, g1, g2, g3, g4, g5).arrange(DOWN, buff=0.38)
    for m in col:
        fit(m, 12.6)
    col.arrange(DOWN, buff=0.38).move_to(np.array([0, 0.15, 0]))
    assert_in_frame(col)
    s.say("把整串走法數打包成生成函數", wait=0.3)
    s.play(FadeIn(g0), run_time=0.8)
    s.play(Write(g1[0]), FadeIn(g1[1]), run_time=2.0)
    s.say("右邊的 ", "$N-1$", " 換個名字，就是鄰居的生成函數", wait=3.0)
    s.play(FadeIn(g2), run_time=0.8)
    s.say("左邊少了 ", "$N=0$", " 那一項：0 步時只有起點算 1", wait=3.0)
    s.play(FadeIn(g3, shift=0.2 * UP), run_time=1.0)
    s.say("得到 (2.9)：生成函數也滿足「四個鄰居」的關係", wait=3.5)
    s.play(FadeIn(g4), run_time=0.8)
    s.say("只剩 ", "$z$", " 一個變數，但每一格都牽扯到鄰居", wait=3.0)
    s.play(FadeIn(g5), run_time=0.6)
    s.say("下一幕：用傅立葉一次解出全部", wait=2.5)
    s.clear_stage()


def part_close(s):
    s.say("每一格＝四個鄰居的和：這是二維的 Pascal 三角形", wait=4.0)


class W203LatticeRecursion(ChuScene):
    def construct(self):
        self.show_header("W203　遞迴：最後一步來自四個鄰居")
        check_math()
        part_question(self)
        part_recap(self)
        part_recursion(self)
        part_grids(self)
        part_generating(self)
        part_close(self)
        self.wait(1.0)
