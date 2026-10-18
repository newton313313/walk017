"""
W902 O(n) 模型與對稱破缺（§9.2.1–9.2.2，pp. 205–211，(9.4)–(9.15)，Fig. 9.2）

問：能量在「所有自旋一起轉」之下不變，所以平均磁化一定是 0——那永久磁鐵是怎麼來的？
  1. O(n)：n = 1 Ising、2 xy、3 Heisenberg、0 自迴避漫步（四個小圖示）；(9.4) 能量、(9.5) 總磁化。
  2. 8 個平面自旋一圈（WC.rotation_spins()，種子 902）全部轉 120°：E = −6.3851、|M| = 6.7943 不變，方向 −13.1° → 106.9°。
  3. 4 個自旋一圈（2×2 方塊）的 16 個組態 → Z = 2e^{4K} + 12 + 2e^{−4K}；P(M) 長條圖（K = 1、0.2）；⟨|m|⟩ = 0.9337、0.4871；2 個自旋 0.8808。
  4. Metropolis 軌跡（L = 8、T = 3.0；L = 4、T = 2.0 翻 30 次；L = 8、T = 2.0 翻 2 次）；平均翻轉等待 133、608、2151 sweep。
  5. 小外場 (9.6)：4 自旋 0.3506、0.0367；4×4 週期 0.9210、0.1584 → M_s = 0.9993；相圖（Fig. 9.2，無數字）。
畫面上的數字都由 utils.walk_crit 算出，check_math() 核對（Metropolis 用模組的預設種子 902）。
模組層級的畫圖函式：on_icons、ring_spins、mini_config、config_columns、pm_chart、trace_panel、phase_diagram。
"""
import math
import sys
import warnings
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_crit as WC
from utils import lattice_draw as L
from utils import walk_draw as D

ROT_DEG = 120
TINY_K = WC.TINY_K                     # (1.0, 0.2)
FIELD_H = WC.FIELD_H                   # (0.1, 0.01)
TRACES = WC.MC_TRACES                  # ((8, 3.0), (4, 2.0), (8, 2.0))
TAU_L = (4, 6, 8)
MS = (-4, -2, 0, 2, 4)
UP_COLOR = "#F2A33A"                   # 朝上：暖色
DOWN_COLOR = "#4A90D9"                 # 朝下：冷色
TRACE_COLORS = (GREY_A, TEAL_B, YELLOW)
GOOD = "#7CE38B"
WEIGHT_TEX = {4: "e^{4K}", 0: "1", -4: "e^{-4K}"}


def f4(x):
    return f"{float(x):.4f}"


def ring4_configs():
    """RING4（4 個自旋圍一圈＝2×2 方塊）的 16 個組態 [(s0..s3), M, Σ_鍵 s_i s_j]，依 M、再依鍵和由大到小。"""
    Ns, edges = WC.RING4
    out = []
    for b in range(2 ** Ns):
        sp = tuple(1 - 2 * ((b >> i) & 1) for i in range(Ns))
        out.append((sp, sum(sp), sum(sp[i] * sp[j] for i, j in edges)))
    return sorted(out, key=lambda c: (-c[1], -c[2]))


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        S, S2 = WC.rotation_spins()
        E, E2, M, M2 = WC.rotation_check()
        nb = {"S": np.array(S), "S2": np.array(S2), "E": (E, E2), "M": (np.array(M), np.array(M2)),
              "ring4": {K: WC.ising_exact(WC.RING4, K) for K in TINY_K},
              "ring4c": {K: WC.ring4_closed(K) for K in TINY_K},
              "pair": WC.pair_closed(1), "pair_ex": WC.ising_exact(WC.PAIR, 1.0),
              "traces": {lt: WC.mc_trace(*lt) for lt in TRACES},
              "flips": {lt: WC.trace_flips(*lt) for lt in TRACES},
              "tau": [WC.flip_tau(n) for n in TAU_L],
              "field": {H: (WC.ising_exact(WC.RING4, 1.0, H)["m"], WC.ising_exact(WC.TORUS4, 1.0, H)["m"])
                        for H in FIELD_H},
              "Ms1": WC.onsager_M(1.0),
              "configs": ring4_configs()}
    return nb


def tau_ratios(tau):
    return [f"{b / a:.1f}" for a, b in zip(tau, tau[1:])]


def check_math():
    nb = numbers()
    S, S2 = nb["S"], nb["S2"]
    assert S.shape == S2.shape == (8, 2) and np.allclose(np.linalg.norm(S, axis=1), 1)
    for a, b in zip(S, S2):
        assert abs((WC.angle_deg(b) - WC.angle_deg(a)) % 360 - ROT_DEG) < 1e-9
    (E, E2), (M, M2) = nb["E"], nb["M"]
    assert f4(E) == f4(E2) == "-6.3851" and abs(WC.on_energy(S, WC.ring_edges(8)) - E) < 1e-12
    assert f4(np.linalg.norm(M)) == f4(np.linalg.norm(M2)) == "6.7943" and np.allclose(M, S.sum(0))
    assert f"{WC.angle_deg(M):.1f}" == "-13.1" and f"{WC.angle_deg(M2):.1f}" == "106.9"
    want = {1.0: ("121.2329", ["0.4504", "0.0330", "0.0333", "0.0330", "0.4504"], "0.9337"),
            0.2: ("17.3497", ["0.1283", "0.2306", "0.2823", "0.2306", "0.1283"], "0.4871")}
    cf = nb["configs"]
    assert len(cf) == 16 and [c[1] for c in cf].count(0) == 6 and {c[2] for c in cf} == set(WEIGHT_TEX)
    for K in TINY_K:
        ex, cl = nb["ring4"][K], nb["ring4c"][K]
        Z, P, am = want[K]
        assert f4(ex["Z"]) == f4(cl[0]) == Z
        assert [f4(ex["PM"][m]) for m in MS] == P and f4(cl[2]) == P[4] and f4(cl[3]) == P[3] and f4(cl[4]) == P[2]
        assert f4(ex["absm"]) == f4(cl[1]) == am and abs(ex["m"]) < 1e-12
        Zc = sum(math.exp(K * c[2]) for c in cf)
        assert abs(Zc - ex["Z"]) < 1e-9 * Zc
        assert abs(Zc - (2 * math.exp(4 * K) + 12 + 2 * math.exp(-4 * K))) < 1e-9 * Zc
        assert abs(sum(math.exp(K * c[2]) for c in cf if c[1] == 4) / Zc - ex["PM"][4]) < 1e-12
    assert f"{1 / TINY_K[0]:g}" == "1" and f"{1 / TINY_K[1]:g}" == "5"
    Zp, amp = nb["pair"]
    assert f4(Zp) == f4(nb["pair_ex"]["Z"]) == "6.1723" and f4(amp) == f4(nb["pair_ex"]["absm"]) == "0.8808"
    for lt in TRACES:
        sw, m = nb["traces"][lt]
        assert len(sw) == len(m) == 400 and 0 <= sw[0] and sw[-1] <= 4000 and np.all(np.abs(m) <= 1)
    assert nb["flips"][TRACES[1]] == 30 and nb["flips"][TRACES[2]] == 2
    assert TRACES[0][1] > WC.TC_2D > TRACES[1][1] == TRACES[2][1] == WC.MC_T == 2.0
    assert (WC.MC_SWEEPS, WC.MC_REPLICAS, WC.MC_M0) == (4000, 50, 0.5)       # 字幕：4000 sweep、50 個複本、±0.5
    assert (len(S), WC.PAIR[0], WC.RING4[0], WC.TORUS4[0]) == (8, 2, 4, 16) and TAU_L == (4, 6, 8)
    assert [f"{v:.0f}" for v in nb["tau"]] == ["133", "608", "2151"] and tau_ratios(nb["tau"]) == ["4.6", "3.5"]
    assert [f4(nb["field"][H][0]) for H in FIELD_H] == ["0.3506", "0.0367"]
    assert [f4(nb["field"][H][1]) for H in FIELD_H] == ["0.9210", "0.1584"]
    assert f4(nb["Ms1"]) == "0.9993"


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


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def row_cells(tab, i, ncols):
    return tab[0][i * ncols:(i + 1) * ncols]


def v3(v):
    return np.array([v[0], v[1], 0.0])


def spin_arrow(center, direction, length, color, thickness=4):
    c, d = np.array(center, dtype=float), v3(direction) / np.linalg.norm(direction)
    return Arrow(c - 0.5 * length * d, c + 0.5 * length * d, buff=0, thickness=thickness,
                 max_width_to_length_ratio=0.16).set_fill(color, 1)


# ---------------------------------------------------------------- 畫圖元件
def on_icons(size=1.3):
    """O(n) 的四個小圖示（n = 1、2、3、0），各自置中在 ORIGIN → list。只是 2D 的示意圖。"""
    r = 0.5 * size
    ising = VGroup(spin_arrow(0.28 * size * LEFT, (0, 1), 0.9 * size, UP_COLOR, 5),
                   spin_arrow(0.28 * size * RIGHT, (0, -1), 0.9 * size, DOWN_COLOR, 5))
    xy = VGroup(Circle(radius=r).set_stroke(GREY_B, 2),
                Arrow(ORIGIN, r * np.array([math.cos(0.7), math.sin(0.7), 0]), buff=0, thickness=5).set_fill(UP_COLOR, 1))
    heis = VGroup(Circle(radius=r).set_stroke(GREY_B, 2),
                  DashedVMobject(Ellipse(width=2 * r, height=0.6 * r), num_dashes=24).set_stroke(GREY_B, 1.5),
                  Arrow(ORIGIN, r * np.array([0.45, 0.8, 0]), buff=0, thickness=5).set_fill(UP_COLOR, 1))
    pts = [(0, 0), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
    saw = L.trail(pts, 0.32 * size, ORIGIN, stroke_width=4, end_radius=0.06)
    return [ising, xy, heis, saw]


def ring_spins(S, center, radius=1.55, length=0.8):
    """Ns 個平面自旋排成一圈 → VGroup(鍵（多邊形）, 格點, 箭頭們)；箭頭中心在圓周上、方向 = S_i。"""
    c = np.array(center, dtype=float)
    Ns = len(S)
    pos = [c + radius * np.array([math.cos(PI / 2 - 2 * PI * k / Ns), math.sin(PI / 2 - 2 * PI * k / Ns), 0])
           for k in range(Ns)]
    bonds = Polygon(*pos).set_stroke(GREY_D, 2)
    sites = VGroup(*[Dot(p, radius=0.05).set_fill(GREY_B, 1) for p in pos])
    arrows = VGroup(*[spin_arrow(p, s, length, UP_COLOR) for p, s in zip(pos, S)])
    return VGroup(bonds, sites, arrows)


def mini_config(sp, center, cell=0.34):
    """RING4 的一個組態：格點 0→左下、1→右下、2→右上、3→左上（沿著圈）；VGroup(鍵的方框, 4 支箭頭)。"""
    c = np.array(center, dtype=float)
    off = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    pos = [c + cell * np.array([a, b, 0]) for a, b in off]
    box = Square(side_length=cell).move_to(c).set_stroke(GREY_D, 1.5)
    arrows = VGroup(*[spin_arrow(p, (0, s), 0.72 * cell, UP_COLOR if s > 0 else DOWN_COLOR, 4) for p, s in zip(pos, sp)])
    return VGroup(box, arrows)


def config_columns(configs, top_y=2.1, pitch_x=1.15, pitch_y=1.25, gap=0.55, cell=0.42):
    """16 個組態依 M 分成 5 欄（M = 4, 2, 0, −2, −4；組態多的欄排成兩行）。
    回傳 (heads, minis, weights)：每欄的 M 標題、每個組態、每個組態下的波茲曼權重（e^{4K}、1、e^{−4K}）。"""
    groups = [(M, [c for c in configs if c[1] == M]) for M in (4, 2, 0, -2, -4)]
    widths = [(1 if len(g) == 1 else 2) * pitch_x for _, g in groups]
    x = -(sum(widths) + gap * (len(groups) - 1)) / 2
    heads, minis, weights = VGroup(), VGroup(), VGroup()
    for (M, g), w in zip(groups, widths):
        ncol = 1 if len(g) == 1 else 2
        heads.add(tx(f"M={M}", 28, YELLOW).move_to([x + w / 2, top_y + 0.75, 0]))
        for k, (sp, _, b) in enumerate(g):
            cx = x + (k % ncol + 0.5) * pitch_x
            cy = top_y - (k // ncol) * pitch_y
            minis.add(mini_config(sp, [cx, cy, 0], cell))
            weights.add(tx(WEIGHT_TEX[b], 24, GREY_A).move_to([cx, cy - 0.6, 0]))
        x += w + gap
    return heads, minis, weights


def pm_chart(PM, center_x, base_y, color, ux=0.42, uy=4.4):
    """P(M) 長條圖（M = −4..4），附 M 刻度與底線 → VGroup(長條 VGroup(長條們, 標籤們), 刻度, 底線)。"""
    ch = D.bar_chart({m: PM[m] for m in MS}, ux=ux, uy=uy, origin=[center_x, base_y, 0], color=color,
                     labels=True, label_fmt=f4, font_size=20)
    ticks = VGroup(*[tx(str(m), 22, GREY_B).move_to([center_x + m * ux, base_y - 0.25, 0]) for m in MS])
    base = Line([center_x - 5.3 * ux, base_y, 0], [center_x + 5.3 * ux, base_y, 0]).set_stroke(GREY_B, 1.5)
    return VGroup(ch, ticks, base)


def trace_panel(origin, width=6.0, height=1.2, xticks=False):
    """m(t) 的小圖：橫軸 sweep 0–4000、縱軸 m ∈ [−1, 1]（畫成 m + 1 ∈ [0, 2]，軸才會在底下）；原點 (0, −1) 放在 origin。"""
    kw = dict(x_ticks=(0, 2000, 4000)) if xticks else {}
    ax = D.plot_axes((0, 4000, 1000), (0, 2, 1), width=width, height=height, y_ticks=(0, 1, 2),
                     y_tick_labels=("-1", "0", "1"), font_size=18, **kw)
    ax.shift(np.array(origin, dtype=float) - ax.c2p(0, 0))
    zero = DashedLine(ax.c2p(0, 1), ax.c2p(4000, 1)).set_stroke(GREY_C, 1)
    return ax, zero


def trace_curve(ax, sw, m, color):
    """一條軌跡＝一個 VMobject 折線（400 點）。"""
    return VMobject().set_points_as_corners([ax.c2p(x, y + 1) for x, y in zip(sw, m)]).set_stroke(color, 1.3)


def phase_diagram(center, width=5.4, height=4.4):
    """Fig. 9.2 的示意：橫軸 T、縱軸 h；h = 0、0 < T < T_c 的粗線（一級相變線）、終點的臨界點、越過粗線的跳躍。
    回傳 VGroup(軸們, 粗線, 臨界點, 跳躍箭頭, 標籤們)。沒有數字。"""
    c = np.array(center, dtype=float)
    o = c + width / 2 * LEFT
    t_axis = Arrow(o, o + width * RIGHT, buff=0, thickness=2.5).set_fill(GREY_B, 1)
    h_axis = Arrow(o + height / 2 * DOWN, o + height / 2 * UP, buff=0, thickness=2.5).set_fill(GREY_B, 1)
    tc = o + 0.58 * width * RIGHT
    line = Line(o, tc).set_stroke(YELLOW, 9)
    cp = Dot(tc, radius=0.11).set_fill(RED, 1)
    jx = o + 0.2 * width * RIGHT
    jump = Arrow(jx + 0.85 * DOWN, jx + 0.85 * UP, buff=0, thickness=3).set_fill(WHITE, 1)
    labs = VGroup(tx("T", 30, GREY_B).next_to(t_axis.get_end(), DOWN, buff=0.2),
                  tx("h", 30, GREY_B).next_to(h_axis.get_end(), LEFT, buff=0.2),
                  tx("T_c", 28, RED_B).next_to(cp, DOWN, buff=0.15),
                  zh("臨界點", font_size=24, color=RED_B).next_to(cp, UR, buff=0.08),
                  mixed("$h>0$", "：", "$M>0$", font_size=24).move_to(o + 0.3 * width * RIGHT + 0.36 * height * UP),
                  mixed("$h<0$", "：", "$M<0$", font_size=24).move_to(o + 0.3 * width * RIGHT + 0.36 * height * DOWN),
                  mixed("$-M_s\\to+M_s$", font_size=24).next_to(jump, RIGHT, buff=0.12).shift(0.45 * UP),
                  mixed("$T>T_c$", "：連續", font_size=24).move_to(o + 0.82 * width * RIGHT + 0.22 * height * UP))
    return VGroup(VGroup(t_axis, h_axis), line, cp, jump, labs)


# ---------------------------------------------------------------- 各段
def part_on_model(s):
    icons = on_icons()
    names = [("$n=1$", " Ising", "上、下"), ("$n=2$", " xy", "平面上任意方向"),
             ("$n=3$", " Heisenberg", "空間任意方向"), ("$n=0$", "：自迴避漫步", "W906 會說明")]
    cards = VGroup()
    for icon, (a, b, c) in zip(icons, names):
        lab = VGroup(mixed(a, b, font_size=26), zh(c, font_size=22, color=GREY_B)).arrange(DOWN, buff=0.1)
        cards.add(VGroup(icon, lab).arrange(DOWN, buff=0.25))
    cards.arrange(RIGHT, buff=0.9)
    e94 = mixed("(9.4)　", "$\\mathcal H=-J\\sum_{\\langle ij\\rangle}\\vec S_i\\cdot\\vec S_j$", font_size=32)
    e95 = mixed("(9.5)　", "$\\vec M=\\sum_i\\vec S_i$", font_size=32)
    enorm = mixed("$|\\vec S_i|=1$", "，", "$\\vec S_i$", " 有 ", "$n$", " 個分量", font_size=28, color=GREY_B)
    eqs = VGroup(VGroup(e94, e95).arrange(RIGHT, buff=1.0), enorm).arrange(DOWN, buff=0.3)
    place_in(cards, -6.7, 6.7, 0.0, 3.1)
    place_in(eqs, -6.7, 6.7, -2.8, -0.5)
    assert_in_frame(cards, eqs)
    assert_apart(cards, eqs)
    for a, b in zip(cards, cards[1:]):
        assert_apart(a, b)
    s.say("問：能量在「所有自旋一起轉」之下不變")
    s.say("所以平均磁化一定是 0——那永久磁鐵是怎麼來的？", wait=2.5)
    s.say("$\\mathrm{O}(n)$", " 模型：每個自旋是 ", "$n$", " 個分量的單位向量")
    s.play(FadeIn(cards[0]), FadeIn(cards[1]), run_time=1.0)
    s.say("$n=1$", " Ising：上或下；", "$n=2$", " xy：平面上任意方向", wait=2.0)
    s.play(FadeIn(cards[2]), FadeIn(cards[3]), run_time=1.0)
    s.say("$n=3$", " Heisenberg：空間任意方向；", "$n=0$", "：自迴避漫步", wait=2.5)
    s.play(FadeIn(eqs), run_time=1.0)
    s.say("(9.4) 能量：相鄰自旋的內積；(9.5) 總磁化是向量和", wait=3.0)
    s.clear_stage()


def part_rotation(s):
    nb = numbers()
    S = nb["S"]
    (E, _), (M, M2) = nb["E"], nb["M"]
    c = np.array([-3.3, 0.15, 0])
    ring = ring_spins(S, c)
    m_arrow = Arrow(c, c + 0.15 * v3(M), buff=0, thickness=6).set_fill(YELLOW, 1)
    leg = mixed("黃色：總磁化 ", "$\\vec M$", "（縮小畫）", font_size=24, color=GREY_B).move_to(c + 2.55 * DOWN)
    vE = mixed("能量 ", f"$E={f4(E)}$", font_size=32)
    vM = mixed("總磁化 ", f"$|\\vec M|={f4(np.linalg.norm(M))}$", font_size=32)
    th0 = f"{WC.angle_deg(M):.1f}^\\circ"
    th1 = f"{WC.angle_deg(M2):.1f}^\\circ"
    vT0 = mixed("方向 ", f"${th0}$", font_size=32)
    vT1 = mixed("方向 ", f"${th0}\\to{th1}$", font_size=32)
    vals = VGroup(vE, vM, vT0).arrange(DOWN, buff=0.5, aligned_edge=LEFT).move_to([3.2, 0.6, 0])
    vT1.move_to(vT0, aligned_edge=LEFT)
    tags = VGroup(zh("不變", font_size=28, color=GOOD).next_to(vE, RIGHT, buff=0.35),
                  zh("不變", font_size=28, color=GOOD).next_to(vM, RIGHT, buff=0.35))
    concl = mixed("各個方向一樣好 ", "$\\Rightarrow\\ \\langle\\vec M\\rangle=0$", font_size=32, color=YELLOW)
    concl.move_to([3.2, -1.9, 0])
    assert_in_frame(ring, m_arrow, leg, vals, vT1, tags, concl)
    assert_apart(ring, vals)
    assert_apart(ring, leg)
    assert_apart(vals, concl)
    assert_apart(vT1, concl)
    assert_column(vE, vM, vT1)
    s.say("8 個平面自旋排成一圈（", "$J=1$", "），箭頭是 ", "$\\vec S_i$")
    s.play(ShowCreation(ring[0]), FadeIn(ring[1]), FadeIn(ring[2], lag_ratio=0.1), run_time=1.5)
    s.play(GrowArrow(m_arrow), FadeIn(leg), FadeIn(vals), run_time=1.2)
    s.say("能量 ", f"$E={f4(E)}$", "，總磁化 ", f"$|\\vec M|={f4(np.linalg.norm(M))}$", wait=2.5)
    s.say(f"全部一起轉 {ROT_DEG} 度……")
    rots = [Rotate(a, ROT_DEG * DEGREES, about_point=a.get_center().copy()) for a in ring[2]]
    s.play(*rots, Rotate(m_arrow, ROT_DEG * DEGREES, about_point=c.copy()), run_time=3.0)
    s.play(FadeTransform(vT0, vT1), FadeIn(tags), run_time=1.0)
    s.say("能量不變、", "$|\\vec M|$", " 不變，方向從 ", f"${th0}$", " 變成 ", f"${th1}$", wait=3.0)
    s.play(FadeIn(concl), run_time=0.8)
    s.say("各個方向一樣好：平均起來 ", "$\\langle\\vec M\\rangle=0$", wait=3.0)
    s.clear_stage()


def part_tiny(s):
    nb = numbers()
    heads, minis, weights = config_columns(nb["configs"])
    eZ = mixed("16 個組態的權重 ", "$e^{K\\sum s_is_j}$", " 加起來：", "$Z=2e^{4K}+12+2e^{-4K}$", font_size=30)
    zs = VGroup(*[mixed(f"$K={K:g}$", "（", f"$T={1 / K:g}$", "）：", f"$Z={f4(nb['ring4'][K]['Z'])}$", font_size=30)
                  for K in TINY_K]).arrange(RIGHT, buff=1.2)
    low = VGroup(eZ, zs).arrange(DOWN, buff=0.3)
    place_in(low, -6.7, 6.7, -2.9, -1.35)
    top = VGroup(heads, minis, weights)
    assert_in_frame(top, low)
    assert_apart(top, low)
    for a, b in zip(minis, minis[1:]):
        assert_apart(a, b)
    for w, mn in zip(weights, minis):
        assert_apart(w, mn, 0.02)
    for h in heads:
        assert_apart(h, minis)
    s.say("小系統可以逐一列舉：4 個自旋圍一圈（", "$2\\times2$", " 方塊）")
    s.play(FadeIn(heads), FadeIn(minis, lag_ratio=0.05), run_time=2.0)
    s.say("$2^4=16$", " 個組態，依總磁化 ", "$M$", " 分成五欄", wait=2.5)
    s.play(FadeIn(weights, lag_ratio=0.05), run_time=1.2)
    s.say("全部同向：", "$e^{4K}$", "；斷兩根鍵：1；對角兩個翻：", "$e^{-4K}$", wait=2.5)
    s.play(FadeIn(eZ), run_time=0.8)
    s.say("$Z=2e^{4K}+12+2e^{-4K}$", wait=2.0)
    s.play(FadeIn(zs), run_time=0.8)
    s.say("$K=J/T$", "：低溫 ", "$K=1$", "、高溫 ", "$K=0.2$", wait=2.0)
    s.clear_stage()

    c1 = pm_chart(nb["ring4"][TINY_K[0]]["PM"], -3.4, -0.9, BLUE_C)
    c2 = pm_chart(nb["ring4"][TINY_K[1]]["PM"], 3.4, -0.9, ORANGE)
    t1 = mixed(f"$K={TINY_K[0]:g}$", "（", f"$T={1 / TINY_K[0]:g}$", "）", font_size=30).move_to([-3.4, 2.15, 0])
    t2 = mixed(f"$K={TINY_K[1]:g}$", "（", f"$T={1 / TINY_K[1]:g}$", "）", font_size=30).move_to([3.4, 2.15, 0])
    a1 = tx(r"\langle|m|\rangle=%s" % f4(nb["ring4"][TINY_K[0]]["absm"]), 30, YELLOW).move_to([-3.4, -1.85, 0])
    a2 = tx(r"\langle|m|\rangle=%s" % f4(nb["ring4"][TINY_K[1]]["absm"]), 30, YELLOW).move_to([3.4, -1.85, 0])
    ylab = mixed("$P(M)$", font_size=26, color=GREY_B).move_to([0, 1.4, 0])
    pair = mixed("2 個自旋（", "$K=1$", "）：", f"$Z={f4(nb['pair'][0])}$", "，",
                 f"$\\langle|m|\\rangle={f4(nb['pair'][1])}$", font_size=28, color=GREY_A)
    place_in(pair, -6.7, 6.7, -2.9, -2.35)
    assert_in_frame(c1, c2, t1, t2, a1, a2, ylab, pair)
    for a, b in ((c1, c2), (t1, c1), (t2, c2), (a1, c1), (a2, c2), (a1, pair), (a2, pair), (ylab, c1), (ylab, c2)):
        assert_apart(a, b)
    for ch in (c1, c2):
        for a, b in zip(ch[0][1], ch[0][1][1:]):
            assert_apart(a, b, 0.03)
    s.say("$P(M)$", "：總磁化是 ", "$M$", " 的機率")
    s.play(FadeIn(t1), FadeIn(t2), FadeIn(ylab), ShowCreation(c1[2]), ShowCreation(c2[2]), FadeIn(c1[1]), FadeIn(c2[1]),
           run_time=1.0)
    s.play(GrowFromEdge(c1[0][0], DOWN), GrowFromEdge(c2[0][0], DOWN), run_time=1.2)
    s.play(FadeIn(c1[0][1]), FadeIn(c2[0][1]), run_time=0.8)
    s.say("低溫：兩端高高的；高溫：中間最高", wait=2.5)
    s.say("兩張圖都左右對稱：", "$\\langle M\\rangle=0$", wait=2.5)
    s.play(FadeIn(a1), FadeIn(a2), run_time=0.8)
    s.say("可是 ", f"$\\langle|m|\\rangle={f4(nb['ring4'][TINY_K[0]]['absm'])}$", " 對 ",
          f"${f4(nb['ring4'][TINY_K[1]]['absm'])}$", "：低溫時不是正就是負", wait=3.0)
    s.play(FadeIn(pair), run_time=0.8)
    s.say("最小的例子：2 個自旋也一樣，", f"$\\langle|m|\\rangle={f4(nb['pair'][1])}$", wait=2.5)
    s.clear_stage()


def part_metropolis(s):
    nb = numbers()
    origins = ([-6.0, 1.7, 0], [-6.0, -0.1, 0], [-6.0, -1.9, 0])
    panels, curves, labels = [], [], []
    words = ("高溫：在 0 附近亂跳", f"翻 {nb['flips'][TRACES[1]]} 次", f"只翻 {nb['flips'][TRACES[2]]} 次")
    for k, (lt, o, col, w) in enumerate(zip(TRACES, origins, TRACE_COLORS, words)):
        ax, zero = trace_panel(o, xticks=(k == 2))
        sw, m = nb["traces"][lt]
        panels.append(VGroup(ax, zero))
        curves.append(trace_curve(ax, sw, m, col))
        lab = VGroup(tx(f"L={lt[0]},\\ T={lt[1]:.1f}", 28, col), zh(w, font_size=24, color=col)).arrange(DOWN, buff=0.12)
        lab.move_to([1.55, o[1] + 0.6, 0])
        labels.append(lab)
    xl = zh("sweep", font_size=22, color=GREY_B).next_to(panels[2][0].c2p(2000, 0), DOWN, buff=0.42)
    mls = VGroup(*[tx("m", 28, GREY_B).move_to([-6.6, o[1] + 0.6, 0]) for o in origins])
    tau = nb["tau"]
    rat = tau_ratios(tau)
    tab = D.tex_table([["L", r"\tau", VectorizedPoint()],
                       [str(TAU_L[0]), f"{tau[0]:.0f}", VectorizedPoint()],
                       [str(TAU_L[1]), f"{tau[1]:.0f}", r"\times%s" % rat[0]],
                       [str(TAU_L[2]), f"{tau[2]:.0f}", r"\times%s" % rat[1]]], [0.9, 1.4, 1.2], row_h=0.7, font_size=32)
    for i in (2, 3):
        D.table_cell(tab, i, 2, 3).set_color(YELLOW)
    ttl = VGroup(mixed("$T=2.0$", "：平均翻轉等待 ", "$\\tau$", font_size=26),
                 zh("（50 個複本，單位 sweep）", font_size=22, color=GREY_B)).arrange(DOWN, buff=0.1)
    right = VGroup(ttl, tab).arrange(DOWN, buff=0.3)
    place_in(right, 3.1, 6.7, -1.6, 2.4)
    allp = VGroup(*panels, *curves)
    assert_in_frame(allp, xl, mls, *labels, right)
    for a, b in zip(panels, panels[1:]):
        assert_apart(a, b)
    for m_lab, p in zip(mls, panels):
        assert_apart(m_lab, p)
    assert_apart(xl, panels[2])
    for lab in labels:
        assert_apart(lab, allp)
        assert_apart(lab, right)
    assert_apart(allp, right)
    assert_table(tab, 3)
    s.say("Metropolis 模擬：", "$L\\times L$", " 週期格子，從全部朝上開始")
    s.play(*[FadeIn(p) for p in panels], FadeIn(mls), FadeIn(xl), run_time=1.0)
    s.say("縱軸 ", "$m$", "：每個自旋的平均磁化；橫軸：sweep 數", wait=1.5)
    s.play(ShowCreation(curves[0]), FadeIn(labels[0]), run_time=2.5, rate_func=linear)
    s.say("高溫 ", f"$T={TRACES[0][1]:.1f}>T_c$", "：", "$m$", " 在 0 附近亂跳", wait=2.0)
    s.play(ShowCreation(curves[1]), FadeIn(labels[1]), run_time=2.5, rate_func=linear)
    s.say("低溫 ", "$T=2.0$", "、", "$L=4$", f"：4000 sweep 裡翻了 {nb['flips'][TRACES[1]]} 次", wait=2.0)
    s.play(ShowCreation(curves[2]), FadeIn(labels[2]), run_time=2.5, rate_func=linear)
    s.say("$L=8$", f"：只翻 {nb['flips'][TRACES[2]]} 次，大部分時間停在同一邊", wait=2.5)
    s.say("翻轉：", "$m$", " 從大於 0.5 跑到小於 ", "$-0.5$", "（或反過來）", wait=1.5)
    s.play(FadeIn(ttl), FadeIn(row_cells(tab, 0, 3)), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][3:]], lag_ratio=0.12), run_time=1.2)
    s.say("平均等待：", "、".join(f"{v:.0f}" for v in tau), " sweep", wait=2.5)
    s.say("$L$", " 每加 2，等待時間大約乘 4", wait=2.0)
    s.say("書上說等待時間是「天文數字」，而且系統越大越久", wait=2.0)
    s.say("無限大的系統：永遠等不到翻轉", wait=2.0)
    s.clear_stage()


def part_field(s):
    nb = numbers()
    e96 = mixed("(9.6)　", "$\\langle\\vec M\\rangle=-\\frac{\\partial F}{\\partial\\vec h}$", font_size=32)
    sub = mixed("$K=1$", "（", "$T=1$", "）時的 ", "$\\langle m\\rangle$", font_size=26, color=GREY_B)
    hcell = mixed("$H=h/T$", font_size=28, color=YELLOW)
    tab = D.tex_table([[hcell] + [f"{H:g}" for H in FIELD_H],
                       [zh("4 個自旋", font_size=28)] + [f4(nb["field"][H][0]) for H in FIELD_H],
                       [mixed("$4\\times4$", " 週期", font_size=28)] + [f4(nb["field"][H][1]) for H in FIELD_H]],
                      [2.4, 1.6, 1.6], row_h=0.72, font_size=32)
    for j in (1, 2):
        D.table_cell(tab, 2, j, 3).set_color(YELLOW)
    note = mixed("無限大格子：", f"$M_s={f4(nb['Ms1'])}$", font_size=30, color=YELLOW)
    left = VGroup(e96, sub, tab, note).arrange(DOWN, buff=0.35)
    place_in(left, -6.7, -0.4, -2.8, 3.1)
    pd = phase_diagram([3.45, 0.2, 0])
    assert_in_frame(left, pd)
    assert_apart(left, pd)
    assert_column(e96, sub, tab, note)
    assert_table(tab, 3)
    labs = pd[4]
    items = list(labs) + [pd[2], pd[3]]
    for k, a in enumerate(items):
        for b in items[k + 1:]:
            assert_apart(a, b, 0.03)
    for lab in labs:
        for other in (pd[0][0], pd[0][1], pd[1]):
            assert_apart(lab, other, 0.03)
    s.say("加一個小外場 ", "$h$", "：(9.6) ", "$\\langle\\vec M\\rangle=-\\partial F/\\partial\\vec h$")
    s.play(FadeIn(e96), run_time=0.8)
    s.say("無因次外場 ", "$H=h/T$", "；溫度 ", "$T=1$", "（", "$K=1$", "）", wait=1.5)
    s.play(FadeIn(sub), FadeIn(row_cells(tab, 0, 3)), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][3:]], lag_ratio=0.12), run_time=1.2)
    s.say("4 個自旋：", f"$H={FIELD_H[0]:g}$", f" 只有 {f4(nb['field'][FIELD_H[0]][0])}；",
          "$4\\times4$", f" 就有 {f4(nb['field'][FIELD_H[0]][1])}", wait=2.5)
    s.play(FadeIn(note), run_time=0.8)
    s.say("格子越大，同一個小外場越接近 ", f"$M_s={f4(nb['Ms1'])}$", wait=2.5)
    s.say("先讓格子變成無限大、再讓 ", "$h\\to0^+$", "：得到 ", "$M_s\\neq0$", wait=2.5)
    s.play(FadeIn(pd[0]), FadeIn(labs[0]), FadeIn(labs[1]), run_time=0.8)
    s.play(ShowCreation(pd[1]), FadeIn(labs[4]), FadeIn(labs[5]), run_time=1.0)
    s.say("相圖：", "$h=0$", "、", "$0<T<T_c$", " 的粗線是一級相變線", wait=2.5)
    s.play(GrowArrow(pd[3]), FadeIn(labs[6]), run_time=1.0)
    s.say("越過粗線，", "$M$", " 從 ", "$-M_s$", " 跳到 ", "$+M_s$", wait=2.5)
    s.play(FadeIn(pd[2]), FadeIn(labs[2]), FadeIn(labs[3]), FadeIn(labs[7]), run_time=0.8)
    s.say("粗線的終點 ", "$(T_c,0)$", " 是臨界點；", "$T>T_c$", " 時 ", "$M$", " 連續變化", wait=2.5)
    s.say("對稱說平均是 0；可是低溫時翻一次要等很久", wait=2.5)
    s.say("格子越大越久，無限大就永遠不翻：這就是對稱破缺", wait=3.0)


class W902BrokenSymmetry(ChuScene):
    def construct(self):
        self.show_header("W902　O(n) 模型與對稱破缺")
        check_math()
        part_on_model(self)
        part_rotation(self)
        part_tiny(self)
        part_metropolis(self)
        part_field(self)
        self.wait(1.0)
