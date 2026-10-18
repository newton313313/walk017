"""
W606 1/d 展開：尖峰變寬（§6.2.2–6.2.5，pp. 142–152，(6.38)–(6.71)，Fig. 6.6–6.14，Ex. 6.7）

問：在 100 維，三個最大的主迴轉半徑是固定的數嗎？20 維、3 維呢？
  1. Wick 配對（Ex. 6.7，SH.wick_pairings、SH.wick_count）：(2m−1)!! = 1、3、15、105、945；m = 2 的三種配對＝⟨x⁴⟩ = 3。
  2. 預解式 (6.38)–(6.42)：Im R(λ − iε) = π Σ δ(λ − λ_j)；筆誤框 (6.52) 漏了 π（SH.resolvent_weight = 3.1416）；
     零階 (6.50)：δ 尖峰在 1/(π²k²) = 0.1013、0.0253、0.0113（除以 N+1）。
  3. 一階 (6.58)–(6.62)：位移 3/(4d)（SH.sum_identity → 3/4）、相對寬度 √(2/d)。
  4. Fig. 6.12–6.14（SH.eig_sim(d)，N = 50）三張直方圖（每色 30 格、以最高點為 1）＋SH.onedee_table(d) 的比值。
畫面上的數字都由 utils.walk_shape 算出，check_math() 核對。
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
from utils import walk_shape as SH
from utils import walk_draw as D

DIMS = (100, 20, 3)
K_COLORS = [BLUE_C, ORANGE, PINK]
XMAX = {100: 0.16, 20: 0.25, 3: 0.5}
XTICKS = {100: [0.05, 0.1, 0.15], 20: [0.1, 0.2], 3: [0.1, 0.2, 0.3, 0.4]}
NBINS = 30
EPS = (0.4, 0.2, 0.1)
EPS_COLORS = [TEAL_C, YELLOW, RED_B]
BOOK = {100: ("var", (0.97, 0.93, 0.95)), 3: ("both", (0.98, 0.48, 0.23), (1.01, 0.91, 0.75))}


def f4(x):
    return f"{x:.4f}"


def f2s(vals):
    return [f"{v:.2f}" for v in vals]


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    nb = {"wick": [SH.wick_count(m) for m in range(1, 6)], "pairs2": SH.wick_pairings(2),
          "x4": SH.gauss_moment(2), "res": SH.resolvent_weight(), "res_eps": [SH.resolvent_weight(e) for e in EPS],
          "lead": [SH.eigs_asym(SH.N_D, k) / (SH.N_D + 1) for k in (1, 2, 3)],
          "sumid": [SH.sum_identity(n) for n in (1, 2, 3)],
          "shift": [SH.shift_1d(d) for d in DIMS], "width": [math.sqrt(2 / d) for d in DIMS]}
    nb["table"] = {d: SH.onedee_table(d) for d in DIMS}
    return nb


@lru_cache(maxsize=None)
def hist_shape(d, k):
    """第 k 個特徵值（除以 N+1）的直方圖：30 格，範圍取 0.1%–99.9% 分位數（不超過 XMAX），高度以最高格為 1。"""
    v = SH.eig_sim(d)[0][:, k - 1]
    lo = max(0.0, float(np.quantile(v, 0.001)))
    hi = min(XMAX[d], float(np.quantile(v, 0.999)))
    cnt, edges = np.histogram(v, bins=NBINS, range=(lo, hi))
    return edges, cnt / cnt.max()


def check_math():
    nb = numbers()
    assert nb["wick"] == [1, 3, 15, 105, 945]
    assert all(len(SH.wick_pairings(m)) == SH.wick_count(m) for m in range(1, 5))
    assert nb["pairs2"] == [[(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]]
    assert abs(nb["x4"] - 3) < 1e-8
    assert f4(nb["res"]) == "3.1416" and all(abs(w - math.pi) < 1e-6 for w in nb["res_eps"])
    assert [f4(v) for v in nb["lead"]] == ["0.1013", "0.0253", "0.0113"]
    assert all(abs(nb["lead"][k - 1] - 1 / (math.pi**2 * k * k)) < 1e-15 for k in (1, 2, 3))
    assert all(f4(v) == "0.7500" for v in nb["sumid"])
    assert [f4(v) for v in nb["shift"]] == ["0.0075", "0.0375", "0.2500"]
    assert [f4(v) for v in nb["width"]] == ["0.1414", "0.3162", "0.8165"]
    tb = nb["table"]
    assert f2s(tb[100]["mean_ratio"]) == ["1.01", "1.00", "1.00"] and f2s(tb[100]["var_ratio"]) == ["0.98", "0.97", "0.94"]
    assert f2s(tb[20]["mean_ratio"]) == ["1.01", "1.01", "1.00"] and f2s(tb[20]["var_ratio"]) == ["1.02", "0.87", "0.68"]
    assert f2s(tb[3]["mean_ratio"]) == ["1.01", "0.90", "0.75"] and f2s(tb[3]["var_ratio"]) == ["1.01", "0.47", "0.23"]
    r100 = sorted(f2s(tb[100]["mean_ratio"]) + f2s(tb[100]["var_ratio"]))
    assert (r100[0], r100[-1]) == ("0.94", "1.01")                                          # 字幕「0.94–1.01」
    assert np.all(np.abs(tb[100]["var_ratio"] - np.array(BOOK[100][1])) <= 0.05)          # 書上 (6.63)–(6.71)
    assert np.all(np.abs(tb[3]["var_ratio"] - np.array(BOOK[3][1])) <= 0.05)
    assert np.all(np.abs(tb[3]["mean_ratio"] - np.array(BOOK[3][2])) <= 0.05)
    for d in DIMS:
        top, _ = SH.eig_sim(d)
        assert top.shape == (SH.M_D[d], 3)
        assert np.allclose(tb[d]["pred_mean"], [SH.mean_1d(k, d) / (SH.N_D + 1) for k in (1, 2, 3)])
        assert np.allclose(tb[d]["lead"], nb["lead"])
        for k in (1, 2, 3):
            edges, h = hist_shape(d, k)
            assert len(h) == NBINS and edges[-1] <= XMAX[d] + 1e-12 and h.max() == 1.0
    assert SH.M_D[100] == SH.M_D[20] == 2000 and SH.M_D[3] == 20000 and SH.N_D == 50


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


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
        assert_apart(a, b)


def spike_axes():
    """d → ∞ 的三根 δ 尖峰（除以 N+1）。"""
    nb = numbers()
    ax = D.plot_axes((0, 0.12, 0.02), (0, 1.3, 0.5), width=7.0, height=3.4, x_ticks=[0.05, 0.1], font_size=20,
                     x_label=r"R_k^2/(N+1)")
    spikes, labs = VGroup(), VGroup()
    for k, x in enumerate(nb["lead"]):
        a = Arrow(ax.c2p(x, 0), ax.c2p(x, 1.0), buff=0, thickness=3).set_color(K_COLORS[k])
        spikes.add(a)
        labs.add(tx(f4(x), 22, K_COLORS[k]).next_to(a, UP, buff=0.1))
    return VGroup(ax, spikes, labs)


def pairing_diagram(pairs, gap=0.75):
    dots = VGroup(*[Dot(i * gap * RIGHT, radius=0.07).set_color(WHITE) for i in range(4)])
    nums = VGroup(*[tx(str(i + 1), 22, GREY_B).next_to(dots[i], DOWN, buff=0.12) for i in range(4)])
    arcs = VGroup(*[ArcBetweenPoints(dots[a].get_center(), dots[b].get_center(), angle=-PI * 0.75)
                    .set_stroke(YELLOW, 3) for a, b in pairs])
    expr = tx("".join(r"\langle\eta_%d\eta_%d\rangle" % (a + 1, b + 1) for a, b in pairs), 24)
    expr.next_to(nums, DOWN, buff=0.2)
    return VGroup(dots, nums, arcs, expr)


def lorentz_plot():
    ax = D.plot_axes((-2, 2, 1), (0, 10, 2), width=5.6, height=3.6, x_ticks=[-2, -1, 1, 2], y_ticks=[5, 10],
                     font_size=20, x_label="x")
    curves, labs = VGroup(), VGroup()
    for e, c in zip(EPS, EPS_COLORS):
        curves.add(ax.get_graph(lambda x, e=e: e / (x * x + e * e), x_range=(-2, 2, 0.01)).set_stroke(c, 3))
        labs.add(tx(r"\epsilon=%.1f" % e, 22, c))
    xs = np.arange(-2, 2.001, 0.02)
    e = EPS[1]
    area = Polygon(*([ax.c2p(-2, 0)] + [ax.c2p(x, e / (x * x + e * e)) for x in xs] + [ax.c2p(2, 0)]))
    area.set_fill(YELLOW, 0.25).set_stroke(width=0)
    labs.arrange(DOWN, buff=0.15, aligned_edge=LEFT).move_to(ax.c2p(1.35, 7.5))
    return VGroup(ax, area, curves, labs)


def hist_panel(d):
    nb = numbers()
    ax = D.plot_axes((0, XMAX[d], 0.01), (0, 1.15, 0.5), width=3.9, height=2.5, x_ticks=XTICKS[d], font_size=18)
    polys = VGroup()
    for k in (1, 2, 3):
        edges, h = hist_shape(d, k)
        pts = [ax.c2p(edges[0], 0)]
        for i in range(NBINS):
            pts += [ax.c2p(edges[i], h[i]), ax.c2p(edges[i + 1], h[i])]
        pts.append(ax.c2p(edges[-1], 0))
        pts = [p for i, p in enumerate(pts) if i == 0 or np.linalg.norm(p - pts[i - 1]) > 1e-9]   # 去掉重複的頂點
        polys.add(Polygon(*pts).set_fill(K_COLORS[k - 1], 0.45).set_stroke(K_COLORS[k - 1], 2))
    lines = VGroup(*[DashedLine(ax.c2p(x, 0), ax.c2p(x, 1.12), dash_length=0.06).set_stroke(WHITE, 1.5)
                     for x in nb["lead"]])
    title = mixed("$d=%d$" % d, "（%d 條）" % SH.M_D[d], font_size=22).next_to(ax, UP, buff=0.12)
    t = nb["table"][d]
    l1 = mixed("平均除以 (6.61)：", "$%s$" % r",\ ".join(f2s(t["mean_ratio"])), font_size=20)
    l2 = mixed("變異數除以 (6.62)：", "$%s$" % r",\ ".join(f2s(t["var_ratio"])), font_size=20)
    txt = VGroup(l1, l2).arrange(DOWN, buff=0.18, aligned_edge=LEFT)
    return VGroup(ax, polys, lines, title), txt


# ---------------------------------------------------------------- 各段
def part_question(s):
    sp = spike_axes()
    q = VGroup(mixed("$d\\to\\infty$", "：三根尖峰（W604）", font_size=28),
               mixed("$d=100$", "：還是尖峰？", font_size=28),
               mixed("$d=20$", "？", font_size=28),
               mixed("$d=3$", "？", font_size=28, color=YELLOW)).arrange(DOWN, buff=0.55, aligned_edge=LEFT)
    place_in(sp, -6.7, 1.2, -2.8, 2.9)
    place_in(q, 1.8, 6.7, -2.5, 2.5)
    assert_in_frame(sp, q)
    assert_apart(sp, q)
    assert_column(*q)
    for a, b in zip(sp[2], sp[2][1:]):
        assert_apart(a, b)
    s.say("問：在 100 維，三個最大的主迴轉半徑是固定的數嗎？", wait=1.0)
    s.play(ShowCreation(sp[0]), run_time=0.8)
    s.play(LaggedStart(*[GrowArrow(a) for a in sp[1]], lag_ratio=0.3), FadeIn(sp[2]), FadeIn(q[0]), run_time=1.2)
    s.say("維度無限大（W604）：每條都是樓梯，", "$R_k^2=\\frac{N+1}{\\pi^2k^2}$", wait=3.0)
    s.say("除以 ", "$N+1$", "：三根尖峰 0.1013、0.0253、0.0113", wait=2.5)
    s.play(FadeIn(q[1:]), run_time=0.8)
    s.say("20 維、3 維呢？尖峰會不會變寬、疊在一起？（先猜）", wait=3.5)
    s.clear_stage()


def part_wick(s):
    nb = numbers()
    rows = [["m", "1", "2", "3", "4", "5"], ["(2m-1)!!"] + [str(v) for v in nb["wick"]]]
    tab = D.tex_table(rows, [2.0, 1.1, 1.1, 1.1, 1.1, 1.1], row_h=0.62, font_size=30)
    diags = VGroup(*[pairing_diagram(p) for p in nb["pairs2"]]).arrange(RIGHT, buff=1.4, aligned_edge=UP)
    bot = mixed("$\\langle\\eta_1\\eta_2\\eta_3\\eta_4\\rangle$", "＝三種配對的和；全都相同時：",
                "$\\langle x^4\\rangle=3$", font_size=26, color=YELLOW)
    place_in(tab, -6.75, 6.75, 1.7, 3.1)
    place_in(diags, -6.75, 6.75, -1.5, 1.3)
    place_in(bot, -6.75, 6.75, -2.9, -2.0)
    assert_in_frame(tab, diags, bot)
    assert_column(tab, diags, bot)
    assert_table(tab, 6)
    for a, b in zip(diags, diags[1:]):
        assert_apart(a, b)
    s.play(FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=0.8)
    s.say("Wick 配對：高斯的高次平均＝所有「兩兩配對」的和", wait=2.0)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:]], lag_ratio=0.25), run_time=1.2)
    s.say("$2m$", " 個 ", "$\\eta$", " 配成 ", "$m$", " 對：共 ", "$(2m-1)!!$", " 種", wait=3.0)
    for dg in diags:
        s.play(FadeIn(dg[0]), FadeIn(dg[1]), ShowCreation(dg[2]), FadeIn(dg[3]), run_time=0.9)
    s.say("$m=2$", "：四個點，三種連法", wait=2.0)
    s.play(FadeIn(bot), run_time=0.8)
    s.say("四個都是同一個 ", "$x$", " 時：標準高斯的 ", "$\\langle x^4\\rangle=3$", wait=2.5)
    s.say("1/d 展開：把這些配對按 ", "$1/d$", " 的次方分類", wait=2.5)
    s.clear_stage()


def part_resolvent(s):
    nb = numbers()
    r1 = eq_label(tx(r"\mathcal R(\lambda)=\mathrm{Tr}\,\frac{1}{\lambda\mathbf I-\mathbf T}"
                     r"=\sum_j\frac{1}{\lambda-\lambda_j}", 32), "(6.38)--(6.40)")
    r2 = eq_label(tx(r"\mathrm{Im}\,\mathcal R(\lambda-i\epsilon)=\sum_j\frac{\epsilon}{(\lambda-\lambda_j)^2+\epsilon^2}", 32),
                  "(6.41)")
    r3 = eq_label(tx(r"\xrightarrow{\ \epsilon\to0\ }\ \pi\sum_j\delta(\lambda-\lambda_j)", 34).set_color(YELLOW), "(6.42)")
    box = D.erratum_box("p. 146，(6.52)", r"\mathrm{Im}\,\mathcal{R}(\lambda-i\epsilon)=\sum_i\delta(\lambda-\lambda_i)",
                        r"\mathrm{Im}\,\mathcal{R}(\lambda-i\epsilon)=\pi\sum_i\delta(\lambda-\lambda_i)", width=6.6)
    left = VGroup(r1, r2, r3, box).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    plot = lorentz_plot()
    z0 = mixed("零階（只配相鄰的 ", "$\\eta$", "，(6.50)）：尖峰在 ", "$\\frac{1}{\\pi^2k^2}=%s$" % ",\\ ".join(
        f4(v) for v in nb["lead"]), "（除以 ", "$N+1$", "）", font_size=24, color=GREEN_B)
    place_in(left, -6.75, -0.1, -1.75, 3.1)
    place_in(plot, 0.3, 6.75, -1.75, 3.1)
    place_in(z0, -6.75, 6.75, -2.9, -2.1)
    assert_in_frame(left, plot, z0)
    assert_apart(left, plot)
    assert_apart(left, z0)
    assert_apart(plot, z0)
    assert_column(*left)
    s.play(FadeIn(r1), run_time=0.8)
    s.say("預解式：", "$\\mathbf{T}$", " 的特徵值變成 ", "$\\mathcal R$", " 的極點", wait=2.5)
    s.play(FadeIn(r2), ShowCreation(plot[0]), run_time=0.8)
    s.say("$\\lambda$", " 往下移一點點（", "$-i\\epsilon$", "）：虛部是一群尖峰", wait=1.5)
    s.play(*[ShowCreation(c) for c in plot[2]], FadeIn(plot[3]), run_time=1.5)
    s.play(FadeIn(plot[1]), run_time=0.6)
    s.say("$\\frac{\\epsilon}{x^2+\\epsilon^2}$", " 底下的面積是 ", "$\\pi=3.1416$", "，不管 ", "$\\epsilon$", " 多小",
          wait=2.5)
    s.play(FadeIn(r3), run_time=0.8)
    s.say("所以 ", "$\\frac1\\pi\\mathrm{Im}\\,\\mathcal R$", " 就是特徵值的分布", wait=2.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (6.52) 漏了 ", "$\\pi$", "：面積對不上", wait=2.5)
    s.play(FadeIn(z0), run_time=0.8)
    s.say("零階：尖峰就在樓梯的值 0.1013、0.0253、0.0113", wait=3.5)
    s.clear_stage()


def part_first_order(s):
    nb = numbers()
    o1 = eq_label(tx(r"\langle\Lambda_k\rangle=\frac{N+1}{\pi^2k^2}\Big(1+\frac{3}{4d}\Big)", 34), "(6.58)--(6.61)")
    o2 = eq_label(tx(r"\sum_{m\ne n}\frac{n^2}{m^2-n^2}=\frac34", 32), "(6.60)")
    o2b = mixed("數值：", "$n=1,2,3$", " 都是 ", "$%s$" % f4(nb["sumid"][0]), font_size=24, color=GREY_A)
    o3 = eq_label(tx(r"\mathrm{Var}\,\Lambda_k=\frac{2}{d}\Big(\frac{N+1}{\pi^2k^2}\Big)^2", 34), "(6.62)")
    left = VGroup(o1, o2, o2b, o3).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    rows = [["d", mixed("相對位移 ", "$\\frac{3}{4d}$", font_size=24, color=YELLOW),
             mixed("相對寬度 ", "$\\sqrt{2/d}$", font_size=24, color=YELLOW)]]
    rows += [[str(d), f4(sh), f4(w)] for d, sh, w in zip(DIMS, nb["shift"], nb["width"])]
    tab = D.tex_table(rows, [1.2, 2.6, 2.6], row_h=0.8, font_size=30)
    D.table_cell(tab, 0, 0, 3).set_color(YELLOW)
    place_in(left, -6.75, -0.2, -2.9, 3.1)
    place_in(tab, 0.3, 6.75, -2.0, 2.0)
    assert_in_frame(left, tab)
    assert_apart(left, tab)
    assert_column(*left)
    assert_table(tab, 3)
    s.play(FadeIn(o1), run_time=0.8)
    s.say("一階：多算一層配對，尖峰往右移", wait=2.0)
    s.play(FadeIn(o2), FadeIn(o2b), run_time=0.8)
    s.say("相對位移 ", "$\\frac{3}{4d}$", "：來自一個漂亮的恆等式，和是 ", "$\\frac34$", wait=2.5)
    s.play(FadeIn(o3), run_time=0.8)
    s.say("尖峰也變寬：標準差除以平均＝", "$\\sqrt{2/d}$", wait=2.5)
    s.play(FadeIn(tab[0][:3]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(tab[0][3 * i:3 * i + 3]) for i in range(1, 4)], lag_ratio=0.4), run_time=1.5)
    s.say("$d=100$", "：位移 0.0075、寬度 0.1414，幾乎還是尖峰", wait=2.5)
    s.say("$d=3$", "：位移 0.25、寬度 0.82，", "$1/d$", " 已經不小", wait=2.5)
    s.clear_stage()


def part_hist(s):
    panels, texts = [], []
    boxes = [(-6.75, -2.35), (-2.15, 2.15), (2.35, 6.75)]
    for d, (x0, x1) in zip(DIMS, boxes):
        p, t = hist_panel(d)
        place_in(p, x0, x1, -0.95, 3.0)
        place_in(t, x0, x1, -2.9, -1.2)
        panels.append(p)
        texts.append(t)
    assert_in_frame(*panels, *texts)
    for a, b in zip(panels, panels[1:]):
        assert_apart(a, b)
    for a, b in zip(texts, texts[1:]):
        assert_apart(a, b)
    for p, t in zip(panels, texts):
        assert_apart(p, t)
    s.say("模擬：高斯步伐、", "$N=50$", "，前三個特徵值除以 ", "$N+1$", wait=1.0)
    for p in panels:
        ax, polys, lines, title = p
        s.play(ShowCreation(ax), FadeIn(title), ShowCreation(lines), run_time=0.8)
        s.play(LaggedStart(*[FadeIn(q) for q in polys], lag_ratio=0.3), run_time=1.2)
    s.say("藍、橘、粉紅：第 1、2、3 大（每色以最高點為 1）", wait=2.5)
    s.say("虛線＝零階的尖峰：100 維幾乎貼著，3 維寬到疊在一起", wait=3.0)
    s.say("下面：模擬除以一階公式 (6.61)(6.62)，1 就是完全吻合", wait=1.0)
    for t in texts:
        s.play(FadeIn(t), run_time=0.6)
    s.wait(1.0)
    s.say("$d=100$", "：平均與變異數都在 0.94–1.01", wait=3.0)
    s.say("$d=20$", "：第 3 個的變異數只剩 0.68 倍", wait=3.0)
    s.say("$d=3$", "：第 2、3 個偏離（變異數 0.47、0.23）", wait=3.0)
    s.say("書上 ", "$d=100$", " 的變異數比：0.97、0.93、0.95", wait=2.5)
    s.say("書上 ", "$d=3$", "：變異數比 0.98、0.48、0.23", wait=2.5)
    s.say("平均比 1.01、0.91、0.75：都和模擬一致", wait=2.5)
    s.say("維度越低，尖峰越寬、越疊在一起", wait=2.5)
    s.say("最大的主半徑一直最準（平均比 1.01）", wait=3.0)


class W606OneOverD(ChuScene):
    def construct(self):
        self.show_header("W606　1/d 展開：尖峰變寬")
        check_math()
        part_question(self)
        part_wick(self)
        part_resolvent(self)
        part_first_order(self)
        part_hist(self)
        self.wait(1.0)
