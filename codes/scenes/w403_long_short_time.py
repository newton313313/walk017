"""
W403 短時間像沒牆、長時間指數衰減（§4.1.3，pp. 76–79，(4.16)–(4.27)，Fig. 4.2）

問：兩面牆的距離加倍，還活著的比例衰減變慢幾倍？
  1. (4.16)–(4.19)：∫ sin(nπx/L) dx 偶數 n 是 0、奇數 n 是 2L/(nπ) → 2^N (8L/π²) Σ e^{−N(2m+1)²π²/2L²}/(2m+1)²。
  2. 短時間 (4.20)–(4.23)：Σ_奇 1/n² = π²/8 = 1.2337（部分和長條）→ 2^N L；筆誤框 p. 77 (4.22)：k = 0 → k = 1。
  3. 長時間 (4.24)：對數圖（L = 10）：BC.box_survival 的點、三條直線（BC.box_survival_long）、斜率（BC.decay_rate）；
     筆誤框 p. 77 (4.24)：指數少了 ½。
  4. 揭曉：L = 20 衰減率 0.01239 → 慢 4.05 倍（連續極限 4 倍）；每 2L²/π² = 20.26 步少 e 倍。
  5. Fig. 4.2：每步分裂 ×2、牆邊漏光 → n_w ≈ 2^N e^{−ακN/L²}，ακ = π²/2 = 4.935。
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

LW = 10
DEAD = D.WRONG_COLOR
KIND_COLOR = {"lattice": YELLOW, "cont": D.EXACT_COLOR, "book": DEAD}
L20_COLOR = TEAL
NPTS = list(range(0, 201, 10))
NB = 8                                   # 長條數


def sci(v, d=3):
    """3.881\\times10^{-5} 形式。"""
    m, e = f"{v:.{d}e}".split("e")
    return m + r"\times10^{" + str(int(e)) + "}"


def partial_sums(M):
    return [math.fsum(1 / (2 * m + 1) ** 2 for m in range(k + 1)) for k in range(M)]


@lru_cache(maxsize=None)
def numbers():
    nb = {"ps": partial_sums(NB), "lim": math.pi ** 2 / 8}
    nb["surv"] = {N: float(BC.box_survival(N, LW)) for N in NPTS}
    nb["long100"] = {k: BC.box_survival_long(100, LW, k) for k in KIND_COLOR}
    nb["rate"] = {k: BC.decay_rate(LW, k) for k in KIND_COLOR}
    nb["r20"] = BC.decay_rate(20, "lattice")
    nb["slow"] = nb["rate"]["lattice"] / nb["r20"]
    nb["efold"] = 2 * LW ** 2 / math.pi ** 2
    nb["ak"] = math.pi ** 2 / 2
    return nb


def check_math():
    nb = numbers()
    assert [f"{v:.4f}" for v in nb["ps"][:4]] == ["1.0000", "1.1111", "1.1511", "1.1715"]
    assert f"{nb['lim']:.4f}" == "1.2337" and abs(BC.odd_inv_squares(60) - nb["lim"]) < 1e-9
    assert all(a < b < nb["lim"] for a, b in zip(nb["ps"], nb["ps"][1:]))
    assert abs(BC.zeta2(200000) - math.pi ** 2 / 6) < 1e-5 and abs(0.75 * math.pi ** 2 / 6 - nb["lim"]) < 1e-15
    for n in (1, 2, 3, 4, 5):                                   # (4.16)–(4.18)：∫_0^L sin(nπx/L) dx
        xs = [(i + 0.5) * LW / 4000 for i in range(4000)]
        integ = math.fsum(math.sin(n * math.pi * x / LW) for x in xs) * LW / 4000
        assert abs(integ - (2 * LW / (n * math.pi) if n % 2 else 0)) < 1e-5
    assert abs(2 / LW * (2 * LW / math.pi) ** 2 - 8 * LW / math.pi ** 2) < 1e-12
    s = nb["surv"]
    assert s[0] == 1 and f"{s[50]:.6f}" == "0.072103" and f"{s[100]:.6f}" == "0.005865"
    assert sci(s[200]) == r"3.881\times10^{-5}"
    assert all(s[a] > s[b] for a, b in zip(NPTS, NPTS[1:])) and math.log10(s[200]) > -5
    lg = nb["long100"]
    assert f"{lg['lattice']:.6f}" == "0.005861" and f"{lg['cont']:.6f}" == "0.005830"
    assert sci(lg["book"], 2) == r"4.19\times10^{-5}" and f"{8 / math.pi ** 2:.4f}" == "0.8106"
    assert abs(BC.box_survival_long(0, LW, "cont") - 8 / math.pi ** 2) < 1e-15
    r = nb["rate"]
    assert f"{r['lattice']:.5f}" == "0.05018" and abs(r["lattice"] + math.log(math.cos(math.pi / LW))) < 1e-15
    assert f"{r['cont']:.5f}" == "0.04935" and abs(r["cont"] - math.pi ** 2 / 200) < 1e-15
    assert f"{r['book']:.4f}" == "0.0987" and abs(r["book"] - math.pi ** 2 / 100) < 1e-15
    assert f"{math.exp(50 * r['lattice']):.0f}" == "12"
    assert f"{nb['r20']:.5f}" == "0.01239" and f"{nb['slow']:.2f}" == "4.05"
    assert abs(BC.decay_rate(10, "cont") / BC.decay_rate(20, "cont") - 4) < 1e-12
    assert f"{nb['efold']:.2f}" == "20.26" and f"{nb['ak']:.3f}" == "4.935"
    assert BC.box_survival_long(114, LW, "book") > 1e-5 > BC.box_survival_long(115, LW, "book")


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


def segment(length, unit, y, label):
    """0..length 的走廊（紅牆）＋下方大括號標籤。"""
    x0 = -length * unit / 2
    base = Line(np.array([x0, y, 0]), np.array([-x0, y, 0])).set_stroke(GREY_B, 3)
    walls = VGroup(*[Line(np.array([x, y - 0.3, 0]), np.array([x, y + 0.3, 0])).set_stroke(DEAD, 6) for x in (x0, -x0)])
    br = Brace(base, DOWN, buff=0.15).set_color(GREY_A)
    lab = tx(label, 30, GREY_A).next_to(br, DOWN, buff=0.1)
    dot = D.walker_dot().move_to(np.array([0, y + 0.3, 0]))
    return VGroup(base, walls, br, lab, dot)


def log_plot():
    """L = 10 的存活比例對數圖：縱軸畫 log10(值) + 5。"""
    ax = D.plot_axes((0, 210, 50), (0, 5.4, 1), width=6.2, height=4.9, x_ticks=[0, 50, 100, 150, 200],
                     y_ticks=[0, 1, 2, 3, 4, 5], y_tick_labels=[r"10^{-5}", r"10^{-4}", r"10^{-3}", r"10^{-2}",
                                                                r"10^{-1}", "1"], font_size=20, x_label="N")
    ax.move_to(np.array([-3.35, -0.05, 0]))
    return ax


def log_line(ax, fn, n0, n1, color, dashed=False):
    pts = [ax.c2p(n, math.log10(fn(n)) + 5) for n in np.linspace(n0, n1, 60)]
    c = VMobject().set_points_as_corners(pts).set_stroke(color, 3)
    return DashedVMobject(c, num_dashes=30).set_stroke(color, 3) if dashed else c


# ---------------------------------------------------------------- 各段
def part_question(s):
    a = segment(10, 0.4, 1.6, "L")
    b = segment(20, 0.4, -0.6, "2L")
    q = at(mixed("衰減變慢幾倍？", "$2$", " 倍？", "$4$", " 倍？", font_size=32, color=YELLOW), 0, -2.4)
    assert_in_frame(a, b, q)
    assert_apart(a, b)
    assert_apart(b, q)
    s.play(FadeIn(a), run_time=0.8)
    s.say("走廊兩端都是懸崖：活著的比例一路往下掉", wait=3.0)
    s.play(FadeIn(b), run_time=0.8)
    s.say("問：兩面牆的距離加倍，衰減變慢幾倍？", wait=2.5)
    s.play(FadeIn(q), run_time=0.6)
    s.say("2 倍？4 倍？先猜", wait=2.5)
    s.clear_stage()


def part_total(s):
    e1 = eq_label(tx(r"\int_0^L\sin\frac{n\pi x}{L}\,dx=\frac{L}{n\pi}\big(1-\cos n\pi\big)"
                     r"=\begin{cases}0, & n\ \text{even}\\[2pt] \dfrac{2L}{n\pi}, & n\ \text{odd}\end{cases}", 32),
                  "(4.16)(4.17)")
    at(e1, 0, 2.2)
    e2 = at(tx(r"\frac2L\Big(\frac{2L}{n\pi}\Big)^2=\frac{8L}{\pi^2n^2}", 32), 0, 0.65)
    e3 = eq_label(VGroup(zh("存活總數", font_size=30), tx(r"\approx2^N\,\frac{8L}{\pi^2}\sum_{m\ge0}"
                                                         r"\frac{e^{-N(2m+1)^2\pi^2/2L^2}}{(2m+1)^2}", 34))
                  .arrange(RIGHT, buff=0.2), "(4.19)")
    at(e3, 0, -1.1)
    box = SurroundingRectangle(e3[0], buff=0.15).set_stroke(YELLOW, 2)
    assert_in_frame(e1, e2, e3, box)
    assert_apart(e1, e2)
    assert_apart(e2, box)
    s.play(FadeIn(e1), run_time=1.2)
    s.say("從每個起點出發、終點不限：活著的總數", wait=2.5)
    s.say("起點與終點都對正弦積分：偶數模態全部抵銷", wait=3.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("奇數模態：兩個積分各給一個 ", "$2L/n\\pi$", wait=3.0)
    s.play(FadeIn(e3, shift=0.2 * UP), ShowCreation(box), run_time=1.2)
    s.say("連續版：每個模態 ", "$e^{-Nq^2/2}$", "，", "$q=n\\pi/L$", wait=3.5)
    s.clear_stage()


def part_short(s):
    nb = numbers()
    e3 = eq_label(tx(r"N\ll L^2:\ \ 2^N\frac{8L}{\pi^2}\sum_{n\ \text{odd}}\frac1{n^2}"
                     r"=2^N\frac{8L}{\pi^2}\cdot\frac{\pi^2}{8}=2^NL", 32), "(4.20)(4.23)")
    at(e3, 0, 2.45)
    z = eq_label(tx(r"\sum_{n\ \text{odd}}\frac1{n^2}=\sum_{k=1}^{\infty}\frac1{k^2}-\frac14\sum_{k=1}^{\infty}\frac1{k^2}"
                    r"=\frac34\cdot\frac{\pi^2}6=\frac{\pi^2}8", 30), "(4.21)(4.22)")
    at(z, 0, 1.1)
    ux, uy = 0.72, 2.0
    org = np.array([-6.2, -2.45, 0])
    vals = {k + 1: v for k, v in enumerate(nb["ps"])}
    bars = D.bar_chart(vals, ux=ux, uy=uy, origin=org, width=0.82)
    labs = VGroup(*[tx(f"{nb['ps'][k]:.4f}" if k else "1", 17, WHITE).move_to(bars[0][k].get_top() + 0.2 * DOWN)
                    for k in range(4)])
    yl = org[1] + nb["lim"] * uy
    lim = DashedLine(np.array([org[0] + 0.3, yl, 0]), np.array([org[0] + (NB + 0.6) * ux, yl, 0])).set_stroke(YELLOW, 2)
    llab = tx(rf"\pi^2/8={nb['lim']:.4f}", 24, YELLOW).next_to(lim, UP, buff=0.08).align_to(lim, RIGHT)
    cap = mixed("奇數項 ", "$1/n^2$", " 的部分和", font_size=22, color=GREY_A).next_to(bars, DOWN, buff=0.08)
    box = D.erratum_box("p. 77，(4.22)", r"\frac14\sum_{k=0}^{\infty}\frac{1}{k^2}",
                        r"\frac14\sum_{k=1}^{\infty}\frac{1}{k^2}", width=5.4)
    at(box, 3.65, -1.35)
    assert_in_frame(e3, z, bars, labs, lim, llab, cap, box)
    assert_apart(e3, z)
    assert_apart(z, llab)
    assert_apart(z, bars)
    assert_apart(bars, box)
    assert_apart(llab, box)
    assert_apart(z, box)
    s.play(FadeIn(e3), run_time=1.2)
    s.say("短時間：", "$N\\ll L^2$", "，每個指數都差不多是 1", wait=3.0)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars[0]], lag_ratio=0.15), FadeIn(cap), run_time=2.0)
    s.play(FadeIn(labs), ShowCreation(lim), FadeIn(llab), run_time=1.0)
    s.say("$1+\\frac19+\\frac1{25}+\\cdots$", " 慢慢爬向 ", "$\\pi^2/8=1.2337$", wait=3.0)
    s.play(FadeIn(z), run_time=1.0)
    s.say("全部的和減掉偶數項：偶數項是全部的 1/4", wait=3.0)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上從 ", "$k=0$", " 加起：那一項是無窮大", wait=3.0)
    s.say("結果 ", "$2^NL$", "：和沒有牆一樣——短時間牆不重要", wait=3.0)
    s.clear_stage()


def part_long(s):
    nb = numbers()
    ax = log_plot()
    pts = VGroup(*[Dot(ax.c2p(N, math.log10(nb["surv"][N]) + 5), radius=0.06).set_fill(WHITE, 1) for N in NPTS])
    lines = {k: log_line(ax, lambda n, k=k: BC.box_survival_long(n, LW, k), 0, 200 if k != "book" else 114,
                         KIND_COLOR[k], dashed=(k == "book")) for k in KIND_COLOR}
    hdr = mixed("存活比例，", "$L=10$", font_size=24, color=GREY_A).next_to(ax, UP, buff=0.05).align_to(ax, LEFT)
    hdr.shift(0.9 * RIGHT)
    s_ = nb["surv"]
    t1 = D.tex_table([["N", zh("存活比例", font_size=22, color=YELLOW)], ["50", f"{s_[50]:.6f}"],
                      ["100", f"{s_[100]:.6f}"], ["200", sci(s_[200])]], [1.0, 2.6], row_h=0.42, font_size=24)
    t1.move_to(np.array([3.6, 3.05 - t1.get_height() / 2, 0]))
    lg, r = nb["long100"], nb["rate"]
    rows = [[zh("直線", font_size=20, color=YELLOW), zh("每步衰減率", font_size=20, color=YELLOW), "N=100"],
            [zh("格子最低模態", font_size=20, color=KIND_COLOR["lattice"]),
             rf"-\ln\cos\frac{{\pi}}{{10}}={r['lattice']:.5f}", f"{lg['lattice']:.6f}"],
            [zh("更正 (4.24)", font_size=20, color=KIND_COLOR["cont"]), rf"\frac{{\pi^2}}{{200}}={r['cont']:.5f}",
             f"{lg['cont']:.6f}"],
            [zh("書上 (4.24)", font_size=20, color=KIND_COLOR["book"]), rf"\frac{{\pi^2}}{{100}}={r['book']:.4f}",
             sci(lg["book"], 2)]]
    t2 = D.tex_table(rows, [1.7, 2.75, 1.75], row_h=0.58, font_size=22)
    t2.move_to(np.array([3.55, 1.0 - t2.get_height() / 2, 0]))
    pref = mixed("更正：", "$\\frac{8}{\\pi^2}e^{-\\pi^2N/2L^2}$", "，", "$\\frac8{\\pi^2}=0.8106$", font_size=22,
                 color=KIND_COLOR["cont"])
    pref.next_to(t2, DOWN, buff=0.18)
    box = D.erratum_box("p. 77，(4.24)", r"e^{-\pi^2\Delta^2N/L^2}", r"e^{-\pi^2\Delta^2N/2L^2}", width=5.4)
    box.move_to(np.array([3.6, 3.05 - box.get_height() / 2, 0]))
    assert_in_frame(ax, pts, hdr, t1, t2, pref, box, *lines.values())
    assert_table(t1, 2)
    assert_table(t2, 3)
    for m in (t1, t2, pref, box):
        assert_apart(ax, m)
    assert_apart(t1, t2)
    assert_apart(t2, pref)
    assert_apart(box, t2)
    s.play(FadeIn(ax), FadeIn(hdr), run_time=1.0)
    s.play(LaggedStart(*[FadeIn(p) for p in pts], lag_ratio=0.05), run_time=1.6)
    s.say("長時間：", "$L=10$", " 的精確存活比例（縱軸取對數）", wait=2.5)
    s.play(FadeIn(t1[0]), ShowCreation(t1[1]), run_time=1.0)
    s.say("每多走 50 步，大約少 12 倍：對數圖上是直線", wait=3.0)
    s.play(ShowCreation(lines["lattice"]), FadeIn(VGroup(*t2[0][:3], *t2[0][3:6])), ShowCreation(t2[1]), run_time=1.2)
    s.say("(4.24)：只剩最低模態 ", "$m=1$", "，其他早就衰減光了", wait=3.0)
    s.play(ShowCreation(lines["cont"]), FadeIn(VGroup(*t2[0][6:9])), FadeIn(pref), run_time=1.2)
    s.say("連續近似的斜率 ", "$\\pi^2/2L^2$", "：和格子只差一點", wait=3.0)
    s.play(ShowCreation(lines["book"]), FadeIn(VGroup(*t2[0][9:12])), run_time=1.2)
    s.say("照書上的指數：斜率大一倍，完全對不上", wait=3.0)
    s.play(FadeOut(t1), FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("(4.19) 的指數有 ", "$\\frac12$", "，(4.24) 把它掉了", wait=3.0)
    s.play(FadeOut(VGroup(t2, pref, box, lines["book"])), run_time=0.8)

    l20 = log_line(ax, lambda n: BC.box_survival_long(n, 20, "lattice"), 0, 200, L20_COLOR)
    tag = tx("L=20", 24, L20_COLOR).next_to(l20.get_end(), UP, buff=0.12).shift(0.2 * LEFT)
    r1 = mixed("$L=10$", "：每步 ", f"${r['lattice']:.5f}$", font_size=26, color=YELLOW)
    r2 = mixed("$L=20$", "：每步 ", f"${nb['r20']:.5f}$", font_size=26, color=L20_COLOR)
    r3 = mixed("慢 ", f"${nb['slow']:.2f}$", " 倍（連續極限正好 4 倍）", font_size=26, color=ORANGE)
    r4 = tx(r"\text{rate}\approx\frac{\pi^2}{2L^2}", 32, ORANGE)
    r5 = mixed("每多走 ", "$2L^2/\\pi^2$", " 步少 ", "$e$", " 倍", font_size=26)
    r6 = mixed("$L=10$", "：", f"${nb['efold']:.2f}$", " 步", font_size=26)
    rev = VGroup(r1, r2, r3, r4, r5, r6).arrange(DOWN, aligned_edge=LEFT, buff=0.32)
    rev.move_to(np.array([3.65, 0.2, 0]))
    assert_in_frame(rev, l20, tag)
    assert_apart(ax, rev)
    s.play(ShowCreation(l20), FadeIn(tag), run_time=1.2)
    s.say("揭曉：牆距加倍成 ", "$L=20$", wait=2.0)
    s.play(FadeIn(r1), FadeIn(r2), run_time=0.8)
    s.play(FadeIn(r3), FadeIn(r4), run_time=0.8)
    s.say("衰減慢 4 倍：衰減率跟 ", "$L^2$", " 成反比", wait=3.0)
    s.play(FadeIn(r5), FadeIn(r6), run_time=0.8)
    s.say("$L=10$", "：每 20 步左右，活著的少 ", "$e$", " 倍", wait=3.0)
    s.clear_stage()


def part_leak(s):
    nb = numbers()
    unit, amp = 0.55, 1.6
    o = np.array([-5.9, -1.9, 0])
    base = Line(o, o + LW * unit * RIGHT).set_stroke(GREY_B, 2)
    walls = VGroup(*[Line(o + w * unit * RIGHT, o + w * unit * RIGHT + 2.0 * UP).set_stroke(DEAD, 5) for w in (0, LW)])
    tent = D.curve(lambda x: math.sin(math.pi * x / LW), 0, LW, ux=unit, uy=amp, origin=o, color=YELLOW, stroke_width=4)
    rng = np.random.default_rng(4)
    xs = rng.uniform(0.6, LW - 0.6, 40)
    ys = [rng.uniform(0.1, 0.92) * math.sin(math.pi * x / LW) for x in xs]
    crowd = VGroup(*[Dot(o + x * unit * RIGHT + y * amp * UP, radius=0.045).set_fill(BLUE_B, 0.9) for x, y in zip(xs, ys)])
    leaks = VGroup(Arrow(o + 0.5 * unit * RIGHT + 0.35 * UP, o + 0.7 * LEFT + 0.35 * UP, buff=0, thickness=3),
                   Arrow(o + (LW - 0.5) * unit * RIGHT + 0.35 * UP, o + LW * unit * RIGHT + 0.7 * RIGHT + 0.35 * UP,
                         buff=0, thickness=3)).set_color(DEAD)
    ltxt = zh("牆邊漏掉", font_size=22, color=DEAD).next_to(leaks[1].get_end(), RIGHT, buff=0.12)
    fig = zh("Fig. 4.2：中間放人、兩邊漏光", font_size=22, color=GREY_A).next_to(base, DOWN, buff=0.2)
    l1 = mixed("每步：分裂 ", "$\\times2$", "；漏掉的比例 ", "$\\propto$", " 斜率 ", "$\\propto n_w/L^2$", font_size=26)
    l2 = eq_label(tx(r"n_w(N)=2^N\Big(1-\frac{\alpha\kappa}{L^2}\Big)^N\approx2^N\,e^{-\alpha\kappa N/L^2}", 32),
                  "(4.25)(4.27)")
    l3 = mixed("精確計算：", f"$\\alpha\\kappa=\\pi^2/2={nb['ak']:.3f}$", font_size=28, color=D.EXACT_COLOR)
    rt = VGroup(l1, l2, l3).arrange(DOWN, buff=0.4)
    rt.move_to(np.array([0, 1.85, 0]))
    assert_in_frame(base, walls, tent, crowd, leaks, ltxt, fig, rt)
    assert_apart(rt, walls)
    assert_apart(rt, ltxt)
    assert_apart(ltxt, walls)
    assert_apart(rt, tent)
    s.play(FadeIn(base), FadeIn(walls), ShowCreation(tent), run_time=1.2)
    s.play(LaggedStart(*[FadeIn(d) for d in crowd], lag_ratio=0.03), FadeIn(fig), run_time=1.2)
    s.say("另一種看法：人群的分布像帳篷，牆上是 0", wait=3.0)
    s.play(GrowArrow(leaks[0]), GrowArrow(leaks[1]), FadeIn(ltxt), run_time=1.0)
    s.play(FadeIn(l1), run_time=0.8)
    s.say("帳篷寬 ", "$L$", "、高 ", "$\\propto n_w/L$", "：牆邊斜率 ", "$\\propto n_w/L^2$", wait=3.5)
    s.play(FadeIn(l2), run_time=1.0)
    s.say("每步都漏掉固定比例 ", "$\\alpha\\kappa/L^2$", "：指數衰減", wait=3.0)
    s.play(FadeIn(l3), run_time=0.8)
    s.say("常數要精確計算：", "$\\alpha\\kappa=\\pi^2/2$", wait=2.5)
    s.clear_stage()


def part_close(s):
    s.say("牆的影響要等 ", "$N\\sim L^2$", " 步才出現", wait=2.5)
    s.say("之後活著的比例像 ", "$e^{-\\pi^2N/2L^2}$", " 那樣消失", wait=3.0)


class W403LongShortTime(ChuScene):
    def construct(self):
        self.show_header("W403　短時間像沒牆、長時間指數衰減")
        check_math()
        part_question(self)
        part_total(self)
        part_short(self)
        part_long(self)
        part_leak(self)
        part_close(self)
        self.wait(1.0)
