"""
W306 三種立方格子的逃脫機率（§3.7，pp. 62–63，(3.34)–(3.41)，Ex. 3.7–3.8）

問：鄰居越多，越容易逃走嗎？（sc 6 個、bcc 8 個、fcc 12 個最近鄰）
  1. 三種格子並排（L.lattice3d）＋各自的 χ（整數座標的鄰居；check_math 用鄰居直接加總核對）；
     (3.35) G(z_c;0) = (2π)^{−3}∫d³q/(1 − χ/χ(0))（Watson 積分）。
  2. 先猜再揭曉：G = 1.516386、1.393204、1.344661 → 逃脫機率 0.659463、0.717770、0.743682（R.G_C、R.escape）；
     Γ 閉式（本檔 closed_form 用 math.gamma 重算）；第二種算法 R.watson_sc_bessel()、R.watson_bcc_series()。
  3. 練習 3.8：R = 0.340537，f_m = R^m(1 − R)，Σf_m = 1，平均 G − 1 = 0.516386 次；有限 N = 1000 的精確值
     （R.return_count_dist）；長條：理論 1000·f_m 對模擬 R.SIM_COUNTS_SC（重算 assert）。
     最後一根是「≥ 4 次」：理論用 Σ_{m≥4} f_m = R⁴ = 0.0134。
畫面上的數字都由 utils.walk_rec 算出，check_math() 核對。
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
from utils import walk_rec as R
from utils import walk_gf as GF
from utils import lattice_draw as L
from utils import walk_draw as D

KINDS = ("sc", "bcc", "fcc")
COLORS = {"sc": TEAL, "bcc": "#F2A33A", "fcc": "#C792EA"}
XS = (-4.5, 0.0, 4.5)
NAMES = {"sc": "簡單立方 sc", "bcc": "體心立方 bcc", "fcc": "面心立方 fcc"}
STEPS = {"sc": GF.SC_STEPS, "bcc": GF.BCC_STEPS, "fcc": GF.FCC_STEPS}
CHI = {"sc": [r"\chi=2(\cos q_x+\cos q_y+\cos q_z)"], "bcc": [r"\chi=8\cos q_x\cos q_y\cos q_z"],
       "fcc": [r"\chi=4(\cos q_x\cos q_y+\cos q_y\cos q_z", r"+\cos q_x\cos q_z)"]}
CLOSED = {"sc": r"\frac{\sqrt6}{32\pi^3}\Gamma(\tfrac{1}{24})\Gamma(\tfrac{5}{24})\Gamma(\tfrac{7}{24})\Gamma(\tfrac{11}{24})",
          "bcc": r"\frac{\Gamma(1/4)^4}{4\pi^3}", "fcc": r"\frac{9\,\Gamma(1/3)^6}{2^{14/3}\pi^4}"}
MS = range(5)


def chi_closed(kind, q):
    x, y, z = q
    c = math.cos
    if kind == "sc":
        return 2 * (c(x) + c(y) + c(z))
    if kind == "bcc":
        return 8 * c(x) * c(y) * c(z)
    return 4 * (c(x) * c(y) + c(y) * c(z) + c(x) * c(z))


def closed_form(kind):
    g, pi = math.gamma, math.pi
    if kind == "sc":
        return math.sqrt(6) / (32 * pi ** 3) * g(1 / 24) * g(5 / 24) * g(7 / 24) * g(11 / 24)
    if kind == "bcc":
        return g(0.25) ** 4 / (4 * pi ** 3)
    return 9 * g(1 / 3) ** 6 / (2 ** (14 / 3) * pi ** 4)


@lru_cache(maxsize=None)
def numbers():
    nb = {"G": {k: f"{R.G_C[k]:.6f}" for k in KINDS}, "esc": {k: f"{R.escape(k):.6f}" for k in KINDS}}
    nb["esc3"] = [f"{R.escape(k):.3f}" for k in KINDS]
    Rr = 1 - R.escape("sc")
    nb["R"] = f"{Rr:.6f}"
    nb["fm"] = [f"{R.f_m('sc', m):.4f}" for m in MS]
    nb["fin"] = [f"{v:.4f}" for v in R.return_count_dist("sc", 1000, 4)]
    nb["mean"] = f"{R.G_C['sc'] - 1:.6f}"
    nb["R4"] = f"{Rr ** 4:.4f}"
    nb["theory"] = [1000 * R.f_m("sc", m) for m in range(4)] + [1000 * Rr ** 4]
    return nb


def check_math():
    nb = numbers()
    assert nb["G"] == {"sc": "1.516386", "bcc": "1.393204", "fcc": "1.344661"}
    assert nb["esc"] == {"sc": "0.659463", "bcc": "0.717770", "fcc": "0.743682"}
    assert nb["esc3"] == ["0.659", "0.718", "0.744"]
    q = (0.37, -1.1, 2.3)
    for k in KINDS:
        assert len(STEPS[k]) == R.Z[k] == round(chi_closed(k, (0, 0, 0))) == R.returns_count(k, 2)
        direct = sum(math.cos(sum(a * b for a, b in zip(q, w))) for w in STEPS[k])
        assert abs(direct - chi_closed(k, q)) < 1e-12
        assert abs(closed_form(k) - R.G_C[k]) < 1e-12 and abs(R.escape(k) * R.G_C[k] - 1) < 1e-15
    assert abs(R.watson_sc_bessel() - R.G_C["sc"]) < 1e-10 and abs(R.watson_bcc_series() - R.G_C["bcc"]) < 1e-10
    assert R.Z["sc"] < R.Z["bcc"] < R.Z["fcc"] and R.escape("sc") < R.escape("bcc") < R.escape("fcc")
    assert nb["R"] == "0.340537" and nb["mean"] == "0.516386" and nb["R4"] == "0.0134"
    assert nb["fm"] == ["0.6595", "0.2246", "0.0765", "0.0260", "0.0089"]
    assert nb["fin"] == ["0.6685", "0.2217", "0.0735", "0.0243", "0.0080"]
    Rr = 1 - R.escape("sc")
    assert abs(sum(R.f_m("sc", m) for m in range(200)) - 1) < 1e-12
    assert abs(sum(m * R.f_m("sc", m) for m in range(400)) - (R.G_C["sc"] - 1)) < 1e-12
    assert abs(Rr / (1 - Rr) - (R.G_C["sc"] - 1)) < 1e-12
    assert abs(sum(R.f_m("sc", m) for m in range(4, 200)) - Rr ** 4) < 1e-12
    fin = R.return_count_dist("sc", 1000, 4)
    assert fin[0] > R.f_m("sc", 0) and all(fin[m] < R.f_m("sc", m) for m in range(1, 5))
    counts = R.simulate_return_counts("sc", 1000, 1000, R.SIM_COUNTS_SEED)
    bins = [sum(1 for c in counts if c == m) for m in range(4)] + [sum(1 for c in counts if c >= 4)]
    assert len(counts) == 1000 and bins == R.SIM_COUNTS_SC == [672, 222, 60, 28, 18]
    assert f"{nb['theory'][0] / 1000:.4f}" == "0.6595" and f"{R.SIM_COUNTS_SC[0] / 1000:.3f}" == "0.672"


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


def assert_row(grp):
    for a, b in zip(grp, grp[1:]):
        assert_apart(a, b)


def assert_table(tab, ncols, gap=0.04):
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


# ---------------------------------------------------------------- 各段
def part_lattices(s):
    titles, lats, chis = VGroup(), VGroup(), VGroup()
    for k, x in zip(KINDS, XS):
        titles.add(at(fit(zh(f"{NAMES[k]}：{R.Z[k]} 個鄰居", font_size=24), 4.2), x, 2.85))
        lats.add(L.lattice3d(k, unit=0.8, origin=np.array([x, 1.35, 0]), color=COLORS[k]))
        g = VGroup(*[tx(t, 28) for t in CHI[k]]).arrange(DOWN, buff=0.12, aligned_edge=RIGHT)
        chis.add(at(fit(g, 4.3), x, -0.5))
    watson = VGroup(eq_label(tx(r"G(z_c;\vec 0)=\frac{1}{(2\pi)^3}\int\frac{d^3q}{1-\chi(\vec q)/\chi(\vec 0)}", 34),
                             "(3.35)"), zh("Watson 積分", font_size=26, color=GREY_A)).arrange(RIGHT, buff=0.5)
    at(watson, 0, -2.0)
    assert_in_frame(titles, lats, chis, watson)
    assert_row(titles)
    assert_row(chis)
    for t, m, c in zip(titles, lats, chis):
        assert_apart(t, m)
        assert_apart(m, c)
        assert_apart(c, watson)
    s.play(*[FadeIn(t) for t in titles], *[ShowCreation(m[0]) for m in lats], run_time=1.2)
    s.play(*[FadeIn(m[3]) for m in lats], run_time=0.5)
    s.say("問：鄰居越多，越容易逃走嗎？", wait=0.5)
    for m in lats:
        s.play(ShowCreation(m[1]), FadeIn(m[2]), run_time=0.8)
    s.say("簡單立方 6 個、體心立方 8 個、面心立方 12 個最近鄰", wait=3.5)
    s.play(LaggedStart(*[FadeIn(c) for c in chis], lag_ratio=0.3), run_time=1.2)
    s.say("鄰居的位置決定 ", "$\\chi(\\vec q)=\\sum_{\\vec w}e^{i\\vec q\\cdot\\vec w}$", "（鄰居座標取整數）", wait=3.5)
    s.say("$\\vec q=\\vec 0$", "：", "$\\chi(\\vec 0)=Z$", "，正好是鄰居的個數", wait=3.0)
    s.say("第 2 步就回家的機率是 ", "$1/Z$", "：鄰居越多越小", wait=3.5)
    s.play(FadeIn(watson, shift=0.2 * UP), run_time=1.0)
    s.say("W305：逃脫機率 ", "$=1/G(z_c;\\vec 0)$", "，要算這個三維積分", wait=3.5)
    s.say("三維的被積函數在 ", "$\\vec q=\\vec 0$", " 附近 ", "$\\propto 1/q^2$", "，積分有限", wait=3.5)
    s.play(FadeOut(chis), FadeOut(watson), run_time=0.7)
    return lats


def part_reveal(s, lats):
    nb = numbers()
    qs = VGroup(*[at(tx("G=?", 34), x, -0.45) for x in XS])
    gs = VGroup(*[at(tx("G=" + nb["G"][k], 34), x, -0.45) for k, x in zip(KINDS, XS)])
    es = VGroup(*[at(tx("1/G=" + nb["esc"][k], 34, COLORS[k]), x, -1.25) for k, x in zip(KINDS, XS)])
    cl = VGroup(*[at(fit(tx(CLOSED[k], 32), 4.3), x, -2.35) for k, x in zip(KINDS, XS)])
    assert_in_frame(qs, gs, es, cl)
    for grp in (gs, es, cl):
        assert_row(grp)
    for m, g, e, c in zip(lats, gs, es, cl):
        assert_apart(m, g)
        assert_apart(g, e)
        assert_apart(e, c)
    s.play(LaggedStart(*[FadeIn(q) for q in qs], lag_ratio=0.3), run_time=1.0)
    s.say("猜猜看：哪一種格子最容易逃走？", wait=3.5)
    s.play(*[FadeOut(q) for q in qs], *[FadeIn(g) for g in gs], run_time=1.0)
    s.say("Watson 積分的值：鄰居越多，", "$G$", " 越小", wait=3.5)
    s.play(LaggedStart(*[FadeIn(e, shift=0.2 * DOWN) for e in es], lag_ratio=0.3), run_time=1.2)
    s.say("逃脫機率 ", "$1/G$", "：" + "、".join(nb["esc3"]), wait=3.5)
    s.play(Indicate(es[2], scale_factor=1.1), run_time=1.0)
    s.say("面心立方最容易逃走：每 4 個人就有 3 個一去不返", wait=3.5)
    s.play(LaggedStart(*[FadeIn(c) for c in cl], lag_ratio=0.3), run_time=1.5)
    s.say("三個積分都有閉式：全是 ", "$\\Gamma$", " 函數", wait=3.5)
    s.say("第二種算法：", "$\\mathrm{sc}=\\int_0^\\infty e^{-t}I_0(t/3)^3dt$", wait=3.5)
    s.say("$\\mathrm{bcc}=\\sum_k\\big(\\binom{2k}{k}/4^k\\big)^3$", "：都和閉式一致到小數第 10 位", wait=3.5)
    s.clear_stage()


def part_counts(s):
    nb = numbers()
    e1 = at(tx(r"R=1-\frac{1}{G}=1-" + nb["esc"]["sc"] + "=" + nb["R"], 34), 0, 2.55)
    e2 = at(mixed("$f_m=R^m(1-R)$", "：剛好回家 ", "$m$", " 次", font_size=30), 0, 1.5)
    e3 = at(tx(r"\sum_{m\ge0}f_m=1,\qquad\sum_{m\ge0}m\,f_m=\frac{R}{1-R}=G-1=" + nb["mean"], 32), 0, 0.45)
    head = [zh("回家次數 m", font_size=26, color=YELLOW)] + [str(m) for m in MS]
    rows = [head, [zh("理論（無窮多步）", font_size=24)] + nb["fm"], [zh("精確（N＝1000）", font_size=24)] + nb["fin"]]
    tab = D.tex_table(rows, [3.4] + [1.5] * 5, row_h=0.62, font_size=30)
    tab.move_to(np.array([0, -0.5 - tab.get_height() / 2, 0]))
    assert_in_frame(e1, e2, e3, tab)
    assert_table(tab, 6)
    for a, b in ((e1, e2), (e2, e3), (e3, tab)):
        assert_apart(a, b)
    s.say("練習 3.8：在簡單立方上，每個人一共回家幾次？", wait=0.5)
    s.play(FadeIn(e1), run_time=0.8)
    s.say("回家的機率 ", "$R=1-1/G$", "，約三分之一", wait=3.5)
    s.say("每次回到家，就像重新出發：再回來的機率還是 ", "$R$", wait=3.5)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("回家 ", "$m$", " 次之後一去不返：", "$f_m=R^m(1-R)$", wait=3.5)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("全部加起來是 1；平均回家 ", "$G-1\\approx0.516$", " 次", wait=3.5)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.2)
    s.say("上一列：走無窮多步；下一列：只走 1000 步的精確值", wait=3.5)
    s.say("1000 步還沒走完：0 次略多，其他略少", wait=3.5)
    s.clear_stage()


def part_bars(s):
    nb = numbers()
    ux, uy, x0, y0 = 1.8, 4.0 / 680, -3.6, -2.2
    org = np.array([x0, y0, 0])
    th = D.bar_chart({m - 0.22: v for m, v in zip(MS, nb["theory"])}, ux=ux, uy=uy, origin=org,
                     color=D.BAR_COLOR, width=0.4)
    sim = D.bar_chart({m + 0.22: v for m, v in zip(MS, R.SIM_COUNTS_SC)}, ux=ux, uy=uy, origin=org,
                      color=ORANGE, width=0.4, labels=True, font_size=26)
    base = Line(org + 0.7 * ux * LEFT, org + 4.7 * ux * RIGHT).set_stroke(GREY_B, 2)
    xl = VGroup(*[tx(t, 28, GREY_A).move_to(org + m * ux * RIGHT + 0.32 * DOWN)
                  for m, t in zip(MS, ["0", "1", "2", "3", r"\ge4"])])
    xname = zh("回家次數", font_size=24, color=GREY_A).next_to(base, RIGHT, buff=0.2)
    title = at(mixed("簡單立方：1000 步內回家 ", "$m$", " 次的人數", font_size=28), -2.4, 2.75)
    sw = lambda c: Square(side_length=0.25).set_fill(c, 0.8).set_stroke(c, 1.5)
    legend = VGroup(VGroup(sw(D.BAR_COLOR), mixed("理論 ", "$1000\\,f_m$", font_size=24)).arrange(RIGHT, buff=0.2),
                    VGroup(sw(ORANGE), zh("模擬：1000 人各走 1000 步", font_size=24)).arrange(RIGHT, buff=0.2))
    legend.arrange(DOWN, aligned_edge=LEFT, buff=0.2)
    at(legend, 6.6 - legend.get_width() / 2, 1.9)
    assert_in_frame(th, sim, base, xl, xname, title, legend)
    assert_apart(title, sim)
    assert_apart(title, legend)
    assert_apart(legend, VGroup(*[b for b in th[0]][1:], *[b for b in sim[0]][1:], *sim[1][1:]))
    s.play(ShowCreation(base), FadeIn(xl), FadeIn(xname), FadeIn(title), run_time=1.0)
    s.say("1000 人各走 1000 步：數每個人回家幾次", wait=0.5)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in sim[0]], lag_ratio=0.15), run_time=1.6)
    s.play(FadeIn(sim[1]), FadeIn(legend[1]), run_time=0.7)
    s.say("672 人一次都沒回來，222 人回來一次", wait=3.5)
    s.play(LaggedStart(*[GrowFromEdge(b, DOWN) for b in th[0]], lag_ratio=0.15), FadeIn(legend[0]), run_time=1.6)
    s.say("藍色是理論 ", "$1000\\,f_m$", "：和模擬幾乎一樣高", wait=3.5)
    s.say("$m\\ge4$", "：理論 ", "$\\sum_{m\\ge4}f_m=R^4=" + nb["R4"] + "$", "，模擬 18 人", wait=3.5)
    s.say("三分之二的人一次也沒回來", wait=3.0)


def part_close(s):
    s.say("鄰居越多逃得越遠：", "$" + "<".join(numbers()["esc3"]) + "$", wait=3.5)
    s.say("但每種三維格子都有一定比例的人永遠回不了家", wait=3.5)


class W306Watson(ChuScene):
    def construct(self):
        self.show_header("W306　三種立方格子的逃脫機率")
        check_math()
        lats = part_lattices(self)
        part_reveal(self, lats)
        part_counts(self)
        part_bars(self)
        part_close(self)
        self.wait(1.0)
