"""
W305 維度的分水嶺（§3.6，pp. 58–62，(3.26)–(3.33)，Ex. 3.5–3.6）

問：為什麼二維和三維之間會分出「一定回家」與「可能永遠回不了家」？
  1. (3.27)(3.28)：q 小時 1 − zχ(q) ≈ c1(z_c − z) + c2 q²/c1（c1 = χ(0) = Z）；筆誤框 c2 = 1/2 → c2 = 1。
  2. 角度積分 d^dq = K_d q^{d−1} dq，K_d = 2π^{d/2}/Γ(d/2)（(3.29)，R.K_d）；練習 3.5：∫e^{−r²}d^dr = π^{d/2}。
  3. z = z_c 時被積函數 ∝ q^{d−3}：三條曲線＋表 ∫_ε^1 q^{d−3} dq（R.tail_integral）→ d ≤ 2 發散、d = 3 有限。
  4. G(z;0) 的表（R.G0，z = x·z_c）與奇異部分 (z_c − z)^{(d−2)/2}（(3.32)）。
  5. (3.33) P_∞ = 1/G(z_c;0) = 0.659463；練習 3.6：P_N ≈ 0.659463 + 0.28696/√N（R.sc_correction_coeff）。
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
from utils import walk_draw as D

DIM_COLORS = {1: BLUE_C, 2: GREEN_C, 3: ORANGE}
KIND = {1: "1d", 2: "sq", 3: "sc"}
XS = (0.9, 0.99, 0.999)
EPS = (0.1, 0.01, 0.001)
NS = (100, 1000)


def g_fmt(v):
    return f"{v:.2f}" if v >= 10 else f"{v:.3f}"


def t_fmt(d, v):
    return str(round(v)) if d == 1 else (f"{v:.3f}" if d == 2 else f"{v:g}")


@lru_cache(maxsize=None)
def numbers():
    nb = {"K": [R.K_d(d) for d in (1, 2, 3, 4)]}
    nb["tail"] = {d: [t_fmt(d, R.tail_integral(d, e)) for e in EPS] for d in (1, 2, 3)}
    nb["G"] = {d: [g_fmt(R.G0(KIND[d], x)) for x in XS] for d in (1, 2, 3)}
    nb["Gc"] = R.G_C["sc"]
    nb["Pinf"] = 1 / R.G_C["sc"]
    nb["A"] = R.sc_correction_coeff()
    P = R.no_return_probs("sc", 1000)
    nb["PN"] = {N: (nb["Pinf"] + nb["A"] / math.sqrt(N), P[N]) for N in NS}
    return nb


def check_math():
    nb = numbers()
    K = nb["K"]
    assert K[0] == 2 and [f"{k:.3f}" for k in K[1:]] == ["6.283", "12.566", "19.739"]
    for d, v in ((2, 2 * math.pi), (3, 4 * math.pi), (4, 2 * math.pi ** 2)):
        assert abs(K[d - 1] - v) < 1e-12
    assert all(abs(R.K_d_gauss(d) - R.K_d(d)) < 1e-8 for d in (1, 2, 3))           # 練習 3.5
    assert nb["tail"] == {1: ["9", "99", "999"], 2: ["2.303", "4.605", "6.908"], 3: ["0.9", "0.99", "0.999"]}
    assert abs(R.tail_integral(3, 1e-12) - 1) < 1e-9 and R.tail_integral(2, 1e-12) > 27
    assert R.tail_integral(1, 1e-12) > 1e11 and abs(R.tail_integral_numeric(2, 0.01) - 4.605170186) < 1e-6
    assert abs(R.tail_integral(2, 0.01) - R.tail_integral(2, 0.1) - math.log(10)) < 1e-12
    assert nb["G"] == {1: ["2.294", "7.089", "22.37"], 2: ["1.452", "2.137", "2.862"], 3: ["1.233", "1.409", "1.480"]}
    assert f"{nb['Gc']:.6f}" == "1.516386" and abs(R.G0("sc", 1) - nb["Gc"]) < 1e-9 and f"{R.G0('sc', 1):.3f}" == "1.516"
    g1 = [R.G0("1d", x) for x in XS]
    g2 = [R.G0("sq", x) for x in XS]
    assert abs(g1[2] / g1[1] - math.sqrt(10)) < 0.01 and abs(g2[2] - g2[1] - math.log(10) / math.pi) < 0.01
    assert f"{math.log(10) / math.pi:.2f}" == "0.73"
    assert R.G0("1d", 1 - 1e-12) > 1e5 and R.G0("sq", 1 - 1e-12) > 9                  # 一、二維在 z_c 發散
    assert f"{nb['Pinf']:.6f}" == "0.659463" and f"{1 - nb['Pinf']:.4f}" == "0.3405"
    assert f"{nb['A']:.5f}" == "0.28696" and abs(R.B1_SC - 9 / math.pi) < 1e-15
    assert abs(nb["A"] - R.B1_SC * math.sqrt(1 / 6) / (nb["Gc"] ** 2 * math.sqrt(math.pi))) < 1e-12
    shown = {N: tuple(f"{v:.6f}" for v in nb["PN"][N]) for N in NS}
    assert shown == {100: ("0.688159", "0.688062"), 1000: ("0.668537", "0.668534")}
    assert f"{nb['PN'][100][0] - nb['PN'][100][1]:.4f}" == "0.0001"
    assert f"{nb['PN'][1000][0] - nb['PN'][1000][1]:.6f}" == "0.000003"


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
    """表格裡左右、上下相鄰的格子不重疊。"""
    cells = tab[0]
    nrows = len(cells) // ncols
    for i in range(nrows):
        for j in range(ncols):
            if j + 1 < ncols:
                assert_apart(cells[i * ncols + j], cells[i * ncols + j + 1], gap)
            if i + 1 < nrows:
                assert_apart(cells[i * ncols + j], cells[(i + 1) * ncols + j], gap)


def assert_clear_of(box_mob, pts, gap=0.05):
    """曲線的點都不在 box_mob 的外框裡。"""
    x0, x1 = box_mob.get_left()[0] - gap, box_mob.get_right()[0] + gap
    y0, y1 = box_mob.get_bottom()[1] - gap, box_mob.get_top()[1] + gap
    assert not any(x0 <= p[0] <= x1 and y0 <= p[1] <= y1 for p in pts), box_mob


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    head = at(zh("走無窮多步，回過家的機率", font_size=30, color=GREY_A), 0, 2.6)
    cols = VGroup()
    for d, x, name, val, fate in ((1, -4.6, "一維", "1", "一定回家"), (2, -1.6, "二維", "1", "一定回家"),
                                  (3, 3.4, "三維", f"{1 - nb['Pinf']:.4f}", "可能永遠回不了家")):
        c = DIM_COLORS[d]
        col = VGroup(zh(name, font_size=36, color=c), tx(val, 52, c), zh(fate, font_size=26, color=GREY_A))
        cols.add(at(col.arrange(DOWN, buff=0.35), x, 0.7))
    line = DashedLine(np.array([0.9, 1.9, 0]), np.array([0.9, -0.9, 0])).set_stroke(YELLOW, 3)
    lab = zh("分水嶺", font_size=30, color=YELLOW).next_to(line, DOWN, buff=0.2)
    assert_in_frame(head, cols, line, lab)
    assert_apart(cols[1], line)
    assert_apart(line, cols[2])
    assert_apart(head, cols)
    s.play(FadeIn(head), run_time=0.6)
    s.play(LaggedStart(*[FadeIn(c, shift=0.2 * UP) for c in cols], lag_ratio=0.3), run_time=1.6)
    s.say("問：一、二維一定回家，三維卻可能永遠回不了家", wait=3.0)
    s.play(ShowCreation(line), FadeIn(lab), run_time=1.0)
    s.say("分水嶺為什麼剛好落在二維和三維之間？", wait=3.0)
    s.clear_stage()


def part_expand(s):
    e1 = at(tx(r"G(z;\vec 0)=\frac{1}{(2\pi)^d}\int\frac{d^dq}{1-z\chi(\vec q)}", 36), 0, 2.45)
    e2 = at(eq_label(tx(r"1-z\chi(\vec q)\approx c_1(z_c-z)+\frac{c_2}{c_1}\,q^2", 36), "(3.27)(3.28)"), 0, 1.15)
    sub = at(tx(r"\chi(\vec q)\approx Z-c_2q^2,\qquad c_1=\chi(\vec 0)=Z,\qquad z_c=1/Z", 32), 0, 0.05)
    box = at(D.erratum_box("p. 59，(3.27) 後", r"c_2=1/2", r"c_2=1\ (\chi=2\sum\cos q_i\approx 2d-q^2)",
                           width=7.6), 0, -1.75)
    assert_in_frame(e1, e2, sub, box)
    assert_apart(e1, e2)
    assert_apart(e2, sub)
    assert_apart(sub, box)
    s.play(Write(e1), run_time=1.5)
    s.say("回到原點的生成函數：第 2 章的傅立葉積分", wait=3.0)
    s.say("$z\\to z_c$", " 時分母在 ", "$\\vec q=\\vec 0$", " 變成 0：麻煩全在那裡", wait=3.0)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=1.0)
    s.say("$q$", " 小時展開：一項管 ", "$z_c-z$", "，一項管 ", "$q^2$", wait=3.0)
    s.play(FadeIn(sub), run_time=0.8)
    s.say("正方、立方格子：", "$\\chi=2\\sum\\cos q_i\\approx2d-q^2$", "，所以 ", "$c_2=1$", wait=3.5)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上寫 ", "$c_2=1/2$", "；照 ", "$\\chi$", " 的展開應是 ", "$c_2=1$", wait=3.5)
    s.play(FadeOut(VGroup(e1, sub, box)), e2.animate.move_to(np.array([0, 2.6, 0])), run_time=0.9)

    rad = at(VGroup(tx(r"d^dq=K_d\,q^{d-1}dq", 34), tx(r"K_d=\frac{2\pi^{d/2}}{\Gamma(d/2)}", 34),
                    tx(r"\text{(3.29)}", 22, GREY_B)).arrange(RIGHT, buff=0.6), 0, 1.4)
    rint = at(tx(r"G(z;\vec 0)\approx\text{const}+\frac{K_d}{(2\pi)^d}\int_0^{q_0}"
                 r"\frac{q^{d-1}\,dq}{c_1(z_c-z)+c_2q^2/c_1}", 34), 0, 0.1)
    rows = [["d", "1", "2", "3", "4"],
            ["K_d", "2", r"2\pi\approx6.283", r"4\pi\approx12.566", r"2\pi^2\approx19.739"]]
    tab = D.tex_table(rows, [1.2, 1.0, 2.6, 2.8, 2.9], row_h=0.55, font_size=30)
    at(tab, 0, -1.3)
    ex = at(mixed("練習 3.5：", "$\\pi^{d/2}=\\int e^{-r^2}d^dr=K_d\\int_0^\\infty r^{d-1}e^{-r^2}dr"
                  "=K_d\\,\\frac{\\Gamma(d/2)}{2}$", font_size=24), 0, -2.45)
    fit(ex, 13.0)
    assert_in_frame(e2, rad, rint, tab, ex)
    assert_table(tab, 5)
    for a, b in ((e2, rad), (rad, rint), (rint, tab), (tab, ex)):
        assert_apart(a, b)
    s.say("只看 ", "$q<q_0$", " 的小球：先做角度積分", wait=0.5)
    s.play(FadeIn(rad), run_time=1.0)
    s.play(FadeIn(rint, shift=0.2 * DOWN), run_time=1.0)
    s.say("$K_d$", "＝", "$d$", " 維單位球面的面積", wait=2.5)
    s.play(FadeIn(tab[0]), ShowCreation(tab[1]), run_time=1.2)
    s.say("$K_1=2$", "（兩個端點）、", "$K_2=2\\pi$", "（圓周）、", "$K_3=4\\pi$", "（球面）", wait=3.5)
    s.play(FadeIn(ex, shift=0.2 * UP), run_time=0.8)
    s.say("練習 3.5：高斯積分算兩次，就得到 ", "$K_d$", wait=3.0)
    s.play(Indicate(rint, scale_factor=1.05), run_time=1.0)
    s.say("剩下一維的積分：分子 ", "$q^{d-1}$", " 對上分母的 ", "$q^2$", wait=3.5)
    s.clear_stage()


def part_divide(s):
    nb = numbers()
    ax = D.plot_axes((0, 1.05, 0.25), (0, 10.5, 2), width=5.4, height=4.3, x_ticks=[0.5, 1],
                     y_ticks=[2, 4, 6, 8, 10], font_size=22, x_label="q")
    ax.move_to(np.array([-3.75, -0.3, 0]))
    title = mixed("$z=z_c$", "：被積函數 ", "$\\dfrac{q^{d-1}}{q^2}=q^{d-3}$", font_size=26)
    title.next_to(ax, UP, buff=0.15).align_to(ax, LEFT)
    curves, legend = VGroup(), VGroup()
    for d in (1, 2, 3):
        q0 = {1: 10.5 ** -0.5, 2: 1 / 10.5, 3: 0.0}[d]
        pts = [ax.c2p(q, q ** (d - 3)) for q in np.linspace(q0, 1.0, 160)]
        curves.add(VMobject().set_points_as_corners(pts).set_stroke(DIM_COLORS[d], 4))
        legend.add(tx(f"d={d}:\\ q^{{{d - 3}}}", 28, DIM_COLORS[d]))
    legend.arrange(DOWN, aligned_edge=LEFT, buff=0.25).move_to(ax.c2p(0.8, 8.2))
    hdr = [tx(r"\int_\varepsilon^1q^{d-3}dq", 26, YELLOW), r"\varepsilon=0.1", "0.01", "0.001", r"\varepsilon\to0"]
    rows = [hdr] + [[tx(f"d={d}", 28, DIM_COLORS[d])] + nb["tail"][d] + [r"\infty" if d < 3 else "1"]
                    for d in (1, 2, 3)]
    tab = D.tex_table(rows, [1.9, 1.3, 1.0, 1.1, 1.1], row_h=0.7, font_size=28)
    tab.move_to(np.array([0.35 + tab.get_width() / 2, 2.25 - tab.get_height() / 2, 0]))
    v1 = mixed("$d\\le2$", "：在 ", "$q=0$", " 發散", font_size=28, color=GREEN_C)
    v2 = mixed("$d=3$", "：有限", font_size=28, color=ORANGE)
    verdict = VGroup(v1, v2).arrange(DOWN, aligned_edge=LEFT, buff=0.3).next_to(tab, DOWN, buff=0.5)
    verdict.align_to(tab, LEFT).shift(0.3 * RIGHT)
    assert_in_frame(ax, title, legend, tab, verdict)
    assert_table(tab, 5)
    assert_apart(ax, tab)
    assert_apart(tab, verdict)
    assert_apart(title, tab)
    for c in curves:
        assert_clear_of(legend, c.get_points())
    late = VGroup(*[D.table_cell(tab, i, 4, 5) for i in range(4)])
    early = VGroup(*[c for c in tab[0] if c not in late.submobjects])
    s.play(FadeIn(ax), FadeIn(title), run_time=1.0)
    s.say("在 ", "$z=z_c$", "：分母只剩 ", "$c_2q^2/c_1$", wait=2.5)
    s.play(LaggedStart(*[ShowCreation(c) for c in curves], lag_ratio=0.4), FadeIn(legend), run_time=2.4)
    s.say("一維 ", "$q^{-2}$", "、二維 ", "$q^{-1}$", "、三維是常數", wait=3.0)
    s.play(FadeIn(early), ShowCreation(tab[1]), run_time=1.4)
    s.say("$\\varepsilon$", " 每小 10 倍：一維約乘 10，二維加 2.303", wait=4.0)
    s.say("三維：0.9、0.99、0.999，收斂到 1", wait=2.5)
    s.play(LaggedStart(*[FadeIn(c) for c in late], lag_ratio=0.2), run_time=1.0)
    s.play(FadeIn(verdict, shift=0.2 * UP), run_time=0.8)
    s.say("$d\\le2$", "：", "$q=0$", " 附近積分發散；", "$d=3$", "：有限", wait=3.5)
    s.clear_stage()


def part_G(s):
    nb = numbers()
    head = [zh("格子", font_size=24, color=YELLOW), "x=0.9", "0.99", "0.999", "1", zh("奇異部分", font_size=24, color=YELLOW)]
    sing = {1: r"(z_c-z)^{-1/2}", 2: r"-\ln(z_c-z)", 3: r"-(z_c-z)^{1/2}"}
    names = {1: "一維", 2: "二維", 3: "三維 sc"}
    rows = [head] + [[zh(names[d], font_size=26, color=DIM_COLORS[d])] + nb["G"][d]
                     + [r"\infty" if d < 3 else f"{nb['Gc']:.3f}", sing[d]] for d in (1, 2, 3)]
    tab = D.tex_table(rows, [1.9, 1.5, 1.5, 1.5, 1.3, 3.0], row_h=0.72, font_size=30)
    tab.move_to(np.array([0, 2.1 - tab.get_height() / 2, 0]))
    title = mixed("$G(z;\\vec 0)$", "，", "$z=x\\,z_c$", font_size=28).next_to(tab, UP, buff=0.12)
    law = at(eq_label(tx(r"G(z;\vec 0)\approx\text{const}+\text{const}\cdot(z_c-z)^{(d-2)/2}", 34), "(3.32)"), 0, -1.55)
    note = at(mixed("$d=2$", "：次方是 0，換成對數", font_size=26, color=GREY_A), 0, -2.45)
    assert_in_frame(tab, title, law, note)
    assert_table(tab, 6)
    assert_apart(title, tab)
    assert_apart(tab, law)
    assert_apart(law, note)
    col = lambda j: VGroup(*[D.table_cell(tab, i, j, 6) for i in range(4)])
    s.play(FadeIn(title), FadeIn(col(0)), ShowCreation(tab[1]), run_time=1.0)
    s.say("$G(z;\\vec 0)$", " 本身：令 ", "$z=x\\,z_c$", "，讓 ", "$x$", " 逼近 1", wait=1.0)
    s.play(LaggedStart(*[FadeIn(col(j)) for j in (1, 2, 3)], lag_ratio=0.5), run_time=2.0)
    s.say("一維：每多一個 9，大約乘上 ", "$\\sqrt{10}$", wait=3.0)
    s.say("二維：每多一個 9，只多加 ", "$\\ln10/\\pi\\approx0.73$", wait=3.0)
    s.play(FadeIn(col(4)), run_time=0.8)
    s.say("三維：停在 1.516；到了 ", "$z_c$", " 也不會變成無限大", wait=3.0)
    s.play(FadeIn(col(5)), FadeIn(law, shift=0.2 * UP), run_time=1.2)
    s.say("(3.32)：奇異部分 ", "$(z_c-z)^{(d-2)/2}$", wait=3.0)
    s.play(FadeIn(note), run_time=0.6)
    s.say("$d=1$", " 是 ", "$-\\tfrac12$", " 次方；", "$d=3$", " 是 ", "$+\\tfrac12$", " 次方：值有限", wait=3.5)
    s.clear_stage()


def part_escape(s):
    nb = numbers()
    e1 = eq_label(tx(r"P_\infty=\frac{1}{G(z_c;\vec 0)}", 36), "(3.33)")
    e2 = tx(r"\mathrm{sc}:\ \frac{1}{" + f"{nb['Gc']:.6f}" + "}=" + f"{nb['Pinf']:.6f}", 36, ORANGE)
    top = at(VGroup(e1, e2).arrange(RIGHT, buff=1.0), 0, 2.5)
    e3 = at(mixed("練習 3.6：", "$G(z;\\vec 0)\\approx B_0-B_1\\sqrt{z_c-z}$", "，", f"$B_0={nb['Gc']:.6f}$",
                  "，", "$B_1=9/\\pi$", font_size=26), 0, 1.35)
    e4 = at(tx(r"P_N\approx\frac{1}{B_0}+\frac{B_1\sqrt{z_c}}{B_0^2\sqrt{\pi N}}=" + f"{nb['Pinf']:.6f}"
               + r"+\frac{" + f"{nb['A']:.5f}" + r"}{\sqrt N}", 36), 0, 0.15)
    rows = [["N", zh("近似", font_size=26, color=YELLOW), zh("精確", font_size=26, color=YELLOW)]]
    rows += [[str(N)] + [f"{v:.6f}" for v in nb["PN"][N]] for N in NS]
    tab = D.tex_table(rows, [1.6, 2.4, 2.4], row_h=0.55, font_size=30)
    tab.move_to(np.array([0, -0.85 - tab.get_height() / 2, 0]))
    assert_in_frame(top, e3, e4, tab)
    assert_table(tab, 3)
    assert_apart(e1, e2)
    for a, b in ((top, e3), (e3, e4), (e4, tab)):
        assert_apart(a, b)
    s.play(FadeIn(e1), run_time=0.8)
    s.say("$z\\to z_c$", " 時只剩 ", "$g$", " 的極點：", "$P_\\infty=1/G(z_c;\\vec 0)$", wait=3.0)
    s.say("一、二維 ", "$G(z_c;\\vec 0)=\\infty$", "，所以 ", "$P_\\infty=0$", "：一定回家", wait=3.0)
    s.play(FadeIn(e2, shift=0.2 * LEFT), run_time=0.8)
    s.say("簡單立方：約 66% 的人永遠回不了家", wait=3.0)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("練習 3.6：三維的 ", "$G$", " 在 ", "$z_c$", " 附近多一項 ", "$\\sqrt{z_c-z}$", wait=3.0)
    s.play(FadeIn(e4, shift=0.2 * DOWN), run_time=1.0)
    s.say("這一項給出 ", "$1/\\sqrt N$", " 的修正", wait=2.5)
    s.play(FadeIn(tab), run_time=1.0)
    s.say("$N=100$", " 差約 0.0001；", "$N=1000$", " 只差 0.000003", wait=3.5)


def part_close(s):
    s.say("$G(z_c;\\vec 0)$", " 有限 ", "$\\Leftrightarrow$", " 有機會永遠逃走", wait=3.0)
    s.say("分水嶺：", "$\\int q^{d-3}dq$", " 在 ", "$d=2$", " 剛好發散", wait=3.0)


class W305DimensionDivide(ChuScene):
    def construct(self):
        self.show_header("W305　維度的分水嶺")
        check_math()
        part_question(self)
        part_expand(self)
        part_divide(self)
        part_G(self)
        part_escape(self)
        part_close(self)
        self.wait(1.0)
