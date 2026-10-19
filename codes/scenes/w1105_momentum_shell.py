"""
W1105 動量殼與連續重整化（§11.2.1–11.2.2，pp. 272–277，(11.64)–(11.87)，Fig. 11.8）

問：在 Fourier 空間「積掉一半」長什麼樣？
  1. Fig. 11.8：格距 c → 2c ⇔ 布里淵區 Q → Q/2；外圈殼層塗色後淡出；拿掉的比例 1 − 2^{−d}（d = 2：3/4、d = 3：7/8）。
  2. (11.64)–(11.67)：Fourier 模態不耦合；α(q) 表（N = 8）；筆誤框 1（(11.66)(11.67) 的 α）。
  3. (11.68)–(11.73)：q → q/2、N → 2^d N、S → 2S ⇒ r → 4r、h → 2^{(d+2)/2} h；與 W1102 的 (11.22) 一樣。
  4. (11.74)–(11.83)：b = e^{δℓ} ⇒ dr/dℓ = 2r、dh/dℓ = (d+2)h/2；(r, h) 平面的流線；表 ℓ = 1、2、ln 2。
  5. (11.84)–(11.87)：m 與 ℓ 無關 ⇒ 偏微分方程；(11.24) 代入三項相消；收束。
畫面上的數字都由 utils.walk_rg 算出，check_math() 核對。
模組層級的畫圖函式：real_lattice（實空間格點）、bz_shell（殼層）、flow_axes／flow_curve（流線）、check_mark。
"""
import math
import sys
import warnings
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rg as RG
from utils import lattice_draw as L
from utils import walk_draw as D

SHELL_COLOR = TEAL
OK_COLOR = "#7CE38B"
FLOW_D = 3
FLOW_STARTS = ((0.01, 0.01), (-0.01, 0.01), (0.03, 0.003), (-0.03, 0.003))
FLOW_RMAX = 1.0
FLOW_HMAX = 1.5
FLOW_COLORS = (YELLOW, ORANGE, TEAL, BLUE_B)


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


def ell_max(r0, h0, d=FLOW_D):
    """流線畫到 |r| ≤ 1 且 h ≤ 1.5 為止。"""
    return min(math.log(FLOW_RMAX / abs(r0)) / 2, math.log(FLOW_HMAX / h0) / ((d + 2) / 2))


def landau_residual(d=3, r=0.4, h=0.3, eps=1e-5):
    """Landau 的組合 m = r^{1/2} f(h/r^{3/2}) 放進 (11.87)（取 d）的殘差 ÷ |m|（f(x) = x/(1+x²)）。"""
    m = lambda rr, hh: rr ** 0.5 * (hh / rr ** 1.5) / (1 + (hh / rr ** 1.5) ** 2)
    dr = (m(r + eps, h) - m(r - eps, h)) / (2 * eps)
    dh = (m(r, h + eps) - m(r, h - eps)) / (2 * eps)
    return abs(-(d - 2) / 2 * m(r, h) + 2 * r * dr + (d + 2) / 2 * h * dh) / abs(m(r, h))


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {"frac": {d: RG.shell_fraction(d) for d in (2, 3)},
                "alpha": RG.alpha_table(),
                "fourier": RG.fourier_check(),
                "shell": RG.shell_rg(1.0, 1.0, 3),
                "steps": RG.shell_rg_steps(1.0, 1.0, 3),
                "flow": {t: RG.flow(0.01, 0.01, FLOW_D, t) for t in (1.0, 2.0, math.log(2))},
                "flow_num": {t: RG.flow_numeric(0.01, 0.01, FLOW_D, t) for t in (1.0, 2.0)},
                "b2": RG.shell_rg(0.01, 0.01, FLOW_D),
                "lines": [RG.flow_line(r0, h0, FLOW_D, ell_max(r0, h0)) for r0, h0 in FLOW_STARTS],
                "pde": {d: RG.pde_residual(d) for d in (1, 3, 4)}}


def check_math():
    nb = numbers()
    # §1 殼層
    assert nb["frac"] == {2: Fraction(3, 4), 3: Fraction(7, 8)}
    # §2 α(q)（N = 8）與 Fourier 核對（種子 1105）
    tab = nb["alpha"]
    assert [k for k, *_ in tab] == [0, 1, 2, 3, 4]
    assert [f4(q) for _, q, *_ in tab] == ["0.0000", "0.7854", "1.5708", "2.3562", "3.1416"]
    assert [f4(a) for _, _, a, _, _ in tab] == ["0.0000", "0.2929", "1.0000", "1.7071", "2.0000"]
    assert [f4(a) for _, _, _, a, _ in tab] == ["0.0000", "0.3084", "1.2337", "2.7758", "4.9348"]
    err, real, book = nb["fourier"]
    assert err < 1e-12 and f4(real) == "3.0394" and f4(real * (1 + err)) == "3.0394" and f4(book) == "-4.3508"
    # §3 b = 2
    assert nb["shell"][0] == 4.0 and f4(nb["shell"][1]) == "5.6569" and f4(2 ** 2.5) == "5.6569"
    assert all(abs(a - b) < 1e-12 for a, b in zip(nb["shell"], nb["steps"]))
    assert all(RG.gauss_rescale_check(d) < 1e-12 for d in (1, 3, 4))
    # §4 流
    r1, h1 = nb["flow"][1.0]
    r2, h2 = nb["flow"][2.0]
    assert [f4(v) for v in (r1, h1, r2, h2)] == ["0.0739", "0.1218", "0.5460", "1.4841"]
    for t in (1.0, 2.0):
        assert all(abs(a - b) < 1e-8 for a, b in zip(nb["flow"][t], nb["flow_num"][t]))
    rl, hl = nb["flow"][math.log(2)]
    assert abs(rl - nb["b2"][0]) < 1e-14 and abs(hl - nb["b2"][1]) < 1e-14 and [f4(rl), f4(hl)] == ["0.0400", "0.0566"]
    for (r0, h0), pts in zip(FLOW_STARTS, nb["lines"]):
        assert len(pts) == 40 and pts[0] == (r0, h0)
        assert all(abs(r) <= FLOW_RMAX + 1e-9 and 0 < h <= FLOW_HMAX + 1e-9 for r, h in pts)
        c = h0 / abs(r0) ** 1.25
        assert all(abs(h / abs(r) ** 1.25 - c) < 1e-9 * c for r, h in pts)          # h ∝ |r|^{5/4}
        assert abs(pts[-1][0]) > 0.5 or pts[-1][1] > 1.0                             # 都離開原點
    # §5 (11.87)
    assert all(nb["pde"][d] < 1e-7 for d in (1, 3, 4))
    assert all(RG.char_check(d) < 1e-9 for d in (1, 3, 4))
    assert landau_residual(3) > 0.1                                                  # Landau 的組合在 d = 3 不成立


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
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


def framed(mob, color, buff=0.25):
    return VGroup(SurroundingRectangle(mob, buff=buff).set_stroke(color, 2), mob)


def check_mark(color=OK_COLOR, size=0.3):
    return VMobject().set_points_as_corners([size * np.array([-0.5, 0.0, 0]), size * np.array([-0.15, -0.4, 0]),
                                             size * np.array([0.55, 0.45, 0])]).set_stroke(color, 5)


# ---------------------------------------------------------------- 畫圖元件
def real_lattice(origin, c=0.7, n=2):
    """(2n+1)² 個格點，格距 c：回傳 (保留的 VGroup, 拿掉的 VGroup)（保留 = 兩個座標都是偶數）。"""
    keep, drop = VGroup(), VGroup()
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            d = Dot(np.array(origin) + c * np.array([i, j, 0]), radius=0.07).set_color(WHITE)
            (keep if i % 2 == 0 and j % 2 == 0 else drop).add(d)
    return keep, drop


def bz_shell(width, origin, color=SHELL_COLOR):
    """外圈殼層（|q| 介於 Q/2 與 Q）：四塊長方形；內框（Q/2）。"""
    o = np.array(origin, dtype=float)
    w, h = width, width / 4
    strips = VGroup(Rectangle(width=w, height=h).move_to(o + (w / 2 - h / 2) * UP),
                    Rectangle(width=w, height=h).move_to(o + (w / 2 - h / 2) * DOWN),
                    Rectangle(width=h, height=w / 2).move_to(o + (w / 2 - h / 2) * LEFT),
                    Rectangle(width=h, height=w / 2).move_to(o + (w / 2 - h / 2) * RIGHT))
    strips.set_stroke(width=0).set_fill(color, 0.45)
    inner = Square(side_length=width / 2).move_to(o).set_stroke(YELLOW, 3)
    return strips, inner


def flow_axes():
    ax = D.plot_axes((-1.1, 1.15, 0.5), (0, 1.65, 0.5), width=5.6, height=4.4, x_ticks=(-1, -0.5, 0.5, 1),
                     x_tick_labels=["-1", "-0.5", "0.5", "1"], x_label="r", y_label="h")
    ax.add(VGroup(*[tx(s, 20, GREY_B).next_to(ax.c2p(0, v), RIGHT, buff=0.1)
                    for s, v in (("0.5", 0.5), ("1", 1.0), ("1.5", 1.5))]))
    return ax


def flow_curve(ax, pts, color):
    """一條流線（40 點）＋終點箭頭。"""
    ps = [ax.c2p(r, h) for r, h in pts]
    line = VMobject().set_points_smoothly(ps).set_stroke(color, 4)
    u = (ps[-1] - ps[-2]) / np.linalg.norm(ps[-1] - ps[-2])
    tip = Arrow(ps[-1] - 0.35 * u, ps[-1], buff=0, thickness=4).set_color(color)
    return VGroup(line, tip)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    keep, drop = real_lattice(np.array([-3.6, 0.8, 0]))
    cl = tx("c", 30, GREY_A).next_to(VGroup(drop[0], keep[0]), LEFT, buff=0.15)
    c2 = tx("2c", 30, YELLOW).move_to(cl)
    rlab = zh("實空間", font_size=26, color=GREY_A).move_to(np.array([-3.6, 2.85, 0]))
    bw = 3.4
    bo = np.array([3.0, 0.8, 0])
    bz = L.bz_frame(bw, bo)
    strips, inner = bz_shell(bw, bo)
    qlab = tx(r"\pm Q/2", 24, YELLOW).next_to(inner, UR, buff=0.02)
    klab = zh("布里淵區", font_size=26, color=GREY_A).move_to(np.array([3.0, 2.95, 0]))
    frac = mixed("拿掉的比例 ", "$1-2^{-d}$", "：", "$d=2$", "：", "$\\tfrac34$", "，", "$d=3$", "：", "$\\tfrac78$",
                 font_size=26)
    frac.move_to(np.array([0, -2.45, 0]))
    assert_in_frame(keep, drop, cl, rlab, bz, klab, frac, qlab)
    assert_apart(VGroup(keep, drop, cl), bz)
    assert_apart(bz, frac)
    assert_apart(VGroup(keep, drop), frac)
    s.say("問：在 Fourier 空間「積掉一半」長什麼樣？", wait=1.8)
    s.play(FadeIn(rlab), FadeIn(keep), FadeIn(drop), FadeIn(cl), run_time=1.0)
    s.play(FadeIn(klab), ShowCreation(bz), run_time=1.0)
    s.say("實空間每邊只留一半格點，", "$N$", " 變成 ", "$N/2^d$", wait=2.0)
    s.say("格距 ", "$c$", " 的格子：波數只到 ", "$|q|\\le Q=\\pi/c$", wait=2.0)
    s.play(FadeOut(drop), FadeIn(strips), ShowCreation(inner), FadeIn(qlab), Transform(cl, c2), run_time=1.4)
    s.say("Fig. 11.8：格距 ", "$c\\to2c$", " ", "$\\Leftrightarrow$", " 布里淵區 ", "$Q\\to Q/2$", wait=2.2)
    s.say("外圈的殼層：短波長的漲落，就是要積掉的部分", wait=2.2)
    s.play(FadeIn(frac), run_time=0.8)
    s.say("$d=2$", " 拿掉 ", "$\\tfrac34$", "，", "$d=3$", " 拿掉 ", "$\\tfrac78$", " 的模態", wait=2.0)
    s.play(FadeOut(strips), run_time=1.0)
    s.say("Gaussian 模型的好處：殼層和內部完全不耦合", wait=2.0)
    s.clear_stage()


def part_modes(s):
    nb = numbers()
    err, real, book = nb["fourier"]
    eqs = VGroup(tx(r"S_i=N^{-1/2}\sum_qS(q)e^{iqR_i}", 32),
                 tx(r"\mathcal H=\sum_q\Big[\frac r2+\alpha(q)\Big]S(q)S(-q)-\sqrt N S(0)h", 32)).arrange(RIGHT, buff=0.7)
    tag = zh("(11.64)–(11.67)　一維鏈 N = 8，a = 1", font_size=22, color=GREY_B)
    rows = [["k"] + [str(k) for k, *_ in nb["alpha"]],
            ["q=2\\pi k/8"] + [f4(q) for _, q, *_ in nb["alpha"]],
            [r"\alpha=a(1-\cos q)"] + [f4(a) for _, _, a, _, _ in nb["alpha"]],
            [r"\tfrac12aq^2"] + [f4(a) for _, _, _, a, _ in nb["alpha"]]]
    tab = D.tex_table(rows, [3.0] + [1.5] * 5, row_h=0.55, font_size=28)
    box = D.erratum_box("p. 273，(11.66)(11.67)", r"\alpha(q)=-2a(1-\cos qc)\sim ac^2q^2",
                        r"\alpha(q)=a(1-\cos qc)\sim\frac12ac^2q^2", width=11.0, font_size=28)
    col = VGroup(tag, eqs, tab, box).arrange(DOWN, buff=0.38)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(tag, eqs, tab, box)
    assert_table(tab, 6)
    s.play(FadeIn(tag), Write(eqs[0]), run_time=1.0)
    s.say("換成 Fourier 模態 ", "$S(q)$", wait=1.6)
    s.play(Write(eqs[1]), run_time=1.2)
    s.say("每個 ", "$q$", " 各自獨立：", "$\\mathcal H$", " 是各模態的和", wait=2.0)
    s.say("外場只碰到 ", "$q=0$", " 的模態", wait=1.6)
    s.say("每個模態都是一個獨立的 Gaussian 積分", wait=2.0)
    s.play(FadeIn(tab), run_time=1.0)
    s.say("$\\alpha(q)=a(1-\\cos q)$", "：小 ", "$q$", " 才像 ", "$\\tfrac12aq^2$", wait=2.2)
    s.say("我們只關心長波長，所以只留 ", "$q^2$", " 項", wait=2.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("隨機 8 個自旋：實空間 ", "$%s$" % f4(real), "＝Fourier ", "$%s$" % f4(real * (1 + err)), wait=2.2)
    s.say("照書上的 ", "$\\alpha$", " 得 ", "$%s$" % f4(book), "：負的能量", wait=2.2)
    s.clear_stage()


def part_shell(s):
    nb = numbers()
    steps = VGroup(mixed("$\\vec q\\to\\vec q/2$", "：", "$Jq^2\\to Jq^2/4$", font_size=28),
                   mixed("$N\\to2^dN$", "：模態數回到原來", font_size=28),
                   mixed("$S\\to2S$", "：把 ", "$J$", " 調回原值", font_size=28)).arrange(DOWN, buff=0.3,
                                                                                     aligned_edge=LEFT)
    tag = zh("(11.68)–(11.73)　b = 2 的一步", font_size=22, color=GREY_B)
    res = tx(r"r\to4r,\qquad h\to2^{(d+2)/2}h", 40).set_color(YELLOW)
    num = mixed("$d=3$", "：", "$h\\times%s$" % f4(nb["shell"][1]), font_size=28)
    gen = tx(r"m(r,h)=b^{(2-d)/2}\,m\big(rb^2,\ hb^{(d+2)/2}\big)", 38)
    genbox = framed(VGroup(gen, zh("和 W1102 的 (11.22) 一樣", font_size=24, color=GREY_A)).arrange(DOWN, buff=0.2),
                    TEAL, 0.25)
    col = VGroup(tag, steps, res, num, genbox).arrange(DOWN, buff=0.38)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(tag, steps, res, num, genbox)
    s.say("積掉殼層之後，再把系統放回原樣：", wait=1.6)
    s.play(FadeIn(tag), FadeIn(steps[0]), run_time=0.8)
    s.say("波數縮一半：", "$q\\to q/2$", "，梯度項變成四分之一", wait=2.0)
    s.play(FadeIn(steps[1]), run_time=0.6)
    s.say("格點變稀了：", "$N\\to2^dN$", " 把個數補回來", wait=1.8)
    s.play(FadeIn(steps[2]), run_time=0.6)
    s.say("自旋放大兩倍，", "$J$", " 就回到原值", wait=1.8)
    s.play(Write(res), run_time=1.0)
    s.play(FadeIn(num), run_time=0.6)
    s.say("$r$", " 乘 4 與維度無關；", "$h$", " 的倍數才看 ", "$d$", wait=2.0)
    s.say("$r\\to4r$", "、", "$h\\to2^{(d+2)/2}h$", "；", "$d=3$", " 時 ", "$h$", " 乘 ", "$%s$" % f4(nb["shell"][1]), wait=2.2)
    s.play(FadeIn(genbox), run_time=0.8)
    s.say("換成任意 ", "$b$", "：正是 W1102 被動重整化的結果", wait=2.2)
    s.clear_stage()


def part_flow(s):
    nb = numbers()
    ax = flow_axes()
    ax.move_to(np.array([-3.4, 0.05, 0]))
    curves = VGroup(*[flow_curve(ax, pts, c) for pts, c in zip(nb["lines"], FLOW_COLORS)])
    origin = Dot(ax.c2p(0, 0), radius=0.09).set_color(WHITE)
    tag = mixed("(11.74)–(11.83)　", "$b=e^{\\delta\\ell}$", font_size=22, color=GREY_B)
    ch = VGroup(tx(r"\vec q\to\vec q\,(1-\delta\ell)", 28), tx(r"S\to(1+\delta\ell)S", 28),
                tx(r"N\to N(1+d\,\delta\ell)", 28)).arrange(DOWN, buff=0.18)
    odes = VGroup(tx(r"\frac{dr}{d\ell}=2r", 34), tx(r"\frac{dh}{d\ell}=\frac{d+2}{2}h", 34)).arrange(RIGHT, buff=0.6)
    for m in odes:
        m.set_color(YELLOW)
    (r1, h1), (r2, h2), (rl, hl) = nb["flow"][1.0], nb["flow"][2.0], nb["flow"][math.log(2)]
    rows = [[r"\ell", "r", "h"], ["1", f4(r1), f4(h1)], ["2", f4(r2), f4(h2)], [r"\ln2", f4(rl), f4(hl)]]
    tab = D.tex_table(rows, [1.0, 1.5, 1.5], row_h=0.5, font_size=28)
    ttl = mixed("$d=3$", "，起點 ", "$(0.01,0.01)$", font_size=22, color=GREY_A)
    right = VGroup(tag, ch, odes, ttl, tab).arrange(DOWN, buff=0.3)
    place_in(right, 0.6, 6.7, -2.85, 3.1)
    assert_in_frame(ax, curves, right)
    assert_apart(ax, right)
    assert_apart(curves, right)
    assert_column(tag, ch, odes, ttl, tab)
    assert_table(tab, 3)
    s.say("讓 ", "$b$", " 連續地變：", "$b=e^{\\delta\\ell}$", "，每次只放大一點點", wait=2.0)
    s.play(FadeIn(tag), FadeIn(ch), run_time=1.0)
    s.play(Write(odes), run_time=1.2)
    s.say("重整化變成微分方程：", "$\\frac{dr}{d\\ell}=2r$", "、", "$\\frac{dh}{d\\ell}=\\frac{d+2}2h$", wait=2.4)
    s.play(ShowCreation(ax), FadeIn(origin), run_time=0.8)
    s.play(*[ShowCreation(c[0]) for c in curves], run_time=2.0)
    s.play(*[GrowArrow(c[1]) for c in curves], run_time=0.6)
    s.say("$(r,h)$", " 平面的流線：", "$h\\propto|r|^{5/4}$", "，都離開原點", wait=2.2)
    s.say("原點 ", "$r=h=0$", " 是不動點——臨界點", wait=1.8)
    s.say("只要稍微偏離，", "$r$", " 和 ", "$h$", " 都被放大：兩個都是相關變數", wait=2.2)
    s.play(FadeIn(ttl), FadeIn(tab), run_time=1.0)
    s.say("$\\ell=1$", "：", "$r=%s$" % f4(r1), "、", "$h=%s$" % f4(h1), "；解析解＝數值積分", wait=2.2)
    s.say("$\\ell=\\ln2$", " 就是一次 ", "$b=2$", "：", "$r$", " 乘 4", wait=2.0)
    s.clear_stage()


def part_pde(s):
    eq1 = tx(r"S_i=e^{\frac{2-d}{2}\ell}S_i^{(\ell)}\ \Rightarrow\ m=e^{\frac{2-d}{2}\ell}f\big(r(\ell),h(\ell)\big)", 34)
    eq2 = tx(r"0=\Big(-\frac{d-2}{2}+2r\partial_r+\frac{d+2}{2}h\partial_h\Big)m", 36).set_color(YELLOW)
    tag = mixed("(11.84)–(11.87)　", "$m$", " 與 ", "$\\ell$", " 無關", font_size=22, color=GREY_B)
    sub = tx(r"m=r^{\frac{d-2}{4}}F(x),\quad x=h/r^{\frac{d+2}{4}}", 32)
    terms = VGroup(tx(r"-\frac{d-2}{2}\,m", 30),
                   tx(r"2r\partial_rm=\frac{d-2}{2}\,m-\frac{d+2}{2}\,r^{\frac{d-2}{4}}xF'(x)", 30),
                   tx(r"\frac{d+2}{2}h\partial_hm=\frac{d+2}{2}\,r^{\frac{d-2}{4}}xF'(x)", 30)
                   ).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
    tick = check_mark().next_to(terms, RIGHT, buff=0.5)
    num = mixed("數值（中央差分）：", "$d=1,3,4$", " 的殘差都 ", "$<10^{-7}$", font_size=24, color=GREY_A)
    col = VGroup(tag, eq1, eq2, sub, VGroup(terms, tick), num).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -2.9, 3.1)
    assert_in_frame(col)
    assert_column(tag, eq1, eq2, sub, terms, num)
    s.say("最後：自旋也連續地重新標度", wait=1.6)
    s.play(FadeIn(tag), Write(eq1), run_time=1.4)
    s.say("同一個磁化 ", "$m$", "，用哪個 ", "$\\ell$", " 去算都一樣", wait=2.0)
    s.play(Write(eq2), run_time=1.2)
    s.say("對 ", "$\\ell$", " 微分＝0：一條一階偏微分方程", wait=2.0)
    s.play(FadeIn(sub), run_time=0.8)
    s.say("把 Gaussian 的解 (11.24) 代進去", wait=1.6)
    s.play(LaggedStart(*[FadeIn(t) for t in terms], lag_ratio=0.4), run_time=1.6)
    s.say("$\\frac{d-2}2m$", " 對消、", "$xF'$", " 對消：三項加起來是 0", wait=2.2)
    s.play(ShowCreation(tick), FadeIn(num), run_time=0.8)
    s.say("Landau 的組合 ", "$r^{1/2}f(h/r^{3/2})$", " 放進 ", "$d=3$", " 就不成立", wait=2.4)
    s.say("讓尺度連續地變，重整化變成一組微分方程", wait=2.2)
    s.say("——參數的「流」", wait=2.0)


class W1105MomentumShell(ChuScene):
    def construct(self):
        self.show_header("W1105　動量殼與連續重整化")
        check_math()
        part_question(self)
        part_modes(self)
        part_shell(self)
        part_flow(self)
        part_pde(self)
        self.wait(1.0)
