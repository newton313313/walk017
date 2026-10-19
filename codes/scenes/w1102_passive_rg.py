"""
W1102 被動的重整化：Landau 與 Gaussian（§11.1.2，pp. 261–264，(11.7)–(11.25)，Ex. 11.4）

問：把整個磁鐵放大 b 倍，自由能的每一項會怎樣？
  1. (11.9) ℋ = V[rm² + um⁴ − hm]，V → b^d V ⇒ m' = m b^{d/4}、r' = r b^{d/2}、h' = h b^{3d/4}。
  2. (11.11)–(11.16) m = r^{1/2} f(h/r^{3/2})；m(−0.25, 0) = 0.5000、m(0, h) = h^{1/3}（0.2000、0.5000）；
     筆誤框 1（p. 263，(11.16) 後：h^{1/3} r → h^{1/3} r^{−1/2}）：r = 1e−4、h = 1e−3：f = 9.9667、更正 10.0000、書上 10^{−5}。
  3. 資料塌縮：m 對 h（r = 0.25、0.5、1）→ m/r^{1/2} 對 h/r^{3/2} 疊成 y + y³ = x；f(1)、f(2)、f(10)；(11.17) γ = 1。
  4. Gaussian (11.18)–(11.25)：m = r^{(d−2)/4} F(h/r^{(d+2)/4})；Ex. 11.4 表；d = 4 時與 Landau 的組合相同。
畫面上的數字都由 utils.walk_rg 算出，check_math() 核對。
"""
import sys
import warnings
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_rg as RG
from utils import walk_draw as D

M_COLOR, R_COLOR, H_COLOR = BLUE_B, GREEN_B, ORANGE
R_COLORS = (BLUE_B, GREEN_B, ORANGE)
H_LIM = (0.008, 0.125)
F_X = (1, 2, 10)
F_VAL = ("0.6823", "1.0000", "2.0000")
PT_M = ("0.5000", "0.7071", "1.0000")


def f4(x):
    return f"{float(x):.4f}"


def ftex(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else rf"\frac{{{q.numerator}}}{{{q.denominator}}}"


def fslash(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def mf_exponents():
    """Landau 的組合 r^{1/2} f(h/r^{3/2})：β = 1/2、γ = 3/2 − 1/2、δ = (3/2)/(1/2)。"""
    a, _, c, _ = RG.scaling_combos(4)
    return {"beta": c, "gamma": a - c, "delta": a / c}


def check_math():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert f4(RG.landau_m(-0.25, 0)) == "0.5000" == f4(np.sqrt(0.25))
        assert [f4(RG.landau_m(0, h)) for h in H_LIM] == ["0.2000", "0.5000"]
        assert all(abs(RG.landau_m(0, h) - h ** (1 / 3)) < 1e-9 for h in H_LIM)
        f, book, right = RG.book1116_limit(1e-4, 1e-3)
        assert f4(f) == "9.9667" and f4(right) == "10.0000" and abs(book - 1e-5) < 1e-15
        assert [f4(RG.landau_scaled(x)) for x in F_X] == list(F_VAL)
        assert RG.COLLAPSE_R == (0.25, 0.5, 1.0)
        for (r, h), m in zip(RG.COLLAPSE_POINTS, PT_M):
            assert f4(RG.landau_m(r, h)) == m
            assert abs(h / r ** 1.5 - 2) < 1e-12 and abs(RG.landau_m(r, h) / r ** 0.5 - 1) < 1e-9
        assert f4(RG.COLLAPSE_POINTS[1][1]) == "0.7071"
        for r in RG.COLLAPSE_R:
            lc = RG.landau_curve(r)
            assert len(lc) <= 40 and max(m for _, m in lc) < 1.5 and abs(lc[-1][0] - 2) < 1e-12
            cc = [(x, y) for x, y in RG.collapse_curve(r) if x <= 4]
            assert all(abs(y + y ** 3 - x) < 1e-8 for x, y in cc) and max(y for _, y in cc) < 1.5
        g3, g4 = RG.gauss_exponents(3), RG.gauss_exponents(4)
        assert (g3["beta"], g3["gamma"], g3["delta"]) == (Fraction(1, 4), 1, 5)
        assert (g4["beta"], g4["gamma"], g4["delta"]) == (Fraction(1, 2), 1, 3)
        mf = mf_exponents()
        assert (mf["beta"], mf["gamma"], mf["delta"]) == (g4["beta"], g4["gamma"], g4["delta"])
        c4, c3 = RG.scaling_combos(4), RG.scaling_combos(3)
        assert c4 == (Fraction(3, 2), Fraction(3, 2), Fraction(1, 2), Fraction(1, 2))
        assert (c3[1], c3[3]) == (Fraction(5, 4), Fraction(1, 4))
        for d in (3, 4):
            g = RG.gauss_exponents(d)
            assert g["beta"] == Fraction(d - 2, 4) and g["delta"] == Fraction(d + 2, d - 2) and g["gamma"] == 1


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_left(m, x0, y, x1=6.7):
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


def plot_line(ax, pts, color, width=3):
    P = [ax.c2p(x, y) for x, y in pts]
    return VMobject().set_points_as_corners(P).set_stroke(color, width)


# ---------------------------------------------------------------- 1. 放大整個磁鐵
def part_scale(s):
    h0 = tx(r"\mathcal H=V\big[rm^2+um^4-hm\big]", 40)
    h0.move_to([0, 2.4, 0])
    tag = tx("(11.9)", 26, GREY_B).next_to(h0, LEFT, buff=0.5)
    v = tx(r"V\to b^dV", 36, YELLOW).move_to([0, 1.4, 0])
    lines = VGroup(tx(r"m\ \longrightarrow\ m'=m\,b^{d/4}", 38, M_COLOR),
                   tx(r"r\ \longrightarrow\ r'=r\,b^{d/2}", 38, R_COLOR),
                   tx(r"h\ \longrightarrow\ h'=h\,b^{3d/4}", 38, H_COLOR))
    for m, yy in zip(lines, (0.35, -0.6, -1.55)):
        m.move_to([0, yy, 0])
    notes = VGroup(*[mixed("由 ", f"${t}$", " 項決定", font_size=24, color=GREY_B) for t in ("um^4", "rm^2", "hm")])
    for n, m in zip(notes, lines):
        n.next_to(m, RIGHT, buff=0.6)
    note9 = mixed("本系列：", r"$\frac r2m^2+\frac u4m^4-hm$", font_size=24, color=GREY_B)
    place_left(note9, -6.3, -2.5)
    assert_in_frame(h0, tag, v, lines, notes, note9)
    for a, b in ((h0, v), (v, lines[0]), (lines[0], lines[1]), (lines[1], lines[2]), (lines[2], note9), (h0, tag)):
        assert_apart(a, b)

    s.say("把整個磁鐵放大 ", "$b$", " 倍，自由能的每一項會怎樣？", wait=2.0)
    s.play(FadeIn(h0), FadeIn(tag), run_time=0.8)
    s.say("Landau 自由能：體積 ", "$V$", " 乘上三項", wait=2.0)
    s.play(FadeIn(v), run_time=0.6)
    s.say("放大 ", "$b$", " 倍：", "$d$", " 維的體積變成 ", "$b^dV$", wait=2.0)
    s.say("要讓 ", r"$\mathcal H$", " 的樣子不變，每個參數都要跟著改", wait=2.0)
    for i in range(3):
        s.play(Write(lines[i]), FadeIn(notes[i]), run_time=0.9)
    s.say("$m$", "、", "$r$", "、", "$h$", " 各自乘上 ", "$b$", " 的某個次方", wait=2.5)
    s.play(FadeIn(note9), run_time=0.6)
    s.say("本系列照第 9 章寫 ", r"$\frac r2m^2+\frac u4m^4$", "——常數不改變標度", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 2. 標度形式與極限
def part_form(s):
    e1 = place_left(mixed("(11.11)  ", r"$m(r,h)=b^{-d/4}F\big(rb^{d/2},\ hb^{3d/4}\big)$", font_size=30), -6.3, 2.55)
    e2 = place_left(mixed("與 ", "$b$", " 無關：取 ", r"$rb^{d/2}=1$", font_size=28), -6.3, 1.7)
    e3 = place_left(tx(r"m=r^{1/2}f\big(h/r^{3/2}\big)", 40, YELLOW), -6.3, 0.75)
    box3 = SurroundingRectangle(e3, buff=0.15).set_stroke(YELLOW, 2)
    l1 = place_left(mixed("$h=0$", "：", r"$m=\sqrt{-r/u}$", "，", f"$m(-0.25,0)={f4(RG.landau_m(-0.25, 0))}$",
                          font_size=26), -6.3, -0.4, 0.0)
    l2 = place_left(mixed("$r=0$", "：", r"$m=h^{1/3}$", font_size=26), -6.3, -1.25, 0.0)
    l3 = place_left(tx(rf"h={H_LIM[0]}\to{f4(RG.landau_m(0, H_LIM[0]))},\quad h={H_LIM[1]}\to{f4(RG.landau_m(0, H_LIM[1]))}",
                       30), -5.8, -2.05, 0.0)
    err = D.erratum_box("p. 263，(11.16) 後", r"f(h/r^{3/2})\to h^{1/3}\,r", r"f(h/r^{3/2})\to h^{1/3}\,r^{-1/2}",
                        width=6.4)
    err.move_to([3.5, -1.2, 0])
    assert_in_frame(e1, e2, e3, box3, l1, l2, l3, err)
    for a, b in ((e1, e2), (e2, box3), (box3, l1), (l1, l2), (l2, l3), (err, l1), (err, l2), (err, l3), (err, box3)):
        assert_apart(a, b)

    s.say("答案與 ", "$b$", " 無關：可以挑一個方便的 ", "$b$", wait=2.0)
    s.play(FadeIn(e1), run_time=0.8)
    s.say("放大前後是同一個磁鐵，所以 ", "$m$", " 不變", wait=2.0)
    s.play(FadeIn(e2), run_time=0.6)
    s.say("挑 ", "$b$", " 讓第一個參數變成 1", wait=2.2)
    s.play(FadeIn(e3), ShowCreation(box3), run_time=0.9)
    s.say("兩個變數的函數，變成一個變數的函數", wait=2.0)
    s.play(FadeIn(l1), run_time=0.6)
    s.say("沒有外場、", "$r<0$", "：自發磁化 ", r"$\sqrt{-r/u}$", wait=2.0)
    s.play(FadeIn(l2), FadeIn(l3), run_time=0.6)
    s.say("臨界點上 ", "$r=0$", "：", r"$m=h^{1/3}$", "，也就是 ", r"$\delta=3$", wait=2.0)
    s.play(FadeIn(err), run_time=0.8)
    s.say("書上寫成 ", r"$h^{1/3}r$", "：", "$r$", " 的次方錯了", wait=2.0)
    f, book, right = RG.book1116_limit(1e-4, 1e-3)
    s.say(r"$r=10^{-4}$", "、", r"$h=10^{-3}$", "：", f"$f={f4(f)}$", "，更正式 ", f"${f4(right)}$", wait=2.0)
    s.say("書上的只有 ", r"$10^{-5}$", "：差了一百萬倍", wait=2.0)
    s.clear_stage()


# ---------------------------------------------------------------- 3. 資料塌縮
def part_collapse(s):
    axL = D.plot_axes((0, 2, 0.5), (0, 1.5, 0.5), width=5.0, height=3.2, x_ticks=(1, 2), y_ticks=(0.5, 1),
                      x_label="h", y_label="m")
    axR = D.plot_axes((0, 4, 1), (0, 1.5, 0.5), width=5.0, height=3.2, x_ticks=(1, 2, 3, 4), y_ticks=(0.5, 1),
                      x_label=r"h/r^{3/2}", y_label=r"m/r^{1/2}")
    axL.move_to([-3.4, 1.0, 0])
    axR.move_to([3.5, 1.0, 0])
    curvesL = VGroup(*[plot_line(axL, RG.landau_curve(r), c) for r, c in zip(RG.COLLAPSE_R, R_COLORS)])
    curvesR = VGroup(*[plot_line(axR, [(x, y) for x, y in RG.collapse_curve(r) if x <= 4], c)
                       for r, c in zip(RG.COLLAPSE_R, R_COLORS)])
    dotsL = VGroup(*[Dot(axL.c2p(h, RG.landau_m(r, h)), radius=0.07).set_color(c)
                     for (r, h), c in zip(RG.COLLAPSE_POINTS, R_COLORS)])
    dotsR = VGroup(*[Dot(axR.c2p(h / r ** 1.5, RG.landau_m(r, h) / r ** 0.5), radius=0.07).set_color(c)
                     for (r, h), c in zip(RG.COLLAPSE_POINTS, R_COLORS)])
    leg = VGroup(*[tx(f"r={r:g}", 26, c) for r, c in zip(RG.COLLAPSE_R, R_COLORS)]).arrange(RIGHT, buff=0.5)
    leg.move_to([-3.4, -1.38, 0])
    pts = place_left(mixed("三個點：", "$m=" + ",\\ ".join(PT_M) + "$", font_size=26), -6.3, -2.0, 1.1)
    yx = tx(r"y+y^3=x", 30, YELLOW).move_to([3.5, -1.38, 0])
    fv = place_left(tx(",\\quad ".join(f"f({x})={F_VAL[k] if k == 0 else int(float(F_VAL[k]))}" for k, x in enumerate(F_X)),
                       28), 1.4, -2.0)
    gam = place_left(mixed("(11.17)  ", r"$\frac{\partial m}{\partial h}=r^{-1}f'(0)$", r"$\Rightarrow$", r"$\gamma=1$",
                           font_size=24), -2.0, -2.63)
    assert_in_frame(axL, axR, leg, pts, yx, fv, gam)
    for a, b in ((axL, axR), (axL, leg), (axR, yx), (leg, pts), (yx, fv), (pts, gam), (fv, gam), (leg, yx)):
        assert_apart(a, b)

    s.say("數值上看看：三個 ", "$r$", "，畫 ", "$m$", " 對 ", "$h$", wait=1.0)
    s.play(FadeIn(axL), FadeIn(leg), run_time=0.8)
    s.play(*[ShowCreation(c) for c in curvesL], run_time=1.5)
    s.say("三條不同的曲線", wait=2.0)
    s.play(FadeIn(axR), run_time=0.8)
    s.say("改用 ", r"$m/r^{1/2}$", " 對 ", r"$h/r^{3/2}$", wait=2.0)
    s.play(*[TransformFromCopy(a, b) for a, b in zip(curvesL, curvesR)], run_time=2.0)
    s.play(FadeIn(yx), run_time=0.5)
    s.say("三條疊成一條：", r"$y+y^3=x$", wait=2.0)
    s.play(FadeIn(fv), run_time=0.6)
    s.say(f"$f(1)={F_VAL[0]}$", "、", "$f(2)=1$", "、", "$f(10)=2$", wait=2.0)
    s.play(FadeIn(dotsL), FadeIn(pts), run_time=0.6)
    s.say("三個點 ", r"$(r,h)=(0.25,0.25),\ (0.5,0.7071),\ (1,2)$", wait=2.5)
    s.say("左邊三個點：", "$m$", " 各不相同", wait=2.0)
    s.play(TransformFromCopy(dotsL, dotsR), run_time=1.2)
    s.say("右邊：全部落在 ", "$x=2$", "、", "$y=1$", wait=2.0)
    s.play(FadeIn(gam), run_time=0.6)
    s.say("對 ", "$h$", " 微分：", r"$\chi\propto r^{-1}$", "，所以 ", r"$\gamma=1$", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 4. Gaussian 模型
def part_gauss(s):
    x0 = -6.3
    g1 = place_left(tx(r"d^dx\to b^d\,d^dx,\qquad \nabla\to b^{-1}\nabla", 30), x0, 2.6, 0.9)
    g2 = place_left(tx(r"S(b\vec x)=b^{(2-d)/2}S'", 30), x0, 1.8, 0.9)
    g3 = place_left(tx(r"r=b^{-2}r',\qquad h=b^{-(d+2)/2}h'", 30), x0, 1.0, 0.9)
    g4 = place_left(tx(r"m=r^{(d-2)/4}F\Big(h/r^{(d+2)/4}\Big)", 36, YELLOW), x0, 0.05, 0.9)
    g5 = place_left(tx(r"\beta=\frac{d-2}4,\quad \gamma=1,\quad \delta=\frac{d+2}{d-2}", 30), x0, -0.95, 0.9)
    g3e, g4e, mf = RG.gauss_exponents(3), RG.gauss_exponents(4), mf_exponents()
    rows = [[zh("Ex. 11.4", font_size=22, color=YELLOW), r"\beta", r"\gamma", r"\delta"],
            ["d=3"] + [fslash(g3e[k]) for k in ("beta", "gamma", "delta")],
            ["d=4"] + [fslash(g4e[k]) for k in ("beta", "gamma", "delta")],
            [zh("平均場", font_size=24)] + [fslash(mf[k]) for k in ("beta", "gamma", "delta")]]
    tab = D.tex_table(rows, [1.8, 1.2, 1.2, 1.2], row_h=0.66, font_size=30).shift(np.array([1.3, 2.9, 0]))
    hl = SurroundingRectangle(VGroup(*[D.table_cell(tab, i, j, 4) for i in (2, 3) for j in range(4)]), buff=0.06)
    hl.set_stroke(GREEN_B, 2)
    c4, c3 = RG.scaling_combos(4), RG.scaling_combos(3)
    k1 = place_left(mixed("$d=4$", "：", rf"$\frac{{d+2}}4={ftex(c4[1])}$", "、", rf"$\frac{{d-2}}4={ftex(c4[3])}$",
                          "＝Landau 的 ", rf"${ftex(c4[0])}$", "、", rf"${ftex(c4[2])}$", font_size=26), x0, -1.85)
    k2 = place_left(mixed("$d=3$", "：", rf"${ftex(c3[1])}$", "、", rf"${ftex(c3[3])}$", "（和 Landau 不同）", font_size=26),
                    x0, -2.6)
    assert_in_frame(g1, g2, g3, g4, g5, tab, hl, k1, k2)
    assert_table(tab, 4)
    for a, b in ((g1, g2), (g2, g3), (g3, g4), (g4, g5), (g5, k1), (k1, k2), (tab, g1), (tab, g3), (hl, g4), (hl, g5)):
        assert_apart(a, b)

    s.say("Gaussian 模型：只有 ", "$m^2$", " 和梯度項", wait=2.0)
    s.play(FadeIn(g1), run_time=0.7)
    s.say("放大 ", "$b$", " 倍：體積元變大、梯度變小", wait=2.0)
    s.play(FadeIn(g2), run_time=0.7)
    s.say("梯度項不變：自旋要乘上 ", r"$b^{(2-d)/2}$", wait=2.0)
    s.play(FadeIn(g3), run_time=0.7)
    s.say("$r$", " 和 ", "$h$", " 也各自重新標度", wait=2.2)
    s.play(FadeIn(g4), run_time=0.8)
    s.say("同樣的論證：", "$m$", " 又是一個變數的函數", wait=2.0)
    s.play(FadeIn(g5), run_time=0.7)
    s.say("從 ", "$m$", " 的形式直接讀出 ", r"$\beta$", "、", r"$\gamma$", "、", r"$\delta$", wait=2.0)
    s.play(FadeIn(tab), run_time=1.0)
    s.say("Ex. 11.4：", "$d=3$", " 時 ", r"$\beta=\frac14$", "、", r"$\delta=5$", wait=2.0)
    s.play(ShowCreation(hl), run_time=0.7)
    s.say("$d=4$", "：", r"$\frac12$", "、1、3，和平均場一樣", wait=2.0)
    s.play(FadeIn(k1), run_time=0.7)
    s.say("$d=4$", " 時兩種組合完全相同", wait=2.2)
    s.play(FadeIn(k2), run_time=0.7)
    s.say("第 8 章說過：", "$d=4$", " 是特殊的維度", wait=2.0)
    s.say("只看「放大之後參數怎麼變」，就讀出了 ", r"$\gamma=1$", wait=2.5)
    s.say("但這還是「被動」的：漲落沒有被積掉", wait=2.5)


class W1102PassiveRG(ChuScene):
    def construct(self):
        self.show_header("W1102　被動的重整化：Landau 與 Gaussian")
        check_math()
        part_scale(self)
        part_form(self)
        part_collapse(self)
        part_gauss(self)
        self.wait(1.0)
