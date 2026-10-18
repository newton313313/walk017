"""
W804 放大鏡：重新標度推出 Flory＋第 8 章回顧（§8.1.3，pp. 199–202，(8.21)–(8.25)，Ex. 8.6–8.8）

問：把漫步「縮小來看」——步數變少、長度變短——自迴避的效果會變大還是變小？
  1. (8.21) 動能＋自迴避；s′ = l₂s、r⃗′ = l₁r⃗（放大鏡）→ (8.22) 動能 × l₂/l₁²、自迴避 × l₁^d/l₂²；數值核對 0.8000、3.1250；筆誤框 1。
  2. Ex. 8.6–8.7：l₁ = √l₂ → l₂^{(d−4)/2}；N = 100：1000、100、10、1、0.1 → 上臨界維度 4；筆誤框 2；Ex. 8.8。
  3. (8.24)(8.25)：l₁ = l₂^{3/(d+2)} → ν = 1、3/4、3/5、1/2；d = 3、l₂ = 0.01：l₁ = 0.0631、1/l₁ = 15.85；拔河 vs 放大鏡。
  4. 第 8 章回顧（五列）。
畫面上的數字都由 utils.walk_scale 算出，check_math() 核對。模組層級的畫圖函式：magnifier（放大鏡示意）。
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
from utils import walk_scale as WS
from utils import walk_saw as SW
from utils import walk_draw as D

RC = (0.5, 0.2, 3)                    # 數值核對用的 (l1, l2, d)
DS = (1, 2, 3, 4, 5)
KIN_COLOR = BLUE_B
INT_COLOR = ORANGE


def f2(x):
    return f"{float(x):.2f}"


def f4(x):
    return f"{float(x):.4f}"


def ftex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return r"%s\frac{%d}{%d}" % ("-" if q < 0 else "", abs(q.numerator), q.denominator)


def fnum(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f"{float(q):g}"


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = {"rc": WS.rescale_check(*RC), "kin": WS.kin_factor(RC[0], RC[1]), "int": WS.int_factor(*RC),
              "amp": [WS.sa_amplitude(WS.SA_N, d) for d in DS],
              "l1": WS.flory_l1(0.01, 3), "l1n": WS.flory_l1_num(0.01, 3),
              "nu": [WS.nu_rescale(d) for d in range(1, 5)],
              "wl": WS.median_width(1000, 1), "wu": WS.median_width(1000, None),
              "fisher": [WS.fisher(Fraction(1, 2)), WS.fisher(Fraction(3, 4))],
              "g2d": WS.gamma_exp(WS.ETA_2D, WS.NU_2D)}
    return nb


def check_math():
    nb = numbers()
    h, q = Fraction(1, 2), Fraction(1, 4)
    assert WS.kin_factor(h, q) == 1 and WS.int_factor(h, q, 3) == 2
    assert f4(nb["kin"]) == "0.8000" and f4(nb["int"]) == "3.1250"
    assert abs(nb["rc"][0] - nb["kin"]) < 1e-9 and abs(nb["rc"][1] - nb["int"]) < 1e-9
    for d in DS:                                                    # l1 = √l2：自迴避 × l2^{(d−4)/2}；書上 l2^{(d−2)/2}
        l2 = 0.01
        assert abs(WS.kin_factor(math.sqrt(l2), l2) - 1) < 1e-12
        assert abs(WS.int_factor(math.sqrt(l2), l2, d) - l2 ** ((d - 4) / 2)) < 1e-9 * l2 ** ((d - 4) / 2)
        assert abs(WS.int_factor_book86(math.sqrt(l2), l2, d) - l2 ** ((d - 2) / 2)) < 1e-9 * l2 ** ((d - 2) / 2)
    assert nb["amp"] == [1000, 100, 10, 1, Fraction(1, 10)] and [fnum(v) for v in nb["amp"]] == ["1000", "100", "10", "1", "0.1"]
    assert WS.UPPER_CRITICAL_D == 4 and WS.sa_amplitude(100, 4) == 1
    assert WS.ex88_relevant(3, 3) and WS.ex88_relevant(5, 3) and WS.ex88_relevant(3, 5) and not WS.ex88_relevant(5, 5)
    assert abs(WS.ex88_factor(0.1, 0.01, 3) - 0.01 ** (-0.5)) < 1e-9
    assert f4(nb["l1"]) == "0.0631" and f2(1 / nb["l1"]) == "15.85" and f2(100 ** 0.6) == "15.85"
    assert abs(nb["l1"] - nb["l1n"]) < 1e-9
    assert nb["nu"] == [1, Fraction(3, 4), Fraction(3, 5), Fraction(1, 2)]
    assert all(WS.nu_rescale(d) == SW.flory_nu(d) for d in range(1, 5))
    assert f"{nb['wl']:.0f}" == "1822" and f2(nb["wu"]) == "28.06"
    assert nb["fisher"] == [2, 4] and nb["g2d"] == Fraction(43, 32)
    assert all(abs(WS.collapse_limit(x) - 2 / (4 + x * x)) < 1e-12 for x in WS.COLLAPSE_X)


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


def magnifier(width=4.0, center=ORIGIN):
    """放大鏡示意：一條細長的足跡（WS.trail(400, None) 的前兩欄）、鏡片圓圈＋把手、鏡片裡放大的一小段折線。
    回傳 VGroup(足跡, 圓圈, 把手, 放大的一段)。"""
    center = np.array(center, dtype=float)
    pts = WS.trail(400, None)[:, :2]
    mid = (pts.min(axis=0) + pts.max(axis=0)) / 2
    k = width / max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
    p3 = [center + np.array([(x - mid[0]) * k, (y - mid[1]) * k, 0.0]) for x, y in pts]
    path = VMobject().set_points_as_corners(p3).set_stroke(GREY_B, 1.2)
    seg = pts[200:212]
    lens_c = center + np.array([width * 0.32, -width * 0.30, 0])
    lens_r = width * 0.26
    sm = (seg.min(axis=0) + seg.max(axis=0)) / 2
    kk = 1.5 * lens_r / max(np.ptp(seg[:, 0]), np.ptp(seg[:, 1]), 1e-9)
    zoom = VMobject().set_points_as_corners([lens_c + np.array([(x - sm[0]) * kk, (y - sm[1]) * kk, 0]) for x, y in seg])
    zoom.set_stroke(YELLOW, 3)
    circ = Circle(radius=lens_r).move_to(lens_c).set_stroke(WHITE, 4).set_fill(BLACK, 0.85)
    d = np.array([1, -1, 0]) / math.sqrt(2)
    handle = Line(lens_c + lens_r * d, lens_c + (lens_r + width * 0.22) * d).set_stroke(GREY_A, 9)
    return VGroup(path, circ, handle, zoom)


# ---------------------------------------------------------------- 各段
def part_821(s):
    mg = magnifier(width=4.2)
    e21 = tx(r"\frac{3}{2a^2}\int\Big|\frac{d\vec r}{ds}\Big|^2ds\ +\ \frac{v_0}{2}\iint\delta\big(\vec r(s_1)-\vec r(s_2)\big)\,ds_1ds_2", 38)
    lab = VGroup(zh("動能（鏈的彈性）", font_size=28, color=KIN_COLOR), zh("自迴避", font_size=28, color=INT_COLOR)).arrange(RIGHT, buff=2.0)
    resc = tx(r"s'=l_2s,\qquad \vec r\,'=l_1\vec r", 38)
    e22 = VGroup(mixed("動能 ", "$\\times\\frac{l_2}{l_1^2}$", font_size=36, color=KIN_COLOR),
                 mixed("自迴避 ", "$\\times\\frac{l_1^d}{l_2^2}$", font_size=36, color=INT_COLOR)).arrange(RIGHT, buff=1.0)
    why = mixed("$\\delta(\\vec r/l_1)=l_1^d\\delta(\\vec r)$", "、", "$ds_1ds_2=ds_1'ds_2'/l_2^2$", font_size=32, color=GREY_B)
    kin, it = numbers()["kin"], numbers()["int"]
    rc = numbers()["rc"]
    num = mixed("數值核對（", "$l_1=%s,\\ l_2=%s,\\ d=%d$" % RC, "）：", "$%s$" % f4(rc[0]), "、", "$%s$" % f4(rc[1]),
                font_size=28, color=GREY_B)
    right = VGroup(e21, lab, resc, e22, why, num).arrange(DOWN, buff=0.3)
    place_in(mg, -6.7, -3.4, -2.4, 2.0)
    place_in(right, -3.1, 6.7, -2.9, 3.1)
    assert_in_frame(mg, right)
    assert_apart(mg, right)
    assert_column(*right)
    assert f4(kin) == f4(rc[0]) and f4(it) == f4(rc[1])
    s.say("問：把漫步「縮小來看」——步數變少、長度變短", wait=0.5)
    s.play(ShowCreation(mg[0]), run_time=2.0)
    s.say("自迴避的效果，會變大還是變小？", wait=1.0)
    s.play(FadeIn(mg[1]), ShowCreation(mg[2]), run_time=0.8)
    s.play(ShowCreation(mg[3]), run_time=1.0)
    s.play(Write(e21), FadeIn(lab), run_time=1.5)
    s.say("(8.21) 指數裡兩項：動能（彈性）與自迴避", wait=3.0)
    s.play(FadeIn(resc), run_time=0.8)
    s.say("重新標度：", "$s'=l_2s$", "、", "$\\vec r\\,'=l_1\\vec r$", "（", "$l_1,l_2<1$", "：縮小來看）", wait=3.0)
    s.play(FadeIn(e22), run_time=0.8)
    s.say("(8.22)：動能乘 ", "$l_2/l_1^2$", "，自迴避乘 ", "$l_1^d/l_2^2$", wait=3.0)
    s.play(FadeIn(why), run_time=0.8)
    s.say("$\\delta$", " 函數在 ", "$d$", " 維裡帶出 ", "$l_1^d$", "；兩個 ", "$ds$", " 帶出 ", "$1/l_2^2$", wait=3.0)
    s.play(FadeIn(num), run_time=0.8)
    s.say("數值核對：", "$l_2/l_1^2=0.8000$", "、", "$l_1^3/l_2^2=3.1250$", wait=3.0)
    s.clear_stage()


def part_ex86(s):
    nb = numbers()
    box1 = D.erratum_box("p. 200，(8.22) 前", r"\vec r\,'(s')=l_1\vec r(l_2s)", r"\vec r\,'(s')=l_1\vec r(s'/l_2)",
                         width=9.0, font_size=44)
    e1 = mixed("Ex. 8.6–8.7：令 ", "$l_1=\\sqrt{l_2}$", "（動能不變、", "$a'=a$", "）→ 自迴避 ", "$\\times l_2^{(d-4)/2}$",
               font_size=34)
    e2 = mixed("$l_2=\\frac1N$", "：只剩一步，", "$R\\sim a\\sqrt N$", "；自迴避被放大 ", "$N^{(4-d)/2}$", " 倍",
               font_size=34)
    rows = [["d"] + [str(d) for d in DS], [r"N=100"] + [fnum(v) for v in nb["amp"]]]
    tab = D.tex_table(rows, [1.6, 1.1, 1.1, 1.1, 1.1, 1.1], row_h=0.7, font_size=36)
    D.table_cell(tab, 1, 4, 6).set_color(YELLOW)
    D.table_cell(tab, 0, 4, 6).set_color(YELLOW)
    e3 = mixed("Ex. 8.8：", "$|\\vec r_1-\\vec r_2|^{-q}$", " → ", "$l_2^{(q-4)/2}$", "：", "$q$", " 取代 ", "$d$",
               "；重要 ", "$\\Leftrightarrow\\min(d,q)<4$", font_size=30)
    box2 = D.erratum_box("p. 200，Ex. 8.6", r"l_1^d/l_2", r"l_1^d/l_2^2", width=5.0, font_size=44)
    place_in(box1, -6.7, 6.7, 1.95, 3.15)
    col = VGroup(e1, e2, tab).arrange(DOWN, buff=0.3)
    place_in(col, -6.7, 6.7, -1.0, 1.75)
    place_in(box2, -6.7, -1.9, -2.9, -1.25)
    place_in(e3, -1.6, 6.7, -2.6, -1.4)
    assert_in_frame(box1, col, box2, e3)
    assert_column(box1, *col)
    assert_apart(col, box2)
    assert_apart(col, e3)
    assert_apart(box2, e3)
    assert_table(tab, 6)
    s.play(FadeIn(box1), run_time=0.8)
    s.say("書上寫 ", "$l_1\\vec r(l_2s)$", "；應該是 ", "$l_1\\vec r(s'/l_2)$", wait=2.5)
    s.say("新曲線在 ", "$s'$", " 的位置，是舊曲線在 ", "$s=s'/l_2$", " 的位置放大 ", "$l_1$", " 倍", wait=3.0)
    s.play(FadeIn(e1), run_time=0.8)
    s.say("令 ", "$l_1=\\sqrt{l_2}$", "：動能不變，自迴避乘 ", "$l_2^{(d-4)/2}$", wait=3.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("$l_2=1/N$", "：整條鏈縮成一步，", "$R\\sim a\\sqrt N$", wait=3.0)
    s.play(FadeIn(tab[0][:6]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][6:]], lag_ratio=0.15), run_time=1.0)
    s.say("$N=100$", "：1000、100、10、1、0.1", wait=3.0)
    s.say("$d<4$", " 越縮越重要；", "$d=4$", " 是分界：上臨界維度 4", wait=3.0)
    s.say("呼應 W705 的冪次表：", "$d=4$", " 那列全是 1", wait=3.0)
    s.play(FadeIn(box2), run_time=0.8)
    s.say("Ex. 8.6 書上寫 ", "$l_1^d/l_2$", "；應該是 ", "$l_1^d/l_2^2$", wait=2.5)
    s.say("照書上得 ", "$l_2^{(d-2)/2}$", "，上臨界維度會變成 2", wait=3.0)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("Ex. 8.8：冪次位能 ", "$|\\vec r_1-\\vec r_2|^{-q}$", "：", "$q$", " 取代了 ", "$d$", wait=3.0)
    s.clear_stage()


def part_flory(s):
    nb = numbers()
    e1 = tx(r"\frac{l_2}{l_1^2}=\frac{l_1^d}{l_2^2}\ \Rightarrow\ l_1=l_2^{3/(d+2)}\ \Rightarrow\ "
            r"R\sim l_1^{-1}=N^{3/(d+2)}", 42)
    rows = [["d", "1", "2", "3", "4"], [r"\nu=\frac{3}{d+2}"] + [ftex(v) for v in nb["nu"]]]
    tab = D.tex_table(rows, [2.0, 1.1, 1.1, 1.1, 1.1], row_h=0.85, font_size=36)
    ex = mixed("例 ", "$d=3,\\ l_2=0.01$", "：", "$l_1=%s$" % f4(nb["l1"]), "，", "$1/l_1=%s=100^{3/5}$" % f2(1 / nb["l1"]),
               font_size=34)
    top = VGroup(e1, tab, ex).arrange(DOWN, buff=0.3)
    place_in(top, -6.7, 6.7, 0.0, 3.1)
    b1 = VGroup(zh("拔河（W706）", font_size=32, color=YELLOW),
                mixed("彈性 ", "$R^2/N$", " 對 排斥 ", "$N^2/R^d$", font_size=28),
                mixed("取最小 → ", "$R\\sim N^{3/(d+2)}$", font_size=28)).arrange(DOWN, buff=0.2)
    b2 = VGroup(zh("放大鏡（W804）", font_size=32, color=YELLOW),
                mixed("動能 ", "$\\frac{l_2}{l_1^2}$", " 和 自迴避 ", "$\\frac{l_1^d}{l_2^2}$", " 同比例", font_size=28),
                mixed("→ ", "$R\\sim N^{3/(d+2)}$", font_size=28)).arrange(DOWN, buff=0.2)
    f1 = SurroundingRectangle(b1, buff=0.25).set_stroke(GREY_B, 2)
    f2_ = SurroundingRectangle(b2, buff=0.25).set_stroke(GREY_B, 2)
    pair = VGroup(VGroup(f1, b1), VGroup(f2_, b2)).arrange(RIGHT, buff=0.8)
    place_in(pair, -6.7, 6.7, -2.9, -0.35)
    assert_in_frame(top, pair)
    assert_column(*top, pair)
    assert_table(tab, 5)
    s.say("(8.24)：讓兩項用同一個比例縮放", wait=0.3)
    s.play(Write(e1), run_time=2.0)
    s.say("$l_2/l_1^2=l_1^d/l_2^2$", "：", "$l_1=l_2^{3/(d+2)}$", wait=3.0)
    s.say("(8.25)：", "$R\\sim l_1^{-1}=N^{3/(d+2)}$", wait=3.0)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(c) for c in tab[0][5:]], lag_ratio=0.15), run_time=1.0)
    s.say("$\\nu=1,\\ \\frac34,\\ \\frac35,\\ \\frac12$", "（", "$d=1..4$", "）", wait=3.0)
    s.play(FadeIn(ex), run_time=0.8)
    s.say("$d=3$", "：", "$l_2=0.01$", " 時 ", "$l_1=0.0631$", "，", "$1/l_1=15.85=100^{3/5}$", wait=3.0)
    s.play(FadeIn(pair[0]), run_time=0.8)
    s.play(FadeIn(pair[1]), run_time=0.8)
    s.say("和 W706 的拔河是同一個答案：Flory 指數", wait=3.0)
    s.clear_stage()


def part_recap(s):
    nb = numbers()
    fz = 24
    rows = [[zh("幕", font_size=26, color=YELLOW), zh("第 8 章重點", font_size=26, color=YELLOW)],
            [tx("W801", fz), mixed("長尾巴寬度 ", "$\\propto N^{1/y}$", "（", "$%s$" % f"{nb['wl']:.0f}", " 對 ",
                                   "$%s$" % f2(nb["wu"]), "）、尾巴 ", "$r^{-5+\\delta}$", font_size=fz)],
            [tx("W802", fz), mixed("塌縮 ", "$\\to\\frac{2}{4+x^2}$", "、", "$\\xi\\propto(z_c-z)^{-1/2}$", font_size=fz)],
            [tx("W803", fz), mixed("$e^{-(r/N^\\nu)^{1/(1-\\nu)}}$", "、Fisher ", "$%s$" % ftex(nb["fisher"][0]), " 與 ",
                                   "$%s$" % ftex(nb["fisher"][1]), "、", "$\\gamma=(2-\\eta)\\nu=%s$" % ftex(nb["g2d"]),
                                   font_size=fz)],
            [tx("W804", fz), mixed("$l_2^{(d-4)/2}$", "、上臨界維度 ", "$%d$" % WS.UPPER_CRITICAL_D, font_size=fz)],
            [tx(r"\text{W804}=\text{W706}", fz), mixed("$\\nu=\\frac{3}{d+2}$", font_size=fz)]]
    tab = D.tex_table(rows, [2.8, 8.6], row_h=0.78, font_size=fz)
    place_in(tab, -6.7, 6.7, -2.6, 2.9)
    assert_in_frame(tab)
    assert_table(tab, 2)
    s.say("第 8 章回顧", wait=0.3)
    s.play(FadeIn(tab[0][:2]), ShowCreation(tab[1]), run_time=0.8)
    for i in range(1, 6):
        s.play(FadeIn(tab[0][2 * i]), FadeIn(tab[0][2 * i + 1]), run_time=0.7)
        s.wait(1.3)
    s.say("長尾巴：寬度 ", "$N^{1/y}$", "，尾巴是冪次——不是高斯", wait=3.0)
    s.say("標度形式：用 ", "$\\xi$", " 當尺，所有曲線疊成一條", wait=3.0)
    s.say("放大鏡底下，四維以下自迴避越看越重要", wait=3.0)
    s.say("讓兩項一起縮放，就得到 Flory 的 ", "$\\frac{3}{d+2}$", wait=3.0)


class W804RescalingFlory(ChuScene):
    def construct(self):
        self.show_header("W804　放大鏡：重新標度推出 Flory")
        check_math()
        part_821(self)
        part_ex86(self)
        part_flory(self)
        part_recap(self)
        self.wait(1.0)
