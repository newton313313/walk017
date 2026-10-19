"""
W1201 加上四次項：怎樣正確地積掉自由度（§12.1–12.2.4 前半，pp. 285–290，(12.1)–(12.11)、(12.20)–(12.34)，Fig. 12.1–12.5）

問：加上四次項 u S⁴ 之後，還能「積掉一半」嗎？
  1. (12.3) 一分量 GLW；Fig. 12.1–12.2：正方晶格與倒晶格並排（無數字）。
  2. Fig. 12.3–12.4：半徑 Λ = e^{−ℓ} 的布里淵區與殼層；(12.24) dq；殼層比例表 (12.20)–(12.23)；格點計數。
  3. (12.5)–(12.10) 流程圖；兩變數玩具：對 S' 積分得 ℋ'(s)，展開出 12、36 與 s⁶ 項；Z = A·Z'。
  4. (12.11) 五項 1, 4, 6, 4, 1；重新標度 (12.29)–(12.34) 的乘數；筆誤框 2 個（(12.32)、(12.33)）。
畫面上的數字都由 utils.walk_eps 算出，check_math() 核對。
模組層級的畫圖函式：lattice_panel、bz_disk、flow_box、toy_axes、toy_curve、vertex_icon。
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
from utils import walk_eps as E
from utils import lattice_draw as L
from utils import walk_draw as D

SHELL_COLOR = TEAL
INNER_COLOR = WHITE
ELL, DELL = 0.5, 0.01
TOY_N = 40


# ---------------------------------------------------------------- 數學
def f4(x):
    return f"{float(x):.4f}"


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r, u = E.TOY
        ss = np.linspace(-2, 2, TOY_N)
        rp, up, c6 = E.toy_coeffs(r, u)
        return {"lam": math.exp(-ELL),
                "dq": E.shell_dq(ELL, DELL),
                "shell": [(d, dl, E.shell_fraction(d, dl), E.shell_linear(d, dl)) for d, dl in E.SHELL_CASES],
                "lattice": [abs(E.shell_fraction_lattice(3, dl) - E.shell_fraction(3, dl)) for dl in (0.1, 0.01)],
                "toy": (rp, up, c6),
                "Z": E.toy_Z(r, u),
                "A": E.toy_A(r),
                "heff": [(s, E.toy_Heff(s, r, u)) for s in ss],
                "quart": [(s, rp / 2 * s * s + up * s ** 4) for s in ss],
                "resc": E.rescale_factors(3, 0.1),
                "book": E.rescale_factors(3, 0.1, book=True)}


def check_math():
    nb = numbers()
    assert f4(nb["lam"]) == "0.6065"
    dq, lin = nb["dq"]
    assert f"{dq:.6f}" == "0.006035" and f"{lin:.6f}" == "0.006065" and abs(lin - nb["lam"] * DELL) < 1e-15
    assert E.SHELL_CASES == ((3, 0.1), (3, 0.01), (4, 0.01))
    assert [(f4(a), f4(b)) for _, _, a, b in nb["shell"]] == [("0.2592", "0.3000"), ("0.0296", "0.0300"),
                                                            ("0.0392", "0.0400")]
    assert all(x < 2e-4 for x in nb["lattice"])
    assert E.TOY == (0.2, 0.01)
    rp, up, c6 = nb["toy"]
    assert f4(rp) == "0.3000" and f4(up) == "0.0075" and f"{c6:.6f}" == "0.000167"
    assert abs(rp - (0.2 + 12 * 0.01 / 1.2)) < 1e-12 and abs(up - (0.01 - 36e-4 / 1.44)) < 1e-12
    z1, z2 = nb["Z"]
    assert abs(z1 - z2) < 1e-6 and f4(z1) == "9.0517" and f4(nb["A"]) == "2.2882" and f4(z1 / nb["A"]) == "3.9558"
    assert len(nb["heff"]) <= 40 and max(abs(a[1] - b[1]) for a, b in zip(nb["heff"], nb["quart"])) < 0.02
    assert max(h for _, h in nb["heff"]) < 0.8
    assert [math.comb(4, k) for k in range(5)] == [1, 4, 6, 4, 1]
    assert [f4(x) for x in nb["resc"]] == ["1.0000", "1.2214", "1.1052"]
    assert f4(nb["book"][0]) == "0.6703"
    assert f4(math.exp(3 * 0.01)) == "1.0305" and f4(1 + 0.01) == "1.0100"


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def place_in(m, x0, x1, y0, y1):
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def assert_column(*mobs):
    for a, b in zip(mobs, mobs[1:]):
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


# ---------------------------------------------------------------- 畫圖元件
def lattice_panel(cell, origin, label_tex, title):
    """5×5 格點＋格距標示＋中文標題。"""
    o = np.array(origin, dtype=float)
    dots = L.lattice_dots(2, cell, o, 0.06)
    p0, p1 = L.lat_point(-2, -2, cell, o), L.lat_point(-1, -2, cell, o)
    brace = Line(p0 + 0.3 * DOWN, p1 + 0.3 * DOWN).set_stroke(YELLOW, 3)
    lab = tx(label_tex, 26, YELLOW).next_to(brace, DOWN, buff=0.1)
    ttl = zh(title, font_size=26, color=GREY_A).next_to(dots, UP, buff=0.3)
    return VGroup(dots, brace, lab, ttl)


def bz_disk(radius, origin, shell_w):
    """布里淵區（圓）＋外圈殼層（環）＋半徑標示。"""
    o = np.array(origin, dtype=float)
    disk = Circle(radius=radius).move_to(o).set_stroke(YELLOW, 3)
    ring = Annulus(inner_radius=radius - shell_w, outer_radius=radius).move_to(o)
    ring.set_fill(SHELL_COLOR, 0.6).set_stroke(width=0)
    rad = Line(o, o + radius * np.array([math.cos(PI / 5), math.sin(PI / 5), 0])).set_stroke(WHITE, 2)
    rlab = tx(r"\Lambda=e^{-\ell}", 26).next_to(rad.get_center(), UL, buff=0.05)
    return VGroup(disk, ring, rad, rlab)


def flow_box(parts, tag, width=2.8, color=BLUE_C):
    body = mixed(*parts, font_size=22)
    if body.get_width() > width - 0.25:
        body.set_width(width - 0.25)
    t = zh(tag, font_size=18, color=GREY_B)
    inner = VGroup(body, t).arrange(DOWN, buff=0.15)
    box = RoundedRectangle(width=width, height=1.25, corner_radius=0.1).set_stroke(color, 2).move_to(inner)
    return VGroup(box, inner)


def toy_axes():
    ax = D.plot_axes((-2.2, 2.3, 1), (0, 0.85, 0.2), width=5.2, height=3.0, x_ticks=(-2, -1, 1, 2),
                     x_tick_labels=["-2", "-1", "1", "2"], x_label="s")
    ax.add(VGroup(*[tx(s, 20, GREY_B).next_to(ax.c2p(0, v), RIGHT, buff=0.1) for s, v in (("0.4", 0.4),
                                                                                           ("0.8", 0.8))]))
    return ax


def toy_curve(ax, pts, color, dashed=False):
    c = VMobject().set_points_smoothly([ax.c2p(x, y) for x, y in pts]).set_stroke(color, 4)
    if dashed:
        c = DashedVMobject(c, num_dashes=30).set_stroke(color, 4)
    return c


def vertex_icon(k_shell, origin, arm=0.45):
    """一個四次頂點：點＋四隻手（前 k 隻在殼層，青色）。"""
    o = np.array(origin, dtype=float)
    g = VGroup()
    for j in range(4):
        ang = PI / 4 + j * PI / 2
        col = SHELL_COLOR if j < k_shell else INNER_COLOR
        g.add(Line(o, o + arm * np.array([math.cos(ang), math.sin(ang), 0])).set_stroke(col, 4))
    g.add(Dot(o, radius=0.08).set_color(YELLOW))
    return g


# ---------------------------------------------------------------- 各段
def part_question(s):
    ham = tx(r"\mathcal H=\sum_{\vec q}\frac12(r+q^2)|S(\vec q)|^2+\frac uN\sum S^4\,\delta", 34)
    tag = zh("(12.3)　一分量 GLW 模型", font_size=22, color=GREY_B)
    top = VGroup(tag, ham).arrange(DOWN, buff=0.2).move_to(np.array([0, 2.45, 0]))
    real = lattice_panel(0.45, np.array([-3.2, -0.6, 0]), "a", "實空間（Fig. 12.1）")
    recip = lattice_panel(0.65, np.array([3.0, -0.6, 0]), r"2\pi/a", "倒晶格（Fig. 12.2）")
    assert_in_frame(top, real, recip)
    assert_apart(top, real)
    assert_apart(top, recip)
    assert_apart(real, recip)
    s.say("問：上一章只有二次項；加上四次項，還能「積掉一半」嗎？", wait=2.2)
    s.play(FadeIn(tag), Write(ham), run_time=1.4)
    s.say("GLW 模型：二次項 ", "$\\frac12(r+q^2)$", "，加上四次項 ", "$u\\,S^4$", wait=2.2)
    s.say("四次項把不同的 ", "$\\vec q$", " 綁在一起：模態不再獨立", wait=2.2)
    s.play(FadeIn(real), run_time=1.0)
    s.say("實空間：格距 ", "$a$", " 的正方晶格", wait=1.6)
    s.play(FadeIn(recip), run_time=1.0)
    s.say("倒晶格：格距 ", "$2\\pi/a$", "；每個點是一個 ", "$\\vec q$", wait=2.0)
    s.say("實空間越大，", "$\\vec q$", " 越密（", "$L\\propto N^{1/d}$", "）", wait=2.2)
    s.clear_stage()


def part_shell(s):
    nb = numbers()
    dq, lin = nb["dq"]
    disk = bz_disk(2.1, np.array([-3.7, 0.15, 0]), 0.45)
    cap = zh("Fig. 12.3–12.4（平面示意）", font_size=20, color=GREY_B).next_to(disk, DOWN, buff=0.15)
    lam = mixed("$\\ell=0.5$", "：", "$\\Lambda=%s$" % f4(nb["lam"]), font_size=26)
    dqt = tx(r"dq=e^{-\ell}-e^{-\ell-\delta\ell}=%s\approx e^{-\ell}\delta\ell=%s" % (f"{dq:.6f}", f"{lin:.6f}"), 26)
    dtag = mixed("(12.24)　", "$\\delta\\ell=0.01$", font_size=20, color=GREY_B)
    rows = [["d", r"\delta\ell", r"1-e^{-d\delta\ell}", r"d\,\delta\ell"]]
    rows += [[str(d), str(dl), f4(a), f4(b)] for d, dl, a, b in nb["shell"]]
    tab = D.tex_table(rows, [0.8, 1.0, 2.0, 1.4], row_h=0.5, font_size=26)
    res = tx(r"\sum_{\rm shell}1\approx N\,d\,\delta\ell", 32).set_color(YELLOW)
    rtag = zh("(12.20)–(12.23)", font_size=20, color=GREY_B)
    right = VGroup(lam, dtag, dqt, tab, VGroup(rtag, res).arrange(RIGHT, buff=0.3)).arrange(DOWN, buff=0.3)
    place_in(right, -1.0, 6.7, -2.85, 3.1)
    assert_in_frame(disk, cap, right)
    assert_apart(VGroup(disk, cap), right)
    assert_column(*right)
    assert_table(tab, 4)
    s.say("把布里淵區換成半徑 ", "$\\Lambda=e^{-\\ell}$", " 的球", wait=2.0)
    s.play(ShowCreation(disk[0]), FadeIn(disk[2]), FadeIn(disk[3]), FadeIn(cap), run_time=1.2)
    s.play(FadeIn(disk[1]), run_time=0.8)
    s.say("外圈的殼層：", "$e^{-\\ell-\\delta\\ell}<q<e^{-\\ell}$", "，就是要積掉的部分", wait=2.2)
    s.play(FadeIn(lam), FadeIn(dtag), Write(dqt), run_time=1.4)
    s.say("殼層的厚度 ", "$dq\\approx e^{-\\ell}\\delta\\ell$", "：", "$\\ell=0.5$", " 時差不到 1%", wait=2.2)
    s.play(FadeIn(tab), run_time=1.0)
    s.say("殼層佔的比例：", "$1-e^{-d\\delta\\ell}$", " 對比 ", "$d\\,\\delta\\ell$", wait=2.0)
    s.say("$\\delta\\ell$", " 越小，兩欄越接近", wait=1.8)
    s.play(FadeIn(right[-1]), run_time=0.8)
    s.say("所以殼層裡的模態數 ", "$\\approx N\\,d\\,\\delta\\ell$", wait=2.0)
    s.say("半徑 200 的球內數格點，結果一樣（差不到 0.0002）", wait=2.2)
    s.clear_stage()


def part_toy(s):
    nb = numbers()
    rp, up, c6 = nb["toy"]
    z1, _ = nb["Z"]
    boxes = VGroup(flow_box(["分出殼層 ", "$S'$"], "(12.5)"),
                   flow_box(["對 ", "$S'$", " 求和：", "$A(K')e^{-\\mathcal H'}$"], "(12.6)(12.7)"),
                   flow_box(["重新標度"], "(12.29)–(12.34)"),
                   flow_box(["$K'=\\mathcal R(K)$"], "(12.10)")).arrange(RIGHT, buff=0.5)
    arrows = VGroup(*[Arrow(a.get_right(), b.get_left(), buff=0.05, thickness=3).set_color(GREY_B)
                      for a, b in zip(boxes, boxes[1:])])
    chart = VGroup(boxes, arrows)
    place_in(chart, -6.7, 6.7, 1.0, 2.9)
    assert_in_frame(chart)
    s.say("怎樣「正確地」積掉自由度：四個步驟", wait=1.8)
    s.play(LaggedStart(*[FadeIn(b) for b in boxes], lag_ratio=0.3), FadeIn(arrows), run_time=2.0)
    s.say("對殼層的 ", "$S'$", " 求和，剩下內區的 ", "$\\mathcal H'$", "，前面多一個常數 ", "$A$", wait=2.2)
    s.say("再重新標度，就得到新的參數 ", "$K'=\\mathcal R(K)$", wait=2.0)
    ham = tx(r"\mathcal H=\frac r2s^2+us^4+\frac{r+1}2S'^2+6u\,s^2S'^2", 30)
    ttag = mixed("兩變數玩具（書外）：", "$(r,u)=(0.2,0.01)$", font_size=22, color=GREY_A)
    top = VGroup(ttag, ham).arrange(DOWN, buff=0.15).move_to(np.array([0, 0.1, 0]))
    assert_in_frame(top)
    assert_apart(chart, top)
    s.play(FadeIn(ttag), Write(ham), run_time=1.4)
    s.say("玩具：", "$s$", " 代表內區、", "$S'$", " 代表殼層；只留 ", "$6u$", " 頂點", wait=2.2)
    s.play(FadeOut(chart), top.animate.move_to(np.array([0, 2.6, 0])), run_time=1.0)
    ax = toy_axes()
    ax.move_to(np.array([-3.6, -1.15, 0]))
    c1 = toy_curve(ax, nb["heff"], YELLOW)
    c2 = toy_curve(ax, nb["quart"], ORANGE, dashed=True)
    l1 = tx(r"\mathcal H'(s)", 24, YELLOW)
    l2 = tx(r"\frac{r'}2s^2+u's^4", 24, ORANGE)
    VGroup(l1, l2).arrange(RIGHT, buff=0.8).next_to(ax, UP, buff=0.15)
    assert_apart(l1, l2)
    heff = tx(r"\mathcal H'(s)=\frac r2s^2+us^4+\frac12\ln\Big(1+\frac{12us^2}{r+1}\Big)", 28)
    r1 = tx(r"r'=r+\frac{12u}{1+r}=%s" % f4(rp), 28).set_color(YELLOW)
    u1 = tx(r"u'=u-\frac{36u^2}{(1+r)^2}=%s" % f4(up), 28).set_color(YELLOW)
    s6 = mixed("新的 ", "$s^6$", " 項：", "$%s$" % f"{c6:.6f}", font_size=24, color=GREY_A)
    zz = tx(r"Z=%s=A\cdot Z'=%s\times%s" % (f4(z1), f4(nb["A"]), f4(z1 / nb["A"])), 28)
    right = VGroup(heff, r1, u1, s6, zz).arrange(DOWN, buff=0.28)
    place_in(right, 0.0, 6.7, -2.85, 2.1)
    assert_in_frame(ax, c1, c2, l1, l2, right)
    assert_apart(VGroup(ax, l1, l2), right)
    assert_apart(top, right)
    assert_apart(top, VGroup(ax, l1, l2))
    assert_column(*right)
    s.play(Write(heff), run_time=1.4)
    s.say("$S'$", " 的 Gaussian 積分做得出來：多一個對數", wait=2.0)
    s.say("把對數展開：", "$\\frac12\\ln(1+x)=\\frac x2-\\frac{x^2}4+\\frac{x^3}6-\\cdots$", wait=2.4)
    s.play(FadeIn(r1), run_time=0.8)
    s.say("$s^2$", " 的係數：", "$r$", " 變成 ", "$%s$" % f4(rp), wait=1.8)
    s.play(FadeIn(u1), run_time=0.8)
    s.say("$s^4$", " 的係數：", "$u$", " 變成 ", "$%s$" % f4(up), wait=1.8)
    s.play(FadeIn(s6), run_time=0.6)
    s.say("還多了 ", "$s^6$", " 項——書上把這種高階項丟掉", wait=2.0)
    s.play(ShowCreation(ax), run_time=0.8)
    s.play(ShowCreation(c1), FadeIn(l1), run_time=1.2)
    s.play(ShowCreation(c2), FadeIn(l2), run_time=1.2)
    s.say("只留 ", "$r'$", "、", "$u'$", " 的曲線和真正的 ", "$\\mathcal H'$", " 幾乎重合", wait=2.2)
    s.play(FadeIn(zz), run_time=0.8)
    s.say("配分函數不變：", "$Z=A\\cdot Z'$", "（(12.8)(12.9)）", wait=2.0)
    s.say("12 和 36 就是下一幕要在殼層裡數出來的數字", wait=2.2)
    s.say("書上 (12.28) 還留著 ", "$d$", "；本系列把 ", "$d$", " 併進 ", "$u$", wait=2.2)
    s.clear_stage()


def part_rescale(s):
    nb = numbers()
    icons = VGroup()
    for k in range(5):
        ic = vertex_icon(k, ORIGIN)
        c = tx(str(math.comb(4, k)), 30, YELLOW).next_to(ic, DOWN, buff=0.2)
        n = mixed("$%d$" % k, " 隻在殼層", font_size=20, color=GREY_B).next_to(c, DOWN, buff=0.12)
        icons.add(VGroup(ic, c, n))
    icons.arrange(RIGHT, buff=0.8)
    tag = mixed("(12.11)　四次項依殼層腳數分五項：", "$\\binom4k$", font_size=24, color=GREY_A)
    top = VGroup(tag, icons).arrange(DOWN, buff=0.25)
    place_in(top, -6.7, 6.7, 0.7, 3.1)
    cross = VGroup(*[Cross(icons[k][0], stroke_width=4).set_color(GREY_B) for k in (3, 4)])
    steps = VGroup(tx(r"N\to N'e^{d\delta\ell}", 28), tx(r"\vec q\to\vec q\,'e^{-\delta\ell}", 28),
                   tx(r"S\to S''e^{\delta\ell}", 28)).arrange(DOWN, buff=0.2)
    rtag = mixed("(12.29)–(12.34)　", "$d=3$", "、", "$\\delta\\ell=0.1$", " 的乘數", font_size=22, color=GREY_B)
    a, b, c = nb["resc"]
    rows = [[zh("項", font_size=24, color=YELLOW), zh("乘數", font_size=24, color=YELLOW)],
            ["q^2", f4(a)], ["r", r"e^{2\delta\ell}=%s" % f4(b)], ["u", r"e^{(4-d)\delta\ell}=%s" % f4(c)]]
    tab = D.tex_table(rows, [1.2, 3.4], row_h=0.48, font_size=26)
    bottom = VGroup(steps, VGroup(rtag, tab).arrange(DOWN, buff=0.2)).arrange(RIGHT, buff=1.0)
    place_in(bottom, -6.7, 6.7, -2.85, 0.4)
    assert_in_frame(top, bottom)
    assert_apart(top, bottom)
    assert_table(tab, 2)
    s.say("回到真正的模型：四次項有四隻腳", wait=1.6)
    s.play(FadeIn(tag), LaggedStart(*[FadeIn(i) for i in icons], lag_ratio=0.25), run_time=2.0)
    s.say("依幾隻腳在殼層分成五項，係數 ", "$1,4,6,4,1$", wait=2.0)
    s.play(ShowCreation(cross), run_time=0.8)
    s.say("三、四隻腳在殼層的兩項：", "$\\delta\\ell\\to0$", " 時可略", wait=2.2)
    s.play(FadeIn(steps), run_time=1.0)
    s.say("積完之後重新標度：模態數、波數、自旋各換一次", wait=2.0)
    s.play(FadeIn(rtag), FadeIn(tab), run_time=1.0)
    s.say("$q^2$", " 項乘 1，模型保持同形；", "$r$", " 乘 ", "$%s$" % f4(b), wait=2.2)
    s.say("$u$", " 乘 ", "$e^{(4-d)\\delta\\ell}=%s$" % f4(c), "：四維以下 ", "$u$", " 會長大", wait=2.2)
    s.play(FadeOut(top), FadeOut(cross), FadeOut(bottom), run_time=0.8)
    box1 = D.erratum_box("p. 295，(12.32)", r"N\to N'e^{d\delta\ell}\approx N'(1+\delta\ell)",
                         r"N\to N'e^{d\delta\ell}\approx N'(1+d\,\delta\ell)", width=11.0, font_size=30)
    box2 = D.erratum_box("p. 295，(12.33)", r"\vec q\to\vec q\,'e^{-d\delta\ell}",
                         r"\vec q\to\vec q\,'e^{-\delta\ell}", width=11.0, font_size=30)
    col = VGroup(box1, box2).arrange(DOWN, buff=0.6)
    place_in(col, -6.7, 6.7, -2.6, 2.9)
    assert_in_frame(col)
    assert_apart(box1, box2)
    s.play(FadeIn(box1), run_time=0.8)
    s.say("$d=3$", "、", "$\\delta\\ell=0.01$", "：", "$e^{0.03}=%s$" % f4(math.exp(0.03)), "，不是 ",
          "$%s$" % f4(1.01), wait=2.4)
    s.play(FadeIn(box2), run_time=0.8)
    s.say("照書上 ", "$q^2$", " 項要乘 ", "$%s$" % f4(nb["book"][0]), "（", "$d=3$", "、", "$\\delta\\ell=0.1$",
          "），模型就不再同形", wait=2.6)
    s.say("積掉殼層、重新標度之後，只剩 ", "$r$", " 與 ", "$u$", " 兩個數在變", wait=2.2)
    s.say("——下一幕數殼層裡的配對", wait=1.8)


class W1201QuarticShell(ChuScene):
    def construct(self):
        self.show_header("W1201　加上四次項：怎樣正確地積掉自由度")
        check_math()
        part_question(self)
        part_shell(self)
        part_toy(self)
        part_rescale(self)
        self.wait(1.0)
