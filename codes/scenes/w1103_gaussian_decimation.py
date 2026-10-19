"""
W1103 實空間重整化：一維 Gaussian 鏈（§11.2 前半，pp. 264–269，(11.26)–(11.49)，Fig. 11.6–11.7，Ex. 11.5–11.6）

問：把一半的自旋積掉，剩下的系統還長得一樣嗎？
  1. 9 個連續自旋（長短不一的箭頭），奇數格點淡出、偶數格點之間換上新彈簧；(11.27) ℋ = Σ[r/2 S_i² + a/2 (S_i − S_{i+1})²]。
  2. (11.28)–(11.30) 完成平方：數值積分 1.224819 ＝ 閉式 1.224819；(11.32)(11.33) r'、a'；
     筆誤框 1（p. 265，(11.31)：(2π/(r+2a))^{N/2} → ^{N/4}）：8 個自旋的環，照書上 ln Z 多 1.8432。
  3. (11.34)–(11.38) 重新標度：r'' = r(r+4a)/a、a'' = a、h''；Fig. 11.6–11.7 蛛網 0.05 → 0.2025 → 0.8510 → 4.1282（→ 33.5553）；
     筆誤框 2（p. 267，Fig. 11.6–11.7：座標 r/√a → r/a）；正文的 t'' 也是 r''（字幕）。
  4. 小 r 表（r = h = 1e−4，k = 1..5）⇒ r^(n) = 4^n r、h^(n) = 2^{3n/2} h、b^(n) = 2^{n/2}；h''/h → 2.8284。
  5. (11.46)–(11.49) m = 2^{n/2} f(…) ⇒ γ = 1；筆誤框 3（p. 269，2^{−n/2} → 2^{+n/2}）：照書上 0.0500、更正 0.1000 = h/r。
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

SPINS = (0.6, 0.3, -0.5, -0.7, 0.4, 0.8, -0.3, 0.5, -0.6)      # S_1..S_9（S_2、S_4 和 square_check 相同）
SPRING_COLOR = GREY_B
NEW_SPRING = YELLOW
ODD_COLOR = BLUE_B
EVEN_COLOR = ORANGE
ITER = ("0.0500", "0.2025", "0.8510", "4.1282", "33.5553")
SMALL_R = (("1.0000", "1.0001", "1.0005", "1.0021", "1.0086"),
           ("1.0000", "1.0000", "1.0000", "1.0000", "1.0000"),
           ("1.0000", "1.0001", "1.0005", "1.0021", "1.0085"))


def f4(x):
    return f"{float(x):.4f}"


def check_math():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        num, closed = RG.square_check()
        assert f"{num:.6f}" == f"{closed:.6f}" == "1.224819"
        assert (SPINS[1], SPINS[3]) == (0.3, -0.7)
        r, a = 0.5, 1.0
        rp, ap, _ = RG.decimate_gauss(r, a)
        assert abs(rp - r * (r + 4 * a) / (r + 2 * a)) < 1e-12 and abs(ap - a * a / (r + 2 * a)) < 1e-12
        assert f4(RG.ring_decimation_check(8, 0.5, 1.0, 0.0, book=True)) == "1.8432"
        assert abs(RG.ring_decimation_check(8, 0.5, 1.0, 0.0)) < 1e-9
        rr, aa, hh = RG.gauss_rescaled(0.5, 1.0, 0.1)
        assert abs(rr - RG.rr_gauss(0.5)) < 1e-12 and abs(aa - 1.0) < 1e-12 and abs(hh - RG.hh_gauss(0.5, 0.1)) < 1e-12
        assert abs(RG.spin_factor(0.0) - np.sqrt(2)) < 1e-12
        it = RG.gauss_iterate(RG.COBWEB_R0, 4)
        assert [f4(x) for x in it] == list(ITER) and RG.COBWEB_R0 == 0.05
        cw = RG.cobweb_points()
        assert len(cw) == 6 and abs(cw[-1][1] - it[3]) < 1e-12 and all(0 <= x <= 10 and 0 <= y <= 10 for x, y in cw)
        assert it[4] > 10
        rc = RG.rr_curve()
        assert len(rc) == 40 and f4(rc[-1][0]) == "1.7417" and abs(rc[-1][1] - 10) < 1e-9 and rc[0] == (0.0, 0.0)
        assert RG.gauss_fixed_points() == (0.0, -3.0)
        eps = 1e-6
        assert abs((RG.rr_gauss(eps) - RG.rr_gauss(0.0)) / eps - 4) < 1e-4            # 斜率 4
        tab = RG.small_r_table()
        assert [k for k, *_ in tab] == [1, 2, 3, 4, 5]
        for col in range(3):
            assert tuple(f4(row[col + 1]) for row in tab) == SMALL_R[col]
        assert f4(RG.hh_gauss(0.0, 1.0)) == "2.8284" == f4(2 ** 1.5)
        lr, lh, _ = RG.chain_lambdas()
        assert f4(lr) == "2.0000" and f4(lh) == "1.5000" and abs(np.log(2) / np.log(4) - 0.5) < 1e-12
        book, right, exact = RG.book1147_m(0.01, 0.001, 1)
        assert (f4(book), f4(right), f4(exact)) == ("0.0500", "0.1000", "0.1000") and abs(exact - 0.001 / 0.01) < 1e-12


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


def spring(p, q, color=SPRING_COLOR, coils=6, amp=0.09, width=2.5):
    p, q = np.array(p, dtype=float), np.array(q, dtype=float)
    d = q - p
    n = np.array([-d[1], d[0], 0]) / np.linalg.norm(d)
    pts = [p, p + 0.15 * d]
    for k in range(2 * coils):
        t = 0.15 + 0.7 * (k + 0.5) / (2 * coils)
        pts.append(p + t * d + (amp if k % 2 == 0 else -amp) * n)
    pts += [p + 0.85 * d, q]
    return VMobject().set_points_as_corners(pts).set_stroke(color, width)


def chain(y0=0.0, unit=1.5):
    xs = [-4.0 + i for i in range(9)]
    dots = VGroup(*[Dot([x, y0, 0], radius=0.07) for x in xs])
    arrows = VGroup()
    for i, (x, sv) in enumerate(zip(xs, SPINS)):
        col = ODD_COLOR if i % 2 == 0 else EVEN_COLOR           # i = 0 是 S_1（奇數）
        arrows.add(Arrow([x, y0, 0], [x, y0 + unit * sv, 0], buff=0, thickness=4).set_color(col))
    springs = VGroup(*[spring([xs[i] + 0.1, y0, 0], [xs[i + 1] - 0.1, y0, 0]) for i in range(8)])
    labels = VGroup(*[tx(f"S_{i + 1}", 26, GREY_B).move_to([x, y0 - 1.45, 0]) for i, x in enumerate(xs)])
    new = VGroup(*[spring([xs[i] + 0.1, y0, 0], [xs[i + 2] - 0.1, y0, 0], NEW_SPRING, coils=9, amp=0.12, width=3)
                   for i in (1, 3, 5)])
    return dots, arrows, springs, labels, new


# ---------------------------------------------------------------- 1. 積掉一半
def part_chain(s):
    y0 = 0.3
    dots, arrows, springs, labels, new = chain(y0)
    h0 = tx(r"\mathcal H=\sum_i\Big[\frac r2S_i^2+\frac a2\big(S_i-S_{i+1}\big)^2\Big]", 38)
    h0.move_to([0, 2.55, 0])
    tag = tx("(11.27)", 24, GREY_B).next_to(h0, LEFT, buff=0.4)
    odd = [0, 2, 4, 6, 8]
    note = mixed("偶數格點之間：新的耦合 ", "$a'$", "、新的 ", "$r'$", font_size=28).move_to([0, -2.3, 0])
    assert_in_frame(h0, tag, dots, arrows, springs, labels, new, note)
    for a, b in ((h0, arrows), (arrows, labels), (labels, note), (h0, tag)):
        assert_apart(a, b)

    s.say("把一半的自旋積掉，剩下的系統還長得一樣嗎？", wait=1.5)
    s.play(FadeIn(dots), ShowCreation(springs), FadeIn(labels), run_time=1.0)
    s.play(*[GrowArrow(a) for a in arrows], run_time=1.0)
    s.say("一維 Gaussian 鏈：自旋是連續的實數，畫成長短不一的箭頭", wait=2.0)
    s.play(FadeIn(h0), FadeIn(tag), run_time=0.8)
    s.say("每個自旋有 ", r"$\frac r2S_i^2$", "，相鄰的用彈簧 ", r"$\frac a2(S_i-S_{i+1})^2$", " 連起來", wait=2.5)
    s.say("奇數格點的自旋：一個一個積分掉", wait=1.5)
    s.play(*[FadeOut(arrows[i]) for i in odd], *[dots[i].animate.set_opacity(0.2) for i in odd],
           *[labels[i].animate.set_opacity(0.3) for i in odd], run_time=1.2)
    s.play(FadeOut(springs), *[ShowCreation(m) for m in new], run_time=1.2)
    s.play(FadeIn(note), run_time=0.6)
    s.say("剩下的偶數自旋之間，換上了新的彈簧", wait=2.0)
    s.say("模型的樣子不變，只是參數變了：", r"$(r,a)\to(r',a')$", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 2. 完成平方
def part_square(s):
    e1 = tx(r"\int dS_3\,\exp\Big[-\Big(\frac r2+a\Big)S_3^2+aS_3(S_2+S_4)-\frac a2\big(S_2^2+S_4^2\big)\Big]", 30)
    e2 = tx(r"=\sqrt{\frac{2\pi}{r+2a}}\,\exp\Big[-\frac{a^2(S_2-S_4)^2}{2(r+2a)}-\frac{ar\big(S_2^2+S_4^2\big)}{2(r+2a)}\Big]", 30)
    place_left(e1, -6.3, 2.55)
    place_left(e2, -5.6, 1.45)
    num, closed = RG.square_check()
    e3 = place_left(mixed(r"$r=0.5,\ a=1,\ S_2=0.3,\ S_4=-0.7$", "：數值積分 ", f"${num:.6f}$", "＝閉式 ",
                          f"${closed:.6f}$", font_size=26), -6.3, 0.4)
    e4 = tx(r"r'=\frac{r(r+4a)}{r+2a},\qquad a'=\frac{a^2}{r+2a}", 36, YELLOW)
    e4.move_to([-0.6, -0.65, 0])
    tag = tx("(11.32)(11.33)", 24, GREY_B).next_to(e4, RIGHT, buff=0.5)
    err = D.erratum_box("p. 265，(11.31)", r"\Big(\frac{2\pi}{r+2a}\Big)^{N/2}", r"\Big(\frac{2\pi}{r+2a}\Big)^{N/4}", width=6.4)
    err.move_to([0, -2.1, 0])
    assert_in_frame(e1, e2, e3, e4, tag, err)
    for a, b in ((e1, e2), (e2, e3), (e3, e4), (e4, err), (tag, err), (e3, tag)):
        assert_apart(a, b)

    s.say("積掉 ", "$S_3$", "：它只和鄰居 ", "$S_2$", "、", "$S_4$", " 有關", wait=1.5)
    s.play(FadeIn(e1), run_time=1.0)
    s.say("指數是 ", "$S_3$", " 的二次式：完成平方", wait=2.0)
    s.play(FadeIn(e2), run_time=1.0)
    s.say("剩下 ", "$S_2$", " 和 ", "$S_4$", " 之間的新彈簧，和各自的新 ", "$r$", wait=2.5)
    s.play(FadeIn(e3), run_time=0.7)
    s.say("數值核對：積分和閉式到小數第六位都相同", wait=2.0)
    s.play(FadeIn(e4), FadeIn(tag), run_time=0.8)
    s.say("讀出新的參數 ", "$r'$", " 和 ", "$a'$", wait=2.0)
    s.play(FadeIn(err), run_time=0.8)
    s.say("$N/2$", " 個奇數格點，每個只給一個平方根：次方是 ", "$N/4$", wait=2.5)
    s.say("8 個自旋的環用行列式算：照書上 ", r"$\ln Z$", " 多 ",
          f"${f4(RG.ring_decimation_check(8, 0.5, 1.0, 0.0, book=True))}$", wait=2.5)
    s.clear_stage()


# ---------------------------------------------------------------- 3. 重新標度與蛛網
def part_cobweb(s):
    x0 = -6.3
    t1 = place_left(tx(r"S_{2i}\to S_{2i}\sqrt{(r+2a)/a}", 32), x0, 2.6, 0.6)
    t2 = place_left(tx(r"r''=\frac{r(r+4a)}a,\qquad a''=a", 34, YELLOW), x0, 1.65, 0.6)
    t3 = place_left(tx(r"h''=h\,\frac{r+4a}{\sqrt{a(r+2a)}}", 32), x0, 0.6, 0.6)
    it = RG.gauss_iterate(RG.COBWEB_R0, 4)
    vals = place_left(tx(r"r_0=" + r"\to".join(f4(v) for v in it[:4]), 30), x0, -0.45, 0.6)
    nxt = place_left(mixed("下一步 ", f"${f4(it[4])}$", "：跑出圖外", font_size=26), x0, -1.15, 0.6)
    fixed = place_left(mixed("不動點 ", "$r=0$", "（斜率 4，不穩定）、", "$r=-3a$", "（不物理）", font_size=26),
                       x0, -1.9, 0.6)
    ax = D.plot_axes((0, 10, 2), (0, 10, 2), width=4.6, height=4.6, x_ticks=(2, 4, 6, 8, 10), y_ticks=(2, 4, 6, 8, 10),
                     x_label="r/a", y_label="r''/a")
    ax.move_to([3.9, 0.3, 0])
    curve = VMobject().set_points_as_corners([ax.c2p(x, y) for x, y in RG.rr_curve()]).set_stroke(TEAL_B, 3)
    diag = DashedLine(ax.c2p(0, 0), ax.c2p(10, 10)).set_stroke(GREY_B, 2)
    cw = RG.cobweb_points() + [(it[3], it[3])]
    segs = VGroup(*[Line(ax.c2p(*cw[k]), ax.c2p(*cw[k + 1])).set_stroke(YELLOW, 2.5) for k in range(len(cw) - 1)])
    out = Arrow(ax.c2p(it[3], it[3]), ax.c2p(it[3], 10.4), buff=0, thickness=3).set_color(RED)
    err = D.erratum_box("p. 267，Fig. 11.6–11.7", r"r/\sqrt{a},\ r''/\sqrt{a}", r"r/a,\ r''/a", width=6.4)
    err.move_to([-3.2, -1.1, 0])
    assert_in_frame(t1, t2, t3, vals, nxt, fixed, ax, segs, out, err)
    for a, b in ((t1, t2), (t2, t3), (t3, vals), (vals, nxt), (nxt, fixed), (t3, err), (ax, t1), (ax, t2), (ax, vals),
                 (ax, fixed), (ax, err)):
        assert_apart(a, b)

    s.say("再把剩下的自旋放大：讓 ", "$a$", " 回到原來的值", wait=1.5)
    s.play(FadeIn(t1), run_time=0.7)
    s.say("自旋乘上 ", r"$\sqrt{(r+2a)/a}$", "：小 ", "$r$", " 時就是 ", r"$\sqrt2$", wait=2.0)
    s.play(FadeIn(t2), run_time=0.7)
    s.say("模型完全回到原樣，只有 ", "$r$", " 和 ", "$h$", " 變了", wait=2.0)
    s.play(FadeIn(t3), run_time=0.7)
    s.play(FadeIn(ax), ShowCreation(diag), run_time=1.0)
    s.play(ShowCreation(curve), run_time=1.0)
    s.say("Fig. 11.6–11.7：", "$r''$", " 對 ", "$r$", "，虛線是 ", "$r''=r$", wait=2.0)
    for seg in segs:
        s.play(ShowCreation(seg), run_time=0.5)
    s.play(FadeIn(vals), run_time=0.6)
    s.say("從 ", "$r_0=0.05$", " 出發：越跑越大", wait=2.0)
    s.play(GrowArrow(out), FadeIn(nxt), run_time=0.8)
    s.say("下一步 ", f"${f4(it[4])}$", "：跑出圖外", wait=1.8)
    s.play(FadeIn(fixed), run_time=0.6)
    s.say("$r=0$", " 是不動點，斜率 4：一碰就跑開，不穩定", wait=2.5)
    s.say("另一個不動點 ", "$r=-3a$", "：", "$r<0$", " 不是物理的", wait=2.0)
    s.play(FadeOut(vals), FadeOut(nxt), FadeOut(fixed), run_time=0.5)
    s.play(FadeIn(err), run_time=0.8)
    s.say("$r$", " 和 ", "$a$", " 同單位：座標是 ", "$r/a$", "，", r"$r''/a=(r/a)(r/a+4)$", wait=2.5)
    s.say("正文的「", "$t''$", "」也就是 ", "$r''$", wait=2.0)
    s.clear_stage()


# ---------------------------------------------------------------- 4. 小 r
def part_small(s):
    tab_rows = RG.small_r_table()
    rows = [["k"] + [str(k) for k, *_ in tab_rows],
            [r"r^{(k)}/(4^kr)"] + [f4(row[1]) for row in tab_rows],
            [r"h^{(k)}/(2^{3k/2}h)"] + [f4(row[2]) for row in tab_rows],
            [mixed("自旋總因子", r"$/2^{k/2}$", font_size=22)] + [f4(row[3]) for row in tab_rows]]
    tab = D.tex_table(rows, [3.4] + [1.3] * 5, row_h=0.62, font_size=28).shift(np.array([-5.6, 2.85, 0]))
    cond = tx(r"r=h=10^{-4}", 28, GREY_B).move_to([5.6, 2.55, 0])
    c1 = tx(r"r^{(n)}=4^nr,\qquad h^{(n)}=2^{3n/2}h,\qquad b^{(n)}=2^{n/2}", 36, YELLOW).move_to([0, -0.4, 0])
    c2 = mixed(r"$r\to0$", "：", rf"$h''/h\to2^{{3/2}}={f4(RG.hh_gauss(0.0, 1.0))}$", font_size=28).move_to([0, -1.5, 0])
    assert_in_frame(tab, cond, c1, c2)
    assert_table(tab, 6)
    for a, b in ((tab, cond), (tab, c1), (c1, c2)):
        assert_apart(a, b)

    s.say("$r$", " 很小時，每做一次會怎樣？", wait=1.5)
    s.play(FadeIn(tab), FadeIn(cond), run_time=1.0)
    s.say("做 ", "$k$", " 次，除以猜的答案：比值都接近 1", wait=2.5)
    s.say("$r$", " 每次乘 4，", "$h$", " 每次乘 ", r"$2^{3/2}$", "，自旋每次乘 ", r"$\sqrt2$", wait=2.5)
    s.play(FadeIn(c1), run_time=0.8)
    s.say("做 ", "$n$", " 次：長度單位變成 ", "$2^n$", " 倍", wait=2.0)
    s.play(FadeIn(c2), run_time=0.6)
    s.say("$h$", " 的放大倍數趨向 ", f"${f4(RG.hh_gauss(0.0, 1.0))}$", wait=2.0)
    s.clear_stage()


# ---------------------------------------------------------------- 5. 標度形式
def part_scaling(s):
    e1 = place_left(mixed("(11.46)  ", r"$m=2^{n/2}f\big(r2^{2n},\ h2^{3n/2}\big)$", "  與 ", "$n$", " 無關", font_size=30),
                    -6.0, 2.55)
    e2 = place_left(tx(r"\Rightarrow\ m=r^{-1/4}F\big(h/r^{3/4}\big)", 36), -6.0, 1.6)
    e3 = place_left(mixed(r"$\frac{\partial m}{\partial h}=r^{-1}F'$", r"$\Rightarrow$", r"$\gamma=1$",
                          "：正是 (11.24) 取 ", "$d=1$", font_size=30), -6.0, 0.5)
    err = D.erratum_box("p. 269，(11.46)(11.47)", r"m=2^{-n/2}f\big(r2^{2n},h2^{3n/2}\big)",
                        r"m=2^{+n/2}f\big(r2^{2n},h2^{3n/2}\big)", width=8.0)
    err.move_to([0, -1.25, 0])
    assert_in_frame(e1, e2, e3, err)
    for a, b in ((e1, e2), (e2, e3), (e3, err)):
        assert_apart(a, b)
    book, right, exact = RG.book1147_m(0.01, 0.001, 1)

    s.say("自旋每次乘 ", r"$\sqrt2$", "：磁化 ", "$m$", " 也跟著標度", wait=1.5)
    s.play(FadeIn(e1), run_time=0.8)
    s.say("同一個系統，做幾次重整化都一樣", wait=2.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("挑 ", "$n$", " 讓 ", r"$r2^{2n}=1$", "：一個變數的函數", wait=2.0)
    s.play(FadeIn(e3), run_time=0.8)
    s.say("又是 ", r"$\gamma=1$", "：和 W1102 的 Gaussian 模型一致", wait=2.0)
    s.play(FadeIn(err), run_time=0.8)
    s.say("精確解 ", "$m=h/r$", "（", "$r=0.01$", "、", "$h=0.001$", "、", "$n=1$", "）", wait=2.0)
    s.say("照書上 ", f"${f4(book)}$", "，更正後 ", f"${f4(right)}$", "＝", "$h/r$", wait=2.5)
    s.say("積掉一半、重新標度：只有 ", "$r$", " 變成 ", "$4r$", wait=2.0)
    s.say("由「每次乘 4」讀出 ", r"$\nu=\frac12$", "、", r"$\gamma=1$", wait=2.5)


class W1103GaussianDecimation(ChuScene):
    def construct(self):
        self.show_header("W1103　實空間重整化：一維 Gaussian 鏈")
        check_math()
        part_chain(self)
        part_square(self)
        part_cobweb(self)
        part_small(self)
        part_scaling(self)
        self.wait(1.0)
