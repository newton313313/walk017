"""
W212 補充 3：從奇點讀出係數的成長＋第 2 章回顧（§2.4，pp. 47–50，(2:S3-1)–(2:S3-13)）

問：知道生成函數長什麼樣，不展開，能直接說出第 N 項大概多大嗎？
  1. 單極點 (2:S3-2)：[z^N] 1/(zc − z) = zc^{−(N+1)}；zc = 1/2、N = 10 → 2048。兩個極點 (2:S3-3)–(2:S3-5)：
     1/(1/2 − z) + 1/(1 − z) → c_N = 2^{N+1} + 1：3, 5, 9, 17, 33, 65；第二項的比重 e^{−δ(N+1)}（δ = ln 2）→ 最近的奇點說了算。
  2. 分支點 (2:S3-6)–(2:S3-10)：(zc − z)^{−α} 的係數 ≈ N^{α−1} zc^{−N−α}/Γ(α)；(1 − 4z)^{−1/2}：184756 vs 187078.97（0.9876）；
     一維回到原點 (1 − 4z²)^{−1/2}：分支點 ±1/2 一樣遠 → 偶數加倍、奇數抵消：252 vs 258.4（0.9754）。
  3. 對數奇點 (2:S3-11)：筆誤框 +ln zc → −ln zc；α = 1、zc = 1/2、N = 10：7418.1、7317.4、4478.3、4715.7。
     主例：正方格子 G(z; 0) ≈ −(1/π) ln(1 − 16z²) → 2·4^N/(πN)：66754 vs 63504（0.9513）…… N = 100：0.995。
  4. 第 2 章回顧：7 格地圖依序亮起。
畫面上的數字都由 utils.walk_gf 算出，check_math() 核對。
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import boxed, mixed, zh
from utils import walk_gf as GF
from utils import lattice_draw as L
from utils import walk_draw as D

RN = (10, 20, 50, 100)


def numbers():
    H10 = sum(1 / k for k in range(1, 11))
    return {"single": f"{GF.coeff_two_poles(10, 1, 0, 0.5, 1):.0f}",
            "two": [f"{GF.coeff_two_poles(N, 1, 1, 0.5, 1):.0f}" for N in range(6)],
            "w2": f"{GF.coeff_two_poles(10, 0, 1, 0.5, 1) / GF.coeff_two_poles(10, 1, 0, 0.5, 1):.5f}",
            "binom": f"{0.5 * GF.coeff_power_exact(10, 0.5, 0.25):.0f}", "basym": f"{0.5 * GF.coeff_power_asym(10, 0.5, 0.25):.2f}",
            "bratio": f"{GF.coeff_power_exact(10, 0.5, 0.25) / GF.coeff_power_asym(10, 0.5, 0.25):.4f}",
            "r1d": str(math.comb(10, 5)), "a1d": f"{GF.asym_returns_1d(10):.1f}",
            "q1d": f"{math.comb(10, 5) / GF.asym_returns_1d(10):.4f}",
            "log": [f"{v:.1f}" for v in (GF.coeff_log_exact(10, 1, 0.5), GF.coeff_log_asym(10, 1, 0.5),
                                         GF.coeff_log_asym(10, 1, 0.5, book=True), GF.coeff_log_asym(10, 1, 0.5, terms=1))],
            "logx": f"{2 ** 11 * (math.log(2) + H10):.1f}",
            "sq": str(GF.returns_sq(10)), "sqa": f"{GF.asym_returns_sq(10):.0f}",
            "sqr": [f"{GF.returns_sq(N) / GF.asym_returns_sq(N):.4f}".rstrip("0") for N in RN]}


def check_math():
    nb = numbers()
    assert nb["single"] == "2048" and nb["two"] == ["3", "5", "9", "17", "33", "65"] and nb["w2"] == "0.00049"
    assert all(GF.coeff_two_poles(N, 1, 1, 0.5, 1) == 2 ** (N + 1) + 1 for N in range(12))
    assert abs(GF.coeff_two_poles(10, 0, 1, 0.5, 1) / GF.coeff_two_poles(10, 1, 0, 0.5, 1) - math.exp(-math.log(2) * 11)) < 1e-15
    assert nb["binom"] == "184756" == str(math.comb(20, 10)) and nb["basym"] == "187078.97" and nb["bratio"] == "0.9876"
    assert abs(0.5 * GF.coeff_power_asym(10, 0.5, 0.25) - 4 ** 10 / math.sqrt(10 * math.pi)) < 1e-6
    assert abs(0.5 * GF.series_coeffs_power(0.5, 0.25, 5)[5] - 252) < 1e-9                 # (1 − 4w)^{−1/2}，w = z²
    assert (nb["r1d"], nb["a1d"], nb["q1d"]) == ("252", "258.4", "0.9754") and GF.asym_returns_1d(9) == 0
    assert nb["log"] == ["7418.1", "7317.4", "4478.3", "4715.7"] and nb["logx"] == nb["log"][0]
    assert (nb["sq"], nb["sqa"]) == ("63504", "66754") and [float(r) for r in nb["sqr"]] == GF.SQ_RATIOS
    assert abs(16 ** 5 / (5 * math.pi) - GF.asym_returns_sq(10)) < 1e-6                    # [z^10] −(1/π)ln(1 − 16z²)


# ---------------------------------------------------------------- 小工具
def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def tagged(tex, tag, font_size=30, color=WHITE):
    return VGroup(tx(tex, font_size, color), tx(r"\text{" + tag + "}", 22, GREY_B)).arrange(RIGHT, buff=0.3)


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


def zplane(center, w=4.4, h=3.5):
    ax = D.plot_axes((-1.5, 1.5, 0.5), (-1.2, 1.2, 0.5), width=w, height=h, font_size=20,
                     x_label=r"\mathrm{Re}\,z", y_label=r"\mathrm{Im}\,z")
    return ax.move_to(np.array(center))


def dashed_circle(ax, r, color=GREY_B):
    rad = np.linalg.norm(ax.c2p(r, 0) - ax.c2p(0, 0))
    return DashedVMobject(Circle(radius=rad).move_to(ax.c2p(0, 0)), num_dashes=30).set_stroke(color, 2)


def stage_card(title, tex, tag, w=3.15, h=1.75):
    t = zh(title, font_size=24, color=YELLOW)
    f = fit(Tex(tex, font_size=26), w - 0.3)
    g = zh(tag, font_size=17, color=GREY_B)
    body = VGroup(t, f, g).arrange(DOWN, buff=0.16)
    frame = RoundedRectangle(width=w, height=h, corner_radius=0.15).set_stroke(GREY_B, 2).set_fill(BLACK, 0.5)
    frame.move_to(body)
    return VGroup(frame, body)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：知道生成函數長什麼樣，不展開", wait=1.5)
    s.say("能不能直接說出第 ", "$N$", " 項大概多大？", wait=3.0)


def part_poles(s):
    nb = numbers()
    ax = zplane((-4.4, -0.2, 0))
    near = D.pole_marker(ax, 0.5 + 0j, r"\tfrac12", color=YELLOW)
    far = D.pole_marker(ax, 1 + 0j, "1", color=GREY_B)
    circ = dashed_circle(ax, 0.5)
    r1 = tagged(r"\frac{1}{z_c-z}=\frac{1}{z_c}\sum_n\Big(\frac{z}{z_c}\Big)^n\ \Rightarrow\ [z^N]=z_c^{-(N+1)}", "(2:S3-2)")
    r2 = mixed("$z_c=\\tfrac12,\\ N=10$", "：", "$2^{11}=" + nb["single"] + "$", font_size=28)
    r3 = tagged(r"\frac{1}{\frac12-z}+\frac{1}{1-z}\ \Rightarrow\ c_N=2^{N+1}+1", "(2:S3-4)")
    r4 = tx(r"c_N=" + r",\ ".join(nb["two"]) + r",\ \dots", 30)
    r5 = tagged(r"c_N=2^{N+1}\big[1+e^{-\delta(N+1)}\big],\quad\delta=\ln2", "(2:S3-5)")
    r6 = mixed("$N=10$", "：第二項只佔 ", "$2^{-11}\\approx" + nb["w2"] + "$", font_size=28)
    col = VGroup(r1, r2, r3, r4, r5, r6).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    for m in col:
        fit(m, 7.9)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.3), -1.25, 3.1)
    assert_in_frame(ax, near, far, col)
    assert_apart(ax, col)
    s.play(FadeIn(r1), run_time=1.2)
    s.say("最簡單的奇點：一個極點；係數是等比數列", wait=3.0)
    s.play(FadeIn(r2), run_time=0.8)
    s.say("$z_c=\\tfrac12$", "：第 10 項是 ", "$2^{11}=" + nb["single"] + "$", wait=2.5)
    s.play(FadeIn(ax), FadeIn(near), FadeIn(far), run_time=1.0)
    s.play(FadeIn(r3), FadeIn(r4), run_time=1.0)
    s.say("兩個極點：", "$z=\\tfrac12$", " 和 ", "$z=1$", wait=2.5)
    s.say("係數 ", "$2^{N+1}+1$", "：", "$+1$", " 來自遠的那個，很快就微不足道", wait=3.5)
    s.play(FadeIn(r5), FadeIn(r6), run_time=1.0)
    s.play(ShowCreation(circ), Indicate(near[0], color=YELLOW, scale_factor=1.8), run_time=1.2)
    s.say("最近的奇點說了算：它決定指數成長 ", "$z_c^{-N}$", wait=3.5)
    s.clear_stage()


def part_branch(s):
    nb = numbers()
    b1 = tagged(r"[z^N](z_c-z)^{-\alpha}=\frac{\Gamma(N+\alpha)}{\Gamma(\alpha)\Gamma(N+1)}z_c^{-N-\alpha}"
                r"\approx\frac{N^{\alpha-1}}{\Gamma(\alpha)}z_c^{-N-\alpha}", "(2:S3-10)")
    b2 = mixed("例：", "$(1-4z)^{-1/2}=\\sum_N\\binom{2N}{N}z^N$", "，", "$N=10$", "：", font_size=28)
    b3 = tx(r"\binom{20}{10}=" + nb["binom"] + r"\quad\text{vs}\quad\frac{4^{10}}{\sqrt{10\pi}}=" + nb["basym"]
            + r"\quad(" + nb["bratio"] + ")", 30)
    b23 = VGroup(b2, b3).arrange(RIGHT, buff=0.2)
    top = VGroup(b1, b23).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    for m in top:
        fit(m, 13.2)
    place_col(top.arrange(DOWN, aligned_edge=LEFT, buff=0.3), -6.6, 3.1)
    ax = zplane((-4.6, -1.0, 0), w=3.8, h=2.9)
    cuts = VGroup(L.branch_cut(ax, 0.5, 1.45), L.branch_cut(ax, -0.5, -1.45))
    bps = VGroup(*[Dot(ax.c2p(x, 0), radius=0.08).set_fill(YELLOW, 1) for x in (0.5, -0.5)])
    bl = VGroup(tx(r"\tfrac12", 24, YELLOW).next_to(ax.c2p(0.5, 0), UP, buff=0.15),
                tx(r"-\tfrac12", 24, YELLOW).next_to(ax.c2p(-0.5, 0), UP, buff=0.15))
    circ = dashed_circle(ax, 0.5)
    c1 = mixed("一維回到原點：", "$G(z;0)=(1-4z^2)^{-1/2}$", font_size=28)
    c2 = mixed("分支點 ", "$z=\\pm\\tfrac12$", "：一樣遠，兩個都要算", font_size=28)
    c3 = mixed("偶數 ", "$N$", "：加倍；奇數：抵消成 0", font_size=28, color=YELLOW)
    c4 = tx(r"\binom{10}{5}=" + nb["r1d"] + r"\quad\text{vs}\quad 2^{10}\sqrt{\frac{2}{10\pi}}=" + nb["a1d"]
            + r"\quad(" + nb["q1d"] + ")", 30)
    right = VGroup(c1, c2, c3, c4).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
    for m in right:
        fit(m, 7.6)
    place_col(right.arrange(DOWN, aligned_edge=LEFT, buff=0.3), -1.6, top.get_bottom()[1] - 0.35)
    assert_in_frame(top, ax, bl, right)
    assert_apart(top, ax)
    assert_apart(top, right)
    assert_apart(ax, right)
    s.say("分支點：", "$(z_c-z)^{-\\alpha}$", "，", "$\\alpha$", " 不是整數", wait=0.3)
    s.play(FadeIn(b1), run_time=1.2)
    s.say("$N$", " 很大時：指數成長 ", "$z_c^{-N}$", " 再乘一個冪次 ", "$N^{\\alpha-1}$", wait=3.5)
    s.play(FadeIn(b2), FadeIn(b3), run_time=1.0)
    s.say("$\\alpha=\\tfrac12$", "：中央二項式係數，差 1.2%", wait=3.0)
    s.play(FadeIn(ax), FadeIn(c1), run_time=1.0)
    s.play(FadeIn(bps), FadeIn(bl), ShowCreation(cuts), run_time=1.2)
    s.say("一維走回原點：分支點在 ", "$\\pm\\tfrac12$", "（鋸齒線是分支切割）", wait=3.0)
    s.play(ShowCreation(circ), FadeIn(c2), run_time=1.0)
    s.say("兩個奇點離原點一樣遠：兩個都要算", wait=3.0)
    s.play(FadeIn(c3), run_time=0.8)
    s.say("$-\\tfrac12$", " 那個多一個 ", "$(-1)^N$", "：偶數加倍、奇數抵消", wait=3.5)
    s.play(FadeIn(c4), run_time=0.8)
    s.say("又是奇偶性的因子 2（W206）", wait=3.0)
    s.clear_stage()


def part_log(s):
    nb = numbers()
    l1 = tagged(r"-(z_c-z)^{-\alpha}\ln(z_c-z)=\frac{d}{d\alpha}(z_c-z)^{-\alpha}", "(2:S3-11)")
    l2 = tx(r"[z^N]\approx\big[-\ln z_c-\psi(\alpha)+\ln N\big]\frac{N^{\alpha-1}}{\Gamma(\alpha)}z_c^{-N-\alpha}", 32)
    top = VGroup(l1, l2).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    fit(top, 13.0)
    place_col(top, -6.6, 3.1)
    rows = [[mixed("$\\alpha=1,\\ z_c=\\tfrac12,\\ N=10$", font_size=24, color=YELLOW), VectorizedPoint()],
            [mixed("精確 ", "$2^{11}(\\ln2+H_{10})$", font_size=24), nb["log"][0]],
            [zh("三項公式（正確）", font_size=24, color=D.RIGHT_COLOR), nb["log"][1]],
            [zh("照書上的號", font_size=24, color=D.WRONG_COLOR), nb["log"][2]],
            [mixed("只留 ", "$\\ln N$", font_size=24), nb["log"][3]]]
    tab = D.tex_table(rows, [4.4, 1.8], row_h=0.56, font_size=28)
    place_col(tab, -6.6, top.get_bottom()[1] - 0.4)
    D.table_cell(tab, 3, 1, 2).set_color(D.WRONG_COLOR)
    box = D.erratum_box("p. 50，(2:S3-11)", r"+\ln z_c", r"-\ln z_c", width=4.4)
    box.move_to(np.array([6.7 - box.get_width() / 2, tab.get_center()[1], 0]))
    assert_in_frame(top, tab, box)
    assert_apart(top, tab)
    assert_apart(tab, box)
    s.say("對數奇點：對 ", "$\\alpha$", " 微分就造得出來", wait=0.3)
    s.play(FadeIn(l1), run_time=1.2)
    s.play(FadeIn(l2), run_time=1.0)
    s.say("係數多了一個 ", "$\\ln N$", wait=2.5)
    s.play(FadeIn(box, shift=0.2 * LEFT), run_time=0.8)
    s.say("書上第一項的號印成 ", "$+\\ln z_c$", "：應是 ", "$-\\ln z_c$", wait=3.5)
    s.play(FadeIn(tab), run_time=1.2)
    s.say("數字核對：正確的三項公式差 1.4%，照書上的差四成", wait=3.5)
    s.say("只留 ", "$\\ln N$", " 也不夠：常數項在 ", "$N=10$", " 還很重要", wait=3.0)
    s.clear_stage()


def part_square(s):
    nb = numbers()
    ax = zplane((-4.4, -0.3, 0))
    cuts = VGroup(L.branch_cut(ax, 0.25, 1.45), L.branch_cut(ax, -0.25, -1.45))
    sps = VGroup(*[Dot(ax.c2p(x, 0), radius=0.08).set_fill(YELLOW, 1) for x in (0.25, -0.25)])
    sl = VGroup(tx(r"\tfrac14", 24, YELLOW).next_to(ax.c2p(0.25, 0), UP, buff=0.15),
                tx(r"-\tfrac14", 24, YELLOW).next_to(ax.c2p(-0.25, 0), UP, buff=0.15))
    s1 = mixed("主例：正方格子的 ", "$G(z;\\vec 0)$", font_size=28)
    s2 = tx(r"G(z;\vec 0)\approx-\frac{1}{\pi}\ln(1-16z^2)+\text{const}", 32)
    s3 = boxed(tx(r"C(N;\vec 0)\approx\frac{2\cdot4^N}{\pi N}", 36), color=YELLOW, buff=0.15)
    s4 = mixed("$N=10$", "：", "$" + nb["sqa"] + "$", " vs 精確 ", "$" + nb["sq"] + "$", font_size=28)
    hdr = [tx("N", 26, YELLOW)] + [str(N) for N in RN]
    tab = D.tex_table([hdr, [zh("精確／近似", font_size=22)] + nb["sqr"]], [2.0, 1.25, 1.25, 1.25, 1.25],
                      row_h=0.55, font_size=26)
    col = VGroup(s1, s2, s3, s4, tab).arrange(DOWN, aligned_edge=LEFT, buff=0.32)
    place_col(col, -1.2, 3.0)
    assert_in_frame(ax, sl, col)
    assert_apart(ax, col)
    s.say("回到主例：正方格子走回原點", wait=0.3)
    s.play(FadeIn(s1), FadeIn(ax), run_time=1.0)
    s.play(FadeIn(sps), FadeIn(sl), ShowCreation(cuts), run_time=1.2)
    s.say("$G(z;\\vec 0)$", " 在 ", "$z=\\pm\\tfrac14$", " 有對數奇點", wait=3.0)
    s.play(FadeIn(s2), run_time=1.0)
    s.play(FadeIn(s3), run_time=1.0)
    s.say("兩個奇點一樣遠：又是因子 2；對數給出 ", "$1/N$", wait=3.5)
    s.play(FadeIn(s4), FadeIn(tab), run_time=1.0)
    s.say("和 W206 的高斯極限一模一樣：", "$N=100$", " 時差 0.5%", wait=3.5)
    s.clear_stage()


def part_review(s):
    cards = [
        stage_card("生成函數", r"G(z)=\sum_N z^NC(N)", "W201、W202"),
        stage_card("遞迴 (2.9)", r"G=z\sum_iG(\vec w_i)+\delta", "W203"),
        stage_card("傅立葉", r"g(z;\vec q)=\frac{1}{1-z\chi(\vec q)}", "W204"),
        stage_card("取係數", r"C(N;\vec r)=\Big(\frac{a}{2\pi}\Big)^d\int_{BZ}e^{-i\vec q\cdot\vec r}\chi^N", "W205"),
        stage_card("高斯與普適性", r"\chi\approx Z\Big(1-\frac{q^2}{6}\Big)", "W206、W207、W208"),
        stage_card("擴散與格林函數", r"\frac{\partial P}{\partial N}=\frac{a^2}{2d}\nabla^2P", "W209"),
        stage_card("奇點讀係數", r"(z_c-z)^{-\alpha}\to\frac{N^{\alpha-1}}{\Gamma(\alpha)}z_c^{-N-\alpha}", "W212（補充 W210、W211）"),
    ]
    pos = [(-5.1, 1.55), (-1.7, 1.55), (1.7, 1.55), (5.1, 1.55), (5.1, -1.35), (1.7, -1.35), (-1.7, -1.35)]
    for c, (x, y) in zip(cards, pos):
        c.move_to(np.array([x, y, 0]))
    arrows = VGroup()
    for a, b in zip(cards, cards[1:]):
        if abs(a.get_y() - b.get_y()) < 0.1:
            st, en = (a[0].get_right(), b[0].get_left()) if b.get_x() > a.get_x() else (a[0].get_left(), b[0].get_right())
        else:
            st, en = a[0].get_bottom(), b[0].get_top()
        arrows.add(Arrow(st, en, buff=0.06, thickness=3.0).set_color(GREY_B))
    title = zh("第 2 章回顧", font_size=30).move_to(np.array([0, 2.95, 0]))
    assert_in_frame(title, *cards, arrows)
    assert_apart(title, VGroup(*cards[:4]))
    for c in cards:
        c.set_opacity(0.25)
    s.say("第 2 章回顧：一張地圖", wait=0.3)
    s.play(FadeIn(title), *[FadeIn(c) for c in cards], run_time=1.0)
    caps = [("生成函數：把整串數列裝進一個函數",), ("遞迴：最後一步來自鄰居",),
            ("傅立葉：鄰居的關係變成乘上 ", "$\\chi(\\vec q)$"), ("取係數＝對布里淵區積分",),
            ("$N$", " 很大：高斯；只看 ", "$q^2$", " 項，所以普適"), ("放大來看是擴散；", "$G$", " 是格林函數"),
            ("奇點的位置與形狀決定係數怎麼長",)]
    for k, c in enumerate(cards):
        anims = [c.animate.set_opacity(1.0)]
        if k:
            anims.append(GrowArrow(arrows[k - 1]))
        s.play(*anims, run_time=0.7)
        s.play(c[0].animate.set_stroke(YELLOW, 3), run_time=0.3)
        s.say(*caps[k], wait=2.0)
        s.play(c[0].animate.set_stroke(GREY_B, 2), run_time=0.3)


def part_close(s):
    s.say("奇點的位置決定指數成長，奇點的形狀決定冪次", wait=3.5)
    s.say("第 2 章完", wait=2.0)


class W212Singularities(ChuScene):
    def construct(self):
        self.show_header("W212　補充 3：從奇點讀出係數的成長＋第 2 章回顧")
        check_math()
        part_question(self)
        part_poles(self)
        part_branch(self)
        part_log(self)
        part_square(self)
        part_review(self)
        part_close(self)
        self.wait(1.0)
