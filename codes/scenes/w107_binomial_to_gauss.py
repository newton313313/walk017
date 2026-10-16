"""
W107 二項分布變成高斯（§1.3，pp. 12–13，(1.37)–(1.41)）

問：N 越來越大，二項分布會變成什麼形狀？
  1. 長條 P(N,d)：N = 10, 20, 40, 100 依序變形（同一個 d 座標）；表格列出峰高 P(N,0) 與 √N·P(N,0)。
  2. 四條輪廓疊在一起，換成 u = d/√N、縱軸乘 √N → 幾乎重合（極限 2/√(2π)·e^{−u²/2}）。
  3. Stirling：n = N/2 + m → ln C(N,n) ≈ N ln2 − 2m²/N − ½ ln(πN/2)；核對 N = 100；筆誤框 (1.38)。
  4. 奇偶性與因子 2：寬 2 的柱子面積 ≈ 曲線在 [d−1, d+1] 下的面積 → 柱高 = 2 × 密度；筆誤框 (1.40)。
  5. 中央極限定理（第 2 章「普適性」預告）。
畫面上的數字都由 utils.walk1d 算出，check_math() 核對（目標值＝場景卡上的數字）。
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
from utils import walk1d as W
from utils import walk_draw as D

NS = (10, 20, 40, 100)
N_COLORS = {10: TEAL_B, 20: GREEN_B, 40: BLUE_B, 100: PURPLE_B}
DMAX = 32                                  # 長條只畫 |d| ≤ 32（更遠的機率合計 < 0.001，看不見）
O2 = np.array([-1.3, -2.3, 0])             # 第 1–2 步的圖：原點
UX, UY = 0.14, 18.0                        # d 座標
UXU, UYU = 1.4, 5.0                        # u = d/√N 座標（縱軸是 √N·P）
O4 = np.array([-3.0, -2.3, 0])             # 第 4 步（N = 10）的圖
UX4, UY4 = 0.33, 16.0


# ---------------------------------------------------------------- 數字（畫面上的字串都從這裡來）
def peak_rows():
    out = []
    for N in NS:
        p0 = float(W.walk_prob(N, 0))
        out.append((N, f"{p0:.4f}", f"{math.sqrt(N) * p0:.3f}"))
    return out


def numbers():
    lnC = W.ln_factorial(100) - 2 * W.ln_factorial(50)
    lnA = W.ln_binom_approx(100, 0)
    return {
        "lim": f"{2 / math.sqrt(2 * math.pi):.3f}",
        "lnC": f"{lnC:.4f}", "lnA": f"{lnA:.4f}", "gap": f"{lnA - lnC:.4f}",
        "P0": f"{float(W.walk_prob(10, 0)):.4f}", "rho0": f"{W.gauss_density(0, 10):.5f}",
        "g0": f"{W.gauss_prob(10, 0):.4f}",
        "P4": f"{float(W.walk_prob(10, 4)):.4f}", "g4": f"{W.gauss_prob(10, 4):.4f}",
    }


def simpson(f, a, b, n=200):
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) + sum((4 if k % 2 else 2) * f(a + k * h) for k in range(1, n)))


def check_math():
    assert peak_rows() == [(10, "0.2461", "0.778"), (20, "0.1762", "0.788"),
                           (40, "0.1254", "0.793"), (100, "0.0796", "0.796")]
    assert numbers() == {"lim": "0.798", "lnC": "66.7838", "lnA": "66.7863", "gap": "0.0025",
                         "P0": "0.2461", "rho0": "0.12616", "g0": "0.2523", "P4": "0.1172", "g4": "0.1134"}
    for N in NS:                                            # 沒畫出來的長條合計 < 0.001
        assert float(sum(W.walk_prob(N, d) for d in W.distribution(N) if abs(d) > DMAX)) < 1e-3
    assert abs(math.sqrt(1000) * float(W.walk_prob(1000, 0)) - 2 / math.sqrt(2 * math.pi)) < 1e-3
    for d in (0, 4):                                        # 寬 2 的柱子面積 ≈ 曲線在 [d−1, d+1] 下的面積
        area = simpson(lambda u: W.gauss_density(u, 10), d - 1, d + 1)
        assert abs(area - float(W.walk_prob(10, d))) < 0.01
    assert abs(simpson(lambda u: 2 * W.gauss_density(u, 1), -8, 8) - 2) < 1e-9   # 曲線 2φ(u) 下面積 2


# ---------------------------------------------------------------- 小工具
def eq(tex, tag=None, font_size=32, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return m
    return VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


def bars_for(N):
    vals = {d: float(W.walk_prob(N, d)) for d in range(-DMAX, DMAX + 1, 2)}
    bars = D.bar_chart(vals, ux=UX, uy=UY, origin=O2, width=1.6)[0]
    for r in bars:
        if r.get_height() < 0.01:          # 到不了或太小：透明（保留位置，Transform 才能一對一）
            r.set_fill(opacity=0).set_stroke(opacity=0)
    return bars


def profile(N, rescaled):
    dm = 2 * int(3.2 * math.sqrt(N) / 2 + 1e-9)            # |u| ≤ 3.2
    pts = []
    for d in range(-dm, dm + 1, 2):
        p = float(W.walk_prob(N, d))
        if rescaled:
            pts.append(O2 + d / math.sqrt(N) * UXU * RIGHT + math.sqrt(N) * p * UYU * UP)
        else:
            pts.append(O2 + d * UX * RIGHT + p * UY * UP)
    line = VMobject().set_points_as_corners(pts).set_stroke(N_COLORS[N], 3)
    dots = VGroup(*[Dot(q, radius=0.045).set_fill(N_COLORS[N], 1) for q in pts])
    return VGroup(line, dots)


def axis_parts():
    line = Line(O2 + (-DMAX - 1.5) * UX * RIGHT, O2 + (DMAX + 1.5) * UX * RIGHT).set_stroke(D.MUTED, 2)
    dl = VGroup(*[Tex(str(t), font_size=22).set_color(D.MUTED).next_to(O2 + t * UX * RIGHT, DOWN, buff=0.12)
                  for t in range(-30, 31, 10)])
    dl.add(Tex("d", font_size=30).set_color(D.MUTED).next_to(line, RIGHT, buff=0.15))
    ul = VGroup(*[Tex(str(t), font_size=22).set_color(D.MUTED).next_to(O2 + t * UXU * RIGHT, DOWN, buff=0.12)
                  for t in range(-3, 4)])
    ul.add(Tex(r"u=d/\sqrt N", font_size=28).set_color(D.MUTED).next_to(line, RIGHT, buff=0.15))
    return line, dl, ul


def peak_table():
    rows = [["N", "P(N,0)", r"\sqrt N\,P(N,0)"]] + [[str(N), p, q] for N, p, q in peak_rows()]
    t = D.tex_table(rows, [0.8, 1.35, 1.75], row_h=0.46, font_size=26)
    t.move_to(np.array([4.8, 1.95, 0]))
    return t


def area_under(fn, a, b, ux, uy, origin, color, n=40):
    xs = np.linspace(a, b, n)
    pts = ([origin + a * ux * RIGHT] + [origin + x * ux * RIGHT + fn(x) * uy * UP for x in xs]
           + [origin + b * ux * RIGHT])
    return Polygon(*pts).set_fill(color, 0.55).set_stroke(color, 0)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：步數 ", "$N$", " 越來越大，二項分布會變成什麼形狀？", wait=2.5)


def part_shapes(s):
    line, dl, ul = axis_parts()
    bars = bars_for(10)
    lab = Tex("N=10", font_size=40).move_to(np.array([-4.8, 2.4, 0]))
    table = peak_table()
    row = lambda i: VGroup(*[D.table_cell(table, i, j, 3) for j in range(3)])
    assert_in_frame(line, dl, ul, bars, lab, table)
    s.play(ShowCreation(line), FadeIn(dl), run_time=0.8)
    s.play(LaggedStartMap(GrowFromEdge, bars, edge=DOWN, lag_ratio=0.03), FadeIn(lab), run_time=1.5)
    s.play(ShowCreation(table[1]), FadeIn(row(0)), FadeIn(row(1)), run_time=1.0)
    s.say("先看 ", "$N=10$", "：峰在 ", "$d=0$", "，高 ", f"$P(10,0)={peak_rows()[0][1]}$", wait=3.0)
    s.say("步數加倍、再加倍：看長條怎麼變", wait=0.5)
    for i, N in enumerate(NS[1:], start=2):
        s.play(Transform(bars, bars_for(N)), Transform(lab, Tex(f"N={N}", font_size=40).move_to(lab)),
               run_time=1.6)
        s.play(FadeIn(row(i)), run_time=0.5)
        s.wait(1.0)
    s.say("同一個 ", "$d$", " 座標：越來越寬、越來越矮", wait=3.0)
    lim = Tex(r"\sqrt N\,P(N,0)\to\frac{2}{\sqrt{2\pi}}=" + numbers()["lim"], font_size=28)
    lim.next_to(table, DOWN, buff=0.3)
    assert_in_frame(lim)
    s.play(Write(lim), run_time=1.2)
    s.say("峰高乘上 ", "$\\sqrt N$", " 幾乎不變：高度 ", "$\\propto 1/\\sqrt N$", wait=3.5)

    profs = {N: profile(N, False) for N in NS}
    legend = VGroup(*[Tex(f"N={N}", font_size=28).set_color(N_COLORS[N]) for N in NS])
    legend.arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to(np.array([-5.75, 2.1, 0]))
    s.play(FadeOut(bars), FadeOut(lab), run_time=0.6)
    s.play(*[ShowCreation(profs[N][0]) for N in NS], *[FadeIn(profs[N][1]) for N in NS], FadeIn(legend),
           run_time=2.0)
    s.say("四個 ", "$N$", " 的輪廓疊在一起：寬窄高矮都不同", wait=3.0)
    s.say("橫軸換成 ", "$u=d/\\sqrt N$", "，縱軸乘上 ", "$\\sqrt N$", wait=0.3)
    s.play(*[Transform(profs[N], profile(N, True)) for N in NS], FadeOut(dl), FadeIn(ul), run_time=2.5)
    s.say("四張圖幾乎重合！", wait=2.5)
    gc = D.curve(lambda u: 2 * W.gauss_density(u, 1), -3.3, 3.3, ux=UXU, uy=UYU, origin=O2, dashed=True)
    glab = Tex(r"\frac{2}{\sqrt{2\pi}}e^{-u^2/2}", font_size=32).set_color(D.GAUSS_COLOR)
    glab.move_to(np.array([-3.7, 1.6, 0]))
    assert_in_frame(gc, glab, legend)
    assert_apart(glab, legend)
    s.play(ShowCreation(gc), FadeIn(glab), run_time=1.5)
    s.say("極限是一條鐘形曲線（前面那個 2 等一下解釋）", wait=3.5)
    s.clear_stage()


def part_stirling(s):
    nb = numbers()
    l1 = mixed("令 ", "$n=\\tfrac N2+m$", "，", "$d=2m$", "（", "$m$", "：正面數比一半多出的量）", font_size=26)
    l2 = eq(r"\ln\Big(1\pm\frac{2m}{N}\Big)\approx\pm\frac{2m}{N}-\frac12\Big(\frac{2m}{N}\Big)^2", "(1.38)")
    l3 = eq(r"\ln\binom{N}{n}\approx N\ln2-\frac{2m^2}{N}-\frac12\ln\frac{\pi N}{2}", "(1.39)")
    l4 = boxed(Tex(r"\binom{N}{n}\approx2^N\sqrt{\frac{2}{\pi N}}\,e^{-d^2/2N}", font_size=36),
               color=D.GAUSS_COLOR, buff=0.18)
    chk = mixed("核對 ", "$N=100,\\ d=0$", "：", f"$\\ln\\tbinom{{100}}{{50}}={nb['lnC']}$", "，近似 ",
                f"${nb['lnA']}$", font_size=26)
    col = place_col(VGroup(l1, l2, l3, l4, chk).arrange(DOWN, aligned_edge=LEFT, buff=0.4), -6.4, 3.0)
    box = D.erratum_box("p. 12，(1.38)", r"+\frac12\Big(\frac{2m}{N}\Big)^2", r"-\frac12\Big(\frac{2m}{N}\Big)^2",
                        width=4.3)
    box.move_to(np.array([6.7 - box.get_width() / 2, l3.get_center()[1], 0]))
    assert_in_frame(col, box)
    for m in (l2, l3, l4, chk):
        assert_apart(box, m)

    s.say("怎麼證明？用 Stirling 公式展開 ", "$\\ln\\tbinom Nn$", wait=0.3)
    s.play(FadeIn(l1, shift=0.2 * DOWN), run_time=0.8)
    s.wait(2.5)
    s.say("關鍵一步：", "$\\ln(1\\pm x)$", " 展開到二次項", wait=0.3)
    s.play(Write(l2[0]), FadeIn(l2[1]), run_time=1.5)
    s.wait(2.5)
    s.say("書上 (1.38) 二次項的正負號印反了", wait=0.3)
    s.play(FadeIn(box, shift=0.2 * LEFT), run_time=0.8)
    s.wait(3.0)
    s.say("一次項互相抵消，只剩 ", "$-2m^2/N$", wait=0.3)
    s.play(Write(l3[0]), FadeIn(l3[1]), run_time=1.5)
    s.wait(3.0)
    s.say("取指數、代回 ", "$m=d/2$", "：一條高斯曲線", wait=0.3)
    s.play(FadeIn(l4, shift=0.2 * UP), run_time=1.0)
    s.wait(3.0)
    s.play(FadeIn(chk), run_time=0.8)
    s.say("核對 ", "$N=100,\\ d=0$", "：只差 ", f"${nb['gap']}$", wait=3.5)
    s.clear_stage()


def part_factor_two(s):
    nb = numbers()
    vals = {d: float(W.walk_prob(10, d)) for d in range(-10, 11, 2)}
    narrow = D.bar_chart(vals, ux=UX4, uy=UY4, origin=O4, width=0.9)[0]
    wide = D.bar_chart({d: v / 2 for d, v in vals.items()}, ux=UX4, uy=UY4, origin=O4, width=2.0,
                       fill_opacity=0.45)[0]
    axis = Line(O4 + 11.4 * UX4 * LEFT, O4 + 11.4 * UX4 * RIGHT).set_stroke(D.MUTED, 2)
    nums = VGroup(*[Tex(str(d), font_size=20).set_color(D.MUTED).next_to(O4 + d * UX4 * RIGHT, DOWN, buff=0.1)
                    for d in range(-10, 11, 2)])
    odd = VGroup(*[Circle(radius=0.055).set_stroke(GREY_C, 2).move_to(O4 + d * UX4 * RIGHT)
                   for d in range(-9, 10, 2)])
    odd_lab = VGroup(Circle(radius=0.07).set_stroke(GREY_C, 2), zh("奇數位置：到不了", font_size=22, color=GREY_B))
    odd_lab.arrange(RIGHT, buff=0.15).move_to(np.array([-5.0, 2.5, 0]))
    rho = D.curve(lambda u: W.gauss_density(u, 10), -11, 11, ux=UX4, uy=UY4, origin=O4)
    rho_lab = Tex(r"\rho(d)", font_size=30).set_color(D.GAUSS_COLOR).move_to(O4 + 5.5 * UX4 * RIGHT + 1.2 * UP)

    t1 = mixed("只有偶數 ", "$d$", " 到得了：柱子間隔 2", font_size=26)
    t2 = eq(r"\rho(d)=\frac{1}{\sqrt{2\pi N}}\,e^{-d^2/2Nl^2}", "(1.41)", font_size=30, color=D.GAUSS_COLOR)
    t3 = Tex(f"P(10,0)={nb['P0']}", font_size=30)
    t3b = Tex(rf"2\rho(0)=2\times{nb['rho0']}={nb['g0']}", font_size=30).set_color(D.GAUSS_COLOR)
    t4 = Tex(f"P(10,4)={nb['P4']}", font_size=30)
    t4b = Tex(rf"2\rho(4)={nb['g4']}", font_size=30).set_color(D.GAUSS_COLOR)
    col = place_col(VGroup(t1, t2, t3, t3b, t4, t4b).arrange(DOWN, aligned_edge=LEFT, buff=0.3), 1.3, 3.0)
    box = D.erratum_box("(1.40)：少了因子 2", r"\tfrac{1}{\sqrt{2\pi N}}\,e^{-d^2/2N}",
                        r"\tfrac{2}{\sqrt{2\pi N}}\,e^{-d^2/2N}", width=5.2, font_size=26)
    box.move_to(np.array([1.3 + box.get_width() / 2, -2.95 + box.get_height() / 2, 0]))
    shade = area_under(lambda u: W.gauss_density(u, 10), 3, 5, UX4, UY4, O4, D.GAUSS_COLOR)
    assert_in_frame(axis, nums, narrow, odd_lab, rho, rho_lab, col, box)
    assert_apart(col, box)
    assert_apart(axis, col)

    s.say("回到 ", "$N=10$", "：極限曲線前面為什麼有個 2？", wait=0.3)
    s.play(ShowCreation(axis), FadeIn(nums), run_time=0.8)
    s.play(LaggedStartMap(GrowFromEdge, narrow, edge=DOWN, lag_ratio=0.05), run_time=1.2)
    s.play(FadeIn(odd), FadeIn(odd_lab), FadeIn(t1), run_time=0.8)
    s.say("奇偶性：", "$d$", " 一定和 ", "$N$", " 同奇偶，奇數位置一個人也沒有", wait=3.0)
    s.play(ShowCreation(rho), FadeIn(rho_lab), FadeIn(t2), run_time=1.5)
    s.say("疊上密度曲線 ", "$\\rho(d)$", "（(1.41)）：柱子大約高了一倍", wait=2.5)
    s.play(FadeIn(t3), FadeIn(t3b), run_time=0.8)
    s.wait(2.5)
    s.say("原因：每根柱子其實佔了寬度 2", wait=0.3)
    s.play(Transform(narrow, wide), run_time=2.0)
    s.say("攤平成寬 2、高度減半：柱子貼住了曲線", wait=2.5)
    s.play(narrow[7].animate.set_stroke(YELLOW, 3), FadeIn(shade), run_time=1.0)
    s.say("寬 2 的柱子面積 ", "$\\approx$", " 曲線在 ", "$[d-1,d+1]$", " 下的面積", wait=3.5)
    s.play(FadeIn(t4), FadeIn(t4b), run_time=0.8)
    s.say("所以 柱高 ", "$=2\\times$", " 密度：", f"${nb['P4']}$", " 對 ", f"${nb['g4']}$", wait=3.5)
    s.play(FadeIn(box, shift=0.2 * UP), run_time=0.8)
    s.say("書上 (1.40) 正好漏了這個因子 2", wait=3.5)
    s.clear_stage()


def part_clt(s):
    title = zh("中央極限定理", font_size=40, color=YELLOW).move_to(np.array([0, 2.6, 0]))
    items = VGroup(
        mixed("每步 ", "$\\pm1$", "（本章的硬幣）", font_size=26),
        mixed("每步 ", "$-1,\\ 0,\\ +1$", "（偶爾停下來）", font_size=26),
        mixed("每步長度隨機（變異數有限）", font_size=26),
    ).arrange(DOWN, aligned_edge=LEFT, buff=0.55)
    place_col(items, -6.3, 1.5)
    o = np.array([3.6, -1.4, 0])
    base = Line(o + 3.0 * LEFT, o + 3.0 * RIGHT).set_stroke(D.MUTED, 2)
    bell = D.curve(lambda u: math.exp(-u * u / 2), -3.2, 3.2, ux=0.85, uy=2.6, origin=o)
    arrows = VGroup(*[Arrow(it.get_right() + 0.25 * RIGHT, o + 1.5 * LEFT + 1.5 * UP, buff=0.1, thickness=2.0)
                      .set_color(GREY_B) for it in items])
    tag = boxed(mixed("第 2 章：", "普適性", font_size=28), color=TEAL, buff=0.2).move_to(np.array([3.6, -2.45, 0]))
    assert_in_frame(title, items, base, bell, arrows, tag)
    s.play(FadeIn(title, shift=0.2 * DOWN), run_time=0.8)
    s.play(LaggedStartMap(FadeIn, items, shift=0.2 * RIGHT, lag_ratio=0.4), run_time=1.8)
    s.play(ShowCreation(base), LaggedStartMap(GrowArrow, arrows, lag_ratio=0.3), run_time=1.2)
    s.play(ShowCreation(bell), run_time=1.2)
    s.say("細節不同的漫步，放大來看都是這條鐘形曲線", wait=3.5)
    s.play(FadeIn(tag), run_time=0.6)
    s.say("這就是中央極限定理：第 2 章「普適性」的預告", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("二項分布 ", "$\\to$", " 高斯：寬度 ", "$\\sqrt N$", "、高度 ", "$1/\\sqrt N$", wait=3.0)
    s.say("而且柱子的高度是密度的兩倍", wait=2.5)


class W107BinomialToGauss(ChuScene):
    def construct(self):
        self.show_header("W107　二項分布變成高斯")
        check_math()
        part_question(self)
        part_shapes(self)
        part_stirling(self)
        part_factor_two(self)
        part_clt(self)
        part_close(self)
        self.wait(1.0)
