"""
W205 從 χ^N 取回走法數（§2.1.3，pp. 33–34，(2.14)–(2.17)）

問：把 (2cos q_x + 2cos q_y)^N 對一個方塊積分，為什麼會數出回到原點的走法？
  1. (2.14) 1/(1−zχ) = Σ (zχ)^n；z^N 的係數 c(N;q) = χ(q)^N（(2.15)(2.16)），c 是 C(N;r) 的傅立葉轉換。
  2. χ = e^{iq_x} + e^{−iq_x} + e^{iq_y} + e^{−iq_y}（E、W、N、S）；χ² 乘開 16 項＝16 種兩步走法，
     常數項 4 個（EW、WE、NS、SN）→ C(2;0) = 4；一般 e^{iq·r} 的係數＝C(N;r)（例 EN、NE → C(2;(1,1)) = 2）。
  3. (2.17) C(N;r) = (a/2π)^d ∫_BZ e^{−iq·r} χ^N d^dq；均分求積：正方 4, 36, 400, 4900, 63504；
     簡單立方 6, 90, 1860, 44730, 1172556；熱圖 (χ/4)^N（N = 2, 4, 10）：峰越來越尖，中心與四個角都是 1。
  4. 2(cos q_x + cos q_y) = 4cos((q_x+q_y)/2)cos((q_x−q_y)/2) → C(N;0) = C_1D(N,0)²：63504 = 252²；C(10;(1,1)) = 210×252。
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
from utils import lattice_draw as L
from utils import walk_draw as D
from utils import walk_gf as GF

NS = [2, 4, 6, 8, 10]
HEAT_N = [2, 4, 10]
FACT = {"E": r"e^{iq_x}", "N": r"e^{iq_y}", "W": r"e^{-iq_x}", "S": r"e^{-iq_y}"}
RET_COLOR = D.EXACT_COLOR


def step(ch):
    return GF.SQ_STEPS[GF.SQ_NAMES.index(ch)]


def numbers():
    return {"sq": [round(GF.fourier_count_sq(N, 0, 0)) for N in NS],
            "sc": [round(GF.fourier_count_sc(N)) for N in NS],
            "c11": GF.count_grid(2)[(1, 1)], "c00": GF.count_grid(2)[(0, 0)],
            "c1": GF.count1(10, 0), "c2": GF.count1(10, 2), "r10": GF.returns_sq(10), "r11": GF.count_sq(10, 1, 1)}


def check_math():
    nb = numbers()
    assert nb["sq"] == [4, 36, 400, 4900, 63504] == [GF.returns_sq(N) for N in NS]
    assert nb["sc"] == [6, 90, 1860, 44730, 1172556] == [GF.returns_sc(N) for N in NS]
    assert all(abs(GF.fourier_count_sq(N, 0, 0) - GF.returns_sq(N)) < 1e-6 for N in NS)
    assert all(abs(GF.fourier_count_sc(N) - GF.returns_sc(N)) < 1e-4 for N in NS)
    terms = GF.chi2_terms()
    assert len(terms) == 16 and sum(1 for t in terms if t[2]) == 4 == nb["c00"]
    assert {a + b for a, b, r in terms if r} == {"EW", "WE", "NS", "SN"}
    assert all(r == (step(a)[0] + step(b)[0] == 0 and step(a)[1] + step(b)[1] == 0) for a, b, r in terms)
    both = sorted(a + b for a, b, _ in terms if (step(a)[0] + step(b)[0], step(a)[1] + step(b)[1]) == (1, 1))
    assert both == ["EN", "NE"] and nb["c11"] == 2
    for qx, qy in [(0.3, -1.2), (2.0, 0.7), (-2.9, 3.1)]:
        assert abs(GF.chi_sq(qx, qy) - 4 * math.cos((qx + qy) / 2) * math.cos((qx - qy) / 2)) < 1e-12
    assert nb["r10"] == 63504 == nb["c1"] ** 2 == 252 ** 2
    assert nb["r11"] == 52920 == nb["c2"] * nb["c1"] == 210 * 252
    assert abs(GF.fourier_count_sq(10, 1, 1) - 52920) < 1e-6
    for N in HEAT_N:                                       # 偶數 N：中心與四個角都是 1
        assert abs((GF.chi_sq(0, 0) / 4) ** N - 1) < 1e-12 and abs((GF.chi_sq(math.pi, math.pi) / 4) ** N - 1) < 1e-12


# ---------------------------------------------------------------- 小工具
def eq(tex, tag=None, font_size=32, color=WHITE):
    m = Tex(tex, font_size=font_size).set_color(color)
    if tag is None:
        return m
    return VGroup(m, Tex(tag, font_size=22).set_color(D.MUTED)).arrange(RIGHT, buff=0.3)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def place_col(col, left, top):
    col.move_to(np.array([left + col.get_width() / 2, top - col.get_height() / 2, 0]))
    return col


# ---------------------------------------------------------------- 各段
def part_question(s):
    e = Tex(r"\frac{1}{(2\pi)^2}\int_{-\pi}^{\pi}\int_{-\pi}^{\pi}(2\cos q_x+2\cos q_y)^N\,dq_x\,dq_y=\ ?",
            font_size=44).move_to(np.array([0, 0.4, 0]))
    assert_in_frame(e)
    s.play(Write(e), run_time=1.8)
    s.say("問：把 ", "$(2\\cos q_x+2\\cos q_y)^N$", " 對一個方塊積分", wait=2.0)
    s.say("為什麼會數出回到原點的走法？", wait=3.5)
    s.clear_stage()


def part_series(s):
    e1 = eq(r"\frac{1}{1-z\chi}=1+z\chi+(z\chi)^2+\cdots+(z\chi)^n+\cdots", "(2.14)", font_size=38)
    e2 = eq(r"g(z;\vec q)=\sum_{N\ge0}z^N\,c(N;\vec q)", "(2.15)", font_size=38)
    e3 = boxed(eq(r"c(N;\vec q)=\chi(\vec q)^N", "(2.16)", font_size=44), color=D.GAUSS_COLOR, buff=0.2)
    e4 = mixed("而 ", "$g=\\sum_{\\vec r}G(z;\\vec r)\\,e^{i\\vec q\\cdot\\vec r}$", "，所以", font_size=26)
    e5 = Tex(r"\chi(\vec q)^N=\sum_{\vec r}C(N;\vec r)\,e^{i\vec q\cdot\vec r}", font_size=40).set_color(YELLOW)
    col = VGroup(e1, e2, e3, e4, e5).arrange(DOWN, buff=0.38).move_to(np.array([0, 0.15, 0]))
    assert_in_frame(col)
    s.say("回到 ", "$g=1/(1-z\\chi)$", "：展開成等比級數", wait=0.3)
    s.play(Write(e1[0]), FadeIn(e1[1]), run_time=1.8)
    s.wait(2.0)
    s.play(FadeIn(e2), run_time=0.8)
    s.say("$z^N$", " 的係數叫 ", "$c(N;\\vec q)$", wait=2.5)
    s.play(FadeIn(e3, shift=0.2 * UP), run_time=1.0)
    s.say("一比就知道：", "$c(N;\\vec q)=\\chi(\\vec q)^N$", wait=3.0)
    s.play(FadeIn(e4), FadeIn(e5), run_time=1.0)
    s.say("$\\chi^N$", " 乘開後，", "$e^{i\\vec q\\cdot\\vec r}$", " 的係數就是走法數", wait=4.0)
    s.clear_stage()


def part_chi2(s):
    pieces = [Tex(r"\chi=", font_size=40)]
    keys = ["E", "W", "N", "S"]
    for k, ch in enumerate(keys):
        if k:
            pieces.append(Tex("+", font_size=40))
        pieces.append(Tex(FACT[ch], font_size=40))
    top = VGroup(*pieces).arrange(RIGHT, buff=0.15).move_to(np.array([-2.6, 2.6, 0]))
    tlabs = VGroup(*[Tex(rf"\text{{{ch}}}", font_size=26).set_color(GREY_B).next_to(pieces[2 * k + 1], DOWN, buff=0.12)
                     for k, ch in enumerate(keys)])

    cs, gx0, gy0 = 0.9, 3.2, 0.75
    order = list(GF.SQ_NAMES)
    cells, texts = VGroup(), VGroup()
    ret_cells, en_cells = VGroup(), VGroup()
    for a, b, ret in GF.chi2_terms():
        i, j = order.index(a), order.index(b)
        c = np.array([gx0 + j * cs, gy0 - i * cs, 0])
        sq = RoundedRectangle(width=cs * 0.9, height=cs * 0.9, corner_radius=0.08).set_stroke(GREY_D, 1.5).move_to(c)
        t = Tex(rf"\text{{{a}{b}}}", font_size=26).move_to(c)
        cells.add(sq)
        texts.add(t)
        if ret:
            ret_cells.add(sq)
        if a + b in ("EN", "NE"):
            en_cells.add(sq)
    rhead = VGroup(*[Tex(FACT[ch], font_size=30).set_color(GREY_B).move_to(np.array([gx0 - 1.0, gy0 - i * cs, 0]))
                     for i, ch in enumerate(order)])
    chead = VGroup(*[Tex(FACT[ch], font_size=30).set_color(GREY_B).move_to(np.array([gx0 + j * cs, gy0 + 0.75, 0]))
                     for j, ch in enumerate(order)])
    gl = VGroup(zh("第一步", font_size=22, color=GREY_B).rotate(PI / 2).next_to(rhead, LEFT, buff=0.15),
                zh("第二步", font_size=22, color=GREY_B).next_to(chead, UP, buff=0.12))
    nb = numbers()
    l1 = mixed("$\\chi^2$", " 乘開：", "$4\\times4=16$", " 項，每項是「兩步」", font_size=24)
    l2 = zh("常數項（回到原點）4 個：EW、WE、NS、SN", font_size=24, color=RET_COLOR)
    l3 = Tex(rf"C(2;\vec 0)={nb['c00']}", font_size=38).set_color(RET_COLOR)
    l4 = mixed("一般：", "$e^{i\\vec q\\cdot\\vec r}$", " 的係數＝", "$C(N;\\vec r)$", font_size=24)
    l5 = mixed("例：", "$e^{i(q_x+q_y)}$", "：EN、NE，", f"$C(2;(1,1))={nb['c11']}$", font_size=24, color=YELLOW)
    col = VGroup(l1, l2, l3, l4, l5)
    for m in col:
        fit(m, 7.0)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.4), -6.5, 1.35)
    grid = VGroup(cells, texts, rhead, chead, gl)
    assert_in_frame(top, tlabs, grid, col)
    assert_apart(col, grid)
    assert_apart(VGroup(top, tlabs), grid)
    assert_apart(VGroup(top, tlabs), col)

    s.say("正方格子：", "$\\chi$", " 的四項＝四個方向", wait=0.3)
    s.play(Write(top), run_time=1.5)
    s.play(FadeIn(tlabs), run_time=0.6)
    s.say("$e^{iq_x}$", " 往東、", "$e^{-iq_x}$", " 往西、", "$e^{\\pm iq_y}$", " 往北、往南", wait=3.0)
    s.play(FadeIn(rhead), FadeIn(chead), FadeIn(gl), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(VGroup(c, t)) for c, t in zip(cells, texts)], lag_ratio=0.08), FadeIn(l1),
           run_time=2.0)
    s.say("$\\chi^2$", " 乘開 16 項：每一項是一種「兩步」走法", wait=3.0)
    s.play(*[c.animate.set_fill(RET_COLOR, 0.35) for c in ret_cells], FadeIn(l2), run_time=1.0)
    s.say("走回原點的是 EW、WE、NS、SN：相位抵消成常數", wait=3.0)
    s.play(FadeIn(l3), run_time=0.8)
    s.say("常數項係數 4 ＝ ", "$C(2;\\vec 0)$", wait=3.0)
    s.play(FadeIn(l4), run_time=0.8)
    s.play(*[c.animate.set_stroke(YELLOW, 3) for c in en_cells], FadeIn(l5), run_time=1.0)
    s.say("一般：", "$e^{i\\vec q\\cdot\\vec r}$", " 的係數就是走到 ", "$\\vec r$", " 的走法數", wait=3.5)
    s.clear_stage()


def part_inversion(s):
    nb = numbers()
    f1 = boxed(eq(r"C(N;\vec r)=\Big(\frac{a}{2\pi}\Big)^d\int_{BZ}e^{-i\vec q\cdot\vec r}\,\chi(\vec q)^N\,d^dq", "(2.17)",
                  font_size=40), color=D.GAUSS_COLOR, buff=0.2)
    f2 = mixed("布里淵區 BZ：", "$-\\pi/a<q_i\\le\\pi/a$", "；乘上 ", "$e^{-i\\vec q\\cdot\\vec r}$", " 再積分＝取出係數",
               font_size=24)
    rows = [[zh("N", font_size=24, color=YELLOW)] + [str(N) for N in NS],
            [zh("正方", font_size=24)] + [str(v) for v in nb["sq"]],
            [zh("簡單立方", font_size=24)] + [str(v) for v in nb["sc"]]]
    tab = D.tex_table(rows, [2.2, 1.2, 1.2, 1.5, 1.7, 2.2], row_h=0.6, font_size=32)
    for j in range(1, 6):
        D.table_cell(tab, 0, j, 6).set_color(YELLOW)
    rowm = [VGroup(*[D.table_cell(tab, i, j, 6) for j in range(6)]) for i in range(3)]
    f3 = zh("數值積分（均分求積，結果是精確的整數）", font_size=22, color=GREY_B)
    col = VGroup(f1, f2, tab, f3).arrange(DOWN, buff=0.45).move_to(np.array([0, 0.3, 0]))
    assert_in_frame(col)
    s.say("所以：取係數＝對布里淵區積分", wait=0.3)
    s.play(FadeIn(f1, shift=0.2 * UP), run_time=1.2)
    s.wait(2.0)
    s.play(FadeIn(f2), run_time=0.8)
    s.say("$e^{-i\\vec q\\cdot\\vec r}$", " 把 ", "$e^{i\\vec q\\cdot\\vec r}$", " 那一項變成常數，其餘積分是 0", wait=3.0)
    s.play(FadeIn(tab[1]), FadeIn(rowm[0]), FadeIn(f3), run_time=0.8)
    s.play(FadeIn(rowm[1]), run_time=0.8)
    s.say("正方格子：4、36、400、4900、63504", wait=3.0)
    s.play(FadeIn(rowm[2]), run_time=0.8)
    s.say("簡單立方：6、90、1860、44730、1172556", wait=3.5)
    s.clear_stage()


def part_heat(s):
    w, org = 4.5, np.array([-3.3, 0.05, 0])
    maps = [L.heatmap(lambda x, y, N=N: (GF.chi_sq(x, y) / 4) ** N, n=32, width=w, origin=org, vmax=1.0) for N in HEAT_N]
    frame = L.bz_frame(width=w, origin=org)
    nlab = [Tex(rf"\Big(\frac{{\chi}}{{4}}\Big)^{{{N}}}", font_size=40).set_color(YELLOW) for N in HEAT_N]
    for m in nlab:
        place_col(m, 1.3, 2.7)
    h1 = zh("峰越來越尖：寬度約 1/√N", font_size=24)
    h2 = mixed("中心 ", "$\\chi=4$", "，四個角 ", "$\\chi=-4$", font_size=24)
    h3 = mixed("偶數 ", "$N$", "：", "$(-4)^N=4^N$", "，角落也是 1", font_size=24, color=YELLOW)
    h4 = zh("（W206：角落的峰貢獻另一半）", font_size=22, color=GREY_B)
    col = VGroup(h1, h2, h3, h4).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    place_col(col, 1.3, nlab[0].get_bottom()[1] - 0.45)
    assert_in_frame(maps[0], frame, col, *nlab)
    assert_apart(frame, col)
    assert_apart(frame, nlab[0])

    hm = maps[0]
    s.say("被積函數長什麼樣？畫 ", "$(\\chi/4)^N$", wait=0.3)
    s.play(FadeIn(hm), FadeIn(frame), FadeIn(nlab[0]), run_time=1.2)
    s.wait(1.5)
    s.play(Transform(hm, maps[1]), FadeTransform(nlab[0], nlab[1]), run_time=1.5)
    s.play(FadeIn(h1), run_time=0.6)
    s.say("峰越來越尖：只剩中心附近有貢獻", wait=2.0)
    s.play(Transform(hm, maps[2]), FadeTransform(nlab[1], nlab[2]), run_time=1.5)
    s.say("$N=2,4,10$", "：峰越來越尖", wait=2.5)
    s.play(FadeIn(h2), FadeIn(h3), run_time=0.8)
    s.say("中心是 1，四個角也是 1：", "$(-4)^N=4^N$", wait=3.0)
    s.play(FadeIn(h4), run_time=0.6)
    s.wait(2.0)
    s.clear_stage()


def part_factor(s):
    nb = numbers()
    k1 = Tex(r"2(\cos q_x+\cos q_y)=4\cos\frac{q_x+q_y}{2}\,\cos\frac{q_x-q_y}{2}", font_size=40)
    k2 = mixed("換變數 ", "$u=x+y,\\ v=x-y$", "：兩枚獨立的硬幣（W203 的彩蛋）", font_size=26)
    k3 = Tex(rf"C(N;\vec 0)=C_{{1D}}(N,0)^2:\quad {nb['r10']}={nb['c1']}^2", font_size=40)
    k4 = Tex(rf"C(10;(1,1))={nb['c2']}\times{nb['c1']}={nb['r11']}", font_size=40)
    k5 = mixed("（", "$u=2$", " 的硬幣 ", f"$C_{{1D}}(10,2)={nb['c2']}$", "，", "$v=0$", f" 的硬幣 {nb['c1']}）",
               font_size=24, color=GREY_B)
    col = VGroup(k1, k2, k3, k4, k5).arrange(DOWN, buff=0.45).move_to(np.array([0, 0.25, 0]))
    k3.set_color(YELLOW)
    assert_in_frame(col)
    s.say("正方格子還有一個小秘密：", "$\\chi$", " 可以因式分解", wait=0.3)
    s.play(Write(k1), run_time=1.8)
    s.wait(1.5)
    s.play(FadeIn(k2), run_time=0.8)
    s.say("轉 45 度看：", "$x+y$", " 與 ", "$x-y$", " 各自是一維的擲幣", wait=3.0)
    s.play(FadeIn(k3), run_time=0.8)
    s.say("所以 ", "$63504=252^2$", "：兩個一維的回家走法相乘", wait=3.0)
    s.play(FadeIn(k4), FadeIn(k5), run_time=0.8)
    s.say("$C(10;(1,1))=210\\times252=52920$", wait=3.5)
    s.clear_stage()


def part_close(s):
    s.say("取係數＝對布里淵區積分", wait=3.0)
    s.say("主例：", "$1,4,36,400,4900,63504$", wait=3.5)


class W205FourierInversion(ChuScene):
    def construct(self):
        self.show_header("W205　從 χ 的 N 次方取回走法數")
        check_math()
        part_question(self)
        part_series(self)
        part_chi2(self)
        part_inversion(self)
        part_heat(self)
        part_factor(self)
        part_close(self)
        self.wait(1.0)
