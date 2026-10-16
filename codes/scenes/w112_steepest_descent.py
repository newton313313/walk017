"""
W112 最陡下降法（鞍點法）＋第 1 章回顧（§1.6，pp. 20–23，(1:S-1)–(1:S-15)，Fig. 1.5–1.6）

問：(1−z)^{−2} 的 z^{10} 係數是 11。能不能用「翻山越嶺」的方法把它算出來？
  1. C_N = (1/2πi)∮ e^{p(z)} dz，p(z) = −α ln(1−z) − (N+1) ln z（(1:S-3)，α = 2，N = 10）。
     兩條剖面：沿實軸 e^p 是谷，最低點 z* = (N+1)/(N+1+α) = 11/13 ≈ 0.846（書上近似 1 − α/(N+1) = 0.818）；
     沿通過 z* 的垂直線是峰 → 鞍點（Fig. 1.5 的馬鞍面只用字幕示意）。
  2. 圓周變形成穿過 z* 的圓（在 z* 垂直穿過：最陡下降路徑，Fig. 1.6）；p ≈ p(z*) + ½p''(z*)(z − z*)²（(1:S-5)）
     → 高斯積分 → C_N ≈ e^{p(z*)}/√(2π p''(z*)) = 265.4/√(2π·99.86) = 10.59（精確 11）。
  3. 一般公式（(1:S-15) 更正版）(1/√(2π)) α^{−α+1/2} e^α (N+1)^{α−1} → 11.46 = 1.0422 × 11；
     筆誤框：p. 23 前因子 1/(2π) → 1/√(2π)（照印的只有 4.57）；1.0422 = 2!/Stirling(2)（呼應 W106）：
     公式＝把精確答案 N^{α−1}/Γ(α) 的 Γ(α) 換成 Stirling。α = ½：精確 C(20,10)/4^10 = 0.1762，公式 0.1983。
  4. 第 1 章回顧：6 格地圖依序亮起。
畫面上的數字都由 utils.walk1d 算出，check_math() 核對。
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

N, A = 10, 2
ZS = W.sd_saddle_exact(N, A)                  # 11/13
X0 = 0.6                                      # 實軸剖面的橫軸從 0.6 開始（畫 u = x − 0.6）
SADDLE_COLOR = D.EXACT_COLOR


def ep_real(x):
    return math.exp(W.sd_p(x, N, A))


def ep_vertical(y):
    return math.exp(W.sd_p(complex(ZS, y), N, A).real)


def numbers():
    ps, p2 = W.sd_p(ZS, N, A), W.sd_p2(ZS, N, A)
    sad, ex = W.sd_full(N, A), W.sd_exact(N, A)
    return {"zs": f"{ZS:.3f}", "zb": f"{W.sd_saddle_book(N, A):.3f}", "ep": f"{math.exp(ps):.1f}",
            "p2": f"{p2:.2f}", "sad": f"{sad:.2f}", "exact": f"{ex:.0f}", "err": f"{(ex - sad) / ex * 100:.0f}",
            "approx": f"{W.sd_approx(N, A):.2f}", "ratio": f"{W.sd_approx(N, A) / ex:.4f}",
            "book": f"{W.sd_book(N, A):.2f}", "h_exact": f"{W.sd_exact(N, 0.5):.4f}",
            "h_approx": f"{W.sd_approx(N, 0.5):.4f}"}


def check_math():
    assert numbers() == {"zs": "0.846", "zb": "0.818", "ep": "265.4", "p2": "99.86", "sad": "10.59", "exact": "11",
                         "err": "4", "approx": "11.46", "ratio": "1.0422", "book": "4.57", "h_exact": "0.1762",
                         "h_approx": "0.1983"}
    assert abs(ZS - 11 / 13) < 1e-15 and abs(W.sd_exact(N, A) - 11) < 1e-9 and math.comb(11, 1) == 11
    assert f"{W.stirling_ratio(2):.4f}" == "1.0422"
    assert abs(W.sd_approx(N, A) / W.sd_exact(N, A) - W.stirling_ratio(2)) < 1e-12        # 正是 2!/Stirling(2)
    assert abs(W.sd_exact(N, 0.5) - math.comb(20, 10) / 4 ** 10) < 1e-12
    assert abs(math.exp(W.sd_p(ZS, N, A)) / math.sqrt(2 * math.pi * W.sd_p2(ZS, N, A)) - W.sd_full(N, A)) < 1e-12
    # 鞍點：實軸上是最低點、垂直方向是最高點
    for h in (0.01, 0.05):
        assert ep_real(ZS - h) > ep_real(ZS) < ep_real(ZS + h)
        assert ep_vertical(h) < ep_real(ZS) and ep_vertical(-h) < ep_real(ZS)
    # 兩張剖面圖的範圍
    assert max(ep_real(x) for x in np.linspace(0.615, 0.965, 200)) < 1600
    assert abs(ep_vertical(0) - ep_real(ZS)) < 1e-9 and ep_real(ZS) < 300
    # 1/Γ(α) 用 Stirling：α^{−α+1/2} e^α/√(2π)
    for a in (0.5, 2):
        assert abs(a / W.stirling(a) - a ** (-a + 0.5) * math.exp(a) / math.sqrt(2 * math.pi)) < 1e-12


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


def stage_card(title, tex, tag, w=3.9, h=1.9):
    t = title if isinstance(title, Mobject) else zh(title, font_size=26, color=YELLOW)
    f = fit(Tex(tex, font_size=28), w - 0.4)
    g = zh(tag, font_size=18, color=GREY_B)
    body = VGroup(t, f, g).arrange(DOWN, buff=0.18)
    frame = RoundedRectangle(width=w, height=h, corner_radius=0.15).set_stroke(GREY_B, 2).set_fill(BLACK, 0.5)
    frame.move_to(body)
    return VGroup(frame, body)


# ---------------------------------------------------------------- 各段
def part_question(s):
    s.say("問：", "$(1-z)^{-2}$", " 的 ", "$z^{10}$", " 係數是 11。能不能「翻山越嶺」算出來？", wait=3.5)


def part_profiles(s):
    nb = numbers()
    f1 = Tex(r"C_N=\frac{1}{2\pi i}\oint e^{p(z)}\,dz", font_size=34)
    f2 = eq(r"p(z)=-\alpha\ln(1-z)-(N+1)\ln z", r"\text{(1:S-3)}", font_size=34)
    f3 = mixed("$\\alpha=2,\\ N=10$", "：", "$(1-z)^{-2}=\\sum(N+1)z^N$", "，係數 ", "$11$", font_size=26)
    top = VGroup(VGroup(f1, f2).arrange(RIGHT, buff=0.8), f3).arrange(DOWN, buff=0.25)
    top.move_to(np.array([0, 3.1 - top.get_height() / 2, 0]))
    ax1 = D.plot_axes((0, 0.4, 0.1), (0, 1600, 500), width=5.0, height=2.2, x_ticks=[0.1, 0.2, 0.3],
                      x_tick_labels=["0.7", "0.8", "0.9"], y_ticks=[500, 1000, 1500], font_size=20, x_label="x")
    ax1.move_to(np.array([-3.3, -0.5, 0]))
    g1 = ax1.get_graph(lambda u: ep_real(X0 + u), x_range=[0.015, 0.365, 0.005]).set_stroke(D.BAR_COLOR, 3)
    dot1 = Dot(ax1.c2p(ZS - X0, ep_real(ZS)), radius=0.08).set_fill(SADDLE_COLOR, 1)
    lab1 = Tex(r"z^*\approx" + nb["zs"], font_size=26).set_color(SADDLE_COLOR).move_to(ax1.c2p(ZS - X0, 800))
    book = DashedLine(ax1.c2p(W.sd_saddle_book(N, A) - X0, 0), ax1.c2p(W.sd_saddle_book(N, A) - X0, 650))
    book.set_stroke(GREY_B, 2)
    ax2 = D.plot_axes((-0.3, 0.3, 0.1), (0, 300, 100), width=3.8, height=2.2, x_ticks=[-0.2, 0, 0.2],
                      y_ticks=[100, 200], font_size=20, x_label="y")
    ax2.move_to(np.array([3.6, -0.5, 0]))
    g2 = ax2.get_graph(ep_vertical, x_range=[-0.3, 0.3, 0.005]).set_stroke(D.GAUSS_COLOR, 3)
    dot2 = Dot(ax2.c2p(0, ep_vertical(0)), radius=0.08).set_fill(SADDLE_COLOR, 1)
    lab2 = Tex(nb["ep"], font_size=24).set_color(SADDLE_COLOR).next_to(dot2, RIGHT, buff=0.15)
    t1 = mixed("沿實軸：", "$e^{p(x)}$", " 是谷", font_size=26, color=D.BAR_COLOR).next_to(ax1, UP, buff=0.12)
    t2 = mixed("沿 ", "$z=z^*+iy$", "：", "$|e^{p}|$", " 是峰", font_size=26, color=D.GAUSS_COLOR)
    t2.next_to(ax2, UP, buff=0.12)
    bot = mixed("最低點 ", "$z^*=\\tfrac{N+1}{N+1+\\alpha}=\\tfrac{11}{13}\\approx" + nb["zs"] + "$",
                "（書上近似 ", "$1-\\tfrac{\\alpha}{N+1}=" + nb["zb"] + "$", "，灰色虛線）", font_size=24)
    fit(bot, 13.0).move_to(np.array([0, -2.55, 0]))
    assert_in_frame(top, ax1, ax2, t1, t2, lab1, lab2, bot)
    assert_apart(top, VGroup(t1, t2))
    assert_apart(ax1, ax2)
    assert_apart(bot, VGroup(ax1, ax2))

    s.play(Write(f1), run_time=1.2)
    s.play(FadeIn(f2), FadeIn(f3), run_time=1.0)
    s.say("把被積函數寫成 ", "$e^{p(z)}$", "：", "$p$", " 是「地形的高度」", wait=3.0)
    s.play(FadeIn(ax1), FadeIn(t1), ShowCreation(g1), run_time=1.5)
    s.say("沿著實軸走：先下坡、再上坡，是一個谷", wait=2.5)
    s.play(FadeIn(dot1), FadeIn(lab1), ShowCreation(book), FadeIn(bot), run_time=1.0)
    s.say("谷底 ", "$z^*=\\tfrac{11}{13}$", "；書上的近似 0.818 稍微偏左", wait=3.5)
    s.play(FadeIn(ax2), FadeIn(t2), ShowCreation(g2), run_time=1.5)
    s.play(FadeIn(dot2), FadeIn(lab2), run_time=0.6)
    s.say("同一點往垂直方向走：卻是一座峰", wait=3.0)
    s.say("一個方向是谷、另一個方向是峰：像馬鞍，叫鞍點", wait=3.5)
    s.clear_stage()


def part_deform(s):
    nb = numbers()
    cp = D.plot_axes((-1.2, 1.4, 0.5), (-1.2, 1.2, 0.5), width=4.55, height=4.2,
                     x_label=r"\mathrm{Re}\,z", y_label=r"\mathrm{Im}\,z")
    cp.move_to(np.array([-4.2, -0.3, 0]))
    small = D.contour_circle(cp, 0.35)
    big = D.contour_circle(cp, ZS)
    pole = D.pole_marker(cp, 1 + 0j, "z=1")
    o = Dot(cp.c2p(0, 0), radius=0.06).set_fill(WHITE, 1)
    sdot = Dot(cp.c2p(ZS, 0), radius=0.08).set_fill(SADDLE_COLOR, 1)
    slab = Tex("z^*", font_size=26).set_color(SADDLE_COLOR).next_to(sdot, LEFT, buff=0.12).shift(0.2 * UP)
    seg = Line(cp.c2p(ZS, -0.32), cp.c2p(ZS, 0.32)).set_stroke(YELLOW, 7)
    s1 = eq(r"p(z)\approx p(z^*)+\tfrac12\,p''(z^*)(z-z^*)^2", r"\text{(1:S-5)}", font_size=32)
    s2 = mixed("令 ", "$z=z^*+iy$", "：", "$e^{p}\\approx e^{p(z^*)}\\,e^{-p''(z^*)y^2/2}$", font_size=26)
    s3 = boxed(Tex(r"C_N\approx\frac{e^{p(z^*)}}{\sqrt{2\pi\,p''(z^*)}}", font_size=40), color=SADDLE_COLOR, buff=0.18)
    s4 = Tex(r"=\frac{" + nb["ep"] + r"}{\sqrt{2\pi\cdot" + nb["p2"] + r"}}=" + nb["sad"], font_size=36)
    s5 = mixed("精確答案 ", "$" + nb["exact"] + "$", "：誤差約 ", "$" + nb["err"] + "\\%$", font_size=28)
    col = VGroup(s1, s2, s3, s4, s5).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
    for m in col:
        fit(m, 7.3)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.35), -0.9, 2.9)
    assert_in_frame(cp, small, big, pole, slab, col)
    assert_apart(cp, col)
    s.say("Cauchy 的圓周可以自由變形：只要不碰到奇點", wait=0.3)
    s.play(FadeIn(cp), FadeIn(o), FadeIn(pole), ShowCreation(small[0]), FadeIn(small[1]), run_time=1.2)
    s.wait(1.5)
    s.play(Transform(small, big), FadeIn(sdot), FadeIn(slab), run_time=2.0)
    s.say("把圓放大到穿過 ", "$z^*$", "：在那裡正好垂直穿過鞍點", wait=3.0)
    s.play(ShowCreation(seg), run_time=0.8)
    s.say("沿這個方向，高度掉得最快：最陡下降路徑", wait=3.0)
    s.play(FadeIn(s1), run_time=0.8)
    s.say("在 ", "$z^*$", " 附近把 ", "$p$", " 展開到二次", wait=2.5)
    s.play(FadeIn(s2), run_time=0.8)
    s.say("垂直方向：", "$(z-z^*)^2=-y^2$", "，變成一個高斯", wait=3.0)
    s.play(FadeIn(s3, shift=0.2 * UP), run_time=1.0)
    s.say("高斯積分做完：只看山口附近", wait=2.5)
    s.play(FadeIn(s4), run_time=0.8)
    s.play(FadeIn(s5), run_time=0.6)
    s.say("得到 ", "$" + nb["sad"] + "$", "，精確答案是 11", wait=3.5)
    s.clear_stage()


def part_general(s):
    nb = numbers()
    g1 = eq(r"C_N\approx\frac{1}{\sqrt{2\pi}}\,\alpha^{-\alpha+1/2}e^{\alpha}(N+1)^{\alpha-1}", r"\text{(1:S-15)}",
            font_size=36)
    g1.move_to(np.array([0, 2.5, 0]))
    g2 = mixed("$\\alpha=2,\\ N=10$", "：", "$" + nb["approx"] + "=" + nb["ratio"] + "\\times11$", font_size=28)
    g3 = mixed("$" + nb["ratio"] + "=2!/\\mathrm{Stirling}(2)$", "：W106 表上的數", font_size=28)
    g4 = Tex(r"\frac{1}{\Gamma(\alpha)}=\frac{\alpha}{\alpha!}\approx\frac{\alpha}{\sqrt{2\pi\alpha}\,\alpha^{\alpha}e^{-\alpha}}"
             r"=\frac{\alpha^{-\alpha+1/2}e^{\alpha}}{\sqrt{2\pi}}", font_size=30)
    g5 = mixed("公式＝精確答案 ", "$\\tfrac{N^{\\alpha-1}}{\\Gamma(\\alpha)}$", " 的 ", "$\\Gamma(\\alpha)$",
               " 換成 Stirling", font_size=26)
    g6 = mixed("$\\alpha=\\tfrac12$", "：精確 ", "$\\tbinom{20}{10}/4^{10}=" + nb["h_exact"] + "$", "，公式 ",
               "$" + nb["h_approx"] + "$", font_size=26)
    col = VGroup(g2, g3, g4, g5, g6).arrange(DOWN, aligned_edge=LEFT, buff=0.32)
    for m in col:
        fit(m, 8.3)
    place_col(col.arrange(DOWN, aligned_edge=LEFT, buff=0.32), -6.4, 1.75)
    box = D.erratum_box("p. 23，(1:S-15)", r"\frac{1}{2\pi}", r"\frac{1}{\sqrt{2\pi}}", width=4.0)
    note = mixed("照印的算：只有 ", "$" + nb["book"] + "$", font_size=24, color=D.WRONG_COLOR)
    err = VGroup(box, note).arrange(DOWN, buff=0.2)
    place_col(err, 6.7 - err.get_width(), 1.6)
    assert_in_frame(g1, col, err)
    assert_apart(col, err)
    assert_apart(g1, err)
    s.say("一般的 ", "$(1-z)^{-\\alpha}$", "：同樣的鞍點法得到 (1:S-15)", wait=0.3)
    s.play(Write(g1[0]), FadeIn(g1[1]), run_time=1.8)
    s.wait(2.0)
    s.play(FadeIn(g2), run_time=0.8)
    s.say("$\\alpha=2$", "：", "$" + nb["approx"] + "$", "，比 11 多了 4%", wait=3.0)
    s.play(FadeIn(err, shift=0.2 * LEFT), run_time=0.8)
    s.say("書上前因子印成 ", "$\\tfrac{1}{2\\pi}$", "：照印的只有 ", "$" + nb["book"] + "$", wait=3.5)
    s.play(FadeIn(g3), run_time=0.8)
    s.say("多出來的 ", "$" + nb["ratio"] + "$", "：正是 W106 表上的 ", "$2!/\\mathrm{Stirling}(2)$", wait=3.5)
    s.play(FadeIn(g4), run_time=1.0)
    s.play(FadeIn(g5), run_time=0.8)
    s.say("公式其實是把 ", "$\\Gamma(\\alpha)$", " 換成了 Stirling 近似", wait=3.5)
    s.play(FadeIn(g6), run_time=0.8)
    s.say("$\\alpha=\\tfrac12$", "：精確 ", "$" + nb["h_exact"] + "$", "，公式 ", "$" + nb["h_approx"] + "$", wait=3.5)
    s.clear_stage()


def part_review(s):
    cards = [
        stage_card("數路線", r"C(N,d)=\binom{N}{(N+d)/2}", "W101、W102"),
        stage_card(mixed("動差與 ", "$\\sqrt N$", font_size=26, color=YELLOW), r"\overline{x^2}=Nl^2",
                   "W103、W104、W105"),
        stage_card("Stirling 與高斯", r"N!\approx\sqrt{2\pi N}\,N^Ne^{-N}", "W106、W107"),
        stage_card("遞迴與擴散", r"\frac{\partial P}{\partial N}=\frac{l^2}{2}\frac{\partial^2P}{\partial y^2}",
                   "W108、W109"),
        stage_card("傅立葉與生成函數", r"C(N,d)=\frac{1}{2\pi}\int_{-\pi}^{\pi}e^{-iqd}\chi^N\,dq", "W110、W111"),
        stage_card("鞍點法", r"C_N\approx\frac{e^{p(z^*)}}{\sqrt{2\pi\,p''(z^*)}}", "W112"),
    ]
    pos = [(-4.4, 1.5), (0, 1.5), (4.4, 1.5), (4.4, -1.3), (0, -1.3), (-4.4, -1.3)]     # 蛇形
    for c, (x, y) in zip(cards, pos):
        c.move_to(np.array([x, y, 0]))
    arrows = VGroup()
    for a, b in zip(cards, cards[1:]):
        va = b[0].get_center() - a[0].get_center()
        u = va / np.linalg.norm(va)
        if abs(u[0]) > 0.5:
            st, en = a[0].get_edge_center(RIGHT if u[0] > 0 else LEFT), b[0].get_edge_center(LEFT if u[0] > 0 else RIGHT)
        else:
            st, en = a[0].get_bottom(), b[0].get_top()
        arrows.add(Arrow(st, en, buff=0.08, thickness=3.0).set_color(GREY_B))
    title = zh("第 1 章回顧", font_size=30, color=WHITE).move_to(np.array([0, 3.0, 0]))
    assert_in_frame(title, *cards, arrows)
    for c in cards:
        c.set_opacity(0.25)
    s.say("第 1 章回顧：一張地圖", wait=0.3)
    s.play(FadeIn(title), *[FadeIn(c) for c in cards], run_time=1.0)
    caps = [("數路線：走法數＝二項式係數",), ("平均位置 0，但距離長成 ", "$\\sqrt N$"),
            ("Stirling 把二項分布變成高斯",), ("Pascal 遞迴放大就是擴散方程",),
            ("傅立葉與生成函數：數走法＝取係數",), ("鞍點法：積分只看山口附近",)]
    for k, c in enumerate(cards):
        anims = [c.animate.set_opacity(1.0)]
        if k:
            anims.append(GrowArrow(arrows[k - 1]))
        s.play(*anims, run_time=0.8)
        s.play(c[0].animate.set_stroke(YELLOW, 3), run_time=0.3)
        s.say(*caps[k], wait=2.2)
        s.play(c[0].animate.set_stroke(GREY_B, 2), run_time=0.3)


def part_close(s):
    s.say("鞍點法：", "$N$", " 很大時，積分只看山口附近", wait=3.0)
    s.say("第 1 章完", wait=2.0)


class W112SteepestDescent(ChuScene):
    def construct(self):
        self.show_header("W112　最陡下降法與第 1 章回顧")
        check_math()
        part_question(self)
        part_profiles(self)
        part_deform(self)
        part_general(self)
        part_review(self)
        part_close(self)
        self.wait(1.0)
