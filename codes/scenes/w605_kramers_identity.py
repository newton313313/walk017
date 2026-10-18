"""
W605 Kramers 恆等式：形狀是步伐的二次式（§6.2.1，pp. 138–142，(6.21)–(6.37)，Ex. 6.5–6.6；補充 2–3 只取結論）

問：T 是 N+1 個點的平均。能不能直接寫成「步伐乘步伐」？
  1. (6.21)(6.22)：T_ij = Σ a_{αβ} η_{α,i} η_{β,j}；N = 4 的 a（× 1/25，SH.kramers）；補充 3：整數＝點對數 m(N+1−n)（SH.kramers_count）；
     代入「東東北西」→ (1/25)[[14,3],[3,6]]（SH.tensor_kramers_frac）；筆誤框 1 (6.22) α<β → α≤β（SH.kramers_book）。
  2. 正弦模態 (6.24)–(6.27)：sin(kπα/(N+1))（SH.sine_mode）；特徵值 vs (N+1)/(π²k²)（SH.eigs_exact、SH.eigs_asym）；Tr a = 4/5。
  3. 高斯步伐＋配對 (6.28)–(6.37)（SH.moments_exact）：N²/20d²、N²/36d²、N²/90d² → 9/5、9/2 → (4+2d)/(4+5d)；筆誤框 2 (6.35b)。
  4. 對照表（d = 3）：SH.A_exact、SH.sim_gauss_moments。
畫面上的數字都由 utils.walk_shape 算出，check_math() 核對。
"""
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_shape as SH
from utils import walk_draw as D

F = Fraction
STEPS = SH.SMALL_STEPS
MODE_COLORS = [BLUE_C, TEAL_C, ORANGE, PINK]
STEP_COLORS = [YELLOW, YELLOW, GREEN_C, YELLOW]


def f4(x):
    return f"{x:.4f}"


# ---------------------------------------------------------------- 數學
@lru_cache(maxsize=None)
def numbers():
    nb = {}
    nb["a25"] = [[int(v * 25) for v in row] for row in SH.kramers(4)]
    nb["count"] = SH.kramers_count(4)
    nb["Tk"] = SH.tensor_kramers_frac(STEPS)
    nb["Tdef"] = SH.small_example()["T"]
    eta = np.array(STEPS, dtype=float)
    nb["Tbook"] = eta.T @ SH.kramers_book(4, "nodiag") @ eta
    nb["eig"] = SH.eigs_exact(4)
    nb["asym"] = [SH.eigs_asym(4, k) for k in range(1, 5)]
    nb["modes"] = [SH.sine_mode(4, k) for k in range(1, 5)]
    nb["tra"] = SH.trace_a(4)
    nb["mex"] = [float(v) for v in SH.moments_exact(100, 3)]
    nb["sim"] = SH.sim_gauss_moments()[:3]
    nb["A4"] = float(SH.A_exact(4, 3))
    nb["A100"] = float(SH.A_exact(100, 3))
    nb["Ainf"] = SH.A_d(3)
    nb["Asim"] = SH.A_from_moments(3, *nb["sim"])
    return nb


def check_math():
    nb = numbers()
    a = nb["a25"]
    assert a == [[4, 3, 2, 1], [3, 6, 4, 2], [2, 4, 6, 3], [1, 2, 3, 4]]
    assert nb["count"] == a                                                 # 補充 3：整數＝點對數
    assert all(nb["count"][m - 1][n - 1] == m * (5 - n) for m in range(1, 5) for n in range(m, 5))
    assert nb["count"][1][2] == 4 and len([j for j in range(5) if j < 2]) * len([k for k in range(5) if k >= 3]) == 4
    assert nb["Tk"] == nb["Tdef"] == [[F(14, 25), F(3, 25)], [F(3, 25), F(6, 25)]]
    assert a[0][0] + a[1][1] + a[3][3] + 2 * a[0][1] - 2 * a[0][3] - 2 * a[1][3] == 14
    assert a[0][2] + a[1][2] - a[3][2] == 3 and a[2][2] == 6
    assert (4 + 6 + 4 + 2 * 3 - 2 * 1 - 2 * 2, 2 + 4 - 3) == (14, 3)
    assert abs(nb["Tbook"][0][0]) < 1e-12 and abs(nb["Tbook"][1][1]) < 1e-12      # 照書上 (6.22)：T11 = T22 = 0
    assert [f4(v) for v in nb["eig"]] == ["0.5236", "0.1447", "0.0764", "0.0553"]
    assert [f4(v) for v in nb["asym"]] == ["0.5066", "0.1267", "0.0563", "0.0317"]
    an = SH.kramers_np(4)
    for k in range(4):                                                      # 正弦向量精確是特徵向量
        v = nb["modes"][k]
        assert np.allclose(an @ v, nb["eig"][k] * v, atol=1e-13)
    assert nb["tra"] == F(4, 5) and abs(sum(nb["eig"]) - 0.8) < 1e-12
    big = 10**6                                                             # 大 N：N/6、N²/90 → N²/20d² …
    assert abs(float(SH.trace_a(big)) / big - 1 / 6) < 1e-6 and abs(float(SH.trace_a2(big)) / big**2 - 1 / 90) < 1e-6
    for d in (2, 3, 4, 10):
        t11, t1122, t12 = SH.moments_exact(big, d)
        assert abs(float(t11) * d * d / big**2 - 1 / 20) < 1e-5 and abs(float(t1122) * d * d / big**2 - 1 / 36) < 1e-5
        assert abs(float(t12) * d * d / big**2 - 1 / 90) < 1e-5
        assert SH.A_from_moments(d, F(1, 20), F(1, 36), F(1, 90)) == SH.A_d(d) == F(4 + 2 * d, 4 + 5 * d)
    assert F(1, 36) + 2 * F(1, 90) == F(1, 20) and F(1, 36) + F(1, 90) == F(7, 180)      # 筆誤框 2
    assert F(1, 20) / F(1, 36) == F(9, 5) and F(1, 20) / F(1, 90) == F(9, 2)
    assert f4(nb["A4"]) == "0.6013" and f4(nb["A100"]) == "0.5265"
    assert [f"{v:.2f}" for v in nb["sim"]] == ["56.09", "31.31", "12.56"]
    assert [f"{v:.2f}" for v in nb["mex"]] == ["56.67", "31.48", "12.60"]
    assert f4(nb["Asim"]) == "0.5261" and nb["Ainf"] == F(10, 19) and f4(float(nb["Ainf"])) == "0.5263"


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def place_in(m, x0, x1, y0, y1):
    """縮到放得進方框 [x0,x1]×[y0,y1]（只縮不放），再置中在方框裡。"""
    if m.get_width() > x1 - x0:
        m.set_width(x1 - x0)
    if m.get_height() > y1 - y0:
        m.set_height(y1 - y0)
    return m.move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))


def eq_label(tex_mob, label, buff=0.4):
    return VGroup(tex_mob, tx(r"\text{" + label + "}", 22, GREY_B)).arrange(RIGHT, buff=buff)


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


def small_walk(unit, origin):
    """東東北西：點、步伐箭頭與 η 標籤。"""
    pts = np.asarray(SH.sites(STEPS))
    org = np.array(origin)
    P = [org + unit * (x * RIGHT + y * UP) for x, y in pts]
    dots = VGroup(*[Dot(p, radius=0.08).set_color(WHITE) for p in P])
    arrows = VGroup(*[Arrow(P[i], P[i + 1], buff=0.1, thickness=4).set_color(c) for i, c in enumerate(STEP_COLORS)])
    offs = [DOWN, DOWN, RIGHT, UP]
    labs = VGroup(*[tx(r"\eta_%d" % (i + 1), 28, STEP_COLORS[i]).next_to(arrows[i], offs[i], buff=0.12)
                    for i in range(4)])
    return VGroup(dots, arrows, labs)


def mode_chart(k, origin):
    """第 k 個正弦模態：4 根長條＋虛線 sin(kπx/5)（同一單位化）。"""
    nb = numbers()
    v = nb["modes"][k - 1]
    norm = float(np.linalg.norm(np.sin(k * np.pi * np.arange(1, 5) / 5)))
    ax = D.plot_axes((0, 5, 1), (-0.8, 0.8, 0.4), width=2.7, height=1.7, x_ticks=[1, 2, 3, 4], font_size=18)
    ax.move_to(origin)
    col = MODE_COLORS[k - 1]
    bars = VGroup()
    for al in range(1, 5):
        p0, p1 = ax.c2p(al - 0.28, 0), ax.c2p(al + 0.28, v[al - 1])
        r = Rectangle(width=abs(p1[0] - p0[0]), height=max(abs(p1[1] - p0[1]), 1e-3))
        r.set_fill(col, 0.7).set_stroke(col, 1.5).move_to((p0 + p1) / 2)
        bars.add(r)
    g = ax.get_graph(lambda x: np.sin(k * np.pi * x / 5) / norm, x_range=(0, 5, 0.05))
    sine = DashedVMobject(g, num_dashes=24).set_stroke(WHITE, 2)
    title = tx(r"k=%d" % k, 26, col).next_to(ax, UP, buff=0.08)
    return VGroup(ax, bars, sine, title)


# ---------------------------------------------------------------- 各段
def part_question(s):
    walk = small_walk(1.3, np.array([-5.6, -0.9, 0]))
    q1 = eq_label(tx(r"\mathbf{T}=\frac{1}{N+1}\sum_{j=0}^{N}(\mathbf r_j-\langle\mathbf r\rangle)"
                     r"(\mathbf r_j-\langle\mathbf r\rangle)^{\mathsf T}", 36), "(6.1)--(6.3)")
    q2 = tx(r"\mathbf r_j=\eta_1+\eta_2+\dots+\eta_j", 34)
    q3 = mixed("$\\mathbf{T}=\\sum_{\\alpha,\\beta}\\ ?\\ \\eta_\\alpha\\,\\eta_\\beta$", "　係數是什麼？",
               font_size=30, color=YELLOW)
    col = VGroup(q1, q2, q3).arrange(DOWN, buff=0.65)
    place_in(walk, -6.7, -1.6, -2.8, 2.8)
    place_in(col, -1.2, 6.7, -2.8, 3.1)
    assert_in_frame(walk, col)
    assert_apart(walk, col)
    assert_column(*col)
    s.play(FadeIn(walk[0]), run_time=0.6)
    s.play(LaggedStart(*[GrowArrow(a) for a in walk[1]], lag_ratio=0.4), FadeIn(walk[2]), run_time=1.6)
    s.say("W602 的小例子：東、東、北、西，", "$N=4$", " 步、5 個點", wait=1.5)
    s.play(FadeIn(q1), run_time=0.8)
    s.say("問：", "$\\mathbf{T}$", " 是 ", "$N+1$", " 個點的平均（減掉重心）", wait=2.0)
    s.play(FadeIn(q2), run_time=0.6)
    s.say("每個點又是步伐的和", wait=1.5)
    s.play(FadeIn(q3), run_time=0.6)
    s.say("能不能直接寫成「步伐乘步伐」？係數是什麼？（先猜）", wait=3.0)
    s.clear_stage()


def part_kramers(s):
    nb = numbers()
    a = nb["a25"]
    k1 = eq_label(tx(r"T_{ij}=\sum_{\alpha,\beta=1}^{N}a_{\alpha\beta}\,\eta_{\alpha,i}\,\eta_{\beta,j}", 34), "(6.21)")
    k2 = eq_label(tx(r"a_{\alpha\beta}=\frac{1}{(N+1)^2}\begin{cases}\alpha\,(N+1-\beta), & \alpha\le\beta\\"
                     r"\beta\,(N+1-\alpha), & \alpha>\beta\end{cases}", 32), "(6.22)")
    top = VGroup(k1, k2).arrange(RIGHT, buff=0.8)
    place_in(top, -6.75, 6.75, 1.75, 3.15)
    rows = [[r"\alpha\backslash\beta", "1", "2", "3", "4"]] + [[str(i + 1)] + [str(v) for v in a[i]] for i in range(4)]
    tab = D.tex_table(rows, [1.0, 0.62, 0.62, 0.62, 0.62], row_h=0.55, font_size=30)
    for i in range(1, 5):
        D.table_cell(tab, i, 0, 5).set_color(YELLOW)
    pre = tx(r"\mathbf a=\frac{1}{25}\times", 32)
    left = VGroup(pre, tab).arrange(RIGHT, buff=0.3)
    place_in(left, -6.75, -0.2, -1.25, 1.6)
    hl = VGroup(SurroundingRectangle(D.table_cell(tab, 2, 3, 5), buff=0.08),
                SurroundingRectangle(D.table_cell(tab, 3, 2, 5), buff=0.08)).set_stroke(ORANGE, 3)
    p1 = tx(r"\mathbf{T}=\frac{1}{(N+1)^2}\sum_{j<k}(\mathbf r_k-\mathbf r_j)(\mathbf r_k-\mathbf r_j)^{\mathsf T}", 30)
    p2 = tx(r"\mathbf r_k-\mathbf r_j=\eta_{j+1}+\dots+\eta_k", 30)
    p3 = mixed("含 ", "$\\eta_2\\eta_3$", " 的點對：", "$j\\in\\{0,1\\},\\ k\\in\\{3,4\\}$", font_size=26)
    p4 = mixed("個數 ", "$m\\,(N+1-n)=2\\times2=4$", font_size=26, color=ORANGE)
    right = VGroup(p1, p2, p3, p4).arrange(DOWN, buff=0.4, aligned_edge=LEFT)
    place_in(right, 0.2, 6.75, -2.9, 1.6)
    assert_in_frame(top, left, right)
    assert_apart(top, left)
    assert_apart(top, right)
    assert_apart(left, right)
    assert_column(*right)
    assert_table(tab, 5)
    s.play(FadeIn(k1), run_time=0.8)
    s.say("Kramers 恆等式：", "$\\mathbf{T}$", " 是步伐的二次式", wait=2.0)
    s.play(FadeIn(k2), run_time=0.8)
    s.say("係數 ", "$a_{\\alpha\\beta}$", " 只看兩步的編號，和怎麼走無關", wait=2.5)
    s.play(FadeIn(pre), FadeIn(tab[0][:5]), ShowCreation(tab[1]), run_time=0.8)
    s.play(LaggedStart(*[FadeIn(tab[0][5 * i:5 * i + 5]) for i in range(1, 5)], lag_ratio=0.3), run_time=1.5)
    s.say("$N=4$", "：乘 ", "$\\frac{1}{25}$", "，對角線最大、離得越遠越小", wait=2.5)
    s.play(FadeIn(p1), FadeIn(p2), run_time=0.8)
    s.say("為什麼？", "$\\mathbf{T}$", " 也是所有點對差的平方平均", wait=2.5)
    s.play(ShowCreation(hl), FadeIn(p3), FadeIn(p4), run_time=0.8)
    s.say("表裡的整數＝含 ", "$\\eta_m\\eta_n$", " 的點對數 ", "$m(N+1-n)$", wait=3.0)
    # 代入東東北西
    c1 = mixed("東東北西：", "$\\eta_{\\alpha,1}=(1,1,0,-1)$", "，", "$\\eta_{\\alpha,2}=(0,0,1,0)$", font_size=26)
    c2 = tx(r"T_{11}=\frac{4+6+4+2\cdot3-2\cdot1-2\cdot2}{25}=\frac{14}{25}", 30)
    c3 = tx(r"T_{12}=\frac{2+4-3}{25}=\frac{3}{25},\qquad T_{22}=a_{33}=\frac{6}{25}", 30)
    c4 = mixed("$\\mathbf{T}=\\frac{1}{25}\\begin{pmatrix}14&3\\\\3&6\\end{pmatrix}$", "＝W602 的定義式",
               font_size=28, color=GREEN_B)
    sub = VGroup(c1, c2, c3, c4).arrange(DOWN, buff=0.32, aligned_edge=LEFT)
    place_in(sub, 0.2, 6.75, -2.9, 1.6)
    box = D.erratum_box("p. 139，(6.22)", r"\alpha<\beta,\ \ \alpha>\beta", r"\alpha\le\beta,\ \ \alpha>\beta", width=5.6)
    place_in(box, -6.75, -0.2, -2.95, -1.4)
    assert_in_frame(sub, box)
    assert_apart(sub, left)
    assert_apart(sub, top)
    assert_apart(sub, box)
    assert_apart(left, box)
    assert_column(*sub)
    s.play(FadeOut(right), FadeOut(hl), run_time=0.5)
    s.play(FadeIn(c1), run_time=0.6)
    s.say("代入「東東北西」：第 1、2 個分量的步伐", wait=2.0)
    s.play(FadeIn(c2), FadeIn(c3), run_time=0.8)
    s.say("只用表裡的整數，三個分量就算出來了", wait=2.0)
    s.play(FadeIn(c4), run_time=0.6)
    s.say("和 W602 用定義式算的一模一樣", wait=2.0)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上只寫 ", "$\\alpha<\\beta$", "、", "$\\alpha>\\beta$", "：少了對角線", wait=2.0)
    s.say("少了對角線，這條漫步會算出 ", "$T_{11}=T_{22}=0$", wait=3.0)
    s.clear_stage()


def part_sine(s):
    nb = numbers()
    charts = VGroup(*[mode_chart(k, np.array([-5.1 + 3.4 * (k - 1), 1.6, 0])) for k in range(1, 5)])
    rows = [["k", "1", "2", "3", "4"],
            [r"\lambda_k"] + [f4(v) for v in nb["eig"]],
            [r"\frac{N+1}{\pi^2k^2}"] + [f4(v) for v in nb["asym"]]]
    tab = D.tex_table(rows, [2.0, 1.6, 1.6, 1.6, 1.6], row_h=0.66, font_size=28)
    for j in range(1, 5):
        D.table_cell(tab, 0, j, 5).set_color(MODE_COLORS[j - 1])
    place_in(charts, -6.75, 6.75, 0.3, 3.15)
    place_in(tab, -6.75, 6.75, -2.15, 0.1)
    tr = mixed("$\\mathrm{Tr}\\,\\mathbf a=\\lambda_1+\\lambda_2+\\lambda_3+\\lambda_4=\\frac45$", "＝平均的 ", "$R_g^2$",
               font_size=26, color=GREEN_B)
    place_in(tr, -6.75, 6.75, -2.95, -2.3)
    assert_in_frame(charts, tab, tr)
    for c1, c2 in zip(charts, charts[1:]):
        assert_apart(c1, c2)
    assert_apart(charts, tab)
    assert_apart(tab, tr)
    assert_table(tab, 5)
    s.say("$\\mathbf a$", " 是對稱矩陣：它的特徵向量長什麼樣子？", wait=1.0)
    for k in range(4):
        ax, bars, sine, title = charts[k]
        s.play(FadeIn(ax), FadeIn(title), LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.15),
               run_time=1.0)
    s.say("長條＝特徵向量的 4 個分量", wait=1.5)
    s.play(*[ShowCreation(c[2]) for c in charts], run_time=1.2)
    s.say("虛線：", "$\\sin\\frac{k\\pi\\alpha}{N+1}$", "，精確地穿過每根長條", wait=2.5)
    s.say("像一條兩端固定的弦：", "$k$", " 個半波", wait=2.0)
    s.play(FadeIn(tab[0][:5]), ShowCreation(tab[1]), FadeIn(tab[0][5:10]), run_time=0.8)
    s.say("特徵值 0.5236、0.1447、0.0764、0.0553", wait=2.0)
    s.play(FadeIn(tab[0][10:]), run_time=0.6)
    s.say("大 ", "$N$", " 的公式 ", "$\\frac{N+1}{\\pi^2k^2}$", "：", "$k$", " 小時接近", wait=2.5)
    s.play(FadeIn(tr), run_time=0.6)
    s.say("特徵值加起來 ", "$\\frac45$", "：步伐長 1 時的平均 ", "$R_g^2$", wait=2.5)
    s.clear_stage()


def part_gauss(s):
    g1 = eq_label(tx(r"\langle\eta_{\alpha,i}\,\eta_{\beta,j}\rangle=\frac1d\,\delta_{\alpha\beta}\,\delta_{ij}", 34),
                  "(6.28)--(6.30)")
    g2 = eq_label(tx(r"\langle\eta_\alpha\eta_\beta\eta_\gamma\eta_\delta\rangle="
                     r"\langle\eta_\alpha\eta_\beta\rangle\langle\eta_\gamma\eta_\delta\rangle"
                     r"+\langle\eta_\alpha\eta_\gamma\rangle\langle\eta_\beta\eta_\delta\rangle"
                     r"+\langle\eta_\alpha\eta_\delta\rangle\langle\eta_\beta\eta_\gamma\rangle", 30), "(6.33)")
    g3 = eq_label(tx(r"\langle T_{11}^2\rangle=\frac{(\mathrm{Tr}\,\mathbf a)^2+2\,\mathrm{Tr}\,\mathbf a^2}{d^2},\quad"
                     r"\langle T_{11}T_{22}\rangle=\frac{(\mathrm{Tr}\,\mathbf a)^2}{d^2},\quad"
                     r"\langle T_{12}^2\rangle=\frac{\mathrm{Tr}\,\mathbf a^2}{d^2}", 32), "(6.36)(6.37)")
    col = VGroup(g1, g2, g3).arrange(DOWN, buff=0.7)
    place_in(col, -6.6, 6.6, -2.8, 3.1)
    assert_in_frame(col)
    assert_column(*col)
    s.play(FadeIn(g1), run_time=0.8)
    s.say("高斯步伐：每個分量獨立、常態、變異數 ", "$1/d$", wait=2.5)
    s.play(FadeIn(g2), run_time=0.8)
    s.say("四個高斯相乘的平均＝三種「兩兩配對」的和", wait=2.5)
    s.play(FadeIn(g3), run_time=0.8)
    s.say("代入 Kramers 恆等式：只剩 ", "$\\mathrm{Tr}\\,\\mathbf a$", " 與 ", "$\\mathrm{Tr}\\,\\mathbf a^2$", wait=3.0)
    g4 = tx(r"\mathrm{Tr}\,\mathbf a\approx\frac N6,\ \ \mathrm{Tr}\,\mathbf a^2\approx\frac{N^2}{90}"
            r"\ \ \Rightarrow\ \ \frac{N^2}{20d^2},\ \ \frac{N^2}{36d^2},\ \ \frac{N^2}{90d^2}", 34)
    g5 = tx(r"\frac{\langle T_{11}^2\rangle}{\langle T_{11}T_{22}\rangle}=\frac95,\quad"
            r"\frac{\langle T_{11}^2\rangle}{\langle T_{12}^2\rangle}=\frac92"
            r"\quad\overset{(6.14)}{\Longrightarrow}\quad A_d=\frac{4+2d}{4+5d}", 36)
    g5.set_color(YELLOW)
    box = D.erratum_box("p. 141，(6.35b)", r"\frac{\pi^4}{36}+\frac{\pi^4}{90}", r"\frac{\pi^4}{36}+\frac{2\pi^4}{90}",
                        width=6.0)
    col2 = VGroup(g4, box, g5).arrange(DOWN, buff=0.45)
    fit(col2, 13.0)
    gy = 3.15 - g3.get_height() / 2                                   # g3 移到最上面
    room = gy - g3.get_height() / 2 - 0.3 - (-2.95)
    if col2.get_height() > room:
        col2.set_height(room)
    at(col2, 0, gy - g3.get_height() / 2 - 0.3 - col2.get_height() / 2)
    s.play(FadeOut(g1), FadeOut(g2), g3.animate.move_to(np.array([0, gy, 0])), run_time=0.8)
    assert_in_frame(g3, col2)
    assert_apart(g3, col2)
    assert_column(*col2)
    s.play(FadeIn(g4), run_time=0.8)
    s.say("大 ", "$N$", "：", "$\\mathrm{Tr}\\,\\mathbf a\\approx N/6$", "、", "$\\mathrm{Tr}\\,\\mathbf a^2\\approx N^2/90$",
          wait=2.5)
    s.play(FadeIn(box), run_time=0.8)
    s.say("書上 (6.35b) 漏了配對的 2 倍：", "$\\frac{\\pi^4}{90}$", " 應為 ", "$\\frac{2\\pi^4}{90}$", wait=2.5)
    s.say("照書上得 ", "$\\frac{7N^2}{180d^2}$", "，不是 ", "$\\frac{N^2}{20d^2}$", wait=2.5)
    s.play(FadeIn(g5), run_time=0.8)
    s.say("比值 ", "$\\frac95$", "、", "$\\frac92$", " 代入 (6.14)：正是 ", "$\\frac{4+2d}{4+5d}$", wait=3.0)
    s.clear_stage()


def part_table(s):
    nb = numbers()
    hz = lambda t, c=YELLOW: zh(t, font_size=24, color=c)
    rows = [[hz("d = 3，N = 100"), hz("高斯模擬（2 萬條）"), hz("精確")],
            [r"\langle T_{11}^2\rangle", f"{nb['sim'][0]:.2f}", f"{nb['mex'][0]:.2f}"],
            [r"\langle T_{11}T_{22}\rangle", f"{nb['sim'][1]:.2f}", f"{nb['mex'][1]:.2f}"],
            [r"\langle T_{12}^2\rangle", f"{nb['sim'][2]:.2f}", f"{nb['mex'][2]:.2f}"],
            ["A", f4(nb["Asim"]), f4(nb["A100"])]]
    tab = D.tex_table(rows, [2.6, 2.8, 1.6], row_h=0.72, font_size=30)
    place_in(tab, -6.75, -0.3, -2.9, 3.1)
    hd = zh("精確的平均非球度（d = 3）", font_size=26, color=YELLOW)
    r1 = mixed("$N=4$", "：", "$%s$" % f4(nb["A4"]), font_size=28)
    r2 = mixed("$N=100$", "：", "$%s$" % f4(nb["A100"]), font_size=28)
    r3 = mixed("$N\\to\\infty$", "：", "$\\frac{10}{19}=%s$" % f4(float(nb["Ainf"])), font_size=28, color=GREEN_B)
    right = VGroup(hd, r1, r2, r3).arrange(DOWN, buff=0.45, aligned_edge=LEFT)
    place_in(right, 0.3, 6.75, -2.9, 3.1)
    assert_in_frame(tab, right)
    assert_apart(tab, right)
    assert_table(tab, 3)
    assert_column(*right)
    s.play(FadeIn(tab[0][:3]), ShowCreation(tab[1]), run_time=0.8)
    s.say("核對：三維高斯步伐 2 萬條、每條 100 步", wait=1.0)
    s.play(LaggedStart(*[FadeIn(tab[0][3 * i:3 * i + 3]) for i in range(1, 5)], lag_ratio=0.35), run_time=1.6)
    s.say("模擬 56.09、31.31、12.56；精確 56.67、31.48、12.60", wait=3.0)
    s.play(FadeIn(right[:3]), run_time=0.8)
    s.say("$N=4$", " 就有精確值 0.6013；", "$N=100$", "：0.5265", wait=2.5)
    s.play(FadeIn(r3), run_time=0.6)
    s.say("$N\\to\\infty$", "：", "$\\frac{10}{19}=0.5263$", "（W603 的答案）", wait=2.5)
    s.say("Kramers 恆等式把形狀變成「步伐的二次式」", wait=2.5)
    s.say("高斯步伐一配對，", "$\\frac{10}{19}$", " 就算出來了", wait=2.5)


class W605KramersIdentity(ChuScene):
    def construct(self):
        self.show_header("W605　Kramers 恆等式：形狀是步伐的二次式")
        check_math()
        part_question(self)
        part_kramers(self)
        part_sine(self)
        part_gauss(self)
        part_table(self)
        self.wait(1.0)
