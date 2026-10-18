"""
W602 迴轉半徑張量：形狀的指紋（§6.1.2–6.1.4，pp. 129–133，(6.1)–(6.9)，Fig. 6.4；補充 1）

問：怎麼用幾個數字說出「這條漫步往哪個方向長、長多少」？（提示：轉動慣量）
  1. (6.1)–(6.3) 的定義；小例子 T = (1/25)[[14,3],[3,6]]，特徵值 3/5、1/5，主軸 (3,1)（18.43°）；不變量 (6.6)–(6.9)。
  2. 補充 1：沿方向 n 的平方寬度 n·T·n（0–180°）：最大 3/5 在 18.43°、最小 1/5 在 108.43°。
  3. Fig. 6.4：二維 1000 步，重心 (−17.9, −8.9)，T11 207.5、T22 82.9、T12 95.0，R1² 258.8、R2² 31.6，主軸 28.4°。
  4. 繞重心轉 30°：T 變成 94.1、196.3、101.4，特徵值不變，主軸 58.4°；對方向平均 ⟨T⟩ = T̄ I（(6.5)）。
畫面上的數字都由 utils.walk_shape 算出，check_math() 核對。
"""
import math
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

THETA = math.pi / 6
COM_COL = ORANGE
AX1_COL = "#7CE38B"
AX2_COL = PINK
CURVE_COL = TEAL_B


# ---------------------------------------------------------------- 數學
def f1(x):
    return f"{float(x):.1f}"


def f2(x):
    return f"{float(x):.2f}"


def frac_tex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return ("-" if q < 0 else "") + r"\frac{%d}{%d}" % (abs(q.numerator), q.denominator)


def deg(v):
    return math.degrees(math.atan2(v[1], v[0]))


@lru_cache(maxsize=None)
def numbers():
    se = SH.small_example()
    Tf = se["T"]
    Tn = np.array(Tf, dtype=float)
    tr, tr2, _, det = SH.invariants(Tn)
    phi1 = math.atan2(se["vec"][1], se["vec"][0])
    fs = SH.fig64_stats(THETA)
    T, Tr = fs["T"], fs["Trot"]
    vals_r, vecs_r = SH.principal(Tr)
    vecs = SH.principal(T)[1]
    p = SH.fig64()
    ang = np.linspace(0, 2 * math.pi, 720, endpoint=False)
    Tavg = np.mean([SH.tensor(SH.rotate(p, a)) for a in ang[::8]], axis=0)
    return dict(se=se, Tn=Tn, inv=(tr, tr2, det), phi1=phi1, phi_deg=f2(math.degrees(phi1)),
                phi2_deg=f2(math.degrees(phi1) + 90),
                ext=(float(SH.extent(Tn, phi1)), float(SH.extent(Tn, phi1 + math.pi / 2))),
                pts=p, com=fs["com"], vec1=vecs[:, 0], vec2=vecs[:, 1], Rraw=fs["R"],
                com_s=(f1(fs["com"][0]), f1(fs["com"][1])),
                T=(f1(T[0, 0]), f1(T[1, 1]), f1(T[0, 1])), eig=(f1(fs["eig"][0]), f1(fs["eig"][1])),
                R=(f2(fs["R"][0]), f2(fs["R"][1])), ratio=f2(fs["ratio"]), angle=f1(fs["angle"]),
                Trot=(f1(Tr[0, 0]), f1(Tr[1, 1]), f1(Tr[0, 1])), eig_rot=(f1(vals_r[0]), f1(vals_r[1])),
                angle_rot=f1(deg(vecs_r[:, 0])), Tbar=f1(np.trace(T) / 2), Tavg=Tavg)


def check_math():
    nb = numbers()
    se = nb["se"]
    F = Fraction
    assert se["T"] == [[F(14, 25), F(3, 25)], [F(3, 25), F(6, 25)]]
    assert SH.tensor_frac(se["pts"]) == se["T"]
    assert se["eig"] == (F(3, 5), F(1, 5)) and se["vec"] == (3, 1)
    T = se["T"]
    # (3,1) 是特徵向量：T·(3,1) = 3/5·(3,1)
    assert T[0][0] * 3 + T[0][1] * 1 == F(3, 5) * 3 and T[1][0] * 3 + T[1][1] * 1 == F(3, 5) * 1
    trf = T[0][0] + T[1][1]
    tr2f = T[0][0] ** 2 + 2 * T[0][1] ** 2 + T[1][1] ** 2
    detf = T[0][0] * T[1][1] - T[0][1] ** 2
    assert (trf, tr2f, detf) == (F(4, 5), F(2, 5), F(3, 25)) and detf == F(3, 5) * F(1, 5)
    tr, tr2, det = nb["inv"]
    assert abs(tr - 0.8) < 1e-12 and abs(tr2 - 0.4) < 1e-12 and abs(det - 0.12) < 1e-12
    assert nb["phi_deg"] == "18.43" and nb["phi2_deg"] == "108.43"
    assert abs(nb["ext"][0] - 0.6) < 1e-12 and abs(nb["ext"][1] - 0.2) < 1e-12
    phis = np.linspace(0, math.pi, 18001)
    w = SH.extent(nb["Tn"], phis)
    assert f2(math.degrees(phis[np.argmax(w)])) == "18.43" and f2(math.degrees(phis[np.argmin(w)])) == "108.43"
    assert abs(w.max() - 0.6) < 1e-6 and abs(w.min() - 0.2) < 1e-6
    p = nb["pts"]
    assert len(p) == 1001 and tuple(p.min(axis=0)) == (-49, -24) and tuple(p.max(axis=0)) == (4, 7)
    assert nb["com_s"] == ("-17.9", "-8.9")
    assert nb["T"] == ("207.5", "82.9", "95.0")
    assert nb["eig"] == ("258.8", "31.6") and nb["R"] == ("16.09", "5.62") and nb["ratio"] == "8.18"
    assert nb["angle"] == "28.4"
    assert nb["Trot"] == ("94.1", "196.3", "101.4") and nb["eig_rot"] == ("258.8", "31.6")
    assert nb["angle_rot"] == "58.4"
    assert nb["Tbar"] == "145.2"
    Ta = nb["Tavg"]                                                          # (6.5)：對方向平均 → T̄ I
    assert abs(Ta[0, 1]) < 1e-9 and abs(Ta[0, 0] - Ta[1, 1]) < 1e-9 and f1(Ta[0, 0]) == "145.2"


# ---------------------------------------------------------------- 小工具
def tx(tex, font_size=30, color=WHITE):
    return Tex(tex, font_size=font_size).set_color(color)


def fit(m, width):
    if m.get_width() > width:
        m.set_width(width)
    return m


def at(m, x, y):
    return m.move_to(np.array([x, y, 0]))


def eq_label(tex_mob, label, buff=0.35):
    return VGroup(tex_mob, tx(r"\text{" + label + "}", 22, GREY_B)).arrange(RIGHT, buff=buff)


def to_screen(P, unit, c, origin):
    P = np.asarray(P, dtype=float)
    out = np.zeros((len(P), 3))
    out[:, 0] = origin[0] + unit * (P[:, 0] - c[0])
    out[:, 1] = origin[1] + unit * (P[:, 1] - c[1])
    return out


def cross_mark(p, size=0.15, color=COM_COL, width=4):
    a = Line(p + size * UL, p + size * DR).set_stroke(color, width)
    b = Line(p + size * UR, p + size * DL).set_stroke(color, width)
    return VGroup(a, b)


def axis_line(center, vec, half, color, width=4):
    v = np.array([vec[0], vec[1], 0.0])
    v /= np.linalg.norm(v)
    return Line(center - half * v, center + half * v).set_stroke(color, width)


def principal_axes(center, R1, R2, vec1, unit):
    """兩條主軸（半長 2R_i）＋淡色橢圓（半軸 2R_i）。"""
    a1 = axis_line(center, vec1, 2 * R1 * unit, AX1_COL)
    a2 = axis_line(center, (-vec1[1], vec1[0]), 2 * R2 * unit, AX2_COL)
    e = Ellipse(width=4 * R1 * unit, height=4 * R2 * unit).rotate(math.atan2(vec1[1], vec1[0])).move_to(center)
    e.set_stroke(GREY_B, 1.5).set_fill(GREY_B, 0.10)
    return VGroup(e, a1, a2)


def small_pic(unit=1.2, com_screen=(-4.1, -0.45)):
    se = numbers()["se"]
    com = np.array([float(v) for v in se["com"]])
    cs = np.array([com_screen[0], com_screen[1], 0.0])
    scr = to_screen(se["pts"], unit, com, cs)
    bg = VGroup(*[Dot(to_screen([(x, y)], unit, com, cs)[0], radius=0.03).set_fill(GREY_D, 1)
                  for x in range(-1, 4) for y in range(-1, 3)])
    path = VGroup(*[Arrow(scr[k], scr[k + 1], buff=0.08, thickness=4).set_color(YELLOW) for k in range(4)])
    dots = VGroup(*[Dot(p, radius=0.07).set_fill(WHITE, 1) for p in scr])
    l1, l2 = (math.sqrt(float(v)) for v in se["eig"])
    axes = principal_axes(cs, l1, l2, np.array(se["vec"], dtype=float), unit)
    return dict(bg=bg, path=path, dots=dots, cross=cross_mark(cs), axes=axes, com=cs, unit=unit)


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    p = nb["pts"]
    unit = 0.10
    c = (p.min(axis=0) + p.max(axis=0)) / 2
    walk = VMobject().set_points_as_corners(to_screen(p, unit, c, (-3.5, 0.2))).set_stroke(YELLOW, 1.8)
    q1 = at(zh("這條漫步往哪個方向長？", font_size=32), 3.5, 1.7)
    q2 = at(zh("長多少？", font_size=32), 3.5, 1.0)
    q3 = at(zh("能不能用幾個數字說出來？", font_size=30, color=YELLOW), 3.5, -0.2)
    q4 = at(zh("提示：轉動慣量", font_size=26, color=GREY_B), 3.5, -1.2)
    assert_in_frame(walk, q1, q2, q3, q4)
    for a, b in ((walk, q1), (walk, q3), (walk, q4), (q1, q2), (q2, q3), (q3, q4)):
        assert_apart(a, b)
    s.play(ShowCreation(walk), run_time=3.0)
    s.say("一條二維 1000 步的漫步", wait=0.3)
    s.play(FadeIn(q1), FadeIn(q2), run_time=0.6)
    s.say("問：怎麼用幾個數字說出它往哪個方向長、長多少？", wait=2.8)
    s.play(FadeIn(q3), FadeIn(q4), run_time=0.6)
    s.say("提示：轉動慣量（W601 的 ", "$R_g^2$", " 是它的跡）", wait=2.8)
    s.clear_stage()


def part_tensor(s):
    nb = numbers()
    se = nb["se"]
    sp = small_pic()
    pic = VGroup(sp["bg"], sp["path"], sp["dots"], sp["cross"], sp["axes"])
    e1 = eq_label(tx(r"T_{kl}=\frac{1}{N+1}\sum_{j=0}^{N}\big(r_{jk}-\langle r_k\rangle\big)"
                     r"\big(r_{jl}-\langle r_l\rangle\big)", 32), "(6.1)--(6.3)")
    at(e1, 0, 2.5)
    e2 = at(tx(r"\mathbf T=\begin{pmatrix}T_{11}&T_{12}\\T_{12}&T_{22}\end{pmatrix}"
               r"=\frac{1}{25}\begin{pmatrix}14&3\\3&6\end{pmatrix}", 32), 2.9, 1.15)
    l1, l2 = se["eig"]
    e3 = at(mixed("$R_1^2=%s,\\quad R_2^2=%s$" % (frac_tex(l1), frac_tex(l2)), "；主軸 ",
                  "$(3,1)$", "：", "$%s^\\circ$" % nb["phi_deg"], font_size=30), 2.9, -0.1)
    e3[0].set_color(AX1_COL)
    e4 = eq_label(tx(r"\mathrm{Tr}\,\mathbf T=\frac{4}{5},\quad \mathrm{Tr}\,\mathbf T^2=\frac{2}{5},\quad "
                     r"\det\mathbf T=\frac{3}{25}=\frac{3}{5}\cdot\frac{1}{5}", 30), "(6.6)--(6.9)")
    at(e4, 2.9, -1.3)
    e5 = at(mixed("$\\mathrm{Tr}\\,\\mathbf T=R_1^2+R_2^2=R_g^2$", "（W601）", font_size=28, color=GREY_A), 2.9, -2.3)
    for m in (e2, e3, e4, e5):
        fit(m, 6.9)
    assert_in_frame(pic, e1, e2, e3, e4, e5)
    for a, b in ((e1, pic), (e1, e2), (e2, e3), (e3, e4), (e4, e5), (pic, e2), (pic, e3), (pic, e4), (pic, e5)):
        assert_apart(a, b)
    s.play(Write(e1), run_time=1.5)
    s.say("定義：每個點離重心的位移，分量兩兩相乘再平均", wait=2.8)
    s.say("除以 ", "$N+1$", "：", "$N$", " 步有 ", "$N+1$", " 個點", wait=2.4)
    s.play(FadeIn(sp["bg"]), FadeIn(sp["dots"]), LaggedStart(*[GrowArrow(a) for a in sp["path"]], lag_ratio=0.3),
           run_time=1.2)
    s.play(ShowCreation(sp["cross"]), run_time=0.5)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("二維是 ", "$2\\times2$", " 對稱矩陣；W601 的小例子手算", wait=2.9)
    s.play(FadeIn(e3), run_time=0.6)
    s.play(FadeIn(sp["axes"][0]), ShowCreation(sp["axes"][1]), ShowCreation(sp["axes"][2]), run_time=1.0)
    s.say("特徵值 ", "$R_1^2=\\frac35$", "、", "$R_2^2=\\frac15$", "：主迴轉半徑的平方", wait=2.9)
    s.say("主軸沿 ", "$(3,1)$", "，", "$%s^\\circ$" % nb["phi_deg"], "（綠線、粉線半長 ", "$2R_i$", "）", wait=2.9)
    s.play(FadeIn(e4), run_time=0.8)
    s.say("不變量：跡、平方的跡、行列式，轉了都不變", wait=2.9)
    s.play(FadeIn(e5), run_time=0.6)
    s.say("$\\mathrm{Tr}\\,\\mathbf T=\\frac45$", "：正是 W601 的 ", "$R_g^2$", wait=2.6)
    s.play(*[FadeOut(m) for m in (e1, e2, e3, e4, e5)], run_time=0.6)
    return sp


def part_extent(s, sp):
    nb = numbers()
    Tn = nb["Tn"]
    ax = D.plot_axes((0, 180, 30), (0, 0.7, 0.1), width=6.2, height=3.3, x_ticks=(0, 45, 90, 135, 180),
                     y_ticks=(0.2, 0.4, 0.6), x_tick_labels=[r"0^\circ", r"45^\circ", r"90^\circ", r"135^\circ",
                                                             r"180^\circ"])
    ax.move_to(np.array([2.9, -0.3, 0]))
    ax.add(tx(r"\varphi", 26, GREY_B).next_to(ax.x_axis.get_end(), RIGHT, buff=0.12))
    graph = ax.get_graph(lambda x: float(SH.extent(Tn, math.radians(x))), x_range=(0, 180)).set_stroke(CURVE_COL, 3)
    head = at(tx(r"\vec n\cdot\mathbf T\cdot\vec n,\qquad \vec n=(\cos\varphi,\sin\varphi)", 32), 3.0, 2.5)
    p1, p2 = math.degrees(nb["phi1"]), math.degrees(nb["phi1"]) + 90
    d1 = Dot(ax.c2p(p1, nb["ext"][0]), radius=0.08).set_fill(AX1_COL, 1)
    d2 = Dot(ax.c2p(p2, nb["ext"][1]), radius=0.08).set_fill(AX2_COL, 1)
    m1 = tx(r"\max=\frac{3}{5}\ \ (%s^\circ)" % nb["phi_deg"], 26, AX1_COL).next_to(d1, UP, buff=0.12)
    m1.shift((ax.c2p(0, 0)[0] + 0.25 - m1.get_left()[0]) * RIGHT * (m1.get_left()[0] < ax.c2p(0, 0)[0] + 0.25))
    m2 = tx(r"\min=\frac{1}{5}\ \ (%s^\circ)" % nb["phi2_deg"], 24, AX2_COL).next_to(d2, DOWN, buff=0.12)
    pic = VGroup(sp["bg"], sp["path"], sp["dots"], sp["cross"], sp["axes"])
    assert_in_frame(ax, graph, head, m1, m2, pic)
    for a, b in ((head, ax), (pic, ax), (pic, head), (m1, m2), (m1, graph), (m1, head), (m2, d2)):
        assert_apart(a, b)
    assert m2.get_bottom()[1] > ax.c2p(0, 0)[1] + 0.05                    # 在曲線最低點與橫軸之間
    phi = ValueTracker(0.0)
    cs = sp["com"]

    def dir_line():
        a = math.radians(phi.get_value())
        v = np.array([math.cos(a), math.sin(a), 0])
        return DashedLine(cs - 1.9 * v, cs + 1.9 * v, dash_length=0.1).set_stroke(CURVE_COL, 3)

    def cur_dot():
        x = phi.get_value()
        return Dot(ax.c2p(x, float(SH.extent(Tn, math.radians(x)))), radius=0.07).set_fill(CURVE_COL, 1)

    line = always_redraw(dir_line)
    dot = always_redraw(cur_dot)
    s.play(FadeIn(ax), FadeIn(head), run_time=0.8)
    s.say("補充 1：沿方向 ", "$\\vec n$", " 量「平方寬度」", "$\\vec n\\cdot\\mathbf T\\cdot\\vec n$", wait=2.6)
    s.add(line, dot)
    s.say("方向轉半圈，寬度跟著變", wait=0.2)
    s.play(phi.animate.set_value(180), ShowCreation(graph), run_time=5.0, rate_func=linear)
    s.play(phi.animate.set_value(p1), run_time=1.2)
    s.play(FadeIn(d1), FadeIn(m1), run_time=0.6)
    s.say("最大 ", "$\\frac35$", " 在 ", "$%s^\\circ$" % nb["phi_deg"], "：最長的方向（綠軸）", wait=2.9)
    s.play(phi.animate.set_value(p2), run_time=1.2)
    s.play(FadeIn(d2), FadeIn(m2), run_time=0.6)
    s.say("最小 ", "$\\frac15$", " 在 ", "$%s^\\circ$" % nb["phi2_deg"], "：最短的方向（粉軸）", wait=2.9)
    s.say("特徵向量＝最長、最短的方向；特徵值＝那個方向的平方寬度", wait=3.2)
    s.remove(line, dot)
    s.clear_stage()


def fig64_group(unit=0.10, center=(-3.5, 0.2)):
    nb = numbers()
    p = nb["pts"]
    c = (p.min(axis=0) + p.max(axis=0)) / 2
    walk = VMobject().set_points_as_corners(to_screen(p, unit, c, center)).set_stroke(YELLOW, 1.6)
    cs = to_screen([nb["com"]], unit, c, center)[0]
    axes = principal_axes(cs, nb["Rraw"][0], nb["Rraw"][1], nb["vec1"], unit)
    return walk, axes, cross_mark(cs), cs


def part_fig64(s):
    nb = numbers()
    walk, axes, cross, cs = fig64_group()
    a, b, c = nb["T"]
    r1 = at(mixed("重心 ", "$(%s,\\ %s)$" % nb["com_s"], font_size=30), 3.6, 2.5)
    r1[1].set_color(COM_COL)
    mat = at(tx(r"\mathbf T=\begin{pmatrix}%s&%s\\%s&%s\end{pmatrix}" % (a, c, c, b), 34), 3.6, 1.35)
    r3 = at(tx(r"R_1^2=%s,\quad R_2^2=%s" % nb["eig"], 32), 3.6, 0.2)
    r4 = at(tx(r"R_1=%s,\ \ R_2=%s,\ \ R_1^2/R_2^2=%s" % (nb["R"][0], nb["R"][1], nb["ratio"]), 28), 3.6, -0.75)
    ang = at(mixed("主軸 ", "$%s^\\circ$" % nb["angle"], font_size=30, color=AX1_COL), 3.6, -1.75)
    for m in (r1, mat, r3, r4, ang):
        fit(m, 5.8)
    left = VGroup(walk, axes, cross)
    assert_in_frame(left, r1, mat, r3, r4, ang)
    for x, y in ((left, r1), (left, mat), (left, r3), (left, r4), (left, ang), (r1, mat), (mat, r3), (r3, r4),
                 (r4, ang)):
        assert_apart(x, y)
    s.play(ShowCreation(walk), run_time=3.0)
    s.say("主例（Fig. 6.4）：二維 1000 步", wait=1.5)
    s.play(ShowCreation(cross), FadeIn(r1), run_time=0.6)
    s.say("重心 ", "$(%s,\\ %s)$" % nb["com_s"], wait=2.4)
    s.play(FadeIn(mat), run_time=0.8)
    s.say("張量：", "$T_{11}=%s$" % a, "、", "$T_{22}=%s$" % b, "、", "$T_{12}=%s$" % c, wait=2.9)
    s.play(FadeIn(r3), FadeIn(r4), run_time=0.8)
    s.say("特徵值 %s 與 %s：長 %s、寬 %s" % (nb["eig"] + nb["R"]), wait=2.9)
    s.play(FadeIn(axes[0]), ShowCreation(axes[1]), ShowCreation(axes[2]), FadeIn(ang), run_time=1.2)
    s.say("兩條主軸（半長 ", "$2R_i$", "）：長軸在 ", "$%s^\\circ$" % nb["angle"], wait=2.9)
    s.say("平方比 %s：又是一條長條" % nb["ratio"], wait=2.4)
    return dict(walk=walk, axes=axes, cross=cross, cs=cs, mat=mat, r3=r3, ang=ang, r1=r1, r4=r4)


def part_rotate(s, g):
    nb = numbers()
    a, b, c = nb["Trot"]
    mat2 = tx(r"\mathbf T'=\begin{pmatrix}%s&%s\\%s&%s\end{pmatrix}" % (a, c, c, b), 34).move_to(g["mat"])
    ang2 = mixed("主軸 ", "$%s^\\circ$" % nb["angle_rot"], font_size=30, color=AX1_COL).move_to(g["ang"])
    e65 = eq_label(tx(r"\langle\mathbf T\rangle=\bar T\,\mathbf I,\qquad \bar T=\tfrac12\,\mathrm{Tr}\,\mathbf T=%s"
                      % nb["Tbar"], 30), "(6.5)")
    at(e65, 3.6, -2.55)
    fit(e65, 5.9)
    fit(mat2, 5.8)
    rot = VGroup(g["walk"], g["axes"], g["cross"])
    rotated = rot.copy().rotate(THETA, about_point=g["cs"])
    assert_in_frame(rotated, mat2, ang2, e65)
    for x, y in ((rotated, mat2), (rotated, g["r3"]), (rotated, ang2), (rotated, g["r1"]), (rotated, g["r4"]),
                 (e65, rotated), (e65, ang2)):
        assert_apart(x, y)
    s.say("把整條漫步繞重心轉 ", "$30^\\circ$", wait=0.2)
    s.play(Rotate(rot, THETA, about_point=g["cs"]), run_time=2.5)
    s.play(FadeTransform(g["mat"], mat2), run_time=0.8)
    s.say("矩陣全變了：", "$T'_{11}=%s$" % a, "、", "$T'_{22}=%s$" % b, "、", "$T'_{12}=%s$" % c, wait=2.9)
    s.play(Indicate(g["r3"], color=YELLOW), run_time=1.0)
    s.say("特徵值不變：還是 %s、%s" % nb["eig_rot"], wait=2.6)
    s.play(FadeTransform(g["ang"], ang2), run_time=0.8)
    s.say("主軸跟著轉：", "$%s^\\circ$" % nb["angle"], " 變成 ", "$%s^\\circ$" % nb["angle_rot"], wait=2.6)
    s.play(FadeIn(e65, shift=0.2 * UP), run_time=0.8)
    s.say("對所有方向平均：", "$\\langle\\mathbf T\\rangle=\\bar T\\,\\mathbf I$", wait=2.6)
    s.say("平均之後只剩一個數 %s：形狀資訊被抹掉" % nb["Tbar"], wait=2.9)
    s.clear_stage()


def part_close(s):
    s.say("張量會跟著轉，特徵值不會：主迴轉半徑就是形狀的指紋", wait=3.4)


class W602GyrationTensor(ChuScene):
    def construct(self):
        self.show_header("W602　迴轉半徑張量：形狀的指紋")
        check_math()
        part_question(self)
        sp = part_tensor(self)
        part_extent(self, sp)
        g = part_fig64(self)
        part_rotate(self, g)
        part_close(self)
        self.wait(1.0)
