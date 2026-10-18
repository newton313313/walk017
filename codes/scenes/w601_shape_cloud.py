"""
W601 漫步長什麼樣子？（§6.1–6.1.1，pp. 127–129，Fig. 6.1–6.3）

問：1000 個醉漢各走 100 步（二維），疊起來是一團圓雲。單獨一條也是圓的嗎？
  1. 小例子：4 步「東、東、北、西」→ 5 個點、重心 (6/5, 2/5)、平方距離 → R_g² = 4/5。
  2. Fig. 6.1 的雲：平均張量 (6.5) 8.35、8.59、0.03（理論 8.42、8.42、0），特徵值 8.60、8.34。
  3. 揭曉 Fig. 6.3：6 條 1000 步的 R_1²/R_2²＝2.1 … 14.5；1000 條 100 步平均 14.14 : 2.80（比 5.05）。
  4. 疊在一起：長軸指向不同方向 → 要用不隨旋轉改變的量。
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

N_SHOW = 200                     # 雲裡畫出的條數
CLOUD_COL = BLUE_B
COM_COL = ORANGE
D2_COL = TEAL_B
TILE_COLS = [YELLOW, TEAL_B, PINK, "#7CE38B", BLUE_B, ORANGE]


# ---------------------------------------------------------------- 數學
def f2(x):
    return f"{float(x):.2f}"


def frac_tex(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return ("-" if q < 0 else "") + r"\frac{%d}{%d}" % (abs(q.numerator), q.denominator)


@lru_cache(maxsize=None)
def numbers():
    se = SH.small_example()
    cs = SH.cloud_stats()
    pts = SH.cloud()
    Tm = cs["Tmean"]
    nb = dict(se=se, cloud=pts,
              Tm=(f2(Tm[0, 0]), f2(Tm[1, 1]), f2(Tm[0, 1])), Tth=f2(cs["exact_T11"]),
              eig=tuple(f2(v) for v in cs["eig_of_mean"]), eig_raw=tuple(float(v) for v in cs["eig_of_mean"]),
              vec_mean=SH.principal(Tm)[1][:, 0],
              R=(f2(cs["R1"]), f2(cs["R2"]), f2(cs["ratio"])),
              ratios=tuple(f"{r:.1f}" for r in SH.fig63_ratios()))
    walks = []
    for p in SH.fig63():
        vals, vecs = SH.principal(SH.tensor(p))
        walks.append(dict(pts=p, com=SH.center_of_mass(p), R=np.sqrt(vals),
                          angle=math.atan2(vecs[1, 0], vecs[0, 0])))
    nb["walks"] = walks
    return nb


def check_math():
    nb = numbers()
    se = nb["se"]
    assert se["pts"] == [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1)]
    assert list(map(tuple, SH.SMALL_STEPS)) == [(1, 0), (1, 0), (0, 1), (-1, 0)]
    assert se["com"] == [Fraction(6, 5), Fraction(2, 5)]
    assert se["d2"] == [Fraction(8, 5), Fraction(1, 5), Fraction(4, 5), Fraction(1), Fraction(2, 5)]
    assert se["Rg2"] == Fraction(4, 5) and sum(se["d2"]) / 5 == Fraction(4, 5)
    T = se["T"]
    assert T[0][0] + T[1][1] == se["Rg2"]                                  # R_g² = Tr T
    pts = nb["cloud"]
    assert pts.shape == (1000, 101, 2) and int(np.abs(pts).max()) <= 30
    assert nb["Tm"] == ("8.35", "8.59", "0.03") and nb["Tth"] == "8.42"
    assert SH.cloud_stats()["exact_T11"] == Fraction(101 ** 2 - 1, 12 * 101)              # ((N+1)²−1)/(12(N+1))
    assert max(nb["ratios"], key=float) == "14.5"
    assert nb["eig"] == ("8.60", "8.34")
    assert nb["R"] == ("14.14", "2.80", "5.05")
    assert nb["ratios"] == ("2.1", "3.0", "6.7", "5.5", "14.5", "4.6")
    assert len(nb["walks"]) == 6 and all(len(w["pts"]) == 1001 for w in nb["walks"])
    assert all(float(r) > 2 for r in nb["ratios"])                          # 沒有一條是圓的


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


def path_mob(P, unit, c, origin, color=YELLOW, width=2.0, opacity=1.0):
    m = VMobject().set_points_as_corners(to_screen(P, unit, c, origin))
    return m.set_stroke(color, width, opacity)


def cross_mark(p, size=0.16, color=COM_COL, width=4):
    a = Line(p + size * (UL), p + size * (DR)).set_stroke(color, width)
    b = Line(p + size * (UR), p + size * (DL)).set_stroke(color, width)
    return VGroup(a, b)


def shape_ellipse(center, R1, R2, angle, unit, color=COM_COL, width=2.5):
    e = Ellipse(width=4 * R1 * unit, height=4 * R2 * unit)          # 半軸 2R_i
    e.rotate(angle).move_to(center)
    return e.set_stroke(color, width).set_fill(color, 0.08)


def cloud_group(unit, origin):
    pts = numbers()["cloud"]
    return VGroup(*[path_mob(pts[m], unit, (0, 0), origin, CLOUD_COL, 1.0, 0.28) for m in range(N_SHOW)])


# ---------------------------------------------------------------- 各段
def part_question(s):
    nb = numbers()
    unit, org = 0.09, np.array([-3.3, 0.15, 0])
    cloud = cloud_group(unit, org)
    one = path_mob(nb["cloud"][0], unit, (0, 0), org, YELLOW, 3.0)
    t1 = at(mixed("1000 個醉漢、各走 100 步", font_size=30), 3.3, 1.9)
    t2 = at(zh("疊起來：一團圓雲", font_size=30, color=CLOUD_COL), 3.3, 1.0)
    t3 = at(zh("單獨一條也是圓的嗎？", font_size=34, color=YELLOW), 3.3, -0.4)
    assert_in_frame(cloud, t1, t2, t3)
    for a, b in ((cloud, t1), (cloud, t2), (cloud, t3), (t1, t2), (t2, t3)):
        assert_apart(a, b)
    s.play(FadeIn(t1), run_time=0.6)
    s.say("問：1000 個醉漢各走 100 步（二維）", wait=2.0)
    s.say("每一步：往東、西、南、北之一走 1 格", wait=0.3)
    s.play(LaggedStart(*[ShowCreation(m) for m in cloud], lag_ratio=0.01), run_time=4.0)
    s.play(FadeIn(t2), run_time=0.5)
    s.say("疊起來是一團圓圓的雲（畫出前 200 條）", wait=2.4)
    s.play(ShowCreation(one), run_time=2.0)
    s.play(FadeIn(t3), run_time=0.5)
    s.say("單獨一條（黃）也是圓的嗎？先猜", wait=3.4)
    s.clear_stage()


def part_small(s):
    nb = numbers()
    se = nb["se"]
    unit, org = 1.45, np.array([-5.0, -0.9, 0])
    com = np.array([float(v) for v in se["com"]])
    scr = to_screen(se["pts"], unit, (0, 0), org)
    com_p = to_screen([com], unit, (0, 0), org)[0]
    dots_bg = VGroup(*[Dot(to_screen([(x, y)], unit, (0, 0), org)[0], radius=0.035).set_fill(GREY_D, 1)
                       for x in range(-1, 4) for y in range(-1, 3)])
    arrows = VGroup(*[Arrow(scr[k], scr[k + 1], buff=0.1, thickness=4).set_color(YELLOW) for k in range(4)])
    pdots = VGroup(*[Dot(p, radius=0.08).set_fill(WHITE, 1) for p in scr])
    cross = cross_mark(com_p)
    dashes = VGroup(*[DashedLine(p, com_p, dash_length=0.08).set_stroke(D2_COL, 2.5) for p in scr])
    labs = VGroup()
    for p, q, d2 in zip(scr, se["pts"], se["d2"]):
        v = np.array([q[0] - com[0], q[1] - com[1], 0])
        v /= np.linalg.norm(v)
        labs.add(tx(frac_tex(d2), 28, D2_COL).move_to(p + 0.48 * v))
    pic = VGroup(dots_bg, arrows, pdots, cross, dashes, labs)
    r1 = at(mixed("4 步：東、東、北、西", font_size=30), 3.0, 2.6)
    r2 = at(tx(r"(0,0),\ (1,0),\ (2,0),\ (2,1),\ (1,1)", 30), 3.0, 1.85)
    r3 = at(mixed("重心 ", "$\\langle\\vec r\\rangle=\\Big(\\frac{6}{5},\\frac{2}{5}\\Big)$", font_size=30,
                  t2c=None), 3.0, 0.95)
    r3[1].set_color(COM_COL)
    r4 = at(tx(r"R_g^2=\frac{1}{5}\Big(\frac{8}{5}+\frac{1}{5}+\frac{4}{5}+1+\frac{2}{5}\Big)=\frac{4}{5}", 32), 3.0,
            -0.35)
    r5 = at(zh("迴轉半徑的平方：轉動慣量除以總質量", font_size=26, color=GREY_A), 3.0, -1.5)
    for m in (r1, r2, r3, r4, r5):
        fit(m, 6.6)
    assert_in_frame(pic, r1, r2, r3, r4, r5)
    for a, b in ((pic, r1), (pic, r2), (pic, r3), (pic, r4), (pic, r5), (r1, r2), (r2, r3), (r3, r4), (r4, r5)):
        assert_apart(a, b)
    for i in range(5):
        for j in range(i + 1, 5):
            assert_apart(labs[i], labs[j])
        assert_apart(labs[i], cross)
    s.play(FadeIn(dots_bg), FadeIn(r1), run_time=0.6)
    s.say("先看一個小例子：4 步，東、東、北、西", wait=0.3)
    s.play(FadeIn(pdots[0]), run_time=0.3)
    for k in range(4):
        s.play(GrowArrow(arrows[k]), FadeIn(pdots[k + 1]), run_time=0.5)
    s.play(FadeIn(r2), run_time=0.6)
    s.say("走過 5 個點（", "$N$", " 步有 ", "$N+1$", " 個點）", wait=2.9)
    s.play(ShowCreation(cross), FadeIn(r3), run_time=0.8)
    s.say("重心：5 個點的平均", wait=2.6)
    s.play(LaggedStart(*[ShowCreation(d) for d in dashes], lag_ratio=0.15), run_time=1.2)
    s.play(LaggedStart(*[FadeIn(m) for m in labs], lag_ratio=0.15), run_time=1.0)
    s.say("每個點到重心的平方距離", wait=2.9)
    s.play(Write(r4), run_time=1.2)
    s.say("平均起來：", "$R_g^2=\\frac45$", "，迴轉半徑的平方", wait=2.8)
    s.say("分母是 5：", "$N+1$", " 個點，不是 ", "$N$", " 步", wait=2.6)
    s.play(FadeIn(r5), run_time=0.6)
    s.say("像轉動慣量：點離重心越遠，貢獻越大", wait=2.8)
    s.clear_stage()


def part_cloud(s):
    nb = numbers()
    unit, org = 0.09, np.array([-3.3, 0.15, 0])
    cloud = cloud_group(unit, org)
    v = nb["vec_mean"]
    l1, l2 = nb["eig_raw"]
    circ = shape_ellipse(org, math.sqrt(l1), math.sqrt(l2), math.atan2(v[1], v[0]), unit, COM_COL, 3)
    a, b, c = nb["Tm"]
    e1 = eq_label(tx(r"\bar{\mathbf T}=\frac{1}{1000}\sum_{m=1}^{1000}\mathbf T^{(m)}", 32), "(6.5)")
    e2 = tx(r"\bar{\mathbf T}=\begin{pmatrix}%s&%s\\%s&%s\end{pmatrix}" % (a, c, c, b), 34)
    th = nb["Tth"]
    e3 = mixed("理論：", "$\\begin{pmatrix}%s&0\\\\0&%s\\end{pmatrix}$" % (th, th), font_size=30)
    e3b = tx(r"\bar T_{11}=\bar T_{22}=\frac{(N+1)^2-1}{12\,(N+1)}=%s" % th, 28)
    e4 = mixed("特徵值 %s、%s：幾乎是圓" % nb["eig"], font_size=30, color=COM_COL)
    for m, y in ((e1, 2.55), (e2, 1.35), (e3, 0.05), (e3b, -0.95), (e4, -1.95)):
        at(m, 3.4, y)
        fit(m, 6.0)
    assert_in_frame(cloud, circ, e1, e2, e3, e3b, e4)
    for p, q in ((cloud, e1), (cloud, e2), (cloud, e3), (cloud, e3b), (cloud, e4), (e1, e2), (e2, e3), (e3, e3b),
                 (e3b, e4)):
        assert_apart(p, q)
    s.play(FadeIn(cloud), run_time=0.8)
    s.say("回到 1000 條：每條算一個張量 ", "$\\mathbf T$", "，再平均", wait=0.3)
    s.play(Write(e1), run_time=1.0)
    s.play(FadeIn(e2, shift=0.2 * DOWN), run_time=0.8)
    s.say("平均張量：對角 %s、%s，非對角只有 %s" % (a, b, c), wait=2.8)
    s.play(FadeIn(e3), run_time=0.6)
    s.say("理論值：對角都是 %s，非對角是 0" % th, wait=2.6)
    s.play(FadeIn(e3b), run_time=0.6)
    s.say("（", "$N=100$", "；W605 的 Kramers 恆等式會導出這個式子）", wait=2.6)
    s.play(FadeIn(e4), ShowCreation(circ), run_time=0.8)
    s.say("兩個特徵值 %s、%s：平均的形狀幾乎是圓" % nb["eig"], wait=2.8)
    s.say("但這是「先平均、再看形狀」", wait=2.6)
    s.clear_stage()


def tiles_layout():
    nb = numbers()
    xs, ys = (-4.35, 0.0, 4.35), (2.1, -0.45)
    W, H = 3.9, 1.95
    out = []
    for i, w in enumerate(nb["walks"]):
        cx, cy = xs[i % 3], ys[i // 3]
        P = w["pts"]
        a, b, th = 2 * w["R"][0], 2 * w["R"][1], w["angle"]
        half = np.array([math.hypot(a * math.cos(th), b * math.sin(th)), math.hypot(a * math.sin(th), b * math.cos(th))])
        lo = np.minimum(P.min(axis=0), w["com"] - half)
        hi = np.maximum(P.max(axis=0), w["com"] + half)
        unit = min((W - 0.25) / (hi[0] - lo[0]), (H - 0.25) / (hi[1] - lo[1]))
        c = (lo + hi) / 2
        org = np.array([cx, cy, 0])
        frame = Rectangle(width=W, height=H).set_stroke(GREY_D, 1.5).move_to(org)
        path = path_mob(P, unit, c, org, TILE_COLS[i], 1.6)
        com_p = to_screen([w["com"]], unit, c, org)[0]
        ell = shape_ellipse(com_p, w["R"][0], w["R"][1], w["angle"], unit, WHITE, 2)
        lab = mixed("$R_1^2/R_2^2=%s$" % nb["ratios"][i], font_size=26)
        lab.next_to(frame, DOWN, buff=0.12)
        out.append(dict(frame=frame, path=path, ell=ell, lab=lab, unit=unit, com=com_p))
    return out


def part_reveal(s):
    nb = numbers()
    tiles = tiles_layout()
    R1, R2, rat = nb["R"]
    summ = at(mixed("1000 條 100 步的平均：", "$\\langle R_1^2\\rangle=%s$" % R1, "，",
                    "$\\langle R_2^2\\rangle=%s$" % R2, "（比 %s）" % rat, font_size=28, color=YELLOW), 0, -2.5)
    fit(summ, 13.0)
    allm = [VGroup(t["frame"], t["path"], t["ell"]) for t in tiles]
    assert_in_frame(*allm, *[t["lab"] for t in tiles], summ)
    for i in range(6):
        for j in range(i + 1, 6):
            assert_apart(allm[i], allm[j])
            assert_apart(tiles[i]["lab"], allm[j]) if j != i else None
        assert_apart(tiles[i]["lab"], allm[i])
        assert_apart(tiles[i]["lab"], summ)
        assert_apart(allm[i], summ)
    s.say("揭曉：單獨看 6 條，各走 1000 步（各自縮放）", wait=0.3)
    s.play(*[FadeIn(t["frame"]) for t in tiles], run_time=0.5)
    s.play(*[ShowCreation(t["path"]) for t in tiles], run_time=3.5)
    s.say("白色橢圓：每條自己的形狀（半軸 ", "$2R_1$", "、", "$2R_2$", "）", wait=0.3)
    s.play(*[ShowCreation(t["ell"]) for t in tiles], run_time=1.2)
    s.wait(2.0)
    s.play(LaggedStart(*[FadeIn(t["lab"]) for t in tiles], lag_ratio=0.15), run_time=1.2)
    s.say("長短軸平方比 ", "$R_1^2/R_2^2$", "：從 2.1 到 14.5", wait=2.8)
    s.say("沒有一條是圓的：典型的一條是長條", wait=2.6)
    s.say("最長的一條：長軸平方是短軸平方的 %s 倍" % max(nb["ratios"], key=float), wait=2.6)
    s.play(FadeIn(summ), run_time=0.6)
    s.say("1000 條 100 步平均：", "$\\langle R_1^2\\rangle$", " 是 ", "$\\langle R_2^2\\rangle$",
          " 的 %s 倍" % rat, wait=3.0)
    return tiles


def part_overlay(s, tiles):
    nb = numbers()
    org = np.array([-3.0, 0.15, 0])
    u = 0.042
    over = VGroup()
    ells = VGroup()
    for i, w in enumerate(nb["walks"]):
        over.add(path_mob(w["pts"], u, w["com"], org, TILE_COLS[i], 1.4, 0.8))
        ells.add(shape_ellipse(org, w["R"][0], w["R"][1], w["angle"], u, TILE_COLS[i], 3).set_fill(opacity=0))
    t1 = at(zh("雲是圓的：", font_size=30), 3.6, 1.6)
    t2 = at(zh("每條往不同方向長", font_size=30), 3.6, 0.95)
    t3 = at(zh("一條的形狀＝不隨旋轉改變的量", font_size=28, color=YELLOW), 3.6, -0.4)
    t4 = at(zh("（下一幕：迴轉半徑張量）", font_size=24, color=GREY_B), 3.6, -1.1)
    assert_in_frame(over, ells, t1, t2, t3, t4)
    for a, b in ((over, t1), (over, t2), (over, t3), (over, t4), (ells, t3), (t1, t2), (t2, t3), (t3, t4)):
        assert_apart(a, b)
    keep = VGroup(*[t["path"] for t in tiles])
    rest = [m for m in s.mobjects if m is not s.frame and m is not s.header_mob and m is not s.caption_mob
            and m not in keep]
    s.play(*[FadeOut(m) for m in rest], run_time=0.6)
    s.play(*[Transform(tiles[i]["path"], over[i]) for i in range(6)], run_time=2.0)
    s.say("把 6 條疊在一起（同一比例、重心對齊）", wait=1.5)
    s.play(ShowCreation(ells), run_time=1.2)
    s.say("各自的長軸指向不同方向", wait=2.8)
    s.play(FadeIn(t1), FadeIn(t2), run_time=0.6)
    s.say("雲是圓的，因為每條往不同方向長", wait=2.8)
    s.say("平均張量把每條的方向都平均掉了", wait=2.6)
    s.play(FadeIn(t3), FadeIn(t4), run_time=0.6)
    s.say("要描述一條的形狀，得用不隨旋轉改變的量", wait=3.0)
    s.clear_stage()


def part_close(s):
    s.say("一群醉漢疊起來是圓的，單獨一個卻幾乎總是長條", wait=3.2)


class W601ShapeCloud(ChuScene):
    def construct(self):
        self.show_header("W601　漫步長什麼樣子？")
        check_math()
        part_question(self)
        part_small(self)
        part_cloud(self)
        tiles = part_reveal(self)
        part_overlay(self, tiles)
        part_close(self)
        self.wait(1.0)
