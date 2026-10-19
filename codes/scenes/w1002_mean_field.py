"""
W1002 O(n) 模型的平均場（§10.2，pp. 236–239，(10.2)–(10.15)，註 2–3）

問：如果每個自旋都一樣（平均場），n → 0 的 O(n) 模型會給出什麼？
  1. (10.2)–(10.5)：能量、G = lim (1/n) Σ⟨S S⟩ = Σ C(N) e^{−rN}、Γ(r) = Σ Γ_N e^{−rN}；這章 e^{−r} = z/z_c（第 9 章 e^{−rN/2}）。
  2. (10.8)–(10.12)：S_i = M ⇒ Γ_MF = N0 lim⟨M1²⟩；插入 δ 函數，n 維積分變一維；表 Γ_n（n = 3, 2, 1, 0.1, 0.01）→ 1.8682。
  3. (10.13)–(10.15)：對 S 積分（註 2，極點在實軸上方）= 2πi／0（ε = 0.1：5.6853i、0）⇒ Γ_MF = (N0/2) ∫_0^∞ f；
     筆誤框 1（p. 238，(10.15)：N0/4 → N0/2）：u = 0 時應得 1/r = 2（書上 1）；表的極限 1.8682（書上 0.9341）。
  4. f(λ) 兩條曲線（r = ±0.5）；r > 0：Γ_MF → 1/r（1.8682、1.9844、1.9984）；r < 0：精確 ÷ 漸近 = 0.9998、1.0000、1.0000，
     ln Γ_MF / N0 = 0.0912、0.0665、0.0630 → r²/4u = 0.0625；極限的順序：Γ_n/Γ_MF = 0.3299、0.8311、0.9801、0.9998。
  5. 註 3：(10.15) 與 (7.57) 一項項換（N0 → N、λ → 2λ、u → f、r → 1 − z/z_c）⇒ Γ_MF = N G_MF；1.8682 = 100 G_MF(1/2)；
     筆誤框 2（p. 238，footnote 3）：書上的換法只對它自己的 N0/4 成立（0.9341 = 0.9341）。
畫面上的數字都由 utils.walk_cond 算出，check_math() 核對。
模組層級的畫圖函式：spin_row、s_plane、f_curve。
"""
import inspect
import math
import sys
import warnings
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import mixed, zh
from utils import walk_cond as WK
from utils import walk_draw as D

NS = (3, 2, 1, 0.1, 0.01)
GN = ("1.7183", "1.7634", "1.8131", "1.8624", "1.8676")
N0S = (100, 1000, 10000)
N0_TEX = ("100", "1000", "10^4")
POS = ("1.8682", "1.9844", "1.9984")
RATIO = ("0.9998", "1.0000", "1.0000")
LNG = ("0.0912", "0.0665", "0.0630")
NS_NEG = (1e-2, 1e-3, 1e-4, 1e-6)
NS_NEG_TEX = ("10^{-2}", "10^{-3}", "10^{-4}", "10^{-6}")
ORDER = ("0.3299", "0.8311", "0.9801", "0.9998")
R_POS = WK.MF_R                            # 0.5
R_NEG = -WK.MF_R                           # −0.5
EPS = 0.1
POS_COLOR = "#4A90D9"
NEG_COLOR = "#F2A33A"
GOOD = "#7CE38B"
SUB_COLORS = (YELLOW, TEAL_B, "#F2A33A", "#7CE38B")
F_TEX = (
    r"\Gamma_{MF}=\frac{N_0}{2}\int_0^\infty e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}d\lambda",
    r"\Gamma_{MF}=\frac{N}{2}\int_0^\infty e^{-N(\frac r2\lambda+\frac u4\lambda^2)}d\lambda",
    r"\Gamma_{MF}=N\int_0^\infty e^{-N(r\lambda+u\lambda^2)}d\lambda",
    r"\Gamma_{MF}=N\int_0^\infty e^{-N(r\lambda+f\lambda^2)}d\lambda",
    r"\Gamma_{MF}=N\int_0^\infty e^{N(-f\lambda^2-\lambda(1-x))}d\lambda=N\,G_{MF}",
)
F_KEYS = ("N", r"\lambda", r"f\lambda", "1-x")       # 每一步換掉的那一項（上色）
RULES = (r"N_0\to N", r"\lambda\to2\lambda", r"u\to f", r"r\to1-z/z_c")
X0 = -1.0                                  # 第 4 段右欄的左緣


def f4(x):
    return f"{float(x):.4f}"


@lru_cache(maxsize=None)
def numbers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return dict(
            gn=tuple(WK.gamma_n(n) for n in NS), g0=WK.gamma_mf(R_POS), book=WK.gamma_mf(R_POS, book=True),
            s_pos=WK.s_integral(1.0), s_neg=WK.s_integral(-1.0), gauss=WK.gamma_gauss_check(),
            pos=tuple(WK.gamma_mf(R_POS, 1, N0) for N0 in N0S),
            ratio=tuple(WK.gamma_mf(R_NEG, 1, N0) / WK.gamma_mf_asym(R_NEG, 1, N0) for N0 in N0S),
            lng=tuple(WK.ln_gamma_mf(R_NEG, 1, N0) / N0 for N0 in N0S),
            lng_asym=tuple(WK.ln_gamma_mf_asym(R_NEG, 1, N0) / N0 for N0 in N0S),
            order=tuple(WK.gamma_n(n, R_NEG) / WK.gamma_mf(R_NEG) for n in NS_NEG),
            fn3=tuple(WK.fn3_check()))


def check_math():
    nb = numbers()
    assert (WK.MF_N0, WK.MF_R, WK.MF_U) == (100, 0.5, 1.0) and N0S[0] == WK.MF_N0
    assert [f4(g) for g in nb["gn"]] == list(GN) and f4(nb["g0"]) == "1.8682" and f4(nb["book"]) == "0.9341"
    assert inspect.signature(WK.s_integral).parameters["eps"].default == EPS
    sp, sn = nb["s_pos"], nb["s_neg"]
    assert f4(sp.imag) == "5.6853" and abs(sp.real) < 1e-6 and abs(sp - 2j * math.pi * math.exp(-EPS)) < 1e-5
    assert abs(sn) < 1e-6 and f"{abs(sn):.0f}" == "0"
    gns, gmf, gbook = nb["gauss"]
    assert all(f4(g) == "2.0000" for g in gns) and f4(gmf) == "2.0000" and f4(gbook) == "1.0000"
    assert f"{1 / R_POS:g}" == "2"
    assert [f4(v) for v in nb["pos"]] == list(POS)
    assert [f4(v) for v in nb["ratio"]] == list(RATIO)
    assert [f4(v) for v in nb["lng"]] == list(LNG) == [f4(v) for v in nb["lng_asym"]]
    assert f4(R_NEG ** 2 / (4 * WK.MF_U)) == "0.0625"
    assert [f4(v) for v in nb["order"]] == list(ORDER)
    assert [f4(v) for v in nb["fn3"]] == ["1.8682"] * 3 + ["0.9341"] * 2
    assert 1 - R_POS == 0.5                                    # x = z/z_c = 1/2
    lam = np.linspace(0, 1, 201)
    fn = [WK.mf_f(x, R_NEG) for x in lam]
    fp = [WK.mf_f(x, R_POS) for x in lam]
    assert abs(lam[int(np.argmax(fn))] - abs(R_NEG) / WK.MF_U) < 1e-9 and f"{abs(R_NEG) / WK.MF_U:g}" == "0.5"
    assert all(a > b for a, b in zip(fp, fp[1:]))


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


def place_left(m, x0, y, x1=6.7):
    """左緣對齊 x0、中心高度 y；太寬就縮到 x1 為止。"""
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


# ---------------------------------------------------------------- 畫圖元件
def spin_row(x0, x1, y, count=8, length=0.75, color=YELLOW):
    """一排完全相同（朝上）的自旋：平均場。"""
    return VGroup(*[Arrow(np.array([x, y, 0]), np.array([x, y + length, 0]), buff=0, thickness=4,
                          max_width_to_length_ratio=0.16).set_fill(color, 1) for x in np.linspace(x0, x1, count)])


def s_plane(center, half_w=2.3, R=1.25, eps_h=0.4):
    """S 平面示意：實軸、極點 S = iε（實軸上方）、λ > 0 往上包、λ < 0 往下包。
    → VGroup(實軸, Re S, 上半圓, 下半圓, ×, 極點標籤, λ>0, λ<0)"""
    c = np.array(center, dtype=float)
    axis = Arrow(c + half_w * LEFT, c + half_w * RIGHT, buff=0, thickness=2.5,
                 max_width_to_length_ratio=0.05).set_fill(GREY_B, 1)
    re_lab = tx(r"\mathrm{Re}\,S", 24, GREY_B).next_to(axis.get_right(), DOWN, buff=0.12)
    pole = tx(r"\times", 36, RED).move_to(c + eps_h * UP)
    pole_lab = tx(r"S=i\epsilon", 24, RED).next_to(pole, UP, buff=0.06)
    up = Arc(start_angle=0, angle=PI, radius=R, arc_center=c).set_stroke(GOOD, 3)
    dn = Arc(start_angle=0, angle=-PI, radius=R, arc_center=c).set_stroke(POS_COLOR, 3)
    d45 = R / math.sqrt(2)
    up_lab = tx(r"\lambda>0", 26, GOOD).next_to(c + d45 * (RIGHT + UP), UR, buff=0.06)
    dn_lab = tx(r"\lambda<0", 26, POS_COLOR).next_to(c + d45 * (RIGHT + DOWN), DR, buff=0.06)
    return VGroup(axis, re_lab, up, dn, pole, pole_lab, up_lab, dn_lab)


def f_curve(ax, r, color, n=60):
    """f(λ)/峰值，n 點折線（r > 0 在 λ = 0 附近取密一點）。"""
    lams = [(k / (n - 1)) ** 2 for k in range(n)] if r > 0 else [k / (n - 1) for k in range(n)]
    peak = WK.mf_f(0.0 if r > 0 else abs(r) / WK.MF_U, r)
    pts = [ax.c2p(x, WK.mf_f(x, r) / peak) for x in lams]
    return VMobject().set_points_as_corners(pts).set_stroke(color, 4)


# ---------------------------------------------------------------- 各段
def part_model(s):
    rows = (("(10.2)", r"\mathcal H=\sum_i\Big\{\frac r2|\vec S_i|^2+\frac12|\nabla\vec S_i|^2"
                       r"+\frac u4|\vec S_i|^4-\vec h\cdot\vec S_i\Big\}"),
            ("(10.3)(10.4)", r"G=\lim_{n\to0}\frac1n\sum_\alpha\langle S_{i\alpha}S_{j\alpha}\rangle"
                             r"=\sum_N C(N)\,e^{-rN}"),
            ("(10.5)", r"\Gamma(r)=\sum_N\Gamma_N\,e^{-rN}"))
    ys = (2.2, 0.75, -0.7)
    tags, eqs = VGroup(), VGroup()
    for (tag, tex), y in zip(rows, ys):
        tags.add(zh(tag, font_size=22, color=GREY_B).move_to([-5.95, y, 0]))
        eqs.add(place_left(tx(tex, 38), -5.0, y, 6.6))
    note = zh("（對終點加總＝磁化率）", font_size=28).next_to(eqs[2], RIGHT, buff=0.4)
    assert_in_frame(tags, eqs, note)
    for a, b in zip(eqs, eqs[1:]):
        assert_apart(a, b)
    for t, e in zip(tags, eqs):
        assert_apart(t, e)

    s.say("如果每個自旋都一樣（平均場）", wait=1.5)
    s.say(r"$n\to0$", " 的 ", r"$\mathrm{O}(n)$", " 模型會給出什麼？", wait=2.0)
    s.play(FadeIn(tags[0]), Write(eqs[0]), run_time=1.5)
    s.say("(10.2)：", r"$\mathrm{O}(n)$", " 模型的能量，", r"$\vec S_i$", " 有 ", "$n$", " 個分量", wait=2.0)
    s.play(FadeIn(tags[1]), Write(eqs[1]), run_time=1.5)
    s.say(r"$n\to0$", " 時關聯函數就是自迴避漫步的加總", wait=2.0)
    s.say("這章 ", r"$e^{-r}=z/z_c$", "（第 9 章是 ", r"$e^{-rN/2}$", "）", wait=2.0)
    s.play(FadeIn(tags[2]), Write(eqs[2]), FadeIn(note), run_time=1.2)
    s.say("(10.5)：對終點加總，", r"$\Gamma(r)$", " 就是磁化率", wait=2.0)
    s.clear_stage()


def part_mean_field(s):
    nb = numbers()
    row = spin_row(-6.2, -3.0, 1.95)
    eq1 = place_in(tx(r"\vec S_i=\vec M\ \Rightarrow\ \Gamma_{MF}=N_0\lim_{n\to0}\langle M_1^2\rangle", 36),
                   -2.4, 6.6, 1.85, 2.95)
    eq2 = place_in(tx(r"1=\int d\lambda\,\delta(\lambda-|\vec M|^2)", 34), -6.4, -0.6, 0.45, 1.45)
    eq3 = place_in(mixed(r"$d^nM\propto\lambda^{n/2-1}d\lambda$", "：", "$n$", " 維積分變成一維", font_size=28),
                   0.0, 6.6, 0.45, 1.45)
    head = ["n"] + [f"{n:g}" for n in NS] + [r"n\to0"]
    vals = [r"N_0\langle M_1^2\rangle_n"] + [f4(g) for g in nb["gn"]] + [f4(nb["g0"])]
    tab = D.tex_table([head, vals], [2.6] + [1.45] * 6, row_h=0.6, font_size=30).shift(np.array([-5.65, -0.1, 0]))
    lim = SurroundingRectangle(VGroup(D.table_cell(tab, 0, 6, 7), D.table_cell(tab, 1, 6, 7)), buff=0.08)
    lim.set_stroke(GOOD, 2.5)
    par = mixed(f"$N_0={WK.MF_N0}$", "、", f"$r={WK.MF_R:g}$", "、", f"$u={WK.MF_U:g}$", font_size=28)
    par.move_to([0, -1.9, 0])
    assert_in_frame(row, eq1, eq2, eq3, tab, lim, par)
    assert_table(tab, 7)
    for a, b in ((row, eq1), (row, eq2), (eq1, eq3), (eq1, eq2), (eq2, eq3), (eq2, tab), (eq3, tab), (tab, par)):
        assert_apart(a, b)
    ns_txt = ",".join(f"{n:g}" for n in NS)

    s.play(FadeIn(row, lag_ratio=0.1), run_time=1.0)
    s.say("平均場：所有自旋都一樣，", r"$\vec S_i=\vec M$", wait=2.0)
    s.play(Write(eq1), run_time=1.5)
    s.say("只剩一個 ", "$n$", " 維向量 ", r"$\vec M$", " 要積分", wait=2.0)
    s.play(FadeIn(eq2), run_time=0.8)
    s.say("插入 ", r"$1=\int d\lambda\,\delta(\lambda-|\vec M|^2)$", wait=2.0)
    s.play(FadeIn(eq3), run_time=0.8)
    s.say("$n$", " 維積分變成對 ", r"$\lambda=|\vec M|^2$", " 的一維積分", wait=2.0)
    s.play(FadeIn(tab), FadeIn(par), run_time=1.0)
    s.say("有限的 ", "$n$", " 也能算：", f"$n={ns_txt}$", wait=2.0)
    s.play(ShowCreation(lim), run_time=0.8)
    s.say(r"$n\to0$", " 收斂到 ", f"${f4(nb['g0'])}$", wait=2.0)
    s.clear_stage()


def part_s_integral(s):
    nb = numbers()
    tag = zh("(10.13)", font_size=22, color=GREY_B).move_to([-5.95, 2.45, 0])
    top = place_left(tx(r"\Gamma_{MF}=\frac{N_0}{2\pi}\iint d\lambda\,dS\,f(\lambda)\,\frac{e^{iS\lambda}}{2iS}", 36),
                     -5.0, 2.45, 6.6)
    fdef = place_left(mixed(r"$f(\lambda)=e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}$", "；註 2：",
                            r"$-iS\to-iS-\epsilon$", "，極點在實軸上方", font_size=26), -6.4, 1.5, 6.6)
    sp = s_plane(np.array([-4.0, -0.35, 0]))
    cases = place_in(tx(r"\int dS\,\frac{e^{iS\lambda}}{S-i\epsilon}=\begin{cases}2\pi i\,e^{-\epsilon\lambda}"
                        r"&\lambda>0\\0&\lambda<0\end{cases}", 34), -0.6, 6.6, -0.05, 1.05)
    sv = f4(nb["s_pos"].imag)
    s0 = f"{abs(nb['s_neg']):.0f}"
    num1 = place_left(mixed(rf"$\epsilon={EPS:g}$", "，", r"$\lambda=1$", "：", rf"${sv}i=2\pi i\,e^{{-{EPS:g}}}$",
                            font_size=26), -0.6, -0.55, 6.6)
    num2 = place_left(mixed(rf"$\epsilon={EPS:g}$", "，", r"$\lambda=-1$", "：", f"${s0}$", font_size=26),
                      -0.6, -1.15, 6.6)
    res = place_in(tx(r"\Rightarrow\ \Gamma_{MF}=\frac{N_0}2\int_0^\infty e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}"
                      r"d\lambda", 36), -6.4, 6.6, -2.9, -1.85)
    assert_in_frame(tag, top, fdef, sp, cases, num1, num2, res)
    for a, b in ((tag, top), (top, fdef), (fdef, sp), (fdef, cases), (sp, cases), (sp, num1), (sp, num2),
                 (cases, num1), (num1, num2), (sp, res), (num2, res)):
        assert_apart(a, b)
    res_top = res.copy().move_to([0, 2.3, 0])
    err = D.erratum_box("p. 238，(10.15)", r"\Gamma_{MF}=\frac{N_0}{4}\int_0^\infty d\lambda\,f(\lambda)",
                        r"\Gamma_{MF}=\frac{N_0}{2}\int_0^\infty d\lambda\,f(\lambda)", width=7.5)
    place_in(err, -3.9, 3.9, -0.4, 1.6)
    gns, gmf, gbook = nb["gauss"]
    gl = mixed("$u=0$", "（高斯）：每個 ", "$n$", " 都得 ", rf"$1/r={gmf:.0f}$", "，書上的式子得 ", f"${gbook:.0f}$",
               font_size=28).move_to([0, -1.2, 0])
    ll = mixed("表的極限 ", f"${f4(nb['g0'])}$", "，書上 ", f"${f4(nb['book'])}$", font_size=28).move_to([0, -1.9, 0])
    assert_in_frame(res_top, err, gl, ll)
    for a, b in ((res_top, err), (err, gl), (gl, ll)):
        assert_apart(a, b)

    s.play(FadeIn(tag), Write(top), run_time=1.5)
    s.say("(10.13)：把 ", r"$\delta$", " 函數寫成對 ", "$S$", " 的積分", wait=2.0)
    s.play(FadeIn(fdef), run_time=0.8)
    s.play(GrowArrow(sp[0]), FadeIn(sp[1]), FadeIn(sp[4]), FadeIn(sp[5]), run_time=1.0)
    s.say("註 2：", r"$-iS\to-iS-\epsilon$", "，極點 ", r"$S=i\epsilon$", " 在實軸上方", wait=2.0)
    s.play(ShowCreation(sp[2]), FadeIn(sp[6]), ShowCreation(sp[3]), FadeIn(sp[7]), FadeIn(cases), run_time=1.2)
    s.say(r"$\lambda>0$", " 往上包住極點：", r"$2\pi i$", "；", r"$\lambda<0$", " 往下：", "$0$", wait=2.0)
    s.play(FadeIn(num1), FadeIn(num2), run_time=0.8)
    s.say("數值積分：", r"$\lambda=1$", " 得 ", rf"${sv}i=2\pi i\,e^{{-{EPS:g}}}$", "，", r"$\lambda=-1$",
          f" 得 {s0}", wait=2.0)
    s.play(Write(res), run_time=1.5)
    s.say("只剩 ", r"$\lambda>0$", "：", r"$\Gamma_{MF}=\frac{N_0}2\int_0^\infty f\,d\lambda$", wait=2.0)
    s.play(*[FadeOut(m) for m in (tag, top, fdef, sp, cases, num1, num2)],
           Transform(res, res_top), run_time=1.0)
    s.play(FadeIn(err), run_time=0.9)
    s.say("書上說對 ", "$S$", " 的積分是 ", r"$i\pi$", "（那是主值），應為 ", r"$2\pi i$", wait=2.0)
    s.play(FadeIn(gl), run_time=0.8)
    s.say("$u=0$", "（高斯）時每個 ", "$n$", " 都得 ", rf"$1/r={gmf:.0f}$", f"，書上的式子得 {gbook:.0f}", wait=2.0)
    s.play(FadeIn(ll), run_time=0.8)
    s.say("表的極限 ", f"${f4(nb['g0'])}$", "，書上 ", f"${f4(nb['book'])}$", wait=2.0)
    s.clear_stage()


def part_asymptotics(s):
    nb = numbers()
    ax = D.plot_axes((0, 1, 0.25), (0, 1.2, 0.5), width=4.8, height=3.6, x_ticks=(0, 0.5, 1),
                     y_ticks=(0.5, 1), x_tick_labels=["0", "0.5", "1"], y_tick_labels=["0.5", "1"],
                     x_label=r"\lambda")
    ax.move_to(np.array([-3.85, -0.1, 0]))
    ylab = mixed(r"$f(\lambda)$", "（各自正規化到峰值 1）", font_size=22).move_to([-3.85, 2.45, 0])
    cp = f_curve(ax, R_POS, POS_COLOR)
    cn = f_curve(ax, R_NEG, NEG_COLOR)
    lp = tx(rf"r={R_POS:g}", 26, POS_COLOR).move_to(ax.c2p(0.2, 0.92))
    lq = tx(rf"r={R_NEG:g}", 26, NEG_COLOR).move_to(ax.c2p(0.84, 0.92))
    lam_s = abs(R_NEG) / WK.MF_U
    vline = DashedLine(ax.c2p(lam_s, 0), ax.c2p(lam_s, 1.0), dash_length=0.08).set_stroke(GREY_B, 2)
    lstar = tx(rf"\lambda^*=|r|/u={lam_s:g}", 22, GREY_A).next_to(ax.c2p(lam_s, 1.0), UP, buff=0.12)
    widths = [1.9, 1.15, 1.15, 1.15]
    tw = sum(widths)
    ta_title = place_left(mixed("$r>0$", "：", r"$\Gamma_{MF}\approx1/r$", "（普通漫步 ", r"$\textstyle\sum_N e^{-rN}$",
                                "）", font_size=26), X0, 2.82, 6.7)
    ta = D.tex_table([[r"N_0", *N0_TEX], [r"\Gamma_{MF}", *[f4(v) for v in nb["pos"]]]], widths,
                     row_h=0.48, font_size=28).shift(np.array([X0, 2.4, 0]))
    ta_lim = place_left(tx(rf"\to1/r={1 / R_POS:g}", 28, GOOD), X0 + tw + 0.2, 2.4 - 1.5 * 0.48, 6.75)
    tb_title = place_left(mixed("$r<0$", "：", r"$\Gamma_{MF}\approx\sqrt{\pi N_0/u}\,e^{N_0r^2/4u}$", font_size=26),
                          X0, 1.0, 6.7)
    tb = D.tex_table([[r"N_0", *N0_TEX], [zh("精確除以漸近", font_size=22), *[f4(v) for v in nb["ratio"]]],
                      [r"\ln\Gamma_{MF}/N_0", *[f4(v) for v in nb["lng"]]]], widths, row_h=0.48,
                     font_size=28).shift(np.array([X0, 0.62, 0]))
    tb_lim = place_left(tx(rf"\to r^2/4u={f4(R_NEG ** 2 / (4 * WK.MF_U))}", 26, GOOD), X0 + tw + 0.15,
                        0.62 - 2.5 * 0.48, 6.75)
    o_title = place_left(mixed("極限的順序（", rf"$r={R_NEG:g}$", "，", f"$N_0={WK.MF_N0}$", "）", font_size=28),
                         X0, 2.6, 6.7)
    ot = D.tex_table([["n", *NS_NEG_TEX], [r"\Gamma_n/\Gamma_{MF}", *[f4(v) for v in nb["order"]]]],
                     [1.9, 1.2, 1.2, 1.2, 1.2], row_h=0.55, font_size=28).shift(np.array([X0, 2.2, 0]))
    o1 = place_left(mixed("要 ", r"$n\,e^{N_0r^2/4u}\ll1$", " 才到極限", font_size=28), X0, 0.45, 6.7)
    o2 = place_left(mixed("必須先取 ", r"$n\to0$", "，再讓盒子變大", font_size=28, color=YELLOW), X0, -0.3, 6.7)
    plot = VGroup(ax, ylab, cp, cn, lp, lq, vline, lstar)
    assert_in_frame(plot, ta_title, ta, ta_lim, tb_title, tb, tb_lim, o_title, ot, o1, o2)
    assert_table(ta, 4)
    assert_table(tb, 4)
    assert_table(ot, 5)
    for m in (ta_title, ta, ta_lim, tb_title, tb, tb_lim, o_title, ot, o1, o2):
        assert_apart(plot, m)
    for a, b in ((ta_title, ta), (ta, ta_lim), (ta, tb_title), (ta_lim, tb_title), (tb_title, tb), (tb, tb_lim),
                 (o_title, ot), (ot, o1), (o1, o2), (lp, lstar), (lq, lstar), (lp, lq)):
        assert_apart(a, b)

    s.play(FadeIn(ax), FadeIn(ylab), run_time=0.8)
    s.say("兩種情況：", r"$f(\lambda)$", " 在 ", "$r>0$", " 和 ", "$r<0$", " 長得很不一樣", wait=2.0)
    s.play(ShowCreation(cp), FadeIn(lp), run_time=1.2)
    s.say("$r>0$", "：", "$f$", " 從 ", r"$\lambda=0$", " 一路遞減", wait=2.0)
    s.play(FadeIn(ta_title), FadeIn(ta), FadeIn(ta_lim), run_time=1.0)
    s.say(rf"$\Gamma_{{MF}}\to1/r={1 / R_POS:g}$", "：就是普通漫步 ", r"$\sum_N e^{-rN}$", wait=2.0)
    s.play(ShowCreation(cn), FadeIn(lq), ShowCreation(vline), FadeIn(lstar), run_time=1.2)
    s.say("$r<0$", "：峰在 ", rf"$\lambda^*=|r|/u={lam_s:g}$", wait=2.0)
    s.play(FadeIn(tb_title), FadeIn(tb), FadeIn(tb_lim), run_time=1.0)
    s.say("鞍點近似越來越準：", r"$\ln\Gamma_{MF}/N_0\to r^2/4u$", wait=2.0)
    s.play(*[FadeOut(m) for m in (ta_title, ta, ta_lim, tb_title, tb, tb_lim)], run_time=0.6)
    s.play(FadeIn(o_title), FadeIn(ot), run_time=1.0)
    s.say("有限的 ", "$n$", "：", "$r<0$", " 時離極限很遠", wait=2.0)
    s.play(FadeIn(o1), run_time=0.6)
    s.say("要 ", r"$n\,e^{N_0r^2/4u}\ll1$", " 才到極限", wait=2.0)
    s.play(FadeIn(o2), run_time=0.6)
    s.say("所以必須先取 ", r"$n\to0$", "，再讓盒子變大", wait=2.0)
    s.clear_stage()


def part_footnote(s):
    nb = numbers()
    g0, gb = f4(nb["fn3"][0]), f4(nb["fn3"][1])
    b1, b2 = f4(nb["fn3"][3]), f4(nb["fn3"][4])
    tag_a = zh("(10.15)", font_size=22, color=GREY_B).move_to([-6.05, 2.45, 0])
    forms = [tx(F_TEX[0], 36)] + [Tex(F_TEX[k], font_size=36, t2c={F_KEYS[k - 1]: SUB_COLORS[k - 1]})
                                  for k in range(1, 5)]
    for f in forms:
        place_left(f, -5.2, 2.45, 6.7)
    tag_b = zh("(7.57)", font_size=22, color=GREY_B).move_to([-6.05, 1.3, 0])
    g757 = place_left(tx(r"G_{MF}=\int_0^\infty e^{N(-f\lambda^2-\lambda(1-x))}d\lambda", 36), -5.2, 1.3, 6.7)
    rules = VGroup(*[tx(t, 34, c) for t, c in zip(RULES, SUB_COLORS)]).arrange(RIGHT, buff=0.7)
    place_in(rules, -6.4, 6.6, -0.05, 0.6)
    res = mixed(f"$r={R_POS:g}$", "（", r"$x=z/z_c=\frac12$", "，第 7 章 Ex. 7.5 的 ", "$x$", "）：",
                rf"$\Gamma_{{MF}}={g0}={WK.MF_N0}\,G_{{MF}}(\tfrac12)$", font_size=26)
    place_in(res, -6.5, 6.6, -1.05, -0.4)
    err = D.erratum_box("p. 238，footnote 3", r"\lambda\to4\lambda,\ 4u\to f,\ 2r\to1-z/z_c",
                        r"\lambda\to2\lambda,\ u\to f,\ r\to1-z/z_c", width=9.0)
    place_in(err, -4.6, 4.6, -2.9, -1.3)
    assert g0 == gb
    assert_in_frame(tag_a, tag_b, g757, rules, res, err, *forms)
    for f in forms:
        assert_apart(tag_a, f)
        assert_apart(f, g757)
    for a, b in ((tag_b, g757), (g757, rules), (rules, res), (res, err)):
        assert_apart(a, b)
    caps = (("第一步：", r"$N_0\to N$"),
            ("第二步：", r"$\lambda\to2\lambda$", "（", r"$d\lambda$", " 也多一個 2）"),
            ("第三步：", r"$u\to f$"),
            ("第四步：", r"$r\to1-z/z_c$", "（", r"$x=z/z_c$", "）"))

    s.play(FadeIn(tag_a), FadeIn(forms[0]), run_time=1.0)
    s.play(FadeIn(tag_b), FadeIn(g757), run_time=1.0)
    s.say("註 3：(10.15) 和第 7 章的 (7.57) 是同一個積分", wait=2.0)
    cur = forms[0]
    for k in range(4):
        s.play(FadeIn(rules[k]), run_time=0.5)
        s.play(ReplacementTransform(cur, forms[k + 1]), run_time=1.2)
        cur = forms[k + 1]
        s.say(*caps[k], wait=1.5)
    s.say(r"$\Gamma_{MF}=N\,G_{MF}$", "：對終點加總，多乘一個 ", "$N$", wait=2.0)
    s.play(FadeIn(res), run_time=0.8)
    s.say(f"$r={R_POS:g}$", " 就是 ", r"$x=\frac12$", "：", rf"${g0}={WK.MF_N0}\,G_{{MF}}(\tfrac12)$", wait=2.0)
    s.play(FadeIn(err), run_time=0.9)
    s.say("書上的換法只對它自己的 ", "$N_0/4$", " 成立", wait=2.0)
    s.say("兩個 2 倍抵掉：", f"${b1}={b2}$", wait=2.0)
    s.say(r"$r\to1-z/z_c$", " 才和 (10.4) 的 ", r"$e^{-r}=z/z_c$", " 一致", wait=2.0)
    s.say("平均場把關在盒子裡的自迴避漫步變成一維積分", wait=2.0)
    s.say("$r<0$", " 多出的 ", r"$e^{N_0r^2/4u}$", " 就是第 7 章的那個因子", wait=2.5)


class W1002MeanField(ChuScene):
    def construct(self):
        self.show_header("W1002　O(n) 模型的平均場")
        check_math()
        part_model(self)
        part_mean_field(self)
        part_s_integral(self)
        part_asymptotics(self)
        part_footnote(self)
        self.wait(1.0)
