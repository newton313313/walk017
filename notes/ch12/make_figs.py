r"""
make_figs.py —— 隨機漫步導讀第 12 章（Rudnick & Gaspari《Elements of the Random Walk》第 12 章，動畫 W1201–W1206）的圖與表

在 notes/ 或 notes/ch12/ 裡執行：
    python ch12/make_figs.py
會先檢查 ch12/part12.tex 寫到的每一個數字（check_numbers_ch12；最後再把 part12.tex 裡出現的每一個小數和檢查過的數比對一遍），
再寫出 ch12/figs.tex（part12.tex 第一行 \input{ch12/figs}）。圖都用 TikZ 畫，座標由這裡算好（不用 pgfplots、不在 TeX 裡算）。
數學一律由 codes/utils/walk_eps.py 算出（全書回顧的數字用 walk_eps.recap_numbers()）；例子與動畫 W1201–W1206 相同。巨集名稱都以 \Twelve 開頭
（不和前幾章的 \Fig…、\Two…、…、\Eleven… 撞名）；part12.tex 的 \label 用 twelve-…。顏色沿用導言區：cBar 精確值、cWalk 路線、
cExact 更正、cGauss 近似、cWrong 筆誤。全部約 2 秒（r_critical 由 walk_eps 快取）。
"""
from __future__ import annotations

import math
import re
import sys
import warnings
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent            # notes/ch12
sys.path.insert(0, str(HERE.parent.parent / "codes"))
warnings.filterwarnings("ignore")                  # scipy 的 IntegrationWarning 不影響結果

from utils import walk_eps as E    # noqa: E402   （套件模組：walk_eps 的相對 import 才找得到 walk_saw、walk_crit、walk_rg）

OUT = HERE / "figs.tex"
TEX = HERE / "part12.tex"
F = Fraction
OK = set()                                         # 檢查過的小數字串（part12.tex 的小數都要在這裡）


def _s(x, d):
    """小數 d 位的字串（和講義上寫的一樣），並記進 OK。"""
    s = f"{float(x):.{d}f}"
    OK.add(s.lstrip("-"))
    return s


def _fr(x):
    x = F(x)
    if x.denominator == 1:
        return str(x.numerator)
    sgn = "-" if x < 0 else ""
    return rf"{sgn}\frac{{{abs(x.numerator)}}}{{{x.denominator}}}"


def _close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) < tol


def _all(vals, d, expect):
    got = tuple(_s(v, d) for v in vals)
    assert got == tuple(expect), (got, expect)


# =============================================================================
# 數字檢查（part12.tex 寫到的每一個數）
# =============================================================================
def check_numbers_ch12():
    OK.update({"0.2", "0.01", "0.1", "0.5", "0.05", "0.005", "0.06", "0.001", "0.588", "0.630", "0.672", "0.711",
               "0.61", "0.563", "0.583", "0.6", "0.600", "0.02", "0.03", "0.45", "0.04", "0.002", "1.0", "2.5", "0.25"})
    # ---- W1201 殼層、兩變數玩具、重新標度
    r, u = E.TOY
    assert (r, u) == (0.2, 0.01)
    rp, up, c6 = E.toy_coeffs(r, u)
    assert (_s(rp, 4), _s(up, 4), _s(c6, 6)) == ("0.3000", "0.0075", "0.000167")
    assert _close(rp, r + 12 * u / (1 + r)) and _close(up, u - 36 * u ** 2 / (1 + r) ** 2) and _close(c6, (12 * u / 1.2) ** 3 / 6)
    assert all(_close(a, b, 1e-5) for a, b in zip(E.toy_coeffs(r, u), E.toy_coeffs_numeric(r, u)))
    assert _close(E.toy_Heff(1.0, r, u), E.toy_Heff_numeric(1.0, r, u), 1e-8) and _s(E.toy_Heff(1.0, r, u), 4) == "0.1577"
    Z2, Zf = E.toy_Z(r, u)
    A = E.toy_A(r)
    assert (_s(Z2, 4), _s(Zf, 4), _s(A, 4), _s(Zf / A, 4)) == ("9.0517", "9.0517", "2.2882", "3.9558")
    assert _close(A, math.sqrt(2 * math.pi / 1.2))
    sh = [(E.shell_fraction(*c), E.shell_linear(*c)) for c in E.SHELL_CASES]
    assert [(_s(a, 4), _s(b, 4)) for a, b in sh] == [("0.2592", "0.3000"), ("0.0296", "0.0300"), ("0.0392", "0.0400")]
    lat = [E.shell_fraction_lattice(3, dl) for dl in (0.1, 0.01)]
    assert abs(lat[0] - sh[0][0]) < 2e-4 and abs(lat[1] - sh[1][0]) < 2e-4
    assert (_s(lat[0], 4), _s(lat[1], 4)) == ("0.2591", "0.0296")
    dq, lin = E.shell_dq(0.5, 0.01)
    assert (_s(dq, 6), _s(lin, 6), _s(math.exp(-0.5), 4)) == ("0.006035", "0.006065", "0.6065")
    assert _all(E.rescale_factors(3, 0.1), 4, ["1.0000", "1.2214", "1.1052"]) is None
    assert _s(E.rescale_factors(3, 0.1, book=True)[0], 4) == "0.6703" and _close(E.rescale_factors(3, 0.1, book=True)[0], math.exp(-0.4))
    assert (_s(math.exp(0.03), 4), _s(1.01, 4)) == ("1.0305", "1.0100")
    assert [math.comb(4, k) for k in range(5)] == [1, 4, 6, 4, 1]
    # ---- W1202 配對、組合因子、連結叢集
    assert [len(E.matchings(n)) for n in (3, 4, 5, 6)] == [15, 105, 945, 10395]
    assert E.type_census(3) == {(1, 1, 1): 1, (2, 1): 6, (3,): 8} and E.diagram_census(3) == {(1, 1, 1): 1, (2, 1): 3, (3,): 1}
    assert E.type_census(4) == {(1, 1, 1, 1): 1, (2, 1, 1): 12, (3, 1): 32, (2, 2): 12, (4,): 48}
    assert E.diagram_census(4) == {(1, 1, 1, 1): 1, (2, 1, 1): 6, (3, 1): 4, (2, 2): 3, (4,): 3}
    assert E.type_census(5)[(2, 1, 1, 1)] == 20 and E.diagram_census(5)[(2, 1, 1, 1)] == 10 == E.pair_factor(3, 1)
    assert math.factorial(5) // (math.factorial(2) * math.factorial(1)) == 60                    # 書上第二次寫的 5!/(2!1!)
    for n in (3, 4, 5):
        assert sum(E.type_census(n).values()) == len(E.matchings(n))
    pf = {(2, 1): 6, (0, 2): 3, (3, 1): 10, (4, 1): 15, (2, 2): 45, (0, 3): 15}
    for (n1, n2), v in pf.items():
        assert E.pair_factor(n1, n2) == v == E.diagram_census(n1 + 2 * n2)[E.sizes_of(n1, n2)]
    assert [E.ring_aut(k) for k in range(1, 6)] == [1, 2, 6, 8, 10] and all(E.ring_aut(k) == 2 * k for k in (4, 5))
    assert [E.ring_weight(k) for k in range(1, 5)] == [1, 2, 8, 48]
    assert all(E.type_census(k)[(k,)] == E.ring_weight(k) for k in range(1, 6))
    assert all(E.diagram_census(k)[(k,)] == math.factorial(k) // E.ring_aut(k) for k in range(1, 6))
    assert E.general_factor({1: 3, 2: 1}) == 10 and E.general_factor({3: 1}) == 1
    assert 6 * 6 * 2 == 72 == 8 * (1 + 8)
    _all(E.egf_check(0.1), 4, ["1.1180"] * 3)
    _all([E.lc_average(*E.LC, m) for m in ("closed", "quad", "series", "linked")], 4, ["0.9129"] * 4)
    assert _s(math.log(E.lc_average(*E.LC)), 4) == "-0.0912"
    _all(E.lc_order(*E.LC, 2), 4, ["0.0150", "0.0100", "0.0050"])
    _all(E.lc_cumulants_rings(*E.LC, 4), 4, ["-0.1000", "0.0200", "-0.0080", "0.0048"])
    assert all(_close(a, b) for a, b in zip(E.lc_cumulants_rings(*E.LC, 4), E.lc_cumulants_moments(*E.LC, 4)))
    _all([E.lc_multi(method=m) for m in ("det", "quad", "linked")], 4, ["0.8087"] * 3)
    # ---- W1203 一迴圈與流
    assert [E.AB(n) for n in range(4)] == [(8, 32), (12, 36), (16, 40), (20, 44)]
    assert [tuple(round(x) for x in E.on_counts(n)) for n in range(1, 5)] == [(6, 36), (8, 40), (10, 44), (12, 48)]
    assert [2 * n + 4 for n in range(4)] == [4, 6, 8, 10]
    fps = {e: E.fixed_point(e) for e in E.NU_EPS}
    assert fps[F(1)] == (F(-1, 7), F(1, 49)) and fps[F(1, 2)] == (F(-1, 13), F(2, 169)) and fps[F(1, 10)] == (F(-1, 61), F(10, 3721))
    _all([x for e in reversed(E.NU_EPS) for x in fps[e]], 4, ["-0.1429", "0.0204", "-0.0769", "0.0118", "-0.0164", "0.0027"])
    _all([x for e in reversed(E.NU_EPS) for x in E.fixed_point_lowest(e)], 4,
         ["-0.1667", "0.0278", "-0.0833", "0.0139", "-0.0167", "0.0028"])
    assert all(_close(a, b, 1e-8) for e in E.NU_EPS for a, b in zip(fps[e], E.fixed_point_numeric(float(e))))
    assert all(max(abs(x) for x in E.beta(float(fps[e][0]), float(fps[e][1]), float(e))) < 1e-12 for e in E.NU_EPS)
    _all(E.eigen(1.0, 1, "gauss"), 4, ["2.0000", "1.0000"])
    for s0 in E.flow_starts():
        L = E.flow_line(*s0)
        assert len(L) == 40
    mids = [E.flow_line(*E.flow_starts()[k]) for k in (1, 4)]
    assert all(abs(m[-1][0] + 1 / 7) < 2e-3 and abs(m[-1][1] - 1 / 49) < 2e-4 for m in mids)
    # ---- W1204 ε 展開
    _all([E.u_analytic(l, 0.005, 1.0) for l in (1, 2, 5)], 4, ["0.0104", "0.0172", "0.0270"])
    _all([E.u_analytic(l, 0.06, 1.0) for l in (1, 2, 5)], 4, ["0.0346", "0.0300", "0.0279"])
    assert all(_close(E.u_analytic(l, u0, 1.0), E.u_numeric(l, u0, 1.0), 1e-7) for l in (1, 2, 5) for u0 in (0.005, 0.06))
    assert _s(1 / 36, 4) == "0.0278" and all(_close(E.u_from_1240(0.005, E.u_analytic(l, 0.005, 1.0), 1.0), l, 1e-7) for l in (1, 2))
    assert (_s(E.u_analytic(2, 1e-6, 1.0) / 1e-6, 4), _s(math.exp(2), 4)) == ("7.3874", "7.3891")
    assert (_s(E.int_1248(0.001, 0.005, 1.0, 3), 4), _s(E.int_1249(0.001, 0.005, 1.0, 3), 4),
            _s(E.int_1248(0.001, 0.005, 1.0, 3, book=True), 4)) == ("0.2207", "0.2207", "81.7804")
    assert (_s(E.t_closed(3, 1e-4, 0.01, 1.0), 4), _s(E.t_numeric(3, 1e-4, 0.01, 1.0), 4)) == ("0.0203", "0.0203")
    assert _close(E.t_closed(3, 1e-4, 0.01, 1.0), E.t_numeric(3, 1e-4, 0.01, 1.0), 1e-9)
    nus = [(E.nu_first(e), E.nu_unexpanded(e), E.nu_eigen(float(e))) for e in E.NU_EPS]
    _all([x[0] for x in nus], 4, ["0.5083", "0.5417", "0.5833"])
    _all([x[1] for x in nus], 4, ["0.5085", "0.5455", "0.6000"])
    _all([x[2] for x in nus], 4, ["0.5076", "0.5256", "0.5272"])
    assert E.nu_first(1) == F(7, 12) and E.nu_unexpanded(1) == F(3, 5)
    assert all(E.nu_unexpanded(e) == F(3) / (6 - e) == E.flory_nu(4 - e) for e in E.NU_EPS)
    assert _close(E.nu_eigen(1.0), E.nu_from_flow(1.0), 1e-4) and _close(1 / E.lambda_r_closed(1.0), E.nu_eigen(1.0))
    J = E.jacobian(*fps[F(1)], F(1))
    assert [[F(x).limit_denominator(100) for x in row] for row in np.asarray(J, dtype=float)] == [[F(5, 3), 14], [F(1, 21), -1]]
    _all(E.eigen(1.0), 4, ["1.8968", "-1.2301"])
    assert _close(E.eigen(1.0)[0], (1 + math.sqrt(22)) / 3) and _close(E.eigen(1.0)[1], (1 - math.sqrt(22)) / 3)
    _all([E.r_critical(1e-4, 0.1), E.r_critical_linear(1e-4, 0.1), -6e-4], 6, ["-0.000631", "-0.000632", "-0.000600"])
    _all([E.r_critical(1e-3, 1.0), E.r_critical_linear(1e-3, 1.0)], 4, ["-0.0111", "-0.0120"])
    # ---- W1205 O(n)
    _all([E.nu_first(1, n) for n in range(4)], 4, ["0.5625", "0.5833", "0.6000", "0.6136"])
    assert [E.nu_first(1, n) for n in range(4)] == [F(9, 16), F(7, 12), F(3, 5), F(27, 44)]
    _all([E.nu_unexpanded(1, n) for n in range(4)], 4, ["0.5714", "0.6000", "0.6250", "0.6471"])
    _all([E.nu_eigen(1.0, n) for n in range(4)], 4, ["0.5205", "0.5272", "0.5325", "0.5368"])
    _all([E.NU_3D_CITED[n] for n in range(4)], 3, ["0.588", "0.630", "0.672", "0.711"])
    assert E.nu_first(F(1), 0) == F(1, 2) + F(1, 16) and E.fixed_point(F(1), 0) == (F(-1, 9), F(2, 81))
    assert [E.flory_series(k) for k in range(4)] == [F(1, 2), F(1, 12), F(1, 72), F(1, 432)] and E.flory_nu(3) == F(3, 5)
    _all([E.saw_size(1000, nu) for nu in (0.5, F(9, 16), 0.588)], 2, ["31.62", "48.70", "58.08"])
    # ---- W1206 圖解法
    assert (_s(math.exp(-5), 4), _s(E.g_numeric(math.exp(-5), 0.005, 1.0), 4)) == ("0.0067", "0.0270")
    assert _close(E.g_numeric(math.exp(-5), 0.005, 1.0), E.u_analytic(5, 0.005, 1.0), 1e-7)
    assert [(E.pc_exponent(d), E.pc_exponent(d, True), E.pc_rhs(d), E.pc_rhs(d, True)) for d in (2, F(5, 2), 3, 4)] == \
        [(-1, 1, -1, -3), (F(-1, 2), 2, F(-1, 2), F(-5, 2)), (0, 3, 0, -2), (1, 5, 1, -1)]
    _all([E.pc_numeric(2), E.pc_numeric(2.5)], 4, ["-1.0000", "-0.5000"])
    # ---- 全書回顧
    rec = E.recap_numbers()
    assert rec == E.RECAP_EXPECT
    for ch, vals in rec.items():
        for v in vals:
            if "." in v:
                OK.add(v)
    assert sum(E.BOOK_SCENES.values()) == 87 and len(E.BOOK_RECAP) == 12
    assert [E.scene_range(c) for c in (1, 12)] == ["W101–W112", "W1201–W1206"]
    assert [math.comb(N, N // 2) ** 2 for N in (0, 2, 4, 6, 8, 10)] == [1, 4, 36, 400, 4900, 63504]
    assert _s(F(10, 19), 4) == "0.5263" and _s(252 / 1024, 4) == "0.2461"


def check_tex_decimals():
    """part12.tex（去掉註解）裡的每一個小數都必須是上面檢查過的數。"""
    if not TEX.exists():
        print("（尚無 part12.tex，略過小數比對）")
        return
    txt = re.sub(r"(?<!\\)%.*", "", TEX.read_text(encoding="utf-8"))
    txt = re.sub(r"[\d.]+\\linewidth|[\d.]+pt\b|[\d.]+em\b|[\d.]+ex\b|[\d.]+cm\b|\\arraystretch\}\{[\d.]+\}", "", txt)   # 版面尺寸
    eqno = re.compile(r"12\.\d{1,3}")                    # 本章的式號、節號、圖號、練習號（12.1–12.89）
    bad = sorted({m for m in re.findall(r"\d+\.\d+", txt) if m not in OK and not eqno.fullmatch(m)})
    assert not bad, f"part12.tex 有未檢查的小數：{bad}"


# =============================================================================
# 圖與表
# =============================================================================
ARROW = r"postaction={decorate},decoration={markings,mark=at position 0.5 with {\arrow{Stealth[length=3pt]}}}"


def _poly(pts, sx=1.0, sy=1.0, x0=0.0, y0=0.0, d=3):
    return " -- ".join(f"({x0 + sx * x:.{d}f},{y0 + sy * y:.{d}f})" for x, y in pts)


def _axes(x0, x1, y0, y1, sx, sy, xt, yt, xl, yl, fx="{:g}", fy="{:g}", ox=None):
    """座標軸：x ∈ [x0, x1]、y ∈ [y0, y1]，(x0, y0) 對到原點；ox 是 y 軸所在的 x（預設 x0）。"""
    ox = x0 if ox is None else ox
    X = lambda x: (x - x0) * sx
    Y = lambda y: (y - y0) * sy
    out = [rf"\draw[->,gray] (0,0) -- ({X(x1) + 0.25:.3f},0) node[right,black] {{{xl}}};",
           rf"\draw[->,gray] ({X(ox):.3f},0) -- ({X(ox):.3f},{Y(y1) + 0.25:.3f}) node[above,black] {{{yl}}};"]
    for t in xt:
        out.append(rf"\draw[gray] ({X(t):.3f},0) -- ++(0,-0.06) node[below,black] {{${fx.format(t)}$}};")
    for t in yt:
        out.append(rf"\draw[gray] ({X(ox):.3f},{Y(t):.3f}) -- ++(-0.06,0) node[left,black] {{${fy.format(t)}$}};")
    return out


def fig_shell():
    """Fig. 12.3–12.4 的平面示意：半徑 Λ 的圓形布里淵區、外圈殼層（橙）；格點 = 允許的 q。"""
    R, dl, h = 1.45, 0.28, 0.2
    Rin = R * math.exp(-dl)
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]",
           rf"\fill[cWalk!22,even odd rule] (0,0) circle ({R:.3f}) (0,0) circle ({Rin:.3f});",
           rf"\draw[cBar,thick] (0,0) circle ({R:.3f});",
           rf"\draw[cBar,densely dashed] (0,0) circle ({Rin:.3f});"]
    k = int(R / h)
    for i in range(-k, k + 1):
        for j in range(-k, k + 1):
            x, y = i * h, j * h
            rr = math.hypot(x, y)
            if rr <= R:
                c = "cWalk" if rr > Rin else "cBar!60"
                out.append(rf"\fill[{c}] ({x:.2f},{y:.2f}) circle (0.7pt);")
    out += [rf"\draw[->] (0,0) -- ({R * math.cos(0.6):.3f},{R * math.sin(0.6):.3f}) node[pos=0.55,above left=-1pt] {{$\Lambda=e^{{-\ell}}$}};",
            rf"\node[align=center] at (0,{-R - 0.45:.3f}) {{殼層（橙）：$e^{{-\ell-\delta\ell}}<|\vec q|<e^{{-\ell}}$}};",
            r"\end{tikzpicture}"]
    return "\n".join(out)


def tab_shell():
    rows = []
    for (d, dl), lat in zip(E.SHELL_CASES, (E.shell_fraction_lattice(3, 0.1), E.shell_fraction_lattice(3, 0.01), None)):
        a, b = E.shell_fraction(d, dl), E.shell_linear(d, dl)
        L = _s(lat, 4) if lat is not None else "--"
        rows.append(rf"{d} & {dl:g} & {_s(a, 4)} & {L} & {_s(b, 4)}\\")
    return "\n".join([r"\begin{tabular}{@{}ccrrr@{}}", r"\toprule",
                      r"$d$ & $\delta\ell$ & $1-e^{-d\,\delta\ell}$ & 數格點 & $d\,\delta\ell$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def fig_toy():
    """兩變數玩具：原來的 r s²/2 + u s⁴（灰虛線）、精確的 ℋ'(s)（藍）與截斷的 r' s²/2 + u' s⁴（橙點線）。"""
    r, u = E.TOY
    rp, up, _ = E.toy_coeffs(r, u)
    sx, sy = 0.95, 3.6
    xs = np.linspace(-2, 2, 40)
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(-2, 2, 0, 0.75, sx, sy, (-2, -1, 0, 1, 2), (0.25, 0.5, 0.75), r"$s$", "", ox=0)
    X = lambda x: (x + 2)
    pts = lambda f: [(X(s), f(s)) for s in xs]
    out.append(r"\draw[gray,densely dashed] " + _poly(pts(lambda s: r / 2 * s * s + u * s ** 4), sx, sy) + ";")
    out.append(r"\draw[cBar,thick] " + _poly(pts(lambda s: E.toy_Heff(s, r, u)), sx, sy) + ";")
    out.append(r"\draw[cWalk,thick,densely dotted] " + _poly(pts(lambda s: rp / 2 * s * s + up * s ** 4), sx, sy) + ";")
    H2 = E.toy_Heff(2.0, r, u)
    out.append(rf"\node[cBar,anchor=east] at ({X(1.95) * sx:.3f},{H2 * sy + 0.18:.3f}) {{$\mathcal H'(s)$}};")
    out.append(rf"\node[gray,anchor=west] at ({X(1.55) * sx:.3f},{0.22 * sy:.3f}) {{$\frac r2s^2+us^4$}};")
    out.append(rf"\node[cWalk,anchor=west] at ({X(-2.0) * sx:.3f},{0.92 * sy:.3f}) {{點線：$\frac{{r'}}2s^2+u's^4$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


# ---- 配對圖 -------------------------------------------------------------------
def _frag(x, y, legs=True):
    out = [rf"\fill ({x:.3f},{y:.3f}) circle (1.6pt);"]
    if legs:
        out.append(rf"\draw[gray] ({x - 0.17:.3f},{y - 0.2:.3f}) -- ({x:.3f},{y:.3f}) -- ({x + 0.17:.3f},{y - 0.2:.3f});")
    return out


def _loop(x, y, rad=0.16, up=True):
    c = y + rad if up else y - rad
    return [rf"\draw[cWalk,thick] ({x:.3f},{c:.3f}) circle ({rad:.3f});"]


def fig_pairs():
    """上：type 1／2／3（黑點＝四次片段，灰＝兩隻內區的手，橙＝殼層的手牽起來）。下：Ex. 12.9 的 10 張片段圖。"""
    out = [r"\begin{tikzpicture}[font=\scriptsize,baseline=(current bounding box.center)]"]
    # type 1
    out += _frag(0, 0) + _loop(0, 0)
    out.append(r"\node at (0,-0.55) {type 1};")
    # type 2
    out += _frag(1.3, 0) + _frag(2.1, 0)
    out.append(r"\draw[cWalk,thick] (1.3,0) to[bend left=45] (2.1,0);")
    out.append(r"\draw[cWalk,thick] (1.3,0) to[bend right=45] (2.1,0);")
    out.append(r"\node at (1.7,-0.55) {type 2};")
    # type 3
    tri = [(3.2, 0.0), (4.0, 0.0), (3.6, 0.6)]
    for x, y in tri:
        out += _frag(x, y)
    out.append(r"\draw[cWalk,thick] (3.2,0) -- (4.0,0) -- (3.6,0.6) -- cycle;")
    out.append(r"\node at (3.6,-0.55) {type 3};")
    # Ex. 12.9：5 個片段 → 3 個 type 1 ＋ 1 個 type 2
    pairs = list(combinations(range(5), 2))
    assert len(pairs) == E.pair_factor(3, 1) == 10
    w, gap, dx = 1.12, 0.3, 0.26
    y0 = -1.45
    for idx, (i, j) in enumerate(pairs):
        row, col = divmod(idx, 5)
        x0 = 5.0 + col * (w + gap)
        y = 0.55 - row * 1.0
        out.append(rf"\draw[gray!40,rounded corners=2pt] ({x0 - 0.12:.3f},{y - 0.32:.3f}) rectangle ({x0 + 4 * dx + 0.12:.3f},{y + 0.4:.3f});")
        for k in range(5):
            x = x0 + k * dx
            out += _frag(x, y, legs=False)
            if k not in (i, j):
                out += _loop(x, y, 0.09, up=False)
        xi, xj = x0 + i * dx, x0 + j * dx
        out.append(rf"\draw[cWalk,thick] ({xi:.3f},{y:.3f}) to[bend left=50] ({xj:.3f},{y:.3f});")
        out.append(rf"\draw[cWalk,thick] ({xi:.3f},{y:.3f}) to[bend left=20] ({xj:.3f},{y:.3f});")
    out.append(rf"\node at ({5.0 + 2 * (w + gap) + 0.6:.3f},{y0 + 0.4:.3f}) {{Ex.~12.9：$\frac{{5!}}{{3!\,2!}}={E.pair_factor(3, 1)}$ 張片段圖}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def _sz(t):
    return "(" + ",".join(str(k) for k in t) + ")"


def tab_census():
    t3, d3 = E.type_census(3), E.diagram_census(3)
    t4, d4 = E.type_census(4), E.diagram_census(4)
    k3, k4 = list(t3), list(t4)
    rows = []
    for i in range(max(len(k3), len(k4))):
        a = rf"${_sz(k3[i])}$ & {t3[k3[i]]} & {d3[k3[i]]}" if i < len(k3) else " & & "
        b = rf"${_sz(k4[i])}$ & {t4[k4[i]]} & {d4[k4[i]]}"
        rows.append(f"{a} & {b}\\\\")
    rows.append(r"\midrule")
    rows.append(rf"共 & {sum(t3.values())} & & 共 & {sum(t4.values())} & \\")
    return "\n".join([r"\begin{tabular}{@{}lrr@{\quad}lrr@{}}", r"\toprule",
                      r"\multicolumn{3}{c}{$n=3$} & \multicolumn{3}{c}{$n=4$}\\",
                      r"環的大小 & 配對 & 圖 & 環的大小 & 配對 & 圖\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_s4():
    rows = []
    for n1, n2 in ((2, 1), (0, 2), (3, 1), (4, 1), (2, 2), (0, 3)):
        n = n1 + 2 * n2
        v = E.pair_factor(n1, n2)
        assert v == E.diagram_census(n)[E.sizes_of(n1, n2)]
        rows.append(rf"$({n1},{n2})$ & {n} & {v} & {E.diagram_census(n)[E.sizes_of(n1, n2)]}\,$\checkmark$\\")
    return "\n".join([r"\begin{tabular}{@{}lrrr@{}}", r"\toprule",
                      r"$(n_1,n_2)$ & $n$ & $\frac{n!}{n_1!\,n_2!\,2^{n_2}}$ & 列舉\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


# ---- (r,u) 流 ------------------------------------------------------------------
def fig_flow():
    """(r, u) 平面（ε = 1）：6 條流線（中間那條是臨界線，流進 WF 點），Gaussian 點與 Wilson–Fisher 點。"""
    sx, sy = 7.6, 92.0
    x0, y0 = -0.45, 0.0
    X = lambda x: (x - x0) * sx
    Y = lambda y: (y - y0) * sy
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(-0.45, 0.45, 0, 0.045, sx, sy, (-0.4, -0.2, 0, 0.2, 0.4), (0.01, 0.02, 0.03, 0.04), r"$r$", r"$u$", ox=0)
    starts = E.flow_starts()
    for k, s0 in enumerate(starts):
        L = E.flow_line(*s0)
        crit = k in (1, 4)
        sty = "cExact,thick" if crit else "cBar"
        P = [(X(a), Y(b)) for a, b in L]
        out.append(rf"\draw[{sty},{'' if crit else '->,'}{ARROW}] " + _poly(P) + ";")
        out.append(rf"\fill[{sty.split(',')[0]}] ({P[0][0]:.3f},{P[0][1]:.3f}) circle (0.9pt);")
    rs, us = E.fixed_point(F(1))
    out.append(rf"\draw[fill=white,thick] ({X(0):.3f},{Y(0):.3f}) circle (2.2pt);")
    out.append(rf"\draw[fill=white,thick] ({X(float(rs)):.3f},{Y(float(us)):.3f}) circle (2.2pt);")
    out.append(rf"\draw[gray] ({X(float(rs)) - 0.06:.3f},{Y(float(us)) + 0.06:.3f}) -- ({X(-0.36):.3f},{Y(0.046):.3f}) "
               rf"node[above,black] {{WF $\big({_fr(rs)},{_fr(us)}\big)$}};")
    out.append(rf"\node[anchor=north west] at ({X(0) + 0.05:.3f},{Y(0) - 0.25:.3f}) {{Gaussian $(0,0)$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_fixed():
    rows = []
    for e in reversed(E.NU_EPS):
        rs, us = E.fixed_point(e)
        lr, lu = E.fixed_point_lowest(e)
        rows.append(rf"${_fr(e)}$ & ${_fr(rs)}={_s(rs, 4)}$ & ${_fr(us)}={_s(us, 4)}$ & ${_s(lr, 4)}$ & ${_s(lu, 4)}$\\")
    return "\n".join([r"\begin{tabular}{@{}lllrr@{}}", r"\toprule",
                      r"$\epsilon$ & $r^*=-\frac\epsilon{6+\epsilon}$ & $u^*=\frac\epsilon{(6+\epsilon)^2}$ & "
                      r"一階 $-\frac\epsilon6$ & $\frac\epsilon{36}$\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


# ---- u(ℓ) 與 ν(ε) ---------------------------------------------------------------
def fig_u():
    """(12.41) 的 u(ℓ)（ε = 1）：u0 = 0.005 與 0.06 兩條、數值積分的點、u* = 1/36；◇ 是 (12.86) 的 g0(q0 = e^{−5})。"""
    sx, sy = 0.82, 62.0
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(0, 6, 0, 0.06, sx, sy, (0, 1, 2, 3, 4, 5, 6), (0.01, 0.02, 0.03, 0.04, 0.05, 0.06), r"$\ell$", r"$u(\ell)$")
    ls = np.linspace(0, 6, 40)
    for (u0, eps), c in zip(E.U_CASES, ("cBar", "cWalk")):
        out.append(rf"\draw[{c},thick] " + _poly([(l, E.u_analytic(l, u0, eps)) for l in ls], sx, sy) + ";")
        for l in (1, 2, 5):
            out.append(rf"\fill[{c}] ({l * sx:.3f},{E.u_numeric(l, u0, eps) * sy:.3f}) circle (1.3pt);")
        out.append(rf"\node[{c},anchor=west] at (0.12,{u0 * sy + (0.12 if u0 > 0.03 else 0.1):.3f}) {{$u_0={u0:g}$}};")
    us = 1 / 36
    out.append(rf"\draw[gray,densely dashed] (0,{us * sy:.3f}) -- ({6 * sx:.3f},{us * sy:.3f}) node[right,black] {{$u^*=\frac1{{36}}$}};")
    g = E.g_numeric(math.exp(-5), 0.005, 1.0)
    out.append(rf"\draw[cExact,thick] ({5 * sx:.3f},{g * sy - 0.1:.3f}) -- ++(0.1,0.1) -- ++(-0.1,0.1) -- ++(-0.1,-0.1) -- cycle;")
    out.append(rf"\node[cExact,anchor=north] at ({5 * sx:.3f},{g * sy - 0.12:.3f}) {{$g_0(e^{{-5}})$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def fig_nu():
    """ν 的三種讀法對 ε（Ising）：一階、不展開、線性化矩陣；ε = 0 相切於 1/2。"""
    sx, sy = 5.2, 26.0
    y0 = 0.5
    out = [r"\begin{tikzpicture}[font=\scriptsize,>=Stealth,baseline=(current bounding box.center)]"]
    out += _axes(0, 1, 0.5, 0.6, sx, sy, (0, 0.5, 1), (0.5, 0.55, 0.6), r"$\epsilon$", r"$\nu$", fy="{:.2f}")
    es = np.linspace(0, 1, 30)
    curves = [("cBar,thick", lambda e: float(E.nu_first(e))),
              ("cWalk,thick,densely dashed", lambda e: float(E.nu_unexpanded(e))),
              ("cExact,thick,densely dotted", lambda e: 0.5 if e == 0 else E.nu_eigen(float(e)))]
    for sty, f in curves:
        out.append(rf"\draw[{sty}] " + _poly([(e, f(e) - y0) for e in es], sx, sy) + ";")
    out.append(rf"\node[cWalk,anchor=west] at ({1.0 * sx + 0.05:.3f},{(0.6 - y0) * sy:.3f}) {{不展開 $\frac1{{2-\epsilon/3}}$}};")
    out.append(rf"\node[cBar,anchor=west] at ({1.0 * sx + 0.05:.3f},{(float(E.nu_first(1)) - y0) * sy - 0.05:.3f}) {{一階 $\frac12+\frac\epsilon{{12}}$}};")
    out.append(rf"\node[cExact,anchor=west] at ({1.0 * sx + 0.05:.3f},{(E.nu_eigen(1.0) - y0) * sy - 0.12:.3f}) {{矩陣 $1/\lambda_r$}};")
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)


def tab_nu():
    rows = []
    for e in E.NU_EPS:
        rows.append(rf"${_fr(e)}$ & {_s(E.nu_first(e), 4)} & {_s(E.nu_unexpanded(e), 4)} & {_s(E.nu_eigen(float(e)), 4)}\\")
    return "\n".join([r"\begin{tabular}{@{}lrrr@{}}", r"\toprule",
                      r"$\epsilon$ & 一階 & 不展開 & 矩陣\\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def tab_on():
    rows = []
    for n in range(4):
        a, b = E.AB(n)
        rows.append(rf"{n} & {2 * n + 4} & {a} & {b} & ${_fr(E.nu_first(1, n))}={_s(E.nu_first(1, n), 4)}$ & "
                    rf"{_s(E.nu_unexpanded(1, n), 4)} & {_s(E.nu_eigen(1.0, n), 4)} & {_s(E.NU_3D_CITED[n], 3)}\\")
    return "\n".join([r"\begin{tabular}{@{}crrrlrrr@{}}", r"\toprule",
                      r"$n$ & $2n+4$ & $A=4(n+2)$ & $B=4(n+8)$ & 一階 (12.68) & 不展開 & 矩陣 & 已知（引用）\\",
                      r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"])


def tab_pc():
    rows = []
    num = {2: E.pc_numeric(2), F(5, 2): E.pc_numeric(2.5)}
    for d in (2, F(5, 2), 3, 4):
        v = f"${_s(num[d], 4)}$" if d in num else "--"
        rows.append(rf"${_fr(d)}$ & ${_fr(E.pc_exponent(d))}$ & ${_fr(E.pc_rhs(d))}$ & {v} & "
                    rf"{{\color{{cWrong}}${_fr(E.pc_exponent(d, True))}$}} & {{\color{{cWrong}}${_fr(E.pc_rhs(d, True))}$}}\\")
    return "\n".join([r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
                      r"$d$ & $d-3$ & 右邊 & 數值 & 書上 $2d-3$ & 書上右邊\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


ERRATA = [
    (r"p.~293，(12.17) 前", r"Figure 12.16(b)", r"Figure 12.11(b)"),
    (r"p.~295，(12.32)", r"$N\to N'e^{d\delta\ell}\approx N'(1+\delta\ell)$", r"$\approx N'(1+d\,\delta\ell)$"),
    (r"p.~295，(12.33)", r"$\vec q\to\vec q\,'e^{-d\delta\ell}$", r"$\vec q\to\vec q\,'e^{-\delta\ell}$"),
    (r"pp.~296--297，(12.36)(12.38)(12.39)", r"$(4-d)u-36u^2$", r"$(4-d)u-36d\,u^2$"),
    (r"p.~298，(12.47)", r"$\int_0^\ell u\,d\ell'=12d\int\cdots$", r"$12d\int_0^\ell u\,d\ell'=12d\int\cdots$"),
    (r"p.~298，(12.48)", r"$12d\int_{u_0}^u\frac{du'}{\epsilon u'-36du'^2}$", r"$12d\int_{u_0}^u\frac{u'\,du'}{\epsilon u'-36du'^2}$"),
    (r"p.~299，Ex.~12.4", r"(12.53) follows from (12.51)", r"(12.52)"),
    (r"p.~299，(12.60) 後", r"3D Ising $\nu\sim0.61$", r"約 0.630（文獻值）"),
    (r"p.~315，(12.89) 後", r"$u^2r^{2d-3}=r^{-1}(ur^{-(4-d)/2})^2$", r"$u^2r^{d-3}=r\,(ur^{-(4-d)/2})^2$"),
    (r"p.~320", r"$5!/(2!\,1!)$（第二次）", r"$5!/(3!\,2!)=10$"),
    (r"p.~322，(12:S-7)", r"$\mathcal T_1^{n_1}\mathcal T_3^{n_3}\mathcal T_3^{n_3}$，$n_1\ge1$",
     r"$\mathcal T_1^{n_1}\mathcal T_2^{n_2}\mathcal T_3^{n_3}$，$n_1\ge0$"),
]


def tab_errata():
    assert len(ERRATA) == 11
    rows = [rf"{a} & {{\color{{cWrong}}{b}}} & {{\color{{cExact}}{c}}}\\" for a, b, c in ERRATA]
    return "\n".join([r"\begin{tabular}{@{}lll@{}}", r"\toprule", r"位置 & 書上 & 應為\\", r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}"])


def tab_book():
    """全書回顧：每章一列（章、主角、動畫、關鍵結果）；數字取 walk_eps.recap_numbers()（＝ W1206 的回顧）。"""
    R = E.recap_numbers()
    key = {
        1: rf"$\langle x^2\rangle=N$；$P_{{10}}(0)=\frac{{252}}{{1024}}={R[1][0]}$，高斯近似 {R[1][1]}",
        2: rf"二維回原點 $\binom N{{N/2}}^2$：{'、'.join(R[2])}",
        3: rf"Pólya：一、二維一定回家；簡單立方晶格的回家機率 {R[3][0]}",
        4: rf"吸收牆與鏡像法；$x_0=2$、$L=5$ 先到左牆的機率 ${_fr(F(R[4][0]))}$",
        5: rf"持續長度 $l_p=\frac Q{{1-Q}}l$；$Q=\frac12$、$N=1000$：$\langle x^2\rangle={R[5][0]}$",
        6: rf"非球度 $A_d=\frac{{4+2d}}{{4+5d}}$，$A_3={_fr(F(R[6][0]))}={R[6][1]}$",
        7: rf"$c_{{14}}={R[7][0]}$；Flory $\nu_F=\frac3{{d+2}}$",
        8: rf"上臨界維度 $d_c=4$；二維的連結常數 $\mu={R[8][0]}$",
        9: rf"O($n$)，$n\to0$＝自迴避；Flory 0.600 對 0.588（差 2\%，{R[9][0]}）",
        10: rf"$\rho=|r|/2u$；$\nu'={_fr(F(R[10][0]))}$",
        11: rf"Koch 曲線 $\frac{{\ln4}}{{\ln3}}={R[11][0]}$；$\nu=1/\Lambda_r$",
        12: rf"$\nu=\frac12+\frac{{(n+2)\epsilon}}{{4(n+8)}}$：Ising {R[12][0]}、自迴避 {R[12][1]}",
    }
    rows = []
    for ch, name, _ in E.BOOK_RECAP:
        name = name.replace("ε", r"$\epsilon$")
        rows.append(rf"{ch} & {name} & {E.scene_range(ch).replace('–', '--')} & \raggedright {key[ch]}\tabularnewline")
    return "\n".join([r"\begin{tabular}{@{}rllp{10.5cm}@{}}", r"\toprule", r"章 & 主角 & 動畫 & 關鍵結果與數字\\", r"\midrule",
                      *rows, r"\midrule",
                      rf" & 全書 & {sum(E.BOOK_SCENES.values())} 幕 & 從 $\sqrt N$ 到 $N^{{0.588}}$（一階 $N^{{9/16}}$）\\",
                      r"\bottomrule", r"\end{tabular}"])


FIGURES = [
    ("TwelveFigShell", fig_shell),
    ("TwelveTabShell", tab_shell),
    ("TwelveFigToy", fig_toy),
    ("TwelveFigPairs", fig_pairs),
    ("TwelveTabCensus", tab_census),
    ("TwelveTabSfour", tab_s4),
    ("TwelveFigFlow", fig_flow),
    ("TwelveTabFixed", tab_fixed),
    ("TwelveFigU", fig_u),
    ("TwelveFigNu", fig_nu),
    ("TwelveTabNu", tab_nu),
    ("TwelveTabOn", tab_on),
    ("TwelveTabPc", tab_pc),
    ("TwelveTabErrata", tab_errata),
    ("TwelveTabBook", tab_book),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    check_numbers_ch12()
    parts = ["% figs.tex —— 由 ch12/make_figs.py 產生，請勿手動修改。\n"
             "% 重新產生：python ch12/make_figs.py（第 12 章，動畫 W1201–W1206）\n"]
    for name, func in FIGURES:
        body = func()
        assert "\n\n" not in body
        parts.append(f"\\newcommand{{\\{name}}}{{%\n{body}}}\n")
    check_tex_decimals()                                   # 圖表裡用到的小數也已記進 OK
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"數字檢查通過；已寫出 {OUT.name}：{len(FIGURES)} 張圖表")


if __name__ == "__main__":
    main()
