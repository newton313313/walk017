"""
check_env.py —— 檢查這台電腦能不能跑 walk017（《Elements of the Random Walk》導讀動畫）

在 codes/ 資料夾內執行：
    python check_env.py            檢查套件、外部程式、字型、設定檔、combkit、數學驗證、講義的數字檢查
    python check_env.py --render   另外實際渲染一小段含 LaTeX 的場景，確認整條管線可用

每一項會顯示 OK / 缺少，最後列出要修的項目與建議指令。
"""
import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMBKIT = Path(r"C:\mypython\manim\projects\combkit")
results = []


def report(name, ok, detail="", fix=""):
    results.append((name, ok, detail, fix))
    mark = "OK  " if ok else "缺少"
    print(f"[{mark}] {name}" + (f"  —  {detail}" if detail else ""))


def run(cmd, cwd=HERE):
    try:
        env = dict(os.environ, PYTHONIOENCODING="utf-8")
        out = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                             timeout=600, cwd=cwd, env=env)
        return out.returncode, ((out.stdout or "") + (out.stderr or "")).strip()
    except Exception as err:          # noqa: BLE001
        return -1, str(err)


def main():
    render = "--render" in sys.argv
    sys.argv = sys.argv[:1]     # manimlib 在 import 時會解析命令列，先清掉本程式自己的參數
    sys.stdout.reconfigure(encoding="utf-8")
    print(f"Python：{sys.version.split()[0]}  （{sys.executable}）\n")

    report("Python 3.10 以上", sys.version_info >= (3, 10), sys.version.split()[0], "安裝 Python 3.10+（建議 Anaconda）")

    try:
        import manimlib  # noqa: F401
        from importlib.metadata import version
        report("ManimGL（manimlib）", True, f"manimgl {version('manimgl')}，位置 {Path(manimlib.__file__).parent}")
    except Exception as err:          # noqa: BLE001
        report("ManimGL（manimlib）", False, str(err).splitlines()[-1],
               "見 C:\\mypython\\manim\\環境與交接_ENVIRONMENT.md 第 2.1 節（不要升級或重裝 ManimGL）")

    for mod in ["av", "manimpango", "numpy", "scipy", "yaml", "pymupdf"]:
        try:
            importlib.import_module(mod)
            report(f"Python 套件 {mod}", True)
        except Exception as err:      # noqa: BLE001
            report(f"Python 套件 {mod}", False, str(err), "pip install -r C:\\mypython\\manim\\requirements.txt"
                   if mod != "pymupdf" else "pip install pymupdf（只有 ref/ 的轉檔與截圖需要）")

    for exe, fix in [
        ("manimgl", "裝好 ManimGL 後應該會有；確認 Python 的 Scripts 資料夾在 PATH 上"),
        ("ffmpeg", "conda install -c conda-forge ffmpeg"),
        ("latex", "安裝 MiKTeX：https://miktex.org/download"),
        ("xelatex", "安裝 MiKTeX（編譯講義用）"),
        ("dvisvgm", "安裝 MiKTeX（manim 把 LaTeX 轉成 SVG 用）"),
    ]:
        path = shutil.which(exe)
        report(f"外部程式 {exe}", path is not None, path or "", fix)

    if shutil.which("kpsewhich"):
        code, out = run(["kpsewhich", "sfrm1000.pfb"])
        found = [line for line in out.splitlines() if line.lower().endswith(".pfb")]
        report("MiKTeX 套件 cm-super", bool(found), found[0] if found else "",
               "miktex packages install cm-super，然後 manimgl --clear-cache")
    else:
        report("MiKTeX 套件 cm-super", False, "找不到 kpsewhich", "安裝 MiKTeX")

    try:
        import manimpango
        fonts = manimpango.list_fonts()
        for font in ["Microsoft JhengHei", "DFKai-SB"]:
            report(f"字型 {font}", font in fonts, "", "Windows 內建；其他系統請改 combkit/draw/styles.py 的 CJK_FONT 與講義的字型設定")
    except Exception as err:          # noqa: BLE001
        report("字型檢查", False, str(err))

    cfg = HERE / "custom_config.yml"
    try:
        cfg.read_bytes().decode("ascii")
        report("custom_config.yml 為純 ASCII", True)
    except Exception:                 # noqa: BLE001
        report("custom_config.yml 為純 ASCII", False, "", "把檔案裡的中文刪掉（ManimGL 用 cp950 讀它）")

    sys.path.insert(0, str(COMBKIT))
    try:
        import combkit
        from combkit.draw import scene_base  # noqa: F401
        report("combkit（場景共用的基底與版面檢查）", True, f"版本 {combkit.__version__}，位置 {COMBKIT}")
    except Exception as err:          # noqa: BLE001
        report("combkit", False, str(err).splitlines()[-1], f"確認 {COMBKIT} 存在（或 pip install -e 它）")

    code, out = run([sys.executable, str(HERE / "utils" / "walk1d.py")])
    ok = code == 0 and "walk1d（第 1 章）：計數、動差、Stirling、高斯近似、遞迴、傅立葉、生成函數、最陡下降全部驗證通過" in out
    report("數學驗證 utils/walk1d.py", ok, "公式、暴力列舉、書上的數字、場景卡的數字全部一致" if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_gf.py")])
    ok = code == 0 and "walk_gf（第 2 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_gf.py（第 2 章）", ok, "遞迴＝兩枚硬幣＝傅立葉、書上的數字、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_rec.py")])
    ok = code == 0 and "walk_rec（第 3 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_rec.py（第 3 章）", ok, "生成函數＝暴力列舉＝逐步遞迴、Watson 積分三種算法、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_bc.py")])
    ok = code == 0 and "walk_bc（第 4 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_bc.py（第 4 章）", ok, "鏡像＝暴力列舉＝逐步遞迴＝模態和、穩定態的精確解、Green 函數、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_var.py")])
    ok = code == 0 and "walk_var（第 5 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_var.py（第 5 章）", ok, "精確＝逐步遞迴＝轉移矩陣、電報方程與二三維的公式＝Laplace 反轉＝模擬、守恆、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_shape.py")])
    ok = code == 0 and "walk_shape（第 6 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_shape.py（第 6 章）", ok, "定義式＝兩兩差＝Kramers、非球度的公式＝精確動差＝暴力列舉＝模擬、1/d 展開與 Fixman＝模擬、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_saw.py")])
    ok = code == 0 and "walk_saw（第 7 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_saw.py（第 7 章）", ok, "切段＝公式＝暴力列舉、自迴避的 DFS＝文獻＝暴力、1−f 展開三種算法＝c_N、Dyson、平均場與 Flory 的數值＝解析、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_scale.py")])
    ok = code == 0 and "walk_scale（第 8 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_scale.py（第 8 章）", ok, "χ(q) 數值＝閉式、模擬寬度＝穩定分布、逆 Fourier＝換變數的尾巴、塌縮、格子 G＝級數＝Yukawa、最陡下降＝數值極值＝精確二項式、比值法、重新標度的因子、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_crit.py")])
    ok = code == 0 and "walk_crit（第 9 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_crit.py（第 9 章）", ok, "平均場＝牛頓＝漸近、Onsager＝低溫展開、小 Ising 列舉＝閉式、Metropolis＝精確列舉、一維 G＝轉移矩陣＝環、標度律（Fraction）、Widom 解析＝差分、Landau 根＝格點、(9.34) quad＝梯形、係數 lgamma＝遞迴、ξ² 兩種列舉、高斯積分＝dblquad＝Gauss–Hermite、Trotter 階數、梯子 n＝1、2 暴力＝閉式、n→0＝自迴避、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_cond.py")])
    ok = code == 0 and "walk_cond（第 10 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_cond.py（第 10 章）", ok, "環面列舉 8 重對稱＝不用對稱＝4^N 暴力＝Gray 碼、平均場 erfcx＝quad、n→0 徑向＝拋物柱＝一維、註 3＝walk_saw 的 G_MF、tanh＝BVP、自旋波 Fraction＝numpy、D 的 FFT＝直接求和＝無限格子＋Ewald、erfc 形 vs 精確步數＝三維遞迴、G(R) Bessel＝G_c、ρ 比值＝差分、指數兩條路（Fraction）、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_rg.py")])
    ok = code == 0 and "walk_rg（第 11 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_rg.py（第 11 章）", ok, "Koch／Cantor 質量－半徑＝盒子計數＝複本公式、⟨R_g²⟩ 公式＝4^N 列舉、Landau 三次方程＝二分法＝walk_crit、塌縮、完成平方＝數值積分、Gaussian 環的行列式＝抽格、Ising 環 暴力＝轉移矩陣＝連乘、x 遞迴＝K 遞迴、ξ(K')＝ξ/2、流 解析＝solve_ivp、(11.87) 殘差、指數公式＝Widom＝walk_crit 指數表、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])
    code, out = run([sys.executable, str(HERE / "utils" / "walk_eps.py")])
    ok = code == 0 and "walk_eps（第 12 章）" in out and "全部驗證通過" in out
    report("數學驗證 utils/walk_eps.py（第 12 章）", ok, "兩變數玩具 閉式＝數值積分＝擬合、殼層 指數＝格點計數、配對列舉＝(12:S-4)(12:S-6)、連結叢集四種算法、O(n) 縮併 2n+4 與 4(n+8)、不動點 Fraction＝fsolve、線性化 閉式＝numpy＝差分＝流的時間、u(ℓ)／t(ℓ) 解析＝數值、ν 一階公式與 Flory 展開、(12.86)＝(12.39)、冪次數算、全書回顧、筆誤的更正、場景卡的數字全部一致"
           if ok else out[-300:])

    notes = HERE.parent / "notes"
    code, out = run([sys.executable, str(notes / "ch1" / "make_figs.py")], cwd=notes)
    ok = code == 0 and "數字檢查通過" in out
    report("講義第 1 章的數字檢查 notes/ch1/make_figs.py", ok, "並重新產生 ch1/figs.tex" if ok else out[-300:])
    code, out = run([sys.executable, str(notes / "ch2" / "make_figs.py")], cwd=notes)
    ok = code == 0 and "數字檢查通過" in out
    report("講義第 2 章的數字檢查 notes/ch2/make_figs.py", ok, "並重新產生 ch2/figs.tex" if ok else out[-300:])
    code, out = run([sys.executable, str(notes / "ch3" / "make_figs.py")], cwd=notes)
    ok = code == 0 and "數字檢查通過" in out
    report("講義第 3 章的數字檢查 notes/ch3/make_figs.py", ok, "並重新產生 ch3/figs.tex" if ok else out[-300:])
    for ch in (4, 5, 6, 7, 8, 9, 10, 11, 12):
        code, out = run([sys.executable, str(notes / f"ch{ch}" / "make_figs.py")], cwd=notes)
        ok = code == 0 and "數字檢查通過" in out
        report(f"講義第 {ch} 章的數字檢查 notes/ch{ch}/make_figs.py", ok, f"並重新產生 ch{ch}/figs.tex" if ok else out[-300:])

    if render:
        print("\n渲染測試中（約 30 秒）……")
        code, out = run(["manimgl", "scenes/w102_counting_walks.py", "W102CountingWalks",
                         "-w", "-l", "-n", "0,3", "--file_name", "_check_env_test"])
        video = HERE / "renders" / "_check_env_test.mp4"
        ok = code == 0 and video.exists() and video.stat().st_size > 0
        report("渲染測試（含 LaTeX 公式）", ok, "成功輸出測試影片" if ok else out[-400:],
               "看上面的錯誤；LaTeX 錯誤先試 manimgl --clear-cache")
        for leftover in (HERE / "renders").glob("_check_env_test*"):
            leftover.unlink()

    bad = [r for r in results if not r[1]]
    print()
    if not bad:
        print("全部通過：可以在 codes/ 裡執行 manimgl，或用 render_all.py／render_all.bat 輸出影片。")
        return 0
    print(f"有 {len(bad)} 項要處理：")
    for name, _, detail, fix in bad:
        print(f"  - {name}：{fix or detail}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
