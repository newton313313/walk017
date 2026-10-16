"""
render_all.py —— 一次把場景輸出成影片（存到 renders/）。在 codes/ 裡執行；一次只跑一個 manimgl。

    python render_all.py                第 1 章 12 幕（W101–W112），1080p
    python render_all.py --ch1          同上（以後的章節是 --ch2、--ch3 …）
    python render_all.py W105 110       只輸出指定的幕（代號；"110" 也可以，等於 W110）
    python render_all.py --ch1-demo     另外輸出連續播放 Ch1Demo（W101–W112）；只要連續播放就再加 --demos-only
    python render_all.py -l             改用 480p（最快，用來預覽）；另有 -m（720p）、--uhd（4K）
    python render_all.py --force        影片已經存在也重新輸出（預設跳過：已完成的影片不覆寫）
"""
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

CH1 = [
    ("W101", "scenes/w101_coin_walk.py", "W101CoinWalk"),
    ("W102", "scenes/w102_counting_walks.py", "W102CountingWalks"),
    ("W103", "scenes/w103_mean_position.py", "W103MeanPosition"),
    ("W104", "scenes/w104_sqrt_n_law.py", "W104SqrtNLaw"),
    ("W105", "scenes/w105_moment_generating.py", "W105MomentGenerating"),
    ("W106", "scenes/w106_stirling.py", "W106Stirling"),
    ("W107", "scenes/w107_binomial_to_gauss.py", "W107BinomialToGauss"),
    ("W108", "scenes/w108_pascal_recursion.py", "W108PascalRecursion"),
    ("W109", "scenes/w109_diffusion_equation.py", "W109DiffusionEquation"),
    ("W110", "scenes/w110_fourier_count.py", "W110FourierCount"),
    ("W111", "scenes/w111_generating_function.py", "W111GeneratingFunction"),
    ("W112", "scenes/w112_steepest_descent.py", "W112SteepestDescent"),
]
CHAPTERS = {"--ch1": CH1}
DEMOS = {"--ch1-demo": ("main.py", "Ch1Demo")}
QUALITY_FLAGS = {"-l", "-m", "--hd", "--uhd"}


def main(argv):
    quality, force, demos_only = "--hd", False, False
    groups, keys, demos = [], [], []
    by_key = {k: (f, s) for k, f, s in sum(CHAPTERS.values(), [])}
    for arg in argv:
        a = arg.upper()
        if arg in QUALITY_FLAGS:
            quality = arg
        elif arg == "--force":
            force = True
        elif arg == "--demos-only":
            demos_only = True
        elif arg in CHAPTERS:
            groups.append(CHAPTERS[arg])
        elif arg in DEMOS:
            demos.append(DEMOS[arg])
        elif a in by_key or f"W{a}" in by_key:
            keys.append(a if a in by_key else f"W{a}")
        else:
            print(f"看不懂的參數：{arg}")
            print(__doc__)
            return 2

    jobs = []
    if not demos_only:
        jobs += [by_key[k] for k in keys]
        for g in groups:
            jobs += [(f, s) for _, f, s in g]
        if not keys and not groups:
            jobs += [(f, s) for _, f, s in sum(CHAPTERS.values(), [])]
    jobs += demos

    failed, skipped = [], []
    for file, scene in jobs:
        out = HERE / "renders" / f"{scene}.mp4"
        if out.exists() and not force:
            skipped.append(scene)
            print(f"略過（影片已存在，要重做請加 --force）：{scene}")
            continue
        print(f"\n=== {scene} ({quality}) ===", flush=True)
        start = time.time()
        result = subprocess.run(["manimgl", file, scene, "-w", quality], cwd=HERE)
        if result.returncode != 0 or not out.exists():
            failed.append(scene)
            print(f"失敗：{scene}")
        else:
            print(f"完成：{out}  （{time.time() - start:.0f} 秒）")

    print(f"\n略過 {len(skipped)} 支。" if skipped else "")
    print("全部完成。" if not failed else f"有 {len(failed)} 個場景失敗：{failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
