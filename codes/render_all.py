"""
render_all.py —— 一次把場景輸出成影片（存到 renders/）。在 codes/ 裡執行；一次只跑一個 manimgl。

    python render_all.py                全部章節的單幕（第 1 章 W101–W112、第 2 章 W201–W212、第 3 章 W301–W308、第 4 章 W401–W408、第 5 章 W501–W507、第 6 章 W601–W607、第 7 章 W701–W706、第 8 章 W801–W804、第 9 章 W901–W906），1080p
    python render_all.py --ch1          只做第 1 章；--ch2 … --ch9 同理（以後的章節是 --ch10 …）
    python render_all.py W105 210       只輸出指定的幕（代號；"210" 也可以，等於 W210）
    python render_all.py --ch2-demo     另外輸出連續播放 Ch2Demo（W201–W212；第 1 章是 --ch1-demo、第 3 章是 --ch3-demo、第 4 章是 --ch4-demo、第 5 章是 --ch5-demo、第 6 章是 --ch6-demo、第 7 章是 --ch7-demo、第 8 章是 --ch8-demo、第 9 章是 --ch9-demo）；只要連續播放就再加 --demos-only
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
CH2 = [
    ("W201", "scenes/w201_what_is_gf.py", "W201WhatIsGF"),
    ("W202", "scenes/w202_gf_is_answer.py", "W202GFIsAnswer"),
    ("W203", "scenes/w203_lattice_recursion.py", "W203LatticeRecursion"),
    ("W204", "scenes/w204_structure_function.py", "W204StructureFunction"),
    ("W205", "scenes/w205_fourier_inversion.py", "W205FourierInversion"),
    ("W206", "scenes/w206_gaussian_limit.py", "W206GaussianLimit"),
    ("W207", "scenes/w207_universality.py", "W207Universality"),
    ("W208", "scenes/w208_off_lattice.py", "W208OffLattice"),
    ("W209", "scenes/w209_diffusion_green.py", "W209DiffusionGreen"),
    ("W210", "scenes/w210_gaussian_integrals.py", "W210GaussianIntegrals"),
    ("W211", "scenes/w211_reciprocal_lattice.py", "W211ReciprocalLattice"),
    ("W212", "scenes/w212_singularities.py", "W212Singularities"),
]
CH3 = [
    ("W301", "scenes/w301_polya_question.py", "W301PolyaQuestion"),
    ("W302", "scenes/w302_visit_gf.py", "W302VisitGF"),
    ("W303", "scenes/w303_one_dimension.py", "W303OneDimension"),
    ("W304", "scenes/w304_two_dimensions.py", "W304TwoDimensions"),
    ("W305", "scenes/w305_dimension_divide.py", "W305DimensionDivide"),
    ("W306", "scenes/w306_watson.py", "W306Watson"),
    ("W307", "scenes/w307_sites_visited.py", "W307SitesVisited"),
    ("W308", "scenes/w308_footprint.py", "W308Footprint"),
]
CH4 = [
    ("W401", "scenes/w401_image_walker.py", "W401ImageWalker"),
    ("W402", "scenes/w402_two_walls.py", "W402TwoWalls"),
    ("W403", "scenes/w403_long_short_time.py", "W403LongShortTime"),
    ("W404", "scenes/w404_sphere_cube.py", "W404SphereCube"),
    ("W405", "scenes/w405_steady_tent.py", "W405SteadyTent"),
    ("W406", "scenes/w406_point_source.py", "W406PointSource"),
    ("W407", "scenes/w407_absorbing_sphere.py", "W407AbsorbingSphere"),
    ("W408", "scenes/w408_receptors.py", "W408Receptors"),
]
CH5 = [
    ("W501", "scenes/w501_biased_walk.py", "W501BiasedWalk"),
    ("W502", "scenes/w502_persistent_walk.py", "W502PersistentWalk"),
    ("W503", "scenes/w503_persistence_length.py", "W503PersistenceLength"),
    ("W504", "scenes/w504_telegrapher.py", "W504Telegrapher"),
    ("W505", "scenes/w505_run_and_tumble.py", "W505RunAndTumble"),
    ("W506", "scenes/w506_waiting_time.py", "W506WaitingTime"),
    ("W507", "scenes/w507_ctrw_gauss.py", "W507CtrwGauss"),
]
CH6 = [
    ("W601", "scenes/w601_shape_cloud.py", "W601ShapeCloud"),
    ("W602", "scenes/w602_gyration_tensor.py", "W602GyrationTensor"),
    ("W603", "scenes/w603_asphericity.py", "W603Asphericity"),
    ("W604", "scenes/w604_infinite_dim.py", "W604InfiniteDim"),
    ("W605", "scenes/w605_kramers_identity.py", "W605KramersIdentity"),
    ("W606", "scenes/w606_one_over_d.py", "W606OneOverD"),
    ("W607", "scenes/w607_fixman_saw.py", "W607FixmanSaw"),
]
CH7 = [
    ("W701", "scenes/w701_path_sum.py", "W701PathSum"),
    ("W702", "scenes/w702_self_avoiding.py", "W702SelfAvoiding"),
    ("W703", "scenes/w703_one_minus_f.py", "W703OneMinusF"),
    ("W704", "scenes/w704_dyson.py", "W704Dyson"),
    ("W705", "scenes/w705_mean_field.py", "W705MeanField"),
    ("W706", "scenes/w706_flory.py", "W706Flory"),
]
CH8 = [
    ("W801", "scenes/w801_non_gaussian.py", "W801NonGaussian"),
    ("W802", "scenes/w802_scaling_form.py", "W802ScalingForm"),
    ("W803", "scenes/w803_steepest_descent.py", "W803SteepestDescent"),
    ("W804", "scenes/w804_rescaling_flory.py", "W804RescalingFlory"),
]
CH9 = [
    ("W901", "scenes/w901_critical_exponent.py", "W901CriticalExponent"),
    ("W902", "scenes/w902_broken_symmetry.py", "W902BrokenSymmetry"),
    ("W903", "scenes/w903_scaling_laws.py", "W903ScalingLaws"),
    ("W904", "scenes/w904_landau_wilson.py", "W904LandauWilson"),
    ("W905", "scenes/w905_walk_from_magnet.py", "W905WalkFromMagnet"),
    ("W906", "scenes/w906_n_to_zero.py", "W906NToZero"),
]
CHAPTERS = {"--ch1": CH1, "--ch2": CH2, "--ch3": CH3, "--ch4": CH4, "--ch5": CH5, "--ch6": CH6, "--ch7": CH7, "--ch8": CH8, "--ch9": CH9}
DEMOS = {"--ch1-demo": ("main.py", "Ch1Demo"), "--ch2-demo": ("main.py", "Ch2Demo"), "--ch3-demo": ("main.py", "Ch3Demo"),
         "--ch4-demo": ("main.py", "Ch4Demo"),
         "--ch5-demo": ("main.py", "Ch5Demo"),
         "--ch6-demo": ("main.py", "Ch6Demo"),
         "--ch7-demo": ("main.py", "Ch7Demo"),
         "--ch8-demo": ("main.py", "Ch8Demo"),
         "--ch9-demo": ("main.py", "Ch9Demo")}
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
