"""
main.py —— walk017 全部場景的入口與連續播放（在 codes/ 裡執行）

    manimgl main.py               列出所有場景，選一個預覽
    manimgl main.py Ch1Demo -w    第 1 章 12 幕連續播放（W101–W112），輸出到 renders/Ch1Demo.mp4
    manimgl main.py Ch2Demo -w    第 2 章 12 幕連續播放（W201–W212），輸出到 renders/Ch2Demo.mp4
    manimgl main.py Ch3Demo -w    第 3 章 8 幕連續播放（W301–W308），輸出到 renders/Ch3Demo.mp4
    manimgl main.py Ch4Demo -w    第 4 章 8 幕連續播放（W401–W408），輸出到 renders/Ch4Demo.mp4
    manimgl main.py Ch5Demo -w    第 5 章 7 幕連續播放（W501–W507），輸出到 renders/Ch5Demo.mp4
    manimgl main.py Ch6Demo -w    第 6 章 7 幕連續播放（W601–W607），輸出到 renders/Ch6Demo.mp4
    manimgl main.py Ch7Demo -w    第 7 章 6 幕連續播放（W701–W706），輸出到 renders/Ch7Demo.mp4
    manimgl main.py Ch8Demo -w    第 8 章 4 幕連續播放（W801–W804），輸出到 renders/Ch8Demo.mp4
    manimgl main.py Ch9Demo -w    第 9 章 6 幕連續播放（W901–W906），輸出到 renders/Ch9Demo.mp4

ManimGL 只列出「定義在這個檔案裡」的場景，所以每個場景在這裡都有一個同名子類別。
連續播放是拿 Ch1Demo 這個物件依序去跑各場景的 construct，所以場景的輔助函式都寫在模組層級。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))          # codes/，讓 scenes.*、utils.* 可以 import
sys.path.insert(0, r"C:\mypython\manim\projects\combkit")

from manimlib import *

from scenes import w101_coin_walk as w101
from scenes import w102_counting_walks as w102
from scenes import w103_mean_position as w103
from scenes import w104_sqrt_n_law as w104
from scenes import w105_moment_generating as w105
from scenes import w106_stirling as w106
from scenes import w107_binomial_to_gauss as w107
from scenes import w108_pascal_recursion as w108
from scenes import w109_diffusion_equation as w109
from scenes import w110_fourier_count as w110
from scenes import w111_generating_function as w111
from scenes import w112_steepest_descent as w112
from scenes import w201_what_is_gf as w201
from scenes import w202_gf_is_answer as w202
from scenes import w203_lattice_recursion as w203
from scenes import w204_structure_function as w204
from scenes import w205_fourier_inversion as w205
from scenes import w206_gaussian_limit as w206
from scenes import w207_universality as w207
from scenes import w208_off_lattice as w208
from scenes import w209_diffusion_green as w209
from scenes import w210_gaussian_integrals as w210
from scenes import w211_reciprocal_lattice as w211
from scenes import w212_singularities as w212
from scenes import w301_polya_question as w301
from scenes import w302_visit_gf as w302
from scenes import w303_one_dimension as w303
from scenes import w304_two_dimensions as w304
from scenes import w305_dimension_divide as w305
from scenes import w306_watson as w306
from scenes import w307_sites_visited as w307
from scenes import w308_footprint as w308
from scenes import w401_image_walker as w401
from scenes import w402_two_walls as w402
from scenes import w403_long_short_time as w403
from scenes import w404_sphere_cube as w404
from scenes import w405_steady_tent as w405
from scenes import w406_point_source as w406
from scenes import w407_absorbing_sphere as w407
from scenes import w408_receptors as w408
from scenes import w501_biased_walk as w501
from scenes import w502_persistent_walk as w502
from scenes import w503_persistence_length as w503
from scenes import w504_telegrapher as w504
from scenes import w505_run_and_tumble as w505
from scenes import w506_waiting_time as w506
from scenes import w507_ctrw_gauss as w507
from scenes import w601_shape_cloud as w601
from scenes import w602_gyration_tensor as w602
from scenes import w603_asphericity as w603
from scenes import w604_infinite_dim as w604
from scenes import w605_kramers_identity as w605
from scenes import w606_one_over_d as w606
from scenes import w607_fixman_saw as w607
from scenes import w701_path_sum as w701
from scenes import w702_self_avoiding as w702
from scenes import w703_one_minus_f as w703
from scenes import w704_dyson as w704
from scenes import w705_mean_field as w705
from scenes import w706_flory as w706
from scenes import w801_non_gaussian as w801
from scenes import w802_scaling_form as w802
from scenes import w803_steepest_descent as w803
from scenes import w804_rescaling_flory as w804
from scenes import w901_critical_exponent as w901
from scenes import w902_broken_symmetry as w902
from scenes import w903_scaling_laws as w903
from scenes import w904_landau_wilson as w904
from scenes import w905_walk_from_magnet as w905
from scenes import w906_n_to_zero as w906


def play_sequence(scene, classes):
    """依序借用各場景的 construct 連續播放；場景之間清空畫面（fade_all 會重設字幕與標題）、鏡頭歸位。"""
    for k, cls in enumerate(classes):
        cls.construct(scene)
        if k < len(classes) - 1:
            scene.fade_all(run_time=0.8)
            scene.frame.to_default_state()
            scene.wait(0.4)


class W101CoinWalk(w101.W101CoinWalk):
    pass


class W102CountingWalks(w102.W102CountingWalks):
    pass


class W103MeanPosition(w103.W103MeanPosition):
    pass


class W104SqrtNLaw(w104.W104SqrtNLaw):
    pass


class W105MomentGenerating(w105.W105MomentGenerating):
    pass


class W106Stirling(w106.W106Stirling):
    pass


class W107BinomialToGauss(w107.W107BinomialToGauss):
    pass


class W108PascalRecursion(w108.W108PascalRecursion):
    pass


class W109DiffusionEquation(w109.W109DiffusionEquation):
    pass


class W110FourierCount(w110.W110FourierCount):
    pass


class W111GeneratingFunction(w111.W111GeneratingFunction):
    pass


class W112SteepestDescent(w112.W112SteepestDescent):
    pass


class W201WhatIsGF(w201.W201WhatIsGF):
    pass


class W202GFIsAnswer(w202.W202GFIsAnswer):
    pass


class W203LatticeRecursion(w203.W203LatticeRecursion):
    pass


class W204StructureFunction(w204.W204StructureFunction):
    pass


class W205FourierInversion(w205.W205FourierInversion):
    pass


class W206GaussianLimit(w206.W206GaussianLimit):
    pass


class W207Universality(w207.W207Universality):
    pass


class W208OffLattice(w208.W208OffLattice):
    pass


class W209DiffusionGreen(w209.W209DiffusionGreen):
    pass


class W210GaussianIntegrals(w210.W210GaussianIntegrals):
    pass


class W211ReciprocalLattice(w211.W211ReciprocalLattice):
    pass


class W212Singularities(w212.W212Singularities):
    pass


class W301PolyaQuestion(w301.W301PolyaQuestion):
    pass


class W302VisitGF(w302.W302VisitGF):
    pass


class W303OneDimension(w303.W303OneDimension):
    pass


class W304TwoDimensions(w304.W304TwoDimensions):
    pass


class W305DimensionDivide(w305.W305DimensionDivide):
    pass


class W306Watson(w306.W306Watson):
    pass


class W307SitesVisited(w307.W307SitesVisited):
    pass


class W308Footprint(w308.W308Footprint):
    pass


class W401ImageWalker(w401.W401ImageWalker):
    pass


class W402TwoWalls(w402.W402TwoWalls):
    pass


class W403LongShortTime(w403.W403LongShortTime):
    pass


class W404SphereCube(w404.W404SphereCube):
    pass


class W405SteadyTent(w405.W405SteadyTent):
    pass


class W406PointSource(w406.W406PointSource):
    pass


class W407AbsorbingSphere(w407.W407AbsorbingSphere):
    pass


class W408Receptors(w408.W408Receptors):
    pass


class W501BiasedWalk(w501.W501BiasedWalk):
    pass


class W502PersistentWalk(w502.W502PersistentWalk):
    pass


class W503PersistenceLength(w503.W503PersistenceLength):
    pass


class W504Telegrapher(w504.W504Telegrapher):
    pass


class W505RunAndTumble(w505.W505RunAndTumble):
    pass


class W506WaitingTime(w506.W506WaitingTime):
    pass


class W507CtrwGauss(w507.W507CtrwGauss):
    pass


class W601ShapeCloud(w601.W601ShapeCloud):
    pass


class W602GyrationTensor(w602.W602GyrationTensor):
    pass


class W603Asphericity(w603.W603Asphericity):
    pass


class W604InfiniteDim(w604.W604InfiniteDim):
    pass


class W605KramersIdentity(w605.W605KramersIdentity):
    pass


class W606OneOverD(w606.W606OneOverD):
    pass


class W607FixmanSaw(w607.W607FixmanSaw):
    pass


class W701PathSum(w701.W701PathSum):
    pass


class W702SelfAvoiding(w702.W702SelfAvoiding):
    pass


class W703OneMinusF(w703.W703OneMinusF):
    pass


class W704Dyson(w704.W704Dyson):
    pass


class W705MeanField(w705.W705MeanField):
    pass


class W706Flory(w706.W706Flory):
    pass


class W801NonGaussian(w801.W801NonGaussian):
    pass


class W802ScalingForm(w802.W802ScalingForm):
    pass


class W803SteepestDescent(w803.W803SteepestDescent):
    pass


class W804RescalingFlory(w804.W804RescalingFlory):
    pass


class W901CriticalExponent(w901.W901CriticalExponent):
    pass


class W902BrokenSymmetry(w902.W902BrokenSymmetry):
    pass


class W903ScalingLaws(w903.W903ScalingLaws):
    pass


class W904LandauWilson(w904.W904LandauWilson):
    pass


class W905WalkFromMagnet(w905.W905WalkFromMagnet):
    pass


class W906NToZero(w906.W906NToZero):
    pass


CH1 = [w101.W101CoinWalk, w102.W102CountingWalks, w103.W103MeanPosition, w104.W104SqrtNLaw,
       w105.W105MomentGenerating, w106.W106Stirling, w107.W107BinomialToGauss, w108.W108PascalRecursion,
       w109.W109DiffusionEquation, w110.W110FourierCount, w111.W111GeneratingFunction, w112.W112SteepestDescent]


class Ch1Demo(w101.W101CoinWalk):
    """第 1 章 12 幕連續播放（W101–W112）。"""

    def construct(self):
        play_sequence(self, CH1)


CH2 = [w201.W201WhatIsGF, w202.W202GFIsAnswer, w203.W203LatticeRecursion, w204.W204StructureFunction,
       w205.W205FourierInversion, w206.W206GaussianLimit, w207.W207Universality, w208.W208OffLattice,
       w209.W209DiffusionGreen, w210.W210GaussianIntegrals, w211.W211ReciprocalLattice, w212.W212Singularities]


class Ch2Demo(w201.W201WhatIsGF):
    """第 2 章 12 幕連續播放（W201–W212）。"""

    def construct(self):
        play_sequence(self, CH2)


CH3 = [w301.W301PolyaQuestion, w302.W302VisitGF, w303.W303OneDimension, w304.W304TwoDimensions,
       w305.W305DimensionDivide, w306.W306Watson, w307.W307SitesVisited, w308.W308Footprint]


class Ch3Demo(w301.W301PolyaQuestion):
    """第 3 章 8 幕連續播放（W301–W308）。"""

    def construct(self):
        play_sequence(self, CH3)


CH4 = [w401.W401ImageWalker, w402.W402TwoWalls, w403.W403LongShortTime, w404.W404SphereCube, w405.W405SteadyTent, w406.W406PointSource, w407.W407AbsorbingSphere, w408.W408Receptors]


class Ch4Demo(w401.W401ImageWalker):
    """第 4 章 8 幕連續播放（W401–W408）。"""

    def construct(self):
        play_sequence(self, CH4)


CH5 = [w501.W501BiasedWalk, w502.W502PersistentWalk, w503.W503PersistenceLength, w504.W504Telegrapher, w505.W505RunAndTumble, w506.W506WaitingTime, w507.W507CtrwGauss]


class Ch5Demo(w501.W501BiasedWalk):
    """第 5 章 7 幕連續播放（W501–W507）。"""

    def construct(self):
        play_sequence(self, CH5)


CH6 = [w601.W601ShapeCloud, w602.W602GyrationTensor, w603.W603Asphericity, w604.W604InfiniteDim, w605.W605KramersIdentity, w606.W606OneOverD, w607.W607FixmanSaw]


class Ch6Demo(w601.W601ShapeCloud):
    """第 6 章 7 幕連續播放（W601–W607）。"""

    def construct(self):
        play_sequence(self, CH6)


CH7 = [w701.W701PathSum, w702.W702SelfAvoiding, w703.W703OneMinusF, w704.W704Dyson, w705.W705MeanField, w706.W706Flory]


class Ch7Demo(w701.W701PathSum):
    """第 7 章 6 幕連續播放（W701–W706）。"""

    def construct(self):
        play_sequence(self, CH7)


CH8 = [w801.W801NonGaussian, w802.W802ScalingForm, w803.W803SteepestDescent, w804.W804RescalingFlory]


class Ch8Demo(w801.W801NonGaussian):
    """第 8 章 4 幕連續播放（W801–W804）。"""

    def construct(self):
        play_sequence(self, CH8)


CH9 = [w901.W901CriticalExponent, w902.W902BrokenSymmetry, w903.W903ScalingLaws, w904.W904LandauWilson, w905.W905WalkFromMagnet, w906.W906NToZero]


class Ch9Demo(w901.W901CriticalExponent):
    """第 9 章 6 幕連續播放（W901–W906）。"""

    def construct(self):
        play_sequence(self, CH9)
