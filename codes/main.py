"""
main.py —— walk017 全部場景的入口與連續播放（在 codes/ 裡執行）

    manimgl main.py               列出所有場景，選一個預覽
    manimgl main.py Ch1Demo -w    第 1 章 12 幕連續播放（W101–W112），輸出到 renders/Ch1Demo.mp4
    manimgl main.py Ch2Demo -w    第 2 章 12 幕連續播放（W201–W212），輸出到 renders/Ch2Demo.mp4
    manimgl main.py Ch3Demo -w    第 3 章 8 幕連續播放（W301–W308），輸出到 renders/Ch3Demo.mp4

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
