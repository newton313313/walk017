"""
main.py —— walk017 全部場景的入口與連續播放（在 codes/ 裡執行）

    manimgl main.py               列出所有場景，選一個預覽
    manimgl main.py Ch1Demo -w    第 1 章 12 幕連續播放（W101–W112），輸出到 renders/Ch1Demo.mp4

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


CH1 = [w101.W101CoinWalk, w102.W102CountingWalks, w103.W103MeanPosition, w104.W104SqrtNLaw,
       w105.W105MomentGenerating, w106.W106Stirling, w107.W107BinomialToGauss, w108.W108PascalRecursion,
       w109.W109DiffusionEquation, w110.W110FourierCount, w111.W111GeneratingFunction, w112.W112SteepestDescent]


class Ch1Demo(w101.W101CoinWalk):
    """第 1 章 12 幕連續播放（W101–W112）。"""

    def construct(self):
        play_sequence(self, CH1)
