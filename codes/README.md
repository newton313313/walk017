# walk017：《Elements of the Random Walk》導讀動畫（ManimGL）

J. Rudnick & G. Gaspari, *Elements of the Random Walk: An Introduction for Advanced Students and Researchers*（Cambridge University Press 2004）
的繁體中文導讀動畫與講義。目前完成**第 1–3 章**：第 1 章（Introduction to techniques）12 幕動畫 W101–W112＋連續播放 `Ch1Demo`；
第 2 章（Generating functions I）12 幕動畫 W201–W212＋連續播放 `Ch2Demo`；第 3 章（Generating functions II：回家、維度、走過幾個點）
8 幕動畫 W301–W308＋連續播放 `Ch3Demo`；第 4 章（邊界條件、穩定態與靜電類比）8 幕動畫 W401–W408＋連續播放 `Ch4Demo`；第 5 章（漫步的變化形：偏向、慣性、連續時間）7 幕動畫 W501–W507＋連續播放 `Ch5Demo`；第 6 章（漫步的形狀：迴轉半徑張量、非球度、1/d 展開）7 幕動畫 W601–W607＋連續播放 `Ch6Demo`；第 7 章（路徑積分與自迴避）6 幕動畫 W701–W706＋連續播放 `Ch7Demo`；第 8 章（標度入門）4 幕動畫 W801–W804＋連續播放 `Ch8Demo`；第 9 章（漫步的標度與臨界現象）6 幕動畫 W901–W906＋連續播放 `Ch9Demo`；第 10 章（O(n) 模型：平均場與自旋波）5 幕動畫 W1001–W1005＋連續播放 `Ch10Demo`；第 11 章（標度、碎形與重整化）6 幕動畫 W1101–W1106＋連續播放 `Ch11Demo`；第 12 章（重整化群的 ε 展開）6 幕動畫 W1201–W1206＋連續播放 `Ch12Demo`；導讀講義 `../notes/walk_導讀.pdf`（全書十二章）。

**路徑的寫法：** 本檔裡沒有寫磁碟機的相對路徑，一律相對於本資料夾 `C:\mypython\manim\projects\walk017\codes`；
寫成 `../notes/…`、`../ref/…` 的在上一層 `C:\mypython\manim\projects\walk017`。權威說明是上一層的 `交接說明_HANDOFF.md`。

## 1. 環境與依賴

- ManimGL **1.7.2**（3Blue1Brown 版，commit fafa083，`manimgl` 指令、`from manimlib import *`），**不是** Manim Community；不要升級或重裝。
  全機環境見 `C:\mypython\manim\環境與交接_ENVIRONMENT.md`、`C:\mypython\manim\新計畫起步指南_NEW_PROJECT_GUIDE.md`。
- 共用工具包 **combkit 1.13.0 以上**（`C:\mypython\manim\projects\combkit`）：場景基底 `ChuScene`、中文字 `zh()`／`mixed()`、版面檢查 `assert_in_frame`。
  場景檔開頭用 `sys.path.insert` 載入它，不必安裝。
- LaTeX：MiKTeX（`latex`＋`dvisvgm`＋`cm-super`；講義用 `xelatex`＋`xeCJK`、`tcolorbox`、`tikz`）。字型：微軟正黑體、標楷體。
- 檢查：`python check_env.py`（加 `--render` 會實際渲染一小段）。

## 2. 指令（都在本資料夾執行）

| 做什麼 | 指令 |
|---|---|
| 預覽一幕 | `manimgl scenes/w101_coin_walk.py W101CoinWalk`（Ctrl+Q 結束、滾輪縮放、按住 f 平移、r 重設） |
| 輸出一幕 1080p | `manimgl scenes/w101_coin_walk.py W101CoinWalk -w --hd`（低畫質把 `--hd` 換 `-l`） |
| 一次輸出全部 | `python render_all.py`（三章 32 幕）；`--ch1`／`--ch2`／`--ch3` 只做一章；`--ch1-demo`／`--ch2-demo`／`--ch3-demo` 加連續播放；`-l` 480p；`W105 210` 只做指定的幕；已存在的影片預設**跳過**，要重做加 `--force` |
| 按兩下的選單 | `render_all.bat`（會自動啟用 conda） |
| 連續播放 | `manimgl main.py Ch1Demo -w --hd`、`manimgl main.py Ch2Demo -w --hd`、`manimgl main.py Ch3Demo -w --hd`、`manimgl main.py Ch4Demo -w --hd`、`manimgl main.py Ch5Demo -w --hd`、`manimgl main.py Ch6Demo -w --hd`、`manimgl main.py Ch7Demo -w --hd`、`manimgl main.py Ch8Demo -w --hd`、`manimgl main.py Ch9Demo -w --hd`、`manimgl main.py Ch10Demo -w --hd`、`manimgl main.py Ch11Demo -w --hd`、`manimgl main.py Ch12Demo -w --hd` |
| 數學自我測試 | `python utils/walk1d.py`（第 1 章）、`python utils/walk_gf.py`（第 2 章）、`python utils/walk_rec.py`（第 3 章，約 2 秒）、`python utils/walk_bc.py`（第 4 章，約 1 秒）、`python utils/walk_var.py`（第 5 章，約 1 秒）、`python utils/walk_shape.py`（第 6 章，約 2 秒）、`python utils/walk_saw.py`（第 7 章，約 1 秒）、`python utils/walk_scale.py`（第 8 章，約 1 秒）、`python utils/walk_crit.py`（第 9 章）、`python utils/walk_cond.py`（第 10 章）、`python utils/walk_rg.py`（第 11 章）、`python utils/walk_eps.py`（第 12 章） |
| 環境檢查 | `python check_env.py`、`python check_env.py --render` |
| 試渲染＋總覽圖 | 在 `C:\mypython\manim\projects\combkit` 執行 `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，看 `runs/<類別>/sheet.png` |
| 講義 | 在 `../notes` 執行 `python ch1/make_figs.py`、`python ch2/make_figs.py`、`python ch3/make_figs.py`，再 `xelatex walk_導讀.tex` 兩次 |
| 參考文獻轉檔 | 在 `../ref` 執行 `python convert_to_md.py`（全書 → `md/`）；`python page_png.py 5`（書上第 5 頁 → `png/p005.png`） |

## 3. 場景（第 1 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W101 | `w101_coin_walk.py` | `W101CoinWalk` | 擲硬幣走路：主例 HTHHTHTTHH、十步的終點分布 | §1.1，pp. 1–2 | 2:09 |
| W102 | `w102_counting_walks.py` | `W102CountingWalks` | 數路線：二項式係數、洗牌論證、$C(10,0)=252$、偏硬幣（練習 1.1） | §1.1，pp. 2–5 | 2:09 |
| W103 | `w103_mean_position.py` | `W103MeanPosition` | 平均位置是 0：20 人的鎖步系綜、對稱抵消 | §1.2，pp. 5–6 | 2:04 |
| W104 | `w104_sqrt_n_law.py` | `W104SqrtNLaw` | 均方位移與 √N 定律：$\Delta_j\Delta_k$ 格子、2000 人的扇形 | §1.2，pp. 6–7 | 2:07 |
| W105 | `w105_moment_generating.py` | `W105MomentGenerating` | 用生成函數算所有動差：$\cosh(y/2)^N$、比值 → 3（練習 1.2、1.3） | §1.2，pp. 7–10 | 2:11 |
| W106 | `w106_stirling.py` | `W106Stirling` | Stirling 公式：$w^{10}e^{-w}$ 與高斯、$\ln10!$（練習 1.4） | §1.3，pp. 10–13 | 2:05 |
| W107 | `w107_binomial_to_gauss.py` | `W107BinomialToGauss` | 二項分布變成高斯、奇偶性的因子 2 | §1.3，pp. 12–13 | 2:18 |
| W108 | `w108_pascal_recursion.py` | `W108PascalRecursion` | 遞迴關係與 Pascal 三角形、Markov | §1.4，p. 13 | 2:05 |
| W109 | `w109_diffusion_equation.py` | `W109DiffusionEquation` | 從遞迴到擴散方程（練習 1.5） | §1.4，pp. 13–14 | 2:03 |
| W110 | `w110_fourier_count.py` | `W110FourierCount` | 傅立葉：用積分數走法、兩個鞍點 | §1.5，pp. 15–16 | 2:04 |
| W111 | `w111_generating_function.py` | `W111GeneratingFunction` | 生成函數的奇點、Cauchy 積分、動差（練習 1.6、1.7） | §1.5.1–1.5.2，pp. 16–20 | 2:02 |
| W112 | `w112_steepest_descent.py` | `W112SteepestDescent` | 最陡下降法（鞍點法）＋第 1 章回顧 | §1.6，pp. 20–23 | 2:08 |

連續播放：`main.py` 的 `Ch1Demo`（W101–W112），長度 25:40。影片都在 `renders/`（不進 git）。

## 3b. 場景（第 2 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W201 | `w201_what_is_gf.py` | `W201WhatIsGF` | 生成函數：一個函數裝下整個數列（$e^x$、Legendre、練習 2.1–2.3） | §2.1.1，pp. 25–27 | 2:13 |
| W202 | `w202_gf_is_answer.py` | `W202GFIsAnswer` | 生成函數本身就是答案：巨配分函數、會消失的走路人 | §2.1.1–2.1.2，pp. 27–29 | 2:06 |
| W203 | `w203_lattice_recursion.py` | `W203LatticeRecursion` | 遞迴：最後一步來自四個鄰居、$C(4;\vec0)=36$、兩枚硬幣 | §2.1.2，pp. 29–31 | 2:09 |
| W204 | `w204_structure_function.py` | `W204StructureFunction` | 平移對稱與傅立葉：$g=1/(1-z\chi)$、五種格子的 $\chi$（練習 2.6） | §2.1.2，pp. 31–33 | 2:06 |
| W205 | `w205_fourier_inversion.py` | `W205FourierInversion` | 從 $\chi^N$ 取回走法數：布里淵區積分、主例 63504 | §2.1.3，pp. 33–34 | 2:05 |
| W206 | `w206_gaussian_limit.py` | `W206GaussianLimit` | 高斯極限與奇偶性的因子 2（布里淵區角落） | §2.1.3，pp. 34–35 | 2:08 |
| W207 | `w207_universality.py` | `W207Universality` | 普適性：sc、bcc、fcc 同一條高斯；$q^4$ 修正 | §2.1.4，pp. 35–36 | 2:16 |
| W208 | `w208_off_lattice.py` | `W208OffLattice` | 離開格子：Cauchy 步長、$\chi=e^{-qr_0}$、距離 ∝ $N$（練習 2.7） | §2.1.4，pp. 36–38 | 2:20 |
| W209 | `w209_diffusion_green.py` | `W209DiffusionGreen` | 擴散方程與格林函數（練習 2.8） | §2.1.5，pp. 38–40 | 2:12 |
| W210 | `w210_gaussian_integrals.py` | `W210GaussianIntegrals` | 補充 1：高斯積分與它的生成函數 | §2.2，pp. 41–42 | 2:06 |
| W211 | `w211_reciprocal_lattice.py` | `W211ReciprocalLattice` | 補充 2：倒晶格、布里淵區、單位根 | §2.3，pp. 42–47 | 2:08 |
| W212 | `w212_singularities.py` | `W212Singularities` | 補充 3：從奇點讀出係數的成長＋第 2 章回顧 | §2.4，pp. 47–50 | 2:19 |

連續播放：`main.py` 的 `Ch2Demo`（W201–W212），長度 26:20。計畫與場景卡：`../notes/plan_ch02.md`。

## 3c. 場景（第 3 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W301 | `w301_polya_question.py` | `W301PolyaQuestion` | 醉漢回得了家，醉鳥不一定：16 條路線、1000 人模擬、Pólya 定理 | §3.1、§3.4，pp. 51, 55 | 2:03 |
| W302 | `w302_visit_gf.py` | `W302VisitGF` | 新的生成函數：用 $(1-t)^M$ 數到訪次數、排容（Fig. 3.1）、幾何級數 (3.9) | §3.2–3.3，pp. 51–54 | 2:14 |
| W303 | `w303_one_dimension.py` | `W303OneDimension` | 一維一定回家：從不回家的走法 252、$\sqrt{2/\pi N}$（(3.24) 的筆誤）、偏硬幣（練習 3.2–3.4） | §3.4，pp. 55–58 | 2:07 |
| W304 | `w304_two_dimensions.py` | `W304TwoDimensions` | 二維也會回家，但慢得驚人：第一次回家、$\pi/\ln N$、$10^{136}$ 步 | §3.5，pp. 58–60 | 2:09 |
| W305 | `w305_dimension_divide.py` | `W305DimensionDivide` | 維度的分水嶺：$q^{d-1}/q^2$、$K_d$、$G(z;0)$ 的奇異部分、練習 3.5–3.6 | §3.6，pp. 58–62 | 2:24 |
| W306 | `w306_watson.py` | `W306Watson` | 三種立方格子的逃脫機率（Watson 積分）、回家 $m$ 次（練習 3.8） | §3.7，pp. 62–63 | 2:09 |
| W307 | `w307_sites_visited.py` | `W307SitesVisited` | 走過幾個點：新點的機率＝還沒回家的機率（倒著走、Fig. 3.3） | §3.8，pp. 63–65 | 2:11 |
| W308 | `w308_footprint.py` | `W308Footprint` | 足跡大小 $\sqrt{8n/\pi}$、$\pi n/\ln n$、$0.6595n$＋第 3 章回顧（練習 3.9–3.12） | §3.8，pp. 65–68 | 2:14 |

連續播放：`main.py` 的 `Ch3Demo`（W301–W308），長度 17:39。計畫與場景卡：`../notes/plan_ch03.md`。數學：`utils/walk_rec.py`（依賴 `walk_gf`）；第 3 章沒有新的畫圖模組（小元件寫在場景檔的模組層級函式）。

## 3d. 場景（第 4 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W401 | `w401_image_walker.py` | `W401ImageWalker` | 鏡像走路人：16 條路線裡 6 條活著、5 對鏡射配對、(4.1)、一維與 sc 的存活機率 | §4.1.1，pp. 70–73 | 2:22 |
| W402 | `w402_two_walls.py` | `W402TwoWalls` | 兩牆之間：逐步計數 324/1024、特徵值 $2\cos(m\pi/6)$、模態和 162 與 81、鏡像（練習 4.3） | §4.1.2，pp. 73–76 | 2:12 |
| W403 | `w403_long_short_time.py` | `W403LongShortTime` | 短時間像沒牆、長時間 $e^{-\pi^2N/2L^2}$（(4.24) 的筆誤）、牆距加倍慢 4.05 倍、漏失的解釋 | §4.1.3，pp. 76–79 | 2:09 |
| W404 | `w404_sphere_cube.py` | `W404SphereCube` | 球裡的 $j_0$ 模態、權重和＝體積、(4.35)(4.37) 的筆誤、立方體（練習 4.5）、同體積球活得久 | §4.1.4，pp. 79–82 | 2:14 |
| W405 | `w405_steady_tent.py` | `W405SteadyTent` | 穩定態＝Poisson 方程（(4.40)(4.41) 的筆誤）、Fick 定律、帳篷 6/5…、賭徒破產 3/5（練習 4.6） | §4.2–4.2.1，pp. 82–86 | 2:18 |
| W406 | `w406_point_source.py` | `W406PointSource` | 點源的 $1/r$（高斯對時間積分）、sc 格子 $nG\to3/2\pi$、吸收牆的鏡像源與偶極 | §4.2.2–4.2.3，pp. 86–89 | 2:14 |
| W407 | `w407_absorbing_sphere.py` | `W407AbsorbingSphere` | 吸收球：流量 $4\pi D r_0c_\infty$、細胞最大半徑（假設值）、命中機率 $r_0/R$、電容 | §4.2.4–4.2.5，pp. 89–91 | 2:05 |
| W408 | `w408_receptors.py` | `W408Receptors` | 布滿受器的球（Berg–Purcell，書外）、吸收邊界 (4:S-1)–(4:S-3)、Catalan、第 4 章回顧 | §4.2.6、§4.3，pp. 91–94 | 2:04 |

連續播放：`main.py` 的 `Ch4Demo`（W401–W408），長度 17:49。計畫與場景卡：`../notes/plan_ch04.md`。數學：`utils/walk_bc.py`（依賴 `walk_gf`、`walk_rec`）；第 4 章沒有新的畫圖模組。

## 3e. 場景（第 5 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W501 | `w501_biased_walk.py` | `W501BiasedWalk` | 有偏向的漫步：不公正硬幣 $p=0.6$、移動的高斯、寬度 $\propto\sqrt{p(1-p)}$、三維的「先吹再走」（(5.3) 的筆誤、練習 5.1） | §5.1–5.1.1，pp. 95–98 | 2:22 |
| W502 | `w502_persistent_walk.py` | `W502PersistentWalk` | 有慣性的漫步：$L$、$R$ 兩個狀態、逐步表、轉移矩陣與特徵值、生成函數 | §5.2，pp. 98–100 | 2:12 |
| W503 | `w503_persistence_length.py` | `W503PersistenceLength` | 持續長度 $l_p=\frac{Q}{1-Q}l$：1000 步的變異數 2996 vs 3000、$Q=1$ 與 $Q=0$、二階遞迴（(5.37) 的筆誤、練習 5.2） | §5.2–5.2.1，pp. 100–103 | 2:14 |
| W504 | `w504_telegrapher.py` | `W504Telegrapher` | 電報方程：兩個尖峰 $e^{-\gamma t}/2$＋拖在後面的分布（Fig. 5.3）、因果性、大 $t$ 的高斯（(5.61) 的筆誤、練習 5.3–5.6） | §5.2.2，pp. 103–109 | 2:17 |
| W505 | `w505_run_and_tumble.py` | `W505RunAndTumble` | 二維、三維的「跑一跑、轉個彎」：算子 $\Gamma$、(5.89)、三維 $P_3$（Fig. 5.6）、10 萬隻模擬、前緣（練習 5.7） | §5.2.3–5.2.5，pp. 109–118 | 2:26 |
| W506 | `w506_waiting_time.py` | `W506WaitingTime` | 連續時間的漫步：等待時間 $a(t)$、Poisson、兩段等待、$A=1-\lambda B$、Bessel 解（練習 5.8） | §5.3，pp. 118–122 | 2:13 |
| W507 | `w507_ctrw_gauss.py` | `W507CtrwGauss` | $B/(1-\chi A)$ 與高斯極限：變異數 $\propto t$（書上寫寬度）、會消失的走路人（練習 5.9）＋第 5 章回顧 | §5.3，pp. 122–125 | 2:08 |

連續播放：`main.py` 的 `Ch5Demo`（W501–W507），長度 16:02。計畫與場景卡：`../notes/plan_ch05.md`。數學：`utils/walk_var.py`（依賴 numpy、scipy.special 與 `walk_bc`）；第 5 章沒有新的畫圖模組。

## 3f. 場景（第 6 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W601 | `w601_shape_cloud.py` | `W601ShapeCloud` | 漫步長什麼樣子？1000 條 100 步疊成圓雲（平均張量對角 8.35、8.59）、單獨一條卻是長條（6 條各帶橢圓）、4 步「東東北西」的重心與 $R_g^2=\frac45$ | §6.1–6.1.1，pp. 127–129 | 2:02 |
| W602 | `w602_gyration_tensor.py` | `W602GyrationTensor` | 迴轉半徑張量：定義、$\vec n\cdot\mathbf T\cdot\vec n$ 隨方向變、種子 643 的 1000 步（$R_1^2=258.8$、$R_2^2=31.6$、主軸 28.4°）、轉 30° 特徵值不變 | §6.1.2–6.1.3，pp. 129–133 | 2:03 |
| W603 | `w603_asphericity.py` | `W603Asphericity` | 非球度 $A_d=\frac{4+2d}{4+5d}$：三維 10/19、2 萬條的直方圖（模擬 0.5231）、個別的平均 $A'=0.3939$（(6.16) 的筆誤） | §6.1.4，pp. 133–135 | 2:04 |
| W604 | `w604_infinite_dim.py` | `W604InfiniteDim` | 無限多維的樓梯：練習 6.2 的 0.6282、$N=3$ 的特徵值 0.4268、0.1250、0.0732、$R_k^2=\frac{N+1}{\pi^2k^2}$、非球度 2/5、9 : 2.25 : 1 vs 11.95 : 2.69 : 1（(6.54) 的筆誤） | §6.2，pp. 135–138 | 2:01 |
| W605 | `w605_kramers_identity.py` | `W605KramersIdentity` | Kramers 恆等式：$a_{\alpha\beta}$ 矩陣、$\sin k\pi x$ 模態、高斯步伐與 Wick 配對、$\langle T_{11}^2\rangle$ 等的模擬 vs 精確（(6.22)、(6.35b) 的筆誤） | §6.2.1，pp. 138–142 | 2:12 |
| W606 | `w606_one_over_d.py` | `W606OneOverD` | 1/d 展開：尖峰變寬、resolvent、位移 $\frac{3}{4d}$、寬度 $\sqrt{2/d}$、$d=100,20,3$ 的特徵值直方圖（(6.52) 的筆誤、練習 6.7） | §6.2.2–6.2.5，pp. 142–152 | 2:09 |
| W607 | `w607_fixman_saw.py` | `W607FixmanSaw` | 主迴轉半徑的分布：gamma 近似、Fixman 的精確分布（平均 1/6、變異數 1/135）、自迴避的 (6.77)＋第 6 章回顧（(6.76) 的筆誤） | §6.2.6–6.3，pp. 152–154 | 2:05 |

連續播放：`main.py` 的 `Ch6Demo`（W601–W607），長度 14:46。計畫與場景卡：`../notes/plan_ch06.md`。數學：`utils/walk_shape.py`（依賴 numpy、scipy.integrate）；第 6 章沒有新的畫圖模組。

## 3g. 場景（第 7 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W701 | `w701_path_sum.py` | `W701PathSum` | 漫步就是路徑的和：切成兩段、對中間點求和（一維 $1\cdot1+2\cdot2+1\cdot1=6$、二維 36）、高斯摺積的變異數相加（練習 7.1、(7.6) 的筆誤）、路徑積分 (7.8) | §7.1，pp. 167–171 | 2:19 |
| W702 | `w702_self_avoiding.py` | `W702SelfAvoiding` | 自迴避漫步：三種定義、權重 $e^{-v}$、二維 $N=2,3$ 的 12、36 條、精確列舉到 $c_{14}=2374444$（存活比例 0.0088）、$\langle R^2\rangle_{14}=42.79$ vs 14 | §7.2–7.2.1，pp. 171–172 | 2:06 |
| W703 | `w703_one_minus_f.py` | `W703OneMinusF` | $1-f$ 展開：$f(x)$、一階的圖、二階 3 種、三階 15 種（練習 7.2）、$N=4$ 逐階 256、28、124、96、100、$N=8$ 的部分和忽正忽負（(7.11) 的筆誤） | §7.2.2，pp. 172–176 | 2:12 |
| W704 | `w704_dyson.py` | `W704Dyson` | 生成函數讓摺積變乘積、費曼規則、幾何級數 → 4（照書上的號得 4/3）、不可約的圖 1、2、10、74 → Dyson 方程（(7.35) 的筆誤、Fig. 7.3(a) 應為 (b)） | §7.2.3–7.2.4，pp. 176–184 | 2:23 |
| W705 | `w705_mean_field.py` | `W705MeanField` | 微擾論失效與平均場：冪次表 $1-l(4-d)/2$、一維 $\Sigma_1=1/\sqrt{1-4z^2}$、$d=4$ 是分界、$C_{MF}\to e^{-fm^2/\mathcal N}$（練習 7.4、7.5；(7.59) 的筆誤） | §7.2.5–7.2.6，pp. 184–191 | 2:20 |
| W706 | `w706_flory.py` | `W706Flory` | Flory 的拔河：$R^*=177.83$、$\nu=3/(d+2)$、二維列舉的局部斜率 1.43→1.46 vs $2\nu=3/2$＋第 7 章回顧 | §7.2.7，pp. 191–192 | 2:11 |

連續播放：`main.py` 的 `Ch7Demo`（W701–W706），長度 13:40。計畫與場景卡：`../notes/plan_ch07.md`。數學：`utils/walk_saw.py`（依賴 numpy、scipy.integrate／optimize／special）；第 7 章沒有新的畫圖模組。

## 3h. 場景（第 8 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W801 | `w801_non_gaussian.py` | `W801NonGaussian` | 長尾巴的步伐：二階動差發散（9、99、999）、$1-\chi\approx B|q|^{2-\delta}$（$B=\pi/4$）、寬度 $\propto N$（1822 對 28.06）、尾巴 $r^{-5+\delta}$、換變數 $w=kr$（練習 8.1；p. 195 的筆誤） | §8.1–8.1.1，pp. 193–195 | 2:06 |
| W802 | `w802_scaling_form.py` | `W802ScalingForm` | 標度形式：一維 $g(z;0)=5,10,50$、塌縮到 $\frac{2}{4+x^2}$、(8.9)(8.10) 的 $\nu,\eta$、實空間 Yukawa 形式 $\xi\propto(z_c-z)^{-1/2}$（2.0412、4.0825；$z_c$ 的筆誤） | §8.1.2 前半，pp. 195–197 | 2:03 |
| W803 | `w803_steepest_descent.py` | `W803SteepestDescent` | 最陡下降：(8.18)(8.19)、$C\sim e^{-(r/N^\nu)^{1/(1-\nu)}}$、一維 $N=20$ 精確 0.2098 vs 更正 0.2179 vs 書上 0.0475、Fisher 指數 2、4、$\gamma=(2-\eta)\nu$、比值法 $\hat\mu\to2.638$（練習 8.2–8.5；(8.18)(8.20) 的筆誤） | §8.1.2 後半，pp. 197–199 | 2:10 |
| W804 | `w804_rescaling_flory.py` | `W804RescalingFlory` | 放大鏡：重新標度 (8.22)、自迴避項乘 $l_2^{(d-4)/2}$、上臨界維度 4、兩項同比例 → $\nu=3/(d+2)$（和 W706 同一個答案）＋第 8 章回顧（練習 8.6–8.8；p. 200 與 Ex. 8.6 的筆誤） | §8.1.3，pp. 199–202 | 2:04 |

連續播放：`main.py` 的 `Ch8Demo`（W801–W804），長度 8:28。計畫與場景卡：`../notes/plan_ch08.md`。數學：`utils/walk_scale.py`（依賴 numpy、scipy 與 `walk_gf`、`walk_saw`；場景用 `from utils import walk_scale as WS`）；第 8 章沒有新的畫圖模組。

## 3i. 場景（第 9 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W901 | `w901_critical_exponent.py` | `W901CriticalExponent` | 漫步的指數就是臨界指數：Flory $6/(d+2)$ vs 已知 $2\nu$（三維差 2%）、二維列舉的 $2\nu_{\rm eff}$ 1.4442→1.4556、平均場 vs Onsager 的 $M_s$（$\beta=\frac12$ 對 $\frac18$）、de Gennes 的對照表 | §9.1–9.2，pp. 203–205 | 2:07 |
| W902 | `w902_broken_symmetry.py` | `W902BrokenSymmetry` | O(n) 模型與對稱破缺：8 個平面自旋一起轉（$E=-6.3851$、$|\vec M|=6.7943$ 不變）、4 自旋的 16 個組態（$Z=121.2329$、$\langle|m|\rangle=0.9337$ 對 0.4871）、Metropolis 翻轉等待 133、608、2151 sweep、外場與相圖 | §9.2.1–9.2.2，pp. 205–211 | 2:25 |
| W903 | `w903_scaling_laws.py` | `W903ScalingLaws` | 臨界指數與標度律：Table 9.1、一維 Ising 鏈 $\xi=3.6719$、$\chi=e^{2K}=7.3891$、平均場／二維 Ising／二維 $n=0$ 驗算四條標度律、Widom 形式（練習 9.1–9.2；Table 9.1 的 $\alpha$ 與第二個關聯長度兩處筆誤） | §9.2.2–9.2.3，pp. 209–215 | 2:24 |
| W904 | `w904_landau_wilson.py` | `W904LandauWilson` | Ginzburg–Landau–Wilson：可伸縮的自旋 (9.28)、8 點環的 Fourier 核對（練習 9.3）、Landau 雙井 $m^*=1.0467$、練習 9.4（$\sqrt3$、$\mathcal M(1)=2.1038$、$\beta=\frac12$、$\Delta=\frac32$）、(9.34) 的有限體積積分 | §9.3，pp. 215–217 | 2:21 |
| W905 | `w905_walk_from_magnet.py` | `W905WalkFromMagnet` | 從磁鐵讀出漫步：(9.36)(9.37)、$c_N\langle R^2\rangle_N=4,32,164,\dots$、$(z_c-z)^{-s}$ 的係數漸近、$\langle R^2\rangle_N\propto N^{2\nu}$、前因子 $\mu^{3/2}=4.2850$（(9.38)(9.39)、(9.41) 的筆誤） | §9.4，pp. 218–219 | 2:03 |
| W906 | `w906_n_to_zero.py` | `W906NToZero` | $n\to0$ 的魔術：輔助場 (9.46)、高斯積分 5.0627 與 Wick、Trotter 誤差 $\propto\epsilon^2$（練習 9.6）、每次相交乘 $e^{-u/2}=0.6065$、$N=8$ 的 65536→5916、$2\times3$ 梯子 $n\to0$ 得 $K+K^3+K^5$＋第 9 章回顧（(9.70)、(9.75) 的筆誤） | §9.5–9.6，pp. 219–231 | 2:25 |

連續播放：`main.py` 的 `Ch9Demo`（W901–W906），長度 13:53。計畫與場景卡：`../notes/plan_ch09.md`。數學：`utils/walk_crit.py`（依賴 numpy、scipy 與 `walk_saw`、`walk_scale`；場景用 `from utils import walk_crit as ...`）；第 9 章沒有新的畫圖模組。

## 3j. 場景（第 10 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W1001 | `w1001_confined_polymer.py` | `W1001ConfinedPolymer` | 關在盒子裡的高分子：和單體池平衡（Fig. 10.1）、$G(z)=\sum_NC(N)z^N$、$4\times4$ 環面：$c_N$ 從 $N=4$ 起比無限平面少（96 對 100）、走滿 16 格的 5712 條、三種 $z$ 的 $\langle N\rangle=1.44$、4.89、11.74（(10.1) 的筆誤） | §10.1，pp. 233–236 | 2:18 |
| W1002 | `w1002_mean_field.py` | `W1002MeanField` | O(n) 的平均場：插入 $\int d\lambda\,\delta(\lambda-|\vec M|^2)$、$n\to0$ 收斂到 1.8682、對 $S$ 的積分是 $2\pi i$（不是 $i\pi$）、$\Gamma_{MF}$ 的兩種漸近（$1/r$；$r^2/4u=0.0625$）、極限的順序、和 (7.57) 對照（(10.15) 與註 3 的筆誤） | §10.2，pp. 236–239 | 2:25 |
| W1003 | `w1003_spin_waves.py` | `W1003SpinWaves` | 自旋波：$\vec S=\vec M+\vec\sigma$、牆邊的 $\tanh$（癒合長度 1.4142，練習 10.1）、盒子裡的模態、$X_{\vec q}$ 的特徵值 $\frac12$（橫向 $n-1$ 重）與 $\frac32$（縱向）、逆矩陣（練習 10.3–10.4）、$n\to0$ 時 $-1$ 個橫向模態、$(\det X)^{-1/2}=1/\sqrt3$（(10.32)、(10.33) 的筆誤） | §10.3，pp. 239–243 | 2:14 |
| W1004 | `w1004_correlation_hole.py` | `W1004CorrelationHole` | 相關洞：(10.37)–(10.45) 的結構、$r\gg0$ 時盒子沒影響、$D=\sum_{q\ne0}\frac{2|r|\psi\psi^*}{q^2(q^2+2|r|)}$、$64^3$ 週期格子的洞（頭尾相鄰少 11%）、Coulomb 形 $1/(4\pi\Delta)$ 與 Ewald 常數、普通漫步離起點 $R$ 的步數 $\propto1/R$（erfc 形；練習 10.5；(10.46)、(10.83) 的筆誤） | §10.3 後半、§10.4，pp. 243–246、251–253 | 2:10 |
| W1005 | `w1005_redner_reynolds.py` | `W1005RednerReynolds` | Redner–Reynolds：$\Omega$、$\langle N\rangle$、密度 $\rho$；環面 $3\times3$、$4\times4$、$5\times5$ 的 $\rho(z)$；平均場 $\rho=|r|/2u$ → $\beta'=1$、$\gamma'=\alpha'=0$、$\nu'=\frac12$、上臨界維度 4；Table 10.1（二維 $\frac12,\frac12$；三維 0.764、0.236）＋第 10 章回顧（(10.53)(10.55)、(10.75) 的筆誤） | §10.3.1，pp. 246–251 | 2:06 |

連續播放：`main.py` 的 `Ch10Demo`（W1001–W1005），長度 11:19。計畫與場景卡：`../notes/plan_ch10.md`。數學：`utils/walk_cond.py`（依賴 numpy、scipy 與 `walk_saw`、`walk_crit` 等；場景用 `from utils import walk_cond as ...`）；第 10 章沒有新的畫圖模組。

## 3k. 場景（第 11 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W1101 | `w1101_fractals.py` | `W1101Fractals` | 碎形：點、直線、圓弧放大；Koch 曲線（$4^k$ 段、長度 $(4/3)^k\to\infty$）、Cantor 集（$(2/3)^k\to0$）；質量－半徑 $\ln4/\ln3=1.2619$、$\ln2/\ln3=0.6309$、盒子計數斜率 1.2362、複本公式（Sierpinski 1.5850）；普通漫步斜率 2.0258、自迴避 $d_f=1/\nu=\frac43$（練習 11.1–11.3） | §11.1–11.1.1，pp. 255–261 | 2:32 |
| W1102 | `w1102_passive_rg.py` | `W1102PassiveRG` | 被動的重整化：放大 $b$ 倍，$m\to mb^{d/4}$、$r\to rb^{d/2}$、$h\to hb^{3d/4}$ → $m=r^{1/2}f(h/r^{3/2})$、三條曲線塌縮成 $y+y^3=x$、$\gamma=1$；Gaussian 模型 $m=r^{(d-2)/4}F(h/r^{(d+2)/4})$、練習 11.4（$d=3$：$\beta=\frac14$、$\delta=5$；$d=4$ 與 Landau 相同；(11.16) 後的筆誤） | §11.1.2，pp. 261–264 | 2:04 |
| W1103 | `w1103_gaussian_decimation.py` | `W1103GaussianDecimation` | 一維 Gaussian 鏈的抽格：完成平方（數值積分 1.224819）、$r''=r(r+4a)/a$、$h''$、自旋乘 $\sqrt{(r+2a)/a}$；蛛網圖 0.05→0.2025→0.8510→4.1282、$r^{(n)}\approx4^nr$、$h$ 每次乘 $2^{3/2}$ → $\gamma=1$（練習 11.5–11.6；(11.31)、Fig. 11.6–11.7、(11.46)(11.47) 的筆誤） | §11.2 前半，pp. 264–269 | 2:06 |
| W1104 | `w1104_ising_decimation.py` | `W1104IsingDecimation` | 一維 Ising 的抽格：$A=2\sqrt{\cosh2K}=3.8793$、$K'=\frac12\ln\cosh2K=0.6625$、$K_n$ 一路趨向 0、加外場（練習 11.8）、4 自旋環 $Z=121.2329$（暴力＝轉移矩陣＝抽格連乘）、$x'=\frac{2x}{1+x^2}$ 的蛛網圖與兩個不動點、$\xi/e^{2K}\to\frac12$（練習 11.7、11.9；(11.53)、(11.63) 的筆誤） | §11.2 例題，pp. 269–272 | 2:07 |
| W1105 | `w1105_momentum_shell.py` | `W1105MomentumShell` | 動量殼與連續重整化：布里淵區縮小一半、積掉殼層、$\alpha(q)=a(1-\cos q)$、$r\to4r$、$h\to2^{(d+2)/2}h$；$\frac{dr}{d\ell}=2r$、$\frac{dh}{d\ell}=\frac{d+2}2h$ 的流線、(11.24) 代入 (11.87) 驗證（殘差 $<10^{-7}$；(11.66)(11.67) 的筆誤） | §11.2.1–11.2.2，pp. 272–277 | 2:07 |
| W1106 | `w1106_fixed_points.py` | `W1106FixedPoints` | 不動點與臨界指數：區塊自旋、臨界曲面（Fig. 11.9 示意）、$2\times2$ 線性化矩陣（示意，特徵值 2、$-3$：相關／無關）、$\nu=1/\Lambda_r$、四個指數公式 → Gaussian（$d=3,4$）、二維 Ising、二維自迴避的指數全部重現第 9 章 W903 的表＋第 11 章回顧 | §11.3，pp. 277–284 | 2:12 |

連續播放：`main.py` 的 `Ch11Demo`（W1101–W1106），長度 13:15。計畫與場景卡：`../notes/plan_ch11.md`。數學：`utils/walk_rg.py`（依賴 numpy、scipy 與前幾章的模組；場景用 `from utils import walk_rg as ...`）；第 11 章沒有新的畫圖模組。

## 3l. 場景（第 12 章，`scenes/`）

| 代號 | 檔名 | 類別 | 內容 | 書上 | 長度 |
|---|---|---|---|---|---|
| W1201 | `w1201_quartic_shell.py` | `W1201QuarticShell` | 加上四次項：GLW 模型、倒晶格、球殼 $e^{-\ell-\delta\ell}<q<e^{-\ell}$、殼層比 $1-e^{-d\delta\ell}$（0.2592 對 0.3000）、四個步驟；兩變數玩具 $r'=0.3000$、$u'=0.0075$、重新標度的乘數（(12.32)(12.33) 的筆誤） | §12.1–12.2.4，pp. 285–289 | 2:07 |
| W1202 | `w1202_linked_clusters.py` | `W1202LinkedClusters` | 配對遊戲：頂點與殼層手、$15=1+6+8$、$f_3=6$、5 個片段的 10 張圖（$\frac{5!}{3!2!}$）、(12:S-4) 的列舉核對、連結叢集＝累積量（四種算法都是 0.9129；p. 320 的筆誤） | §12.2.4、§12.7，pp. 289–292、317–322 | 2:00 |
| W1203 | `w1203_one_loop_flow.py` | `W1203OneLoopFlow` | 一迴圈與 $(r,u)$ 的流：兩張圖、動量守恆、$\frac{dr}{d\ell}=2r+\frac{12u}{1+r}$、$\frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}$、Gaussian 點與 Wilson–Fisher 點 $(-\frac17,\frac1{49})$、臨界線都流向 WF 點（Fig. 12.16(b)、(12.36)(12.38)(12.39) 的筆誤） | §12.2.4–12.2.5，pp. 292–297 | 2:04 |
| W1204 | `w1204_epsilon_nu.py` | `W1204EpsilonNu` | ε 展開求 $\nu$：$u(\ell)$ 的閉式與 $u^*=\frac1{36}$（練習 12.2–12.3）、$t$ 的流、$\nu=\frac1{2-\epsilon/3}\approx\frac12+\frac{\epsilon}{12}$、$\epsilon=1$ 得 $\frac7{12}=0.5833$、三種讀法的表、已知 0.630 只引用（(12.48)、練習 12.4、「0.61」的筆誤） | §12.2.5–12.2.6，pp. 296–300 | 2:03 |
| W1205 | `w1205_on_saw.py` | `W1205OnSaw` | O(n) 與自迴避：$2n+4$、$4(n+8)$、$\nu=\frac12+\frac{(n+2)\epsilon}{4(n+8)}$、$n=0$ 得 0.5625、$n=0..3$ 的表對已知值、Flory 的展開 $\frac12+\frac{\epsilon}{12}$、$N=1000$ 的大小 31.62／48.70／58.08、ε 展開的四個隱憂（練習 12.5–12.8） | §12.3–12.4，pp. 300–304 | 2:08 |
| W1206 | `w1206_diagrammatic_recap.py` | `W1206DiagrammaticRecap` | 圖解法的重整化群：質量算子與 Dyson 形、$g_0(q_0)$ 的方程＝(12.39)、冪次數算（數值 $-1$、$-\frac12$；(12.89) 的筆誤）；第 12 章回顧＋**全書回顧**（12 列，兩個畫面） | §12.5–12.6，pp. 304–317 | 2:07 |

連續播放：`main.py` 的 `Ch12Demo`（W1201–W1206），長度 12:38。計畫與場景卡：`../notes/plan_ch12.md`。數學：`utils/walk_eps.py`（依賴 numpy、scipy 與前幾章的模組；場景用 `from utils import walk_eps as ...`）；第 12 章沒有新的畫圖模組。全書 12 章到此全部完成。

## 4. 程式結構

```text
codes/
├── README.md             本檔
├── custom_config.yml     ManimGL 設定（只能寫英文）
├── main.py               每個場景的同名子類別＋Ch1Demo、Ch2Demo、Ch3Demo、Ch4Demo、Ch5Demo、Ch6Demo、Ch7Demo、Ch8Demo、Ch9Demo、Ch10Demo、Ch11Demo、Ch12Demo（play_sequence 依序借用各場景的 construct）
├── render_all.py / .bat  批次輸出（.bat 必須純 ASCII＋CRLF）
├── check_env.py          環境檢查
├── scenes/               w101_…py … w112_…py、w201_…py … w212_…py、w301_…py … w308_…py（一個檔案一幕；輔助函式都是模組層級）
├── utils/
│   ├── walk1d.py         第 1 章的數學（純 Python；check() 核對書上的數字與每張場景卡的數字）
│   ├── walk_draw.py      第 1 章的畫圖元件（ManimGL）：硬幣、數線、走法圖、扇形、長條、Pascal、Δ 格子、表格、筆誤框、圍道
│   ├── walk_gf.py        第 2 章的數學（純 Python，不依賴 walk1d）：遞迴計數、結構函數、布里淵區積分、高斯與普適性、Cauchy、擴散、補充 1–3
│   ├── lattice_draw.py   第 2 章的畫圖元件：格線、數字格、鄰居、足跡、熱圖、布里淵區、單位根、三維晶格、分支切割
│   ├── walk_rec.py       第 3 章的數學（純 Python，依賴 walk_gf）：到訪次數的生成函數、還沒回家的機率、Watson 積分、足跡、固定種子的模擬
│   ├── walk_bc.py        第 4 章的數學（純 Python，依賴 walk_gf、walk_rec）：反射原理、兩牆之間、球與立方體、穩定態、點源與鏡像、吸收球、Berg–Purcell
│   ├── walk_var.py       第 5 章的數學（純 Python，依賴 numpy、scipy.special、walk_bc）：偏向漫步、有慣性的漫步、電報方程、二三維的跑與轉、連續時間漫步
│   ├── walk_shape.py     第 6 章的數學（純 Python，依賴 numpy、scipy.integrate）：迴轉半徑張量、Kramers 恆等式、非球度、無限多維的樓梯、1/d 展開、Fixman 分布
│   ├── walk_saw.py       第 7 章的數學（純 Python，依賴 numpy、scipy）：路徑和的切段、自迴避漫步的精確列舉、1−f 展開、弧線圖與 Dyson 方程、冪次數算、平均場、Flory 指數
│   ├── walk_scale.py     第 8 章的數學（純 Python，依賴 numpy、scipy、walk_gf、walk_saw）：長尾巴步伐、標度形式與塌縮、最陡下降與 Fisher 指數、比值法、重新標度與上臨界維度
│   ├── walk_crit.py      第 9 章的數學（純 Python，依賴 numpy、scipy、walk_saw、walk_scale）：自發磁化、小 Ising 系統、一維 Ising 鏈、指數表與標度律、Landau 自由能、高斯積分與 Wick、Trotter 分解
│   ├── walk_cond.py      第 10 章的數學（純 Python，依賴 numpy、scipy 與前幾章的模組）：環面上的自迴避漫步、O(n) 平均場、自旋波的橫向與縱向模態、相關洞、Redner–Reynolds 的指數
│   ├── walk_rg.py        第 11 章的數學（純 Python，依賴 numpy、scipy 與前幾章的模組）：Koch 曲線與 Cantor 集、碎形維度、被動重整化、Gaussian 與 Ising 的抽格、動量殼與連續重整化、不動點與指數
│   └── walk_eps.py       第 12 章的數學（純 Python，依賴 numpy、scipy 與前幾章的模組）：積掉殼層、連結叢集的組合學、(r,u) 的流與 Wilson–Fisher 不動點、ε 展開的 ν、O(n) 與自迴避、冪次數算
└── renders/              影片（不進 git）
```

`walk1d.py`、`walk_draw.py`、`walk_gf.py`、`lattice_draw.py`、`walk_rec.py` 的複本在 combkit（`combkit.walk1d`、`combkit.walk_gf`、`combkit.walk_rec`、
`combkit.draw.walk_draw`、`combkit.draw.lattice_draw`），用 combkit 的 `tools/sync_from_walk017.py` 同步；
改這裡的正本，再到 combkit 執行同步、升版、打包。

## 5. 每個數字怎麼核對

- 場景：每幕的 `check_math()` 用 `utils.walk1d`（第 1 章）、`utils.walk_gf`（第 2 章）、`utils.walk_rec`（第 3 章）、`utils.walk_bc`（第 4 章）、`utils.walk_var`（第 5 章）、`utils.walk_shape`（第 6 章）、`utils.walk_saw`（第 7 章）、`utils.walk_scale`（第 8 章）、`utils.walk_crit`（第 9 章）、`utils.walk_cond`（第 10 章）、`utils.walk_rg`（第 11 章）或 `utils.walk_eps`（第 12 章）重算畫面上的每個數並 assert；版面用 `assert_in_frame`／`assert_apart`。
- 數學模組：`python utils/walk1d.py`（公式 vs 暴力列舉、三種動差算法、書上的 15.1044／15.0961 等、`../notes/plan_ch01.md` 場景卡上的全部數字）；
  `python utils/walk_gf.py`（遞迴＝兩枚硬幣＝傅立葉積分、書上的數字、筆誤的更正、`../notes/plan_ch02.md` 場景卡上的全部數字）；
  `python utils/walk_rec.py`（生成函數＝暴力列舉＝逐步遞迴、Watson 積分的 Gamma 閉式＝Bessel 積分＝級數、模擬的固定種子、`../notes/plan_ch03.md` 場景卡上的全部數字）。
- 講義：`../notes/ch1/make_figs.py` 的 `check_numbers_ch1()`、`../notes/ch2/make_figs.py` 的 `check_numbers_ch2()`、`../notes/ch3/make_figs.py` 的 `check_numbers_ch3()`。
- 書上的筆誤：`../ref/md/ch01_formulas.md`（13 處）、`../ref/md/ch02_formulas.md`（16 處）、`../ref/md/ch03_formulas.md`（12 處，含 2 條〔注〕），動畫裡用紅色「書上筆誤」框標出。
