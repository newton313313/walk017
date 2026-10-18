# walk017：《Elements of the Random Walk》導讀動畫（ManimGL）

J. Rudnick & G. Gaspari, *Elements of the Random Walk: An Introduction for Advanced Students and Researchers*（Cambridge University Press 2004）
的繁體中文導讀動畫與講義。目前完成**第 1–3 章**：第 1 章（Introduction to techniques）12 幕動畫 W101–W112＋連續播放 `Ch1Demo`；
第 2 章（Generating functions I）12 幕動畫 W201–W212＋連續播放 `Ch2Demo`；第 3 章（Generating functions II：回家、維度、走過幾個點）
8 幕動畫 W301–W308＋連續播放 `Ch3Demo`；第 4 章（邊界條件、穩定態與靜電類比）8 幕動畫 W401–W408＋連續播放 `Ch4Demo`；第 5 章（漫步的變化形：偏向、慣性、連續時間）7 幕動畫 W501–W507＋連續播放 `Ch5Demo`；第 6 章（漫步的形狀：迴轉半徑張量、非球度、1/d 展開）7 幕動畫 W601–W607＋連續播放 `Ch6Demo`；導讀講義 `../notes/walk_導讀.pdf`（六章）。

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
| 連續播放 | `manimgl main.py Ch1Demo -w --hd`、`manimgl main.py Ch2Demo -w --hd`、`manimgl main.py Ch3Demo -w --hd`、`manimgl main.py Ch4Demo -w --hd`、`manimgl main.py Ch5Demo -w --hd`、`manimgl main.py Ch6Demo -w --hd` |
| 數學自我測試 | `python utils/walk1d.py`（第 1 章）、`python utils/walk_gf.py`（第 2 章）、`python utils/walk_rec.py`（第 3 章，約 2 秒）、`python utils/walk_bc.py`（第 4 章，約 1 秒）、`python utils/walk_var.py`（第 5 章，約 1 秒）、`python utils/walk_shape.py`（第 6 章，約 2 秒） |
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

## 4. 程式結構

```text
codes/
├── README.md             本檔
├── custom_config.yml     ManimGL 設定（只能寫英文）
├── main.py               每個場景的同名子類別＋Ch1Demo、Ch2Demo、Ch3Demo、Ch4Demo、Ch5Demo、Ch6Demo（play_sequence 依序借用各場景的 construct）
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
│   └── walk_shape.py     第 6 章的數學（純 Python，依賴 numpy、scipy.integrate）：迴轉半徑張量、Kramers 恆等式、非球度、無限多維的樓梯、1/d 展開、Fixman 分布
└── renders/              影片（不進 git）
```

`walk1d.py`、`walk_draw.py`、`walk_gf.py`、`lattice_draw.py`、`walk_rec.py` 的複本在 combkit（`combkit.walk1d`、`combkit.walk_gf`、`combkit.walk_rec`、
`combkit.draw.walk_draw`、`combkit.draw.lattice_draw`），用 combkit 的 `tools/sync_from_walk017.py` 同步；
改這裡的正本，再到 combkit 執行同步、升版、打包。

## 5. 每個數字怎麼核對

- 場景：每幕的 `check_math()` 用 `utils.walk1d`（第 1 章）、`utils.walk_gf`（第 2 章）、`utils.walk_rec`（第 3 章）、`utils.walk_bc`（第 4 章）、`utils.walk_var`（第 5 章）或 `utils.walk_shape`（第 6 章）重算畫面上的每個數並 assert；版面用 `assert_in_frame`／`assert_apart`。
- 數學模組：`python utils/walk1d.py`（公式 vs 暴力列舉、三種動差算法、書上的 15.1044／15.0961 等、`../notes/plan_ch01.md` 場景卡上的全部數字）；
  `python utils/walk_gf.py`（遞迴＝兩枚硬幣＝傅立葉積分、書上的數字、筆誤的更正、`../notes/plan_ch02.md` 場景卡上的全部數字）；
  `python utils/walk_rec.py`（生成函數＝暴力列舉＝逐步遞迴、Watson 積分的 Gamma 閉式＝Bessel 積分＝級數、模擬的固定種子、`../notes/plan_ch03.md` 場景卡上的全部數字）。
- 講義：`../notes/ch1/make_figs.py` 的 `check_numbers_ch1()`、`../notes/ch2/make_figs.py` 的 `check_numbers_ch2()`、`../notes/ch3/make_figs.py` 的 `check_numbers_ch3()`。
- 書上的筆誤：`../ref/md/ch01_formulas.md`（13 處）、`../ref/md/ch02_formulas.md`（16 處）、`../ref/md/ch03_formulas.md`（12 處，含 2 條〔注〕），動畫裡用紅色「書上筆誤」框標出。
