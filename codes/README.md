# walk017：《Elements of the Random Walk》導讀動畫（ManimGL）

J. Rudnick & G. Gaspari, *Elements of the Random Walk: An Introduction for Advanced Students and Researchers*（Cambridge University Press 2004）
的繁體中文導讀動畫與講義。目前完成**第 1–2 章**：第 1 章（Introduction to techniques）12 幕動畫 W101–W112＋連續播放 `Ch1Demo`；
第 2 章（Generating functions I）12 幕動畫 W201–W212＋連續播放 `Ch2Demo`；導讀講義 `../notes/walk_導讀.pdf`（兩章）。

**路徑的寫法：** 本檔裡沒有寫磁碟機的相對路徑，一律相對於本資料夾 `C:\mypython\manim\projects\walk017\codes`；
寫成 `../notes/…`、`../ref/…` 的在上一層 `C:\mypython\manim\projects\walk017`。權威說明是上一層的 `交接說明_HANDOFF.md`。

## 1. 環境與依賴

- ManimGL **1.7.2**（3Blue1Brown 版，commit fafa083，`manimgl` 指令、`from manimlib import *`），**不是** Manim Community；不要升級或重裝。
  全機環境見 `C:\mypython\manim\環境與交接_ENVIRONMENT.md`、`C:\mypython\manim\新計畫起步指南_NEW_PROJECT_GUIDE.md`。
- 共用工具包 **combkit 1.12.0 以上**（`C:\mypython\manim\projects\combkit`）：場景基底 `ChuScene`、中文字 `zh()`／`mixed()`、版面檢查 `assert_in_frame`。
  場景檔開頭用 `sys.path.insert` 載入它，不必安裝。
- LaTeX：MiKTeX（`latex`＋`dvisvgm`＋`cm-super`；講義用 `xelatex`＋`xeCJK`、`tcolorbox`、`tikz`）。字型：微軟正黑體、標楷體。
- 檢查：`python check_env.py`（加 `--render` 會實際渲染一小段）。

## 2. 指令（都在本資料夾執行）

| 做什麼 | 指令 |
|---|---|
| 預覽一幕 | `manimgl scenes/w101_coin_walk.py W101CoinWalk`（Ctrl+Q 結束、滾輪縮放、按住 f 平移、r 重設） |
| 輸出一幕 1080p | `manimgl scenes/w101_coin_walk.py W101CoinWalk -w --hd`（低畫質把 `--hd` 換 `-l`） |
| 一次輸出全部 | `python render_all.py`（兩章 24 幕）；`--ch1`／`--ch2` 只做一章；`--ch1-demo`／`--ch2-demo` 加連續播放；`-l` 480p；`W105 210` 只做指定的幕；已存在的影片預設**跳過**，要重做加 `--force` |
| 按兩下的選單 | `render_all.bat`（會自動啟用 conda） |
| 連續播放 | `manimgl main.py Ch1Demo -w --hd`、`manimgl main.py Ch2Demo -w --hd` |
| 數學自我測試 | `python utils/walk1d.py`（第 1 章）、`python utils/walk_gf.py`（第 2 章） |
| 環境檢查 | `python check_env.py`、`python check_env.py --render` |
| 試渲染＋總覽圖 | 在 `C:\mypython\manim\projects\combkit` 執行 `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，看 `runs/<類別>/sheet.png` |
| 講義 | 在 `../notes` 執行 `python ch1/make_figs.py` 與 `python ch2/make_figs.py`，再 `xelatex walk_導讀.tex` 兩次 |
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

## 4. 程式結構

```text
codes/
├── README.md             本檔
├── custom_config.yml     ManimGL 設定（只能寫英文）
├── main.py               每個場景的同名子類別＋Ch1Demo、Ch2Demo（play_sequence 依序借用各場景的 construct）
├── render_all.py / .bat  批次輸出（.bat 必須純 ASCII＋CRLF）
├── check_env.py          環境檢查
├── scenes/               w101_…py … w112_…py、w201_…py … w212_…py（一個檔案一幕；輔助函式都是模組層級）
├── utils/
│   ├── walk1d.py         第 1 章的數學（純 Python；check() 核對書上的數字與每張場景卡的數字）
│   ├── walk_draw.py      第 1 章的畫圖元件（ManimGL）：硬幣、數線、走法圖、扇形、長條、Pascal、Δ 格子、表格、筆誤框、圍道
│   ├── walk_gf.py        第 2 章的數學（純 Python，不依賴 walk1d）：遞迴計數、結構函數、布里淵區積分、高斯與普適性、Cauchy、擴散、補充 1–3
│   └── lattice_draw.py   第 2 章的畫圖元件：格線、數字格、鄰居、足跡、熱圖、布里淵區、單位根、三維晶格、分支切割
└── renders/              影片（不進 git）
```

`walk1d.py`、`walk_draw.py`、`walk_gf.py`、`lattice_draw.py` 的複本在 combkit（`combkit.walk1d`、`combkit.walk_gf`、`combkit.draw.walk_draw`、
`combkit.draw.lattice_draw`），用 combkit 的 `tools/sync_from_walk017.py` 同步；
改這裡的正本，再到 combkit 執行同步、升版、打包。

## 5. 每個數字怎麼核對

- 場景：每幕的 `check_math()` 用 `utils.walk1d`（第 1 章）或 `utils.walk_gf`（第 2 章）重算畫面上的每個數並 assert；版面用 `assert_in_frame`／`assert_apart`。
- 數學模組：`python utils/walk1d.py`（公式 vs 暴力列舉、三種動差算法、書上的 15.1044／15.0961 等、`../notes/plan_ch01.md` 場景卡上的全部數字）；
  `python utils/walk_gf.py`（遞迴＝兩枚硬幣＝傅立葉積分、書上的數字、筆誤的更正、`../notes/plan_ch02.md` 場景卡上的全部數字）。
- 講義：`../notes/ch1/make_figs.py` 的 `check_numbers_ch1()`、`../notes/ch2/make_figs.py` 的 `check_numbers_ch2()`。
- 書上的筆誤：`../ref/md/ch01_formulas.md`（13 處）、`../ref/md/ch02_formulas.md`（16 處），動畫裡用紅色「書上筆誤」框標出。
