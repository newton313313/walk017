# walk017：《Elements of the Random Walk》導讀動畫（ManimGL）

J. Rudnick & G. Gaspari, *Elements of the Random Walk: An Introduction for Advanced Students and Researchers*（Cambridge University Press 2004）
的繁體中文導讀動畫與講義。目前完成**第 1 章**（Introduction to techniques）：12 幕動畫 W101–W112、連續播放 `Ch1Demo`、
導讀講義 `notes/walk_導讀.pdf`。

**路徑的寫法：** 本檔裡沒有寫磁碟機的相對路徑，一律相對於本資料夾 `C:\mypython\manim\projects\walk017\codes`；
寫成 `../notes/…`、`../ref/…` 的在上一層 `C:\mypython\manim\projects\walk017`。權威說明是上一層的 `交接說明_HANDOFF.md`。

## 1. 環境與依賴

- ManimGL **1.7.2**（3Blue1Brown 版，commit fafa083，`manimgl` 指令、`from manimlib import *`），**不是** Manim Community；不要升級或重裝。
  全機環境見 `C:\mypython\manim\環境與交接_ENVIRONMENT.md`、`C:\mypython\manim\新計畫起步指南_NEW_PROJECT_GUIDE.md`。
- 共用工具包 **combkit 1.11.0 以上**（`C:\mypython\manim\projects\combkit`）：場景基底 `ChuScene`、中文字 `zh()`／`mixed()`、版面檢查 `assert_in_frame`。
  場景檔開頭用 `sys.path.insert` 載入它，不必安裝。
- LaTeX：MiKTeX（`latex`＋`dvisvgm`＋`cm-super`；講義用 `xelatex`＋`xeCJK`、`tcolorbox`、`tikz`）。字型：微軟正黑體、標楷體。
- 檢查：`python check_env.py`（加 `--render` 會實際渲染一小段）。

## 2. 指令（都在本資料夾執行）

| 做什麼 | 指令 |
|---|---|
| 預覽一幕 | `manimgl scenes/w101_coin_walk.py W101CoinWalk`（Ctrl+Q 結束、滾輪縮放、按住 f 平移、r 重設） |
| 輸出一幕 1080p | `manimgl scenes/w101_coin_walk.py W101CoinWalk -w --hd`（低畫質把 `--hd` 換 `-l`） |
| 一次輸出全部 | `python render_all.py`（第 1 章 12 幕）；`--ch1-demo` 加連續播放；`-l` 480p；`W105 110` 只做指定的幕；已存在的影片預設**跳過**，要重做加 `--force` |
| 按兩下的選單 | `render_all.bat`（會自動啟用 conda） |
| 連續播放 | `manimgl main.py Ch1Demo -w --hd` |
| 數學自我測試 | `python utils/walk1d.py` |
| 環境檢查 | `python check_env.py`、`python check_env.py --render` |
| 試渲染＋總覽圖 | 在 `C:\mypython\manim\projects\combkit` 執行 `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，看 `runs/<類別>/sheet.png` |
| 講義 | 在 `../notes` 執行 `python ch1/make_figs.py`，再 `xelatex walk_導讀.tex` 兩次 |
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

## 4. 程式結構

```text
codes/
├── README.md             本檔
├── custom_config.yml     ManimGL 設定（只能寫英文）
├── main.py               每個場景的同名子類別＋Ch1Demo（play_sequence 依序借用各場景的 construct）
├── render_all.py / .bat  批次輸出（.bat 必須純 ASCII＋CRLF）
├── check_env.py          環境檢查
├── scenes/               w101_…py … w112_…py（一個檔案一幕；輔助函式都是模組層級）
├── utils/
│   ├── walk1d.py         第 1 章的數學（純 Python；check() 核對書上的數字與每張場景卡的數字）
│   └── walk_draw.py      第 1 章的畫圖元件（ManimGL）：硬幣、數線、走法圖、扇形、長條、Pascal、Δ 格子、表格、筆誤框、圍道
└── renders/              影片（不進 git）
```

`walk1d.py`、`walk_draw.py` 的複本在 combkit（`combkit.walk1d`、`combkit.draw.walk_draw`），用 combkit 的 `tools/sync_from_walk017.py` 同步；
改這裡的正本，再到 combkit 執行同步、升版、打包。

## 5. 每個數字怎麼核對

- 場景：每幕的 `check_math()` 用 `utils.walk1d` 重算畫面上的每個數並 assert；版面用 `assert_in_frame`／`assert_apart`。
- 數學模組：`python utils/walk1d.py`（公式 vs 暴力列舉、三種動差算法、書上的 15.1044／15.0961 等、`notes/plan_ch01.md` 場景卡上的全部數字）。
- 講義：`../notes/ch1/make_figs.py` 的 `check_numbers_ch1()`。
- 書上的筆誤：`../ref/md/ch01_formulas.md`（逐式校對＋筆誤表），動畫裡用紅色「書上筆誤」框標出。
