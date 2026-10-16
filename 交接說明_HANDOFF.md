# walk017 交接說明（權威說明，給繼任者與下一個 AI）

整理日期：2026-10-17。**路徑的寫法：** 本文件裡沒有寫磁碟機的相對路徑，一律相對於 `C:\mypython\manim\projects\walk017`（本文件所在的資料夾）。

## 1. 一句話

把 J. Rudnick & G. Gaspari《Elements of the Random Walk》（Cambridge University Press 2004，12 章）逐章做成繁體中文的
**ManimGL 動畫＋導讀講義**。目前完成：全書分章轉成 md（`ref/md/`）、**第 1 章** 12 幕動畫 W101–W112（＋連續播放 `Ch1Demo`）、
導讀講義 `notes/walk_導讀.pdf` 第 1 章、第 1 章逐式校對與筆誤表。

## 2. 現況

| 項目 | 狀態 |
|---|---|
| 參考文獻 | `ref/ELEMENTS OF THE RANDOM WALK.pdf`（345 頁，有文字層；正文 PDF 頁序＝書頁＋16）。僅供個人研讀，不可散布 |
| md 轉檔 | `ref/convert_to_md.py` → `ref/md/00_front_preface.md`…`13_references_index.md`（14 檔）＋`00_INDEX.md`（每章頁碼與全部小節標題）。數學字型依 /Differences 還原成 Unicode，獨立式排成二維文字圖 |
| 第 1 章校對 | `ref/md/ch01_formulas.md`：全部式子改寫成 LaTeX（依頁面影像核對）＋13 處筆誤表 |
| 計畫 | `notes/plan_ch01.md`（使用者指定的六步：大綱 → 生動度 → 例子 → 程式 → 省 token → 定稿場景卡） |
| 數學 | `codes/utils/walk1d.py`（`python utils/walk1d.py` 自我測試；也在 `check_env.py` 裡） |
| 畫圖 | `codes/utils/walk_draw.py` |
| 動畫 | W101–W112，1080p 在 `codes/renders/`；`Ch1Demo` 25:40 |
| 講義 | `notes/walk_導讀.tex/.pdf`（17 頁）；內文 `notes/ch1/part1.tex`，圖表與數字檢查 `notes/ch1/make_figs.py` → `notes/ch1/figs.tex` |
| combkit | 1.11.0：`combkit.walk1d`、`combkit.draw.walk_draw`、`examples/ex14_random_walk.py`、`tools/sync_from_walk017.py`；`projects/combkit_v1.11.0.zip`（路徑相對於 `C:\mypython\manim`） |
| git | 本資料夾是獨立的 repo（master，無遠端）；`codes/renders/`、LaTeX 中間檔、`ref/png/` 不進 git |

## 3. 檔案結構

```text
walk017/
├── 交接說明_HANDOFF.md        本文件（權威）
├── codes/                      manimgl 一律在這裡執行；說明見 codes/README.md
│   ├── custom_config.yml       只能寫英文
│   ├── main.py                 同名子類別＋Ch1Demo
│   ├── render_all.py / .bat    批次輸出（已存在的影片預設跳過）
│   ├── check_env.py            環境檢查
│   ├── scenes/                 w101_coin_walk.py … w112_steepest_descent.py
│   ├── utils/                  walk1d.py（數學）、walk_draw.py（畫圖）
│   └── renders/                影片（不進 git）
├── notes/
│   ├── plan_ch01.md            第 1 章計畫（§6 定稿場景卡、§7 製作規則、§8 分工）
│   ├── walk_導讀.tex / .pdf    導讀講義（主檔；以後每章一個 \section）
│   └── ch1/                    part1.tex（第 1 章內文）、make_figs.py、figs.tex（產生的）
└── ref/                        僅供個人研讀，不可散布
    ├── ELEMENTS OF THE RANDOM WALK.pdf
    ├── convert_to_md.py        全書 → md/
    ├── page_png.py             書頁截圖 → png/（不進 git）
    └── md/                     00_INDEX.md、00–13 各章、ch01_formulas.md
```

## 4. 怎麼跑

- 動畫、講義、轉檔的指令表：`codes/README.md` 第 2 節。
- 環境：`C:\mypython\manim\環境與交接_ENVIRONMENT.md`；新計畫的通用流程與坑：`C:\mypython\manim\新計畫起步指南_NEW_PROJECT_GUIDE.md`。
- 讀原文省 token：先看 `ref/md/00_INDEX.md` 找到小節與頁碼，再讀那一章的 md 的那幾段；式子看不清楚就 `python page_png.py 頁碼`（在 `ref/`）截圖看。

## 5. 場景與影片（第 1 章）

| 代號 | 主題 | 長度 |
|---|---|---|
| W101 | 擲硬幣走路 | 2:09 |
| W102 | 數路線：二項式係數 | 2:09 |
| W103 | 平均位置是 0 | 2:04 |
| W104 | 均方位移與 √N 定律 | 2:07 |
| W105 | 用生成函數算所有動差 | 2:11 |
| W106 | Stirling 公式 | 2:05 |
| W107 | 二項分布變成高斯 | 2:18 |
| W108 | 遞迴關係與 Pascal 三角形 | 2:05 |
| W109 | 從遞迴到擴散方程 | 2:03 |
| W110 | 傅立葉：用積分數走法 | 2:04 |
| W111 | 生成函數的奇點與動差 | 2:02 |
| W112 | 最陡下降法＋第 1 章回顧 | 2:08 |
| Ch1Demo | 第 1 章連續播放 | 25:40 |

12 幕合計 25:27；每幕都通過 combkit 的嚴格模式試渲染（`tools/render_check.py`），1080p 由 `render_all.py` 一次一支輸出。
和計畫卡（`notes/plan_ch01.md` §6）不同的地方（都是代理人為了版面或正確性做的小調整）：
- W102：字幕的「÷」在 `Text` 裡變成缺字，改寫成「除以 $2^N$」；開頭多畫 3 條回到原點的路線。
- W103：收束改成「大多數人都不在 0」（20 人裡有 6 人剛好在 0）；多了單人示範、即時平均與 $x^2$ 一列。
- W104：多一張對數圖（斜率 1/2）。W106：多了 1!–10! 的乘積、$\Gamma(N+1)=N\Gamma(N)$ 與練習 1.4 的台階圖。
- W107：$N=10,20,40,100$ 重新縮放後用折線疊在一起（不是四張長條圖）；只畫 $|d|\le32$ 的柱子。
- W108：Markov 用兩條路線（HHTT、THTH）示範。W110：書上 (1.54) 後「誤差可忽略」的說法用字幕點出（不是紅框）。
- W111：多一張 $\sum_dd^2C$、$\sum_dd^4C$ 與平均的表（$N=1,2,3,4,10$）。W112：馬鞍面（Fig. 1.5）不畫，只畫兩條剖面；變形後的圍道取 $|z|=z^*$ 的圓。

## 6. 下一章怎麼做（照這個流程，省 token）

1. 讀 `ref/md/00_INDEX.md` 與那一章的 md；式子用 `ref/page_png.py` 截圖核對，必要時照 `ref/md/ch01_formulas.md` 寫一份 `chNN_formulas.md`（含筆誤表）。
2. 寫 `notes/plan_chNN.md`（六步，格式照 `notes/plan_ch01.md`）；場景代號 `WN01…`（第 2 章是 W201…），類別名同代號。
3. 先寫數學模組 `codes/utils/<主題>.py`（純 Python、`check()` 核對書上的例子與所有場景卡數字），需要新圖形就擴充或新增畫圖模組並試畫一次。
4. 派 2 個代理人（Agent 工具，背景）各做一半場景，提示詞照 `C:\mypython\manim\新計畫起步指南_NEW_PROJECT_GUIDE.md` 第 4.3 節；自己同時寫講義的新一節
   （`notes/chNN/partNN.tex`＋`make_figs.py`，在 `walk_導讀.tex` 加 `\input`）。
5. 全部回報後看總覽圖、接進 `main.py`（加 `ChNDemo`）與 `render_all.py`（`CHNN` 清單、`--chN`、`--chN-demo`）與 `.bat` 的 DEMOS，一次一支輸出 1080p。
6. combkit：把新模組名稱加進 combkit 的 `tools/sync_from_walk017.py` 的 FILES，執行它，加範例、README、升版、`make_zip.py`。
7. 更新本文件與 `codes/README.md`、`C:\mypython\manim\install.log` 附一行進度（混合編碼：讀 bytes、接上 UTF-8 bytes 再寫回）、git commit。

## 7. 本專案的注意事項（通用的坑見起步指南第 6 節）

1. **書上的筆誤很多**（第 1 章 13 處，見 `ref/md/ch01_formulas.md`）：照抄會算錯，例如 (1.34) 的指數、(1.48) 的負號、(1:S-15) 的 $1/\sqrt{2\pi}$。
   做動畫與講義前，每個式子都先用程式核對；在畫面上用 `walk_draw.erratum_box` 紅框標出。
2. **奇偶性的因子 2**：一維格子上 $d$ 與 $N$ 同奇偶，離散機率 $P(N,d)\approx2\times$ 高斯密度；傅立葉積分裡是 $q=0$ 與 $q=\pm\pi$ 兩個鞍點。書上的 (1.40)(1.54) 沒處理好。
3. **PDF 轉檔**：這本 PDF 是 Cambridge Books Online 的各章合併，每章的數學字型子集不同，字碼要用各自字型物件的 /Differences 對應（`convert_to_md.py` 已處理）；
   頁首的 y 座標各章不同（用「本頁最上面一行」判斷）；有些行排得太緊沒有空白字元（依字距補）。用 `python -I` 跑腳本時 `PYTHONIOENCODING` 會被忽略，
   要在程式裡 `sys.stdout.reconfigure(encoding="utf-8")`。
4. 已完成的影片不刪不覆寫：`render_all.py` 預設跳過已存在的影片，要重做才加 `--force`（舊影片先改名保存）。
5. 參考文獻 PDF、md 轉檔、頁面截圖都不可散布（不放進公開的壓縮檔或網路）。
6. 字幕（`zh()`／`Text`）裡的「÷」會變成缺字：寫成「除以」，或放進 `$…$` 用 `\div`。
7. 代理人用 `-s` 想存單張畫面時，ManimGL 開了預覽視窗卡住（這台的 Pillow 9.4 本來就不支援 `-s`）：看畫面一律用 render_check 的總覽圖或 ffmpeg 抽格。

## 8. 進度紀錄

- 2026-10-17：建立專案；全書分章 md 轉檔與索引；第 1 章逐式校對（13 處筆誤）；計畫檔；`walk1d`／`walk_draw`；W101–W112 由 2 個代理人平行製作；
  導讀講義第 1 章（17 頁）；combkit 1.11.0。
