# 《Elements of the Random Walk》第 7 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 7 章 Path integrals and self-avoidance，pp. 167–192（PDF 頁序 183–208）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_saw.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch06.md` §6–§9。

## 1. 草稿大綱（6 幕；本章只有 26 頁，比第 5、6 章少一幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W701 | 漫步就是「路徑的和」：把 $N$ 步切成兩段、對中間點求和 (7.3)；切成很多段＋高斯 → 路徑積分 (7.8)，$|d\vec r/ds|^2$ 像動能 | §7.1，pp. 167–171，(7.1)–(7.8)，Ex. 7.1 |
| W702 | 自迴避漫步：會記得走過的地方；嚴格（踩到就停）、溫和、「真」自迴避三種；權重 $e^{-v}$ (7.9)(7.10)；二維數一數有幾條、走多遠 | §7.2–7.2.1，pp. 171–172 |
| W703 | 把自迴避展開：$e^{-v}=1-f$ (7.17)、δ 型交互作用 (7.15)(7.16)；一階、二階（3 種）、三階的圖（Fig. 7.1–7.3，Ex. 7.2）；小 $N$ 的逐階部分和 | §7.2.2，pp. 172–176，(7.11)–(7.19) |
| W704 | 生成函數讓摺積變乘積 (7.21)、Fourier 也一樣 (7.24)；費曼規則；幾何級數 → 自能 $\Sigma$ → Dyson 方程 $g=1/(g_0^{-1}-\Sigma)$ | §7.2.3–7.2.4，pp. 176–184，(7.20)–(7.37)，Fig. 7.4–7.10，Ex. 7.3 |
| W705 | 微擾論為什麼失效：數冪次 $(z_c-z)^{1-l(4-d)/2}$，$d=4$ 是分界；改用平均場：只留一個模態 (7.51)、$C_{MF}=\frac{1}{\mathcal N}z_c^{-m}e^{-fm^2/\mathcal N}$ | §7.2.5–7.2.6，pp. 184–191，(7.38)–(7.61)，Ex. 7.4–7.5 |
| W706 | Flory 的論證：熵 $e^{-R^2/m}$ vs 排斥 $e^{-um^2/R^d}$ → $R\propto m^{3/(d+2)}$；和二維精確計數、已知指數比較＋第 7 章回顧 | §7.2.7，pp. 191–192，(7.62)–(7.64) |

## 2. 生動度檢討

- W701 開場：在格子上畫兩條路線「從 A 到 B 的 4 步」，問「一共幾條？」——先走 2 步到中間點、再走 2 步，把每個中間點的「進來方法數 × 出去方法數」加起來（像轉機：經過哪個城市）。
- W702 的對比：同樣 10 步，普通漫步常常繞回原地，自迴避漫步被「逼」得越走越遠；計數表 $c_N$ 和 $4^N$ 並排，存活比例一路掉。
- W703 把「$1-f$」變成橡皮擦：展開後每一項是「一條線＋幾條弧線（碰撞）」；二階三種圖（並排、巢狀、交叉），三階 15 種，和第 6 章的 Wick 配對數 $(2k-1)!!$ 呼應。小 $N$ 逐階加上去：部分和忽大忽小，全部加完才剛好等於自迴避的條數——直接看到「低階不可信」。
- W704 用「火車車廂」比喻幾何級數：$g_0$、$g_0\Sigma g_0$、$g_0\Sigma g_0\Sigma g_0$……接起來就是 Dyson 方程。
- W705 一張表：$d=3$ 時每多一階，奇點的冪次少 $\frac12$；$d=5$ 時反而變好——四維是分界。
- W706 拔河：熵把漫步往內拉、排斥把它往外推，平衡點就是 Flory 指數；$d=2$ 的 3/4 恰好是精確值，$d=3$ 的 0.6 和數值 0.588 很接近。
- 找錯遊戲照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W701 | 一維 $N=4$ 回原點：$1\cdot1+2\cdot2+1\cdot1=6$；二維 $N=4$ 回原點 $=\sum_r C(2;0,r)^2=36$ | 高斯的摺積還是高斯：變異數相加（Ex. 7.1） |
| W702 | 二維 $N=2,3$ 的自迴避條數 12、36（逐條畫出） | 精確計數 $c_N$（$N\le12$ 左右）、$\langle R^2\rangle$ 與普通漫步的 $N$ 比較 |
| W703 | 二維 $N=4$：全部 256 條、自迴避 100 條；用 $1-f$ 逐階展開 | 二、三階圖的數目 3、15；較大 $N$ 的部分和 |
| W704 | 兩個數列的摺積 → 生成函數相乘（具體係數） | 幾何級數的和＝Dyson 方程（代數值核對） |
| W705 | 冪次表 $1-l(4-d)/2$（$d=1..5$，$l=0..4$） | (7.59) 的積分 vs (7.60) 的漸近式（Ex. 7.4）；$z<z_c$ 時沒有效果（Ex. 7.5） |
| W706 | 最大化 $-R^2/m-um^2/R^d$：$d=1,2,3,4$ 得 $1,\frac34,\frac35,\frac12$ | 二維精確計數的局部斜率 vs 3/4；三維已知值（只引用） |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_saw.py`（純 Python＋numpy，可 import 同資料夾的 `walk1d`、`walk_gf` 等）：
  路徑和的切段（(7.3)、(7.5)，一維與二維的精確計數 vs 直接計數）、高斯摺積（Ex. 7.1）；二維自迴避漫步的精確列舉（$c_N$、$\langle R^2\rangle$，有快取，控制在 `check()` 約 10 秒內）、
  嚴格自迴避的 $1-f$ 展開（逐階部分和＝排容原理，全部加完＝$c_N$）、弧線圖的數目 $(2k-1)!!$ 與二階三種圖的列表；生成函數的摺積 (7.21) 與幾何級數／Dyson 方程的數值核對；
  冪次表 (7.39)；平均場 (7.57)–(7.61)：數值積分 vs 漸近式（Ex. 7.4、7.5）；Flory：數值最大化 vs $3/(d+2)$、精確計數的局部斜率。
  `check()` 用兩種以上獨立算法核對全部場景卡數字。
- **重用（不改）**：`walk_draw`（長條、`curve`、`plot_axes`、表格、筆誤框）、`lattice_draw`（格子、足跡、`trail`）。
- **畫圖**：不新增模組；弧線圖（費曼圖）、車廂、拔河圖寫在場景檔的模組層級函式。
- combkit：`walk_saw.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex20_saw_flory.py`；升 1.17.0。

## 5. 省 token 檢討

1. **26 頁 → 6 幕**：§7.2.3 的費曼規則（6 條）在 W704 只畫一張圖說明，細節放講義；(7.43)–(7.50) 的模態展開只取「只留一個模態」的結論；(7.31)(7.34) 的多重積分只寫不算。
   場景代理人 A 做 W701–W703、B 做 W704–W706。
2. **重活交給上下文短的代理人**（同第 4–6 章）：基礎代理人（校對、模組、定稿卡、combkit）、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁；`ch07_formulas.md` 只收場景與講義用得到的式子＋筆誤表（我讀正文時已看到可疑處：(7.11) 指數裡的 $3M_1a^2$ 應該是 $2M_1a^2$？(7.59) 指數的寫法、(7.43dd)、腳註 1 的冪次 $p$）。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_saw.check()` assert；列舉與模擬在模組裡算好（有快取），場景只呼叫。
5. **不新增畫圖模組、不做 3D**；自迴避漫步只畫二維。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 7 節目標 7–9 頁**；練習解答 Ex. 7.1–7.5。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_saw.py` 已完成（`python utils/walk_saw.py` 自我測試通過，約 1 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch07_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$（平均場與 Flory 用 $m$）、走法數 $C(N;\vec x,\vec y)$、格距 $a$、
二維自迴避漫步的條數 $c_N$、終點的平方距離 $\langle R^2\rangle$（普通漫步 $=N$）、逸度 $z$（$z_c=\frac1{2d}$）、傳播子 $g_0=\frac1{1-z\chi(\vec q)}$、自能 $\Sigma$、
格點數 $\mathcal N$（Tex `\mathcal N`；書上 §7.2.6 也寫成 $N$，本系列一律 $\mathcal N$）。二維格子每步 E、N、W、S 之一；「嚴格」自迴避＝$N+1$ 個點兩兩不同。
程式代號：`SW` = `utils.walk_saw`、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。

**W701 漫步就是「路徑的和」**（§7.1，pp. 167–171，(7.1)–(7.8)，Ex. 7.1）
1. 問：從原點走 4 步回到原點，一維有幾條？二維呢？（先猜。）
2. (7.3) $C(N;\vec x,\vec y)=\sum_{\vec r}C(M;\vec x,\vec r)\,C(N-M;\vec r,\vec y)$：先走 $M$ 步到中間點、再走 $N-M$ 步（像轉機：經過哪個城市）。
   一維 $N=4$、$M=2$（`SW.split_terms(4, 2, 0, 1)`）：中間點 $r=-2,0,2$，進來×出去 $=1\cdot1+2\cdot2+1\cdot1=6$（＝暴力列舉 `SW.count_brute(4, 0, 1)`）；數線＋三組小路線圖。
3. 二維 $N=4$（`SW.split_terms(4, 2, (0, 0), 2)`；`L.count_grid` 在 9 個中間點寫 $C(2;\vec r)$）：中心 $4\cdot4=16$、四個 $(\pm2,0),(0,\pm2)$ 各 $1\cdot1$、四個 $(\pm1,\pm1)$ 各 $2\cdot2$
   → $16+4\cdot1+4\cdot4=36=\binom42^2$。$M=N-1$ 退化成 (2.7)：一維 $3+3=6$、二維 $4\cdot9=36$。切 $2+2+2$（(7.5)，`SW.multi_split`）：一維 20、二維 400。
4. 每段很長 → 高斯 (7.6)，Ex. 7.1（`SW.gauss1`、`SW.conv_ratio`）：三維 $M_1=3$、$M_2=6$（$a=1$），一個分量的變異數 $1+2=3$；
   三條曲線（`D.plot_axes`，$x\in[-6,6]$）峰高 0.3989、0.2821 → 摺積 0.2303：高斯的摺積還是高斯，變異數相加。
   筆誤框（where「p. 169，(7.6)、Ex. 7.1」）：錯 `\big(2\pi M a^2\big)^{-3/2}`、對 `\big(2\pi M a^2/3\big)^{-3/2}`；字幕：照書上，每摺積一次多出 $3^{-3/2}=0.1925$（一個分量 0.5774）。
5. 切成很多段 (7.7)(7.8)：指數 $-\frac{3}{2a^2}\frac{|\Delta\vec r_i|^2}{\Delta N_i}$ → $C\propto\prod_s\int d\vec r(s)\exp\Big[-\frac{3}{2a^2}\int_0^N\Big|\frac{d\vec r}{ds}\Big|^2ds\Big]$：像「速度平方」（動能）。
   折線逐漸變成光滑曲線的示意；一行：$\Delta N_i$ 其實要 $\gg1$ 高斯才成立，$d\vec r/ds$ 只是記號。無其他數字。
收束：漫步的條數＝所有路徑的和；切段、求和、再切下去，就是路徑積分。

**W702 自迴避漫步：會記得走過的地方**（§7.2–7.2.1，pp. 171–172，(7.9)(7.10)）
1. 問：如果醉漢不准踩到自己走過的點，10 步之後離起點會比較遠還是比較近？
2. 三種自迴避（三個小圖示＋字幕）：嚴格（踩到就停）、溫和（越靠近越可能停）、「真」自迴避（偏向不踩，但一直走）；本章用嚴格的。
   權重 (7.9) $\prod_{i>j}e^{-v(|\vec w_i-\vec w_j|)}$：重合時 $v=\infty$ → 權重 0 → 那條路徑被刪掉。
3. 逐條畫（`SW.saw_list(N)`，每條一個小路線圖）：$N=2$：16 條中 12 條（4 條原地折返，灰色打叉）；$N=3$：64 條中 36 條，畫第一步向東的 9 條（×4）；$N=4$：256 條中 100 條（只寫數字）。
4. 精確列舉（`SW.saw_c`、`SW.saw_frac`、`SW.saw_R2`；表四列九欄，或表＋`D.plot_axes` 的 $\langle R^2\rangle$ 對 $N$ 點圖與直線 $N$）：
   $N=1,2,3,4,6,8,10,12,14$；$c_N=4,12,36,100,780,5916,44100,324932,2374444$；$c_N/4^N=1.0000,0.7500,0.5625,0.3906,0.1904,0.0903,0.0421,0.0194,0.0088$；
   $\langle R^2\rangle=1.0000,2.6667,4.5556,7.0400,12.5744,19.0128,26.2425,34.1870,42.7864$（普通漫步 $=N$）。
5. 對比（`SW.sample_saw(10)`、`SW.sample_rw(10)`，種子 702，`L.trail`）：自迴避的 10 步終點 $(-3,3)$、$R^2=18$；普通漫步終點 $(1,-1)$、$R^2=2$（中間踩回 3 次）。無筆誤框。
收束：自迴避漫步被自己的足跡推著走：條數越來越少，卻走得越來越遠。

**W703 把自迴避展開：$1-f$**（§7.2.2，pp. 172–176，(7.11)–(7.19)，Fig. 7.1–7.4，Ex. 7.2）
1. 問：「不准踩到自己」很難算。能不能先算普通漫步，再一點一點扣掉？
2. (7.17) $e^{-v}=1-f$：$f$ 只在 $x\approx0$ 附近不是 0（Fig. 7.1 的小圖）；嚴格自迴避：重合時 $f=1$。每一對點乘 $(1-f)$ 再展開：
   0 階＝普通漫步；一階 (7.19) $-\sum C(M_1;\vec x,\vec u_1)C(M_2;\vec u_1,\vec u_2)C(M_3;\vec u_2,\vec y)f(|\vec u_1-\vec u_2|)$＝一條線＋一個圈＋一條虛線（Fig. 7.2）。
   筆誤框（where「p. 173，(7.11)」）：錯 `e^{-3|\vec x-\vec r_1|^2/3M_1a^2}`、對 `e^{-3|\vec x-\vec r_1|^2/2M_1a^2}`；字幕：照書上，一個方向的變異數變成 $M/2$（$M=6$：3），格子漫步是 $M/3$（2）。
3. 弧線圖（模組層級函式：一條水平線上 $2k$ 個點、上方的虛線弧＝交互作用；`SW.matchings(k)`、`SW.SECOND_ORDER`）：二階 3 種＝Fig. 7.3 (a) 巢狀、(b) 並排、(c) 交叉；
   三階 15 種（Ex. 7.2，$3\times5$ 小圖）；$(2k-1)!!=1,3,15,105$——和第 6 章的 Wick 配對一樣。
4. 二維 $N=4$ 逐階算（`SW.ie_terms(4)`、`SW.ie_partial(4)`）：第 $k$ 階 $T_k=256,228,96,28,4$（一階 $=3\cdot4^2\cdot4+36=228$：三對隔 2 步的點、一對隔 4 步，`SW.order1_formula(4)`）；
   部分和 $256\to28\to124\to96\to100$（`D.bar_chart`），最後剛好＝自迴避的 100 條（排容原理）。二階的 96＝巢狀 16＋交叉 32＋共用一點 48（`SW.order2_types(4)`）。
5. $N=8$（`SW.ie_partial(8)`）：$65536\to-119332\to190940\to-222036\to\cdots\to5916$（共 17 項）：低階的部分和連正負都不對——要全部加起來才可信（伏筆：W705）。
收束：$1-f$ 展開把自迴避變成「一條線加很多弧」的圖；只取前幾階，數字可能錯得離譜。

**W704 生成函數與 Dyson 方程**（§7.2.3–7.2.4，pp. 176–184，(7.20)–(7.37)，Fig. 7.5–7.10，Ex. 7.3）
1. 問：展開的每一項都是「切段再求和」的摺積，越高階越難算。有沒有辦法把摺積變成乘法？
2. (7.20)–(7.22) $\sum_Nz^N\sum_{M}H_1(M)H_2(N-M)=h_1(z)h_2(z)$。例（`SW.EX_CONV`、`SW.convolve`）：$H_1=H_2=2^K=1,2,4,8,16$（一維 $K$ 步的走法數）
   → 摺積 $1,4,12,32,80=(N+1)2^N$ ↔ $\frac1{(1-2z)^2}$；$z=0.1$：1.5625（`SW.gf_check(0.1)`）。Fourier 也一樣 (7.24)（一行）。
   筆誤框 1（where「p. 176，(7.21)」）：錯 `z^{M_1}H_1(M_1)\,Z^{N-M_1}H_2(N-M_1)`、對 `z^{M_1}H_1(M_1)\,z^{N-M_1}H_2(N-M_1)`；字幕：大寫 $Z$ 成了另一個變數，乘積就拆不開。
3. 費曼規則（一張圖，模組層級函式）：實線＝傳播子 $g_0=\frac1{1-z\chi(\vec q)}$、虛線＝$-F(\vec Q)$、頂點動量守恆；一階 (7.26) $g_0\Big[\sum_{\vec Q}F(\vec q-\vec Q)g_0(\vec Q)\Big]g_0$。
   並排的兩個圈（Fig. 7.3(b)）＝$g_0\Sigma_1g_0\Sigma_1g_0$（字幕：書上 p. 178 寫成 (a)，(a) 其實是巢狀）。
4. 火車車廂（(7.30)(7.36)）：$g_0+g_0\Sigma g_0+g_0\Sigma g_0\Sigma g_0+\cdots=\frac{1}{g_0^{-1}-\Sigma}$。數值（`SW.DYSON_TOY`：$g_0=2$、$\Sigma=\frac14$；`SW.dyson_partial(2, Fraction(1, 4), 5)`）：
   部分和 $2,\ 3,\ \frac72,\ \frac{15}4,\ \frac{31}8$ → $4=\frac1{1/2-1/4}$（車廂一節節接上、長條逼近虛線 4）。
   筆誤框 2（where「p. 183，(7.35)」）：錯 `\frac{g_0}{1+(\Sigma_1+\Sigma_2)g_0}`、對 `\frac{g_0}{1-(\Sigma_1+\Sigma_2)g_0}`；字幕：照書上得 $\frac43$，不是 4（同頁 (7.33) 少了中間的 $g_0$，照抄得 $\frac{10}3$）。
5. 自能＝不可約（1PI，剪一條實線分不開）的圖；二階有巢狀 (7.34) 與交叉（Fig. 7.9）兩種。計數版的 Dyson（`SW.dyson_counts(4)`）：
   全部弧線圖 $1,3,15,105$＝不可約的 $1,2,10,74$ 排成一串，即 $\frac{1}{1-(x+2x^2+10x^3+74x^4)}$ 的係數。
收束：生成函數把摺積變乘積，幾何級數把無限多張圖收成一個分母：$g=1/(g_0^{-1}-\Sigma)$。

**W705 微擾論為什麼失效、平均場**（§7.2.5–7.2.6，pp. 184–191，(7.38)–(7.61)，Ex. 7.4–7.5）
1. 問：自能的每一階都帶一個小小的 $f$。是不是只要算前幾階就夠了？
2. 數冪次 (7.38)(7.39)：$l$ 階＝$2l-1$ 條傳播子 $\frac1{(z_c-z)+q^2}$、$l$ 個 $d^dq$ → 令 $q\to(z_c-z)^{1/2}q$ → $(z_c-z)^{1-l(4-d)/2}$。表（`SW.power_table()`，列 $d=1..5$、欄 $l=0..4$）：
   $d=1$：$1,-\frac12,-2,-\frac72,-5$；$d=2$：$1,0,-1,-2,-3$；$d=3$：$1,\frac12,0,-\frac12,-1$；$d=4$：$1,1,1,1,1$；$d=5$：$1,\frac32,2,\frac52,3$。
   數值積分量出來（`SW.loop_slope`）：$l=1,d=1$：$-0.50$；$l=1,d=3$：$0.50$；$l=2,d=1$：$-2.00$。
3. 一維一階自能 $\int\frac{dq}{2\pi}\frac1{1-2z\cos q}=\frac1{\sqrt{1-4z^2}}$（`SW.sigma1_1d`）：$z=0.4,\ 0.49,\ 0.499$ → 1.6667、5.0252、15.8193（離 $z_c=\frac12$ 每近 10 倍就大 $\sqrt{10}$ 倍）。
   (7.42)：截到 $n_h$ 階，奇點 $(z_c-z)^{-p}$、$p=1-n_h(4-d)/2$：$d<4$ 每多一階 $p$ 就變小（$d=3$ 每階少 $\frac12$），$f$ 再小也壓不住；$d=4$ 是分界，$d>4$ 高階反而不要緊（腳註 1）。
4. 平均場：只留一個模態 (7.51) $G\approx z_c\frac{\psi_0(\vec x)\psi_0^*(\vec y)}{z_0-z}$、$\psi_0=\frac1{\sqrt{\mathcal N}}$ → (7.57) $G_{MF}=\int_0^\infty\exp\big[\mathcal N\big(-f\lambda^2-\lambda(1-z/z_c)\big)\big]d\lambda$
   → (7.59)(7.60) $C_{MF}(m)=\frac{\mathcal N^m}{z_c^mm!}\int_0^\infty e^{-\mathcal N(f\lambda^2+\lambda)+m\ln\lambda}d\lambda\approx\frac1{\mathcal N}z_c^{-m}e^{-fm^2/\mathcal N}$。
   筆誤框（where「p. 190，(7.59)」）：錯 `e^{-\mathcal N(f\lambda^2+\mathcal N\lambda)+m\ln\lambda}`、對 `e^{-\mathcal N(f\lambda^2+\lambda)+m\ln\lambda}`；字幕：照書上，$f=0$ 也回不到普通漫步的 $z_c^{-m}/\mathcal N$。
5. Ex. 7.4（`SW.cmf_ratio`、`SW.cmf_asym`、`SW.MF_CASES`）：$\mathcal Nz_c^mC_{MF}(m)$ 的精確積分 vs $e^{-fm^2/\mathcal N}$（$fm^2/\mathcal N=1$）：$(m,\mathcal N)=(10,100)$：0.3385；$(100,10^4)$：0.3643；$(1000,10^6)$：0.3675 → $e^{-1}=0.3679$。
   Ex. 7.5（`SW.ex75`）：$z/z_c=\frac12$：$\mathcal N(1-z/z_c)G_{MF}=0.6899,\ 0.9341,\ 0.9922,\ 0.9992$（$\mathcal N=10,10^2,10^3,10^4$）→ 1：條數就是 $z_c^{-m}/\mathcal N$，自迴避沒有效果；
   對照 $z/z_c=\frac32$（`SW.ex75_big`）：$\ln G_{MF}/\mathcal N=0.0452,\ 0.0596,\ 0.0621$（$\mathcal N=10^2,10^3,10^4$）→ $0.0625$。
收束：四維以下，低階微擾一定失效；平均場只留一個模態，得到 $e^{-fm^2/\mathcal N}$。

**W706 Flory 的論證＋第 7 章回顧**（§7.2.7，pp. 191–192，(7.61)–(7.64)）
1. 問：$m$ 個單體的高分子（自迴避漫步），大小 $R$ 怎麼隨 $m$ 變大？
2. (7.61) 體積 $\mathcal N=\alpha R^d$、$u=f/\alpha$：$C_{MF}=\frac{(2d)^m}{\mathcal N}e^{-um^2/R^d}$（排斥：$R$ 越大越好）；(7.62) 熵 $e^{-R^2/m}$（$R$ 越小越好）。拔河圖（模組層級函式）。
3. (7.63)(7.64)：最大化 $-\frac{R^2}{m}-\frac{um^2}{R^d}$（`SW.flory_F`；`D.plot_axes` 畫 $d=2$、$m=1000$、$u=1$ 的兩項與和，$R\in[50,600]$）：最大在 $R^*=177.83$（`SW.flory_Rstar(1000, 2)`），
   此時兩項都是 $-31.62$、和 $-63.25$；一般 $R^{d+2}=\frac{du}{2}m^3$ → $R\propto m^{3/(d+2)}$。
4. 表（`SW.flory_nu`、`SW.flory_slope`、`SW.NU_QUOTED`）：$d=1,2,3,4$：Flory $1,\frac34,\frac35,\frac12$（數值最大化的斜率 1.0000、0.7500、0.6000、0.5000）；已知值 1、$\frac34$（精確）、0.588（數值，只引用）、$\frac12$。
   二維精確列舉的局部斜率（`SW.two_nu_eff`，$\log\langle R^2\rangle$ 對 $\log N$ 的斜率＝$2\nu$）：$N=6,8,10,12,14$：1.4306、1.4372、1.4442、1.4505、1.4556 → $\frac32$（橫虛線）。
5. 第 7 章回顧（表六列，`D.tex_table`＋中文標籤）：路徑和＝切段求和、高斯變異數相加（W701：6、36）；自迴避的精確計數 $c_{14}=2374444$、存活 0.0088、$\langle R^2\rangle=42.79$ vs 14（W702）；
   $1-f$ 展開＝弧線圖 $(2k-1)!!$，部分和 256, 28, 124, 96, 100（W703）；Dyson $g=1/(g_0^{-1}-\Sigma)$（W704）；冪次 $1-l(4-d)/2$、平均場 $e^{-fm^2/\mathcal N}$（W705）；Flory $\nu=\frac3{d+2}$（W706）。無筆誤框。
收束：熵和排斥拔河，拔出 $R\propto m^{3/(d+2)}$：二維剛好是精確的 $\frac34$，三維的 0.6 也很接近 0.588。

## 7. 製作規則（代理人必讀）

- ManimGL 1.7.2（`from manimlib import *`），**不是** Manim Community（用 `ShowCreation` 不用 `Create`）。繼承 `ChuScene`；`construct` 只呼叫**模組層級函式**
  （連續播放會借用 `construct`，寫成方法會 AttributeError）。開頭 `self.show_header("…")`＋`check_math()`，結尾收束字幕後 `self.wait(1.0)`。
- 場景骨架：
  ```python
  import sys
  from fractions import Fraction
  from pathlib import Path
  sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
  sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
  from manimlib import *
  from combkit.draw.fulton_draw import assert_apart, assert_in_frame
  from combkit.draw.scene_base import ChuScene
  from combkit.draw.styles import mixed, zh
  from utils import walk_saw as SW
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_saw` 算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串，分數用 `Fraction`）；
  版面用 `assert_in_frame`／`assert_apart`。不要在場景裡呼叫 `SW.check()` 或 `SW.summary()`（整套自我測試由 `check_env.py` 跑）。隨機路線一律用模組的預設種子。
  最慢的呼叫：`SW.saw_enum(14)`（`saw_c`、`saw_R2`、`two_nu_eff` 第一次呼叫時觸發）約 0.2 秒、`SW.ie_partial(8)`／`SW.order2_types(8)` 約 0.1 秒（有快取）；其餘都很快。
  模組有 numpy／scipy 的 RuntimeWarning 時不影響結果（`check_math` 裡可用 `warnings.catch_warnings()` 包起來）。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W701 1 個、W702 無、W703 1 個、W704 2 個、W705 1 個、W706 無。
  框裡的 Tex 不放中文（英文字用 `\text{}`）。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；字幕裡的「÷」寫成「除以」。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「R 平方」或省略）。
  `\mathcal N` 在 Tex 裡可用；負數的分數寫 `-\frac12`。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——不要拿來當變數；`SW` 也是。
- 格子路線：一條路線＝一個 `VMobject().set_points_as_corners(...)` 或 `L.trail(points, cell, origin)`；點座標 $(x,y)\mapsto$ `origin + cell*x*RIGHT + cell*y*UP`。
  W702 的小路線圖：每張是「淡格點＋一條折線＋起點綠點」，cell ≈ 0.22，**同一畫面最多 12 張**（N=2 的 12 條、N=3 的 9 條）。
- 弧線圖（W703、W704）：模組層級函式 `chord_diagram(m, width, origin)`：水平實線（漫步）上 $2k$ 個點，每對 $(a,b)$ 一條上方的虛線弧（`ArcBetweenPoints(..., angle=-PI)` 再 `DashedVMobject`，
  或直接畫實線弧、顏色區分）；巢狀／並排／交叉用同一函式。三階 15 張排成 $3\times5$，每張寬 ≤ 2.2。費曼圖（實線箭頭＋虛線）也寫成模組層級函式。
- 長條與曲線：`D.bar_chart`（部分和有正有負時，長條往下畫；N=8 的 17 項數值差很大，畫成表或只畫前 5 項＋最後一項的長條，不要硬塞）；`D.plot_axes` 的刻度數字自己加。
- 大量物件：長條最多 60 根（同一畫面）、表格最多 6 列 × 10 欄、小路線圖最多 12 張、點最多 60 個；避免數千個小物件。不做 3D。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式（弧線圖、費曼圖、車廂、拔河圖、矩陣框），回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W701 | `w701_path_sum.py` | `W701PathSum` | A |
| W702 | `w702_self_avoiding.py` | `W702SelfAvoiding` | A |
| W703 | `w703_one_minus_f.py` | `W703OneMinusF` | A |
| W704 | `w704_dyson.py` | `W704Dyson` | B |
| W705 | `w705_mean_field.py` | `W705MeanField` | B |
| W706 | `w706_flory.py` | `W706Flory` | B |

連續播放：`main.py` 的 `Ch7Demo`（W701–W706）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch7/part7}`；內文 `notes/ch7/part7.tex`，圖表與數字檢查 `notes/ch7/make_figs.py`
（`check_numbers_ch7()` 逐一 assert，所有數字都從 `codes/utils/walk_saw.py` 算；表格與 TikZ／pgfplots 圖產生成 `notes/ch7/figs.tex`；巨集都以 `\Seven` 開頭、`\label` 用 `seven-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–6 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch07_formulas.md`、`walk_saw.py` 的 docstring、
`notes/ch6/part6.tex` 與 `notes/ch6/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 7–9 頁、圖表 4–5 張：(1) 路徑和：一維、二維 $N=4$ 的切段表（中間點、進來、出去、乘積；6、36）＋切 $2+2+2$ 的 20、400＋Ex. 7.1 的三條高斯（pgfplots，變異數 1、2、3，峰高 0.3989、0.2821、0.2303）；
(2) 自迴避的精確計數：§6 W702 的表（$c_N$、$c_N/4^N$、$\langle R^2\rangle$，$N=1..14$ 可全列）＋$\log\langle R^2\rangle$ 對 $\log N$ 的圖（自迴避的點、普通漫步的直線、斜率 $\frac32$ 的參考線）＋12 條兩步自迴避的 TikZ 小圖；
(3) $1-f$ 展開：二階 3 種弧線圖（TikZ）＋$N=4,6,8$ 的逐階部分和表（W703 的數字；$N=8$ 列全部 17 項）＋二階的分類（巢狀 16、交叉 32、共用一點 48）；
(4) Dyson：(7.30)(7.37)、費曼規則的摘要、`DYSON_TOY` 的部分和、計數版 $(2k-1)!!$ vs 不可約 $1,2,10,74$ 的表；
(5) 冪次表（§6 W705，5×5）＋平均場 Ex. 7.4、7.5 的數字表＋Flory 表（$\nu_F$、數值斜率、已知值、二維的 $2\nu_{\rm eff}$）。
小節對應書上 7.1–7.2.7，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch07_formulas.md` 的筆誤表，13 列；
另附〔注〕：(7.6) 省略 $6^{M}$ 與格點和的 $a^{-3}$、(7.13)(7.14) 的 $v_0$、規則 3 的 $-F$ 沒寫進 (7.26)–(7.34)、Fig. 7.10 是交叉圖而 (7.34) 是巢狀圖、
幾何級數只在 $|\Sigma g_0|<1$ 收斂、(7.52) 少 $z_c^2$ 與 (7.57) 併掉 $\frac12$、書上格點數也寫 $N$、Flory 的熵因子正確是 $e^{-dR^2/2m}$）。
最後練習 7.1–7.5 的解答：7.1 每個分量配方：$\int e^{-\alpha(x-r)^2-\beta(r-y)^2}dr=\sqrt{\frac{\pi}{\alpha+\beta}}e^{-\frac{\alpha\beta}{\alpha+\beta}(x-y)^2}$，$\alpha=\frac{3}{2M_1a^2}$、$\beta=\frac3{2M_2a^2}$ →
$\frac{\alpha\beta}{\alpha+\beta}=\frac{3}{2(M_1+M_2)a^2}$，前因子合成 $(2\pi Na^2/3)^{-1/2}$，歸納到 $k$ 段（前因子要用 $(2\pi Ma^2/3)^{-3/2}$；例 $1+2=3$）；
7.2 三階 $(2\cdot3-1)!!=15$ 種（列出），其中不可約 10 種、可約 5 種（$1\cdot1\cdot1$ 一種、$1\cdot2$ 與 $2\cdot1$ 各兩種）；頭尾倒過來算同一種則 11 種；
7.3 (a) 巢狀 $g_0(\vec q)\Big[\sum_{\vec Q_1,\vec Q_2}F(\vec q-\vec Q_1)g_0(\vec Q_1)F(\vec Q_1-\vec Q_2)g_0(\vec Q_2)g_0(\vec Q_1)\Big]g_0(\vec q)$；(b) 並排 $g_0\Big[\sum_{\vec Q}F(\vec q-\vec Q)g_0(\vec Q)\Big]g_0\Big[\sum_{\vec Q}F(\vec q-\vec Q)g_0(\vec Q)\Big]g_0$；
(c) 交叉 $g_0(\vec q)\Big[\sum_{\vec Q_1,\vec Q_2}F(\vec q-\vec Q_1)F(\vec Q_1-\vec Q_2)g_0(\vec Q_1)g_0(\vec q+\vec Q_2-\vec Q_1)g_0(\vec Q_2)\Big]g_0(\vec q)$（每條虛線 $-F$，兩條 → 正號）；
7.4 令 $t=\mathcal N\lambda$：$\mathcal Nz_c^mC_{MF}(m)=\langle e^{-fT^2/\mathcal N}\rangle$，$T\sim\mathrm{Gamma}(m+1)$（平均 $m+1$、標準差 $\sqrt{m+1}$）→ $e^{-fm^2/\mathcal N}$，修正約 $f^2m^3/\mathcal N^2$；表 0.3385、0.3643、0.3675 → 0.3679；
7.5 $z<z_c$ 時指數在 $\lambda=0$ 最大：$G_{MF}\approx\frac{1}{\mathcal N(1-z/z_c)}-\frac{2f}{\mathcal N^2(1-z/z_c)^3}+\cdots$ → 條數 $z_c^{-m}/\mathcal N$（數字 0.6899、0.9341、0.9922、0.9992）；
$z>z_c$ 時最大在 $\lambda^*=\frac{z/z_c-1}{2f}$，$G_{MF}\sim e^{\mathcal N(z/z_c-1)^2/4f}$（0.0452、0.0596、0.0621 → 0.0625）。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；TikZ 從一點畫回同一點不用 `to[looseness]`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點、長條與表的數字由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；$\mathcal N$ 與步數 $N$、$m$ 不要混用。
