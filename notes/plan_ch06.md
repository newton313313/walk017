# 《Elements of the Random Walk》第 6 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 6 章 The shape of a random walk，pp. 127–166（PDF 頁序 143–182）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_shape.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch05.md` §6–§9。

## 1. 草稿大綱（7 幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W601 | 漫步長什麼樣子？平均起來是一團圓雲，單獨一條卻是長條 | §6.1–6.1.1，pp. 127–129，Fig. 6.1–6.3 |
| W602 | 迴轉半徑張量：像轉動慣量的矩陣；特徵向量＝最長、最短的方向，特徵值＝主迴轉半徑 | §6.1.2–6.1.4，pp. 129–133，(6.1)–(6.9)，Fig. 6.4；補充 1（§6.4、6.4.1） |
| W603 | 非球度：0 是圓球、1 是一條線；三維平均 10/19，但個別值散佈在 0 到 1 | §6.1.4，pp. 133–135，(6.10)–(6.18)，Fig. 6.5，Ex. 6.1 |
| W604 | 無限多維：每一步都往新方向走，每條漫步形狀都一樣；本徵值 $\propto1/k^2$、非球度 2/5 | §6.2 開頭，pp. 135–138，(6.19)–(6.20)，Ex. 6.2–6.4 |
| W605 | Kramers 恆等式：$T$ 寫成步伐的二次式；$\sin k\pi x$ 模態；高斯步伐與配對（Wick）算出 $A_d$ | §6.2.1，pp. 138–142，(6.21)–(6.37)，Ex. 6.5–6.6；補充 2–3 只取結論 |
| W606 | 1/d 展開：本徵值從 δ 尖峰變成有寬度的峰；$d=100,20,3$ 的模擬對照 | §6.2.2–6.2.5，pp. 142–152，(6.38)–(6.71)，Fig. 6.6–6.14，Ex. 6.7 |
| W607 | 主迴轉半徑的分布、Fixman 的精確結果、自迴避漫步的形狀＋第 6 章回顧 | §6.2.6–6.3，pp. 152–154，(6.72)–(6.77)，Fig. 6.15，Ex. 6.8–6.9 |

## 2. 生動度檢討

- 開場的「錯覺」：先畫 1000 條漫步疊成的圓雲（Fig. 6.1），問「一條漫步也是圓的嗎？」再一條條抽出來——幾乎都是長條（Fig. 6.3）。
- W602 用二維的主例（固定種子的一條 1000 步漫步，Fig. 6.4）：畫出重心與兩條主軸，數字就是兩個特徵值；再旋轉整條漫步，特徵值不變、矩陣會變。
- W603 先猜「三維漫步平均有多扁？」再揭曉 10/19；然後用 2 萬條漫步的直方圖（Fig. 6.5）說明「平均」不代表典型。
- W604 的「樓梯」：維度比步數還多時，每一步都垂直於之前所有步——所有漫步都是同一個形狀；本徵值比 $1:\frac14:\frac19$ 和三維模擬的比值很接近（(6.53) vs (6.54)）。
- W606 把 $d=100$、$20$、$3$ 的本徵值直方圖並排（Fig. 6.12–6.14），垂直線是 $d=\infty$ 的預測；看峰怎麼變寬、互相重疊。
- 找錯遊戲照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W601 | 4 步的二維漫步：重心與到重心的平方距離 | 1000 條 100 步的雲（Fig. 6.1） |
| W602 | 一條很短的二維漫步：$2\times2$ 的 $T$ 與兩個特徵值（手算） | 固定種子 1000 步（Fig. 6.4）的主軸；旋轉不變 |
| W603 | 三種形狀的非球度：圓 0、線 1、中間的例子 | $A_3=10/19$、$d\to\infty$ 時 $2/5$；2 萬條的直方圖；$A'_3$（(6.17)） |
| W604 | $N=3$ 步的「樓梯」：重心 $\langle x_n\rangle=l(N+1-n)/(N+1)$（Ex. 6.3(a)） | 本徵值 $(N+1)/\pi^2k^2$；$\sum1/k^4/(\sum1/k^2)^2=2/5$ |
| W605 | 小 $N$ 的 Kramers 矩陣 $a_{\alpha\beta}$ 與它的本徵值 vs $\sin$ 模態 | 配對算 $\langle T_{11}\rangle$、$\langle T_{11}^2\rangle$ → $A_d$ |
| W606 | Wick 配對數（Ex. 6.7）$(2m-1)!!$ | $d=100,20,3$ 的本徵值模擬 vs (6.61)(6.62) |
| W607 | Fixman 分布的平均與變異數 | 最大本徵值的分布（Fig. 6.15）；(6.77) 自迴避 |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_shape.py`（純 Python＋numpy，可 import 同資料夾的模組）：迴轉半徑張量與主軸（(6.1)–(6.4)，兩種算法：定義式 vs Kramers (6.21)）、
  不變量與非球度（(6.6)–(6.15)；$A_d$ 的公式 vs 固定種子模擬；$A'_d$ (6.17)）、$d=\infty$ 的樓梯（本徵值、2/5）、Kramers 矩陣的本徵值 vs $(N+1)/\pi^2k^2$、Wick 配對數、
  1/d 展開的位移與寬度（(6.59)–(6.62)）vs $d=100,20,3$ 的模擬、Fixman (6.74)–(6.76)、自迴避的 (6.77) 只引用。模擬都用固定種子、控制在 `check()` 約 10 秒內。
  `check()` 用兩種以上獨立算法核對全部場景卡數字。
- **重用（不改）**：`walk_draw`（長條、`curve`、`plot_axes`、表格、筆誤框）、`lattice_draw`（足跡、`trail`）。
- **畫圖**：不新增模組；橢圓／主軸、樓梯、直方圖寫在場景檔的模組層級函式。
- combkit：`walk_shape.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex19_walk_shape.py`；升 1.16.0。

## 5. 省 token 檢討

1. **40 頁 → 7 幕**：補充 1（轉動）只在 W602 用一句話；補充 2–3（非球度的長計算、(6.21) 的推導）只取結論，細節放講義；§6.2.2–6.2.5 的圖解（Fig. 6.6–6.11）合成 W606 的一段示意，不逐張重畫。
   場景代理人 A 做 W601–W604、B 做 W605–W607。
2. **重活交給上下文短的代理人**（同第 4、5 章）：基礎代理人（校對、模組、定稿卡、combkit）、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁；`ch06_formulas.md` 只收場景與講義用得到的式子＋筆誤表；補充 2 的逐式計算不校對，只核對最後的 (6.15) 與比值。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_shape.check()` assert；模擬一律固定種子、在模組裡算好，場景只呼叫。
5. **不新增畫圖模組、不做 3D 動畫**（三維漫步用二維投影；Fig. 6.2 的立體對不重畫）。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 6 節目標 8–10 頁**；練習解答 Ex. 6.1–6.9。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_shape.py` 已完成（`python utils/walk_shape.py` 自我測試通過，約 2 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch06_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$（走過 $N+1$ 個點）、維度 $d$、迴轉半徑張量 $\overleftrightarrow T$（Tex 寫 `\mathbf{T}`）、
主迴轉半徑的平方 $R_1^2\ge R_2^2\ge R_3^2$（$R_i$ 是長度）、步伐 $\eta_\alpha$、Kramers 矩陣 $a_{\alpha\beta}$、非球度 $A_d$（平均的比）與 $A'_d$（個別值的平均）。
格子漫步＝每步沿 $2d$ 個方向之一走 1；「高斯步伐」＝每個分量常態、變異數 $1/d$（(6.28)）。1/d 與 Fixman 的模擬：高斯步伐、$N=50$、數值都**除以 $N+1$**（書上 Fig. 6.12–6.15 的橫軸）。
程式代號：`SH` = `utils.walk_shape`、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。
〔注〕書上 (6.1) 寫「$N$ 個足跡、除以 $N$」；本系列一律「$N$ 步、$N+1$ 個點、除以 $N+1$」（和 (6.22) 一致）。Fig. 6.4 說主軸線長「正比於特徵值」：我們畫成半長 $2R_i$（長度單位）。

**W601 漫步長什麼樣子？**（§6.1–6.1.1，pp. 127–129，Fig. 6.1–6.3）
1. 問：1000 個醉漢各走 100 步（二維），疊起來是一團圓雲。單獨一條也是圓的嗎？（先猜。）
2. 小例子（`SH.small_example()`，`SH.SMALL_STEPS`）：4 步「東、東、北、西」→ 點 $(0,0),(1,0),(2,0),(2,1),(1,1)$；重心 $(\frac65,\frac25)$（十字）；
   到重心的平方距離 $\frac85,\frac15,\frac45,1,\frac25$（虛線＋數字）→ 平均 $R_g^2=\frac45$（迴轉半徑的平方：像轉動慣量）。
3. Fig. 6.1 的雲（`SH.cloud()`：1000 條 100 步，種子 601，座標 $|x|,|y|\le30$；畫前 200 條細、半透明的折線，數字用全部 1000 條）：
   平均張量 (6.5)（`SH.cloud_stats()["Tmean"]`）：$\bar T_{11}=8.35$、$\bar T_{22}=8.59$、$\bar T_{12}=0.03$（理論 $8.42,8.42,0$）；平均張量的兩個特徵值 8.60、8.34 → 幾乎是圓。
4. 揭曉（Fig. 6.3，`SH.fig63()`：6 條 1000 步，種子 603，$2\times3$ 小圖，各自置中縮放）：每條的 $R_1^2/R_2^2$（`SH.fig63_ratios()`）＝2.1、3.0、6.7、5.5、14.5、4.6；
   1000 條 100 步的平均（`cloud_stats()`）：$\langle R_1^2\rangle=14.14$、$\langle R_2^2\rangle=2.80$（比 5.05）——典型的一條是長條。
5. 一行：雲是圓的，因為每條往不同方向長；要描述「一條」的形狀，得用不隨旋轉改變的量（下一幕）。無筆誤框。
收束：一群醉漢疊起來是圓的，單獨一個卻幾乎總是長條。

**W602 迴轉半徑張量：形狀的指紋**（§6.1.2–6.1.4，pp. 129–133，(6.1)–(6.9)，Fig. 6.4；補充 1）
1. 問：怎麼用幾個數字說出「這條漫步往哪個方向長、長多少」？（提示：轉動慣量。）
2. (6.1)–(6.3)：$T_{kl}=\frac1{N+1}\sum_j(r_{jk}-\langle r_k\rangle)(r_{jl}-\langle r_l\rangle)$；二維 $2\times2$。小例子手算（`SH.tensor_frac`）：$\overleftrightarrow T=\frac1{25}\begin{pmatrix}14&3\\3&6\end{pmatrix}$；
   特徵值 $R_1^2=\frac35$、$R_2^2=\frac15$，主軸 $(3,1)$（$18.43^\circ$）；不變量 (6.6)–(6.9)：$\mathrm{Tr}\,T=\frac45$、$\mathrm{Tr}\,T^2=\frac25$、$\det T=\frac3{25}$（＝$\frac35\cdot\frac15$）。
3. 補充 1（一張曲線圖 `D.plot_axes`，橫軸角度 $0$–$180^\circ$）：沿方向 $\vec n$ 的平方寬度 $\vec n\cdot\overleftrightarrow T\cdot\vec n$（`SH.extent(T, phi)`）：最大 $\frac35$ 在 $18.43^\circ$、最小 $\frac15$ 在 $108.43^\circ$
   → 特徵向量＝最長、最短的方向，特徵值＝那個方向的平方寬度。
4. 主例 Fig. 6.4（`SH.fig64()`：二維 1000 步，種子 643，$x\in[-49,4]$、$y\in[-24,7]$；`SH.fig64_stats()`）：重心 $(-17.9,-8.9)$；$T_{11}=207.5$、$T_{22}=82.9$、$T_{12}=95.0$；
   $R_1^2=258.8$、$R_2^2=31.6$（$R_1=16.09$、$R_2=5.62$，比 8.18），主軸 $28.4^\circ$。畫兩條主軸（半長 $2R_i$）＋重心。
5. 旋轉（(6:S1-17)，`SH.rotate(p, π/6)`、`fig64_stats()["Trot"]`）：繞重心轉 $30^\circ$ → $T_{11}=94.1$、$T_{22}=196.3$、$T_{12}=101.4$（矩陣變了），特徵值仍 258.8、31.6，主軸變 $58.4^\circ$。
   一行：對所有方向平均，$\langle\overleftrightarrow T\rangle=\bar T\,\overleftrightarrow I$（(6.5)）——形狀資訊被抹掉。無筆誤框。
收束：張量會跟著轉，特徵值不會：主迴轉半徑就是形狀的指紋。

**W603 非球度：0 是圓球、1 是一條線**（§6.1.4，pp. 133–135，(6.10)–(6.18)，Fig. 6.5，Ex. 6.1）
1. 問：三維漫步平均起來有多扁？在 0（圓球）到 1（一條線）之間猜一個數。
2. (6.10)(6.13)：$\sum_{i<j}(R_i^2-R_j^2)^2=d\,\mathrm{Tr}\,T^2-(\mathrm{Tr}\,T)^2$；個別非球度 $\frac{\sum_{i<j}(R_i^2-R_j^2)^2}{(d-1)(\sum R_i^2)^2}$（`SH.asph`）。
   形狀例子（三個示意形狀＋數字）：圓球 $(1,1,1)$：0；線 $(1,0,0)$：1；雪茄 $(4,1,1)$ 與扁盤 $(1,1,0)$：都是 $\frac14$；W602 小例子 $(\frac35,\frac15)$：$\frac14$；Fig. 6.4 主例：0.6117。
3. 揭曉 (6.12)(6.14)(6.15)：$A_d=\frac{\sum\langle(R_i^2-R_j^2)^2\rangle}{(d-1)\langle(\sum R_i^2)^2\rangle}=\frac{4+2d}{4+5d}$（`SH.A_d`）。表：$d=2$：$\frac47=0.5714$；3：$\frac{10}{19}=0.5263$；4：$\frac12$；10：$\frac49=0.4444$；$\infty$：$\frac25$。
4. Fig. 6.5（`SH.sim3d()`：三維格子 2 萬條 100 步，種子 605；`SH.sim3d_stats()`）：個別非球度的直方圖 20 根（`hist`，寬 0.05）：
   119, 616, 1183, 1631, 1890, 1795, 1808, 1621, 1669, 1529, 1505, 1305, 1065, 1000, 655, 389, 162, 56, 2, 0；
   平均非球度 0.5231（$N=100$ 格子的精確值 0.5235，`SH.A_exact(100, 3, "lattice")`；$N\to\infty$：0.5263）；小於 0.5 的占 69.3%——從 0 到 1 都有，平均不代表典型。
5. (6.16)–(6.18)：個別值的平均 $A'_3=0.394274$（`SH.Aprime(3)`；二維 0.396400、四維 0.393829）vs 模擬 0.3939——比 $\frac{10}{19}$ 小（比值的平均 ≠ 平均的比）。
   筆誤框（where「p. 135，(6.16)」）：錯 `(d-1)\sum_{k=1}^{N}R_k^2`、對 `(d-1)\Big(\sum_{k=1}^{d}R_k^2\Big)^2`；字幕：照書上，整條漫步放大 2 倍，值從 0.2 變 0.8（`SH.asph_book16`），就不是形狀了。
收束：三維漫步平均「一半扁」（10/19），但個別的形狀從近乎圓球到近乎直線都有。

**W604 無限多維：每條漫步都是同一座樓梯**（§6.2 開頭，pp. 135–138，(6.19)–(6.20)，(6.53)–(6.54)，Ex. 6.2–6.4）
1. 問：如果空間有 100 維，醉漢走 10 步，會不會走回已經走過的方向？
2. Ex. 6.2（`SH.orth_fraction(10, 100)`）：每一步都換新方向的比例 0.6282（近似 $e^{-N^2/2d}=0.6065$，`SH.orth_approx`）；$d\gg N^2$ 時幾乎每一步都垂直於之前所有步
   → 轉一轉、照鏡子，每條漫步都是同一個形狀（樓梯）。
3. $N=3$ 的樓梯（`SH.staircase(3)`；三根互相垂直的箭頭用 `L.iso_point` 斜投影）：重心 $\langle x_n\rangle=\frac34,\frac12,\frac14$（Ex. 6.3(a)，`SH.staircase_com(3)`）；
   $\overleftrightarrow T=\frac1{16}\begin{pmatrix}3&2&1\\2&4&2\\1&2&3\end{pmatrix}$；特徵值 0.4268、0.1250、0.0732（精確 $\frac{2+\sqrt2}8,\frac18,\frac{2-\sqrt2}8$，`SH.eigs_exact(3)`）
   vs 大 $N$ 公式 $\frac{N+1}{\pi^2k^2}$：0.4053、0.1013、0.0450（`SH.eigs_asym`）；特徵向量是正弦的形狀（三組長條 `D.bar_chart`：$(1,\sqrt2,1)$、$(1,0,-1)$、$(1,-\sqrt2,1)$，`SH.sine_mode(3, k)`）。
4. (6.19)(6.20)：$R_n^2=\frac{N+1}{\pi^2n^2}$；$A_\infty=\frac{\sum n^{-4}}{(\sum n^{-2})^2}=\frac{\pi^4/90}{(\pi^2/6)^2}$＝$\frac{1.0823}{2.7058}=\frac25$（$\sum n^{-2}=1.6449$，`SH.zeta2`、`SH.zeta4`）；
   有限 $N$ 的樓梯（`SH.staircase_asph`）：$N=3$：0.5200（$\frac{13}{25}$）；10：0.4150；100：0.4002 → 0.4。
5. 和三維比：(6.53) $9:2.25:1$ vs 三維模擬（`SH.sim3d_stats()["ratios"]`）$11.95:2.69:1$——差不多。
   筆誤框（where「p. 146，(6.54)」）：錯 `\langle R_1^2\rangle:\langle R_2^2\rangle:\langle R_3^2\rangle=13.76:3.03:1`、對 `\langle R_1^2\rangle:\langle R_2^2\rangle:\langle R_3^2\rangle\approx11.9:2.7:1`；
   字幕：書上的數字接近「比值的平均」$\langle R_1^2/R_3^2\rangle:\langle R_2^2/R_3^2\rangle$（模擬 13.69 : 2.98，`sim3d_stats()["ratio_of_ratios"]`）。
收束：維度夠高時，每條漫步都是同一座樓梯，非球度正好 2/5。

**W605 Kramers 恆等式：形狀是步伐的二次式**（§6.2.1，pp. 138–142，(6.21)–(6.37)，Ex. 6.5–6.6；補充 2–3 只取結論）
1. 問：$\overleftrightarrow T$ 是 $N+1$ 個點的平均。能不能直接寫成「步伐乘步伐」？
2. (6.21)(6.22)：$T_{ij}=\sum_{\alpha,\beta}a_{\alpha\beta}\eta_{\alpha,i}\eta_{\beta,j}$，$a_{\alpha\beta}=\frac{\min(\alpha,\beta)(N+1-\max(\alpha,\beta))}{(N+1)^2}$。$N=4$（`SH.kramers(4)`，$4\times4$ 表 `D.tex_table`，乘 $\frac1{25}$）：
   $\begin{pmatrix}4&3&2&1\\3&6&4&2\\2&4&6&3\\1&2&3&4\end{pmatrix}$；補充 3：表裡的整數＝含 $\eta_m\eta_n$ 的點對數 $m(N+1-n)$（`SH.kramers_count(4)`）。
   代入「東東北西」→ $\frac1{25}\begin{pmatrix}14&3\\3&6\end{pmatrix}$，和 W602 的定義式一樣（`SH.tensor_kramers_frac`）。
   筆誤框 1（where「p. 139，(6.22)」）：錯 `\alpha<\beta,\ \ \alpha>\beta`、對 `\alpha\le\beta,\ \ \alpha>\beta`；字幕：少了對角線，這條漫步會算出 $T_{11}=T_{22}=0$。
3. 正弦模態 (6.24)–(6.27)：$a$ 的特徵向量精確是 $\sin\frac{k\pi\alpha}{N+1}$（`SH.sine_mode(4, k)`：四組長條＋虛線正弦曲線）；特徵值 0.5236、0.1447、0.0764、0.0553
   vs $\frac{N+1}{\pi^2k^2}$：0.5066、0.1267、0.0563、0.0317（$k$ 小時接近）；$\mathrm{Tr}\,a=\frac45$（＝平均的 $R_g^2$）。
4. 高斯步伐 (6.28)–(6.30)＋配對 (6.33)(6.36)(6.37)（`SH.moments_exact`）：$\langle T_{11}^2\rangle=\frac{(\mathrm{Tr}a)^2+2\mathrm{Tr}a^2}{d^2}$、$\langle T_{11}T_{22}\rangle=\frac{(\mathrm{Tr}a)^2}{d^2}$、
   $\langle T_{12}^2\rangle=\frac{\mathrm{Tr}a^2}{d^2}$ → 大 $N$：$\frac{N^2}{20d^2},\frac{N^2}{36d^2},\frac{N^2}{90d^2}$ → 比值 $\frac95$、$\frac92$ → 代入 (6.14) 得 $\frac{4+2d}{4+5d}$。
   筆誤框 2（where「p. 141，(6.35b)」）：錯 `\frac{\pi^4}{36}+\frac{\pi^4}{90}`、對 `\frac{\pi^4}{36}+\frac{2\pi^4}{90}`；字幕：照書上得 $\frac{7N^2}{180d^2}$，不是 $\frac{N^2}{20d^2}$。
5. 對照表（$d=3$）：$N=4$ 精確 0.6013；$N=100$ 精確 0.5265（`SH.A_exact(N, 3)`）；高斯 2 萬條 $N=100$ 模擬（`SH.sim_gauss_moments()`）：
   $\langle T_{11}^2\rangle$ 56.09（精確 56.67）、$\langle T_{11}T_{22}\rangle$ 31.31（31.48）、$\langle T_{12}^2\rangle$ 12.56（12.60）、$A=0.5261$；$N\to\infty$：$\frac{10}{19}=0.5263$。
收束：Kramers 恆等式把形狀變成「步伐的二次式」；高斯步伐一配對，10/19 就算出來了。

**W606 1/d 展開：尖峰變寬**（§6.2.2–6.2.5，pp. 142–152，(6.38)–(6.71)，Fig. 6.6–6.14，Ex. 6.7）
1. 問：在 100 維，三個最大的主迴轉半徑是固定的數嗎？20 維、3 維呢？
2. Wick 配對（Ex. 6.7，`SH.wick_pairings(m)`、`SH.wick_count`）：$2m$ 個 $\eta$ 配成 $m$ 對的方法數 $(2m-1)!!$：1、3、15、105、945（$m=1..5$）；
   $m=2$ 畫出 3 種配對（四個點、三種連線）＝ 標準高斯的 $\langle x^4\rangle=3$。
3. 預解式 (6.38)–(6.42)：$\mathcal R(\lambda)=\mathrm{Tr}\frac1{\lambda I-T}=\sum_j\frac1{\lambda-\lambda_j}$；$\mathrm{Im}\,\mathcal R(\lambda-i\epsilon)=\pi\sum_j\delta(\lambda-\lambda_j)$。
   筆誤框（where「p. 146，(6.52)」）：錯 `\mathrm{Im}\,\mathcal{R}(\lambda-i\epsilon)=\sum_i\delta(\lambda-\lambda_i)`、對 `\mathrm{Im}\,\mathcal{R}(\lambda-i\epsilon)=\pi\sum_i\delta(\lambda-\lambda_i)`；
   字幕：$\frac{\epsilon}{x^2+\epsilon^2}$ 底下的面積是 $\pi=3.1416$（`SH.resolvent_weight()`）。零階（只配相鄰的 $\eta$，(6.50)）：δ 尖峰在 $\frac{1}{\pi^2k^2}=0.1013,0.0253,0.0113$（除以 $N+1$）。
4. 一階 (6.58)–(6.62)：$\langle\Lambda_k\rangle=\frac{N+1}{\pi^2k^2}(1+\frac3{4d})$（$\sum_{m\ne n}\frac{n^2}{m^2-n^2}=\frac34$，`SH.sum_identity`）、變異數 $\frac2d(\frac{N+1}{\pi^2k^2})^2$：
   相對位移 $\frac3{4d}$：$d=100$：0.0075；20：0.0375；3：0.2500；相對寬度 $\sqrt{2/d}$：0.1414、0.3162、0.8165。
5. Fig. 6.12–6.14（`SH.eig_sim(d)`：$d=100,20$ 各 2000 條、$d=3$ 2 萬條，$N=50$；三張直方圖並排，每張三色、每色 ≤ 30 根，垂直虛線＝0.1013、0.0253、0.0113）＋表（`SH.onedee_table(d)`，「平均÷(6.61)；變異數÷(6.62)」）：
   $d=100$：1.01、1.00、1.00；0.98、0.97、0.94。$d=20$：1.01、1.01、1.00；1.02、0.87、0.68。$d=3$：1.01、0.90、0.75；1.01、0.47、0.23
   （書上 (6.63)–(6.71)：$d=100$ 變異數 0.97、0.93、0.95；$d=3$：0.98、0.48、0.23 與 1.01、0.91、0.75）。
收束：維度越低，尖峰越寬、越疊在一起；最大的主半徑一直最準。

**W607 主半徑的分布、Fixman、自迴避＋第 6 章回顧**（§6.2.6–6.3，pp. 152–154，(6.72)–(6.77)，Fig. 6.15，Ex. 6.8–6.9）
1. 問：三維漫步最大的 $R_1^2$ 分布長什麼樣子？三個加起來（迴轉半徑的平方）呢？
2. (6.72)(6.73)：每個 $R_i^2$ 近似為 $(R_i^2)^{d/2-1}e^{-dR_i^2/2\alpha_i}$（`SH.gamma_density`；Ex. 6.8：變異數 $\frac{2\alpha_i^2}d$）。Fig. 6.15（`D.plot_axes`，$x\in[0,0.6]$）：
   模擬 $R_1^2/(N+1)$ 的直方圖（`SH.eig_sim(3)[0][:, 0]`，30 根，寬 0.02）vs 曲線（$\alpha_1=\frac1{\pi^2}$，峰在 0.0338、高 7.16）；實際平均 0.1280 vs $\alpha_1=0.1013$（(6.61) 的 0.1267 比較準）。
3. Fixman (6.74)–(6.76)（`SH.fixman_density(x, 3)`，≤ 120 點）＋模擬 $\mathrm{Tr}\,T/(N+1)$ 的直方圖（`SH.eig_sim(3)[1]`）：峰在 0.111（高 6.20）；$x=0.1,0.2,0.3$：6.078、3.281、1.011；
   平均 $\frac16=0.1667$、變異數 $\frac1{135}=0.0074$（標準差 0.0861）；數值積分（`SH.fixman_moments_numeric(3)`）1.0000、0.1667、0.0074；模擬 0.1671、0.0076（$N=50$ 精確 0.1666、0.0074）。
   筆誤框（where「p. 153，(6.76)」）：錯 `K(R)=\prod_{n}\Big(1-\frac{2is\alpha_n}{d}\Big)^{-d/2}`、對 `K(s)=\prod_{n}\Big(1+\frac{2is\alpha_n}{d}\Big)^{-d/2}`；
   字幕：照書上的號，$R^2>0$ 的地方算出來全是 0（`SH.fixman_density(x, 3, book=True)`）。Ex. 6.9 一行：獨立的量相加 → $K$ 是各自的乘積。
4. 自迴避 (6.77)（只引用）：$A_d=\frac{2d+4}{5d+4}+0.008(4-d)$，三維 $0.5263+0.008=0.5343$——自迴避讓漫步稍微更扁，但影響不大。
5. 第 6 章回顧（表六列，`D.tex_table`＋中文標籤）：形狀 → $\overleftrightarrow T$ 的特徵值 $R_i^2$（W601–W602）；非球度 $\frac{4+2d}{4+5d}$，三維 $\frac{10}{19}$、個別平均 0.3943（W603）；
   $d\to\infty$ 樓梯 $R_n^2=\frac{N+1}{\pi^2n^2}$、$\frac25$（W604）；Kramers＋Wick（W605）；1/d：位移 $\frac3{4d}$、寬度 $\sqrt{2/d}$（W606）；Fixman 的精確分布（W607）。
收束：漫步沒有「標準形狀」，但形狀有統計：三維平均一半扁，維度越高越確定。

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
  from utils import walk_shape as SH
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_shape` 算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串）；版面用 `assert_in_frame`／`assert_apart`。
  不要在場景裡呼叫 `SH.check()` 或 `SH.summary()`（整套自我測試由 `check_env.py` 跑）。模擬一律用模組的預設參數（種子、條數、步數都在模組裡），數字才對得上。
  最慢的呼叫：`SH.sim3d()` 約 0.4 秒、`SH.eig_sim(d)` 三個 $d$ 合計約 1 秒（都有快取）、`SH.fixman_density` 120 點約 0.08 秒；其餘都很快。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W601 無、W602 無、W603 1 個、W604 1 個、W605 2 個、W606 1 個、W607 1 個。
  框裡的 Tex 不放中文（英文字用 `\text{}`）。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；字幕裡的「÷」寫成「除以」。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「R1 平方」或省略）。
  張量的雙箭頭不要用 `\overleftrightarrow`（Tex 裡寫 `\mathbf{T}`、`\mathbf{a}`）。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——不要拿來當變數；分數用 `Fraction`，Tex 裡寫成 `\frac{a}{b}`。
- 漫步路線：一條路線＝一個 `VMobject().set_points_as_corners(...)`（1001 點沒問題）；點座標 $(x,y)\mapsto$ `origin + unit*(x-x_c)*RIGHT + unit*(y-y_c)*UP`，
  unit 由寬度決定（Fig. 6.4：寬 53 格 → 約 0.11；雲：$\pm30$ → 約 0.09）。主軸：重心出發、方向＝特徵向量、半長 $2R_i\cdot$unit，可加淡色橢圓（半軸 $2R_i$）。
- 直方圖：`D.plot_axes`＋長條（`Rectangle`），每組 ≤ 30 根；高度用密度（計數 ÷ (總數 × 寬度)）時與曲線同一套座標；三張並排時縮小 `plot_axes` 的 width、height，不要 `.scale()`。
- 大量物件：長條最多 60 根（同一畫面）、表格最多 6×6、小路線圖最多 6 張、路線最多 200 條（雲）、點最多 60 個；避免數千個小物件。不做 3D 動畫（三維只畫投影與示意圖）。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式（例如路線、主軸、直方圖、矩陣框），回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W601 | `w601_shape_cloud.py` | `W601ShapeCloud` | A |
| W602 | `w602_gyration_tensor.py` | `W602GyrationTensor` | A |
| W603 | `w603_asphericity.py` | `W603Asphericity` | A |
| W604 | `w604_infinite_dim.py` | `W604InfiniteDim` | A |
| W605 | `w605_kramers_identity.py` | `W605KramersIdentity` | B |
| W606 | `w606_one_over_d.py` | `W606OneOverD` | B |
| W607 | `w607_fixman_saw.py` | `W607FixmanSaw` | B |

連續播放：`main.py` 的 `Ch6Demo`（W601–W607）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch6/part6}`；內文 `notes/ch6/part6.tex`，圖表與數字檢查 `notes/ch6/make_figs.py`
（`check_numbers_ch6()` 逐一 assert，所有數字都從 `codes/utils/walk_shape.py` 算；表格與 TikZ／pgfplots 圖產生成 `notes/ch6/figs.tex`；巨集都以 `\Six` 開頭、`\label` 用 `six-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–5 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch06_formulas.md`、`walk_shape.py` 的 docstring、
`notes/ch5/part5.tex` 與 `notes/ch5/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 8–10 頁、圖表 4–5 張：(1) 迴轉半徑張量：「東東北西」的 TikZ 圖（點、重心、主軸 $(3,1)$）＋手算（$T$、$\frac35$、$\frac15$、$\frac14$）＋Fig. 6.4 主例的數字表（重心、$T$、轉 $30^\circ$ 後的 $T$、$R_i^2$、主軸角）；
(2) 非球度：$A_d$ 表（$d=2,3,4,10,\infty$）＋$A'_d$（0.396400、0.394274、0.393829）＋Fig. 6.5 的直方圖（pgfplots `ybar`，`SH.sim3d_stats()["hist"]` 的 20 個計數）與 0.5231／0.5235／0.5263、0.3939；
(3) 樓梯與 Kramers：$N=3$ 的 $T$、$N=4$ 的 $a$、特徵值精確 vs $\frac{N+1}{\pi^2k^2}$ 的表（$N=3,4$）、樓梯非球度 0.5200／0.4150／0.4002、比值 9 : 2.25 : 1 vs 11.95 : 2.69 : 1 vs 13.69 : 2.98；
(4) 1/d 展開：§6 W606 的三列表（平均、(6.61)、標準差、兩個比值）＋三張小直方圖或一張 $d=3$ 的圖（取樣由 `make_figs.py` 算好寫成座標）；
(5) Fixman：$d=3$ 的密度曲線（`SH.fixman_density`，≤ 120 點）＋Fig. 6.15 的 Gamma 曲線 vs 模擬直方圖，數字 0.1667、0.0074、0.111／6.20、0.0338／7.16、0.1280。
小節對應書上 6.1–6.3（含補充 1–3 的大要），每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch06_formulas.md` 的筆誤表，12 列；
另附〔注〕：$N$ vs $N+1$、(6.19) 的精確 $\sin^2$ 公式、(6.73) 其實是精確的 Gamma、Fixman 要用精確的 $\alpha_n$）。
最後練習 6.1–6.9 的解答：6.1 展開 $\sum_i\sum_j(\vec r_i-\vec r_j)(\vec r_i-\vec r_j)^{\mathsf T}=2n\sum_i\vec r_i\vec r_i^{\mathsf T}-2(\sum\vec r_i)(\sum\vec r_j)^{\mathsf T}$，除以 $2n^2$ 即定義式；
6.2 $(2d)(2d-2)\cdots(2d-2N+2)=(2d)^N\prod_{k<N}(1-\frac kd)\approx(2d)^Ne^{-N^2/2d}$（$d=100$、$N=10$：0.6282 vs 0.6065）；
6.3 (a) $\langle x_n\rangle=\frac{l(N+1-n)}{N+1}$，(c) 共變異 $l^2(N+1)\frac n{N+1}(1-\frac m{N+1})$，(d) $\sin k\pi y$ 的特徵值 $\frac{l^2(N+1)}{\pi^2k^2}$（精確 $\frac{l^2}{4(N+1)\sin^2\frac{k\pi}{2(N+1)}}$），(e) $R_k^2=\frac{l^2(N+1)}{\pi^2k^2}$；
6.4 $-\partial_x^2[x_<(1-x_>)]=\delta(x-y)$，$L(0,y)=L(1,y)=0$（離散版：$(N+1)a$ 是三對角 $(-1,2,-1)$ 的反矩陣）；
6.5 $\sum_\beta\langle\eta_{\alpha,i}\eta_{\beta,j}\rangle=0$ → 非對角 $-\frac{\delta_{ij}}{Nd}$（對角變成 $\frac{1}{d}(1-\frac1N)$）；6.6 配對 → $\frac{(\mathrm{Tr}a)^2}{d^2}=\frac{N^2}{36d^2}$、$\frac{\mathrm{Tr}a^2}{d^2}=\frac{N^2}{90d^2}$；
6.7 $(2m-1)!!=\frac{(2m)!}{2^mm!}$（1、3、15、105、945）；6.8 $\frac{2\alpha_i^2}d$；6.9 獨立的 $R_i^2$ 相加 → 特徵函數相乘，$K(s)=\prod_n(1+\frac{2is\alpha_n}d)^{-d/2}$。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；TikZ 從一點畫回同一點不用 `to[looseness]`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點與直方圖計數由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；`\overleftrightarrow{T}` 可用（講義是 xelatex），或統一寫 $\mathbf T$。
