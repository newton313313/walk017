# 《Elements of the Random Walk》第 10 章 → Manim：大綱與省 token 檢討（2026-10-19）

來源：Rudnick & Gaspari 第 10 章 Walks and the O(n) model: mean field theory and spin waves，pp. 233–254（PDF 頁序 249–270）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_cond.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch09.md` §6–§9。

## 1. 草稿大綱（5 幕；本章 22 頁）

| 幕 | 主題 | 書上 |
|---|---|---|
| W1001 | 關在盒子裡的高分子：和單體池平衡的高分子（Fig. 10.1）、巨配分函數 $G(z)=\sum_NC(N)z^N$ (10.1)；週期邊界（走出去從對面回來）；先取 $n\to0$ 再讓盒子變大；三種情況：漫步比盒子小、差不多、比盒子大（擠滿盒子的「凝聚相」） | §10.1，pp. 233–236，(10.1) |
| W1002 | O(n) 模型的平均場：均勻的 $\vec M$ (10.8)–(10.10)、用 $\delta(\lambda-|\vec M|^2)$ 把 $n$ 維積分變成一維 (10.11)–(10.13)、$n\to0$ → (10.15)，和第 7 章的 (7.57) 一模一樣（註 3 的對照）；$r>0$ 是普通漫步，$r<0$ 多出 $e^{N_0r^2/4u}$ | §10.2，pp. 236–239，(10.2)–(10.15) |
| W1003 | 自旋波：$\vec S=\vec M+\vec\sigma$ (10.18)、一階項因狀態方程消失、二階項 (10.19)–(10.21)；盒子裡的 Fourier 模態 $q_i=2\pi n_i/L$ (10.22)–(10.24)；$n\times n$ 矩陣的兩個特徵值：橫向（$n-1$ 重）與縱向（1 重）(10.28)–(10.31)；$n\to0$ 時橫向模態的重數是 $-1$；自由邊界的狀態方程 (10.17) 與牆邊的 $\tanh$ 形（Ex. 10.1–10.4） | §10.3，pp. 239–243，(10.16)–(10.36) |
| W1004 | 相關洞：生成函數 (10.37)–(10.43)；$r\gg0$：盒子裡的自由漫步 (10.44)；$r\ll0$：平均場 × $[1-\cdots D(\vec x,\vec x')]$ (10.45)，$-D\to$ 常數 $-c/|\vec x-\vec x'|^{d-2}$（Coulomb）(10.46)–(10.48)（Ex. 10.5、Fig. 10.2）；§10.4 的解釋：普通漫步「離起點 $R$ 的步數」$\propto\frac1R\,\mathrm{erfc}(\cdots)$ (10.76)–(10.84) → 修正 $-u/R$ | §10.3 後半、§10.4，pp. 243–246、251–253 |
| W1005 | Redner–Reynolds：把漫步看成巨正則系綜 (10.49)–(10.52)，密度 $\rho=\langle N\rangle/V$ 是序參量；平均場 $\beta'=1$、$\gamma'=\alpha'=0$、$\nu'=\frac12$ (10.53)–(10.64)，標度律 → 上臨界維度 4 (10.65)(10.66)；一般 $\beta'=\nu d-1$、$\alpha'=\gamma'=2-\nu d$ (10.67)–(10.75)、Table 10.1＋第 10 章回顧 | §10.3.1，pp. 246–251 |

## 2. 生動度檢討

- 全章的主角是一個看得見的例子：**小的週期格子（環面）上的自迴避漫步**。W1001 先在 $4\times4$（或基礎代理人判斷可行的大小）環面上畫幾條漫步：短的像在無限大平面上，長的繞過邊界、最後把格子擠滿（「一盤義大利麵」）。
- W1002 讓第 7 章的式子「再登場」：把 (10.15) 和 (7.57) 上下並排，按註 3 的對照表一項項換，兩式重合。
- W1003 的「負一個模態」：$n=3$ 時兩個橫向一個縱向；$n=1$ 只有縱向；$n\to0$ 時橫向的重數是 $-1$——數學上允許，物理上奇妙。
- W1004 用三維週期格子上算出的 $D(\vec x,\vec x')$ 畫成「洞」的剖面，旁邊放普通漫步的「離起點 $R$ 的步數」$\propto1/R$，兩條曲線形狀相同。
- W1005 用環面列舉直接算 $\rho(z)=\langle N\rangle/V$：$z$ 小時接近 0，$z$ 過了 $z_c$ 迅速上升到接近 1；格子越大轉折越陡——有限大小的「相變」。
- 找錯遊戲照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W1001 | 環面上 $N\le3$ 的自迴避漫步條數和無限平面一樣（4、12、36）；從哪一步開始不一樣 | 環面的 $c_N$ 一路到填滿（$N=V-1$）；和無限平面的 $c_N$ 並排 |
| W1002 | $\int_0^\infty e^{-N_0(r\lambda+u\lambda^2)}d\lambda$ 在 $r>0$、$r<0$ 的兩種漸近 | 數值積分 vs 漸近式；與 (7.57) 的對照（註 3） |
| W1003 | $n=3$ 的 $3\times3$ 矩陣：特徵值與重數、逆矩陣（Ex. 10.3–10.4 的數值核對） | 自由邊界的 $M(x)$（$\tanh$ 形，Ex. 10.1）；Ex. 10.2 的三階項 |
| W1004 | 三維格子 Green 函數的 $\frac{3}{2\pi R}$（第 3 章 `walk_rec`／`walk_gf` 可重用） | 週期盒子裡的 $D(\vec x,\vec x')$ 數值 vs Coulomb 形；(10.83) 的 erfc 形 vs 模擬或精確計數 |
| W1005 | 平均場的 $\rho$：$r>0$ 為 0、$r<0$ 正比於 $|r|$ | 環面列舉的 $\rho(z)$ 曲線（兩種大小）；Table 10.1 在 $d=2$（$\nu=\frac34$）、$d=3$（0.588，只引用）、平均場 |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_cond.py`（純 Python＋numpy，scipy 可用；可 import 同資料夾的 `walk_saw`、`walk_crit`、`walk_gf`、`walk_rec`，照 `walk_crit` 的相對 import 寫法）：
  環面上的自迴避漫步精確列舉（位元遮罩 DFS、快取；大小以 `check()` 約 10 秒內為準）與 $\langle N\rangle(z)/V$；平均場積分 (10.13)–(10.15) 的數值 vs 漸近、與 (7.57) 的對照；
  自旋波矩陣的特徵值、重數與逆矩陣（任意 $n$、隨機 $\vec M$ 的數值核對）、自由邊界的 $\tanh$ 解；週期盒子裡的 $D(\vec x,\vec x')$ 與 Coulomb 形；(10.83) 的 erfc 形 vs 普通漫步的精確計數（或固定種子模擬）；
  Redner–Reynolds 的平均場與一般指數（Fraction）、標度律驗算。`check()` 用兩種以上獨立算法核對全部場景卡數字，約 10 秒內。
- **重用（不改）**：`walk_draw`、`lattice_draw`（格子、足跡）、`walk_saw`、`walk_crit`、`walk_gf`／`walk_rec`（三維 Green 函數）。
- **畫圖**：不新增模組；環面上的漫步（越界時從對面畫回來）、洞的剖面寫在場景檔的模組層級函式。
- combkit：`walk_cond.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex23_torus_saw.py`；升 1.20.0。

## 5. 省 token 檢討

1. **22 頁 → 5 幕**：(10.32)–(10.43) 的長式只取結構（「平均場 × 自旋波的行列式 × 兩種積分」），細節放講義；§10.4 併進 W1004。
   場景代理人 A 做 W1001–W1003、B 做 W1004–W1005（5 幕給一個代理人時，它的上下文會隨每次 render 越長越大，總用量反而不比兩個代理人省）。
2. **重活交給上下文短的代理人**：基礎代理人、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁；`ch10_formulas.md` 只收場景與講義用得到的式子＋Table 10.1＋筆誤表（我讀正文時看到的可疑處：(10.4) 的 $r=(T-T_c)/T_c$ 與第 9 章 $r$ 的 2 倍差、(10.13)(10.14) 的「$2\pi f(0)=2\pi$」與「對 $S$ 的積分是 $i\pi$」、註 3 的對照、Case 2 的 $\lambda\approx|r|/u$ 與極大值位置、(10.58) 的 $\rho_{MF}\to|r|$、$\gamma'$ 的平均場值（文字層斷掉）、(10.73)、(10.83) 的 $\sqrt{2N/z_c}$）。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_cond.check()` assert；列舉與模擬在模組裡算好（有快取、固定種子）。
5. **不新增畫圖模組、不做 3D 動畫**（三維的 $D$ 只畫剖面曲線）。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 10 節目標 7–8 頁**；練習解答 Ex. 10.1–10.5。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log、記憶；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_cond.py` 已完成（`python utils/walk_cond.py` 自我測試通過，約 4 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch10_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$、每步逸度 $z$、二維 $z_c=1/\mu=0.3791$（`WC.ZC_2D`）；這一章 $e^{-r}=z/z_c$（第一次出現時一行字幕：第 9 章是 $e^{-rN/2}$）；
O($n$) 的分量數 $n$；格點數 $N_0$（$=V$，格距 1）；平均場 $f(\lambda)=e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}$；自旋波 $\lambda_1$（橫向）、$\lambda_2$（縱向）；相關洞 $D$；漫步的新指數 $\beta',\gamma',\alpha',\nu'$（磁鐵的不加撇號）。
程式代號：`WK` = `utils.walk_cond`、`WC` = `utils.walk_crit`（只用 `WC.ZC_2D`）、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。無限平面的 $c_N$ 用 `WK.SAW_INF`。

**W1001 關在盒子裡的高分子**（§10.1，pp. 233–236，(10.1)，Fig. 10.1）
1. 問：一條高分子關在小盒子裡、和一池單體平衡——盒子的牆會怎樣改變它能走的路？Fig. 10.1 示意（模組層級函式：方框裡一條彎曲的長鏈＋約 30 個散落的短線段「單體」；無數字）；
   (10.1) $G(z;\vec x,\vec y)=\sum_NC(N;\vec x,\vec y)z^N$：$z$ 是單體的逸度（fugacity），單體可以接上鏈的兩端、也可以蒸發回池裡。
   筆誤框 1（where「p. 235，(10.1)」）：錯 `\sum_N C(N;\vec x,\vec y)\,z^n`、對 `\sum_N C(N;\vec x,\vec y)\,z^N`；字幕：$n$ 是自旋的分量數；照書上 $n\to0$ 時和 $z$ 無關（4×4 環面永遠是 90677 條，$\langle N\rangle=0$，`WK.book101_mean_N(4, 0.3)`）。
2. 週期邊界：$4\times4$ 環面（16 個格點畫在格子中心，盒子 $[-\frac12,\frac72]^2$；左右邊、上下邊各畫一對相同記號的小箭頭表示「黏在一起」；字幕：走出右邊就從左邊回來）。
   三條例子（`WK.torus_example(name)`，起點 (1, 1)；線段用 `ex["segments"]`，越界的一步拆成兩個半段，越界的半段用另一個顏色）：'short' $N=3$ `WNE`（不越界，像在無限平面上）、
   'mid' $N=8$ `WSSEENNE`（越界 2 次）、'full' $N=15$ `WSSEEENNNEEESSW`（越界 3 次，16 格全部走過——「一盤義大利麵」）。旁邊小圖：'mid' 攤開到無限平面（`ex["lift"]`，`L.trail`，終點 (3, 1)），背景淡淡畫出盒子的複本。
3. 表（`D.tex_table`，分兩段各 3 列 × 9 欄：$N$｜環面 $c_N$｜無限平面 $c_N$；`WK.torus_counts(4)`、`WK.SAW_INF`）：$N=0..7$：1, 4, 12, 36, 96, 264, 624, 1536 vs 1, 4, 12, 36, 100, 284, 780, 2172；
   $N=8..15$：3072, 6240, 9840, 15912, 17472, 19200, 10656, 5712 vs 5916, 16268, 44100, 120292, 324932, 881500, 2374444, 6416596。$N\le3$ 相同（高亮）；$N=4$：96 vs 100——少了 4 條直線（EEEE 繞一圈回到起點）；
   $N=15$（走滿）5712 條；全部 90677 條。一行字：$N=L$ 時剛好少 4 條：3×3 在 $N=3$（32 vs 36）、5×5 在 $N=5$（280 vs 284）（`WK.torus_first_diff(L)` = (L, 4)）。
   字幕：4×4 環面就是四維超立方體：$16\times5712=91392$，正是四維 Gray 碼的條數（`WK.GRAY_Q4`）。
4. 三種情況（三個 4×4 小格子，各畫一條例子走法；`WK.torus_mean_N(4, z)`）：$z=0.2$：$\langle N\rangle=1.4400$（漫步比盒子小）、$z=z_c=0.3791$：4.8899（差不多）、$z=1.0$：11.7409（擠滿，最多 15）。
   字幕：先取 $n\to0$，再讓盒子變大——順序對調會出錯（W1002 會看到為什麼）。
收束：盒子小、漫步長，漫步就被擠成一盤義大利麵——本章要用 O($n$) 模型的平均場與自旋波來描述這個「凝聚相」。筆誤框 1 個。

**W1002 O($n$) 模型的平均場**（§10.2，pp. 236–239，(10.2)–(10.15)，註 2–3）
1. 問：如果每個自旋都一樣（平均場），$n\to0$ 的 O($n$) 模型會給出什麼？(10.2)–(10.4)：$\mathcal H=\sum_i\{\frac r2|\vec S_i|^2+\frac12|\nabla\vec S_i|^2+\frac u4|\vec S_i|^4-\vec h\cdot\vec S_i\}$、
   $G=\lim_{n\to0}\frac1n\sum_\alpha\langle S_{i\alpha}S_{j\alpha}\rangle=\sum_NC(N)e^{-rN}$；(10.5) $\Gamma(r)=\sum_N\Gamma_Ne^{-rN}$（終點加總＝磁化率）。字幕：這章 $e^{-r}=z/z_c$（第 9 章是 $e^{-rN/2}$）。
2. (10.8)–(10.12)：$\vec S_i=\vec M$ ⇒ $\Gamma_{MF}=N_0\lim_{n\to0}\langle M_1^2\rangle$；插入 $1=\int d\lambda\,\delta(\lambda-|\vec M|^2)$：$n$ 維積分變成一維（$d^nM\propto\lambda^{n/2-1}d\lambda$）。
   表（`WK.gamma_n(n)`，$N_0=100$、$r=0.5$、$u=1$）：$n$＝3、2、1、0.1、0.01：1.7183、1.7634、1.8131、1.8624、1.8676 → $n\to0$：1.8682（`WK.gamma_mf(0.5)`）。
3. (10.13)–(10.15)：$\frac{N_0}{2\pi}\iint d\lambda\,dS\,f(\lambda)\frac{e^{iS\lambda}}{2iS}$；對 $S$ 積分（註 2：$-iS\to-iS-\epsilon$，極點在實軸上方）＝ $2\pi i$（$\lambda>0$）、0（$\lambda<0$）
   （`WK.s_integral(1.0)`、`WK.s_integral(-1.0)`，$\epsilon=0.1$：$5.6853i=2\pi ie^{-0.1}$、0）⇒ $\Gamma_{MF}=\frac{N_0}2\int_0^\infty e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}d\lambda$。
   筆誤框 1（where「p. 238，(10.15)」）：錯 `\Gamma_{MF}=\frac{N_0}{4}\int_0^\infty d\lambda\,f(\lambda)`、對 `\Gamma_{MF}=\frac{N_0}{2}\int_0^\infty d\lambda\,f(\lambda)`；
   字幕：書上說對 $S$ 的積分是 $i\pi$（那是主值），應為 $2\pi i$；$u=0$（高斯）時每個 $n$ 都得 $1/r=2$，書上的式子得 1（`WK.gamma_gauss_check()`）；表的極限 1.8682，書上 0.9341。
4. 兩種漸近（`D.plot_axes`：$f(\lambda)$ 兩條曲線，$r=0.5$ 從 $\lambda=0$ 遞減、$r=-0.5$ 在 $\lambda^*=|r|/u=0.5$ 有峰；各 ≤ 60 點，縱軸各自正規化到峰值 1）：
   $r>0$：$\Gamma_{MF}\approx\frac1r$（普通漫步 $\sum e^{-rN}$）：$N_0=100,1000,10^4$：1.8682、1.9844、1.9984 → 2（`WK.gamma_mf(0.5, 1, N0)`）；
   $r<0$：$\Gamma_{MF}\approx\sqrt{\pi N_0/u}\,e^{N_0r^2/4u}$：精確 ÷ 漸近 0.9998、1.0000、1.0000；$\ln\Gamma_{MF}/N_0$＝0.0912、0.0665、0.0630 → $r^2/4u=0.0625$（`WK.ln_gamma_mf`、`WK.ln_gamma_mf_asym`）。
   極限的順序（一行表，`WK.gamma_n(n, -0.5) / WK.gamma_mf(-0.5)`）：$n=10^{-2},10^{-3},10^{-4},10^{-6}$：0.3299、0.8311、0.9801、0.9998——要 $n\,e^{N_0r^2/4u}\ll1$ 才到極限：必須先取 $n\to0$。
5. 註 3：(10.15) 與第 7 章 (7.57) $G_{MF}=\int_0^\infty e^{N(-f\lambda^2-\lambda(1-x))}d\lambda$ 上下並排，一項項換（每換一項變色）：$N_0\to N$、$\lambda\to2\lambda$、$u\to f$、$r\to1-z/z_c$ ⇒ $\Gamma_{MF}=N\,G_{MF}$（終點加總乘 $N$）；
   數值（`WK.fn3_check()`）：$r=0.5$（$x=z/z_c=\frac12$，就是第 7 章 Ex. 7.5 的 $x$）：$\Gamma_{MF}=1.8682=100\,G_{MF}(\frac12)$。
   筆誤框 2（where「p. 238，footnote 3」）：錯 `\lambda\to4\lambda,\ 4u\to f,\ 2r\to1-z/z_c`、對 `\lambda\to2\lambda,\ u\to f,\ r\to1-z/z_c`；
   字幕：書上的換法只對它自己的 $N_0/4$ 成立（兩個 2 倍抵掉：0.9341＝0.9341）；$r\to1-z/z_c$ 才和 (10.4) 的 $e^{-r}=z/z_c$ 一致。
收束：平均場把關在盒子裡的自迴避漫步變成一維積分；$r<0$ 多出的 $e^{N_0r^2/4u}$ 就是第 7 章的那個因子。筆誤框 2 個。

**W1003 自旋波：在平均場上加漣漪**（§10.3，pp. 239–243，(10.16)–(10.36)，Ex. 10.1–10.4）
1. 問：平均場把所有自旋看成一樣——在它上面加一點漣漪（自旋波），會多出什麼？(10.18) $\vec S=\vec M+\vec\sigma$；(10.19)–(10.21)：零階 $\mathcal H_{MF}$、
   一階 $\int(r\vec M-\nabla^2\vec M+u\vec M|\vec M|^2-\vec h)\cdot\vec\sigma$（＝0，因為狀態方程 (10.17)）、二階 $\int\{\frac12|\nabla\vec\sigma|^2+\frac r2|\vec\sigma|^2+\frac u2|\vec M|^2|\vec\sigma|^2+u(\vec M\cdot\vec\sigma)^2\}$。
   Ex. 10.2 一行：三階 $u\int(\vec M\cdot\vec\sigma)|\vec\sigma|^2$、四階 $\frac u4\int|\vec\sigma|^4$（`WK.ex102_coeffs()`：數值內插＝公式，誤差 $<10^{-11}$）。
2. 自由邊界 (10.17)、Ex. 10.1：$-\nabla^2\vec M+r\vec M+u\vec M|\vec M|^2-\vec h=0$，牆上 $M=0$：一維 $M(x)=\sqrt{|r|/u}\tanh(x\sqrt{|r|/2})$（`WK.wall_M(x)`；`D.plot_axes` 橫軸 $x\in[0,6]$、縱軸 $M\in[0,1.2]$，≤ 60 點；虛線 $M=1$；左邊畫一道牆）。
   $r=-1$、$u=1$：$M(1),M(2),M(4)$＝0.6089、0.8884、0.9930；癒合長度 $\sqrt{2/|r|}=1.4142$；有限差分＋牛頓解邊值問題與 tanh 相差 $2.8\times10^{-6}$（`WK.wall_errors()[0]`）。
   字幕：週期邊界沒有牆，$\vec M$ 就是常數（本章接下來都用週期邊界）。
3. Fourier 模態 (10.22)–(10.24)：$\psi_{\vec q}=\prod_ie^{iq_ix_i}/\sqrt{L_i}$、$q_i=2\pi n_i/L$（一維 $L=8$ 的 $\cos(2\pi nx/8)$，$n=0,1,2$ 三條小曲線示意；無數字）。(10.27) $X_{\vec q}=(\frac r2+\frac12q^2+\frac12u|\vec M|^2)I+u|\vec M\rangle\langle\vec M|$。
   例子（Fraction，`WK.sw_X()`、`WK.sw_lambdas()`）：$r=-1$、$u=1$、$q^2=1$、$\vec M=(\frac23,\frac23,\frac13)$（$|\vec M|^2=1$ 正是平均場的鞍點 $|r|/u$）：$X=\frac1{18}\begin{pmatrix}17&8&4\\8&17&4\\4&4&11\end{pmatrix}$；
   特徵值 $\frac12,\frac12,\frac32$：橫向 $\lambda_1=\frac r2+\frac12q^2+\frac12u|\vec M|^2=\frac12$（$n-1=2$ 重：垂直於 $\vec M$ 的平面，畫兩支箭頭）、縱向 $\lambda_2=\lambda_1+u|\vec M|^2=\frac32$（1 重：沿 $\vec M$）。
   字幕：$q=0$ 時 $\lambda_1=0$——整體轉動 $\vec M$ 不花能量。
4. (10.30)(10.31)、Ex. 10.3–10.4：$X^{-1}=\frac I{\lambda_1}-\frac{u|\vec M\rangle\langle\vec M|}{\lambda_1\lambda_2}=2I-\frac43|\vec M\rangle\langle\vec M|=\frac1{27}\begin{pmatrix}38&-16&-8\\-16&38&-8\\-8&-8&50\end{pmatrix}$（`WK.sw_Xinv()`）；$XX^{-1}=X^{-1}X=I$（精確）。
   隨機 $\vec M$、$n=1..6$（`WK.sw_check()`）：特徵值、重數 $(n-1,1)$、左右逆全對（誤差 $<10^{-11}$）。高斯積分 (10.32)：$\prod_{\vec q}\frac{\pi^{n/2}}{\sqrt{\lambda_2}\sqrt{\lambda_1^{n-1}}}\,e^{+\frac14\vec h\cdot X^{-1}\cdot\vec h}$。
   筆誤框 1（where「p. 242，(10.32)」）：錯 `\exp\Big[-\frac14\,\vec h\cdot X_{\vec q}^{-1}\cdot\vec h\Big]`、對 `\exp\Big[+\frac14\,\vec h\cdot X_{\vec q}^{-1}\cdot\vec h\Big]`；字幕：$n=2$ 的數值積分 3.9891，負號得 3.2988（`WK.gauss_n2_check()`）。
5. $n\to0$：表（`WK.mode_table()`）$n$＝3、2、1、0：橫向重數 2、1、0、$-1$；縱向 1；$\det X=\lambda_1^{n-1}\lambda_2$＝$\frac38,\frac34,\frac32,3$；$(\det X)^{-1/2}$＝1.6330、1.1547、0.8165、0.5774 ＝ $\sqrt{\lambda_1/\lambda_2}=1/\sqrt3$。
   字幕：負一個模態——數學上允許，物理上奇妙（$n=1$ 只有縱向、$n=3$ 兩個橫向一個縱向）。
   筆誤框 2（where「p. 242，(10.33)」）：錯 `\sqrt{\frac{\lambda_2}{\lambda_1}}=\sqrt{\frac{r+q^2+u|\vec M|^2}{r+q^2+3u|\vec M|^2}}`、對 `\sqrt{\frac{\lambda_1}{\lambda_2}}=\sqrt{\frac{r+q^2+u|\vec M|^2}{r+q^2+3u|\vec M|^2}}`；
   字幕：右邊對、左邊顛倒（書上的左邊是 $\sqrt3=1.7321$：`float(l2 / l1) ** 0.5`，`l1, l2 = WK.sw_lambdas()`）。
收束：每個 $\vec q$ 有 $n-1$ 個橫向、1 個縱向模態；$n\to0$ 時橫向變成 $-1$ 個，留下 $\sqrt{\lambda_1/\lambda_2}$。筆誤框 2 個。

**W1004 相關洞：頭尾互相排斥**（§10.3 後半、§10.4，pp. 243–246、251–253，(10.37)–(10.48)、(10.76)–(10.84)，Fig. 10.2–10.3，Ex. 10.5）
1. 問：擠滿盒子的漫步，頭和尾喜不喜歡靠在一起？(10.37) $G(z;\vec x,\vec x\,')=\lim_{n\to0}\frac1n\frac{\partial^2\ln Z}{\partial h(\vec x)\partial h(\vec x\,')}$；(10.38)–(10.43) 只取結構
   （三個方塊：平均場 × 自旋波的行列式 × 兩種 $\lambda$ 積分 → $\frac{\vec h\cdot\vec h'}4\{G_0-(\int_0^\infty f\,d\lambda)D\}$）。
2. 兩種情況：$r\gg0$ (10.44) $G_0=\sum_{\vec q}\frac{\psi\psi'}{r+q^2}$——漫步比盒子小，盒子沒影響：$32^3$ 週期格子與無限大格子（$r=1$）相差 $1.2\times10^{-11}$（`WK.case1_check()`）；
   $r\ll0$ (10.45)：$\int_0^\infty e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}d\lambda\,(1-D)$——平均場 ×（$1-D$）。
   字幕（不畫框）：書上 (10.45) 的指數印成 $-N_0(r\lambda+u\lambda^2)$，和 (10.15) 不一致（極大在 $|r|/2u=0.5$，不是正文的 $|r|/u=1$；`WK.book45_argmax()`）。
3. (10.36)→(10.46)：鞍點 $u|\vec M|^2=|r|$ ⇒ $\frac1{q^2}-\frac1{q^2+2|r|}=\frac{2|r|}{q^2(q^2+2|r|)}$。$64^3$ 週期格子（$r=-1$，`WK.D_axis()`）沿軸 $\Delta=0,1,2,4,8,16,32$：$D$＝0.1083、0.0613、0.0359、0.0167、0.0065、0.0016、$-0.0001$。
   畫 $-D(\Delta)$（Fig. 10.2 的樣子；`D.plot_axes` 橫軸 $\Delta\in[0,32]$，33 點）與 Coulomb 形 $c_0-c/\Delta$（`WK.coulomb_fit()`：$c_0=0.0035$、$c=0.0801\approx\frac1{4\pi}=0.0796$，$\Delta\ge2$ 才畫）；
   一行：$1-D$ 從 0.8917（頭尾相鄰）到 1.0001——頭尾靠近時條數少 11%：相關洞。
   筆誤框 1（where「p. 245，(10.46)」）：錯 `D=\sum_{\vec q}\frac{\psi_{\vec q}(\vec x)\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}`、對 `D=\sum_{\vec q\neq0}\frac{2|r|\,\psi_{\vec q}(\vec x)\psi_{\vec q}(\vec x\,')}{q^2(2|r|+q^2)}`；
   字幕：照書上 $D(0)=0.0541$，只有一半（`WK.D_book46()[0]`）。字幕（不畫框）：書上說 $D$「隨距離增加」——其實 $D$ 遞減，增加的是 $-D$（Fig. 10.2 畫的就是 $-D$）。
4. Ex. 10.5、(10.47)(10.48)：連續三維 $D(\rho)=\frac1{4\pi\rho}(1-e^{-\sqrt{2|r|}\rho})\to\frac1{4\pi\rho}$（與 $r$ 無關）：$\rho=1,2,4,8$：0.0602、0.0374、0.0198、0.0099（`WK.D_cont`）；
   無限大格子 0.0648、0.0394、0.0202、0.0100（`WK.D_inf((p, 0, 0))`）——遠處就是 Coulomb：$D\propto|\vec x-\vec x\,'|^{2-d}$。
   字幕：週期盒子只讓 $D$ 整體平移一個常數 $\frac{2.837297}{4\pi L}=0.0035$（Ewald 常數，就是上面的 $c_0$）。
5. §10.4、Fig. 10.3：普通漫步「離起點 $R$ 的步數」(10.76)–(10.84)：sc 格子 $N=100$（`WK.sc_visits(100, R)`）vs $\frac3{2\pi R}\mathrm{erfc}\big(R\sqrt{3/(2(N+1))}\big)$（`WK.erfc_form(R, 100)`）：
   表 $R=1,2,4,8,16$：精確 0.4509、0.1931、0.0610、0.0103、0.000183；erfc 形 0.4121、0.1744、0.0586、0.0100、0.000174；$N\to\infty$ 的 $G(R)\times\frac{2\pi R}3$＝1.0815、1.0779、1.0198、1.0041、1.0009 → 1（`WK.green_sc(R)`）：
   $\propto1/R$——和洞的 Coulomb 形一樣（左右兩條曲線並排）。
   筆誤框 2（where「p. 253，(10.83)」）：錯 `\frac{1}{4\pi z_cR}\,\mathrm{erfc}\Big(\frac{z_cR^2}{2(N+1)}\Big)`、對 `\frac{1}{4\pi z_cR}\,\mathrm{erfc}\Big(\sqrt{\frac{z_cR^2}{4(N+1)}}\Big)`；
   字幕：照書上 $R=8$ 得 $4.3\times10^{-4}$（精確 0.0103）、$R=16$ 得 $1.7\times10^{-28}$（`WK.erfc_form_book`）。
收束：步數的密度 $\propto1/R$，自迴避的懲罰 $-u/R$ 就讓頭尾互相排斥——一個 Coulomb 形的洞。筆誤框 2 個（(10.45) 與「increases」用字幕）。

**W1005 Redner–Reynolds：漫步的密度是序參量**（§10.3.1，pp. 246–251，(10.49)–(10.75)，Table 10.1）
1. 問：把一條漫步的步數當成粒子數——它的「密度」能不能像磁化一樣當序參量？(10.49)–(10.52)：$\Gamma=\sum_N\Gamma_Ne^{-rN}$（巨正則，$r$ 是化學勢）、$\Omega=\frac1V\ln\Gamma$、
   $\langle N\rangle=-\frac{\partial}{\partial r}\ln\Gamma$、$\rho=\frac{\langle N\rangle}V=-\frac{\partial\Omega}{\partial r}$。
2. 環面精確列舉（`WK.torus_curve(L)`，$L=3,4,5$ 三條曲線各 60 點，每條一個 VMobject；`D.plot_axes` 橫軸 $z\in[0,1.2]$、縱軸 $\rho\in[0,1]$；$z=z_c$ 畫鉛直虛線）：
   表 $z=0.2,\ z_c,\ 0.6,\ 1.0$（`WK.torus_rho(L, z)`）：3×3：0.1391、0.3374、0.5284、0.6823；4×4：0.0900、0.3056、0.5585、0.7338；5×5：0.0603、0.2769、0.5764、0.7565（$z\to\infty$：$\frac89,\frac{15}{16},\frac{24}{25}$）。
   (10.59) 的漲落 $\mathrm{Var}\,N/V$ 的峰（`WK.torus_peak(L)`）：$z^*$＝0.4443、0.4380、0.4272，高 0.4289、0.5846、0.7242——格子越大越陡、峰往 $z_c$ 移：有限大小的相變。
3. 平均場 (10.53)–(10.58)：$\rho_{MF}=\frac1V\frac1r\to0$（$r>0$）、$\frac{|r|}{2u}$（$r<0$）（`WK.mf_rho(r)`，$V=10^6$、$u=1$）：$r=-0.5,-0.1,-0.01$：0.2500、0.0500、0.0050；$r=0.1,0.5$：$1.0\times10^{-5}$、$2.0\times10^{-6}$；
   斜率 1.0000 ⇒ $\beta'=1$（磁鐵的 $\beta=\frac12$：$\rho\propto|\vec M|^2$）。$\chi'=\partial^2\Omega/\partial r^2$（`WK.mf_chi(r)`）：0.5000（$r=-0.1$）、$1.0\times10^{-4}$、$4.0\times10^{-6}$（$r=0.1,0.5$）⇒ $\gamma'=\alpha'=0$；$\nu'=\frac12$（(10.62)–(10.64)）。
   筆誤框 1（where「pp. 247–248，(10.53)(10.55)」）：錯 `e^{-\frac{V}{a^d}\frac{r^2}{4u}}`、對 `e^{+\frac{V}{a^d}\frac{r^2}{4u}}`；字幕：照書上的負號 $\rho=-0.2500<0$（$r=-0.5$，`WK.book55_rho(-0.5)`）。
4. 標度律 (10.65)–(10.75)、Table 10.1（`D.tex_table`，4 列 × 4 欄：指數｜平均場 $d=4$｜二維 $\nu=\frac34$｜三維 $\nu\approx0.588$（引用）；`WK.rr_exponents(*WK.RR_CASES[k])`）：
   $\beta'=d\nu-1=1-\alpha$：1、$\frac12$、0.764；$\gamma'=\alpha'=2-d\nu=\alpha$：0、$\frac12$、0.236；$\nu'=\nu$：$\frac12$、$\frac34$、0.588。
   驗算（打勾）：$\alpha'+2\beta'+\gamma'=2$ 全成立；$2-\alpha'=d\nu'$：2、$\frac32$、1.764；平均場 $2-0=d_c\cdot\frac12$ ⇒ $d_c=4$（`WK.dc_mean_field()`）。
   筆誤框 2（where「p. 251，(10.75)」）：錯 `2-\alpha'=d_c\,\nu'`、對 `2-\alpha'=d\,\nu'`；字幕：書上說「一般都成立」，但 $d=2$ 時 $d_c\nu'=3\ne\frac32$（`WK.rr_laws`）。
5. 第 10 章回顧（換頁；表五列，`D.tex_table`＋中文標籤）：W1001 4×4 環面 $c_N$ 從 $N=4$ 起變少（96 vs 100）、走滿 5712 條；W1002 $\Gamma_{MF}=\frac{N_0}2\int f$（書上 $\frac{N_0}4$）、先取 $n\to0$；
   W1003 $\lambda_1$（$n-1$ 重）、$\lambda_2$（1 重）、$n\to0$ 時 $-1$ 個模態、$\sqrt{\lambda_1/\lambda_2}$；W1004 洞 $1-D$：0.8917 → 1、$D\approx\frac1{4\pi\Delta}$、erfc 形；W1005 $\beta'=1$、$\gamma'=\alpha'=0$、$d_c=4$。
收束：漫步的密度在 $z_c$ 兩邊從 0 變成有限值——擠滿盒子的漫步就是 O($n$) 磁鐵的有序相。筆誤框 2 個。

## 7. 製作規則（代理人必讀）

- ManimGL 1.7.2（`from manimlib import *`），**不是** Manim Community（用 `ShowCreation` 不用 `Create`）。繼承 `ChuScene`；`construct` 只呼叫**模組層級函式**
  （連續播放會借用 `construct`，寫成方法會 AttributeError）。開頭 `self.show_header("…")`＋`check_math()`，結尾收束字幕後 `self.wait(1.0)`。
- 場景骨架：
  ```python
  import sys
  import warnings
  from fractions import Fraction
  from pathlib import Path
  sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
  sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
  from manimlib import *
  from combkit.draw.fulton_draw import assert_apart, assert_in_frame
  from combkit.draw.scene_base import ChuScene
  from combkit.draw.styles import mixed, zh
  from utils import walk_cond as WK
  from utils import walk_crit as WC
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_cond`（卡上註明的 `WC.ZC_2D`）算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串，
  科學記號用 `f"{x:.1e}"`（得 `1.0e-05`），分數用 `Fraction`；Fraction 不能直接 `:.4f`，要 `float()`）；版面用 `assert_in_frame`／`assert_apart`。
  不要在場景裡呼叫 `WK.check()`（整套自我測試由 `check_env.py` 跑）；`WK.summary()` 約 2 秒，可以不用（直接呼叫卡上的函式）。
  隨機一律用模組的預設種子（環面例子 1004、自旋波 1003）。最慢的呼叫：`WK.torus_counts(5)`（W1001 的 5×5 一行字、W1005）第一次約 1 秒（同一個行程內有快取）、
  `WK.D_box()` 約 0.1 秒（快取）、`WK.green_sc`／`WK.D_inf` 每點約 5 毫秒、`WK.sc_visits(100, R)` 約 10 毫秒、`WK.torus_peak` 約 10 毫秒；其餘都很快。
  scipy 的 IntegrationWarning 不影響結果（`check_math` 裡用 `warnings.catch_warnings()` 包起來）。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W1001 1 個、W1002 2 個、W1003 2 個、W1004 2 個、W1005 2 個。框裡的 Tex 不放中文（英文字用 `\text{}`）。
  卡上標「不畫框」的書上問題只用一行字幕。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。字幕裡的「÷」寫成「除以」。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；矩陣用 `\begin{pmatrix}…\end{pmatrix}`。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「n 趨近 0」之類）。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——不要拿來當變數（格子邊長用 `Lside`、`size` 之類；**相關洞的 $D$ 在程式裡叫 `Dv`、`hole`**）；`WK`、`WC` 也是。
- 環面（W1001）：格點、盒子、對邊記號都是模組層級函式；漫步的每段用 `ex["segments"]` 的 `((x0, y0), (x1, y1), wrapped)`，座標乘 cell 加 origin（同 `L.lat_point`）；
  一條漫步＝**一個** `VMobject().set_points_as_corners(...)` 不行（越界處斷開），所以每段一條 `Line`（最多 18 段）包成 VGroup，越界的段換色；攤開圖用 `L.trail(ex["lift"], …)`。
  曲線（W1002 的 $f(\lambda)$、W1003 的 tanh、W1004 的 $-D$ 與 erfc、W1005 的 $\rho(z)$）一條＝一個 VMobject（`D.curve` 或 `set_points_as_corners`），每條 ≤ 60 點；刻度數字自己加。
  W1003 的自旋（$\vec M$ 與兩個橫向方向）用平面上的箭頭示意（不做 3D）；W1004 的三個方塊與 W1005 的回顧表沒有新數字。
- 大量物件：長條最多 60 根（同一畫面）、表格最多 6 列 × 10 欄、點最多 60 個、箭頭最多 64 支；避免數千個小物件。不做 3D。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式，回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W1001 | `w1001_confined_polymer.py` | `W1001ConfinedPolymer` | A |
| W1002 | `w1002_mean_field.py` | `W1002MeanField` | A |
| W1003 | `w1003_spin_waves.py` | `W1003SpinWaves` | A |
| W1004 | `w1004_correlation_hole.py` | `W1004CorrelationHole` | B |
| W1005 | `w1005_redner_reynolds.py` | `W1005RednerReynolds` | B |

連續播放：`main.py` 的 `Ch10Demo`（W1001–W1005）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch10/part10}`；內文 `notes/ch10/part10.tex`，圖表與數字檢查 `notes/ch10/make_figs.py`
（`check_numbers_ch10()` 逐一 assert，所有數字都從 `codes/utils/walk_cond.py` 算；表格與 TikZ／pgfplots 圖產生成 `notes/ch10/figs.tex`；巨集都以 `\Ten` 開頭、`\label` 用 `ten-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–9 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch10_formulas.md`、`walk_cond.py` 的 docstring、
`notes/ch9/part9.tex` 與 `notes/ch9/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 7–8 頁、圖表 5 組：(1) §10.1：4×4 環面上的三條例子（TikZ，座標由 `WK.torus_example(name)["segments"]` 產生）＋環面 vs 無限平面的 $c_N$ 表（$N=0..15$，分兩段）＋3×3、5×5 的「$N=L$ 少 4 條」與 Gray 碼 91392；
(2) §10.2：$n$ 的表（1.7183 … → 1.8682）＋$u=0$ 的核對＋漸近表（$N_0=100,1000,10^4$）＋極限順序的表（0.3299 … 0.9998）＋註 3 的對照（兩欄並排）＋$f(\lambda)$ 兩條曲線（pgfplots，各 ≤ 60 點）；
(3) §10.3：牆邊的 tanh（pgfplots）＋$X$、$X^{-1}$ 的例子（Fraction）＋模態表（$n=3,2,1,0$）＋(10.32) 的 $n=2$ 數值；
(4) §10.3 後半＋§10.4：$-D(\Delta)$ 剖面與 Coulomb 擬合（pgfplots 33 點）＋連續 vs 格子的表＋Ewald 常數＋sc 的步數 vs erfc 形的表（含書上引數的值）＋$G(R)\times2\pi R/3$；
(5) §10.3.1：$\rho(z)$ 三條曲線（pgfplots，各 60 點）＋$\rho$ 表＋漲落峰的表＋平均場 $\rho$、$\chi'$ 表＋Table 10.1（照抄）與三欄指數＋標度律驗算（含 (10.75) 的 $d_c\nu'$ 對照）。
小節對應書上 10.1、10.2、10.3、10.3 後半＋10.4、10.3.1，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch10_formulas.md` 的筆誤表，11 列；
另附〔注〕：(10.4) 的 $r$ 與第 9 章差 2 倍、(10.13) 的 $2\pi f(0)=2\pi$（$\frac1{2\pi}$ 相消）、(10.53) 多一個 $V/a^d$ 且隨 (10.15) 再乘 2、(10.55) 的 $\frac32\ln V$、
(10.58) 的「$\rho_{MF}\to|r|$」是 $\propto|r|$、p. 249 的 $\gamma'_{MF}=0$、(10.71) 少負號與 (10.72) 的 $|r|$、p. 246 的 $\gg1/|r|$ 應為 $\gg1/\sqrt{2|r|}$、(10.79) 的 $d^kk_1$、(10.80) 省了 $z_c^{-N}$、
p. 253 的轉折 $R\approx2\sqrt{N/z_c}$ 與前因子依正規化而定、(10.32)–(10.43) 只核對結構）。
最後練習 10.1–10.5 的解答：
10.1 $\vec S\to\vec S+\delta\vec S$：$\delta\mathcal H=\int d^dx\{(r\vec S+u\vec S|\vec S|^2-\vec h)\cdot\delta\vec S+\nabla\vec S\cdot\nabla\delta\vec S\}$，分部積分（週期邊界或牆上 $\delta\vec S=0$，邊界項為 0）⇒ $-\nabla^2\vec M+r\vec M+u\vec M|\vec M|^2-\vec h=0$；
離散化的數值梯度核對（相對誤差 $<10^{-7}$）；牆邊的解 $M=\sqrt{|r|/u}\tanh(x\sqrt{|r|/2})$（代回：$M''=-2k^2M(1-M^2/M_0^2)$、$k^2=|r|/2$）：0.6089、0.8884、0.9930，邊值數值解差 $2.8\times10^{-6}$；
10.2 $|\vec M+\vec\sigma|^4=(|\vec M|^2+2\vec M\cdot\vec\sigma+|\vec\sigma|^2)^2$ 的三階項 $4(\vec M\cdot\vec\sigma)|\vec\sigma|^2$ ⇒ $\mathcal H_3=u\int(\vec M\cdot\vec\sigma)|\vec\sigma|^2d^dx$；四階 $\frac u4\int|\vec\sigma|^4$；梯度項與 $r$ 項只到二階；數值內插誤差 $<10^{-11}$；
10.3 $X=\lambda_1I+uP$、$P=|\vec M\rangle\langle\vec M|$：$P\vec M=|\vec M|^2\vec M$ ⇒ 沿 $\vec M$ 的特徵值 $\lambda_1+u|\vec M|^2=\lambda_2$（1 重）；$\vec v\perp\vec M$ ⇒ $P\vec v=0$、特徵值 $\lambda_1$（$n-1$ 維）；
例子 $\frac12,\frac12,\frac32$；$\det X=\lambda_1^{n-1}\lambda_2$、$n\to0$ 時重數 $-1$；
10.4 $(\lambda_1I+uP)\big(\frac I{\lambda_1}-\frac{uP}{\lambda_1\lambda_2}\big)=I+uP\,\frac{\lambda_2-\lambda_1-u|\vec M|^2}{\lambda_1\lambda_2}=I$（用 $P^2=|\vec M|^2P$）；$P$ 與 $I$ 可交換所以左右逆相同；例子 $2I-\frac43P=\frac1{27}(\cdots)$，Fraction 精確；
10.5 令 $\vec q=\vec p/|\Delta|$：$\int d^dq\frac{e^{i\vec q\cdot\vec\Delta}}{q^2(2|r|+q^2)}=|\Delta|^{4-d}\int d^dp\frac{e^{i\vec p\cdot\hat\Delta}}{p^2(2|r||\Delta|^2+p^2)}\approx\frac{|\Delta|^{2-d}}{2|r|}\int d^dp\frac{e^{i\vec p\cdot\hat\Delta}}{p^2}$（$|\Delta|^2|r|\gg1$）⇒ (10.48)；
乘上 (10.46) 漏掉的 $2|r|$ 後與 $r$ 無關；三維精確 $D=\frac1{4\pi\rho}(1-e^{-\sqrt{2|r|}\rho})$：0.0602、0.0374、0.0198、0.0099（$r=-1$）；格子與週期盒子的數值見圖 (4)。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點、表的數字由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；講義裡「÷」可以用 `$\div$`；TikZ 從一點畫回同一點不要用 `to[looseness]`。
