# 《Elements of the Random Walk》第 11 章 → Manim：大綱與省 token 檢討（2026-10-19）

來源：Rudnick & Gaspari 第 11 章 Scaling, fractals, and renormalization，pp. 255–284（PDF 頁序 271–300）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_rg.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch10.md` §6–§9。

## 1. 草稿大綱（6 幕；本章 30 頁）

| 幕 | 主題 | 書上 |
|---|---|---|
| W1101 | 碎形：點和直線放大後不變、圓有自己的尺度；Koch 曲線（每一代長度 ×4/3 → 無限長）、Cantor 集（長度 → 0）；以端點為圓心的「質量－半徑」$n(R/3)=n(R)/4$ → 維度 $\ln4/\ln3$ (11.1)–(11.5)、複本公式 $d_f=\ln\mathcal R/\ln\mathcal M$ (11.6)（Ex. 11.1–11.3）；普通漫步的碎形維度是 2，自迴避漫步是 $1/\nu$ | §11.1–11.1.1，pp. 255–261，Fig. 11.1–11.5 |
| W1102 | 「被動」的重整化：把系統放大 $b$ 倍，Landau 自由能 (11.9) 各項乘 $b^d$ → $m\to mb^{d/4}$、$r\to rb^{d/2}$、$h\to hb^{3d/4}$ → $m=r^{1/2}f(h/r^{3/2})$ (11.10)–(11.17)；Gaussian 模型 (11.18)–(11.25) 的另一種組合，$d=4$ 時兩者相同（Ex. 11.4） | §11.1.2，pp. 261–264 |
| W1103 | 實空間重整化（一維 Gaussian 鏈）：固定偶數格點、把奇數格點積掉 (11.26)–(11.33)，重新標度自旋 (11.34)–(11.37)，$r$ 的遞迴（Fig. 11.6）與蛛網圖（Fig. 11.7）；小 $r$ 時 $r^{(n)}\approx4^nr$、$h^{(n)}$、自旋的總因子 (11.41)–(11.45) → $m$ 的標度形式、$\gamma=1$ (11.46)–(11.49)（Ex. 11.5–11.6） | §11.2 前半，pp. 264–269 |
| W1104 | 一維 Ising 的抽格：$2\cosh K(S_2+S_4)=Ae^{K'S_2S_4}$、$K'=\frac12\ln\cosh2K$ (11.50)–(11.55)（Ex. 11.7）、加外場（Ex. 11.8）、配分函數變成連乘 (11.57)–(11.63)；$x=e^{-2K}$ 的遞迴、不動點 $x=0,1$、$\xi\propto e^{2J/k_BT}$、$F\to-NJ$（Ex. 11.9） | §11.2 例題，pp. 269–272 |
| W1105 | 動量殼與連續重整化：布里淵區縮小一半（Fig. 11.8）、Gaussian 模型在 Fourier 空間各模態不耦合 (11.64)–(11.67)、$\vec q\to\vec q/2$、$N\to2^dN$、自旋重新標度 → $r'=b^2r$、$h'=b^{(d+2)/2}h$ (11.68)–(11.73)；讓 $b=e^{\delta\ell}$ 連續變化 → $\frac{dr}{d\ell}=2r$、$\frac{dh}{d\ell}=\frac{d+2}2h$ (11.74)–(11.83)，$m$ 的偏微分方程與特徵線 (11.84)–(11.87) | §11.2.1–11.2.2，pp. 272–277 |
| W1106 | 不動點與臨界指數：區塊自旋、參數的流 (11.88)–(11.92)、臨界超曲面（Fig. 11.9、11.10）、線性化與特徵值：相關／無關／邊際 (11.93)–(11.96)；$\nu=1/\Lambda_r$ (11.97)–(11.103)、自由能的標度 (11.104) 與四個指數公式 (11.105)–(11.108)；Gaussian 模型 $\Lambda_r=2$、$\Lambda_h=\frac{d+2}2$ (11.109)–(11.114)＋第 11 章回顧 | §11.3，pp. 277–284 |

## 2. 生動度檢討

- W1101 的「無限放大」：畫 Koch 曲線第 0–5 代，再把鏡頭放大到其中一段，看到同樣的形狀；Cantor 集一層層挖掉中間三分之一。質量－半徑用書上 Fig. 11.5 的兩個圓（半徑 $L$ 與 $L/3$）直接數段數。最後一段把第 7–9 章的自迴避漫步接回來：「漫步的碎形維度＝$1/\nu$」，普通漫步 2、二維自迴避 4/3。
- W1102 用「資料塌縮」（同 W802 的手法）：幾個 $r$ 的 Landau 解 $m(r,h)$ 畫成不同曲線，把軸換成 $m/r^{1/2}$ 對 $h/r^{3/2}$ 後疊成一條。
- W1103、W1104 是全章最直觀的兩幕：一排自旋，奇數格點「被積掉」（變淡消失），偶數格點之間的耦合變了；蛛網圖一步步畫，參數往不動點跑（Gaussian：$r$ 越來越大；Ising：$K\to0$）。
- W1104 的精確核對：小 $N$ 的 Ising 環用暴力加總、轉移矩陣 $(2\cosh K)^N+(2\sinh K)^N$、抽格的連乘三種算法得同一個 $Z$。
- W1105 用二維布里淵區的正方形一圈圈縮小，殼層塗色後拿掉；$(r,h)$ 平面上畫出 $r e^{2\ell}$、$he^{(d+2)\ell/2}$ 的流線。
- W1106 的流向圖：臨界超曲面上的點流向不動點、離開超曲面的方向被推開（相關）；最後用 (11.105)–(11.108) 把第 9 章 W903 的二維 Ising 指數（$\Lambda_r=1$、$\Lambda_h=\frac{15}8$）再算一次。
- 找錯遊戲照舊（基礎代理人校對後列出）。不要在畫面或講義引用 Swift／De Morgan 的詩句原文，只說「跳蚤身上還有更小的跳蚤」的意思。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W1101 | Koch 曲線前幾代的段數 $4^k$、長度 $(4/3)^k$；Cantor 集 $2^k$ 段、總長 $(2/3)^k$ | 盒子計數法（box counting）估 $\ln4/\ln3=1.2619$、$\ln2/\ln3=0.6309$；普通漫步的質量－半徑斜率 ≈ 2（固定種子） |
| W1102 | Landau 的 $m(r,0)=\sqrt{-r/u}$、$m(0,h)\propto h^{1/3}$ | 塌縮圖；Gaussian 模型的指數（Ex. 11.4）在 $d=3$：$\beta=\frac14$、$\gamma=1$、$\delta=5$ |
| W1103 | 三個自旋積掉中間那個：完成平方 (11.28)–(11.30) | 遞迴 $r''(r)$ 的數值與蛛網圖；小 $r$ 時 $4^n$ 的成長 |
| W1104 | 4 個自旋的 Ising 環：暴力加總 vs 轉移矩陣 vs 抽格 | $K'=\frac12\ln\cosh2K$ 迭代幾次的表、$x'=\frac{2x}{1+x^2}$、$\xi=-1/\ln\tanh K$ 對 $e^{2K}/2$ |
| W1105 | 一維 Gaussian 鏈的 $\alpha(q)=a(1-\cos q)$（書上 $-2a(1-\cos q)$ 是筆誤）與小 $q$ 近似 | $(r,h)$ 的流線；(11.24) 代進 (11.87) 驗證 |
| W1106 | $2\times2$ 的線性化矩陣與它的特徵值（相關／無關） | 指數表：Gaussian（$d=3,4$）、二維 Ising（$\Lambda_r=1$、$\Lambda_h=\frac{15}8$）對照第 9 章 |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_rg.py`（純 Python＋numpy，scipy 可用；可 import 同資料夾的 `walk_crit`、`walk_saw`、`walk_scale`，照 `walk_crit` 的相對 import 寫法）：
  Koch 曲線與 Cantor 集的頂點／區段（給場景畫圖用）、盒子計數與質量－半徑的維度估計、普通漫步的質量－半徑（固定種子）、自迴避的 $1/\nu$（引用第 7–9 章的數字）；
  Landau 解與塌縮、Gaussian 模型的標度與指數（Fraction）；一維 Gaussian 鏈的抽格（完成平方的數值核對、$r''(r)$、$h''$、自旋因子、蛛網圖的點列）；
  一維 Ising 的抽格（$A$、$K'$、加外場的 Ex. 11.8、連乘的 $Z$ 對暴力加總與轉移矩陣、$x=e^{-2K}$ 的遞迴與不動點、$\xi$、$F\to-NJ$）；
  動量殼與連續重整化的流（解析解 vs 數值積分）、(11.87) 的驗證；線性化矩陣的特徵值與指數公式 (11.105)–(11.108)（Gaussian、二維 Ising）。`check()` 用兩種以上獨立算法核對全部場景卡數字，約 10 秒內。
- **重用（不改）**：`walk_draw`、`lattice_draw`、`walk_crit`（指數表、二維 Ising 的精確值）。
- **畫圖**：不新增模組；Koch 曲線、自旋鏈、蛛網圖、布里淵區、流向圖寫在場景檔的模組層級函式。
- combkit：`walk_rg.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex24_decimation_flow.py`；升 1.21.0。

## 5. 省 token 檢討

1. **30 頁 → 6 幕**：§11.3 的一般理論只取主軸（不動點 → 線性化 → 特徵值 → 指數），(11.97)–(11.103) 的推導細節放講義；(11.84)–(11.87) 只做「代入驗證」，特徵線法放講義。
   場景代理人 A 做 W1101–W1103、B 做 W1104–W1106。
2. **重活交給上下文短的代理人**：基礎代理人、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁；`ch11_formulas.md` 只收場景與講義用得到的式子＋筆誤表（我讀正文時看到的可疑處：p. 267 的「$t''$」、(11.36)(11.37) 與 Fig. 11.6 的座標 $r/a$、(11.41) 的小 $r$ 近似、(11.43) 與 (11.45) 的因子、(11.54)(11.55)、(11.71)(11.72) 的 4 與 $b^2$、(11.73) 的 $b^{(2-d)/2}$、(11.82)(11.83)、(11.105)–(11.108) 與 (11.111)–(11.114)；Ex. 11.9 的題目文字斷掉）。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_rg.check()` assert；模擬固定種子、在模組裡算好（有快取）。
5. **不新增畫圖模組、不做 3D**。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 11 節目標 9–10 頁**；練習解答 Ex. 11.1–11.9。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log、記憶；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_rg.py` 已完成（`python utils/walk_rg.py` 自我測試通過，約 1 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch11_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：碎形維度 $p$（＝$d_f$）；Landau 的 $r,h,m$（本系列用第 9 章的 $\frac r2m^2+\frac u4m^4-hm$，$u=1$）；
一維 Gaussian 鏈的耦合 $a$（$=1$）；Ising 的 $K=J/k_BT$、外場 $H=h/k_BT$；抽格 $b=2$、連續 $b=e^\ell$；特徵值 $\Lambda_r,\Lambda_h$。
不要引用 Swift／De Morgan 的詩句原文，只說「跳蚤身上還有更小的跳蚤」的意思。
程式代號：`RG` = `utils.walk_rg`、`WC` = `utils.walk_crit`（只在 W1106 用 `WC.EXPONENTS`）、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。

**W1101 碎形：放大之後還一樣**（§11.1–11.1.1，pp. 255–261，(11.1)–(11.6)，Fig. 11.1–11.5，Ex. 11.1–11.3）
1. 問：什麼東西放大之後還長得一樣？點、直線（不變，平凡）；圓弧越放大越直（有自己的尺度 = 半徑）。字幕：跳蚤身上還有更小的跳蚤——一層又一層（不引用原詩）。無數字。
2. Koch 曲線第 0–5 代（`RG.koch_points(g)`，每代一個 VMobject，第 5 代 1025 點；突起在左邊＝畫成向上）：表（`RG.koch_counts(g)`）代 $k=0..5$：段數 1、4、16、64、256、1024；
   總長 $(\frac43)^k$：1、1.3333、1.7778、2.3704、3.1605、4.2140 → ∞。放大：取第 5 代的前四分之一（`koch_points(5)[:257]`）放大 3 倍，和整條疊在一起——一模一樣。
   Cantor 集第 0–5 代（`RG.cantor_intervals(g)`，最多 32 段，一列一代）：段數 $2^k$、總長 $(\frac23)^k$：1、0.6667、0.4444、0.2963、0.1975、0.1317 → 0。
3. 質量－半徑（Fig. 11.5）：以 Koch 曲線左端為圓心畫半徑 $L$、$L/3$ 的圓；(11.1) $n(R/3)=\frac14n(R)$（`RG.mass_table('koch')`）：$R=1,\frac13,\frac19,\frac1{27},\frac1{81}$：
   $n$＝1、0.2500、0.0625、0.0156、0.0039；(11.2)–(11.5) $n=KR^p$ ⇒ $3^{-p}=\frac14$ ⇒ $p=\frac{\ln4}{\ln3}=1.2619$（`RG.LN4_LN3`）。Cantor（`mass_table('cantor')`）：1、0.5000、0.2500、0.1250、0.0625 ⇒ $\frac{\ln2}{\ln3}=0.6309$（Ex. 11.2）。
   一行：盒子計數（`RG.box_counts('koch', RG.BOX_K, 3)`）$\varepsilon=\frac19,\frac1{27},\frac1{81},\frac1{243}$：22、83、325、1291 個盒子 → 斜率 1.2362（`RG.box_dim('koch')`）；Cantor 0.6254。
4. 複本公式 (11.6) $p=\ln\mathcal R/\ln\mathcal M$（`D.tex_table`，5 列 × 4 欄，`RG.REPLICAS`、`RG.replica_dim`）：線段 (3, 3) 1.0000、正方形 (9, 3) 2.0000（Ex. 11.1）、
   Koch (4, 3) 1.2619、Cantor (2, 3) 0.6309（Ex. 11.3）、Sierpinski (3, 2) 1.5850（本系列加的）。字幕：維度介於 1 和 2 之間——介於線和面之間。
5. 漫步也是碎形：一條普通漫步（`RG.rw_trail()`，1501 點，一個 VMobject，種子 1101）；$N$ 步走 $R\sim N^{1/2}$ ⇒ $N\propto R^2$。$\langle R_g^2\rangle=\frac{N(N+2)}{6(N+1)}$（`RG.rg2_exact`）；
   400 條（`RG.rw_df()`）$\ln N$ 對 $\ln R_g$ 的斜率 2.0258、公式 2.0069 → 2。自迴避 $d_f=1/\nu$（`RG.saw_df(d)`）：二維 $\frac43$（精確列舉的有效值 $N=10,12,14$：1.3848、1.3788、1.3740，`RG.saw_df_eff`）、
   三維 1.70（引用 $\nu=0.588$）、四維 2。
收束：碎形維度＝「放大 $\mathcal M$ 倍、份量變 $\mathcal R$ 倍」的比的對數——重整化群就是把這個想法用在漲落的系統上。筆誤框 0 個。

**W1102 被動的重整化：Landau 與 Gaussian**（§11.1.2，pp. 261–264，(11.7)–(11.25)，Ex. 11.4）
1. 問：把整個磁鐵放大 $b$ 倍，自由能的每一項會怎樣？(11.9) $\mathcal H=V[rm^2+um^4-hm]$，$V\to b^dV$ ⇒ $m'=mb^{d/4}$、$r'=rb^{d/2}$、$h'=hb^{3d/4}$（三個箭頭各自變色）。
   字幕：本系列照第 9 章寫 $\frac r2m^2+\frac u4m^4$——常數不改變標度。
2. (11.11)–(11.16)：$m(r,h)=b^{-d/4}F(rb^{d/2},hb^{3d/4})$ 與 $b$ 無關 ⇒ $m=r^{1/2}f(h/r^{3/2})$；兩個極限（`RG.landau_m`）：$m(-0.25,0)=0.5000=\sqrt{-r/u}$；
   $m(0,h)=h^{1/3}$：$h=0.008\to0.2000$、$h=0.125\to0.5000$。
   筆誤框 1（where「p. 263，(11.16) 後」）：錯 `f(h/r^{3/2})\to h^{1/3}\,r`、對 `f(h/r^{3/2})\to h^{1/3}\,r^{-1/2}`；
   字幕：$r=10^{-4}$、$h=10^{-3}$：$f=9.9667$，更正式 10.0000，書上的只有 $10^{-5}$（`RG.book1116_limit(1e-4, 1e-3)`）。
3. 資料塌縮（`D.plot_axes` 兩張並排，各 3 條 ≤ 40 點）：左 $m$ 對 $h$（`RG.landau_curve(r)`，$r=0.25,0.5,1$，$h\in[0,2]$）；右 $m/r^{1/2}$ 對 $h/r^{3/2}$（`RG.collapse_curve(r)`，只畫 $x\le4$）——疊成一條 $y+y^3=x$。
   $f(1)=0.6823$、$f(2)=1$、$f(10)=2$（`RG.landau_scaled`）；三個點 $(r,h)=(0.25,0.25),(0.5,0.7071),(1,2)$（`RG.COLLAPSE_POINTS`）：$m=0.5000,0.7071,1.0000$——都在 $x=2$、$y=1$。
   (11.17) $\partial m/\partial h=r^{-1}f'(0)$ ⇒ $\gamma=1$。
4. Gaussian 模型 (11.18)–(11.25)：$d^dx\to b^dd^dx$、$\nabla\to b^{-1}\nabla$ ⇒ $S(b\vec x)=b^{(2-d)/2}S'$、$r=b^{-2}r'$、$h=b^{-(d+2)/2}h'$ ⇒ $m=r^{(d-2)/4}F\big(h/r^{(d+2)/4}\big)$。
   Ex. 11.4 表（`D.tex_table`，`RG.gauss_exponents(d)`，Fraction）：$d=3$：$\beta=\frac14$、$\gamma=1$、$\delta=5$；$d=4$：$\frac12$、1、3（＝平均場）。
   $d=4$ 時兩種組合相同：$\frac{d+2}4=\frac32$、$\frac{d-2}4=\frac12$（`RG.scaling_combos(4)`；$d=3$ 是 $\frac54$、$\frac14$）——第 8 章說過 $d=4$ 的特殊地位。
收束：只看「放大之後參數怎麼變」，就讀出了 $\gamma=1$——但這兩個模型還是「被動」的：漲落沒有被積掉。筆誤框 1 個。

**W1103 實空間重整化：一維 Gaussian 鏈**（§11.2 前半，pp. 264–269，(11.26)–(11.49)，Fig. 11.6–11.7，Ex. 11.5–11.6）
1. 問：把一半的自旋積掉，剩下的系統還長得一樣嗎？一排 9 個自旋（連續值，畫成長短不一的箭頭），奇數格點變淡消失、偶數格點之間的彈簧變了。
   (11.27) $\mathcal H=\sum_i[\frac r2S_i^2+\frac a2(S_i-S_{i+1})^2]$。
2. (11.28)–(11.30) 完成平方：$\int dS_3\exp[-(\frac r2+a)S_3^2+aS_3(S_2+S_4)-\frac a2(S_2^2+S_4^2)]=\sqrt{\frac{2\pi}{r+2a}}\exp[-\frac{a^2(S_2-S_4)^2}{2(r+2a)}-\frac{ar(S_2^2+S_4^2)}{2(r+2a)}]$；
   數值（`RG.square_check()`，$r=0.5$、$a=1$、$S_2=0.3$、$S_4=-0.7$）：數值積分 1.224819 ＝ 閉式 1.224819。(11.32)(11.33)：$r'=\frac{r(r+4a)}{r+2a}$、$a'=\frac{a^2}{r+2a}$。
   筆誤框 1（where「p. 265，(11.31)」）：錯 `\Big(\frac{2\pi}{r+2a}\Big)^{N/2}`、對 `\Big(\frac{2\pi}{r+2a}\Big)^{N/4}`；
   字幕：$N/2$ 個奇數格點，每個只給一個平方根；8 個自旋的環用行列式算：照書上 $\ln Z$ 多 1.8432（`RG.ring_decimation_check(8, 0.5, 1.0, 0.0, book=True)`）。
3. (11.34)–(11.38)：$S_{2i}\to S_{2i}\sqrt{(r+2a)/a}$ ⇒ $r''=\frac{r(r+4a)}a$、$a''=a$、$h''=h\frac{r+4a}{\sqrt{a(r+2a)}}$。Fig. 11.6–11.7（`D.plot_axes` $[0,10]^2$；曲線 `RG.rr_curve()` 40 點，$r\in[0,1.7417]$；
   虛線 $r''=r$；蛛網 `RG.cobweb_points()` 6 點逐段畫）：$r_0=0.05\to0.2025\to0.8510\to4.1282$（下一步 33.5553 出圖）——越跑越大。不動點 $r=0$（斜率 4，不穩定）與 $-3a$（不物理）。
   筆誤框 2（where「p. 267，Fig. 11.6–11.7」）：錯 `r/\sqrt{a},\ r''/\sqrt{a}`、對 `r/a,\ r''/a`；字幕：$r$ 和 $a$ 同單位，$r''/a=(r/a)(r/a+4)$；正文的「$t''$」也是 $r''$。
4. 小 $r$（(11.41)–(11.45)）：表（`RG.small_r_table()`，$r=h=10^{-4}$，$k=1..5$）$r^{(k)}/(4^kr)$：1.0000、1.0001、1.0005、1.0021、1.0086；$h^{(k)}/(2^{3k/2}h)$：全是 1.0000；
   自旋總因子 ÷ $2^{k/2}$：1.0000、1.0001、1.0005、1.0021、1.0085 ⇒ $r^{(n)}=4^nr$、$h^{(n)}=2^{3n/2}h$、$b^{(n)}=2^{n/2}$（$h''/h\to2^{3/2}=2.8284$，`RG.hh_gauss(0.0, 1.0)`）。
5. (11.46)–(11.49)：$m=2^{n/2}f(r2^{2n},h2^{3n/2})$ 與 $n$ 無關 ⇒ $f=r^{-1/4}F(h/r^{3/4})$ ⇒ $\partial m/\partial h=r^{-1}F'$：$\gamma=1$，正是 (11.24) 取 $d=1$。
   筆誤框 3（where「p. 269，(11.46)(11.47)」）：錯 `m=2^{-n/2}f\big(r2^{2n},h2^{3n/2}\big)`、對 `m=2^{+n/2}f\big(r2^{2n},h2^{3n/2}\big)`；
   字幕：精確解 $m=h/r$（$r=0.01$、$h=0.001$、$n=1$）：照書上 0.0500，更正 0.1000＝$h/r$（`RG.book1147_m(0.01, 0.001, 1)`）。
收束：積掉一半、重新標度，模型的樣子不變，只有 $r$ 變成 $4r$——由「每次乘 4」讀出 $\nu=\frac12$、$\gamma=1$。筆誤框 3 個（$t''$ 用字幕）。

**W1104 一維 Ising 的抽格**（§11.2 例題，pp. 269–272，(11.50)–(11.63)，Ex. 11.7–11.9）
1. 問：Ising 鏈的奇數自旋加總掉，剩下的自旋之間是什麼耦合？一排 9 個 ±1 自旋（上下箭頭），奇數格點淡出。(11.51)(11.52) $\sum_{S_3=\pm1}e^{KS_3(S_2+S_4)}=2\cosh K(S_2+S_4)$：
   兩種情況（$K=1$，`RG.book1153(1.0)[2]`）：$S_2=S_4$：$2\cosh2K=7.5244$；$S_2=-S_4$：2。
   筆誤框 1（where「p. 270，(11.53)」）：錯 `\cosh K(S_2+S_4)=2\cosh 2K\ (S_2=S_4)`、對 `2\cosh K(S_2+S_4)=2\cosh 2K\ (S_2=S_4)`；字幕：$K=1$：$\cosh2=3.7622$，不是 7.5244。
2. Ex. 11.7 (11.54)(11.55)：$Ae^{K'}=2\cosh2K$、$Ae^{-K'}=2$ ⇒ $A=2\sqrt{\cosh2K}=3.8793$、$K'=\frac12\ln\cosh2K=0.6625$（$K=1$，`RG.ising_A`、`RG.ising_Kp`）。
   (11.62) 迭代表（`RG.ising_iterate(1.0, 5)`）：$K_n$＝1.0000、0.6625、0.3501、0.1137、0.0128、0.0002 → 0（越看越像高溫）。
   Ex. 11.8 一行（`RG.ex118(1.0, 0.1)`）：$K=1$、$H=0.1$：$A=3.8897$、$K_1'=0.6602$、$K_2'=0.0482$、新外場 $H'=H+2K_2'=0.1964$。
3. (11.57)–(11.63) 配分函數變成連乘：表（`D.tex_table`，3 列）$K=1$：4 個自旋的環：暴力 $2^4$ 項（`RG.ring_Z_brute(4, 1.0)`）＝ 轉移矩陣 $(2\cosh K)^4+(2\sinh K)^4$ ＝ 抽格連乘 ＝ 121.2329
   （字幕：和第 9 章 W902 的 4 自旋環同一個數）；8 個自旋：9160.4387。
   筆誤框 2（where「p. 271，(11.63)」）：錯 `\cdots A(K_j)^{N/2^{j+1}}\cdots`、對 `\cdots A(K_j)^{N/2^{j}}\cdots`；字幕：照書上 18.5498、143.2705（`RG.ring_Z_decimation(N, 1.0, book=True)`）。
4. Ex. 11.9(a)：$x=e^{-2K}$ ⇒ $x'=\frac{2x}{1+x^2}$；蛛網（`D.plot_axes` $[0,1]^2$，曲線 ≤ 40 點、對角虛線；點列由 `RG.x_iterate(RG.x_of_K(1.0), 5)` 在場景裡組）：
   0.1353、0.2658、0.4965、0.7966、0.9747、0.9997 → 1。不動點 $x=0$（$T=0$，斜率 2，不穩定）、$x=1$（$T=\infty$，斜率 0，穩定）（`RG.x_slope`）。
5. Ex. 11.9(b)(c)：$x$ 小時 $x'\approx2x$、$\xi'=\xi/2$ ⇒ $\xi\propto1/x=e^{2J/k_BT}$。表（`RG.xi_exact(K)`、`RG.xi_ratio(K)`）$K=1,2,3$：$\xi=3.6719,27.2960,201.7140$；
   $e^{2K}/2=3.6945,27.2991,201.7144$；比 0.4969、0.4999、0.5000。字幕：其實 $\tanh K'=\tanh^2K$，$\xi$ 每次剛好減半。
   $\ln Z/(NK)$（`RG.free_energy_ratio(K)`）$K=1,2,5$：1.1269、1.0091、1.0000 ⇒ $F\to-NJ$。
收束：一維 Ising 的流只有兩個不動點，$T=0$ 不穩定、$T=\infty$ 穩定——所以一維沒有有限溫度的相變。筆誤框 2 個。

**W1105 動量殼與連續重整化**（§11.2.1–11.2.2，pp. 272–277，(11.64)–(11.87)，Fig. 11.8）
1. 問：在 Fourier 空間「積掉一半」長什麼樣？Fig. 11.8：格距 $c\to2c$ ⇔ 布里淵區 $Q\to Q/2$（`L.bz_frame` 或自畫正方形；外圈殼層塗色後淡出）；
   拿掉的比例 $1-2^{-d}$（`RG.shell_fraction(d)`）：$d=2$：$\frac34$、$d=3$：$\frac78$。
2. (11.64)–(11.67)：$S_i=N^{-1/2}\sum_qS(q)e^{iqR_i}$、$\mathcal H=\sum_q[\frac r2+\alpha(q)]S(q)S(-q)-\sqrt NS(0)h$，各模態不耦合。
   表（`RG.alpha_table()`，$N=8$，$k=0..4$）：$q$＝0、0.7854、1.5708、2.3562、3.1416；$\alpha=a(1-\cos q)$＝0、0.2929、1.0000、1.7071、2.0000；$\frac12aq^2$＝0、0.3084、1.2337、2.7758、4.9348（小 $q$ 才像）。
   筆誤框 1（where「p. 273，(11.66)(11.67)」）：錯 `\alpha(q)=-2a(1-\cos qc)\sim ac^2q^2`、對 `\alpha(q)=a(1-\cos qc)\sim\frac12ac^2q^2`；
   字幕：隨機 8 個自旋：實空間 3.0394＝Fourier 3.0394；照書上的 $\alpha$ 得 $-4.3508$（負的能量，`RG.fourier_check()`）。
3. (11.68)–(11.73)：$\vec q\to\vec q/2$（$Jq^2\to Jq^2/4$）、$N\to2^dN$、$S\to2S$ ⇒ $r\to4r$、$h\to2^{(d+2)/2}h$（$d=3$：×5.6569，`RG.shell_rg(1.0, 1.0, 3)`）；換成 $b$：$m(r,h)=b^{(2-d)/2}m(rb^2,hb^{(d+2)/2})$——和 W1102 的 (11.22) 一樣。
4. (11.74)–(11.83) 連續：$b=e^{\delta\ell}$：$\vec q\to\vec q(1-\delta\ell)$、$S\to(1+\delta\ell)S$、$N\to N(1+d\,\delta\ell)$ ⇒ $\frac{dr}{d\ell}=2r$、$\frac{dh}{d\ell}=\frac{d+2}2h$。
   $(r,h)$ 平面的流線（`RG.flow_line(r0, h0, 3, ell_max)`，起點 $(\pm0.01,0.01)$、$(\pm0.03,0.003)$，$\ell_{\max}$ 取到 $|r|\le1$ 或 $h\le1.5$；每條 ≤ 40 點，加箭頭）：$h\propto|r|^{5/4}$，都離開原點。
   表（`RG.flow(0.01, 0.01, 3, ell)`）：$\ell=1$：$r=0.0739$、$h=0.1218$；$\ell=2$：0.5460、1.4841（解析＝數值積分）；$\ell=\ln2$ 就是一次 $b=2$。
5. (11.84)–(11.87)：$S_i=e^{\frac{2-d}2\ell}S_i^{(\ell)}$ ⇒ $m=e^{\frac{2-d}2\ell}f(r(\ell),h(\ell))$ 與 $\ell$ 無關 ⇒ $0=\big(-\frac{d-2}2+2r\partial_r+\frac{d+2}2h\partial_h\big)m$；
   把 (11.24) 代進去：三項一一相消（打勾；`RG.pde_residual(d)` 在 $d=1,3,4$ 都 $<10^{-7}$，卡上只寫「$<10^{-7}$」）；字幕：Landau 的組合 $r^{1/2}f(h/r^{3/2})$ 放進 $d=3$ 就不成立。
收束：讓尺度連續地變，重整化變成一組微分方程——參數的「流」。筆誤框 1 個。

**W1106 不動點與臨界指數**（§11.3，pp. 277–284，(11.88)–(11.114)，Fig. 11.9–11.10）
1. 問：參數一直流——會流到哪裡停下來？區塊自旋 (11.88)–(11.92)：$\kappa_i\to R_i(\kappa)$、$\frac{d\kappa_i}{d\ell}=W_i(\kappa)$、$s=e^\ell$、$\frac{d\xi}{d\ell}=-\xi$；
   臨界超曲面（$\xi=\infty$）被映到自己；不動點 $W_i(\kappa^*)=0$（Fig. 11.9 示意：曲面上的點流向一個黑點；無數字）。
2. 線性化 (11.93)–(11.96)：$\frac{d\delta\kappa_i}{d\ell}=\sum_jA_{ij}\delta\kappa_j$ ⇒ $\delta\kappa=\sum_\alpha C_\alpha e^{\Lambda_\alpha\ell}$。示意矩陣（`RG.TOY_A`、`RG.toy_eigen()`，Fraction）
   $A=\begin{pmatrix}1&2\\2&-2\end{pmatrix}$：特徵值 2（$(2,1)$，相關）、$-3$（$(1,-2)$，無關）；流向圖（`RG.toy_flow_line(c1, c2)`，`RG.TOY_STARTS` 5 條，各 30 點；
   $c_1=0$ 那條沿臨界線流進不動點，其餘先靠近再沿 $(2,1)$ 被推開）。註 6：$\Lambda>0$ 相關、$<0$ 無關、$=0$ 邊際。
3. (11.97)–(11.103)：偏離 $\propto re^{\Lambda_r\ell}$、$\xi^{(\ell)}=\xi e^{-\ell}$；到 $\xi^{(\ell^*)}=1$ 時離臨界面距離 $D$：$Kre^{\Lambda_r\ell^*}=D$ ⇒ $\xi=(Kr/D)^{-1/\Lambda_r}$ ⇒ $\nu=1/\Lambda_r$（Fig. 11.10 示意）。
   數值（`RG.xi_from_flow(0.01, Lr)`，$K=D=1$）：$\Lambda_r=2$：$\ell^*=2.3026$、$\xi=10.0000$；$\Lambda_r=1$：4.6052、100.0000。
4. (11.104)–(11.108)：$f(r,h)=e^{-\ell d}f(re^{\Lambda_r\ell},he^{\Lambda_h\ell})$ ⇒ $\beta=\frac{d-\Lambda_h}{\Lambda_r}$、$-\gamma=\frac{d-2\Lambda_h}{\Lambda_r}$、$2-\alpha=\frac d{\Lambda_r}$、$\delta=\frac{\Lambda_h}{d-\Lambda_h}$、$\nu=\frac1{\Lambda_r}$。
   表（`D.tex_table`，6 列 × 5 欄，`RG.rg_exponents(*RG.RG_CASES[k])`，Fraction）欄：Gaussian $d=3$（$\Lambda_r=2$、$\Lambda_h=\frac52$）｜Gaussian $d=4$（2、3）｜二維 Ising（1、$\frac{15}8$）｜二維自迴避（$\frac43$、$\frac{91}{48}$）：
   $\alpha$：$\frac12$、0、0、$\frac12$；$\beta$：$\frac14$、$\frac12$、$\frac18$、$\frac5{64}$；$\gamma$：1、1、$\frac74$、$\frac{43}{32}$；$\delta$：5、3、15、$\frac{91}5$；$\nu$：$\frac12$、$\frac12$、1、$\frac34$。
   打勾：後三欄＝第 9 章 W903 的指數表（`WC.EXPONENTS` 的 'mf'、'ising2d'、'saw2d'）；每欄 $\alpha+2\beta+\gamma=2$。一行：W1103 的鏈 $\log_2 4=2$、$\log_2 2^{3/2}=\frac32$（`RG.chain_lambdas()`：2.0000、1.5000）＝ $d=1$ 的 $\Lambda_r$、$\Lambda_h$。
5. 第 11 章回顧（換頁；表六列，`D.tex_table`＋中文標籤）：W1101 Koch $\frac{\ln4}{\ln3}=1.2619$、Cantor 0.6309、漫步 2、自迴避 $\frac43$；W1102 $m=r^{1/2}f(h/r^{3/2})$、$d=4$ 與 Gaussian 一致；
   W1103 $r''=r(r+4a)/a\approx4r$、$h''\approx2^{3/2}h$；W1104 $K'=\frac12\ln\cosh2K\to0$、$x'=\frac{2x}{1+x^2}$、$\xi\propto e^{2K}$；W1105 $\frac{dr}{d\ell}=2r$、$\frac{dh}{d\ell}=\frac{d+2}2h$；W1106 $\nu=1/\Lambda_r$、$\beta=(d-\Lambda_h)/\Lambda_r$。
收束：臨界指數只由不動點附近兩個相關特徵值決定——下一章把這套用到 O($n$) 模型與自迴避漫步。筆誤框 0 個。

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
  from utils import walk_rg as RG
  from utils import walk_crit as WC          # 只有 W1106 需要
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_rg`（W1106 的對照用 `WC.EXPONENTS`）算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串，
  分數用 `Fraction`；Fraction 不能直接 `:.4f`，要 `float()`）；版面用 `assert_in_frame`／`assert_apart`。
  不要在場景裡呼叫 `RG.check()`（整套自我測試由 `check_env.py` 跑）；`RG.summary()` 約 1 秒，可以不用（直接呼叫卡上的函式）。
  隨機一律用模組的預設種子（普通漫步 1101、Fourier 核對 1105）。最慢的呼叫：`RG.box_counts('koch', …)` 約 0.3 秒、`RG.rw_df()` 第一次約 0.2 秒（快取）；其餘都很快。
  scipy 的 IntegrationWarning 不影響結果（`check_math` 裡用 `warnings.catch_warnings()` 包起來）。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W1101 0 個、W1102 1 個、W1103 3 個、W1104 2 個、W1105 1 個、W1106 0 個。
  框裡的 Tex 不放中文（英文字用 `\text{}`）。卡上標「用字幕」的書上問題只用一行字幕。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。字幕裡的「÷」寫成「除以」。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；矩陣用 `\begin{pmatrix}…\end{pmatrix}`。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「1 除以 ν」之類）。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——不要拿來當變數（Koch 的端點距離用 `Lk`、`size`；**(11.101) 的距離 $D$ 在程式裡叫 `Dist`**）；`RG`、`WC` 也是。
- 曲線一條＝一個 VMobject（`D.curve` 或 `set_points_as_corners`）：Koch 第 5 代 1025 點、漫步 1501 點是例外（各一個 VMobject，可以）；其他曲線每條 ≤ 60 點，刻度數字自己加。
  Cantor 集每代一列 `Line`（最多 32 段）。蛛網圖：轉角點逐段 `ShowCreation`（W1103 用 `RG.cobweb_points()`；W1104 的 $x$ 蛛網在場景裡由 `RG.x_iterate` 組點）。
  自旋鏈：9 個格點、箭頭長度（Gaussian）或上下（Ising）；奇數格點淡出（`FadeOut` 或 `set_opacity`）。布里淵區的殼層、流向圖、Fig. 11.9–11.10 都是平面示意，不做 3D。
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
| W1101 | `w1101_fractals.py` | `W1101Fractals` | A |
| W1102 | `w1102_passive_rg.py` | `W1102PassiveRG` | A |
| W1103 | `w1103_gaussian_decimation.py` | `W1103GaussianDecimation` | A |
| W1104 | `w1104_ising_decimation.py` | `W1104IsingDecimation` | B |
| W1105 | `w1105_momentum_shell.py` | `W1105MomentumShell` | B |
| W1106 | `w1106_fixed_points.py` | `W1106FixedPoints` | B |

連續播放：`main.py` 的 `Ch11Demo`（W1101–W1106）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch11/part11}`；內文 `notes/ch11/part11.tex`，圖表與數字檢查 `notes/ch11/make_figs.py`
（`check_numbers_ch11()` 逐一 assert，所有數字都從 `codes/utils/walk_rg.py` 算（W1106 的對照可用 `walk_crit.EXPONENTS`）；表格與 TikZ／pgfplots 圖產生成 `notes/ch11/figs.tex`；
巨集都以 `\Eleven` 開頭、`\label` 用 `eleven-…`；執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–10 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch11_formulas.md`、
`walk_rg.py` 的 docstring、`notes/ch10/part10.tex` 與 `notes/ch10/make_figs.py` 的開頭（格式範例）；不讀原文；不引用 Swift／De Morgan 的詩句。
目標 9–10 頁、圖表 6 組：(1) §11.1.1：Koch 第 0–4 代與 Cantor 第 0–5 代（TikZ，座標由 `RG.koch_points(g)`、`RG.cantor_intervals(g)` 產生；第 4 代 257 點）＋段數／總長表＋
質量－半徑表（Koch、Cantor）＋盒子計數（22、83、325、1291 → 1.2362；Cantor 0.6254；$R=2^{-j}$ 的 1.2571）＋複本公式表（5 列）＋普通漫步的 $\langle R_g^2\rangle$ 表（`RG.rw_rg2(N)` vs `RG.rg2_exact(N)`，
$N=64..4096$）與斜率 2.0258／2.0069＋自迴避 $d_f$（$\frac43$、有效值、1.70、2）；
(2) §11.1.2：塌縮兩張（pgfplots，`RG.landau_curve`、`RG.collapse_curve`，各 3 條 ≤ 40 點）＋$f(1),f(2),f(10)$＋(11.16) 後的筆誤數值＋Ex. 11.4 指數表（$d=3,4$）與 `RG.scaling_combos`；
(3) §11.2 前半：$r''(r)$ 曲線與蛛網（pgfplots）＋小 $r$ 表（$k=1..5$）＋(11.31) 的行列式核對（多 1.8432）＋(11.47) 的數值（0.0500 vs 0.1000）＋`RG.chain_lambdas()`；
(4) §11.2 例題：$K_n$ 迭代表＋Ex. 11.8 的數值＋$Z$ 的表（$N=4,8,16$：暴力＝轉移矩陣＝連乘，書上指數的值）＋$x$ 的蛛網（pgfplots）＋$\xi$ 表＋$\ln Z/(NK)$ 表（含 $N=16$ 的環 1.1277）；
(5) §11.2.1–11.2.2：布里淵區的殼層（TikZ）＋$\alpha(q)$ 表（含書上的值）＋`RG.fourier_check()`＋$(r,h)$ 流線（pgfplots，$d=3$）＋流的數值表；
(6) §11.3：示意矩陣的流向圖（pgfplots，`RG.toy_flow_line`）＋$\xi$ 的 $\ell^*$ 例子＋指數表（4 欄，同 W1106 卡）＋與 `WC.EXPONENTS` 的對照。
小節對應書上 11.1.1、11.1.2、11.2 前半、11.2 例題、11.2.1–11.2.2、11.3，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch11_formulas.md` 的筆誤表，8 列；
另附〔注〕：(11.9) 沒有 $\frac12,\frac14$、(11.13) 是 $|r|^{1/2}$、(11.24) 少「$m=$」、(11.18) 梯度項的約定、$r''=r$ 的第二根 $-3a$、Ex. 11.8 的 $K_1,K_2$ 撞名、
Ex. 11.9 的「such that that」與「$E=-NJ$」、(11.78) 的 $N$、(11.98) 的 $C_{i,t}$ 與 $y$、Koch 的 $n(L/3)=\frac14n(L)$ 精確（30° 錐））。
最後練習 11.1–11.9 的解答：
11.1 平面上以面內一點為圓心：$n(R)=\pi R^2$（在邊上 $\frac\pi2R^2$、角上 $\frac\pi4R^2$）⇒ $n(R/3)=\frac19n(R)$ ⇒ $p=\ln9/\ln3=2$；
11.2 Cantor：以左端為中心，$R=\frac13$ 內正好是左邊那一份 ⇒ $n(R/3)=\frac12n(R)$ ⇒ $3^{-p}=\frac12$、$p=\ln2/\ln3=0.6309$（數值表 1、0.5000、0.2500、0.1250、0.0625）；
11.3 Koch：$\mathcal R=4$、$\mathcal M=3$；Cantor：$\mathcal R=2$、$\mathcal M=3$ ⇒ (11.6) 給同樣的值；盒子計數 1.2362、0.6254 是第三種估計（有限尺寸偏低）；
11.4 由 (11.24)：$h=0$、$r<0$ 形式上 $m\propto|r|^{(d-2)/4}$ ⇒ $\beta=\frac{d-2}4$；$\chi=r^{-1}F'(0)$ ⇒ $\gamma=1$；$r=0$ 時 $m$ 與 $r$ 無關要求 $F(x)\propto x^{(d-2)/(d+2)}$ ⇒ $m\propto h^{(d-2)/(d+2)}$、
$\delta=\frac{d+2}{d-2}$；$d=3$：$\frac14,1,5$；$d=4$：$\frac12,1,3$（〔注〕Gaussian 模型在 $r<0$ 不穩定、$m=h/r$，$\beta$ 只是標度形式的讀法）；
11.5 (11.28)→(11.29) 展開；$\int e^{-AS^2+BS}dS=\sqrt{\pi/A}e^{B^2/4A}$，$A=\frac{r+2a}2$、$B=a(S_2+S_4)$；用 $(S_2+S_4)^2=2(S_2^2+S_4^2)-(S_2-S_4)^2$ 得 (11.30)；
每個偶數格點從兩邊各得 $-\frac{ar}{2(r+2a)}S^2$ ⇒ $\frac{r'}2=\frac r2+\frac{ar}{r+2a}$ ⇒ $r'=\frac{r(r+4a)}{r+2a}$、$a'=\frac{a^2}{r+2a}$；前因子 $(\frac{2\pi}{r+2a})^{N/4}$（書上 $N/2$）；數值 1.224819；
11.6 $B=a(S_2+S_4)+h$：$\frac{B^2}{2(r+2a)}$ 的一次項 $\frac{ah(S_2+S_4)}{r+2a}$ ⇒ $h'=h+\frac{2ah}{r+2a}=h\frac{r+4a}{r+2a}$，乘自旋因子 $\sqrt{(r+2a)/a}$ 得 (11.38)；$S^2$ 的係數與 $h$ 無關 ⇒ (11.36) 不變；
常數項 $\frac{h^2}{2(r+2a)}$ 只進自由能；核對：$h/r$ ＝ 自旋因子 × $h''/r''$（精確）；
11.7 $Ae^{K'}=2\cosh2K$、$Ae^{-K'}=2$ ⇒ 相乘 $A^2=4\cosh2K$、相除 $e^{2K'}=\cosh2K$ ⇒ (11.54)(11.55)；$K=1$：3.8793、0.6625；
11.8 三種組合 $2\cosh(2K\pm H)$、$2\cosh H$ ⇒ $K_2'=\frac14\ln\frac{\cosh(2K+H)}{\cosh(2K-H)}$、$K_1'=\frac14\ln\frac{\cosh(2K+H)\cosh(2K-H)}{\cosh^2H}$、$A=2[\cosh(2K+H)\cosh(2K-H)\cosh^2H]^{1/4}$；
偶數格點的新外場 $H'=H+2K_2'$；$K=1$、$H=0.1$：3.8897、0.6602、0.0482、0.1964；8 自旋環 $Z_8(K,H)=A^4Z_4(K_1',H')$；
11.9 (a) $x'=\frac{2x}{1+x^2}$，不動點 $x=0$（$T=0$，不穩定）、$x=1$（$T=\infty$，穩定）——$\frac12\ln\cosh2K=K$ 只在 $K=0,\infty$；(b) $x'\approx2x$、$\xi(x')=\xi(x)/2$ ⇒ $\xi\propto x^{-1}=e^{2J/k_BT}$；
精確 $\tanh K'=\tanh^2K$ ⇒ $\xi=-1/\ln\tanh K\approx\frac12e^{2K}$（比 0.4969、0.4999、0.5000）；(c) $A(K)\approx\sqrt2e^K$、$K'\approx K-\frac12\ln2$ ⇒ $\ln Z/N=\sum_j2^{-j}\ln A(K_j)\to K$ ⇒ $F\to-NJ$
（$\ln Z/(NK)=1.1269,1.0091,1.0000$）。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點、表的數字由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；講義裡「÷」可以用 `$\div$`；TikZ 從一點畫回同一點不要用 `to[looseness]`。
