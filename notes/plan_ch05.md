# 《Elements of the Random Walk》第 5 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 5 章 Variations on the random walk，pp. 95–126（PDF 頁序 111–142）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_var.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch04.md` §6–§9。

## 1. 草稿大綱（7 幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W501 | 有偏向的漫步：會移動的高斯，寬度 $\propto\sqrt{p(1-p)}$；三維的「先走一步再亂走」 | §5.1–5.1.1，pp. 95–98，(5.1)–(5.10)，Fig. 5.1，Ex. 5.1 |
| W502 | 有慣性的漫步：記得上一步，左右兩個狀態與 $2\times2$ 轉移矩陣 | §5.2，pp. 98–100，(5.11)–(5.22) |
| W503 | 持續長度：大 $N$ 時又是高斯，步長換成 $l_p$；另一種推法（二階遞迴） | §5.2–5.2.1，pp. 100–103，(5.23)–(5.38)，Ex. 5.2 |
| W504 | 電報方程：等速 $v$、以速率 $\gamma$ 反轉；兩個往外跑的 δ 尖峰＋拖在後面的分布 | §5.2.2，pp. 103–109，(5.39)–(5.64)，Fig. 5.2–5.3，Ex. 5.3–5.6 |
| W505 | 二維、三維的「跑一跑、轉個彎」（細菌的 run-and-tumble）：因果性、前緣、算子 $\Gamma$ | §5.2.3–5.2.5，pp. 109–118，(5.65)–(5.106)，Fig. 5.4–5.6，Ex. 5.7 |
| W506 | 連續時間的漫步：等一下再走——等待時間、走了幾步的機率、Laplace 轉換把卷積變乘法 | §5.3，pp. 118–122，(5.107)–(5.123)，Ex. 5.8 |
| W507 | 等待＋步伐：$A(\lambda)/(1-B(\lambda)\chi(q))$ 與高斯極限；會消失的走路人（Ex. 5.9）＋第 5 章回顧 | §5.3，pp. 122–125，(5.124)–(5.133)，Ex. 5.9 |

## 2. 生動度檢討

- 主例延續第 1 章的擲硬幣：W501 換成「不公正的硬幣」$p=0.6$，同樣的十步醉漢一邊走一邊被風往右吹；先猜「平均位置」與「寬度」再揭曉。
- W502–W503：「習慣往同一邊走」的醉漢（$Q$ 大）與「老是回頭」的醉漢（$Q$ 小）並排走；兩個極端 $Q=1$（直線前進）、$Q=0$（原地來回）一眼看出寬度為何會到無限大與 0。
- W504 重現書上 Fig. 5.3（$v=\gamma=1$，$t=1,3,5,8$）：兩個尖峰往外跑、越來越矮（$e^{-\gamma t}$），後面的分布慢慢變成高斯；先問「有沒有人跑得比光（$v$）快？」——因果性。
- W505 用真實的例子：大腸桿菌「游一段、翻滾、換方向」；二維、三維的分布用曲線（書上 Fig. 5.6 的 $r^2\Sigma$）與模擬點雲對照，不做 3D 動畫。
- W506–W507：等公車的比喻（等待時間的分布）；固定速率 → 走了幾步是 Poisson 分布；最後回到高斯，寬度 $\propto t$。
- 筆誤紅框照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W501 | $p=0.6$、$N=10$：平均 $2$、變異數 $9.6$，精確二項分布 | Fig. 5.1 的移動高斯；「直接寫擴散方程」會把 $\sqrt{p(1-p)}$ 換成 $1/2$（Ex. 5.1） |
| W502 | $Q$ 給定、$N\le4$：逐步數 $L$、$R$ 的機率 vs 矩陣乘冪 | 特徵值 (5.17)(5.18)、生成函數 (5.20)–(5.22) |
| W503 | 精確變異數 vs $N\,Q/(1-Q)\,l^2$（大 $N$） | 持續長度 $l_p$；$Q=1$、$Q=0$ 兩個極端 |
| W504 | $t=1,3,5,8$ 的 $\Sigma(x,t)$ 與尖峰重量 $e^{-\gamma t}$ | 大 $t$ 的高斯，$D=v^2/2\gamma$（Ex. 5.5–5.6） |
| W505 | 二維 (5.89) 的數值：總數守恆（積分＋尖峰＝1） | 三維 $r^2\Sigma(r,t)$（Fig. 5.6）、前緣的對數發散 |
| W506 | 固定速率 $p_0$：走了 $m$ 步的機率 $(p_0t)^me^{-p_0t}/m!$ | 一般等待分布的 Laplace 轉換 (5.119)–(5.123) |
| W507 | Ex. 5.8 的閉式；Ex. 5.9 的衰減率 $p_1$ | (5.131)–(5.133) 的高斯，寬度 $\propto t$ |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_var.py`（純 Python，可 import 同資料夾的 `walk_gf`、`walk_rec`、`walk_bc`）：偏向漫步（精確二項 vs 高斯，(5.7) 的寬度，Ex. 5.1）、
  有慣性的漫步（逐步 DP vs 轉移矩陣特徵值，(5.17)(5.18)，精確變異數 vs $NQ/(1-Q)$，持續長度）、電報方程 (5.61)（Bessel $I_0,I_1$；與 $v,\gamma$ 的粒子模擬或格子極限交叉驗證；尖峰重量、總數守恆、大 $t$ 高斯）、
  二維 (5.89) 與三維 (5.98)(5.101)(5.103) 的數值（總數守恆、Fig. 5.6 的曲線；可用固定種子的模擬核對）、連續時間漫步（Poisson、Laplace 轉換的數值反演或閉式、(5.128) 的高斯極限、Ex. 5.8–5.9）。`check()` 用兩種以上獨立算法核對全部場景卡數字。
- **重用（不改）**：`walk_draw`（數線、長條、`curve`、`plot_axes`、表格、筆誤框）、`lattice_draw`（足跡、熱圖）。
- **畫圖**：不新增模組；缺的小元件寫在場景檔的模組層級函式。
- combkit：`walk_var.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex18_persistent_walk.py`；升 1.15.0。

## 5. 省 token 檢討

1. **32 頁 → 7 幕**：§5.1.1 併進 W501；§5.2.1（二階遞迴）併進 W503；§5.2.3–5.2.5（二維、三維與算子 $\Gamma$）合成一幕 W505，只展示結果與數值曲線，(5.85)–(5.106) 的圍道細節只放講義；回顧併進 W507。
   場景代理人 A 做 W501–W504、B 做 W505–W507（B 的三幕較重）。
2. **我不在這個長對話裡做重活**：照第 4 章——基礎代理人（校對、模組、定稿卡、combkit）、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁（110 dpi）；`ch05_formulas.md` 只收場景與講義用得到的式子＋筆誤表。三維 (5.95)–(5.106) 的長推導不逐式校對，只核對最後的數值曲線與守恆。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_var.check()` assert；場景與講義代理人都不讀原文。
5. **不新增畫圖模組、不做 3D 動畫**。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 5 節目標 8–10 頁**（第 4 節寫到 11 頁）；練習解答 Ex. 5.1–5.9。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log；walkit 是否重新複製由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_var.py` 已完成（`python utils/walk_var.py` 自我測試通過，約 1 秒）：下面每個數字都在它的 `check()` 裡 assert 過（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。
書上的式子與筆誤見 `ref/md/ch05_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$、向右機率 $p$、「和上一步同方向」的機率 $Q$（第一步左右各半）、$L$／$R$＝最後一步向左／向右的機率、
電報方程的速率 $v$、反轉率 $\gamma$（$v=\gamma=1$）、二三維的轉向率 $\gamma$（轉成完全隨機的新方向）、連續時間漫步的等待：「固定速率」$q(t)=p_0=1$（Poisson）與「兩段等待」（等待時間密度 $te^{-t}$，$q(t)=\frac t{1+t}$，平均等 2）。
密度一律標準歸一化：二維徑向密度 $P_2=2\pi r\rho$、三維 $P_3=4\pi r^2\rho$（$\int P\,dr$＋尖峰＝1）；書上 Fig. 5.6 的縱軸 $r^2\Sigma=\frac\pi2P_3$。
程式代號：`V` = `utils.walk_var`、`BC` = `utils.walk_bc`、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7；`L`、`D` 是模組，左右狀態寫 `Lw`、`Rw`，擴散係數寫 `Dc`）。
**大綱 §1 的 W507 寫的 $A(\lambda)/(1-B(\lambda)\chi(q))$ 是筆誤**：照書上 (5.128) 用 $B(\lambda)/(1-\chi(q)A(\lambda))$。§2–§3 的「寬度 $\propto t$」照抄了書上的筆誤：應為變異數 $\propto t$（W507 放筆誤框）。

**W501 有偏向的漫步**（§5.1–5.1.1，pp. 95–98，(5.1)–(5.10)，Fig. 5.1，Ex. 5.1）
1. 問：不公正的硬幣（向右的機率 $p=0.6$）。醉漢走 10 步，平均到哪裡？分布比公正硬幣寬還是窄？（先猜。）
2. 精確二項 (5.1)＝逐步遞迴 (5.2)（`V.biased_pmf(10, 0.6)`＝`V.biased_dp`）：長條圖 $x=-10..10$（`D.bar_chart`）；平均 $(2p-1)N=2$、變異數 $4Np(1-p)=9.6$（公正硬幣 10）。
   三根標數字：$P(10;0)=0.2007$、$P(10;2)=0.2508$、$P(10;4)=0.2150$；疊上 $2\times$(5.7) 的高斯（`D.curve`，`2*V.biased_gauss(x, 10, 0.6)`）：0.2091、0.2575、0.2091（一行字幕：只有偶數格有人，所以 ×2）。
3. Fig. 5.1：$N=10,20,40$ 三組往右移、變寬：平均 2、4、8；變異數 9.6、19.2、38.4；峰值 0.2508、0.1797、0.1279（$2\times$高斯 0.2575、0.1821、0.1288）。
   (5.7)：寬度 $\propto\sqrt{p(1-p)}$，$p\to1$ 時寬度 $\to0$（沒有隨機性了）。
4. 泰勒展開 (5.3)＋筆誤框（where「p. 96，(5.3)」）：錯 `l(2p-1)\,\partial_y P`、對 `l(1-2p)\,\partial_y P`；字幕：往右偏、峰往右移，$\partial P/\partial y$ 前面要是負號（(5.4) 已經寫對）。
   Ex. 5.1 的陷阱：直接寫成擴散方程 → $D=l^2/2=0.5$、變異數 10、寬度 3.162；多保留一項 → $D=2p(1-p)=0.48$、變異數 9.6、寬度 3.098（`V.diffusion_D`）。
   寬度比 $\sqrt{4p(1-p)}$：$p=0.6$：0.9798、$p=0.9$：0.6（`V.width_ratio`）→ 擴散方程只在 $p\approx\frac12$ 時對。
5. (5.9)(5.10) 三維的「先走一步再亂走」：$\chi(\vec q)=e^{-i\vec q\cdot\vec r_B}\chi_0(\vec q)$ → $P(N;\vec r-N\vec r_B)$：整團平移、寬度不變。一維對照表（`D.tex_table`）：
   「被風吹 $b=0.2$ 再擲公正硬幣」（`V.shift_pmf(10, 0.2)`）平均 2、變異數 10 vs「不公正硬幣 $p=0.6$」平均 2、變異數 9.6。例：在流動的水裡擴散的染料。
收束：偏向讓高斯往前走、寬度乘上 $\sqrt{4p(1-p)}$；「先走一步再亂走」則只平移、不變窄。

**W502 有慣性的漫步：左右兩個狀態**（§5.2，pp. 98–100，(5.11)–(5.22)）
1. 問：醉漢有 $\frac34$ 的機率照上一步的方向走（$Q=\frac34$）。走 4 步之後，變異數比普通醉漢的 4 大還是小？（先猜。）
2. 兩個狀態 (5.11)(5.12)：$L$＝最後一步向左、$R$＝最後一步向右。逐步表（`V.persistent_dp(N, Fraction(3, 4))`，每格寫 $(L,R)$；＝暴力列舉 $2^N$ 條路線）：
   $N=0$：$0{:}(\frac12,\frac12)$；$N=1$：$-1{:}(\frac12,0)$、$1{:}(0,\frac12)$；$N=2$：$-2{:}(\frac38,0)$、$0{:}(\frac18,\frac18)$、$2{:}(0,\frac38)$；
   $N=3$：$-3{:}(\frac9{32},0)$、$-1{:}(\frac18,\frac3{32})$、$1{:}(\frac3{32},\frac18)$、$3{:}(0,\frac9{32})$；
   $N=4$：$-4{:}(\frac{27}{128},0)$、$-2{:}(\frac{15}{128},\frac9{128})$、$0{:}(\frac{13}{128},\frac{13}{128})$、$2{:}(\frac9{128},\frac{15}{128})$、$4{:}(0,\frac{27}{128})$。
3. 揭曉（表，`V.persistent_var`）：變異數 $N=1..4$：$Q=\frac34$：1、3、$\frac{11}2$、$\frac{33}4$；$Q=\frac12$（普通）：1、2、3、4；$Q=\frac14$（老是回頭）：1、1、$\frac32$、$\frac74$。
4. 轉移矩陣 (5.15)(5.16)：`T = [[Q e^{-iq}, (1-Q) e^{-iq}], [(1-Q) e^{iq}, Q e^{iq}]]`；特徵值 (5.17)(5.18)（`V.eigs`）：$q=0.1$：$\lambda_1=0.9848$、$\lambda_2=0.5077$；
   $q=0.5$：$0.6582\pm0.2584i$（$Q|\sin q|>1-Q$ 時成共軛複數）。傅立葉核對（`V.total_fourier`＝`V.T_total`＝`V.total_eig`）：$N=5$、$q=0.3$：0.5699；$N=20$：0.0388。
5. 生成函數 (5.19)–(5.22)：$(1-\mathbf Tz)^{-1}$＝伴隨矩陣 $/D$，$D=1-2zQ\cos q+z^2(2Q-1)=(1-z\lambda_1)(1-z\lambda_2)$；(5.24)(5.25)：$[z^N]\frac1D=\frac{\lambda_1^{N+1}-\lambda_2^{N+1}}{\lambda_1-\lambda_2}$；
   從 $L=R=\frac12$ 出發，總數 $=c_N-(2Q-1)\cos q\,c_{N-1}$（$c_k=[z^k]\frac1D$）。一行：大 $N$ 由最大的 $\lambda_1$ 主導。無筆誤框。
收束：記得上一步的醉漢，用兩個狀態就又是「乘矩陣」；走久了由最大的特徵值說了算。

**W503 持續長度**（§5.2–5.2.1，pp. 100–103，(5.23)–(5.38)，Ex. 5.2）
1. 問：$Q=\frac34$ 的醉漢走 1000 步，寬度相當於普通醉漢走幾步？（先猜。）
2. (5.23)：$\lambda_1\approx1-\frac{Q}{1-Q}\frac{q^2}2$（$q=0.1$：0.985 vs 精確 0.9848）→ (5.26)(5.27) 高斯，變異數 $NQ/(1-Q)$。表（`V.persistent_var`、`V.persistent_var_asym`）：
   $Q=\frac34$：$N=10,100,1000$：26.00、296.00、2996.00 vs 30、300、3000；$Q=\frac14$：3.777、33.778、333.778 vs 3.333、33.333、333.333（差趨近常數 4 與 $-\frac49$）。
   $N=100$、$x=0$：精確 0.04587 vs $2\times$(5.27) 0.04607（奇偶性的因子 2）。
3. 持續長度 (5.28)–(5.31)：$n=\frac{Q}{1-Q}$：$Q=\frac34$ → 3、$Q=\frac12$ → 1、$Q=\frac14$ → $\frac13$（`V.persistence_length`），$l_p=nl$。
   揭曉：1000 步 ≈ 普通醉漢 3000 步（變異數 3000，精確 2996）＝ 每步長 3、走 333 步。
4. 兩個極端（各畫一條路線 `D.walk_path`）：$Q=1$：變異數 $N^2$（$N=10$：100），直線前進、永遠不是高斯；$Q=0$：$N=10$：0、$N=11$：1（原地來回）。
5. 二階遞迴 (5.34)–(5.38)：$R(N+2)=2Q\cos q\,R(N+1)+(1-2Q)R(N)$（`V.recursion2_residual` = 0）、矩陣 (5.36)、特徵方程 (5.37)。
   筆誤框（where「p. 103，(5.37)」）：錯 `\lambda^2-2Q\cos ql\,\lambda+(1-2Q)=0`、對 `\lambda^2-2Q\cos ql\,\lambda-(1-2Q)=0`；
   字幕：把 $\lambda_1=0.9848$ 代入書上的式子差 1.000（`V.eig_residual(l1, 0.1, 0.75, book=True)`）。Ex. 5.2：$R=A\lambda^N$ 的兩個解相加，回到 W502 的總數。
收束：慣性只是把步長放大成持續長度 $l_p=\frac{Q}{1-Q}l$；走得夠久，又是高斯。

**W504 電報方程：兩個尖峰＋拖在後面的分布**（§5.2.2，pp. 103–109，(5.39)–(5.64)，Fig. 5.2–5.3，Ex. 5.3–5.6）
1. 問：走路人以速率 $v=1$ 前進，平均每單位時間反轉 $\gamma=1$ 次。$t=3$ 時有沒有人在 $x=3$ 之外？一次都沒轉身的有多少？（先猜。）
2. (5.39)–(5.44)：$R$、$L$ 隨波前進＋以速率 $\gamma$ 互換；$\Sigma=R+L$、$\Delta=R-L$（(5.45)–(5.48)）；Ex. 5.4：$\Sigma_{tt}+2\gamma\Sigma_t=v^2\Sigma_{xx}$（電報方程）。
   Ex. 5.3 格子極限（`V.tele_lattice(1.0, dt)`：$Q=1-\gamma\delta t$、$l=v\delta t$ 的有慣性漫步，密度＝機率 $/2l$）：$t=1$、$x=0$：$\delta t=0.01$：0.3392；$\delta t=0.005$：0.3380；外插 0.3368＝公式。
3. 解 (5.61)（更正後）：$\Sigma=\frac{\gamma e^{-\gamma t}}{2v}\big[I_0(\gamma s)+\frac tsI_1(\gamma s)\big]\Theta(vt-|x|)+\frac{e^{-\gamma t}}2[\delta(x-vt)+\delta(x+vt)]$，$s=\sqrt{t^2-x^2/v^2}$（`V.tele_cont`、`V.tele_spike`）。
   筆誤框（where「p. 107，(5.60)(5.61)」）：錯 `\frac{1}{2}\Big(2\gamma+\frac{\partial}{\partial t}\Big)`、對 `\frac{1}{2v}\Big(2\gamma+\frac{\partial}{\partial t}\Big)`；
   字幕：$v=2$、$t=3$ 照書上積分得 1.9502，不是 1（(5.58) 有 $1/v$，`V.tele_mass(3, v=2.0, book=True)`）。
4. Fig. 5.3（`D.plot_axes`，$x\in[0,9]$，$t=1,3,5,8$ 四條曲線，尖峰用垂直箭頭）：$\Sigma(0,t)=0.3368$、0.2199、0.1738、0.1388；前緣後方 0.2759、0.06223、0.01179、0.0008387；
   每個尖峰 $\frac12e^{-t}=0.1839$、0.02489、0.003369、0.0001677（因果性：$|x|>vt$ 處是 0）。守恆：連續部分 $1-e^{-t}=0.6321$、0.9502、0.9933、0.9997，加上尖峰＝1（`V.tele_mass`）。
   揭曉：$t=3$ 沒轉身的 $e^{-3}=0.0498$；10 萬人模擬（`V.rt_mc(1, 3.0, 100000, 5)`，`V.mc_fraction`）：$|x|<1.5$ 的比例 0.6060（公式 0.6051，`V.radial_fraction(1, 0, 1.5, 3)`），沒轉身 0.0488。
5. 大 $t$（Ex. 5.5–5.6）：高斯 $D=v^2/2\gamma=0.5$（`V.tele_gauss`）：$\Sigma(0,t)$ vs $1/\sqrt{4\pi Dt}$：0.3368 vs 0.3989、0.2199 vs 0.2303、0.1738 vs 0.1784、0.1388 vs 0.1410；
   精確變異數 $\frac{v^2}\gamma\big[t-\frac{1-e^{-2\gamma t}}{2\gamma}\big]$（`V.tele_var`）：0.5677、2.5012、4.5000、7.5000（→ $2Dt-0.5$）；持續長度 $l_p=v/\gamma=1$（$D=\frac12vl_p$）。
收束：電報方程是有慣性的擴散：前緣以 $v$ 跑、尖峰以 $e^{-\gamma t}$ 變小，後面的人慢慢變成 $D=v^2/2\gamma$ 的高斯。

**W505 二維、三維的「跑一跑、轉個彎」**（§5.2.3–5.2.5，pp. 109–118，(5.65)–(5.106)，Fig. 5.4–5.6，Ex. 5.7）
1. 問：大腸桿菌以速率 $v$ 游、平均每單位時間翻滾 $\gamma$ 次換一個隨機方向（$v=\gamma=1$）。$t=3$ 時有多少還在第一條直線上？分布在前緣長什麼樣子？（先猜。）
2. (5.65)(5.69)：$\Gamma=\gamma-\gamma|0\rangle\langle0|$（翻滾後方向完全隨機）；(5.73)–(5.75) 的反算子（Ex. 5.7：角度離散成 64 格核對，誤差 $<10^{-12}$，`V.sm_residual()`）；
   $\langle0|A^{-1}|0\rangle$：二維 $1/\sqrt{(\lambda+\gamma)^2+q^2v^2}$、三維 $\arctan\frac{vq}{\lambda+\gamma}/vq$（$\lambda=q=1$：0.4472、0.4636，`V.angle_avg_2d_closed`＝數值角度平均）。圍道細節不放（講義才放）。
3. 二維（更正後、極點＋切割合起來）：$P_2=2\pi r\rho=\frac{\gamma r}{v}\frac{e^{-\gamma t+\gamma s/v}}{s}$，$s=\sqrt{v^2t^2-r^2}$，外加 $r=vt$ 的圓環 $e^{-\gamma t}$（`V.rt2_radial`）。
   $t=3$：$r=1,2,2.5,2.9$：0.2978、0.4167、0.3941、0.4052（前緣 $1/\sqrt{\ }$ 發散、可積）；圓環 0.0498、其餘 0.9502。
   筆誤框 1（where「p. 112，(5.84)」）：錯 `2\gamma\sqrt{\pi/2}\,e^{-\gamma t}\cosh(\gamma\sqrt{v^2t^2-r^2})`、對 `2\gamma\,e^{-\gamma t}\cosh(\gamma\sqrt{v^2t^2-r^2}/v)`；
   筆誤框 2（where「p. 113，(5.88)(5.89)」）：錯 `\frac{2\pi}{x}\delta(x-y)`、對 `\frac{1}{x}\delta(x-y)`；字幕：照書上的係數積分，$t=3$ 得 1.2658，不是 1（`V.rt2_mass(3, book=True)`）。
4. 三維 (5.98)+(5.101)+(5.103)：$P_3=4\pi r^2\rho$（`V.rt3_radial`，每條曲線 ≤ 120 點）。重畫 Fig. 5.6（`D.plot_axes`，$r\in[0,8]$，$t=1,3,5,8$）；表 $r=\frac t4,\frac t2,\frac{3t}4,0.9t$：
   $t=1$：0.1472、0.5334、1.0275、1.3699；$t=3$：0.1729、0.4558、0.4828、0.3594；$t=5$：0.1869、0.3630、0.2169、0.0926；$t=8$：0.1948、0.2411、0.0622、0.0118。
   $t=3$ 三類人：沒翻滾 $e^{-3}=0.0498$（球殼）、翻一次 $3e^{-3}=0.1494$（(5.103) 最後一項，前緣 $\ln$ 發散）、其餘 0.8009。
   因果性：$r=3.5>vt$ 時 (5.101) $=0.0113$、(5.103) $=-0.0113$，相加為 0（`V.rt3_pole_book`、`V.rt3_cut_book`）。
5. 模擬對照（`V.rt_mc(d, 3.0, 100000, 5)`，直方圖＋公式曲線）：$r<1.5$ 的比例：二維 0.3335（公式 0.3310）、三維 0.2880（公式 0.2881，`V.radial_fraction`）；圓環／球殼 0.0493、0.0497（公式 0.0498）。
   前緣一行表：一維有限（$\Sigma(vt^-)=0.06223$）、二維 $1/\sqrt{\ }$ 發散、三維 $\ln$ 發散——都可積，總數都守恆。
收束：維度越高，前緣越尖；但沒有人跑得比 $v$ 快，也一個都不少。

**W506 連續時間的漫步：等一下再走**（§5.3，pp. 118–122，(5.107)–(5.123)，Ex. 5.8）
1. 問：走路人每落腳一次就要「等公車」才走下一步。3 分鐘內走了幾步？固定速率和「越等越可能來」哪個走得多？（先猜。）
2. (5.107)(5.108)：$a(\tau)=\exp(-\int_0^\tau q\,dt)$（`V.survival`＝`V.survival_numeric`；畫 $q(t)$ 與 $a(t)$ 兩條曲線）：固定速率 $p_0=1$：$a(3)=e^{-3}=0.0498$；
   兩段等待（$q(t)=\frac t{1+t}$）：$a(3)=4e^{-3}=0.1991$。(5.109)：$\int aq\,dt=1-a(\infty)=1$（一定會走）。
3. (5.110)(5.111)：走 $m$ 步＝卷積；Laplace 轉換 (5.112)–(5.121) 把卷積變乘法：$A^mB$，$A=1-\lambda B$（(5.119)）。
   $\lambda=0.7$：固定速率 $A=B=0.5882$；兩段等待 $A=0.3460$、$B=0.9343$（`V.AB`＝`V.AB_numeric`）。
4. $t=3$ 走了 $m$ 步的機率（兩組長條，`V.steps_prob`＝`V.steps_prob_laplace`＝`V.steps_prob_conv`）：固定速率 $m=0..6$：0.0498、0.1494、0.2240、0.2240、0.1680、0.1008、0.0504；
   兩段等待 $m=0..4$：0.1991、0.4481、0.2689、0.0720、0.0108。平均 3 與 1.2506（`V.mean_steps`）。揭曉：固定速率走得多。(5.122)：$\sum A^mB=\frac B{1-A}=\frac1\lambda$ → 機率加起來是 1。
5. Ex. 5.8 加上步伐（一維格子 $\chi=\cos q$）：$\frac1{\lambda+p_0(1-\cos q)}$ → $P(x,t)=e^{-p_0t}I_{|x|}(p_0t)$（`V.ctrw_bessel`＝`V.ctrw_lattice`）：$t=3$：$x=0,\pm1,\pm2$：0.2430、0.1968、0.1118；
   字幕：步數不固定，所以每一格都可能有人（沒有奇偶性）。無筆誤框。
收束：等待時間的卷積，Laplace 轉換後只是乘法：走 $m$ 步的轉換就是 $A^mB$。

**W507 等待＋步伐＝高斯；會消失的走路人；第 5 章回顧**（§5.3，pp. 122–125，(5.124)–(5.133)，Ex. 5.9）
1. 問：把等待和步伐一起算，走很久以後分布是什麼？時間變 4 倍，寬度變幾倍？（先猜。）
2. (5.124)–(5.128)：$n$ 步位移的轉換 $\chi(q)^n$ → 總和 $\frac{B(\lambda)}{1-\chi(q)A(\lambda)}$（`V.ctrw_transform`＝Ex. 5.8 的 `V.ex58`）。
3. 高斯極限 (5.129)–(5.133)：$\chi\approx1-Cq^2$、$A\approx1-\lambda B(0)$ → $\frac1{\lambda+(C/B(0))q^2}$ → 變異數 $2Ct/B(0)$；一維格子 $C=\frac12$，$B(0)$＝平均等待。
   $P(0,10)$：固定速率（$B(0)=1$）0.1278 vs 高斯 0.1262；兩段等待（$B(0)=2$）0.1805 vs 0.1784（`V.ctrw_lattice`、`V.ctrw_gauss`；長條＋曲線）。
   筆誤框（where「p. 124，(5.133) 後」）：錯 `\text{width}\propto t`、對 `\text{variance}\propto t,\ \ \text{width}\propto\sqrt{t}`；
   數字：固定速率變異數 $t=10,40$：10.00、40.00 → 寬度 3.162 → 6.325（時間 ×4、寬度 ×2）；兩段等待 4.75、19.75（$\approx t/2$，`V.ctrw_var`）。
4. Ex. 5.9（$p_0=1$、$p_1=0.1$）：$A=\frac{p_0}{p_0+p_1+\lambda}$、$B=\frac1{p_0+p_1+\lambda}$ → $\frac1{\lambda+p_1+p_0(1-\chi)}$ → 總數 $e^{-p_1t}$，衰減率 $p_1$；
   $t=10$：0.3679（`V.mortal_total`＝`V.mortal_series`＝`V.mortal_laplace`），半衰期 6.931。
5. 第 5 章回顧（表五列，`D.tex_table`＋中文標籤）：偏向 → 移動的高斯、寬度 $\times\sqrt{4p(1-p)}$（W501）；慣性 → $l_p=\frac{Q}{1-Q}l$（W502–W503）；電報 → 尖峰 $e^{-\gamma t}$＋$D=\frac{v^2}{2\gamma}$（W504）；
   二三維 → 前緣、因果性、守恆（W505）；等待 → $\frac B{1-\chi A}$、變異數 $2Ct/B(0)$（W506–W507）。
收束：改了規則，走路人最後還是回到高斯——只是平均、寬度與時間尺度換了樣子。

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
  from utils import walk_var as V
  from utils import walk_bc as BC
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_var` 算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串）；版面用 `assert_in_frame`／`assert_apart`。
  不要在場景裡呼叫 `V.check()`（整套自我測試由 `check_env.py` 跑）。模擬一律用卡上的參數（`V.rt_mc(d, 3.0, 100000, 5)`，種子 5、10 萬個），數字才對得上。
  最慢的呼叫：`V.rt3_radial` 每點約 0.1 毫秒（一條 120 點的曲線約 0.01 秒）、`V.rt_mc(3, 3.0, 100000, 5)` 約 0.05 秒（$t=8$ 約 0.1 秒）、`V.rt3_mass` 約 0.03 秒；其餘都很快。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W501 1 個、W502 無、W503 1 個、W504 1 個、W505 2 個、W506 無、W507 1 個。
  框裡的 Tex 不放中文（英文字用 `\text{}`）。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；字幕裡的「÷」寫成「除以」。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「e 的負 γt 次方」或省略）。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——左右狀態的機率寫 `Lw`、`Rw`，擴散係數寫 `Dc`；分數用 `Fraction`，Tex 裡寫成 `\frac{a}{b}`。
- 曲線：`D.plot_axes` 上用 `axes.get_graph` 或 `D.curve`；有尖峰（δ）的地方畫垂直箭頭並在旁邊寫重量；前緣發散的曲線（W505）在 $y$ 軸上限處截斷（`min(值, 上限)`），並用一行字幕說明「發散但可積」。
- 大量物件：長條最多 60 根、表格最多 6×6、小路線圖最多 16 張、點最多 60 個；直方圖最多 30 根；避免數千個小物件。不做 3D 動畫（三維只畫徑向曲線與示意圖）。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式（例如矩陣框、尖峰箭頭、直方圖），回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W501 | `w501_biased_walk.py` | `W501BiasedWalk` | A |
| W502 | `w502_persistent_walk.py` | `W502PersistentWalk` | A |
| W503 | `w503_persistence_length.py` | `W503PersistenceLength` | A |
| W504 | `w504_telegrapher.py` | `W504Telegrapher` | A |
| W505 | `w505_run_and_tumble.py` | `W505RunAndTumble` | B |
| W506 | `w506_waiting_time.py` | `W506WaitingTime` | B |
| W507 | `w507_ctrw_gauss.py` | `W507CtrwGauss` | B |

連續播放：`main.py` 的 `Ch5Demo`（W501–W507）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch5/part5}`；內文 `notes/ch5/part5.tex`，圖表與數字檢查 `notes/ch5/make_figs.py`
（`check_numbers_ch5()` 逐一 assert，所有數字都從 `codes/utils/walk_var.py` 算；表格與 TikZ 圖產生成 `notes/ch5/figs.tex`；巨集都以 `\Five` 開頭、`\label` 用 `five-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–4 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch05_formulas.md`、`walk_var.py` 的 docstring、
`notes/ch4/part4.tex` 與 `notes/ch4/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 8–10 頁、圖表 4–5 張：(1) 偏向漫步：$p=0.6$、$N=10$ 的精確 vs $2\times$高斯（0.2007/0.2091、0.2508/0.2575、0.2150/0.2091）＋Fig. 5.1 的平均、變異數、峰值表；Ex. 5.1 的寬度表（$p=0.6$：3.098 vs 3.162、比 0.9798；$p=0.9$：比 0.6）；
(2) 有慣性的漫步：變異數表（$Q=\frac34,\frac14$；$N=10,100,1000$；精確 vs $NQ/(1-Q)$）、$Q=\frac34$ 的 $N\le4$ 小表（$(L,R)$）、持續長度 3、1、$\frac13$；
(3) 電報方程：Fig. 5.3 的四條曲線（TikZ／pgfplots，由 `V.tele_cont` 取樣，尖峰畫箭頭標重量 0.1839、0.02489、0.003369、0.0001677）＋表（$\Sigma(0,t)$、前緣、$1-e^{-t}$、變異數、高斯 $1/\sqrt{4\pi Dt}$）；
(4) 二三維：$P_3$ 的 Fig. 5.6 重畫（四條，前緣截斷）＋表（§6 W505 的 4×4 數字、二維 $P_2$ 在 $t=3$ 的四個值、模擬 vs 公式 0.3335/0.3310、0.2880/0.2881、三類人 0.0498/0.1494/0.8009）；
(5) 連續時間漫步：$t=3$ 的 $P_m$ 兩組長條或表＋$P(0,10)$ 與變異數表（0.1278/0.1262、0.1805/0.1784；10、40；4.75、19.75）＋Ex. 5.9 的 0.3679。
小節對應書上 5.1–5.3，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch05_formulas.md` 的筆誤表，10 列）。
二三維一節可以放 (5.65)–(5.106) 的推導大要（極點 (5.99)、切割兩邊的 $\pm\frac\pi2$、「只轉一次彎」的 $\gamma te^{-\gamma t}$ 與 $\int_0^1u\ln\frac{1+u}{1-u}du=1$），並註明書上的歸一化（二維 $2\pi\rho$、三維 $2\pi^2\rho$，所以「$r^2\Sigma$ 積分 $\pi/2$」）。
最後練習 5.1–5.9 的解答：5.1 $D=\frac{l^2}2-\frac{v^2}2=2p(1-p)l^2$（只有 $p\to\frac12$ 時等於 $\frac{l^2}2$）；5.2 $R=A_1\lambda_1^N+A_2\lambda_2^N$ → 總數 $c_N-(2Q-1)\cos q\,c_{N-1}$；
5.3 $R(x,t+\delta t)=(1-\gamma\delta t)R(x-v\delta t,t)+\gamma\delta t\,L(x-v\delta t,t)$ 展開到 $\delta t$ 一階；5.4 $\Sigma_{tt}+2\gamma\Sigma_t=v^2\Sigma_{xx}$；5.5 $I_0(z)\approx e^z/\sqrt{2\pi z}$ → $D=\frac{v^2}{2\gamma}$ 的高斯；
5.6 變異數 $\frac{v^2}\gamma[t-\frac{1-e^{-2\gamma t}}{2\gamma}]\approx\frac{v^2t}\gamma$，持續長度 $v/\gamma$；5.7 直接相乘（Sherman–Morrison）；5.8 $\frac1{\lambda+p_0(1-\chi(q))}$，一維格子 $e^{-p_0t}I_{|x|}(p_0t)$；
5.9 $\frac1{\lambda+p_1+p_0(1-\chi)}$，總數 $e^{-p_1t}$、衰減率 $p_1$。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；TikZ 從一點畫回同一點不用 `to[looseness]`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算 Bessel 函數）。
