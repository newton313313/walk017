# 《Elements of the Random Walk》第 2 章 → Manim：大綱、檢討與定稿（2026-10-17）

來源：J. Rudnick & G. Gaspari, *Elements of the Random Walk*（CUP 2004）第 2 章 Generating functions I，pp. 25–50。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。原文轉檔 `ref/md/02_generating_functions_I.md`，
逐式校對與筆誤表 `ref/md/ch02_formulas.md`（代理人**不必**讀，本檔 §6 已含全部需要的數學）。

## 1. 草稿大綱（12 幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W201 | 生成函數：一個函數裝下整個數列（e^x、Legendre） | §2.1.1，pp. 25–27，(2.1)，Ex. 2.1–2.3 |
| W202 | 生成函數本身就是答案（巨配分函數、會消失的走路人） | §2.1.1–2.1.2，pp. 27–29，(2.2)(2.3)，Ex. 2.5，Fig. 2.1 |
| W203 | 遞迴：最後一步來自四個鄰居 | §2.1.2，pp. 29–31，(2.4)–(2.9)，Fig. 2.2–2.3 |
| W204 | 平移對稱與傅立葉：一行代數解出 G | §2.1.2，pp. 31–33，(2.10)–(2.13)，Fig. 2.4，Ex. 2.6 |
| W205 | 從 χ^N 取回走法數（布里淵區積分） | §2.1.3，pp. 33–34，(2.14)–(2.17) |
| W206 | 高斯極限與奇偶性的因子 2 | §2.1.3，pp. 34–35，(2.18)–(2.20) |
| W207 | 普適性：三種晶格、同一條高斯 | §2.1.4，pp. 35–36，(2.21)–(2.22) |
| W208 | 離開格子：長尾步長（Cauchy）讓 √N 失效 | §2.1.4，pp. 36–38，(2.23)–(2.31)，Ex. 2.7 |
| W209 | 擴散方程與格林函數 | §2.1.5，pp. 38–40，(2.32)–(2.38)，Ex. 2.8–2.9 |
| W210 | 補充 1：高斯積分與它的生成函數 | §2.2，pp. 41–42，(2:S1-1)–(2:S1-6) |
| W211 | 補充 2：倒晶格與布里淵區 | §2.3，pp. 42–47，(2:S2-1)–(2:S2-34) |
| W212 | 補充 3：從奇點讀出係數的成長＋第 2 章回顧 | §2.4，pp. 47–50，(2:S3-1)–(2:S3-13) |

## 2. 生動度檢討（哪裡太乾 → 怎麼改）

- 全系列主例：**正方格子（Fig. 2.2）上回到原點的走法數** $C(N;\vec 0)=1,4,36,400,4900,63504$（$N=0,2,\dots,10$，$=\binom{N}{N/2}^2$）。
  W203 用遞迴數出 36、W205 用積分算出 63504、W206 用高斯估 66754、W212 用對數奇點再估一次；立方格子的 1172556 當「較大的例子」。
- 「兩枚硬幣」彩蛋接回第 1 章：轉 45° 後 $u=x+y$、$v=x-y$ 各自是一維擲幣 → $C(N;\vec 0)=C_{1D}(N,0)^2$，$63504=252^2$。
- 每幕一個問題開場、結尾一句收束；先猜再揭曉：W202（停在原點的機率 0.145）、W203（36）、W206（高斯只拿到一半）、
  W207（三種晶格其實一樣）、W208（中位距離 236 而不是 15）、W210（$3\sqrt\pi/4$）。
- 「找錯遊戲」紅框（`walk_draw.erratum_box`）：W201 p. 26 的「右邊、對 $y$」、W207 (2.22) 的指數、W208 Ex. 2.7 的歸一化、
  W209 (2.38) 的正負號、W210 (2:S1-4) 的 $b^2/2a$、W212 (2:S3-11) 的 $+\ln z_c$。小筆誤（$\chi(k)$、$d\vec q$、$l_2,l_2$、「左邊」）只用字幕或講義。
- 抽象的幕改成畫面：W204 平面波在格點上「平移＝乘相位」；W205 熱圖 $(\chi/4)^N$ 越來越尖、四個角也亮；W209 機率在格子上擴散（棋盤格）；
  W211 混疊（兩條曲線過同一組格點）與單位根圍成封閉多邊形；W212 複平面上「最近的奇點說了算」。
- 跨幕呼應：奇偶性的因子 2 出現三次（W206 布里淵區角落、W209 棋盤格、W212 兩個等距奇點 $\pm z_c$）；W210 的峰度 3 呼應第 1 章 W105。

## 3. 例子檢討（每幕至少一個小例子算到底＋一個書上的例子或較大的例子）

| 幕 | 從頭算到尾的小例子 | 書上／較大的例子 |
|---|---|---|
| W201 | $(1-t)^{-1/2}$ 展開到 $t^3$ → $P_0..P_3$；$P_n(1)=1$ | Ex. 2.2 $\int P_n=2,0,0,0$；Ex. 2.3 $2/(2n+1)$；$P_5$ |
| W202 | 理想氣體 $x=2$：$Z(N)=1,2,2,\frac43,\frac23,\frac4{15}$，Poisson | 2000 個會消失的人：平均 9.34、停在原點 296 人 vs 0.1452 |
| W203 | $N=2$：中心 4、斜角 2、直走 1 | $N=4$ 整張格子，中心 36，總和 256 |
| W204 | $q=\pi/3$ 的平面波：鄰居和＝自己 | 五種格子的 $\chi$；$\chi(0)=2,4,6,8,12$ |
| W205 | $\chi^2$ 的 16 項，常數項 4 | 正方 4…63504；立方 6…1172556 |
| W206 | 立方 $N=10$：630848 vs 1172556 | 加上角落：1261696（0.929）；正方比值 → 1 |
| W207 | 擲幣 $N=10$：0.246094 / 0.252313 / 0.246005 | sc、bcc、fcc 的 $N^{3/2}P/(g\Omega)$ → 0.3299 |
| W208 | $\chi=e^{-qr_0}$，$q=1$：0.3679 | 中位距離 15.79 vs 236.06（$N=100$） |
| W209 | 二階差分；$P(20;\vec0)=0.031045$ vs 0.031831 | Ex. 2.8 代入；(2.38) 的 4.24213 / 4.23729 / 6.09756 |
| W210 | $x^2-2x=(x-1)^2-1$；$e\sqrt\pi=4.818$ | $n=0,2,4,6$ 的動差表 |
| W211 | 8 個單位根的和 0 或 8 | sc、bcc、fcc 的單胞體積 1、0.7698、0.7071 |
| W212 | $2^{N+1}+1$；$C(20,10)=184756$ vs 187078.97 | 對數奇點 7418.1 / 7317.4 / 4478.3；正方 $2\cdot4^N/(\pi N)$ |

## 4. 程式大綱（已完成，代理人直接用）

- 數學 `codes/utils/walk_gf.py`（純 Python，自我測試 `python utils/walk_gf.py` 通過；不依賴 walk1d）：Legendre 展開、理想氣體、
  會消失的走路人（模擬＋橢圓積分精確值）、格子遞迴計數（正方、sc、bcc、fcc）、兩枚硬幣公式、結構函數、均分求積的傅立葉逆變換、
  高斯近似（含奇偶因子與單胞體積）、(2.22) 的 $q^4$ 修正、Cauchy 步長（$\chi$、$P_N$、取樣、中位距離）、二維擴散、Ex. 2.8、
  (2.38) 的傅立葉核對、高斯積分、倒晶格、單位根、冪級數係數的漸近（極點、分支點、對數）。模組開頭的說明列出每個函式。
- 畫圖 `codes/utils/lattice_draw.py`（已用試畫場景跑過 render_check，全部元件正常）：方格線、格點、數字格、Fig. 2.3 的鄰居星、足跡、
  熱圖與色塊、布里淵區框、單位根與頭接尾的鏈、三維晶格斜投影（sc／bcc／fcc）、分支切割。
- 第 1 章的 `codes/utils/walk_draw.py` 照用：`plot_axes`、`curve`、`bar_chart`、`tex_table`、`erratum_box`、`pole_marker`、`contour_circle`。
- combkit 現成：`ChuScene`（`show_header`、`say`、`clear_stage`、`fade_all`）、`styles.zh/mixed/boxed`、`fulton_draw.assert_in_frame/assert_apart`。

## 5. 省 token 檢討（改了什麼）

- 不做 3D 渲染：Fig. 2.4 與 bcc／fcc 用斜投影的線框（`lattice3d`）；三維的走法只算數字、足跡畫 $(x,y)$ 投影。不另建資料檔：所有數字由 `walk_gf` 現算並 assert。
- 原本想把 §2.1.1 的 Bessel（Ex. 2.4）與偏向擴散（Ex. 2.9）各做一段：只放講義。§2.1.6 Summing up 併進 W212 的回顧。
- 熱圖最多 40×40 色塊、數字格最多 9×9；足跡最多 400 步。
- 代理人**不讀原文**，只讀本檔 §6–§8 與模組開頭的說明。2 個代理人平行（各 6 幕）；每幕最多 2 次 render_check；只看 `--tile 480` 總覽圖。
- 所有代理人回報完才輸出 1080p（一次一支，背景執行）；文件、接線、combkit、git 由我一次做完。講義在代理人做場景時由我同時寫。

## 6. 定稿場景卡（代理人讀這裡）

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$、格距 $a$（畫面上取 $a=1$）、
走法數 $C(N;\vec x,\vec y)$（起點 $\vec x$、終點 $\vec y$；從原點出發時寫 $C(N;\vec r)$）、生成函數 $G(z;\vec x,\vec y)=\sum_N z^NC(N;\vec x,\vec y)$、
結構函數 $\chi(\vec q)$、傅立葉 $g(z;\vec q)$、$c(N;\vec q)=\chi(\vec q)^N$、布里淵區 BZ。
主例：正方格子回到原點 $C(N;\vec 0)=1,4,36,400,4900,63504$（`GF.returns_sq(N)`，$N=0,2,\dots,10$）。
程式代號：`GF` = `utils.walk_gf`、`L` = `utils.lattice_draw`、`D` = `utils.walk_draw`（見 §7）。

**W201 生成函數：一個函數裝下整個數列**（§2.1.1，pp. 25–27，(2.1)，Ex. 2.1–2.3）
1. 問：無窮多個數，能不能裝進一個式子裡？
2. $e^x=1+\frac{x}{1!}+\frac{x^2}{2!}+\cdots$：係數 $1,1,\frac12,\frac16,\frac1{24}$（`GF.exp_coeffs(4)`）——「$e^x$ 生成數列 $1/n!$」。
   一行回顧：第 1 章已見過兩個，$\frac1{2^N}(1+e^y)^N$（W105）與 $\frac1{1-2z}=1+2z+4z^2+\cdots$（W111）。
3. (2.1)：$\dfrac1{\sqrt{z^2-2zx+1}}=\sum_n P_n(x)z^n$（Legendre 多項式）。令 $t=2xz-z^2$：$(1-t)^{-1/2}=1+\frac12t+\frac38t^2+\frac5{16}t^3+\cdots$
   （`GF.half_binom_series(3)`），收集 $z$ 的次方：$P_0=1$、$P_1=x$、$P_2=\frac{3x^2-1}2$（練習 2.1）、$P_3=\frac{5x^3-3x}2$（`GF.legendre_gf(3)`）。
   `D.plot_axes` 畫 $P_0..P_3$，$x\in[-1,1]$。先猜再揭曉：$P_n(1)=?$ → 全是 1（$x=1$ 時左邊 $=\frac1{1-z}=1+z+z^2+\cdots$）。
   筆誤框（`D.erratum_box`，where「p. 26，(2.1) 後」）：書上說「把**右邊**對 $y$ 展開到五次」得 $P_5$ → 應是「把**左邊**對 $z$ 展開」；
   TeX：錯 `\text{right side},\ y^5`、對 `\text{left side},\ z^5`；旁邊 $P_5=\frac{63x^5-70x^3+15x}8$。
4. 練習 2.2：(2.1) 兩邊對 $x$ 從 $-1$ 積到 $1$：左邊 $=\frac{1+z-|1-z|}{z}=2$（每個 $|z|<1$ 都是 2，`GF.ex22_integral`）
   → $\int_{-1}^1P_0\,dx=2$，其餘 $0$：表 $n=0,1,2,3$ → $2,0,0,0$（`GF.poly_int11`）。
5. 練習 2.3：兩個生成函數相乘再積分 $=\frac1{\sqrt{yz}}\ln\frac{1+\sqrt{yz}}{1-\sqrt{yz}}=2+\frac23yz+\frac25(yz)^2+\frac27(yz)^3+\cdots$（`GF.ex23_coeffs(3)`）
   → 只有 $y^nz^n$ 的項 → $\int_{-1}^1P_nP_m\,dx=\frac{2\delta_{nm}}{2n+1}$：$2,\frac23,\frac25,\frac27$（正交）。
收束：生成函數把整串數列打包成一個函數：要哪一項就展開取係數，還能一次證明整串數的性質。

**W202 生成函數本身就是答案**（§2.1.1–2.1.2，pp. 27–29，(2.2)(2.3)，Ex. 2.5，Fig. 2.1）
1. 問：走路的人每走一步有 10% 的機會消失（存活率 $z=0.9$）。平均能走幾步？最後停在原點的機率多大？
2. 統計力學 (2.2)：$Q(\alpha)=\sum_NZ(N)e^{-\alpha N}$；fugacity $z=e^{-\alpha}$（$\alpha=-\mu/kT$）→ $Q(z)=\sum_NZ(N)z^N$：巨配分函數＝正則配分函數的生成函數。
   練習 2.5（理想氣體）：$Z(N)=\frac1{N!}(\mathcal CT^{3/2}V)^N$ → $Q(z)=\exp(z\,\mathcal CT^{3/2}V)$。小例子 $x=\mathcal CT^{3/2}V=2$：
   $Z(N)=1,2,2,\frac43,\frac23,\frac4{15}$（$N=0..5$，`GF.ideal_gas_Z`）；$z=1$：$Q=e^2\approx7.389$；
   粒子數的機率 $Z(N)z^N/Q=0.1353,0.2707,0.2707,0.1804,0.0902,0.0361$（Poisson，`GF.grand_probs(2,1.0,5)`），平均 $2=zx$。
3. 會消失的走路人（p. 28）：剛好走 $n$ 步的機率 $z^n(1-z)$：$0.1,0.09,0.081,\dots$；平均 $\frac z{1-z}=9$。
   模擬 2000 人（`GF.mortal_walkers(2000, 0.9, seed=2)`）：步數 $n=0..9$ 的人數長條 214, 167, 156, 138, 119, 123, 98, 94, 88, 71，
   疊上理論 $2000\cdot0.1\cdot0.9^n$ = 200, 180, 162, 146, 131, 118, 106, 96, 86, 77（四捨五入）；樣本平均 9.34（`GF.MORTAL_MEAN`；理論 9）。
4. 格子上 (2.3)：$C(N;\vec x,\vec y)$＝從 $\vec x$ 走 $N$ 步到 $\vec y$ 的走法數，$G(z;\vec x,\vec y)=\sum_Nz^NC(N;\vec x,\vec y)$。
   主例：$G(z;\vec0,\vec0)=1+4z^2+36z^4+400z^6+4900z^8+63504z^{10}+\cdots$。畫一個會消失的人的足跡（Fig. 2.1 的樣子；
   `GF.mortal_trail(0.9, seed=25)`：24 步，停在 $(3,-3)$；`L.trail` + `L.lattice_lines`）。
5. 揭曉：停在原點的機率 $=\sum_Nz^N(1-z)\frac{C(N;\vec0)}{4^N}=(1-z)\,G(\tfrac z4;\vec0,\vec0)=0.1452$（`GF.mortal_origin_exact(0.9)`）；
   模擬 296/2000 = 0.148（`GF.MORTAL_AT_ORIGIN`）。把 $z$ 換成存活率，生成函數本身就是答案。
收束：$z$ 是「每多一個（粒子、單體、一步）就打的折扣」；生成函數＝所有大小加權後的總和，本身就有物理意義。

**W203 遞迴：最後一步來自四個鄰居**（§2.1.2，pp. 29–31，(2.4)–(2.9)，Fig. 2.2–2.3）
1. 問：正方格子上走 4 步回到原點，有幾種走法？（先猜；答案 36）
2. 回顧一維（第 1 章 W110、W111）：$C(N;x,y)=\frac1{2\pi}\int_{-\pi}^{\pi}e^{-iq(x-y)}\chi(q)^Ndq$（(2.4)）、
   $G(z;x,y)=\frac1{2\pi}\int_{-\pi}^{\pi}\frac{e^{-iq(x-y)}}{1-\chi(q)z}dq$（(2.5)）、$\chi(q)=2\cos q$（(2.6)；書上寫成 $\chi(k)$，字幕帶過）。怎麼推廣到二維、三維？
3. 遞迴 (2.7)（Fig. 2.3，`L.neighbor_star`）：最後一步一定來自四個鄰居 $\vec w_1..\vec w_4$ 之一 →
   $C(N;\vec x,\vec y)=\sum_iC(N-1;\vec x,\vec w_i)$，$C(0;\vec x,\vec y)=\delta_{\vec x,\vec y}$。
4. 數字格（`L.count_grid(GF.count_grid(N), …)`，$-4..4$）：$N=0,1,2,3,4$ 一張張長出來，每格＝上一張四個鄰居的和。
   原點 1, 0, 4, 0, 36（揭曉 36）；總和 1, 4, 16, 64, 256 $=4^N$。$N=2$：中心 4、$(\pm1,\pm1)$ 2、$(\pm2,0),(0,\pm2)$ 1。
   $N=4$ 整張：中心 36；$(\pm1,\pm1)$ 24；$(\pm2,0),(0,\pm2)$ 16；$(\pm2,\pm2)$ 6；$(\pm3,\pm1),(\pm1,\pm3)$ 4；$(\pm4,0),(0,\pm4)$ 1。
   棋盤格：$x+y$ 與 $N$ 不同奇偶的格子都是 0。彩蛋（兩枚硬幣：$u=x+y$、$v=x-y$ 各擲一次，`GF.count_sq`）：$36=6\times6$、$16=4\times4$。
5. 乘上 $z^N$ 對 $N$ 加總（(2.8)），補上 $N=0$ 的 $\delta$：$G(z;\vec x,\vec y)=z\sum_iG(z;\vec x,\vec w_i)+\delta_{\vec x,\vec y}$（(2.9)）。
   只剩 $z$ 一個變數，但右邊要用到鄰居（非局部）→ 下一幕用傅立葉解。
收束：每一格＝四個鄰居的和：這是二維的 Pascal 三角形。

**W204 平移對稱與傅立葉：一行代數解出 G**（§2.1.2，pp. 31–33，(2.10)–(2.13)，Fig. 2.4，Ex. 2.6）
1. 問：(2.9) 每一格都牽扯到鄰居，怎麼一次解出全部？
2. 平面波（一維格點）：$f(x)=\cos(qx)$、$q=\pi/3$：$x=-3..3$ 的值 $-1,-\frac12,\frac12,1,\frac12,-\frac12,-1$；
   左右鄰居的和 $f(x-1)+f(x+1)=2\cos q\cdot f(x)=f(x)$（$\chi(\pi/3)=1$，`GF.plane_wave_check`）：平移只乘上相位，
   「和鄰居的關係」變成「乘上 $\chi(q)$」。
3. 平移對稱：$G$ 只依賴 $\vec x-\vec y$。(2.10) $g(z;\vec q)=\sum_{\vec x}G(z;\vec x-\vec y)e^{i\vec q\cdot(\vec x-\vec y)}$；把 (2.9) 乘上 $e^{i\vec q\cdot(\vec x-\vec y)}$ 對 $\vec x$ 加總
   （(2.11)，相位拆成 $e^{i\vec q\cdot(\vec x-\vec w_i)}e^{i\vec q\cdot(\vec w_i-\vec y)}$）：$g=z\,g\,\chi(\vec q)+1$ → $g(z;\vec q)=\dfrac1{1-z\chi(\vec q)}$（(2.13)），
   $\chi(\vec q)=\sum_ie^{i\vec q\cdot(\vec w_i-\vec y)}$（結構函數）。
4. 各種格子的 $\chi$（Fig. 2.4，`L.lattice3d("sc"/"bcc"/"fcc")`）：一維 $2\cos q$；正方 $2(\cos q_x+\cos q_y)$；簡單立方 (2.12) $2(\cos q_xa+\cos q_ya+\cos q_za)$；
   練習 2.6：bcc $8\cos\frac{q_xa}{\sqrt3}\cos\frac{q_ya}{\sqrt3}\cos\frac{q_za}{\sqrt3}$、fcc $4\big[\cos\frac{q_xa}{\sqrt2}\cos\frac{q_ya}{\sqrt2}+\cos\frac{q_ya}{\sqrt2}\cos\frac{q_za}{\sqrt2}+\cos\frac{q_za}{\sqrt2}\cos\frac{q_xa}{\sqrt2}\big]$
   （$a$＝最近鄰距離）。$\chi(\vec0)$＝鄰居數：2, 4, 6, 8, 12（`GF.chi_1d/chi_sq/chi_sc/chi_bcc/chi_fcc`）。
5. 正方格子 $\chi$ 的熱圖（`L.heatmap(lambda x, y: GF.chi_sq(x, y) / 4, …)` + `L.bz_frame`）：中心 4 最大、四個角 $(\pm\pi,\pm\pi)$ 是 $-4$（伏筆）。
   $\vec q=\vec0$：$g(z;\vec0)=\frac1{1-4z}=1+4z+16z^2+64z^3+256z^4+\cdots$：所有 $N$ 步走法的總數 $4^N$。
收束：傅立葉把「和鄰居的關係」變成「乘上 $\chi(\vec q)$」：$G$ 一行代數就解出來。

**W205 從 χ^N 取回走法數**（§2.1.3，pp. 33–34，(2.14)–(2.17)）
1. 問：把 $(2\cos q_x+2\cos q_y)^N$ 對一個方塊積分，為什麼會數出回到原點的走法？
2. (2.14)：$\frac1{1-z\chi}=1+z\chi+\cdots+(z\chi)^n+\cdots$；$z^N$ 的係數 $c(N;\vec q)=\chi(\vec q)^N$（(2.15)(2.16)）。
3. $\chi=e^{iq_x}+e^{-iq_x}+e^{iq_y}+e^{-iq_y}$（E、W、N、S 四個方向）；$\chi^2$ 展開成 16 項（`GF.chi2_terms()`），每一項是「兩步」；
   常數項 4 個：EW、WE、NS、SN（回到原點）→ $C(2;\vec0)=4$。一般：$e^{i\vec q\cdot\vec r}$ 的係數就是 $C(N;\vec r)$。
4. (2.17)：$C(N;\vec r)=\big(\frac a{2\pi}\big)^d\int_{BZ}e^{-i\vec q\cdot\vec r}\chi(\vec q)^Nd^d q$，BZ：$-\pi/a<q_i\le\pi/a$。
   數值（均分求積，精確，`GF.fourier_count_sq`、`GF.fourier_count_sc`）：正方 $N=2,4,6,8,10$ → 4, 36, 400, 4900, 63504；
   簡單立方 → 6, 90, 1860, 44730, 1172556。熱圖 $(\chi/4)^N$，$N=2,4,10$（`L.heatmap`，vmax = 1）：峰越來越尖，中心與四個角都是 1。
5. 因式分解：$2(\cos q_x+\cos q_y)=4\cos\frac{q_x+q_y}2\cos\frac{q_x-q_y}2$ → $C(N;\vec0)=C_{1D}(N,0)^2$：$63504=252^2$（兩枚硬幣）；
   $C(10;(1,1))=210\times252=52920$（`GF.count_sq(10,1,1)`）。
收束：取係數＝對布里淵區積分；主例 1, 4, 36, 400, 4900, 63504。

**W206 高斯極限與奇偶性的因子 2**（§2.1.3，pp. 34–35，(2.18)–(2.20)）
1. 問：$N=10$，用高斯近似估立方格子回到原點的走法數（精確 1172556），會差多少？
2. (2.18)(2.19)：$(\cos q_xa+\cos q_ya+\cos q_za)^N=\exp[N\ln(3-\frac12a^2q^2+O(q^4))]=3^Ne^{-\frac N6a^2|\vec q|^2+\cdots}$。
   圖（`D.plot_axes`，$a=1$）：沿 $q_x$ 軸 $\big(\frac{\cos q+2}3\big)^{10}$ 與 $e^{-10q^2/6}$，$q\in[-\pi,\pi]$，幾乎重合。
3. (2.20)：$C(N;\vec r)\approx6^N\big(\frac3{2\pi N}\big)^{3/2}e^{-3|\vec r|^2/2Na^2}$；$N=10$、$\vec r=0$：630848（`GF.gauss_sc(10)`）。
   先猜再揭曉：精確 1172556 是它的 1.86 倍——少了一半！
4. 另一半在布里淵區的角落 $\vec q=(\pi,\pi,\pi)/a$：$\chi=-6$，$(-6)^{10}=6^{10}$；沿體對角線畫 $\cos^{10}q$，$q\in[0,\pi]$：兩端都是 1。
   角落的高斯也貢獻一份，相位 $(-1)^{N+x+y+z}$：同奇偶 → 乘 2，不同 → 0（棋盤格）。注記框（`D.erratum_box`，where「p. 34，(2.20) 前」）：
   錯 `\text{only the zone center}`、對 `\text{center}+\text{corner}:\ \times2`（同第 1 章 W110 的 $q=\pm\pi$）。
   更正後 $2\times630848=1261696$，精確／近似 0.929；$N=20$：0.964。
5. 正方格子（主例）：$\frac{2\cdot4^N}{\pi N}$：$N=10$ → 66754（精確 63504；沒乘 2 的 $\frac{4^N}{\pi N}$ = 33377）；
   比值 $N=10,20,50,100$：0.9513, 0.9753, 0.9901, 0.995 → 1（`GF.SQ_RATIOS`）。
收束：$N$ 很大時走法數是高斯：寬度 $\sqrt N$；格子的細節只剩奇偶性的因子 2。

**W207 普適性：三種晶格、同一條高斯**（§2.1.4，pp. 35–36，(2.21)–(2.22)）
1. 問：簡單立方、體心立方、面心立方——走很多步之後，回到原點的機率會不會不一樣？
2. 三種晶格（`L.lattice3d`：鄰居在立方體的面心 6 個、角 8 個、稜中點 12 個；最近鄰距離 $a=1$）：$\chi\approx Z-c_2q^2$：
   $6-q^2$、$8-\frac43q^2$、$12-2q^2$（`GF.CHI_SMALL_Q`）→ $\chi/Z\approx1-\frac{q^2}6$，三個一樣！
3. 回到原點的機率 $P(N;\vec0)=C(N;\vec0)/Z^N$（`GF.UNIV_TABLE`），$N=2,4,6,8,10$：
   sc 0.16667, 0.06944, 0.03987, 0.02663, 0.01939；bcc 0.125, 0.05273, 0.03052, 0.02044, 0.0149；fcc 0.08333, 0.02604, 0.01415, 0.00939, 0.00684——看起來不一樣。
   除以單胞體積 $\Omega$（1、0.7698、0.7071，`GF.CELL_VOLUME`）與奇偶因子 $g$（2、2、1）再乘 $N^{3/2}$：$N=10$ 時 0.307、0.306、0.306（`GF.UNIV_SCALED10`），
   $N=40$ 時三個都是 0.3238 → 極限 $(3/2\pi)^{3/2}=0.3299$。點列圖：三種顏色的點收斂到同一條虛線。
4. 為什麼只看 $q^2$：(2.21) $\int e^{-iqr-\frac{a^2N}6q^2-bNq^4}dq$，把 $q^4$ 項展開到一次 (2.22)：
   $\sqrt{\frac{6\pi}{Na^2}}e^{-3r^2/2a^2N}\big[1-\frac{27b}{Na^4}\big(1-6\frac{r^2}{Na^2}+3(\frac{r^2}{Na^2})^2\big)\big]$，修正 $\propto1/N\to0$。
   一維擲幣的實例（$\cos^Nq=e^{-Nq^2/2-Nq^4/12-\cdots}$，即 $a^2=3$、$b=\frac1{12}$；`GF.coin_quartic`）：
   $N=10,d=0$：精確 0.246094、高斯 0.252313、加 $q^4$ 修正 0.246005；$d=2$：0.205078、0.206577、0.205268；$N=20,d=0$：0.176197、0.178412、0.176182。
   筆誤框（where「p. 35，(2.22)」）：錯 `e^{-3r^2/2aN}`、對 `e^{-3r^2/2a^2N}`（$N=10$、$r=2$、$a=\sqrt3$：照印 0.7072，正確 0.8187）。
5. 普適性的兩個來源：① $\chi(\vec q)=c_0+c_2q^2+\cdots$（最低修正是 $q^2$）；② $N$ 大時只有最低兩項重要。反例：$\chi=c_0+c'|q|^{1.3}+\cdots$ → 不是高斯（下一幕）。
收束：只要 $\chi(\vec q)$ 從 $q^2$ 開始變，走很多步之後一律是同一條高斯——這就是普適性。

**W208 離開格子：長尾步長讓 √N 失效**（§2.1.4，pp. 36–38，(2.23)–(2.31)，Ex. 2.7）
1. 問：如果偶爾會跳很遠很遠，走 $N$ 步後離家多遠？還是 $\sqrt N$ 嗎？
2. (2.23)：$\chi(\vec q)=\int e^{i\vec q\cdot\vec r}p(r)d^3r$；(2.24) 動差展開（書上左邊的 $d\vec q$ 應為 $d^3r$，字幕帶過）；(2.25) 歸一化 → $\chi(0)=1$。
   練習 2.7（高斯步長）：$p(r)=(2\pi r_0^2)^{-3/2}e^{-r^2/2r_0^2}$ → $\chi=e^{-q^2r_0^2/2}$ → $N$ 步仍是高斯，每軸變異數 $Nr_0^2$。
   筆誤框（where「p. 37，Ex. 2.7」）：錯 `\sqrt{2/\pi}\,r_0^{-3}`、對 `(2\pi)^{-3/2}r_0^{-3}`（照印的 $\int p\,d^3r=4\pi\approx12.566$，`GF.radial_integral`）。
3. 長尾的例子 (2.26)：$p(r)=\frac{r_0}{\pi^2(r^2+r_0^2)^2}$（$\int p\,d^3r=1$，但 $\int r^2p\,d^3r=\infty$）。(2.27)–(2.29)：先做角度積分 (2.28)，
   再在上半平面繞 $r=ir_0$ 的二階極點（`D.plot_axes` + `D.pole_marker` 示意）→ (2.30) $\chi(q)=e^{-qr_0}$（$q=1,r_0=1$：0.3679，`GF.chi_cauchy`）。
   關鍵：$e^{-qr_0}\approx1-r_0|q|$：從 $|q|$ 開始變，不是 $q^2$。
4. (2.31)：$\chi^N=e^{-Nqr_0}$ → $P(N;\vec r)=\frac{Nr_0}{\pi^2(r^2+N^2r_0^2)^2}$：形狀和一步一樣，寬度乘 $N$。
   中位距離（`GF.MEDIANS`：400 人、seed 5、$r_0=1$）：高斯步長 $N=1,10,100$ → 1.45, 4.92, 15.79（理論 $1.538\sqrt N$）；
   Cauchy 步長 → 2.5, 23.28, 236.06（理論 $2.264N$，`GF.cauchy_radial_median()`）。先猜再揭曉：$N=100$ 是 236，不是 15。
5. 畫面：兩條 200 步足跡的 $(x,y)$ 投影（`GF.trail_3d(GF.steps_gauss_3d(200, 1.0, random.Random(7)))` 與 `steps_cauchy_3d`，
   各自縮放到 4.5×4.5 的方框；`L.trail`）：左邊一團、右邊幾次大跳躍（Lévy 飛行）。種子不好看可以換，回報時說明。
收束：步長的二階動差發散時，$\chi(q)$ 從 $|q|$ 開始變：普適性失效，距離 $\propto N$ 而不是 $\sqrt N$。

**W209 擴散方程與格林函數**（§2.1.5，pp. 38–40，(2.32)–(2.38)，Ex. 2.8–2.9）
1. 問：一滴墨水滴在格子中央，機率怎麼散開？
2. (2.32)：$P(N;\vec x,\vec y)=\frac1{2d}\sum_{nns}P(N-1;\vec x,\vec y-\vec w_i)$（$2d$＝配位數）；偏向的走法 (2.33) 把 $\frac1{2d}$ 換成 $p_i$。
   動畫：正方格子 `L.heat_cells(GF.prob_grid_sq(N))`，$N=0..20$：棋盤格（奇偶性）一路擴散；變異數 $\sum r^2P=N$（`GF.variance`，$N=20$ → 20）；
   $P(20;\vec0)=0.031045$ vs $\frac2{\pi N}=0.031831$。
3. (2.34)(2.35)：加減 $P(N-1;\vec x,\vec y)$ → $P(N)-P(N-1)=\frac16\sum_j[P(N-1;\vec y-\hat e_ja)+P(N-1;\vec y+\hat e_ja)-2P(N-1;\vec y)]$（二階差分）
   →（(2.36) 左邊 $\approx\partial P/\partial N$）→ (2.37) $\frac{\partial P}{\partial N}=\frac{a^2}6\nabla^2P$（擴散方程；$d$ 維是 $\frac{a^2}{2d}$）。
   書上說解是「(2.20) 的**左邊**去掉 $6^N$」→ 應為右邊（字幕帶過）。
4. 練習 2.8：$P=\big(\frac3{2\pi Na^2}\big)^{3/2}e^{-3r^2/2Na^2}$：$\frac{\partial P}{\partial N}=\big(-\frac3{2N}+\frac{3r^2}{2N^2a^2}\big)P=\frac{a^2}6\nabla^2P$（`GF.ex28_sides`）；$\int P\,d^3x=1$。
5. 格林函數 (2.38)：(2.9) 同樣取連續極限：$(1-2dz)G-za^2\nabla^2G=\delta$。筆誤框（where「p. 40，(2.38)」）：
   錯 `za^2\nabla^2G-(2dz-1)G=\delta`、對 `za^2\nabla^2G+(2dz-1)G=-\delta`（$\nabla^2$ 項的正負號）。
   核對（一維、傅立葉）：精確 $\frac1{1-z\chi(q)}$、正確版 $\frac1{1-2z+zq^2}$、照書上的號 $\frac1{1-2z-zq^2}$；$z=0.4,q=0.3$：4.24213、4.23729、6.09756（`GF.g_fourier/g_cont/g_cont_book`）。
收束：隨機漫步放大來看就是擴散；生成函數 $G$ 是擴散算子的格林函數。

**W210 補充 1：高斯積分與它的生成函數**（§2.2，pp. 41–42，(2:S1-1)–(2:S1-6)）
1. 問：$\int_{-\infty}^\infty x^4e^{-x^2}dx=?$ 能不能一次算出所有 $\int x^ne^{-ax^2}dx$？
2. (2:S1-1)：$\int e^{-ax^2}dx=\sqrt{\pi/a}$（$a=1$：$\sqrt\pi=1.7725$；一句話：平方後換極座標，$I^2=\iint e^{-(x^2+y^2)}=\pi$）。
3. 生成函數 (2:S1-3) $\int e^{-ax^2+bx}dx$：配方 (2:S1-4) $ax^2-bx=a(x-\frac b{2a})^2-\frac{b^2}{4a}$。筆誤框（where「p. 41，(2:S1-4)」）：
   錯 `-\frac{b^2}{2a}`、對 `-\frac{b^2}{4a}`（$a=1,b=2$：$x^2-2x=(x-1)^2-1$；照印會是 $(x-1)^2-2$）。
   動畫（`D.plot_axes` + `D.curve`）：$e^{-x^2+2x}$＝$e^{-(x-1)^2}$ 往右移 1 再乘 $e$ → (2:S1-5) $e^{b^2/4a}\sqrt{\pi/a}$；$a=1,b=2$：$e\sqrt\pi=4.818$（`GF.gauss_gen(1,2)`）。
4. 展開比較係數 (2:S1-6)：$e^{b^2/4a}=\sum_k\frac{b^{2k}}{k!(4a)^k}$ → $\int x^ne^{-ax^2}dx=\sqrt{\frac\pi a}\frac{n!}{(n/2)!(4a)^{n/2}}$（$n$ 偶）、0（$n$ 奇）。
   $a=1$：$n=0,2,4,6$ → $\sqrt\pi\cdot(1,\frac12,\frac34,\frac{15}8)$ = 1.7725, 0.8862, 1.3293, 3.3234（`GF.gauss_moment`，數值積分核對）。揭曉 $\frac{3\sqrt\pi}4\approx1.3293$。
5. 呼應：四階／二階² $=\frac{3/4}{(1/2)^2}=3$（第 1 章 W105 的峰度 3）；(2.20)、(2.22) 用的正是這些積分（W207 的 $27b/Na^4$）。
收束：一個生成函數 $e^{b^2/4a}\sqrt{\pi/a}$，一次給出高斯的所有動差。

**W211 補充 2：倒晶格與布里淵區**（§2.3，pp. 42–47，(2:S2-1)–(2:S2-34)）
1. 問：格子上的函數做傅立葉，為什麼 $\vec k$ 只要在一個方塊（布里淵區）裡取就夠了？
2. 混疊：一維格點 $n=-4..4$ 上，$e^{ikn}$ 與 $e^{i(k+2\pi)n}$ 完全一樣：畫 $\cos(\pi x/4)$ 與 $\cos(9\pi x/4)$ 兩條連續曲線，在 9 個整數點重合。
   所以 $f(\vec k)=\sum_{\vec R}F(\vec R)e^{i\vec k\cdot\vec R}$（(2:S2-21)）是週期函數 $f(\vec k+\vec K)=f(\vec k)$（(2:S2-22)），因為 $\vec K\cdot\vec R=2\pi\times$整數（(2:S2-8)）。
3. 倒晶格 (2:S2-5)–(2:S2-8)：$\vec R=n_1\vec a_1+n_2\vec a_2+n_3\vec a_3$；$\vec b_i=\frac{2\pi}\Omega\vec a_j\times\vec a_k$（(2:S2-7)），$\vec a_i\cdot\vec b_j=2\pi\delta_{ij}$。
   二維示意：正方格子（格距 1）與倒晶格（間距 $2\pi$），布里淵區 $[-\pi,\pi]^2$（`L.bz_frame`）。單胞體積（最近鄰距離 $a$；(2:S2-9)–(2:S2-20)）：
   sc $a^3$、bcc $\frac4{3\sqrt3}a^3=0.7698a^3$、fcc $\frac{a^3}{\sqrt2}=0.7071a^3$（`GF.CELL_VOLUME`）；布里淵區體積 $(2\pi)^3/\Omega$。
   書上小筆誤 p. 43「$l_2,l_2$ and $l_3$」→ $l_1,l_2,l_3$（字幕帶過）。
4. 反演 (2:S2-23)–(2:S2-25)：$\int_{BZ}e^{i\vec k\cdot(\vec R-\vec R')}d^3k=\frac{(2\pi)^3}\Omega\delta_{\vec R,\vec R'}$ → $F(\vec R)=\frac\Omega{(2\pi)^3}\int_{BZ}f(\vec k)e^{-i\vec k\cdot\vec R}d^3k$——就是 (2.17) 的 $(a/2\pi)^3$。
5. 有限格子（週期邊界，(2:S2-27)–(2:S2-34)）：$\vec k=\sum_i\frac{m_i}{N_i}\vec b_i$，布里淵區裡剛好 $N_0=N_1N_2N_3$ 個；$\sum_{\vec R}e^{i\vec k\cdot\vec R}=N_0\delta_{\vec k,0}$。
   $N=8$ 的單位根（`L.roots_of_unity(8)`）：$m=0$ → 8 個箭頭排成一直線（和 8，`L.root_chain(8, 0, …)`）；$m=1$、$m=3$ → 頭接尾圍成封閉多邊形（和 0，`GF.roots_sum`）。
收束：格子上的波，$\vec k$ 和 $\vec k+\vec K$ 分不出來：所以只需要一個布里淵區。

**W212 補充 3：從奇點讀出係數的成長＋第 2 章回顧**（§2.4，pp. 47–50，(2:S3-1)–(2:S3-13)）
1. 問：知道生成函數長什麼樣，不展開，能直接說出第 $N$ 項大概多大嗎？
2. 單極點 (2:S3-2)：$\frac1{z_c-z}=\frac1{z_c}\sum(z/z_c)^n$ → $z^N$ 的係數 $z_c^{-(N+1)}$；$z_c=\frac12$、$N=10$ → 2048。
   兩個極點 (2:S3-3)–(2:S3-5)：$\frac1{1/2-z}+\frac1{1-z}$ → $c_N=2^{N+1}+1$：3, 5, 9, 17, 33, 65（$N=0..5$，`GF.coeff_two_poles`）；
   第二項的比重 $e^{-\delta(N+1)}$、$\delta=\ln2$ → 0：最近的奇點說了算（複平面 `D.plot_axes` + 兩個 `D.pole_marker`，近的那個亮起）。
3. 分支點 (2:S3-6)–(2:S3-10)：$(z_c-z)^{-\alpha}$ 的係數 $\frac{\Gamma(N+\alpha)}{\Gamma(\alpha)\Gamma(N+1)}z_c^{-N-\alpha}\approx\frac{N^{\alpha-1}}{\Gamma(\alpha)}z_c^{-N-\alpha}$。
   例：$(1-4z)^{-1/2}=\sum\binom{2N}Nz^N$：$N=10$：184756 vs $\frac{4^{10}}{\sqrt{10\pi}}=187078.97$（比值 0.9876）。
   一維回到原點 $G(z;0)=(1-4z^2)^{-1/2}$：兩個分支點 $z=\pm\frac12$（`L.branch_cut`）離原點一樣遠 → 兩個都要算：偶數 $N$ 加倍、奇數 $N$ 抵消成 0——
   又是奇偶性的因子 2（呼應 W206）：$\binom{10}5=252$ vs $2^{10}\sqrt{2/(10\pi)}=258.4$（比值 0.9754，`GF.asym_returns_1d`）。
4. 對數奇點 (2:S3-11)：$-(z_c-z)^{-\alpha}\ln(z_c-z)=\frac d{d\alpha}(z_c-z)^{-\alpha}$ → 係數多一個 $\ln N$。筆誤框（where「p. 50，(2:S3-11)」）：
   錯 `+\ln z_c`、對 `-\ln z_c`。例 $\alpha=1$、$z_c=\frac12$、$N=10$（`GF.coeff_log_exact/coeff_log_asym`）：精確 $2^{11}(\ln2+H_{10})=7418.1$；
   正確的三項公式 7317.4；照書上的號 4478.3；只留 $\ln N$ 4715.7。主例：正方格子的 $G(z;\vec0)$ 在 $z=\pm\frac14$ 有對數奇點，
   $\approx-\frac1\pi\ln(1-16z^2)$ → $C(N;\vec0)\approx\frac{2\cdot4^N}{\pi N}$：$N=10$：66754 vs 63504（0.9513）；$N=100$：0.995。
5. 第 2 章回顧（地圖 7 格依序亮起）：生成函數 → 遞迴 (2.9) → 傅立葉 $g=1/(1-z\chi)$ → 取係數 $\chi^N$ → 高斯與普適性 → 擴散與格林函數 → 奇點讀係數。
收束：奇點的位置決定指數成長，奇點的形狀決定冪次。第 2 章完。

## 7. 製作規則（代理人必讀）

- ManimGL 1.7.2（`from manimlib import *`），**不是** Manim Community（用 `ShowCreation` 不用 `Create`）。繼承 `ChuScene`；`construct` 只呼叫**模組層級函式**
  （連續播放會借用 `construct`，寫成方法會 AttributeError）。開頭 `self.show_header("…")`＋`check_math()`，結尾收束字幕後 `self.wait(1.0)`。
- 場景骨架：
  ```python
  import random
  import sys
  from pathlib import Path
  sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
  sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
  from manimlib import *
  from combkit.draw.fulton_draw import assert_apart, assert_in_frame
  from combkit.draw.scene_base import ChuScene
  from combkit.draw.styles import mixed, zh
  from utils import walk_gf as GF
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_gf` 算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標）；版面用 `assert_in_frame`／`assert_apart`。
  不要在場景裡呼叫 `GF.check()`（整套自我測試由 `check_env.py` 跑）。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；字幕裡的「÷」寫成「除以」。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。
- 大量物件：熱圖最多 40×40 色塊、數字格最多 9×9、足跡最多 400 步、長條最多 60 根；避免數千個小物件。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式，回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；每幕最多 2 次。不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W201 | `w201_what_is_gf.py` | `W201WhatIsGF` | A |
| W202 | `w202_gf_is_answer.py` | `W202GFIsAnswer` | A |
| W203 | `w203_lattice_recursion.py` | `W203LatticeRecursion` | A |
| W204 | `w204_structure_function.py` | `W204StructureFunction` | A |
| W205 | `w205_fourier_inversion.py` | `W205FourierInversion` | A |
| W206 | `w206_gaussian_limit.py` | `W206GaussianLimit` | A |
| W207 | `w207_universality.py` | `W207Universality` | B |
| W208 | `w208_off_lattice.py` | `W208OffLattice` | B |
| W209 | `w209_diffusion_green.py` | `W209DiffusionGreen` | B |
| W210 | `w210_gaussian_integrals.py` | `W210GaussianIntegrals` | B |
| W211 | `w211_reciprocal_lattice.py` | `W211ReciprocalLattice` | B |
| W212 | `w212_singularities.py` | `W212Singularities` | B |

連續播放：`main.py` 的 `Ch2Demo`（W201–W212）。

## 9. 講義（我寫，與代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch2/part2}`；內文 `notes/ch2/part2.tex`，圖表與數字檢查 `notes/ch2/make_figs.py`
（`check_numbers_ch2()` 逐一 assert，表格與 TikZ 圖產生成 `notes/ch2/figs.tex`）。小節對應書上 2.1.1–2.1.6 與三個補充，每節：
先看這裡框 → 直觀 → 定義 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒 → 練習與答案（含書上 Ex. 2.1–2.9 的完整解答）。
