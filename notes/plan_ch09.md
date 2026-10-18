# 《Elements of the Random Walk》第 9 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 9 章 Scaling of walks and critical phenomena，pp. 203–232（PDF 頁序 219–248）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_crit.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch08.md` §6–§9。

## 1. 草稿大綱（6 幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W901 | 漫步的指數就是臨界指數：$\langle R^2\rangle\propto N^{2\nu}$ (9.1)(9.2)，Flory 在 1、2、4 維精確、三維差約 2%（書上說一個百分點左右）；de Gennes：$\nu$ 是一個假想磁鐵的指數；磁鐵的自發磁化 $M_s\propto(T_c-T)^\beta$（Fig. 9.1） | §9.1–9.2 開頭，pp. 203–205，(9.1)–(9.3) |
| W902 | O(n) 模型與對稱破缺：$n=1$ Ising、2 xy、3 Heisenberg、$n=0$ 漫步；能量 (9.4) 在整體旋轉下不變 → $\langle\vec M\rangle=0$；可是低溫時翻轉要等「天文數字」那麼久；外場 (9.6)–(9.11) 與相圖（Fig. 9.2） | §9.2.1–9.2.2，pp. 205–211，(9.4)–(9.15) |
| W903 | 臨界指數與標度律：$\chi\propto|t|^{-\gamma}$、關聯函數與關聯長度 $\xi$ (9.14)(9.15)、Table 9.1；Rushbrooke (9.16)、Fisher $\gamma=(2-\eta)\nu$ (9.17)、hyperscaling (9.18)(9.19)；Widom 的齊次形式 (9.20)–(9.23)（Ex. 9.1–9.2）；$\chi=\sum_jG_{ij}$ (9.24)–(9.27) | §9.2.2 後半–9.2.3，pp. 209–215 |
| W904 | Ginzburg–Landau–Wilson：讓自旋長度可變的有效 Hamiltonian (9.28)(9.29)、Fourier (9.30)–(9.32)、只留 $q=0$ 的平均場 (9.33)(9.34)；Landau 自由能 $\frac r2m^2+\frac u4m^4-hm$ 的雙井（Ex. 9.4：$\beta=\frac12$、$\delta=3$） | §9.3，pp. 215–217，(9.28)–(9.35)，Ex. 9.3–9.4 |
| W905 | 從磁鐵讀出漫步：$\xi^2=\sum r^2G/\sum G$ (9.36)；分子是 $\sum_N\langle R^2\rangle\Gamma_Nz^N$、分母是 $\sum_N\Gamma_Nz^N$ (9.37)；取 $z^N$ 的係數 (9.38)–(9.41) → $\langle R^2\rangle\propto N^{2\nu}$ | §9.4，pp. 218–219，(9.36)–(9.41) |
| W906 | $n\to0$ 的魔術：把四次項換成輔助場 (9.46)、高斯積分（補充 9.6）、$Z_n\to1$；Trotter 分解 → 路徑積分 (9.66)–(9.71)（Ex. 9.6）；每次相交乘 $e^{-u/2}$ (9.72)–(9.74) → 關聯函數＝自迴避漫步的生成函數 (9.75)–(9.77)＋第 9 章回顧 | §9.5–9.6，pp. 219–231，(9.42)–(9.77)，(9:S-1)–(9:S-22)，Ex. 9.5–9.6 |

## 2. 生動度檢討

- W901 開場把兩件看似無關的事並排：左邊是一條自迴避漫步越走越大，右邊是一塊磁鐵加熱後磁化消失；字幕「這兩個是同一個數學問題」。
- W902 用一個小磁鐵（例如 $8\times8$ 的二維 Ising，固定種子的 Metropolis）演示：高溫時磁化忽正忽負，低溫時卡在一邊很久不翻——格子越大卡越久（「天文數字」的縮小版）。
- W903 用可以完全算出來的一維 Ising 鏈：關聯 $G(r)=\tanh^rK$、$\xi=-1/\ln\tanh K$、$\chi=\sum_rG(r)$ 兩種算法一致；再把平均場、二維 Ising（精確）、二維自迴避（$n=0$）的指數填進同一張表，逐一驗算標度律。
- W904 的雙井：$r>0$ 一個谷、$r<0$ 兩個谷，小球滾到其中一邊＝對稱破缺；谷底位置 $\propto\sqrt{-r}$ 就是 $\beta=\frac12$。
- W905 回到第 7、8 章的自迴避計數：用 `walk_saw` 的 $c_N$、$\langle R^2\rangle_N$ 組出兩個生成函數，看比值怎麼給出 $\xi^2$。
- W906 是全章最抽象的一幕：用「每個自旋分量的高斯積分給出 $(\cdots)^{n/2}$，$n\to0$ 時整個分母變 1」這一句話當主軸；Trotter 分解用小矩陣的數值誤差 $\propto\epsilon^2$ 示範。
- 找錯遊戲照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W901 | Flory $3/(d+2)$ vs 已知 $\nu$：$d=1,2,3,4$（三維 0.588 只引用） | 平均場 $m=\tanh(m/T)$ 的 $M_s(T)$（$\beta=\frac12$）vs 二維 Ising 的 Onsager 精確式（$\beta=\frac18$） |
| W902 | 2 個或 4 個自旋的 Ising：逐一列出所有組態、算 $Z$ 與 $\langle M\rangle=0$、$\langle|M|\rangle\ne0$ | 小格子 Metropolis 的翻轉等待時間隨 $L$ 變長（固定種子；數字在模組裡算好） |
| W903 | 一維 Ising 鏈的 $G(r)$、$\xi$、$\chi$ | 指數表＋四條標度律的驗算（平均場、二維 Ising、二維 $n=0$）；Ex. 9.1–9.2 |
| W904 | Landau 自由能的極小值：$m^*=\sqrt{-r/u}$ | (9.34) 的平均場積分（有限體積）vs 直接取極小值；Ex. 9.4 的 $\beta,\Delta,\alpha$ |
| W905 | $(z_c-z)^{-s}$ 的 $z^N$ 係數 vs $N^{s-1}z_c^{-N-s}/\Gamma(s)$ | 二維自迴避：兩個生成函數的係數比＝$\langle R^2\rangle_N$ |
| W906 | $2\times2$ 高斯積分 (9:S-11) 與 Wick 配對數 $(2m-1)!!$ | Trotter 誤差 $\propto\epsilon^2$（Ex. 9.6 的交換子項）；在很小的圖上驗證「$n\to0$ 的高溫展開只剩自迴避路徑」（基礎代理人判斷可行的大小） |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_crit.py`（純 Python＋numpy，scipy 可用；可 import 同資料夾的 `walk_saw`、`walk_scale`）：
  平均場與 Onsager 的 $M_s(T)$；小 Ising 系統的精確列舉（$Z$、$\langle M\rangle$、$\langle|M|\rangle$、外場）；固定種子的 Metropolis 翻轉時間（小格子、控制時間）；
  一維 Ising 的 $G(r)$、$\xi$、$\chi$（兩種算法）；指數表與標度律的驗算（Fraction）；Widom 形式（Ex. 9.1–9.2）；Landau 自由能與 (9.34) 的積分（Ex. 9.4）；
  $(z_c-z)^{-s}$ 的係數漸近；用 `walk_saw` 的 $c_N$、$\langle R^2\rangle_N$ 組 (9.36)–(9.41)；高斯積分 (9:S-11)、Wick；Trotter 誤差（Ex. 9.6）；小圖上的 $n\to0$ 驗證（如果可行）。
  `check()` 用兩種以上獨立算法核對全部場景卡數字，約 10 秒內。
- **重用（不改）**：`walk_draw`、`lattice_draw`（格子、自旋箭頭可用格點＋箭頭自己畫）、`walk_saw`、`walk_scale`。
- **畫圖**：不新增模組；自旋格子、雙井、相圖寫在場景檔的模組層級函式。
- combkit：`walk_crit.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex22_ising_scaling.py`；升 1.19.0。

## 5. 省 token 檢討

1. **30 頁 → 6 幕**：(9.47)–(9.65) 的長推導只取主軸（輔助場 → 高斯積分 → $n\to0$ → 指數算子的矩陣元素），細節放講義；補充 9.6 併進 W906。
   場景代理人 A 做 W901–W903、B 做 W904–W906（兩個代理人：6 幕給一個代理人太長）。
2. **重活交給上下文短的代理人**：基礎代理人、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁（Table 9.1 一定要看影像）；`ch09_formulas.md` 只收場景與講義用得到的式子＋筆誤表（我讀正文時看到的可疑處：$\beta\approx1/3$ 的說法、(9.21) 的 $r$ 與 p. 227 的「$r/2=1-T_c/T$」、$z_c=z_0/\mu$、(9.39) 的係數、(9.72)–(9.74) 的 $e^{-u/4}$、$e^{-u}$、$e^{-u/2}$）。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_crit.check()` assert；Metropolis 一律固定種子、在模組裡算好（有快取）。
5. **不新增畫圖模組、不做 3D**。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 9 節目標 9–10 頁**；練習解答 Ex. 9.1–9.6。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log、記憶；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_crit.py` 已完成（`python utils/walk_crit.py` 自我測試通過，約 4 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch09_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$、$\langle R^2\rangle\propto N^{2\nu}$；溫度 $T$、臨界溫度 $T_c$、自發磁化 $M_s$；$K=J/T$（$J=k_B=1$）、無因次外場 $H=h/T$；
O($n$) 的分量數 $n$；臨界指數 $\alpha,\beta,\gamma,\delta,\nu,\eta$；GLW 的約化溫度 $r$、四次項係數 $u$；自迴避條數 $c_N$（書上 §9.4 叫 $\Gamma_N$，畫面第一次出現時一行字幕說明）、每步逸度 $z$、$z_c=1/\mu$。
程式代號：`WC` = `utils.walk_crit`、`SW` = `utils.walk_saw`（只在 W901 用 `SW.sample_saw`、`SW.two_nu_eff`，W905 用 `SW.two_nu_eff`）、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。

**W901 漫步的指數就是臨界指數**（§9.1–9.2 開頭，pp. 203–205，(9.1)–(9.3)，Fig. 9.1）
1. 問：自迴避漫步越走越大、磁鐵加熱就失去磁性——為什麼說這是同一個數學問題？左框：一條 20 步的二維自迴避漫步（`SW.sample_saw(20)`，種子 702）一步一步長出來（`L.trail`）；
   右框：$8\times8$ 的自旋箭頭（上／下，模組層級函式畫）由 `WC.mc_snapshot(8, 1.5)`（64 個全朝上）換成 `WC.mc_snapshot(8, 3.5)`（33 上、31 下）：加熱，磁化消失。
2. (9.1)(9.2)：$\langle R^2\rangle\propto N^p$，Flory $p=\frac6{d+2}$（`WC.flory_p(d)`）。表 $d=1,2,3,4$：Flory $p$＝2、$\frac32$、$\frac65$、1；已知 $2\nu$＝2、$\frac32$、1.176、1（三維 $\nu=0.588$ 只引用，`WC.NU_KNOWN`）。
   一、二、四維精確；三維 0.600 vs 0.588 差 2%（`WC.flory_gap_3d()`＝0.0204；字幕：書上說「一個百分點左右」）。二維精確列舉（第 7 章）的局部指數 $2\nu_{\rm eff}$＝1.4442、1.4505、1.4556（$N=10,12,14$，`SW.two_nu_eff`）→ $\frac32$。
3. (9.3)、Fig. 9.1：$M_s\propto(T_c-T)^\beta$。兩條曲線（`D.plot_axes`，橫軸 $T/T_c\in[0,1.2]$、縱軸 $M_s\in[0,1]$；各取樣 ≤ 60 點）：平均場 $m=\tanh(m/T)$（`WC.mf_m(x)`，$T_c=1$）與二維 Ising 的 Onsager 精確式
   $M_s=(1-\sinh^{-4}(2/T))^{1/8}$（`WC.onsager_M(x*WC.TC_2D)`，$T_c=2.2692$）。表 $T/T_c=0.5,0.9,0.99$：平均場 0.9575、0.5254、0.1725；Onsager 0.9980、0.8959、0.6859。
   log–log 斜率（`WC.beta_slope`）0.5000、0.1250 → $\beta=\frac12$、$\frac18$。字幕：書上說 $\beta\approx\frac13$——那是三維（Ising 約 0.326，只引用）。
4. de Gennes：$p=2\nu$，$\nu$ 是一個假想磁鐵（$n=0$ 的 O($n$) 模型）的指數。對照表（`D.tex_table`＋中文標籤，四列；無數字）：步數 $N$ ↔ $1/(z_c-z)$（離臨界點多遠）；$\langle R^2\rangle\propto N^{2\nu}$ ↔ $\xi^2\propto(z_c-z)^{-2\nu}$；
   條數 $c_N\propto\mu^NN^{\gamma-1}$ ↔ $\chi\propto(z_c-z)^{-\gamma}$；自迴避 ↔ $n=0$ 個分量。
收束：漫步的 $\nu$ 和磁鐵的 $\beta$ 都是臨界指數；第 9 章要說明為什麼「$n=0$ 的磁鐵」會數自迴避漫步。無筆誤框（三維的 2% 用字幕）。

**W902 O($n$) 模型與對稱破缺**（§9.2.1–9.2.2，pp. 205–211，(9.4)–(9.15)，Fig. 9.2）
1. 問：能量在「所有自旋一起轉」之下不變，所以平均磁化一定是 0——那永久磁鐵是怎麼來的？
2. O($n$)：$n=1$ Ising（上下）、2 xy（平面）、3 Heisenberg（空間）、$n=0$ 自迴避漫步（四個小圖示）。(9.4) $H=-J\sum\vec S_i\cdot\vec S_j$、(9.5) $\vec M=\sum\vec S_i$：
   8 個平面箭頭排一圈（`WC.rotation_spins()`，種子 902），全部轉 120°：能量 $E=-6.3851$ 不變、$|\vec M|=6.7943$ 不變，方向 $-13.1^\circ\to106.9^\circ$（`WC.rotation_check()`、`WC.angle_deg`）→ 各個方向一樣好，平均 $\langle\vec M\rangle=0$。
3. 小系統逐一列舉（`WC.ising_exact(WC.RING4, K)`、`WC.ring4_closed(K)`）：4 個自旋圍一圈（$2\times2$ 方塊，16 個組態排成小圖）：$Z=2e^{4K}+12+2e^{-4K}$。
   長條圖 $P(M)$（$M=-4,-2,0,2,4$，`D.bar_chart`）：$K=1$（$T=1$）：0.4504、0.0330、0.0333、0.0330、0.4504（$Z=121.2329$）；$K=0.2$（$T=5$）：0.1283、0.2306、0.2823、0.2306、0.1283（$Z=17.3497$）。
   $\langle M\rangle=0$（左右對稱），但 $\langle|m|\rangle=0.9337$ vs 0.4871：低溫時磁化「不是正就是負」。2 個自旋（$K=1$）：$Z=6.1723$、$\langle|m|\rangle=0.8808$（一行字，`WC.pair_closed(1)`）。
4. Metropolis（`WC.mc_trace(L, T)` → (sweep, m) 各 400 點，一條＝一個 VMobject；`D.plot_axes` 橫軸 sweep 0–4000、縱軸 $m\in[-1,1]$）三條：$L=8$、$T=3.0$（高溫，在 0 附近亂跳）；$L=4$、$T=2.0$（翻 30 次）；
   $L=8$、$T=2.0$（只翻 2 次）（`WC.trace_flips(L, T)`）。表（`WC.flip_tau(L)`，$T=2.0$、50 個複本）：平均翻轉等待 $L=4,6,8$：133、608、2151 sweep——$L$ 每加 2 約乘 4。
   字幕：書上說等待時間是「天文數字」，而且系統越大越久；無限大的系統永遠等不到。
5. 外場與相圖 (9.6)–(9.11)、Fig. 9.2：$\langle\vec M\rangle=-\partial F/\partial\vec h$。小外場（`WC.ising_exact(g, 1.0, H)["m"]`，$K=1$）表：$H=0.1$、0.01：4 個自旋 0.3506、0.0367；$4\times4$ 週期（`WC.TORUS4`）0.9210、0.1584
   ——格子越大，同一個小外場越接近無限大格子的 $M_s=0.9993$（`WC.onsager_M(1.0)`）。相圖（模組層級函式）：$h=0$、$0<T<T_c$ 的粗線是一級相變線（越過去 $M$ 從 $-M_s$ 跳到 $+M_s$），終點是臨界點 $(T_c,0)$。
收束：對稱說平均是 0；可是低溫時翻一次要等很久——格子越大越久，無限大就永遠不翻：這就是對稱破缺。無筆誤框。

**W903 臨界指數與標度律**（§9.2.2 後半–9.2.3，pp. 209–215，(9.13)–(9.27)，Table 9.1，Ex. 9.1–9.2）
1. 問：一個磁鐵有六個臨界指數 $\alpha,\beta,\gamma,\delta,\nu,\eta$——它們彼此獨立嗎？
2. Table 9.1（`D.tex_table`，6 列：量｜行為｜平均場｜三維 $n=0$（書上），`WC.BOOK_3D`）：比熱 $\propto|T-T_c|^{-\alpha}$：0｜0.12；$M\propto|T_c-T|^\beta$：$\frac12$｜0.308；$\chi\propto|T-T_c|^{-\gamma}$：1｜1.19；
   $M(T_c,h)\propto|h|^{1/\delta}$：3｜4.865；$\xi\propto|T-T_c|^{-\nu}$：$\frac12$｜0.602；$C(T_c,\rho)\propto\rho^{-d+2-\eta}$：0｜0.023。
   筆誤框 1（where「p. 210，Table 9.1」）：錯 `\alpha=\alpha'=0.12`、對 `\alpha\approx0.19`；字幕：同一欄的 $\beta,\gamma$ 代 (9.16) 得 $2-2\beta-\gamma=0.194$，$\nu$ 代 (9.18) $2-3\nu$ 也是 0.194（`WC.book3d_alpha()`）。
3. 一維 Ising 鏈（可以完全算出來；$K=1$、$t=\tanh K=0.7616$）：長條圖 $G(r)=t^{|r|}$（$r=-8..8$，`WC.ising1d_G`）：1、0.7616、0.5800、0.4417、0.3364、…；$\xi=-1/\ln t=3.6719$（`WC.ising1d_xi`；虛線 $e^{-|r|/\xi}$）。
   (9.24)–(9.27) $\chi=\sum_rG(r)$：一根一根加（`WC.ising1d_chi_partial(1, R)`）$|r|\le5,10,20$：5.7520、6.9696、7.3615 → $\frac{1+t}{1-t}=e^2=7.3891$（`WC.ising1d_chi`）；16 個自旋的環逐一列舉 7.2021（有限大小，`WC.ising1d_ring_closed(16, 1)`）。
   $\sum r^2G/\sum G$ 除以 $\xi^2$（`WC.ising1d_xi2_moment(K) / WC.ising1d_xi(K)**2`）＝1.9877、1.9998、2.0000（$K=1,2,3$）：是常數 → 那是長度的平方。
   筆誤框 2（where「p. 210，Table 9.1 第二個關聯長度」）：錯 `\xi=\frac{\int C\rho^2d^d\rho}{\int C\,d^d\rho}`、對 `\xi^2=\frac{\int C\rho^2d^d\rho}{\int C\,d^d\rho}`；字幕：比值是長度的平方（第一列寫對了）。
   一行字：$\chi/\xi$＝2.0123、2.0002、2.0000（$K=1,2,3$）→ $\chi\propto\xi^{2-\eta}$，一維 $\eta=1$。
4. 標度律（`WC.EXPONENTS`、`WC.laws`；表三欄：平均場｜二維 Ising（精確）｜二維自迴避 $n=0$）：$\alpha$：0、0、$\frac12$；$\beta$：$\frac12$、$\frac18$、$\frac5{64}$；$\gamma$：1、$\frac74$、$\frac{43}{32}$；$\delta$：3、15、$\frac{91}5$；$\nu$：$\frac12$、1、$\frac34$；$\eta$：0、$\frac14$、$\frac5{24}$。
   驗算（每格打勾）：(9.16) $\alpha+2\beta+\gamma=2$；(9.17) $(2-\eta)\nu=\gamma$；(9.18) $2-d\nu=\alpha$（$d=4,2,2$）；(9.19) $\frac{d+2-\eta}{d-2+\eta}=\delta$。
   字幕：平均場的 hyperscaling 只在 $d=4$ 成立（$d=3$ 得 $\frac12\ne0$）；書上三維欄 $\alpha+2\beta+\gamma=1.926$（框 1 的由來）。
5. Widom 的齊次形式 (9.20)–(9.23)：$G=|r|^{2-\alpha}f\big(h/|r|^\Delta,\ r/|r|\big)$；Ex. 9.1：$\beta=2-\alpha-\Delta$、$\gamma=\Delta-\beta$；Ex. 9.2：$M\propto h^{1/\delta}$、$\delta=\Delta/\beta$ ⇒ $\beta(\delta-1)=\gamma$（`WC.widom_exponents`）。
   二維 Ising $\alpha=0$、$\Delta=\frac{15}8$：$\beta=\frac18$、$\gamma=\frac74$、$\delta=15$；一個具體的 $f$ 量出的斜率（`WC.widom_slopes(*WC.WIDOM_CASE)`）0.1250、1.7500、0.0667（$=\frac1{15}$）；$\frac18\times14=\frac74$。
收束：六個指數只有兩個獨立：知道 $\nu$ 和 $\eta$，其餘都由標度律決定。

**W904 Ginzburg–Landau–Wilson：讓自旋可以伸縮**（§9.3，pp. 215–217，(9.28)–(9.35)，Ex. 9.3–9.4）
1. 問：如果自旋不是固定長度的箭頭，而是可以伸縮的「場」，磁鐵的自由能長什麼樣子？
2. (9.28) $\mathcal H_{GL}=\sum_i\{\frac r2|\vec S_i|^2+\frac12|\nabla\vec S_i|^2+\frac u4|\vec S_i|^4-\vec h\cdot\vec S_i\}$、(9.29) $Z=\prod_i\int d^nS_i\,e^{-\mathcal H}$（沒有 $k_BT$：「有效」Hamiltonian；$r$ 是約化溫度）。
   (9.30)(9.31) Fourier $\vec S(\vec q)=\frac1N\sum_i\vec S_ie^{i\vec q\cdot\vec r_i}$；Ex. 9.3 的離散版（`WC.fourier_check()`，8 點環、種子 904）：$\sum S_i^2=N\sum_q|S(q)|^2$、$\sum(\nabla S)^2=N\sum_q4\sin^2\frac q2|S(q)|^2$（小 $q$ 就是 $q^2$）、
   四次項變成 $q_1+q_2+q_3+q_4=0$ 的和：三個誤差都 $<10^{-12}$。(9.32)(9.33)：只留 $\vec q=0$ → $\mathcal H_{MF}=V[\frac r2m^2+\frac u4m^4-hm]$。字幕：(9.32) 的 $V$、$(2\pi)^d$ 因子書上沒對好，不影響 (9.33)。
3. 雙井（`D.plot_axes`，$m\in[-2,2]$、$F\in[-0.5,1.5]$；曲線 `WC.landau_F(m, r)`，$u=1$）：$r=1,0,-1$ 三條；$r=-1$ 兩個谷 $m^*=\pm1$、$F_{\min}=-0.25$（`WC.landau_closed`）；小球滾到其中一邊＝對稱破缺；
   $h=0.1$ 把雙井歪向右：$m^*=1.0467$（`WC.landau_min(-1, 1, 0.1)`）。谷底 $m^*=\sqrt{-r/u}\propto(-r)^{1/2}$ → $\beta=\frac12$。
4. Ex. 9.4（書上 (9.35) $F=\frac r2m^2+\frac{m^4}{12}-hm$，即 $u=\frac13$；字幕說明）：$\partial F/\partial m=rm+\frac{m^3}3-h=0$；$r=-1$：$m^*=\sqrt3=1.7321$、$F_{\min}=-0.75$（`WC.ex94_closed`）。
   塌縮：$m=|r|^{1/2}\mathcal M(h/|r|^{3/2})$，$x=1$ 時 $\mathcal M=2.1038$（$r=-0.01$ 與 $-0.0001$ 相同，`WC.ex94_scaling`）⇒ $\beta=\frac12$、$\Delta=\frac32$；$F_{\min}\propto|r|^2$ ⇒ $\alpha=0$；$\delta=\Delta/\beta=3$、$\gamma=1$。
   數值斜率（`WC.ex94_slopes()`）：$\beta,\Delta,\alpha,1/\delta,\gamma$＝0.5000、1.5000、0.0000、0.3333、1.0000。字幕：書上 (b) 寫 $r^{2-\alpha}$，應為 $|r|^{2-\alpha}$（不畫框）。
5. (9.34) 有限體積的平均場積分 $Z_{MF}=\int dm\,e^{-VF(m)}$（$n=1$、$r=-1$、$u=1$；`WC.zmf_avg`）表：$h=0.01$：$\langle m\rangle$＝0.0869、0.7576、1.0042（$V=10,100,1000$）→ 谷底 $m^*=1.0050$；
   $h=0$：$\langle m\rangle=0$（任何 $V$）、$\langle|m|\rangle$＝0.8888、0.9922、0.9992；$-\ln Z/V$＝$-0.2719$、$-0.2397$、$-0.2478$ → $F_{\min}=-0.25$。
   字幕：先 $V\to\infty$、再 $h\to0$ 才得到 $m^*\ne0$——和 W902 的翻轉等待是同一件事。
收束：GLW 把自旋變成可伸縮的場；只留 $\vec q=0$ 就是平均場：雙井給出 $\beta=\frac12$、$\delta=3$、$\alpha=0$。無筆誤框。

**W905 從磁鐵讀出漫步：$\langle R^2\rangle\propto N^{2\nu}$**（§9.4，pp. 218–219，(9.36)–(9.41)）
1. 問：只知道磁鐵的 $\xi\propto(z_c-z)^{-\nu}$、$\chi\propto(z_c-z)^{-\gamma}$，能不能推出漫步的 $\langle R^2\rangle\propto N^{2\nu}$？
2. (9.36) $\xi^2=\sum_rr^2G/\sum_rG$；(9.37) 分子＝$\sum_Nc_N\langle R^2\rangle_Nz^N$、分母＝$\chi=\sum_Nc_Nz^N$（$z$ 是每步的逸度，$z_c=1/\mu=0.3791$，`WC.ZC_2D`）。
   表（`WC.saw_table()`）$N=1..6$：$c_N$＝4、12、36、100、284、780；$c_N\langle R^2\rangle_N$＝4、32、164、704、2716、9808。部分和（$N\le14$，`WC.saw_gf(z)[2]`）$\xi^2$＝0.6347、1.2799、2.4593（$z=0.1,0.15,0.2$）：$z$ 越接近 $z_c$，$\xi$ 越大。
   字幕：書上寫 $(z/z_0)^N$，所以說 $z_c=z_0/\mu$；用每步的逸度就是 $z_c=1/\mu$。
3. 係數的漸近：$[z^N](z_c-z)^{-s}\approx\frac{N^{s-1}}{\Gamma(s)}z_c^{-N-s}$（第 2 章補充 3）；精確 ÷ 漸近（`WC.coef_ratio(N, s)`）$N=10,100,1000$：$s=\frac{43}{32}$（$\gamma$）：1.0227、1.0023、1.0002；$s=\frac{91}{32}$（$\gamma+2\nu$）：1.2760、1.0264、1.0026。
   筆誤框 1（where「p. 218，(9.38)→(9.39)」）：錯 `z_c^{-(N+s)}N^{\delta-1}/\Gamma(s)`、對 `z_c^{-(N+s)}N^{s-1}/\Gamma(s)`；字幕：照書上（$\delta=\frac{91}5$）比值是 $1.4\times10^{-17}$（$N=10$），完全不對（`WC.coef_ratio_book(10, WC.S_GAMMA, Fraction(91, 5))`）。
4. (9.38)–(9.41)：$\xi^2\chi=K_1K_2(z_c-z)^{-(\gamma+2\nu)}$ 的係數 ÷ $\chi=K_2(z_c-z)^{-\gamma}$ 的係數 ⇒ $\langle R^2\rangle_N=\frac{K_1}{z_c^{2\nu}}\frac{\Gamma(\gamma)}{\Gamma(\gamma+2\nu)}N^{2\nu}$；
   兩個精確係數之比 ÷ 這個式子（`WC.r2_ratio(N, WC.S_GAMMA, Fraction(3, 4), WC.ZC_2D)`）：1.2476、1.0240、1.0024（$N=10,100,1000$）→ 1。
   筆誤框 2（where「p. 219，(9.41)」）：錯 `\frac{K_1}{z_c^{\gamma}}\frac{\Gamma(\gamma)}{\Gamma(\gamma+2\nu)}N^{2\nu}`、對 `\frac{K_1}{z_c^{2\nu}}\frac{\Gamma(\gamma)}{\Gamma(\gamma+2\nu)}N^{2\nu}`；
   字幕：$z_c^{-(N+\gamma+2\nu)}\div z_c^{-(N+\gamma)}=z_c^{-2\nu}$；二維 $\mu^{3/2}=4.2850$，書上 $\mu^{43/32}=3.6823$（`WC.zc_pow()`）。
   字幕一行：書上說關聯長度 $\propto K_1/(z_c-z)^\nu$；要得到 $K_1K_2$，應是 $\xi^2\approx K_1/(z_c-z)^{2\nu}$（不畫框）。
5. 回到第 7 章的自迴避數據（`SW.two_nu_eff(N)`、`WC.saw_amp(N)`）$N=10,12,14$：$2\nu_{\rm eff}$＝1.4442、1.4505、1.4556 → $\frac32$；$\langle R^2\rangle_N/N^{3/2}$＝0.8299、0.8224、0.8168（慢慢趨於常數）：前因子與 $N$ 無關 ⇒ $\langle R^2\rangle\propto N^{2\nu}$——de Gennes。
收束：磁鐵的 $\xi$ 和 $\chi$ 只要是冪次，取 $z^N$ 的係數一除，就得到漫步的 $\langle R^2\rangle\propto N^{2\nu}$。

**W906 $n\to0$ 的魔術：O($n$) 模型就是自迴避漫步＋第 9 章回顧**（§9.5–9.6，pp. 219–231，(9.42)–(9.77)，(9:S-1)–(9:S-22)，Ex. 9.5–9.6）
1. 問：「$n=0$ 個分量的磁鐵」聽起來荒謬——它為什麼會數自迴避漫步？目標 (9.45)：$\lim_{n\to0}\frac1n\sum_\alpha\langle S_{i\alpha}S_{j\alpha}\rangle=\sum_NC(N;\vec R_i,\vec R_j)z^N$。
2. (9.46) 輔助場：$e^{-\frac u4|\vec S|^4}=\frac1{\sqrt{\pi u}}\int d\phi\,e^{-\phi^2/u-i\phi|\vec S|^2}$（數值 `WC.hs_rhs(s)`，$u=1$：$|\vec S|^2=1,2$ → 0.7788、0.3679 $=e^{-1/4},e^{-1}$）。四次項變二次 → 每個分量一個一樣的高斯積分 → $(\cdots)^n$ → $n\to0$ 時 $Z_n\to1$。
   高斯積分 (9:S-11)（`WC.gauss_I()`）：$\hat T=\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$、$\vec h=(1,0)$：$\frac{2\pi}{\sqrt3}e^{1/3}=5.0627$（數值積分相同）；$\langle x_1x_2\rangle=(\hat T^{-1})_{12}=\frac13$、$\langle x_1^2\rangle=\frac23$（(9:S-18)，`WC.gauss_moment`）；
   Wick (9:S-20)：$\langle x_1^2x_2^2\rangle=\frac49+\frac29=\frac23$、$\langle x_1^4\rangle=3\cdot(\frac23)^2=\frac43$；配對數 $\frac{(2m)!}{2^mm!}$＝1、3、15、105（$m=1..4$，`WC.wick_count`；第 7 章的弧線圖）。
3. Trotter (9.66)、Ex. 9.6：$e^{(\frac12\nabla^2-i\phi)\epsilon}=e^{-i\phi\epsilon}e^{\frac12\nabla^2\epsilon}e^{O(\epsilon^2)}$，下一項 $e^{\frac{i\epsilon^2}4[\phi,\nabla^2]}$（`WC.trotter_err(eps, order)`，6 點環、$\phi$ 種子 906）：
   $\epsilon=0.1,0.05,0.025$：一階誤差 0.005784、0.001502、0.0003827（每次 ÷3.85、÷3.92 → 4：$\propto\epsilon^2$）；乘上交換子項 0.0009465、0.000124、$1.589\times10^{-5}$（÷7.63、÷7.81 → 8：$\propto\epsilon^3$）。
   切成 $M$ 小段、中間插入 $\sum|\vec R\rangle\langle\vec R|$ → 路徑積分 (9.67)–(9.71)（示意：一條折線穿過幾個時間片；無數字）。
4. 相交的權重 (9.72)–(9.74)：一個格點走過 $k$ 次 → $\frac1{\sqrt{\pi u}}\int d\phi\,e^{-\phi^2/u-ik\phi}=e^{-uk^2/4}$（`WC.visit_weight`）：$k=1$：$e^{-u/4}$、$k=2$：$e^{-u}$ ⇒ 每相交一次 $e^{-u/2}=0.6065$（$u=1$，`WC.intersection_factor()`）；
   $\exp(-\frac u4\sum_{l,l'}\delta)$：$k^2=k+2\binom k2$。筆誤框 1（where「p. 226，(9.72)」）：錯 `\int_0^\infty`、對 `\int_{-\infty}^{\infty}`；字幕：從 0 積得 $0.3894-0.2395i$，不是 $e^{-1/4}=0.7788$（`WC.hs_book72()`）。
   $N=8$ 的全部 $4^8$ 條漫步，每相交一次乘 $e^{-u/2}$（`WC.weak_saw_count(8, u)`）：$u=0,1,4,\infty$：65536、23695.11、7988.15、5916（$u\to\infty$ 只剩自迴避 $c_8=5916$）。
   筆誤框 2（where「p. 226，(9.75)」）：錯 `-\int_0^t\frac d2\Big|\frac{d\vec R}{ds}\Big|^2ds\ \frac u4\int\!\!\int\delta`、對 `-\int_0^t\frac d2\Big|\frac{d\vec R}{ds}\Big|^2ds-\frac u4\int\!\!\int\delta`；
   字幕：少了負號就變成每相交一次乘 $e^{+u/2}$：$u=1$ 得 $6.2\times10^5$，越相交越多（`WC.weak_saw_count(8, 1, sign=+1)`）。
5. 小梯子上的 $n\to0$（$2\times3$ 梯子、6 點 7 鍵，格點 $(x,y)\to x+3y$；`WC.ladder_G`、`WC.ladder_saws`；畫格子＋依序高亮自迴避路徑）：權重 $\prod(1+K\vec S_i\cdot\vec S_j)$、$|\vec S|^2=n$：每個閉圈乘 $n$，
   ⇒ $G_{03}=\frac{K+K^3+K^5+nK^5}{1+2nK^4+nK^6}$。表 $K=0.5$：$n=1$（Ising）0.6027、$n=2$（xy）0.5610、$n\to0$：0.65625＝$K+K^3+K^5$（三條自迴避路徑，長 1、3、5）；
   對角 $G_{05}$：0.3562、0.3171、0.40625＝$3K^3+K^5$（長 3、3、3、5）。(9.76)(9.77)：$G=\sum_NC(N)e^{-rN/2}$、$e^{-r/2}=z/z_c$：關聯函數＝自迴避漫步的生成函數。
   字幕：書上 p. 227 的 $r/2=1-T_c/T$ 和 (9.21) 差 2 倍（只差常數倍，不畫框）。
   第 9 章回顧（換頁；表六列，`D.tex_table`＋中文標籤）：W901 $p=2\nu=\frac6{d+2}$、$\beta=\frac12$ vs $\frac18$；W902 $\langle M\rangle=0$ 但翻轉等待 133、608、2151；W903 $\chi=e^{2K}$、四條標度律、$\alpha=0.12\to0.19$；
   W904 $m^*=\sqrt{-r/u}$、$\beta=\frac12$、$\Delta=\frac32$、$\alpha=0$；W905 $\langle R^2\rangle\propto N^{2\nu}$、$z_c^{2\nu}$；W906 每相交 $e^{-u/2}$、$n\to0$ 只剩自迴避。
收束：$n\to0$ 時每個閉圈都帶一個 $n$ 而消失，剩下的就是一條從 $i$ 到 $j$、不碰自己的路——O($n$) 的關聯函數就是自迴避漫步的生成函數。

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
  from utils import walk_crit as WC
  from utils import walk_saw as SW
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_crit`（或卡上註明的 `walk_saw`）算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串，分數用 `Fraction`；
  Fraction 不能直接 `:.4f`，要 `float()`）；版面用 `assert_in_frame`／`assert_apart`。不要在場景裡呼叫 `WC.check()` 或 `WC.summary()`（整套自我測試由 `check_env.py` 跑）。
  隨機一律用模組的預設種子（Metropolis 902、旋轉 902、Fourier 904、Trotter 906；`SW.sample_saw` 702）。最慢的呼叫：`WC.flip_tau(L)`／`WC.mc_trace(L, T)`／`WC.trace_flips` 每個 $(L,T)$ 第一次約 0.3–0.5 秒
  （同一個行程內有快取；W902 全部約 2 秒）、`WC.mc_snapshot` 約 0.03 秒、`WC.weak_saw_count(8, u)` 約 0.01 秒、`WC.ladder_G_xy` 約 5 毫秒、`WC.zmf_avg` 約 1 毫秒；其餘都很快。
  scipy 的 IntegrationWarning 不影響結果（`check_math` 裡用 `warnings.catch_warnings()` 包起來）。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W901 0 個、W902 0 個、W903 2 個、W904 0 個、W905 2 個、W906 2 個。框裡的 Tex 不放中文（英文字用 `\text{}`）。
  卡上標「不畫框」的書上問題只用一行字幕。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。字幕裡的「÷」寫成「除以」。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；矩陣用 `\begin{pmatrix}…\end{pmatrix}`。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「N 的 2ν 次方」或省略）。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——不要拿來當變數（格子邊長用 `Lside`、`size` 之類）；`WC`、`SW` 也是。
- 自旋箭頭（W901、W902）：模組層級函式，一個格點一支短箭頭（`Arrow(…, buff=0)` 或 `Line`＋小三角形），$8\times8$ 最多 64 支、包成一個 VGroup；上＝暖色、下＝冷色。
  Metropolis 軌跡（W902）：一條＝**一個** `VMobject().set_points_as_corners(...)`（400 點），線寬 1–1.5；不要每段一個 `Line`。
  長條圖、雙井、log–log 圖都用 `D.plot_axes`／`D.bar_chart`／`D.curve`，刻度數字自己加。梯子（W906）：`L.lat_point` 定座標、鍵用 `Line`、高亮的路徑用一條 `VMobject` 折線。
  相圖（W902）與路徑積分示意（W906）：模組層級函式，沒有數字。
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
| W901 | `w901_critical_exponent.py` | `W901CriticalExponent` | A |
| W902 | `w902_broken_symmetry.py` | `W902BrokenSymmetry` | A |
| W903 | `w903_scaling_laws.py` | `W903ScalingLaws` | A |
| W904 | `w904_landau_wilson.py` | `W904LandauWilson` | B |
| W905 | `w905_walk_from_magnet.py` | `W905WalkFromMagnet` | B |
| W906 | `w906_n_to_zero.py` | `W906NToZero` | B |

連續播放：`main.py` 的 `Ch9Demo`（W901–W906）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch9/part9}`；內文 `notes/ch9/part9.tex`，圖表與數字檢查 `notes/ch9/make_figs.py`
（`check_numbers_ch9()` 逐一 assert，所有數字都從 `codes/utils/walk_crit.py`（與卡上註明的 `walk_saw`）算；表格與 TikZ／pgfplots 圖產生成 `notes/ch9/figs.tex`；巨集都以 `\Nine` 開頭、`\label` 用 `nine-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–8 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch09_formulas.md`、`walk_crit.py` 的 docstring、
`notes/ch8/part8.tex` 與 `notes/ch8/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 9–10 頁、圖表 6 組：(1) $M_s(T)$：平均場 vs Onsager（pgfplots，各 ≤ 60 點）＋Flory 表（$d=1..4$）＋$\beta$ 的比較（1/2、1/8、三維 0.326 只引用）；
(2) 對稱破缺：4 個自旋的 $P(M)$ 長條（$K=1$、0.2）＋三條 Metropolis 軌跡（`WC.mc_trace`，每條取 200 點即可）＋翻轉等待表（133、608、2151）＋外場表（0.3506、0.0367 vs 0.9210、0.1584）；
(3) 指數：Table 9.1（照抄並標出兩處更正）＋三欄指數與四條標度律的驗算表＋一維 Ising 的 $G(r)$ 長條與 $\chi$ 部分和＋Widom 斜率（0.1250、1.7500、0.0667）；
(4) GLW：雙井（$r=1,0,-1$）＋Ex. 9.4（$\mathcal M(1)=2.1038$、五個斜率）＋(9.34) 的表（$V=10,100,1000$）；
(5) §9.4：$c_N$、$c_N\langle R^2\rangle_N$ 表＋$\xi^2(z)$ 表＋係數比表（兩個 $s$）＋$r_2$ 比值表＋$2\nu_{\rm eff}$ 與 $\langle R^2\rangle_N/N^{3/2}$；
(6) §9.5–9.6：(9.46) 的數值、高斯積分與 Wick、Trotter 表（一階／二階／總誤差）、弱自迴避計數（65536、23695.11、7988.15、5916）、梯子圖（TikZ）＋$G_{03}$、$G_{05}$ 表。
小節對應書上 9.1–9.2、9.2.1–9.2.2、9.2.3、9.3、9.4、9.5–9.6，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch09_formulas.md` 的筆誤表，10 列；
另附〔注〕：三維 Flory 差 2%（不是 1%）、$\beta\approx\frac13$ 是三維、Table 9.1 的 $\xi^2\propto|T-T_c|^{-\nu}$ 要讀成 $\xi$、三維欄是舊估計（文獻 $\nu\approx0.5876$、$\gamma\approx1.157$，只引用）、
(9.27) 被積函數少 $|\vec x|^{-(d-2+\eta)}$、(9.32) 的正規化、(9.35) 是 $m^4/12$、Ex. 9.4(b) 的 $|r|$、「Chapter 1」應為第 2 章補充 3、$z_0$ 與 $z_c=z_0/\mu$ 的記號、(9.44) 的 $\partial^2\ln Z$、
(9.65) 沒寫出 $e^{-\sum\phi_i^2/u}$、$C(t;\cdots)$ 是弱自迴避（$u\to\infty$ 才嚴格）、p. 227 的 $r/2=1-T_c/T$ 與「$=J/T$」、(9:S-16) 的 $\langle x_rx_k\rangle$）。
最後練習 9.1–9.6 的解答：
9.1 $M=-\partial G/\partial h=-|r|^{2-\alpha-\Delta}f_1(0,\pm1)$ ⇒ $\beta=2-\alpha-\Delta$；$\chi=\partial M/\partial h\propto|r|^{2-\alpha-2\Delta}$ ⇒ $\gamma=2\Delta-2+\alpha=\Delta-\beta$；兩式合起來就是 Rushbrooke (9.16)；
二維 Ising（$\alpha=0$、$\Delta=\frac{15}8$）：$\beta=\frac18$、$\gamma=\frac74$，數值斜率 0.1250、1.7500；
9.2 $r=0$ 時 $G$ 必須與 $r$ 無關 ⇒ $f(x)\sim x^{(2-\alpha)/\Delta}$ ⇒ $G\propto h^{(2-\alpha)/\Delta}$、$M\propto h^{(2-\alpha-\Delta)/\Delta}=h^{\beta/\Delta}$ ⇒ $\delta=\Delta/\beta$ ⇒ $\beta(\delta-1)=\Delta-\beta=\gamma$（Widom）；平均場 3、二維 Ising 15、二維自迴避 $\frac{91}5$；
9.3 代入 (9.31)、用 $\sum_ie^{i(\vec q+\vec q\,')\cdot\vec r_i}=(2\pi)^d\delta(\vec q+\vec q\,')$（格距 1）：二次項 $\frac{V^2}{2(2\pi)^d}\int(r+|\vec q|^2)|\vec S(\vec q)|^2d^dq$（格子上 $|\vec q|^2\to\sum4\sin^2\frac{q_k}2$）、四次項帶 $\delta(\vec q_1+\cdots+\vec q_4)$；
書上 (9.32) 的前因子只差常數；離散版（8 點環）三個恆等式誤差 $<10^{-12}$；
9.4 $rm+\frac{m^3}3-h=0$；$m=|r|^{1/2}\mathcal M(h/|r|^{3/2},\pm1)$：$\beta=\frac12$、$\Delta=\frac32$；$F_{\min}=|r|^2\mathcal K$ ⇒ $\alpha=0$；例 $\mathcal M(1)=2.1038$、$m^*=\sqrt{-3r}$、$F_{\min}=-\frac34r^2$；另 $\delta=3$、$\gamma=\gamma'=1$（振幅比 2）；
9.5 $\overleftrightarrow V\,\overleftrightarrow V^{-1}=\int_0^\infty\overleftrightarrow Ve^{-t\overleftrightarrow V}dt=\int_0^\infty\big(-\frac{d}{dt}e^{-t\overleftrightarrow V}\big)dt=\overleftrightarrow I-\lim_{t\to\infty}e^{-t\overleftrightarrow V}=\overleftrightarrow I$（$r>0$ 時 $\frac r2-\frac12\nabla^2$ 正定、$i\phi$ 只轉相位，所以 $e^{-t\overleftrightarrow V}\to0$）；
$\overleftrightarrow V$ 與 $e^{-t\overleftrightarrow V}$ 對易 ⇒ 左右逆相同；數值（6 點環，`WC.ex95_check()`）誤差 $<10^{-10}$；
9.6 $e^{\epsilon(X+Y)}=e^{\epsilon X}e^{\epsilon Y}e^{-\frac{\epsilon^2}2[X,Y]+O(\epsilon^3)}$（BCH／Zassenhaus），$X=-i\phi$、$Y=\frac12\nabla^2$：下一項 $e^{\frac{i\epsilon^2}4[\phi,\nabla^2]}$；
數值：一階誤差比 3.85、3.92 → 4、加上這一項後 7.63、7.81 → 8；$t=1$ 切 $M=10,20,40$ 段的總誤差 0.02366、0.0117、0.005821（$\propto1/M$，路徑積分在 $M\to\infty$ 時才精確）。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點、長條與表的數字由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；講義裡「÷」可以用 `$\div$`；TikZ 從一點畫回同一點不要用 `to[looseness]`。
