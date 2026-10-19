# 《Elements of the Random Walk》第 12 章 → Manim：大綱與省 token 檢討（2026-10-19）

來源：Rudnick & Gaspari 第 12 章 More on the renormalization group，pp. 285–322（PDF 頁序 301–338）。全書最後一章。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_eps.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch11.md` §6–§9。

## 1. 草稿大綱（6 幕；本章 38 頁）

| 幕 | 主題 | 書上 |
|---|---|---|
| W1201 | 加上四次項之後：Fourier 空間的 GLW Hamiltonian (12.1)–(12.3)、倒晶格（Fig. 12.1–12.2）、球形布里淵區 $\Lambda=e^{-\ell}$ 與殼層（Fig. 12.3–12.4）；「怎樣正確地積掉自由度」(12.5)–(12.10)：先對 $S'$ 求和得 $K'$，再重新標度回到原形——用一個算得完的小例子示範 | §12.1–12.2.4 前半，pp. 285–289 |
| W1202 | 積掉殼層：把四次項分成殼內／殼外 (12.11)、對殼層的 Gaussian 平均＝配對 (12.12)–(12.15)、只剩連通的圖（連結叢集，Fig. 12.6–12.11、補充 12.7）；組合因子 $\frac{5!}{3!\,2!}=10$（Ex. 12.9）、$\frac{n!}{n_1!\,n_2!\,(2!)^{n_2}}$ (12:S-4)、$f_3=6$（Ex. 12.11）；Gaussian 變數的 $\langle e^X\rangle=e^{\text{連通部分}}$ 數值核對（Ex. 12.1） | §12.2.4 後半、§12.7，pp. 289–292、317–322 |
| W1203 | 一迴圈修正與遞迴：兩張圖的值 (12.16)–(12.19)、殼層體積 (12.20)–(12.25)、動量守恆 (12.26)(12.27)、重新標度 (12.28)–(12.36) → $\frac{dr}{d\ell}=2r+\frac{12u}{1+r}$、$\frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}$ (12.37)(12.38)；$(r,u)$ 平面的流：Gaussian 不動點與 Wilson–Fisher 不動點 | §12.2.4–12.2.5，pp. 292–297 |
| W1204 | $\epsilon$ 展開求 $\nu$：$u(\ell)$ 的解與不動點 (12.39)–(12.41)（Ex. 12.2–12.3，Fig. 12.13）、$r=t-\cdots$ 的換變數與線性化 (12.42)–(12.51)、$\xi(\ell^*)=1$ → $\nu=\frac12+\frac{\epsilon}{12}$ (12.52)–(12.59)、$\epsilon=1$：0.583 對三維 Ising 的已知值 (12.60)（Ex. 12.4） | §12.2.5–12.2.6，pp. 296–300 |
| W1205 | O(n) 與自迴避漫步：$n$ 分量的四次項與圖 (12.61)–(12.65)（Ex. 12.5–12.7）、12 → $4(n+2)$、36 → $4(n+8)$ (12.66)(12.67)、$\nu=\frac12+\frac{(n+2)\epsilon}{4(n+8)}$ (12.68)；$n=0$：$\nu=\frac12+\frac{\epsilon}{16}$、$\epsilon=1$ 得 0.5625 (12.69)(12.70)（Ex. 12.8）；Flory $\frac{3}{d+2}$ 展開成 $\frac12+\frac{\epsilon}{12}$ (12.71)(12.72)；ε 展開的四個隱憂 | §12.3–12.4，pp. 300–304 |
| W1206 | 圖解法的重整化群（概念）：質量算子與重整化傳播子 (12.73)–(12.80)、四點函數 $\Gamma_4$ 與 $g_0(q_0)$ 的微分方程 (12.81)–(12.86)＝(12.39)、冪次數算 $(ur^{-(4-d)/2})^n$ 被 $u_R\propto q^{4-d}$ 救回 (12.87)–(12.89)、反項與 Callan–Symanzik 方程；**第 12 章與全書回顧** | §12.5–12.6，pp. 304–317 |

## 2. 生動度檢討

- W1201 用一張「積掉自由度」的流程圖貫穿最後兩章：分殼層 → 平均 → 重新標度 → 參數變了。小例子沿用第 11 章 W1104 的手法（先對一部分變數求和，看耦合常數怎麼變），但這次多一個四次項。
- W1202 的「配對遊戲」：四次頂點像有四隻手，殼內的手兩兩牽起來；畫出 type 1、2、3 的圖，數組合因子；指數化只留下連通圖——和第 6 章的 Wick 配對、第 7 章的 Dyson 方程互相呼應。
- W1203、W1204 是全書數學的高潮：$(r,u)$ 平面的流向圖（數值積分畫流線），兩個不動點一目了然；$\epsilon$ 越小 Wilson–Fisher 點越靠近 Gaussian 點，$\epsilon=0$ 時重合。$\nu$ 用兩種方法算：線性化矩陣的特徵值（數值）與一階公式 $\frac12+\frac{\epsilon}{12}$，畫出兩者隨 $\epsilon$ 的差。
- W1205 的總表：$n=0,1,2,3$ 在 $\epsilon=1$ 的一階 $\nu$ vs 已知值（只引用）vs Flory；自迴避漫步從第 1 章的 $\sqrt N$ 一路走到這裡的 $N^{0.5625}$（已知 0.588）。
- W1206 最後一段是全書回顧：12 章各一行（主角＋關鍵數字），配上全書 87 幕的代號（用文字表即可，不載入影片）。
- 找錯遊戲照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W1201 | 兩個變數的玩具：$e^{KSS'}$ 對 $S'$ 積掉（Gaussian 或 ±1）得 $K'$ | 二維正方晶格的倒晶格、球形布里淵區殼層的體積比 $1-e^{-d\,\delta\ell}\approx d\,\delta\ell$ |
| W1202 | 5 個四次片段配成 3 個 type 1＋1 個 type 2：$\frac{5!}{3!\,2!}=10$（列舉驗證） | $\frac{n!}{n_1!n_2!(2!)^{n_2}}$ 的一般式（列舉 vs 公式）、$f_3=6$；連結叢集的數值核對 |
| W1203 | 一個殼層積分的值 $\frac{12u}{1+r}\delta\ell$（Ising） | $(r,u)$ 流線（數值積分）、兩個不動點的位置（$\epsilon=1$、0.5、0.1） |
| W1204 | $u(\ell)$ 的解析解 vs 數值積分；不穩定的不動點 $u=0$（Ex. 12.3） | 線性化特徵值 → $\nu$（數值）vs $\frac12+\frac{\epsilon}{12}$；$\epsilon=1$ 的 0.583 |
| W1205 | $n=0$：$\nu=\frac12+\frac{\epsilon}{16}$，$\epsilon=1$ → 0.5625 | $n=0..3$ 的表；Flory 的展開 $\frac{3}{6-\epsilon}=\frac12+\frac{\epsilon}{12}+\cdots$；已知值只引用 |
| W1206 | 冪次數算：書上 $u^2r^{2d-3}=r^{-1}(ur^{-(4-d)/2})^2$ 的指數核對（應為 $u^2r^{d-3}=r\,(ur^{-(4-d)/2})^2$，基礎代理人更正） | $g_0(q_0)$ 方程與 (12.39) 同形（換變數 $q_0=e^{-\ell}$）；全書回顧表 |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_eps.py`（純 Python＋numpy，scipy 可用；可 import 同資料夾的 `walk_rg`、`walk_crit`、`walk_saw`，照 `walk_rg` 的相對 import 寫法）：
  兩變數玩具的耦合重整化；殼層體積比；組合因子的列舉（5 個片段、一般式、$f_3$）與連結叢集的數值核對；
  遞迴方程 (12.37)(12.38)／(12.66)(12.67) 的數值積分（流線點列、不動點、線性化矩陣與特徵值）、$u(\ell)$ 的解析解 (12.41)、$\nu$ 的一階公式 vs 數值特徵值（Ising 與 O(n)，$n=0..3$）、Flory 的 ε 展開、冪次數算的指數。
  已知的三維指數（$n=0$：0.588；Ising 0.630 等）只當引用常數。`check()` 用兩種以上獨立算法核對全部場景卡數字，約 10 秒內。
- **重用（不改）**：`walk_draw`、`lattice_draw`、`walk_rg`（第 11 章的流與指數公式）、`walk_crit`（指數表）。
- **畫圖**：不新增模組；流向圖、配對圖、殼層、全書回顧表寫在場景檔的模組層級函式。
- combkit：`walk_eps.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex25_wilson_fisher.py`；升 1.22.0。

## 5. 省 token 檢討

1. **38 頁 → 6 幕**：§12.5–12.6 的圖解法只取概念與兩個可核對的式子（(12.86) 與冪次數算），反項與 Callan–Symanzik 只用文字；補充 12.7 的組合學併進 W1202。
   場景代理人 A 做 W1201–W1203、B 做 W1204–W1206。
2. **重活交給上下文短的代理人**：基礎代理人、2 個場景代理人、1 個講義代理人；我只寫大綱、看總覽圖、接線、輸出、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁；`ch12_formulas.md` 只收場景與講義用得到的式子＋筆誤表（我讀正文時看到的可疑處：p. 293 的「Figure 12.16(b)」應是 12.11(b)？、p. 297 的「$u^*\to\epsilon/36d$」、p. 298 的「$r=t-6ud$」、(12.59)(12.60) 與「三維 Ising $\nu\sim0.61$」（已知約 0.630）、Ex. 12.4 的文字斷掉、p. 320 前後兩次組合因子寫成 $\frac{5!}{3!2!}$ 與 $\frac{5!}{2!1!}$、(12.68)–(12.72)、(12.89) 的 $u^2r^{2d-3}$）。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_eps.check()` assert；流線與特徵值在模組裡算好（有快取）。
5. **不新增畫圖模組、不做 3D**。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 12 節目標 9–10 頁**，最後加一頁「全書回顧」；練習解答 Ex. 12.1–12.11。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log、記憶；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_eps.py` 已完成（`python utils/walk_eps.py` 自我測試通過，約 1.2 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch12_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：$\epsilon=4-d$；GLW 的 $r,u$；布里淵區半徑 $\Lambda=e^{-\ell}$；殼層裡 $r+q^2\to r+1$；
**遞迴方程把書上 (12.37) 的 $d$ 吸收進 $u$**：$\frac{dr}{d\ell}=2r+\frac{12u}{1+r}$、$\frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}$（O($n$)：12 → $4(n+2)$、36 → $4(n+8)$；畫面上第一次出現時用字幕說「書上留著 $d$，本系列併進 $u$」）。
程式代號：`E` = `utils.walk_eps`、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。三維的已知指數只當引用值（`E.NU_3D_CITED`），畫面標「已知（引用）」。

**W1201 加上四次項：怎樣正確地積掉自由度**（§12.1–12.2.4 前半，pp. 285–290，(12.1)–(12.11)、(12.20)–(12.34)，Fig. 12.1–12.5）
1. 問：上一章只有二次項；加上四次項 $u\,S^4$ 之後，還能「積掉一半」嗎？(12.3) 一分量 GLW $\mathcal H=\sum_q\frac12(r+q^2)|S(\vec q)|^2+\frac uN\sum S^4\delta$。
   Fig. 12.1–12.2：正方晶格（格距 $a$，5×5 點）與倒晶格（格距 $2\pi/a$）並排；字幕：實空間越大，$\vec q$ 越密（$L\propto N^{1/d}$）。無數字。
2. 球形布里淵區（Fig. 12.3–12.4，平面示意：半徑 $\Lambda=e^{-\ell}$ 的圓、外圈殼層塗色）：$\ell=0.5$：$\Lambda=0.6065$；(12.24) $\delta\ell=0.01$：$dq=e^{-\ell}-e^{-\ell-\delta\ell}=0.006035\approx e^{-\ell}\delta\ell=0.006065$
   （`E.shell_dq(0.5, 0.01)`）。殼層比例表（`D.tex_table`，3 列；`E.SHELL_CASES`、`E.shell_fraction`、`E.shell_linear`）：$d=3$、$\delta\ell=0.1$：$1-e^{-d\delta\ell}=0.2592$ vs $d\,\delta\ell=0.3000$；
   $d=3$、0.01：0.0296 vs 0.0300；$d=4$、0.01：0.0392 vs 0.0400 ⇒ (12.20)–(12.23) $\sum_{\rm shell}1\approx N\,d\,\delta\ell$。字幕：半徑 200 的球內數格點，結果一樣（差不到 0.0002）。
3. (12.5)–(12.10) 流程圖（四格：分出殼層 → 對 $S'$ 求和得 $A(K')e^{-H'}$ → 重新標度 → $K'=\mathcal R(K)$）。兩變數玩具（書外；`E.TOY` $=(r,u)=(0.2,0.01)$）：
   $\mathcal H=\frac r2s^2+us^4+\frac{r+1}2S'^2+6u\,s^2S'^2$（$s$ 在內區、$S'$ 在殼層）；對 $S'$ 積分（`E.toy_Heff`）：$\mathcal H'(s)=\frac r2s^2+us^4+\frac12\ln\big(1+\frac{12us^2}{r+1}\big)$；
   展開（`E.toy_coeffs(0.2, 0.01)`）：$r'=r+\frac{12u}{1+r}=0.3000$、$u'=u-\frac{36u^2}{(1+r)^2}=0.0075$、新的 $s^6$ 項 0.000167（書上丟掉的高階項）。
   小圖（`D.plot_axes`，$s\in[-2,2]$）：$\mathcal H'(s)$ 與 $\frac{r'}2s^2+u's^4$ 兩條（各 ≤ 40 點，幾乎重合）；$Z=9.0517=A\cdot Z'=2.2882\times3.9558$（`E.toy_Z`、`E.toy_A`）。
   字幕：12 和 36 就是下一幕要在殼層裡數出來的數字。
4. (12.11) 四次項按「幾隻腳在殼層」分五項，係數 $1,4,6,4,1$（$\binom4k$）；三、四隻腳在殼層的兩項在 $\delta\ell\to0$ 時可略。重新標度 (12.29)–(12.34)：$N\to N'e^{d\delta\ell}$、
   $\vec q\to\vec q\,'e^{-\delta\ell}$、$S\to S''e^{\delta\ell}$ ⇒ 乘數（`E.rescale_factors(3, 0.1)`）：$q^2$ 項 1.0000、$r$ 項 $e^{2\delta\ell}=1.2214$、$u$ 項 $e^{(4-d)\delta\ell}=1.1052$。
   筆誤框 1（where「p. 295，(12.32)」）：錯 `N\to N'e^{d\delta\ell}\approx N'(1+\delta\ell)`、對 `N\to N'e^{d\delta\ell}\approx N'(1+d\,\delta\ell)`；字幕：$d=3$、$\delta\ell=0.01$：$e^{0.03}=1.0305$，不是 1.0100。
   筆誤框 2（where「p. 295，(12.33)」）：錯 `\vec q\to\vec q\,'e^{-d\delta\ell}`、對 `\vec q\to\vec q\,'e^{-\delta\ell}`；字幕：照書上 $q^2$ 項要乘 0.6703（$d=3$、$\delta\ell=0.1$），模型就不再同形（`E.rescale_factors(3, 0.1, book=True)[0]`）。
收束：積掉殼層、重新標度之後，只剩 $r$ 與 $u$ 兩個數在變——下一幕數殼層裡的配對。筆誤框 2 個。

**W1202 配對遊戲：只剩連通的圖**（§12.2.4 後半、§12.7，pp. 289–292、317–322，(12.12)–(12.15)、(12:S-1)–(12:S-7)，Fig. 12.6–12.11、12:S-2–12:S-5，Ex. 12.1、12.9–12.11）
1. 問：殼層的變數積掉之後，四次頂點會變成什麼？頂點畫成一個點＋四隻手（Fig. 12.6）：(a) 一隻手在殼層 $\frac{4u}N$、(b) 兩隻手在殼層 $\frac{6u}N$（(12.14)(12.15)）；
   對殼層做 Gaussian 平均＝殼層的手兩兩牽起來（(12:S-2) $\langle S(\vec Q_1)S(\vec Q_2)\rangle=\frac{\delta(\vec Q_1+\vec Q_2)}{r+Q_1^2}$）。字幕：和第 6 章的 Wick 配對同一件事。無數字。
2. 片段（Fig. 12:S-2）：每個 (b) 片段有兩隻殼層手；全部牽好之後連成一個個環：type 1（自己的兩手相牽）、type 2（兩個片段雙線）、type 3（三角環）。
   3 個片段（`E.type_census(3)`）：15 種牽法 $=1+6+8$（三個 type 1／type 1＋type 2／三角環）；不同的片段圖 1、3、1（`E.diagram_census(3)`）。
   $f_3=6$（Ex. 12.11，`E.ring_aut(3)`）：三角形的 $3!$ 種換頂點都不破壞配對；字幕：4 個片段的環只有 8 種換法（`E.ring_aut(4)`），不是 $4!=24$。
3. Fig. 12:S-5（Ex. 12.9）：5 個片段 → 3 個 type 1＋1 個 type 2。列舉：945 種牽法（`len(E.matchings(5))`）裡有 20 種（`E.type_census(5)[(2, 1, 1, 1)]`）＝10 張片段圖×2 種雙線接法；
   $\frac{5!}{3!\,2!}=10$（`E.pair_factor(3, 1)`）。畫 10 張小圖（2×5 排列，每張 5 個點、標出連雙線的那一對）。
   筆誤框 1（where「p. 320」）：錯 `\frac{5!}{2!\,1!}`、對 `\frac{5!}{3!\,2!}=10`；字幕：照書上是 60；同一頁前面寫對了。
4. 一般式 (12:S-4) $\frac{n!}{n_1!\,n_2!\,(2!)^{n_2}}$（Ex. 12.10）表（`D.tex_table`，6 列；`E.pair_factor` 與 `E.diagram_census` 打勾）：$(n_1,n_2)$＝(2,1) 6、(0,2) 3、(3,1) 10、(4,1) 15、(2,2) 45、(0,3) 15，全部＝列舉。
   (12:S-3)(12:S-5)：$\sum\frac{\mathcal T_1^{n_1}\mathcal T_2^{n_2}}{n_1!\,n_2!}=e^{\mathcal T_1+\mathcal T_2}$。
   筆誤框 2（where「p. 322，(12:S-7)」）：錯 `\mathcal T_1^{n_1}\mathcal T_3^{n_3}\mathcal T_3^{n_3}`、對 `\mathcal T_1^{n_1}\mathcal T_2^{n_2}\mathcal T_3^{n_3}`；字幕：$n_1$ 也要從 0 開始加。
5. 連結叢集的數值核對（Ex. 12.1；`E.LC` $=(a,\sigma^2)=(0.1,1)$，$X=-as^2$）：$\langle e^X\rangle$：數值積分＝所有圖的和＝$e^{\text{連通圖}}$＝$(1+2a\sigma^2)^{-1/2}$＝0.9129
   （`E.lc_average(0.1, 1.0, method)`，四種 method）；二階：全部 0.0150 ＝ 連通 0.0100 ＋ 不連通 0.0050（`E.lc_order(0.1, 1.0, 2)`）；連通圖 $\kappa_k$＝−0.1000、0.0200、−0.0080、0.0048（`E.lc_cumulants_rings(0.1, 1.0, 4)`）。
   一行：$\sum_n(2n-1)!!\frac{x^n}{n!}=\exp\big(\sum_kw_k\frac{x^k}{k!}\big)=(1-2x)^{-1/2}$，$w_k=1,2,8,48$（`E.ring_weight`）；$x=0.1$：1.1180（`E.egf_check(0.1)`）；三個殼層變數也一樣：0.8087（`E.lc_multi()`）。
收束：指數化自動把不連通的圖算好——只要畫連通的圖；Fig. 12.11 的兩張一迴圈圖就是下一幕的主角。筆誤框 2 個。

**W1203 一迴圈與 $(r,u)$ 平面的流**（§12.2.4–12.2.5，pp. 292–297，(12.16)–(12.38)，Fig. 12.11–12.12）
1. 問：兩張一迴圈圖各給多少？Fig. 12.11：(a) 一個頂點、殼層線自己接成圈 → 改 $r$；(b) 兩個頂點、兩條殼層線 → 改 $u$。
   (12.18)(12.19)：$\text{(a)}=\frac{6u}{N(r+1)}\sum_{\vec Q}\sum_{\vec q}SS$、$\text{(b)}=\frac{36u^2}{N^2(r+1)^2}\sum_{\vec Q,\vec Q'}\sum S\cdots S\delta$（6 = W1202 的 6u 頂點；36 = 6×6）。
   筆誤框 1（where「p. 293，(12.17) 前」）：錯 `\text{Figure 12.16(b)}`、對 `\text{Figure 12.11(b)}`。
2. $\sum_{\vec Q\in\rm shell}1=N\,d\,\delta\ell$（(12.20)–(12.25)；回扣 W1201 的 0.0392 vs 0.0400）；動量守恆 (12.26) $\vec Q+\vec Q'+\vec q_1+\vec q_2=0$，取 $\vec q_1+\vec q_2=0$ 得 (12.27)；
   (12.28) $r+\frac{12u}{r+1}d\,\delta\ell$、$u-\frac{36u^2}{(r+1)^2}d\,\delta\ell$；乘上 W1201 的重新標度乘數 ⇒ (12.37)(12.38)。
   筆誤框 2（where「p. 296–297，(12.36)(12.38)(12.39)」）：錯 `(4-d)u-36u^2`、對 `(4-d)u-36d\,u^2`；字幕：(12.28)、(12.40)、(12.41) 都有 $d$；本系列把 $d$ 併進 $u$：
   $\frac{dr}{d\ell}=2r+\frac{12u}{1+r}$、$\frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}$（`E.beta(r, u, eps)`）。
3. 不動點：Gaussian $(0,0)$；Wilson–Fisher（`E.fixed_point(eps)`，Fraction）$r^*=-\frac\epsilon{6+\epsilon}$、$u^*=\frac\epsilon{(6+\epsilon)^2}$。表（`D.tex_table`，3 列 × 5 欄，`E.NU_EPS`；一階 `E.fixed_point_lowest`）：
   $\epsilon=1$：$-\frac17=-0.1429$、$\frac1{49}=0.0204$｜一階 $-\frac\epsilon6=-0.1667$、$\frac\epsilon{36}=0.0278$；$\epsilon=\frac12$：$-\frac1{13}=-0.0769$、$\frac2{169}=0.0118$｜−0.0833、0.0139；
   $\epsilon=0.1$：$-\frac1{61}=-0.0164$、$\frac{10}{3721}=0.0027$｜−0.0167、0.0028。字幕：$\epsilon\to0$ 時兩個不動點合而為一（書上 $u^*\to\epsilon/36d$）。
4. $(r,u)$ 平面的流（`D.plot_axes` $r\in[-0.45,0.45]$、$u\in[0,0.05]$；`E.flow_starts()` 6 條，每條 `E.flow_line(r0, u0)` 40 點、末端加箭頭；兩個不動點畫白點）：
   每組中間那條（臨界線，$r_0=r_c(u_0)$）流進 WF 點 $(-0.1429,0.0204)$，旁邊兩條先靠近、再往 $r\to\pm$ 跑掉（停在 $|r|=0.45$）。Gaussian 點的特徵值 2.0000、1.0000（`E.eigen(1.0, 1, "gauss")`）：
   $\epsilon>0$ 時 $u$ 方向也是相關的——Gaussian 點不穩定。
收束：四維以下，臨界的系統都流向 Wilson–Fisher 點——臨界指數由它決定。筆誤框 2 個。

**W1204 ε 展開求 $\nu$**（§12.2.5–12.2.6，pp. 296–300，(12.39)–(12.60)，Fig. 12.13，Ex. 12.2–12.4）
1. 問：耦合 $u$ 最後流到哪裡？(12.39)（$r=0$）$\frac{du}{d\ell}=\epsilon u-36u^2$ ⇒ (12.41) $u(\ell)=\frac{\epsilon u_0e^{\epsilon\ell}}{\epsilon+36u_0(e^{\epsilon\ell}-1)}$（Ex. 12.2）。
   圖（`D.plot_axes` $\ell\in[0,6]$；`E.u_analytic` 兩條各 40 點＋`E.u_numeric` 的 6 個點疊上去；`E.U_CASES`、$\epsilon=1$）：$u_0=0.005$：$\ell=1,2,5$ → 0.0104、0.0172、0.0270；
   $u_0=0.06$：0.0346、0.0300、0.0279 → 虛線 $u^*=\frac1{36}=0.0278$。Fig. 12.13 的數線：$u=0$ 斜率 $+\epsilon$（不穩定，Ex. 12.3；$u_0=10^{-6}$：$u(2)/u_0=7.3874\approx e^2=7.3891$）、$u^*$ 斜率 $-\epsilon$（穩定）。
2. $r=t-6u$（書上 $t-6ud$）⇒ (12.45) $\frac{dt}{d\ell}=2t-12u\,t$ ⇒ (12.49)–(12.51) $t(\ell)=t(0)e^{2\ell}\big[\frac{\epsilon}{\epsilon+36u_0(e^{\epsilon\ell}-1)}\big]^{1/3}$；
   數值 $t_0=10^{-4}$、$u_0=0.01$、$\ell=3$：0.0203（閉式＝數值積分，`E.t_closed(3, 1e-4, 0.01, 1.0)`、`E.t_numeric`）；大 $\ell$：$t\propto e^{(2-\epsilon/3)\ell}$。
   筆誤框 1（where「p. 298，(12.48)」）：錯 `12d\int_{u_0}^{u}\frac{du'}{\epsilon u'-36du'^2}`、對 `12d\int_{u_0}^{u}\frac{u'\,du'}{\epsilon u'-36du'^2}`；
   字幕：$u_0=0.001\to0.005$（$\epsilon=1$、$d=3$）：照書上 81.7804，應為 0.2207＝(12.49)（`E.int_1248`、`E.int_1249`）；(12.47) 左邊也少了 $12d$。
3. (12.52)–(12.59)：$\xi(\ell^*)=1$、$r(\ell^*)=D$ ⇒ $e^{\ell^*}=(T-T_c^*)^{-1/(2-\epsilon/3)}$ ⇒ $\nu=\frac1{2-\epsilon/3}\approx\frac12+\frac\epsilon{12}$；$\epsilon=1$：$\frac7{12}=0.5833$（(12.60)）。
   Ex. 12.4：臨界點移到 $r_c\approx-6u_0$（書上 $-6u_0d$）：$u_0=10^{-4}$、$\epsilon=0.1$：二分法 −0.000631、特徵向量 $-\frac{12u_0}{2-\epsilon}=-0.000632$、一階 $-6u_0=-0.000600$
   （`E.r_critical(1e-4, 0.1)`、`E.r_critical_linear(1e-4, 0.1)`）。
   筆誤框 2（where「p. 299，Ex. 12.4」）：錯 `\text{(12.53) follows from (12.51)}`、對 `\text{(12.52) follows from (12.51)}`；字幕：(12.53) 在練習之後才定義。
4. $\nu$ 的三種讀法（`D.tex_table`，4 列 × 4 欄，Ising；`E.nu_first`、`E.nu_unexpanded`、`E.nu_eigen`）：$\epsilon=0.1$、0.5、1：一階 $\frac12+\frac\epsilon{12}$：0.5083、0.5417、0.5833；
   不展開 $\frac1{2-\epsilon/3}$：0.5085、0.5455、0.6000；線性化矩陣 $1/\lambda_r$：0.5076、0.5256、0.5272。$\epsilon=1$ 的矩陣 $\begin{pmatrix}5/3&14\\1/21&-1\end{pmatrix}$，特徵值 1.8968、−1.2301（`E.eigen(1.0)`）。
   小圖：三條 $\nu(\epsilon)$，$\epsilon\in[0,1]$（各 ≤ 30 點）在 $\epsilon=0$ 相切。字幕：三者只在 $\epsilon$ 的一次相同——一迴圈只準到 $O(\epsilon)$。
5. 和已知值比：平均場 0.5 → 一階 0.5833 → 三維 Ising 已知約 0.630（引用）：大幅改善。
   筆誤框 3（where「p. 299，(12.60) 後」）：錯 `\nu\sim0.61`、對 `\nu\approx0.630`；字幕：0.630 是文獻值（只引用）。
收束：Wilson–Fisher 點把 $\nu$ 從平均場的 $\frac12$ 推到 0.58——第一次「算」出不是 $\frac12$ 的臨界指數。筆誤框 3 個。

**W1205 O($n$) 與自迴避漫步**（§12.3–12.4，pp. 300–304，(12.61)–(12.72)，Fig. 12.14–12.17，Ex. 12.5–12.8）
1. 問：自旋有 $n$ 個分量時，12 和 36 會變成什麼？(12.61) $\frac uN\sum(\vec S\cdot\vec S)(\vec S\cdot\vec S)\delta$，波浪線表示 $u$（Fig. 12.14，兩個點積各畫一條實線）；
   Ex. 12.5 (12.63)：$\langle S^\alpha(\vec Q_1)S^\beta(\vec Q_2)\rangle=\frac{\delta_{\alpha\beta}}{1+r}\delta(\vec Q_1+\vec Q_2)$。
2. $r$ 的兩張圖（Fig. 12.15）：殼層的兩隻手在同一個點積 → 繞一圈 $\sum_\beta\delta_{\beta\beta}=n$（2 種）；在不同點積 → 1（4 種）⇒ $2n+4$（(12.64)，Ex. 12.6）；$u$ 的圖（Fig. 12.16–12.17，Ex. 12.7）→ $4(n+8)$。
   表（`D.tex_table`，4 列；`E.on_counts(n)`、`E.AB(n)`，$n=0..3$）：$2n+4$：4、6、8、10；$A=4(n+2)$：8、12、16、20；$B=4(n+8)$：32、36、40、44。字幕：程式真的做了 $n=1..4$ 的張量縮併，$n=1$ 回到 12、36。
3. (12.66)(12.67) ⇒ (12.68) $\nu=\frac12+\frac{(n+2)\epsilon}{4(n+8)}$；$n=0$（自迴避）：$\frac12+\frac\epsilon{16}$（(12.69)），$\epsilon=1$：$\frac9{16}=0.5625$（(12.70)；書上的 0.563 是四捨五入）。
   表（`D.tex_table`，4 列 × 4 欄，$\epsilon=1$）：$n=0,1,2,3$：一階 0.5625、0.5833、0.6000、0.6136（`E.nu_first(1, n)`）；不展開 0.5714、0.6000、0.6250、0.6471（`E.nu_unexpanded(1, n)`）；
   已知（引用，`E.NU_3D_CITED`）0.588、0.630、0.672、0.711。WF 點 $n=0$：$(-\frac19,\frac2{81})$（`E.fixed_point(Fraction(1), 0)`）。
4. Flory（(12.71)(12.72)）：$\nu_F=\frac3{d+2}=\frac35=0.6$；$\frac3{6-\epsilon}=\frac12+\frac\epsilon{12}+\frac{\epsilon^2}{72}+\frac{\epsilon^3}{432}+\cdots$（`E.flory_series(k)`）；一次項 $\frac1{12}$ vs $n=0$ 的 $\frac1{16}$ ⇒ Flory 不是精確的
   （但數值接近：0.5833 vs 0.5625）。字幕：Ising 的 $\frac1{2-\epsilon/3}$ 恰好就是 $\frac3{6-\epsilon}$。三根長條（`E.saw_size(1000, nu)`）：$N=1000$ 步：$\sqrt N=31.62$、$N^{9/16}=48.70$、$N^{0.588}=58.08$。
5. ε 展開的四個隱憂（p. 304，四行字）：$\epsilon=1$ 不小；級數是漸近的（收斂半徑 0）；可能有 $e^{-A/\epsilon}$ 這種冪級數抓不到的項；重整化群本身沒有被嚴格證明。無數字。
收束：同一套遞迴，$n\to0$ 就是自迴避漫步——第 1 章的 $\sqrt N$ 在這裡變成 $N^{0.5625}$（已知 $N^{0.588}$）。筆誤框 0 個。

**W1206 圖解法的重整化群＋全書回顧**（§12.5–12.6，pp. 304–317，(12.73)–(12.89)，Fig. 12.18–12.28）
1. 問：不切殼層、直接用費曼圖，也能做重整化嗎？(12.73) $\frac1{r+q^2}+\frac1{r+q^2}\Sigma\frac1{r+q^2}+\cdots=\frac1{r+q^2-\Sigma}$（Fig. 12.18–12.20：質量算子 $\Sigma$、雙線＝重整化傳播子）；
   動量依大小排序、截斷 $q<q_0<q_1,q_2,q_3$（Fig. 12.21 西瓜圖示意）。字幕：和第 7 章 W704 的 Dyson 方程同一個形狀。無數字。
2. 四點函數 $\Gamma_4$ 與 $u_R=q_0^{4-d+2\eta}g$（(12.81)(12.82)）；(12.86) $q_0\frac{dg_0}{dq_0}=(-4+d-2\eta)g_0+\mathcal Bg_0^2$；換變數 $q_0=e^{-\ell}$ ⇒ $\frac{dg_0}{d\ell}=\epsilon g_0-\mathcal Bg_0^2$ ＝ (12.39)（$\eta=0$）。
   數值（`E.g_numeric(math.exp(-5), 0.005, 1.0)`）：$g_0(1)=0.005$ 積到 $q_0=e^{-5}=0.0067$：0.0270 ＝ W1204 的 $u(\ell=5)$；不動點 $g^*=\frac1{36}=0.0278$。
3. 冪次數算（Fig. 12.27、(12.89)）：第 $n$ 階含 $(ur^{-(4-d)/2})^n$，$d<4$、$r\to0$ 時發散；$u_R\propto r^{(4-d)/2}$ 剛好抵消。$Q\propto r^{1/2}$：三個 $d^dQ$ 給 $r^{3d/2}$、$\delta$ 給 $r^{-d/2}$、三個 $G$ 給 $r^{-3}$ ⇒ $u^2r^{d-3}=r\,(ur^{-(4-d)/2})^2$（`E.pc_exponent(d)`、`E.pc_rhs(d)`）。
   筆誤框 1（where「p. 315，(12.89) 後」）：錯 `u^2r^{2d-3}=r^{-1}(ur^{-(4-d)/2})^2`、對 `u^2r^{d-3}=r\,(ur^{-(4-d)/2})^2`；
   字幕：數值積分量出的冪次 $d=2$：−1.0000、$d=2.5$：−0.5000（`E.pc_numeric(2)`、`E.pc_numeric(2.5)`）；照書上是 +1、+2。
   一行字：反項（Fig. 12.28）與 Callan–Symanzik 方程——截斷 $q_0$ 是算法的產物，答案不能依賴它（只用文字）。
4. 第 12 章回顧（換頁；表 6 列，沒有表頭）：W1201 $r'=r+\frac{12u}{1+r}$、$u'=u-\frac{36u^2}{(1+r)^2}$（玩具 0.3000、0.0075）；W1202 只剩連通圖、$\frac{5!}{3!2!}=10$、$f_3=6$；W1203 WF 點 $(-\frac17,\frac1{49})$；
   W1204 $\nu=\frac12+\frac\epsilon{12}=0.5833$；W1205 $\nu=\frac12+\frac{(n+2)\epsilon}{4(n+8)}$、$n=0$：0.5625；W1206 (12.86)＝(12.39)、$u^2r^{d-3}$。
5. 全書回顧（換頁，可分兩頁各 6 列；`E.BOOK_RECAP`（章, 主角, Tex）＋`E.recap_numbers()`＋`E.scene_range(ch)`；左欄中文用 `zh()`、右欄 Tex；不載入影片）：
   1 一維漫步 W101–W112：$\langle x^2\rangle=N$；$P_{10}(0)=\frac{252}{1024}=0.2461$ vs 高斯 0.2523；2 生成函數與格子 W201–W212：回原點 1、4、36、400、4900、63504；
   3 回家問題 W301–W308：一、二維必回家，sc 回家機率 0.3405；4 牆與邊界 W401–W408：鏡像；$x_0=2$、$L=5$ 先到左牆 3/5；
   5 偏向與持續 W501–W507：$l_p=\frac Q{1-Q}l$；$Q=\frac12$、$N=1000$：2996；6 漫步的形狀 W601–W607：$A_3=\frac{10}{19}=0.5263$；
   7 自迴避 W701–W706：$c_{14}=2374444$、$\nu_F=\frac3{d+2}$；8 長尾與標度 W801–W804：上臨界維度 4、$\mu=2.638$；
   9 臨界現象 W901–W906：$n\to0$ 的 O($n$)＝自迴避，Flory 0.600 vs 0.588（差 2%，0.0204）；10 受限的高分子 W1001–W1005：$\rho=|r|/2u$、$\nu'=\frac12$；
   11 重整化群 W1101–W1106：$\frac{\ln4}{\ln3}=1.2619$、$\nu=1/\Lambda_r$；12 ε 展開 W1201–W1206：0.5833、0.5625。字幕：全書 87 幕（`sum(E.BOOK_SCENES.values())`）。
收束：從一枚硬幣到一個 ε——漫步、磁鐵與重整化群說的是同一個故事。筆誤框 1 個。

## 7. 製作規則（代理人必讀）

- ManimGL 1.7.2（`from manimlib import *`），**不是** Manim Community（用 `ShowCreation` 不用 `Create`）。繼承 `ChuScene`；`construct` 只呼叫**模組層級函式**
  （連續播放會借用 `construct`，寫成方法會 AttributeError）。開頭 `self.show_header("…")`＋`check_math()`，結尾收束字幕後 `self.wait(1.0)`。
- 場景骨架：
  ```python
  import math
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
  from utils import walk_eps as E
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_eps` 算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串，
  分數用 `Fraction`；Fraction 不能直接 `:.4f`，要 `float()`）；版面用 `assert_in_frame`／`assert_apart`。
  不要在場景裡呼叫 `E.check()`（整套自我測試由 `check_env.py` 跑）；`E.summary()` 約 1 秒，可以不用（直接呼叫卡上的函式）。
  最慢的呼叫：`E.r_critical` 每個約 0.2 秒（`E.flow_starts()` 叫兩次，有快取）、`E.toy_Z` 約 0.1 秒；其餘都很快。
  scipy 的 IntegrationWarning 不影響結果（`check_math` 裡用 `warnings.catch_warnings()` 包起來）。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W1201 2 個、W1202 2 個、W1203 2 個、W1204 3 個、W1205 0 個、W1206 1 個。
  框裡的 Tex 不放中文（英文字用 `\text{}`）。卡上標「用字幕」的書上問題只用一行字幕。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。字幕裡的「÷」寫成「除以」。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；矩陣用 `\begin{pmatrix}…\end{pmatrix}`；
  ManimGL 的 TeX 模板載入 tipa，避免 `\!\!`。`get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。
  標題的 `Text` 不能有上標或希臘字母以外的數學（寫「ε 展開」「n 分量」之類）。
- 名稱：`L`（lattice_draw）、`D`（walk_draw）、`E`（walk_eps）是模組——不要拿來當變數；**(12.53) 的常數 $D$ 在程式裡叫 `Dist`**，能量、期望值不要叫 `E`。
- 曲線一條＝一個 VMobject（`D.curve` 或 `set_points_as_corners`），每條 ≤ 60 點，刻度數字自己加。流線用 `E.flow_line` 的 40 點、末端一個小箭頭。
  頂點／片段／環的小圖：點＋短線（每張 ≤ 12 個物件），10 張小圖排成 2×5；布里淵區、殼層、費曼圖都是平面示意，不做 3D。
- 大量物件：長條最多 60 根（同一畫面）、表格最多 6 列 × 10 欄、點最多 60 個、箭頭最多 64 支；避免數千個小物件。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式，回報時說明。每幕 250–450 行（W1206 可到 550）。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W1201 | `w1201_quartic_shell.py` | `W1201QuarticShell` | A |
| W1202 | `w1202_linked_clusters.py` | `W1202LinkedClusters` | A |
| W1203 | `w1203_one_loop_flow.py` | `W1203OneLoopFlow` | A |
| W1204 | `w1204_epsilon_nu.py` | `W1204EpsilonNu` | B |
| W1205 | `w1205_on_saw.py` | `W1205OnSaw` | B |
| W1206 | `w1206_diagrammatic_recap.py` | `W1206DiagrammaticRecap` | B |

連續播放：`main.py` 的 `Ch12Demo`（W1201–W1206）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch12/part12}`；內文 `notes/ch12/part12.tex`，圖表與數字檢查 `notes/ch12/make_figs.py`
（`check_numbers_ch12()` 逐一 assert，所有數字都從 `codes/utils/walk_eps.py` 算；表格與 TikZ／pgfplots 圖產生成 `notes/ch12/figs.tex`；
巨集都以 `\Twelve` 開頭、`\label` 用 `twelve-…`；執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–11 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch12_formulas.md`、
`walk_eps.py` 的 docstring、`notes/ch11/part11.tex` 與 `notes/ch11/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 9–10 頁＋最後一頁「全書回顧」，圖表 6 組：
(1) §12.1–12.2.4 前半：布里淵區殼層（TikZ）＋殼層比例表（`E.SHELL_CASES`，含格點計數 `E.shell_fraction_lattice`）＋兩變數玩具：$\mathcal H'(s)$ 與 $\frac{r'}2s^2+u's^4$（pgfplots，各 ≤ 40 點）、
$r'=0.3000$、$u'=0.0075$、$s^6$ 係數 0.000167、$Z=9.0517=2.2882\times3.9558$＋重新標度乘數（含書上 (12.33) 的 0.6703）；
(2) §12.2.4 後半與補充 12.7：type 1／2／3 的示意（TikZ）＋配對普查表（`E.type_census(n)`、`E.diagram_census(n)`，$n=3,4,5$）＋(12:S-4) 表（6 列）＋$f_k$、$w_k$＋連結叢集的數值
（0.9129 四種算法、二階 0.0150＝0.0100＋0.0050、$\kappa_k$、1.1180、0.8087）；
(3) §12.2.4–12.2.5：不動點表（3 個 $\epsilon$，精確分數＋一階）＋$(r,u)$ 流線（pgfplots，`E.flow_starts()`、`E.flow_line`）＋Gaussian 點的特徵值；
(4) §12.2.5–12.2.6：$u(\ell)$ 兩條（pgfplots）＋$t(\ell)$ 的數值＋(12.48) 的筆誤數值＋$\nu$ 三種讀法對 $\epsilon$ 的曲線（pgfplots，≤ 30 點）與表＋$\epsilon=1$ 的矩陣與特徵值＋Ex. 12.4 的 $r_c$；
(5) §12.3–12.4：$2n+4$、$A$、$B$ 表＋$\nu$ 的 O($n$) 表（一階／不展開／特徵值 `E.nu_eigen(1.0, n)`／引用值）＋Flory 的係數＋$N^\nu$ 三個數；
(6) §12.5–12.6：$g_0(q_0)$ 與 $u(\ell)$ 疊在一起的曲線（pgfplots，橫軸 $\ell=-\ln q_0$）＋冪次表（$d=2,\frac52,3,4$：$d-3$、書上 $2d-3$、數值 −1.0000、−0.5000）。
最後一頁「全書回顧」：12 列的表（章、主角、場景代號 `E.scene_range(ch)`、關鍵結果與數字，數字取 `E.recap_numbers()`），附一句：87 幕、從 $\sqrt N$ 到 $N^{0.588}$。
小節對應書上 12.1–12.2.3、12.2.4、12.7（補充）、12.2.5、12.2.6、12.3–12.4、12.5–12.6，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch12_formulas.md` 的筆誤表，11 列；
另附〔注〕：本系列把 $d$ 吸收進 $u$、(12.35) 少負號、(12.63) 多一個 $S($、(12.64) 少 $\frac uN$、$r=t-6ud$ 與 $u^*\to\epsilon/36d$ 是對的、(12.70) 的四捨五入、$f_k=2k$（$k\ge4$）、
(12.54) 不展開等於 Flory（Ising）、完整一迴圈矩陣的 0.5272 只有一次項可信、(12.86) 的 $q_0=e^{-\ell}$）。
最後練習 12.1–12.11 的解答：
12.1 $n$ 個片段的 Gaussian 平均＝所有配對的和；每個配對把片段分成環（type $k$），同一種分法的數目 $\frac{n!}{\prod_kn_k!\,f_k^{n_k}}$，配上原本的 $\frac1{n!}$ 就是 $\prod_k\frac{\mathcal T_k^{n_k}}{n_k!}$，
對所有 $n_k$ 求和 ⇒ $\exp(\sum_k\mathcal T_k)$；數值：$X=-as^2$ 時 0.9129（四種算法）、二階 0.0150＝0.0100＋0.0050；三個變數 0.8087；
12.2 部分分式 $\frac1{\epsilon u-36du^2}=\frac1\epsilon\big(\frac1u+\frac{36d}{\epsilon-36du}\big)$ ⇒ $\frac1\epsilon\ln\frac{u(\epsilon-36du_0)}{u_0(\epsilon-36du)}=\ell$ ⇒ (12.41)；數值 $u_0=0.005$：0.0104、0.0172、0.0270；
12.3 另一個不動點 $u=0$：$\delta u'=\epsilon\,\delta u$ ⇒ $\delta u\propto e^{\epsilon\ell}$ 變大（不穩定）；$u^*$ 處斜率 $-\epsilon$（穩定）⇒ Fig. 12.13 的箭頭；$u_0=10^{-6}$：$u(2)/u_0=7.3874\approx e^2$；
12.4 大 $\ell$ 時 $\big[\frac\epsilon{\epsilon+36du_0(e^{\epsilon\ell}-1)}\big]^{1/3}\to\big(\frac\epsilon{36du_0}\big)^{1/3}e^{-\epsilon\ell/3}$ ⇒ $t\propto t(0)e^{(2-\epsilon/3)\ell}$（題目的 (12.53) 應為 (12.52)）；
$r^*(0)=t(0)=r(0)+6u_0d\propto T-T_c^*$ ⇒ 臨界點移到 $r_c\approx-6u_0d$（本系列 $-6u_0$；數值 $-0.000631$ vs 特徵向量 $-0.000632$ vs $-0.000600$）；
12.5 殼層的權重 $\exp[-\frac12(1+r)\sum_\alpha S^\alpha(\vec Q)S^\alpha(-\vec Q)]$ 對各分量獨立 ⇒ $\langle S^\alpha S^\beta\rangle=\frac{\delta_{\alpha\beta}}{1+r}\delta(\vec Q_1+\vec Q_2)$；
12.6 四隻手選兩隻在殼層（6 種）：同一個點積（2 種）$\sum_\beta\delta_{\beta\beta}=n$，不同點積（4 種）各 1 ⇒ $2n+4$；乘上 $\frac12$ 前的 2 ⇒ $A=4(n+2)$（$n=1$：12）；
12.7 兩個頂點各選兩隻手、互接（$6\times6\times2=72$ 種）；對外面四隻手對稱化後是 $(\vec S\cdot\vec S)^2$ 的 $8(n+8)$ 倍，乘二階的 $\frac12$ ⇒ $4(n+8)$（$n=1$：36；einsum 核對 $n=1..4$）；
12.8 把 12 → $4(n+2)$、36 → $4(n+8)$：$t\propto e^{(2-\frac{(n+2)\epsilon}{n+8})\ell}$ ⇒ $\nu=\frac1{2-(n+2)\epsilon/(n+8)}\approx\frac12+\frac{(n+2)\epsilon}{4(n+8)}$；$n=0$、$\epsilon=1$：$\frac9{16}=0.5625$；
12.9 先排 5 個片段（$5!$），前三個各自成 type 1、後兩個成 type 2；前三個互換（$3!$）、後兩個互換（$2!$）不產生新的圖 ⇒ $\frac{5!}{3!\,2!}=10$（列舉：945 種配對裡 20 種＝10 張圖×2 種雙線接法；書上的 $\frac{5!}{2!1!}=60$ 是筆誤）；
12.10 $n=n_1+2n_2$ 個片段排成一列：$n!$；type 1 之間互換 $n_1!$、type 2 之間互換 $n_2!$、每個 type 2 內部互換 $2!$ ⇒ $\frac{n!}{n_1!\,n_2!\,(2!)^{n_2}}$（$n\le6$ 全部＝列舉）；
12.11 三角環的 $3!=6$ 種換頂點都保持「每個頂點和另外兩個各連一線」⇒ $f_3=6$（三角環的腳配對 8 種＝$2^3$；〔注〕$k\ge4$ 的環 $f_k=2k$，例如 $f_4=8$）。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點、表的數字由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；講義裡「÷」可以用 `$\div$`；TikZ 從一點畫回同一點不要用 `to[looseness]`。
