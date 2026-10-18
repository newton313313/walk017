# 《Elements of the Random Walk》第 9 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 9 章 Scaling of walks and critical phenomena，pp. 203–232（PDF 頁序 219–248）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_crit.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch08.md` §6–§9。

## 1. 草稿大綱（6 幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W901 | 漫步的指數就是臨界指數：$\langle R^2\rangle\propto N^{2\nu}$ (9.1)(9.2)，Flory 在 1、2、4 維精確、三維差約 1%；de Gennes：$\nu$ 是一個假想磁鐵的指數；磁鐵的自發磁化 $M_s\propto(T_c-T)^\beta$（Fig. 9.1） | §9.1–9.2 開頭，pp. 203–205，(9.1)–(9.3) |
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
