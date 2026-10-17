# 《Elements of the Random Walk》第 4 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 4 章 Boundary conditions, steady state, and the electrostatic analogy，pp. 69–94（PDF 頁序 85–110）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_bc.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch03.md` §6–§9。

## 1. 草稿大綱（8 幕）

| 幕 | 主題 | 書上 |
|---|---|---|
| W401 | 鏡像走路人：碰到牆的路線 ↔ 從鏡像點出發的路線（反射原理） | §4.1.1，pp. 70–73，(4.1)–(4.8)，Fig. 4.1 |
| W402 | 夾在兩面吸收牆之間：正弦模態取代平面波 | §4.1.2，pp. 73–76，(4.9)–(4.15)，Ex. 4.1–4.4 |
| W403 | 短時間像沒牆、長時間指數衰減 $e^{-\pi^2N/2L^2}$；熵與「漏失」的解釋 | §4.1.3，pp. 76–79，(4.16)–(4.27)，Fig. 4.2 |
| W404 | 關在球裡（$j_0=\sin x/x$）與立方體（Ex. 4.5） | §4.1.4，pp. 79–82，(4.28)–(4.37) |
| W405 | 穩定態＝靜電學：源＋漫步 → Poisson 方程；Fick 定律；兩牆之間的「帳篷」 | §4.2–4.2.1，pp. 82–86，(4.38)–(4.49)，Ex. 4.6 |
| W406 | 點源的 $1/r$（把高斯對時間積分）與吸收牆的鏡像源（偶極、最後全被吸收） | §4.2.2–4.2.3，pp. 86–89，(4.50)–(4.60) |
| W407 | 吸收球：$c=c_\infty(1-r_0/r)$、流量 $\propto r_0$ → 細胞大小的上限；電容類比 | §4.2.4–4.2.5，pp. 89–91，(4.61)–(4.66) |
| W408 | 布滿受器的球（Berg–Purcell）＋吸收邊界條件的推導＋第 4 章回顧 | §4.2.6、§4.3，pp. 91–94，(4:S-1)–(4:S-3) |

## 2. 生動度檢討

- 主例延續第 1 章的一維格子漫步：W401 從 $x_0=3$ 出發、牆在 0；鏡像點 $-3$。「碰過牆的路線」與「從 $-3$ 出發的路線」逐條配對（第一次碰牆後鏡射），先猜再數。
- 和第 3 章呼應：W401 一維「從不碰牆」機率 $\approx x_0\sqrt{2/(\pi N)}$ 趨近 0（一維一定會碰牆）；W406「三維點源旁的吸收平面最後吃掉所有人」——因為 $x$ 座標本身是一維漫步。
- W403 用對數圖：存活路線比例 vs $N$，長時間是一條直線（斜率 $\ln\cos(\pi/(L+1))$，連續極限 $-\pi^2/2L^2$）；先問「牆距加倍，衰減變慢幾倍？」（4 倍）。
- W405 生活化：源源不絕放人（像一直噴香水），穩定後的濃度＝靜電位；一維兩牆之間是「帳篷」形；Ex. 4.6 就是賭徒破產：被左牆吸收的機率 $1-x_0/L$。
- W407 細胞問題：營養流入 $\propto r_0$、需求 $\propto r_0^3$ → 存在最大半徑；W408 Berg–Purcell 的驚喜：受器只占表面極小比例就能接近全吸收（書外補充，標明出處 Berg & Purcell 1977）。
- 筆誤紅框照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W401 | $x_0=1$、$N=4$：16 條路線裡不碰牆的條數＝無限制－鏡像 | $x_0=3$、$N=10,100$：精確存活機率 vs $x_0\sqrt{2/\pi N}$；書上三維 $\sqrt{6x_0^2/\pi N}$ |
| W402 | $L$ 小（例如牆在 0 與 $L=6$）：轉移矩陣的特徵值 $2\cos(m\pi/L)$ | 模態和＝精確計數（某個 $N$） |
| W403 | $\sum_{\text{奇數}}1/m^2=\pi^2/8$ | 存活比例的對數圖；短時間 $\approx L$、長時間 $\approx(8L/\pi^2)e^{-\pi^2N/2L^2}$ |
| W404 | $j_0$ 的零點 $k_n=n\pi/R$ | 球與立方體的衰減率比較 |
| W405 | 帳篷 $c(x)$；$x_0=2$、$L=5$：被左牆吸收 $3/5$（精確格子計算） | Poisson ↔ 靜電；Fick 定律 |
| W406 | $\int_0^\infty(4\pi Dt)^{-3/2}e^{-r^2/4Dt}dt=1/(4\pi Dr)$ | sc 格子 Green 函數 $\times r$ 趨近常數；吸收牆的偶極場 |
| W407 | 流量 $4\pi Dr_0c_\infty$；從 $R$ 出發撞到球的機率 $r_0/R$ | 細胞最大半徑的數量級（標明假設） |
| W408 | Berg–Purcell 半飽和受器數 $n=\pi r_0/a$ | 吸收邊界 $G(-1)=0$ 的格子推導 |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_bc.py`（純 Python，可 import 同資料夾的 `walk_gf`、`walk_rec`）：反射原理（暴力列舉 vs 鏡像公式）、
  兩牆之間的精確計數（轉移矩陣／動態規劃）vs 正弦模態和、長短時間的極限、球與立方體的模態（$j_0$ 零點）、一維穩定態帳篷與吸收機率（線性方程精確解 vs 公式）、
  點源時間積分、sc 格子 Green 函數（可用第 3 章的 Watson 計算法）、鏡像源與偶極、吸收球流量與命中機率 $r_0/R$、Berg–Purcell。`check()` 核對書上的數字與全部場景卡數字，至少兩種獨立算法。
- **重用（不改）**：`walk_draw`（數線、走法圖、長條、`plot_axes`、表格、筆誤框）、`lattice_draw`（格點、`heat_cells`、`lattice3d`、`trail`）、`walk_gf`、`walk_rec`。
- **畫圖**：不新增模組；缺的小元件（鏡射路線、等位線、圓形吸收體）寫在場景檔的模組層級函式。
- combkit：`walk_bc.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex17_image_walker.py`；升 1.14.0。

## 5. 省 token 檢討

1. **26 頁 → 8 幕**（第 1、2 章 12 幕）：§4.2.2＋§4.2.3 合一幕、§4.2.4＋§4.2.5 合一幕、§4.2.6＋§4.3＋回顧合一幕；(4.4)(4.6)(4.7) 的換變數細節只放講義。2 個場景代理人 × 4 幕。
2. **我不在這個長對話裡做重活**：這個對話的上下文已很長，所以「讀原文細節、截圖校對、寫模組、補定稿卡、combkit」整包交給**一個基礎代理人**（新的短上下文）；
   場景交給 2 個代理人；講義第 4 節交給 1 個代理人（和場景代理人同時跑）。我只做：大綱（本節）、看總覽圖、接線、1080p、用腳本改文件、git。
3. **只截有疑問的頁面**：基礎代理人最多截 8 頁（110 dpi）；`ch04_formulas.md` 只收場景與講義用得到的式子＋筆誤表。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_bc.check()` assert；場景與講義代理人都**不讀原文**。
5. **不新增畫圖模組、不做 3D 動畫**（球、立方體用 2D 剖面與 `lattice3d` 斜投影）。
6. **少一輪重渲染**：代理人的最後一次 render_check 必須是最終版（改過就重跑，每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 4 節目標 8–10 頁**，圖 4–5 張，練習解答 Ex. 4.1–4.6。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_bc.py` 已完成（`python utils/walk_bc.py` 自我測試通過，約 1 秒）：下面每個數字都在它的 `check()` 裡 assert 過（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。
書上的式子與筆誤見 `ref/md/ch04_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$、起點 $x_0$、吸收牆「踩到就消失」、鏡像點 $-x_0$、
兩面牆在 $0$ 與 $L$（內點 $1..L-1$；格子特徵值 $2\cos(m\pi/L)$——**本章定稿用這個約定**，§2 寫的 $\ln\cos(\pi/(L+1))$ 是「$L$ 個內點」的另一種約定，不用）、
穩定態濃度 $w(\vec r)$（格子）／$c(\vec r)$（連續）、擴散係數 $\beta=D$、源強度 $s_0$。主例：第 1 章的一維醉漢；「0 是懸崖」。
程式代號：`BC` = `utils.walk_bc`、`R` = `utils.walk_rec`、`GF` = `utils.walk_gf`、`L` = `utils.lattice_draw`、`D` = `utils.walk_draw`（見 §7；注意 `L` 是模組，牆距寫 `Lw` 之類的變數名）。
書外的內容（Berg–Purcell、Faber–Krahn、圓盤電容）在畫面上標「書外」並寫出處。

**W401 鏡像走路人**（§4.1.1，pp. 69–73，(4.1)–(4.8)，Fig. 4.1）
1. 問：醉漢從 $x_0=1$ 出發，0 是懸崖（踩到就掉下去）。走 4 步還活著的路線有幾條？（先猜。）
2. 16 條小走法圖（`BC.walks(4, 1)`，4×4 排，`D.walk_axes`＋`D.walk_path`，0 那條水平線畫成紅色懸崖）：碰過 0 的 10 條標紅、沒碰過的 6 條標綠
   （`BC.touches_wall`；`BC.survivors_brute(4, 1)` = {1: 2, 3: 3, 5: 1}）。16 條裡終點 $>0$ 的 11 條，其中 5 條碰過牆（`BC.reflection_pairs(4, 1)`）。
3. 鏡射動畫（Fig. 4.1）：路線 $1\to0\to1\to2\to3$ 的「第一次碰牆之前」那段對 0 鏡射 → 從 $-1$ 出發的 $-1\to0\to1\to2\to3$（`BC.reflect_first`）；之後兩條重合。
   一一對應：碰過牆、終點 $>0$ 的 5 條 ↔ 從 $-1$ 出發、終點 $>0$ 的 5 條。(4.1) $C'(N;x,x_0)=C(N;x-x_0)-C(N;x+x_0)$（`BC.survivors_image`），表三列：
   $x=1$：$6-4=2$；$x=3$：$4-1=3$；$x=5$：$1-0=1$ → 共 6（$=11-5$）。
4. 主例 $x_0=3$：存活機率（`BC.survival_prob`）vs $x_0\sqrt{2/\pi N}$（`BC.survival_asym`）：$N=10$：$21/32=0.65625$ vs 0.7569；$N=100$：0.2356 vs 0.2394；$N=1000$：0.07557 vs 0.07569（→ 0：一維一定掉下去）。
   書上的三維 (4.5)：$\sqrt{6x_0^2/\pi N}$；sc 格子、吸收平面精確（`BC.survival_3d(N, 3)`）：$N=100$：0.3958 vs 0.4146；$N=1000$：0.1305 vs 0.1311
   （一行字幕：$x$ 座標只在 1/3 的步數裡動，所以三維碰平面還是一維的問題）。表格 `D.tex_table`。
5. (4.8) 一行：牆吃掉的走法數 $\mu^N\sqrt{2N/3\pi}$：「離牆 $\sqrt N$ 之內出發的，大多活不久」。無筆誤框（本節的式子都對；$x_0^2\gg2N$ 只是數量級，只放講義）。
收束：碰過牆的路線，正好被從鏡像點出發的路線抵銷——數活著的路線只要做一次減法。

**W402 夾在兩面吸收牆之間：正弦模態**（§4.1.2，pp. 73–76，(4.9)–(4.15)，Ex. 4.1–4.4）
1. 問：懸崖在 0 和 6，醉漢從正中間 3 出發。走 10 步都沒掉下去的路線有幾條？（先猜。）
2. 逐步遞迴（數線 0..6，`D.number_line`，每格上方寫數字；`BC.box_dp(N, 6, 3)`）：$N=0..4$：$\{3{:}1\}$、$\{2{:}1,4{:}1\}$、$\{1{:}1,3{:}2,5{:}1\}$、$\{2{:}3,4{:}3\}$、$\{1{:}3,3{:}6,5{:}3\}$（共 12）
   → 每兩步乘 3 → $N=10$：$81,\ 0,\ 162,\ 0,\ 81$（共 324 條，總共 1024 條）。
3. 為什麼是 3：轉移矩陣（5×5，鄰居相加）的特徵值 $2\cos(m\pi/6)=\sqrt3,\ 1,\ 0,\ -1,\ -\sqrt3$（`BC.box_eigs(6)`：1.7321、1.0000、0.0000、−1.0000、−1.7321），
   特徵向量 $\sin(m\pi x/6)$（(4.9)(4.10)：牆上自動是 0；畫五條正弦曲線在格點上的值）。從中間出發只有 $m=1,5$（和特徵值 0 的 $m=3$）有份：
   $(\sqrt3)^{10}=243$，$x=3$：$\frac26(243+243)=162$；$x=1$：$\frac26(\frac12\cdot243+\frac12\cdot243)=81$（`BC.box_modes`）。
   格子版 (4.15)：$C=\sum_m\frac2L\sin\frac{m\pi x_0}{L}\sin\frac{m\pi x}{L}\big(2\cos\frac{m\pi}{L}\big)^N$；連續極限 $2\cos q\approx2-q^2$ → (4.14)(4.15) 的分母 $1-2z+zq^2$。
4. 第三種算法（Ex. 4.3）：無窮多個鏡像點，正的在 $3+12k$（…, −9, 3, 15, …）、負的在 $-3+12k$（…, −15, −3, 9, …）（`BC.box_image_points(6, 3)`、`BC.box_images`）——也得 81、162、81。
5. 筆誤框（`D.erratum_box`，where「p. 75，(4.12)(4.13)」）：錯 `\sqrt{2/L}\,\sin q x_1`、對 `\sqrt{2/L}\,\sin q x_2`；
   字幕：遞迴移動的是 $x_1$，展開係數裡留下的是另一端 $x_2$（(4.14) 已經寫對）。
收束：有牆的時候，正弦波取代平面波：每個模態每走一步就乘上 $2\cos(m\pi/L)$。

**W403 短時間像沒牆、長時間指數衰減**（§4.1.3，pp. 76–79，(4.16)–(4.27)，Fig. 4.2）
1. 問：兩面牆的距離加倍，「還活著的比例」衰減變慢幾倍？2 倍？4 倍？（先猜。）
2. 總數 (4.16)–(4.19)：$\int_0^L\sin\frac{n\pi x}{L}dx$ 偶數 $n$ 是 0、奇數 $n$ 是 $\frac{2L}{n\pi}$ → $2^N\frac{8L}{\pi^2}\sum_m\frac{e^{-N(2m+1)^2\pi^2/2L^2}}{(2m+1)^2}$。
   短時間 (4.20)–(4.23)：$\sum_{\text{奇}}\frac1{m^2}=\frac{\pi^2}8=1.2337$：長條累加 $1,\ 1.1111,\ 1.1511,\ 1.1715,\dots\to1.2337$（`D.bar_chart`；`BC.odd_inv_squares`）→ $2^NL$：短時間牆不重要。
   筆誤框 1（where「p. 77，(4.22)」）：錯 `\frac14\sum_{k=0}^{\infty}\frac{1}{k^2}`、對 `\frac14\sum_{k=1}^{\infty}\frac{1}{k^2}`（$k=0$ 那項是無窮大）。
3. 長時間 (4.24)：只剩第一項 $2^N\frac{8L}{\pi^2}e^{-\pi^2N/2L^2}$。對數圖（`D.plot_axes`，$N\in[0,200]$，縱軸 $\log_{10}$，$L=10$）：存活比例（`BC.box_survival(N, 10)`，$N=0,10,\dots,200$ 的點）：
   $N=50$：0.072103、$N=100$：0.005865、$N=200$：$3.881\times10^{-5}$；三條直線（`BC.box_survival_long(N, 10, kind)`）：格子最低模態（$N=100$：0.005861，穿過點）、
   更正 $\frac8{\pi^2}e^{-\pi^2N/2L^2}$（綠，$\frac8{\pi^2}=0.8106$，$N=100$：0.005830）、書上 $\frac8{\pi^2}e^{-\pi^2N/L^2}$（紅虛線，$N=100$：$4.19\times10^{-5}$，兩倍陡）。
   斜率（每步，`BC.decay_rate`）：格子 $-\ln\cos(\pi/10)=0.05018$、更正 $\pi^2/200=0.04935$、書上 $\pi^2/100=0.0987$。
   筆誤框 2（where「p. 77，(4.24)」）：錯 `e^{-\pi^2\Delta^2N/L^2}`、對 `e^{-\pi^2\Delta^2N/2L^2}`；字幕：(4.19) 的指數有 ½，(4.24) 把它掉了。
4. 揭曉：$L=20$：格子 0.01239 → 慢 4.05 倍（連續極限正好 4 倍：衰減率 $\pi^2/2L^2$）。每多走 $2L^2/\pi^2$ 步（$L=10$：20.26 步）少 $e$ 倍。
5. 熵／漏失的解釋 (4.25)–(4.27)（Fig. 4.2：中間放人、兩邊漏光的帳篷示意）：每步「分裂」×2，牆邊漏掉的比例 ∝ 斜率 ∝ $n_w/L^2$ →
   $n_w(N)=2^N(1-\frac{\alpha\kappa}{L^2})^N\approx2^Ne^{-\alpha\kappa N/L^2}$；精確計算給 $\alpha\kappa=\pi^2/2=4.935$。
收束：牆的影響要等 $N\sim L^2$ 才出現；之後活著的比例像 $e^{-\pi^2N/2L^2}$ 那樣消失。

**W404 關在球裡與立方體裡**（§4.1.4，pp. 79–82，(4.28)–(4.37)，Ex. 4.5）
1. 問：同樣的體積，關在球裡和關在立方體裡，哪一個活得久？（先猜。）
2. 球的模態（(4.28)–(4.33)）：$j_0(x)=\frac{\sin x}{x}$，牆上 $j_0(kR)=0$ → $k_nR=\pi,\ 2\pi,\ 3\pi$（`BC.j0_zero`）。圖（`D.plot_axes`，$x\in[0,10]$）：$\sin x/x$ 曲線與三個零點；
   一行字幕：角度積分只留 $l=0$（(4.30)）。左邊畫球的剖面（`Circle`）。
3. 權重 (4.34)：$\frac{8R^3}{\pi n^2}$；$R=1$：$n=1$：2.546、$n=2$：0.637（`BC.sphere_weight`＝數值積分 `BC.sphere_weight_numeric`）；短時間 (4.36)：$\frac{8R^3}{\pi}\cdot\frac{\pi^2}6=\frac43\pi R^3=4.189$——球的體積。
   筆誤框 1（where「p. 80，(4.31)」）：錯 `\frac{2}{R}\,\frac{j_0\,j_0}{1-2z+z\Delta^2k^2}`、對 `k_{n_0}^2\,\frac{j_0\,j_0}{1-6z+z\Delta^2k^2}`；字幕：照書上算 $n=1$ 的權重得 0.516，不是 2.546。
4. 長時間 (4.37)：$(6-\Delta^2k^2)^N\approx6^Ne^{-N\Delta^2k^2/6}$ → $6^N\frac{8R^3}{\pi}e^{-\pi^2N/6R^2}$。筆誤框 2（where「p. 81，(4.35)(4.37)」）：錯 `e^{-\Delta^2\pi^2N/R^2}`、
   對 `e^{-\Delta^2\pi^2n^2N/6R^2}`。格子核對（`BC.ball_rate(8)`，半徑 8 的格子球，2103 個格點）：每步衰減 0.0242；更正 $\pi^2/(6\cdot64)=0.0257$；書上 $\pi^2/64=0.154$（差 6 倍）。
5. 立方體（Ex. 4.5，`L.lattice3d("sc", …)` 當示意）：三個方向各取奇數模態 → 短時間 $6^NL^3$、長時間 $6^N\frac{512L^3}{\pi^6}e^{-\pi^2N/2L^2}$（衰減率和一維線段一樣）；
   格子 $L=6$：$N=0$：125 條、$N=2$：2970 條（`BC.cube_total`＝模態和 `BC.cube_total_modes`）。揭曉：同體積時「立方體的衰減率 ÷ 球」$=3(3/4\pi)^{2/3}=1.155$ → 球活得久
   （書外：同體積裡球的最低模態最慢，Faber–Krahn）。表（長時間前因子 ÷ 體積，`BC.PREFACTOR`）：線段 $8/\pi^2=0.811$、立方體 $512/\pi^6=0.533$、球 $6/\pi^2=0.608$。
收束：形狀只改兩件事：衰減率（最低模態）與前因子；同樣的體積，球最耐關。

**W405 穩定態＝靜電學；兩牆之間的帳篷**（§4.2–4.2.1，pp. 82–86，(4.38)–(4.49)，Ex. 4.6）
1. 問：走廊 0..5，兩端是懸崖；每一步都在 $x_0=2$ 放進一個醉漢。等很久以後，每一格平均站著幾個人？（先猜。）
2. (4.38)(4.39)：$w(\vec r)=\frac16\sum_iw(\vec\rho_i)+S(\vec r)$；泰勒展開 (4.40)：$\frac16\sum_iw(\vec r+a\hat e_i)\approx w+\frac{a^2}6\nabla^2w$ → (4.41) $-\nabla^2w=\frac6{a^2}S$＝Poisson 方程。
   數值核對（`BC.laplacian_coeff()`）：0.1667（＝1/6）。筆誤框（where「p. 83，(4.40)(4.41)」）：錯 `\frac{a^2}{3}\nabla^2 w,\ \ -\nabla^2 w=\frac{3}{a^2}S`、對 `\frac{a^2}{6}\nabla^2 w,\ \ -\nabla^2 w=\frac{6}{a^2}S`。
3. Fick 定律 (4.42)–(4.45)：$\partial_tc=-\nabla\cdot\vec j+S$、$\vec j=-\beta\nabla c$ → $-\beta\nabla^2c=S$。對照表：濃度 ↔ 電位、源 ↔ 電荷、吸收牆 ↔ 接地導體（$c=0$）、反射牆 ↔ $\hat n\cdot\nabla c=0$。
4. 帳篷 (4.46)–(4.49)：格子精確解（`BC.tent_lattice(5, 2)`）$w=0,\ \frac65,\ \frac{12}5,\ \frac85,\ \frac45,\ 0$（1.2、2.4、1.6、0.8；長條 `D.bar_chart`）＝(4.49) 取 $s_0=1$、$\beta=\frac12$
   （`BC.tent_formula`；格子上剛好精確）；斜率 $\frac65$、$-\frac45$，折角處斜率跳 $-s_0/\beta=-2$（(4.48)）。
5. Ex. 4.6（賭徒破產）：流進左牆 $w(1)/2=\frac35$、右牆 $w(4)/2=\frac25$（`BC.wall_fluxes(5, 2)`）＝從 2 出發先掉進左邊懸崖的機率 $\frac{L-x_0}{L}=\frac35$
   （`BC.absorb_left`＝解方程 `BC.absorb_left_solve`＝逐步遞迴 `BC.absorb_left_dp`）；$L\to\infty$：$\to1$（一維走不到無窮遠）。
收束：源源不絕＋亂走＝靜電：穩定態的濃度就是 Poisson 方程的解。

**W406 點源的 $1/r$ 與吸收牆的鏡像源**（§4.2.2–4.2.3，pp. 86–89，(4.50)–(4.60)，Fig. 4.4–4.5）
1. 問：在原點不斷放出走路人，穩定以後，離源 $r$ 處的濃度怎麼隨 $r$ 變？像高斯那樣急速下降嗎？
2. (4.53)–(4.56)：一次放出的高斯對時間積分：$\int_0^\infty(4\pi Dt)^{-3/2}e^{-r^2/4Dt}dt=\frac1{4\pi Dr}$（$t=Tr^2$ 把 $r$ 提出來，$K=\int_0^\infty T^{-3/2}e^{-1/aT}dT=\sqrt{\pi a}$）。
   動畫：幾條不同 $t$ 的高斯（`D.plot_axes`）疊起來變成 $1/r$；數值 $D=1$、$r=2$：0.03979＝$\frac1{8\pi}$（`BC.time_integral(2.0)`）。
3. 格子上（sc，每步在原點放 1 人）：$w(\vec r)=G(z_c;\vec r)$（第 3 章的 Green 函數，`BC.sc_green(n, 0, 0)`）。表：$n=0..5$：1.5164、0.5164、0.2573、0.1653、0.1217、0.0966；
   $n\,G$：$n=1$：0.5164、2：0.5147、3：0.4958、4：0.4869、5：0.4830、10：0.4787 → $\frac3{2\pi}=0.4775$（更正後的 (4.41)：$\frac{6}{4\pi r}$）。
   一行：$G(1)=G(0)-1$（原點的方程）。照書上的 $3/a^2$ 會得 $\frac3{4\pi}=0.2387$（差 2 倍；呼應 W405 的筆誤框，不再放框）。
4. 吸收牆旁的點源 (4.57)(4.58)：源在 $(l,0,0)$、鏡像源在 $(-l,0,0)$、牆是 $x=0$ 平面：$c\propto\frac1{|\vec r-\vec r_s|}-\frac1{|\vec r-\vec r_i|}\to\frac{2lx}{r^3}$（偶極，比 $1/r$ 快）。
   熱圖（`L.heatmap`，$z=0$ 剖面，$l=1$，`BC.conc_wall(x, y, 0, 1)`，色塊 ≤ 40×40，vmax 截斷源附近）。
   筆誤框（where「p. 87，(4.58) 前」）：錯 `x\text{-}z\ \text{plane}`、對 `y\text{-}z\ \text{plane}\ (x=0)`；字幕：源在 $x$ 軸上，牆要垂直 $x$ 軸。
5. 穿過遠處平面的總流量 (4.59)(4.60)：$\int_0^\infty\big[\frac{3x^2}{(x^2+R^2)^{5/2}}-\frac1{(x^2+R^2)^{3/2}}\big]RdR=0$（`BC.dipole_flux_integral`）；精確：每個點源穿過無限大平面的流量是 $\pm\frac12$：
   源 $+\frac12$、鏡像 $-\frac12$ → 0（任何 $X>l$，`BC.plane_flux(3, 1)`）；在牆與源之間（$X=0.5$）：$-1$——全部流進牆。對照 W401：三維走路人的 $x$ 座標本身是一維漫步，終究碰到牆。
收束：持續的點源給 $1/r$（像電荷）；一面吸收牆把它變成偶極——最後全被牆吃掉。

**W407 吸收球：細胞長不大；電容**（§4.2.4–4.2.5，pp. 89–91，(4.61)–(4.66)，Fig. 4.6–4.7）
1. 問：細胞靠擴散吃營養。細胞半徑加倍，吃進來的營養變幾倍？需要的呢？
2. (4.61)–(4.63)：$c=c_\infty(1-\frac{r_0}r)$（曲線 `D.plot_axes`，$r\in[1,8]$，$r_0=1$）；$\vec j=-\beta\nabla c=-\hat r\frac{\beta c_\infty r_0}{r^2}$；
   流量 $4\pi\beta c_\infty r_0$——任何半徑的球面都一樣（`BC.sphere_flux_numeric(r, 1)` = 12.566 = $4\pi$）。
3. 細胞的上限：吃進 $\propto r_0$、需要 $\frac43\pi r_0^3q$ → 半徑加倍：吃進 ×2、需要 ×8；$r_{\max}=\sqrt{3Dc_\infty/q}$。例（標「假設值，只看數量級」）：
   $D=2\times10^{-9}\,\mathrm{m^2/s}$、$c_\infty=0.25\,\mathrm{mol/m^3}$、$q=1\,\mathrm{mol/(m^3\,s)}$ → $r_{\max}\approx3.9\times10^{-5}$ m（約 39 µm；`BC.cell_rmax(**BC.CELL_EXAMPLE)`）。
4. 球外的點源（Fig. 4.6）：鏡像源在 $\rho=r_0^2/R$、強度 $-r_0/R$；$r_0=1$、$R=4$：鏡像在 0.25、強度 −0.25（`BC.kelvin_image(4, 1)`）；球面上濃度處處是 0；
   被球吃掉的比例＝球面流量 0.25（`BC.hit_prob_numeric`）$=r_0/R$；逃到無窮遠 $\frac{R-r_0}{R}=0.75$（高斯定律：總強度 $1-r_0/R$）。圖：圓（球）、源點、鏡像點、幾條示意路線。
5. 電容類比 (4.64)–(4.66)：$\frac{-\oint\vec j\cdot d\vec S}{c_\infty}=4\pi\beta C$；球 $C=r_0$ → 回到 $4\pi\beta r_0c_\infty$（Fig. 4.7：內導體＋遠處球殼）。
   筆誤框（where「p. 91，(4.64) 前」）：錯 `\rho=-4\pi\nabla^2\phi`、對 `\nabla^2\phi=-4\pi\rho`；字幕：(4.64) 本身是對的。
收束：擴散送來的營養跟半徑成正比、需要的跟體積成正比——只靠擴散的細胞長不大。

**W408 布滿受器的球（Berg–Purcell）＋吸收邊界＋第 4 章回顧**（§4.2.6、§4.3，pp. 91–94，(4:S-1)–(4:S-3)）
1. 問：細胞表面只有一小部分是受器（吸收點）。吃營養的速率會打幾折？
2. 書上的論證（p. 92）：$n$ 個半徑 $a$ 的受器，每個帶 $Q/n$：電位 $\approx\frac{Q}{na}$（自己）$+\frac{Q}{r_0}$（其他受器，像均勻帶電球）；$na\gg r_0$ 時第二項主導 → 和整顆球一樣好。
   筆誤框（where「p. 92–93」）：錯 `na\ll r_0,\ \ a\ll r/\sqrt{n}`、對 `na^2\ll r_0^2,\ \ a\ll r_0/\sqrt{n}`；字幕：受器只占小面積是 $na^2\ll r_0^2$，和 $na\gg r_0$ 不矛盾。
3. 書外補充（標「Berg & Purcell 1977」）：$\frac{J}{J_{\max}}=\frac{na}{na+\pi r_0}$（電容串聯 $\frac1C=\frac1{r_0}+\frac{\pi}{na}$；`BC.berg_purcell`＝`BC.berg_purcell_series`）；半飽和 $n=\pi r_0/a$。
   例 $r_0=1\,\mu$m、$a=1$ nm（`BC.BP_EXAMPLE`）：$n=3142$，只占表面 0.079%。表：$n=100,\ 1000,\ 3142,\ 10^4,\ 10^5$：0.031、0.241、0.500、0.761、0.970（曲線＋點，橫軸 $\log_{10}n$）。
4. §4.3 吸收邊界（(4:S-1)–(4:S-3)）：格點 $n\ge0$、懸崖在 $-1$：在 $y=0$ 少了左鄰居 → 加一個虛擬格點 $G(z;x,-\Delta)=0$，整條遞迴就統一了；鏡像點在 $-x-2$。
   從 0 出發、$2k$ 步回到 0 的走法數 $C(2k;0)-C(2k;-2\to0)$＝Catalan 數 1, 1, 2, 5, 14, 42（$k=0..5$；`BC.half_counts`＝`BC.half_counts_image`＝`BC.catalan`）。
   筆誤框（where「p. 93，(4:S-2)」）：錯 `G(z;0,y)=zG(z;x,\Delta)+\delta_{x,\Delta}`、對 `G(z;x,0)=zG(z;x,\Delta)+\delta_{x,0}`。
5. 第 4 章回顧（表三列，`D.tex_table`＋中文標籤）：吸收牆 → 鏡像（W401、W406）；有限區域 → 正弦／$j_0$ 模態、$e^{-\pi^2N/2L^2}$（W402–W404）；
   穩定態 → Poisson：帳篷、$1/r$、$4\pi Dr_0c_\infty$（W405–W407）。底下一行：吸收牆＝接地導體。
收束：吸收牆就是接地導體：鏡像、模態、電容——靜電學的工具全都能拿來數走法。

## 7. 製作規則（代理人必讀）

- ManimGL 1.7.2（`from manimlib import *`），**不是** Manim Community（用 `ShowCreation` 不用 `Create`）。繼承 `ChuScene`；`construct` 只呼叫**模組層級函式**
  （連續播放會借用 `construct`，寫成方法會 AttributeError）。開頭 `self.show_header("…")`＋`check_math()`，結尾收束字幕後 `self.wait(1.0)`。
- 場景骨架：
  ```python
  import sys
  from pathlib import Path
  sys.path.insert(0, str(Path(__file__).resolve().parent.parent))      # codes/
  sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
  from manimlib import *
  from combkit.draw.fulton_draw import assert_apart, assert_in_frame
  from combkit.draw.scene_base import ChuScene
  from combkit.draw.styles import mixed, zh
  from utils import walk_bc as BC
  from utils import walk_rec as R
  from utils import walk_gf as GF
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_bc`（必要時 `walk_rec`／`walk_gf`）算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串）；
  版面用 `assert_in_frame`／`assert_apart`。不要在場景裡呼叫 `BC.check()`（整套自我測試由 `check_env.py` 跑）。
  最慢的呼叫：`BC.sc_green` 第一次約 0.15 秒（之後有快取）、`BC.ball_rate(8)` 約 0.15 秒、`BC.survival_3d_dp(1000, 3)` 約 0.3 秒；其餘都很快。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W401 無、W402 1 個、W403 2 個、W404 2 個、W405 1 個、W406 1 個、W407 1 個、W408 2 個。
  框裡的 Tex 不放中文（英文字用 `\text{}`）。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）；字幕裡的「÷」寫成「除以」；µm 寫成 `\mu\mathrm{m}`。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「1/r」「e 的負 N 次方」或省略）。
- 對數圖：`D.plot_axes` 的縱軸直接畫 $\log_{10}$ 值，刻度標籤用 Tex 寫 `10^{-2}` 之類；橫軸 $\log_{10}n$ 同理。
- 大量物件：熱圖／色塊最多 40×40、數字格最多 9×9、長條最多 60 根、小走法圖最多 16 張、點最多 60 個；避免數千個小物件。不做 3D 動畫（球、立方體用 2D 剖面、`Circle`、`L.lattice3d` 斜投影）。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式（例如鏡射路線、等位線、圓形吸收體、Catalan 小路線），回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W401 | `w401_image_walker.py` | `W401ImageWalker` | A |
| W402 | `w402_two_walls.py` | `W402TwoWalls` | A |
| W403 | `w403_long_short_time.py` | `W403LongShortTime` | A |
| W404 | `w404_sphere_cube.py` | `W404SphereCube` | A |
| W405 | `w405_steady_tent.py` | `W405SteadyTent` | B |
| W406 | `w406_point_source.py` | `W406PointSource` | B |
| W407 | `w407_absorbing_sphere.py` | `W407AbsorbingSphere` | B |
| W408 | `w408_receptors.py` | `W408Receptors` | B |

連續播放：`main.py` 的 `Ch4Demo`（W401–W408）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch4/part4}`；內文 `notes/ch4/part4.tex`，圖表與數字檢查 `notes/ch4/make_figs.py`
（`check_numbers_ch4()` 逐一 assert，所有數字都從 `codes/utils/walk_bc.py` 算；表格與 TikZ 圖產生成 `notes/ch4/figs.tex`；巨集都以 `\Four` 開頭、`\label` 用 `four-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–3 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch04_formulas.md`、`walk_bc.py` 的 docstring、
`notes/ch3/part3.tex` 與 `notes/ch3/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 8–10 頁、圖表 4–5 張：(1) W401 的存活機率表（一維與 sc，精確 vs $x_0\sqrt{2/\pi N}$、$\sqrt{6x_0^2/\pi N}$；加一列一維 (4.8)：$N=10$：2.084 vs 2.523、$N=100$：7.499 vs 7.979）；
(2) $L=10$ 存活比例的對數圖（點＋格子／更正／書上三條直線）與 $L=10,20$ 的衰減率表（0.05018、0.04935、0.0987；0.01239、0.01234、0.0247）；
(3) 線段／立方體／球的表（前因子 0.811、0.533、0.608；衰減率；格子球 0.0242 vs 0.0257 vs 0.154；同體積比 1.155）；
(4) 帳篷圖（$L=5$、$x_0=2$：$\frac65,\frac{12}5,\frac85,\frac45$，標出流進兩牆的 $\frac35$、$\frac25$）；(5) sc Green 函數 $n\,G$ 的表＋Berg–Purcell 表。
小節對應書上 4.1.1–4.3，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch04_formulas.md` 的筆誤表，16 列）；
書外內容（Berg–Purcell、Faber–Krahn、圓盤電容 $2a/\pi$）標明出處。最後練習 4.1–4.6 的解答（答案見 `ch04_formulas.md`：4.1 鏡像 $G(z;x_1-x_2)-G(z;x_1+x_2)$、
4.2 取 $[z^N]$ 再積分得 (4.3)、4.3 鏡像點 $2kL\pm x_0$、4.4 $\sin q(x+\Delta)+\sin q(x-\Delta)=2\cos q\Delta\sin qx$、4.5 立方體、4.6 $(L-x_0)/L$）。
注意（第 3 章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；TikZ 從一點畫回同一點不用 `to[looseness]`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc。
