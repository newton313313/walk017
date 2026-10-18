# 《Elements of the Random Walk》第 8 章 → Manim：大綱與省 token 檢討（2026-10-18）

來源：Rudnick & Gaspari 第 8 章 Properties of the random walk: introduction to scaling，pp. 193–202（PDF 頁序 209–218）。
本檔的相對路徑都相對於 `C:\mypython\manim\projects\walk017`。§1–§5 由我寫；§6 定稿場景卡、§7 製作規則、§8 分工、§9 講義
由「基礎代理人」在寫好數學模組 `codes/utils/walk_scale.py`（每個數字算出來）之後補上，格式照 `notes/plan_ch07.md` §6–§9。

## 1. 草稿大綱（4 幕；本章只有 10 頁）

| 幕 | 主題 | 書上 |
|---|---|---|
| W801 | 普適性與它的例外：步長分布的二階動差存在 → 高斯（第 2 章）；步長的尾巴太肥（$p(r)\propto r^{-\alpha}$，(8.1)）→ $\chi(\vec q)$ 在 $q=0$ 不可展開 (8.3)(8.4) → 不是高斯 (8.5)；換變數 $w=kr$ 讀出 $r$ 的冪次 (8.6)(8.7) | §8.1–8.1.1，pp. 193–195，(8.1)–(8.7)，Ex. 8.1 |
| W802 | 生成函數的標度形式：$g(z;\vec q)\approx\frac{1}{(z_c-z)+\dots q^2}$ (8.8)、$g=(z_c-z)^{-1}f\big(q(z_c-z)^{-1/2}\big)$ (8.9)、推廣成 $\nu,\eta$ (8.10)；回到實空間 (8.11)–(8.14)：$G=r^{-(d-2+\eta)}F\big(r(z_c-z)^\nu\big)$ | §8.1.2 前半，pp. 195–197，(8.8)–(8.14) |
| W803 | 最陡下降法：用 Cauchy 公式取 $z^N$ 的係數 (8.15)(8.16)、極值條件 (8.18)(8.19) → $C(N;r)\sim\exp[-(r/N^\nu)^{1/(1-\nu)}]$ (8.20)：$\nu=\frac12$ 回到高斯，$\nu\ne\frac12$ 就不是；$\gamma=(2-\eta)\nu$、總條數 $\propto x^NN^{p}$（Ex. 8.2–8.5） | §8.1.2 後半，pp. 197–199，(8.15)–(8.20) |
| W804 | 用標度重新推出 Flory：$s'=l_2s$、$\vec r'=l_1\vec r$ (8.21)(8.22)；保住動能項 → 自迴避項乘上 $l_2^{(d-4)/2}$，上臨界維度 $d=4$（Ex. 8.6–8.8）；讓兩項同比例 (8.24)(8.25) → $\nu=3/(d+2)$，和 W706 的拔河同一個答案＋第 8 章回顧 | §8.1.3，pp. 199–202，(8.21)–(8.25) |

## 2. 生動度檢討

- W801 開場對比：同樣 1000 步，普通的步伐（每步長 1）vs 長尾巴的步伐（偶爾一大跳，像 Lévy 飛行）；前者是一團霧，後者是幾團霧被長跳連起來——寬度隨 $N$ 長得比 $\sqrt N$ 快。
- W802 的「資料塌縮」：幾個不同的 $z$ 畫 $g(z;q)$ 是幾條不同的曲線；把縱軸乘 $(z_c-z)$、橫軸除 $(z_c-z)^{1/2}$，全部疊成一條——這就是標度形式。
- W803 用一維的精確二項式係數對照最陡下降法的近似（差多少一目了然）；再把指數 $1/(1-\nu)$ 列表：$\nu=\frac12$ 得 2（高斯），$\nu=\frac34$ 得 4。
- W804 呼應 W706：一次是「拔河」、一次是「放大鏡」（重新標度），兩條路得到同一個 $3/(d+2)$；上臨界維度 4 也和 W705 的冪次表互相印證。
- 找錯遊戲照舊（基礎代理人校對後列出）。

## 3. 例子檢討

| 幕 | 小例子（算到底） | 較大／書上的例子 |
|---|---|---|
| W801 | $\int e^{ikr}k^{y}dk$ 換變數 $w=kr$ 得 $r^{-(y+1)}$（具體 $y$ 的數值核對） | 長尾巴步伐的 $\chi(q)$ vs $1-c|q|^{y}$；固定種子的模擬：寬度 $\propto N^{1/y}$ 而不是 $\sqrt N$；Ex. 8.1 |
| W802 | 一維格子 $g(z;q)=1/(1-z\chi(q))$ 在三個 $z$ 的值 | 塌縮圖；實空間 $d=3$：$G\propto e^{-r/\xi}/r$、$\xi\propto(z_c-z)^{-1/2}$ |
| W803 | 一維 $N=20$（或 100）的二項式係數 vs 最陡下降法 | $\gamma=(2-\eta)\nu$、Ex. 8.4 的 $x=1/z_c$、$p=\gamma-1$；二維自迴避計數（第 7 章 `walk_saw` 的 $c_N$）的比值法估 $\mu$、$\gamma$（已知 $\gamma=43/32$ 只引用） |
| W804 | $l_2^{(d-4)/2}$ 在 $d=1..5$ 的表 | $l_1=l_2^{3/(d+2)}$ → $\nu$；Ex. 8.8 冪次型交互作用的上臨界維度 |

## 4. 程式大綱

- **數學（基礎代理人寫）**：`codes/utils/walk_scale.py`（純 Python＋numpy，scipy 可用；可 import 同資料夾的 `walk1d`、`walk_gf`、`walk_saw`）：
  長尾巴步伐的 $\chi(q)$（數值）vs 小 $q$ 的奇異展開 (8.3)(8.4)、模擬寬度的冪次（固定種子）、$C(N;r)$ 的尾巴 (8.5) vs 數值逆 Fourier、換變數技巧 (8.6)(8.7)；
  格子生成函數與標度形式 (8.8)(8.9) 的塌縮、實空間 (8.13)(8.14)；最陡下降法 (8.15)–(8.20) vs 精確係數；$\gamma=(2-\eta)\nu$ 與 Ex. 8.3–8.5；自迴避 $c_N$ 的比值法；
  重新標度 (8.22)(8.24)(8.25) 的指數、上臨界維度（Ex. 8.6–8.8）。`check()` 用兩種以上獨立算法核對全部場景卡數字，約 10 秒內。
- **重用（不改）**：`walk_draw`、`lattice_draw`、`walk_saw`（$c_N$）。
- **畫圖**：不新增模組；塌縮圖、放大鏡示意寫在場景檔的模組層級函式。
- combkit：`walk_scale.py` 加進 `tools/sync_from_walk017.py` 的 FILES；範例 `examples/ex21_scaling_collapse.py`；升 1.18.0。

## 5. 省 token 檢討

1. **10 頁 → 4 幕**；(8.21)(8.22) 的路徑積分改寫只取「兩項各乘什麼因子」，細節放講義。
2. **場景代理人只用 1 個**（4 幕，W801–W804），和講義代理人同時進行；4 幕分給兩個代理人時，各自讀規則與模組文件的固定成本會重複一次，不划算。基礎代理人照舊。
3. **只截有疑問的頁面**：基礎代理人最多截 6 頁；`ch08_formulas.md` 只收場景與講義用得到的式子＋筆誤表（我讀正文時已看到可疑處：(8.8c) 的 $z_c=1/a=1/2d$ 和第 2 章 $\chi$ 的正規化是否一致、(8.5) 的冪次、(8.20) 的指數、Ex. 8.6 的 $l_1^d/l_2^2$）。
4. **數字一次算好**：場景卡上的每個數字先在 `walk_scale.check()` assert；模擬一律固定種子、在模組裡算好；自迴避的 $c_N$ 直接用 `walk_saw` 的快取。
5. **不新增畫圖模組、不做 3D**。
6. 代理人的最後一次 render_check 必須是最終版（每幕最多 3 次）；我只看一次 `--tile 480` 總覽圖。
7. **講義第 8 節目標 5–6 頁**；練習解答 Ex. 8.1–8.8。
8. 1080p 背景一次一支、不輪詢；收尾用腳本一次改 README、交接說明、install.log；walkit 是否更新由使用者決定。

## 6. 定稿場景卡（代理人讀這裡）

數學模組 `codes/utils/walk_scale.py` 已完成（`python utils/walk_scale.py` 自我測試通過，約 1 秒；`--summary` 印出全部數字）：下面每個數字都在它的 `check()` 裡 assert 過
（多數用字串比對小數位，卡上怎麼寫就怎麼 assert）。書上的式子與筆誤見 `ref/md/ch08_formulas.md`（場景代理人不必讀；需要的式子都抄在卡上）。

共同：每幕 2–3 分鐘；開頭一個問題（字幕），結尾一句收束；字幕繁中。記號照書：步數 $N$、條數 $C(N;\vec r)$、結構函數 $\chi(\vec q)$、生成函數 $g(z;\vec q)$（傅立葉）與 $G(z;\vec r)$（實空間）、
逸度 $z$（$z_c=\frac1{2d}$）、標度指數 $\nu,\eta$（普通漫步 $\frac12$、0）、$\gamma$、$\mu$；長尾巴步伐：步長 $r\ge1$、$P(\text{步長}>r)=r^{-y}$（三維空間密度 $p\propto r^{-3-y}$，書上寫 $r^{-5+\delta}$，$y=2-\delta$；主例 $y=1$、$\delta=1$）。
§8.1.2 的 $a$ 是 $\chi(\vec 0)$（簡單立方 $a=2d$、$b=1$），不是格距——畫面上第一次出現時一行字幕說明。重新標度：$s'=l_2s$、$\vec r\,'=l_1\vec r$。
程式代號：`WS` = `utils.walk_scale`、`SW` = `utils.walk_saw`（只在 W803 用 `SW.saw_c`）、`D` = `utils.walk_draw`、`L` = `utils.lattice_draw`（見 §7）。

**W801 長尾巴的步伐：不再是高斯**（§8.1–8.1.1，pp. 193–195，(8.1)–(8.7)，Ex. 8.1）
1. 問：如果醉漢偶爾來一個超大步，走 1000 步之後，還是一團高斯的霧嗎？
2. 對比足跡（`WS.trail(1000, None)` 與 `WS.trail(1000, 1)`，種子 801；取第 0、1 欄畫投影，左右兩框各自縮放、各畫比例尺「10」與「1000」）：左＝步長 1；右＝$P(\text{步長}>r)=r^{-1}$，
   幾團霧被長跳連起來。(8.1)(8.2)：$p\propto r^{-\alpha}$ → $\int p\,r^2d^3r\propto\int r^{4-\alpha}dr$；$\alpha=4$（$y=1$）截到 $R$（`WS.moment2_trunc(1, R)`）：$R=10,100,1000$ → $\langle r^2\rangle_R=9,\ 99,\ 999$：二階動差不存在。
3. (8.3)(8.4) $\chi(\vec q)\approx A-B|\vec q|^{2-\delta}$：主例 $\delta=1$：$B=\frac\pi4=0.7854$（`WS.B_pareto(1)`）；數值（`WS.chi_ratio(q, 1)`）$\frac{1-\chi(q)}{q}=0.7687,\ 0.7837,\ 0.7852$（$q=0.1,0.01,0.001$）。
   log–log 圖（`D.plot_axes`，橫軸 $\log_{10}q\in[-3,-1]$，縱軸 $\log_{10}(1-\chi)$；點由 `WS.chi_pareto(q, y)` 算）：$\delta=1$ 斜率 1.00、$\delta=\frac12$（$y=1.5$）斜率 1.49（`WS.chi_slope`）、高斯步伐的參考線斜率 2。
4. 寬度表（`WS.median_width(N, y)`、`WS.width_pred`、`WS.width_slope`；各 2000 人，種子 801）：$|R_N|$ 的中位數，$N=10,100,1000$：
   長尾巴 18、172、1822（除以 $N$：1.77、1.72、1.82 → 理論 1.7785）；步長 1：2.86、8.98、28.06（除以 $\sqrt N$：0.90、0.90、0.89 → 理論 0.8881）。log–log 斜率 1.01 vs 0.50：寬度 $\propto N^{1/y}$，不是 $\sqrt N$。
5. (8.5) 尾巴 $C(N;\vec r)\propto|\vec r|^{-5+\delta}$（`WS.levy_density` 數值逆 Fourier ÷ `WS.levy_tail` 換變數的漸近式，$s=1$）：$\delta=1$：0.9803、0.9978、0.9998；$\delta=\frac12$：1.1755、1.0316、1.0051（$r=10,30,100$）；
   局部斜率（$r=100$，`WS.levy_slope`）$-4.00$、$-4.51$ → $-4$、$-4.5$。換變數 (8.6)(8.7)：$\int e^{ikr}k^ydk$，令 $w=kr$ → $r^{-(y+1)}\int e^{iw}w^ydw$；$y=\frac12$：$|I(2)|/|I(1)|=0.3536=2^{-3/2}$（`WS.subst_integral`）。
   Ex. 8.1 一行：$\frac{4\pi}r\int q\sin(qr)\,q^{2-\delta}dq$ → $r^{-1}\cdot r^{-(4-\delta)}=r^{-5+\delta}$。
   筆誤框（where「p. 195，(8.6)→(8.7)」）：錯 `w=k/r`、對 `w=kr`；字幕：照書上 $k=rw$，積分變成 $r^{y+1}\int e^{ir^2w}w^ydw$，$r$ 提不出來。
收束：二階動差一發散，中央極限定理就不管用：寬度長得像 $N^{1/y}$，尾巴是冪次 $r^{-5+\delta}$。

**W802 標度形式：所有曲線疊成一條**（§8.1.2 前半，pp. 195–197，(8.8)–(8.14)）
1. 問：$z$ 越靠近 $z_c$，$g(z;q)$ 越尖越高；不同 $z$ 的曲線，能不能換個座標疊成一條？
2. (8.8)：$g=\frac1{1-z\chi}$、$\chi=a-bq^2$ → $g\approx\frac1a\frac1{z_c-z+b/a^2\,q^2}$，$z_c=\frac1a=\frac1{2d}$。一維 $g=\frac1{1-2z\cos q}$（$a=2$、$b=1$，`WS.g_lattice`）：$z=0.4,0.45,0.49$ → $g(z;0)=5,\ 10,\ 50$。
   筆誤框（where「p. 196，第一段」）：錯 `z_c=\text{coordination number}=2d`、對 `z_c=1/(2d)`；字幕：照書上一維 $z_c=2$，可是 $z>\frac12$ 時 $\sum(2z)^N$ 早就發散了。
3. 塌縮（`D.plot_axes` 兩張；先畫三條 $g$ 對 $q\in[-1,1]$，再變形成 $Y=(z_c-z)g$ 對 $x=q/\sqrt{z_c-z}\in[-4,4]$；曲線由 `WS.collapse(eps, x)` 取樣）：
   表 $x=0,1,2,4$：$z_c-z=0.1$：0.5000、0.4172、0.2819、0.1317；$0.05$：0.5000、0.4085、0.2652、0.1145；$0.01$：0.5000、0.4017、0.2529、0.1027 → 極限 $\frac2{4+x^2}$（`WS.collapse_limit`）：0.5000、0.4000、0.2500、0.1000（虛線）。
   (8.9) $g=|\vec q|^{-2}f\big(q(z_c-z)^{-1/2}\big)$ → 推廣 (8.10) $g=|\vec q|^{-2+\eta}f\big(q(z_c-z)^{-\nu}\big)$（普通漫步 $\nu=\frac12$、$\eta=0$；自迴避漫步不同）。
4. 回實空間 (8.11)–(8.13)：令 $\vec q=\vec\kappa/|\vec r|$ → $G(z;\vec r)=|\vec r|^{-d+2-\eta}F\big(|\vec r|(z_c-z)^\nu\big)$。三維（簡單立方）：$G\approx\frac3{2\pi}\frac{e^{-r/\xi}}r$，$\xi=\frac1{\sqrt{6(1-z/z_c)}}\propto(z_c-z)^{-1/2}$
   （`WS.xi_88c`、`WS.xi_lattice`）：$z/z_c=0.96$：$\xi=2.0412$（格子精確 2.0205）；$0.99$：$4.0825$（4.0722）：$z_c-z$ 變 $\frac14$，$\xi$ 變 2 倍。格子精確 $G$ ÷ Yukawa（$r=8$，`WS.G_ratio(x, 8)`）：1.0231、1.0178。
   $rG=F(u)=\frac3{2\pi}e^{-6u}$、$u=r\sqrt{z_c-z}$，$F(0)=0.4775$（`WS.F3`）。
5. 書上對 $F$ 的兩個假設（兩行字＋小示意曲線）：$x\to0$ 趨於常數、$x\to\infty$ 指數衰減（三維的 $\frac3{2\pi}e^{-6u}$ 兩條都滿足）。無其他數字。
收束：靠近 $z_c$ 只剩一個長度 $\xi\propto(z_c-z)^{-\nu}$；用 $\xi$ 當尺，所有曲線疊成一條——這就是標度形式。

**W803 最陡下降：從標度形式讀出分布**（§8.1.2 後半，pp. 197–199，(8.15)–(8.20)，Ex. 8.2–8.5）
1. 問：知道 $G(z;\vec r)$ 的標度形式，怎麼拿回「$N$ 步走到 $\vec r$」的條數？
2. (8.15)–(8.17)：$C=\frac1{2\pi i}\oint\frac{G}{z^{N+1}}dz$ → 極值化 $A\big(|\vec r|(z_c-z)^\nu\big)-(N+1)\ln z$，$A=\ln F\approx-Bx$；(8.18) $-\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}z=0$ → (8.19) $\delta=z_c-z\approx\Big(\frac{N+1}{\nu Bz_c|\vec r|}\Big)^{1/(\nu-1)}$。
   筆誤框 1（where「p. 198，(8.18)」）：錯 `\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}{z}=0`、對 `-\nu|\vec r|(z_c-z)^{\nu-1}A'-\frac{N+1}{z}=0`；字幕：照書上 $A'=-B$ 時兩項都是負的，沒有解。
3. (8.20) $C(N;\vec r)\propto\exp\Big[-(1-\nu)\big((\nu z_c)^\nu B\big)^{\frac1{1-\nu}}(N+1)^{-\frac\nu{1-\nu}}|\vec r|^{\frac1{1-\nu}}\Big]$，即 $\exp\big[-c\,(r/N^\nu)^{1/(1-\nu)}\big]$。
   檢驗：簡單立方 $\nu=\frac12$、$B=2d$、$z_c=\frac1{2d}$ → 係數 $\frac d2$（`WS.sd_K_dim`；$d=3$：$\frac32$）→ $e^{-dr^2/2N}$＝第 2 章的高斯。
   筆誤框 2（where「p. 198，(8.20)」）：錯 `e^{-((\nu z_c)^\nu B)^{1/(1-\nu)}\cdots}`、對 `e^{-(1-\nu)((\nu z_c)^\nu B)^{1/(1-\nu)}\cdots}`；字幕：照書上一維得 $e^{-x^2/N}$，變異數少一半（第一行的 $-(N+1)/z_c$ 也應為 $+(N+1)\delta/z_c$）。
4. 一維 $N=20$ 長條圖（`D.bar_chart`，$x=-20..20$ 偶數，高 `WS.ratio_exact(20, x)`）＋兩條曲線 $e^{-x^2/42}$（更正，綠）與 $e^{-x^2/21}$（書上，紅虛線）；表（精確｜更正｜書上，`WS.sd_ratio`）：
   $N=20$：$x=4$：0.6818｜0.6832｜0.4668；$x=8$：0.2098｜0.2179｜0.0475。$N=100$：$x=8$：0.7282｜0.7285｜0.5306；$x=12$：0.4894｜0.4902｜0.2403。
   字幕一行：一維在 $z=\pm\frac12$ 各有一個鞍點，完整的條數要乘 2（取比值時消掉）。
5. 指數 $\frac1{1-\nu}$（Fisher，`WS.fisher`）：$\nu=\frac12\to2$（高斯）、$\frac35\to\frac52$、$\frac34\to4$（二維自迴避）；三維 $0.588$（只引用）→ 2.43。Ex. 8.3–8.4（`WS.gamma_exp`、`WS.growth`）：總條數 $\propto(z_c-z)^{-\gamma}$，$\gamma=(2-\eta)\nu$：
   普通 $(2-0)\cdot\frac12=1$；二維自迴避 $\eta=\frac5{24}$、$\nu=\frac34$ → $\frac{43}{32}$（已知值，只引用）；條數 $\propto x^NN^p$，$x=\frac1{z_c}$、$p=\gamma-1$：二維 $x=\mu\approx2.638$、$p=\frac{11}{32}=0.34375$。
   比值法（`WS.ratio_method(N)`，第 7 章的 $c_N$；$c_{15},c_{16}$ 是文獻值）：$N=10,12,14,16$：$\hat\mu=2.6592,\ 2.6529,\ 2.6483,\ 2.6454\to2.638$；$\hat\gamma=1.3552,\ 1.3519,\ 1.3496,\ 1.3483\to1.34375$（表或兩張點圖＋虛線）。
   Ex. 8.5 一行：小 $r$ 的條數 $\propto r^\theta$，$\theta=\frac{\gamma-1}\nu$（二維自迴避 $\frac{11}{24}$，`WS.theta`）。
收束：$\nu=\frac12$ 才是高斯；$\nu\ne\frac12$ 時 $C\sim\exp[-(r/N^\nu)^{1/(1-\nu)}]$——自迴避漫步的分布不是高斯。

**W804 放大鏡：重新標度推出 Flory＋第 8 章回顧**（§8.1.3，pp. 199–202，(8.21)–(8.25)，Ex. 8.6–8.8）
1. 問：把漫步「縮小來看」——步數變少、長度變短——自迴避的效果會變大還是變小？
2. (8.21) 指數兩項：動能 $\frac3{2a^2}\int|\frac{d\vec r}{ds}|^2ds$、自迴避 $\frac{v_0}2\iint\delta(\vec r(s_1)-\vec r(s_2))$；$s'=l_2s$、$\vec r\,'=l_1\vec r$（放大鏡示意，模組層級函式）→ (8.22)：動能 $\times\frac{l_2}{l_1^2}$、自迴避 $\times\frac{l_1^d}{l_2^2}$
   （$\delta(\vec r/l_1)=l_1^d\delta(\vec r)$、$ds_1ds_2=ds_1'ds_2'/l_2^2$；`WS.kin_factor`、`WS.int_factor`，數值核對 `WS.rescale_check`）。
   筆誤框 1（where「p. 200，(8.22) 前」）：錯 `\vec r\,'(s')=l_1\vec r(l_2s)`、對 `\vec r\,'(s')=l_1\vec r(s'/l_2)`；字幕：新曲線在 $s'$ 的位置，是舊曲線在 $s=s'/l_2$ 的位置放大 $l_1$ 倍。
3. Ex. 8.6–8.7：令 $l_1=\sqrt{l_2}$（動能項不變、$a'=a$）→ 自迴避項 $\times l_2^{(d-4)/2}$；$l_2=\frac1N$ 只剩一步 → $R\sim a\sqrt N$。表（`WS.sa_amplitude(100, d)`）$N=100$：
   $d=1,2,3,4,5$：1000、100、10、1、0.1 → 上臨界維度 4（呼應 W705 冪次表 $d=4$ 那列全是 1）。
   筆誤框 2（where「p. 200，Ex. 8.6」）：錯 `l_1^d/l_2`、對 `l_1^d/l_2^2`；字幕：照書上得 $l_2^{(d-2)/2}$，上臨界維度會變成 2。
   Ex. 8.8 一行：$|\vec r_1-\vec r_2|^{-q}$ → $l_2^{(q-4)/2}$，$q$ 取代了 $d$；重要 ⇔ $\min(d,q)<4$（`WS.ex88_relevant`）。
4. (8.24)(8.25)：讓兩項同比例 $\frac{l_2}{l_1^2}=\frac{l_1^d}{l_2^2}$ → $l_1=l_2^{3/(d+2)}$ → $R\sim l_1^{-1}=N^{3/(d+2)}$：$\nu=1,\frac34,\frac35,\frac12$（$d=1..4$，`WS.nu_rescale`＝`SW.flory_nu`）；
   例 $d=3$、$l_2=0.01$：$l_1=0.0631$、$1/l_1=15.85=100^{3/5}$（`WS.flory_l1`）。和 W706 的拔河同一個答案（左右兩欄並列：拔河 vs 放大鏡）。
5. 第 8 章回顧（表五列，`D.tex_table`＋中文標籤）：長尾巴寬度 $\propto N^{1/y}$（1822 vs 28.06）、尾巴 $r^{-5+\delta}$（W801）；塌縮 $\to\frac2{4+x^2}$、$\xi\propto(z_c-z)^{-1/2}$（W802）；
   $e^{-(r/N^\nu)^{1/(1-\nu)}}$、Fisher 2 與 4、$\gamma=(2-\eta)\nu=\frac{43}{32}$（W803）；$l_2^{(d-4)/2}$、上臨界維度 4（W804）；$\nu=\frac3{d+2}$（W804＝W706）。無筆誤框。
收束：放大鏡底下，四維以下自迴避越看越重要；讓兩項一起縮放，就得到 Flory 的 $\frac3{d+2}$。

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
  from utils import walk_scale as WS
  from utils import walk_saw as SW
  from utils import lattice_draw as L
  from utils import walk_draw as D
  ```
- 每個畫面上的數字都由 `walk_scale` 算出、在 `check_math()` 裡 assert（卡片上的數字就是 assert 的目標；照卡上的小數位用 `f"{x:.4f}"` 之類比對字串，分數用 `Fraction`；
  Fraction 不能直接 `:.4f`，要 `float()`）；版面用 `assert_in_frame`／`assert_apart`。不要在場景裡呼叫 `WS.check()` 或 `WS.summary()`（整套自我測試由 `check_env.py` 跑）。
  隨機足跡與模擬一律用模組的預設種子（801）。最慢的呼叫：`WS.median_width(1000, 1)`／`(1000, None)` 各約 0.1 秒（有快取）、`WS.G_lattice_sc` 每次約 0.01 秒、
  `WS.chi_pareto`／`WS.levy_density` 每次約 1–5 毫秒（曲線取樣 ≤ 40 點）；其餘都很快。scipy 的 IntegrationWarning 不影響結果（`check_math` 裡用 `warnings.catch_warnings()` 包起來）。
- 筆誤框一律 `D.erratum_box(where, wrong_tex, right_tex)`（卡上寫了 where 與兩段 Tex）；每幕的筆誤框：W801 1 個、W802 1 個、W803 2 個、W804 2 個。框裡的 Tex 不放中文（英文字用 `\text{}`）。
- 字幕 `self.say("中文", "$數學$", "中文")`（以 `$` 開頭的段落是 LaTeX），每句 ≤ 30 字；字幕在 $y=-3.45$，**內容放在 $y\in[-2.9,3.2]$**，左上角留給標題。
  以 `$` 開頭的段落整段當 LaTeX：「`$N$ 步都…`」要拆成 `"$N$", " 步都…"`。字幕裡的「÷」寫成「除以」。
- `Tex(..., color=X)` 無效，用 `.set_color(X)`；同一個 `play` 不要同時動畫父物件與它的子物件；Tex 開頭結尾不放 `\,`、`\;`；
  用 `t2c` 時不要 `\left(\right)`（用 `\Big(`）；`\iff` 改 `\Leftrightarrow`；Tex 裡不放中文（中文用 `zh()`／`mixed()`）。
  `get_center()` 回傳活的陣列，要存就 `.copy()`；不要 `.scale()` 元件（用參數 cell、unit、width 決定大小）。標題的 `Text` 不能有上標（寫「N 的 1/y 次方」或省略）。
- 名稱：`L`（lattice_draw）與 `D`（walk_draw）是模組——不要拿來當變數；`WS`、`SW` 也是。
- 足跡（W801）：一條足跡＝**一個** `VMobject().set_points_as_corners(...)`（1001 個點也只是一個物件），線寬 1–1.5；兩框各自縮放（長尾巴的框按最大範圍縮，步長 1 的框按自己的範圍），
  框外畫比例尺；不要每一步一個 `Line`。log–log 圖、塌縮圖、長條圖都用 `D.plot_axes`／`D.bar_chart`／`D.curve`，刻度數字自己加。
- 塌縮動畫（W802）：同一組取樣點，`Transform` 舊曲線到新座標的曲線（三條顏色固定，極限曲線虛線）；點數每條 ≤ 60。
- 放大鏡（W804）：模組層級函式（圓圈＋把手＋框內放大的一小段折線即可）；「拔河 vs 放大鏡」並列用兩個小框＋文字，不必重畫 W706。
- 大量物件：長條最多 60 根（同一畫面）、表格最多 6 列 × 10 欄、點最多 60 個；避免數千個小物件。不做 3D。
- 只建立自己的場景檔；共用檔（`utils/`、`main.py`、`render_all.py`）不改，缺什麼寫成場景檔裡的模組層級函式，回報時說明。每幕 250–450 行。
- 自查：在 `C:\mypython\manim\projects\combkit` 執行
  `python tools/render_check.py C:/mypython/manim/projects/walk017/codes/scenes/<檔名> <類別> --config C:/mypython/manim/projects/walk017/codes/custom_config.yml --tile 480`，
  看 `runs/<類別>/sheet.png`（重疊、出框、缺字、數字）；**每幕最多 3 次，且最後一次 render_check 必須是最終版**：檔案改過（哪怕只改一行）就要再跑一次，不能停在沒驗證過的版本。
  不要直接在 codes/ 跑 manimgl，也不要用 `-s`。不要 commit。
- 回報 ≤ 8 行：檔名、長度、STATUS（最後一次 render_check 的結果）、偏離本檔之處。

## 8. 檔名與分工（`codes/scenes/`）

| 代號 | 檔名 | 類別 | 代理人 |
|---|---|---|---|
| W801 | `w801_non_gaussian.py` | `W801NonGaussian` | A |
| W802 | `w802_scaling_form.py` | `W802ScalingForm` | A |
| W803 | `w803_steepest_descent.py` | `W803SteepestDescent` | A |
| W804 | `w804_rescaling_flory.py` | `W804RescalingFlory` | A |

連續播放：`main.py` 的 `Ch8Demo`（W801–W804）。

## 9. 講義（1 個代理人，與場景代理人同時進行）

`notes/walk_導讀.tex`（主檔）加 `\input{ch8/part8}`；內文 `notes/ch8/part8.tex`，圖表與數字檢查 `notes/ch8/make_figs.py`
（`check_numbers_ch8()` 逐一 assert，所有數字都從 `codes/utils/walk_scale.py` 算；表格與 TikZ／pgfplots 圖產生成 `notes/ch8/figs.tex`；巨集都以 `\Eight` 開頭、`\label` 用 `eight-…`；
執行後印出含「數字檢查通過」的一行，`check_env.py` 之後會照第 1–7 章的樣子加一項）。只讀：本檔 §6、`ref/md/ch08_formulas.md`、`walk_scale.py` 的 docstring、
`notes/ch7/part7.tex` 與 `notes/ch7/make_figs.py` 的開頭（格式範例）；不讀原文。
目標 5–6 頁、圖表 4 張：(1) 長尾巴：$1-\chi(q)$ 的 log–log 圖（pgfplots，$\delta=1$、$\frac12$ 的數值點＋斜率 2 的參考線）＋寬度表（W801 的 18、172、1822 vs 2.86、8.98、28.06，理論 1.7785、0.8881）
＋尾巴比值表（0.9803… 與 1.1755…）；(2) 塌縮圖（三條 $(z_c-z)g$ 對 $x$＋極限 $\frac2{4+x^2}$）＋W802 的表＋三維 $\xi$ 表（2.0412／2.0205、4.0825／4.0722）；
(3) 最陡下降：一維 $N=20$ 的長條＋$e^{-x^2/42}$ 與 $e^{-x^2/21}$ 兩條曲線＋W803 的表（$N=20$、100）＋Fisher 指數表＋比值法表；
(4) 重新標度：$N^{(4-d)/2}$ 表（$d=1..5$：1000、100、10、1、0.1）＋$\nu=\frac3{d+2}$ 表（並列第 7 章 Flory 與已知值 1、$\frac34$、0.588、$\frac12$）。
小節對應書上 8.1.1、8.1.2（標度形式／最陡下降兩小節）、8.1.3，每節：先看這裡框 → 直觀 → 例子（算到底）→ 書上的式子（標頁碼）→ 筆誤提醒（照 `ch08_formulas.md` 的筆誤表，7 列；
另附〔注〕：(8.2) 是 $\alpha\le5$（$\alpha=5$ 對數發散）、(8.5) 是 $r\gg N^{1/(2-\delta)}$ 的尾巴且正比於 $N$、(8.11) 省略 $(2\pi)^{-d}$、§8.1.2 的 $a$ 是 $\chi(0)$ 不是格距、
$F(x)\to$ 常數只在 $d>2$ 成立（一維 $F(u)=e^{-2u}/2u$）、一維最陡下降的奇偶性因子 2、Ex. 8.6 的「first term in the double integral」指動能項）。
最後練習 8.1–8.8 的解答：
8.1 三維 $\int d^3q\,e^{i\vec q\cdot\vec r}f(q)=\frac{4\pi}r\int_0^\infty q\sin(qr)f(q)dq$；$A^Ne^{-N(B/A)q^{2-\delta}}\approx A^N(1-N\frac BAq^{2-\delta})$，第一項只在 $\vec r=0$；第二項令 $w=qr$ → $r^{-1}r^{-(4-\delta)}=r^{-5+\delta}$，
係數 $\frac{NB}{A}A^N\frac{\Gamma(4-\delta)\sin(\pi(2-\delta)/2)}{2\pi^2}$（$\delta=1$ 即第 2 章 Cauchy 的 $\frac{N}{\pi^2r^4}$；數值比 0.9803、0.9978、0.9998）；
8.2 $\frac{d}{dz}$：$-\nu|\vec r|\delta^{\nu-1}A'-\frac{N+1}z=0$、$A'=-B$ → (8.19)；代回：$-B r\delta^\nu=-\frac{(N+1)\delta}{\nu z_c}$，加上 $+\frac{(N+1)\delta}{z_c}$ → 指數 $-(1-\nu)\big((\nu z_c)^\nu B\big)^{1/(1-\nu)}(N+1)^{-\nu/(1-\nu)}r^{1/(1-\nu)}$；
$\nu=\frac12$、簡單立方 → $-\frac{dr^2}{2N}$（書上的係數少 $1-\nu$，得 $-\frac{dr^2}N$；一維 $N=20$、$x=8$：0.0475 vs 精確 0.2098）；
8.3 $\sum_{\vec r}G(z;\vec r)=g(z;\vec 0)$；假設 (8.10) 一直用到 $q\to0$、$g(z;0)$ 在 $z<z_c$ 有限 ⇒ $f(x)\approx c\,x^{2-\eta}$（$x\to0$）⇒ $g(z;0)=c(z_c-z)^{-(2-\eta)\nu}$，$\gamma=(2-\eta)\nu$（普通 1；二維 $\frac{43}{32}$）；
8.4 (2:S3-10)：$(z_c-z)^{-\gamma}$ 的係數 $\approx\frac{N^{\gamma-1}}{\Gamma(\gamma)}z_c^{-N-\gamma}$ → $x=\frac1{z_c}$、$p=\gamma-1$（普通 $2d$、0；二維自迴避 $\mu\approx2.638$、$\frac{11}{32}$；比值法 $\hat\mu$、$\hat\gamma$ 的表）；
8.5 固定 $r$、$N\to\infty$ 要 $N^{-d\nu}$ ⇒ $F(x)=F_0+F_1x^{s}+\cdots$ 的非解析項 $s=d-\frac1\nu$（$(z_c-z)^{\nu s}$ 的係數 $\propto N^{-\nu s-1}=N^{-d\nu}$）⇒ 條數 $\propto z_c^{-N}N^{-d\nu}r^{\theta}$，$\theta=\frac{\gamma-1}\nu$
（普通 0、$C(N;0)\propto(2d)^NN^{-d/2}$；二維自迴避 $\frac{11}{24}$）；
8.6 $a'=al_1/\sqrt{l_2}$ 把動能項變回 $\frac3{2a'^2}$；自迴避項 $\frac{l_1^d}{l_2^2}$；$l_1=\sqrt{l_2}$ → $l_2^{(d-4)/2}$；
8.7 $l_2=\frac1N$：一步長 $a'=a$，$R=R'/l_1=a\sqrt N$；自迴避放大 $N^{(4-d)/2}$（$N=100$：1000、100、10、1、0.1）→ 上臨界維度 4；
8.8 $|\vec r|^{-q}$ 給 $\frac{l_1^q}{l_2^2}$ → $l_2^{(q-4)/2}$（與 $d$ 無關）：$q<d$ 時長程、$q<4$ 就在任何維度都重要；$q>d$ 時短程部分主導、回到 $d_c=4$；合起來重要 ⇔ $\min(d,q)<4$；Flory 版 $\nu=\frac3{q+2}$。
注意（前幾章的坑）：`\subsubsection` 標題裡有數學用 `\texorpdfstring`；含反斜線的 Python 用 Write 工具寫檔，不用 heredoc；
pgfplots 取樣點、長條與表的數字由 `make_figs.py` 算好寫成座標（不要在 TeX 裡算）；講義裡「÷」可以用 `$\div$`。
