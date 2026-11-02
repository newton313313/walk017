# walk017：Rudnick & Gaspari《Elements of the Random Walk》全書 12 章的 ManimGL 動畫與中文導讀

（本文的相對路徑都以本 repo 的根目錄為基準，開發機上是 `C:\mypython\manim\projects\walk017`。）

把 J. Rudnick & G. Gaspari, *Elements of the Random Walk*（Cambridge University Press 2004）逐章做成
繁體中文的 ManimGL 動畫（W101–W1206，每章另有連續播放的 `Ch1Demo`–`Ch12Demo`）與
導讀講義 `notes/walk_導讀.pdf`（140 頁，全書第 1–12 章）。

| 位置 | 內容 |
|---|---|
| `codes/` | ManimGL 場景與各章數學模組 `codes/utils/walk_*.py`；用法見 `codes/README.md`，一鍵輸出全部影片用 `codes/render_all.bat` |
| `notes/` | 導讀講義 `walk_導讀.tex`／`.pdf`、各章內文 `notes/chN/`、各章計畫 `notes/plan_chNN.md` |
| `交接說明_HANDOFF.md` | 現況、檔案說明與交接（權威說明） |

- 環境：3Blue1Brown 的 **ManimGL 1.7.2**（3b1b/manim commit fafa083），不是 Manim Community。
- 原書 PDF、全書 md 轉檔與逐式校對（本機的 `ref/`）受著作權保護，**不在這個 repo 裡**，請自行取得原書。
- 影片（`codes/renders/`）不進 git，可用 `codes/render_all.bat` 重新產生。
