# 第 23 版證書重播紀錄

於 2026-11-23 UTC（臺北 11-24）從指定原檔解壓至新暫存目錄，只在暫存副本執行。ZIP 長度為 34,246,542 bytes，含 5,891 筆成員。逐一檢查 CRC 和路徑。SHA-256 為

```text
eb9731587a519742de65d3d8a951c005f435680a0444dbdbff9f1ff1c46eb196
```

環境為 macOS 26.5.1 / arm64、CPython 3.12.14、NumPy 2.3.5、Node v24.22.3、Apple clang 21.0.0。重播需要 Python 3.11 或更新版本、NumPy、Node.js 及 C++17 編譯器；附加的三支 Python 程式僅用標準函式庫。請讓 python3、node、c++ 在 PATH，啟用 Python 斷言。若執行未定義行為檢查，編譯器還須支援 -fsanitize=undefined。範例假設解壓後的 EL18_58_search 與這份說明所在的 audit 放在同一個目錄；從該目錄執行以下命令。

```sh
audit_dir=$(cd audit && pwd)
cd EL18_58_search
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 -B -u joint23/replay.py --full --manifest --node "$(command -v node)"
python3 -B joint23/verification/verify_host_list.py joint23
python3 -B joint23/verification/integer-width/replay.py .
python3 -B joint23/verification/exceptional-arithmetic/check_exceptional_arithmetic.py joint23
```

四項均以 0 退出。完整重播紀錄為 replay.log，計 36 行，含清單內 585 份檔案。另三項分別見 host-list.log、integer-width.log、exceptional-arithmetic.log。獨立枚舉有根宿主型得到 23，扣 A₁²⁴ 後確為 18 種普通宿主與 4 種特殊宿主。十八種鏈證書最大的有向數上界為 112。秩 17 有 16 筆平方超過有號 64 位，運算型別重新檢得為 Python 任意精度整數。四個 C++ 計算在未定義行為偵測下均無診斷，輸出比對含 24 份主法向輸入、15 份第二算法輸入、兩種奇素數計算。15 份和 24 份中的相同資料不可加總為新排除。

## 從原始鏈另核宿主清單

```sh
python3 -B "$audit_dir/host_lists.py" joint23
```

程式由 Coxeter 數回推各 ADE 分量並生成秩 24 的分拆。對三批原始 gzip 證書逐條檢查其索引確實分成覆蓋所有點的互不相交路徑，從路徑數算回每型最大上界，再檢查四個特殊型與 A₁²⁴ 為精確補集。逐型輸出見 host-chains.out。鏈比較的根座標條件與標記完整性由原始各 verify.py 的完整重播另行驗證。

## 來自 Gram 的另一種下界檢查

```sh
python3 -B "$audit_dir/gram_polynomial.py" joint23/a2/nonad57.json
```

程式僅取輸入的 gram 欄。逐項驗對稱性、對角 5 及非對角 ±1，用 57 個整數冪跡、Newton 恆等式和整多項式乘法得到

```text
det(tI+G) = t^39 (t+14)^8 (t+18)^8 (t^2+29t+192).
```

因特徵值皆為實數，正特徵值有 18 個、零特徵值有 39 個，且無重複軸。完整係數見 gram-polynomial.out。此項獨立於原檔所附的整數座標分解和有理 Schur 消去。

## 法向篩的正反控制

```sh
sieve_build=$(mktemp -d)
c++ -O3 -std=c++17 joint23/normal_sieve/sieve.cpp -o "$sieve_build/sieve"
python3 -B "$audit_dir/normal_controls.py" "$sieve_build/sieve"
```

控制輸入是抽象標籤圖。每筆原始二元列空間維數為 d。選定 58 個落在前 18 維的不同標籤，單胞內先取 58 團；多出的 d−18 個標籤各取獨立新基方向並設為孤立點。負控制刪除前兩點的相容邊，容量恰降到 57。對 d=18,19,20,21 各試兩個控制，全部與解析答案一致，在 d=20,21 的兩份肯定結果另逐字驗回所報法向的獨立性、正交性及容量。結果在 normal-controls.out。另造 4 筆包含額外商空間標籤的負控制。對核心標籤的 262,143 個非零二元方程逐一計數，每個至少刪去 11 點，超過最多 7 個外點。因此 58 點落入 18 空間的唯一可能是核心空間；其中一組僅有 57 點，另一組有 58 點但恰缺一條相容邊。最深的控制實際經 6 個雙法向測試（全通過），再由 3 個三法向測試拒絕。總計 12 筆皆與已知答案相符。

審讀輸出處理時，主 C++ 在輸出可行見證和排除時均可返回 0。封存重播對每份資料另逐行比對既存輸出，並檢查排除字串；前階段、見證和窮盡三種終止互斥。另一路法向驗證器對可行必要條件返回非零。A₅⁴D₄ 的有節點預算團搜尋只在完成時把團上界收緊；補圖搜尋提前停止時，已找到的反團僅用作頂點傳遞容量不等式的分母下界。這些是對特定終止條件的檢查。重播與控制的效力取決於正文中由任意實配置推到允許輸入的數學約化。
