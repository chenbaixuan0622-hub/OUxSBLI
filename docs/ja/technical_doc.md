# OUxSBLI 技術ドキュメント

## 1. プロジェクト概要

**OUxSBLI** は、**CUDA Fortran** で記述された GPU 加速型計算流体力学（CFD）ソルバーです。構造化直線格子上で、陽的高次有限差分法を用いて、圧縮性ナビエ・ストークス方程式およびオイラー方程式を解きます。

### 主な機能
- **GPU アクセラレーション**: NVIDIA GPU 上で動作する CUDA Fortran カーネル
- **高次スキーム**: 2 次、4 次、6 次、8 次の空間離散化オプション
- **複数のスキーム**: KEEP (運動エネルギー保存型)、SLAU (低散逸 AUSM)、Roe
- **並列スケーラビリティ**：GPUを意識した通信機能を備えたMPIドメイン分解
- **LESサポート**：選択的混合スケールモデルを用いた大渦シミュレーション（LES）
- **参照状態制御**：衝撃波・乱流相互作用（SBLI）シミュレーション用の特別な再スケーリングモジュール

### 依存関係
- HPC SDK（バージョン24.*または25.*）
- CUDAツールキット（HPC SDKと互換性あり）
- MPIライブラリ
- ParaView（VTK出力可視化）

---

## 2. コードアーキテクチャ


### ディレクトリ構造

~~~
OUxSBLI/
├── 3D_solver/          # メインの3Dソルバー（最も完成度が高い）
│   ├── ETGV/           # オイラー・テイラー・グリーン渦
│   ├── IVST/           # 等エントロピー渦
│   ├── KHI/            # ケルビン・ヘルムホルツ不安定
│   ├── NSTGV/          # ナビエ・ストークス・テイラー・グリーン渦
│   ├── SBLI/           # 衝撃波-境界層相互作用
│   ├── TBL/            # 乱流境界層
│   └── src/            # 3D コアルーチンのソースコード
├── src/                # トップレベルのユーティリティ
│   ├── main.f90        # メインエントリポイント
│   ├── sbli.f90        # SBLI専用の代替メイン関数
│   ├── mod_constant.f90       # 物理定数
│   ├── cpu_gpu_mpi.f90        # CPU-GPU MPI抽象化
│   ├── print.f90              # 出力/印刷ユーティリティ
│   ├── calc_muscl.f90         # MUSCL再構成
│   ├── calc_physical_quantities.f90  # プリミティブ保存型変換
│   ├── set_compressible_bl.f90      # 境界層の設定
│   ├── set_coordinate.f90           # 座標変換およびヤコビアン
│   └── (その他のフラックス/粘性モジュール)    # calc_keep_kernel, calc_visc2/4, など
└── docs/
    ├── api.md          # 自動生成されたAPIリファレンス
    └── technical_doc.md # このファイル
~~~

### ケースごとの設定

各ケースフォルダ（3D_solver/ETGV、SBLIなど）には以下が含まれます：

~~~
case_name/
├── mod_globals.f90     # グローバルパラメータ（方程式タイプ、スキーム、格子サイズ、ブロック寸法）
├── set.f90             # 格子、初期条件、境界条件
├── Makefile            # ビルド設定
├── calc.sh             # 実行スクリプト
├── profile.sh          # プロファイリングスクリプト（利用可能な場合）
├── data/               # 出力ディレクトリ
└── recal/              # リスタート/チェックポイントファイル
~~~

---

## 3. 空間離散化

### 3.1 対流スキーム

#### KEEP (運動エネルギーおよびエントロピー保存型)
- **ファイル**: src/ 内の `calc_keep_kernel.f90`
- **精度**: 最大6次
- **特徴**:
  - 運動エネルギーとエントロピーを保存するために対流形式を分割
  - 滑らかな流れにおいて高精度
  - 従来のスキームよりも散逸が小さい
- **カーネル**:
  - `calc_keep_x[2,4,6]()`: X方向のフラックス
  - `calc_keep_y[2,4,6]()`: Y方向のフラックス
  - `calc_keep_z[2,4,6]()`: Z方向のフラックス

#### SLAU (Simple Low Dissipation AUSM)
- **ファイル**: src/ 内の `calc_slau_kernel.f90`
- **次数**: 衝撃波センサー付きで最大6次
- **特徴**:
  - 衝撃波から離れた領域での低散逸
  - センサー関数による衝撃波適応性
  - 圧縮性乱流においてRoe法よりも滑らか
- **カーネル**:
  - `calc_slau_x[2,4,6]()`: センサー付き X 方向
  - `calc_slau_y[2,4,6]()`: センサー付き Y 方向
  - `calc_slau_z[2,4,6]()`: センサー付き Z 方向

#### Roe スキーム
- **ファイル**: src/ 内の `calc_roe_kernel.f90`
- **特徴**:
  - 古典的な Roe フラックス分割
  - センサー付き衝撃波適応型
  - 超音速/極超音速流に有効

#### ハイブリッドスキーム
- **ファイル**: `calc_hybrid_kernel.f90`
- **特徴**: 流れの条件に基づいてスキームをブレンド

### 3.2 粘性項

#### 2次スキーム
- **ファイル**: `calc_visc2.f90`
- **カーネル**:
  - `calc_Ev2()`: X方向の粘性フラックス
  - `calc_Fv2()`: Y方向の粘性フラックス
  - `calc_Gv2()`: Z方向の粘性フラックス
  - `calc_Ev_LES2()`、`calc_Fv_LES2()`、`calc_Gv_LES2()`: LES乱流モデルを使用

#### 4次スキーム（サンドハムのラプラシアン形式）
- **ファイル**: `calc_visc4.f90`
- **特徴**: 最近追加されたもので、検証中
- **カーネル**:
  - `calc_Ev4()`, `calc_Fv4()`, `calc_Gv4()`: LESなし
  - `calc_Ev_LES4()`, `calc_Fv_LES4()`, `calc_Gv_LES4()`: LESあり

#### 粘性境界条件
- **ファイル**: `set_bc_common.f90`
- 複数の精度階数（2次、4次、6次）を持つ周期的境界条件
- LES変異の同期: `set_bc_mut_common()`

### 3.3 乱流モデル

#### LES（大渦シミュレーション）
- **ファイル**: `calc_les.f90`
- **SGSモデル**: 選択的混合スケールモデル（開発中）
- `mut`（乱流粘度）および `qc2`（フィルタリングされた運動エネルギー）フィールドを通じて追加

---

## 4. 時間離散化

### 時間積分法

`calc_time_dev.f90` に実装：

#### 明示的ルンゲ・クッタ 3 次（3-3 TVD）
```fortran
subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
```
- 3段階TVDルンゲ・クッタ法
- 高次空間スキームとの安定性確保に使用

#### 明示的ルンゲ・クッタ法 4次 (4-4)
```fortran
subroutine RungeKutta_4th(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
```
- 4段階の古典的RK4
- より高精度な時間積分

#### リスケール変種
- `RungeKutta_3rd_rescale()`: 参照状態リスケールを伴う3次
- `RungeKutta_4th_rescale()`: 基準状態の再スケーリングを伴う4次
- 衝撃波-境界層相互作用（SBLI）のケースに使用

---

## 5. CUDA Fortran GPU実装

### 5.1 GPUカーネルパターン

#### スレッド構成
- **3Dスレッドブロック**: 各カーネルは多次元ステンシル用に3Dスレッドブロックを使用する
  - E-fluxカーネル: `threadsE = (x, y, z)`
  - F-fluxカーネル: `threadsF = (x, y, z)`
  - G-fluxカーネル: `threadsG = (x, y, z)`
  - 粘性カーネル: 構造は同様だが、タイルの命名規則が異なる

- **共有メモリ**:
  - 内側のループでのグローバルメモリ読み込みを避けるため、スレッドブロックごとに読み込まれる
  - ステンシルの幅に基づいて余分な要素でパディングされる
  - 6次の場合の例: `sx = threadsE%x + 5` (各辺に3点 + 中心)

#### カーネル構造の例 (calc_keep_x6)

~~~fortran
attributes(global) subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E)
  ! 宣言: スレッドインデックス (it, jt, kt)、グローバル座標 (i, j, k)
  integer, parameter :: sx = threadsE%x + 5  ! パディング付きタイルサイズ
  real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp

  ! 共有メモリへの協調読み込み
  do ii = it-2, threadsE%x+3, blockDim%x
    i = i_base + ii
    if (within_bounds) then
      idx = ii + offset_yz
      rho(idx) = Q(1,i,j,k)     ! パディング付きで読み込み
      ...
    endif
  enddo
  call syncthreads()  ! すべてのスレッドを同期

  ! 6点ステンシルを用いた局所フラックス計算
  if (3 <= i かつ i <= nx-3) then
    E(:,i,j-1,k-1) = KEEP6(rho(idx-2:idx+3), u(...), ...)
  elseif (2 <= i かつ i <= nx-2) then
    E(:,i,j-1,k-1) = KEEP4(rho(idx-1:idx+2), ...)  ! 境界では4点法に簡略化
  else
    E(:,i,j-1,k-1) = KEEP2(rho(idx:idx+1), ...)    ! 端点では2点法
  endif
end subroutine
~~~

#### 主な最適化
- **共有メモリ**: グローバルメモリ帯域幅を 10～100 倍削減
- **アクセス結合**: メモリトランザクションの効率化を確保
- **境界フォールバック**: 低次スキームへのシームレスな降格
- **レジスタ使用**: 占有率を最大化するために慎重に管理

### 5.2 デバイスとホストのメモリ

#### 状態ベクトル Q
- **ホスト (CPU)**: `real(8), allocatable :: Q(5,nx,ny,nz)`
- **デバイス (GPU)**: カーネル内で `intent(in/out), device` を通じて暗黙的に割り当てられる
- **転送**: MPIランク間（cpu_gpu_mpiモジュールを参照）

#### 中間フィールド
- **E, F, G フラックス**: デバイス配列、カーネル内で計算
- **T (温度)**: デバイス配列、Qから計算
- **mu (動粘度)**: デバイス配列
- **mut, qc2**: LES乱流用のデバイス配列

#### 永続的なデバイスデータ
- 座標微分: デバイス上の `dx, dy, dz`
- メトリック: `xix, etay, zetaz` (座標変換)
- ヤコビアン: デバイス上の `Jacobian(nx,ny)`

---

## 6. MPI並列化

### 6.1 ドメイン分解

- **分解戦略**: 3次元空間におけるz方向に沿った1次元ドメイン分割
- **ランク割り当て**:
  - GPUあたり2つのMPIランク: `mygpu = myrank / 2`
  - 偶数ランク (0, 2, 4, ...): GPU上で計算
  - 奇数ランク (1, 3, 5, ...): 出力のためのデータ収集

### 6.2 CPU-GPU MPI通信 (cpu_gpu_mpiモジュール)

GPUメモリを用いたMPIの抽象化を提供します：

~~~fortran
interface CPUGPU_MPI_SEND
  module procedure CPU_MPI_SEND, GPU_MPI_SEND
end interface

interface CPUGPU_MPI_ISEND
  module procedure CPU_MPI_ISEND, GPU_MPI_ISEND
end interface
~~~

#### CPU_MPI_* (固定ホストメモリ)
- デバイス → ホストメモリへのコピー
- MPI操作を実行
- 追加のメモリコピーが必要だが、すべてのシステムで動作

#### GPU_MPI_* (CUDA対応MPI)
- GPUメモリからの直接MPI
- CUDA対応MPIライブラリ（NVIDIA OpenMPI）が必要
- 低レイテンシ、高帯域幅

#### ノンブロッキング操作
- `CPU_MPI_ISEND/IRECV`: 非同期操作
- `GPU_MPI_ISEND/IRECV`: GPUダイレクト非同期
- 通信と計算のオーバーラップを可能にする

### 6.3 同期

境界データ交換は各時間ステップで行われる：
1. 内部フラックスの計算
2. ゴースト領域の非同期送信/受信
3. 通信完了を待機
4. 境界状態を更新

---

## 7. 状態の表現

### 保存変数ベクトル Q

```
Q(1, i, j, k) = ρ      (密度)
Q(2, i, j, k) = ρu     (x方向の運動量)
Q(3, i, j, k) = ρv     (y方向の運動量)
Q(4, i, j, k) = ρw     (z方向の運動量)
Q(5, i, j, k) = E      (全エネルギー) または p (圧力、構成に依存)
```

### フラックス出力

**E-フラックス** (x方向):
- 形状: ステンシル計算後の `(5, nx-1, ny-2, nz-2)`
- インデックスの範囲: x [2, nx-2], y [2, ny-1], z [2, nz-1]

**F-フラックス** (y方向):
- 形状: ステンシル計算後の `(5, nx-2, ny-1, nz-2)`
- インデックス範囲: x [2, nx-1], y [2, ny-2], z [2, nz-1]

**G-flux** (z方向):
- 形状: ステンシル計算後の `(5, nx-2, ny-1, nz-2)`
- インデックス範囲: x [2, nx-1], y [2, ny-1], z [2, nz-2]

### 派生フィールド

**T (温度)**: `real(8), device :: T(nx,ny,nz)`
- 状態方程式を用いて Q から計算

**mu (動粘度)**: `real(8), device :: mu(nx,ny,nz)`
- T を用いてサザーランドの法則または類似の式から計算される

**mut, qc2 (LES)**: `real(8), device :: mut(nx,ny,nz), qc2(nx,ny,nz)`
- 乱流粘度およびフィルタリングされた運動エネルギー
- LES 粘性項の計算に使用される

---

## 8. コアモジュールリファレンス

### 8.1 mod_globals (パラメータモジュール)

**典型的な内容** (各ケースの mod_globals.f90 内):

```fortran
implicit none
! 格子パラメータ
integer, parameter :: nx = 256, ny = 128, nz = 256
real(8), parameter :: Lx = 2.d0*pi, Ly = 1.d0*pi, Lz = 2.d0*pi

! 方程式タイプ: 0=オイラー, 1=NS
integer(kind=2), parameter :: id_visc = 1

! スキームの選択
real(2), parameter         :: id_scheme   = 1      ! 1=KEEP, 2=SLAU, 3=Roe
integer(kind=8), parameter :: id_accuracy = 8      ! 2,4,6,8 (ポイント)
integer(kind=2), parameter :: id_tvd      = 2      ! 時間ステップの種類
integer(kind=4), parameter :: id_slau     = 0      ! SLAU バリアント

! GPU スレッド/ブロック構成
type(dim3) :: threads, threadsE, threadsF, threadsG
type(dim3) :: blocks, blocksE, blocksF, blocksG
! 粘性流体の場合も同様: threadsEv, threadsFv, threadsGv など

! リスケーリング (SBLI の場合)
integer(kind=2), parameter :: id_rescale = 0      ! 再スケーリング間隔、または 0

! 再起動フラグ
integer(kind=2), parameter :: id_recal = 2         ! 2=新規、4=再起動
```

### 8.2 set.f90 (問題設定)

**典型的な実装**:

```fortran
module set
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    ! 座標配列と格子間隔を作成
  end subroutine

  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    ! 初期条件を設定
    ! 各格子点における Q(:,:,:,:) = IC
  end subroutine

  subroutine set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)
    ! 座標変換のヤコビアン行列を計算
  end subroutine
end module set
```

### 8.3 calc_flux_base.f90 (フラックス・ラッパー)

フラックス計算を調整する:

```fortran
subroutine calc_EFG_Euler(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, ...)
  ! 適切な対流スキーム（KEEP/SLAU/Roe）へディスパッチ
  ! オイラーフラックスを物理フラックスに変換
end subroutine

subroutine calc_EFG_visc(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, ...)
  ! 粘性フラックスを計算し、総和に加算
end subroutine

subroutine calc_EFG_LES(id_visc, nx, ny, nz, dx, dy, dz, ...)
  ! 乱流散逸を含むLES変種
end subroutine
```

### 8.4 set_bc_common.f90 (境界条件)

**機能**:
- **周期的境界条件** (周期的): `set_bc_cyclic2/4/6()`
- **精度バリエーション**: 内部コピー用の2点、4点、6点ステンシル
- **カーネルバリエーション**: CPU版およびGPU版
- **LES 変異**: 境界で `mut` と `qc2` を同期

**フロー**:
```
set_bc_cyclic2_init(Q) -> 初期条件用の CPU バージョン
set_bc_cyclic2(Q)      -> 時間ステップ実行中の GPU カーネル
```

### 8.5 calc_time_dev.f90 (時間ステップ)

**メインインターフェース**:

```fortran
subroutine RungeKutta(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
  ! 特定のRK法へディスパッチ
  if (id_RungeKutta == 2) then
    call RungeKutta_3rd(...)
  else if (id_RungeKutta == 4) then
    call RungeKutta_4th(...)
  endif
end subroutine

subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
  ! 3段階TVD RK3:
  ! Q1 = Q^n + Δt * RHS(Q^n)
  ! Q2 = (3*Q^n + Q1 + Δt*RHS(Q1)) / 4
  ! Q^{n+1} = (Q^n + 2*Q2 + 2*Δt*RHS(Q2)) / 3

  ! 各段階内:
  do step = 1, nt
    call calc_physical_quantities(Q, T, mu)          ! 温度、粘度
    call calc_EFG_Euler(..., E, F, G)                ! 対流フラックス
    call calc_EFG_visc(..., E, F, G)                 ! 粘性項の追加
    call compute_RHS(...)                            ! フラックスの有限差分微分
    call UpdateQ(...)                                ! Qの更新
    call set_bc_cyclic(Q)                            ! 境界条件
    call MPI_exchange(Q)                             ! ゴースト領域の同期
    if (mod(step, id_rescale) == 0) then
      call step_rescale(...)                         ! 必要に応じて再スケーリング
    endif
  enddo
end subroutine
```

### 8.6 set_coordinate.f90 (メトリクス)

**関数**:
```fortran
subroutine set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)
  ! 2D/3D 曲線座標グリッドの場合: J = ∂(ξ,η)/∂(x,y)
end subroutine

subroutine set_coordinate_derivative(dx, dy, dz, xix, etay, zetaz)
  ! ∂ξ/∂x = 1/dx, ∂η/∂y = 1/dy, など
end subroutine
```

### 8.7 calc_physical_quantities.f90 (物性)

**機能**:
- 保存変数と原始変数の相互変換
- 内部エネルギーから温度を計算
- 温度から動粘度を計算（サザーランドの法則）
- 音速を計算

### 8.8 cpu_gpu_mpi.f90 (通信)

**インターフェース**:
- `CPUGPU_MPI_SEND(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)`
- `CPUGPU_MPI_RECV(id_gpu_mpi, buf, count, src, tag, comm, ireq, ierr)`
- `CPUGPU_MPI_ISEND(...)`
- `CPUGPU_MPI_IRECV(...)`

**ディスパッチロジック**:
- `id_gpu_mpi == 2`: ピン留めメモリ経路 (CPU_MPI_*)
- `id_gpu_mpi == 4`: CUDA対応経路 (GPU_MPI_*)

### 8.9 calc_rescale.f90 (参照状態)

乱流を発生させつつ背景流を保持するためのSBLIシミュレーションに使用される:

**ルーチン**:
- `calc_mean(step, ireq, flag_re, nx, ny, nz, Jacobian, QJ, Qm)`: 平均値を蓄積
- `step_rescale(...)`: 再スケーリング変換を適用
- `rescale_recv_send(...)`: 平均データのMPI通信
- `set_rescale(flag_re, step, nx, ny, nz, y, Jacobian, Qm, bltre, Qre)`: Qの逆変換

### 8.10 calc_keep_kernel_internal.f90 (最適化されたKEEPカーネル)

**目的**: 共有メモリ管理を最適化した、統合型多精度KEEPスキームカーネルの実装。

**主な機能**:
- 2次、4次、6次精度をサポートする統一インターフェース
- `kind(id_accuracy)`パラメータによる精度の選択
- ステンシルパディングの自動調整: `sx = threadsE%x + 2*io + 1`
- GPUキャッシュ効率に最適化された共有メモリアクセスパターン
- KEEPフラックス計算のために `include 『calc_keep_3d.f90』` を呼び出す

**カーネル**:
```fortran
attributes(global) subroutine calc_keep_x_in(nx, ny, nz, Q, T, E)
attributes(global) subroutine calc_keep_y_in(nx, ny, nz, Q, T, F)
attributes(global) subroutine calc_keep_z_in(nx, ny, nz, Q, T, G)
```

### 8.11 calc_slau_kernel_internal.f90 (最適化されたSLAUカーネル)

**目的**: 効率的なGPU計算のために、適応型補間と衝撃検知機能を備えたSLAUスキームをリファクタリングしたもの。

**主な最適化**:
- **多精度補間**: 柔軟な精度設定のためのデバイス関数 `interp2`, `interp4`, `interp6`
- **衝撃適応型散逸**: 自動的な衝撃検出のためにウィグル検出器を統合
- **メモリ効率**: 補間とフラックス計算を分離し、共有メモリの再利用を可能化
- **スタックフレームの回避**: 大規模なスタック割り当ての代わりに、共有メモリ配列を再利用
- **呼び出し**: SLAU1/HR-SLAU2フラックス関数については `include 『calc_slau_3d.f90』`

**カーネル**:
```fortran
attributes(global) subroutine calc_slau_x_in(nx, ny, nz, Q, sensor, E)
attributes(global) subroutine calc_slau_y_in(nx, ny, nz, Q, sensor, F)
attributes(global) subroutine calc_slau_z_in(nx, ny, nz, Q, sensor, G)
```

**補間ディスパッチ**:
```fortran
interface interp
  module procedure interp2, interp4, interp6  ! id_accuracy に基づいて選択
end interface interp
```

---

## 12. メイン実行フロー

### 12.1 プログラムのエントリポイント (src/main.f90)

```fortran
program main
  ! 1. MPIの初期化
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  mygpu = myrank / 2

  ! 2. GPUブロック/スレッドの設定（偶数ランクのみ）
  if (mod(myrank,2) == 0) then
    call set_block3(nx, ny, nz, threads, threadsE, ..., blocks, blocksE, ...)
  endif

  ! 3. メモリの割り当て
  allocate(Q(5,nx,ny,nz), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))

  ! 4. グリッドの設定
  call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
  call set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)

  ! 5. 初期条件（または再開）
  if (id_recal == 4) then
    ! チェックポイントファイルから読み込み
    open(10, file="recal/Q" // formatted_rank // 「.dat」, ...)
    read(10) header, Q
  elseif (id_recal == 2) then
    call set_init(myrank, nx, ny, nz, x, y, z, Q)
  endif

  ! 6. 時間ステップ
  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)

  ! 7. 出力の保存（偶数ランク）
  if (mod(myrank,2) == 0) then
    do l = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(:,i,j,l) = Jacobian(i,j) * Q(:,i,j,l)    ! ヤコビアンを適用
        enddo
      enddo
    enddo
    write(filename, 「(a, i5.5, a)」) 「recal/Q」, int(myrank/2+1), 「.dat」
    open(10, file=filename, ...)
    write(10) 『SEQFMT01』
    write(10) Q
    close(10)
  endif

  call MPI_FINALIZE(ierr)
end program main
```

### 12.2 SBLI バリアント (src/sbli.f90)

- デュアル領域（境界層と衝撃波の相互作用）に対応
- 領域ごとに異なる格子サイズ：`nx1/ny1/nz1` 対 `nx2/ny2/nz2`
- ランク 0-1 は領域 1 を、ランク 2 以上は領域 2 を処理
- フラックス整合のためにリスケーリングモジュールを使用

---

## 13. ビルドおよび実行ワークフロー

### 13.1 コンパイル

```bash
cd 3D_solver/ETGV    # ケースディレクトリに移動
vi mod_globals.f90   # パラメータを編集（格子、スキーム、精度、ブロック/スレッド）
vi set.f90           # 格子、初期条件、境界条件を編集
make                 # HPC SDKコンパイラでコンパイル
```

**Makefileでは通常、以下を使用します**:
- `pgfortran` または `nvfortran` (HPC SDKコンパイラ)
- CUDA Fortranコンパイルフラグ
- MPIライブラリへのリンク

### 13.2 実行

```bash
bash calc.sh                 # シミュレーションを実行
# 出力は data/ および recal/ ディレクトリに保存されます

# オプションのプロファイリング
bash profile.sh             # nsys/ncu プロファイルを生成
```

### 13.3 チェックポイントファイル

- **入力**: `recal/Q00001.dat`, `recal/Q00002.dat`, ... （ランクペアごとに1つ）
- **出力**: 同じ構造。出力間隔ごとに上書きされる
- **形式**: シーケンシャルな非フォーマットFortranまたはストリームバイナリ
- **データ**: ヤコビアンが適用されたQ配列

### 13.4 パフォーマンスの監視

- **nsys** (システムプロファイラ): GPU使用率、メモリ帯域幅
- **ncu** (カーネルプロファイラ): カーネルごとのメトリクス (占有率、帯域幅など)
- `profile.sh` を実行すると、`.nsys-rep` および `.ncu-rep` ファイルが生成されます

---

## 付録 A: クイックリファレンス

### 環境設定
```bash
module load hpc-sdk/24.x     # または適切な HPC SDK モジュール
source activate hpc_env      # conda を使用している場合
export CUDA_VISIBLE_DEVICES=0  # GPU の選択
```

### 必須コマンド
```bash
make clean && make           # 完全な再ビルド
make -j4                     # 並列コンパイル
bash calc.sh                 # スクリプト設定で実行
mpirun -np 4 ./a.out         # MPIの直接起動 (4ランク、2つのGPU)
```

### デバッグ
```bash
pgdbg ./a.out                # NVIDIAデバッガー
cuda-memcheck ./a.out        # メモリエラー検出 (古いCUDA)
compute-sanitizer ./a.out    # メモリエラー検出 (新しいCUDA)
```

### プロファイリング
```bash
nsys profile ./a.out         # システムプロファイラ
ncu ./a.out                  # カーネルプロファイラ
ncu --config full ./a.out    # 完全なメトリクス
```

---

## 付録B: ファイルタイプの規約

| 拡張子 | 目的 | 例 |
|-----------|---------|---------|
| `.f90` | Fortran 90 ソース | `calc_keep_kernel.f90` |
| `.mod` | コンパイル済みモジュール (生成済み) | `calc_keep_kernel.mod` |
| `.dat` | バイナリチェックポイントデータ | `recal/Q00001.dat` |
| `.sh` | Bashスクリプト | `calc.sh`, `profile.sh` |
| `.md` | Markdownドキュメント | `api.md`, `technical_doc.md` |
| `Makefile` | ビルド自動化 | コンパイル設定 |

---

## 付録 C: 用語集

**対流**: オイラー方程式／ナヴィエ・ストークス方程式の双曲部（流束発散）
**KEEP**: 運動エネルギー・エントロピー保存スキーム
**SLAU**: 単純低散逸AUSM（近似リーマンソルバー）
**TVD**: 総変分減少（安定性特性）
**LES**: 大渦シミュレーション（乱流モデル）
**SGS**: サブグリッドスケール（LESにおけるモデル化された乱流）
**ヤコビアン**: 座標変換の行列式；幾何学的因子
**ステンシル**: 有限差分法で使用される近傍の格子点の集合
**共有メモリ**: スレッドブロック間で共有される高速なオンチップCUDAメモリ
**SBLI**: 衝撃波-境界層相互作用
**ETGV**: オイラー・テイラー・グリーン渦
**TBL**: 乱流境界層
**リスケーリング**: 基準平均状態を維持するためのフラックス・リスケーリング

---

**ドキュメントバージョン**: 1.2 (2026年4月5日)
**生成元**: OUxSBLI ソースコード解析
**関連ソースドキュメント**: `api.md`, `README.md`, `CODE_ANALYSIS.md`
**最終更新**: 最適化されたカーネル実装のドキュメント (calc_keep_kernel_internal, calc_slau_kernel_internal)

