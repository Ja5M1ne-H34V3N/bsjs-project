# TSP 串行 vs MPI 岛模型：30 次重复实验 + Welch t 检验

## 实验设置

- 数据集: `pcb442/pcb442.tsp`（pcb442）
- 串行程序: `tsp.c` 编译为 `task1`
- MPI 程序: `mpitsp.c` 编译为 `task2`
- 共同设置: `N_COLONY=200`
- 迭代代数: 默认 `maxGen=200000`（两者均跑满）
- 重复次数: 串行 30 次；MPI 30 次
- 指标: 最终 best distance（串行取最后一条 `Gen:best`；MPI 取 `Final best distance`）

## MPI 运行参数

- `mpirun -np 4 ./task2 2000 2 chain best`

## 统计结果

| 项目 | 串行 | MPI |
|---|---:|---:|
| 平均值 mean | 51568.433 | 51327.700 |
| 方差 var (sample) | 16309.564 | 13794.769 |

## Welch 双样本 t 检验（双侧）

- t = 7.599453
- df = 57.598
- p-value = 2.99944e-10

判定（常用阈值 α=0.05）：

- **存在显著差异**（p < 0.05）

## 原始数据

- `results_raw.txt` 保存了 30 个串行结果与 30 个 MPI 结果（按顺序）。
