#!/usr/bin/env python3
import re
import subprocess
import math
from typing import List, Tuple

# =============================
# Config (adjust here if needed)
# =============================
N_RUNS = 30

SERIAL_EXE = "./task1"
MPI_EXE = "./task2"

# MPI run config (mpitsp.c requires 4 args)
MPI_NP = 4
MIG_INTERVAL = 2000
MIG_COUNT = 2
MIG_TOPO = "chain"     # chain | island
MIG_STRATEGY = "best"  # best | worst
# NOTE: do NOT pass maxGen -> use default 200000 to match requirement

DATASET_PATH = "pcb442/pcb442.tsp"

# =============================
# Helpers: t distribution CDF (pure python)
# Using regularized incomplete beta (Numerical Recipes style)
# =============================

def _betacf(a: float, b: float, x: float) -> float:
    # Continued fraction for incomplete beta
    # Reference: Numerical Recipes in C
    MAXIT = 200
    EPS = 3e-14
    FPMIN = 1e-300

    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d

    for m in range(1, MAXIT + 1):
        m2 = 2 * m

        # even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c

        # odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delh = d * c
        h *= delh

        if abs(delh - 1.0) < EPS:
            break

    return h

def betai(a: float, b: float, x: float) -> float:
    # Regularized incomplete beta I_x(a,b)
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0

    ln_bt = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    bt = math.exp(ln_bt)

    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    else:
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b

def student_t_cdf(t: float, df: float) -> float:
    # CDF for Student's t distribution with df degrees of freedom
    if df <= 0:
        raise ValueError("df must be > 0")
    if t == 0:
        return 0.5

    x = df / (df + t * t)
    a = df / 2.0
    b = 0.5
    ib = betai(a, b, x)

    if t > 0:
        return 1.0 - 0.5 * ib
    else:
        return 0.5 * ib

def welch_ttest(x: List[float], y: List[float]) -> Tuple[float, float, float]:
    # Returns (t_stat, df, p_value_two_sided)
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        raise ValueError("Need at least 2 samples per group")

    mean1 = sum(x) / n1
    mean2 = sum(y) / n2

    var1 = sum((v - mean1) ** 2 for v in x) / (n1 - 1)
    var2 = sum((v - mean2) ** 2 for v in y) / (n2 - 1)

    se = math.sqrt(var1 / n1 + var2 / n2)
    if se == 0:
        return float("nan"), float("nan"), float("nan")

    t = (mean1 - mean2) / se

    # Welch–Satterthwaite df
    num = (var1 / n1 + var2 / n2) ** 2
    den = (var1 * var1) / (n1 * n1 * (n1 - 1)) + (var2 * var2) / (n2 * n2 * (n2 - 1))
    df = num / den if den != 0 else float("inf")

    # two-sided p-value
    cdf = student_t_cdf(abs(t), df)
    p = 2.0 * (1.0 - cdf)
    return t, df, p

# =============================
# Running and parsing programs
# =============================

SERIAL_LINE_RE = re.compile(r"^(\d+):(\d+(?:\.\d+)?)\s*$")
MPI_FINAL_RE = re.compile(r"Final best distance:\s*([0-9]+(?:\.[0-9]+)?)")

def run_cmd(cmd: List[str]) -> str:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.stdout

def parse_serial_best(output: str) -> float:
    last = None
    for line in output.splitlines():
        m = SERIAL_LINE_RE.match(line.strip())
        if m:
            last = float(m.group(2))
    if last is None:
        raise RuntimeError("Cannot parse serial best distance from output")
    return last

def parse_mpi_best(output: str) -> float:
    last = None
    for line in output.splitlines():
        m = MPI_FINAL_RE.search(line)
        if m:
            last = float(m.group(1))
    if last is None:
        raise RuntimeError("Cannot parse MPI final best distance from output")
    return last

def compile_all():
    subprocess.check_call(["gcc", "-O2", "-Wall", "-Wextra", "-o", "task1", "./tsp.c", "-lm"])
    subprocess.check_call(["mpicc", "-O2", "-Wall", "-Wextra", "-o", "task2", "./mpitsp.c", "-lm"])

def mean(xs: List[float]) -> float:
    return sum(xs) / len(xs)

def var_sample(xs: List[float]) -> float:
    mu = mean(xs)
    return sum((v - mu) ** 2 for v in xs) / (len(xs) - 1)

def main():
    print("Compiling...")
    compile_all()
    print("OK")

    serial_vals: List[float] = []
    mpi_vals: List[float] = []

    print(f"\nRunning serial {N_RUNS} times (maxGen=200000)...")
    for i in range(N_RUNS):
        out = run_cmd([SERIAL_EXE])
        v = parse_serial_best(out)
        serial_vals.append(v)
        print(f"serial {i+1:02d}/{N_RUNS}: {v}")

    print(
        f"\nRunning MPI {N_RUNS} times (np={MPI_NP}, interval={MIG_INTERVAL}, count={MIG_COUNT}, "
        f"topo={MIG_TOPO}, strategy={MIG_STRATEGY}, maxGen=200000)..."
    )
    for i in range(N_RUNS):
        out = run_cmd(
            ["mpirun", "-np", str(MPI_NP), MPI_EXE, str(MIG_INTERVAL), str(MIG_COUNT), MIG_TOPO, MIG_STRATEGY]
        )
        v = parse_mpi_best(out)
        mpi_vals.append(v)
        print(f"mpi    {i+1:02d}/{N_RUNS}: {v}")

    serial_mean, mpi_mean = mean(serial_vals), mean(mpi_vals)
    serial_var, mpi_var = var_sample(serial_vals), var_sample(mpi_vals)

    t, df, p = welch_ttest(serial_vals, mpi_vals)

    with open("results_raw.txt", "w", encoding="utf-8") as f:
        f.write("# serial\n")
        for v in serial_vals:
            f.write(f"{v}\n")
        f.write("# mpi\n")
        for v in mpi_vals:
            f.write(f"{v}\n")

    with open("RESULTS.md", "w", encoding="utf-8") as f:
        f.write("# TSP 串行 vs MPI 岛模型：30 次重复实验 + Welch t 检验\n\n")
        f.write("## 实验设置\n\n")
        f.write(f"- 数据集: `{DATASET_PATH}`（pcb442）\n")
        f.write("- 串行程序: `tsp.c` 编译为 `task1`\n")
        f.write("- MPI 程序: `mpitsp.c` 编译为 `task2`\n")
        f.write("- 共同设置: `N_COLONY=200`\n")
        f.write("- 迭代代数: 默认 `maxGen=200000`（两者均跑满）\n")
        f.write(f"- 重复次数: 串行 {N_RUNS} 次；MPI {N_RUNS} 次\n")
        f.write("- 指标: 最终 best distance（串行取最后一条 `Gen:best`；MPI 取 `Final best distance`）\n\n")

        f.write("## MPI 运行参数\n\n")
        f.write(f"- `mpirun -np {MPI_NP} ./task2 {MIG_INTERVAL} {MIG_COUNT} {MIG_TOPO} {MIG_STRATEGY}`\n\n")

        f.write("## 统计结果\n\n")
        f.write("| 项目 | 串行 | MPI |\n|---|---:|---:|\n")
        f.write(f"| 平均值 mean | {serial_mean:.3f} | {mpi_mean:.3f} |\n")
        f.write(f"| 方差 var (sample) | {serial_var:.3f} | {mpi_var:.3f} |\n\n")

        f.write("## Welch 双样本 t 检验（双侧）\n\n")
        f.write(f"- t = {t:.6f}\n")
        f.write(f"- df = {df:.3f}\n")
        f.write(f"- p-value = {p:.6g}\n\n")
        f.write("判定（常用阈值 α=0.05）：\n\n")
        if p < 0.05:
            f.write("- **存在显著差异**（p < 0.05）\n")
        else:
            f.write("- **未观察到显著差异**（p ≥ 0.05）\n")
        f.write("\n## 原始数据\n\n")
        f.write("- `results_raw.txt` 保存了 30 个串行结果与 30 个 MPI 结果（按顺序）。\n")

    print("\nDone. Wrote RESULTS.md and results_raw.txt")
    print(f"Welch t-test: t={t}, df={df}, p={p}")

if __name__ == "__main__":
    main()
