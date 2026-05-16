#define _GNU_SOURCE
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define FILE_PATH  "pcb442/pcb442.tsp"
#define N_COLONY 200
#define CITY 442

/* 输出文件名：与 tsp.c 的 tsp0.txt 区分 */
#define OUTPUT_TXT "mpitsp0.txt"

/* ===== 参数默认值（可通过 argv 覆盖） =====
 * argv[1] 迁移间隔 migration_interval (代)
 * argv[2] 迁移数量 migration_count (每次迁移个体数)
 * argv[3] 迁移拓扑 topology: island | chain
 *         island: 全连接(all-to-all) 通过 MPI_Allgather
 *         chain : 环形链(ring) 只和相邻rank迁移（send to next, recv from prev）
 * argv[4] 迁移策略 strategy: best | worst （迁移发送最好/最坏个体）
 */
static long migration_interval = 2000;
static int migration_count = 2;
typedef enum { TOPO_ISLAND = 0, TOPO_CHAIN_RING = 1 } Topology;
typedef enum { MIG_BEST = 0, MIG_WORST = 1 } MigStrategy;
static Topology topo = TOPO_ISLAND;
static MigStrategy mig_strategy = MIG_BEST;

/* ===== GA 参数（沿用 tsp.c） ===== */
static int xColony = N_COLONY;
static int xCity = CITY;
static double probab1 = 0.02;
static long maxGen = 200000;

static int colony[N_COLONY * 2][CITY];
static double cityXY[CITY][2];
static double city_dis[CITY][CITY];
static double dis_p[N_COLONY * 2];
static int temp[CITY];

static long GenNum = 0;

static int position_in(int *tmp, int C);
static void invert_segment(int pos_start, int pos_end);
static void init_data_and_population(int rank);
static void select1_local();
static double compute_distance_of_path(const int *path);
static void local_improve_one(int idx);
static void migrate_if_needed(int rank, int size);
static void parse_args_or_die(int rank, int argc, char **argv);
static void usage_and_exit(int rank, const char *prog);

static void usage_and_exit(int rank, const char *prog) // 如果参数输入的不对就直接打印参数输入方法，然后结束程序
{
    if (rank == 0)
    {
        fprintf(stderr,
                "Usage: mpirun -np <P> %s <migration_interval> <migration_count> <topology> <strategy> [maxGen]\n"
                "  topology: island | chain\n"
                "    island = 全连接(all-to-all)迁移\n"
                "    chain  = 环形链(ring)相邻迁移\n"
                "  strategy: best | worst\n"
                "  maxGen  : 可选，最大迭代代数（默认 200000）\n"
                "Example : mpirun -np 4 %s 2000 4 island best 20000\n",
                prog, prog);
    }
    MPI_Abort(MPI_COMM_WORLD, 2); //异常就退出
    exit(2);
}

static void parse_args_or_die(int rank, int argc, char **argv)  //接收启动程序时的参数
{
    if (argc < 5)
        usage_and_exit(rank, argv[0]);

    migration_interval = atol(argv[1]);
    migration_count = atoi(argv[2]);
    if (migration_interval <= 0 || migration_count <= 0 || migration_count > xColony)
    {
        if (rank == 0)
            fprintf(stderr, "Invalid migration_interval/migration_count. interval>0, 1<=count<=%d\n", xColony);
        usage_and_exit(rank, argv[0]);
    } //遗传频率和数量填写发生问题时报错退出

    if (strcmp(argv[3], "island") == 0 || strcmp(argv[3], "all") == 0)
        topo = TOPO_ISLAND;
    else if (strcmp(argv[3], "chain") == 0 || strcmp(argv[3], "ring") == 0 || strcmp(argv[3], "line") == 0)
        topo = TOPO_CHAIN_RING;//遗传拓扑发生问题时报错退出
    else
    {
        if (rank == 0)
            fprintf(stderr, "Invalid topology: %s\n", argv[3]);
        usage_and_exit(rank, argv[0]);
    }

    if (strcmp(argv[4], "best") == 0)
        mig_strategy = MIG_BEST;
    else if (strcmp(argv[4], "worst") == 0)
        mig_strategy = MIG_WORST;
    else
    {
        if (rank == 0)
            fprintf(stderr, "Invalid strategy: %s\n", argv[4]);
        usage_and_exit(rank, argv[0]); // 接收遗传策略参数，并将正确的参数传到全局变量
    }

    /* 可选：maxGen 控制最大迭代次数 ，可填可不填*/
    if (argc >= 6)
    {
        long v = atol(argv[5]);
        if (v > 0)
            maxGen = v;
        else if (rank == 0)
            fprintf(stderr, "Warning: ignore invalid maxGen=%s, keep default %ld\n", argv[5], maxGen);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // 获取总进程数和当前进程编号

    parse_args_or_die(rank, argc, argv);

    /* 每个rank用不同随机种子，保证岛之间变化的多样性 */
    init_data_and_population(rank); // 根据rank初始化不同的种群起点

    /* rank0 输出到 txt（名字与 tsp.c 的 tsp0.txt 不同） */
    FILE *fpout = NULL;
    clock_t timeStart = 0;
    if (rank == 0)
    {
        fpout = fopen(OUTPUT_TXT, "a"); // 打开文件准备记录实验数据
        timeStart = clock(); //用来计算运行时间的，不重要

        fprintf(stdout,
                "MPI TSP island-model GA start: P=%d, interval=%ld, count=%d, topo=%s, strategy=%s\n",
                size, migration_interval, migration_count,
                (topo == TOPO_ISLAND ? "island(all-to-all)" : "chain(ring)"),
                (mig_strategy == MIG_BEST ? "best" : "worst"));

        if (fpout)
        {
            fprintf(fpout,
                    "# P=%d interval=%ld count=%d topo=%s strategy=%s maxGen=%ld\n",
                    size, migration_interval, migration_count,
                    (topo == TOPO_ISLAND ? "island" : "chain"),
                    (mig_strategy == MIG_BEST ? "best" : "worst"),
                    maxGen);
            fprintf(fpout, "# Gen\tTime(s)\tBestDistance\n");
            fflush(fpout);
        }
    }

    for (GenNum = 0; GenNum < maxGen; GenNum++)
    {
        /* 一代：对每个个体做一次局部改进 */
        for (int i = 0; i < xColony; i++)
            local_improve_one(i); // 执行局部搜索优化个体路径

        /* 选择：保留改进后的更优解 */
        select1_local(); // 比较新老个体，保留距离更短的解

        /* 迁移（岛模型） */
        migrate_if_needed(rank, size); // 到达迁移代数时，在进程间交换优良个体

        /* 打印全局最优（每隔若干代） */
        if ((GenNum % 2000) == 0)
        {
            double local_best = dis_p[0];
            for (int i = 1; i < xColony; i++)
                if (dis_p[i] < local_best)
                    local_best = dis_p[i];

            double global_best = 0.0;
            MPI_Allreduce(&local_best, &global_best, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); // 汇总所有进程的最优值

            if (rank == 0)
            {
                printf("%ld:%f\n", GenNum, global_best); // 只有主进程负责在屏幕上打印进度

                if (fpout)
                {
                    clock_t timeNow = clock();
                    fprintf(fpout, "%ld\t%4.2f\t%.0f\n",
                            GenNum,
                            (double)(timeNow - timeStart) / CLOCKS_PER_SEC,
                            global_best);
                    fflush(fpout);
                }
            }
        }
    }

    //最终输出全局最优
    double local_best = dis_p[0];
    for (int i = 1; i < xColony; i++)
        if (dis_p[i] < local_best)
            local_best = dis_p[i];
    double global_best = 0.0;
    MPI_Allreduce(&local_best, &global_best, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); // 结束前最后一次全进程同步
    if (rank == 0)
    {
        printf("Final best distance: %f\n", global_best);
        if (fpout)
        {
            clock_t timeNow = clock();
            fprintf(fpout, "%ld\t%4.2f\t%.0f\n",
                    GenNum,
                    (double)(timeNow - timeStart) / CLOCKS_PER_SEC,
                    global_best);
            fclose(fpout);
            fpout = NULL;
        }
    }// 主进程输出最好的成绩

    MPI_Finalize();
    return 0;
}

static void init_data_and_population(int rank)
//初始化种群和个体的信息
{
    FILE *fp = fopen(FILE_PATH, "r");
    if (xCity > CITY)
    {
        fprintf(stderr, "Rank %d: xCity(%d) > CITY(%d) compile-time limit\n", rank, xCity, CITY);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }

    for (int i = 0; i < xCity; i++)
    {
        double x, y;
        if (fscanf(fp, "%*d%lf%lf", &x, &y) != 2) // 读取文件的城市ID和坐标
        {
            fprintf(stderr, "Rank %d: failed to read city coord at line %d\n", rank, i);
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1);
        }
        cityXY[i][0] = x;
        cityXY[i][1] = y;
    }
    fclose(fp); // 坐标读入完成，关闭文件

    for (int i = 0; i < xCity; i++)
    {
        for (int j = 0; j < xCity; j++)
        {
            if (j > i)
            {
                double d = (cityXY[i][0] - cityXY[j][0]) * (cityXY[i][0] - cityXY[j][0]) +
                           (cityXY[i][1] - cityXY[j][1]) * (cityXY[i][1] - cityXY[j][1]);
                city_dis[i][j] = (int)(sqrt(d) + 0.5); // 计算城市间距离并四舍五入
            }
            else if (j == i)
            {
                city_dis[i][j] = 0;
            }
            else
            {
                city_dis[i][j] = city_dis[j][i];
            }
        }
    }

    // 设置rank相关随机种子
    srand((unsigned)time(NULL) ^ (unsigned)(rank * 0x9e3779b1));

    int array[CITY];
    for (int i = 0; i < xCity; i++)
        array[i] = i;

    /* 初始化种群：随机排列 */
    for (int i = 0; i < xColony; i++)
    {
        int mod = xCity;
        for (int j = 0; j < xCity; j++)
        {
            int sign = rand() % mod;
            colony[i][j] = array[sign];
            int t = array[mod - 1];
            array[mod - 1] = array[sign];
            array[sign] = t;
            mod--;
            if (mod == 1)
            {
                colony[i][++j] = array[0];
                break;
            }
        }
        dis_p[i] = compute_distance_of_path(colony[i]);
    }
}

static double compute_distance_of_path(const int *path)
// 计算给定路径的总距离
{
    double d = 0;
    for (int j = 0; j < xCity - 1; j++)
        d += city_dis[path[j]][path[j + 1]];
    d += city_dis[path[0]][path[xCity - 1]];
    return d;
}

static int position_in(int *tmp, int C)
// 查找城市C在路径数组中的位置
{
    for (int j = 0; j < xCity; j++)
        if (tmp[j] == C)
            return j;
    return -1;
}

static void invert_segment(int pos_start, int pos_end)
{
    int j, k, t;
    if (pos_start < pos_end)
    {
        j = pos_start + 1;
        k = pos_end;
        for (; j <= k; j++, k--)
        {
            t = temp[j];
            temp[j] = temp[k];
            temp[k] = t;
        }
    }
    else
    {
        if (xCity - 1 - pos_start <= pos_end + 1)
        {
            j = pos_end;
            k = pos_start + 1;
            for (; k < xCity; j--, k++)
            {
                t = temp[j];
                temp[j] = temp[k];
                temp[k] = t;
            }
            k = 0;
            for (; k <= j; k++, j--)
            {
                t = temp[j];
                temp[j] = temp[k];
                temp[k] = t;
            }
        }
        else
        {
            j = pos_end;
            k = pos_start + 1;
            for (; j >= 0; j--, k++)
            {
                t = temp[j];
                temp[j] = temp[k];
                temp[k] = t;
            }
            j = xCity - 1;
            for (; k <= j; k++, j--)
            {
                t = temp[j];
                temp[j] = temp[k];
                temp[k] = t;
            }
        }
    }
}

static void local_improve_one(int idx)
// 对种群中索引为idx的个体进行局部优化（2-opt/变异逻辑）
{
    int C1, j, k, pos_C, pos_C1;
    int k1, k2, l1, l2, pos_flag;
    double disChange = 0;

    for (j = 0; j < xCity; j++)
        temp[j] = colony[idx][j];

    pos_flag = 0;
    pos_C = rand() % xCity;

    for (;;)
    {
        if ((rand() / 32768.0) < probab1)
        {
            do
            {
                pos_C1 = rand() % xCity;
            } while (pos_C1 == pos_C);
            C1 = colony[idx][pos_C1];
        }
        else
        {
            do
            {
                j = rand() % xColony;
            } while (j == idx);

            k = position_in(colony[j], temp[pos_C]);
            C1 = colony[j][(k + 1) % xCity];
            pos_C1 = position_in(temp, C1);
        }

        if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)
            break;

        k1 = temp[pos_C];
        k2 = temp[(pos_C + 1) % xCity];
        l1 = temp[pos_C1];
        l2 = temp[(pos_C1 + 1) % xCity];

        disChange += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
        invert_segment(pos_C, pos_C1);

        pos_flag++;
        if (pos_flag > xCity - 1)
            break;

        pos_C++;
        if (pos_C >= xCity)
            pos_C = 0;
    }

    dis_p[N_COLONY + idx] = dis_p[idx] + disChange;
    for (j = 0; j < xCity; j++)
        colony[N_COLONY + idx][j] = temp[j];
}

static void select1_local()
// 本地选择过程：如果新生成的个体更好，则替换旧个体
{
    for (int j = 0; j < xColony; j++)
    {
        if (dis_p[N_COLONY + j] < dis_p[j])
        {
            dis_p[j] = dis_p[N_COLONY + j];
            for (int k = 0; k < xCity; k++)
                colony[j][k] = colony[N_COLONY + j][k];
        }
    }
}

/* ===== 迁移（岛模型） ===== */
static void pick_k_indices(int *out_idx, int k, int pick_best)
// 根据策略选出k个最好或最坏个体的索引
{
    for (int t = 0; t < k; t++)
    {
        int chosen = -1;
        for (int i = 0; i < xColony; i++)
        {
            int used = 0;
            for (int u = 0; u < t; u++)
                if (out_idx[u] == i)
                    used = 1;
            if (used)
                continue;

            if (chosen == -1)
            {
                chosen = i;
                continue;
            }

            if (pick_best)
            {
                if (dis_p[i] < dis_p[chosen])
                    chosen = i;
            }
            else
            {
                if (dis_p[i] > dis_p[chosen])
                    chosen = i;
            }
        }
        out_idx[t] = chosen;
    }
}

static void replace_k_individuals(const int *incoming, int incoming_cnt, int replace_worst)
// 用接收到的外部个体替换本地较差的个体
{
    /* 选择要被替换的本地个体：通常替换最差的 */
    int repl_idx[N_COLONY];
    pick_k_indices(repl_idx, incoming_cnt, replace_worst ? 0 : 1);

    for (int m = 0; m < incoming_cnt; m++)
    {
        int dst = repl_idx[m];

        /* 拷贝路径 */
        const int *src_path = incoming + (m * xCity);
        for (int j = 0; j < xCity; j++)
            colony[dst][j] = src_path[j];

        dis_p[dst] = compute_distance_of_path(colony[dst]);
    }
}

static void migrate_if_needed(int rank, int size)
// 迁移逻辑主函数：按代数间隔进行进程间的数据交换
{
    if (migration_interval <= 0)
        return;
    if (GenNum == 0)
        return;
    if ((GenNum % migration_interval) != 0)
        return;

    /* 选出要发送的个体 */
    int send_idx[N_COLONY];
    int pick_best = (mig_strategy == MIG_BEST);
    pick_k_indices(send_idx, migration_count, pick_best);

    /* 打包：migration_count * xCity */
    int send_buf[N_COLONY * CITY]; /* 上限足够：100*442=44200 ints */
    for (int m = 0; m < migration_count; m++)
    {
        int src = send_idx[m];
        for (int j = 0; j < xCity; j++)
            send_buf[m * xCity + j] = colony[src][j];
    }

    if (topo == TOPO_CHAIN_RING)
    {
        int prev = (rank - 1 + size) % size;
        int next = (rank + 1) % size;

        int recv_buf[N_COLONY * CITY];

        MPI_Status st;
        MPI_Sendrecv(send_buf, migration_count * xCity, MPI_INT, next, 100,
                     recv_buf, migration_count * xCity, MPI_INT, prev, 100,
                     MPI_COMM_WORLD, &st);

        /* 注入：统一替换本地最差个体 */
        replace_k_individuals(recv_buf, migration_count, /*replace_worst=*/1);
    }
    else
    {
        /* 全连接：Allgather 收集所有rank的 migrants，然后逐批注入（忽略本rank自己的那批） */
        int total_cnt = migration_count * xCity * size;
        int *all_buf = (int *)malloc(sizeof(int) * (size_t)total_cnt);
        if (!all_buf)
        {
            fprintf(stderr, "Rank %d: malloc failed in migrate\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 3);
            exit(3);
        }

        MPI_Allgather(send_buf, migration_count * xCity, MPI_INT,
                      all_buf, migration_count * xCity, MPI_INT,
                      MPI_COMM_WORLD);

        for (int r = 0; r < size; r++)
        {
            if (r == rank)
                continue;
            const int *incoming = all_buf + r * migration_count * xCity;
            replace_k_individuals(incoming, migration_count, /*replace_worst=*/1);
        }

        free(all_buf);
    }
}
