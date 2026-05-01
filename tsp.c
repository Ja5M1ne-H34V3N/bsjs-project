#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#define FILE_PATH  "./pcb442.tsp","r"  //�����ļ���
#define N_COLONY 100  // N_COLONY>=xColony
#define CITY     442  // CITY>=xCity
int     xColony = 100;     //##//  ������
int     xCity = CITY;
double  probab1 = 0.02;    //##//  �������
long    NOCHANGE = 200000;  //##//  ���ֹͣ�ı����
long    maxGen = 200000;    //##//  ͣ������
int     colony[N_COLONY * 2][CITY], colony2[N_COLONY][CITY]; //zhongqun 
double  cityXY[CITY][2];
double  city_dis[CITY][CITY];
double  dis_p[N_COLONY * 2]; //��Ӧֵ
double  sumbest, sumTemp;
int     temp[CITY], ibest;
clock_t timeStart, timeNow, timeTemp;
long    GenNum, Ni;
void    init();
int     position(int *tmp, int C);
void    invert(int pos_start, int pos_end);
void    printBest(long GenNum);
void    tempTest(int i);
void    mapped();
void    LastCP();
double  path(int tmp[], int k1, int k2);
void select1();
void select2();
double SPAD_compute();
FILE *fpme;
int main()
{
	register int C1, j, k, pos_C, pos_C1; int k1, k2, l1, l2, pos_flag, icount, mem;
	register double disChange;
	static int i = 0;
	timeStart = timeNow = timeTemp = clock();
	init();
	for (;;)//进行无限次循环
	{
		for (j = 0; j<xCity; j++)temp[j] = colony[i][j];
		disChange = 0; pos_flag = 0;
		pos_C = rand() % xCity;
		for (;;)
		{
			if ((rand() / 32768.0)<probab1)     //�ڱ�������
			{
				do
				pos_C1 = rand() % xCity;
				while (pos_C1 == pos_C);
				C1 = colony[i][pos_C1];
			}
			else
			{
				do
				j = rand() % xColony;
				while (j == i);
				k = position(colony[j], temp[pos_C]);
				C1 = colony[j][(k + 1) % xCity];
				pos_C1 = position(temp, C1);
			}
			if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)break;
			k1 = temp[pos_C];
			k2 = temp[(pos_C + 1) % xCity];
			l1 = temp[pos_C1];
			l2 = temp[(pos_C1 + 1) % xCity];
			disChange += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
			invert(pos_C, pos_C1);
			pos_flag++;
			if (pos_flag>xCity - 1)
				break;  ////////////
			pos_C++;
			if (pos_C >= xCity)pos_C = 0;                 /**********************/
		}
		dis_p[N_COLONY + i] = dis_p[i] + disChange;
		disChange = 0;
		for (j = 0; j<xCity; j++)
			colony[N_COLONY + i][j] = temp[j];
		i++;
		if (i >= xColony)//�˴��Ǳ���+�Ӵ��� 
		{
			select1();
			Ni++; GenNum++; i = 0;
			sumbest = dis_p[0];
			for (j = 0; j<N_COLONY; j++)
			if (sumbest>dis_p[j])
				sumbest = dis_p[j];
			printf("%d:%f\n", GenNum, sumbest);
			if (GenNum % 2000 == 0 && GenNum<maxGen)
				printBest(GenNum);
			if (GenNum >= maxGen)
			{
				timeNow = clock();
				printf("Finial solution:");
				printBest(GenNum);
				exit(1);
			}
		}
	}

	return 1;
}
void select1()
{
	int j, k;
	for (j = 0; j<N_COLONY; j++)
	if (dis_p[N_COLONY + j]<dis_p[j])
	{
		dis_p[j] = dis_p[N_COLONY + j];
		for (k = 0; k<CITY; k++)
			colony[j][k] = colony[N_COLONY + j][k];
	}
}
void init()//初始化种群的信息
{
	int i, j, t, sign, mod, array[CITY];
	double x, y;
	double d;
	FILE *fp;
	srand((unsigned)time(NULL));//设置随机数
	if ((fp = fopen(FILE_PATH)) == NULL)exit(0);
	fscanf(fp, "%d", &xCity);
	for (i = 0; i<xCity; i++)      /*  init cityXY[][]  */
	{
		fscanf(fp, "%*d%Lf%Lf", &x, &y);
		cityXY[i][0] = x;
		cityXY[i][1] = y;
	}
	fclose(fp); //读取城市信息，并将其坐标存在数组中

	for (i = 0; i<xCity; i++)    /*  init city_dis[] */
	for (j = 0; j<xCity; j++)
	{
		if (j>i)
		{
			d = (cityXY[i][0] - cityXY[j][0])*(cityXY[i][0] - cityXY[j][0])*1.0 +
				(cityXY[i][1] - cityXY[j][1])*(cityXY[i][1] - cityXY[j][1])*1.0;
			city_dis[i][j] = (int)(sqrt(d) + 0.5);
			continue;
		}
		if (j == i) { city_dis[i][j] = 0; continue; }
		if (j<i)  city_dis[i][j] = city_dis[j][i];
	}// 初始化城市的距离矩阵，表示城市之间的距离

	mod = xCity;//  我们将xcity 也就是城市的总数作为后续计算的初始取模标准
	for (i = 0; i<xCity; i++)array[i] = i;     //初始化城市选取列表array     
	for (i = 0; i<xColony; i++, mod = xCity) // 接下来定义当前种群里的每一个个体 根据变量xColony的信息，我们知道本次种群初始设置为100个，同时在每定义一个新的个体的时候，都重置模系数的值为城市总数
	for (j = 0; j<xCity; j++)  // 接下来进行xcity次循环 旨在填好某个特定的个体
	{
		sign = rand() % mod; //定义sign值，意为随机选取一个城市作为填入这个个体这一位的城市
		colony[i][j] = array[sign]; //i指的是特定的个体 j指的是个体中第j个城市，这里被填入sign 也就是刚刚随机出来的城市
		t = array[mod - 1];  
		array[mod - 1] = array[sign];
		array[sign] = t;//这三个是换位的代码，将array 也就是城市列表第mod - 1 个城市和刚刚选中的随机城市互换
		mod--; //将mod见一，也就是下一轮循环中，只在array的前mod - 1个城市里进行随机抽取，这样做的目的是将已经被选择的城市移到后面就能实现其不会再被选择（比如mod系数是5肯定不会出现结果6）
		if (mod == 1) colony[i][++j] = array[0]; //当mod是1的时候，直接选剩下的最后一个城市就好了，方便快捷
	}// 这样 经过i次循环 整个种群的个体就定义好了

	for (i = 0; i<xColony; i++)	// 接下来遍历整个种群	    /*    init dis_p[]       */
	{
		dis_p[i] = 0;//dis_p保存着每个种群的全部距离 ，在算法里就是适应值
		for (j = 0; j<xCity - 1; j++)
			dis_p[i] = dis_p[i] + city_dis[*(*(colony + i) + j)][*(*(colony + i) + j + 1)];
		dis_p[i] = dis_p[i] + city_dis[**(colony + i)][*(*(colony + i) + xCity - 1)];
	} //计算一个个体的具体的适应值，这里就是按照个体保存的城市顺序走一圈的总距离

	ibest = 0; sumbest = dis_p[0];	    /*  init ibest & sumbest */
	sumTemp = sumbest * 5;
	GenNum = 0;	Ni = 0;               /*   initialize GunNum & Ni    */
	printf("init success!!!\n");//将最强个体ibest初始化为0号个体 最好长度初始化为第一个个体的长度 初始化代数计数器为0
}

void invert(int pos_start, int pos_end)
{
	int j, k, t;
	if (pos_start<pos_end)
	{
		j = pos_start + 1; k = pos_end;
		for (; j <= k; j++, k--)
		{
			t = temp[j]; temp[j] = temp[k]; temp[k] = t;
		}
	}
	else
	{
		if (xCity - 1 - pos_start <= pos_end + 1)
		{
			j = pos_end; k = pos_start + 1;
			for (; k<xCity; j--, k++)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
			k = 0;
			for (; k <= j; k++, j--)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
		}
		else
		{
			j = pos_end; k = pos_start + 1;
			for (; j >= 0; j--, k++)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
			j = xCity - 1;
			for (; k <= j; k++, j--)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
		}
	}
}
int position(int *tmp, int C)
{
	int j;
	for (j = 0; j<xCity; j++)
	if (*(tmp + j) == C)break;
	return(j);
}
void printBest(long GenNum)
{
	int i;
	if ((fpme = fopen("e:\\tsp0.txt", "a")) == NULL)exit(0);
	//fprintf(fpme,"\n   CITY      %d\t\tN_COLONY  %d",CITY,N_COLONY);
	//fprintf(fpme,"\ntime     %4.2f",(double)(timeNow-timeStart)/CLOCKS_PER_SEC);
	//fprintf(fpme,"\n   distance  %f",sumbest);
	//fprintf(fpme,"\n   GenNum    %d\n\n",GenNum);
	fprintf(fpme, "%d\t%4.2f\t%d\n", GenNum, (double)(timeNow - timeStart) / CLOCKS_PER_SEC, (int)sumbest);
	fclose(fpme);

}
double path(int tmp[], int k1, int k2)
{
	int j, t1, t2; double temp_dis = 0;
	if (k2>k1)
	for (j = k1; j<k2; j++)
		temp_dis += city_dis[tmp[j]][tmp[j + 1]];
	else
	for (j = k1; j<k2 + xCity; j++)
	{
		t1 = j%xCity; t2 = (j + 1) % xCity;
		temp_dis += city_dis[tmp[t1]][tmp[t2]];
	}
	return temp_dis;
}

