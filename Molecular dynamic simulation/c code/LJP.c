#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define N (6 * 6 * 6)    // 原子的个数
#define L 30             // 盒子的边长
#define mass 1           // 原子的质量
#define dt 0.00005       // 步长 单位 秒
#define T 50             // 温度
#define stepLimit 10000  // 步数
#define epsilon 1        // LJ势参数
#define sigma 1          // LJ势参数
#define rc (2.5 * sigma) // 截止距离
#define V_0 1.0          // 高斯势阱深度
#define a1 0.0025
#define s -0.02               // 高斯势阱变化率
#define rc_g (3.0 / sqrt(a1)) // 高斯势阱的截断距离
#define kb 1.0                // 玻尔兹曼常数
#define PI 3.1415926          //

typedef struct
{
    double pos[3];
    double vel[3];
    double a[3];
} Mol;

typedef struct
{
    double val;
    double sum;
    double sum_squ;
} Prop;

Mol mol[N];
int color[N];
double timeNow = 0;
double vSum[3] = {0.0, 0.0, 0.0};
double kinEnergy;
double temps = 0.0;
int number = N; 
int stepCount = 0;

double norm(double v[3])
{
    if (v == NULL) {
        fprintf(stderr, "Error: NULL pointer passed to norm function.\n");
        return -1.0; // Return a negative value to indicate an error
    }
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void GaussianForce(double v[3], double force[3])
{
    double r = norm(v);
    double V_0_new = V_0 + s * timeNow;
    for (int i = 0; i < 3; i++)
    {
        force[i] = -V_0_new * 2 * a1 * exp(-a1 * r * r) * v[i];
    }
}

void LeapfrogStep(int part)
{
    if (part == 1)
    {
        for (int n = 0; n < N; n++)
        {
            for (int i = 0; i < 3; i++)
            {
                mol[n].vel[i] += 0.5 * dt * mol[n].a[i];
                mol[n].pos[i] += dt * mol[n].vel[i];
            }
        }
    }
    else
    {
        for (int n = 0; n < N; n++)
        {
            for (int i = 0; i < 3; i++)
                mol[n].vel[i] += 0.5 * dt * mol[n].a[i];
        }
    }
}

void ComputeForces()
{
    double rrCut = rc * rc;
    for (int n = 0; n < N; n++)
    {
        for (int i = 0; i < 3; i++)
            mol[n].a[i] = 0.0;
    }

    for (int j1 = 0; j1 < N - 1; j1++)
    {
        for (int j2 = j1 + 1; j2 < N; j2++)
        {
            double dr[3];
            for (int i = 0; i < 3; i++)
                dr[i] = mol[j1].pos[i] - mol[j2].pos[i];

            double rr = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (rr < rrCut)
            {
                double r = sqrt(rr);
                double fcVal = 48 * epsilon * pow(sigma, 12) / pow(r, 14) - 24 * epsilon * pow(sigma, 6) / pow(r, 8);
                for (int i = 0; i < 3; i++)
                {
                    mol[j1].a[i] += fcVal * dr[i] / mass;
                    mol[j2].a[i] -= fcVal * dr[i] / mass;
                }
                if (norm(mol[j1].pos) < rc_g)
                {
                    double force[3];
                    GaussianForce(mol[j1].pos, force);
                    for (int i = 0; i < 3; i++)
                        mol[j1].a[i] += force[i] / mass;
                }
                if (norm(mol[j2].pos) < rc_g)
                {
                    double force[3];
                    GaussianForce(mol[j2].pos, force);
                    for (int i = 0; i < 3; i++)
                        mol[j1].a[i] += force[i] / mass;
                }
            }
        }
    }
}
void EvalProps()
{
    double vvSum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        vSum[i] = 0.0;
    }

    for (int n = 0; n < N; n++)
    {
        if (color[n] == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                vSum[i] += mol[n].vel[i];
            }
            double vv = mol[n].vel[0] * mol[n].vel[0] + mol[n].vel[1] * mol[n].vel[1] + mol[n].vel[2] * mol[n].vel[2];
            vvSum += vv;
        }
    }

    kinEnergy = (0.5 * vvSum) / number;
    temps = (kinEnergy * 2) / (3 * kb);
}

void singlestep()
{
    stepCount += 1;
    timeNow = stepCount * dt;
    LeapfrogStep(1);
    ComputeForces();
    LeapfrogStep(2);
    EvalProps();
}

double rand_normal()
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    return z0;
}

int main()
{
    int n = cbrt(N), x = 0, y = 0, z = 0;
    double d = L / n;
    // 初始化分子位置
    for (int i = 0; i < N; i++)
    {
        x = (i + 1) / (n * n);
        y = ((i + 1) % (n * n)) / n;
        z = ((i + 1) % (n * n)) - y * n;
        mol[i].pos[0] = x * d - L / 2 + d / 2;
        mol[i].pos[1] = y * d - L / 2 + d / 2;
        mol[i].pos[2] = z * d - L / 2 + d / 2;
    }
    for (int i = 0; i < 3; i++)
    {
        mol[N - 1].pos[i] = -L / 2 + d / 2;
    }
    // 初始化分子速度
    srand(time(NULL));
    int moreCycles = 1;
    double scale = sqrt(kb * T / mass);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mol[i].vel[j] = (1.0 / sqrt(3.0)) * rand_normal() * scale;
            vSum[j] += mol[i].vel[j];
        }
    }
    
    //调整速度分布，使总速度为0
    for (int j = 0; j < 3; j++)
    {
        vSum[j] /= N;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mol[i].vel[j] -= vSum[j];
        }
    }
    //初始化color数列
    for (int i = 0; i < N; i++)
    {
        color[i]=0;
    }
    
    // 输出文件
    FILE *fp = fopen("LJP_c/data_c.xyz", "w");
    FILE *fp_temp = fopen("LJP_c/temperature_c.txt", "w");
    FILE *fp_number = fopen("LJP_c/number_c.txt", "w");
    while (moreCycles)
    {
        for (int n; n<N;n++){
            if (norm(mol[n].pos)>=rc_g && color[n]==0 ){
                number -= 1;
                color[n] =1;   
            }          }    
        singlestep();//求解每一步的运动方程
        // 输出 XYZ 格式的数据
        fprintf(fp, "%d\nFrame %d\n", N, stepCount);
        for (int i = 0; i < N; i++)
        {
            fprintf(fp, "Ar %.4f %.4f %.4f\n", mol[i].pos[0], mol[i].pos[1], mol[i].pos[2]);
        }
        
        fprintf(fp_temp, "Temperature: %.4f K\n", temps);
        fprintf(fp_number, "number: %d K\n", number);
        if (stepCount>=stepLimit || number == 0 )
        {
            moreCycles=0;
        }
        
    }
    
    fclose(fp);
    fclose(fp_temp);
    fclose(fp_number);
}
