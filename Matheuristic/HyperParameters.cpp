#include"HyperParameters.h"

std::mt19937 gen(10);
const double Noise = 0.05;//噪声比例
const double TLimit = 1800;//无人机飞行时间限制
const double QDd = 4;//配送无人机的载重量
const double QDr = 10;//补货无人机的载重量
const double SpeedT = 10;//卡车速度
const double SpeedD = 10;//无人机速度
const double OPERATOR_D = 0.9;
const int kCluster = 1;
const double Max_dis = 5000;//补货无人机的最大飞行距离
const double  PercentageRemove = 0.2;//移除比率
const float P = 2;
const int Max_time = 10000;
const double T= 100;
const double K = 0.95;
const double DELTA1 = 9, DELTA2 = 8, DELTA3 = 5;
const double R = 0.5;

const std::string address = "D:\\Document\\VS_code\\TSPD-RD\\u_100.txt";