#include"HyperParameters.h"

std::mt19937 gen(10);//随机种子

//与算法有关的
const double W = 0.5;//比初始解差W的解，以0.5的概率采纳
const double RATE_T = 0.997;//温度下降速率
const double NOISE = 0.05;//噪声比例5%
const int kCluster = 1;//初始簇个数
const double  PercentageRemove = 0.2;//移除比率
const int Max_time = 10000;//最大允许时长
const double T_end = 0.001;//冷却后的温度

//移除的随机性
const const float P = 2;;//y的p次方的p

//算子的更新
const int SEG = 10;//算子得分的更新步长
const double R = 0.1;//算子得分衰减率
const double DELTA1 = 15;//best
const double DELTA2 = 10;//依概率接受
const double DELTA3 = 5;//其他情况
const double OPERATOR_D = 1;

//无人机参数
const double SPEED_D = 11.5;//无人机的飞行速度（计算时间）
const double TLimit = 1800;//无人机飞行时间限制
const double QDd = 4;//配送无人机的载重量
const double QDr = 10;//补货无人机的载重量
const double Max_dis = 7000;//补货无人机的最大飞行距离

//卡车参数
const double SPEED_T = 11.5;//卡车速度

const std::string address = "D:\\Document\\VS_code\\TSPD-RD\\u_100.txt";