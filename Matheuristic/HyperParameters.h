#pragma once

#include <random>
#include <string>

extern std::mt19937 gen;

//与算法有关的
extern const double W;//比初始解差W的解，以0.5的概率采纳
extern const double RATE_T;//温度下降速率
extern const double NOISE;//噪声比例5%
extern const int kCluster;
extern const double  PercentageRemove;//移除比率
extern const int Max_time;//最大允许时长
extern const double T_end;//冷却后的温度

//移除的随机性
extern const const float P;//y的p次方的p

//算子的更新
extern const int SEG;//算子得分的更新步长
extern const double R;//算子得分衰减率
extern const double DELTA1;//best
extern const double DELTA2;//依概率接受
extern const double DELTA3;//其他情况
extern const double OPERATOR_D;

//无人机参数
extern const double SPEED_D;//无人机的飞行速度（计算时间）
extern const double TLimit;//无人机飞行时间限制
extern const double QDd;//配送无人机的载重量
extern const double QDr;//补货无人机的载重量
extern const double Max_dis;//补货无人机的最大飞行距离

//卡车参数
extern const double SPEED_T;//卡车速度

extern const std::string address;//文件地址
