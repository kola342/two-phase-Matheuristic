#pragma once

#include <random>
#include <string>

extern std::mt19937 gen;
extern const double Noise;//噪声比例
extern const double TLimit;//配送无人机时长限制
extern const double	QDd;//配送无人机最大载重
extern const double	QDr;//补货无人机最大载重
extern const double SpeedD;//无人机速度
extern const double SpeedT;//卡车速度
extern const double OPERATOR_D;
extern const int kCluster;
extern const double Max_dis;//补货无人机的最大飞行距离
extern const double  PercentageRemove;
extern const float P;

extern const std::string address;
extern const int Max_time;
extern const double T;
extern const double K;
extern const double DELTA1, DELTA2, DELTA3;
extern const double R;