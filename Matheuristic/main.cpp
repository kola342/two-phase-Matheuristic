
//两阶段数学启发式


#include <iostream>
#include "CustInfo.h"
#include "HyperParameters.h"
#include "Algorithm.h"
#include "Solution.h"
#include <chrono>
#include "Functions.h"

using namespace std;

int main() {
	auto start = std::chrono::high_resolution_clock::now(); // 开始计时
	int outer = 1;//外部循环

	uniform_real_distribution<double> prob(0, 1);//0~1均匀分布
	Solution BEST('B');//全局最优解

	long long duration = 0;//初始化运行时长

	while (duration < Max_time) {//多起点框架
		Solution s_current;//构建当前解
		Solution::initial_LS();//初始化局部搜索
		s_current.solveResupply();//求解补货方案
		Solution s_best = s_current;//局部最优解

		localSearch(s_current, s_best, outer);//局部搜索

		//更新构造算子得分
		if (s_best.get_makespan() < BEST.get_makespan()) {
			Solution::algorithm.updatePoint(s_current.construct_index, DELTA1);
			BEST = s_best;
		}
		else {
			Solution::algorithm.updatePoint(s_current.construct_index, DELTA2);
		}
		
		if (outer % SEG == 0) Solution::algorithm.updateWeight();//更新权重-外层

		auto now = std::chrono::high_resolution_clock::now();   // 结束计时

		duration = duration_cast<std::chrono::seconds>(now - start).count();//计算运行时长
		
		++outer;
	}

	BEST.cout();//输出最优结果

	return 0;
}