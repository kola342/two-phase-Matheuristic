
//两阶段数学启发式


#include <iostream>
#include "CustInfo.h"
#include "HyperParameters.h"
#include "Algorithm.h"
#include "Solution.h"
#include <chrono>

using namespace std;

int main() {
	auto start = std::chrono::high_resolution_clock::now(); // 开始计时
	int i = 0;

	uniform_real_distribution<double> prob(0, 1);//0~1均匀分布
	Solution BEST('B');//最优解

	int n = 0;
	int irra = 1;
	while (true) {
		if (irra == 32)
			cout << "find" << endl;
		cout << irra++ << endl;
		Solution s_current;//构建当前解
		double best_delivery = s_current.get_makespan();//最优配送得分
		double current_delivery = best_delivery;//当前解配送得分
		s_current.solveResupply();//求解补货方案
		Solution s_best = s_current;//更新最优解

		while (T < 0.001) {
			Solution s = s_current;
			int index = s.localSearch();

			if (s.get_makespan() < best_delivery) {
				s_current = move(s);
				s_current.solveResupply();
				if (s_current.get_makespan() < s_best.get_makespan()) s_best = s_current;
				s_current.updatePoint(index, DELTA1);
			}
			else if (prob(gen) <= exp((current_delivery - s.get_makespan()) / T)) {
				s_current = move(s);
				s_current.solveResupply();
				if (s_current.get_makespan() < s_best.get_makespan()) s_best = s_current;
				s_current.updatePoint(index, DELTA2);
			}
			else {
				s_current.updatePoint(index, DELTA3);
			}

			n++;

			if (n % 10 == 0) s_current.updateWeight();//更新权重
		}

		if (s_best.get_makespan() < BEST.get_makespan()) {
			Solution::algorithm.updatePoint(s_current.construct_index, DELTA1);
			BEST = s_best;
		}

		i++;

		if (n % 10 == 0) Solution::algorithm.updateWeight();//更新权重-外层

		auto now = std::chrono::high_resolution_clock::now();   // 结束计时
		auto duration = duration_cast<std::chrono::seconds>(now - start).count();
		if (duration > Max_time) break;
	}

	BEST.cout();

	return 0;
}