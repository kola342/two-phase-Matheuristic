#include "Functions.h"
#include <algorithm>
#include <iostream>

using namespace std;

void setT(double& T, double cost) {
	T = -cost * W / log(0.5);
	cout << "初始温度：" << T << endl;
}

void localSearch(Solution& s_current, Solution& s_best, int outer) {
	int inner = 1;
	uniform_real_distribution<double> prob(0, 1);//0~1均匀分布

	double T;	
	setT(T, s_best.get_makespan());

	while (T > T_end) {
		Solution s = s_current;
		int index = s.localSearch();

		cout << outer << "." << inner << ":" << endl;
		cout << "	新解:" << endl;
		cout << "		配送成本：" << s.get_delivery_makespan() << endl;

		if (s.get_delivery_makespan() < s_best.get_delivery_makespan()) {
			s_current = move(s);
			s_current.solveResupply();
			if (s_current.get_makespan() < s_best.get_makespan()) 
				s_best = s_current;
			s_current.updatePoint(index, DELTA1);
		}
		else if (prob(gen) <= exp((s_current.get_delivery_makespan() - s.get_delivery_makespan()) / T)) {
			s_current = move(s);
			s_current.solveResupply();
			if (s_current.get_makespan() < s_best.get_makespan()) 
				s_best = s_current;
			s_current.updatePoint(index, DELTA2);
		}
		else {
			s_current.updatePoint(index, DELTA3);
		}

		if (inner % SEG == 0) 
			s_current.updateWeight();//更新局部搜索权重
		cout << "	当前解:" << endl;
		cout << "		配送成本：" << s_current.get_delivery_makespan();
		cout << "		总成本：" << s_current.get_makespan() << endl;
		cout << "	局部最优解:" << endl;
		cout << "		配送成本：" << s_best.get_delivery_makespan();
		cout << "		总成本：" << s_best.get_makespan() << endl;

		T *= RATE_T;
		++inner;
	}
}