#pragma once

#include <vector>
#include "HyperParameters.h"
#include "CustInfo.h"

using namespace std;

//算法权重类
class Algorithm
{//记录共享数据
public:
	//获取构造算子权重
	const vector<double>& get_Construct() const {
		return this->weightConstruct;
	}
	//获取局部搜索算子的权重的引用
	vector<double> get_LS() {
		return this->weightLS;
	}
	//算法初始化
	Algorithm(CustInfo& info) {
		if (info.get_k() == 1) {
			this->weightConstruct[2] = 0.0;
		}
	}
	//更新构造算子的的得分
	void updatePoint(int index, double delta) {
		this->pointConstruct[index] += delta;
		this->numConstruct[index] += 1;
	}
	//更新构造算子的权重
	void updateWeight() {
		for (size_t i = 0; i < this->weightConstruct.size(); i++) {
			if (this->numConstruct[i] != 0) {
				this->weightConstruct[i] = this->weightConstruct[i] * (1 - R) + R * this->pointConstruct[i] / this->numConstruct[i];
			}
		}
	}

private:
	vector<double> weightConstruct = vector<double>(3, OPERATOR_D);//构造启发式的权重
	vector<double> weightLS = vector<double>(6, OPERATOR_D);//局部搜索的权重
	vector<double> pointConstruct = vector<double>(3, 0);//构造启发式的得分
	vector<int> numConstruct = vector<int>(3, 0);//构造启发式的的使用次数
};

