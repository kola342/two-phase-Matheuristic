#pragma once
#include <unordered_map>
#include <set>
#include "Solution.h"

using namespace std;

//第二阶段的C_D
struct SeriesTime
{
	int id = 0;
	double t = 0;//卡车出发时间
	double e = 0;//卡车等待时间
	set<int> AP;//可以接收的包裹
	SeriesTime(int id, double t, double e, const set<int> AP) :t(t), e(e), AP(AP), id(id) {};
};

struct PairHash {
	size_t operator()(const std::pair<int, int>& p) const noexcept {
		auto h1 = std::hash<int>{}(p.first);
		auto h2 = std::hash<int>{}(p.second);
		return h1 ^ (h2 << 1); // 简单且常用的组合方式
	}
};

class ResupplyInfo
{
public:
	ResupplyInfo(list<TruckCust>& solu, CustInfo& c) :info(c) {
		set<int> package;
		//vector<SeriesTime> temp;
		for (auto it = solu.rbegin(); it != solu.rend(); ++it) {
			if (it->id == 0) continue;
			//收集可补货的货物集合
			if (c.at_custDelay(it->id) != 0) package.insert(it->id);
			if (!it->droneRoute.empty()) {
				auto last = prev(it->droneRoute.end());
				for (auto it1 = next(it->droneRoute.begin()); it1 != last; it1++) {
					if (c.at_custDelay(it1->id) != 0) package.insert(it1->id);
				}
			}

			if (c.in_CD(it->id)) this->C_D.emplace_back(it->id, it->t, it->e, package);
		}		

		reverse(this->C_D.begin(), this->C_D.end());

		C_delay = vector<int>(C_D[0].AP.begin(), C_D[0].AP.end());

		//处理P
		for (auto it = solu.begin(); it != solu.end(); ++it) {
			if (it->droneRoute.empty()) continue;
			auto temp_i = it;
			while (temp_i != solu.begin() && !c.in_CD(temp_i->id)) { --temp_i; }
			auto temp_j = prev(it->end);
			while (temp_j != solu.begin() && !c.in_CD(temp_j->id)) { --temp_j; }
			if (temp_i->id != temp_j->id) {
				double tm = TLimit - (it->end->t - it->end->e - it->t);
				if (tm <= 0) {
					int y = 0;
				}
				this->P[{ temp_i->id, temp_j->id }] = TLimit - (it->end->t - it->end->e - it->t);
			}
		}

		this->n_1 = c.get_num() + 1;
	}

	void solve(map<int, DroneResupply>& obj) {
		this->branch_benders_cut(obj);
	}

	double get_dt(unsigned int a) {
		return this->dt.at(a);
	}

private:
	CustInfo& info;
	vector<SeriesTime> C_D;//可补货节点
	vector<int> C_delay;
	int n_1;
	vector<double> dt;
	unordered_map<pair<int, int>, double, PairHash> P;
	void branch_benders_cut(map<int, DroneResupply>&);
	void gurobi_solve(map<int, DroneResupply>&);
};

