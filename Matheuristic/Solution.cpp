#include "Solution.h"
#include "HyperParameters.h"
#include <numeric>
#include <algorithm>
#include "ResupplyInfo.h"

//初始化Solution
CustInfo Solution::info = move(CustInfo(address));//客户信息
Algorithm Solution::algorithm = move(Algorithm(Solution::info));//权重信息

//选择构造启发式
int Solution::intinialize(const vector<double>& weightConstruct) {//构造初始解
	discrete_distribution<> prob(weightConstruct.begin(), weightConstruct.end());
	int t = prob(gen);

	vector<function<void()>> construct_functions = {
		[this]()->void { this->greedyConstructNoise(); },
		[this]()->void { this->randomGreedyConstruct(); },
		[this]()->void { this->greedyConstructNoiseCluster(); } };//构造函数的数组

	construct_functions[t]();
	return t;
}

//局部搜索
int Solution::localSearch() {
	discrete_distribution<> prob(this->weightLS.begin(), weightLS.end());
	int t = prob(gen);

	vector<function<void()>> operatorLS = {
			[this]()->void {
				this->randomRemove();
				list<int> temp;
				this->seperate(temp);
				
				this->greedyInsert(temp);
				this->greedyInsert(this->custInsert);
				this->greedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->randomRemove();
				list<int> temp;
				this->seperate(temp);

				this->randomGreedyInsert(temp);
				this->randomGreedyInsert(this->custInsert);
				this->randomGreedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->shawRemove();
				list<int> temp;
				this->seperate(temp);

				this->greedyInsert(temp);
				this->greedyInsert(this->custInsert);
				this->greedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->shawRemove();
				list<int> temp;
				this->seperate(temp);

				this->randomGreedyInsert(temp);
				this->randomGreedyInsert(this->custInsert);
				this->randomGreedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->worstRemove();
				list<int> temp;
				this->seperate(temp);

				this->greedyInsert(temp);
				this->greedyInsert(this->custInsert);
				this->greedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->worstRemove();
				list<int> temp;
				this->seperate(temp);

				this->randomGreedyInsert(temp);
				this->randomGreedyInsert(this->custInsert);
				this->randomGreedyInsert(this->custInsert_delay);
			}
	};//局部搜索函数的数组

	operatorLS[t]();
	return t;
}

//贪婪构造
void Solution::greedyConstruct(list<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		//计算各个待插入客户的最优插入位置
		for (auto it = a.begin(); it != a.end(); it++) {
			this->calBestInsertion(it, best, this->solu.end(), Noise);
		}

		this->insert(best);//执行插入操作
		a.erase(best.id);//更新待插入集合
	}
}

void Solution::greedyConstruct(list<int>& a, list<TruckCust>::iterator solu_end) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		for (auto it = a.begin(); it != a.end(); it++) {
			this->calBestInsertion(it, best, solu_end, Noise);
		}

		this->insert(best);
		a.erase(best.id);//更新待插入集合
	}
}

//随机贪婪构造
void Solution::randomGreedyConstruct(list<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		//随机选择客户
		uniform_int_distribution<int> prob(0, static_cast<int>(a.size()) - 1);
		auto it = a.begin();
		advance(it, prob(gen));

		this->calBestInsertion(it, best, this->solu.end());//最优插入客户);

		this->insert(best);
		a.erase(best.id);//更新待插入集合
	}
}

//贪婪构造_簇情况下
void Solution::greedyConstructNoiseCluster() {
	list<list<TruckCust>> solu_part;//部分解集合

	for (auto& cluster : this->info.at_clusters()) {
		this->custInsert.clear();
		this->custInsert_delay.clear();
		for (auto customer : cluster) {
			if (this->info.at_custDelay(customer) == -1) this->custInsert.push_back(customer);
			else this->custInsert_delay.push_back(customer);
		}
		this->greedyConstructNoise();//构造启发式
		solu_part.push_back(move(this->solu));
		this->initial_solu();//解的初始化
	}

	//节约里程法
	while (solu_part.size() != 1) {
		pair<list<list<TruckCust>>::iterator, list<list<TruckCust>>::iterator> best;
		double best_value = static_cast<double>(INT_MAX);
		for (auto i = solu_part.begin(); i != solu_part.end(); i++) {
			int a = prev(i->end(), 2)->id;
			for (auto j = solu_part.begin(); j != solu_part.end(); j++) {
				if (i == j) continue;
				int b = next(j->begin())->id;
				double temp = this->info.at_TT(a, b) - this->info.at_TT(a, 0) - this->info.at_TT(0, b);
				if (temp < best_value) {
					best_value = temp;
					best = { i,j };
				}
			}
		}

		best.first->erase(prev(best.first->end()));
		best.first->insert(best.first->end(), next(best.second->begin()), best.second->end());
		solu_part.erase(best.second);
	}
	this->solu = move(solu_part.front());
	this->constraintCD();
	this->calMakespan(this->solu.begin());
}


/*
* 给定客户序号obj，查找在路线初始到solu_end的最优插入位置信息best，给定噪声max_noise
*/
void Solution::calBestInsertion(list<int>::iterator obj, Insertion& best, list<TruckCust>::iterator solu_end, double max_noise) {
	
	uniform_real_distribution<double> noise(-max_noise, max_noise);

	auto td = this->info.at_custDelay(*obj);//延迟时间

	for (auto it = next(this->solu.begin()); it != solu_end; ++it) {
		auto pre = prev(it);//记录it的前一个卡车客户

		//插入卡车的情况
		double delta = this->info.at_ST(*obj) + this->info.at_TT(pre->id, *obj) + this->info.at_TT(*obj, it->id) - this->info.at_TT(pre->id, it->id);//插入代价
		double increment = max(td + this->info.at_TD(0, *obj) - (pre->t + this->info.at_TT(pre->id, *obj)), 0.0);//时间窗造成的延迟时间，
		delta += increment;//综合插入代价
		if (pre->numDrone == 0) {//插入在无人机跨度的卡车路线
			auto start = pre;
			for (; start != this->solu.begin(); --start) {
				if (!start->droneRoute.empty()) break;
			}

			//修改无人机路线的降落点
			auto end = start->end;
			for (; end != it; --end) {
				if (delta + end->t - end->e - start->droneRoute.front().t <= TLimit) break;
			}//插入点的前半部分
			if (end == it) {
				for (end = it; end != start; --end) {
					if (end->t - end->e - start->droneRoute.front().t <= TLimit) break;
				}
				if (end == start) delta = static_cast<double>(INT_MAX);
			}//插入点的后半部分
			delta += delta * noise(gen);
			if (delta < best.delta)	best = Insertion(delta, 'x', obj, it, start, end);
		}
		else {
			delta += delta * noise(gen);
			if (delta < best.delta)	best = Insertion(delta, 't', obj, it);
		}


		//无人机的情况
		if (this->info.at_demand(*obj) > QDd) continue;

		if (pre->numDrone) {//新增无人机的情况，以pre为起点
			increment = max(td + this->info.at_TD(0, it->id) - (pre->t + this->info.at_TT(pre->id, *obj)), 0.0);//时间窗造成的延迟时间
			for (auto it1 = it; it1 != this->solu.end(); ++it1) {//寻找降落点
				if (it1->id == pre->id) break;
				double span = this->info.at_TD(pre->id, *obj) + this->info.at_TD(*obj, it1->id) + this->info.at_SD(*obj);//无人机需要飞行的最小时长
				if (span <= TLimit) {//满足无人机飞行限制
					if (it1->t - it1->e - pre->t > TLimit) break;

					delta = (span + pre->t - pre->e) - it1->t;//为多个客户的情况提供可能
					delta = (delta + increment) * (1 + noise(gen));
					if (delta < best.delta) best = Insertion(delta, 'n', obj, pre, it1);
				}
				if (!it1->droneRoute.empty()) break;
			}
		}
		else if (!pre->droneRoute.empty()
			&& pre->droneRoute.front().q + this->info.at_demand(*obj) <= QDd) {//插入无人机路线的情况，pre作起点
			for (auto it1 = next(pre->droneRoute.begin()); it1 != pre->droneRoute.end(); ++it1) {
				increment = max(td + this->info.at_TD(0, pre->id) - pre->t, 0.0);//时间窗造成的延迟时间

				auto pre1 = prev(it1);//插入位置的前一个客户
				delta = this->info.at_SD(*obj) + this->info.at_TD(pre1->id, *obj) + this->info.at_TD(*obj, it1->id) - this->info.at_TD(pre1->id, it1->id);
				double span = pre->droneRoute.back().t - pre->droneRoute.back().e + delta;//插入后，到达降落点的时间
				if (span - pre->droneRoute.front().t > TLimit) continue;

				delta = (span - pre->end->t);//插入后，无人机路线的代价
				bool up = false;//是否修改无人机降落点
				if (delta > 0 && pre->end->numDrone != 0) {//尝试修改降落点
					auto end = prev(pre->droneRoute.end(), 2);
					for (auto it2 = next(pre->end); it2 != this->solu.end(); ++it2) {
						if (it2->numDrone == 0 && it2->droneRoute.empty()) break;
						if (it2->id == pre->id) continue;
						if (it2->t - it2->e - pre->t > TLimit) break;
						double new_span = span - this->info.at_TD(end->id, pre->end->id) + this->info.at_TD(end->id, it2->id);
						if (new_span < TLimit && new_span < it2->t) {
							delta = increment * (1 + noise(gen));
							if (delta < best.delta) {
								best = Insertion(delta, 'o', obj, pre, it1, it2);
								up = true;
							}
							break;
						}
					}
				}
				delta += increment;
				delta *= 1 + noise(gen);
				if (!up && delta < best.delta) best = Insertion(delta, 'd', obj, pre, it1);
			}
		}
	}
}


//插入并更新决策变量，计算目标函数
void Solution::insert(Insertion& best) {
	switch (best.kind)
	{
	case 't'://插入卡车路线
		this->solu.insert(best.loc_T, TruckCust(*best.id));
		this->calMakespan(prev(best.loc_T));//定位到插入点
		break;
	case 'n'://新增无人机路线
		best.loc_T->droneRoute.emplace_back(best.loc_T->id);
		best.loc_T->droneRoute.emplace_back(*best.id);
		best.loc_T->droneRoute.emplace_back(best.land->id);
		best.loc_T->end = best.land;
		best.land->drone_end = &(best.loc_T->droneRoute.back());
		best.loc_T->droneRoute.front().q += this->info.at_demand(*(best.id));
		this->calMakespan(best.loc_T);//定位到插入点
		break;
	case 'd'://插入无人机路线
		best.loc_T->droneRoute.insert(best.loc_D, DroneCust(*best.id));
		best.loc_T->droneRoute.front().q += this->info.at_demand(*(best.id));
		this->calMakespan(best.loc_T);//定位到插入点
		break;
	case 'o'://改变降落点
		best.loc_T->droneRoute.insert(best.loc_D, DroneCust(*best.id));
		best.loc_T->droneRoute.back().id = best.land->id;
		best.loc_T->end->drone_end = nullptr;
		best.loc_T->end = best.land;
		best.land->drone_end = &(best.loc_T->droneRoute.back());
		best.loc_T->droneRoute.front().q += this->info.at_demand(*(best.id));
		this->calMakespan(best.loc_T);//定位到插入点
		break;
	case 'x'://插入卡车，修改降落点
		this->solu.insert(best.loc_T, TruckCust(*best.id));
		best.start->droneRoute.back().id = best.land->id;
		best.start->end->drone_end = nullptr;
		best.start->end = best.land;
		best.land->drone_end = &(best.start->droneRoute.back());
		best.land->numDrone = 1;
		this->calMakespan(best.start);//定位到插入点
		break;
	}
}

/*
* 计算目标函数，从obj点开始计算后续makespan
*/
void Solution::calMakespan(list<TruckCust>::iterator obj) {
	list<TruckCust>::iterator start;
	bool s = false;//表示是否存在start

	//寻找start点
	if (obj != this->solu.begin() && obj->droneRoute.empty() && prev(obj)->numDrone == 0) {
		for (start = prev(obj); start != this->solu.begin(); --start) {
			if (!start->droneRoute.empty()) break;
		}
		s = true;
		obj->numDrone = 0;
	}

	//遍历更新决策变量
	for (auto it = obj; it != this->solu.end(); it++) {
		double td = this->info.at_custDelay(it->id);//延迟时间
		double arrival = this->info.at_TD(0, it->id);//假设，补货无人机飞行时间		

		if (s && it == start->end) {//计算降落点
			if (!it->droneRoute.empty()) it->numDrone = 0;
			else it->numDrone = 1;
			auto last = prev(start->droneRoute.end(), 2);//无人机服务的最后一个客户
			double span = last->t + this->info.at_TD(last->id, start->droneRoute.back().id);
			double estimate = prev(it)->t + this->info.at_TT(prev(it)->id, it->id);//到达时间
			start->droneRoute.back().t = max(span, estimate);
			start->droneRoute.back().e = start->droneRoute.back().t - span;
			it->t = max(estimate + this->info.at_ST(it->id), span);//预测离开时间
			it->t = max(it->t, td + arrival);//考虑延迟时间的离开时间
			it->e = it->t - estimate;//等待时间
		}
		else {//计算卡车路线
			if (it != this->solu.begin()) {
				if (!it->droneRoute.empty()) it->numDrone = 0;
				else it->numDrone = prev(it)->numDrone;
				double estimate = prev(it)->t + this->info.at_TT(prev(it)->id, it->id);//到达时间
				it->t = estimate + this->info.at_ST(it->id);//服务完成时间
				it->t = max(it->t, td + arrival);
				it->e = it->t - estimate;
			}
			else {//起点
				if (!it->droneRoute.empty()) it->numDrone = 0;
				else it->numDrone = 1;
				it->t = 0;
				it->e = 0;
			}
		}

		//计算无人机路线
		if (!it->droneRoute.empty()) {
			start = it;//定位无人机起点
			s = true;

			it->droneRoute.front().t = it->t;
			double demand = 0;
			for (auto it1 = next(it->droneRoute.begin()); it1 != it->droneRoute.end(); it1++) {
				it1->e = it1->id != it->droneRoute.back().id ? this->info.at_SD(it1->id) : 0;
				it1->t = prev(it1)->t + this->info.at_TD(prev(it1)->id, it1->id) + it1->e;
				it1->q = prev(it1)->q - demand;
				demand = this->info.at_demand(it1->id);
				if (it1->id != it->end->id) td = max(td, this->info.at_custDelay(it1->id));//延迟时间计算
			}
		}
	}
	if (s) this->makespan = max(this->solu.back().t, start->droneRoute.back().t);
	else this->makespan = this->solu.back().t;
}

//延迟客户约束
void Solution::constraintCD() {
	int num = 0;//记录目标数量
	double d_sum = 0.0;//记录目标重量
	for (auto it = this->solu.begin(); it != this->solu.end(); it++) {
		if (it->id == 0) continue;
		if (this->info.in_CD(it->id)) num++;
		if (this->info.at_custDelay(it->id) != 0) d_sum += this->info.at_demand(it->id);
		if (!it->droneRoute.empty()) {
			auto last = prev(it->droneRoute.end());
			for (auto it1 = next(it->droneRoute.begin()); it1 != last; it1++) {
				if (this->info.at_custDelay(it1->id) != 0) d_sum += this->info.at_demand(it1->id);
			}
		}
		int obj = static_cast<int>(ceil(d_sum / QDr));
		if (num < obj) {//违反约束的情况
			vector<pair<double, list<TruckCust>::iterator>> custs;
			for (auto it1 = next(it); it1 != this->solu.end(); it1++) {
				if (this->info.in_CD(it1->id) && it1->droneRoute.empty() 
					&& it1->drone_end == nullptr && this->info.at_custDelay(it1->id) == 0) {
					auto pre = prev(it1);
					auto nex = next(it1);
					double temp = this->info.at_TT(pre->id, it1->id) + this->info.at_TT(it1->id, nex->id) - this->info.at_TT(pre->id, nex->id);
					custs.emplace_back(temp, it1);
				}
			}

			std::sort(custs.begin(), custs.end(), [](pair<double, list<TruckCust>::iterator>& a, pair<double, list<TruckCust>::iterator>& b)->bool
				{ return a.first > b.first; });

			for (int i = 0; i < obj - num; i++) {
				int c = custs[i].second->id;
				this->custInsert.push_back(c);
				this->solu.erase(custs[i].second);
				num++;
			}

			this->greedyConstruct(this->custInsert, it);
		}
	}
}


void Solution::calMakespan(map<int, DroneResupply>& res, ResupplyInfo& result) {
	list<TruckCust>::iterator start;
	bool s = false;//表示是否存在start
	
	double rd = 0;//前一个补货造成的延迟时间

	for (auto it = next(this->solu.begin()); it != this->solu.end(); it++) {
		//计算补货无人机造成的延迟时间

		//添加补货方案
		auto resupply = res.find(it->id);
		if (resupply != res.end()) {
			it->resupply = resupply->second;
			rd = result.get_dt(it->id);
			it->e += rd;
			this->tw_resupply[it->id] = { it->t - it->e,it->t + rd };//记录time window
		}
		else it->resupply.clear();

		if (rd == 0) continue;
		it->t += rd;

		//计算无人机路线
		if (!it->droneRoute.empty()) {
			start = it;//定位无人机起点
			s = true;

			it->droneRoute.front().t = it->t;
			for (auto it1 = next(it->droneRoute.begin()); it1 != it->droneRoute.end(); it1++) {
				it1->t = prev(it1)->t + this->info.at_TD(prev(it1)->id, it1->id) + it1->e;
			}
		}

		if (s && it == start->end) {//计算降落点
			start->droneRoute.back().t = max(it->t - it->e, start->droneRoute.back().t);
		}
	}

	if (s) this->makespan = max(this->solu.back().t, start->droneRoute.back().t);
	else this->makespan = this->solu.back().t;
}

void Solution::solveResupply() {
	this->tw_resupply.clear();
	this->seq = vector<int>(this->info.get_num() + 1, 0);

	ResupplyInfo resupply(this->solu, this->info);
	map<int, DroneResupply> res;
	resupply.solve(res);
	this->calMakespan(res, resupply);

	this->tw_resupply.clear();
	for (auto& [key, value] : res) {
		for (auto& j : value.goods) {
			this->seq[j] = key;
		}
	}
}

Location Solution::find(int obj) {
	for (auto it = this->solu.begin(); it != this->solu.end(); it++) {
		if (it->id == obj) return Location('t', it);

		if (!it->droneRoute.empty()) {
			for (auto it1 = it->droneRoute.begin(); it1 != it->droneRoute.end(); it1++) {
				if (it1->id == obj) return Location('d', it, it1);
			}
		}
	}
	return Location('t', this->solu.end());
}

/*
查找目标并移除，同时分类移除的目标
*/
bool Solution::remove(Location& obj, vector<pair<int, double>>& cust, int& num) {
	if (obj.t == this->solu.end()) return false;

	if (obj.kind == 't') {
		if (!obj.t->droneRoute.empty()) {
			for (auto& i : obj.t->droneRoute) {
				if (i.id != obj.t->droneRoute.back().id) {
					if (this->info.at_custDelay(i.id) != 0 || this->tw_resupply.find(i.id) != this->tw_resupply.end())
						this->custInsert_delay.push_back(i.id);
					else 
						this->custInsert.push_back(i.id);
					cust[i.id] = { 0,0 };
					num--;
				}
			}

			obj.t->end->drone_end = nullptr;
			this->solu.erase(obj.t);
		}
		else {
			if (this->info.at_custDelay(obj.t->id) != 0 || this->tw_resupply.find(obj.t->id) != this->tw_resupply.end())
				this->custInsert_delay.push_back(obj.t->id);
			else 
				this->custInsert.push_back(obj.t->id);
			cust[obj.t->id] = { 0,0 };
			this->solu.erase(obj.t);
			num--;
		}
	}
	else {
		if (this->info.at_custDelay(obj.d->id) != 0 || this->tw_resupply.find(obj.d->id) != this->tw_resupply.end()) 
			this->custInsert_delay.push_back(obj.d->id);
		else 
			this->custInsert.push_back(obj.d->id);
		cust[obj.d->id] = { 0,0 };

		obj.t->droneRoute.erase(obj.d);
		num--;
	}
	return true;
}

bool Solution::remove(Location& obj, vector<int>& cust, int& num) {
	if (obj.t == this->solu.end()) return false;

	if (obj.kind == 't') {
		if (!obj.t->droneRoute.empty()) {
			for (auto& i : obj.t->droneRoute) {
				if (i.id != obj.t->droneRoute.back().id) {
					if (this->info.at_custDelay(i.id) != 0 || this->tw_resupply.find(i.id) != this->tw_resupply.end()) 
						this->custInsert_delay.push_back(i.id);
					else 
						this->custInsert.push_back(i.id);
					cust[i.id] = 0;
					num--;
				}
			}

			obj.t->end->drone_end = nullptr;
			this->solu.erase(obj.t);
		}
		else {
			if (this->info.at_custDelay(obj.t->id) != 0 || this->tw_resupply.find(obj.t->id) != this->tw_resupply.end()) 
				this->custInsert_delay.push_back(obj.t->id);
			else
				this->custInsert.push_back(obj.t->id);
			cust[obj.t->id] = 0;
			this->solu.erase(obj.t);
			num--;
		}
	}
	else {
		if (this->info.at_custDelay(obj.d->id) != 0 || this->tw_resupply.find(obj.d->id) != this->tw_resupply.end()) 
			this->custInsert_delay.push_back(obj.d->id);
		else
			this->custInsert.push_back(obj.d->id);
		cust[obj.d->id] = 0;

		obj.t->droneRoute.erase(obj.d);
		num--;
	}
	return true;
}

void Solution::randomRemove() {
	int num = this->info.get_num();

	vector<int> cust;
	int n = this->info.get_num();
	cust.push_back(0);
	for (int i = 1; i <= n; i++) {
		cust.push_back(i);
	}

	while (this->custInsert.size() + this->custInsert_delay.size() < this->numRemove) {
		uniform_int_distribution<int> q(0, num - 1);
		int obj = q(gen);
		int a = 0;
		for (size_t i = 1; i <= n; i++) {
			if (cust[i] == 0) continue;
			if (obj == a) {
				auto loc = find(cust[i]);
				this->remove(loc, cust, num);
				break;
			}
			a++;
		}
	}

	this->calMakespan(this->solu.begin());
}

void Solution::shawRemove() {//距离作为相似性计算
	uniform_int_distribution<int> prob(1, this->info.get_num());
	int obj = prob(gen);

	//计算相似性
	int num = this->info.get_num();
	vector<pair<int, double>> cust;
	cust.emplace_back(0, 0.0);
	int n = this->info.get_num();
	for (int i = 1; i <= n; i++) {
		cust.emplace_back(i, this->info.at_dis(obj, i));
	}

	//移除
	auto loc = find(obj);
	this->remove(loc, cust, num);

	sort(cust.begin(), cust.end(), [](pair<int, double>& a, pair<int, double>& b)->bool
		{ return a.second < b.second; });

	while (this->custInsert.size() + this->custInsert_delay.size() < this->numRemove) {
		uniform_real_distribution<float> y(0, 1);
		int q = static_cast<int>(floor(pow(y(gen), P) * num));
		int a = 0;
		for (size_t i = 1; i <= n; i++) {
			if (cust[i].first == 0) continue;
			if (q == a) {
				loc = find(cust[i].first);
				this->remove(loc, cust, num);
				break;
			}
			a++;
		}
	}

	this->calMakespan(this->solu.begin());
}

void Solution::worstRemove() {//上个客户离开时间至下个客户到达时间，作为cost
	int n = this->info.get_num();
	vector<pair<int, double>> cost(static_cast<size_t>(n) + 1, { 0,0 });
	auto end = prev(this->solu.end());
	for (auto it = next(this->solu.begin()); it != end; it++) {
		auto nex = next(it);
		double temp = nex->t - nex->e - prev(it)->t;
		cost[it->id] = { it->id,temp };
		if (!it->droneRoute.empty()) {
			auto end1 = prev(it->droneRoute.end());
			for (auto it1 = next(it->droneRoute.begin()); it1 != end1; it1++) {
				auto nex1 = next(it1);
				temp = nex1->t - nex1->e - prev(it1)->t;
				cost[it1->id] = { it1->id,temp };
			}
		}
	}

	sort(cost.begin(), cost.end(), [](pair<int, double>& a, pair<int, double>& b)->bool {
		return a.second > b.second; });

	int num = n;

	while (this->custInsert.size() + this->custInsert_delay.size() < this->numRemove) {
		uniform_real_distribution<float> y(0, 1);
		int q = static_cast<int>(floor(pow(y(gen), P) * num));
		int a = 0;
		for (size_t i = 1; i <= n; i++) {
			if (cost[i].first == 0) continue;
			if (q == a) {
				auto loc = find(cost[i].first);
				this->remove(loc, cost, num);
				break;
			}
			a++;
		}
	}

	this->calMakespan(this->solu.begin());
}

void Solution::seperate(list<int>& a) {
	for (auto it = this->custInsert_delay.begin(); it != this->custInsert_delay.end(); it++) {
		if (this->tw_resupply.find(*it) != this->tw_resupply.end()) {
			a.push_back(*it);
			it = this->custInsert_delay.erase(it);
		}
	}
}

//插入后半段
void Solution::calBestInsertionLS(list<int>::iterator obj, Insertion& best, list<TruckCust>::iterator solu_start, double max_noise) {
	uniform_real_distribution<double> noise(-max_noise, max_noise);

	auto td = this->info.at_custDelay(*obj);//延迟时间

	for (auto it = next(solu_start); it != this->solu.end(); it++) {
		auto pre = prev(it);//记录前一个卡车客户

		//插入卡车的情况
		double delta = this->info.at_ST(*obj) + this->info.at_TT(pre->id, *obj) + this->info.at_TT(*obj, it->id) - this->info.at_TT(pre->id, it->id);
		double increment = max(td + this->info.at_TD(0, *obj) - (pre->t + this->info.at_TT(pre->id, *obj)), 0.0);//时间窗延迟时间
		delta += increment;
		if (pre->numDrone == 0) {//插入在无人机跨度的卡车路线
			auto start = pre;
			for (; start != this->solu.begin(); start--) {
				if (!start->droneRoute.empty()) break;
			}
			auto end = start->end;
			if (delta + end->t - end->e - start->droneRoute.front().t > TLimit) delta = static_cast<double>(INT_MAX);
		}
		delta += delta * noise(gen);
		if (delta < best.delta)	best = Insertion(delta, 't', obj, it);


		//无人机的情况
		if (this->info.at_demand(*obj) > QDd) continue;

		increment = max(td + this->info.at_TD(0, pre->id) - pre->t, 0.0);//时间窗延迟时间

		if (pre->numDrone) {//新增无人机的情况
			increment = max(td + this->info.at_TD(0, it->id) - (pre->t + this->info.at_TT(pre->id, *obj)), 0.0);
			for (auto it1 = it; it1 != this->solu.end(); it1++) {
				double span = this->info.at_TD(pre->id, *obj) + this->info.at_TD(*obj, it1->id) + this->info.at_SD(*obj);
				if (span > TLimit) continue;
				if (it1->t - it1->e - pre->t > TLimit) break;

				delta = (span + pre->t - pre->e) - it1->t;//为多个客户的情况提供可能
				delta = (delta + increment) * (1 + noise(gen));
				if (delta < best.delta) best = Insertion(delta, 'n', obj, pre, it1);

			}
		}
		else {//插入无人机路线的情况
			for (auto it1 = next(pre->droneRoute.begin()); it1 != pre->droneRoute.end(); it1++) {
				auto pre1 = prev(it1);
				delta = this->info.at_SD(*obj) + this->info.at_TD(pre1->id, *obj) + this->info.at_TD(*obj, it1->id) - this->info.at_TD(pre1->id, it1->id);
				double span = pre->droneRoute.back().t - pre->droneRoute.back().e + delta;
				if (span - pre->droneRoute.front().t > TLimit) continue;

				delta = (span - pre->end->t) * (1 + noise(gen));
				bool up = false;
				if (delta > 0) {//尝试修改降落点
					auto end = prev(pre->droneRoute.begin(), 2);
					for (auto it2 = next(pre->end); it2 != this->solu.end(); it2++) {
						if (it2->t - it2->e - pre->t > TLimit) break;
						double new_span = span - this->info.at_TD(end->id, pre->end->id) + this->info.at_TD(end->id, it2->id);
						if (new_span < TLimit && new_span < it2->t) {
							delta = increment;
							if (delta < best.delta) {
								best = Insertion(delta, 'o', obj, pre, it1, it2);
								up = true;
							}
							break;
						}
					}
				}
				delta += increment;
				if (!up && delta < best.delta) best = Insertion(delta, 'd', obj, pre, it1);
			}
		}
	}
}


void Solution::greedyInsert(list<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		for (auto it = a.begin(); it != a.end(); it++) {
			if (this->seq[*it] == 0) this->calBestInsertionLS(it, best, this->solu.begin(), Noise);
			else {
				auto start = this->find(this->seq[*it]);
				if (start.t == this->solu.end()) this->calBestInsertionLS(it, best, this->solu.begin(), Noise);
				else this->calBestInsertionLS(it, best, start.t, Noise);
			}
		}

		this->insert(best);
		a.erase(best.id);//更新待插入集合
	}
};

void Solution::randomGreedyInsert(list<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		//随机选择客户
		uniform_int_distribution<int> prob(0, static_cast<int>(a.size()) - 1);
		auto it = a.begin();
		advance(it, prob(gen));

		if (this->seq[*it] == 0) this->calBestInsertionLS(it, best, this->solu.begin());
		else {
			auto start = this->find(this->seq[*it]);
			if (start.t == this->solu.end()) this->calBestInsertionLS(it, best, this->solu.begin());
			else this->calBestInsertionLS(it, best, start.t);
		}

		this->insert(best);
		a.erase(best.id);//更新待插入集合
	}
}

void Solution::updatePoint(int index, double delta) {
	this->pointLS[index] += delta;
	this->numLS[index] += 1;
}

void Solution::updateWeight() {
	for (size_t i = 0; i < this->weightLS.size(); i++) {
		if (this->numLS[i] != 0) {
			this->weightLS[i] = this->weightLS[i] * (1 - R) + R * this->pointLS[i] / this->numLS[i];
		}
	}
}

void Solution::cout() {

}