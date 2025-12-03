#include "Solution.h"
#include "HyperParameters.h"
#include <numeric>
#include <algorithm>
#include "ResupplyInfo.h"

//初始化Solution的静态变量
CustInfo Solution::info = move(CustInfo(address));//客户信息
Algorithm Solution::algorithm = move(Algorithm(Solution::info));//权重信息
int Solution::numRemove = static_cast<int>(ceil(Solution::info.get_num() * PercentageRemove));//定义移除的个数

vector<double> Solution::weightLS;//局部搜索算子的权重
double Solution::relaxation;//松弛系数
vector<double> Solution::pointLS;//局部搜索算子的得分
vector<int> Solution::numLS;//局部搜索算子的使用次数

/*
选择构造启发式
*/
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

/*
选择局部搜索算子
*/
int Solution::localSearch() {
	discrete_distribution<> prob(this->weightLS.begin(), weightLS.end());
	int t = prob(gen);

	vector<function<void()>> operatorLS = {
			[this]()->void {
				this->randomRemove();
				this->greedyInsert(this->custInsert);
				this->greedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->randomRemove();
				this->randomGreedyInsert(this->custInsert);
				this->randomGreedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->shawRemove();
				this->greedyInsert(this->custInsert);
				this->greedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->shawRemove();
				this->randomGreedyInsert(this->custInsert);
				this->randomGreedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->worstRemove();
				this->greedyInsert(this->custInsert);
				this->greedyInsert(this->custInsert_delay);
			},
			[this]()->void {
				this->worstRemove();
				this->randomGreedyInsert(this->custInsert);
				this->randomGreedyInsert(this->custInsert_delay);
			}
	};//局部搜索函数的数组

	operatorLS[t]();
	this->delivery_makespan = makespan;
	return t;
}


/*
* 贪婪构造，计算到solu_end的插入位置
*/
void Solution::greedyConstruct(unordered_set<int>& a, list<TruckCust>::iterator solu_end) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		for (auto it : a) {
			this->calBestInsertion(it, best, this->solu.begin(), solu_end, NOISE);
		}

		this->insert(best);
		a.erase(best.id);//更新待插入集合
		int yy = call();
		if (yy != 100) {
			call();
			throw runtime_error("111");
		}
	}
}

/*
* 随机贪婪构造
*/
void Solution::randomGreedyConstruct(unordered_set<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		//随机选择客户
		uniform_int_distribution<int> prob(0, static_cast<int>(a.size()) - 1);

		auto it = a.begin();
		advance(it, prob(gen));

		this->calBestInsertion(*it, best, this->solu.begin(), this->solu.end());//最优插入客户;

		this->insert(best);
		a.erase(best.id);//更新待插入集合
		int yy = call();
		if (yy != 100) {
			call();
			throw runtime_error("111");
		}
	}
}

/*
贪婪构造_簇情况下
*/
void Solution::greedyConstructNoiseCluster() {
	list<list<TruckCust>> solu_part;//部分解集合

	for (auto& cluster : Solution::info.at_clusters()) {
		this->custInsert.clear();
		this->custInsert_delay.clear();
		for (auto customer : cluster) {
			if (Solution::info.at_custDelay(customer) == -1) this->custInsert.insert(customer);
			else this->custInsert_delay.insert(customer);
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
				double temp = Solution::info.at_TT(a, b) - Solution::info.at_TT(a, 0) - Solution::info.at_TT(0, b);
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
void Solution::calBestInsertion(int obj, Insertion& best, list<TruckCust>::iterator solu_start,
	list<TruckCust>::iterator solu_end, double max_noise) {
	
	uniform_real_distribution<double> noise(-max_noise, max_noise);

	auto td = Solution::info.at_custDelay(obj);//延迟时间

	for (auto it = next(solu_start); it != solu_end; ++it) {
		auto pre = prev(it);//记录it的前一个卡车客户,插入it之前

		//插入卡车的情况
		double delta = Solution::info.at_ST(obj) + Solution::info.at_TT(pre->id, obj) + Solution::info.at_TT(obj, it->id) - Solution::info.at_TT(pre->id, it->id);//插入代价
		double increment = max(td + Solution::info.at_TD(0, obj) - (pre->t + Solution::info.at_TT(pre->id, obj)), 0.0);//时间窗造成的延迟时间，
		delta += increment;//综合插入代价
		if (pre->numDrone == 0) {//插入在无人机跨度的卡车路线
			auto start = pre;
			for (; start != this->solu.begin(); --start) {
				if (!start->droneRoute.empty()) break;
			}

			//修改无人机路线的降落点
			auto end = next(start);
			//记录最优降落点
			double best_delta = static_cast<double>(INT_MAX);
			list<TruckCust>::iterator best_end;
			bool is_find = false;

			for (; end != it; ++end) {//降落点it前的
				double arrival_end = end->t - end->e;//卡车到达end的时间
				auto last = prev(start->droneRoute.end(), 2);//无人机服务的最后一个客户
				double span_drone = last->t + Solution::info.at_TD(last->id, end->id);//无人机到达降落点的时间
				if (arrival_end - start->t <= TLimit && span_drone - start->t <= TLimit) {
					double delta_cost = span_drone - end->t;//等待无人机的时间
					if (abs(delta_cost) < abs(best_delta)) {
						best_delta = delta_cost;
						best_end = end;
						is_find = true;
					}
				}
			}

			for (; end != next(start->end); ++end) {//降落点it及之后
				double arrival_end = delta + end->t - end->e;//卡车到达end的时间
				auto last = prev(start->droneRoute.end(), 2);//无人机服务的最后一个客户
				double span_drone = last->t + Solution::info.at_TD(last->id, end->id);//无人机到达降落点的时间
				if (arrival_end - start->t <= TLimit && span_drone - start->t <= TLimit) {
					double delta_cost = span_drone - end->t;//等待无人机的时间
					if (abs(delta_cost) < abs(best_delta)) {
						best_delta = delta_cost;
						best_end = end;
						is_find = true;
					}
				}
			}

			delta = is_find ? delta * (1 + noise(gen)) : static_cast<double>(INT_MAX);
			if (delta < best.delta)	best = Insertion(delta, 'x', obj, it, start, best_end);
		}
		else {
			delta += delta * noise(gen);
			if (delta < best.delta)	best = Insertion(delta, 't', obj, it);
		}


		//无人机的情况
		if (Solution::info.at_demand(obj) > QDd) continue;

		if (pre->numDrone) {//新增无人机的情况，以pre为起点
			increment = max(td + Solution::info.at_TD(0, it->id) - (pre->t + Solution::info.at_TT(pre->id, obj)), 0.0);//时间窗造成的延迟时间
			for (auto it1 = it; it1 != this->solu.end(); ++it1) {//寻找降落点
				if (it1->id == pre->id) break;
				double span = Solution::info.at_TD(pre->id, obj) + Solution::info.at_TD(obj, it1->id) + Solution::info.at_SD(obj);//无人机需要飞行的最小时长
				if (span <= TLimit) {//满足无人机飞行限制
					if (it1->t - it1->e - pre->t > TLimit) break;

					delta = (span + pre->t - pre->e) - it1->t;//为多个客户的情况提供可能
					delta = (delta + increment) * (1 + noise(gen));
					if (delta < best.delta) best = Insertion(delta, 'n', obj, pre, it1);
				}
				if (!it1->droneRoute.empty()) break;
			}
		}
		else if (!pre->droneRoute.empty() && pre->droneRoute.front().q + Solution::info.at_demand(obj) <= QDd) {//插入无人机路线的情况，pre作起点
			for (auto it1 = next(pre->droneRoute.begin()); it1 != pre->droneRoute.end(); ++it1) {
				increment = max(td + Solution::info.at_TD(0, pre->id) - pre->t, 0.0);//时间窗造成的延迟时间

				auto pre1 = prev(it1);//插入位置的前一个客户
				delta = Solution::info.at_SD(obj) + Solution::info.at_TD(pre1->id, obj) + Solution::info.at_TD(obj, it1->id) - Solution::info.at_TD(pre1->id, it1->id);
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
						double new_span = span - Solution::info.at_TD(end->id, pre->end->id) + Solution::info.at_TD(end->id, it2->id);
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
		this->solu.insert(best.loc_T, TruckCust(best.id));
		this->calMakespan(prev(best.loc_T));//定位到插入点
		break;
	case 'n'://新增无人机路线
		best.loc_T->droneRoute.emplace_back(best.loc_T->id);
		best.loc_T->droneRoute.emplace_back(best.id);
		best.loc_T->droneRoute.emplace_back(best.land->id);
		best.loc_T->end = best.land;
		best.land->drone_end = &(best.loc_T->droneRoute.back());
		best.loc_T->droneRoute.front().q += Solution::info.at_demand(best.id);
		this->calMakespan(best.loc_T);//定位到插入点
		break;
	case 'd'://插入无人机路线
		best.loc_T->droneRoute.insert(best.loc_D, DroneCust(best.id));
		best.loc_T->droneRoute.front().q += Solution::info.at_demand(best.id);
		this->calMakespan(best.loc_T);//定位到插入点
		break;
	case 'o'://改变降落点
		best.loc_T->droneRoute.insert(best.loc_D, DroneCust(best.id));
		best.loc_T->droneRoute.back().id = best.land->id;
		best.loc_T->end->drone_end = nullptr;
		best.loc_T->end = best.land;
		best.land->drone_end = &(best.loc_T->droneRoute.back());
		best.loc_T->droneRoute.front().q += Solution::info.at_demand(best.id);
		this->calMakespan(best.loc_T);//定位到插入点
		break;
	case 'x'://插入卡车，修改降落点
		this->solu.insert(best.loc_T, TruckCust(best.id));
		best.start->droneRoute.back().id = best.land->id;
		best.start->end->drone_end = nullptr;
		best.start->end = best.land;
		best.land->drone_end = &(best.start->droneRoute.back());
		best.land->numDrone = 1;
		this->calMakespan(best.start);//定位到插入点
		break;
	}
}


//延迟客户约束
void Solution::constraintCD() {
	int num = 0;//记录目标数量
	double d_sum = 0.0;//记录目标重量
	for (auto it = this->solu.begin(); it != this->solu.end(); it++) {
		if (it->id == 0) continue;
		if (Solution::info.in_CD(it->id)) num++;
		if (Solution::info.at_custDelay(it->id) != 0) d_sum += Solution::info.at_demand(it->id);
		if (!it->droneRoute.empty()) {
			auto last = prev(it->droneRoute.end());
			for (auto it1 = next(it->droneRoute.begin()); it1 != last; it1++) {
				if (Solution::info.at_custDelay(it1->id) != 0) d_sum += Solution::info.at_demand(it1->id);
			}
		}
		int obj = static_cast<int>(ceil(d_sum / QDr));
		if (num < obj) {//违反约束的情况
			vector<pair<double, list<TruckCust>::iterator>> custs;
			for (auto it1 = next(it); it1 != this->solu.end(); it1++) {
				if (Solution::info.in_CD(it1->id) && it1->droneRoute.empty()
					&& it1->drone_end == nullptr && Solution::info.at_custDelay(it1->id) == 0) {//需要至少满足非起降点的要求
					auto pre = prev(it1);
					auto nex = next(it1);
					double temp = Solution::info.at_TT(pre->id, it1->id) + Solution::info.at_TT(it1->id, nex->id) - Solution::info.at_TT(pre->id, nex->id);
					custs.emplace_back(temp, it1);
				}
			}

			std::sort(custs.begin(), custs.end(), [](pair<double, list<TruckCust>::iterator>& a, pair<double, list<TruckCust>::iterator>& b)->bool
				{ return a.first > b.first; });

			for (int i = 0; i < obj - num; i++) {
				int c = custs[i].second->id;
				this->custInsert.insert(c);
				this->solu.erase(custs[i].second);
				num++;
			}

			this->calMakespan(this->solu.begin());//计算移除后的makespan
			this->greedyConstruct(this->custInsert, it);
		}
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
		double td = Solution::info.at_custDelay(it->id);//延迟时间
		double arrival = Solution::info.at_TD(0, it->id);//假设，补货无人机飞行时间

		if (!it->droneRoute.empty()) {//该点是起始点，计算延迟时间			
			for (auto it1 = next(it->droneRoute.begin()); it1 != prev(it->droneRoute.end()); it1++) {
				td = max(td, Solution::info.at_custDelay(it1->id));//延迟时间计算
			}
		}

		if (s && it == start->end) {//计算降落点
			if (!it->droneRoute.empty()) it->numDrone = 0;
			else it->numDrone = 1;
			auto last = prev(start->droneRoute.end(), 2);//无人机服务的最后一个客户
			double span = last->t + Solution::info.at_TD(last->id, start->droneRoute.back().id);
			double estimate = prev(it)->t + Solution::info.at_TT(prev(it)->id, it->id);//到达时间
			start->droneRoute.back().t = max(span, estimate);
			start->droneRoute.back().e = start->droneRoute.back().t - span;
			it->t = max(estimate + Solution::info.at_ST(it->id), span);//预测离开时间
			it->t = max(it->t, td + arrival);//考虑延迟时间的离开时间
			it->e = it->t - estimate;//等待时间
		}
		else {//计算卡车路线
			if (it != this->solu.begin()) {
				if (!it->droneRoute.empty()) it->numDrone = 0;
				else it->numDrone = prev(it)->numDrone;
				double estimate = prev(it)->t + Solution::info.at_TT(prev(it)->id, it->id);//到达时间
				it->t = estimate + Solution::info.at_ST(it->id);//服务完成时间
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
				it1->e = it1->id != it->droneRoute.back().id ? Solution::info.at_SD(it1->id) : 0;
				it1->t = prev(it1)->t + Solution::info.at_TD(prev(it1)->id, it1->id) + it1->e;
				it1->q = prev(it1)->q - demand;
				demand = Solution::info.at_demand(it1->id);
				if (it1->id != it->end->id) td = max(td, Solution::info.at_custDelay(it1->id));//延迟时间计算
			}
		}
	}
	if (s) this->makespan = max(this->solu.back().t, start->droneRoute.back().t);
	else this->makespan = this->solu.back().t;
}


void Solution::recalculate(list<TruckCust>::iterator it, double& rd) {
	if (it->droneRoute.back().t < it->end->t - it->end->e) {
		double a_s = it->t - it->e;//卡车到达start的时间
		double l_s = it->t;//卡车离开start的时间

		double best = static_cast<double>(INT_MAX);
		list<TruckCust>::iterator best_end;
		char type;

		for (auto end = prev(it->end); end != it; --end) {//降落点it前的
			double a_e = end->t - end->e;//卡车到达end的时间
			double l_e = end->t;//卡车离开end的时间
			auto last = prev(it->droneRoute.end(), 2);//无人机服务的最后一个客户
			double span_drone = last->t + Solution::info.at_TD(last->id, end->id);//无人机到达降落点的时间
			double duration = span_drone - it->droneRoute.front().t;//无人机飞行时长

			if (a_e - l_s <= TLimit && duration <= TLimit) {
				double in = min(duration - (a_e - l_s), 0.0);
				double ex = max(duration - (l_e - a_s), 0.0);
				double temp = duration - (a_e - l_s);
				if (in == 0 && ex == 0) {
					if (temp < best) {
						best = temp;
						best_end = end;
						type = 'b';
					}
				}
				else if (in == 0 && ex != 0) {
					temp *= 1.5;
					if (temp < best) {
						best = temp;
						best_end = end;
						type = 'w';
					}
				}
				else if (in != 0 && ex == 0) {
					temp *= -1.2;
					if (temp < best) {
						best = temp;
						best_end = end;
						type = 'm';
					}
				}
			}
		}

		//修改降落点
		double arrive_e = best_end->t - best_end->e;//卡车到达end的时间
		double leave_e = best_end->t;//卡车离开end的时间

		it->droneRoute.back().id = best_end->id;
		for (auto it1 = best_end; it1 != it->end; ++it) {//更新无人机的可用个数
			it1->numDrone = 1;
		}
		it->end->drone_end = nullptr;
		it->end = best_end;
		best_end->drone_end = &(it->droneRoute.back());

		if (type == 'b') {
			double delta = it->droneRoute.front().t;//起飞的时间变化量
			it->droneRoute.front().t = it->t - best * (l_s - a_s) / (l_s - a_s + leave_e - arrive_e);
			delta = it->droneRoute.front().t - delta;
			for (auto it1 = next(it->droneRoute.begin()); it1 != it->droneRoute.end(); it1++) {
				if (it1->id != best_end->id)
					it1->t += delta;
				else
					it1->t = prev(it1)->t + Solution::info.at_TD(prev(it1)->id, it1->id);
			}
		}
		else if (type == 'w') {
			double delta = it->droneRoute.front().t;//起飞的时间变化量
			it->droneRoute.front().t = it->t - it->e;
			delta = it->droneRoute.front().t - delta;

			double t_delta = 0;

			for (auto it1 = next(it->droneRoute.begin()); it1 != it->droneRoute.end(); it1++) {
				if (it1->id != best_end->id)
					it1->t += delta;
				else {
					it1->t = prev(it1)->t + Solution::info.at_TD(prev(it1)->id, it1->id);
					t_delta = it1->t - best_end->t;//增加降落点的等待时长
				}
			}

			for (auto it1 = best_end; it1 != it->end; ++it) {//更新无人机的可用个数
				it1->t += t_delta;
			}

			rd += t_delta;
		}
		else if (type == 'm') {
			double delta = it->droneRoute.front().t;//起飞的时间变化量
			it->droneRoute.front().t = it->t;
			delta = it->droneRoute.front().t - delta;
			for (auto it1 = next(it->droneRoute.begin()); it1 != it->droneRoute.end(); it1++) {
				if (it1->id != best_end->id)
					it1->t += delta;
				else
					it1->t = best_end->t - best_end->e;
			}
		}
	}
}

/*
* res表示补货方案，result表示补货方案的信息，将补货方案整合到解中
*/
void Solution::calMakespan_resupply(map<int, DroneResupply>& res, ResupplyInfo& result) {
	list<TruckCust>::iterator start;
	bool s = false;//表示是否存在start
	
	double rd = 0;//前一个补货造成的延迟时间

	for (auto it = this->solu.begin(); it != this->solu.end(); ++it) {
		//计算补货无人机造成的延迟时间

		if (!it->droneRoute.empty()) {//记录起飞点
			start = it;
			s = true;
		}

		//添加补货方案
		auto resupply_it = res.find(it->id);
		if (resupply_it != res.end()) {//增加补货方案，针对于补货点
			it->resupply = move(resupply_it->second);
			rd = result.get_dt(it->id);
			it->e += rd;//更新等待时间
		}
		else it->resupply.clear();

		if (rd == 0) continue;

		if (s && it == start->end) recalculate(start, rd);//修改降落点
		it->t += rd;//更新离开时间

		//计算无人机路线
		if (!it->droneRoute.empty()) {
			it->droneRoute.front().t += rd;
			for (auto it1 = next(it->droneRoute.begin()); it1 != it->droneRoute.end(); it1++) {
				it1->t += rd;
			}
		}
	}

	if (s) this->makespan = max(this->solu.back().t, start->droneRoute.back().t);
	else this->makespan = this->solu.back().t;
}

void Solution::solveResupply() {
	this->seq = vector<int>(Solution::info.get_num() + 1, 0);

	ResupplyInfo resupply(this->solu, Solution::info);
	map<int, DroneResupply> res;//补货方案的解
	resupply.solve(res);

	//赋值seq
	for (auto& [key, value] : res) {
		for (auto& j : value.goods) {
			this->seq[j] = key;
		}
		this->node_resupply.insert(key);
	}

	this->calMakespan_resupply(res, resupply);//重计算makespan
}

Location Solution::find(int obj) {
	for (auto it = this->solu.begin(); it != this->solu.end(); ++it) {
		if (it->id == obj) return Location('t', it);
	
		if (!it->droneRoute.empty()) {
			for (auto it1 = next(it->droneRoute.begin()); it1 != prev(it->droneRoute.end()); ++it1) {
				if (it1->id == obj) return Location('d', it, it1);
			}
		}
	}
	return Location('o', this->solu.end());
}


void Solution::handleRemoveCustomer(int id) {
	if (Solution::info.at_custDelay(id) != 0) {
		this->custInsert_delay.insert(id);
	}
	else {
		this->custInsert.insert(id);
	}
}


/*
* 移除obj,
*/
void Solution::remove(Location& obj) {
	if (obj.kind == 'o') throw runtime_error("not find!");

	if (obj.kind == 't') {
		if (!obj.t->droneRoute.empty()) {//移除作为起点的一系列客户
			for (auto& i : obj.t->droneRoute) {
				if (i.id != obj.t->droneRoute.back().id && i.id != obj.t->droneRoute.front().id) {
					handleRemoveCustomer(i.id);
				}
			}

			obj.t->end->drone_end = nullptr;
			obj.t->numDrone = 1;
		}
		if (obj.t->drone_end != nullptr) {//移除作为终点的一系列客户
			auto temp_it = prev(obj.t);
			for (; temp_it != this->solu.begin(); --temp_it) {
				if (!temp_it->droneRoute.empty()) break;
			}

			for (auto& i : temp_it->droneRoute) {
				if (i.id != temp_it->droneRoute.front().id && i.id != temp_it->droneRoute.back().id) {
					handleRemoveCustomer(i.id);
				}
			}

			temp_it->droneRoute.clear();
			temp_it->end = this->solu.begin();
			temp_it->numDrone = 1;
		}

		handleRemoveCustomer(obj.t->id);	
		this->solu.erase(obj.t);//同一删除
	}
	else {
		handleRemoveCustomer(obj.d->id);
		obj.t->droneRoute.erase(obj.d);
		if (obj.t->droneRoute.size() == 2) {
			obj.t->end->drone_end = nullptr;
			obj.t->numDrone = 1;
			obj.t->droneRoute.clear();
		}
	}
}

/*
* 随机移除，直到满足移除个数
*/
void Solution::randomRemove() {
	int num = Solution::info.get_num();//客户总数
	uniform_int_distribution<int> q(1, num);

	while (this->custInsert.size() + this->custInsert_delay.size() < this->numRemove) {
		int obj = q(gen);//目标客户
		if (this->custInsert.find(obj) != this->custInsert.end()
			|| this->custInsert_delay.find(obj) != this->custInsert_delay.end())
			continue;

		auto loc = find(obj);
		this->remove(loc);
		int yy = call();
		if (yy != 100) {
			call();
			throw runtime_error("111");
		}
	}

	this->calMakespan(this->solu.begin());
}

int Solution::call() {
	set<int> temp;
	for (auto& i : custInsert_delay) {
		temp.insert(i);
	}
	for (auto& i : custInsert) {
		temp.insert(i);
	}
	for (auto it = this->solu.begin(); it != solu.end(); ++it) {
		if (it->id != 0) temp.insert(it->id);
		if (!it->droneRoute.empty()) {
			for (auto& i : it->droneRoute) {
				if (i.id != it->droneRoute.front().id && i.id != it->droneRoute.back().id) temp.insert(i.id);
			}
		}
	}
	return temp.size();
}


/*
* 相似性移除
*/
void Solution::shawRemove() {//距离作为相似性计算
	int num = Solution::info.get_num();//客户总数
	uniform_int_distribution<int> prob(1, num);

	int obj = prob(gen);
	//移除
	auto loc = find(obj);
	this->remove(loc);

	uniform_real_distribution<float> y(0, 1);
	while (this->custInsert.size() + this->custInsert_delay.size() < this->numRemove) {
		prob = uniform_int_distribution<int>(0, this->custInsert.size() + this->custInsert_delay.size() - 1);
		obj = prob(gen);

		auto it = (obj < this->custInsert.size()) ? next(this->custInsert.begin(), obj) :
			next(this->custInsert_delay.begin(), obj - this->custInsert.size());//获取随机被移除的客户

		//计算相似性
		vector<pair<int, double>> shaw;
		for (int i = 1; i <= num; i++) {
			if (this->custInsert.find(i) != this->custInsert.end()
				|| this->custInsert_delay.find(i) != this->custInsert_delay.end())
				continue;
			shaw.emplace_back(i, Solution::info.at_dis(obj, i));
		}

		std::sort(shaw.begin(), shaw.end(), [](pair<int, double>& a, pair<int, double>& b)->bool
			{ return a.second < b.second; });

		int q = static_cast<int>(floor(pow(y(gen), P) * shaw.size()));
		loc = find(shaw[q].first);
		this->remove(loc);
		int yy = call();
		if (yy != 100) {
			throw runtime_error("111");
		}
	}

	this->calMakespan(this->solu.begin());
}

/*
* 最差移除，
*/
void Solution::worstRemove() {//上个客户离开时间至下个客户到达时间，作为cost
	int num = this->info.get_num();
	list<pair<int, double>> cost;
	auto end = prev(this->solu.end());
	for (auto it = next(this->solu.begin()); it != end; it++) {
		auto nex = next(it);
		double temp = nex->t - nex->e - prev(it)->t;
		cost.emplace_back(it->id, temp);
		if (!it->droneRoute.empty()) {
			auto end1 = prev(it->droneRoute.end());
			for (auto it1 = next(it->droneRoute.begin()); it1 != end1; it1++) {
				auto nex1 = next(it1);
				temp = nex1->t - nex1->e - prev(it1)->t;
				cost.emplace_back(it1->id, temp);
			}
		}
	}

	cost.sort([](pair<int, double>& a, pair<int, double>& b)->bool {
		return a.second > b.second; });

	uniform_real_distribution<float> y(0, 1);
	while (this->custInsert.size() + this->custInsert_delay.size() < this->numRemove) {
		int q = static_cast<int>(floor(pow(y(gen), P) * cost.size()));
		auto obj = next(cost.begin(), q);
		if (this->custInsert.find(obj->first) != this->custInsert.end()
			|| this->custInsert_delay.find(obj->first) != this->custInsert_delay.end()) {
			cost.erase(obj);
			continue;
		}
		auto loc = find(obj->first);
		this->remove(loc);
		cost.erase(obj);
		int yy = call();
		if (yy != 100) {
			throw runtime_error("111");
		}
	}

	this->calMakespan(this->solu.begin());
}

void Solution::seperate(unordered_set<int>& a) {
	for (auto it : a) {
		if (this->seq[it] != 0) {//补货后存在顺序约束
			a.insert(it);
			it = this->custInsert_delay.erase(it);
		}
	}
}


void Solution::greedyInsert(unordered_set<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		for (auto it : a) {
			if (this->seq[it] == 0) this->calBestInsertion(it, best, this->solu.begin(), this->solu.end(), NOISE);
			else {
				auto start = this->find(this->seq[it]);
				if (start.t == this->solu.end()) this->calBestInsertion(it, best, this->solu.begin(), this->solu.end(), NOISE);
				else this->calBestInsertion(it, best, start.t, this->solu.end(), NOISE);
			}
		}

		this->insert(best);
		a.erase(best.id);//更新待插入集合
	}
};

void Solution::randomGreedyInsert(unordered_set<int>& a) {
	while (!a.empty()) {
		double temp = static_cast<double> (INT_MAX);
		Insertion best(temp);//记录最优插入信息

		//随机选择客户
		uniform_int_distribution<int> prob(0, static_cast<int>(a.size()) - 1);
		auto it = a.begin();
		advance(it, prob(gen));

		if (this->seq[*it] == 0) this->calBestInsertion(*it, best, this->solu.begin(), this->solu.end());
		else {
			auto start = this->find(this->seq[*it]);
			if (start.t == this->solu.end()) this->calBestInsertion(*it, best, this->solu.begin(), this->solu.end());
			else this->calBestInsertion(*it, best, start.t, this->solu.end());
		}

		this->insert(best);
		a.erase(best.id);//更新待插入集合
	}
}

void Solution::updatePoint(int index, double delta) {
	Solution::pointLS[index] += delta;
	Solution::numLS[index] += 1;
}

void Solution::updateWeight() {
	for (size_t i = 0; i < Solution::weightLS.size(); i++) {
		if (Solution::numLS[i] != 0) {
			Solution::weightLS[i] = Solution::weightLS[i] * (1 - R) + R * Solution::pointLS[i] / Solution::numLS[i];
			Solution::pointLS[i] = 0;//重置
			Solution::numLS[i] = 0;//重置
		}
	}
}

void Solution::copy_iterator(const Solution& other) {
	for (auto it = solu.begin(); it != solu.end(); ++it) {
		if (it->droneRoute.empty()) continue;

		it->drone_end = &(it->droneRoute.back());
		for (auto it1 = next(it); it1 != solu.end(); ++it1) {
			if (it1->id == it->end->id) {
				it->end = it1;
				break;
			}
		}
	}
}

void Solution::cout() const {

}