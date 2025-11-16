#pragma once
#include "Algorithm.h"
#include <list>
#include <unordered_set>
#include "CustInfo.h"
#include <functional>
#include <map>


class ResupplyInfo;


struct DroneResupply {//补货信息
	double time;//补货出发时间
	unordered_set<int> goods;//补给的货物集合
	DroneResupply(double a, unordered_set<int>& b) :time(a), goods(b) {};
	DroneResupply() { this->time = -1; }
	void clear() {
		this->time = 0;
		this->goods.clear();
	}
};

struct DroneCust {//无人机配送的客户信息
	int id = -1;//客户id
	double t = 0;//离开时间
	double e = 0;//等待时间
	double q = 0;//到达该点前的载重量
	DroneCust(int a) :id(a) {};
};

struct TruckCust {//卡车配送的客户信息
	int id = -1;//客户id
	int numDrone = 1;//当前点可用的无人机数量
	double t = 0;//离开时间
	double e = 0;//等待时间
	list<DroneCust> droneRoute;//以该点为起点的无人机路线
	DroneCust* drone_end = nullptr;//无人机路线的降落点(drone)
	list<TruckCust>::iterator end;//无人机路线的降落点(truck)
	DroneResupply resupply;//该点的补给信息
	TruckCust(int a) :id(a) {};
};

struct Insertion//插入位置的信息
{
	list<int>::iterator id;
	double delta = 0;
	char kind = 'a';//t,d,n分别代表卡车，无人机，新增无人机, o代表改变降落点，x插入卡车，修改降落点
	list<TruckCust>::iterator loc_T;//插入的卡车位置
	list<TruckCust>::iterator start;//无人机起飞点
	list<DroneCust>::iterator loc_D;//插入的无人机位置
	list<TruckCust>::iterator land;//降落位置
	Insertion(double delta) : delta(delta) {};

	//插入卡车
	Insertion(double delta, char kind, list<int>::iterator id, list<TruckCust>::iterator loc_T) :
		delta(delta), kind(kind), id(id), loc_T(loc_T) {};

	//插入卡车，修改降落点
	Insertion(double delta, char kind, list<int>::iterator id, list<TruckCust>::iterator loc_T, list<TruckCust>::iterator start, list<TruckCust>::iterator land) :
		delta(delta), kind(kind), id(id), loc_T(loc_T), start(start), land(land) {};

	//新增无人机
	Insertion(double delta, char kind, list<int>::iterator id, list<TruckCust>::iterator loc_T, list<TruckCust>::iterator land) :
		delta(delta), kind(kind), id(id), loc_T(loc_T), land(land) {};

	//插入无人机
	Insertion(double delta, char kind, list<int>::iterator id, list<TruckCust>::iterator loc_T, list<DroneCust>::iterator loc_D) :
		delta(delta), kind(kind), id(id), loc_D(loc_D), loc_T(loc_T) {};

	//改变降落点
	Insertion(double delta, char kind, list<int>::iterator id, list<TruckCust>::iterator loc_T, list<DroneCust>::iterator loc_D, list<TruckCust>::iterator land) :
		delta(delta), kind(kind), id(id), loc_D(loc_D), land(land), loc_T(loc_T) {};
};

struct Location
{
	char kind;
	list<TruckCust>::iterator t;
	list<DroneCust>::iterator d;
	Location(char kind_, list<TruckCust>::iterator t_) :t(t_), kind(kind_) {};
	Location(char kind_, list<TruckCust>::iterator t_, list<DroneCust>::iterator d_) :t(t_), d(d_), kind(kind_) {};
};


//一个解
class Solution
{//一个解及其操作
public:
	/*
	* 输入构造启发式的权重，获得配送方案
	*/
	Solution() : weightLS(Solution::algorithm.get_LS()) {
		//初始化待插入的两类客户
		int n = info.get_num();
		for (int i = 1; i <= n; i++) {
			if (info.at_custDelay(i) == 0) this->custInsert.push_back(i);
			else this->custInsert_delay.push_back(i);
		}

		//初始化解
		initial_solu();
		this->construct_index = this->intinialize(algorithm.get_Construct());//构造初始解
		this->numRemove = static_cast<int>(ceil(info.get_num() * PercentageRemove));//定义移除的个数
	}

	Solution(char b) : weightLS(Solution::algorithm.get_LS()) {
		this->makespan = static_cast<double>(INT_MAX);
	}

	Solution& operator=(const Solution& other) {
		if (this != &other) { // 将 data 绑定到其他对象的引用
			this->numRemove = other.numRemove;
			this->solu = other.solu;
			this->makespan = other.makespan;
			this->pointLS = other.pointLS;
			this->numLS = other.numLS;
			this->custInsert = other.custInsert;
			this->custInsert_delay = other.custInsert_delay;
			this->tw_resupply = other.tw_resupply;
			this->seq = other.seq;
		}
		return *this;
	}

	static CustInfo info;//客户信息类的引用
	static Algorithm algorithm;//权重信息
	int construct_index;

	int localSearch();//局部搜索，适用于第一阶段
	void solveResupply();//补货方案求解

	double get_makespan() const {
		return this->makespan;
	}

	void updatePoint(int, double);
	void updateWeight();

	void cout();

private:
	int numRemove;//移除个数
	list<TruckCust> solu;//解
	double makespan;//目标函数
	vector<double>& weightLS;//局部搜索算子的权重
	vector<double> pointLS = vector<double>(6, 0);//局部搜索算子的得分
	vector<int> numLS = vector<int>(6, 0);//局部搜索算子的使用次数

	list<int> custInsert;//待插入客户
	list<int> custInsert_delay;//待插入客户(延迟)

	map<int, pair<double, double>> tw_resupply;//补货方案造成的时间窗约束
	vector<int> seq;//补货方案造成的顺序约束//wu:0

	int intinialize(const vector<double>&);//初始化解

	//贪婪构造
	void greedyConstructNoise() {
		this->greedyConstruct(this->custInsert);//先构造非延迟客户
		this->greedyConstruct(this->custInsert_delay);//再构造延迟客户
		this->constraintCD();//约束检查
	};

	//随机贪婪构造
	void randomGreedyConstruct() {
		this->randomGreedyConstruct(this->custInsert);//先构造非延迟客户
		this->randomGreedyConstruct(this->custInsert_delay);//再构造延迟客户
		this->constraintCD();//约束检查
	};

	//待聚类的随机贪婪构造
	void greedyConstructNoiseCluster();

	void greedyConstruct(list<int>&);
	void greedyConstruct(list<int>&, list<TruckCust>::iterator);
	void randomGreedyConstruct(list<int>&);
	void calBestInsertion(list<int>::iterator, Insertion&, list<TruckCust>::iterator solu_end, double max_noise = 0);
	void calBestInsertionLS(list<int>::iterator, Insertion&, list<TruckCust>::iterator solu_start, double max_noise = 0);
	void insert(Insertion&);
	void calMakespan(list<TruckCust>::iterator);
	void constraintCD();

	//局部搜索
	void randomRemove();
	void shawRemove();
	void worstRemove();

	void greedyInsert(list<int>&);
	void randomGreedyInsert(list<int>&);

	void seperate(list<int>&);

	Location find(int);
	bool remove(Location&, vector<pair<int, double>>&, int&);
	bool remove(Location&, vector<int>&, int&);

	//补货方案求解
	void calMakespan(map<int, DroneResupply>&, ResupplyInfo&);

	
	//构造仓库
	void initial_solu() {
		//初始化解
		this->solu.clear();
		this->solu.push_back(TruckCust(0));
		this->solu.push_back(TruckCust(0));
	}
};
