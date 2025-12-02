#pragma once
#include "Algorithm.h"
#include <list>
#include <unordered_set>
#include "CustInfo.h"
#include <functional>
#include <unordered_map>
#include <map>

class ResupplyInfo;//前置声明


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
	int id;
	double delta = 0;
	char kind = 'a';//t,d,n分别代表卡车，无人机，新增无人机, o代表改变降落点，x插入卡车，修改降落点
	list<TruckCust>::iterator loc_T;//插入的卡车位置
	list<TruckCust>::iterator start;//无人机起飞点
	list<DroneCust>::iterator loc_D;//插入的无人机位置
	list<TruckCust>::iterator land;//降落位置
	Insertion(double delta) : delta(delta) {};

	//插入卡车
	Insertion(double delta, char kind, int id, list<TruckCust>::iterator loc_T) :
		delta(delta), kind(kind), id(id), loc_T(loc_T) {};

	//插入卡车，修改降落点
	Insertion(double delta, char kind, int id, list<TruckCust>::iterator loc_T, list<TruckCust>::iterator start, list<TruckCust>::iterator land) :
		delta(delta), kind(kind), id(id), loc_T(loc_T), start(start), land(land) {};

	//新增无人机
	Insertion(double delta, char kind, int id, list<TruckCust>::iterator loc_T, list<TruckCust>::iterator land) :
		delta(delta), kind(kind), id(id), loc_T(loc_T), land(land) {};

	//插入无人机
	Insertion(double delta, char kind, int id, list<TruckCust>::iterator loc_T, list<DroneCust>::iterator loc_D) :
		delta(delta), kind(kind), id(id), loc_D(loc_D), loc_T(loc_T) {};

	//改变降落点
	Insertion(double delta, char kind, int id, list<TruckCust>::iterator loc_T, list<DroneCust>::iterator loc_D, list<TruckCust>::iterator land) :
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
	Solution() {
		//初始化待插入的两类客户
		int n = info.get_num();
		for (int i = 1; i <= n; i++) {
			if (info.at_custDelay(i) == 0) this->custInsert.insert(i);
			else this->custInsert_delay.insert(i);
		}

		//初始化解
		initial_solu();
		this->construct_index = this->intinialize(algorithm.get_Construct());//构造初始解
		this->numRemove = static_cast<int>(ceil(info.get_num() * PercentageRemove));//定义移除的个数
		this->delivery_makespan = this->makespan;
	}

	Solution(char b) {
		this->makespan = static_cast<double>(INT_MAX);
	}

	/*
	* 生成新变量
	*/
	Solution(const Solution& other) : numRemove(other.numRemove),
		makespan(other.makespan),
		pointLS(other.pointLS),
		numLS(other.numLS),
		custInsert(other.custInsert),
		custInsert_delay(other.custInsert_delay),
		tw_resupply(other.tw_resupply),
		seq(other.seq),
		delivery_makespan(other.delivery_makespan) {
		
		copy_iterator(other);//深拷贝迭代器
	}

	/*
	* 赋值给一个已存在的变量
	*/
	Solution& operator=(const Solution& other) {
		if (this != &other) {
			this->numRemove = other.numRemove;
			this->makespan = other.makespan;
			this->pointLS = other.pointLS;
			this->numLS = other.numLS;
			this->custInsert = other.custInsert;
			this->custInsert_delay = other.custInsert_delay;
			this->tw_resupply = other.tw_resupply;
			this->seq = other.seq;
			this->delivery_makespan = other.delivery_makespan;

			copy_iterator(other);//深拷贝迭代器
		}
		return *this;
	}	

	/*
	* 用于move，传入右值引用
	*/
	Solution& operator=(Solution&& other) noexcept {
		if (this != &other) {
			this->numRemove = std::move(other.numRemove);
			this->makespan = std::move(other.makespan);
			this->pointLS = std::move(other.pointLS);
			this->numLS = std::move(other.numLS);
			this->custInsert = std::move(other.custInsert);
			this->custInsert_delay = std::move(other.custInsert_delay);
			this->tw_resupply = std::move(other.tw_resupply);
			this->seq = std::move(other.seq);
			this->delivery_makespan = std::move(other.delivery_makespan);

			copy_iterator(other);//深拷贝迭代器
		}
		return *this;
	}


	static CustInfo info;//客户信息类的引用
	static Algorithm algorithm;//权重信息
	int construct_index;//记录构造算子的索引

	int localSearch();//局部搜索，适用于第一阶段
	void solveResupply();//补货方案求解

	double get_makespan() const {
		return this->makespan;
	}

	double get_delivery_makespan() const {
		return this->delivery_makespan;
	}

	void updatePoint(int, double);//更新局部算子的得分
	void updateWeight();//更新局部算子的权重

	void cout();//输出结果

	static void initial_LS() {
		Solution::weightLS = algorithm.get_LS();
		Solution::relaxation = 0;
	}

	static void relax() {

	}

private:
	int numRemove;//移除个数
	list<TruckCust> solu;//解
	double makespan;//目标函数
	double delivery_makespan;//配送目标函数
	static vector<double> weightLS;//局部搜索算子的权重
	static double relaxation;//松弛系数
	vector<double> pointLS = vector<double>(6, 0);//局部搜索算子的得分
	vector<int> numLS = vector<int>(6, 0);//局部搜索算子的使用次数

	unordered_set<int> custInsert;//待插入客户
	unordered_set<int> custInsert_delay;//待插入客户(延迟)

	unordered_map<int, pair<double, double>> tw_resupply;//补货方案造成的时间窗约束
	vector<int> seq;//补货方案造成的顺序约束, wu:0
	//vector<int> seq;//补货方案造成的顺序约束, wu:0

	int intinialize(const vector<double>&);//初始化解

	//贪婪构造
	void greedyConstructNoise() {
		this->greedyConstruct(this->custInsert, this->solu.end());//先构造非延迟客户
		this->greedyConstruct(this->custInsert_delay, this->solu.end());//再构造延迟客户
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

	void greedyConstruct(unordered_set<int>&, list<TruckCust>::iterator);
	void randomGreedyConstruct(unordered_set<int>&);
	void calBestInsertion(int, Insertion&, list<TruckCust>::iterator solu_start, list<TruckCust>::iterator solu_end, double max_noise = 0);
	void insert(Insertion&);
	void calMakespan(list<TruckCust>::iterator);
	void constraintCD();

	//局部搜索
	void randomRemove();
	void shawRemove();
	void worstRemove();

	void greedyInsert(unordered_set<int>&);
	void randomGreedyInsert(unordered_set<int>&);

	void seperate(unordered_set<int>&);

	Location find(int);

	void handleRemoveCustomer(int id);

	void remove(Location&);

	//补货方案求解
	void calMakespan_resupply(map<int, DroneResupply>&, ResupplyInfo&);
	void recalculate(list<TruckCust>::iterator, double&);

	
	//构造仓库
	void initial_solu() {
		//初始化解
		this->solu.clear();
		this->solu.push_back(TruckCust(0));
		this->solu.push_back(TruckCust(0));
	}

	void copy_iterator(const Solution&);

	int cal();
	int call();
};
