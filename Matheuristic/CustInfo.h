#pragma once

#include<string>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<iostream>

using namespace std;

//客户信息类
class CustInfo
{//记录客户信息
public:
	//初始化客户信息
	CustInfo(const string& address) {
		this->read(address);//读取文件
		this->calculate_dis();//计算所需参数
		this->cluster();//聚类
	}
	
	//获取距离
	double at_dis(int a, int b) const{
		return this->dis.at(a).at(b);
	}
	//获取无人机时长
	double at_TD(unsigned int a, unsigned int b) const {
		if (a >= TD.size() || b >= TD[0].size())
			std::cout << a << ' ' << b << endl;
		return this->TD.at(a).at(b);
	}
	//获取卡车时长
	double at_TT(unsigned int a, unsigned int b) const {
		return this->TT.at(a).at(b);
	}
	//获取客户点的坐标
	pair<double, double> at_xy(unsigned int a) const {
		return this->xy.at(a);
	}
	//获取客户点的延迟时间（默认初始时间为0）
	double at_custDelay(unsigned int a) const {
		if (this->cust_delay.find(a) == this->cust_delay.end()) {
			return 0;
		}
		return this->cust_delay.at(a);
	}
	//获取总客户数量
	unsigned int get_num() const {
		return this->num_cust;
	}
	//获取客户需求量
	double at_demand(unsigned int a) const {
		return this->demand.at(a);
	}
	//获取聚类结果
	const vector<vector<int>>& at_clusters() const {
		return this->clusters;
	}
	//获取簇的个数
	unsigned int get_k() const {
		return this->best_k;
	}
	//获取卡车服务时长
	double at_ST(unsigned int a) const {
		return this->ST.at(a);
	}
	//获取无人机服务时长
	double at_SD(unsigned int a) const {
		return this->SD.at(a);
	}
	//判断客户是否可作为补货点
	bool in_CD(unsigned int a) const {
		if (this->CD.count(a) > 0) return true;
		else return false;
	}

private:
	unsigned int best_k = 1;//初始化最优簇的个数
	vector<vector<double>> TT;//卡车时间矩阵
	vector<vector<double>> TD;//无人机时间矩阵
	vector<vector<double>> dis;//距离矩阵
	vector<double> ST;//卡车服务时间
	vector<double> SD;//无人机服务时间
	vector<vector<int>> clusters;//簇集合
	vector<pair<double, double>> xy;//坐标
	vector<double> demand;//需求
	unordered_map<int, double> cust_delay;//延迟客户
	unordered_set<int> CD;//补货无人机可达客户
	int num_cust = 0;//客户数量
	bool is_read = false;//读取是否正确

	void read(const string&);//读文件
	void calculate_dis();//计算距离
	void cluster();//聚类操作
	void kmeans(int k, vector<pair<double, double>>& centers, vector<int>& labels, int max_iter = 100);
	void initializeCenters(vector<pair<double, double>>&, int k);
	double silhouetteScore(const vector<int>& labels, int k);
};

