#include "CustInfo.h"
#include"HyperParameters.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <numeric>


inline double euclideanDistance(const pair<double, double>&, const pair<double, double>&);

//读取文件（客户坐标、需求、延迟时间）
void CustInfo::read(const string& filename) {
	ifstream file_xy(filename);  //打开文件，构造 ifstream 对象

	if (!file_xy) {  // 检查文件是否成功打开
		cerr << "Error opening file_xy!" << endl;
	}

	string line;
	int id;
	double d;
	double x, y, st, sd;
	int sum1 = 0;
	getline(file_xy, line);//忽略第一行

	while (getline(file_xy, line)) {//按空字符分割
		istringstream stream(line);
		stream >> id >> x >> y >> d >> st >> sd;
		//时间矩阵
		this->ST.push_back(st);
		this->SD.push_back(sd);
		this->demand.push_back(d);//需求
		this->xy.emplace_back(x, y);//坐标
		if (stream >> d) {//获得延迟时间
			this->cust_delay[id] = d;
		}
		sum1++;
	}

	file_xy.close();  // 关闭文件
	this->num_cust = sum1 - 1;
	this->is_read = true;
}

//计算距离，时间矩阵、计算补货无人机可达客户
void CustInfo::calculate_dis() {
	if (!this->is_read) {
		cerr << "读取错误" << endl;
		return;
	}

	//初始化矩阵
	this->dis = vector<vector<double>>(static_cast<size_t>(this->num_cust) + 1, vector<double>(static_cast<size_t>(this->num_cust) + 1, 0));
	this->TT = vector<vector<double>>(static_cast<size_t>(this->num_cust) + 1, vector<double>(static_cast<size_t>(this->num_cust) + 1, 0));
	this->TD = vector<vector<double>>(static_cast<size_t>(this->num_cust) + 1, vector<double>(static_cast<size_t>(this->num_cust) + 1, 0));

	//计算矩阵元素
	for (size_t i = 0; i <= this->num_cust; i++) {
		this->dis[i][i] = 0;//距离
		this->TT[i][i] = 0;//卡车时长
		this->TD[i][i] = 0;//无人机时长
		for (size_t j = i + 1; j <= this->num_cust; j++) {
			this->dis[i][j] = euclideanDistance(this->xy[i], this->xy[j]);
			this->TT[i][j] = this->dis[i][j] / SpeedT;
			this->TD[i][j] = this->dis[i][j] / SpeedD;
			
			//对称矩阵
			this->dis[j][i] = this->dis[i][j];
			this->TT[j][i] = this->TT[i][j];
			this->TD[j][i] = this->TD[i][j];
		}
	}

	//计算补货无人机的可达客户
	for (int i = 1; i <= this->num_cust; i++) {
		if (this->dis[0][i] <= Max_dis) this->CD.insert(i);
	}
}

//自动聚类操作
void CustInfo::cluster() {
	double best_score = -1.0;//初始化轮廓系数
	vector<pair<double,double>> best_centers;//存储最优簇中心
	vector<int> best_labels;//存储最优簇标签

	//寻找合适的簇个数
	for (int k = 1; k <= kCluster; ++k) {
		vector<pair<double, double>> centers;
		vector<int> labels;
		kmeans(k, centers, labels);
		double score = silhouetteScore(labels, k);
		if (score > best_score) {
			best_score = score;
			this->best_k = k;
			best_centers = centers;
			best_labels = labels;
		}
	}

	//构建簇集合
	this->clusters = vector<vector<int>>(this->best_k, vector<int>());
	for (int i = 1; i <= this->num_cust; i++) {
		this->clusters[best_labels[i]].push_back(i);
	}
}

//计算欧式距离
inline double euclideanDistance(const pair<double,double>& a, const pair<double, double>& b) {
	return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
}

//k-means，给定k的个数，簇中心集合（空），标签（空），最大迭代数
void CustInfo::kmeans(int k, vector<pair<double,double>>& centers, vector<int>& labels, int max_iter) {
	size_t n = this->num_cust;
	labels = vector<int>(n + 1, -1);
	initializeCenters(centers, k);//初始化簇中心

	//总的迭代框架
	for (int iter = 0; iter < max_iter; ++iter) {
		// 分配每个点到最近中心
		bool changed = false;
		for (size_t i = 1; i <= n; ++i) {
			double minDist = numeric_limits<double>::max();
			int bestCenter = -1;//记录近的簇中心
			for (int j = 0; j < k; ++j) {
				double dist = euclideanDistance(this->xy[i], centers[j]);
				if (dist < minDist) {
					minDist = dist;
					bestCenter = j;
				}
			}

			//更新分类
			if (labels[i] != bestCenter) {
				labels[i] = bestCenter;
				changed = true;
			}
		}

		if (!changed) break;//收敛标准

		// 更新中心
		fill(centers.begin(), centers.end(), pair<double,double>(0.0, 0.0));//新的簇中心
		vector<int> count(k, 0);

		for (size_t i = 1; i <= n; ++i) {
			int cluster = labels[i];
			centers[cluster].first += this->xy[i].first;
			centers[cluster].second += this->xy[i].second;
			count[cluster]++;
		}

		for (int j = 0; j < k; ++j) {
			if (count[j] == 0) continue;
			centers[j].first /= count[j];
			centers[j].second /= count[j];
		}
	}
}

//k-means++初始化
void CustInfo::initializeCenters(vector<pair<double, double>>& centers, int k) {
	//random_device rd;
	//mt19937 gen(rd());

	uniform_int_distribution<int> dis(1, this->num_cust);
	centers.push_back(xy[dis(gen)]);  // 随机选一个点作为第一个中心

	while (centers.size() < (size_t)k) {
		vector<double> distances(static_cast<size_t>(this->num_cust) + 1);

		for (size_t i = 1; i <= this->num_cust; ++i) {
			double minDist = numeric_limits<double>::max();
			for (const auto& center : centers)
				minDist = min(minDist, euclideanDistance(this->xy[i], center));
			distances[i] = pow(minDist,2);
		}

		// 计算概率分布
		discrete_distribution<> prob(distances.begin(), distances.end());//离散分布
		centers.push_back(this->xy[prob(gen)]);
	}
}

// 计算 Silhouette Score
double CustInfo::silhouetteScore(const vector<int>& labels, int k) {
	size_t n = this->num_cust;
	vector<double> silhouette(n + 1);

	for (size_t i = 1; i <= n; ++i) {
		int label_i = labels[i];

		double a = 0.0;
		double b = numeric_limits<double>::max();
		int a_count = 0;

		for (size_t j = 1; j <= n; ++j) {
			if (i == j) continue;
			double dist = euclideanDistance(this->xy[i], this->xy[j]);
			if (labels[j] == label_i) {
				a += dist;
				a_count++;
			}
		}

		if (a_count > 0) a /= a_count;

		// 计算 b（最近其他簇的平均距离）
		for (int l = 0; l < k; ++l) {
			if (l == label_i) continue;

			double b_sum = 0.0;
			int b_count = 0;
			for (size_t j = 1; j <= n; ++j) {
				if (labels[j] == l) {
					b_sum += euclideanDistance(this->xy[i], this->xy[j]);
					b_count++;
				}
			}
			if (b_count > 0)
				b = min(b, b_sum / b_count);
		}

		double s = 0.0;
		if (a < b)
			s = 1.0 - a / b;
		else if (a > b)
			s = b / a - 1.0;
		silhouette[i] = s;
	}

	double sum = accumulate(silhouette.begin(), silhouette.end(), 0.0);
	return sum / n;
}

