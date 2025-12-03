#include "ResupplyInfo.h"
#include "gurobi_c++.h"
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <string>

using namespace std;

class MyLazyCallback : public GRBCallback
{
public:
    MyLazyCallback(unordered_map<pair<int, int>, GRBVar, PairHash> r_,
        unordered_map<int, GRBVar> h_,
        GRBVar eta_,
        vector<SeriesTime>& C_D_,
        vector<int>& C_delay_,
        int n_1_,
        CustInfo& info_,
        unordered_map<pair<int, int>, double, PairHash>& P_)
        :r(move(r_)), h(move(h_)), C_D(C_D_), C_delay(C_delay_), n_1(n_1_), info(info_), eta(move(eta_)), P(P_) {
        this->t = vector<double>(n_1_, 0);
        this->dt = vector<double>(static_cast<size_t>(n_1_) + 1, 0);
        this->dt[this->n_1] = static_cast<double>(INT_MAX);
    };

    double get_t(unsigned int a) const{
        return this->t.at(a);
    }

    vector<double> get_dt() const {
        return this->dt;
    }

    bool feasible() const {
        return this->num != 0;
    }

private:
    unordered_map<pair<int, int>, GRBVar, PairHash> r;  // 访问原模型中的变量
    unordered_map<int, GRBVar> h;
    GRBVar eta;
    vector<SeriesTime>& C_D;
    vector<int>& C_delay;
    CustInfo& info;
    int n_1;
    vector<double> t;
    vector<double> dt;
    unordered_map<pair<int, int>, double, PairHash>& P;
    int num = 0;
    GRBEnv subEnv = GRBEnv(true);

protected:
    void callback() override;
};




void MyLazyCallback::callback() {
    if (where != GRB_CB_MIPSOL) return;

    // 1. 当前整数解
    unordered_map<std::pair<int, int>, double, PairHash> r_bar;
    for (auto& kv : this->r)
        r_bar[kv.first] = getSolution(kv.second);

    unordered_map<int, double> h_bar;
    for (auto& kv : this->h)
        h_bar[kv.first] = getSolution(kv.second);

    try {
        subEnv.set(GRB_IntParam_OutputFlag, 0);
        subEnv.set(GRB_IntParam_InfUnbdInfo, 1);
        subEnv.start();

        GRBModel subModel = GRBModel(subEnv);

        // 子问题变量：t, dt
        unordered_map<int, GRBVar> t;
        unordered_map<int, GRBVar> dt;

        for (auto& i : this->C_D) {
            t[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                "t_" + std::to_string(i.id));
            dt[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                "dt_" + std::to_string(i.id));
        }
        dt[0] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_0");
        dt[n_1] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_last");

        // 目标：min dt[n_1]
        subModel.setObjective(dt.at(n_1) + 0, GRB_MINIMIZE);

        std::vector<GRBConstr> constr;  // 用于取 Pi / FarkasDual
        int n = static_cast<int>(this->C_D.size());

        // ---- 对应原模型约束 39: t[i] >= custDelay(j) * r_ij
        for (auto& i : this->C_D) {
            for (auto& j : this->C_delay) {
                double rhs = this->info.at_custDelay(j) * r_bar.at({ i.id, j });
                constr.push_back(subModel.addConstr(t.at(i.id) >= rhs));
            }
        }

        // ---- 对应原模型 2-6: t[ci1] >= t[ci] + 2*TD(0,ci)*h_ci
        for (int k = 0; k < n - 1; ++k) {
            int ci = this->C_D[k].id;
            int ci1 = this->C_D[k + 1].id;
            double rhs = 2.0 * this->info.at_TD(0, ci) * h_bar.at(ci);
            // t[ci1] - t[ci] >= 2*TD*h
            constr.push_back(subModel.addConstr(t.at(ci1) - t.at(ci) >= rhs));
        }

        // ---- 对应原模型 2-7:
        // t[i] + TD(0,i) - 100000(1-h[i]) - i.t - dt[i] <= 0
        // 等价：dt[i] >= t[i] + TD(0,i) - i.t - 100000(1-h[i])
        for (auto& i : this->C_D) {
            double rhs = this->info.at_TD(0, i.id) - i.t - 100000.0 * (1.0 - h_bar.at(i.id));
            constr.push_back(subModel.addConstr(dt.at(i.id) - t.at(i.id) >= rhs));
        }

        // ---- 对应原模型 2-8:
        // t[i] + TD(0,i) - i.t + i.e + 100000(1-h[i]) - dt[prev] >= 0
        // 等价：dt[prev] <= t[i] + TD(0,i) - i.t + i.e + 100000(1-h[i])
        for (int idx = 0; idx < n; ++idx) {
            auto ci = this->C_D[idx];
            int ci_id = ci.id;
            int prevId = (idx == 0 ? 0 : this->C_D[idx - 1].id);
            double rhs = ci.t - ci.e - this->info.at_TD(0, ci_id)
                - 100000.0 * (1.0 - h_bar.at(ci_id));
            // 原式 rearranged 成 t[ci] - dt[prev] >= rhs
            constr.push_back(subModel.addConstr(t.at(ci_id) - dt.at(prevId) >= rhs));
        }

        // ---- 对应原模型 2-9: dt 单调非减
        // dt[C_D[0]] >= dt[0]
        subModel.addConstr(dt.at(this->C_D[0].id) - dt.at(0) >= 0.0);
        // dt[C_D[i]] >= dt[C_D[i-1]]
        for (int i = 1; i < n; ++i) {
            int ci = this->C_D[i].id;
            int cim1 = this->C_D[i - 1].id;
            subModel.addConstr(dt.at(ci) - dt.at(cim1) >= 0.0);
        }
        // dt[n_1] >= dt[last]
        subModel.addConstr(dt.at(this->n_1) - dt.at(this->C_D[n - 1].id) >= 0.0);

        // ---- 对应原模型 2-10: dt[0] == 0
        subModel.addConstr(dt.at(0) == 0.0);

        // 6. 求解子问题
        subModel.optimize();
        int subStatus = subModel.get(GRB_IntAttr_Status);

        if (subStatus == GRB_OPTIMAL) {
            double z_sub = dt[this->n_1].get(GRB_DoubleAttr_X);
            double eta_val = getSolution(this->eta);

            // 顺便保存最优 t, dt（你原先的逻辑）
            if (z_sub < this->dt.at(this->n_1)) {
                for (auto& i : this->C_D) {
                    this->t[i.id] = t.at(i.id).get(GRB_DoubleAttr_X);
                    this->dt[i.id] = dt.at(i.id).get(GRB_DoubleAttr_X);
                }
                this->dt[0] = 0.0;
                this->dt[this->n_1] = z_sub;
            }

            const double TOL = 1e-6;
            if (eta_val + TOL >= z_sub) {
                return; // 当前 eta 已经不低估子问题值，不加最优割
            }

            // ===== 构造最优性割：eta >= cut =====
            GRBLinExpr cut = 0;
            int index = 0;

            // 约束 39: t[i] >= delay(j)*r_ij → 对 r 的项
            for (auto& i : this->C_D) {
                for (auto& j : this->C_delay) {
                    double pi = constr[index].get(GRB_DoubleAttr_Pi);
                    cut += pi * this->info.at_custDelay(j) * this->r.at({ i.id, j });
                    ++index;
                }
            }

            // 约束 2-6: t[ci1]-t[ci] >= 2*TD(0,ci)*h_ci
            for (int k = 0; k < n - 1; ++k) {
                int ci = this->C_D[k].id;
                double pi = constr[index].get(GRB_DoubleAttr_Pi);
                cut += pi * 2.0 * this->info.at_TD(0, ci) * this->h.at(ci);
                ++index;
            }

            // 约束 2-7: dt[i] - t[i] >= TD(0,i) - i.t - M(1-h[i])
            for (auto& i : this->C_D) {
                double pi = constr[index].get(GRB_DoubleAttr_Pi);
                cut += pi * (this->info.at_TD(0, i.id) - i.t
                    - 100000.0 * (1.0 - this->h.at(i.id)));
                ++index;
            }

            // 约束 2-8: t[i] - dt[prev] >= t_i - e_i - TD(0,i) - M(1-h[i])
            for (int idx = 0; idx < n; ++idx) {
                auto ci = this->C_D[idx];
                double pi = constr[index].get(GRB_DoubleAttr_Pi);
                cut += pi * (ci.t - ci.e - this->info.at_TD(0, ci.id)
                    - 100000.0 * (1.0 - this->h.at(ci.id)));
                ++index;
            }

            // 最优性割：eta >= cut
            addLazy(this->eta >= cut);
        }
        else if (subStatus == GRB_INFEASIBLE) {
            // ===== 子问题不可行 → 可行性割 =====
            GRBLinExpr cut = 0;
            int index = 0;

            // 结构和最优性割一样，只是用 FarkasDual 并写成 <= 0
            for (auto& i : this->C_D) {
                for (auto& j : this->C_delay) {
                    double piF = constr[index].get(GRB_DoubleAttr_FarkasDual);
                    cut += piF * this->info.at_custDelay(j) * this->r.at({ i.id, j });
                    ++index;
                }
            }

            for (int k = 0; k < n - 1; ++k) {
                int ci = this->C_D[k].id;
                double piF = constr[index].get(GRB_DoubleAttr_FarkasDual);
                cut += piF * 2.0 * this->info.at_TD(0, ci) * this->h.at(ci);
                ++index;
            }

            for (auto& i : this->C_D) {
                double piF = constr[index].get(GRB_DoubleAttr_FarkasDual);
                cut += piF * (this->info.at_TD(0, i.id) - i.t
                    - 100000.0 * (1.0 - this->h.at(i.id)));
                ++index;
            }

            for (int idx = 0; idx < n; ++idx) {
                auto ci = this->C_D[idx];
                double piF = constr[index].get(GRB_DoubleAttr_FarkasDual);
                cut += piF * (ci.t - ci.e - this->info.at_TD(0, ci.id)
                    - 100000.0 * (1.0 - this->h.at(ci.id)));
                ++index;
            }

            addLazy(cut <= 0.0);
        }
    }
    catch (GRBException& e) {
        std::cout << "Subproblem error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...) {
        std::cout << "Unknown error during subproblem optimization." << std::endl;
    }
}



void ResupplyInfo::branch_benders_cut(map<int, DroneResupply>& obj) {
    try {
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_OutputFlag, 0);
        env.start();

        GRBModel model = GRBModel(env);

        // 主问题变量：r, h, eta
        unordered_map<pair<int, int>, GRBVar, PairHash> r;
        unordered_map<int, GRBVar> h;
        // 如果确定 dt[n_1] >= 0，可以用 0；否则更安全的是 -INF
        GRBVar eta = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "eta");

        for (auto& i : this->C_D) {
            h[i.id] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                "h_" + to_string(i.id));
            for (auto& j : this->C_delay) {
                r[{i.id, j}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                    "r_" + to_string(i.id) + "_" + to_string(j));
            }
        }

        // 目标：min eta
        model.setObjective(eta + 0, GRB_MINIMIZE);

        // 约束 37: 容量约束
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay) {
                expr += this->info.at_demand(j) * r.at({ i.id, j });
            }
            model.addConstr(expr <= QDr, "constr37_" + to_string(i.id));
        }

        // 约束 38: 每个 j 被恰好一个 i 选择
        for (auto& j : C_delay) {
            GRBLinExpr expr = 0;
            for (auto& i : this->C_D) {
                expr += r.at({ i.id, j });
            }
            model.addConstr(expr == 1, "constr38_" + to_string(j));
        }

        // 约束 2-4: r 和 h 的逻辑关系（上界）
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay)
                expr += r.at({ i.id, j });
            // 你原始模型里这里是 1000*h[i]（我直接沿用）
            model.addConstr(expr <= 1000.0 * h.at(i.id),
                "constr2-4_" + to_string(i.id));
        }

        // 约束 2-5: r 和 h 的逻辑关系（下界）
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay)
                expr += r.at({ i.id, j });
            model.addConstr(expr >= h.at(i.id),
                "constr2-5_" + to_string(i.id));
        }

        // 约束 2-11: 不允许的 r[i,j] 强制为 0
        for (auto& i : this->C_D) {
            for (auto& j : C_delay) {
                if (i.AP.count(j) == 0) {
                    model.addConstr(r.at({ i.id, j }) == 0,
                        "constr2-11_" + to_string(i.id) + "_" + to_string(j));
                }
            }
        }

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

        // 绑定 callback
        MyLazyCallback cb(r, h, eta, this->C_D, this->C_delay,
            this->n_1, this->info, this->P);
        model.setCallback(&cb);

        model.optimize();

        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            // 回收主问题解
            for (auto& i : this->C_D) {
                for (auto& j : this->C_delay) {
                    if (r.at({ i.id, j }).get(GRB_DoubleAttr_X) >= 0.99)
                        obj[i.id].goods.insert(j);
                }
                if (obj.find(i.id) != obj.end())
                    obj[i.id].time = cb.get_t(i.id);
            }
            this->dt = cb.get_dt();
            this->dt.back() = model.get(GRB_DoubleAttr_ObjVal);
            std::cout << "      Branch-and-Benders-cut 求解 = "
                << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
        }
        else {
            std::cout << "No optimal solution found in BBC." << std::endl;
            model.computeIIS();
            model.write("infeasible_bbc.ilp");
        }
    }
    catch (GRBException& e) {
        std::cout << "Gurobi exception: " << e.getMessage() << std::endl;
    }
    catch (...) {
        std::cout << "Unknown error during optimization." << std::endl;
    }
}


void ResupplyInfo::gurobi_solve(map<int, DroneResupply>& obj) {
    try {
        // 创建环境
        GRBEnv env = GRBEnv(true);
        //env.set("LogFile", "gurobi.log");
        env.set(GRB_IntParam_OutputFlag, 0); // 不输出日志
        env.start();

        // 创建模型
        GRBModel model = GRBModel(env);

        // 添加变量
        map<pair<int, int>, GRBVar> r;
        map<int, GRBVar> h;
        map<int, GRBVar> t;
        map<int, GRBVar> dt;    

        for (auto& i : this->C_D) {
            h[i.id] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "h_" + to_string(i.id));
            for (auto& j : this->C_delay) {
                r[{i.id, j}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "r_" + to_string(i.id) + "_" + to_string(j));
            }
        }

        for (auto& i : this->C_D) {
            t[i.id] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t_" + to_string(i.id));
            dt[i.id] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(i.id));
        }
        dt[0] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(0));
        dt[n_1] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(n_1));

        // 设置目标函数: maximize x + y
        model.setObjective(dt[n_1] + 0, GRB_MINIMIZE);

        // 添加约束: 37
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay) {
                expr += this->info.at_demand(j) * r.at({i.id, j});
            }
            model.addConstr(expr <= QDr, "constr37_" + to_string(i.id));
        }

        // 添加约束: 38
        for (auto& j : C_delay) {
            GRBLinExpr expr = 0;
            for (auto& i : this->C_D) {
                expr += r.at({i.id, j});
            }
            model.addConstr(expr == 1, "constr38_" + to_string(j));
        }

        //添加约束：39
        for (auto& i : this->C_D) {
            for (auto& j : C_delay) {
                model.addConstr(t.at(i.id) >= this->info.at_custDelay(j) * r.at({i.id, j}), "constr39_" + to_string(i.id) + to_string(j));
            }
        }

        // 添加约束: 2-4
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay) {
                expr += r.at({i.id, j});
            }
            model.addConstr(expr <=  1000 * h.at(i.id), "constr2-4_" + to_string(i.id));
        }

        // 添加约束: 2-5
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay) {
                expr += r.at({i.id, j});
            }
            model.addConstr(expr >= h.at(i.id), "constr2-5_" + to_string(i.id));
        }

        //添加约束2-6
        int n = this->C_D.size();
        for (int i = 0; i < n - 1; ++i) {
            int ci = this->C_D[i].id;
            int ci1 = this->C_D[i + 1].id;
            GRBLinExpr expr = t.at(ci) + 2 * this->info.at_TD(0, ci) * h.at(ci);
            model.addConstr(expr <= t.at(ci1), "constr2-6_" + to_string(ci));
        }

        //添加约束2-7
        for (auto& i : this->C_D) {
            GRBLinExpr expr = t.at(i.id) + this->info.at_TD(0, i.id) - 100000 * (1 - h.at(i.id)) - i.t - dt.at(i.id);
            model.addConstr(expr <= 0, "constr2-7_" + to_string(i.id));
        }

        //添加约束2-8
        for (int i = 0; i < n; ++i) {
            auto ci = this->C_D[i];
            GRBLinExpr expr = t.at(ci.id) + this->info.at_TD(0, ci.id) - ci.t + ci.e + 100000 * (1 - h.at(ci.id));
            expr += i == 0 ? -dt.at(0) : -dt.at(this->C_D[i - 1].id);
            model.addConstr(expr >= 0, "constr2-8_" + to_string(ci.id));
        }

        //添加约束2-9
        model.addConstr(dt.at(this->C_D[0].id) >= dt.at(0), "constr2-9_" + to_string(this->C_D[0].id));
        for (int i = 1; i < n; ++i) {
            model.addConstr(dt.at(this->C_D[i].id) >= dt.at(this->C_D[i - 1].id), "constr2-9_" + to_string(this->C_D[i].id));;
        }
        model.addConstr(dt.at(n_1) >= dt.at(this->C_D[n - 1].id), "constr2-9_" + to_string(this->C_D[n - 1].id));

        //添加约束2-10
        model.addConstr(dt.at(0) == 0, "constr2-10");

        // 添加约束: 2-11
        for (auto& i : this->C_D) {
            for (auto& j : C_delay) {
                if (i.AP.count(j) == 0)
                    model.addConstr(r.at({i.id, j}) == 0, "constr2-11_" + to_string(i.id));
            }
        }

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

        // 求解模型
        model.optimize();

        // 输出结果
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            cout << "Gurobi求解 = " << model.get(GRB_DoubleAttr_ObjVal) << endl;
        }
        else {
            cout << "No optimal solution found." << endl;
            model.computeIIS();
            model.write("infeasible.ilp");
        }
        int yy = 0;
    }
    catch (GRBException& e) {
        cout << "Gurobi exception: " << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Unknown error during optimization." << endl;
    }

}

//std::vector<GRBVar> vars;  // 已定义并加入模型
//std::vector<double> start_vals;  // 热启动用的值

//for (int i = 0; i < vars.size(); ++i) {
//    vars[i].set(GRB_DoubleAttr_Start, start_vals[i]);
//}