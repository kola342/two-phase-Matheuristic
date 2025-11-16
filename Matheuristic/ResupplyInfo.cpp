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
    MyLazyCallback(map<pair<int, int>, GRBVar> r_,
        map<int, GRBVar> h_,
        GRBVar eta_,
        vector<SeriesTime>& C_D_,
        vector<int>& C_delay_,
        int n_1_,
        CustInfo& info_,
        unordered_map<pair<int, int>, double, PairHash>& P_)
        :r(move(r_)), h(move(h_)), C_D(C_D_), C_delay(C_delay_), n_1(n_1_), info(info_), eta(move(eta_)), P(P_) {
        this->t = vector<double>(n_1_, 0);
        this->dt = vector<double>(static_cast<size_t>(n_1_) + 1, 0);
    };

    double get_t(unsigned int a) const{
        return this->t.at(a);
    }

    vector<double> get_dt() const {
        return this->dt;
    }

private:
    map<pair<int, int>, GRBVar> r;  // 访问原模型中的变量
    map<int, GRBVar> h;
    GRBVar eta;
    vector<SeriesTime>& C_D;
    vector<int>& C_delay;
    CustInfo& info;
    int n_1;
    vector<double> t;
    vector<double> dt;
    unordered_map<pair<int, int>, double, PairHash>& P;
    double sp();

protected:
    void callback() override;
};

double MyLazyCallback::sp() {
    // 获取当前解的变量值
    map<pair<int, int>, int> r_bar;
    for (auto& kv : this->r) {
        r_bar[kv.first] = round(getSolution(kv.second));
    }
    map<int, int> h_bar;
    for (auto& kv : this->h) {
        h_bar[kv.first] = round(getSolution(kv.second));
    }

    try {
        // 创建一个独立的环境和模型用于子问题求解
        GRBEnv subEnv1 = GRBEnv(true);
        subEnv1.set(GRB_IntParam_OutputFlag, 0); // 不输出日志
        subEnv1.start();

        GRBModel subModel1 = GRBModel(subEnv1);

        // 创建子问题变量（如 y[i] ∈ [0,1]）
        map<int, GRBVar> t;
        map<int, GRBVar> dt;

        for (auto& i : this->C_D) {
            t[i.id] = subModel1.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t_" + to_string(i.id));
            dt[i.id] = subModel1.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(i.id));
        }
        dt[0] = subModel1.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(0));
        dt[this->n_1] = subModel1.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(this->n_1));


        // 设置目标函数
        subModel1.setObjective(dt[this->n_1] + 0, GRB_MINIMIZE);


        //vector<GRBConstr> constr_t;//t的对偶约束
        //约束1
        for (auto& i : this->C_D) {
            for (auto& j : C_delay) {
                if (r_bar[{i.id, j}] == 1)
                    subModel1.addConstr(t.at(i.id) >= this->info.at_custDelay(j));
            }
        }

        //约束2
        size_t n = this->C_D.size();
        for (size_t k = 0; k < n - 1; k++) {
            int ci = this->C_D[k].id;
            int ci1 = this->C_D[k + 1].id;
            subModel1.addConstr(t.at(ci1) - t.at(ci) >= 2 * this->info.at_TD(0, ci) * h_bar.at(ci));
        }

        //约束3
        for (auto& i : this->C_D) {
            if (h_bar[i.id] == 1)
                subModel1.addConstr(dt.at(i.id) - t.at(i.id) >= this->info.at_TD(0, i.id) - i.t);
        }

        //vector<GRBConstr> constr_dt;//dt的对偶约束

        //约束4
        for (int i = 0; i < n; ++i) {
            int ci = this->C_D[i].id;
            int c1i = i == 0 ? 0 : this->C_D[i - 1].id;
            if (h_bar[ci] == 1)
                subModel1.addConstr(t.at(ci) - dt.at(c1i) >= this->C_D[i].t - this->C_D[i].e - this->info.at_TD(0, ci));
            subModel1.addConstr(dt.at(ci) - dt.at(c1i) >= 0);
        }


        //约束5
        subModel1.addConstr(dt.at(this->n_1) - dt.at(this->C_D.back().id) >= 0);
        subModel1.addConstr(dt.at(0) == 0);

        //约束6
        for (auto& [i, j] : this->P) {
            subModel1.addConstr(dt.at(i.first) - dt.at(i.second) >= j);
        }

        // 求解
        subModel1.optimize();
        return subModel1.get(GRB_IntAttr_Status) == GRB_OPTIMAL ? subModel1.get(GRB_DoubleAttr_ObjVal) : -1000;
    }
    catch (GRBException& e) {
        cout << "Gurobi exception: " << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Unknown error during optimization." << endl;
    }
}


/*
void MyLazyCallback::callback() {   //对偶子问题
    if (where == GRB_CB_MIPSOL) {  // 当前是可行整数解节点
        // 获取当前解的变量值
        map<pair<int, int>, int> r_bar;
        for (auto& kv : this->r) {
            r_bar[kv.first] = round(getSolution(kv.second));
        }
        map<int, int> h_bar;
        for (auto& kv : this->h) {
            h_bar[kv.first] = round(getSolution(kv.second));
        }

        try {
            // 创建一个独立的环境和模型用于子问题求解
            GRBEnv subEnv = GRBEnv(true);
            subEnv.set(GRB_IntParam_OutputFlag, 0); // 不输出日志
            subEnv.start();

            GRBModel subModel = GRBModel(subEnv);

            // 创建子问题变量（如 y[i] ∈ [0,1]）
            map<pair<int, int>, GRBVar> alpha;//==1
            map<int, GRBVar> beta;
            map<int, GRBVar> gamma;//==1
            map<int, GRBVar> delta;//==1
            map<int, GRBVar> epsilon;
            map<int, GRBVar> lambda;//==1
            for (auto& i : this->C_D) {
                if (h_bar[i.id] == 1) {
                    gamma[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "gamma_" + to_string(i.id));
                    delta[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "delta_" + to_string(i.id));
                }
                epsilon[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "epsilon_" + to_string(i.id));
                if (i.id != this->C_D.back().id) beta[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "beta_" + to_string(i.id));
                for (auto& j : this->C_delay) {
                    if (r_bar[{i.id, j}] == 1)
                        alpha[{i.id, j}] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "alpha_" + to_string(i.id) + "_" + to_string(j));
                }
            }
            epsilon[this->n_1] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "epsilon_" + to_string(this->n_1));
            GRBVar zeta = subModel.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "zeta_" + to_string(1));

            for (auto& [i, j] : this->P) {
                lambda[i.first] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "lambda" + to_string(i.first));
            }


            // 设置目标函数
            GRBLinExpr expr0 = 0;
            for (auto& i : this->C_D) {
                double TD = this->info.at_TD(0, i.id);
                for (auto& j : C_delay) {
                    if (r_bar.at({ i.id, j }) == 1)
                        expr0 += this->info.at_custDelay(j) * alpha.at({ i.id, j });
                }
                if (h_bar.at(i.id) == 1) {
                    expr0 += (TD - i.t) * gamma.at(i.id);
                    expr0 += (i.t - TD - i.e) * delta.at(i.id);
                }
                if (i.id != this->C_D.back().id) expr0 += 2 * TD * h_bar.at(i.id) * beta.at(i.id);
            }

            for (auto& [i, j] : this->P) {
                expr0 += j * lambda.at(i.first);
            }

            subModel.setObjective(expr0, GRB_MAXIMIZE);


            //vector<GRBConstr> constr_t;//t的对偶约束
            //约束1
            int c1 = this->C_D.front().id;
            GRBLinExpr expr1 = h_bar.at(c1) == 1 ? -beta.at(c1) + (delta.at(c1) - gamma.at(c1)) : -beta.at(c1);
            for (auto& j : C_delay) {
                if (r_bar.at({ c1, j }) == 1)
                    expr1 += alpha.at({ c1, j });
            }
            subModel.addConstr(expr1 <= 0, "constr1_");

            //约束2
            size_t n = this->C_D.size();
            for (size_t k = 1; k < n - 1; k++) {
                int ci = this->C_D[k].id;
                GRBLinExpr expr2 = 0;
                for (auto& j : C_delay) {
                    if (r_bar.at({ ci, j }) == 1)
                        expr2 += alpha.at({ ci, j });
                }
                expr2 += h_bar.at(ci) == 1 ? -beta.at(ci) + beta.at(this->C_D[k - 1].id) + (delta.at(ci) - gamma.at(ci)) : -beta.at(ci) + beta.at(this->C_D[k - 1].id);
                subModel.addConstr(expr2 <= 0, "constr2_" + to_string(ci));
            }

            //约束3
            int cn = this->C_D[n - 1].id;
            GRBLinExpr expr3 = 0;
            for (auto& j : C_delay) {
                if (r_bar.at({ cn, j }) == 1)
                    expr3 += alpha.at({ cn, j });
            }
            expr3 += h_bar.at(cn) == 1 ? beta.at(this->C_D[n - 2].id) + (delta.at(cn) - gamma.at(cn)) : beta.at(this->C_D[n - 2].id);
            subModel.addConstr(expr3 <= 0, "constr2_" + to_string(cn));

            //vector<GRBConstr> constr_dt;//dt的对偶约束

            //约束4
            GRBLinExpr expr4 = (h_bar.at(this->C_D[1].id) == 1 && h_bar.at(this->C_D[0].id) == 1) ? zeta - delta.at(this->C_D[0].id) - epsilon.at(this->C_D[0].id)
                : zeta - epsilon.at(this->C_D[0].id);
            for (auto& j : this->P) {
                if (j.first.first == 0) {
                    expr4 += lambda.at(0);
                    break;
                }
            }
            subModel.addConstr(expr4 <= 0, "constr3");


            //约束5
            for (size_t k = 0; k < n - 1; k++) {
                int ci = this->C_D[k].id;
                int ci1 = this->C_D[k + 1].id;
                GRBLinExpr expr5 = h_bar.at(ci) == 1 ? gamma.at(ci) + epsilon.at(ci) : epsilon.at(ci);
                expr5 += h_bar.at(ci1) == 1 ? -delta.at(ci1) - epsilon.at(ci1) : -epsilon.at(ci1);

                for (auto& j : this->P) {
                    if (j.first.first == ci) {
                        expr5 += lambda.at(ci);
                        break;
                    }
                }
                for (auto& j : this->P) {
                    if (j.first.second == ci) {
                        expr5 -= lambda.at(j.first.first);
                        break;
                    }
                }

                subModel.addConstr(expr5 <= 0, "constr4_" + to_string(ci));
            }

            //约束6
            GRBLinExpr expr6 = h_bar.at(cn) == 1 ? gamma.at(cn) - epsilon.at(this->n_1) + epsilon.at(cn) : -epsilon.at(this->n_1) + epsilon.at(cn);
            for (auto& j : this->P) {
                if (j.first.first == cn) {
                    expr6 += lambda.at(cn);
                    break;
                }
            }
            for (auto& j : this->P) {
                if (j.first.second == cn) {
                    expr6 -= lambda.at(j.first.first);
                    break;
                }
            }
            subModel.addConstr(expr6 <= 0, "constr5");

            //约束7
            GRBLinExpr expr7 = epsilon.at(this->n_1);
            for (auto& j : this->P) {
                if (j.first.second == this->n_1) {
                    expr7 -= lambda.at(j.first.first);
                    break;
                }
            }
            subModel.addConstr(expr7 <= 1, "constr6");


            // 求解
            subModel.optimize();

            if (subModel.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
                //计算子问题目标函数表达式
                cout << "对偶子问题的求解结果" << subModel.get(GRB_DoubleAttr_ObjVal) << endl;
                cout << "子问题的求解结果" << sp() << endl;

                GRBLinExpr cut = 0;
                for (auto& i : this->C_D) {
                    double TD = this->info.at_TD(0, i.id);
                    for (auto& j : C_delay) {
                        if (r_bar.at({ i.id, j }) == 1) {
                            cut += this->info.at_custDelay(j) * this->r.at({ i.id, j }) * alpha.at({ i.id, j }).get(GRB_DoubleAttr_X);
                        }
                        if (h_bar.at(i.id) == 1) {
                            cut += (TD - i.t) * gamma.at(i.id).get(GRB_DoubleAttr_X);
                            cut += (i.t - TD - i.e) * delta.at(i.id).get(GRB_DoubleAttr_X);
                        }
                        if (i.id != this->C_D.back().id) cut += 2 * TD * this->h.at(i.id) * beta.at(i.id).get(GRB_DoubleAttr_X);
                    }

                    for (auto& [i, j] : this->P) {
                        cut += j * lambda.at(i.first).get(GRB_DoubleAttr_X);
                    }

                    addLazy(cut <= this->eta);
                }
            }
            else if (subModel.get(GRB_IntAttr_Status) == GRB_UNBOUNDED) {
                //cout << "可行割" << endl;

                GRBLinExpr cut = 0;
                for (auto& i : this->C_D) {
                    double TD = this->info.at_TD(0, i.id);
                    for (auto& j : C_delay) {
                        if (r_bar.at({ i.id, j }) == 1)
                            cut += this->info.at_custDelay(j) * this->r.at({ i.id, j }) * alpha.at({ i.id, j }).get(GRB_DoubleAttr_UnbdRay);
                    }
                    if (h_bar.at(i.id) == 1) {
                        cut += (TD - i.t) * gamma.at(i.id).get(GRB_DoubleAttr_UnbdRay);
                        cut += (i.t - TD - i.e) * delta.at(i.id).get(GRB_DoubleAttr_UnbdRay);
                    }
                    if (i.id != this->C_D.back().id) cut += 2 * TD * this->h.at(i.id) * beta.at(i.id).get(GRB_DoubleAttr_UnbdRay);
                }

                for (auto& [i, j] : this->P) {
                    cut += j * lambda.at(i.first).get(GRB_DoubleAttr_UnbdRay);
                }

                addLazy(cut <= 0);
            }
        }
        catch (GRBException& e) {
            cout << "Gurobi exception: " << e.getMessage() << endl;
        }
        catch (...) {
            cout << "Unknown error during optimization." << endl;
        }
    }
}
*/



void MyLazyCallback::callback() {//子问题
    if (where == GRB_CB_MIPSOL) {  // 当前是可行整数解节点
        // 获取当前解的变量值
        map<pair<int, int>, double> r_bar;
        for (auto& kv : this->r) {
            r_bar[kv.first] = getSolution(kv.second);
        }
        map<int, double> h_bar;
        for (auto& kv : this->h) {
            h_bar[kv.first] = getSolution(kv.second);
        }

        try {
            // 创建一个独立的环境和模型用于子问题求解
            GRBEnv subEnv = GRBEnv(true);
            subEnv.set(GRB_IntParam_OutputFlag, 0); // 不输出日志
            subEnv.start();

            GRBModel subModel = GRBModel(subEnv);

            // 创建子问题变量（如 y[i] ∈ [0,1]）
            map<int, GRBVar> t;
            map<int, GRBVar> dt;

            for (auto& i : this->C_D) {
                t[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t_" + to_string(i.id));
                dt[i.id] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(i.id));
            }
            dt[0] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(0));
            dt[n_1] = subModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "dt_" + to_string(n_1));


            // 设置目标函数
            subModel.setObjective(dt[n_1] + 0, GRB_MINIMIZE);


            vector<GRBConstr> constr;//t的对偶约束
            //约束1
            for (auto& i : this->C_D) {
                for (auto& j : C_delay) {
                    if (r_bar.at({ i.id, j }) == 1)
                        constr.push_back(
                            subModel.addConstr(t[i.id] >= this->info.at_custDelay(j)));
                }
            }

            //约束2
            size_t n = this->C_D.size();
            for (size_t k = 0; k < n - 1; k++) {
                int ci = this->C_D[k].id;
                int ci1 = this->C_D[k + 1].id;
                constr.push_back(
                    subModel.addConstr(t[ci1] - t[ci] >= 2 * this->info.at_TD(0, ci) * h_bar.at(ci)));
            }

            //约束3
            for (auto& i : this->C_D) {
                if (h_bar.at(i.id) == 1)
                    constr.push_back(
                        subModel.addConstr(dt[i.id] - t[i.id] >= this->info.at_TD(0, i.id) - i.t));
            }

            //约束4
            for (size_t i = 0; i < n; ++i) {
                int ci = this->C_D[i].id;
                int c1i = i == 0 ? 0 : this->C_D[i - 1].id;
                if (h_bar.at(ci) == 1)
                    constr.push_back(
                        subModel.addConstr(t[ci] - dt[c1i] >= this->C_D[i].t - this->C_D[i].e - this->info.at_TD(0, ci)));
            }

            //约束5
            for (size_t i = 0; i < n; ++i) {
                int ci = this->C_D[i].id;
                int c1i = i == 0 ? 0 : this->C_D[i - 1].id;
                subModel.addConstr(dt[ci] - dt[c1i] >= 0);
            }
            subModel.addConstr(dt[this->n_1] - dt[this->C_D.back().id] >= 0);

            //约束6
            subModel.addConstr(dt[0] == 0);

            //约束7
            for (auto& [i, j] : this->P) {
                constr.push_back(
                    subModel.addConstr(dt[i.first] - dt[i.second] >= j));
            }

            // 求解
            subModel.optimize();

            if (subModel.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {

                //获得t、dt
                for (auto& i : this->C_D) {
                    this->t[i.id] = t[i.id].get(GRB_DoubleAttr_X);
                    this->dt[i.id] = dt[i.id].get(GRB_DoubleAttr_X);
                }
                this->dt[0] = 0;
                this->dt[this->n_1] = dt[this->n_1].get(GRB_DoubleAttr_X);

                //计算对偶子问题目标函数表达式
                GRBLinExpr cut = 0;
                int index = 0;
                for (auto& i : this->C_D) {
                    for (auto& j : C_delay) {
                        if (r_bar.at({i.id, j}) == 1) {
                            cut += constr[index++].get(GRB_DoubleAttr_Pi) * this->info.at_custDelay(j) * this->r.at({ i.id, j });
                        }
                    }
                }

                //约束2
                size_t n = this->C_D.size();
                for (size_t k = 0; k < n - 1; k++) {
                    int ci = this->C_D[k].id;
                    cut += constr[index++].get(GRB_DoubleAttr_Pi) * 2 * this->info.at_TD(0, ci) * this->h.at(ci);
                }

                //约束3
                for (auto& i : this->C_D) {
                    if (h_bar.at(i.id) == 1)
                        cut += constr[index].get(GRB_DoubleAttr_Pi) * constr[index++].get(GRB_DoubleAttr_RHS);
                }

                //约束4
                for (size_t i = 0; i < n; ++i) {
                    int ci = this->C_D[i].id;
                    int c1i = i == 0 ? 0 : this->C_D[i - 1].id;
                    if (h_bar.at(ci) == 1)
                        cut += constr[index].get(GRB_DoubleAttr_Pi) * constr[index++].get(GRB_DoubleAttr_RHS);
                }

                //约束7
                for (auto& [i, j] : this->P) {
                    cut += constr[index].get(GRB_DoubleAttr_Pi) * constr[index++].get(GRB_DoubleAttr_RHS);
                }

                addLazy(cut <= this->eta);
            }
            else if (subModel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {

                //计算对偶子问题目标函数表达式
                GRBLinExpr cut = 0;
                int index = 0;
                for (auto& i : this->C_D) {
                    for (auto& j : C_delay) {
                        if (r_bar.at({ i.id, j }) == 1) {
                            cut += constr[index++].get(GRB_DoubleAttr_FarkasDual) * this->info.at_custDelay(j) * this->r.at({ i.id, j });
                        }
                    }
                }

                //约束2
                size_t n = this->C_D.size();
                for (size_t k = 0; k < n - 1; k++) {
                    int ci = this->C_D[k].id;
                    cut += constr[index++].get(GRB_DoubleAttr_FarkasDual) * 2 * this->info.at_TD(0, ci) * this->h.at(ci);
                }

                //约束3
                for (auto& i : this->C_D) {
                    if (h_bar.at(i.id) == 1)
                        cut += constr[index].get(GRB_DoubleAttr_FarkasDual) * constr[index++].get(GRB_DoubleAttr_RHS);
                }

                //约束4
                for (size_t i = 0; i < n; ++i) {
                    int ci = this->C_D[i].id;
                    int c1i = i == 0 ? 0 : this->C_D[i - 1].id;
                    if (h_bar.at(ci) == 1)
                        cut += constr[index].get(GRB_DoubleAttr_FarkasDual) * constr[index++].get(GRB_DoubleAttr_RHS);
                }

                //约束7
                for (auto& [i, j] : this->P) {
                    cut += constr[index].get(GRB_DoubleAttr_FarkasDual) * constr[index++].get(GRB_DoubleAttr_RHS);
                }

                addLazy(cut <= 0);
            }
        }
        catch (GRBException& e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        }
        catch (...) {
            std::cout << "Unknown error during optimization." << std::endl;
        }
    }
}


void ResupplyInfo::branch_benders_cut(map<int, DroneResupply>& obj) {
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
        GRBVar eta = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "eta");

        for (auto& i:this->C_D) {
            h[i.id] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "h_" + to_string(i.id));
            for (auto& j : this->C_delay) {
                r[{i.id, j}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "r_" + to_string(i.id) + "_" + to_string(j));
            }
        }

        // 设置目标函数: maximize x + y
        model.setObjective(eta + 0, GRB_MINIMIZE);

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

        // 添加约束: 2-4
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay) {
                expr += r.at({i.id, j});
            }
            model.addConstr(expr <= static_cast<double>(this->C_D.size()) * h.at(i.id), "constr2-4_" + to_string(i.id));
        }

        // 添加约束: 2-5
        for (auto& i : this->C_D) {
            GRBLinExpr expr = 0;
            for (auto& j : C_delay) {
                expr += r.at({i.id, j});
            }
            model.addConstr(expr >= h.at(i.id), "constr2-4_" + to_string(i.id));
        }

        // 添加约束: 2-11
        for (auto& i : this->C_D) {
            for (auto& j : C_delay) {
                if (i.AP.count(j) == 0)
                    model.addConstr(r.at({i.id, j}) == 0);
            }
        }

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

        // 求解模型
        MyLazyCallback cb(r, h, eta, this->C_D, this->C_delay, this->n_1, this->info, this->P);
        model.setCallback(&cb);
        //model.write("model.lip");
        model.optimize();

        // 输出结果
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {

            //收集结果
            for (auto& i : this->C_D) {
                for (auto& j : this->C_delay) {
                    if (r.at({i.id, j}).get(GRB_DoubleAttr_X) >= 0.99) obj[i.id].goods.insert(j);
                }
                obj[i.id].time = cb.get_t(i.id);
            }
            this->dt = cb.get_dt();
            cout << "精确算法 = " << model.get(GRB_DoubleAttr_ObjVal) << endl;    
        }
        else {
            cout << "No optimal solution found." << endl;
        }

    }
    catch (GRBException& e) {
        cout << "Gurobi exception: " << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Unknown error during optimization." << endl;
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
            model.addConstr(expr <= static_cast<double>(this->C_D.size()) * h.at(i.id), "constr2-4_" + to_string(i.id));
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
                    model.addConstr(r.at({i.id, j}) == 0);
            }
        }

        //添加约束2-12
        for (auto& [i, j] : this->P) {
            model.addConstr(dt.at(i.first) - dt.at(i.second) <= j, "constr2-10" + to_string(i.first));
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