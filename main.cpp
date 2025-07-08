/*
 * 行程规划程序
 * 使用图结构 + 贪心选择 + 分治策略 + 动态规划优化实现最优旅游路线规划
 * 核心功能：
 * 1. 多天行程规划
 * 2. 考虑景点开放时间限制
 * 3. 预算约束
 * 4. 距离和时间优化
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <climits>
#include <bitset>
using namespace std;

// 全局常量定义
const double SPEED = 30.0;         // 旅行速度，单位：公里/小时
const double COST_PER_KM = 2.3;    // 每公里交通费用，单位：元/公里
const double INF = 1e9;            // 无穷大，用于初始化

// 景点结构体：存储景点的基本信息
struct Spot {
    string name;           // 景点名称
    double ticket;        // 门票价格
    double play_time;     // 游玩所需时间（小时）
    double interest;      // 景点的期待值/吸引力
    double open_time;     // 开放时间（24小时制）
    double close_time;    // 关闭时间（24小时制）
    double x, y;          // 景点的地理坐标
    bool visited = false;  // 是否已经游玩过
};

// 每日时间窗口结构体
struct Day {
    double start_time;    // 一天的开始时间
    double end_time;      // 一天的结束时间
};

// 路线结构体：存储规划好的行程信息
struct Route {
    vector<int> spot_indices;  // 景点访问顺序（存储索引）
    double start_time;         // 行程开始时间
    double end_time;          // 行程结束时间
    double play_duration;      // 总游玩时间
    double travel_duration;    // 总交通时间
    double cost;              // 总花费
    double interest;          // 总期待值
};

/*
 * 计算两点之间的欧几里得距离
 * @param x1, y1: 第一个点的坐标
 * @param x2, y2: 第二个点的坐标
 * @return: 两点间的直线距离
 */
double euclidean_distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

/*
 * 根据距离计算交通时间
 * @param dist: 距离（公里）
 * @return: 所需时间（小时）
 */
double travel_time(double dist) {
    return dist / SPEED;
}

/*
 * 计算交通费用
 * @param dist: 距离（公里）
 * @return: 费用（元）
 */
double travel_cost(double dist) {
    return dist * COST_PER_KM;
}

// 全局距离矩阵
vector<vector<double>> global_dist;

/*
 * 初始化全局距离矩阵
 * 计算所有景点之间以及起点（酒店）到各景点的距离
 * @param spots: 景点列表
 */
void init_global_dist(const vector<Spot>& spots) {
    int N = spots.size();
    global_dist.assign(N + 1, vector<double>(N + 1));

    // 计算起点（索引0）到所有景点的距离
    for (int i = 0; i < N; ++i) {
        global_dist[0][i + 1] = global_dist[i + 1][0] =
            euclidean_distance(0, 0, spots[i].x, spots[i].y);
    }

    // 计算景点之间的距离
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            global_dist[i + 1][j + 1] =
                euclidean_distance(spots[i].x, spots[i].y, spots[j].x, spots[j].y);
        }
    }
}

/*
 * 检查给定路线是否满足时间窗口约束
 * @param spots: 景点列表
 * @param order: 访问顺序
 * @param day_start: 一天的开始时间
 * @param day_end: 一天的结束时间
 * @return: 是否可行
 */
bool check_time_windows(const vector<Spot>& spots, const vector<int>& order,
    double day_start, double day_end) {
    double current_time = day_start;
    double prev_x = 0, prev_y = 0;

    // 遍历每个景点，检查时间约束
    for (int idx : order) {
        double dist = euclidean_distance(prev_x, prev_y, spots[idx].x, spots[idx].y);
        current_time += travel_time(dist);

        // 检查是否能在景点关闭前到达
        if (current_time > spots[idx].close_time) {
            return false;
        }

        // 如果早到，需要等待开放
        if (current_time < spots[idx].open_time) {
            current_time = spots[idx].open_time;
        }

        current_time += spots[idx].play_time;

        // 检查是否能在景点关闭前游玩完
        if (current_time > spots[idx].close_time) {
            return false;
        }

        prev_x = spots[idx].x;
        prev_y = spots[idx].y;
    }

    // 检查是否能在规定时间内返回
    double dist_back = euclidean_distance(prev_x, prev_y, 0, 0);
    current_time += travel_time(dist_back);
    return current_time <= day_end;
}

/*
 * 使用最近邻算法构建初始路线
 * @param spots: 景点列表
 * @param indices: 可选景点的索引
 * @param day_start: 开始时间
 * @param day_end: 结束时间
 * @return: 可行的访问顺序
 */
vector<int> nearest_neighbor(const vector<Spot>& spots, const vector<int>& indices,
    double day_start, double day_end) {
    int n = indices.size();
    vector<bool> visited(n, false);
    vector<int> order;
    double current_time = day_start;
    double prev_x = 0, prev_y = 0;

    while (order.size() < n) {
        double min_dist = INF;
        int next = -1;
        double earliest_close = INF;

        // 寻找最近且满足时间窗口的未访问景点
        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                double dist = euclidean_distance(prev_x, prev_y,
                    spots[indices[i]].x, spots[indices[i]].y);
                double arrive = current_time + travel_time(dist);
                double finish = max(arrive, spots[indices[i]].open_time) +
                    spots[indices[i]].play_time;

                if (arrive <= spots[indices[i]].close_time && finish <= day_end) {
                    if (spots[indices[i]].close_time < earliest_close ||
                        (spots[indices[i]].close_time == earliest_close && dist < min_dist)) {
                        min_dist = dist;
                        earliest_close = spots[indices[i]].close_time;
                        next = i;
                    }
                }
            }
        }

        if (next == -1) break;

        visited[next] = true;
        order.push_back(indices[next]);

        double dist = euclidean_distance(prev_x, prev_y,
            spots[indices[next]].x, spots[indices[next]].y);
        current_time = max(current_time + travel_time(dist),
            spots[indices[next]].open_time) + spots[indices[next]].play_time;
        prev_x = spots[indices[next]].x;
        prev_y = spots[indices[next]].y;
    }

    // 验证并优化路线
    if (!order.empty() && !check_time_windows(spots, order, day_start, day_end)) {
        for (int i = order.size(); i > 0; --i) {
            vector<int> partial_order(order.begin(), order.begin() + i);
            if (check_time_windows(spots, partial_order, day_start, day_end)) {
                return partial_order;
            }
        }
        return {};
    }

    return order;
}

/*
 * 核心路线规划函数
 * 使用动态规划和贪心策略进行多天行程规划
 * @param N: 景点数量
 * @param D: 旅行天数
 * @param total_budget: 总预算
 * @param spots: 景点列表
 * @param days: 每天的时间窗口
 * @return: 规划好的路线列表
 */
vector<Route> plan_routes(int N, int D, double total_budget,
    vector<Spot>& spots, vector<Day>& days) {
    vector<Route> routes(D);
    double budget_used = 0.0;

    // 初始化距离矩阵
    init_global_dist(spots);

    // 对每一天进行规划
    for (int d = 0; d < D; ++d) {
        double day_start = days[d].start_time, day_end = days[d].end_time;

        // 筛选当天可行的景点
        vector<int> valid_indices;
        for (int i = 0; i < N; ++i) {
            if (!spots[i].visited &&
                spots[i].ticket <= total_budget - budget_used &&
                spots[i].open_time <= day_end &&
                spots[i].close_time >= day_start) {
                valid_indices.push_back(i);
            }
        }

        int M = valid_indices.size();
        if (M == 0) continue;

        // 预处理状态信息
        vector<double> ticket_sum(1 << M, 0), interest_sum(1 << M, 0),
            play_sum(1 << M, 0);
        for (int mask = 1; mask < (1 << M); ++mask) {
            for (int i = 0; i < M; ++i) {
                if (mask & (1 << i)) {
                    ticket_sum[mask] += spots[valid_indices[i]].ticket;
                    interest_sum[mask] += spots[valid_indices[i]].interest;
                    play_sum[mask] += spots[valid_indices[i]].play_time;
                }
            }
            if (ticket_sum[mask] + budget_used > total_budget) {
                ticket_sum[mask] = INF;
                interest_sum[mask] = -1;
                play_sum[mask] = INF;
            }
        }

        // 动态规划状态数组
        vector<vector<double>> dp(1 << M, vector<double>(M, INF));
        vector<vector<int>> parent(1 << M, vector<int>(M, -1));

        // 计算每个景点的最早到达和最晚离开时间
        vector<double> earliest_arrival(M), latest_departure(M);
        for (int i = 0; i < M; ++i) {
            earliest_arrival[i] = day_start +
                travel_time(global_dist[0][valid_indices[i] + 1]);
            latest_departure[i] = spots[valid_indices[i]].close_time -
                spots[valid_indices[i]].play_time;
        }

        // 初始化单个景点的状态
        for (int i = 0; i < M; ++i) {
            double arrive = earliest_arrival[i];
            double start_play = max(arrive, spots[valid_indices[i]].open_time);
            double finish_play = start_play + spots[valid_indices[i]].play_time;
            double travel_cost_to_spot = travel_cost(global_dist[0][valid_indices[i] + 1]);
            double travel_cost_back = travel_cost(global_dist[valid_indices[i] + 1][0]);

            if (finish_play <= spots[valid_indices[i]].close_time &&
                finish_play + travel_time(global_dist[valid_indices[i] + 1][0]) <= day_end &&
                ticket_sum[1 << i] + travel_cost_to_spot + travel_cost_back <=
                total_budget - budget_used) {
                dp[1 << i][i] = finish_play;
                parent[1 << i][i] = -1;
            }
        }

        // 动态规划主循环
        const int MAX_K = 5;  // 每天最多访问的景点数
        for (int k = 1; k <= min(M, MAX_K); ++k) {
            for (int mask = 1; mask < (1 << M); ++mask) {
                if (bitset<32>(mask).count() != k || interest_sum[mask] < 0) continue;

                // 检查当前组合是否可行
                vector<int> subset;
                for (int i = 0; i < M; ++i) {
                    if (mask & (1 << i)) subset.push_back(valid_indices[i]);
                }
                if (!check_time_windows(spots, subset, day_start, day_end)) {
                    continue;
                }

                // 尝试添加新的景点
                for (int u = 0; u < M; ++u) {
                    if (!(mask & (1 << u)) || dp[mask][u] == INF) continue;

                    for (int v = 0; v < M; ++v) {
                        if (mask & (1 << v)) continue;

                        int next_mask = mask | (1 << v);
                        double travel_t = travel_time(
                            global_dist[valid_indices[u] + 1][valid_indices[v] + 1]);
                        double arrive = dp[mask][u] + travel_t;

                        if (arrive > spots[valid_indices[v]].close_time) continue;

                        double start_play = max(arrive, spots[valid_indices[v]].open_time);
                        double finish_play = start_play + spots[valid_indices[v]].play_time;
                        double travel_cost_to_spot = travel_cost(
                            global_dist[0][valid_indices[u] + 1]);
                        double travel_cost_next = travel_cost(
                            global_dist[valid_indices[u] + 1][valid_indices[v] + 1]);
                        double travel_cost_back = travel_cost(
                            global_dist[valid_indices[v] + 1][0]);
                        double cost = ticket_sum[next_mask] + travel_cost_to_spot +
                            travel_cost_next + travel_cost_back;

                        if (finish_play <= day_end && cost + budget_used <= total_budget &&
                            finish_play < dp[next_mask][v]) {
                            dp[next_mask][v] = finish_play;
                            parent[next_mask][v] = u;
                        }
                    }
                }
            }
        }

        // 选择最优解
        double best_interest = -1, best_cost = INF, best_end_time = INF;
        int best_mask = 0, best_last = -1;

        for (int mask = 1; mask < (1 << M); ++mask) {
            if (interest_sum[mask] < 0 || bitset<32>(mask).count() > MAX_K) continue;

            vector<int> subset;
            for (int i = 0; i < M; ++i) {
                if (mask & (1 << i)) subset.push_back(valid_indices[i]);
            }

            if (!check_time_windows(spots, subset, day_start, day_end)) continue;

            for (int last = 0; last < M; ++last) {
                if (dp[mask][last] == INF) continue;

                double return_travel_t = travel_time(
                    global_dist[valid_indices[last] + 1][0]);
                double finish_time = dp[mask][last] + return_travel_t;

                if (finish_time <= day_end) {
                    double travel_cost_to_spot = travel_cost(
                        global_dist[0][valid_indices[0] + 1]);
                    double travel_cost_back = travel_cost(
                        global_dist[valid_indices[last] + 1][0]);
                    double cost = ticket_sum[mask];

                    for (int i = 1; i < subset.size(); ++i) {
                        cost += travel_cost(
                            global_dist[subset[i - 1] + 1][subset[i] + 1]);
                    }
                    cost += travel_cost_to_spot + travel_cost_back;

                    if (cost + budget_used <= total_budget &&
                        interest_sum[mask] > best_interest) {
                        best_interest = interest_sum[mask];
                        best_cost = cost;
                        best_end_time = finish_time;
                        best_mask = mask;
                        best_last = last;
                    }
                }
            }
        }

        // 如果找不到可行解，继续下一天
        if (best_interest < 0) {
            continue;
        }

        // 重建最优路径
        vector<int> route_spots;
        int cur_mask = best_mask, cur_pos = best_last;
        while (cur_pos != -1) {
            route_spots.push_back(valid_indices[cur_pos]);
            int prev_pos = parent[cur_mask][cur_pos];
            cur_mask ^= (1 << cur_pos);
            cur_pos = prev_pos;
        }
        reverse(route_spots.begin(), route_spots.end());

        // 使用最近邻算法优化路线
        route_spots = nearest_neighbor(spots, route_spots, day_start, day_end);
        if (route_spots.empty()) {
            continue;
        }

        // 标记已访问的景点
        for (int idx : route_spots) {
            spots[idx].visited = true;
        }

        // 计算路线详细信息
        double play_duration = 0, travel_duration = 0;
        double current_time = day_start, prev_x = 0, prev_y = 0;
        double cost = 0;

        for (int idx : route_spots) {
            double dist = euclidean_distance(prev_x, prev_y, spots[idx].x, spots[idx].y);
            travel_duration += travel_time(dist);
            current_time = max(current_time + travel_time(dist), spots[idx].open_time);
            play_duration += spots[idx].play_time;
            current_time += spots[idx].play_time;
            cost += spots[idx].ticket;
            if (prev_x != 0 || prev_y != 0) {
                cost += travel_cost(dist);
            }
            prev_x = spots[idx].x;
            prev_y = spots[idx].y;
        }

        travel_duration += travel_time(euclidean_distance(prev_x, prev_y, 0, 0));
        cost += travel_cost(euclidean_distance(0, 0,
            spots[route_spots[0]].x, spots[route_spots[0]].y));
        cost += travel_cost(euclidean_distance(prev_x, prev_y, 0, 0));
        best_end_time = current_time +
            travel_time(euclidean_distance(prev_x, prev_y, 0, 0));

        // 保存当天的路线
        routes[d] = { route_spots, day_start, best_end_time, play_duration,
                     travel_duration, cost, best_interest };
        budget_used += cost;
    }

    return routes;
}

/*
 * 打印规划结果
 * @param routes: 规划好的路线列表
 * @param spots: 景点列表
 */
void print_routes(const vector<Route>& routes, const vector<Spot>& spots) {
    double total_cost = 0, total_interest = 0;

    cout << fixed << setprecision(2);
    cout << "\n========== 行程规划结果 ==========" << endl;

    for (int d = 0; d < routes.size(); ++d) {
        cout << "第" << d + 1 << "天行程：" << endl;

        if (routes[d].spot_indices.empty()) {
            cout << "无可行景点行程" << endl;
        }
        else {
            cout << "行程路线：酒店 -> ";
            for (int i = 0; i < routes[d].spot_indices.size(); ++i) {
                cout << spots[routes[d].spot_indices[i]].name;
                if (i != routes[d].spot_indices.size() - 1) {
                    cout << " -> ";
                }
            }
            cout << " -> 酒店" << endl;

            cout << "起始时间：" << routes[d].start_time << " 时" << endl;
            cout << "结束时间：" << routes[d].end_time << " 时" << endl;
            cout << "游玩时间：" << routes[d].play_duration << " 小时" << endl;
            cout << "交通时间：" << routes[d].travel_duration << " 小时" << endl;
            cout << "当日花费：" << routes[d].cost << " 元" << endl;
            cout << "期待总值：" << routes[d].interest << endl;

            total_cost += routes[d].cost;
            total_interest += routes[d].interest;
        }
        cout << "----------------------------------" << endl;
    }

    cout << "\n总预算消耗：" << total_cost << " 元" << endl;
    cout << "总期待程度：" << total_interest << endl;
    cout << "==================================" << endl;
}

/*
 * 主函数
 * 处理输入并调用规划函数
 */
int main() {
    int N, D;
    double B;

    cout << "请输入景点数量 N、旅行天数 D、总预算 B（元）：\n";
    cin >> N >> D >> B;

    vector<Spot> spots(N);
    cout << "请输入每个景点的信息（名称 票价 游玩时间 期待 开始时间 结束时间 x y）：\n";
    for (int i = 0; i < N; ++i) {
        cin >> spots[i].name >> spots[i].ticket >> spots[i].play_time
            >> spots[i].interest >> spots[i].open_time >> spots[i].close_time
            >> spots[i].x >> spots[i].y;
    }

    vector<Day> days(D);
    cout << "请输入每一天的可用时间（开始时间 结束时间）：\n";
    for (int i = 0; i < D; ++i) {
        cin >> days[i].start_time >> days[i].end_time;
    }

    // 调用路线规划函数
    vector<Route> routes = plan_routes(N, D, B, spots, days);

    // 打印规划结果
    print_routes(routes, spots);

    return 0;
}
