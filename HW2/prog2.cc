#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <tuple>
#include <utility>
#include <climits>
#include <set>
#include <map>
#include <random>
#include <chrono>
#include <fstream>

using namespace std;

// #define DEBUG // 取消註解以啟用調試輸出

// Node represents a coordinate (x, y)
typedef pair<int, int> Node;

// Structure to hold information about a connection (net)
struct Connect {
    int idx;
    Node src;
    Node dst;
    int mh_dist; // Manhattan distance
    long long bbox_size; // Bounding box area
    bool routed;
    int dist;
    int bend;
    vector<Node> path;

    Connect(int i, int a, int b, int x, int y) : 
        idx(i), src(a, b), dst(x, y), 
        mh_dist(abs(x - a) + abs(y - b)), 
        bbox_size((long long)(abs(x-a)+1) * (abs(y-b)+1)),
        routed(false), dist(0), bend(0) {}

    Connect() : idx(-1), src(0, 0), dst(0, 0), mh_dist(0), bbox_size(0), routed(false), dist(0), bend(0) {}
};

// Comparison for initial routing order: by MH distance, then bbox size
bool cmp_initial_order(const Connect& a, const Connect& b) {
    if (a.mh_dist != b.mh_dist) return a.mh_dist < b.mh_dist;
    if (a.bbox_size != b.bbox_size) return a.bbox_size < b.bbox_size;
    return a.idx < b.idx;
}

bool cmp_idx(const Connect& a, const Connect& b) {
    return a.idx < b.idx;
}

// State for the A* search algorithm
struct State {
    Node pos;
    int dir;
    int g_cost;
    int f_cost;
    vector<Node> path;

    State(Node p, int d, int g, int f, vector<Node> pa) : 
        pos(p), dir(d), g_cost(g), f_cost(f), path(pa) {}

    bool operator>(const State& other) const {
        if (f_cost != other.f_cost) return f_cost > other.f_cost;
        if (g_cost != other.g_cost) return g_cost > other.g_cost;
        return path.size() > other.path.size();
    }
};

const vector<pair<int, int>> DIR = {{ -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 }};

int heuristic(const Node& a, const Node& b) {
    return abs(a.first - b.first) + abs(a.second - b.second);
}

// A* search algorithm with simplified cost
void A_star_search(vector<vector<int>>& grid, Connect& c, int W, int H, int BEND_PENALTY, int EXPANSION_PENALTY, const map<Node, vector<int>>& pin_map) {
    vector<vector<vector<int>>> cost(H + 1, vector<vector<int>>(W + 1, vector<int>(4, INT_MAX)));
    priority_queue<State, vector<State>, greater<State>> pq;

    vector<Node> initial_path;
    initial_path.push_back(c.src);
    pq.push(State(c.src, -1, 0, heuristic(c.src, c.dst), initial_path));

    Connect best_result;
    best_result.dist = INT_MAX;
    best_result.bend = INT_MAX;

    while (!pq.empty()) {
        State current = pq.top();
        pq.pop();
        Node pos = current.pos;

        if (current.dir != -1 && current.g_cost >= cost[pos.second][pos.first][current.dir]) continue;
        if (current.dir != -1) cost[pos.second][pos.first][current.dir] = current.g_cost;

        if (pos == c.dst) {
            int current_dist = current.path.size() - 1;
            int current_bends = 0;
            for (size_t i = 1; i < current.path.size() - 1; ++i) {
                bool prev_horizontal = current.path[i-1].second == current.path[i].second;
                bool next_horizontal = current.path[i].second == current.path[i+1].second;
                if (prev_horizontal != next_horizontal) {
                    current_bends++;
                }
            }
            
            if (current_bends < best_result.bend || 
                (current_bends == best_result.bend && current_dist < best_result.dist)) {
                best_result.routed = true;
                best_result.dist = current_dist;
                best_result.bend = current_bends;
                best_result.path = current.path;
            }
            continue;
        }
        
        for (int d = 0; d < 4; d++) {
            Node next_pos = {pos.first + DIR[d].first, pos.second + DIR[d].second};

            if (next_pos.first < 1 || next_pos.first > W || next_pos.second < 1 || next_pos.second > H) continue;
            
            int cell_val = grid[next_pos.second][next_pos.first];
            if (cell_val == -1 || (cell_val > 0 && !(next_pos == c.src && grid[next_pos.second][next_pos.first] == c.idx) && !(next_pos == c.dst && grid[next_pos.second][next_pos.first] == c.idx))) continue;

            int is_bend = (current.dir != -1 && current.dir != d);
            int expansion_cost = 0;
            if (next_pos.first < min(c.src.first, c.dst.first) || 
                next_pos.first > max(c.src.first, c.dst.first) || 
                next_pos.second < min(c.src.second, c.dst.second) || 
                next_pos.second > max(c.src.second, c.dst.second)) {
                expansion_cost = EXPANSION_PENALTY;
            }

            int congestion_cost = 0;
            for (const auto& d2 : DIR) {
                Node adj = {next_pos.first + d2.first, next_pos.second + d2.second};
                if (adj.first >= 1 && adj.first <= W && adj.second >= 1 && adj.second <= H && grid[adj.second][adj.first] > 0) {
                    congestion_cost += 5;
                }
            }

            int new_g_cost = current.g_cost + 1 + is_bend * BEND_PENALTY + expansion_cost + congestion_cost;

            if (new_g_cost < cost[next_pos.second][next_pos.first][d]) {
                int h_cost = heuristic(next_pos, c.dst);
                vector<Node> new_path = current.path;
                new_path.push_back(next_pos);
                pq.push(State(next_pos, d, new_g_cost, new_g_cost + h_cost, new_path));
            }
        }
    }
    
    if (best_result.dist != INT_MAX) {
        c.routed = true;
        c.dist = best_result.dist;
        c.bend = best_result.bend;
        c.path = best_result.path;
    } else {
        c.routed = false;
        c.path.clear();
#ifdef DEBUG
        cerr << "Failed to route net " << c.idx << ": src=(" << c.src.first << "," << c.src.second 
             << "), dst=(" << c.dst.first << "," << c.dst.second << ")\n";
#endif
    }
}

void clear_route_from_grid(vector<vector<int>>& grid, Connect& c) {
    if (!c.routed) return;
    for (const auto& p : c.path) {
        if (grid[p.second][p.first] == c.idx) {
            grid[p.second][p.first] = 0;
        }
    }
    c.routed = false;
    c.path.clear();
}

void apply_route_to_grid(vector<vector<int>>& grid, Connect& c) {
    if (!c.routed) return;
    // Check for conflicts
    for (const auto& p : c.path) {
        if (grid[p.second][p.first] != 0 && grid[p.second][p.first] != c.idx) {
#ifdef DEBUG
            cerr << "Conflict at (" << p.first << "," << p.second << ") between net " << c.idx 
                 << " and net " << grid[p.second][p.first] << endl;
#endif
            c.routed = false;
            c.path.clear();
            return;
        }
    }
    // Apply route
    for (const auto& p : c.path) {
        grid[p.second][p.first] = c.idx;
    }
}

// Compute solution score
double compute_score(const vector<Connect>& connects) {
    int total_route = 0;
    long long total_dist = 0;
    long long total_bend = 0;
    for (const auto& c : connects) {
        if (c.routed) {
            total_route++;
            total_dist += c.dist;
            total_bend += c.bend;
        }
    }
    return 1000.0 * total_route - total_dist - 10.0 * total_bend;
}

// Simulated Annealing for Rip-up and Reroute
void simulated_annealing(vector<vector<int>>& grid, vector<Connect>& connects, vector<pair<Node, Node>> blockages, int W, int H, int BEND_PENALTY, int EXPANSION_PENALTY, const map<Node, vector<int>>& pin_map, int max_inner_iterations) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    uniform_int_distribution<> net_dis(0, connects.size() - 1);

    double T = 1000.0; // Initial temperature
    double T_min = 0.1; // Minimum temperature
    double alpha = 0.95; // Cooling rate

    vector<Connect> best_connects = connects;
    double best_score = compute_score(best_connects);

    // Identify high-risk nets with shared endpoints
    vector<int> high_risk_nets;
    for (int i = 0; i < connects.size(); ++i) {
        if (pin_map.at(connects[i].src).size() > 1 || pin_map.at(connects[i].dst).size() > 1) {
            high_risk_nets.push_back(i);
        }
    }

    while (T > T_min) {
        for (int inner = 0; inner < max_inner_iterations; ++inner) {
            vector<Connect> new_connects = connects;
            vector<vector<int>> new_grid = grid;

            // Randomly select 1-5 nets to rip-up
            int num_to_rip = uniform_int_distribution<>(1, min(5, (int)connects.size()))(gen);
            vector<int> rip_indices;
            for (int i = 0; i < num_to_rip; ++i) {
                int idx;
                if (!high_risk_nets.empty() && dis(gen) < 0.5) {
                    idx = high_risk_nets[net_dis(gen) % high_risk_nets.size()];
                } else {
                    idx = net_dis(gen);
                }
                if (find(rip_indices.begin(), rip_indices.end(), idx) == rip_indices.end()) {
                    rip_indices.push_back(idx);
                }
            }

            // Clear selected nets
            for (int idx : rip_indices) {
                clear_route_from_grid(new_grid, new_connects[idx]);
            }

            // Reapply pins for unrouted nets
            map<Node, int> used_pins;
            for (auto& c : new_connects) {
                if (!c.routed) {
                    if (!used_pins.count(c.src)) {
                        new_grid[c.src.second][c.src.first] = c.idx;
                        used_pins[c.src] = c.idx;
                    }
                    if (!used_pins.count(c.dst)) {
                        new_grid[c.dst.second][c.dst.first] = c.idx;
                        used_pins[c.dst] = c.idx;
                    }
                }
            }

            // Shuffle rip_indices for rerouting order
            shuffle(rip_indices.begin(), rip_indices.end(), gen);

            // Reroute selected nets
            for (int idx : rip_indices) {
                A_star_search(new_grid, new_connects[idx], W, H, BEND_PENALTY, EXPANSION_PENALTY, pin_map);
                if (new_connects[idx].routed) {
                    apply_route_to_grid(new_grid, new_connects[idx]);
                }
            }

            // Compute new score
            double new_score = compute_score(new_connects);
            double delta = new_score - compute_score(connects);

            // Accept new solution
            if (delta > 0 || dis(gen) < exp(delta / T)) {
                connects = new_connects;
                grid = new_grid;
                if (new_score > best_score) {
                    best_connects = new_connects;
                    best_score = new_score;
                }
            }
        }
        T *= alpha;
    }

    // Restore best solution
    connects = best_connects;
    grid = vector<vector<int>>(H + 1, vector<int>(W + 1, 0));
    // Reapply blockages
    for (const auto& block : blockages) {
        for (int i = block.first.first; i <= block.second.first; i++) {
            for (int j = block.first.second; j <= block.second.second; j++) {
                if (i >= 1 && i <= W && j >= 1 && j <= H) grid[j][i] = -1;
            }
        }
    }
    // Reapply routes
    map<Node, int> used_pins;
    for (auto& c : connects) {
        if (c.routed) {
            apply_route_to_grid(grid, c);
            used_pins[c.src] = c.idx;
            used_pins[c.dst] = c.idx;
        }
    }
}

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);

    // Check command-line arguments
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " [input_file] [output_file]\n";
        return 1;
    }

    // Open input file
    ifstream in_file(argv[1]);
    if (!in_file.is_open()) {
        cerr << "Error: Cannot open input file " << argv[1] << "\n";
        return 1;
    }

    // Open output file
    ofstream out_file(argv[2]);
    if (!out_file.is_open()) {
        cerr << "Error: Cannot open output file " << argv[2] << "\n";
        in_file.close();
        return 1;
    }

    int W, H, num_block, num_connect;
    string dump;

    // Read grid dimensions
    in_file >> W >> H;
    vector<vector<int>> grid(H + 1, vector<int>(W + 1, 0));
    vector<pair<Node, Node>> blockages;

    // Dynamic parameter tuning
    int BEND_PENALTY = 3;
    int EXPANSION_PENALTY = 0;
    int MAX_ITERATIONS = 20;
    int max_inner_iterations = 50;

    // Read blockages
    in_file >> dump >> num_block;
#ifdef DEBUG
    cerr << "Reading " << num_block << " blockages\n";
#endif
    while (num_block--) {
        int a, b, x, y;
        in_file >> a >> b >> x >> y;
        blockages.emplace_back(Node(a, b), Node(x, y));
        for (int i = a; i <= x; i++) {
            for (int j = b; j <= y; j++) {
                if (i >= 1 && i <= W && j >= 1 && j <= H) grid[j][i] = -1;
            }
        }
    }

    // Read connections
    in_file >> dump >> num_connect;
#ifdef DEBUG
    cerr << "Reading " << num_connect << " connections\n";
#endif
    vector<Connect> connects(num_connect);
    map<Node, vector<int>> pin_map;
    for (int i = 0; i < num_connect; i++) {
        int a, b, x, y;
        in_file >> a >> b >> x >> y;
        connects[i] = Connect(i + 1, a, b, x, y);
        pin_map[Node(a, b)].push_back(i + 1);
        pin_map[Node(x, y)].push_back(i + 1);
    }

    // Adjust parameters based on connections
    if (num_connect > 50) {
        BEND_PENALTY = 1;
        EXPANSION_PENALTY = 2;
        MAX_ITERATIONS = 30;
        max_inner_iterations = 200;
    } else if (num_connect > 20) {
        BEND_PENALTY = 1;
        EXPANSION_PENALTY = 2;
        MAX_ITERATIONS = 25;
        max_inner_iterations = 100;
    }

    // Initial routing
    sort(connects.begin(), connects.end(), cmp_initial_order);
    map<Node, int> used_pins;
    for (auto& c : connects) {
        if (!used_pins.count(c.src) && !used_pins.count(c.dst)) {
            A_star_search(grid, c, W, H, BEND_PENALTY, EXPANSION_PENALTY, pin_map);
            if (c.routed) {
                apply_route_to_grid(grid, c);
                used_pins[c.src] = c.idx;
                used_pins[c.dst] = c.idx;
            }
        }
    }

    // Simulated Annealing
    simulated_annealing(grid, connects, blockages, W, H, BEND_PENALTY, EXPANSION_PENALTY, pin_map, max_inner_iterations);

    // Validate routes
    map<Node, vector<int>> point_to_nets;
    for (auto& c : connects) {
        if (c.routed) {
            for (const auto& p : c.path) {
                point_to_nets[p].push_back(c.idx);
            }
        }
    }
    for (const auto& [p, nets] : point_to_nets) {
        if (nets.size() > 1) {
#ifdef DEBUG
            cerr << "Overlap at (" << p.first << "," << p.second << ") among nets: ";
            for (int idx : nets) cerr << idx << " ";
            cerr << endl;
#endif
            for (int idx : nets) {
                connects[idx - 1].routed = false;
                connects[idx - 1].path.clear();
            }
        }
    }

    // Final Output
    sort(connects.begin(), connects.end(), cmp_idx);

    int total_route = 0;
    long long total_dist = 0;
    long long total_bend = 0;
    int max_idx = 0;
    int max_dist = 0;

    for (const auto& c : connects) {
        if (c.routed) {
            total_route++;
            total_dist += c.dist;
            total_bend += c.bend;
            if (c.dist > max_dist) {
                max_dist = c.dist;
                max_idx = c.idx;
            }
        }
    }

    out_file << "#interconnections routed = " << total_route << "\n";
    out_file << "Total interconnection length = " << total_dist << "\n";
    
    out_file << "The longest interconnection = " << max_idx << "; length = " << max_dist << "\n";
    
    out_file << "Total number of bends = " << total_bend << "\n";

    for (const auto& c : connects) {
        if (c.routed) {
            out_file << "Interconnection " << c.idx << ": length = " << c.dist
                     << ", #bends = " << c.bend << "\n";
            for (size_t j = 0; j < c.path.size(); j++) {
                if (j != 0) out_file << ", ";
                out_file << "(" << c.path[j].first << ", " << c.path[j].second << ")";
            }
            out_file << "\n";
        } else {
            out_file << "Interconnection " << c.idx << ": fails.\n";
        }
    }

    // Close files
    in_file.close();
    out_file.close();

    return 0;
}