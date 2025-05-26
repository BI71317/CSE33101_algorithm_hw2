#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <chrono> 
#include <cmath>

#define MAX_N 100000LL
#define INF 1e300
#define MAX_HEAP_SIZE MAX_N * MAX_N

using namespace std;

int path[MAX_N + 20];

struct City {
    int id;
    double x, y;
};

struct HeapNode {
    int u, v;
    int weight;
};

bool comp(const HeapNode& a, const HeapNode& b) {
    return a.weight < b.weight;
}

struct MinHeap {
    HeapNode* data;
    int capacity;
    int size;

    void init(int max_size) {
        data = new (nothrow) HeapNode[max_size];
        if (!data) {
            fprintf(stderr, "[Heap] allocation failed for size: %d\n", max_size);
            exit(1);
        }
        capacity = max_size;
        size = 0;
    }
    void destroy() {
        delete[] data;
    }

    void push(int u, int v, int weight) {
        int i = size++;
        data[i] = {u, v, weight};

        while (i > 0) {
            int parent = (i - 1) / 2;
            if (comp(data[i], data[parent])) {
                HeapNode tmp = data[i];
                data[i] = data[parent];
                data[parent] = tmp;
                i = parent;
            } else break;
        }
    }

    HeapNode pop() {
        HeapNode ret = data[0];
        data[0] = data[--size];

        int i = 0;
        while (true) {
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            int smallest = i;

            if (left < size && comp(data[left], data[smallest])) smallest = left;
            if (right < size && comp(data[right], data[smallest])) smallest = right;

            if (smallest == i) break;

            HeapNode tmp = data[i];
            data[i] = data[smallest];
            data[smallest] = tmp;
            i = smallest;
        }
        return ret;
    }

    bool empty() {
        return size == 0;
    }
};

City cities[MAX_N + 20];
int city_count = 0;

int dist(int i, int j) {
    double dx = cities[i].x - cities[j].x;
    double dy = cities[i].y - cities[j].y;
    return (int)(round(sqrt(dx * dx + dy * dy)));
}


bool starts_with(const std::string& line, const char* prefix) {
    return line.find(prefix) == 0;
}

bool load_tsp_file(const char* filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "file open fail: " << filename << std::endl;
        return false;
    }

    int MAXDEPT = 10;
    string line;
    while (getline(fin, line)) {
        if (starts_with(line, "DIMENSION")) {
            sscanf(line.c_str(), "DIMENSION : %d", &city_count);
        }
        else if (starts_with(line, "NODE_COORD_SECTION")) {
            break;
        }
        MAXDEPT--;
        if (!MAXDEPT) {
            cerr << "file structure err\n";
            return false;
        }
    }

    if (city_count == 0) {
        cerr << "No dimension given\n";
        return false;
    }

    for (int i = 0; i < city_count; i++) {
        int id;
        double x, y;
        fin >> id >> x >> y;
        if (fin.fail()) {
            cerr << "parsing err\n";
            return false;
        }
        cities[i] = {id, x, y};
    }

    return true;
}

template <typename Func>
void run_with_timer(const string& name, Func algorithm) {
    cout << "[" << name << "] ...\n";
    auto start = chrono::steady_clock::now();

    algorithm(); 

    auto end = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "[" << name << "] done (run time: " << duration << "ms)\n";
}

struct Edge {
    int u, v;      
    int weight;     
};

void compute_mst(int** dist_matrix, int n, int *mst_parent){
    // Prim algorithm
    // starting from 0, using heap, find path
    bool *in_mst = new bool[n]();
    MinHeap heap;
    heap.init(n * n + 500);
    in_mst[0] = true;
    mst_parent[0] = -1;
    for (int v = 1; v < n; ++v) {
        heap.push(0, v, dist_matrix[0][v]);
    }
    while (!heap.empty()) {
        HeapNode e = heap.pop();
        int u = e.u;
        int v = e.v;
        if (in_mst[v]) continue;
        mst_parent[v] = u;
        in_mst[v] = true;
        for (int w = 0; w < n; ++w) {
            if (!in_mst[w]) {
                heap.push(v, w, dist_matrix[v][w]);
            }
        }
    }
}

void find_odd_vertices(int *mst_parent, int n, int *degree, int *odd, int *odd_cnt){
    // from given MST graph, find odd vertices by traversing mst tree
    for(int i = 1; i < n; ++i){
        int u = mst_parent[i];
        int v = i;
        degree[u]++;
        degree[v]++;
    }
    *odd_cnt = 0;
    for (int i = 0; i < n; ++i) {
        if (degree[i] % 2 == 1) {
            odd[(*odd_cnt)++] = i;
        }
    }
}

inline int min(int a, int b) { return a < b ? a : b; }

int find_lca(int* match, int* base, int* parent, int u, int v, int n, bool* in_path) {
    for (int i = 0; i < n; ++i) in_path[i] = false;
    while (true) {
        u = base[u];
        in_path[u] = true;
        if (match[u] == -1) break;
        u = parent[match[u]];
    }
    while (true) {
        v = base[v];
        if (in_path[v]) return v;
        if (match[v] == -1) break;
        v = parent[match[v]];
    }
    return -1;
}

void mark_blossom(int* base, int* parent, bool* in_blossom, int* label, int* queue, int& qt, int lca, int u, int v, int* match, int n) {
    auto blossom_trace = [&](int x) {
        while (base[x] != lca) {
            int b = base[x];
            int m = match[b];

            in_blossom[b] = true;
            in_blossom[m] = true;

            base[b] = lca;
            base[m] = lca;

            if (parent[b] == -1) parent[b] = m;

            x = parent[m];
            if (x == -1) break;
        }
    };

    blossom_trace(u);
    blossom_trace(v);

    for (int i = 0; i < n; ++i) {
        if (in_blossom[base[i]]) {
            base[i] = lca;
            if (label[i] != 1) {  
                label[i] = 1;
                queue[qt++] = i;
            }
        }
    }
}
void expand_and_flip(int* match, int* parent, int u, int v, int* base, int n) {
    int cur = v;
    bool* visited = new bool[n](); // defensive
    while (cur != -1) {
        int next = parent[cur];
        if (next == -1 || visited[cur]) break;
        visited[cur] = true;
        match[cur] = next;
        match[next] = cur;
        cur = parent[next];
    }
    delete[] visited;
}

void edmond_matching(int** dist_matrix, int* odd, int odd_cnt, Edge* matching) {
    const int n = odd_cnt;
    int* match = new int[n];
    int* base = new int[n];
    int* parent = new int[n];
    int* label = new int[n];
    int* queue = new int[n];
    int* dual = new int[n];
    bool* in_path = new bool[n];
    bool* in_blossom = new bool[n];

    for (int i = 0; i < n; ++i) match[i] = -1;
    for (int i = 0; i < n; ++i) {
        dual[i] = 1000000000;
        for (int j = 0; j < n; ++j)
            if (i != j && dist_matrix[odd[i]][odd[j]] < dual[i])
                dual[i] = dist_matrix[odd[i]][odd[j]];
    }

    cout << "[Christofides] BFS\n";
    for (int root = 0; root < n; ++root) {
        if (match[root] != -1) continue;
        cout << "[Christofides] BFSing from " << root << " to " << n << '\n';
    
        for (int i = 0; i < n; ++i) {
            parent[i] = -1;
            base[i] = i;
            label[i] = 0;
        }

        int qh = 0, qt = 0;
        queue[qt++] = root;
        label[root] = 1;

        bool found = false;
        while (qh < qt && !found) {
            int u = queue[qh++];
            for (int v = 0; v < n; ++v) {
                if (base[u] == base[v] || match[u] == v) continue;
                if (dist_matrix[odd[u]][odd[v]] > dual[u] + dual[v]) continue;

                if (label[v] == 0) {

                    label[v] = 2;
                    parent[v] = u;
                    if (match[v] == -1) {
                        cout << "flippin\n";
                        expand_and_flip(match, parent, root, v, base, n);
                        found = true;
                        break;
                    }
                    label[match[v]] = 1;
                    queue[qt++] = match[v];
                    parent[match[v]] = v;
                } 
                else if (label[v] == 1) {
                        cout << "cntractin\n";
                    int lca = find_lca(match, base, parent, u, v, n, in_path);
                    for (int i = 0; i < n; ++i) in_blossom[i] = false;
                    mark_blossom(base, parent, in_blossom, label, queue, qt, lca, u, v, match, n);
                }
            }
        }

        if (!found) {
            int delta = 1000000000;
            for (int i = 0; i < n; ++i) if (label[i] == 1)
                for (int j = 0; j < n; ++j) if (label[j] == 0)
                    delta = min(delta, dist_matrix[odd[i]][odd[j]] - dual[i] - dual[j]);
            for (int i = 0; i < n; ++i) {
                if (label[i] == 1) dual[i] += delta;
                else if (label[i] == 2) dual[i] -= delta;
            }
        }
    }

    cout << "[Christofides] BFS done\n";
    int match_cnt = 0;
    for (int i = 0; i < n; ++i) {
        if (match[i] != -1 && i < match[i]) {
            matching[match_cnt++] = { odd[i], odd[match[i]], dist_matrix[odd[i]][odd[match[i]]] };
        }
    }

    delete[] match;
    delete[] base;
    delete[] parent;
    delete[] label;
    delete[] queue;
    delete[] dual;
    delete[] in_path;
    delete[] in_blossom;
}

void build_multigraph(int n, int* mst_parent, Edge* matching, int match_cnt, int** multigraph) {
    // using adjacent matrix, represent connected edge
    for (int i = 1; i < n; ++i) {
        int u = i;
        int v = mst_parent[i];
        multigraph[u][v]++;
        multigraph[v][u]++;
    }

    for (int i = 0; i < match_cnt; ++i) {
        int u = matching[i].u;
        int v = matching[i].v;
        multigraph[u][v]++;
        multigraph[v][u]++;
    }
}

void find_eulerian_tour(int n, int** multigraph, int* path, int* path_len) {
    // using given adjacent matrix, find path greedy

    int** temp_graph = new int*[n];
    for (int i = 0; i < n; ++i) {
        temp_graph[i] = new int[n];
        for (int j = 0; j < n; ++j)
            temp_graph[i][j] = multigraph[i][j];
    }

    int stack[MAX_N * 2];
    int top = 0, idx = 0;

    stack[top++] = 0; 

    while (top > 0) {
        int u = stack[top - 1];
        int found = 0;
        for (int v = 0; v < n; ++v) {
            if (temp_graph[u][v] > 0) {
                temp_graph[u][v]--;
                temp_graph[v][u]--;
                stack[top++] = v;
                found = 1;
                break;
            }
        }
        if (!found) {
            path[idx++] = u;
            top--;
        }
    }
    *path_len = idx;
    for (int i = 0; i < n; ++i)
        delete[] temp_graph[i];
    delete[] temp_graph;
}

void shortcut_tour(int* euler_path, int euler_len, int* visited, int* tsp_path, int* tsp_len) {
    // remove duplicated node from euler path.
    int idx = 0;
    for (int i = 0; i < euler_len; ++i) {
        int v = euler_path[i];
        if (!visited[v]) {
            visited[v] = 1;
            tsp_path[idx++] = v;
        }
    }
    *tsp_len = idx;
}

void build_2appgraph(int n, int* mst_parent, int** multigraph){
    for (int i = 1; i < n; ++i) {
        int u = i;
        int v = mst_parent[i];
        multigraph[u][v] += 2;
        multigraph[v][u] += 2;
    }
}

void run_held_karp() {
    int n = city_count;
    size_t subsets = 1 << n;
    size_t total_bytes = subsets * n * sizeof(double);

    const size_t MEMORY_LIMIT = 1ULL << 30; // 1 GB
    if (n > 30 || total_bytes > MEMORY_LIMIT) {
        cerr << "[Held-Karp] required too much memory \n";
        return;
    }
    double** dp = new (nothrow) double*[subsets];
    int** parent = new (nothrow) int*[subsets];

    if (!dp || !parent) {
        cerr << "[Held-Karp] allocation fail\n";
        return;
    }
    for (size_t i = 0; i < subsets; ++i) {
        dp[i] = new (nothrow) double[n];
        parent[i] = new (nothrow) int[n];
        if (!dp[i] || !parent[i]) {
            cerr << "[Held-Karp] allocation fail at row " << i << endl;
            for (size_t j = 0; j <= i; ++j) {
                delete[] dp[j];
                delete[] parent[j];
            }
            delete[] dp;
            delete[] parent;
            return;
        }
    }
    cout << "[Held-Karp] memory allocation (" << (total_bytes >> 20) << " MB)\n";
 
    
    for(int mask = 0; mask < (1 << n); ++mask){
        for (int i = 0; i < n; ++i){
            dp[mask][i] = INF;
            parent[mask][i] = -1;
        }
    }
    dp[1][0] = 0; // starting at zero

    for(int mask = 0; mask < (1 << n); ++mask){
        // starting node i
        for(int i = 0; i < n; i++){
            int prev_mask = mask ^ (1 << i);
            // node visited right before Last node
            if(!(mask & (1 << i)))  continue;
            for(int j = 0; j < n; j++){
                // Last visited node 
                if(i == j || !(mask & (1 << i))) continue;
                double new_cost = dp[prev_mask][j] + dist(j, i);                
                if (new_cost < dp[mask][i]) {
                    dp[mask][i] = new_cost;
                    parent[mask][i] = j;
                }
            }
        }
    }
    double cost = INF;
    int full_mask = (1 << n) - 1;
    int last = -1;
   for (int i = 1; i < n; ++i) {
        double cur = dp[full_mask][i] + dist(i, 0);
        if (cur < cost) {
            cost = cur;
            last = i;
        }
    }

    cout << "[Held-Karp] cost: " << cost << '\n';
    int mask = full_mask;
    int path_idx = 0;
    while(last != -1){
        path[path_idx++] = last;
        int prev = parent[mask][last];
        mask ^= (1 << last);
        last = prev;
    }
    cout << "[Held-Karp] route: ";
    for(int i = path_idx - 1; i >= 0; i--){
        cout << path[i] << ' ';
    }
    cout << "0\n";


    for (size_t i = 0; i < subsets; ++i) {
        delete[] dp[i];
        delete[] parent[i];
    }
    delete[] dp;
    delete[] parent;
}

void run_2app(){
    const int n = city_count;
    const size_t MAX_MEM_BYTES = 1ULL << 32; 
    const size_t graph_bytes = n * n * sizeof(int);
    
    if (graph_bytes > MAX_MEM_BYTES) {
        cerr << "[2-app] distance matrix too large for memory limit\n";
        return;
    }

    int** dist_matrix = new (nothrow) int*[n];
    for (int i = 0; i < n; ++i) {
        dist_matrix[i] = new (nothrow) int[n];
        if (!dist_matrix[i]) {
            cerr << "[2-app] failed to allocate row " << i << " for distance matrix\n";
            for (int j = 0; j < i; ++j) delete[] dist_matrix[j];
            delete[] dist_matrix;
            return;
        }
    }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            dist_matrix[i][j] = dist(i, j);

    int* mst_parent = new int[n];
    int** multigraph = new int*[n];
    for (int i = 0; i < n; ++i) {
        multigraph[i] = new int[n];
        memset(multigraph[i], 0, sizeof(int) * n);
    }
    int* euler_path = new int[2 * n];
    int euler_len = 0;
    int* visited = new int[n]();
    int* tsp_path = new int[n];
    int tsp_len = 0;

    cout << "[2-app] distance matrix allocated: " << (graph_bytes >> 20) << " MB\n";

    compute_mst(dist_matrix, n, mst_parent);
    cout << "[2-app] MST computed\n";

    for (int i = 0; i < n; ++i)
        delete[] dist_matrix[i];
    delete[] dist_matrix;

    build_2appgraph(n, mst_parent, multigraph);
    cout << "[2-app] Eulerian multigraph built\n";

    find_eulerian_tour(n, multigraph, euler_path, &euler_len);
    cout << "[2-app] Eulerian tour length: " << euler_len << '\n';

    shortcut_tour(euler_path, euler_len, visited, tsp_path, &tsp_len);
    cout << "[2-app] Hamiltonian TSP path length: " << tsp_len << '\n';

    cout << "[2-app] route: ";
    for (int i = 0; i < tsp_len; ++i) cout << tsp_path[i] << ' ';
    cout << tsp_path[0] << '\n';

    int tsp_cost = 0;
    for (int i = 0; i < tsp_len; ++i) {
        int u = tsp_path[i];
        int v = tsp_path[(i + 1) % tsp_len];
        tsp_cost += dist(u, v);
    }
    cout << "[2-app] total path cost: " << tsp_cost << '\n';

    
    for (int i = 0; i < n; ++i) {
        delete[] multigraph[i];
    }
    delete[] multigraph;
    delete[] mst_parent;
    delete[] euler_path;
    delete[] visited;
    delete[] tsp_path;
}

void run_christofides() {
    const int n = city_count;
    const size_t MAX_MEM_BYTES = 1ULL << 32; // 1 GB
    const size_t graph_bytes = n * n * sizeof(int);
    
    if (graph_bytes > MAX_MEM_BYTES) {
        cerr << "[Christofides] distance matrix too large for memory limit\n";
        return;
    }

    int** dist_matrix = new (nothrow) int*[n];
    for (int i = 0; i < n; ++i) {
        dist_matrix[i] = new (nothrow) int[n];
        for (int j = 0; j < n; ++j)
            dist_matrix[i][j] = dist(i, j);
    }

    int* mst_parent = new int[n];
    int* degree = new int[n]();
    int* odd = new int[n];
    int odd_cnt = 0;
    Edge* matching = new Edge[n]; 
    int** multigraph = new int*[n];
    for (int i = 0; i < n; ++i) {
        multigraph[i] = new int[n];
        memset(multigraph[i], 0, sizeof(int) * n);
    }
    int* euler_path = new int[2 * n];
    int euler_len = 0;
    int* visited = new int[n]();
    int* tsp_path = new int[n];
    int tsp_len = 0;

    cout << "[Christofides] distance matrix allocated: " << (graph_bytes >> 20) << " MB\n";

    compute_mst(dist_matrix, n, mst_parent);
    cout << "[Christofides] MST computed\n";

    find_odd_vertices(mst_parent, n, degree, odd, &odd_cnt);
    cout << "[Christofides] odd degree vertices found: " << odd_cnt << '\n';

    delete[] degree;

    edmond_matching(dist_matrix, odd, odd_cnt, matching);
    cout << "[Christofides] minimum matching complete\n";

    for (int i = 0; i < n; ++i) 
        delete[] dist_matrix[i];
    delete[] dist_matrix;
    delete[] odd;

    build_multigraph(n, mst_parent, matching, odd_cnt / 2, multigraph);
    cout << "[Christofides] Eulerian multigraph built\n";

    find_eulerian_tour(n, multigraph, euler_path, &euler_len);
    cout << "[Christofides] Eulerian tour length: " << euler_len << '\n';

    shortcut_tour(euler_path, euler_len, visited, tsp_path, &tsp_len);
    cout << "[Christofides] Hamiltonian TSP path length: " << tsp_len << '\n';

    // cout << "[Christofides] route: ";
    // for (int i = 0; i < tsp_len; ++i) cout << tsp_path[i] << ' ';
    // cout << tsp_path[0] << '\n';

    int tsp_cost = 0;
    for (int i = 0; i < tsp_len; ++i) {
        int u = tsp_path[i];
        int v = tsp_path[(i + 1) % tsp_len];
        tsp_cost += dist(u, v);
    }
    cout << "[Christofides] total path cost: " << tsp_cost << '\n';

    
    for (int i = 0; i < n; ++i) {
        delete[] multigraph[i];
    }
    delete[] multigraph;
    delete[] mst_parent;
    delete[] matching;
    delete[] euler_path;
    delete[] visited;
    delete[] tsp_path;
}

void compute_tsp_held_karp(int* vertices, int k, int** dist_matrix, int** multigraph) {
    const int INFE = 1e9;

    int full = (1 << k);
    int** dp = new int*[full];
    int** parent = new int*[full];
    for (int i = 0; i < full; ++i) {
        dp[i] = new int[k];
        parent[i] = new int[k];
        for (int j = 0; j < k; ++j) {
            dp[i][j] = INFE;
            parent[i][j] = -1;
        }
    }
    dp[1][0] = 0;

    for(int mask = 0; mask < (1 << k); ++mask){
        // starting node i
        for(int i = 0; i < k; i++){
            int prev_mask = mask ^ (1 << i);
            // node visited right before Last node
            if(!(mask & (1 << i)))  continue;
            for(int j = 0; j < k; j++){
                // Last visited node 
                if(i == j || !(mask & (1 << i))) continue;
                double new_cost = dp[prev_mask][j] + dist(j, i);                
                if (new_cost < dp[mask][i]) {
                    dp[mask][i] = new_cost;
                    parent[mask][i] = j;
                }
            }
        }
    }
    cout << "[my-algo] dp done\n";

    // 최소 TSP 완성 비용 및 마지막 정점 찾기
    int best_cost = INFE, last = -1;
    for (int i = 1; i < k; ++i) {
        int cost = dp[full - 1][i] + dist_matrix[vertices[i]][vertices[0]];
        if (cost < best_cost) {
            best_cost = cost;
            last = i;
        }
    }
    cout << "[my-algo] found last ele\n";

    // 경로 복원
    int* path = new int[k + 1];
    int mask = full - 1;
    int cur = last;
    for (int i = k - 1; i >= 0; --i) {
        path[i] = cur;
        int temp = parent[mask][cur];
        mask ^= (1 << cur);
        cur = temp;
    }
    path[k] = 0; // 돌아오는 edge (0으로)
    cout << "[my-algo] backtrack trace\n";

    // multigraph에 간선 추가
    for (int i = 0; i < k; ++i) {
        int u = vertices[path[i]];
        int v = vertices[path[i + 1]];
        multigraph[u][v]++;
        multigraph[v][u]++;
    }

    // 메모리 해제
    for (int i = 0; i < full; ++i) {
        delete[] dp[i];
        delete[] parent[i];
    }
    delete[] dp;
    delete[] parent;
    delete[] path;
}

void run_my_algorithm1() {
    const int n = city_count;
    const int MAX_TSP_SIZE = 20;

    // (1) dist matrix set
    int** dist_matrix = new int*[n];
    for (int i = 0; i < n; ++i) {
        dist_matrix[i] = new int[n];
        for (int j = 0; j < n; ++j)
            dist_matrix[i][j] = dist(i, j);
    }

    // (2) temporary array
    int* current = new int[n];
    for (int i = 0; i < n; ++i) current[i] = i;
    int current_cnt = n;

    // (3) degree accumulate array
    int* degree = new int[n];
    for (int i = 0; i < n; ++i) degree[i] = 0;

    // (4) multigraph set
    int** multigraph = new int*[n];
    for (int i = 0; i < n; ++i) {
        multigraph[i] = new int[n];
        for (int j = 0; j < n; ++j) multigraph[i][j] = 0;
    }

    cout << "[my-algo] MST loop\n";
    // (5) multi-level MST
    while (true) {
        if (current_cnt <= MAX_TSP_SIZE) {
            compute_tsp_held_karp(current, current_cnt, dist_matrix, multigraph);
            cout << "[my-algo] Held-Karp applied at size " << current_cnt << '\n';
            break;
        }
        
        cout << "[my-algo] small than TSP SIZE" << endl;

        // (5-1) sub distance re-calculate
        int** sub_dist = new int*[current_cnt];
        for (int i = 0; i < current_cnt; ++i) {
            sub_dist[i] = new int[current_cnt];
            for (int j = 0; j < current_cnt; ++j)
                sub_dist[i][j] = dist_matrix[current[i]][current[j]];
        }
        cout << "[my-algo] sub matrix done" << endl;

        // (5-2) MST 
        int* mst_parent = new int[current_cnt];
        compute_mst(sub_dist, current_cnt, mst_parent);

        cout << "[my-algo] MST done\n";

        // (5-3) degree accumulate
        for (int i = 1; i < current_cnt; ++i) {
            int u = current[mst_parent[i]];
            int v = current[i];
            multigraph[u][v]++;
            multigraph[v][u]++;
            degree[u]++; 
            degree[v]++;
        }

        // (5-4) extract odd degree node
        int* next = new int[current_cnt];
        int next_cnt = 0;
        for (int i = 0; i < current_cnt; ++i) {
            int v = current[i];
            if (degree[v] % 2 == 1)
                next[next_cnt++] = v;
        }
        cout << "[my-algo] after iteration, next cnt is " << next_cnt << endl;
        // free temp memory
        for (int i = 0; i < current_cnt; ++i) delete[] sub_dist[i];
        delete[] sub_dist;
        delete[] mst_parent;
        delete[] current;
        current = next;
        current_cnt = next_cnt;
    }

    // (6) Euler tour 
    int* euler_path = new int[n * n + 500];
    int euler_len = 0;
    find_eulerian_tour(n, multigraph, euler_path, &euler_len);
    cout << "[my-algo] Eulerian tour length: " << euler_len << '\n';

    // (7) Euler tour → Hamiltonian TSP approximation by removing duplicate node in route
    int* visited = new int[n];
    for (int i = 0; i < n; ++i) visited[i] = 0;

    int* tsp_path = new int[n];
    int tsp_len = 0;
    shortcut_tour(euler_path, euler_len, visited, tsp_path, &tsp_len);
    cout << "[my-algo] Hamiltonian TSP path length: " << tsp_len << '\n';

    // (8) route computation
    cout << "[my-algo] route: ";
    for (int i = 0; i < tsp_len; ++i) cout << tsp_path[i] << ' ';
    cout << tsp_path[0] << '\n';

    int tsp_cost = 0;
    for (int i = 0; i < tsp_len; ++i) {
        int u = tsp_path[i];
        int v = tsp_path[(i + 1) % tsp_len];
        tsp_cost += dist(u, v);
    }
    cout << "[my-algo] total path cost: " << tsp_cost << '\n';

    // (9) memory free
    for (int i = 0; i < n; ++i) {
        delete[] dist_matrix[i];
        delete[] multigraph[i];
    }
    delete[] dist_matrix;
    delete[] multigraph;
    delete[] current;
    delete[] degree;
    delete[] euler_path;
    delete[] tsp_path;
    delete[] visited;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "call: " << argv[0] << " [dataset] [algorithm]\n";
        return 1;
    }

    const char* dataset = argv[1];
    const char* algorithm = argv[2];

    char filename[128];
    if (strcmp(dataset, "a280") == 0) {
        strcpy(filename, "a280.tsp");
    }
    else if (strcmp(dataset, "kz9976") == 0) {
        strcpy(filename, "kz9976.tsp");
    }
    else if (strcmp(dataset, "xql662") == 0) {
        strcpy(filename, "xql662.tsp");
    }
    else if (strcmp(dataset, "mona") == 0) {
        strcpy(filename, "mona-lisa100k.tsp");
    }
    else if (strcmp(dataset, "ulys") == 0) {
        strcpy(filename, "ulysses22.tsp");
    }
    else if (strcmp(dataset, "ital") == 0) {
        strcpy(filename, "it16862.tsp");
    }
    else {
        cerr << "unknown dataset: " << dataset << std::endl;
        return 1;
    }

    if (!load_tsp_file(filename)) {
        cerr << "load err\n";
        return 1;
    }

    cout << "loaded " << city_count << " cities from " << filename << endl;

    if (strcmp(algorithm, "held") == 0) {
        run_with_timer("Held-Karp", run_held_karp);
    }
    else if (strcmp(algorithm, "chris") == 0) {
        run_with_timer("Christofides", run_christofides);
    }
    else if(strcmp(algorithm, "2-app") == 0) {
        run_with_timer("2-app", run_2app);
    }
    else if (strcmp(algorithm, "myalgo") == 0) {
        run_with_timer("my-algo", run_my_algorithm1);
    }
    else {
        cerr << "unknown algorithm: " << algorithm << endl;
        return 1;
    }

    return 0;
}