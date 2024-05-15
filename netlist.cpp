#include "netlist.hpp"
#include "Grid_Graph.hpp"
#include "batch.hpp"
#include "struct.hpp"

#ifndef ROUTE_MEM_ESTIMATE
#define ROUTE_MEM_ESTIMATE 10
#endif

#define costfile "pattern_costs.txt"

extern int MAZE_ROUTE_ITER;
extern int BOX_MIN_DIM;
extern int NUM_THREADS;

// Constructor
Netlist::Netlist(int NUM_THREADS,Grid_Graph G, const std::vector<int>& v1, const std::vector<int>& v2, const std::vector<int>& v3, const std::vector<int>& v4) {
    
    int N = v1.size();
    nets.resize(N);
    for (auto& net : nets) {
        net.route.reserve(ROUTE_MEM_ESTIMATE); // Initialize route vector with capacity ROUTE_MEM_ESTIMATE, so as to minimize later time wasted in allocations
    }
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < N; ++i) {
        Point point = {v1[i], v4[i]};
        for (int y = v2[i]; y <= point.y; y++) {
            G.Gy[v1[i] * (G.M + 1) + y] += 1;
        }
        for (int y = point.y; y <= v2[i]; y++) {
            G.Gy[v1[i] * (G.M + 1) + y] += 1;
        }
        for (int x = point.x; x < v3[i]; x++) {
            G.Gx[v4[i] * (G.N + 1) + x] += 1;
        }
        for (int x = v3[i]; x < point.x; x++) {
            G.Gx[v4[i] * (G.N + 1) + x] += 1;
        }
        nets[i].x1 = v1[i];
        nets[i].y1 = v2[i];
        nets[i].x2 = v3[i];
        nets[i].y2 = v4[i];
        if (!(((point.x==nets[i].x1)&&(point.y==nets[i].y1))||((point.x==nets[i].x2)&&(point.y==nets[i].y2)))){
            nets[i].route.push_back(point);
        }
    }
}

// Function to schedule patterns
void Netlist::pattern_schedule() {
    int N = nets.size();
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    for (int i = N - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        std::swap(nets[i], nets[j]);
    }

    std::vector<int> clique(N, 0);
    std::vector<int> batch_count;
    int id = 1;

    for (int i = 0; i < N; ++i) {
        if (clique[i] == 0) {
            clique[i] = id;
            int count = 1;
            for (int j = i + 1; j < N; ++j) {
                if (clique[j]==0){
                    clique[j] = id;
                    for (int l = i ; l < j; ++l){
                        if (clique[l] == id && !checkRectangleIntersection(nets[l], nets[j])) {
                        }
                        else if (clique[l] == id && checkRectangleIntersection(nets[l], nets[j])){
                            clique[j] = 0;
                            break;
                        }
                    }
                }
                if (clique[j]){
                    count++;
                }
            }

            batch_count.push_back(count);
            id++;
        }
    }
    id -= 1;
    batches.clear();
    std::cerr << "batches are cleared" << std::endl;
    batches.resize(id);
    for(int j=1; j<id+1; j++){
        for (int l=1; l<clique.size();l++){
            if (clique[l]==j){
                batches[j].nets.emplace_back(nets[l]);
            }
        }
    }
    std::cerr << "Done pushing " << std::endl;
    for (int kk = 0; kk<batches.size();kk++){
        for (int m = 0; m<batches[kk].N;m++){
            for (int o = 0; o<m;o++){
                if (checkRectangleIntersection(batches[kk].nets[m], batches[kk].nets[o])){
                    std::cerr << " overlapping nets " << o << " " << m << " " << kk << std::endl;
                    exit(1);
                }
            }
        }
    }
    std::cerr << "Done checking " << std::endl;
    return;
}


void Netlist::maze_schedule(Grid_Graph G,float k, int BOX_MIN_DIM) {
    int N = nets.size();
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    for (int i = N - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        std::swap(nets[i], nets[j]);
    }

    std::vector<int> clique(N, 0);
    std::vector<int> batch_count;
    int id = 1;

    for (int i = 0; i < N; ++i) {
        if (clique[i] == 0) {
            clique[i] = id;
            int count = 1;
            for (int j = i + 1; j < N; ++j) {
                if (clique[j]==0){
                    clique[j] = id;
                    for (int l = i ; l < j; ++l){
                        if (clique[l] == id && !overlap(nets[l], nets[j], k, BOX_MIN_DIM)) {
                        }
                        else if (clique[l] == id && overlap(nets[l], nets[j], k, BOX_MIN_DIM)){
                            clique[j] = 0;
                            break;
                        }
                    }
                }
                if (clique[j]){
                    count++;
                }
            }

            batch_count.push_back(count);
            id++;
        }
    }
    id -= 1;
    batches.clear();
    std::cerr << "batches are cleared" << std::endl;
    batches.resize(id);
    for(int j=1; j<id; j++){
        for (int l=1; l<clique.size();l++){
            if (clique[l]==j){
                batches[j].nets.emplace_back(nets[l]);
            }
        }
        std::cerr << "pshed back " << j << std::endl;
    }
    std::cerr << "Done pushing " << std::endl;
    for (int kk = 0; kk<batches.size();kk++){
        for (int m = 0; m<batches[kk].N;m++){
            for (int o = 0; o<m;o++){
                if (overlap(batches[kk].nets[m], batches[kk].nets[o], k, BOX_MIN_DIM)){
                    std::cerr << " overlapping nets " << o << " " << m << " " << kk << std::endl;
                    exit(1);
                }
            }
        }
    }
    std::cerr << "Done checking " << std::endl;
    return;
}

// Function for simulated annealing pattern routing
float Netlist::SA_patternroute(Grid_Graph G) {

    float T = 1000;
    int N = nets.size();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<float> costs;

    float tot_cost = std::numeric_limits<float>::max();
    float new_tot_cost = 0;
  
    for (auto& batch : batches) {
        new_tot_cost += batch.pattern_route(G, 0, T, 1);
    }
    tot_cost = new_tot_cost;

    // Start the simulated annealing process
    while (T > 0.01) {
        new_tot_cost = 0;
        for (auto& batch : batches) {
            new_tot_cost += batch.pattern_route(G, 0, T, dist(gen));
        }
        tot_cost = new_tot_cost;
        T *= 0.995;
        costs.push_back(tot_cost);
    }

    // Save patterns for each batch
    int k = 0;
    for (auto& batch : batches) {
        k += batch.N;
    }

    // Write costs to output file
    std::ofstream outputFile(costfile);
    if (!outputFile.is_open()) {
        return -1.0f;
    }

    for (size_t i = 0; i < costs.size(); ++i) {
        outputFile << costs[i];
        if (i != costs.size() - 1) {
            outputFile << " ";
        }
    }

    outputFile.close();
    std::cerr << "Finished writing costs." << std::endl;

    std::cerr << "Total final cost: " << tot_cost << std::endl;
    std::cerr << "Finished SA_patternroute." << std::endl;
    return tot_cost;
}

// Function to perform maze routing
void Netlist::mazer(Grid_Graph G,float k, int NUM_THREADS,int BOX_MIN_DIM, int MAZE_ROUTE_ITER) {
    std::vector<float> Sdist1(G.M * G.N);
    std::vector<char> Sdir1(G.M * G.N);
    std::vector<float> Sdist2(G.M * G.N);
    std::vector<char> Sdir2(G.M * G.N);
    #pragma omp parallel num_threads(NUM_THREADS) 
    {
        for (int i=0; i<MAZE_ROUTE_ITER; i++){
            std::cerr << batches.size();
            for (int j=0; j<batches.size(); j++) {
                std::cerr << "Maze routing a batch " << std::endl;
                #pragma omp barrier
                batches[j].maze_route(G, k, 2, Sdist1, Sdir1, Sdist2, Sdir2, NUM_THREADS, BOX_MIN_DIM);
            }
            std::cerr << "Iteration " << i << std::endl;
        }
    }
}

// Function to check if two nets overlap
inline bool Netlist::overlap(const Net& net1, const Net& net2, float k, int BOX_MIN_DIM) {
    int centerX1 = (net1.x1 + net1.x2) / 2;
    int centerY1 = (net1.y1 + net1.y2) / 2;
    int centerX2 = (net2.x1 + net2.x2) / 2;
    int centerY2 = (net2.y1 + net2.y2) / 2;

    int width1 = std::max((int)std::ceil(k * std::abs(net1.x1 - net1.x2)),BOX_MIN_DIM);
    int height1 = std::max((int)std::ceil(k * std::abs(net1.y1 - net1.y2)),BOX_MIN_DIM);
    int width2 = std::max((int)std::ceil(k * std::abs(net2.x1 - net2.x2)),BOX_MIN_DIM);
    int height2 = std::max((int)std::ceil(k * std::abs(net2.y1 - net2.y2)),BOX_MIN_DIM);

    int minX1 = centerX1 - width1 / 2;
    int maxX1 = centerX1 + width1 / 2;
    int minY1 = centerY1 - height1 / 2;
    int maxY1 = centerY1 + height1 / 2;

    int minX2 = centerX2 - width2 / 2;
    int maxX2 = centerX2 + width2 / 2;
    int minY2 = centerY2 - height2 / 2;
    int maxY2 = centerY2 + height2 / 2;

    return !(minX1 > maxX2 || maxX1 < minX2 || minY1 > maxY2 || maxY1 < minY2);
}

// Function to check if two nets' rectangles intersect
inline bool Netlist::checkRectangleIntersection(const Net& net1, const Net& net2) {
    bool xOverlap = net1.x1 <= net2.x2 && net1.x2 >= net2.x1;
    bool yOverlap = net1.y1 <= net2.y2 && net1.y2 >= net2.y1;

    bool fullyInside = (net1.x1 >= net2.x1 && net1.x2 <= net2.x2 && net1.y1 >= net2.y1 && net1.y2 <= net2.y2) ||
                       (net2.x1 >= net1.x1 && net2.x2 <= net1.x2 && net2.y1 >= net1.y1 && net2.y2 <= net1.y2);

    return xOverlap && yOverlap && !fullyInside;
}
