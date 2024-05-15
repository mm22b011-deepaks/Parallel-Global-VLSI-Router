#include "Grid_Graph.hpp"
#include "batch.hpp"
#include "netlist.hpp"
#include "main.hpp"
#include "struct.hpp"
#include <omp.h>

void readDataFromFile(const char* filename, int* M, int* N, int* C, int* v, std::vector<int>& x1, std::vector<int>& y1, std::vector<int>& x2, std::vector<int>& y2) ;
void StoreToFile(const char* filename,Grid_Graph G, Netlist Netlist); 
void disp_routelist(Netlist netlist) {
    for (const auto& net : netlist.nets) {
        for (const auto& point : net.route) {
            std::cout << "(" << point.x << "," << point.y << ") ";
        }
        std::cout << std::endl;
    }
}

void disp_netlist(Netlist netlist) {
    for (const auto& net : netlist.nets) {
        std::cout << net.x1 << " " << net.y1 << " " << net.x2 << " " << net.y2 << std::endl;
    }
}

void disp_batches(Netlist netlist) {
    for (size_t i = 0; i < netlist.batches.size(); ++i) {
        std::cout << "Batch number " << i << std::endl;
        const Batch& batch = netlist.batches[i];
        for (const auto& net : batch.nets) {
            std::cout << net.x1 << " " << net.y1 << " " << net.x2 << " " << net.y2 << std::endl;
        }
    }
}

void disp_routes(Netlist netlist) {
    for (size_t i = 0; i < netlist.batches.size(); ++i) {
        std::cout << "Batch number " << i << std::endl;
        const Batch& batch = netlist.batches[i];
        for (const auto& net : batch.nets) {
            std::cout << "Net " << net.x1 << "," << net.y1 << " to " << net.x2 << "," << net.y2 << " Route:" << std::endl;
            for (const auto& point : net.route) {
                std::cout << point.x << "," << point.y << std::endl;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    /*
--- variables initialized- 
    M rows
    N cols
    C capacity
    V cost of bending
    x1,x2,y1,y2
    Gx,Gy grid graph edge weights
    net_x stores x coordinates of corners in each path
    net_y stores y coordinates of edges in each path
--- functions called-
    ReadDataFromFile
    Router
    */
    if (argc != 6) {
        std::cerr << "Correct format is " << argv[0] << " <input filename> <output filename> <Bounding Box Dimensions> <Number of Iterations for Maze Route> <Number of Threads for Parallelising>\n";
        return 1;
    }
     
    char* infilename = argv[1];
    char* outfilename = argv[2];
    int BOX_MIN_DIM = atoi(argv[3]);
    int MAZE_ROUTE_ITER = atoi(argv[4]);
    int NUM_THREADS = atoi(argv[5]);
    
    // read data from file
    int M, N, C, v;
    std::vector<int> x1, y1, x2, y2;
    readDataFromFile(infilename, &M, &N, &C, &v, x1, y1, x2, y2);
    // create the return data structures, Gx, Gy are dynamically allocated with pointers due to large size, and as it will be passed to GPU
    Grid_Graph G(M,N,C,v);
    //std::cout << "Grid graph successfully made" << std::endl;
    Netlist Netlist(NUM_THREADS,G,x1,y1,x2,y2);
    std::cerr << "Netlist successfully made" << std::endl;
    disp_routelist(Netlist);
    disp_netlist(Netlist);

    double T = omp_get_wtime();
    std::cout << "Pattern routing successfully scheduled" << std::endl;

    std::cerr << "Pattern routing successfully performed" << std::endl;
    disp_routes(Netlist);
    std::cerr << "Batches for maze routing " << std::endl;
    disp_batches(Netlist);
    float k = 4;
    Netlist.maze_schedule(G,k, BOX_MIN_DIM);
    std::cerr << "Maze is scheduled " << std::endl;
    double T__ = omp_get_wtime();
    Netlist.mazer(G,k, NUM_THREADS, BOX_MIN_DIM, MAZE_ROUTE_ITER);
    double T___ = omp_get_wtime();
    std::cerr << "Maze routing takes  " << T___ -T__  << " seconds" << std::endl; 
    // Store the result
    StoreToFile(outfilename, G, Netlist);
    std:: cerr << "succesfully stored to file " << std::endl;
    return 0;
}

void readDataFromFile(const char* filename, int* M, int* N, int* C, int* v, std::vector<int>& x1, std::vector<int>& y1, std::vector<int>& x2, std::vector<int>& y2) {
    
    std::ifstream file(filename);
    std::string line;
    int lineNumber = 0;

    if (file.is_open()) {
        while (getline(file, line)) {
            if (line.empty() || line[0] == ';') {
                continue;
            }
            
            std::istringstream iss(line);
            int a, b, c, d;
            iss >> a >> b >> c >> d;

            if (lineNumber == 0) {
                // First line contains M, N, C, v
                *M = a;
                *N = b;
                *C = c;
                *v = d;
            } else {
                // Remaining lines contain x1, y1, x2, y2
                x1.push_back(a);
                y1.push_back(b);
                x2.push_back(c);
                y2.push_back(d);
            }

            lineNumber++;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}


void StoreToFile(const char* filename,Grid_Graph G, Netlist Netlist) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Write M and N to the file
        file << G.M << " " << G.N << std::endl;

        // Write Gx to the file
        
        for (int i = 0; i < G.M; ++i) {
            file << 0 << " ";
            for (int j = 1; j < G.N ; ++j) {
                file << G.Gx[i * (G.N + 1) + j];
                if (j < G.N) {
                    file << " ";
                }
            }
            file << 0 << " ";
            file << std::endl;
        }

        // Write Gy to the file
        for (int i = 0; i < G.M ; ++i) {
            for (int j = 0; j < G.N; ++j) {
                file << G.Gy[i * (G.N+1) + j];
                if (j < G.N - 1) {
                    file << " ";
                }
            }
            file << std::endl;
        }
        for (int j = 0; j < G.N; ++j) {
            file << 0;
            if (j < G.N - 1) {
                file << " ";
            }
        }
        file << std::endl;
        // Write net_x and net_y pairs to the file
        for (int k=0; k<Netlist.batches.size();k++){
            std::vector<Net> nets = Netlist.batches[k].nets;
            for (size_t i = 0; i < nets.size(); ++i) {
                file << nets[i].x1 << " " << nets[i].y1 << " ";
                for (size_t j = 0; j < nets[i].route.size(); ++j) {
                    file << nets[i].route[j].x << " " << nets[i].route[j].y << " ";
                }
                file << nets[i].x2 << " " << nets[i].y2 << " ";
                file << std::endl;
            }
        }
        file.close();
    } 
    
    else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
