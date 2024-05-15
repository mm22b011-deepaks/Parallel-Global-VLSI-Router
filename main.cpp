#include "Grid_Graph.hpp"
#include "batch.hpp"
#include "netlist.hpp"
#include "main.hpp"
#include "struct.hpp"
#include <omp.h>

void readDataFromFile(const char* filename, int* M, int* N, int* C, int* v, std::vector<int>& x1, std::vector<int>& y1, std::vector<int>& x2, std::vector<int>& y2) ;
void StoreToFile(const char* filename,Grid_Graph G, Netlist Netlist); 
/*
write a void function named disp_routelist(Netlist Netlist); whose argument is of type class Netlist {
public:
    std::vector<Net> nets;      // Vector of nets} and this is net struct Net {
    int x1, y1, x2, y2;
    std::vector<Point> route;

    // Default constructor
    Net() : x1(0), y1(0), x2(0), y2(0) {}

    // Parameterized constructor
    Net(int v1, int v2, int v3, int v4, size_t size)
        : x1(v1), y1(v2), x2(v3), y2(v4), route(size) {}
};
it has to pint each route in a seperate line over all the nets in the netlist
*/


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

/*
Gen AI prompt-
Write a function void disp_batches(Netlist Netlist); that does the following: Netlist id the following class  class Netlist {
public:
    std::vector<Net> nets;      // Vector of nets
    std::vector<Batch> batches; // Vector of batches} and nets is the following struct struct Net {
    int x1, y1, x2, y2;
    std::vector<Point> route;

    // Default constructor
    Net() : x1(0), y1(0), x2(0), y2(0) {}

    // Parameterized constructor
    Net(int v1, int v2, int v3, int v4, size_t size)
        : x1(v1), y1(v2), x2(v3), y2(v4) {
        route.reserve(size);
    }
}; and batches is the following class  class Batch {
public:
    std::vector<Net> nets;  // Vector of nets
    int N; // size of batch
}  now print all the batches in the followinng format batch number (index of batch in batches vector) - 
and in the following lines print whitespace seperated x1,y1,x2,y2 for each of the nets in that batch, 1 net 
per line
*/
void disp_batches(Netlist netlist) {
    for (size_t i = 0; i < netlist.batches.size(); ++i) {
        std::cout << "Batch number " << i << std::endl;
        const Batch& batch = netlist.batches[i];
        for (const auto& net : batch.nets) {
            std::cout << net.x1 << " " << net.y1 << " " << net.x2 << " " << net.y2 << std::endl;
        }
    }
}

/*
Write a function void disp_routes(Netlist Netlist); that does the following: Netlist id the following class  class Netlist {
public:
    std::vector<Net> nets;      // Vector of nets
    std::vector<Batch> batches; // Vector of batches} and nets is the following struct struct Net {
    int x1, y1, x2, y2;
    std::vector<Point> route;

    // Default constructor
    Net() : x1(0), y1(0), x2(0), y2(0) {}

    // Parameterized constructor
    Net(int v1, int v2, int v3, int v4, size_t size)
        : x1(v1), y1(v2), x2(v3), y2(v4) {
        route.reserve(size);
    }
}; and batches is the following class  class Batch {
public:
    std::vector<Net> nets;  // Vector of nets
    int N; // size of batch
}  now print all the batches in the followinng format batch number (index of batch in batches vector) 
- and in the following lines print route for each of the nets in that batch, 1 net per line
*/
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
        std::cerr << "Correct format is ./" << argv[0] << " <input filename> <output filename> <Bounding Box Dimensions> <Number of Iterations for Maze Route> <Number of Threads for Parallelising>\n";
        return 1;
    }
     
    char* infilename = argv[1];
    char* outfilename = argv[2];
    int BOX_MIN_DIM = atoi(argv[3]);
    int MAZE_ROUTE_ITER = atoi(argv[4]);
    int NUM_THREADS = atoi(argv[5]);
    

    // Get the filename from command line argument
//    char* infilename = argv[1];
//    char* outfilename = argv[2];
    // read data from file
    int M, N, C, v;
    std::vector<int> x1, y1, x2, y2;
    readDataFromFile(infilename, &M, &N, &C, &v, x1, y1, x2, y2);
    // create the return data structures, Gx, Gy are dynamically allocated with pointers due to large size, and as it will be passed to GPU
    Grid_Graph G(M,N,C,v);
    //std::cout << "Grid graph successfully made" << std::endl;
    Netlist Netlist(NUM_THREADS,G,x1,y1,x2,y2);
    std::cerr << "Netlist successfully made" << std::endl;
    //std::cout << "displaying generated route lists and then netlists" << std::endl;
    disp_routelist(Netlist);
    disp_netlist(Netlist);
    //std::cout << std::endl;


    //Do the routing
    //float cost;
    double T = omp_get_wtime();
    //Netlist.pattern_schedule();
    std::cout << "Pattern routing successfully scheduled" << std::endl;

    //cost = Netlist.SA_patternroute(G);
    //double T_ = omp_get_wtime();
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
    //std::cerr << "The final batches are " << std::endl;
    //disp_batches(Netlist);
    //std::cout << "The final routes are " << std::endl;
    //disp_routes(Netlist);
    //G.boundary_cond();
    //std::cout << "Final cost is " << cost << std::endl; 
    //std::cout << "Pattern routing time is " << T_-T << std::endl; 
    std::cerr << "Maze routing takes  " << T___ -T__  << " seconds" << std::endl; 
    // Store the result
    StoreToFile(outfilename, G, Netlist);
    std:: cerr << "succesfully stored to file " << std::endl;
    return 0;
}


// Gen-AI prompt
/*
User
I have a text file which stores data in the following format- All lines starting with a ";" character is to be ignored, the first line which 
is ignored has 4 integers seperated by a whitespace, which are M, N, C, and v. Every line after that which is not ignored consists of 4 integers
seperated by a white space, which are the ith elements of 4 vectors x1,y1,x2,y2. Write a cpp function that creates these 4 arrays as well as 
M,N,C and v. The input to the function are pointers to M,N,C,v and vectors x1,y1,x2,y2 and the function stores the values read in them.
*/
void readDataFromFile(const char* filename, int* M, int* N, int* C, int* v, std::vector<int>& x1, std::vector<int>& y1, std::vector<int>& x2, std::vector<int>& y2) {
    
    std::ifstream file(filename);
    std::string line;
    int lineNumber = 0;

    if (file.is_open()) {
        while (getline(file, line)) {
            if (line.empty() || line[0] == ';') {
                // Ignore comment lines
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

/*
Gen-AI prompt-
Write a cpp function StoreToFile(a string filename,int M, int N,int* Gx, int* GY, vector<vector<int>> net_x, vector<vector<int>> net_y); 
 Gx is a flatenned array with M rows and N+1 columns, while Gy is a flattened array with M+1 rows and N columns. net_x and net_y have the 
 same number of sub-vectors, with corresponding sub-vectors having equall length The function should store the following data in the following format- 
 The first lines should have M N seperated by a whitespace, next Gx should be stored with elements in each row seperated by white spaces and rows in different 
 lines, and then Gy should be stored in the same format. Then for each pair of subvectors in net_x and net_y, store the subvectors as net_x[0] net_y[0] net_x[1] 
 net_y[1] ... net_x[...] net_y[...]. seperated by whitespaces, 1 pair of subvectors per row

*/

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
        /*for (int j = 0; j < G.N; ++j) {
            file << 0;
            if (j < G.N - 1) {
                file << " ";
            }
        }*/
        //file << std::endl;
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
        //std::cerr << "Netlist size is " << Netlist.batches.size() << std::endl;
        // Write net_x and net_y pairs to the file
        for (int k=0; k<Netlist.batches.size();k++){
            std::vector<Net> nets = Netlist.batches[k].nets;
            //std::cout << "k is " << k << std::endl;
            //std::cout << "Batch size is " << Netlist.batches[k].N << std::endl;
            for (size_t i = 0; i < nets.size(); ++i) {
                file << nets[i].x1 << " " << nets[i].y1 << " ";
                for (size_t j = 0; j < nets[i].route.size(); ++j) {
                    //std::cerr << "Net size is " << nets[i].route.size() << std::endl;
                    file << nets[i].route[j].x << " " << nets[i].route[j].y << " ";
                }
                file << nets[i].x2 << " " << nets[i].y2 << " ";
                file << std::endl;
            }
            //std::cerr << "k is " << k << std::endl;
        }
        //std::cerr << "Written to file " << std::endl;
        file.close();
        //std::cerr << "Closed file " << std::endl;
    } 
    
    else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}