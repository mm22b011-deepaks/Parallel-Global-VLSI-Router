#ifndef BATCH_HPP
#define BATCH_HPP

#include "main.hpp"

struct Net;
struct Point;
class Grid_Graph;

// Class representing a batch
class Batch {
public:
    std::vector<Net> nets;  // Vector of nets
    int N; // size of batch

    // Constructor
    Batch(const std::vector<Net> netVector, int i, int n);
    Batch(); 

    // Function to perform pattern route
    float pattern_route(Grid_Graph G, int k, float T, float p);

    // Function to perform maze route
    void maze_route(Grid_Graph G, float k, float c, std::vector<float>& Sdist1, std::vector<char>& Sdir1, std::vector<float>& Sdist2, std::vector<char>& Sdir2, int NUM_THREADS, int BOX_MIN_DIM);

    void save_patterns(std::vector<Point>& L, int k);
    void rip_wire(Grid_Graph G, Point Src, Point Dest, std::vector<Point> path);
    void route_wire(Grid_Graph G, Point Src, Point Dest, std::vector<Point> path);

    inline float weight(float demand, float capacity);
    inline float weight_pr(float demand, float capacity);
    float survey(Grid_Graph G, Net net, int orientation);
    void routeL(Grid_Graph G, Net net, int orientation);
    float ripL(Grid_Graph G, Net net, int orientation);
    void rip_line(Grid_Graph G, Point Src, Point Dest);
    void route_line(Grid_Graph G, Point Src, Point Dest);
    void disp_s(std::vector<float> Sdist1,std::vector<float> Sdist2,std::vector<char> Sdir1,std::vector<char> Sdir2,int M,int N);
    void dispG(Grid_Graph G);
};

#endif // BATCH_HPP

