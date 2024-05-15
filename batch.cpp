#include "batch.hpp"
#include "Grid_Graph.hpp"
#include "struct.hpp"
#include "main.hpp"

extern int NUM_THREADS;
extern int BOX_MIN_DIM;

Batch::Batch(std::vector<Net> netVector, int i, int n) : N(n) {
    // Store elements from netVector starting from index i up to index i+n into a vector named nets
    for (int j = i; j < i + n && j < netVector.size(); ++j) {

        nets.push_back(netVector[j]);
    }
}
Batch::Batch() : N(0) {} // Initializes N to 0

void Batch::save_patterns(std::vector<Point>& L, int k){
    for (int i=0; i<N; i++){
        Point temp = {L[i+k].x,L[i+k].y};
        nets[i].route[0] = temp;
    }
    return;
}

// Function to perform pattern route
float Batch::pattern_route(Grid_Graph G, int k, float T, float p) {
    float tot_cost = 0;
    float cost;
    float cost_new;

    for (int i = 0; i < N; i++) {
        if (nets[i].route.size()){
            Point point = nets[i].route[0];
            int orientation;
            if (point.x == nets[i].x1 && point.y == nets[i].y2) {
                orientation = 0;
            } else {
                orientation = 1;
            }
            cost = ripL(G, nets[i], orientation);
            cost_new = survey(G, nets[i], 1 - orientation);
            if (cost_new < cost) {
                routeL(G, nets[i], 1 - orientation);
                point = {nets[i].x2, nets[i].y1};
                nets[i].route[0] = point;
            } else {
                if (exp((cost - cost_new) / T) > p) {
                    routeL(G, nets[i], 1 - orientation);
                    point = {nets[i].x2, nets[i].y1};
                    nets[i].route[0] = point;
                } else {
                    routeL(G, nets[i], orientation);
                }
            }

            tot_cost += cost_new;
        }
    }

    return tot_cost;
}


float Batch::ripL(Grid_Graph G, Net net, int orientation) {
    float costnew = 0;
    if (orientation == 0) {
        // First traverse vertically, then horizontally (turn at x1, y2)
        // Vertical traversal
        if (net.y1 < net.y2) {
            for (int y = net.y1 + 1; y < net.y2; ++y) {
                int gy_index = net.x1 * (G.M + 1) + y+1;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
            }
        } else {
            for (int y = net.y1 - 1; y > net.y2; --y) {
                int gy_index = net.x1 * (G.M + 1) + y;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
            }
        }

        // Horizontal traversal
        if (net.x1 < net.x2) {
            for (int x = net.x1 + 1; x <= net.x2; ++x) {
                int gx_index = net.y2 * (G.N + 1) + x+1;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
            }
        } 
        else {
            for (int x = net.x1 - 1; x >= net.x2; --x) {
                int gx_index = net.y2 * (G.N + 1) + x;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
            }
        }
    } else {
        // First traverse horizontally, then vertically (turn at x2, y1)
        // Horizontal traversal
        if (net.x1 < net.x2) {
            for (int x = net.x1 + 1; x < net.x2; ++x) {
                int gx_index = net.y1 * (G.N + 1) + x+1;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
            }
        }
        else {
            for (int x = net.x1 - 1; x > net.x2; --x) {
                int gx_index = net.y1 * (G.N + 1) + x;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
            }
        }

        // Vertical traversal
        if (net.y1 < net.y2) {
            for (int y = net.y1 + 1; y <= net.y2; ++y) {
                int gy_index = net.x2 * (G.M + 1) + y+1;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
            }
        }
        else {
            for (int y = net.y1 - 1; y >= net.y2; --y) {
                int gy_index = net.x2 * (G.M + 1) + y;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
                //std::cout << costnew << std::endl;
            }
        }
    }

    return costnew; 
}

void Batch::routeL(Grid_Graph G,Net net,int orientation){
    float cost = 0;
    if (orientation == 0) {
        // Traverse the path with orientation 0 (turn at x1, y2)
        for (int y = net.y1 + 1; y < net.y2; ++y) {
            G.Gy[net.x1*(G.M+1)+y] += 1;
        }
        for (int x = net.x1 + 1; x <= net.x2; ++x) {
            G.Gx[net.y2*(G.N+1)+x] += 1;
        }
    } else {
        // Traverse the path with orientation 1 (turn at x2, y1)
        for (int x = net.x1 + 1; x < net.x2; ++x) {
            G.Gx[net.y1*(G.N+1)+x] += 1;
        }
        for (int y = net.y1; y <= net.y2; ++y) {
            G.Gy[net.x2*(G.M+1)+y] += 1;
        }
    }
}

float Batch::survey(Grid_Graph G,Net net,int orientation){
    float cost = 0;
    if (orientation == 0) {
        // Traverse the path with orientation 0 (turn at x1, y2)
        for (int y = net.y1 + 1; y < net.y2; ++y) {
            cost += 1;
        }
        for (int x = net.x1 + 1; x <= net.x2; ++x) {
            cost += 1;
        }
    } else {
        // Traverse the path with orientation 1 (turn at x2, y1)
        for (int x = net.x1 + 1; x < net.x2; ++x) {
            cost += 1;
        }
        for (int y = net.y1; y <= net.y2; ++y) {
            cost += 1;
        }
    }
    return cost;
}
// Function to perform maze route
void Batch::maze_route(Grid_Graph G, float k, float c,std::vector<float>& Sdist1,std::vector<char>&  Sdir1,std::vector<float>& Sdist2,std::vector<char>&  Sdir2, int NUM_THREADS, int BOX_MIN_DIM) {
    // rip up the full batch
    if (nets.size()<999999999999)
    {
        Point source;
        Point dest;
        int centerX1;
        int centerY1;
        int width1;
        int height1;
        Point cornerl;
        Point cornerh;

        double t_start, t_end;
        int Thread_num = omp_get_thread_num();
        int size = omp_get_num_threads();
	int N = nets.size();
        t_start = omp_get_wtime();
        #pragma omp for
        for (int i=0; i<N; i++){
            source = {nets[i].x1,nets[i].y1};
            dest = {nets[i].x2,nets[i].y2};
            rip_wire(G, source,dest,nets[i].route);
            }
    

        #pragma omp for 
        for (int i=0; i<N; i++){
            source = {nets[i].x1,nets[i].y1};
            dest = {nets[i].x2,nets[i].y2};
            // calculate corners of bounding box
            centerX1 = (source.x + dest.x) / 2;
            centerY1 = (source.y + dest.y) / 2;
            width1 = std::max((int)std::ceil(k * std::abs(source.x- dest.x)),BOX_MIN_DIM);
            height1 = std::max((int)std::ceil(k * std::abs(source.y - dest.y)),BOX_MIN_DIM);
            cornerl = {std::max(0,(int)(centerX1 - width1 / 2)),std::max(0,(int)(centerY1 - height1 / 2))};
            cornerh = {std::min(G.N-1,(int)(centerX1 + width1 / 2)),std::min(G.M-1,(int)(centerY1 + height1 / 2))};
            for (int i = cornerl.y; i<cornerh.y+1;i++){
                for (int j = cornerl.x; j<cornerh.x+1;j++){
                    Sdist1[i*G.N+j] = std::numeric_limits<int>::max();
                    Sdir1[i*G.N+j] = 'x';
                }
            }
            for (int i = cornerl.x; i<cornerh.x+1;i++){
                for (int j = cornerl.y; j<cornerh.y+1;j++){
                    Sdist2[i*G.M+j] = std::numeric_limits<int>::max();
                    Sdir2[i*G.M+j] = 'x';
                }
            }
            Sdist1[source.y*G.N+source.x] = 0;
            // Then do Bellman-Ford
            bool flag = 1; // to keep track of if the relaxation step has caused any change to the distances and routes or not
            while (flag){
                flag = 0;
                // Relaxing Sdir1 (rows)
                for (int k = cornerl.y; k<cornerh.y+1; k++){
                    // left to right
                    for (int j = cornerl.x+1; j<cornerh.x+1; j++){
                        if (Sdist1[G.N*k+j] > Sdist1[G.N*k+j-1]+weight_pr(G.Gx[(G.N+1)*k+j],G.C)){
                            Sdist1[G.N*k+j] = Sdist1[G.N*k+j-1]+weight_pr(G.Gx[(G.N+1)*k+j],G.C);
                            Sdir1[G.N*k+j] = 'l';
                            flag = 1;
                        }
                    }
                    // right to left
                    for (int j = cornerh.x-1; j>cornerl.x-1; j--){
                        if (Sdist1[G.N*k+j] > Sdist1[G.N*k+j+1]+weight_pr(G.Gx[(G.N+1)*k+j+1],G.C)){
                            Sdist1[G.N*k+j] = Sdist1[G.N*k+j+1]+weight_pr(G.Gx[(G.N+1)*k+j+1],G.C);
                            Sdir1[G.N*k+j] = 'r';
                            flag = 1;
                        }
                    }
                }
                // via sweep
                for (int k = cornerl.y; k<cornerh.y+1; k++){
                    for (int j = cornerl.x; j<cornerh.x+1; j++){
                        if (Sdist1[G.N*k+j] > Sdist2[G.M*j+k]+G.v){
                            Sdist1[G.N*k+j] = Sdist2[G.M*j+k]+G.v;
                            Sdir1[G.N*k+j] =  'u';
                            flag = 1;
                        }
                        if (Sdist1[G.N*k+j] + G.v < Sdist2[G.M*j+k]){
                            Sdist2[G.M*j+k] = Sdist1[G.N*k+j]+G.v;
                            Sdir2[G.M*j+k] =  'd';
                            flag = 1;
                        }
                    }
                }
                // Relaxing Sdir2 (cols)
                for (int j = cornerl.x; j<cornerh.x+1; j++){
                    // south to  north
                    for (int k = cornerl.y+1; k<cornerh.y+1; k++){
                        if (Sdist2[G.M*j+k] > Sdist2[G.M*j+k-1]+weight_pr(G.Gy[(G.M+1)*j+k],G.C)){
                            Sdist2[G.M*j+k] = Sdist2[G.M*j+k-1]+weight_pr(G.Gy[(G.M+1)*j+k],G.C);
                            Sdir2[G.M*j+k] = 's';
                            flag = 1;
                        }
                    }
                    // right to left
                    for (int k = cornerh.y-1; k>cornerl.y-1; k--){
                        if (Sdist2[G.M*j+k] > Sdist2[G.M*j+k+1]+weight_pr(G.Gy[(G.M+1)*j+k+1],G.C)){
                            Sdist2[G.M*j+k] = Sdist2[G.M*j+k+1]+weight_pr(G.Gy[(G.M+1)*j+k+1],G.C);
                            Sdir2[G.M*j+k] = 'n';
                            flag = 1;
                        }
                    }
                }
                // via sweep
                for (int k = cornerl.y; k<cornerh.y+1; k++){
                    for (int j = cornerl.x; j<cornerh.x+1; j++){
                        if (Sdist1[G.N*k+j] > Sdist2[G.M*j+k]+G.v){
                            Sdist1[G.N*k+j] = Sdist2[G.M*j+k]+G.v;
                            Sdir1[G.N*k+j] =  'u';
                            flag = 1;
                        }
                        if (Sdist1[G.N*k+j] + G.v < Sdist2[G.M*j+k]){
                            Sdist2[G.M*j+k] = Sdist1[G.N*k+j]+G.v;
                            Sdir2[G.M*j+k] =  'd';
                            flag = 1;
                        }
                    }
                }
                
            }
        }


        #pragma omp single //This part of the code cannot be parallelised
        {
        for(int i=0; i<N;i++){
                Point source = {nets[i].x1,nets[i].y1};
                Point dest = {nets[i].x2,nets[i].y2};
                Point here = dest;
                int layer = 0;
                char dir = 'x';
                char diro = Sdir1[G.N*here.y+here.x]; // so that dest is not saved in route
                nets[i].route.clear();
                int k = 0;
                while (!((here.x == source.x) && (here.y == source.y))){
                    if (layer){
                        dir = Sdir2[G.M*here.x+here.y];
                    }
                    else{
                        dir = Sdir1[G.N*here.y+here.x];
                    }
                    if ((dir != diro)&&(dir!='d')&&(dir!='u')&&!((here.x == dest.x) && (here.y == dest.y))){
                        nets[i].route.push_back(here);
                    }
                    diro = dir;
                    if (dir=='x'){
                        k++;
                    }
                    else{k=0;}
                    switch(dir){
                        case 'd':
                            layer = 0; break;
                        case 'u':
                            layer = 1; break;
                        case 'l':
                            here.x -= 1; break;
                        case 'r':
                            here.x += 1; break;
                        case 'n':
                            here.y += 1; break;
                        case 's':
                            here.y -= 1; break;
                        default: break;
                    }
                    
                }
            }
        
        }

        #pragma omp for 
        for(int i=0;i<N; i++){
            std::reverse(nets[i].route.begin(),nets[i].route.end());

        // update grid graph
            Point source = {nets[i].x1,nets[i].y1};
            Point dest = {nets[i].x2,nets[i].y2};
            route_wire(G, source,dest,nets[i].route);
        }
    }
}


void Batch::rip_line(Grid_Graph G, Point Src, Point Dest)
{
    if (!((Src.x == Dest.x)&&(Src.y == Dest.y))){}
    else if (Src.x == Dest.x) {
        // only one of the for loops below will execute
        for (int y = Src.y + 1; y < Dest.y; ++y) {
            G.Gy[Src.x*(G.M+1)+y] -= 1;
        }
        for (int y = Dest.y+1; y < Src.y ; ++y) {
            G.Gy[Src.x*(G.M+1)+y] -= 1;
        }
    } 
    else {
        // only one of the for loops below will execute
        for (int x = Src.x + 1; x < Dest.x; ++x) {
            G.Gx[Src.y*(G.N+1)+x] -= 1;
        }
        for (int x = Dest.x+1; x < Src.x; ++x) {
            G.Gx[Src.y*(G.N+1)+x] -= 1;
        }
    }
}

void Batch::route_line(Grid_Graph G, Point Src, Point Dest)
{
    if (!((Src.x == Dest.x)&&(Src.y == Dest.y))){}
    else if (Src.x == Dest.x) {
        // only one of the for loops below will execute
        for (int y = Src.y + 1; y < Dest.y; ++y) {
            G.Gy[Src.x*(G.M+1)+y] += 1;
        }
        for (int y = Dest.y+1; y < Src.y ; ++y) {
            G.Gy[Src.x*(G.M+1)+y] += 1;
        }
    } 
    else {
        // only one of the for loops below will execute
        for (int x = Src.x + 1; x < Dest.x; ++x) {
            G.Gx[Src.y*(G.N+1)+x] += 1;
        }
        for (int x = Dest.x+1; x < Src.x; ++x) {
            G.Gx[Src.y*(G.N+1)+x] += 1;
        }
    }
}

void Batch::route_wire(Grid_Graph G, Point Src, Point Dest, std::vector<Point> path)
{
    int i = 0;
    if (!path.size()){
        route_line(G,Src,Dest);
    }
    else{
        route_line(G,Src,path[0]);
        for (i=1; i<path.size()-1;i++)
        {
            route_line(G,path[i-1],path[i]);
        }
        i = path.size()-1;
        route_line(G,path[i],Dest);
    }
    return;
}

void Batch::rip_wire(Grid_Graph G, Point Src, Point Dest, std::vector<Point> path)
{
    int i = 0;
    if (!path.size()){
        rip_line(G,Src,Dest);
    }
    else{
        rip_line(G,Src,path[0]);
        for (i=1; i<path.size()-1;i++)
        {
            rip_line(G,path[i-1],path[i]);
        }
        i = path.size()-1;
        rip_line(G,path[i],Dest);
    }
    return;
}

inline float Batch::weight(float demand, float capacity){
    // the return values of this function can be memoised for further speedup
    return (demand*2)/capacity;
}

inline float Batch::weight_pr(float demand, float capacity){
    // the return values of this function can be memoised for further speedup
    return 3*(std::exp(demand/capacity)-1);
}

void Batch::disp_s(std::vector<float> Sdist1,std::vector<float> Sdist2,std::vector<char> Sdir1,std::vector<char> Sdir2,int M,int N){
    std::cout << "Sdir1" << std::endl ;
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            if (Sdir1[i*M+j]){
                std::cout << Sdir1[i*M+j] << " " ;
            }
            else{
                std::cout << 'y' << " " ;
            }
        }
        std::cout << std::endl; 
    }
    std::cout << std::endl;
    std::cout << "Sdir2" << std::endl ;
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            if (Sdir2[j*N+i]){
                std::cout << Sdir2[j*N+i] << " " ;
            }
            else{
                std::cout << 'y' << " " ;
            }
        }
        std::cout << std::endl; 
    }
    std::cout << std::endl; 
    std::cout << "Sdist1" << std::endl ;
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            if (Sdist1[i*M+j]){
                std::cout << Sdist1[i*M+j] << " " ;
            }
            else{
                std::cout << 'y' << " " ;
            }
        }
        std::cout << std::endl; 
    }
    std::cout << std::endl;
    std::cout << "Sdist2" << std::endl ;
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            if (Sdist2[j*N+i]){
                std::cout << Sdist2[j*N+i] << " " ;
            }
            else{
                std::cout << 'y' << " " ;
            }
        }
        std::cout << std::endl; 
    }
}

