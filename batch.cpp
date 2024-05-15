#include "batch.hpp"
#include "Grid_Graph.hpp"
#include "struct.hpp"
#include "main.hpp"
// Constructor
/*
vector<route> is bsed on old way  of storing data, has to be fixed
*/

extern int NUM_THREADS;
extern int BOX_MIN_DIM;

Batch::Batch(std::vector<Net> netVector, int i, int n) : N(n) {
    // Store elements from netVector starting from index i up to index i+n into a vector named nets
    for (int j = i; j < i + n && j < netVector.size(); ++j) {

        nets.push_back(netVector[j]);
    }
}
Batch::Batch() : N(0) {} // Initializes N to 0

/*
The below function is likely a small mess up, this has to be corrected to account for how data is stored in struct Net, it is called by Netlist.SA_patternroute();
*/
void Batch::save_patterns(std::vector<Point>& L, int k){
    for (int i=0; i<N; i++){
        Point temp = {L[i+k].x,L[i+k].y};
        nets[i].route[0] = temp;
    }
    return;
}

// Function to perform pattern route
/*
L is basically storing the set of routes for pattern routing, is this needed? There may be some issue here
*/
float Batch::pattern_route(Grid_Graph G, int k, float T, float p) {
    //std::cerr << "Performing pattern routing..." << std::endl;
    float tot_cost = 0;
    float cost;
    float cost_new;

    for (int i = 0; i < N; i++) {
    // Go over all 
        if (nets[i].route.size()){
            Point point = nets[i].route[0];
            int orientation;
            if (point.x == nets[i].x1 && point.y == nets[i].y2) {
                orientation = 0;
            } else {
                orientation = 1;
            }

            //std::cerr << "Net #" << i << ": (" << nets[i].x1 << "," << nets[i].y1 << ") to ("
                    //<< nets[i].x2 << "," << nets[i].y2 << "), Orientation: " << orientation << std::endl;

            cost = ripL(G, nets[i], orientation);
            //std::cout << "Cost of ripping is " << cost << std::endl;
            cost_new = survey(G, nets[i], 1 - orientation);
            //std::cout << "Cost of other orientation is " << cost_new << std::endl;
            if (cost_new < cost) {
                routeL(G, nets[i], 1 - orientation);
                point = {nets[i].x2, nets[i].y1};
                nets[i].route[0] = point;
                //std::cerr << "New orientation for Net #" << i << " is " << (1 - orientation) << std::endl;
            } else {
                if (exp((cost - cost_new) / T) > p) {
                    routeL(G, nets[i], 1 - orientation);
                    point = {nets[i].x2, nets[i].y1};
                    nets[i].route[0] = point;
                    //std::cerr << "Switching orientation for Net #" << i << std::endl;
                } else {
                    routeL(G, nets[i], orientation);
                    //std::cerr << "Retaining current orientation for Net #" << i << std::endl;
                }
            }

            tot_cost += cost_new;
            //std::cerr << "Total cost for Net #" << i << ": " << tot_cost << std::endl;
        }
    }

    //std::cerr << "Pattern routing completed. Total cost: " << tot_cost << std::endl;
    return tot_cost;
}


/*
Gen-AI prompt-
User
A net is a struct of 4 int's , x1,y1,x2,y2, along with a vector<Point> route, where (x1,y1) form a point on a grid, and (x2,y2) forms another, and Point is a struct Point{ int x; int y;};. . An L shaped connection is the set of points 
that lie on an L shaped path that goes from the first point in net to the second, excluding the start and end points. Int Orientation = 0  for an L shaped path
 that takes a turn at (x1,y2), and Orientation = 1 if it turns at (x2,y1). write a code that traverses first the path with a given orientation, then the path 
 with the other orientation, and then updates the orientation according to some function, and then traverses according to that orientation. In traversing a path with a given 
 orientation, take care that the traversal from 1 point to the other may need going downwards or upwards, an rightwards or leftwards.
 Also, if G is an object of class Grid_Graph {
public:
    int* Gx;
    int* Gy;
    int C;
    int M;
    int N;
    float v;
}
and the point (i,j) is associated with G.Gx[j*(G.N+1)+i] and G.Gy[i*(G.M+1)+j] if traversing in decreasing i or j,
and G.Gx[j*(G.N+1)+i+1] and G.Gy[i*(G.M+1)+j+1 if traversing in increasong i or j
while traversing vertically, decrement G.Gy of the corresponding point by 1 for each step, and while traversing
horizontally, decrement G.Gx of the corresponding point by 1 for each step
*/

float Batch::ripL(Grid_Graph G, Net net, int orientation) {
    //std::cerr << "Entering ripL function with orientation: " << orientation << std::endl;
    float costnew = 0;
    //std::cout << costnew << std::endl;
    if (orientation == 0) {
        // First traverse vertically, then horizontally (turn at x1, y2)
        // Vertical traversal
        if (net.y1 < net.y2) {
            for (int y = net.y1 + 1; y < net.y2; ++y) {
                int gy_index = net.x1 * (G.M + 1) + y+1;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
                //std::cout << "grid graph" << G.Gy[gy_index] << std::endl;
                //std::cout << costnew << std::endl;
            }
        } else {
            for (int y = net.y1 - 1; y > net.y2; --y) {
                int gy_index = net.x1 * (G.M + 1) + y;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
                //std::cout << "grid graph" << G.Gy[gy_index] << std::endl;
                //std::cout << costnew << std::endl;
            }
        }

        // Horizontal traversal
        if (net.x1 < net.x2) {
            for (int x = net.x1 + 1; x <= net.x2; ++x) {
                int gx_index = net.y2 * (G.N + 1) + x+1;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
                //std::cout << "grid graph" << G.Gx[gx_index] << std::endl;
                //std::cout << costnew << std::endl;
            }
        } 
        else {
            for (int x = net.x1 - 1; x >= net.x2; --x) {
                int gx_index = net.y2 * (G.N + 1) + x;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
                //std::cout << "grid graph" << G.Gx[gx_index] << std::endl;
                //std::cout << costnew << std::endl;
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
                //std::cout << "grid graph" << G.Gx[gx_index] << std::endl;
                //std::cout << costnew << std::endl;
            }
        }
        else {
            for (int x = net.x1 - 1; x > net.x2; --x) {
                int gx_index = net.y1 * (G.N + 1) + x;
                G.Gx[gx_index] -= 1;
                costnew += weight_pr(G.Gx[gx_index],G.C);
                //std::cout << "grid graph" << G.Gx[gx_index] << std::endl;
                //std::cout << costnew << std::endl;
            }
        }

        // Vertical traversal
        if (net.y1 < net.y2) {
            for (int y = net.y1 + 1; y <= net.y2; ++y) {
                int gy_index = net.x2 * (G.M + 1) + y+1;
                G.Gy[gy_index] -= 1;
                costnew += weight_pr(G.Gy[gy_index],G.C);
                //std::cout << "grid graph" << G.Gy[gy_index] << std::endl;
                //std::cout << costnew << std::endl;
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
    //std::cerr << "Exiting ripL with total cost: " << costnew << std::endl;
    //std::cout << costnew << std::endl;
    return costnew;
    //std::cout << "Hi anurag" << std::endl; 
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
            //G.Gy[net.x1*(G.M+1)+y] -= 1;
            cost += 1;
        }
        for (int x = net.x1 + 1; x <= net.x2; ++x) {
            //G.Gx[net.y2*(G.N+1)+x] -= 1;
            cost += 1;
        }
    } else {
        // Traverse the path with orientation 1 (turn at x2, y1)
        for (int x = net.x1 + 1; x < net.x2; ++x) {
            //G.Gx[net.y1*(G.N+1)+x] -= 1;
            cost += 1;
        }
        for (int y = net.y1; y <= net.y2; ++y) {
            //G.Gy[net.x2*(G.M+1)+y] -= 1;
            cost += 1;
        }
    }
    return cost;
}
// Function to perform maze route
void Batch::maze_route(Grid_Graph G, float k, float c,std::vector<float>& Sdist1,std::vector<char>&  Sdir1,std::vector<float>& Sdist2,std::vector<char>&  Sdir2, int NUM_THREADS, int BOX_MIN_DIM) {
    /*
    Look more carefully at if the below double allocation is really necessary
    */
   // First all variables necessary for maze routing are made
   /*
    x,y -> source
    x1,y1 -> destination
    xl,yl -> lower left corner of bounding box
    xu, yu -> upper left corner of bounding box
    Sdist -> stores distance from source to grid point
    Sdir -> stores direction to enter a cell
    Sdist, Sdir are allocated outside the function
    */ 
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

        //The following variables are created to benchmark timings in each #pragma
        double t_start, t_end;
        int Thread_num = omp_get_thread_num();
        int size = omp_get_num_threads();
	int N = nets.size();
        t_start = omp_get_wtime();
        #pragma omp for
        //schedule(dynamic,1)
        
        for (int i=0; i<N; i++){
            source = {nets[i].x1,nets[i].y1};
            dest = {nets[i].x2,nets[i].y2};
            rip_wire(G, source,dest,nets[i].route);
            }
    
        t_end = omp_get_wtime();

        #pragma omp critical
        std::cerr << "Time taken to rip nets " << t_end - t_start << "seconds by thread number " << Thread_num << std::endl;

        //std::cerr << "ripped up " << std::endl;
        t_start = omp_get_wtime();
        #pragma omp for //schedule(dynamic,1)
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
            //std::cout << "corner l " << cornerl.x << " " << cornerl.y << std::endl;
            //std::cerr << "cornerh " << cornerh.x << " " << cornerh.y << std::endl;
            //initializing Sdists, Sdirs
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
            //disp_s(Sdist1,Sdist2,Sdir1,Sdir2,G.M,G.N);
            // Then do Bellman-Ford
            bool flag = 1; // to keep track of if the relaxation step has caused any change to the distances and routes or not
            while (flag){
                //if (source.x==0&&source.y==0&&dest.x==3&&dest.y==3){
                //  disp_s(Sdist1,Sdist2,Sdir1,Sdir2,G.M,G.N);
                //}
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
                //disp_s(Sdist1,Sdist2,Sdir1,Sdir2,G.M,G.N);
                //std::cout << "i " << i << std::endl;
                //std::cout << "i " << i << std::endl;
                //std::cout << "i " << i << std::endl;
                //std::cout << "i " << i << std::endl;
                //std::cout << "i " << i << std::endl;
            }
            //std::cerr << "Bellmann ford relaxations are completed " << std::endl;
            

            // Backtrack and store the results
            // this has to be done completely in serial
        }
        t_end = omp_get_wtime();

        #pragma omp critical
        std::cerr << "Time taken to do Bellmann ford relaxations " << t_end - t_start << "seconds by thread number " << Thread_num << std::endl;

        t_start = omp_get_wtime();
        #pragma omp single
        {
        double t_start_single = omp_get_wtime();
	double hi = omp_get_wtime();
        for(int i=0; i<N;i++){
            //#pragma omp single
                // This save time as each core has its Sdist, etc in cache
                Point source = {nets[i].x1,nets[i].y1};
                Point dest = {nets[i].x2,nets[i].y2};
                Point here = dest;
                //std::cerr << "Dest is " << dest.x << " " << dest.y << std::endl;
                //exit(1);
                //std::cerr << "Source is " << source.x << " " << source.y << std::endl;
                //exit(1);
                int layer = 0;
                char dir = 'x';
                char diro = Sdir1[G.N*here.y+here.x]; // so that dest is not saved in route
                //std::cout << "Located at " << here.x << " " << here.y << std::endl;
                nets[i].route.clear();
                //std::cerr << "Variables for backtracking are initialized " << std::endl;
                //std::cerr << i << std::endl;
                int k = 0;
                while (!((here.x == source.x) && (here.y == source.y))){
                    // first check which direction to enter "here" from
                    //std::cout << "At " << here.x << " " << here.y << std::endl;
                    //std::cerr << "Source is " << source.x << " " << source.y << std::endl;
                    if (layer){
                        //std::cerr << "Reached a layer = 1 " << std::endl;
                        dir = Sdir2[G.M*here.x+here.y];
                        //std::cerr << "Crossed a layer = 1 " << std::endl;
                    }
                    else{
                        ///std::cerr << "Reached a layer = 0 " << std::endl;
                        dir = Sdir1[G.N*here.y+here.x];
                        //std::cerr << "Crossed a layer = 0 " << std::endl;
                    }
                    //std::cerr << "Checking if a bend is there  " << std::endl;
                    if ((dir != diro)&&(dir!='d')&&(dir!='u')&&!((here.x == dest.x) && (here.y == dest.y))){
                        //std::cerr << "Ready to print output " << std::endl;
                        nets[i].route.push_back(here);
                        //std::cerr << "Printed output " << std::endl;
                        //std:: cout << "dir " << dir << std::endl;
                        //std:: cout << "loc " << here.x << here.y << std::endl;
                    }
                    //std::cerr << "dir is " << dir << std::endl;
                    diro = dir;
                    if (dir=='x'){
                        k++;
                    }
                    else{k=0;}
                ///   std::cerr << "Dest is " << here.x << " " << dest.y << std::endl;
                //exit(1);
                //std::cerr << "Source is " << source.x << " " << source.y << std::endl;
                // exit(1);
                // now update "here" or "layer"
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
                //   std::cerr << "Dest is " << here.x << " " << dest.y << std::endl;
                //exit(1);
                //std::cerr << "Source is " << source.x << " " << source.y << std::endl;
                //exit(1);
                    
                    if (k==5){
                        //std:: cerr << "source " << source.x << " " << source.y << std::endl;
                        //std:: cerr << "dest " << dest.x << " " << dest.y << std::endl;
                        //std:: cerr << "at " << here.x << " " << here.y << " " << layer << std::endl;
                        //std:: cerr << " i " << i << std::endl;
                        //disp_s(Sdist1,Sdist2,Sdir1,Sdir2,G.M,G.N);
                        //exit(1);
                    }
                }
            
            //exit(1);
            }
        
        }

        double t_end_single = omp_get_wtime();
        //std::cerr << "Time taken to backtrack under #pragma omp single " << t_end_single - hi << "seconds by thread number " << Thread_num << std::endl;
        #pragma omp barrier
        t_end = omp_get_wtime();
        #pragma omp critical
        std::cerr << "Time taken to backtrack (very imp) " << t_end - t_start << "seconds by thread number " << Thread_num << std::endl;

        t_start = omp_get_wtime();
        #pragma omp for //schedule(dynamic,1)
        for(int i=0;i<N; i++){
            std::reverse(nets[i].route.begin(),nets[i].route.end());


            //std::cerr << "Reversed nets " << std::endl;
            /*
            
            if (1){//source.x==0&&source.y==0&&dest.x==3&&dest.y==3){
                for (int j = 0; j<nets[i].route.size();j++){
                    std::cerr << nets[i].route[j].x << " " << nets[i].route[j].y << std::endl;
                }
                std::cerr << std::endl;
            }
            */
            t_end = omp_get_wtime();
            std::cerr << "Time taken to reverse " << t_end - t_start << "seconds by thread number " << Thread_num << std::endl;
             t_start = omp_get_wtime();
        // update grid graph
            Point source = {nets[i].x1,nets[i].y1};
            Point dest = {nets[i].x2,nets[i].y2};
            route_wire(G, source,dest,nets[i].route);
        }
        t_end = omp_get_wtime();
        #pragma omp critical
        std::cerr << "Time taken to reroute " << t_end - t_start << "seconds by thread number " << Thread_num << std::endl;

    }

}


/*
we need to rip multiple bends
*/

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

void Batch::dispG(Grid_Graph G){
    //std::cout << "G.Gy" <<std::endl ;
    for (int i=0; i<G.M; i++){
        for (int j=0; j<G.N; j++){
            if (G.Gy[j*(G.M+1)+i]){
                //std::cout << G.Gy[j*(G.M+1)+i] << " " ;
            }
            else{
                //std::cout << 'y' << " " ;
            }
        }
        //std::cout << std::endl; 
    }
    //std::cout << "G.Gx" <<std::endl ;
    for (int i=0; i<G.M; i++){
        for (int j=0; j<G.N; j++){
            if (G.Gy[i*(G.N+1)+j]){
                //std::cout << G.Gy[i*(G.N+1)+j] << " " ;
            }
            else{
                //std::cout << 'y' << " " ;
            }
        }
        //std::cout << std::endl; 
    }
}
