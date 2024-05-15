#include "main.hpp"
#include "Grid_Graph.hpp"
#include "struct.hpp"

// Constructor
Grid_Graph::Grid_Graph(int M_val, int N_val, int C_val, float v_val) : M(M_val), N(N_val), C(C_val), v(v_val) {
    Gx = (int*)calloc((N + 1) * M, sizeof(int));
    Gy = (int*)calloc((M + 1) * N, sizeof(int));
}


void Grid_Graph::boundary_cond(){
    for (int i = 0; i<N;i++){
        Gy[(M+1)*i] = 0;
        Gy[(M+1)*(i+1)-1] = 0;
    }
    for (int i = 0; i<M;i++){
        Gx[(N+1)*i] = 0;
        Gx[(N+1)*(i+1)-1] = 0;
    }
}
// calculate number of electrical shorts in the chip
int Grid_Graph::overflows() {
    int count = 0;
    for (int i = 0; i < (N + 1) * M; i++) {
        if (Gx[i] > C) {
            count+=Gx[i]-C;
        }
    }
    for (int i = 0; i < (M + 1) * N; i++) {
        if (Gy[i] > C) {
            count+=Gx[i]-C;
        }
    }
    return count;
}
