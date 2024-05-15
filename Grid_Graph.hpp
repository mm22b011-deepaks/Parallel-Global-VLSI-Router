
class Grid_Graph {
public:
    int* Gx;
    int* Gy;
    int C;
    int M;
    int N;
    float v;
    // Constructor
    Grid_Graph(int M_val, int N_val, int C_val, float v_val);
    void boundary_cond();
    // Function to calculate number of elements greater than C
    int overflows();
};