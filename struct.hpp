#ifndef STRUCT_HPP
#define STRUCT_HPP

// struct Point
struct Point {
    int x;
    int y;
};



// Struct representing a net
struct Net {
    int x1, y1, x2, y2;
    std::vector<Point> route;

    // Default constructor
    Net() : x1(0), y1(0), x2(0), y2(0) {}

    // Parameterized constructor
    Net(int v1, int v2, int v3, int v4, size_t size)
        : x1(v1), y1(v2), x2(v3), y2(v4) {
        route.reserve(size);
    }
};
#endif // STRUCT_HPP


