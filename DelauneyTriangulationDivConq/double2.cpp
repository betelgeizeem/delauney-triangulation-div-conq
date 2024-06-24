#include "double2.h"

// double2 operations
double2 min(double2 a, double2 b) { return { std::min(a.x, b.x), std::min(a.y, b.y) }; }
double2 max(double2 a, double2 b) { return { std::max(a.x, b.x), std::max(a.y, b.y) }; }

double2 operator-(double2 a) { return { -a.x , -a.y }; }
double2 operator-(double2 a, double2 b) { return { a.x - b.x , a.y - b.y }; }
double2 operator+(double2 a, double2 b) { return { a.x + b.x , a.y + b.y }; }

double crossProductZ(double2 vec1, double2 vec2) { return vec1.x * vec2.y - vec2.x * vec1.y; }

std::ostream& operator<<(std::ostream& s, double2 a) { return (s << a.x << " " << a.y); }

double length(double2 a) { return std::sqrt(a.x * a.x + a.y * a.y); }
double length2(double2 a) { return a.x * a.x + a.y * a.y; }

double dotProduct(double2 a, double2 b) { return a.x * b.x + a.y * b.y; }