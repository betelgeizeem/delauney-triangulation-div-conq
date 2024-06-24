#pragma once

#include <iostream>

struct double2
{
    double2(double x, double y) : x(x), y(y) {}
    double2(double a) : x(a), y(a) {}
    double2() : x(0), y() {}
    double x;
    double y;
};

double2 min(double2 a, double2 b);
double2 max(double2 a, double2 b);

double2 operator-(double2 a);
double2 operator-(double2 a, double2 b);
double2 operator+(double2 a, double2 b);

double crossProductZ(double2 vec1, double2 vec2);

std::ostream& operator<<(std::ostream& s, double2 a);

double length(double2 a);
double length2(double2 a);

double dotProduct(double2 a, double2 b);