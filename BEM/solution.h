#pragma once

#include <vector>
#include <cmath>

#include "constants.h"
#include "point.h"
#include "BEM.h"

double V(double x, double Ax, double Ay)
{
    double Rx = Ax - x, Ry = Ay;
    double logR = log(Rx * Rx + Ry * Ry);

    double A = -2.0 * x;
    double B = -Rx;
    double C = -2.0 * Ay;

    double res = A;

    if (fabs(B) > Eps)
        res += B * logR;

    if (fabs(Ry) > Eps)
        res += C * atan(Rx / Ry);

    return res;
}

double K1(double x, double Ax, double Ay)
{
    double Rx = Ax - x, Ry = Ay;

    double C = -1.0;

    double res = 0.0;

    if (fabs(Ry) > Eps)
        res += C * atan(Rx / Ry);

    return res;
}

double Kx(double x, double Ax, double Ay)
{
    double Rx = Ax - x, Ry = Ay;
    double logR = log(Rx * Rx + Ry * Ry);

    double B = 1.0 / 2.0 * Ay;
    double C = -Ax;

    double res = 0.0;

    if (fabs(B) > Eps)
        res += B * logR;

    if (fabs(Ry) > Eps)
        res += C * atan(Rx / Ry);

    return res;
}

double GetValue(
    std::vector<Point>& points, 
    std::vector<BEMElement>& elements,
    std::vector<double>& q,
    std::vector<double>& p, 
    Point& Y)
{
    double result = 0;

    for (auto& el : elements)
    {
        double pvalue = p[el.P];
        double q1 = q[el.V1];
        double q2 = q[el.V2];


        Point A = points[el.V1];
        Point B = points[el.V2];

        Point vec = B - A;
        double h = vec.Norm();
        vec = vec / h;

        double Ax = vec.Dot(Y - A);
        double Ay = -vec.Orthogonal().Dot(Y - A);

        double v = -0.5 * PI_2 * (V(h, Ax, Ay) - V(0, Ax, Ay));
        double k1 = -PI_2 * (K1(h, Ax, Ay) - K1(0, Ax, Ay));
        double kx = -PI_2 / h * (Kx(h, Ax, Ay) - Kx(0, Ax, Ay));

        result += v * pvalue - kx * q1 - (k1 - kx) * q2;
    }

    return result;
}
