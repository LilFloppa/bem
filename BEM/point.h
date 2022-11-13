#pragma once

#include <cmath>

struct Point
{
	double X, Y;

	Point() : X(0), Y(0) {}
	Point(double x, double y) : X(x), Y(y) {}

	Point operator+(Point other)
	{
		return Point(X + other.X, Y + other.Y);
	}

	Point operator-(Point other)
	{
		return Point(X - other.X, Y - other.Y);
	}

	Point operator*(double a)
	{
		return Point(X * a, Y * a);
	}

	Point operator / (double a)
	{
		return Point(X / a, Y / a);
	}

	double Dot(Point other)
	{
		return X * other.X + Y * other.Y;
	}

	Point Normailze()
	{
		double n = Norm();
		return Point(X / n, Y / n);
	}

	double Norm()
	{
		return std::sqrt(X * X + Y * Y);
	}

	Point Orthogonal()
	{
		return Point(Y, -X);
	}
};