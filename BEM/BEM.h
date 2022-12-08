#pragma once

#include <vector>

#include "point.h"
#include "matrix.h"
#include "integrals.h"

struct BEMElement
{
	int V1, V2, P;

	BEMElement(int v1, int v2, int p) : V1(v1), V2(v2), P(p) {}
};

Point ElementMiddle(std::vector<Point>& points, BEMElement& el)
{
	Point a = points[el.V1];
	Point b = points[el.V2];
	return (a + b) / 2;
}

Point ElementNormal(std::vector<Point>& points, BEMElement& el)
{
	Point a = points[el.V1];
	Point b = points[el.V2];

	Point v = (b - a).Normailze();
	return v.Orthogonal();
}

void BuildMatrix(std::vector<Point>& points, std::vector<BEMElement>& elements, matrix<double>& V, matrix<double>& K, matrix<double>& D)
{
	int elementCount = elements.size();

	for (auto& e1 : elements)
	{
		for (auto& e2 : elements)
		{
			Point p1 = points[e1.V1];
			Point p2 = points[e2.V1];

			Point v1 = (points[e1.V2] - points[e1.V1]);
			Point v2 = (points[e2.V2] - points[e2.V1]);

			double h1 = v1.Norm();
			double h2 = v2.Norm();

			v1 = v1 / h1;
			v2 = v2 / h2;

			double cos = v1.Dot(v2);
			double sin = v1.Dot(v2.Orthogonal());

			double Ax = v1.Dot(p2 - p1);
			double Ay = -v1.Orthogonal().Dot(p2 - p1);

			auto f = Fv(Ax, Ay, sin, cos);
			double v = -0.5 * PI_2 * (f(h1, h2) + f(0, 0) - f(h1, 0) - f(0, h2));
			V[e1.P][e2.P] = v;

			double d = v / (h1 * h2);
			D[e1.V1][e2.V1] += d;
			D[e1.V1][e2.V2] -= d;
			D[e1.V2][e2.V1] -= d;
			D[e1.V2][e2.V2] += d;

			if (std::abs(Ay) > Eps || std::abs(sin) > Eps)
			{
				f = Fkx(Ax, Ay, sin, cos);
				double kx = -PI_2 / h1 * (f(h1, h2) + f(0, 0) - f(h1, 0) - f(0, h2));

				f = Fk1(Ax, Ay, sin, cos);
				double k1 = -PI_2 * (f(h1, h2) + f(0, 0) - f(h1, 0) - f(0, h2));

				K[e2.P][e1.V1] += k1 - kx;
				K[e2.P][e1.V2] += kx;
			}

			if (e1.V1 == e2.V1 && e1.V2 == e2.V2)
			{
				K[e2.P][e1.V1] += h1 / 4.0;
				K[e2.P][e1.V2] += h1 / 4.0;
			}
		}
	}
}

void BuildMatrixM(std::vector<Point>& points, std::vector<BEMElement>& elements, matrix<double>& M)
{
	int elementCount = elements.size();

	for (auto& e : elements)
	{
		Point a = points[e.V1];
		Point b = points[e.V2];

		double h = (b - a).Norm();

		M[e.P][e.V1] += h / 2.0;
		M[e.P][e.V2] += h / 2.0;
	}
}
