#pragma once

#include <functional>

#include "point.h"
#include "matrix.h"
#include "mesh.h"
#include "BEM.h"
#include "solver.h"
#include "solution.h"

void BasicDirichletProblem()
{
	auto segments = std::vector<Segment>
	{
		Segment{ { 0.0, 0.0 }, { 1.0, 0.0 },  2 },
		Segment{ { 1.0, 0.0 }, { 1.0, 1.0 },  2 },
		Segment{ { 1.0, 1.0 }, { 0.0, 1.0 },  2 },
		Segment{ { 0.0, 1.0 }, { 0.0, 0.0 },  2 },
	};

	auto& points = *GenerateMesh(segments);
	int nodeCount = points.size();

	std::vector<BEMElement> elements;

	for (int i = 0; i < nodeCount; i++)
		elements.emplace_back(i, i + 1, i);

	elements[elements.size() - 1].V2 = 0;
	int elementCount = elements.size();

	matrix<double> V(elementCount, std::vector<double>(elementCount));
	matrix<double> K(elementCount, std::vector<double>(nodeCount));
	matrix<double> D(nodeCount, std::vector<double>(nodeCount));

	BuildMatrix(points, elements, V, K, D);

	std::vector<double> q = { 0.0, 1.0, 1.0, 0.0 };

	std::vector<double> b(nodeCount);

	std::vector<double> p(elementCount);

	MultVector(K, q, b);
	Gauss(V, p, b);

	for (auto pi : p)
		std::cout << std::setprecision(6) << std::scientific << pi << "\t";

	std::cout << std::endl;
}

void DirichletProblemWithComplexDomain()
{
	auto segments = std::vector<Segment>
	{
		Segment{ { 0.0, 0.0 }, { 0.3, 0.0 },  2 },
		Segment{ { 0.3, 0.0 }, { 0.3, 0.8 },  2 },
		Segment{ { 0.3, 0.8 }, { 0.7, 0.8 },  2 },
		Segment{ { 0.7, 0.8 }, { 0.7, 0.0 },  2 },
		Segment{ { 0.7, 0.0 }, { 1.3, 0.0 },  2 },
		Segment{ { 1.3, 0.0 }, { 1.3, 0.8 },  2 },
		Segment{ { 1.3, 0.8 }, { 1.7, 0.8 },  2 },
		Segment{ { 1.7, 0.8 }, { 1.7, 0.0 },  2 },
		Segment{ { 1.7, 0.0 }, { 2.0, 0.0 },  2 },
		Segment{ { 2.0, 0.0 }, { 2.0, 2.0 },  2 },
		Segment{ { 2.0, 2.0 }, { 0.0, 2.0 },  2 },
		Segment{ { 0.0, 2.0 }, { 0.0, 0.0 },  2 },
	};

	auto& points = *GenerateMesh(segments);
	int nodeCount = points.size();

	std::vector<BEMElement> elements;

	for (int i = 0; i < nodeCount; i++)
		elements.emplace_back(i, i + 1, i);

	elements[elements.size() - 1].V2 = 0;
	int elementCount = elements.size();

	matrix<double> V(elementCount, std::vector<double>(elementCount));
	matrix<double> K(elementCount, std::vector<double>(nodeCount));
	matrix<double> D(nodeCount, std::vector<double>(nodeCount));


	std::function<double(double, double)> u = [](double x, double y) {
		return x + y;
	};

	BuildMatrix(points, elements, V, K, D);

	std::vector<double> q(nodeCount);
	for (int i = 0; i < nodeCount; i++)
		q[i] = u(points[i].X, points[i].Y);

	std::vector<double> b(nodeCount);

	std::vector<double> p(elementCount);

	MultVector(K, q, b);
	Gauss(V, p, b);

	for (auto pi : p)
		std::cout << std::setprecision(6)  << std::scientific << pi << "\n";

	std::cout << std::endl;
}

void DirichletSquareFunction()
{
	int segmentPointCount = 64;
	auto segments = std::vector<Segment>
	{
		Segment{ { 0.0, 0.0 }, { 0.3, 0.0 },  segmentPointCount },
		Segment{ { 0.3, 0.0 }, { 0.3, 0.8 },  segmentPointCount },
		Segment{ { 0.3, 0.8 }, { 0.7, 0.8 },  segmentPointCount },
		Segment{ { 0.7, 0.8 }, { 0.7, 0.0 },  segmentPointCount },
		Segment{ { 0.7, 0.0 }, { 1.3, 0.0 },  segmentPointCount },
		Segment{ { 1.3, 0.0 }, { 1.3, 0.8 },  segmentPointCount },
		Segment{ { 1.3, 0.8 }, { 1.7, 0.8 },  segmentPointCount },
		Segment{ { 1.7, 0.8 }, { 1.7, 0.0 },  segmentPointCount },
		Segment{ { 1.7, 0.0 }, { 2.0, 0.0 },  segmentPointCount },
		Segment{ { 2.0, 0.0 }, { 2.0, 2.0 },  segmentPointCount },
		Segment{ { 2.0, 2.0 }, { 0.0, 2.0 },  segmentPointCount },
		Segment{ { 0.0, 2.0 }, { 0.0, 0.0 },  segmentPointCount },
	};

	auto& points = *GenerateMesh(segments);
	int nodeCount = points.size();

	std::vector<BEMElement> elements;

	for (int i = 0; i < nodeCount; i++)
		elements.emplace_back(i, i + 1, i);

	elements[elements.size() - 1].V2 = 0;
	int elementCount = elements.size();

	matrix<double> V(elementCount, std::vector<double>(elementCount));
	matrix<double> K(elementCount, std::vector<double>(nodeCount));
	matrix<double> D(nodeCount, std::vector<double>(nodeCount));


	std::function<double(double, double)> u = [](double x, double y) {
		return x * x - y * y;
	};
	std::function<Point(double, double)> gradu = [](double x, double y) {
		return Point(2 * x, -2 * y);
	};

	BuildMatrix(points, elements, V, K, D);

	std::vector<double> realq(nodeCount);
	std::vector<double> realp(elementCount);
	for (int i = 0; i < nodeCount; i++)
	{
		realq[i] = u(points[i].X, points[i].Y);
		Point m = ElementMiddle(points, elements[i]);
		Point n = ElementNormal(points, elements[i]);
		realp[i] = gradu(m.X, m.Y).Dot(n);
	}

	std::vector<double> b(nodeCount);

	std::vector<double> p(elementCount);

	MultVector(K, realq, b);
	Gauss(V, p, b);

	double sum = 0;
	for (int i = 0; i < nodeCount; i++)
		sum += (realp[i] - p[i]) * (realp[i] - p[i]);

	std::cout << std::setw(20) << std::left << "Element Count: " << elementCount << std::endl;
	std::cout << std::setw(20) << std::left << "R: " << std::sqrt(sum) << std::endl;
}