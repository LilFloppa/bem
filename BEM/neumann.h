#pragma once

#include <functional>

#include "point.h"
#include "matrix.h"
#include "mesh.h"
#include "BEM.h"
#include "solver.h"
#include "solution.h"

void BasicNeumannProblem()
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

	matrix<double> M(nodeCount, std::vector<double>(nodeCount));
	matrix<double> Kt(nodeCount, std::vector<double>(nodeCount));
	BuildMatrixM(points, elements, M);

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = K[j][i];

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = M[j][i] - K[j][i];


	std::vector<double> q = { 0.0, 1.0, 1.0, 0.0 };
	std::vector<double> p = { 0.0, 1.0, 0.0, -1.0 };

	std::vector<double> Ktp(nodeCount);
	std::vector<double> Dq(nodeCount);

	MultVector(Kt, p, Ktp);
	MultVector(D, q, Dq);

	double sum = 0;
	for (int i = 0; i < nodeCount; i++)
		sum += (Ktp[i] - Dq[i]) * (Ktp[i] - Dq[i]);

	std::cout << "R: " << std::sqrt(sum) << std::endl;

	for (auto qi : q)
		std::cout << std::setprecision(6) << std::scientific << qi << "\t";
}

void NeumannProblemWithComplexDomain()
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
	BuildMatrix(points, elements, V, K, D);

	matrix<double> M(nodeCount, std::vector<double>(nodeCount));
	matrix<double> Kt(nodeCount, std::vector<double>(nodeCount));
	BuildMatrixM(points, elements, M);

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = K[j][i];

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = M[j][i] - K[j][i];


	std::function<double(double, double)> u = [](double x, double y) {
		return x + y;
	};

	std::function<Point(double, double)> gradu = [](double x, double y) {
		return Point(1, 1);
	};

	std::vector<double> q(nodeCount);
	std::vector<double> p(elementCount);
	for (int i = 0; i < nodeCount; i++)
	{
		q[i] = u(points[i].X, points[i].Y);
		Point m = ElementMiddle(points, elements[i]);
		Point n = ElementNormal(points, elements[i]);
		p[i] = gradu(m.X, m.Y).Dot(n);
	}

	std::vector<double> Ktp(nodeCount);
	std::vector<double> Dq(nodeCount);

	MultVector(Kt, p, Ktp);
	MultVector(D, q, Dq);

	double sum = 0;
	for (int i = 0; i < nodeCount; i++)
		sum += (Ktp[i] - Dq[i]) * (Ktp[i] - Dq[i]);

	std::cout << "R: " << std::sqrt(sum) << std::endl;

	for (auto qi : q)
		std::cout << std::setprecision(6) << std::scientific << qi << "\t";
}

void NeumannSquareFunction()
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
	BuildMatrix(points, elements, V, K, D);

	matrix<double> M(nodeCount, std::vector<double>(nodeCount));
	matrix<double> Kt(nodeCount, std::vector<double>(nodeCount));
	BuildMatrixM(points, elements, M);

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = K[j][i];

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = M[j][i] - K[j][i];


	std::function<double(double, double)> u = [](double x, double y) {
		return x * x - y * y;
	};

	std::function<Point(double, double)> gradu = [](double x, double y) {
		return Point(2 * x, -2 * y);
	};

	std::vector<double> q(nodeCount);
	std::vector<double> p(elementCount);
	for (int i = 0; i < nodeCount; i++)
	{
		q[i] = u(points[i].X, points[i].Y);
		Point m = ElementMiddle(points, elements[i]);
		Point n = ElementNormal(points, elements[i]);
		p[i] = gradu(m.X, m.Y).Dot(n);
	}

	std::vector<double> Ktp(nodeCount);
	std::vector<double> Dq(nodeCount);

	MultVector(Kt, p, Ktp);
	MultVector(D, q, Dq);

	double sum = 0;
	for (int i = 0; i < nodeCount; i++)
		sum += (Ktp[i] - Dq[i]) * (Ktp[i] - Dq[i]);

	std::cout << "R: " << std::sqrt(sum) << std::endl;
}