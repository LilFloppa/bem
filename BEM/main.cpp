#include <iostream>
#include <iomanip>
#include <vector>

#include "integrals.h"
#include "point.h"
#include "matrix.h"
#include "mesh.h"
#include "BEM.h"
#include "solver.h"


double U(double x, double y)
{
	return x;
}

void BasicDirichletProblem()
{
	int elementCount = 4;
	int nodeCount = 4;
	std::vector<Point> points;
	points.emplace_back(0, 0);
	points.emplace_back(1, 0);
	points.emplace_back(1, 1);
	points.emplace_back(0, 1);

	std::vector<BEMElement> elements;
	elements.emplace_back(0, 1, 0);
	elements.emplace_back(1, 2, 1);
	elements.emplace_back(2, 3, 2);
	elements.emplace_back(3, 0, 3);

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
		std::cout << pi << "\t";

	std::cout << std::endl;
}

void DirichletProblemWithNonRectDomain()
{
	int elementCount = 6;
	int nodeCount = 6;
	std::vector<Point> points;
	points.emplace_back(0, 0);
	points.emplace_back(0.8, 0);
	points.emplace_back(0.8, 0.7);
	points.emplace_back(1.0, 0.7);
	points.emplace_back(1.0, 1.0);
	points.emplace_back(0.0, 1.0);

	std::vector<BEMElement> elements;
	elements.emplace_back(0, 1, 0);
	elements.emplace_back(1, 2, 1);
	elements.emplace_back(2, 3, 2);
	elements.emplace_back(3, 4, 3);
	elements.emplace_back(4, 5, 4);
	elements.emplace_back(5, 0, 5);

	matrix<double> V(elementCount, std::vector<double>(elementCount));
	matrix<double> K(elementCount, std::vector<double>(nodeCount));
	matrix<double> D(nodeCount, std::vector<double>(nodeCount));

	BuildMatrix(points, elements, V, K, D);

	std::vector<double> q = { 0.0, 0.8, 0.8, 1.0, 1.0, 0.0 };

	std::vector<double> b(nodeCount);

	std::vector<double> p(elementCount);

	MultVector(K, q, b);
	Gauss(V, p, b);
}

void DirichetProblemWithDetailedMesh()
{
	int elementCount = 10;
	int nodeCount = 10;
	std::vector<Point> points;
	points.emplace_back(0, 0);
	points.emplace_back(0.4, 0);
	points.emplace_back(0.8, 0);
	points.emplace_back(0.8, 0.4);
	points.emplace_back(0.8, 0.7);
	points.emplace_back(1.0, 0.7);
	points.emplace_back(1.0, 1.0);
	points.emplace_back(0.4, 1.0);
	points.emplace_back(0.0, 1.0);
	points.emplace_back(0.0, 0.5);

	std::vector<BEMElement> elements;
	elements.emplace_back(0, 1, 0);
	elements.emplace_back(1, 2, 1);
	elements.emplace_back(2, 3, 2);
	elements.emplace_back(3, 4, 3);
	elements.emplace_back(4, 5, 4);
	elements.emplace_back(5, 6, 5);
	elements.emplace_back(6, 7, 6);
	elements.emplace_back(7, 8, 7);
	elements.emplace_back(8, 9, 8);
	elements.emplace_back(9, 0, 9);

	matrix<double> V(elementCount, std::vector<double>(elementCount));
	matrix<double> K(elementCount, std::vector<double>(nodeCount));
	matrix<double> D(nodeCount, std::vector<double>(nodeCount));

	BuildMatrix(points, elements, V, K, D);

	std::vector<double> q = { 0.0, 0.4, 0.8, 0.8, 0.8, 1.0, 1.0, 0.4, 0.0, 0.0 };

	std::vector<double> b(nodeCount);

	std::vector<double> p(elementCount);

	MultVector(K, q, b);
	Gauss(V, p, b);
}


// План
// 1. Сделать задачу Неймана
// 2. Протестировать задачу Неймана
// 3. Протестить функцию x^2 - y^2
// 4. Перенести построение матриц V, K, D на видеокарту

int main()
{
	auto segments = std::vector<Segment>
	{
		Segment{ { 0.0, 0.0 }, { 1.0, 0.0 },  2 },
		Segment{ { 1.0, 0.0 }, { 1.0, 1.0 },  2 },
		Segment{ { 1.0, 1.0 }, { 0.0, 1.0 },  2 },
		Segment{ { 0.0, 1.0 }, { 0.0, 0.0 },  2 },
	};

	auto points = GenerateMesh(segments);
	int nodeCount = points->size();

	std::vector<BEMElement> elements;

	for (int i = 0; i < nodeCount; i++)
		elements.emplace_back(i, i + 1, i);

	elements[elements.size() - 1].V2 = 0;

	int elementCount = elements.size();

	matrix<double> V(elementCount, std::vector<double>(elementCount));
	matrix<double> K(elementCount, std::vector<double>(nodeCount));
	matrix<double> D(nodeCount, std::vector<double>(nodeCount));
	BuildMatrix(*points, elements, V, K, D);

	matrix<double> M(nodeCount, std::vector<double>(nodeCount));
	matrix<double> Kt(nodeCount, std::vector<double>(nodeCount));
	BuildMatrixM(*points, elements, M);

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

	return 0;
}