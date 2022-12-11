#include <iostream>
#include <iomanip>
#include <vector>

#include "integrals.h"
#include "point.h"
#include "matrix.h"
#include "mesh.h"
#include "BEM.h"
#include "solver.h"
#include "solution.h"

#include "dirichlet.h"
#include "neumann.h"

int main()
{
	// BasicDirichletProblem();
	// BasicNeumannProblem();

	// DirichletProblemWithComplexDomain();
	// NeumannProblemWithComplexDomain();
	
	int multiplier = 2;
	for (int i = 0; i < 6; i++)
	{
		DirichletSquareFunction(multiplier);
		multiplier *= 2;
	}

	std::cout << std::endl << std::endl;

	multiplier = 2;
	for (int i = 0; i < 6; i++)
	{
		NeumannSquareFunction(multiplier);
		multiplier *= 2;
	}

	return 0;
}