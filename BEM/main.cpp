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
	// NeumannSquareFunction();
	DirichletSquareFunction();
	return 0;
}