#pragma once

#include <vector>

#include "matrix.h"

double DotProduct(std::vector<double>& a, std::vector<double>& b)
{
	int n = a.size();
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += a[i] * b[i];

	return sum;
}

void Gauss(matrix<double>& A, std::vector<double>& x, std::vector<double>& b)
{
	int N = A.size();

	for (int k = 0; k < N - 1; k++)
	{
		double max = abs(A[k][k]);
		int m = k;
		for (int i = k + 1; i < N; i++)
			if (abs(A[k][i]) > max)
			{
				max = abs(A[k][i]);
				m = i;
			}

		std::swap(b[m], b[k]);
		for (int j = k; j < N; j++)
			std::swap(A[k][j], A[m][j]);

		for (int i = k + 1; i < N; i++)
		{
			double t = A[i][k] / A[k][k];
			b[i] -= t * b[k];
			for (int j = k + 1; j < N; j++)
				A[i][j] -= t * A[k][j];
		}
	}

	x[N - 1] = b[N - 1] / A[N - 1][N - 1];
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
			sum += A[k][j] * x[j];

		x[k] = (b[k] - sum) / A[k][k];
	}
}


void LOS(matrix<double>& A, std::vector<double>& x, std::vector<double>&  f, int maxiter, double eps)
{
	int N = A.size();
	std::vector<double> Ax(A.size());
	std::vector<double> r(A.size());
	std::vector<double> z(A.size());
	std::vector<double> p(A.size());

	// Calculate r0, z0
	MultVector(A, x, Ax);
	for (int i = 0; i < N; i++)
	{
		r[i] = f[i] - Ax[i];
		z[i] = r[i];
	}

	// Calculate p0
	MultVector(A, z, p);

	double diff = DotProduct(r, r);
	for (int k = 0; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(p, p);
		double a = DotProduct(p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		MultVector(A, r, Ax);
		double b = -DotProduct(p, Ax) / dotP;

		// Calculate zk, pk
		for (int i = 0; i < N; i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(r, r);
	}
}