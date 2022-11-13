#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

template <typename T>
using matrix = std::vector<std::vector<double>>;

void MultVector(matrix<double>& A, std::vector<double>& v, std::vector<double>& result)
{
	int n = A.size();

	for (int i = 0; i < n; i++)
		result[i] = 0.0;

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			result[i] += A[i][j] * v[j];
}

void Copy(matrix<double>& src, matrix<double>& dst)
{
	dst.resize(src.size());

	for (int i = 0; i < src.size(); i++)
	{
		dst[i].resize(src[i].size());
		for (int j = 0; j < src[i].size(); j++)
			dst[i][j] = src[i][j];
	}
}

void Print(matrix<double>& A)
{
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
			std::cout << std::setprecision(4) << std::setw(10) << A[i][j];

		std::cout << std::endl;
	}
}

void Transpose(matrix<double>& A)
{
	int n = A.size();
	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
			std::swap(A[i][j], A[j][i]);
}