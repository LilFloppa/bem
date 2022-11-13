#pragma once

#include <vector>

#include "point.h"

struct Segment
{
	Point Begin, End;
	int PointCount;
};

std::vector<Point>* GenerateInterval(Segment& segment)
{
	Point a = segment.Begin;
	Point b = segment.End;

	if (segment.PointCount < 2)
		throw new std::invalid_argument("pointCount can't be less than 2");

	std::vector<Point>* result = new std::vector<Point>(segment.PointCount);

	Point step = (b - a) / (segment.PointCount - 1);

	for (int i = 0; i < segment.PointCount; i++)
		(*result)[i] = a + step * i;

	return result;
}


std::vector<Point>* GenerateMesh(std::vector<Segment>& segments)
{
	int totalPointCount = 0;
	for (auto& s : segments)
		totalPointCount += s.PointCount - 1;

	std::vector<Point>* result = new std::vector<Point>(totalPointCount);

	int index = 0;
	for (auto& s : segments)
	{
		auto points = GenerateInterval(s);
		for (int i = 0; i < points->size() - 1; i++)
			(*result)[index++] = (*points)[i];
	}

	return result;
}