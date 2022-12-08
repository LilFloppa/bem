#pragma once

#include <functional>
#include <cmath>

#include "constants.h"

std::function<double(double, double)>  Fv(double Ax, double Ay, double sin, double cos) {
	return [=](double x, double t)
	{
		double Rx = (Ax + cos * t - x);
		double Ry = (Ay + sin * t);
		double R = Rx * Rx + Ry * Ry;
		double result = -3 * x * t;
		double temp = 0.5 * (-x * x + 2 * Ax * x - t * t - Ax * Ax + Ay * Ay) * cos + (Ay * sin + t) * (x - Ax);

		if (std::abs(temp) > Eps)
			result += temp * std::log(R);

		if (std::abs(sin) < Eps)
		{
			if (std::abs(Ry) > Eps)
			{
				result += 2 * ((-Ax + x) * cos - t) * Ay * std::atan(Rx / Ry);
			}
		}
		else
		{
			double a = (cos * Ry - sin * Rx);
			double temp = a * a / sin;
			if (std::abs(temp) > Eps)
				result += temp * std::atan((cos * Rx + sin * Ry) / a);

			if (std::abs(Ry) > Eps)
				result += -Ry * Ry * std::atan(Rx / Ry) / sin;
		}

		return result;
	};
}

std::function<double(double, double)>  Fk1(double Ax, double Ay, double sin, double cos) {
	return [=](double x, double t)
	{
		double Rx = (Ax + cos * t - x);
		double Ry = (Ay + sin * t);
		double R = Rx * Rx + Ry * Ry;
		double result = (0.5 * cos * Ay + 0.5 * sin * (x - Ax));

		double arctan = std::atan(Rx / Ry);

		if (std::abs(result) > Eps)
			result *= std::log(R);

		if (std::abs(sin) < Eps)
		{
			double temp = cos * (x - Ax) - t;
			if (std::abs(temp) > Eps)
				result += temp * arctan;
		}
		else
		{
			double a = (cos * Ry - sin * Rx);

			if (std::abs(Ry) > Eps)
				result += (-Ry * arctan) / sin;

			if (std::abs(cos * a) > Eps)
				result += cos * a * std::atan((cos * Rx + sin * Ry) / a) / sin;
		}

		return result;
	};
}

std::function<double(double, double)>  Fkx(double Ax, double Ay, double sin, double cos) {
	return [=](double x, double t)
	{
		double Rx = (Ax + cos * t - x);
		double Ry = (Ay + sin * t);
		double R = Rx * Rx + Ry * Ry;

		double result = (0.25 * (-Ax * Ax + Ay * Ay + t * t + x * x) * sin + 0.5 * Ay * (t + Ax * cos));

		if (std::abs(result) > Eps)
			result *= std::log(R);

		double arctan = std::atan(Rx / Ry);

		if (std::abs(sin) < Eps)
		{
			double temp = (0.5 * (-t * t - Ax * Ax + x * x + Ay * Ay) * cos - Ax * t);
			if (std::abs(temp) > Eps)
				result += temp * arctan;
		}
		else
		{
			double a = (cos * Ry - sin * Rx);
			double b = (sin * Ax - cos * Ay);

			double temp = (cos * b + 0.5 * cos * a) * a / (sin * sin);
			if (std::abs(temp) > Eps)
				result += temp * std::atan((cos * Rx + sin * Ry) / a);
			
			temp = (b + 0.5 * cos * Ry) * Ry / (sin * sin);
			if (std::abs(temp) > Eps)
				result += -temp * arctan;
		}

		return result;
	};
}