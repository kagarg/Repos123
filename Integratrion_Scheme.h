#pragma once
#include "Point.h"
#include <vector>

namespace Com_Methods 
{
	//класс схемы интегрирования
	class Integration_Scheme
	{
	protected:
		// узлы интегрирования
		std::vector<Point> Points;
		//веса квадратурных формул
		std::vector<double> Weight;
	public:
		// типы квадратурных формул
		enum Integration_Scheme_Type
		{
			Gauss1 = 1,
			Gauss2,
			Gauss3,
			Gauss4,
			Gauss5,
			Trap,
			Parab
		};
	};
}