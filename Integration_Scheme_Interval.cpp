#include "pch.h"
#include "Integration_Scheme_Interval.h"
#define A 1
namespace Com_Methods
{
	//конуструктору на вход подается тип квадратурной формулы
	Integration_Scheme_Interval::Integration_Scheme_Interval(Integration_Scheme_Type Type)
	{
		// заполнение массивов точек и весов интегриования
		switch (Type)
		{
			case Gauss1: 
			{
				Weight = { 2 };
				Points = { Point(0,0,0) };
				break;
			}
			
			case Gauss4:
			{
				Weight = { (18.0 - sqrt(30.0)) / 36,
						   (18.0 + sqrt(30.0)) / 36,
						   (18.0 - sqrt(30.0)) / 36,
						   (18.0 + sqrt(30.0)) / 36 };
				Points = { Point(-sqrt((3 + 2 * sqrt(6.0 / 5.0)) / 7.0), 0, 0),
						   Point(-sqrt((3 - 2 * sqrt(6.0 / 5.0)) / 7.0), 0, 0),
						   Point(sqrt((3 + 2 * sqrt(6.0 / 5.0)) / 7.0), 0, 0),
						   Point(sqrt((3 - 2 * sqrt(6.0 / 5.0)) / 7.0), 0, 0) };
				break;
			}
			case Parab:
			{
				Weight = { 1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 };
				Points = { Point(-1.0, 0, 0), Point(0,0,0), Point(1.0, 0, 0) };
				break;
			}			
		}
	}

	double Integration_Scheme_Interval::Calculate_Integral(
		const Point& Begin,
		const Point& End,
		int Number_Segments,
		const std::function<double(const Point& P)>& Func) const
	{
		//Результат (квадратурная сумма)
		double Result = 0.0;
		//Начальная точка сегмента
		double X0;
		//шаг на отрезке
		double h = ((End.x() - Begin.x()) / A) / Number_Segments;
		//сумма по всем сегментам разбиения
		for (int i = 0; i < Number_Segments; i++)
		{
			//начальная точка сегмента
			X0 = Begin.x() + i * h;
			//сумма по узлам интегрирования
			for (int Integ_Point = 0; Integ_Point < Points.size(); Integ_Point++)
			{
				//переход с мастер-элемента [-1, 1]
				auto P = Point(X0 + (1 + Points[Integ_Point].x()) * h / 2.0, 0, 0);
				Result += Weight[Integ_Point] * Func(P);
			}
		}
		//формируем результат с учетом якобиана на отрезке [-1, 1]
		return Result * (h / 2.0);
	}
}