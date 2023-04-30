#pragma once
#include "Integratrion_Scheme.h"
#include <functional>

namespace Com_Methods
{
	class Integration_Scheme_Interval : protected Integration_Scheme
	{
	public:
		// конструктору на вход подается тип квадратурной формулы
		Integration_Scheme_Interval(Integration_Scheme_Type Type);
		// метод для выч. определенного интерала. begin и end - промежуток интегрирования
		// num segments - число сегментов,func - подынтегральная функция
		double Calculate_Integral(const Point& begin, const Point& end, int Number_segments,
								  const std:: function<double(const Point& P)>& Func) const;
	};
}