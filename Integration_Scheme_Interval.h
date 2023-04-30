#pragma once
#include "Integratrion_Scheme.h"
#include <functional>

namespace Com_Methods
{
	class Integration_Scheme_Interval : protected Integration_Scheme
	{
	public:
		// ������������ �� ���� �������� ��� ������������ �������
		Integration_Scheme_Interval(Integration_Scheme_Type Type);
		// ����� ��� ���. ������������� ��������. begin � end - ���������� ��������������
		// num segments - ����� ���������,func - ��������������� �������
		double Calculate_Integral(const Point& begin, const Point& end, int Number_segments,
								  const std:: function<double(const Point& P)>& Func) const;
	};
}