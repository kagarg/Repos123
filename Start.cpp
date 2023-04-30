#include <iostream>
#include <iomanip>
#include <functional>
#include "Point.h"
#include "Integration_Scheme_Interval.h"
#define Pi 3.14
#define h 1
int main()
{
	//подынтегральная функция f(x) = sin x + cos x
	std::function<double(const Com_Methods::Point& P)> f =
		[](const Com_Methods::Point& P) { return (sin(P.x()) + cos(P.x())); };
	//первообразная F(x) = -cos x + sin x
	std::function<double(const Com_Methods::Point& P)> F =
		[](const Com_Methods::Point& P) { return (sin(P.x()) - cos(P.x())); };

	//квадратурная формула Гаусс-4
	Com_Methods::Integration_Scheme_Interval Quadrature_Formula1(Com_Methods::Integration_Scheme::Parab);
	Com_Methods::Integration_Scheme_Interval Quadrature_Formula2(Com_Methods::Integration_Scheme::Gauss4);

	//начало и конец отрезка интегрирования
	auto Begin = Com_Methods::Point(0, 0, 0);
	auto End = Com_Methods::Point(Pi/2, 0, 0);
	//число сегментов
	const int Num_Segments = 1;

	//точное значение интеграла (ф. Ньютона-Лейбница)
	double I_True = F(End) - F(Begin);

	//численное значение интеграла
	double I_gauss4[3];
	double I_gauss4h[3];
	std::cout << " Gauss 4 > " << std::endl;
	for (int i = 0; i < 3; i++) {
		I_gauss4[i] = Quadrature_Formula2.Calculate_Integral(Begin, End, Num_Segments *pow(2,i), f);
		I_gauss4h[i] = Quadrature_Formula2.Calculate_Integral(Begin, End,( Num_Segments * pow(2,i) )*2, f);

		std::cout << std::scientific;
		std::cout << "h = " << (End.x() - Begin.x()) /( Num_Segments * pow(2, i) )<< std::endl;
		/*std::cout << "Ih = " << I_gauss4[i] << std::endl;*/
		/*std::cout << "|I* - Ih| = " << fabs(I_gauss4[i] - I_True) << std::endl;*/
		int k = round(log2(fabs((1 + (I_gauss4h[i] - I_gauss4[i]) / (I_True - I_gauss4h[i])))));
		/*std::cout << "k = " << k << std::endl;*/
		
		double IR = I_gauss4h[i] + ((I_gauss4h[i] - I_gauss4[i]) / (pow(2, k) - 1));
		std::cout << "(I* - Ih) /((I* - Ih/2)" << ((I_True - I_gauss4[i]) / (I_True - I_gauss4h[i])) << std::endl;
		std::cout << "((I_gauss4h[i] - I_gauss4[i]) / (pow(2, k) - 1)) = " << ((I_gauss4h[i] - I_gauss4[i]) / (pow(2, k) - 1)) <<std::endl;
		std::cout << "Ir = " << IR << std::endl;
		std::cout << "I* - Ir = " << I_True - IR<< std:: endl;
		std::cout << "I* - Ih = " << I_True - I_gauss4[i] << std::endl;
		std::cout << std::endl;
	}
	
	std::cout << " Simpson > " << std::endl;
	double I_parab[3];
	double I_parabh[3];

	std::cout << I_True << std::endl;

	for (int i = 0; i < 3; i++) {
		I_parab[i] = Quadrature_Formula1.Calculate_Integral(Begin, End, Num_Segments * pow(2, i), f);
		I_parabh[i] = Quadrature_Formula1.Calculate_Integral(Begin, End, (Num_Segments * pow(2, i)) * 2, f);

		//std::cout << std::scientific;
		std::cout << "h = " << (End.x() - Begin.x()) / (Num_Segments * pow(2, i)) << std::endl;
		//std::cout << "I = " << I_parab[i] << std::endl;
		////std::cout << "Ih = " << I_parabh[i] << std::endl;
		//std::cout << "|I* - I| = " << fabs(I_parab[i] - I_True) << std::endl;

		int k = round(log2(fabs(1 + (I_parabh[i] - I_parab[i]) / (I_True - I_parabh[i]))));
		//std::cout << "k = " << k << std::endl;

		double IR = I_parabh[i] + ((I_parabh[i] - I_parab[i]) / (pow(2, k) - 1));
		//std::cout << "((I_parabh[i] - I_parab[i]) / (pow(2, k) - 1)) = " << ((I_parabh[i] - I_parab[i]) / (pow(2, k) - 1)) << std::endl;
		std::cout << "Ir = " << IR << std::endl;
		std::cout << "I* - Ir = " << fabs(I_True - IR) << std::endl;
		//std::cout << std::endl;
		std::cout <<"I* - Ih = " << I_True - I_parab[i]  << std::endl;
		std::cout << "(I* - Ih) /((I* - Ih/2)" << ((I_True - I_parab[i])  / (I_True - I_parabh[i])) << std::endl;
		std::cout << "(Ih/2 - Ih) /(2k-1)" << ((I_parabh[i] - I_parab[i]) / (2 * k - 1)) << std::endl;

	}
}