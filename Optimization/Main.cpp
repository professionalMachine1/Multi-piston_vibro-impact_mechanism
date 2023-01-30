#include "SDE.h"

int main()
{
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(0);

	/*

		Обязательно должно быть выполнено следующее:
		gamma[0] = 1, fi[0] = 0, eps[cur_index] << 1, mu << 1, gamma[cur_index] > 1, sum(lambda[cur_index], cur_index=1,...,N) < 1

	*/

	Study_DE var(2, 500, 100);
	float_type pMin, pMax;
	float_type pMin_Out, pMax_Out, pMin_In, pMax_In;

	var.t_start = 0, var.dx_start = 0.2;
	var.R = 0.5, var.p = 0.13, var.mu = 0.12;
	switch (var.num_of_pists)
	{
	case 1:
		var.fi[0] = 0.0;
		var.gamma[0] = 1;
		var.eps[0] = 0.005;
		var.lambda[0] = 0.2;
		break;
	case 2:
		var.fi[0] = 0, var.fi[1] = 0.2;
		var.gamma[0] = 1, var.gamma[1] = 2;
		var.eps[0] = 0.005, var.eps[1] = 0.005;
		var.lambda[0] = 0.2, var.lambda[1] = 0.2;
		break;
	case 3:
		var.fi[0] = 0, var.fi[1] = 0.5 * PI, var.fi[2] = PI;
		var.gamma[0] = 1, var.gamma[1] = 1.5, var.gamma[2] = 2;
		var.eps[0] = 0.005, var.eps[1] = 0.005, var.eps[2] = 0.005;
		var.lambda[0] = 0.2, var.lambda[1] = 0.2, var.lambda[2] = 0.2;
		break;
	}

	/*

		По умолчанию h_one = 0.001, min_speed = 0.0001
		При желании эти параметры можно менять

	*/

	// Ищем время старта с поверхности f_1
	// start_time();

	// Версия с многократным решением дифференциального уравнения
	var.h_one = pow(10, -3);
	var.h_out = pow(10, -2), var.h_in = pow(10, -3);
	pMin = 0.12, pMax = 0.15;
	pMin_Out = 0.1, pMax_Out = 0.8, pMin_In = 0.1, pMax_In = 0.4;

	// var.oscillogram_MSV(0.01, 50);
	// var.bifurcation_diagramm_MSV(pMin, pMax);
	// var.dynamic_mode_map_MSV(pMin_Out, pMax_Out, pMin_In, pMax_In);

	// Версия с точечным отображением
	var.h_one = pow(10, -3);
	var.h_out = pow(10, -2), var.h_in = pow(10, -3);
	pMin = 0.12, pMax = 0.15;
	pMin_Out = 0.1, pMax_Out = 0.8, pMin_In = 0.1, pMax_In = 0.4;

	// var.point_map_DMV();
	// var.bifurcation_diagramm_DMV(pMin, pMax);
	// var.lyapunov_diagramm_DMV(pMin, pMax);
	// var.dynamic_mode_map_DMV(pMin_Out, pMax_Out, pMin_In, pMax_In);

	return 0;
}
