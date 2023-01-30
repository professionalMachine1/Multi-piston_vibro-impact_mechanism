#pragma once

#ifndef SDE_H
#define SDE_H

#include "ODE.h"
#include "Directives.h"

const std::string settings_term = "set terminal pngcairo size 1080,720 enhanced font 'Times New Roman,20'\n";
const std::string settings_palette = "set palette rgbformulae 23, 13, 10\n";

vector ode_function(const vector& x, float_type t, const vector& var)
{
	// var: {A: 0, FI: 1, p: 2}
	vector product(x.size());
	product[0] = x[1];
	product[1] = var[0] * cos(t - var[1]) - var[2];
	return product;
}

void Plot(const std::string& str)
{
	FILE* gpipe = _popen("gnuplot\\bin\\gnuplot.exe -persist", "w");
	fprintf(gpipe, str.c_str());
	_pclose(gpipe);
}

/*

	----------------------------------

*/

// P_ONE - одиночный цикл, P_OUT - внешний цикл, P_IN - внутренний цикл
#define P_ONE p
std::string SP_ONE = "p";

#define P_OUT R
std::string SP_OUT = "R";
#define P_IN p
std::string SP_IN = "p";

/*

	----------------------------------

*/

class Study_DE
{
public:
	int_type num_of_pists;
	float_type h_one, h_in, h_out, min_speed;
	size_t max_iterations, add_iterations;
	vector eps, gamma, fi, lambda, tau, xder;
	float_type R, p, A, FI, a, b, mu, dx_start, t_start;
public:
	Study_DE(int_type num_pists = 2, size_t max_iter = 500, size_t add_iter = 100);

	void reset_parametrs();
	void start_time_search();
	// ---
	float_type max_eigv();
	float_type dpm_der_value(int_type, int_type);
	float_type surface_der_val(int_type, float_type);
	std::tuple<float_type, int_type> surface_val_ind(float_type);
	// Методы, связанные с многократным решение дифференциального уравнения
	void oscillogram_MSV(float_type step, float_type interval);
	void bifurcation_diagramm_MSV(float_type pMin, float_type pMax);
	void dynamic_mode_map_MSV(float_type pMin_out, float_type pMax_out,
		float_type pMin_in, float_type pMax_in);
	// ---
	std::tuple<int_type, bool> next_point_de(DE_SOLVE& solve);
	// Методы, связанные с точечным отображением
	void point_map_DMV();
	void bifurcation_diagramm_DMV(float_type pMin, float_type pMax);
	void lyapunov_diagramm_DMV(float_type pMin, float_type pMax);
	void dynamic_mode_map_DMV(float_type pMin_out, float_type pMax_out,
		float_type pMin_in, float_type pMax_in);
	// ---
	bool next_point();
	void next_point_fast();
	float_type pm_value(int_type number, float_type t_);
	float_type pm_der_value(int_type number, float_type t_);
	// ---
	std::tuple<float_type, bool> nonlin_solver(float_type t_, float_type Const, int_type index);
	float_type newton_method(float_type a, float_type b, float_type Const, int_type index);
};

Study_DE::Study_DE(int_type num_pists, size_t max_iter, size_t add_iter)
{
	num_of_pists = num_pists;
	max_iterations = max_iter, add_iterations = add_iter;
	R = 0, p = 0, A = 0, FI = 0, a = 0, b = 0, mu = 0, dx_start = 0, t_start = 0;
	lambda = std::vector<float_type>(num_of_pists);
	tau = std::vector<float_type>(num_of_pists + 1);
	eps = tau, fi = tau, gamma = tau, xder = tau;

	h_one = pow(10, -3), h_in = pow(10, -3), h_out = pow(10, -2);
	min_speed = pow(10, -4);
}

/*

	Общие вспомогательные функции

*/

void Study_DE::reset_parametrs()
{
	a = 0, b = 0;
	for (int_type i = 0; i < num_of_pists; ++i)
	{
		a += lambda[i] * gamma[i] * sin(fi[i]);
		b += lambda[i] * gamma[i] * cos(fi[i]);
	}
	a *= mu, b *= mu;
	A = sqrt(pow(a, 2) + pow(b, 2));
	FI = atan(a / b);
	gamma[num_of_pists] = gamma[0];
	eps[num_of_pists] = eps[0];
	fi[num_of_pists] = fi[0];
}

void Study_DE::start_time_search()
{
	size_t find;
	float_type temp;
	for (float_type tau = 0; tau <= 2 * PI; tau += pow(10, -3))
	{
		std::tie(temp, find) = surface_val_ind(tau);
		if (find == 0)
			std::cout << "yes! " << tau << "\n";
	}
}

float_type Study_DE::max_eigv()
{
	matrix JacobiM(2, vector(2)), NextJ(2, vector(2));

	JacobiM[0][0] = dpm_der_value(1, 0), JacobiM[0][1] = dpm_der_value(2, 0);
	JacobiM[1][0] = dpm_der_value(3, 0), JacobiM[1][1] = dpm_der_value(4, 0);

	for (int_type i = 1; i < num_of_pists; ++i)
	{
		NextJ[0][0] = dpm_der_value(1, i), NextJ[0][1] = dpm_der_value(2, i);
		NextJ[1][0] = dpm_der_value(3, i), NextJ[1][1] = dpm_der_value(4, i);
		JacobiM = NextJ * JacobiM;
	}

	vector eigv_abs(2);
	float_type trace = JacobiM[0][0] + JacobiM[1][1], det;
	det = JacobiM[0][0] * JacobiM[1][1] - JacobiM[0][1] * JacobiM[1][0];

	std::complex<float_type> D = pow(trace, 2) - 4 * det;
	eigv_abs[0] = abs((trace + sqrt(D)) / 2.0), eigv_abs[1] = abs((trace - sqrt(D)) / 2.0);

	return (eigv_abs[0] >= eigv_abs[1] ? eigv_abs[0] : eigv_abs[1]);
}

float_type Study_DE::dpm_der_value(int_type elNum, int_type ind)
{
	float_type Value = 0, temp1, temp2;
	switch (elNum)
	{
	case 1:
		temp1 = -R * (A * cos(tau[ind + 1] - FI) - p)
			+ (1 + R) * mu * gamma[ind + 1] * cos(tau[ind + 1] - fi[ind + 1]);
		Value = dpm_der_value(3, ind) * temp1 - R;
		break;
	case 2:
		temp1 = -R * (A * cos(tau[ind + 1] - FI) - p)
			+ (1 + R) * mu * gamma[ind + 1] * cos(tau[ind + 1] - fi[ind + 1]);
		Value = dpm_der_value(4, ind) * temp1 - R * (p - A * cos(tau[ind] - FI));
		break;
	case 3:
		temp1 = mu * gamma[ind + 1] * sin(tau[ind + 1] - fi[ind + 1]) - A * sin(tau[ind + 1] - FI)
			+ p * (tau[ind + 1] - tau[ind]) - xder[ind] + A * sin(tau[ind] - FI);
		Value = (tau[ind + 1] - tau[ind]) / temp1;
		break;
	case 4:
		temp1 = mu * gamma[ind + 1] * sin(tau[ind + 1] - fi[ind + 1]) - A * sin(tau[ind + 1] - FI)
			+ p * (tau[ind + 1] - tau[ind]) - xder[ind] + A * sin(tau[ind] - FI);
		temp2 = mu * gamma[ind] * sin(tau[ind] - fi[ind]) + p * (tau[ind + 1] - tau[ind]) - xder[ind]
			- A * cos(tau[ind] - FI) * (tau[ind + 1] - tau[ind]);
		Value = temp2 / temp1;
		break;
	}

	return Value;
}

std::tuple<float_type, int_type> Study_DE::surface_val_ind(float_type t_)
{
	int_type max_index = 0;
	float_type max_value = eps[0] - mu * gamma[0] * cos(t_ - fi[0]), cur_value;
	for (int_type cur_index = 1; cur_index < num_of_pists; ++cur_index)
	{
		cur_value = eps[cur_index] - mu * gamma[cur_index] * cos(t_ - fi[cur_index]);
		if (cur_value > max_value)
		{
			max_index = cur_index;
			max_value = cur_value;
		}
	}

	return std::make_tuple(max_value, max_index);
}

float_type Study_DE::surface_der_val(int_type surf_index, float_type t_)
{
	return mu * gamma[surf_index] * sin(t_ - fi[surf_index]);
}

/*

	Используем многократное решение дифференциального уравнения

*/

// Осцилограмма
void Study_DE::oscillogram_MSV(float_type step, float_type interval)
{
	int_type surf_index;
	float_type t, tmax = t_start + interval, surf_value;
	std::ofstream surf("data\\surface.txt"), traj("data\\trajectory.txt");

	// Обновляем формулы после смены параметров
	reset_parametrs();

	// Заполнили файл для графика поверхности
	for (t = t_start; t < tmax; t += step)
	{
		std::tie(surf_value, surf_index) = surface_val_ind(t);
		surf << t << " " << surf_value << " " << surf_index << "\n";
	}

	// Задаём начальные условия и отдаём в решатель ДУ
	t = t_start;
	std::tie(surf_value, surf_index) = surface_val_ind(t);
	vector var = { surf_value, dx_start };
	// ---
	DE_SOLVE solve(ode_function, var, t, pow(10, -3), pow(10, -12), { A, FI, p });

	for (size_t i = 0; (i < 10000) && (t < tmax); ++i)
	{
		traj << t << " " << var[0] << "\n";
		std::tie(var, t) = solve.next_point();
		std::tie(surf_value, surf_index) = surface_val_ind(t);

		if (var[0] < surf_value)
		{
			solve.step_back();
			std::tie(var, t) = solve.get_data();
			std::tie(surf_value, surf_index) = surface_val_ind(t);

			while (fabs(var[0] - surf_value) > pow(10, -12))
			{
				std::tie(var, t) = solve.next_point(false);
				std::tie(surf_value, surf_index) = surface_val_ind(t);
				if (var[0] < surf_value)
				{
					solve.step_back();
					std::tie(var, t) = solve.get_data();
					std::tie(surf_value, surf_index) = surface_val_ind(t);
				}
			}
			var[1] = -R * var[1] + (1 + R) * surface_der_val(surf_index, t);

			// Проверка на скорость
			if (fabs(var[1] - surface_der_val(surf_index, t)) < min_speed) {
				break;
			}
			solve.reset_values(var);
			solve.reset_step(pow(10, -3));
		}
	}
	surf.close();
	traj.close();

	// Вывод
	std::string str;
	str = settings_term;
	str += settings_palette;
	str += "set grid\n";
	str += "set title 'oscillogramm'\n";
	str += "set output 'pictures\\oscillogram.png'\n";
	str += "plot 'data\\surface.txt' u 1:2:3 with lines lt 1 lw 1 lc palette title 'surface', ";
	str += "'data\\trajectory.txt' with lines lc 'black' title 'x'\n";
	Plot(str);
}

// Бифуркационная диаграмма
void Study_DE::bifurcation_diagramm_MSV(float_type pMin, float_type pMax)
{
	bool Correct;
	vector var(2);
	int_type surf_index;
	float_type t, surf_value;
	std::ofstream points("data\\BD_MSV.txt");

	DE_SOLVE solve;
	for (P_ONE = pMin; P_ONE <= pMax; P_ONE += h_one)
	{
		// Обновляем формулы после смены параметров
		reset_parametrs();

		// Задаём начальные условия и отдаём в решатель ДУ
		t = t_start;
		std::tie(surf_value, surf_index) = surface_val_ind(t);
		var = { surf_value, dx_start };
		// ---
		solve = DE_SOLVE(ode_function, var, t, pow(10, -3), pow(10, -12), { A, FI, p });

		Correct = true;
		for (size_t i = 0; (i < max_iterations) && Correct; ++i) {
			std::tie(surf_index, Correct) = next_point_de(solve);
		}

		if (Correct)
		{
			points << std::fixed << std::setprecision(10);
			for (size_t i = 0; i < add_iterations; ++i)
			{
				std::tie(surf_index, Correct) = next_point_de(solve);
				std::tie(var, t) = solve.get_data();
				points << P_ONE << " " << var[1] << " " << surf_index << "\n";
			}
		}

	}
	points.close();

	// Вывод
	std::string str;
	str = settings_term;
	str += settings_palette;
	str += "set grid\n";
	str += "set xlabel '" + SP_ONE + "'\n";
	str += "set ylabel 'u'\n";
	str += "set title 'bifurcation diagramm'\n";
	str += "set output 'pictures\\BD_MSV.png'\n";
	str += "plot 'data\\BD_MSV.txt' u 1:2:3 w points pt 7 ps 0.1 lc palette notitle\n";
	Plot(str);
}

//Карта динамических режимов
void Study_DE::dynamic_mode_map_MSV(float_type pMin_out, float_type pMax_out,
	float_type pMin_in, float_type pMax_in)
{
	bool Correct;
	vector var(2);
	float_type t, surf_value;
	std::vector<int_type> count(num_of_pists);
	int_type surf_index, num_cycles = 16, cycle_size;
	matrix dk(num_of_pists, vector(num_cycles, 100.0));

	DE_SOLVE solve;

	// ---

	std::vector<std::ofstream> points(num_of_pists);
	for (int_type i = 0; i < num_of_pists; i++)
		points[i].open("data\\DMM_MSV" + std::to_string(i) + ".txt", std::ios_base::out);

	// ---

	for (P_OUT = pMin_out; P_OUT <= pMax_out; P_OUT += h_out)
	{
		for (P_IN = pMin_in; P_IN <= pMax_in; P_IN += h_in)
		{
			// Обновляем формулы после смены параметров
			reset_parametrs();

			// Задаём начальные условия и отдаём в решатель ДУ
			t = t_start;
			std::tie(surf_value, surf_index) = surface_val_ind(t);
			var = { surf_value, dx_start };
			// ---
			solve = DE_SOLVE(ode_function, var, t, pow(10, -3), pow(10, -12), { A, FI, p });

			Correct = true;
			for (size_t i = 0; (i < max_iterations) && Correct; ++i) {
				std::tie(surf_index, Correct) = next_point_de(solve);
			}

			if (Correct)
			{
				// Берём значения скоростей для всех поршней за несколько итераций
				for (int_type i = 0; i < num_of_pists * num_cycles; i++)
				{
					std::tie(surf_index, Correct) = next_point_de(solve);
					std::tie(var, t) = solve.get_data();

					if (count[surf_index] < num_cycles)
						dk[surf_index][count[surf_index]++] = var[1];
				}

				for (int_type row = 0; row < num_of_pists; ++row)
				{
					cycle_size = 0;
					for (int_type cell = num_cycles - 2; cell >= 0; --cell)
					{
						if (fabs(dk[row][num_cycles - 1] - dk[row][cell]) < pow(10, -6))
						{
							cycle_size = num_cycles - 1 - cell;
							break;
						}
					}
					// Записываем в отдельные файлы размеры циклов для каждого поршня
					points[row] << P_OUT << " " << P_IN << " " << cycle_size << "\n";
				}
			}
			else {
				for (int_type i = 0; i < num_of_pists; ++i) {
					points[i] << P_OUT << " " << P_IN << " " << -1 << "\n";
				}
			}
		}

		// Делаем разделение перед новым слоем
		for (int_type i = 0; i < num_of_pists; ++i) {
			points[i] << "\n";
		}
	}

	// Вывод
	std::string str, name;
	str = settings_term;
	str += settings_palette;
	str += "set xlabel '" + SP_OUT + "'\n";
	str += "set ylabel '" + SP_IN + "'\n";

	for (int_type i = 0; i < num_of_pists; i++)
	{
		points[i].close();
		name = str + "set title 'dynamic mode map" + std::to_string(i + 1) + "'\n";
		name += "set output 'pictures\\DMM_MSV" + std::to_string(i + 1) + ".png'\n";
		name += "plot 'data\\DMM_MSV" + std::to_string(i) + ".txt' u 1:2:3 with image notitle\n";
		Plot(name);
	}
}

/*

	----------------------------------

*/

std::tuple<int_type, bool> Study_DE::next_point_de(DE_SOLVE& solve)
{
	vector var(2);
	int_type surf_index;
	float_type t, surf_value;

	do
	{
		std::tie(var, t) = solve.next_point();
		std::tie(surf_value, surf_index) = surface_val_ind(t);
	} while (var[0] > surf_value);

	if (var[0] < surf_value)
	{
		solve.step_back();
		std::tie(var, t) = solve.get_data();
		std::tie(surf_value, surf_index) = surface_val_ind(t);
		while (fabs(var[0] - surf_value) > pow(10, -12))
		{
			std::tie(var, t) = solve.next_point(false);
			std::tie(surf_value, surf_index) = surface_val_ind(t);
			if (var[0] < surf_value)
			{
				solve.step_back();
				std::tie(var, t) = solve.get_data();
				std::tie(surf_value, surf_index) = surface_val_ind(t);
			}
		}
	}
	var[1] = -R * var[1] + (1 + R) * surface_der_val(surf_index, t); // Переводим в послеударную скорость

	// Проверка на скорость
	if (fabs(var[1] - surface_der_val(surf_index, t)) < min_speed)
		return std::make_tuple(0, false);

	solve.reset_values(var);
	solve.reset_step(pow(10, -3));
	return std::make_tuple(surf_index, true);
}

/*

	Используем точечное отображение

*/

// Точечное отображение
void Study_DE::point_map_DMV()
{
	bool Correct = true;
	std::ofstream points("data\\PM_DMV.txt");

	// Обновляем формулы после смены параметров
	reset_parametrs();

	tau[num_of_pists] = t_start;
	xder[num_of_pists] = dx_start;

	for (size_t i = 0; (i < add_iterations) && Correct; ++i)
	{
		points << xder[0] << " " << xder[num_of_pists] << "\n";
		points << xder[num_of_pists] << " " << xder[num_of_pists] << "\n";
		Correct = next_point();
	}
	points.close();

	// Вывод
	std::string str;
	str = settings_term;
	str += settings_palette;
	str += "set grid\n";
	str += "set title 'point map'\n";
	str += "set output 'pictures\\PM_DMV.png'\n";
	str += "plot 'data\\PM_DMV.txt' with lines lc 'black' notitle\n";
	Plot(str);
}

// Бифуркационная диаграмма
void Study_DE::bifurcation_diagramm_DMV(float_type pMin, float_type pMax)
{
	bool Correct = true;
	std::ofstream points("data\\BD_DMV.txt");

	for (P_ONE = pMin; P_ONE <= pMax; P_ONE += h_one)
	{
		// Обновляем формулы после смены параметров
		reset_parametrs();

		tau[num_of_pists] = t_start;
		xder[num_of_pists] = dx_start;
		
		Correct = true;
		for (size_t i = 0; (i < max_iterations) && Correct; ++i) {
			Correct = next_point();
		}

		if (Correct)
		{
			points << std::fixed << std::setprecision(10);
			for (size_t i = 0; i < add_iterations; ++i)
			{
				for (int_type j = 1; j <= num_of_pists; ++j) {
					points << P_ONE << " " << xder[j] << " " << j % num_of_pists << "\n";
				}
				next_point_fast();
			}
		}
	}
	points.close();

	// Вывод
	std::string str;
	str = settings_term;
	str += settings_palette;
	str += "set grid\n";
	str += "set xlabel '" + SP_ONE + "'\n";
	str += "set ylabel 'u'\n";
	str += "set title 'bifurcation diagramm'\n";
	str += "set output 'pictures\\BD_DMV.png'\n";
	str += "plot 'data\\BD_DMV.txt' u 1:2:3 w points pt 7 ps 0.1 lc palette notitle\n";
	Plot(str);
}

// Диаграмма Ляпунова
void Study_DE::lyapunov_diagramm_DMV(float_type pMin, float_type pMax)
{
	bool Correct = true;
	float_type df = 1, sum = 0;
	std::ofstream points("data\\LD_DMV.txt");

	for (P_ONE = pMin; P_ONE <= pMax; P_ONE += h_one)
	{
		// Обновляем формулы после смены параметров
		reset_parametrs();

		tau[num_of_pists] = t_start;
		xder[num_of_pists] = dx_start;

		sum = 0;
		Correct = true;
		for (size_t i = 0; (i < max_iterations) && Correct; ++i) {
			Correct = next_point();
		}

		if (Correct)
		{
			points << std::fixed << std::setprecision(10);
			for (size_t i = 0; i < add_iterations; ++i)
			{
				df = max_eigv();
				sum += log(df);

				next_point_fast();
			}
			points << P_ONE << " " << sum / (1.0 * add_iterations) << "\n";
		}
	}
	points.close();

	// Вывод
	std::string str;
	str = settings_term;
	str += "set grid\n";
	str += "set xlabel '" + SP_ONE + "'\n";
	str += "set ylabel '{/Symbol l}'\n";
	str += "set title 'lyapunov diagramm'\n";
	str += "set output 'pictures\\LD_DMV.png'\n";
	str += "plot 'data\\LD_DMV.txt' u 1:2 w points pt 7 ps 0.1 lc 'black' notitle\n";
	Plot(str);
}

// Карта динамических режимов
void Study_DE::dynamic_mode_map_DMV(float_type pMin_out, float_type pMax_out, 
	float_type pMin_in, float_type pMax_in)
{
	bool Correct;
	int_type num_cycles = 16, cycle_size;
	matrix dk(num_of_pists, vector(num_cycles));

	// ---

	std::vector<std::ofstream> points(num_of_pists);
	for (int_type i = 0; i < num_of_pists; i++)
		points[i].open("data\\DMM_DMV" + std::to_string(i) + ".txt", std::ios_base::out);

	// ---

	for (P_OUT = pMin_out; P_OUT <= pMax_out; P_OUT += h_out)
	{
		for (P_IN = pMin_in; P_IN <= pMax_in; P_IN += h_in)
		{
			// Обновляем формулы после смены параметров
			reset_parametrs();

			tau[num_of_pists] = t_start;
			xder[num_of_pists] = dx_start;

			Correct = true;
			for (size_t i = 0; (i < max_iterations) && Correct; ++i) {
				Correct = next_point();
			}

			if (Correct)
			{
				// Берём значения скоростей для всех поршней за несколько итераций
				for (int_type cell = 0; cell < num_cycles; ++cell)
				{
					for (int_type row = 1; row <= num_of_pists; ++row) {
						dk[row % num_of_pists][cell] = xder[row];
					}

					next_point_fast();
				}

				// Ищем циклы в полученных последовательностях
				for (int_type row = 0; row < num_of_pists; ++row) 
				{
					cycle_size = 0;
					for (int_type cell = num_cycles - 2; cell >= 0; --cell)
					{
						if (fabs(dk[row][num_cycles - 1] - dk[row][cell]) < pow(10, -6))
						{
							cycle_size = num_cycles - 1 - cell;
							break;
						}
					}
					// Записываем в отдельные файлы размеры циклов для каждого поршня
					points[row] << P_OUT << " " << P_IN << " " << cycle_size << "\n";
				}
			}
			else {
				for (int_type i = 0; i < num_of_pists; i++) {
					points[i] << P_OUT << " " << P_IN << " " << -1 << "\n";
				}
			}
		}

		// Делаем разделение перед новым слоем
		for (int_type i = 0; i < num_of_pists; i++) {
			points[i] << "\n";
		}
	}

	// Вывод
	std::string str, name;
	str = settings_term;
	str += settings_palette;
	str += "set xlabel '" + SP_OUT + "'\n";
	str += "set ylabel '" + SP_IN + "'\n";

	for (int_type i = 0; i < num_of_pists; i++)
	{
		points[i].close();
		name = str + "set title 'dynamic mode map" + std::to_string(i + 1) + "'\n";
		name += "set output 'pictures\\DMM_DMV" + std::to_string(i + 1) + ".png'\n";
		name += "plot 'data\\DMM_DMV" + std::to_string(i) + ".txt' u 1:2:3 with image notitle\n";
		Plot(name);
	}
}

/*

	----------------------------------

*/

bool Study_DE::next_point()
{
	bool Correct = true;
	int_type surf_index;
	float_type t_cur = 0, surf_value;

	tau[0] = tau[num_of_pists];
	xder[0] = xder[num_of_pists];

	//Проверка на принадлежность t_0 поверхности f_1
	std::tie(surf_value, surf_index) = surface_val_ind(tau[0]);
	if (surf_index != 0)
		return false;

	std::tie(t_cur, Correct) = nonlin_solver(tau[0] + pow(10, -2), xder[0], 0);
	for (int_type i = 1; (i < num_of_pists) && Correct; ++i)
	{
		tau[i] = t_cur, xder[i] = pm_value(i, tau[i]);

		//Проверка на скорость и на принадлежность t_i поверхности f_i
		std::tie(surf_value, surf_index) = surface_val_ind(tau[i]);
		if ((fabs(xder[i] - surface_der_val(i, t_cur)) < min_speed) || (surf_index != i))
			return false;

		std::tie(t_cur, Correct) = nonlin_solver(tau[i] + pow(10, -2), xder[i], -i);
	}

	if (Correct)
	{
		xder[num_of_pists] = pm_value(num_of_pists, t_cur);
		std::tie(surf_value, surf_index) = surface_val_ind(t_cur);
		//Проверка на скорость и на принадлежность t_N поверхности f_0
		if ((fabs(xder[num_of_pists] - surface_der_val(0, t_cur)) < min_speed) || (surf_index != 0))
			return false;
	}

	tau[num_of_pists] = t_cur;
	return Correct;
}

void Study_DE::next_point_fast()
{
	bool Correct = true;
	float_type t_cur = 0;

	tau[0] = tau[num_of_pists];
	xder[0] = xder[num_of_pists];

	std::tie(t_cur, Correct) = nonlin_solver(tau[0] + pow(10, -2), xder[0], 0);

	for (int_type i = 1; (i < num_of_pists) && Correct; ++i)
	{
		tau[i] = t_cur, xder[i] = pm_value(i, tau[i]);
		std::tie(t_cur, Correct) = nonlin_solver(tau[i] + pow(10, -2), xder[i], -i);
	}

	tau[num_of_pists] = t_cur;
	xder[num_of_pists] = pm_value(num_of_pists, t_cur);
}

// Функции точечного отображения
float_type Study_DE::pm_value(int_type number, float_type t_)
{
	int_type index;
	float_type value = 0;

	if (number <= 0)
	{
		index = abs(number);
		value = (mu * gamma[index] * cos(tau[index] - fi[index])
			- mu * gamma[index + 1] * cos(t_ - fi[index + 1])
			+ A * (cos(t_ - FI) - cos(tau[index] - FI))
			+ eps[index + 1] - eps[index]) / (t_ - tau[index])
			+ p * (t_ - tau[index]) / 2.0 + A * sin(tau[index] - FI);
	}
	else
	{
		index = number - 1;
		value = -R * (A * (sin(t_ - FI) - sin(tau[index] - FI)) - p * (t_ - tau[index])
			+ pm_value(-index, t_)) + (1 + R) * mu * gamma[index + 1] * sin(t_ - fi[index + 1]);
	}

	return value;
}

// Производные точечного отображения
float_type Study_DE::pm_der_value(int_type number, float_type t_)
{
	int_type index;
	float_type value = 0;

	if (number <= 0)
	{
		index = abs(number);
		value = -(mu * gamma[index] * cos(tau[index] - fi[index])
			- mu * gamma[index + 1] * cos(t_ - fi[index + 1])
			+ A * (cos(t_ - FI) - cos(tau[index] - FI))
			+ eps[index + 1] - eps[index]) / pow((t_ - tau[index]), 2)
			+ (mu * gamma[index + 1] * sin(t_ - fi[index + 1])
				- A * sin(t_ - FI)) / (t_ - tau[index]) + p / 2.0;
	}
	else
	{
		index = number - 1;
		value = -R * (A * cos(t_ - FI) - p + pm_der_value(-index, t_))
			+ (1 + R) * mu * gamma[index + 1] * cos(t_ - fi[index + 1]);
	}

	return value;
}

// Решение нелинейных уравнений
std::tuple<float_type, bool> Study_DE::nonlin_solver(float_type t_, float_type Const, int_type index)
{
	int_type SupInd = abs(index) + 1;
	float_type root = newton_method(t_, t_ + 1.0, Const, index);

	if (pm_value(SupInd, root) - surface_der_val(SupInd, root) > 0)
		return std::make_tuple(root, true);

	return std::make_tuple(t_, false);
}

// Метод Ньютона
float_type Study_DE::newton_method(float_type left_border,
	float_type right_border, float_type Const, int_type index)
{
	size_t max_iter = 1000, count = 0;
	float_type x, f = 1, df, eps = pow(10, -10);
	x = (left_border + right_border) / 2.0;

	f = pm_value(index, x) - Const;
	df = pm_der_value(index, x);
	while (fabs(f) > eps && count < max_iter)
	{
		x = x - f / df;
		f = pm_value(index, x) - Const;
		df = pm_der_value(index, x);
		++count;
	}
	return x;
}

#endif

