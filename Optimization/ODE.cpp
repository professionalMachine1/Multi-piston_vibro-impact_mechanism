#include "ODE.h"

DE_SOLVE::DE_SOLVE(vector(*function_) (const vector&, float_type, const vector&), const vector& x0,
	float_type t0, float_type step_, float_type eps_, const vector& parametrs_)
{
	t = t0;
	step = step_, eps = eps_;
	// ---
	x = x0, xprev = x0;
	parametrs = parametrs_;
	X = vector(x.size()), Xprev = vector(x.size());
	// ---
	function = function_;
}

vector DE_SOLVE::RK4(const vector& v_, float_type t_, float_type h_)
{
	size_t k = x.size();
	vector k1(k), k2(k), k3(k), k4(k), next(k);
	k1 = function(v_, t_, parametrs);
	k2 = v_ + k1 * (h_ / 2), k2 = function(k2, t_ + h_ / 2, parametrs);
	k3 = v_ + k2 * (h_ / 2), k3 = function(k3, t_ + h_ / 2, parametrs);
	k4 = v_ + k3 * h_, k4 = function(k4, t_ + h_, parametrs);

	next = v_ + (k1 + 2 * k2 + 2 * k3 + k4) * (h_ / 6);
	return next;
}

std::tuple<vector, float_type> DE_SOLVE::next_point(bool step_change_access)
{
	float_type MaxS;
	vector S(x.size());

	xprev = x, Xprev = X;
	// X - Vn+1 вычисление с шагом step/2
	X = RK4(x, t, step / 2.0), X = RK4(X, t + step / 2.0, step / 2.0);
	// x - Vn+1 вычисленное с шагом step
	x = RK4(x, t, step);

	// Ѕерем модуль погрешности
	S = fabs((X - x) / 15.0);
	// Ќаходим максимальную погрешность
	MaxS = *std::max_element(S.begin(), S.end());

	tprev = t, t = t + step;
	if ((step_change_access) && (MaxS < eps / 32.0))
		step = step * 2.0;
	else
		if (MaxS > eps)
		{
			step_back();
			next_point();
		}

	return std::make_tuple(x, t);
}

void DE_SOLVE::step_back(bool step_change_access)
{
	t = tprev;
	x = xprev, X = Xprev;
	if (step_change_access)
		step = step / 2.0;
}