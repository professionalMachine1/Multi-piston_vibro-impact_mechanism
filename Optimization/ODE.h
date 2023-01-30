#pragma once

#ifndef ODE_H
#define ODE_H

#include "Directives.h"

class DE_SOLVE
{
private:
	vector parametrs;
	vector x, X, xprev, Xprev;
	float_type eps, tprev, t, step;
	vector(*function) (const vector&, float_type, const vector&);
public:
	DE_SOLVE() 
	{ 
		function = NULL;
		eps = 1, t = tprev = 0, step = 0; 
	}
	DE_SOLVE(vector(*function_) (const vector&, float_type, const vector&), const vector& x0,
		float_type t0, float_type step_, float_type eps_, const vector& parametrs_);
	std::tuple<vector, float_type> next_point(bool step_change_access = true);
	vector RK4(const vector& v_, float_type t_, float_type h_);
	std::tuple<vector, float_type> get_data() { return std::make_tuple(x, t); }
	void reset_values(const vector& x_new) { x = x_new; }
	void reset_step(float_type new_step) { step = new_step; }
	void reset_parametrs(const vector& parametrs_) { parametrs = parametrs_; }
	void step_back(bool step_change_access = true);
};

#endif