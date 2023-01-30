#pragma once

#ifndef DIRECTIVES_H
#define DIRECTIVES_H

#include <algorithm>
#include <iostream>
#include <complex>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>

typedef int32_t int_type;
typedef double float_type;
typedef std::vector<float_type> vector;
typedef std::vector<std::vector<float_type>> matrix;

const float_type PI = 2 * asin(1);

matrix operator*(const matrix& left, const matrix& right);
vector operator*(const vector& vec, float_type num);
vector operator*(float_type num, const vector& vec);
vector operator/(const vector& vec, float_type num);
vector operator+(const vector& vec1, const vector& vec2);
vector operator-(const vector& first_vec, const vector& second_vec);
vector fabs(const vector& vec);

#endif