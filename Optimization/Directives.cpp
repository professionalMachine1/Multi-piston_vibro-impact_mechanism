#include "Directives.h"

matrix operator*(const matrix& left, const matrix& right)
{
	size_t size = left.size();
	matrix res_matrix = matrix(size, vector(size));
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			for (size_t k = 0; k < size; ++k) {
				res_matrix[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return res_matrix;
}

vector operator*(const vector& vec, float_type num)
{
	size_t size = vec.size();
	vector product(size);
	for (size_t i = 0; i < size; ++i) {
		product[i] = vec[i] * num;
	}
	return product;
}

vector operator*(float_type num, const vector& vec)
{
	return vec * num;
}

vector operator/(const vector& vec, float_type num)
{
	return vec * (1 / num);
}

vector operator+(const vector& vec1, const vector& vec2)
{
	size_t size = std::min(vec1.size(), vec2.size());
	vector product(size);
	for (size_t i = 0; i < size; ++i) {
		product[i] = vec1[i] + vec2[i];
	}
	return product;
}

vector operator-(const vector& first_vec, const vector& second_vec)
{
	size_t size = std::min(first_vec.size(), second_vec.size());
	vector product(size);
	for (size_t i = 0; i < size; ++i) {
		product[i] = first_vec[i] - second_vec[i];
	}
	return product;
}

vector fabs(const vector& vec)
{
	size_t size = vec.size();
	vector product(size);
	for (size_t i = 0; i < size; ++i) {
		product[i] = fabs(vec[i]);
	}
	return product;
}
