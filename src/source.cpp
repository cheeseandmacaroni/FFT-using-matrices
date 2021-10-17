#include <complex>
#include <iostream>
#include <vector>
#include "my_matrix.hpp"

using namespace std;

#define el_type double
#define Pi 3.14159265359
#define w(N,k) (complex<el_type>(cos(2*Pi*k/N), -sin(2*Pi*k/N)))
#define dft_w(N,n,k) (complex<el_type>(cos(2*Pi*k*n/N), -sin(2*Pi*k*n/N)))



vector<complex<el_type>> operator*(m_square_matrix<complex<el_type>> &a, vector<complex<el_type>> &b)
{
	int n = a.get_size();
	vector<complex<el_type>> res(n, 0);
	#pragma omp parallel for
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			res[i] += a.get_elem(i,j) * b[j];
	return res;
}

template<typename T>
void matrix_block_multiplication(m_square_matrix<T> &a, m_square_matrix<T> &b, m_square_matrix<T> &res, int sub_size)
{

}

template <typename T>
T m_reverse(T a, int bit_len)
{
	T res = 0;
	for (int i = 0; i < bit_len; ++i)
	{
		bool bit = (a >> i) % 2;
		res |= bit << (bit_len - i - 1);
	}
	return res;
}

int my_log_2(int a)
{
	int res = 0;
	int tmp = 1;
	while (tmp != a)
	{
		//if (tmp > INT32_MAX / 2)
		//	throw exception("not a power of 2");
		tmp *= 2;
		res++;
	}
	return res;
}

vector<complex<el_type>> pre_permutation_algorithm(vector<complex<el_type>> &src)//src.size() must be power of 2
{
	int n = src.size();
	int bit_length = my_log_2(n);
	vector<complex<el_type>> res(n);
	for (int i = 0; i < n; ++i)
	{
		int new_pos = m_reverse(i, bit_length);
		res[new_pos] = src[i];
	}
	return res;
}

void FFT(vector<complex<el_type>> &x)
{
	int n = x.size();
	x = pre_permutation_algorithm(x);
	m_square_matrix<complex<el_type>> tmp(n, complex<el_type>(0, 0));
	vector<complex<el_type>> pre_calculated_wnk(n/2);
	for (int i = 0; i < n / 2; ++i)
	{
		pre_calculated_wnk[i] = w(n, i);
	}
	int iterations = my_log_2(n);
	for (int i = 1; i <= iterations; ++i)
	{
		int sub_sequence_size = pow(2, i);
		for (int j = 0; j < sub_sequence_size / 2; ++j)
		{
			for (int t = 0; t < n / sub_sequence_size; ++t)
			{
				tmp.get_elem(j + t * sub_sequence_size , j + t * sub_sequence_size) = 1;
				tmp.get_elem(j + t * sub_sequence_size , j + sub_sequence_size / 2 + t * sub_sequence_size) = pre_calculated_wnk[j*n/sub_sequence_size];
				tmp.get_elem(j + sub_sequence_size / 2 + t * sub_sequence_size , j  + t * sub_sequence_size) = 1;
				tmp.get_elem(j + sub_sequence_size / 2 + t * sub_sequence_size , j + sub_sequence_size / 2 + t * sub_sequence_size) = -pre_calculated_wnk[j*n / sub_sequence_size];
			}
		}
		x = tmp*x;
		tmp.set_all_elements(complex<el_type>(0, 0));
	}
}

void DFT(vector<complex<el_type>> &x)
{
	int N = x.size();
	vector<complex<el_type>> res(N);
	for (int k = 0; k < N; ++k)
	{
		for (int n = 0; n < N; ++n)
		{
			res[k] += dft_w(N, n, k) * x[n];
		}
	}
	x = res;
}

bool check(el_type epsilon, vector<complex<el_type>> &x, vector<complex<el_type>> &y)
{
	for (int i = 0; i < x.size(); ++i)
	{
		if (abs(x[i].real() - y[i].real()) > epsilon || (x[i].imag() - y[i].imag()) > epsilon)
			return false;
	}
	return true;
}

int main()
{
	int n = 8192;
	double start;
	double finish;
	vector<complex<el_type>> x(n, complex<el_type>());
	vector<complex<el_type>> y(n, complex<el_type>());
	for (int i = 0; i < n; ++i)
	{
		y[i] = x[i] = complex<el_type>(rand() % 20 / 10. - 1, rand() % 20 / 10. - 1);
	}
	start = omp_get_wtime();
	FFT(x);
	finish = omp_get_wtime();
	/*DFT(y);
	for (int i = 0; i < n; ++i)
	{
		cout << x[i] << "!!!!!" << y[i] << endl;
	}*/
	cout << finish - start << "forward fft" << endl; 
	start = omp_get_wtime();
	reverse(++begin(x), end(x));//
	FFT(x);						//inverse FFT using forward FFT
	for (int i = 0; i < n; ++i) //
		x[i] /= n;				//
	finish = omp_get_wtime();
	cout << finish - start << "inverse fft" << endl;


	return 0;
}