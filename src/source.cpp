#include <complex>
#include <iostream>
#include <vector>

using namespace std;

#define el_type double
#define Pi 3.14159265359
#define m_matrix vector<vector<complex<el_type>>> 
#define w(N,k) (complex<el_type>(cos(2*Pi*k/N), -sin(2*Pi*k/N)))
#define dft_w(N,n,k) (complex<el_type>(cos(2*Pi*k*n/N), -sin(2*Pi*k*n/N)))

vector<complex<el_type>> operator*(m_matrix &a, vector<complex<el_type>> &b)
{
	int n = a.size();
	vector<complex<el_type>> res(n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			res[i] += a[i][j] * b[j];
	return res;
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
	vector<complex<el_type>> pre_calculated_wnk(n/2);
	for (int i = 0; i < n / 2; ++i)
	{
		pre_calculated_wnk[i] = w(n, i);
	}
	int iterations = my_log_2(n);
	for (int i = 1; i <= iterations; ++i)
	{
		int sub_sequence_size = pow(2, i);
		m_matrix tmp(n, vector<complex<el_type>>(n, complex<el_type>(0, 0)));
		for (int j = 0; j < sub_sequence_size / 2; ++j)
		{
			for (int t = 0; t < n / sub_sequence_size; ++t)
			{
				tmp[j + t * sub_sequence_size][j + t * sub_sequence_size] = 1;
				tmp[j + t * sub_sequence_size][j + sub_sequence_size / 2 + t * sub_sequence_size] = pre_calculated_wnk[j*n/sub_sequence_size];
				tmp[j + sub_sequence_size / 2 + t * sub_sequence_size][j  + t * sub_sequence_size] = 1;
				tmp[j + sub_sequence_size / 2 + t * sub_sequence_size][j + sub_sequence_size / 2 + t * sub_sequence_size] = -pre_calculated_wnk[j*n / sub_sequence_size];
			}
		}
		x = tmp*x;
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
	int n;
	cin >> n;
	vector<complex<el_type>> x(n, complex<el_type>());
	vector<complex<el_type>> y(n, complex<el_type>());
	for (int i = 0; i < n; ++i)
	{
		y[i] = x[i] = complex<el_type>(rand() % 1000 - 500, rand() % 1000 - 500);
	}
	FFT(x);

	reverse(++begin(x), end(x));//
	FFT(x);						//inverse FFT using forward FFT
	for (int i = 0; i < n; ++i) //
		x[i] /= n;				//
	cout << check(1.e-10, x, y);
	return 0;
}