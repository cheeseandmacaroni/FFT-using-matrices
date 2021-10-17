#ifndef MY_MATRIX_HPP
#define MY_MATRIX_HPP
#include <omp.h>

template <class T>
class m_square_matrix
{
public:
	m_square_matrix(int n) :size(n)
	{
		arr = new T[n * n];
	}
	m_square_matrix(const m_square_matrix &obj):size(0), arr(nullptr)
	{
		*this = obj;
	}
	m_square_matrix(int n, T value):size(n)
	{
		arr = new T[n * n];
		set_all_elements(value);
	}
	~m_square_matrix()
	{
		delete[] arr;
	}
	m_square_matrix& operator=(const m_square_matrix &obj)
	{
		if (this->size != obj.size && this->arr != obj.arr)
		{
			size = obj.size;
			if (arr != nullptr)
				delete[] arr;
			arr = new T[size * size];
			for (int i = 0; i < size*size; ++i)
			{
				arr[i] = obj.arr[i];
			}
		}
		return *this;
	}
	inline T& get_elem(int i, int j)
	{
		return arr[i * size + j];
	}
	int get_size() const { return size; }
	inline void set_all_elements(T value)
	{
		#pragma omp parallel for num_threads(2)
		for (int i = 0; i < size * size; ++i)
		{
			arr[i] = value;
		}
	}
private:
	T* arr;
	int size;
};

#endif
