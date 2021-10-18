#ifndef MY_MATRIX_HPP
#define MY_MATRIX_HPP
template <class T>
class m_square_matrix
{
public:
	m_square_matrix() : size(0), arr(nullptr) {};
	m_square_matrix(int n) : size(n)
	{
		arr = new T[n * n];
	}
	m_square_matrix(const m_square_matrix &other) :size(0), arr(nullptr)
	{
		*this = other;
	}
	m_square_matrix(int n, T value) :size(n)
	{
		arr = new T[n * n];
		set_all_elements(value);
	}
	~m_square_matrix()
	{
		delete[] arr;
	}
	m_square_matrix& operator=(const m_square_matrix &other)
	{
		if (this->size != other.size && this->arr != other.arr)
		{
			size = other.size;
			if (arr != nullptr)
				delete[] arr;
			arr = new T[size * size];
#pragma omp parallel for num_threads(2)
			for (int i = 0; i < size*size; ++i)
			{
				arr[i] = other.arr[i];
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
