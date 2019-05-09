#include <iostream>
#include "Matrix.h"
#include "ImprovedMatrix.h"
#include "ComplexNumber.h"
#include <random>
#include <ctime>
#include <chrono>

class Timer
{
public:
	void start()
	{
		m_StartTime = std::chrono::system_clock::now();
		m_bRunning = true;
	}

	void stop()
	{
		m_EndTime = std::chrono::system_clock::now();
		m_bRunning = false;
	}

	double elapsedMilliseconds()
	{
		std::chrono::time_point<std::chrono::system_clock> endTime;

		if (m_bRunning)
		{
			endTime = std::chrono::system_clock::now();
		}
		else
		{
			endTime = m_EndTime;
		}

		return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
	}

	double elapsedSeconds()
	{
		return elapsedMilliseconds() / 1000.0;
	}

private:
	std::chrono::time_point<std::chrono::system_clock> m_StartTime;
	std::chrono::time_point<std::chrono::system_clock> m_EndTime;
	bool                                               m_bRunning = false;
};

template <class T>
void printMatrix(Matrix<T>& mat) {
	for (int i = 0; i < mat.getRows(); i++)
	{
		for (int j = 0; j < mat.getCols(); j++) {
			std::cout << mat[i][j] << ",";
		}
		std::cout << '\n';
	}
	std::cout << std::endl << std::endl;
}
template <class T>
void printMatrix(ImprovedMatrix<T>& mat) {
	for (int i = 0; i < mat.getRows(); i++)
	{
		for (int j = 0; j < mat.getCols(); j++) {
			std::cout << mat[i][j] << ",";
		}
		std::cout << '\n';
	}
	std::cout << std::endl << std::endl;
}
int main() {
	srand(time(NULL));
	while (true) {
		std::cout << "Insert m and n: ";
		int m, n;
		std::cin >> m >> n;
		Timer t;
		double** array = new double*[m];
		for (register int i = 0; i < m; i++) {
			array[i] = new double[n];
		}

		double temp;
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++) {
				std::cout << "Insert the element (" << i << ',' << j << "): ";
				std::cin >> temp;
				//temp = rand() % 600;
				array[i][j] = temp;
			}
		}
		double* temp_arr = new double[m*n];
		for (register int i = 0; i < m*n; i++) {
			temp_arr[i] = 1;
		}
		ImprovedMatrix<double> b(array, m, n);
		ImprovedMatrix<double> temp_mat(temp_arr, m, n);
		printMatrix(b);
		b.GOD(temp_mat);

		//clearing up mess
		for (int i = 0; i < m; i++)
		{
			delete array[i];
		}
		delete[] array;
	}

	system("pause");

	return 0;
}