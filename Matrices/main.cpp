#include <iostream>
#include "Matrix.h"
#include "ImprovedMatrix.h"
#include "ComplexNumber.h"
#include "SIMDMatrix.h"
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
void printMatrix(SIMDMatrix& mat) {
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
		Timer t2;
		float** array = new float*[m];
		float** arr = new float*[m];
		for (register int i = 0; i < m; i++) {
			array[i] = new float[n];
			arr[i] = new float[n];
		}

		float temp;
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++) {
				//std::cout << "Insert the element (" << i << ',' << j << "): ";
				//std::cin >> temp;
				temp = rand() % 20;
				array[i][j] = temp;
				temp = rand() % 20;
				arr[i][j] = temp;
			}
		}

		SIMDMatrix mat(array, m, n);
		SIMDMatrix mat2(arr, m, n);

		ImprovedMatrix<float> mat1(array, m, n);
		ImprovedMatrix<float> mat3(arr, m, n);
		t.start();
		mat += mat2;
		t.stop();
		t2.start();
		mat1 += mat3;
		t2.stop();
		//printMatrix(mat);
		//printMatrix(mat1);
		std::cout << "SIMD took: " << t.elapsedMilliseconds() << '\n';
		std::cout << "NON-SIMD took: " << t2.elapsedMilliseconds() << '\n';

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