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

void applyKernel(SIMDMatrix &m, SIMDMatrix& ker);

int main() {
	/*
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
				//temp = rand() % 20;
				arr[i][j] = temp;
			}
		}

		SIMDMatrix mat(array, m, n);

		printMatrix(mat);
		std::cout << "\n\n\n";

		t2.start();
		double det2 = mat.det_2();
		t2.stop();

		std::cout << " Det2: " << det2 << '\n';
		std::cout << "time for Det2: " << t2.elapsedMilliseconds() << '\n';

		//clearing up mess
		for (int i = 0; i < m; i++)
		{
			delete array[i];
		}
		delete[] array;
	}
	*/
	float k[] = { 2,3,-2,4,5,-4,-4,-4,2 };
	SIMDMatrix kernel(k, 3, 3);

	float m[] = { 5,4,25,35,20,60,45,40,30,10,11,60,80,15,5,80,22,20,10,1,200,140,50,0,8 };
	SIMDMatrix mat(m, 5, 5);

	printMatrix(mat);

	applyKernel(mat, kernel);

	printMatrix(mat);

	system("pause");

	return 0;
}

void applyKernel(SIMDMatrix &m, SIMDMatrix& ker) {
	int size = ker.getCols();
	int costant = int(size / 2);
	float res = 0;
	for (register int i = 0; i < m.getRows(); i++) {
		for (register int j = 0; j < m.getCols(); j++) {
			for (register int k = -costant; k <= costant; k++) {
				for (register int l = -costant; l <= costant; l++) {
					if (i + k < 0 || i + k >= m.getRows())continue;
					if (j + l < 0 || j + l >= m.getCols())continue;
					res += (m[i + k][j + l] * ker[costant + k][costant + l]);
				}
			}
			m[i][j] = res;
			res = 0;
		}
	}
}