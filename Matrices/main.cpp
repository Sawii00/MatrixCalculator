#include <iostream>
#include "Matrix.h"

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
int main() {
	std::cout << "Insert m and n: ";
	int m, n;
	std::cin >> m >> n;

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
			array[i][j] = temp;
		}
	}

	Matrix<double> a(array, m, n);
	printMatrix(a);

	a.reduceToEchelonForm();
	printMatrix(a);

	system("pause");

	return 0;
}