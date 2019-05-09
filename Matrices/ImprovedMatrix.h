#pragma once

#pragma once
#include <iostream>
#include <cmath>

enum ImprovedMatrixType
{
	Nothingg, Symmetricall, AntiSymmetricall, Diagonall, UpperTriangularr, LowerTriangularr
};

void REQUIRE(bool condition, const char* string) {
	if (!condition)throw string;
}

template <class T>
class ImprovedMatrix
{
private:
	T* m_arr = nullptr;
	int m_rows, m_cols;

	//removes a certain row and a cetain column, and returns a new matrix identical to the previuos one, expet for the row and the col that have been removed
	ImprovedMatrix<T> removeColAndRow(int r, int c) {
		REQUIRE(r >= 0 && r < m_rows && c >= 0 && c < m_cols, "Out of Bound Access");
		REQUIRE(m_cols == m_rows, "This is not a square Matrix");
		T* array = new T[(m_rows - 1)*(m_cols - 1)];
		size_t x = 0, y = 0;
		for (int row = 0; row < m_rows; row++)
		{
			y = 0;
			if (row == r)continue;
			for (int col = 0; col < m_cols; col++) {
				if (col == c)continue;
				array[(m_cols - 1)*x + y] = m_arr[m_cols*row + col];
				y++;
			}
			x++;
		}
		ImprovedMatrix<T> res(array, m_rows - 1, m_cols - 1);

		delete[] array;
		return res;
	}
	//returns a copy of the matrix which is reduced in echelon form
	ImprovedMatrix<T> returnEchelonForm() {
		ImprovedMatrix<T> res(*this);
		res.reduceToEchelonForm();
		return res;
	}

	//Swaps two rows within the matrix
	void swapRows(size_t row1, size_t row2) {
		REQUIRE(row1 > 0 && row1 <= m_rows && row2 > 0 && row2 <= m_rows, "Out of Bounds");
		for (register int col = 0; col < m_cols; col++) {
			T temp = m_arr[m_cols*row1 + col];
			m_arr[m_cols * row1 + col] = m_arr[m_cols * row2 + col];
			m_arr[m_cols * row2 + col] = temp;
		}
	}

	//COUNTING THE NUMBER OF ZEROS BEFORE THE FIRST NON-NULL ELEMENT FOR A SPECIFIED ROW
	size_t n_of_leading_zeros(size_t row) {
		size_t counter = 0;
		for (register int j = 0; j < m_cols; ++j) {
			if (m_arr[m_cols*row + j] == 0)counter++;
			else return counter;
		}
		return counter;
	}

	//returns true in case of a floating point calculation approximation where a 0 is left as a very small floating point number
	bool false_zero(T n) {
		//only works with standard types
		if (std::abs(n) < 0.0000001) return true;
		return false;
	}
	//it reorders the rows depending on the number of leading zeros (more zeros on the bottom)
	void orderByLeadingZeros() {
		int most_zeros = -1;
		size_t n_of_zeros = 0;
		size_t counter = 0;
		for (register int k = m_rows; k > 0; --k) {
			for (register int i = 0; i < k; ++i) {
				counter = n_of_leading_zeros(i);
				if (counter >= n_of_zeros) {
					most_zeros = i;
					n_of_zeros = counter;
				}

				counter = 0;
			}
			//added a modification (k instead of m_rows);
			if (most_zeros != k - 1 && most_zeros >= 0)swapRows(most_zeros, k - 1);
			n_of_zeros = 0;
			most_zeros = -1;
		}
	}

public:
	void print_array() {
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++) {
				std::cout << m_arr[m_cols*i + j] << ',';
			}
			std::cout << '\n';
		}
		std::cout << std::endl;
	}
	void GOD(ImprovedMatrix<T> & temp) {
		add(temp);
		print_array();
		subtract(temp);
		print_array();
		multiplyByScalar(10.0);
		print_array();
		std::cout << trace();
		ImprovedMatrix<T> trans = transposedMatrix();
		trans.print_array();
		std::cout << determinant() << '\n';
		std::cout << rank() << '\n';
		ImprovedMatrix<T> inv = inverseMatrix();
		inv.print_array();
		std::cout << categorizeMatrix();
	}

	//SEARCHES FOR PIVOT IN A SPECIFIED ROW
	T findPivot(size_t row) {
		T val;
		for (register int j = 0; j < m_cols; j++) {
			val = m_arr[m_cols*row + j];
			if (val != 0)return val;
		}
		return 0;
	}

	//searches for very small elements that can be aproximated to 0 and chanes their value to 0
	void check_and_clear_almost_zeros() {
		for (int i = 0; i < m_rows*m_cols; i++)
		{
			if (false_zero(m_arr[i]))m_arr[i] = 0;
		}
	}

	//////////////////////////////////////////////////////////////////////////CONSTRUCTORS AND DESTRUCTORS////////////////////////////////////////7///////////////////
	ImprovedMatrix() = delete;
	ImprovedMatrix(int m, int n)
		:m_rows(m), m_cols(n)
	{
		m_arr = new T[m * n];
		memset(m_arr, 0, sizeof(T)*(m * n));
	}
	ImprovedMatrix(T** arr, int m, int n)
		:m_rows(m), m_cols(n)
	{
		m_arr = new T[m*n];
		//@TODO: profile the cache miss/hit for this setup
		for (int row = 0; row < m_rows; row++)
		{
			for (int col = 0; col < m_cols; col++) {
				m_arr[m_cols*row + col] = arr[row][col];
			}
		}
	}
	ImprovedMatrix(T* arr, int m, int n)
		:m_rows(m), m_cols(n)
	{
		m_arr = new T[m*n];
		memcpy(m_arr, arr, sizeof(T)*(m*n));
	}
	ImprovedMatrix(const ImprovedMatrix& mat)
		:m_rows(mat.getRows()), m_cols(mat.getCols())
	{
		m_arr = new T[m_rows*m_cols];

		memcpy(m_arr, mat.getMatrixPointer(), sizeof(T)*(m_rows*m_cols));
	}
	ImprovedMatrix(ImprovedMatrix&& mat)
		:m_rows(mat.getRows()), m_cols(mat.getCols())
	{
		m_arr = mat.getMatrixPointer();
		mat.nullify();
	}
	~ImprovedMatrix() {
		if (m_arr) {
			delete[] m_arr;
		}
	}

	//////////////////////////////////////////////////////////////////////////CLASSES//////////////////////////////////////////////////////////////////////////

	class MITM {
	public:
		MITM(T* arr, int n_) : arr(arr), n(n_) { }

		T& operator[](int index) {
			REQUIRE(index >= 0 && index < n, "Out Of Bounds");
			return arr[index];
		}

	private:
		T* arr = nullptr;
		int n = 0;
	};

	//////////////////////////////////////////////////////////////////////////UTILITY METHODS//////////////////////////////////////////////////////////////////////////

	//returns number of rows
	inline int getRows() const { return m_rows; }

	//returns number of columns
	inline int getCols() const { return m_cols; }

	//returns the pointer to the first element of the array, so the array itself
	inline T* getMatrixPointer() const { return m_arr; }

	//sets the pointer to the array at null value
	inline void nullify() { m_arr = nullptr; }

	//gets the selected cell and puts the object in it
	void setElement(size_t i, size_t j, T el) {
		REQUIRE(i >= 0 || i < m_rows || j >= 0 || j < m_cols, "OutOfBoundAccess");
		m_arr[m_cols*i + j] = el;
	}
	//////////////////////////////////////////////////////////////////////////OPERATION FUNCTIONS//////////////////////////////////////////////////////////////////////////
	//adds a matrix to the current object
	void add(ImprovedMatrix& mat) {
		if (mat.getRows() != m_rows || mat.getCols() != m_cols)return;
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				m_arr[m_cols* i + j] += mat[i][j];
			}
		}
	}
	//subtracts a matrix from the current object
	void subtract(ImprovedMatrix& mat) {
		if (mat.getRows() != m_rows || mat.getCols() != m_cols)return;
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				m_arr[m_cols* i + j] -= mat[i][j];
			}
		}
	}

	// multiplies the elements of a specified row for a scalar
	void multiplyRowByScalar(size_t row, T scalar) {
		if (row < 0 || row >= m_rows)return;
		for (register int i = 0; i < m_cols; ++i) {
			m_arr[m_cols*row + i] *= scalar;
		}
	}
	//multiplies the whole matrix by a scalar
	void multiplyByScalar(T scalar) {
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				m_arr[m_cols* i + j] *= scalar;
			}
		}
	}
	//checks whether the two matrices are equal
	bool equal(ImprovedMatrix& mat) {
		if (mat.getRows() != m_rows || mat.getCols() != m_cols)return false;
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				if (m_arr[m_cols* i + j] != mat[i][j])return false;
			}
		}
		return true;
	}

	//selected a row, it sums to it another selected row
	void addRowToRow(size_t destination, size_t source, T multiplication_factor = 1) {
		if (destination < 0 || destination >= m_rows || source < 0 || source >= m_rows)return;
		for (int j = 0; j < m_cols; j++) {
			m_arr[m_cols*destination + j] = (m_arr[m_cols*destination + j] + m_arr[m_cols*source + j] * multiplication_factor);
		}
	}

	T trace() {
		REQUIRE(m_rows == m_cols, "Cannot calculate trace of a non-square matrix");
		T res = 0;
		for (int i = 0; i < m_rows; i++)
		{
			res += m_arr[m_cols* i + i];
		}

		return res;
	}

	//creates a new matrix that is the transposed of the original one
	ImprovedMatrix<T> transposedMatrix() {
		ImprovedMatrix<T> res(m_cols, m_rows);

		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				res[j][i] = m_arr[m_cols*i + j];
			}
		}
		return res;
	}
	//it returns the matrix result of the Row by Column product between this and another matrix
	ImprovedMatrix<T> rowByColProduct(ImprovedMatrix& mat) {
		REQUIRE(m_cols == mat.getRows(), "Matrices' sizes don't match the requirements");
		ImprovedMatrix<T> result(m_rows, mat.getCols());
		int counter = 0;
		T temp = 0;
		for (int i = 0; i < m_rows; i++)
		{
			counter = 0;
			while (counter < mat.getCols()) {
				for (int j = 0; j < m_cols; j++)
				{
					temp += m_arr[m_cols*i + j] * mat.getMatrixPointer()[mat.getCols()*j + counter];
				}
				result[i][counter] = temp;
				temp = 0;
				counter++;
			}
		}
		return result;
	}
	//returns the result of the scalar product between this and another matrix
	T scalarProduct(ImprovedMatrix& mat) {
		ImprovedMatrix<T> temp = mat.transposedMatrix();
		ImprovedMatrix<T> res = rowByColProduct(temp);
		return res.trace();
	}
	//it returns the DETERMINANT of the current matrix
	T determinant() {
		REQUIRE(m_rows == m_cols, "Cannot calculate determinant of non-square matrix");
		if (m_rows == 1 && m_cols == 1) return m_arr[0];
		if (m_rows == 2 && m_cols == 2) {
			return m_arr[0] * m_arr[3] - m_arr[1] * m_arr[2];
		}

		T res = 0;
		ImprovedMatrix<T> temp(m_rows - 1, m_cols - 1);
		for (int i = 0; i < m_cols; i++) {
			temp = removeColAndRow(0, i);
			res += (std::pow(-1, i)*m_arr[i] * temp.determinant());
		}
		return res;
	}
	//it returns the RANK of the current matrix
	int rank() {
		ImprovedMatrix<T> tmp = returnEchelonForm();
		int counter = 0;
		for (register int i = 0; i < m_rows; ++i) {
			if (tmp.findPivot(i) != 0)counter++;
			else return counter;
		}
		return counter;
	}
	//it returns the inverse matrix of the current one
	ImprovedMatrix<T> inverseMatrix() {
		T det = determinant();
		ImprovedMatrix<T> res(m_rows, m_cols);
		if (det != 0) {
			ImprovedMatrix<T> temp(m_rows - 1, m_cols - 1);
			for (int i = 0; i < m_rows; i++) {
				for (int j = 0; j < m_cols; j++) {
					temp = removeColAndRow(j, i);
					res[i][j] = std::pow(-1, i + j)*temp.determinant() / det;
				}
			}
		}

		return res;
	}
	//returns true if the current matrix is orthogonal (i.e if its inverse and its transposed are equal)
	bool isOrthogonal() {
		ImprovedMatrix<T> inverse = inverseMatrix();
		ImprovedMatrix<T> transpose = transposedMatrix();
		return inverse == transpose;
	}
	//not the most efficient....
	//complexity : O(n^4)
	//reduces the current matrix to echelon form
	void reduceToEchelonForm() {
		T pivot;
		orderByLeadingZeros();
		for (int i = 1; i < m_rows; i++) {
			if (n_of_leading_zeros(i - 1) < n_of_leading_zeros(i))continue;
			pivot = findPivot(i - 1);
			if (pivot == 0)continue; //OR BREAK, I DON'T KNOW WHAT IS ADVISABLE
			for (int j = i; j < m_rows; j++) {
				if (n_of_leading_zeros(i - 1) < n_of_leading_zeros(j))continue;
				addRowToRow(j, i - 1, -(findPivot(j) / pivot));
				//print_array();
			}
			orderByLeadingZeros();
		}
		check_and_clear_almost_zeros();
	}
	//returns true if the current matrix is symmetrical
	bool isSymmetrical() {
		if (m_cols != m_rows)return false;
		ImprovedMatrix<T> temp = transposedMatrix();
		return temp == *this;
	}
	//returns true if the current Matrix is AntiSymmetrical
	bool isAntisymmetrical() {
		if (m_cols != m_rows)return false;
		ImprovedMatrix<T> temp = transposedMatrix();
		temp *= (-1.0);
		return temp == *this;
	}
	//checks if the matrix is upper triangular
	bool isUpperTriangular() {
		if (m_cols != m_rows)return false;
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_cols; j++) {
				if (i <= j)continue;
				else {
					if (m_arr[m_cols*i + j] != 0)return false;
				}
			}
		}
		return true;
	}
	//returns true if the current matrix is LowerTriangular
	bool isLowerTriangular() {
		if (m_cols != m_rows)return false;
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_cols; j++) {
				if (i >= j)continue;
				else {
					if (m_arr[i*m_cols + j] != 0)return false;
				}
			}
		}
		return true;
	}
	//returns true if the current matrix is Diagonal
	bool isDiagonal() {
		if (m_cols != m_rows)return false;
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_cols; j++) {
				if (i == j)continue;
				else {
					if (m_arr[i*m_cols + j] != 0)return false;
				}
			}
		}
		return true;
	}

	ImprovedMatrixType categorizeMatrix() {
		if (isDiagonal())return ImprovedMatrixType::Diagonall;
		else if (isLowerTriangular())return ImprovedMatrixType::LowerTriangularr;
		else if (isUpperTriangular())return ImprovedMatrixType::UpperTriangularr;
		else if (isSymmetrical())return ImprovedMatrixType::Symmetricall;
		else if (isAntisymmetrical())return ImprovedMatrixType::AntiSymmetricall;
		else return ImprovedMatrixType::Nothingg;
	}

	void pow(int n) {
		while (--n > 0) {
			(*this) *= (*this);
		}
	}

	//////////////////////////////////////////////////////////////////////////OPERATORSSSSSSS///////////////////////////////////////////////////////////////////

	inline MITM operator[](int index) const {
		REQUIRE(index >= 0 && index < m_rows, "Out of Bounds");
		return MITM(m_arr + (m_cols*index), m_cols);
	}

	inline void operator=(const ImprovedMatrix& mat) {
		m_rows = mat.getRows();
		m_cols = mat.getCols();
		if (m_arr)
			delete[] m_arr;
		m_arr = new T[m_rows * m_cols];
		memcpy(m_arr, mat.getMatrixPointer(), sizeof(T)*(m_rows * m_cols));
	}
	inline void operator=(ImprovedMatrix&& mat) {
		m_rows = mat.getRows();
		m_cols = mat.getCols();
		if (m_arr)
			delete[] m_arr;
		m_arr = mat.getMatrixPointer();
		mat.nullify();
	}
	inline bool operator==(ImprovedMatrix& rhs) {
		return equal(rhs);
	}
	inline bool operator!=(ImprovedMatrix& rhs) {
		return !equal(rhs);
	}
	inline ImprovedMatrix<T> operator+(const ImprovedMatrix& rhs) {
		ImprovedMatrix<T> res(*this);
		res.add(rhs);
		return res;
	}
	inline ImprovedMatrix<T> operator-(const ImprovedMatrix& rhs) {
		ImprovedMatrix<T> res(*this);
		res.subtract(rhs);
		return res;
	}
	inline ImprovedMatrix<T>& operator+=(const ImprovedMatrix& rhs) {
		add(rhs);
		return *this;
	}
	inline ImprovedMatrix<T>& operator-=(const ImprovedMatrix& rhs) {
		subtract(rhs);
		return *this;
	}
	inline ImprovedMatrix<T> operator*(const T& scalar) {
		ImprovedMatrix<T> res(*this);
		res.multiplyByScalar(scalar);
		return res;
	}
	inline ImprovedMatrix<T> operator/(const T& scalar) {
		ImprovedMatrix<T> res(*this);
		res.multiplyByScalar(1.0 / scalar);
		return res;
	}
	inline ImprovedMatrix<T>& operator*=(const T& scalar) {
		multiplyByScalar(scalar);
		return *this;
	}
	inline ImprovedMatrix<T>& operator*=(ImprovedMatrix& rhs) {
		*this = rowByColProduct(rhs);
		return *this;
	}
	inline ImprovedMatrix<T>& operator/=(const T& scalar) {
		multiplyByScalar(1 / scalar);
		return *this;
	}
};
