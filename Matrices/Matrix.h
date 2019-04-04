#pragma once
#include <cmath>

enum MatrixType
{
	Nothing, Symmetrical, AntiSymmetrical, Diagonal, UpperTriangular, LowerTriangular
};

template <class T>
class Matrix
{
private:
	T** m_arr = nullptr;
	int m_rows, m_cols;
	Matrix<T> removeColAndRow(int r, int c) {
		if (r < 0 || r > m_rows - 1 || c < 0 || c > m_cols - 1)throw "Out of Bound Access";
		T** array = new T*[m_rows - 1];
		for (register int i = 0; i < m_rows - 1; i++) {
			array[i] = new T[m_cols - 1];
		}
		size_t x = 0, y = 0;
		for (int i = 0; i < m_rows; i++)
		{
			y = 0;
			if (i == r)continue;
			for (int j = 0; j < m_cols; j++) {
				if (j == c)continue;
				array[x][y] = m_arr[i][j];
				y++;
			}
			x++;
		}
		Matrix<T> res(array, m_rows - 1, m_cols - 1);
		for (register int i = 0; i < m_rows - 1; ++i) {
			delete[] array[i];
		}
		delete[] array;
		return res;
	}

	T findPivot(size_t row) {
		T val;
		for (register int j = 0; j < m_cols; j++) {
			val = m_arr[row][j];
			if (val != 0)return val;
		}
		return 0;
	}

public:

	//////////////////////////////////////////////////////////////////////////CONSTRUCTORS AND DESTRUCTORS////////////////////////////////////////7///////////////////

	Matrix() = delete;
	Matrix(int m, int n)
		:m_rows(m), m_cols(n)
	{
		m_arr = new T*[m];
		for (register int i = 0; i < m; i++) {
			m_arr[i] = new T[n];
		}

		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++) {
				m_arr[i][j] = 0;
			}
		}
	}
	Matrix(T** arr, int m, int n)
		:m_rows(m), m_cols(n)
	{
		m_arr = new T*[m];
		for (register int i = 0; i < m; i++) {
			m_arr[i] = new T[n];
		}
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++) {
				m_arr[i][j] = arr[i][j];
			}
		}
	}
	Matrix(const Matrix& mat)
		:m_rows(mat.getRows()), m_cols(mat.getCols())
	{
		m_arr = new T*[m_rows];
		for (register int i = 0; i < m_rows; i++) {
			m_arr[i] = new T[m_cols];
		}
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++) {
				m_arr[i][j] = mat[i][j];
			}
		}
	}
	Matrix(Matrix&& mat)
		:m_rows(mat.getRows()), m_cols(mat.getCols())
	{
		m_arr = mat.getMatrixPointer();
		mat.nullify();
	}
	~Matrix() {
		if (m_arr) {
			for (register int i = 0; i < m_rows; ++i) {
				delete[] m_arr[i];
			}
			delete[] m_arr;
		}
	}

	//////////////////////////////////////////////////////////////////////////CLASSES//////////////////////////////////////////////////////////////////////////

	class MITM {
	public:
		MITM(T* arr, int n_) : arr(arr), n(n_) { }

		T& operator[](int index) {
			if (index >= 0 && index < n) {
				return arr[index];
			}
			else
			{
				throw "Out of Bounds";
			}
		}
	private:
		T* arr = nullptr;
		int n = 0;
	};

	//////////////////////////////////////////////////////////////////////////UTILITY METHODS//////////////////////////////////////////////////////////////////////////

	inline int getRows() { return m_rows; }
	inline int getCols() { return m_cols; }
	inline T** getMatrixPointer() { return m_arr; }
	inline void nullify() { m_arr = nullptr; }
	void setElement(size_t i, size_t j, T el) {
		if (i < 0 || i > m_rows - 1 || j < 0 || j > m_cols - 1) throw "OutOfBoundAccess";
		m_arr[i][j] = el;
	}
	//////////////////////////////////////////////////////////////////////////OPERATION FUNCTIONS//////////////////////////////////////////////////////////////////////////

	void add(Matrix& mat) {
		if (mat.getRows() != m_rows || mat.getCols() != m_cols)return;
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				m_arr[i][j] += mat[i][j];
			}
		}
	}
	void subtract(Matrix& mat) {
		if (mat.getRows() != m_rows || mat.getCols() != m_cols)return;
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				m_arr[i][j] -= mat[i][j];
			}
		}
	}
	void multiplyRowByScalar(size_t row, T scalar) {
		if (row < 0 || row >= m_rows)return;
		for (register int i = 0; i < m_cols; ++i) {
			m_arr[row][i] *= scalar;
		}
	}

	void multiplyByScalar(T scalar) {
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				m_arr[i][j] *= scalar;
			}
		}
	}
	bool equal(Matrix& mat) {
		if (mat.getRows() != m_rows || mat.getCols() != m_cols)return false;
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				if (m_arr[i][j] != mat[i][j])return false;
			}
		}
		return true;
	}
	void addRowToRow(size_t destination, size_t source, T multiplication_factor = 1) {
		if (destination < 0 || destination >= m_rows || source < 0 || source >= m_rows)return;
		for (int j = 0; j < m_cols; j++) {
			m_arr[destination][j] += m_arr[source][j] * multiplication_factor;
		}
	}

	T trace() {
		T res = 0;
		if (m_rows == m_cols) {
			for (int i = 0; i < m_rows; i++)
			{
				res += m_arr[i][i];
			}
		}
		return res;
	}

	Matrix<T> transposedMatrix() {
		Matrix<T> res(m_cols, m_rows);

		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++)
			{
				res[j][i] = m_arr[i][j];
			}
		}
		return res;
	}
	//seems to be working
	Matrix<T> rowByColProduct(Matrix& mat) {
		if (m_cols != mat.getRows())throw "Matrix sizes don't match the requirements";
		Matrix<T> result(m_rows, mat.getCols());
		int counter = 0;
		T temp = 0;
		for (int i = 0; i < m_rows; i++)
		{
			counter = 0;
			while (counter < mat.getCols()) {
				for (int j = 0; j < m_cols; j++)
				{
					temp += m_arr[i][j] * mat.m_arr[j][counter];
				}
				result[i][counter] = temp;
				temp = 0;
				counter++;
			}
		}
		return result;
	}

	T scalarProduct(Matrix& mat) {
		Matrix<T> temp = mat.transposedMatrix();
		Matrix<T> res = rowByColProduct(temp);
		return res.trace();
	}
	T determinant() {
		if (m_rows == 2 && m_cols == 2) {
			return m_arr[0][0] * m_arr[1][1] - m_arr[0][1] * m_arr[1][0];
		}

		T res = 0;
		Matrix<T> temp(m_rows - 1, m_cols - 1);
		for (int i = 0; i < m_cols; i++) {
			temp = removeColAndRow(0, i);
			res += (std::pow(-1, i)*m_arr[0][i] * temp.determinant());
		}
		return res;
	}
	Matrix<T> inverseMatrix() {
		Matrix<T> res(m_rows, m_cols);
		Matrix<T> temp(m_rows - 1, m_cols - 1);
		T det = determinant();
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_cols; j++) {
				temp = removeColAndRow(j, i);
				res[i][j] = std::pow(-1, i + j)*temp.determinant() / det;
			}
		}

		return res;
	}

	bool isOrthogonal() {
		Matrix<T> inverse = inverseMatrix();
		Matrix<T> transpose = transposedMatrix();
		return inverse == transpose;
	}

	void reduceToEchelonForm() {
		T pivot;
		for (int i = 1; i < m_rows; i++) {
			pivot = findPivot(i - 1);
			for (int j = i; j < m_rows; j++) {
				addRowToRow(j, i - 1, -(findPivot(j) / pivot));
			}
		}
	}

	bool isSymmetrical() {
		Matrix<T> temp = transposedMatrix();
		return temp == *this;
	}
	bool isAntisymmetrical() {
		Matrix<T> temp = transposedMatrix();
		temp *= (-1.0);
		return temp == *this;
	}
	bool isUpperTriangular() {
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_rows; j++) {
				if (i <= j)continue;
				else {
					if (m_arr[i][j] != 0)return false;
				}
			}
		}
		return true;
	}
	bool isLowerTriangular() {
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_rows; j++) {
				if (i >= j)continue;
				else {
					if (m_arr[i][j] != 0)return false;
				}
			}
		}
		return true;
	}
	bool isDiagonal() {
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_rows; j++) {
				if (i == j)continue;
				else {
					if (m_arr[i][j] != 0)return false;
				}
			}
		}
		return true;
	}

	MatrixType categorizeMatrix() {
		if (isDiagonal())return MatrixType::Diagonal;
		else if (isLowerTriangular())return MatrixType::LowerTriangular;
		else if (isUpperTriangular())return MatrixType::UpperTriangular;
		else if (isSymmetrical())return MatrixType::Symmetrical;
		else if (isAntisymmetrical())return MatrixType::AntiSymmetrical;
		else return MatrixType::Nothing;
	}

	//////////////////////////////////////////////////////////////////////////OPERATORSSSSSSS///////////////////////////////////////////////////////////////////

	inline MITM operator[](int index) {
		if (index >= 0 && index < m_rows) {
			return MITM(m_arr[index], m_cols);
		}
		else {
			throw "Out of bounds";
		}
	}

	inline void operator=(const Matrix& mat) {
		m_rows = mat.getRows();
		m_cols = mat.getCols();
		if (m_arr) {
			for (register int i = 0; i < m_rows; ++i) {
				delete[] m_arr[i];
			}
			delete[] m_arr;
		}
		m_arr = new T*[m_rows];
		for (register int i = 0; i < m_rows; i++) {
			m_arr[i] = new T[m_cols];
		}
		for (int i = 0; i < m_rows; i++)
		{
			for (int j = 0; j < m_cols; j++) {
				m_arr[i][j] = mat[i][j];
			}
		}
	}
	inline void operator=(Matrix&& mat) {
		m_rows = mat.getRows();
		m_cols = mat.getCols();
		if (m_arr) {
			for (register int i = 0; i < m_rows; ++i) {
				delete[] m_arr[i];
			}
			delete[] m_arr;
		}
		m_arr = mat.getMatrixPointer();
		mat.nullify();
	}
	inline bool operator==(Matrix& rhs) {
		return equal(rhs);
	}
	inline bool operator!=(Matrix& rhs) {
		return !equal(rhs);
	}
	inline Matrix<T> operator+(const Matrix& rhs) {
		Matrix<T> res(*this);
		res.add(rhs);
		return res;
	}
	inline Matrix<T> operator-(const Matrix& rhs) {
		Matrix<T> res(*this);
		res.subtract(rhs);
		return res;
	}
	inline Matrix<T>& operator+=(const Matrix& rhs) {
		add(rhs);
		return *this;
	}
	inline Matrix<T>& operator-=(const Matrix& rhs) {
		subtract(rhs);
		return *this;
	}
	inline Matrix<T> operator*(const T& scalar) {
		Matrix<T> res(*this);
		res.multiplyByScalar(scalar);
		return res;
	}
	inline Matrix<T> operator/(const T& scalar) {
		Matrix<T> res(*this);
		res.multiplyByScalar(1.0 / scalar);
		return res;
	}
	inline Matrix<T>& operator*=(const T& scalar) {
		multiplyByScalar(scalar);
		return *this;
	}
	inline Matrix<T>& operator/=(const T& scalar) {
		multiplyByScalar(1 / scalar);
		return *this;
	}
};
