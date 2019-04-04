#pragma once
class ComplexNumber
{
private:
	double m_real = 0;
	double m_i = 0;

public:
	ComplexNumber() {};
	ComplexNumber(const ComplexNumber& other)
		:m_real(other.getReal()), m_i(other.getImaginary()) {}
	ComplexNumber(int n)
		:m_real(n) {}

	~ComplexNumber() {};
	double getReal() const { return m_real; }
	double getImaginary() const { return m_i; }
	void add(const ComplexNumber& n) {
		m_real += n.getReal();
		m_i += n.getImaginary();
	}
	void subtract(const ComplexNumber& n) {
		m_real -= n.getReal();
		m_i -= n.getImaginary();
	}
	void add(const int& n) {
		m_real += n;
	}
	void subtract(const int& n) {
		m_real -= n;
	}
	void add(const double& n) {
		m_real += n;
	}
	void subtract(const double& n) {
		m_real -= n;
	}

	//////////////////////////////////////////////////////////////////////////OPERATORS
	void operator=(const ComplexNumber& n) {
		m_real = n.getReal();
		m_i = n.getImaginary();
	}
	void operator=(int n) {
		m_real = n;
		m_i = 0;
	}

	inline ComplexNumber operator+(const ComplexNumber& rhs) {
		ComplexNumber res(*this);
		res.add(rhs);
		return res;
	}
	inline ComplexNumber operator+(const int& rhs) {
		ComplexNumber res(*this);
		res.add(rhs);
		return res;
	}
	inline ComplexNumber operator+(const double& rhs) {
		ComplexNumber res(*this);
		res.add(rhs);
		return res;
	}
	inline ComplexNumber operator-(const ComplexNumber& rhs) {
		ComplexNumber res(*this);
		res.subtract(rhs);
		return res;
	}
	inline ComplexNumber operator-(const int& rhs) {
		ComplexNumber res(*this);
		res.subtract(rhs);
		return res;
	}
	inline ComplexNumber operator-(const double& rhs) {
		ComplexNumber res(*this);
		res.subtract(rhs);
		return res;
	}
	inline ComplexNumber& operator+=(const ComplexNumber& rhs) {
		add(rhs);
		return *this;
	}
	inline ComplexNumber& operator-=(const ComplexNumber& rhs) {
		subtract(rhs);
		return *this;
	}
	inline ComplexNumber& operator+=(const int& rhs) {
		add(rhs);
		return *this;
	}
	inline ComplexNumber& operator-=(const int& rhs) {
		subtract(rhs);
		return *this;
	}
	inline ComplexNumber& operator+=(const double& rhs) {
		add(rhs);
		return *this;
	}
	inline ComplexNumber& operator-=(const double& rhs) {
		subtract(rhs);
		return *this;
	}
};
