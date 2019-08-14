#pragma once

#include <initializer_list>
#include <immintrin.h>
#include <malloc.h>

void REQUIRE(bool condition, const char* string) {
	if (!condition)throw string;
}

template <size_t T>
class Vector {
private:
	float* vec;
	size_t m_size;
public:

	Vector(std::initializer_list<float> list) {
		if (list.size() != 1) {
			REQUIRE(list.size() == T, "Wrong number of parameters");

			vec = static_cast<float*>(_aligned_malloc(sizeof(float)*T, sizeof(float)));
			m_size = T;
			for (register int i = 0; i < T; ++i) {
				vec[i] = *(list.begin() + i);
			}
		}
		else {
			vec = static_cast<float*>(_aligned_malloc(sizeof(float)*T, sizeof(float)));
			m_size = T;
			for (register int i = 0; i < T; i++) {
				vec[i] = *list.begin();
			}
		}
	}
	Vector() {
		vec = nullptr;
		m_size = T;
	}
	Vector(const Vector& vec2) {
		if (vec)_aligned_free(vec);
		vec = static_cast<float*>(_aligned_malloc(sizeof(float)*vec2.size(), sizeof(float)));
		m_size = vec2.size();
		memcpy(vec, vec2.getPointer(), sizeof(float)*vec2.size());
	}
	Vector(Vector&& vec2) {
		m_size = vec2.size();
		vec = vec2.getPointer();
		vec2.nullify();
	}

	~Vector() {
		if (vec)_aligned_free(vec);
	};

	inline size_t size() const {
		return m_size;
	}

	inline void nullify() {
		vec = nullptr;
	}
	const float* const getPointer() const { return vec; }

	void print(std::ostream& o) {
		for (register int i = 0; i < size(); i++) {
			o << vec[i] << ' ';
		}
		o << '\n';
	}

	inline float& operator[](int index) {
		REQUIRE(index >= 0 && index < size(), "Out of Bounds");
		return vec[index];
	}
	inline float operator[] (int index) const {
		REQUIRE(index >= 0 && index < size(), "Out of Bounds");
		return vec[index];
	}

	void operator+=(const Vector& vec2) {
		REQUIRE(vec2.size() == size(), "Sizes don't match");
		int remaining = size() % 8;
		int steps = size() / 8;

		for (register int i = 0; i < steps; i++) {
			__m256 first = _mm256_load_ps(vec + i * 8);
			__m256 second = _mm256_load_ps(vec2.getPointer() + i * 8);
			first = _mm256_add_ps(first, second);
			memcpy(vec + i * 8, (float*)&first, sizeof(float) * 8);
		}

		for (register int i = size() - remaining; i < remaining; i++) {
			vec[i] += *(vec2.getPointer() + i);
		}
	}

	void operator+=(const float scalar) {
		int remaining = size() % 8;
		int steps = size() / 8;

		for (register int i = 0; i < steps; i++) {
			__m256 first = _mm256_load_ps(vec + i * 8);
			__m256 second = _mm256_set1_ps(scalar);
			first = _mm256_add_ps(first, second);
			memcpy(vec + i * 8, (float*)&first, sizeof(float) * 8);
		}

		for (register int i = size() - remaining; i < remaining; i++) {
			vec[i] += scalar;
		}
	}

	void operator-=(const Vector& vec2) {
		REQUIRE(vec2.size() == size(), "Sizes don't match");
		int remaining = size() % 8;
		int steps = size() / 8;

		for (register int i = 0; i < steps; i++) {
			__m256 first = _mm256_load_ps(vec + i * 8);
			__m256 second = _mm256_load_ps(vec2.getPointer() + i * 8);
			first = _mm256_sub_ps(first, second);
			memcpy(vec + i * 8, (float*)&first, sizeof(float) * 8);
		}

		for (register int i = size() - remaining; i < remaining; i++) {
			vec[i] -= *(vec2.getPointer() + i);
		}
	}
	void operator-=(const float scalar) {
		int remaining = size() % 8;
		int steps = size() / 8;

		for (register int i = 0; i < steps; i++) {
			__m256 first = _mm256_load_ps(vec + i * 8);
			__m256 second = _mm256_set1_ps(scalar);
			first = _mm256_sub_ps(first, second);
			memcpy(vec + i * 8, (float*)&first, sizeof(float) * 8);
		}

		for (register int i = size() - remaining; i < remaining; i++) {
			vec[i] -= scalar;
		}
	}

	void operator *=(const float scalar) {
		int remaining = size() % 8;
		int steps = size() / 8;

		for (register int i = 0; i < steps; i++) {
			__m256 first = _mm256_load_ps(vec + i * 8);
			__m256 second = _mm256_set1_ps(scalar);
			first = _mm256_mul_ps(first, second);
			memcpy(vec + i * 8, (float*)&first, sizeof(float) * 8);
		}

		for (register int i = size() - remaining; i < remaining; i++) {
			vec[i] *= scalar;
		}
	}
	void operator /=(const float scalar) {
		int remaining = size() % 8;
		int steps = size() / 8;

		for (register int i = 0; i < steps; i++) {
			__m256 first = _mm256_load_ps(vec + i * 8);
			__m256 second = _mm256_set1_ps(scalar);
			first = _mm256_div_ps(first, second);
			memcpy(vec + i * 8, (float*)&first, sizeof(float) * 8);
		}

		for (register int i = size() - remaining; i < remaining; i++) {
			vec[i] /= scalar;
		}
	}

	int operator*(const Vector& vec2) {
		REQUIRE(size() == vec2.size(), "Sizes don't match");
		int res = 0;
		for (register int i = 0; i < size(); i++) {
			res += (vec[i] * vec2[i]);
		}
		return res;
	}

	int scalar(const Vector& vec2) {
		return (*this)*vec2;
	}
};