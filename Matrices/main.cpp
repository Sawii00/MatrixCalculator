#include <iostream>
#include "Vector.h"

int main() {
	Vector<4> vec = { 2 };
	Vector<4> vec2({ 1, 2, 3, 4 });
	vec.print(std::cout);
	vec2.print(std::cout);

	std::cout << vec2.scalar(vec);
	system("pause");
}