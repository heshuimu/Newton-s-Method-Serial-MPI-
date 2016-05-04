//
//  main.cpp
//  Rosenbrock
//
//  Created by 雪竜 on 16/4/27.
//
//

#include <iostream>
#include <cmath>

#define __ITER__ 12

double fx(double x, double y)
{
	return 2*(200*x*x*x - 200*x*y + x - 1);
}

double fy(double x, double y)
{
	return 200*(y-x*x);
}

double fxx(double x, double y)
{
	return 1200*x*x - 400*y + 2;
}

double fxy(double x, double y)
{
	return -400*x;
}

double fyx(double x, double y)
{
	return -400*x;
}

double fyy(double x, double y)
{
	return 200;
}

int main(int argc, const char * argv[]) {
	double x = -1.2;
	double y  = 1;
	
	for (int i = 0; i < __ITER__; i++) {
	
		double
			A = fxx(x, y),
			B = fxy(x, y),
			C = fyx(x, y),
			D = fyy(x, y),
			P = -fx(x, y),
			Q = -fy(x, y);
		
		double delta_x = (P*D - Q*B) / (A*D - C*B);
		double delta_y = (Q - delta_x*C) / D;
		
		x+=delta_x;
		y+=delta_y;
		
		std::cout << i << "," << sqrt(delta_x*delta_x + delta_y*delta_y) << "," << x << "," << y << "\n";
	}
	
	std::cout << "Result: " << x <<", "<< y << "\n";
	
	return 0;
}
