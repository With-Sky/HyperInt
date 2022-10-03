#include <iostream>
#include "hint.hpp"
#include "stopwatch.hpp"

using namespace std;

int main(int argc, char **argv)
{
	StopWatch w(1000); //初始化计时器，单位tick为1000us
	w.start();
	HyperInt a = factorial(50000);
	w.stop();
	double t1 = w.duration();
	w.start();
	a = factorial(100000);
	w.stop();
	double t2 = w.duration();
	w.start();
	a = factorial(1000000);
	w.stop();
	double t3 = w.duration();
	cout << t1 << "ms\t" << t2 << "ms\t" << t3 << "ms";
	return 0;
}