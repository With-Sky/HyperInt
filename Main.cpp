#include <iostream>
#include "hint.hpp"
#include "stopwatch.hpp"

using namespace std;

int main(int argc, char **argv)
{
	StopWatch w(1000);				//创建一个单位tick为1000us的时钟
	w.start();						//开始计时
	HyperInt a = factorial(500000); //计算500000!
	w.stop();						//结束计时
	double t1 = w.duration();
	w.start();
	string s = a.to_string();       //转10进制字符串 to_string(a) is fine
	w.stop();
	double t2 = w.duration();
	cout << s << endl;
	cout << t1 << "ms\t" << t2 << "ms";
	return 0;
}