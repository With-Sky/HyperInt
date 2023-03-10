#pragma once
#include <chrono>
class StopWatch
{
private:
	bool is_stop = false, is_start = false;
	uint64_t tick = 0;
	double rate = 1.0f;
	std::chrono::system_clock::time_point begin;
	std::chrono::system_clock::time_point end;

public:
	StopWatch()
	{
		is_start = false;
		is_stop = false;
		rate = 1.0f;
		tick = 0;
	}
	StopWatch(double rate_in)
	{
		is_start = false;
		is_stop = false;
		rate = rate_in;
		tick = 0;
	}
	void start()
	{
		reset();
		is_start = true;
		is_stop = false;
		begin = std::chrono::system_clock::now();
	}
	void stop()
	{
		if (is_start)
		{
			is_stop = true;
			end = std::chrono::system_clock::now();
			auto delta_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
			tick = delta_time.count();
		}
	}
	void reset()
	{
		is_start = false;
		is_stop = false;
		tick = 0;
	}
	double duration()
	{
		if (!is_stop)
		{
			stop();
			is_stop = false;
		}
		return static_cast<double>(tick / (rate * 1000));
	}
};
StopWatch watch_default(1);