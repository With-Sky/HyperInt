#ifndef HINT_HPP
#define HINT_HPP

#include <iostream>
#include <future>
#include <cstring>
#include <cmath>

#define _MAX_(x, y) (x > y ? x : y)
#define _MIN_(x, y) (x < y ? x : y)
#define _NEG_(x) (x < 0)
#define _ODD_(x) (x & 1)
#define _ABS_(x) (x < 0 ? (-x) : x)
#define _TWICE_(x) (x << 1)
#define _HALF_(x) (x >> 1)
#define _SELFTWICE_(x) (x <<= 1)
#define _SELFHALF_(x) (x >>= 1)
#define _ARYCOPY_(target, source, len) std::memcpy(target, source, len * sizeof(*target))

#define SHORT_BIT 16
#define INT_BIT 32
#define LOG_INT_BIT 5
#define INT16_0XFF 0xFFFF
#define INT32_0XFF 0xFFFFFFFF
#define INT32_0X7F 0x7FFFFFFF
#define INT32_0X01 0x00000001
#define INT32_0X80 0X80000000
#define FFT_MIN 128
#define FFT_MAX 512000
#define PI 3.1415926535897932

#define UINT_8 unsigned char
#define UINT_16 unsigned short
#define UINT_32 unsigned int
#define UINT_64 unsigned long long
#define INT_32 int
#define INT_64 long long
#define REG register

namespace hint
{
    struct h_int
    {
        bool neg = false;
        UINT_32 *array = nullptr;
        size_t len, size;
    };
    //对复数的定义
    struct Complex
    {
        double real, imaginary;
        Complex()
        {
            real = imaginary = 0;
        }
        Complex(const Complex &input)
        {
            real = input.real;
            imaginary = input.imaginary;
        }
        Complex(double r, double i)
        {
            real = r;
            imaginary = i;
        }
        Complex(UINT_32 n)
        {
            real = cos(2 * PI / n);
            imaginary = sin(2 * PI / n);
        }
        Complex operator=(const Complex &input)
        {
            real = input.real;
            imaginary = input.imaginary;
            return *this;
        }
        Complex operator+(Complex input)
        {
            return Complex(real + input.real, imaginary + input.imaginary);
        }
        Complex operator-(Complex input)
        {
            return Complex(real - input.real, imaginary - input.imaginary);
        }
        Complex operator*(Complex input)
        {
            return Complex(real * input.real - imaginary * input.imaginary, real * input.imaginary + imaginary * input.real);
        }
        void console_out()
        {
            std::cout << real;
            if (_NEG_(imaginary))
                std::cout << imaginary << "i";
            else
                std::cout << "+" << imaginary << "i";
        }
    };
    void fft(Complex *input, size_t n, bool is_ifft)
    {
        size_t *rev = new size_t[n];
        REG size_t log_n = 1;
        REG Complex tmp, unit_omega, omega, tmp1, tmp2;
        while ((1 << log_n) < n)
        {
            log_n++;
        }
        rev[0] = 0;
        for (REG size_t i = 1; i < n; i++)
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log_n - 1));
        for (REG size_t i = 0; i < n; i++)
        {
            if (i < rev[i])
            {
                tmp = input[i];
                input[i] = input[rev[i]];
                input[rev[i]] = tmp;
            }
        }
        delete[] rev;
        for (REG size_t rank = 1, gap; rank < n; rank <<= 1)
        {
            gap = rank << 1;
            unit_omega = Complex(gap);
            if (is_ifft)
            {
                unit_omega.imaginary = -unit_omega.imaginary;
            }
            for (REG size_t begin = 0; begin < n; begin += gap)
            {
                omega = Complex(1, 0);
                for (REG size_t pos = begin; pos < begin + rank; pos++)
                {
                    tmp1 = input[pos];
                    tmp2 = input[pos + rank] * omega;
                    input[pos] = tmp1 + tmp2;
                    input[pos + rank] = tmp1 - tmp2;
                    omega = omega * unit_omega;
                }
            }
        }
        if (is_ifft)
        {
            for (REG size_t i = 0; i < n; i++)
                input[i].real /= 1.0 * n;
        }
    } //快速傅里叶变换
}
class HyperInt
{
private:
    hint::h_int data;
    short shorten_count, size_shorten_max = 64; //设定array长度不缩短次数阈值
    size_t true_len() const
    {
        REG size_t t_len = data.len;
        while (t_len > 0 && data.array[t_len - 1] == 0)
        {
            t_len--;
        }
        return t_len;
    }
    void reset_size(const size_t &new_size)
    {
        data.size = 2;
        while (data.size < new_size)
        {
            _SELFTWICE_(data.size);
        }
        data.len = _MIN_(data.len, data.size);
        shorten_count = 0;
        UINT_32 *new_array = new UINT_32[data.size];
        delete[] data.array;
        data.array = new_array;
    } //设定长度为2^n且不小于new_size
    void clear()
    {
        memset(data.array, 0, data.size * sizeof(*data.array));
    } //清空
    void quick_self_twice()
    {
        if (data.len)
        {
            REG UINT_32 tmp;
            tmp = data.array[data.len - 1];
            data.array[data.len - 1] = tmp << 1;
            if (data.size > data.len && (tmp >> (INT_BIT - 1)))
            {
                data.array[data.len] = INT32_0X01;
                data.len++;
            }
        }
    } //快速自增为二倍
    void quick_self_half()
    {
        if (data.len)
        {
            REG UINT_32 tmp = data.array[data.len - 1];
            data.array[data.len - 1] = tmp >> 1;
            if (data.len > 1)
            {
                data.array[data.len - 2] = tmp << (INT_BIT - 1);
            }
        }
    } //快速自减为一半
    HyperInt quick_twice() const
    {
        HyperInt result(*this);
        result.quick_self_twice();
        result.data.len = result.true_len();
        return result;
    }
    HyperInt quick_half() const
    {
        HyperInt result(*this);
        result.quick_self_half();
        result.data.len = result.true_len();
        return result;
    }
    HyperInt add(const HyperInt &input) const
    {
        REG HyperInt result;
        REG size_t len1 = true_len(), len2 = input.true_len();
        REG size_t result_len = _MAX_(data.len, input.data.len) + 1;
        result.reset_size(result_len);
        result.data.len = result_len;
        REG INT_64 tmp = 0;
        REG size_t pos1 = 0, pos2 = 0, count = 0;
        while (count < result_len)
        {
            if (pos1 < len1)
            {
                tmp += data.array[pos1];
                pos1++;
            }
            if (pos2 < len2)
            {
                tmp += input.data.array[pos2];
                pos2++;
            }
            result.data.array[count] = tmp;
            tmp >>= INT_BIT;
            count++;
        }
        result.data.len = result.true_len();
        return result;
    } //基础加法
    HyperInt sub(const HyperInt &input) const
    {
        REG HyperInt result;
        REG size_t len1 = true_len(), len2 = input.true_len();
        REG size_t result_len = _MAX_(data.len, input.data.len) + 1;
        result.reset_size(result_len);
        result.data.len = result_len;
        REG INT_64 tmp = 0;
        REG size_t pos1 = 0, pos2 = 0, count = 0;
        while (count < result_len)
        {
            if (pos1 < len1)
            {
                tmp += data.array[pos1];
                pos1++;
            }
            if (pos2 < len2)
            {
                tmp -= input.data.array[pos2];
                pos2++;
            }
            result.data.array[count] = tmp;
            tmp >>= INT_BIT;
            count++;
        }
        result.data.len = result.true_len();
        return result;
    } //基础减法
    HyperInt multiply(const HyperInt &input) const
    {
        REG HyperInt result;
        if (abs_equal(0) || input.abs_equal(0))
        {
            return result;
        }
        REG size_t len1 = true_len(), len2 = input.true_len();
        REG size_t result_len = data.len + input.data.len;
        result.reset_size(result_len);
        result.clear();
        result.data.len = result_len;
        result.data.neg = data.neg ^ input.data.neg;
        REG UINT_64 tmp, sum;
        for (REG size_t pos1 = 0; pos1 < len1; pos1++)
        {
            for (REG size_t pos2 = 0; pos2 < len2; pos2++)
            {
                tmp = (UINT_64)data.array[pos1] * input.data.array[pos2];
                for (REG size_t pos3 = pos1 + pos2; pos3 < result_len; pos3++)
                {
                    sum = tmp + result.data.array[pos3];
                    result.data.array[pos3] = sum;
                    if ((sum >> INT_BIT) > 0)
                    {
                        tmp = sum >> INT_BIT;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        result.data.len = result.true_len();
        return result;
    } //普通乘法
    HyperInt fft_multiply(HyperInt input) const
    {
        REG HyperInt result;
        if (abs_equal(0) || input.abs_equal(0))
        {
            return result;
        }
        REG size_t len1 = true_len(), len2 = input.true_len(), out_len = len1 + len2, fft_len = 1;
        result.reset_size(out_len);
        result.data.neg = data.neg ^ input.data.neg;
        while (fft_len < _TWICE_(out_len))
        {
            _SELFTWICE_(fft_len);
        }
        hint::Complex *fft_in1 = new hint::Complex[fft_len];
        hint::Complex *fft_in2 = new hint::Complex[fft_len];
        hint::Complex *fft_out = new hint::Complex[fft_len];
        REG UINT_32 ary_tmp;
        REG size_t half_pos, pos = 0;
        //每一位分解为短整数后存入复数组
        while (pos < fft_len)
        {
            half_pos = _HALF_(pos);
            if (half_pos < len1)
            {
                ary_tmp = data.array[half_pos];
                fft_in1[pos].real = ary_tmp & INT16_0XFF;
                fft_in1[pos + 1].real = ary_tmp >> SHORT_BIT;
            }
            if (half_pos < len2)
            {
                ary_tmp = input.data.array[half_pos];
                fft_in2[pos].real = ary_tmp & INT16_0XFF;
                fft_in2[pos + 1].real = ary_tmp >> SHORT_BIT;
            }
            if (half_pos >= _MAX_(len1, len2))
            {
                break;
            }
            pos += 2;
        }
        fft(fft_in1, fft_len, false); //快速傅里叶变换
        fft(fft_in2, fft_len, false);
        //#pragma omp parallel for
        //每一位相乘
        for (REG size_t i = 0; i < fft_len; i += 2)
        {
            fft_out[i] = fft_in1[i] * fft_in2[i];
            fft_out[i + 1] = fft_in1[i + 1] * fft_in2[i + 1];
        }
        fft(fft_out, fft_len, true); //逆变换
        REG size_t fft_tlen = fft_len;
        REG UINT_64 tmp = 0;
        while (!(UINT_8)(fft_out[fft_tlen - 1].real + 0.5))
        {
            fft_tlen--;
        }
        pos = 0;
        //整理每一位
        while (half_pos < out_len)
        {
            half_pos = _HALF_(pos);
            if (pos < fft_tlen)
            {
                tmp += fft_out[pos].real + 0.5;
            }
            if (_ODD_(pos))
            {
                ary_tmp += (UINT_32)tmp << SHORT_BIT;
                result.data.array[half_pos] = ary_tmp;
            }
            else
            {
                ary_tmp = tmp & INT16_0XFF;
            }
            tmp >>= SHORT_BIT;
            pos++;
        }
        result.data.len = half_pos;
        delete[] fft_in1;
        delete[] fft_in2;
        delete[] fft_out;
        return result;
    } //快速傅里叶变换乘法
public:
    void test()
    {
        while (1)
        {
            console_out_hex();
            getchar();
            printf("\n");
            quick_self_twice();
        }
    }
    ~HyperInt()
    {
        if (data.array != nullptr)
        {
            delete[] data.array;
            data.array = nullptr;
        }
    } //析构函数
    HyperInt()
    {
        data.len = 0;
        data.size = 2;
        data.neg = false;
        shorten_count = 0;
        data.array = new UINT_32[data.size];
        data.array[1] = data.array[0] = 0;
    } //无参数构造
    HyperInt(const HyperInt &input)
    {
        if (this != &input)
        {
            data.len = input.true_len();
            data.size = 2;
            while (data.size < data.len)
            {
                _SELFTWICE_(data.size);
            }
            data.neg = input.data.neg;
            shorten_count = 0;
            data.array = new UINT_32[data.size];
            _ARYCOPY_(data.array, input.data.array, data.len);
        }
    } // HyperInt拷贝构造
    HyperInt(HyperInt &&input)
    {
        if (this != &input)
        {
            data.len = input.true_len();
            data.size = input.data.size;
            shorten_count = 0;
            data.neg = input.data.neg;
            data.array = input.data.array;
            input.data.array = nullptr;
        }
    } // HyperInt移动构造
    HyperInt(const hint::h_int &input)
    {
        data.len = input.len;
        data.size = 2;
        while (data.size < data.len)
        {
            _SELFTWICE_(data.size);
        }
        data.neg = input.neg;
        shorten_count = 0;
        data.array = new UINT_32[data.size];
        _ARYCOPY_(data.array, input.array, data.len);
        data.len = true_len();
    } // h_int参数拷贝构造
    HyperInt(INT_64 input)
    {
        data.len = data.size = 2;
        data.neg = _NEG_(input);
        input = _ABS_(input);
        shorten_count = 0;
        data.array = new UINT_32[data.size];
        data.array[0] = input & INT32_0XFF;
        data.array[1] = input >> INT_BIT;
        data.len = true_len();
        //printf("_5_\n");
    } // INT_64参数构造
    HyperInt(const std::string input)
    {
        string_in(input);
    } // string参数构造
    HyperInt operator=(const HyperInt &input)
    {
        if (this != &input)
        {
            data.len = input.true_len();
            if (data.size < input.data.len || shorten_count >= size_shorten_max)
            {
                data.size = 2;
                while (data.size < input.data.len)
                {
                    _SELFTWICE_(data.size);
                }
                delete[] data.array;
                data.array = new UINT_32[data.size];
                if (shorten_count >= size_shorten_max)
                {
                    shorten_count = 0;
                }
            }
            else if (_HALF_(data.size) > data.len)
            {
                shorten_count++;
            }
            data.neg = input.data.neg;
            _ARYCOPY_(data.array, input.data.array, data.len);
        }
        return *this;
    } // HyperInt拷贝赋值
    HyperInt operator=(HyperInt &&input)
    {
        if (this != &input)
        {
            data.len = input.true_len();
            data.size = input.data.size;
            data.neg = input.data.neg;
            delete[] data.array;
            data.array = input.data.array;
            input.data.array = nullptr;
        }
        return *this;
    } // HyperInt移动赋值
    HyperInt operator=(const hint::h_int &input)
    {
        data.len = true_len();
        if (data.size < input.len || shorten_count >= size_shorten_max)
        {
            data.size = 2;
            while (data.size < input.len)
            {
                _SELFTWICE_(data.size);
            }
            delete[] data.array;
            data.array = new UINT_32[data.size];
            if (shorten_count >= size_shorten_max)
            {
                shorten_count = 0;
            }
        }
        else if ((data.size >> 1) > data.len)
        {
            shorten_count++;
        }
        data.neg = input.neg;
        _ARYCOPY_(data.array, input.array, data.len);
        return *this;
    } // h_int拷贝赋值
    HyperInt operator=(hint::h_int &&input)
    {
        data.len = _MAX_(input.len, 2);
        data.size = input.size;
        data.neg = input.neg;
        data.array = input.array;
        input.array = nullptr;
        data.len = true_len();
        return *this;
    } // h_int移动赋值
    HyperInt operator=(INT_64 input)
    {
        if (data.len > 2)
        {
            shorten_count++;
        }
        if (shorten_count > size_shorten_max)
        {
            shorten_count = 0;
            data.size = 2;
            delete[] data.array;
            data.array = new UINT_32[data.size];
        }
        data.len = 2;
        data.neg = _NEG_(input);
        input = _ABS_(input);
        shorten_count = 0;
        data.array[0] = input & INT32_0XFF;
        data.array[1] = input >> INT_BIT;
        data.len = true_len();
        return *this;
    } // longlong赋值
    HyperInt operator=(const std::string input)
    {
        string_in(input);
        return *this;
    } // string赋值

    //基本操作
    bool is_neg() const
    {
        return data.neg;
    }
    size_t length() const
    {
        return data.len;
    }
    size_t size() const
    {
        return data.size;
    }
    INT_64 to_int64() const
    {
        INT_64 out = 0;
        out = (data.array[1] & INT32_0X7F);
        out <<= INT_BIT;
        out += data.array[0];
        if (data.neg)
            out = -out;
        return out;
    }
    std::string to_string() const
    {
        REG std::string result;
        if (abs_equal(0))
            return std::string("0");
        REG HyperInt tmp1(*this), tmp2;
        REG size_t pos = 0;
        char *tmp = new char[10 * data.len];
        while (tmp1 > 0)
        {
            tmp2 = tmp1 / 10;
            tmp[pos] = (tmp1 - tmp2 * 10).to_int64() + '0';
            tmp1 = std::move(tmp2);
            pos++;
        }
        result[pos] = 0;
        while (pos > 0)
        {
            result += tmp[pos - 1];
            pos--;
        }
        return result;
    }
    void console_in()
    {
        clear();
        char tmp = '0';
        while (tmp != '\n' && '0' <= tmp <= '9')
        {
            if (tmp != '\n')
            {
                *this *= 10;
                *this += ((INT_64)(tmp - '0'));
            }
            tmp = getchar();
        }
    }
    void string_in(std::string s)
    {
        clear();
        char tmp = '0';
        for (std::string::iterator i = s.begin(); i != s.end(); i++)
        {
            if ('0' <= *i <= '9')
            {
                *this *= 10;
                *this += ((INT_64)(*i - '0'));
            }
        }
    }
    void console_out_hex() const
    {
        if (data.neg)
        {
            printf("-");
        }
        REG size_t pos = data.len;
        while (pos)
        {
            pos--;
            printf("%X ", data.array[pos]);
        }
    }
    void console_out_dec() const
    {
        REG HyperInt tmp1(*this), tmp2;
        REG size_t pos = 0;
        char *result = new char[10 * data.len];
        while (tmp1 > 0)
        {
            tmp2 = tmp1 / 10;
            result[pos] = (tmp1 - tmp2 * 10).to_int64() + '0';
            tmp1 = std::move(tmp2);
            pos++;
        }
        result[pos] = 0;
        while (pos > 0)
        {
            putchar(result[pos - 1]);
            pos--;
        }
    }
    HyperInt &self_half()
    {
        REG UINT_32 tmp1, tmp2 = 0;
        REG size_t pos = data.len;
        while (pos)
        {
            pos--;
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (tmp1 << (INT_BIT - 1)) + (tmp2 >> 1);
            data.array[pos] = tmp1;
        }
        data.len = true_len();
        return *this;
    } // /2
    HyperInt &self_twice()
    {
        REG UINT_32 tmp1, tmp2 = 0;
        for (REG size_t pos = 0; pos < data.len; pos++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (tmp1 >> (INT_BIT - 1)) + (tmp2 << 1);
            data.array[pos] = tmp1;
        }
        if (data.size > data.len)
        {
            data.array[data.len] = tmp2 >> (INT_BIT - 1);
            data.len++;
        }
        return *this;
    } // *2
    HyperInt half() const
    {
        return r_shift(1);
    } // /2
    HyperInt twice() const
    {
        return l_shift(1);
    } // *2
    HyperInt r_shift(size_t n) const
    {
        if (!n)
            return *this;
        HyperInt result;
        REG size_t shift = n >> LOG_INT_BIT;
        REG size_t offset = n % 32;
        REG size_t len = data.len - shift;
        REG UINT_32 tmp1, tmp2 = 0;
        result.reset_size(len);
        result.data.len = len;
        while (len > 0)
        {
            len--;
            tmp1 = tmp2;
            tmp2 = data.array[len + shift];
            tmp1 = (offset ? (tmp1 << (INT_BIT - offset)) : 0) + (tmp2 >> offset);
            result.data.array[len] = tmp1;
        }
        result.data.len = result.true_len();
        return result;
    }
    HyperInt l_shift(size_t n) const
    {
        if (!n)
            return *this;
        HyperInt result;
        REG size_t shift = n >> LOG_INT_BIT;
        REG size_t offset = n % 32;
        REG UINT_32 tmp1, tmp2 = 0;
        result.reset_size(data.len + shift + 1);
        result.clear();
        result.data.len = data.len + shift + 1;
        for (REG size_t pos = 0; pos < data.len; pos++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (offset ? (tmp1 >> (INT_BIT - offset)) : 0) + (tmp2 << offset);
            result.data.array[pos + shift] = tmp1;
        }
        if (tmp2 >> (INT_BIT - offset) & offset)
            result.data.array[data.len + shift] = tmp2 >> (INT_BIT - offset);
        else
            result.data.len--;
        return result;
    }
    HyperInt abs() const
    {
        HyperInt result(*this);
        result.data.neg = false;
        return result;
    }
    bool abs_larger(const HyperInt &input) const
    {
        REG size_t t_len1 = true_len(), t_len2 = input.true_len();
        if (t_len1 > t_len2)
            return true;
        else if (t_len1 < t_len2)
            return false;
        else
        {
            REG UINT_32 num1, num2;
            while (t_len1 > 0)
            {
                t_len1--;
                num1 = data.array[t_len1];
                num2 = input.data.array[t_len1];
                if (num1 > num2)
                {
                    return true;
                }
                else if (num1 < num2)
                {
                    return false;
                }
            }
            return false;
        }
    }
    bool abs_smaller(const HyperInt &input) const
    {
        REG size_t t_len1 = true_len(), t_len2 = input.true_len();
        if (t_len1 < t_len2)
            return true;
        else if (t_len1 > t_len2)
            return false;
        else
        {
            REG UINT_32 num1, num2;
            while (t_len1 > 0)
            {
                t_len1--;
                num1 = data.array[t_len1];
                num2 = input.data.array[t_len1];
                if (num1 < num2)
                {
                    return true;
                }
                else if (num1 > num2)
                {
                    return false;
                }
            }
            return false;
        }
    }
    bool abs_equal(const HyperInt &input) const
    {
        REG size_t t_len1 = true_len(), t_len2 = input.true_len();
        if (t_len1 != t_len2)
            return false;
        else
        {
            while (t_len1 > 0)
            {
                t_len1--;
                if (data.array[t_len1] != input.data.array[t_len1])
                    return false;
            }
            return true;
        }
    } //比较是否为相反数

    //逻辑运算
    bool operator==(const HyperInt &input) const;
    bool operator==(const INT_64 &input) const;
    bool operator>(const HyperInt &input) const;
    bool operator>(const INT_64 &input) const;
    bool operator>=(const HyperInt &input) const;
    bool operator>=(const INT_64 &input) const;
    bool operator<(const HyperInt &input) const;
    bool operator<(const INT_64 &input) const;
    bool operator<=(const HyperInt &input) const;
    bool operator<=(const INT_64 &input) const;
    // bool operator!();
    bool operator!() const;
    // operator bool();

    //友元函数
    friend HyperInt add(hint::h_int in1, hint::h_int in2);
    friend HyperInt quick_multiply(hint::h_int in1, hint::h_int in2);
    friend HyperInt abs(HyperInt input);
    friend void print(const HyperInt &input);
    friend bool operator==(const INT_64 &input1, const HyperInt &input2);
    friend bool operator>(const INT_64 &input1, const HyperInt &input2);
    friend bool operator>=(const INT_64 &input1, const HyperInt &input2);
    friend bool operator<(const INT_64 &input1, const HyperInt &input2);
    friend bool operator<=(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator+(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator-(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator*(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator/(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator%(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator+=(INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator-=(INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator*=(INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator/=(INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator%=(INT_64 &input1, const HyperInt &input2);
    friend std::ostream &operator<<(std::ostream &output, HyperInt &input);
    friend std::istream &operator>>(std::istream &input, HyperInt &output);

    //算术运算
    HyperInt operator+(const HyperInt &input) const;
    HyperInt operator+() const;
    HyperInt operator-(const HyperInt &input) const;
    HyperInt operator-() const;
    HyperInt operator*(const HyperInt &input) const;
    HyperInt operator*(INT_64 input) const;
    HyperInt operator/(const HyperInt &input) const;
    HyperInt operator%(const HyperInt &input) const;
    HyperInt operator+=(const HyperInt &input);
    HyperInt operator-=(const HyperInt &input);
    HyperInt operator*=(const HyperInt &input);
    HyperInt operator*=(INT_64 input);
    HyperInt operator/=(const HyperInt &input);
    HyperInt operator%=(const HyperInt &input);
    HyperInt operator++(int);
    HyperInt &operator++();
    HyperInt operator--(int);
    HyperInt &operator--();
    HyperInt power(UINT_32 n) const;
};
//逻辑运算
inline bool HyperInt::operator==(const HyperInt &input) const
{
    if (data.neg ^ input.data.neg)
        return false;
    else
        return abs_equal(input);
}
inline bool HyperInt::operator==(const INT_64 &input) const
{
    if (data.neg ^ _NEG_(input))
        return false;
    else
        return abs_equal(input);
}
bool HyperInt::operator>(const HyperInt &input) const
{
    if (data.neg ^ input.data.neg)
        return input.data.neg;
    else
        return data.neg ^ abs_larger(input);
}
inline bool HyperInt::operator>(const INT_64 &input) const
{
    return *this > HyperInt(input);
}
inline bool HyperInt::operator>=(const HyperInt &input) const
{
    return !(*this < input);
}
inline bool HyperInt::operator>=(const INT_64 &input) const
{
    return !(*this < HyperInt(input));
}
bool HyperInt::operator<(const HyperInt &input) const
{
    if (data.neg ^ input.data.neg)
        return data.neg;
    else
        return data.neg ^ abs_smaller(input);
}
inline bool HyperInt::operator<(const INT_64 &input) const
{
    return *this < HyperInt(input);
}
inline bool HyperInt::operator<=(const HyperInt &input) const
{
    return !(*this > input);
}
inline bool HyperInt::operator<=(const INT_64 &input) const
{
    return !(*this > HyperInt(input));
}
inline bool HyperInt::operator!() const
{
    return !(bool)true_len();
}
// inline HyperInt::operator bool()
// {
//     if(data.len>2)
//         return true;
//     return(data.array[0]!=0||data.array[1]!=0);
// }

//友元函数
HyperInt abs(HyperInt input)
{
    HyperInt result(input);
    result.data.neg = false;
    return result;
}
void print(HyperInt &input)
{
    input.console_out_dec();
}
bool operator==(const INT_64 &input1, const HyperInt &input2)
{
    return input2 == input1;
}
bool operator>(const INT_64 &input1, const HyperInt &input2)
{
    return input2 > input1;
}
bool operator>=(const INT_64 &input1, const HyperInt &input2)
{
    return input2 >= input1;
}
bool operator<(const INT_64 &input1, const HyperInt &input2)
{
    return input2 < input1;
}
bool operator<=(const INT_64 &input1, const HyperInt &input2)
{
    return input2 <= input1;
}
HyperInt operator+(const INT_64 &input1, const HyperInt &input2)
{
    return input2 + input1;
}
HyperInt operator-(const INT_64 &input1, const HyperInt &input2)
{
    return input2 - input1;
}
HyperInt operator*(const INT_64 &input1, const HyperInt &input2)
{
    return input2 * input1;
}
HyperInt operator/(const INT_64 &input1, const HyperInt &input2)
{
    return input2 / input1;
}
HyperInt operator%(const INT_64 &input1, const HyperInt &input2)
{
    return input2 % input1;
}
HyperInt operator+=(INT_64 &input1, const HyperInt &input2)
{
    return input1 += input2.to_int64();
}
HyperInt operator-=(INT_64 &input1, const HyperInt &input2)
{
    return input1 -= input2.to_int64();
}
HyperInt operator*=(INT_64 &input1, const HyperInt &input2)
{
    return input1 *= input2.to_int64();
}
HyperInt operator/=(INT_64 &input1, const HyperInt &input2)
{
    return input1 /= input2.to_int64();
}
HyperInt operator%=(INT_64 &input1, const HyperInt &input2)
{
    return input1 %= input2.to_int64();
}
std::ostream &operator<<(std::ostream &output, HyperInt &input)
{
    output << input.to_string();
    return output;
}
std::istream &operator>>(std::istream &input, HyperInt &output)
{
    std::string in;
    input >> in;
    output.string_in(in);
    return input;
}

//算术运算
HyperInt HyperInt::operator+(const HyperInt &input) const
{
    REG HyperInt result;
    if (!data.neg ^ input.data.neg)
    {
        result = add(input);
        result.data.neg = data.neg;
        return result;
    }
    else
    {
        if (abs_equal(input))
        {
            return result;
        }
        else if (abs_larger(input))
        {
            result = sub(input);
            result.data.neg = data.neg;
            return result;
        }
        else
        {
            result = input.sub(*this);
            result.data.neg = input.data.neg;
            return result;
        }
    }
}
inline HyperInt HyperInt::operator+() const
{
    return *this;
}
HyperInt HyperInt::operator-(const HyperInt &input) const
{
    REG HyperInt result;
    if (data.neg ^ input.data.neg)
    {
        result = add(input);
        result.data.neg = data.neg;
        return result;
    }
    else
    {
        if (abs_equal(input))
        {
            return result;
        }
        else if (abs_larger(input))
        {
            result = sub(input);
            result.data.neg = data.neg;
            return result;
        }
        else
        {
            result = input.sub(*this);
            result.data.neg = input.data.neg;
            return result;
        }
    }
}
inline HyperInt HyperInt::operator-() const
{
    HyperInt out(*this);
    out.data.neg = !data.neg;
    return out;
}
HyperInt HyperInt::operator*(const HyperInt &input) const
{
    REG size_t max_len = _MAX_(data.len, input.data.len);
    if (max_len < FFT_MIN || max_len > FFT_MAX)
    {
        return multiply(input);
    }
    else
    {
        return fft_multiply(input);
    }
}
HyperInt HyperInt::operator*(INT_64 input) const
{
    REG HyperInt result;
    if (abs_equal(0) || input == 0)
    {
        return result;
    }
    REG size_t len = true_len();
    REG size_t result_len = len + 2;
    result.reset_size(result_len);
    result.data.len = result_len;
    result.data.neg = data.neg ^ _NEG_(input);
    input = _ABS_(input);
    REG UINT_64 tmp1, tmp2, sum, input_num1, input_num2;
    input_num1 = input & INT32_0XFF;
    input_num2 = input >> INT_BIT;
    tmp2 = data.array[0];
    sum = tmp2 * input_num1;
    result.data.array[0] = sum & INT32_0XFF;
    sum >>= INT_BIT;
    for (REG size_t pos = 1; pos < result_len; pos++)
    {
        tmp1 = tmp2;
        if (pos < len)
        {
            tmp2 = data.array[pos];
        }
        else
        {
            tmp2 = 0;
        }
        sum += input_num2 * tmp1 + input_num1 * tmp2;
        result.data.array[pos] = sum & INT32_0XFF;
        sum >>= INT_BIT;
    }
    result.data.len = result.true_len();
    return result;
}
HyperInt HyperInt::operator/(const HyperInt &input) const
{
    if (!input)
        throw "Can't divided by zero";
    REG HyperInt result, mid, tmp;
    if (abs_smaller(input) || abs_equal(0))
        return result;
    else
    {
        REG size_t len1 = true_len(), len2 = input.true_len(), i = 0;
        tmp.reset_size(len1 - len1 + 1);
        result.data.neg = data.neg ^ input.data.neg;
        tmp=1;
        while (!abs_smaller(input.l_shift(i)))
        {
            tmp.quick_self_twice();
            i++;
        }
        result = tmp.quick_half();
        mid = (result + tmp).half();
        while (!result.abs_equal(mid))
        {
            if (abs_smaller(mid.multiply(input)))
            {
                tmp = mid;
            }
            else
            {
                result = mid;
            }
            mid = (result + tmp).half();
        }
        result.data.len = result.true_len();
        return result;
    }
}
HyperInt HyperInt::operator%(const HyperInt &input) const
{
    HyperInt tmp = *this / input;
    tmp *= input;
    return *this - tmp;
}
HyperInt HyperInt::operator+=(const HyperInt &input)
{
    *this = *this + input;
    return *this;
}
HyperInt HyperInt::operator-=(const HyperInt &input)
{
    *this = *this - input;
    return *this;
}
HyperInt HyperInt::operator*=(const HyperInt &input)
{
    *this = *this * input;
    return *this;
}
HyperInt HyperInt::operator*=(INT_64 input)
{
    *this = *this * input;
    return *this;
}
HyperInt HyperInt::operator/=(const HyperInt &input)
{
    *this = *this / input;
    return *this;
}
HyperInt HyperInt::operator%=(const HyperInt &input)
{
    *this = *this % input;
    return *this;
}
HyperInt HyperInt::operator++(int)
{
    HyperInt tmp(*this);
    *this = *this + 1;
    return tmp;
}
HyperInt &HyperInt::operator++()
{
    *this = *this + 1;
    return *this;
}
HyperInt HyperInt::operator--(int)
{
    HyperInt tmp(*this);
    *this = *this - 1;
    return tmp;
}
HyperInt &HyperInt::operator--()
{
    *this = *this - 1;
    return *this;
}
HyperInt HyperInt::power(UINT_32 n) const
{
    REG HyperInt tmp(*this), result = 1;
    while (n)
    {
        if (n & 1)
        {
            result = result * tmp;
        }
        tmp = tmp * tmp;
        n >>= 1;
    }
    return result;
}
#endif
