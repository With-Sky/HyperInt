#ifndef HINT_HPP
#define HINT_HPP

#include <iostream>
#include <future>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdlib>

#define UINT_8 unsigned char
#define UINT_16 unsigned short
#define UINT_32 unsigned int
#define UINT_64 unsigned long long
#define INT_32 int
#define INT_64 long long

#define size_t_bits 64 //如果size_t不是64位需要修改程序
#define SHORT_BIT 16
#define INT_BIT 32
#define LOG_INT_BIT 5
#define INT16_0XFF UINT16_MAX
#define INT32_0XFF UINT32_MAX
#define INT32_0X01 (UINT32)1
#define INT32_0X80 INT32_MIN
#define INT32_0X7F INT32_MAX
#define INT64_0x80 INT64_MIN
#define INT64_0x7F INT64_MAX
#define FFT_MIN 128
#define FFT_MAX 512000
#define PI 3.1415926535897932

#define _MAX_(x, y) ((x) > (y) ? (x) : (y))
#define _MIN_(x, y) ((x) < (y) ? (x) : (y))
#define _NEG_(x) ((x) < 0)
#define _ODD_(x) ((x)&1)
#define _ABS_(x) ((x) < 0 ? (-(x)) : (x))
#define _TWICE_(x) ((x) << 1)
#define _HALF_(x) ((x) >> 1)
#define _SELFTWICE_(x) ((x) <<= 1)
#define _SELFHALF_(x) ((x) >>= 1)
#define _INT32REARY_(ptr, len) ptr = (UINT_32 *)realloc((ptr), (len) * sizeof(UINT_32))
#define _ARYCOPY_(target, source, len) std::memcpy((target), (source), (len) * sizeof(*target))

namespace hint
{
    struct h_int
    {
        UINT_32 *array = nullptr;
        INT_64 neg_n_len = 0;
        size_t size = 0;
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
        size_t log_n = 1;
        Complex tmp, unit_omega, omega, tmp1, tmp2;
        while ((1 << log_n) < n)
        {
            log_n++;
        }
        rev[0] = 0;
        for (size_t i = 1; i < n; i++)
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log_n - 1));
        for (size_t i = 0; i < n; i++)
        {
            if (i < rev[i])
            {
                tmp = input[i];
                input[i] = input[rev[i]];
                input[rev[i]] = tmp;
            }
        }
        delete[] rev;
        for (size_t rank = 1, gap; rank < n; _SELFTWICE_(rank))
        {
            gap = rank << 1;
            unit_omega = Complex(gap);
            if (is_ifft)
            {
                unit_omega.imaginary = -unit_omega.imaginary;
            }
            for (size_t begin = 0; begin < n; begin += gap)
            {
                omega = Complex(1, 0);
                for (size_t pos = begin; pos < begin + rank; pos++)
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
            for (size_t i = 0; i < n; i++)
                input[i].real /= 1.0 * n;
        }
    } //快速傅里叶变换
    // HyperInt factori(UINT_32 input);
}
class HyperInt
{
private:
    hint::h_int data;
    void change_length(size_t new_length)
    {
        data.neg_n_len = (data.neg_n_len & INT32_0X80) | new_length;
    } //设置新的长度
    void set_true_len()
    {
        size_t t_len = length();
        while (t_len > 0 && data.array[t_len - 1] == 0)
        {
            t_len--;
        }
        change_length(t_len);
    }
    void clear()
    {
        memset(data.array, 0, data.size * sizeof(*data.array));
    } //清空
    void quick_self_twice()
    {
        size_t len_tmp = length();
        if (len_tmp > 0)
        {
            UINT_32 tmp = data.array[len_tmp - 1];
            _SELFTWICE_(tmp);
            if (tmp == 0)
            {
                len_tmp++;
                data.array[len_tmp - 1] = 1;
            }
            else
            {
                data.array[len_tmp - 1] = tmp;
            }
        }
    } //快速自增为二倍
    void quick_self_half()
    {
        size_t len_tmp = length();
        if (len_tmp > 0)
        {
            UINT_32 tmp = data.array[len_tmp - 1];
            _SELFHALF_(tmp);
            if (tmp == 0)
            {
                len_tmp--;
                data.array[len_tmp - 1] = INT32_0X80;
            }
            else
            {
                data.array[len_tmp - 1] = tmp;
            }
        }
    } //快速自减为一半
    HyperInt quick_twice()
    {
        HyperInt result(*this);
        result.quick_self_twice();
        result.set_true_len();
        return result;
    }
    HyperInt quick_half()
    {
        HyperInt result(*this);
        result.quick_self_half();
        result.set_true_len();
        return result;
    }
    HyperInt add(const HyperInt &input) const
    {
        HyperInt result;
        size_t len1 = length(), len2 = input.length();
        size_t result_len = _MAX_(len1, len2) + 1;
        result.reset_size(result_len);
        result.change_length(result_len);
        INT_64 tmp = 0;
        size_t pos1 = 0, pos2 = 0, count = 0;
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
        result.set_true_len();
        return result;
    } //基础加法
    HyperInt sub(const HyperInt &input) const
    {
        HyperInt result;
        size_t len1 = length(), len2 = input.length();
        size_t result_len = _MAX_(len1, len2);
        result.reset_size(result_len);
        result.change_length(result_len);
        INT_64 tmp = 0;
        size_t pos1 = 0, pos2 = 0, count = 0;
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
        result.set_true_len();
        return result;
    } //基础减法
    HyperInt multiply(const HyperInt &input) const
    {
        HyperInt result;
        if (equal_to_z() || input.equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        size_t result_len = len1 + len2;
        result.reset_size(result_len);
        result.clear();
        result.change_length(result_len);
        UINT_64 tmp = 0, sum = 0;
        for (size_t pos1 = 0; pos1 < len1; pos1++)
        {
            for (size_t pos2 = 0; pos2 < len2; pos2++)
            {
                tmp = (UINT_64)data.array[pos1] * input.data.array[pos2];
                for (size_t pos3 = pos1 + pos2; pos3 < result_len; pos3++)
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
        bool result_neg = is_neg() ^ input.is_neg();
        result.neg_sign(result_neg);
        result.set_true_len();
        return result;
    } //普通乘法
    HyperInt fft_multiply(HyperInt input) const
    {
        HyperInt result;
        if (equal_to_z() || input.equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        size_t out_len = len1 + len2, fft_len = 1;
        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);
        while (fft_len < _TWICE_(out_len))
        {
            _SELFTWICE_(fft_len);
        }
        hint::Complex *fft_in1 = new hint::Complex[fft_len];
        hint::Complex *fft_in2 = new hint::Complex[fft_len];
        hint::Complex *fft_out = new hint::Complex[fft_len];
        UINT_32 ary_tmp = 0;
        size_t half_pos = 0, pos = 0;
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
#pragma omp parallel for
        for (int i = 0; i < fft_len; i += 2)
        {
            fft_out[i] = fft_in1[i] * fft_in2[i];
            fft_out[i + 1] = fft_in1[i + 1] * fft_in2[i + 1];
        }                            //每一位相乘
        fft(fft_out, fft_len, true); //逆变换
        size_t fft_tlen = fft_len;
        UINT_64 tmp = 0;
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
        delete[] fft_in1;
        delete[] fft_in2;
        delete[] fft_out;
        result.set_true_len();
        bool result_neg = is_neg() ^ input.is_neg();
        result.neg_sign(result_neg);
        return result;
    } //快速傅里叶变换乘法
public:
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
        change_length(0);
        data.size = 2;
        neg_sign(false);
        if (data.array != nullptr)
            delete[] data.array;
        data.array = new UINT_32[data.size];
        data.array[1] = data.array[0] = 0;
    } //无参数构造
    HyperInt(const HyperInt &input)
    {
        if (this != &input)
        {
            change_length(input.length());
            data.size = 2;
            reset_size(length());
            neg_sign(input.is_neg());
            if (data.array != nullptr)
                delete[] data.array;
            data.array = new UINT_32[data.size];
            _ARYCOPY_(data.array, input.data.array, length());
            set_true_len();
        }
    } // HyperInt拷贝构造
    HyperInt(HyperInt &&input)
    {
        if (this != &input)
        {
            change_length(input.length());
            data.size = input.data.size;
            neg_sign(input.is_neg());
            delete[] data.array;
            data.array = input.data.array;
            input.data.array = nullptr;
            set_true_len();
        }
    } // HyperInt移动构造
    HyperInt(INT_64 input)
    {
        change_length(2);
        data.size = 2;
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        if (data.array != nullptr)
            delete[] data.array;
        data.array = new UINT_32[data.size];
        data.array[0] = input & INT32_0XFF;
        data.array[1] = input >> INT_BIT;
        set_true_len();
    } // INT_64参数构造
    HyperInt(const std::string input)
    {
        string_in(input);
    } // string参数构造
    HyperInt(const char input[])
    {
        char tmp = '0';
        for (size_t pos = 0; input[pos] != 0; pos++)
        {
            tmp = input[pos];
            if ('0' <= tmp && tmp <= '9')
            {
                *this *= 10;
                *this += ((INT_64)(tmp - '0'));
            }
        }
        if (input[0] == '-')
        {
            neg_sign(true);
        }
    }
    HyperInt &operator=(const HyperInt &input)
    {
        if (this != &input)
        {
            change_length(input.length());
            data.size = 2;
            reset_size(length());
            neg_sign(input.is_neg());
            if (data.array != nullptr)
                delete[] data.array;
            data.array = new UINT_32[data.size];
            _ARYCOPY_(data.array, input.data.array, length());
            set_true_len();
        }
        return *this;
    } // HyperInt拷贝赋值
    HyperInt &operator=(HyperInt &&input)
    {
        if (this != &input)
        {
            change_length(input.length());
            data.size = input.data.size;
            neg_sign(input.is_neg());
            delete[] data.array;
            data.array = input.data.array;
            input.data.array = nullptr;
            set_true_len();
        }
        return *this;
    } // HyperInt移动赋值
    HyperInt &operator=(INT_64 input)
    {
        change_length(2);
        data.size = 2;
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        if (data.array != nullptr)
            delete[] data.array;
        data.array = new UINT_32[data.size];
        data.array[0] = input & INT32_0XFF;
        data.array[1] = input >> INT_BIT;
        set_true_len();
        return *this;
    } // longlong赋值
    HyperInt &operator=(const std::string input)
    {
        string_in(input);
        return *this;
    } // string赋值

    //基本操作
    void neg_sign(bool neg)
    {
        if (neg)
        {
            data.neg_n_len = data.neg_n_len | INT64_0x80;
        }
        else
        {
            data.neg_n_len = data.neg_n_len & INT64_0x7F;
        }
    }
    INT_64 div_mod(INT_64 divisor)
    {
        //  lldiv_t div_tmp;
        UINT_64 last_rem = 0, tmp = 0, rem_num = 0;
        bool result_neg = is_neg() ^ _NEG_(divisor);
        neg_sign(result_neg);
        size_t pos = length();
        while (pos > 1)
        {
            pos--;
            tmp = (last_rem << INT_BIT) + data.array[pos];
            auto div_tmp = lldiv(tmp, divisor);      //一次性得到商和余数
            data.array[pos] = (UINT_32)div_tmp.quot; //当前数为变商
            last_rem = div_tmp.rem;                  //得到余数
        }
        tmp = (last_rem << INT_BIT) + data.array[0];
        auto div_tmp = lldiv(tmp, divisor);
        data.array[0] = (UINT_32)div_tmp.quot;
        rem_num = div_tmp.rem;
        set_true_len();
        return rem_num;
    }
    INT_64 mod(INT_64 divisor) const
    {
        //  lldiv_t div_tmp;
        UINT_64 last_rem = 0, tmp = 0, rem_num = 0;
        bool result_neg = is_neg() ^ _NEG_(divisor);
        size_t pos = length();
        while (pos > 1)
        {
            pos--;
            tmp = (last_rem << INT_BIT) + data.array[pos];
            auto div_tmp = lldiv(tmp, divisor); //一次性得到商和余数
            last_rem = div_tmp.rem;             //得到余数
        }
        tmp = (last_rem << INT_BIT) + data.array[0];
        auto div_tmp = lldiv(tmp, divisor);
        rem_num = div_tmp.rem;
        if (is_neg())
        {
            rem_num = -rem_num;
        }
        return rem_num;
    }
    HyperInt power(UINT_64 n) const
    {
        HyperInt tmp(*this), result = 1;
        if (!_ODD_(n))
        {
            result.neg_sign(false);
        }
        while (n)
        {
            if (n & 1)
            {
                result = result * tmp;
            }
            tmp = tmp * tmp;
            _SELFHALF_(n);
        }
        return result;
    }
    bool is_neg() const
    {
        return data.neg_n_len & INT64_0x80;
    }
    size_t length() const
    {
        return data.neg_n_len & INT64_0x7F;
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
        if (is_neg())
            out = -out;
        return out;
    }
    std::string to_string() const
    {
        if (equal_to_z())
            return std::string("0");
        std::string result;
        HyperInt tmp(*this);
        size_t pos = 0;
        INT_64 rem_n = 0;
        while (tmp.abs_larger(0ll))
        {
            rem_n = tmp.div_mod(10);
            result += rem_n + '0';
            pos++;
        }
        if (is_neg())
        {
            result += '-';
        }
        reverse(result.begin(), result.end());
        return result;
    }
    void console_in()
    {
        clear();
        char tmp = '0';
        bool head = true, neg = false;
        while (tmp != '\n' && '0' <= tmp && tmp <= '9')
        {
            if (tmp != '\n')
            {
                *this *= 10;
                *this += ((INT_64)(tmp - '0'));
            }
            tmp = getchar();
            if (head)
            {
                if (tmp == '-')
                    neg = true;
                head = false;
            }
        }
        neg_sign(neg);
    }
    void string_in(std::string s)
    {
        clear();
        char tmp = '0';
        for (std::string::iterator i = s.begin(); i != s.end(); i++)
        {
            if ('0' <= *i && *i <= '9')
            {
                *this *= 10;
                *this += ((INT_64)(*i - '0'));
            }
        }
        if (s[0] == '-')
        {
            neg_sign(true);
        }
    }
    void console_out_hex() const
    {
        if (equal_to_z())
            printf("0");
        else
        {
            if (is_neg())
            {
                printf("-");
            }
            size_t pos = length();
            while (pos)
            {
                pos--;
                printf("%X ", data.array[pos]);
            }
        }
    }
    void console_out_dec() const
    {
        if (equal_to_z())
        {
            putchar('0');
            return;
        }
        HyperInt tmp(*this);
        size_t pos = 0;
        UINT_8 *tmp_str = new UINT_8[5 * length()];
        INT_64 rem_n = 0;
        while (tmp.abs_larger(0ll))
        {
            rem_n = tmp.div_mod(100);
            tmp_str[pos] = rem_n;
            pos++;
        }
        if (is_neg())
        {
            putchar('-');
        }
        pos--;
        printf("%d", tmp_str[pos]);
        while (pos > 0)
        {
            pos--;
            printf("%02d", tmp_str[pos]);
        }
        delete[] tmp_str;
    }
    int quick_mod() const
    {
        size_t len = length();
        UINT_64 sum = 0;
        for (size_t i = 1; i < len; i++)
        {
            UINT_32 tmp = data.array[i];
            sum += tmp & INT32_0XFF;
            sum += ((UINT_64)tmp >> INT_BIT);
        }
        sum = (sum % 10) * 6 % 10;
        sum += data.array[0];
        sum = sum % 10;
        return sum;
    }
    HyperInt &self_half()
    {
        UINT_32 tmp1, tmp2 = 0;
        size_t pos = length();
        while (pos)
        {
            pos--;
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (tmp1 << (INT_BIT - 1)) + (tmp2 >> 1);
            data.array[pos] = tmp1;
        }
        set_true_len();
        return *this;
    } // /2
    HyperInt &self_twice()
    {
        UINT_32 tmp1, tmp2 = 0;
        for (size_t pos = 0; pos < length(); pos++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (tmp1 >> (INT_BIT - 1)) + (tmp2 << 1);
            data.array[pos] = tmp1;
        }
        if (data.size > length())
        {
            data.array[length()] = tmp2 >> (INT_BIT - 1);
            change_length(length() + 1);
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
        if (n == 0)
            return *this;
        HyperInt result;
        size_t shift = n >> LOG_INT_BIT;
        size_t offset = n % 32;
        size_t len = length() - shift;
        UINT_32 tmp1, tmp2 = 0;
        result.reset_size(len);
        result.change_length(len);
        while (len > 0)
        {
            len--;
            tmp1 = tmp2;
            tmp2 = data.array[len + shift];
            tmp1 = (offset ? (tmp1 << (INT_BIT - offset)) : 0) + (tmp2 >> offset);
            result.data.array[len] = tmp1;
        }
        result.set_true_len();
        return result;
    } //右移函数
    HyperInt l_shift(size_t n) const
    {
        if (n == 0)
            return *this;
        HyperInt result;
        size_t shift = n >> LOG_INT_BIT;
        size_t offset = n % 32;
        UINT_32 tmp1 = 0, tmp2 = 0;
        result.reset_size(length() + shift + 1);
        result.clear();
        result.change_length(length() + shift + 1);
        for (size_t pos = 0; pos < length(); pos++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (offset ? (tmp1 >> (INT_BIT - offset)) : 0) + (tmp2 << offset);
            result.data.array[pos + shift] = tmp1;
        }
        if (tmp2 >> (INT_BIT - offset) & offset)
            result.data.array[length() + shift] = tmp2 >> (INT_BIT - offset);
        else
            result.change_length(result.length() - 1);
        return result;
    } //左移函数
    void reset_size(size_t new_size)
    {
        if (_TWICE_(new_size) < new_size)
        {
            data.size = new_size;
        }
        else if (new_size < 2)
        {
            data.size = 2;
        }
        else
        {
            size_t size_tmp = 2, half = 0, mid = 0;
            while (size_tmp < new_size)
            {
                _SELFTWICE_(size_tmp);
            }
            half = _HALF_(size_tmp);
            mid = _HALF_(size_tmp + half);
            if (new_size > mid)
            {
                data.size = size_tmp;
            }
            else
            {
                data.size = mid;
            }
        }
        _INT32REARY_(data.array, data.size);
        change_length(_MIN_(data.size, length()));
    } //重新设定长度不小于new_size，1.5倍长度算法
    HyperInt abs() const
    {
        HyperInt result(*this);
        result.neg_sign(false);
        return result;
    }
    bool abs_larger(const HyperInt &input) const
    {
        size_t t_len1 = length(), t_len2 = input.length();
        if (t_len1 > t_len2)
            return true;
        else if (t_len1 < t_len2)
            return false;
        else
        {
            UINT_32 num1 = 0, num2 = 0;
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
        size_t t_len1 = length(), t_len2 = input.length();
        if (t_len1 < t_len2)
            return true;
        else if (t_len1 > t_len2)
            return false;
        else
        {
            UINT_32 num1 = 0, num2 = 0;
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
        size_t t_len1 = length(), t_len2 = input.length();
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
    } //比较绝对值是否相等
    bool equal_to_z() const
    {
        return (0 == length());
    } //判定是否为零
    bool is_even() const
    {
        if (length() == 0)
        {
            return true;
        }
        return !(bool)data.array[0] & 1;
    }
    bool is_odd() const
    {
        if (length() == 0)
        {
            return false;
        }
        return (bool)data.array[0] & 1;
    }
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
    HyperInt operator/(INT_64 input) const;

    HyperInt operator%(const HyperInt &input) const;
    INT_64 operator%(INT_64 input) const;

    HyperInt &operator+=(const HyperInt &input);
    HyperInt &operator-=(const HyperInt &input);

    HyperInt &operator*=(const HyperInt &input);
    HyperInt &operator*=(INT_64 input);

    HyperInt &operator/=(const HyperInt &input);
    HyperInt &operator%=(const HyperInt &input);
    HyperInt operator++(int);
    HyperInt &operator++();
    HyperInt operator--(int);
    HyperInt &operator--();
    // HyperInt &div_ten(INT_64 &rem_num);
    // HyperInt &div_h(INT_64 div_num, INT_64 &rem_num);
    // HyperInt power(UINT_64 n) const;
};
//逻辑运算
inline bool HyperInt::operator==(const HyperInt &input) const
{
    if (is_neg() ^ input.is_neg())
        return false;
    else
        return abs_equal(input);
}
inline bool HyperInt::operator==(const INT_64 &input) const
{
    if (is_neg() ^ _NEG_(input))
        return false;
    else
        return abs_equal(input);
}
bool HyperInt::operator>(const HyperInt &input) const
{
    if (is_neg() ^ input.is_neg())
        return input.is_neg();
    else
        return is_neg() ^ abs_larger(input);
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
    if (is_neg() ^ input.is_neg())
        return is_neg();
    else
        return is_neg() ^ abs_smaller(input);
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
    return !(bool)length();
}
// inline HyperInt::operator bool()
// {
//     if(length()>2)
//         return true;
//     return(data.array[0]!=0||data.array[1]!=0);
// }

//友元函数
HyperInt abs(HyperInt input)
{
    HyperInt result(input);
    result.neg_sign(false);
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
    HyperInt result;
    if (!is_neg() ^ input.is_neg())
    {
        result = add(input);
        result.neg_sign(is_neg());
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
            result.neg_sign(is_neg());
            return result;
        }
        else
        {
            result = input.sub(*this);
            result.neg_sign(input.is_neg());
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
    HyperInt result;
    if (is_neg() ^ input.is_neg())
    {
        result = add(input);
        result.neg_sign(is_neg());
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
            result.neg_sign(is_neg());
            return result;
        }
        else
        {
            result = input.sub(*this);
            result.neg_sign(input.is_neg());
            return result;
        }
    }
}
inline HyperInt HyperInt::operator-() const
{
    HyperInt result(*this);
    result.neg_sign(!is_neg());
    return result;
}
HyperInt HyperInt::operator*(const HyperInt &input) const
{
    size_t min_len = _MIN_(length(), input.length());
    if (min_len > FFT_MIN && min_len < FFT_MAX)
    //if(true)
    {
        return fft_multiply(input);
    }
    else
    {
        return multiply(input);
    }
}
HyperInt HyperInt::operator*(INT_64 input) const
{
    HyperInt result;
    if (equal_to_z() || input == 0)
    {
        return result;
    }
    size_t len = length();
    size_t result_len = len + 2;
    result.reset_size(result_len);
    result.change_length(result_len);
    bool result_neg = is_neg() ^ _NEG_(input);
    result.neg_sign(result_neg);
    input = _ABS_(input);
    UINT_32 tmp1, tmp2;
    UINT_64 sum, input_num1, input_num2;
    input_num1 = input & INT32_0XFF;
    input_num2 = input >> INT_BIT;
    tmp2 = data.array[0];
    sum = tmp2 * input_num1;
    result.data.array[0] = sum & INT32_0XFF;
    sum >>= INT_BIT;
    for (size_t pos = 1; pos < result_len; pos++)
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
    result.set_true_len();
    return result;
}
HyperInt HyperInt::operator/(const HyperInt &input) const
{
    if (!input)
        throw "Can't divided by zero";
    HyperInt result, mid, tmp;
    if (abs_smaller(input) || equal_to_z())
        return result;
    else
    {
        size_t len1 = length(), len2 = input.length(), i = 0;
        tmp.reset_size(len1 - len1 + 1);
        bool result_neg = is_neg() ^ input.is_neg();
        result.neg_sign(result_neg);
        tmp = 1;
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
        result.set_true_len();
        return result;
    }
}
HyperInt HyperInt::operator/(INT_64 input) const
{
    HyperInt result(*this);
    result.div_mod(input);
    return result;
}
HyperInt HyperInt::operator%(const HyperInt &input) const
{
    HyperInt tmp = *this / input;
    tmp *= input;
    return *this - tmp;
}
INT_64 HyperInt::operator%(INT_64 input) const
{
    INT_64 rem = mod(input);
    return rem;
}
HyperInt &HyperInt::operator+=(const HyperInt &input)
{
    *this = *this + input;
    return *this;
}
HyperInt &HyperInt::operator-=(const HyperInt &input)
{
    *this = *this - input;
    return *this;
}
HyperInt &HyperInt::operator*=(const HyperInt &input)
{
    *this = *this * input;
    return *this;
}
HyperInt &HyperInt::operator*=(INT_64 input)
{
    if (equal_to_z() || input == 0)
    {
        HyperInt();
        return *this;
    }
    size_t len = length();
    if (data.size < len + 2)
    {
        reset_size(len + 2);
    }
    change_length(len + 2);
    bool result_neg = is_neg() ^ _NEG_(input);
    neg_sign(result_neg);
    input = _ABS_(input);
    UINT_32 tmp1, tmp2;
    UINT_64 sum, input_num1, input_num2;
    input_num1 = input & INT32_0XFF;
    input_num2 = input >> INT_BIT;
    tmp2 = data.array[0];
    sum = tmp2 * input_num1;
    data.array[0] = sum & INT32_0XFF;
    sum >>= INT_BIT;
    for (size_t pos = 1; pos < len + 2; pos++)
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
        data.array[pos] = sum & INT32_0XFF;
        sum >>= INT_BIT;
    }
    set_true_len();
    return *this;
}
HyperInt &HyperInt::operator/=(const HyperInt &input)
{
    *this = *this / input;
    return *this;
}
HyperInt &HyperInt::operator%=(const HyperInt &input)
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
#endif