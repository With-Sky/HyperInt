#include <iostream>
#include <future>
#include <thread>
#include <atomic>
#include <algorithm>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cassert>

#ifndef HINT_HPP
#define HINT_HPP

#define UINT_8 uint8_t
#define UINT_16 uint16_t
#define UINT_32 uint32_t
#define UINT_64 uint64_t
#define INT_32 int32_t
#define INT_64 int64_t
#define ULONG unsigned long
#define LONG long

#if SIZE_MAX == 18446744073709551615ull
#define SIZE_T_BITS 64
#elif SIZE_MAX == 4294967295
#define SIZE_T_BITS 32
#else
#error "unknown"
#endif

#define MULTITHREAD
#define HINT_CHAR_BIT 8
#define HINT_SHORT_BIT 16
#define HINT_INT_BIT 32
#define HINT_INT8_0XFF UINT8_MAX
#define HINT_INT16_0XFF UINT16_MAX
#define HINT_INT16_0X10 0X10000ull
#define HINT_INT32_0XFF UINT32_MAX
#define HINT_INT32_0X01 (UINT_32)1
#define HINT_INT32_0X80 0X80000000u
#define HINT_INT32_0X7F INT32_MAX
#define HINT_INT32_0X10 0X100000000ull
#define HINT_INT64_0X80 INT64_MIN
#define HINT_INT64_0X7F INT64_MAX
#define HINT_FFT_MIN 128ull
#define HINT_FFT_MAX 4096ull        // 4096为2分fft结果最大长度
#define HINT_QUAL_FFT_MAX 262144ull // 2^18为4分fft最大长度
#define HINT_NTT_MIN 640ull
#define HINT_NTT_MAX 16777216ull  // 2^24为2分ntt结果最大长度
#define HINT_NTT_MULTHLEN 1024ull //设置触发多线程的ntt长度
#define HINT_PI 3.1415926535897932384626433832795
#define HINT_2PI 6.283185307179586476925286766559

#define _MAX_(x, y) std::max(x, y) //((x) > (y) ? (x) : (y))
#define _MIN_(x, y) std::min(x, y) //((x) < (y) ? (x) : (y))
#define _NEG_(x) ((x) < 0)
#define _ODD_(x) ((x)&1)
#define _ABS_(x) std::abs(x) //((x) < 0 ? (-(x)) : (x))
#define _TWICE_(x) ((x) << 1)
#define _HALF_(x) ((x) >> 1)
#define _SELFTWICE_(x) ((x) <<= 1)
#define _SELFHALF_(x) ((x) >>= 1)

namespace hint
{
    struct h_int
    {
        UINT_32 *array = nullptr;
        INT_64 neg_n_len = 0;
        size_t size = 0;
    };
    struct Complex //对复数的定义
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
        Complex(UINT_64 n) // n等分圆周的复数,即x^n=1的解中除x=1以外辐角最小的一个解
        {
            real = cos(HINT_2PI / n);
            imaginary = sin(HINT_2PI / n);
        }
        Complex operator=(Complex input)
        {
            real = input.real;
            imaginary = input.imaginary;
            return *this;
        }
        Complex operator+(Complex input) //复数加法
        {
            return Complex(real + input.real, imaginary + input.imaginary);
        }
        Complex operator-(Complex input) //复数减法
        {
            return Complex(real - input.real, imaginary - input.imaginary);
        }
        Complex operator*(Complex input) //复数乘法
        {
            return Complex(real * input.real - imaginary * input.imaginary, real * input.imaginary + imaginary * input.real);
        }
        Complex &operator*=(Complex input) //复数乘法
        {
            double tmp_real = real * input.real - imaginary * input.imaginary;
            imaginary = real * input.imaginary + imaginary * input.real;
            real = tmp_real;
            return *this;
        }
        Complex operator/(Complex input) //复数除法
        {
            double tmp = input.real * input.real + input.imaginary * input.imaginary;
            double re = real * input.real + imaginary * input.imaginary;
            double img = imaginary * input.real - real * input.imaginary;
            return Complex(re / tmp, img / tmp);
        }
        void console_out() //打印复数
        {
            std::cout << real;
            if (_NEG_(imaginary))
                std::cout << imaginary << "i";
            else
                std::cout << "+" << imaginary << "i";
        }
    };
    const UINT_32 hint_threads = std::thread::hardware_concurrency();
    const UINT_32 log2_threads = std::ceil(std::log2(hint_threads));
    const Complex unit_omega_ary[] = {Complex(1), Complex(2), Complex(4), Complex(8), Complex(16), Complex(32), Complex(64), Complex(128),
                                      Complex(256), Complex(512), Complex(1024), Complex(2048), Complex(4096), Complex(8192), Complex(16384),
                                      Complex(32768), Complex(131072), Complex(262144), Complex(524288), Complex(1048576), Complex(2097152)};
    UINT_8 *to_base100_unit_ary = nullptr;
    size_t to_base100_unit_ary_len = 0;
    template <typename T>
    T qpow(T m, UINT_64 n) //模板快速幂
    {
        T result = 1;
        while (n > 0)
        {
            if (n & 1 != 0)
            {
                result = result * m;
            }
            m = m * m;
            n >>= 1;
        }
        return result;
    }
    UINT_64 qpow(UINT_64 m, UINT_64 n, UINT_64 mod) //取模快速幂
    {
        UINT_64 result = 1;
        while (n > 0)
        {
            if ((n & 1) != 0)
            {
                result = result * m % mod;
            }
            m = m * m % mod;
            n >>= 1;
        }
        return result;
    }
    UINT_64 gcd(UINT_64 a, UINT_64 b) //最大公因数
    {
        if (b == 0)
        {
            return a;
        }
        UINT_64 tmp = b;
        b = a % b;
        a = tmp;
        while (b > 0)
        {
            tmp = b;
            b = a % b;
            a = tmp;
        }
        return a;
    }
    UINT_64 crt(UINT_64 *mods, UINT_64 *nums, size_t n) //中国剩余定理
    {
        UINT_64 result = 0, mod_product = 1;
        for (size_t i = 0; i < n; i++)
        {
            mod_product *= mods[i];
        }
        for (size_t i = 0; i < n; i++)
        {
            UINT_64 mod = mods[i];
            UINT_64 tmp = mod_product / mod;
            UINT_64 inv = qpow(tmp, mod - 2, mod);
            result += nums[i] * tmp * inv % mod_product;
        }
        return result % mod_product;
    }
    inline UINT_64 qcrt(UINT_64 num1, UINT_64 num2,
                        UINT_64 mod1 = 167772161, UINT_64 mod2 = 469762049,
                        UINT_64 inv1 = 104391568, UINT_64 inv2 = 130489458) //快速计算两模数的中国剩余定理
    {
        if (num1 > num2)
        {
            return ((num1 - num2) * inv2 % mod1) * mod2 + num2;
        }
        else
        {
            return ((num2 - num1) * inv1 % mod2) * mod1 + num1;
        }
    }
    template <typename T>
    inline T *ary_copy(T *target, const T *source, size_t len) //模板数组拷贝
    {
        return static_cast<T *>(std::memcpy(target, source, len * sizeof(T)));
    }
    template <typename T>
    inline void com_ary_copy(Complex *target, const T *source, size_t len) //从其他类型数组拷贝到复数组
    {
        for (size_t i = 0; i < len; i++)
        {
            target[i].real = source[i];
        }
    }
    template <typename T>
    inline void ary_calloc(T *&ptr, size_t len) //模版数组分配内存且清零
    {
        ptr = static_cast<T *>(calloc(len, sizeof(T)));
    }
    template <typename T>
    inline void ary_clr(T *ptr, size_t len) //模版数组清零
    {
        memset(ptr, 0, len * sizeof(T));
    }
    template <typename T>
    inline T *ary_realloc(T *ptr, size_t len) //重分配空间
    {
        return static_cast<T *>(realloc(ptr, len * sizeof(T)));
    }
    void fft(Complex *input, size_t fft_len, bool is_ifft) //快速傅里叶(逆)变换
    {
        size_t log_n = static_cast<size_t>(log2(fft_len));
        fft_len = 1ull << log_n;
        size_t *rev = new size_t[fft_len];
        rev[0] = 0;
        for (size_t i = 1; i < fft_len; i++)
        {
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log_n - 1)); //求rev交换数组
        }
        for (size_t i = 0; i < fft_len; i++)
        {
            if (i < rev[i])
            {
                Complex tmp = input[i];
                input[i] = input[rev[i]];
                input[rev[i]] = tmp;
            }
        }
        delete[] rev;
        Complex unit_omega, omega, tmp1, tmp2;
        for (size_t rank = 1, gap; rank < fft_len; _SELFTWICE_(rank))
        {
            gap = _TWICE_(rank);
            // unit_omega = unit_omega_ary[i];
            unit_omega = Complex(gap);
            if (is_ifft)
            {
                unit_omega.imaginary = -unit_omega.imaginary;
            }
            for (size_t begin = 0; begin < fft_len; begin += gap)
            {
                omega = Complex(1, 0);
                for (size_t pos = begin; pos < begin + rank; pos++)
                {
                    tmp1 = input[pos];
                    tmp2 = input[pos + rank] * omega;
                    input[pos] = tmp1 + tmp2;
                    input[pos + rank] = tmp1 - tmp2;
                    omega *= unit_omega;
                }
            }
        }
        if (is_ifft) //逆变换需除以n
        {
            double inv = (1.0 / fft_len);
            for (size_t i = 0; i < fft_len; i++)
            {
                input[i].real *= inv;
            }
        }
    }
    void ntt(UINT_64 *input, size_t ntt_len, bool is_intt, UINT_64 mod = 998244353, UINT_64 g_root = 3) //快速数论变换
    {
        size_t log_n = static_cast<size_t>(log2(ntt_len));
        size_t *rev = new size_t[ntt_len];
        ntt_len = 1ull << log_n;
        rev[0] = 0;
        for (size_t i = 1; i < ntt_len; i++)
        {
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log_n - 1));
        }
        for (size_t i = 0; i < ntt_len; i++)
        {
            if (i < rev[i])
            {
                UINT_64 tmp = input[i];
                input[i] = input[rev[i]];
                input[rev[i]] = tmp;
            }
        }
        delete[] rev;
        if (is_intt)
        {
            g_root = qpow(g_root, mod - 2, mod);
        }
        UINT_64 unit_omega, omega;
        for (size_t rank = 1, gap; rank < ntt_len; _SELFTWICE_(rank))
        {
            gap = _TWICE_(rank);
            unit_omega = qpow(g_root, (mod - 1) / gap, mod);
            for (size_t begin = 0; begin < ntt_len; begin += gap)
            {
                omega = 1;
                for (size_t pos = begin; pos < begin + rank; pos++)
                {
                    UINT_64 tmp1 = input[pos];
                    UINT_64 tmp2 = (input[pos + rank] % mod) * omega % mod;
                    input[pos] = (tmp1 + tmp2) % mod;
                    input[pos + rank] = (mod + tmp1 - tmp2) % mod;
                    omega = omega * unit_omega % mod;
                }
            }
        }
        if (is_intt)
        {
            UINT_64 inv = qpow(ntt_len, mod - 2, mod);
            for (size_t i = 0; i < ntt_len; ++i)
            {
                input[i] = input[i] * inv % mod;
            }
        }
    }
    void fft_convolution(Complex *const fft_ary1, Complex *const fft_ary2, Complex *const out, size_t fft_len) //快速傅里叶变换卷积分
    {
#ifdef MULTITHREAD
        // UINT_32 cores = std::thread::hardware_concurrency();
        bool multi_threads = hint::hint_threads >= 2 && fft_len > 2 * HINT_FFT_MAX;
        if (multi_threads)
        {
            std::future<void> fft_th = std::async(fft, fft_ary1, fft_len, false); //快速傅里叶变换
            if (fft_ary1 != fft_ary2)
            {
                fft(fft_ary2, fft_len, false);
            }
            fft_th.wait();
        }
        else
#endif
        {
            fft(fft_ary1, fft_len, false); //快速傅里叶变换
            if (fft_ary1 != fft_ary2)
            {
                fft(fft_ary2, fft_len, false);
            }
            //每一位相乘
        }
        for (size_t i = 0; i < fft_len; i++)
        {
            out[i] = fft_ary1[i] * fft_ary2[i];
        }
        fft(out, fft_len, true); //逆变换
    }
    void ntt_convolution(UINT_64 *const ntt_ary1, UINT_64 *const ntt_ary2, UINT_64 *const out, size_t ntt_len) //数论变换卷积分
    {
        const UINT_64 mod1 = 167772161ull, mod2 = 469762049ull;
        const UINT_64 root1 = 3, root2 = 3;
        UINT_64 *ntt_ary3 = new UINT_64[ntt_len];
        ary_copy(ntt_ary3, ntt_ary1, ntt_len);
        UINT_64 *ntt_ary4 = ntt_ary3;
#ifdef MULTITHREAD
        // UINT_32 cores = std::thread::hardware_concurrency();
        bool multi_threads = hint::hint_threads >= 4 && ntt_len > HINT_NTT_MULTHLEN;
        if (multi_threads)
        {
            std::future<void> ntt_th1 = std::async(ntt, ntt_ary1, ntt_len, false, mod1, root1); // 多线程快速数论变换
            std::future<void> ntt_th2 = std::async(ntt, ntt_ary3, ntt_len, false, mod2, root2);
            if (ntt_ary1 != ntt_ary2)
            {
                ntt_ary4 = new UINT_64[ntt_len];
                ary_copy(ntt_ary4, ntt_ary2, ntt_len);
                std::future<void> ntt_th3 = std::async(ntt, ntt_ary2, ntt_len, false, mod1, root1);
                std::future<void> ntt_th4 = std::async(ntt, ntt_ary4, ntt_len, false, mod2, root2);
                ntt_th3.wait();
                ntt_th4.wait();
            }
            ntt_th1.wait();
            ntt_th2.wait();
            std::function<void(UINT_64 *, UINT_64 *)> mul_func = [=](UINT_64 *ary1, UINT_64 *ary2)
            {
                for (size_t i = 0; i < ntt_len; i++)
                {
                    ary1[i] = ary1[i] * ary2[i];
                } //每一位相乘
            };
            std::future<void> mul_th = std::async(mul_func, ntt_ary1, ntt_ary2);
            mul_func(ntt_ary3, ntt_ary4);
            mul_th.wait();

            std::future<void> intt_th = std::async(ntt, ntt_ary1, ntt_len, true, mod1, root1);
            ntt(ntt_ary3, ntt_len, true, mod2, root2);
            intt_th.wait();
        }
        else
#endif
        {
            ntt(ntt_ary1, ntt_len, false, mod1, root1); //快速数论变换
            ntt(ntt_ary3, ntt_len, false, mod2, root2);
            if (ntt_ary1 != ntt_ary2)
            {
                ntt_ary4 = new UINT_64[ntt_len];
                ary_copy(ntt_ary4, ntt_ary2, ntt_len);
                ntt(ntt_ary2, ntt_len, false, mod1, root1); //快速数论变换
                ntt(ntt_ary4, ntt_len, false, mod2, root2);
            }
            for (size_t i = 0; i < ntt_len; i++)
            {
                ntt_ary1[i] = ntt_ary1[i] * ntt_ary2[i];
                ntt_ary3[i] = ntt_ary3[i] * ntt_ary4[i];
            }                                          //每一位相乘
            ntt(ntt_ary1, ntt_len, true, mod1, root1); //逆变换
            ntt(ntt_ary3, ntt_len, true, mod2, root2);
        }

        for (size_t i = 0; i < ntt_len; i++)
        {
            out[i] = qcrt(ntt_ary1[i], ntt_ary3[i]);
        } //使用中国剩余定理变换
        delete[] ntt_ary3;
        if (ntt_ary4 != ntt_ary3)
        {
            delete[] ntt_ary4;
        }
    }
    template <typename T>
    void trans_add(const T *in1, const T *in2, T *out, size_t len1, size_t len2, const INT_64 base = 100) //可计算多项式的加法,默认为100进制
    {
        size_t result_len = _MAX_(len1, len2) + 1;
        INT_64 tmp = 0;
        size_t count = 0;
        while (count < result_len)
        {
            if (count < len1)
            {
                tmp += in1[count];
            }
            if (count < len2)
            {
                tmp += in2[count];
            }
            out[count] = static_cast<T>(tmp % base);
            tmp /= base;
            count++;
        }
    }
    template <typename T>
    void trans_mul(const T *in1, const T *in2, T *out, size_t len1, size_t len2, const INT_64 base = 100) //计算多项式的乘法
    {
        size_t out_len = len1 + len2;
        size_t fft_len = 1ull << static_cast<UINT_16>(ceil(log2(out_len)));
        hint::Complex *fft_in1 = new hint::Complex[fft_len];
        hint::Complex *fft_in2 = new hint::Complex[fft_len];
        size_t pos = 0, max_len = _MAX_(len1, len2);
        while (pos < max_len) //放入复数组
        {
            if (pos < len1)
            {
                fft_in1[pos].real = in1[pos];
            }
            if (pos < len2)
            {
                fft_in2[pos].real = in2[pos];
            }
            pos++;
        }
        fft_convolution(fft_in1, fft_in2, fft_in1, fft_len);
        UINT_64 tmp = 0;
        pos = 0;
        while (pos < out_len)
        {
            tmp += static_cast<UINT_64>(fft_in1[pos].real + 0.5);
            out[pos] = static_cast<T>(tmp % base);
            tmp /= base;
            pos++;
        } //整理每一位
        delete[] fft_in1;
        delete[] fft_in2;
    }
    template <typename T>
    void trans_square(const T *in, T *out, size_t len, const INT_64 base = 100) //平方
    {
        size_t out_len = _TWICE_(len);
        size_t fft_len = 1ull << static_cast<UINT_16>(ceil(log2(out_len)));
        hint::Complex *fft_ary = new hint::Complex[fft_len];
        size_t pos = 0;
        while (pos < len) //放入复数组
        {
            fft_ary[pos].real = in[pos];
            pos++;
        }
        fft_convolution(fft_ary, fft_ary, fft_ary, fft_len);
        UINT_64 tmp = 0;
        pos = 0;
        while (pos < out_len)
        {
            tmp += static_cast<UINT_64>(fft_ary[pos].real + 0.5);
            out[pos] = static_cast<T>(tmp % base);
            tmp /= base;
            pos++;
        } //整理每一位
        delete[] fft_ary;
    }
    template <typename T>
    inline T hint_mul(const T &input1, const T &input2)
    {
        return input1 * input2;
    }
    template <typename T>
    inline T hint_square(const T &input)
    {
        return input.square();
    }
#ifdef MULTITHREAD
    void to_base100_sub_prod(UINT_8 *unit_ary, UINT_8 *result, UINT_8 *tmp_product, size_t len, size_t pos, size_t offset, const size_t result_len)
    {
        static std::atomic<UINT_32> ths;
        // const UINT_32 cores = std::thread::hardware_concurrency();
        if (pos < result_len)
        {
            size_t gap = _TWICE_(len);
            if (ths < hint::hint_threads)
            {
                ths++;
                std::future<void> cur_thread = std::async(to_base100_sub_prod, unit_ary, result, tmp_product, len, pos + gap, offset, result_len);
                hint::trans_mul(unit_ary + offset, result + pos + len, tmp_product, len, len, 100);
                hint::trans_add(tmp_product, result + pos, tmp_product, gap, len, 100);
                hint::ary_copy(result + pos, tmp_product, gap);
                cur_thread.wait();
                ths--;
            }
            else
            {
                hint::trans_mul(unit_ary + offset, result + pos + len, tmp_product, len, len, 100);
                hint::trans_add(tmp_product, result + pos, tmp_product, gap, len, 100);
                hint::ary_copy(result + pos, tmp_product, gap);
                to_base100_sub_prod(unit_ary, result, tmp_product, len, pos + gap, offset, result_len);
            }
        }
    }
#endif
}
class HyperInt
{
private:
    hint::h_int data;
    void change_length(size_t new_length) //设置新的长度
    {
        if (new_length > data.size)
        {
            new_length = data.size;
        }
        data.neg_n_len = (data.neg_n_len & HINT_INT64_0X80) | new_length;
    }
    size_t generate_size(size_t new_size_input) const //生成1.5倍数组空间
    {
        // const size_t max = HINT_NTT_MAX; // 2^24 长度 //占用内存64MB
        // if (new_size_input > max)
        // {
        //     return max;
        // }
        if (new_size_input <= 2)
        {
            return 2;
        }
        size_t size1 = 1ull << static_cast<UINT_16>(ceil(log2(new_size_input)));
        size_t size2 = _HALF_(size1);
        size2 = size2 + _HALF_(size2);
        if (new_size_input <= size2)
        {
            return size2;
        }
        else
        {
            return size1;
        }
    }
    void clear() //清空
    {
        hint::ary_clr(data.array, data.size);
    }
    void fill_element(UINT_32 element) //填充相同的元素
    {
        std::fill(data.array, data.array + length(), element);
    }
    void quick_self_twice() //快速左移一位,只能操作大小为2^n的数
    {
        size_t len_tmp = length();
        if (len_tmp > 0)
        {
            UINT_32 tmp = data.array[len_tmp - 1];
            _SELFTWICE_(tmp);
            if (tmp == 0)
            {
                len_tmp++;
                if (len_tmp > data.size)
                {
                    return;
                }
                change_length(len_tmp);
                data.array[len_tmp - 2] = 0;
                data.array[len_tmp - 1] = 1;
            }
            else
            {
                data.array[len_tmp - 1] = tmp;
            }
        }
    }                      //快速自增为二倍
    void quick_self_half() //快速右移一位,只能操作大小为2^n的数
    {
        size_t len_tmp = length();
        if (len_tmp > 0)
        {
            UINT_32 tmp = data.array[len_tmp - 1];
            _SELFHALF_(tmp);
            if (tmp == 0)
            {
                len_tmp--;
                if (len_tmp == 0)
                {
                    return;
                }
                change_length(len_tmp);
                data.array[len_tmp - 1] = HINT_INT32_0X80;
            }
            else
            {
                data.array[len_tmp - 1] = tmp;
            }
        }
    }
    HyperInt quick_twice() //将自身的二倍返回
    {
        HyperInt result(*this);
        result.quick_self_twice();
        result.set_true_len();
        return result;
    }
    HyperInt quick_half() //将自身的一半返回
    {
        HyperInt result(*this);
        result.quick_self_half();
        result.set_true_len();
        return result;
    }
    HyperInt normal_multiply(const HyperInt &input) const //基础乘法
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
                tmp = static_cast<UINT_64>(data.array[pos1]) * input.data.array[pos2];
                for (size_t pos3 = pos1 + pos2; pos3 < result_len; pos3++)
                {
                    sum = tmp + result.data.array[pos3];
                    result.data.array[pos3] = static_cast<UINT_32>(sum & HINT_INT32_0XFF);
                    if ((sum >> HINT_INT_BIT) > 0)
                    {
                        tmp = sum >> HINT_INT_BIT;
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
    }
    HyperInt normal_square() const //基础平方
    {
        HyperInt result;
        if (equal_to_z())
        {
            return result;
        }
        size_t len = length();
        size_t result_len = _TWICE_(len);
        result.reset_size(result_len);
        result.clear();
        result.change_length(result_len);
        UINT_64 tmp = 0, sum = 0;
        for (size_t pos1 = 0; pos1 < len; pos1++)
        {
            for (size_t pos2 = 0; pos2 < len; pos2++)
            {
                tmp = static_cast<UINT_64>(data.array[pos1]) * data.array[pos2];
                for (size_t pos3 = pos1 + pos2; pos3 < result_len; pos3++)
                {
                    sum = tmp + result.data.array[pos3];
                    result.data.array[pos3] = static_cast<UINT_32>(sum & HINT_INT32_0XFF);
                    if ((sum >> HINT_INT_BIT) > 0)
                    {
                        tmp = sum >> HINT_INT_BIT;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        result.neg_sign(false);
        result.set_true_len();
        return result;
    }
    HyperInt fft_multiply(const HyperInt &input) const //快速傅里叶变换乘法
    {
        HyperInt result;
        if (equal_to_z() || input.equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        size_t out_len = len1 + len2;
        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);
        size_t fft_len = 2ull << static_cast<UINT_16>(ceil(log2(out_len)));
        bool is_qual_div = out_len > HINT_FFT_MAX;
        if (is_qual_div)
        {
            _SELFTWICE_(fft_len);
        }
        hint::Complex *fft_ary1 = new hint::Complex[fft_len * 2];
        hint::Complex *fft_ary2 = fft_ary1 + fft_len; // new hint::Complex[fft_len];

        //每一位分解为短整数后存入复数组
        if (is_qual_div)
        {
            UINT_8 *ary1_8 = reinterpret_cast<UINT_8 *>(data.array);
            UINT_8 *ary2_8 = reinterpret_cast<UINT_8 *>(input.data.array);
            size_t data_len1 = len1 * 4;
            size_t data_len2 = len2 * 4;

            hint::com_ary_copy(fft_ary1, ary1_8, data_len1);
            hint::com_ary_copy(fft_ary2, ary2_8, data_len2);
        }
        else
        {
            UINT_16 *ary1_16 = reinterpret_cast<UINT_16 *>(data.array);
            UINT_16 *ary2_16 = reinterpret_cast<UINT_16 *>(input.data.array);
            size_t data_len1 = len1 * 2;
            size_t data_len2 = len2 * 2;

            hint::com_ary_copy(fft_ary1, ary1_16, data_len1);
            hint::com_ary_copy(fft_ary2, ary2_16, data_len2);
        }

        hint::fft_convolution(fft_ary1, fft_ary2, fft_ary1, fft_len);

        UINT_64 tmp = 0;
        if (is_qual_div)
        {
            UINT_8 *ary_8 = reinterpret_cast<UINT_8 *>(result.data.array);
            size_t data_len = out_len * 4;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<UINT_64>(fft_ary1[i].real + 0.5);
                ary_8[i] = static_cast<UINT_8>(tmp & HINT_INT8_0XFF);
                tmp >>= HINT_CHAR_BIT;
            }
        }
        else
        {
            UINT_16 *ary_16 = reinterpret_cast<UINT_16 *>(result.data.array);
            size_t data_len = out_len * 2;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<UINT_64>(fft_ary1[i].real + 0.5);
                ary_16[i] = static_cast<UINT_16>(tmp & HINT_INT16_0XFF);
                tmp >>= HINT_SHORT_BIT;
            }
        } //整理每一位

        delete[] fft_ary1;
        // delete[] fft_ary2;

        result.set_true_len();
        bool result_neg = is_neg() ^ input.is_neg();
        result.neg_sign(result_neg);
        return result;
    }
    HyperInt fft_square() const // fft平方计算
    {
        HyperInt result;
        if (equal_to_z())
        {
            return result;
        }
        size_t len = length();
        size_t out_len = _TWICE_(len);

        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);
        size_t fft_len = 2ull << static_cast<UINT_16>(ceil(log2(out_len)));

        bool is_qual_div = out_len > HINT_FFT_MAX;
        if (is_qual_div)
        {
            _SELFTWICE_(fft_len);
        }
        hint::Complex *fft_ary = new hint::Complex[fft_len];
        //每一位分解为短整数后存入复数组
        if (is_qual_div)
        {
            UINT_8 *ary_8 = reinterpret_cast<UINT_8 *>(data.array);
            size_t data_len = len * 4;

            hint::com_ary_copy(fft_ary, ary_8, data_len);
        }
        else
        {
            UINT_16 *ary_16 = reinterpret_cast<UINT_16 *>(data.array);
            size_t data_len = len * 2;

            hint::com_ary_copy(fft_ary, ary_16, data_len);
        }

        hint::fft_convolution(fft_ary, fft_ary, fft_ary, fft_len); //卷积

        UINT_64 tmp = 0;
        if (is_qual_div)
        {
            UINT_8 *ary_8 = reinterpret_cast<UINT_8 *>(result.data.array);
            size_t data_len = out_len * 4;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<UINT_64>(fft_ary[i].real + 0.5);
                ary_8[i] = static_cast<UINT_8>(tmp & HINT_INT8_0XFF);
                tmp >>= HINT_CHAR_BIT;
            }
        }
        else
        {
            UINT_16 *ary_16 = reinterpret_cast<UINT_16 *>(result.data.array);
            size_t data_len = out_len * 2;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<UINT_64>(fft_ary[i].real + 0.5);
                ary_16[i] = static_cast<UINT_16>(tmp & HINT_INT16_0XFF);
                tmp >>= HINT_SHORT_BIT;
            }
        } //整理每一位

        delete[] fft_ary;
        result.set_true_len();
        return result;
    }
    HyperInt ntt_multiply(const HyperInt &input) const //快速数论变换乘法
    {
        HyperInt result;
        if (equal_to_z() || input.equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        size_t out_len = len1 + len2;
        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);

        size_t ntt_len = 2ull << static_cast<UINT_16>(ceil(log2(out_len)));
        UINT_64 *ntt_ary1 = new UINT_64[ntt_len * 2];
        UINT_64 *ntt_ary2 = ntt_ary1 + ntt_len; // new UINT_64[ntt_len];

        hint::ary_clr(ntt_ary1, ntt_len * 2);
        // hint::ary_clr(ntt_ary2, ntt_len);

        UINT_16 *ary1_16 = reinterpret_cast<UINT_16 *>(data.array); //每一位分解为短整数后存入数组
        UINT_16 *ary2_16 = reinterpret_cast<UINT_16 *>(input.data.array);
        size_t data_len1 = len1 * 2;
        size_t data_len2 = len2 * 2;
        for (size_t i = 0; i < data_len1; i++)
        {
            ntt_ary1[i] = ary1_16[i];
        }
        for (size_t i = 0; i < data_len2; i++)
        {
            ntt_ary2[i] = ary2_16[i];
        }

        hint::ntt_convolution(ntt_ary1, ntt_ary2, ntt_ary1, ntt_len);

        UINT_64 tmp = 0;
        UINT_16 *ary_16 = reinterpret_cast<UINT_16 *>(result.data.array);
        size_t data_len = out_len * 2;
        for (size_t i = 0; i < data_len; i++)
        {
            tmp += ntt_ary1[i];
            ary_16[i] = static_cast<UINT_16>(tmp & HINT_INT16_0XFF);
            tmp >>= HINT_SHORT_BIT;
        }
        //整理每一位
        delete[] ntt_ary1;
        // delete[] ntt_ary2;
        result.set_true_len();
        bool result_neg = is_neg() ^ input.is_neg();
        result.neg_sign(result_neg);
        return result;
    }
    HyperInt ntt_square() const //快速数论变换平方
    {
        HyperInt result;
        if (equal_to_z())
        {
            return result;
        }
        size_t len = length();
        size_t out_len = _TWICE_(len);
        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);

        size_t ntt_len = 2ull << static_cast<UINT_16>(ceil(log2(out_len)));
        UINT_64 *ntt_ary = new UINT_64[ntt_len];
        hint::ary_clr(ntt_ary, ntt_len);

        //每一位分解为短整数后存入数组
        UINT_16 *ary_16 = reinterpret_cast<UINT_16 *>(data.array);
        size_t data_len = len * 2;
        for (size_t i = 0; i < data_len; i++)
        {
            ntt_ary[i] = ary_16[i];
        }

        hint::ntt_convolution(ntt_ary, ntt_ary, ntt_ary, ntt_len);

        UINT_64 tmp = 0;
        ary_16 = reinterpret_cast<UINT_16 *>(result.data.array);
        data_len = out_len * 2;
        for (size_t i = 0; i < data_len; i++)
        {
            tmp += ntt_ary[i];
            ary_16[i] = static_cast<UINT_16>(tmp & HINT_INT16_0XFF);
            tmp >>= HINT_SHORT_BIT;
        } //整理每一位
        delete[] ntt_ary;
        result.set_true_len();
        return result;
    }
    HyperInt karatsuba_multiply(const HyperInt &input) const // karatsuba乘法,速度较慢
    {
        HyperInt result;
        if (equal_to_z() || input.equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        size_t result_len = len1 + len2;
        if (result_len < 4)
        {
            return normal_multiply(input);
        }
        result.reset_size(result_len);
        result.clear();
        result.change_length(result_len);
        size_t sub_len1, sub_len2, sub_len_public;
        if (len1 > len2)
        {
            sub_len1 = _HALF_(len1);
            sub_len_public = len1 - sub_len1;
            sub_len2 = 0;
            if (len2 > sub_len_public)
            {
                sub_len2 = len2 - sub_len_public;
            }
        }
        else
        {
            sub_len2 = _HALF_(len2);
            sub_len_public = len2 - sub_len2;
            sub_len1 = 0;
            if (len1 > sub_len_public)
            {
                sub_len1 = len1 - sub_len_public;
            }
        }
        HyperInt sub_a = split(0, sub_len_public);
        HyperInt sub_b = split(sub_len_public, sub_len1);
        HyperInt sub_c = input.split(0, sub_len_public);
        HyperInt sub_d = input.split(sub_len_public, sub_len2);

        HyperInt sub_e, sub_f, sub_g;
#ifdef MULTITHREAD
        if (hint::hint_threads > 1)
        {
            std::future<HyperInt> sub_e_th = std::async(hint::hint_mul<HyperInt>, sub_a, sub_c);
            std::future<HyperInt> sub_f_th = std::async(hint::hint_mul<HyperInt>, sub_b, sub_d);
            sub_g = (sub_a + sub_b) * (sub_c + sub_d);
            sub_e = sub_e_th.get();
            sub_f = sub_f_th.get();
        }
        else
#endif
        {
            sub_e = sub_a * sub_c;
            sub_f = sub_b * sub_d;
            sub_g = (sub_a + sub_b) * (sub_c + sub_d);
        }

        size_t count = 0, pos_1 = 0, pos_2 = 0, pos_3 = 0, pos_4 = 0, len_pub = sub_len_public;
        size_t len_e = sub_e.length(), len_f = sub_f.length(), len_g = sub_g.length();
        INT_64 tmp = 0;
        while (count < len_pub)
        {
            if (pos_1 < len_e)
            {
                tmp += sub_e.data.array[pos_1];
                pos_1++;
            }
            result.data.array[count] = static_cast<UINT_32>(tmp + HINT_INT32_0X10);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        _SELFTWICE_(len_pub);
        while (count < len_pub)
        {
            if (pos_1 < len_e)
            {
                tmp += sub_e.data.array[pos_1];
                pos_1++;
            }
            if (pos_2 < len_g)
            {
                tmp += sub_g.data.array[pos_2];
                pos_2++;
            }
            if (pos_3 < len_e)
            {
                tmp -= sub_e.data.array[pos_3];
                pos_3++;
            }
            if (pos_4 < len_f)
            {
                tmp -= sub_f.data.array[pos_4];
                pos_4++;
            }
            result.data.array[count] = static_cast<UINT_32>(tmp + HINT_INT32_0X10);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        pos_1 = 0;
        while (count < result_len)
        {
            if (pos_1 < len_f)
            {
                tmp += sub_f.data.array[pos_1];
                pos_1++;
            }
            if (pos_2 < len_g)
            {
                tmp += sub_g.data.array[pos_2];
                pos_2++;
            }
            if (pos_3 < len_e)
            {
                tmp -= sub_e.data.array[pos_3];
                pos_3++;
            }
            if (pos_4 < len_f)
            {
                tmp -= sub_f.data.array[pos_4];
                pos_4++;
            }
            result.data.array[count] = static_cast<UINT_32>(tmp + HINT_INT32_0X10);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        result.set_true_len();
        return result;
    }
    HyperInt karatsuba_square() const // karatsuba平方,速度较慢
    {
        HyperInt result;
        if (equal_to_z())
        {
            return result;
        }
        size_t len = length();
        size_t result_len = _TWICE_(len);
        if (result_len < 4)
        {
            return normal_square();
        }
        result.reset_size(result_len);
        result.clear();
        result.change_length(result_len);
        size_t sub_len, sub_len_public;
        sub_len = _HALF_(len);
        sub_len_public = len - sub_len;

        HyperInt sub_a = split(0, sub_len_public);
        HyperInt sub_b = split(sub_len_public, sub_len);

        HyperInt sub_ab;
#ifdef MULTITHREAD
        if (hint::hint_threads > 1)
        {
            std::future<HyperInt> sub_a_th = std::async(hint::hint_square<HyperInt>, sub_a);
            std::future<HyperInt> sub_b_th = std::async(hint::hint_square<HyperInt>, sub_b);
            sub_ab = sub_a * sub_b;
            sub_a = sub_a_th.get();
            sub_b = sub_b_th.get();
        }
        else
#endif
        {
            sub_ab = sub_a * sub_b;
            sub_a = sub_a.square();
            sub_b = sub_b.square();
        }

        sub_ab.self_twice();

        size_t count = 0, pos_1 = 0, pos_2 = 0, len_pub = sub_len_public;
        size_t len_aa = sub_a.length(), len_bb = sub_b.length(), len_ab = sub_ab.length();
        INT_64 tmp = 0;
        while (count < len_pub)
        {
            if (pos_1 < len_aa)
            {
                tmp += sub_a.data.array[pos_1];
                pos_1++;
            }
            result.data.array[count] = static_cast<UINT_32>(tmp & HINT_INT32_0XFF);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        _SELFTWICE_(len_pub);
        while (count < len_pub)
        {
            if (pos_1 < len_aa)
            {
                tmp += sub_a.data.array[pos_1];
                pos_1++;
            }
            if (pos_2 < len_ab)
            {
                tmp += sub_ab.data.array[pos_2];
                pos_2++;
            }
            result.data.array[count] = static_cast<UINT_32>(tmp & HINT_INT32_0XFF);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        pos_1 = 0;
        while (count < result_len)
        {
            if (pos_2 < len_ab)
            {
                tmp += sub_ab.data.array[pos_2];
                pos_2++;
            }
            if (pos_1 < len_bb)
            {
                tmp += sub_b.data.array[pos_1];
                pos_1++;
            }
            result.data.array[count] = static_cast<UINT_32>(tmp & HINT_INT32_0XFF);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        result.set_true_len();
        return result;
    }
    HyperInt newton_divide(const HyperInt &input) const //牛顿迭代法除法
    {
        assert(!input.equal_to_z());
        HyperInt result;
        if (abs_smaller(input) || equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        const size_t precise_len = _TWICE_(len1);

        HyperInt input_inv = 1;      //除数的倒数
        size_t input_inv_len = len2; //小数部分长度
        const size_t times = 320;
        for (size_t i = 0; i < times; i++)
        {
            HyperInt input_inv_tmp = input_inv;
            HyperInt input_inv_sq = input_inv.square(); //除数倒数的平方,小数部分为两倍input_inv_len
            input_inv_sq = input * input_inv_sq;
            size_t input_inv_sq_len = input_inv_len * 2;
            input_inv = input_inv.l_shift((input_inv_sq_len - input_inv_len) * HINT_INT_BIT + 1);
            input_inv.add_sub_inplace(input_inv_sq, false);
            input_inv_len = input_inv_sq_len;
            if (input_inv.length() > precise_len)
            {
                size_t cut_len = input_inv.length() - precise_len;
                input_inv = input_inv.split(cut_len, precise_len);
                input_inv_len -= cut_len;
            }
            if (!input_inv.abs_larger(input_inv_tmp))
            {
                break;
            }
            // std::cout << i << "\n";
        }
        input_inv++;
        result = (*this * input_inv);
        result = result.split(input_inv_len, result.length() - input_inv_len);
        if (abs_smaller(result * input))
        {
            result--;
        }
        result.neg_sign(is_neg() ^ input.is_neg());
        return result;
    }
    HyperInt normal_divide(const HyperInt &input) const //模拟手算的除法
    {
        assert(!input.equal_to_z());
        HyperInt result;
        if (abs_smaller(input) || equal_to_z())
        {
            return result;
        }
        size_t len1 = length(), len2 = input.length();
        if (len1 <= 2)
        {
            result = first_int64() / input.first_int64();
            result.neg_sign(static_cast<bool>(is_neg() ^ input.is_neg()));
            return result;
        }
        else if (len2 < 2)
        {
            result = *this;
            result.div_mod(static_cast<INT_64>(input.first_int32()));
            result.neg_sign(static_cast<bool>(is_neg() ^ input.is_neg()));
            return result;
        }
        size_t out_len = len1 - len2 + 1;
        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);
        HyperInt dividend(*this);
        HyperInt sub;
        size_t shift = 0, pos = 0;
        UINT_64 tmp = 0, first_num2 = input.first_int64();
        UINT_64 try_num = 0;
        while (!dividend.abs_smaller(input))
        {
            shift = dividend.length() - len2;
            UINT_64 first_num1 = dividend.first_int64();
            if (first_num1 > first_num2)
            {
                tmp = first_num1 / first_num2;
                sub = input * tmp;
                if (dividend.abs_compare(sub, shift) < 0)
                {
                    tmp--;
                    sub.add_sub_inplace(input, false);
                }
                dividend.add_sub_inplace(sub, false, shift);
            }
            else if (dividend.abs_compare(input, shift) < 0)
            {
                shift--;
                double digit_3 = static_cast<double>(first_num1) * static_cast<double>(HINT_INT32_0X10);
                tmp = static_cast<UINT_64>(digit_3 / first_num2);
                sub = input * tmp;
                if (dividend.abs_compare(sub, shift) < 0)
                {
                    sub.add_sub_inplace(input, false);
                    tmp--;
                }
                dividend.add_sub_inplace(sub, false, shift);
            }
            else
            {
                dividend.add_sub_inplace(input, false, shift);
                tmp = 1;
            }
            result.data.array[shift] = static_cast<UINT_32>(tmp & HINT_INT32_0XFF);
        }
        result.neg_sign(static_cast<bool>(is_neg() ^ input.is_neg()));
        result.set_true_len();
        return result;
    }
    void base10000_in(UINT_16 *input_ary, size_t in_len, bool is_neg = false) //由10000进制的数转换，方便10进制字符串转换
    {
        const UINT_64 target_base = HINT_INT16_0X10;
        if (in_len == 1)
        {
            return;
        }
        size_t max_rank = 1ull << static_cast<UINT_16>(ceil(log2(in_len)) - 1);
        size_t unit_ary_num = static_cast<UINT_16>(log2(max_rank)) + 1; // 2^32,(2^32)^这样的个数
        size_t unit_ary_len = (_TWICE_(max_rank) - 1);
        UINT_16 *unit_ary = new UINT_16[unit_ary_len]; //用一个数组存储(10000)^1,(10000)^2,(10000)^4...
        unit_ary[0] = 10000;
        for (size_t i = 0, offset = 0; i < unit_ary_num - 1; i++)
        {
            size_t len = 1ull << i;
            hint::trans_square(unit_ary + offset, unit_ary + offset + len, len, target_base);
            offset += len;
        }
        size_t result_len = _TWICE_(max_rank);
        UINT_16 *result = input_ary;
        // hint::ary_copy(result, input_ary, in_len);
        // hint::ary_clr(result + in_len, result_len - in_len); //清零
        UINT_16 *tmp_product = new UINT_16[2];
        for (size_t i = 0, offset = 0; i < unit_ary_num; i++)
        {
            size_t len = 1ull << i;
            size_t gap = _TWICE_(len);
            size_t pos = 0;
            tmp_product = hint::ary_realloc(tmp_product, gap + 1);
            hint::ary_clr(tmp_product, gap);
            while (pos < result_len)
            {
                hint::trans_mul(unit_ary + offset, result + pos + len, tmp_product, len, len, target_base);
                hint::trans_add(tmp_product, result + pos, tmp_product, gap, len, target_base);
                hint::ary_copy(result + pos, tmp_product, gap);
                pos += gap;
            }
            offset += len;
        }
        // delete[] result;
        delete[] tmp_product;
        delete[] unit_ary;
    }
    void to_base100(UINT_8 *result, size_t &result_len) const //转为100进制的数组方便打印和转字符串
    {
        size_t in_len = length();
        size_t max_rank = 1ull << static_cast<UINT_16>(ceil(log2(in_len)) - 1);
        const UINT_16 base_len = 5;
        constexpr UINT_16 target_base = 100;                            // 2^32转为100进制的长度为5
        size_t unit_ary_num = static_cast<UINT_16>(log2(max_rank)) + 1; // 2^32,(2^32)^这样的个数

        result_len = _TWICE_(max_rank) * base_len;
        hint::ary_clr(result, result_len); //清零

        size_t unit_ary_len = (_TWICE_(max_rank) - 1) * base_len;
        // UINT_8 *unit_ary = new UINT_8[unit_ary_len]; //用一个数组存储(2^32)^1,(2^32)^2,(2^32)^4...
        UINT_8 *unit_ary = hint::to_base100_unit_ary;
        UINT_64 tmp = HINT_INT32_0X10; // 2^32
        size_t pos = 0;
        if (unit_ary_len > hint::to_base100_unit_ary_len)
        {
            hint::to_base100_unit_ary = hint::ary_realloc(hint::to_base100_unit_ary, unit_ary_len);
            unit_ary = hint::to_base100_unit_ary;
            if (hint::to_base100_unit_ary_len < 5)
            {
                while (pos < 5) //转100进制
                {
                    unit_ary[pos] = tmp % target_base;
                    tmp /= target_base;
                    pos++;
                }
                hint::to_base100_unit_ary_len = 5;
            }
            size_t index_begin = static_cast<UINT_16>(log2((hint::to_base100_unit_ary_len / 5) + 1) - 1);
            size_t offset_begin = ((1ull << index_begin) - 1) * 5;
            for (size_t i = index_begin, offset = offset_begin; i < unit_ary_num - 1; i++)
            {
                size_t len = (1ull << i) * base_len;
                hint::trans_square(unit_ary + offset, unit_ary + offset + len, len, target_base);
                offset += len;
            }
            hint::to_base100_unit_ary_len = unit_ary_len;
        }
        for (size_t i = 0; i < in_len; i++)
        {
            tmp = data.array[i];
            pos = 0;
            while (pos < 5)
            {
                result[i * base_len + pos] = tmp % target_base;
                tmp /= target_base;
                pos++;
            }
        }
        pos = 0;
        UINT_8 *tmp_product = new UINT_8[5];
        for (size_t i = 0, offset = 0; i < unit_ary_num; i++)
        {
            size_t len = (1ull << i) * base_len;
            size_t gap = _TWICE_(len);
            pos = 0;
            tmp_product = hint::ary_realloc(tmp_product, gap + 1);
            hint::ary_clr(tmp_product, gap);
            while (pos < result_len)
            {
                hint::trans_mul(unit_ary + offset, result + pos + len, tmp_product, len, len, target_base);
                hint::trans_add(tmp_product, result + pos, tmp_product, gap, len, target_base);
                hint::ary_copy(result + pos, tmp_product, gap);
                pos += gap;
            }
            // hint::to_base100_sub_prod(unit_ary, result, tmp_product, len, 0, offset, result_len);
            offset += len;
        }
        while (!result[result_len - 1])
        {
            result_len--;
        }
        delete[] tmp_product;
        // delete[] unit_ary;
    }

public:
    HyperInt fftsq() const
    {
        return fft_square();
    }
    HyperInt normsq() const
    {
        return normal_square();
    }
    HyperInt nttsq() const
    {
        return ntt_square();
    }
    HyperInt karsq() const
    {
        return karatsuba_square();
    }
    HyperInt test(const HyperInt &in) const // nor
    {
        return normal_multiply(in);
    }
    HyperInt test1(const HyperInt &in) const // fft
    {
        return fft_multiply(in);
    }
    HyperInt test2(const HyperInt &in) const // ntt
    {
        return ntt_multiply(in);
    }
    HyperInt test3(const HyperInt &in) const // kar
    {
        return newton_divide(in);
    }
    ~HyperInt() //析构函数
    {
        if (data.array != nullptr)
        {
            delete[] data.array;
            data.array = nullptr;
        }
    }
    HyperInt() //无参数构造
    {
        neg_sign(false);
        reset_size(2);
        change_length(0);
        data.array[1] = data.array[0] = 0;
    }
    HyperInt(size_t new_length, UINT_32 num) //填充new_length个num进行构造
    {
        data.size = generate_size(new_length);
        change_length(new_length);
        neg_sign(false);
        if (data.array != nullptr)
        {
            delete[] data.array;
        }
        data.array = new UINT_32[data.size];
        fill_element(num);
        set_true_len();
    }
    HyperInt(const HyperInt &input) // HyperInt 拷贝构造
    {
        if (this != &input)
        {
            size_t len = input.length();
            reset_size(len);
            change_length(len);
            neg_sign(input.is_neg());
            hint::ary_copy(data.array, input.data.array, len);
            set_true_len();
        }
    }
    HyperInt(HyperInt &&input) noexcept // HyperInt 移动构造
    {
        if (this != &input)
        {
            data.size = input.data.size;
            change_length(input.length());
            neg_sign(input.is_neg());
            delete[] data.array;
            data.array = input.data.array;
            input.data.array = nullptr;
            set_true_len();
        }
    }
    HyperInt(UINT_64 input) // UINT_64 参数构造
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(false);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = input >> HINT_INT_BIT;
        set_true_len();
    }
    HyperInt(INT_64 input) // INT_64 参数构造
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = input >> HINT_INT_BIT;
        set_true_len();
    }
    HyperInt(ULONG input) // unsigned long 参数构造
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(false);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = static_cast<UINT_64>(input) >> HINT_INT_BIT;
        set_true_len();
    }
    HyperInt(LONG input) // long 参数构造
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        INT_64 tmp = static_cast<INT_64>(input);
        neg_sign(_NEG_(tmp));
        tmp = _ABS_(tmp);
        data.array[0] = tmp & HINT_INT32_0XFF;
        data.array[1] = static_cast<INT_64>(tmp) >> HINT_INT_BIT;
        set_true_len();
    }
    HyperInt(UINT_32 input) // UINT_32 参数构造
    {
        data.size = 2;
        reset_size(2);
        change_length(1);
        neg_sign(false);
        data.array[0] = input;
        data.array[1] = 0;
        set_true_len();
    }
    HyperInt(INT_32 input) // INT_32 参数构造
    {
        data.size = 2;
        reset_size(2);
        change_length(1);
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        data.array[0] = input;
        data.array[1] = 0;
        set_true_len();
    }
    HyperInt(const std::string &input) // string 参数构造
    {
        string_in(input);
    }
    HyperInt(const char input[]) //字符串构造
    {
        string_in(input);
    }

    HyperInt &operator=(const HyperInt &input) // HyperInt 拷贝赋值
    {
        if (this != &input)
        {
            size_t len = input.length();
            reset_size(len);
            change_length(len);
            neg_sign(input.is_neg());
            hint::ary_copy(data.array, input.data.array, len);
            set_true_len();
        }
        return *this;
    }
    HyperInt &operator=(HyperInt &&input) noexcept // HyperInt 移动赋值
    {
        if (this != &input)
        {
            data.size = input.data.size;
            change_length(input.length());
            neg_sign(input.is_neg());
            delete[] data.array;
            data.array = input.data.array;
            input.data.array = nullptr;
            set_true_len();
        }
        return *this;
    }
    HyperInt &operator=(UINT_64 input) // UINT_64 赋值
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(false);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = input >> HINT_INT_BIT;
        set_true_len();
        return *this;
    }
    HyperInt &operator=(INT_64 input) // INT_64 赋值
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = input >> HINT_INT_BIT;
        set_true_len();
        return *this;
    }
    HyperInt &operator=(ULONG input) // unsigned long 赋值
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(false);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = static_cast<UINT_64>(input) >> HINT_INT_BIT;
        set_true_len();
        return *this;
    }
    HyperInt &operator=(LONG input) // long 赋值
    {
        data.size = 2;
        reset_size(2);
        change_length(2);
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        data.array[0] = input & HINT_INT32_0XFF;
        data.array[1] = static_cast<INT_64>(input) >> HINT_INT_BIT;
        set_true_len();
        return *this;
    }
    HyperInt &operator=(UINT_32 input) // UINT_32 赋值
    {
        data.size = 2;
        reset_size(2);
        change_length(1);
        neg_sign(false);
        data.array[0] = input;
        data.array[1] = 0;
        set_true_len();
        return *this;
    }
    HyperInt &operator=(INT_32 input) // INT_32 赋值
    {
        data.size = 2;
        reset_size(2);
        change_length(1);
        neg_sign(_NEG_(input));
        input = _ABS_(input);
        data.array[0] = input;
        data.array[1] = 0;
        set_true_len();
        return *this;
    }
    HyperInt &operator=(const std::string &input) // string 赋值
    {
        string_in(input);
        return *this;
    }
    HyperInt &operator=(const char input[]) //字符串赋值
    {
        string_in(input);
        return *this;
    }
    //基本操作
    void set_true_len() //去除前导0
    {
        size_t t_len = length();
        while (t_len > 0 && data.array[t_len - 1] == 0)
        {
            t_len--;
        }
        change_length(t_len);
    }
    void neg_sign(bool neg) //设置符号是否为负
    {
        if (neg)
        {
            data.neg_n_len = data.neg_n_len | HINT_INT64_0X80;
        }
        else
        {
            data.neg_n_len = data.neg_n_len & HINT_INT64_0X7F;
        }
    }
    INT_64 div_mod(INT_64 divisor) //自身除以divisor的同时返回余数
    {
        //  lldiv_t div_tmp;
        if (divisor == 0)
        {
            return to_int64();
        }
        UINT_64 last_rem = 0, tmp = 0, rem_num = 0;
        bool result_neg = is_neg() ^ _NEG_(divisor);
        neg_sign(result_neg);
        size_t pos = length();
        while (pos > 1)
        {
            pos--;
            tmp = (last_rem << HINT_INT_BIT) + data.array[pos];
            auto div_tmp = lldiv(tmp, divisor);                   //一次性得到商和余数
            data.array[pos] = static_cast<UINT_32>(div_tmp.quot); //当前数为变商
            last_rem = div_tmp.rem;                               //得到余数
        }
        tmp = (last_rem << HINT_INT_BIT) + data.array[0];
        auto div_tmp = lldiv(tmp, divisor);
        data.array[0] = static_cast<UINT_32>(div_tmp.quot);
        rem_num = div_tmp.rem;
        set_true_len();
        return rem_num;
    }
    INT_64 mod(INT_64 divisor) const //返回对divisor的余数
    {
        //  lldiv_t div_tmp;
        if (divisor == 0)
        {
            return to_int64();
        }
        INT_64 last_rem = 0, tmp = 0, rem_num = 0;
        size_t pos = length();
        while (pos > 1)
        {
            pos--;
            tmp = (last_rem << HINT_INT_BIT) + data.array[pos];
            auto div_tmp = lldiv(tmp, divisor); //一次性得到商和余数
            last_rem = div_tmp.rem;             //得到余数
        }
        tmp = (last_rem << HINT_INT_BIT) + data.array[0];
        auto div_tmp = lldiv(tmp, divisor);
        rem_num = div_tmp.rem;
        if (is_neg())
        {
            rem_num = -rem_num;
        }
        return rem_num;
    }
    HyperInt power(UINT_64 n) const //快速幂
    {
        HyperInt tmp(*this), result = 1;
        if (!_ODD_(n))
        {
            result.neg_sign(false);
        }
        else
        {
            result.neg_sign(is_neg());
        }
        if (abs_equal(result))
        {
            return result;
        }
        while (n)
        {
            if (n & 1)
            {
                result = tmp * result;
            }
            tmp = tmp.square();
            _SELFHALF_(n);
        }
        return result;
    }
    HyperInt square() const //求自身的平方
    {
        size_t len = length();
        if (len <= 48)
        {
            return normal_square();
        }
        else if (len <= HINT_QUAL_FFT_MAX / 2)
        {
            return fft_square();
        }
        else if (len < HINT_NTT_MAX / 16)
        {
            return karatsuba_square();
        }
        else if (len <= HINT_NTT_MAX / 2)
        {
            return ntt_square();
        }
        else
        {
            return karatsuba_square();
        }
    }
    HyperInt square_root() const
    {
        size_t len = length();
        if (len < 2)
        {
            return HyperInt(static_cast<UINT_32>(sqrt(first_int32())));
        }
        HyperInt left, right;
        left.reset_size(_HALF_(len) + 1);
        left.clear();
        left.change_length(_HALF_(len) + 1);
        left.data.array[_HALF_(len) - 1] = 1;
        left.set_true_len();
        right.reset_size(_HALF_(len) + 1);
        right.clear();
        right.change_length(_HALF_(len) + 1);
        right.data.array[_HALF_(len) - 1] = 2;
        right.set_true_len();
        while (!abs_smaller(right.square()))
        {
            // std::cout << left << "\n";
            // left.print_hex();
            // std::cin.get();
            right.quick_self_twice();
            left.quick_self_twice();
        }
        while (left.abs_smaller(right))
        {
            // std::cout <<"l"<< left << "\n";
            // std::cout <<"r"<< right << "\n";
            // std::cin.get();
            HyperInt mid = (left + right).half();
            if (left.abs_equal(mid))
            {
                return left;
            }
            INT_32 cmp = abs_compare(mid.square());
            if (cmp > 0)
            {
                left = std::move(mid);
            }
            else if (cmp < 0)
            {
                right = std::move(mid);
            }
            else
            {
                return mid;
            }
        }
        return left;
    }
    bool is_neg() const //返回符号是否为为负号
    {
        return data.neg_n_len & HINT_INT64_0X80;
    }
    template <typename T>
    static bool is_neg(T input)
    {
        return (input < 0);
    }
    size_t length() const //返回长度
    {
        return data.neg_n_len & HINT_INT64_0X7F;
    }
    size_t size() const //返回分配的数组空间
    {
        return data.size;
    }
    HyperInt split(size_t begin, size_t len) const //返回从下标begi开始长度为len的子数组
    {
        if (len == 0)
        {
            return HyperInt(0);
        }
        size_t data_len = length();
        if (begin >= data_len)
        {
            INT_64 num = data.array[data_len - 1];
            return HyperInt(num);
        }
        else if (len > data_len - begin)
        {
            len = data_len - begin;
        }
        HyperInt result;
        result.reset_size(len);
        result.change_length(len);
        hint::ary_copy(result.data.array, data.array + begin, len);
        result.set_true_len();
        result.neg_sign(is_neg());
        return result;
    }
    INT_64 to_int64() const //转INT_64
    {
        INT_64 out = 0;
        out = (data.array[1] & HINT_INT32_0X7F);
        out <<= HINT_INT_BIT;
        out += data.array[0];
        if (is_neg())
        {
            out = -out;
        }
        return out;
    }
    UINT_64 to_uint64() const //转UINT_64
    {
        UINT_64 out = 0;
        out = data.array[1];
        out <<= HINT_INT_BIT;
        out += data.array[0];
        return out;
    }
    std::string to_string() const //转string,用10进制表示的字符串
    {
        if (equal_to_z())
        {
            return std::string("0");
        }
        size_t in_len = length();
        if (in_len == 1)
        {
            if (is_neg())
            {
                return std::string("-") + std::to_string(data.array[0]);
            }
            return std::to_string(data.array[0]);
        }
        if (in_len < 128)
        {
            char digits[3] = {};
            std::string result_str;
            HyperInt tmp(*this);
            size_t pos = 0;
            INT_64 rem_n = 0;
            HyperInt zero;
            while (tmp.abs_larger(zero))
            {
                rem_n = tmp.div_mod(100);
                digits[0] = static_cast<char>(rem_n % 10) + '0';
                digits[1] = static_cast<char>(rem_n / 10) + '0';
                result_str += digits;
                pos++;
            }
            if (result_str[result_str.size() - 1] == '0')
            {
                result_str.erase(result_str.end() - 1);
            }
            if (is_neg())
            {
                result_str += '-';
            }
            std::reverse(result_str.begin(), result_str.end());
            return result_str;
        }
        size_t max_rank = 1ull << static_cast<UINT_16>(ceil(log2(in_len)) - 1);
        const UINT_16 base_len = 5; // 2^32转为100进制的长度为5
        size_t result_len = _TWICE_(max_rank) * base_len;
        UINT_8 *result = new UINT_8[result_len]; //结果数组
        to_base100(result, result_len);
        char digits[3] = {};
        std::string result_str;
        if (is_neg())
        {
            result_str += '-';
        }
        result_str += std::to_string(result[result_len - 1]);
        result_len--;
        while (result_len > 0) //整理每一位将100进制转为十进制
        {
            result_len--;
            digits[0] = result[result_len] / 10 + '0';
            digits[1] = result[result_len] % 10 + '0';
            result_str += digits;
        }
        delete[] result;
        return result_str;
    }
    void string_in(const std::string &str)
    {
        size_t len = str.size();
        if (len > 1500 && len <= 39000)
        {
            quick_string_in(str);
        }
        else
        {
            normal_string_in(str);
        }
    }
    void string_in(const char str[])
    {
        size_t len = strlen(str);
        if (len > 1500 && len <= 39000)
        {
            quick_string_in(str);
        }
        else
        {
            normal_string_in(str);
        }
    }
    void normal_string_in(const std::string &str) //输入十进制字符串
    {
        clear();
        for (auto &i : str)
        {
            if ('0' <= i && i <= '9')
            {
                *this *= static_cast<UINT_64>(10);
                *this += (i - '0');
            }
        }
        if (str[0] == '-')
        {
            neg_sign(true);
        }
    }
    void normal_string_in(const char str[]) //输入十进制字符串
    {
        clear();
        size_t in_len = strlen(str);
        for (size_t i = 0; i < in_len; i++)
        {
            char c = str[i];
            if ('0' <= c && c <= '9')
            {
                *this *= static_cast<UINT_64>(10);
                *this += (c - '0');
            }
        }
        if (str[0] == '-')
        {
            neg_sign(true);
        }
    }
    void quick_string_in(const std::string &str)
    {
        size_t in_len = str.size();
        size_t trans_len = in_len / 4 + 1;
        size_t ary_len = 1ull << static_cast<UINT_16>(ceil(log2(trans_len)));
        UINT_16 *trans_ary = new UINT_16[ary_len];
        hint::ary_clr(trans_ary, ary_len);
        bool is_neg = (str[0] == '-');
        for (size_t pos = 0; pos < trans_len; pos++)
        {
            UINT_16 tmp = 0, rank = 1;
            for (size_t i = 0; i < 4 && in_len > 0; i++)
            {
                in_len--;
                char c = str[in_len];
                if ('0' <= c && c <= '9')
                {
                    tmp += rank * (c - '0');
                }
                rank *= 10;
            }
            trans_ary[pos] = tmp;
        }
        base10000_in(trans_ary, trans_len, is_neg);

        reset_size((trans_len + 1) / 2);
        change_length((trans_len + 1) / 2);
        neg_sign(is_neg);
        memcpy(data.array, trans_ary, trans_len * sizeof(*trans_ary));
        set_true_len();
        delete[] trans_ary;
    }
    void quick_string_in(const char str[])
    {
        size_t in_len = strlen(str);
        size_t trans_len = in_len / 4 + 1;
        size_t ary_len = 1ull << static_cast<UINT_16>(ceil(log2(trans_len)));
        UINT_16 *trans_ary = new UINT_16[ary_len];
        hint::ary_clr(trans_ary, ary_len);
        bool is_neg = (str[0] == '-');
        for (size_t pos = 0; pos < trans_len; pos++)
        {
            UINT_16 tmp = 0, rank = 1;
            for (size_t i = 0; i < 4 && in_len > 0; i++)
            {
                in_len--;
                char c = str[in_len];
                if ('0' <= c && c <= '9')
                {
                    tmp += rank * (c - '0');
                }
                rank *= 10;
            }
            trans_ary[pos] = tmp;
        }
        base10000_in(trans_ary, trans_len, is_neg);

        reset_size((trans_len + 1) / 2);
        change_length((trans_len + 1) / 2);
        neg_sign(is_neg);
        memcpy(data.array, trans_ary, trans_len * sizeof(*trans_ary));
        set_true_len();
        delete[] trans_ary;
    }
    void console_in() //从控制台读入十进制值
    {
        clear();
        char tmp = '0';
        bool head = true, neg = false;
        while (tmp != '\n' && '0' <= tmp && tmp <= '9')
        {
            if (tmp != '\n')
            {
                *this *= static_cast<UINT_64>(10);
                *this += (tmp - '0');
            }
            tmp = getchar();
            if (head)
            {
                if (tmp == '-')
                {
                    neg = true;
                }
                head = false;
            }
        }
        neg_sign(neg);
    }
    void print_hex() const //向控制台打印十六进制值
    {
        if (equal_to_z())
        {
            printf("0");
        }
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
    void classic_print() const //向控制台打印十进制值,传统算法
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
        while (tmp.abs_larger(HyperInt()))
        {
            rem_n = tmp.div_mod(100);
            tmp_str[pos] = static_cast<UINT_8>(rem_n);
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
    void quick_print() const //向控制台打印十进制值,快速算法
    {
        if (equal_to_z())
        {
            putchar('0');
            return;
        }
        if (is_neg())
        {
            putchar('-');
        }
        if (length() == 1)
        {
            printf("%u", data.array[0]);
            return;
        }
        size_t max_rank = 1ull << static_cast<UINT_16>(ceil(log2(length())) - 1);
        const UINT_16 base_len = 5; // 2^32转为100进制的长度为5
        size_t result_len = _TWICE_(max_rank) * base_len;
        UINT_8 *result = new UINT_8[result_len]; //结果数组
        to_base100(result, result_len);
        printf("%d", result[result_len - 1]);
        result_len--;
        while (result_len > 0)
        {
            result_len--;
            printf("%02d", result[result_len]);
        }
        delete[] result;
    }
    void print_dec() const //向控制台打印十进制值,取两种算法最快者
    {
        if (length() < 3300)
        {
            classic_print();
        }
        else
        {
            quick_print();
        }
    }
    HyperInt add_sub(const HyperInt &input, bool is_add) const //基础加减法a=b.add_sub(c,ture)->a=b+c;a=b.add_sub(c,fasle)->a=b-c,(b>c);
    {
        HyperInt result;
        size_t len1 = length(), len2 = input.length();
        if (is_add)
        {
            size_t result_len = _MAX_(len1, len2) + 1;
            result.reset_size(result_len);
            result.change_length(result_len);
            INT_64 tmp = 0;
            size_t count = 0;
            while (count < result_len)
            {
                if (count < len1)
                {
                    tmp += data.array[count];
                }
                if (count < len2)
                {
                    tmp += input.data.array[count];
                }
                result.data.array[count] = static_cast<UINT_32>(tmp & HINT_INT32_0XFF);
                tmp >>= HINT_INT_BIT;
                count++;
            }
        }
        else
        {
            size_t result_len = _MAX_(len1, len2);
            result.reset_size(result_len);
            result.change_length(result_len);
            INT_64 tmp = 0;
            size_t count = 0;
            while (count < result_len)
            {
                if (count < len1)
                {
                    tmp += data.array[count];
                }
                if (count < len2)
                {
                    tmp -= input.data.array[count];
                }
                result.data.array[count] = static_cast<UINT_32>(tmp + HINT_INT32_0X10);
                tmp >>= HINT_INT_BIT;
                count++;
            }
        }
        result.set_true_len();
        return result;
    }
    void add_sub_inplace(const HyperInt &input, bool is_add, const size_t shift = 0) //就地加减 a+=b;a-=b,a加/减去左移位后的b，默认不移位
    {
        size_t len1 = length(), len2 = input.length();
        if (is_add)
        {
            size_t result_len = _MAX_(len1, len2) + 1;
            reset_size(result_len);
            hint::ary_clr(data.array + len1, result_len - len1);
            change_length(result_len);
            INT_64 tmp = 0;
            size_t count = shift;
            while (count < result_len)
            {
                UINT_32 first_num = 0;
                if (count < len1)
                {
                    first_num = data.array[count];
                    tmp += first_num;
                }
                if (count - shift < len2)
                {
                    tmp += input.data.array[count - shift];
                }
                else if (tmp == static_cast<INT_64>(first_num))
                {
                    break;
                }
                data.array[count] = static_cast<UINT_32>(tmp & HINT_INT32_0XFF);
                tmp >>= HINT_INT_BIT;
                count++;
            }
        }
        else
        {

            size_t result_len = _MAX_(len1, len2);
            reset_size(result_len);
            change_length(result_len);
            INT_64 tmp = 0;
            size_t count = shift;
            while (count < result_len)
            {
                UINT_32 first_num = 0;
                if (count < len1)
                {
                    first_num = data.array[count];
                    tmp += first_num;
                }
                if (count - shift < len2)
                {
                    tmp -= input.data.array[count - shift];
                }
                else if (tmp == static_cast<INT_64>(first_num))
                {
                    break;
                }
                data.array[count] = static_cast<UINT_32>(tmp + HINT_INT32_0X10);
                tmp >>= HINT_INT_BIT;
                count++;
            }
        }
        set_true_len();
    }
    void sub_inplace(const HyperInt &input) //由减数调用,a.sub_inplace(b)->a=b-a;
    {
        size_t len1 = length(), len2 = input.length();
        size_t result_len = _MAX_(len1, len2);
        reset_size(result_len);
        change_length(result_len);
        INT_64 tmp = 0;
        size_t count = 0;
        while (count < result_len)
        {
            if (count < len2)
            {
                tmp += input.data.array[count];
            }
            if (count < len1)
            {
                tmp -= data.array[count];
            }
            data.array[count] = static_cast<UINT_32>(tmp + HINT_INT32_0X10);
            tmp >>= HINT_INT_BIT;
            count++;
        }
        set_true_len();
    }
    HyperInt &self_half() //自身右移一位
    {
        UINT_32 tmp1, tmp2 = 0;
        size_t pos = length();
        while (pos)
        {
            pos--;
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (tmp1 << (HINT_INT_BIT - 1)) + (tmp2 >> 1);
            data.array[pos] = tmp1;
        }
        set_true_len();
        return *this;
    }
    HyperInt &self_twice() //自身左移一位
    {
        UINT_64 tmp = 0;
        size_t len = length();
        for (size_t pos = 0; pos < len; pos++)
        {
            tmp += (static_cast<UINT_64>(data.array[pos]) << 1);
            data.array[pos] = static_cast<UINT_32>(tmp);
            tmp >>= HINT_INT_BIT;
        }
        if (tmp > 0)
        {
            len++;
            reset_size(len);
            change_length(len);
            data.array[len - 1] = HINT_INT32_0X01;
        }
        return *this;
    }
    HyperInt half() const //返回右移后的值
    {
        return r_shift(1);
    }
    HyperInt twice() const //返回左移后的值
    {
        return l_shift(1);
    }
    HyperInt r_shift(size_t n) const //右移n位
    {
        if (n == 0)
        {
            return *this;
        }
        HyperInt result;
        size_t shift = n / 32;
        size_t offset = n % 32;
        size_t len = length();
        shift %= len;
        len -= shift;
        UINT_32 tmp1, tmp2 = 0;
        result.reset_size(len);
        result.change_length(len);
        while (len > 0)
        {
            len--;
            tmp1 = tmp2;
            tmp2 = data.array[len + shift];
            tmp1 = (offset ? (tmp1 << (HINT_INT_BIT - offset)) : 0) + (tmp2 >> offset);
            result.data.array[len] = tmp1;
        }
        result.set_true_len();
        return result;
    }
    HyperInt l_shift(size_t n) const //左移n位
    {
        if (n == 0)
        {
            return *this;
        }
        HyperInt result;
        size_t shift = n / 32;
        size_t offset = n % 32;
        size_t len = length();
        UINT_32 tmp1 = 0, tmp2 = 0;
        result.reset_size(len + shift + 1);
        result.clear();
        result.change_length(len + shift + 1);
        for (size_t pos = 0; pos < length(); pos++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[pos];
            tmp1 = (offset ? (tmp1 >> (HINT_INT_BIT - offset)) : 0) + (tmp2 << offset);
            result.data.array[pos + shift] = tmp1;
        }
        if (tmp2 >> (HINT_INT_BIT - offset) & offset)
        {
            result.data.array[len + shift] = tmp2 >> (HINT_INT_BIT - offset);
        }
        else
        {
            result.change_length(result.length() - 1);
        }
        return result;
    }
    void reset_size(size_t new_size_input) //重新设定长度不小于new_size,1.5倍长度算法,在change_len()之前调用
    {
        size_t size_tmp = generate_size(new_size_input);
        if (data.array == nullptr)
        {
            data.size = size_tmp;
            data.array = new UINT_32[data.size];
            change_length(0);
        }
        else if (data.size != size_tmp)
        {
            data.size = size_tmp;
            data.array = hint::ary_realloc(data.array, data.size);
            change_length(_MIN_(length(), data.size));
        }
    }
    HyperInt abs() const //返回绝对值
    {
        HyperInt result(*this);
        result.neg_sign(false);
        return result;
    }
    UINT_32 first_int32() const
    {
        if (length() >= 1)
        {
            return data.array[length() - 1];
        }
        else
        {
            return 0;
        }
    }
    UINT_64 first_int64() const
    {
        if (length() >= 2)
        {
            UINT_64 result = data.array[length() - 1];
            result <<= HINT_INT_BIT;
            result += data.array[length() - 2];
            return result;
        }
        else
        {
            return static_cast<UINT_64>(first_int32());
        }
    }
    INT_32 abs_compare(const HyperInt &input, INT_64 shift = 0) const //自身和input移位shift比较，大于返回1，小于返回-1，等于返回0
    {
#if SIZE_T_BITS == 64
        INT_64 len1 = static_cast<INT_64>(length());
        INT_64 len2 = static_cast<INT_64>(input.length()) + shift;
#elif SIZE_T_BITS == 32
        INT_32 len1 = static_cast<INT_32>(length());
        INT_32 len2 = static_cast<INT_32>(input.length()) + shift;
#else
#error "unknown"
#endif
        if (len1 < len2)
        {
            return -1;
        }
        else if (len1 > len2)
        {
            return 1;
        }
        UINT_32 num1 = 0, num2 = 0;
        size_t len = len1;
        while (len > 0)
        {
            len--;
            num1 = data.array[len];
            num2 = input.data.array[len - shift];
            if (num1 > num2)
            {
                return 1;
            }
            else if (num1 < num2)
            {
                return -1;
            }
        }
        return 0;
    }
    bool abs_larger(const HyperInt &input) const //绝对值是否大于input
    {
        size_t t_len1 = length(), t_len2 = input.length();
        if (t_len1 > t_len2)
        {
            return true;
        }
        else if (t_len1 < t_len2)
        {
            return false;
        }

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
    bool abs_smaller(const HyperInt &input) const //绝对值是否小于input
    {
        size_t t_len1 = length(), t_len2 = input.length();
        if (t_len1 < t_len2)
        {
            return true;
        }
        else if (t_len1 > t_len2)
        {
            return false;
        }
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
    bool abs_equal(const HyperInt &input) const //绝对值是否等于input
    {
        if (this == &input)
        {
            return true;
        }
        size_t t_len1 = length(), t_len2 = input.length();
        if (t_len1 != t_len2)
        {
            return false;
        }
        while (t_len1 > 0)
        {
            t_len1--;
            if (data.array[t_len1] != input.data.array[t_len1])
            {
                return false;
            }
        }
        return true;
    }
    bool equal_to_z() const //判定是否为零
    {
        return (0 == length());
    }
    bool is_even() const //判断是否为偶数
    {
        if (length() == 0)
        {
            return true;
        }
        return !static_cast<bool>(data.array[0] & 1);
    }
    bool is_odd() const //判断是否为奇数
    {
        if (length() == 0)
        {
            return false;
        }
        return (bool)data.array[0] & 1;
    }
    //逻辑运算
    template <typename T>
    bool operator==(const T &input) const;
    template <typename T>
    bool operator!=(const T &input) const;
    template <typename T>
    bool operator>(const T &input) const;
    template <typename T>
    bool operator>=(const T &input) const;
    template <typename T>
    bool operator<(const T &input) const;
    template <typename T>
    bool operator<=(const T &input) const;
    bool operator!() const;
    // operator bool();

    //友元函数
    friend HyperInt abs(const HyperInt &input);
    friend void print(const HyperInt &input);
    friend bool operator==(const INT_64 &input1, const HyperInt &input2);
    friend bool operator!=(const INT_64 &input1, const HyperInt &input2);
    friend bool operator>(const INT_64 &input1, const HyperInt &input2);
    friend bool operator>=(const INT_64 &input1, const HyperInt &input2);
    friend bool operator<(const INT_64 &input1, const HyperInt &input2);
    friend bool operator<=(const INT_64 &input1, const HyperInt &input2);

    friend HyperInt operator+(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator-(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator*(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator/(const INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator%(const INT_64 &input1, const HyperInt &input2);

    friend INT_64 &operator+=(INT_64 &input1, const HyperInt &input2);
    friend INT_64 &operator-=(INT_64 &input1, const HyperInt &input2);
    friend INT_64 &operator*=(INT_64 &input1, const HyperInt &input2);
    friend INT_64 &operator/=(INT_64 &input1, const HyperInt &input2);
    friend INT_64 &operator%=(INT_64 &input1, const HyperInt &input2);

    friend std::string to_string(const HyperInt &input);
    friend std::ostream &operator<<(std::ostream &output, const HyperInt &input);
    friend std::istream &operator>>(std::istream &input, HyperInt &output);

    friend HyperInt combination(UINT_64 n, UINT_64 m);

    //算术运算
    HyperInt operator+(const HyperInt &input) const;
    HyperInt operator+(INT_64 input) const;
    HyperInt operator+() const;

    HyperInt operator-(const HyperInt &input) const;
    HyperInt operator-(INT_64 input) const;
    HyperInt operator-() const;

    HyperInt operator*(const HyperInt &input) const;
    HyperInt operator*(UINT_64 input) const;
    HyperInt operator*(INT_64 input) const;

    HyperInt operator/(const HyperInt &input) const;
    HyperInt operator/(INT_64 input) const;

    HyperInt operator%(const HyperInt &input) const;
    INT_64 operator%(INT_64 input) const;

    HyperInt &operator+=(const HyperInt &input);

    HyperInt &operator-=(const HyperInt &input);

    HyperInt &operator*=(const HyperInt &input);
    HyperInt &operator*=(UINT_64 input);
    HyperInt &operator*=(INT_64 input);

    HyperInt &operator/=(const HyperInt &input);
    HyperInt &operator/=(INT_64 input);

    HyperInt &operator%=(const HyperInt &input);
    HyperInt &operator%=(INT_64 input);

    HyperInt operator++(int);
    HyperInt &operator++();
    HyperInt operator--(int);
    HyperInt &operator--();
    HyperInt operator~() const;
    HyperInt operator|(const HyperInt &input) const;
    HyperInt operator&(const HyperInt &input) const;
    HyperInt operator^(const HyperInt &input) const;
    HyperInt &operator|=(const HyperInt &input);
    HyperInt &operator&=(const HyperInt &input);
    HyperInt &operator^=(const HyperInt &input);
};
//逻辑运算

template <typename T>
inline bool HyperInt::operator==(const T &input) const
{
    if (is_neg() ^ is_neg(input))
    {
        return false;
    }
    else
    {
        return abs_equal(input);
    }
}

template <typename T>
inline bool HyperInt::operator!=(const T &input) const
{
    if (is_neg() ^ is_neg(input))
    {
        return true;
    }
    else
    {
        return !abs_equal(input);
    }
}

template <typename T>
inline bool HyperInt::operator>(const T &input) const
{
    if (is_neg() ^ is_neg(input))
    {
        return is_neg(input);
    }
    else
    {
        return is_neg() ^ abs_larger(input);
    }
}

template <typename T>
inline bool HyperInt::operator>=(const T &input) const
{
    return !(*this < input);
}

template <typename T>
inline bool HyperInt::operator<(const T &input) const
{
    if (is_neg() ^ is_neg(input))
    {
        return is_neg();
    }
    else
    {
        return is_neg() ^ abs_smaller(input);
    }
}

template <typename T>
inline bool HyperInt::operator<=(const T &input) const
{
    return !(*this > input);
}

inline bool HyperInt::operator!() const
{
    return !(bool)length();
}
// inline HyperInt::operator bool();

//算术运算
HyperInt HyperInt::operator+(const HyperInt &input) const
{
    HyperInt result;
    if (!is_neg() ^ input.is_neg()) //是否同号
    {
        result = add_sub(input, true);
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
            result = add_sub(input, false);
            result.neg_sign(is_neg());
            return result;
        }
        else
        {
            result = input.add_sub(*this, false);
            result.neg_sign(!is_neg());
            return result;
        }
    }
}
HyperInt HyperInt::operator+(INT_64 input) const
{
    return *this + HyperInt(input);
}
inline HyperInt HyperInt::operator+() const
{
    return *this;
}
HyperInt HyperInt::operator-(const HyperInt &input) const
{
    HyperInt result;
    if (is_neg() ^ input.is_neg()) //是否异号
    {
        result = add_sub(input, true);
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
            result = add_sub(input, false);
            result.neg_sign(is_neg());
            return result;
        }
        else
        {
            result = input.add_sub(*this, false);
            result.neg_sign(!is_neg());
            return result;
        }
    }
}
HyperInt HyperInt::operator-(INT_64 input) const
{
    return *this - HyperInt(input);
}
inline HyperInt HyperInt::operator-() const
{
    HyperInt result(*this);
    result.neg_sign(!is_neg());
    return result;
}
inline HyperInt HyperInt::operator*(const HyperInt &input) const
{
    if (this == &input)
    {
        return square();
    }
    else
    {
        size_t min_len = _MIN_(length(), input.length());
        size_t sum_len = length() + input.length();
        size_t factor = 4;
        if (sum_len <= HINT_FFT_MAX)
        {
            factor = 4;
        }
        else if (sum_len <= HINT_QUAL_FFT_MAX)
        {
            factor = 9;
        }
        else
        {
            factor = 43;
        }
        if (sum_len <= 120 || (min_len <= factor * std::ceil(std::log2(sum_len))))
        {
            // printf("normal\n");
            return normal_multiply(input);
        }
        else if (sum_len <= HINT_QUAL_FFT_MAX)
        {
            // printf("fft\n");
            return fft_multiply(input);
        }
        else if (sum_len < HINT_NTT_MAX / 8)
        {
            return karatsuba_multiply(input);
        }
        else if (sum_len <= HINT_NTT_MAX)
        {
            // printf("ntt\n");
            return ntt_multiply(input);
        }
        else
        {
            // printf("over\n");
            return karatsuba_multiply(input);
        }
    }
}
HyperInt HyperInt::operator*(UINT_64 input) const
{
    if (input == 1)
    {
        return *this;
    }
    HyperInt result;
    if (equal_to_z() || input == 0)
    {
        return result;
    }
    size_t len = length();
    size_t result_len = len + 2;
    result.reset_size(result_len);
    result.change_length(result_len);
    result.neg_sign(is_neg());
    UINT_32 tmp = 0, tmp1 = 0, tmp2 = 0;
    UINT_64 sum1 = 0, sum2 = 0;

    const UINT_64 input_num1 = input & HINT_INT32_0XFF;
    const UINT_64 input_num2 = input >> HINT_INT_BIT;

    tmp2 = data.array[0];
    sum1 = tmp2 * input_num1;
    result.data.array[0] = sum1 & HINT_INT32_0XFF;
    sum1 >>= HINT_INT_BIT;
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
        sum1 += input_num2 * tmp1;
        tmp = sum1 & HINT_INT32_0XFF;
        sum2 += input_num1 * tmp2 + tmp;
        tmp = sum2 & HINT_INT32_0XFF;
        result.data.array[pos] = tmp;
        sum1 >>= HINT_INT_BIT;
        sum2 >>= HINT_INT_BIT;
    }
    result.set_true_len();
    return result;
}
inline HyperInt HyperInt::operator*(INT_64 input) const
{
    bool result_neg = is_neg() ^ _NEG_(input);
    input = std::abs(input);
    HyperInt result;
    if (equal_to_z() || input == 0)
    {
        return result;
    }
    result = *this * (static_cast<UINT_64>(input));
    result.neg_sign(result_neg);
    return result;
}
HyperInt HyperInt::operator/(const HyperInt &input) const
{
    if (this == &input || abs_equal(input))
    {
        return HyperInt(1);
    }
    size_t len1 = length(), len2 = input.length();
    if (len1 <= 4096 || len2 <= 4096)
    {
        return normal_divide(input);
    }
    else
    {
        return newton_divide(input);
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
    if (!is_neg() ^ input.is_neg()) //是否同号
    {
        add_sub_inplace(input, true);
        neg_sign(is_neg());
    }
    else
    {
        if (abs_equal(input))
        {
            *this = HyperInt(0);
        }
        else if (abs_larger(input))
        {
            add_sub_inplace(input, false);
            neg_sign(is_neg());
        }
        else
        {
            sub_inplace(input);
            neg_sign(!is_neg());
        }
    }
    return *this;
}
HyperInt &HyperInt::operator-=(const HyperInt &input)
{
    if (is_neg() ^ input.is_neg()) //是否异号
    {
        add_sub_inplace(input, true);
        neg_sign(is_neg());
    }
    else
    {
        if (abs_equal(input))
        {
            *this = HyperInt(0);
        }
        else if (abs_larger(input))
        {
            add_sub_inplace(input, false);
            neg_sign(is_neg());
        }
        else
        {
            sub_inplace(input);
            neg_sign(!is_neg());
        }
    }
    return *this;
}
HyperInt &HyperInt::operator*=(const HyperInt &input)
{
    if (this == &input)
    {
        *this = square();
    }
    else
    {
        *this = *this * input;
    }
    return *this;
}
HyperInt &HyperInt::operator*=(INT_64 input)
{
    bool result_neg = is_neg() ^ _NEG_(input);
    input = std::abs(input);
    if (equal_to_z() || input == 0)
    {
        *this = 0;
        return *this;
    }
    *this *= (static_cast<UINT_64>(input));
    neg_sign(result_neg);
    return *this;
}
HyperInt &HyperInt::operator*=(UINT_64 input)
{
    if (input == 1)
    {
        return *this;
    }
    if (equal_to_z() || input == 0)
    {
        *this = 0;
        return *this;
    }
    size_t len = length();
    size_t result_len = len + 2;
    reset_size(result_len);
    change_length(result_len);
    UINT_32 tmp = 0, tmp1 = 0, tmp2 = 0;
    UINT_64 sum1 = 0, sum2 = 0;

    const UINT_64 input_num1 = input & HINT_INT32_0XFF;
    const UINT_64 input_num2 = input >> HINT_INT_BIT;

    tmp2 = data.array[0];
    sum1 = tmp2 * input_num1;
    data.array[0] = sum1 & HINT_INT32_0XFF;
    sum1 >>= HINT_INT_BIT;
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
        sum1 += input_num2 * tmp1;
        tmp = sum1 & HINT_INT32_0XFF;
        sum2 += input_num1 * tmp2 + tmp;
        tmp = sum2 & HINT_INT32_0XFF;
        data.array[pos] = tmp;
        sum1 >>= HINT_INT_BIT;
        sum2 >>= HINT_INT_BIT;
    }
    set_true_len();
    return *this;
}
HyperInt &HyperInt::operator/=(const HyperInt &input)
{
    *this = *this / input;
    return *this;
}
HyperInt &HyperInt::operator/=(INT_64 input)
{
    div_mod(input);
    return *this;
}
HyperInt &HyperInt::operator%=(const HyperInt &input)
{
    HyperInt tmp = *this / input;
    *this = *this - (tmp * input);
    return *this;
}
HyperInt &HyperInt::operator%=(INT_64 input)
{
    *this = HyperInt(div_mod(input));
    return *this;
}
HyperInt HyperInt::operator++(int)
{
    HyperInt tmp(*this);
    *this += 1;
    return tmp;
}
HyperInt &HyperInt::operator++()
{
    *this += 1;
    return *this;
}
HyperInt HyperInt::operator--(int)
{
    HyperInt tmp(*this);
    *this -= 1;
    return tmp;
}
HyperInt &HyperInt::operator--()
{
    *this -= 1;
    return *this;
}
HyperInt HyperInt::operator~() const
{
    size_t len = length();
    HyperInt result;
    result.reset_size(len);
    result.change_length(len);
    for (size_t i = 0; i < len; i++)
    {
        result.data.array[i] = ~data.array[i];
    }
    result.set_true_len();
    return result;
}
HyperInt HyperInt::operator|(const HyperInt &input) const
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = _MAX_(len1, len2);
    size_t min_len = _MIN_(len1, len2);
    HyperInt result;
    result.reset_size(result_len);
    result.clear();
    result.change_length(result_len);
    size_t pos = 0;
    while (pos < min_len)
    {
        result.data.array[pos] = data.array[pos] | input.data.array[pos];
        pos++;
    }
    if (len1 > len2)
    {
        hint::ary_copy(result.data.array + min_len, data.array + min_len, len1 - len2);
    }
    else
    {
        hint::ary_copy(result.data.array + min_len, input.data.array + min_len, len2 - len1);
    }
    result.set_true_len();
    return result;
}
HyperInt HyperInt::operator&(const HyperInt &input) const
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = _MIN_(len1, len2);
    HyperInt result;
    result.reset_size(result_len);
    result.clear();
    result.change_length(result_len);
    size_t pos = 0;
    while (pos < result_len)
    {
        result.data.array[pos] = data.array[pos] & input.data.array[pos];
        pos++;
    }
    result.set_true_len();
    return result;
}
HyperInt HyperInt::operator^(const HyperInt &input) const
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = _MAX_(len1, len2);
    size_t min_len = _MIN_(len1, len2);
    HyperInt result;
    result.reset_size(result_len);
    result.clear();
    result.change_length(result_len);
    size_t pos = 0;
    while (pos < min_len)
    {
        result.data.array[pos] = data.array[pos] ^ input.data.array[pos];
        pos++;
    }
    if (len1 > len2)
    {
        hint::ary_copy(result.data.array + min_len, data.array + min_len, len1 - len2);
    }
    else
    {
        hint::ary_copy(result.data.array + min_len, input.data.array + min_len, len2 - len1);
    }
    result.set_true_len();
    return result;
}
HyperInt &HyperInt::operator|=(const HyperInt &input)
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = _MAX_(len1, len2);
    size_t min_len = _MIN_(len1, len2);
    reset_size(result_len);
    change_length(result_len);
    size_t pos = 0;
    while (pos < min_len)
    {
        data.array[pos] |= input.data.array[pos];
        pos++;
    }
    if (len1 < len2)
    {
        hint::ary_copy(data.array + len1, input.data.array + len1, len2 - len1);
    }
    set_true_len();
    return *this;
}
HyperInt &HyperInt::operator&=(const HyperInt &input)
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = _MIN_(len1, len2);
    if (len1 > len2)
    {
        reset_size(result_len);
        change_length(result_len);
    }
    size_t pos = 0;
    while (pos < result_len)
    {
        data.array[pos] &= input.data.array[pos];
        pos++;
    }
    set_true_len();
    return *this;
}
HyperInt &HyperInt::operator^=(const HyperInt &input)
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = _MAX_(len1, len2);
    size_t min_len = _MIN_(len1, len2);
    reset_size(result_len);
    change_length(result_len);
    size_t pos = 0;
    while (pos < min_len)
    {
        data.array[pos] ^= input.data.array[pos];
        pos++;
    }
    if (len1 < len2)
    {
        hint::ary_copy(data.array + min_len, input.data.array + min_len, len2 - len1);
    }
    set_true_len();
    return *this;
}
//友元函数
HyperInt abs(const HyperInt &input) //返回input的绝对值
{
    HyperInt result(input);
    result.neg_sign(false);
    return result;
}
void print(const HyperInt &input) //打印input
{
    if (input.length() <= 1000000)
    {
        input.print_dec();
    }
    else
    {
        input.print_hex();
    }
}
bool operator!=(const INT_64 &input1, const HyperInt &input2)
{
    return input2 != input1;
}
bool operator==(const INT_64 &input1, const HyperInt &input2)
{
    return input2 == input1;
}
bool operator>(const INT_64 &input1, const HyperInt &input2)
{
    return input2 < input1;
}
bool operator>=(const INT_64 &input1, const HyperInt &input2)
{
    return input2 <= input1;
}
bool operator<(const INT_64 &input1, const HyperInt &input2)
{
    return input2 > input1;
}
bool operator<=(const INT_64 &input1, const HyperInt &input2)
{
    return input2 >= input1;
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
    return HyperInt(input1) / input2;
}
HyperInt operator%(const INT_64 &input1, const HyperInt &input2)
{
    return HyperInt(input1) % input2;
}
INT_64 &operator+=(INT_64 &input1, const HyperInt &input2)
{
    return input1 += input2.to_int64();
}
INT_64 &operator-=(INT_64 &input1, const HyperInt &input2)
{
    return input1 -= input2.to_int64();
}
INT_64 &operator*=(INT_64 &input1, const HyperInt &input2)
{
    return input1 *= input2.to_int64();
}
INT_64 &operator/=(INT_64 &input1, const HyperInt &input2)
{
    return input1 /= input2.to_int64();
}
INT_64 &operator%=(INT_64 &input1, const HyperInt &input2)
{
    return input1 %= input2.to_int64();
}
inline std::string to_string(const HyperInt &input)
{
    return input.to_string();
}
inline HyperInt sqrt(const HyperInt &input)
{
    return input.square_root();
}
std::ostream &operator<<(std::ostream &output, const HyperInt &input)
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

HyperInt classic_factorial(UINT_64 end, UINT_64 start = 1) //累乘阶乘,可计算排列数A(n,m) n!/(n-m)!
{
    HyperInt result = 1;
    if (end < start)
    {
        return result;
    }
    for (UINT_64 i = start; i <= end; i++)
    {
        result *= i;
    }
    result.set_true_len();
    return result;
}
HyperInt factorial(UINT_64 end, UINT_64 start = 1, const UINT_32 rec_level = 0) //递归拆分,可计算排列数A(n,m) n!/(n-m)!
{
    if (end < start)
    {
        return HyperInt(1);
    }
    UINT_64 len = end - start;
    if (len < 120)
    {
        HyperInt result = 1;
        for (UINT_64 i = start; i <= end; i++)
        {
            result *= i;
        }
        return result;
    }
    UINT_64 mid = start + (len * 3 / 5);
#ifdef MULTITHREAD
    if (rec_level < hint::log2_threads)
    {
        std::future<HyperInt> first = std::async(factorial, mid, start, rec_level + 1);
        std::future<HyperInt> second = std::async(factorial, end, mid + 1, rec_level + 1);
        HyperInt result = first.get() * second.get();
        return result;
    }
    else
#endif
    {
        return factorial(mid, start, rec_level) * factorial(end, mid + 1, rec_level);
    }
}
HyperInt combination(UINT_64 n, UINT_64 m) // return C(n,m)组合数公式n!/((n-m)!m!)
{
    if (m > _HALF_(n))
    {
        return combination(n, n - m);
    }
    else if (n == 0 || m == 0)
    {
        return HyperInt(n);
    }
    else if (m >= n)
    {
        return HyperInt(1);
    }
    return factorial(n, (n - m + 1)) / factorial(m);
}
#endif