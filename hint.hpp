#include <iostream>
#include <future>
#include <thread>
#include <atomic>
#include <algorithm>
#include <random>
#include <stack>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <climits>

#ifndef HINT_HPP
#define HINT_HPP

//取消对宏MULTITHREAD的注释即可开启多线程
// #define MULTITHREAD

#if SIZE_MAX == 18446744073709551615ull
#define SIZE_T_BITS 64
#elif SIZE_MAX == 4294967295
#define SIZE_T_BITS 32
#else
#error "Unknown sys bits"
#endif

namespace hint
{
    using UINT_8 = uint8_t;
    using UINT_16 = uint16_t;
    using UINT_32 = uint32_t;
    using UINT_64 = uint64_t;
    using INT_32 = int32_t;
    using INT_64 = int64_t;
    using ULONG = unsigned long;
    using LONG = long;

    constexpr UINT_64 HINT_CHAR_BIT = 8;
    constexpr UINT_64 HINT_SHORT_BIT = 16;
    constexpr UINT_64 HINT_INT_BIT = 32;
    constexpr UINT_64 HINT_INT8_0XFF = UINT8_MAX;
    constexpr UINT_64 HINT_INT8_0X10 = (UINT8_MAX + 1ull);
    constexpr UINT_64 HINT_INT16_0XFF = UINT16_MAX;
    constexpr UINT_64 HINT_INT16_0X10 = (UINT16_MAX + 1ull);
    constexpr UINT_64 HINT_INT32_0XFF = UINT32_MAX;
    constexpr UINT_64 HINT_INT32_0X01 = 1;
    constexpr UINT_64 HINT_INT32_0X80 = 0X80000000ull;
    constexpr UINT_64 HINT_INT32_0X7F = INT32_MAX;
    constexpr UINT_64 HINT_INT32_0X10 = (UINT32_MAX + 1ull);
    constexpr UINT_64 HINT_INT64_0X80 = INT64_MIN;
    constexpr UINT_64 HINT_INT64_0X7F = INT64_MAX;

    constexpr size_t HINT_FFT_MIN = 128ull;
    constexpr size_t HINT_FFT_MAX = 4096ull;        // 4096为2分fft结果最大长度
    constexpr size_t HINT_QUAL_FFT_MAX = 262144ull; // 2^18为4分fft最大长度
    constexpr size_t HINT_NTT_MIN = 640ull;
    constexpr size_t HINT_NTT_MAX = 67108864ull;  // 2^26为2分ntt结果最大长度
    constexpr size_t HINT_NTT_MULTHLEN = 1024ull; //设置触发多线程的ntt长度
    constexpr double HINT_PI = 3.1415926535897932384626433832795;
    constexpr double HINT_2PI = 6.283185307179586476925286766559;

    struct h_int
    {
        hint::UINT_32 *array = nullptr;
        hint::INT_64 neg_n_len = 0;
        size_t size = 0;
    };
    struct Complex //对复数的定义
    {
        double real = 0, imaginary = 0;
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
        Complex(hint::UINT_64 n) // n等分圆周的复数,即x^n=1的解中除x=1以外辐角最小的一个解
        {
            real = std::cos(HINT_2PI / n);
            imaginary = std::sin(HINT_2PI / n);
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
            if (imaginary < 0)
                std::cout << imaginary << "i";
            else
                std::cout << "+" << imaginary << "i";
        }
    };
    const hint::UINT_32 hint_threads = std::thread::hardware_concurrency();
    const hint::UINT_32 log2_threads = std::ceil(std::log2(hint_threads));
    const Complex unit_omega_ary[] = {Complex(1), Complex(1 << 1), Complex(1 << 2), Complex(1 << 3), Complex(1 << 4), Complex(1 << 5), Complex(1 << 6), Complex(1 << 7),
                                      Complex(1 << 8), Complex(1 << 9), Complex(1 << 10), Complex(1 << 11), Complex(1 << 12), Complex(1 << 13), Complex(1 << 14),
                                      Complex(1 << 15), Complex(1 << 16), Complex(1 << 17), Complex(1 << 18), Complex(1 << 19), Complex(1 << 20)};

    template <typename T>
    constexpr bool is_neg(T x)
    {
        return x < 0;
    }
    template <typename T>
    constexpr bool is_odd(T x)
    {
        return static_cast<bool>(x & 1);
    }
    template <typename T>
    constexpr T twice(T x)
    {
        return x * 2;
    }
    template <typename T>
    constexpr T half(T x)
    {
        return x / 2;
    }
    template <typename T>
    constexpr void self_twice(T &x)
    {
        x *= 2;
    }
    template <typename T>
    constexpr void self_half(T &x)
    {
        x /= 2;
    }
    template <typename T>
    constexpr T qpow(T m, hint::UINT_64 n) //模板快速幂
    {
        T result = 1;
        while (n > 0)
        {
            if ((n & 1) != 0)
            {
                result = result * m;
            }
            m = m * m;
            n >>= 1;
        }
        return result;
    }
    constexpr hint::UINT_64 qpow(hint::UINT_64 m, hint::UINT_64 n, hint::UINT_64 mod) //取模快速幂
    {
        hint::UINT_64 result = 1;
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
    constexpr hint::UINT_64 gcd(hint::UINT_64 a, hint::UINT_64 b) //最大公因数
    {
        if (b == 0)
        {
            return a;
        }
        hint::UINT_64 tmp = b;
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
    hint::UINT_64 crt(hint::UINT_64 *mods, hint::UINT_64 *nums, size_t n) //中国剩余定理
    {
        hint::UINT_64 result = 0, mod_product = 1;
        for (size_t i = 0; i < n; i++)
        {
            mod_product *= mods[i];
        }
        for (size_t i = 0; i < n; i++)
        {
            hint::UINT_64 mod = mods[i];
            hint::UINT_64 tmp = mod_product / mod;
            hint::UINT_64 inv = qpow(tmp, mod - 2, mod);
            result += nums[i] * tmp * inv % mod_product;
        }
        return result % mod_product;
    }
    constexpr hint::UINT_64 qcrt(hint::UINT_64 num1, hint::UINT_64 num2,
                                 hint::UINT_64 mod1 = 167772161, hint::UINT_64 mod2 = 469762049,
                                 hint::UINT_64 inv1 = 104391568, hint::UINT_64 inv2 = 130489458) //快速计算两模数的中国剩余定理
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
        for (size_t rank = 1, gap; rank < fft_len; hint::self_twice(rank))
        {
            gap = hint::twice(rank);
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
    void ntt(hint::UINT_64 *input, size_t ntt_len, bool is_intt, const hint::UINT_64 mod = 998244353, hint::UINT_64 g_root = 3) //快速数论变换
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
                hint::UINT_64 tmp = input[i];
                input[i] = input[rev[i]];
                input[rev[i]] = tmp;
            }
        }
        delete[] rev;
        if (is_intt)
        {
            g_root = qpow(g_root, mod - 2, mod);
        }
        hint::UINT_64 unit_omega, omega;
        for (size_t rank = 1, gap; rank < ntt_len; hint::self_twice(rank))
        {
            gap = hint::twice(rank);
            unit_omega = qpow(g_root, (mod - 1) / gap, mod);
            for (size_t begin = 0; begin < ntt_len; begin += gap)
            {
                omega = 1;
                for (size_t pos = begin; pos < begin + rank; pos++)
                {
                    hint::UINT_64 tmp1 = input[pos];
                    hint::UINT_64 tmp2 = (input[pos + rank] % mod) * omega % mod;
                    input[pos] = (tmp1 + tmp2) % mod;
                    input[pos + rank] = (mod + tmp1 - tmp2) % mod;
                    omega = omega * unit_omega % mod;
                }
            }
        }
        if (is_intt)
        {
            hint::UINT_64 inv = qpow(ntt_len, mod - 2, mod);
            for (size_t i = 0; i < ntt_len; ++i)
            {
                input[i] = input[i] * inv % mod;
            }
        }
    }
    template <typename T1, typename T2, typename T3>
    void normal_convolution(T1 *const ary1, T2 *const ary2, T3 *const out, size_t len1, size_t len2)
    {
        ary_clr(out, len1 + len2);
        for (size_t i = 0; i < len1; i++)
        {
            for (size_t j = 0; j < len2; j++)
            {
                out[i + j] += static_cast<T2>(ary1[i]) * ary2[j];
            }
        }
    }
    void fft_convolution(Complex *const fft_ary1, Complex *const fft_ary2, Complex *const out, size_t fft_len) //快速傅里叶变换卷积分
    {
#ifdef MULTITHREAD
        bool multi_threads = hint::hint_threads >= 2 && fft_len >= 2 * hint::HINT_FFT_MAX;
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
    void ntt_convolution(hint::UINT_64 *const ntt_ary1, hint::UINT_64 *const ntt_ary2, hint::UINT_64 *const out, size_t ntt_len) //数论变换卷积分
    {
        constexpr hint::UINT_64 mod1 = 2013265921, mod2 = 2281701377;
        constexpr hint::UINT_64 root1 = 31, root2 = 3;
        hint::UINT_64 *ntt_ary3 = new hint::UINT_64[ntt_len];
        ary_copy(ntt_ary3, ntt_ary1, ntt_len);
        hint::UINT_64 *ntt_ary4 = ntt_ary3;
        std::function<void(hint::UINT_64 *, hint::UINT_64 *)> mul_func = [=](hint::UINT_64 *ary1, hint::UINT_64 *ary2)
        {
            for (size_t i = 0; i < ntt_len; i++)
            {
                ary1[i] = ary1[i] * ary2[i];
            } //每一位相乘
        };
#ifdef MULTITHREAD
        bool multi_threads = hint::hint_threads >= 2 && ntt_len > hint::HINT_NTT_MULTHLEN;
        if (multi_threads)
        {
            std::future<void> ntt_th1 = std::async(ntt, ntt_ary1, ntt_len, false, mod1, root1); // 多线程快速数论变换
            std::future<void> ntt_th2 = std::async(ntt, ntt_ary3, ntt_len, false, mod2, root2);
            if (ntt_ary1 != ntt_ary2)
            {
                ntt_ary4 = new hint::UINT_64[ntt_len];
                ary_copy(ntt_ary4, ntt_ary2, ntt_len);
                std::future<void> ntt_th3 = std::async(ntt, ntt_ary2, ntt_len, false, mod1, root1);
                std::future<void> ntt_th4 = std::async(ntt, ntt_ary4, ntt_len, false, mod2, root2);
                ntt_th3.wait();
                ntt_th4.wait();
            }
            ntt_th1.wait();
            ntt_th2.wait();

            mul_func(ntt_ary1, ntt_ary2);
            mul_func(ntt_ary3, ntt_ary4); //每一位相乘

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
                ntt_ary4 = new hint::UINT_64[ntt_len];
                ary_copy(ntt_ary4, ntt_ary2, ntt_len);
                ntt(ntt_ary2, ntt_len, false, mod1, root1); //快速数论变换
                ntt(ntt_ary4, ntt_len, false, mod2, root2);
            }

            mul_func(ntt_ary1, ntt_ary2);
            mul_func(ntt_ary3, ntt_ary4); //每一位相乘

            ntt(ntt_ary1, ntt_len, true, mod1, root1); //逆变换
            ntt(ntt_ary3, ntt_len, true, mod2, root2);
        }
        constexpr hint::UINT_64 inv1 = qpow(mod1, mod2 - 2, mod2);
        constexpr hint::UINT_64 inv2 = qpow(mod2, mod1 - 2, mod1);
        for (size_t i = 0; i < ntt_len; i++)
        {
            out[i] = qcrt(ntt_ary1[i], ntt_ary3[i], mod1, mod2, inv1, inv2);
        } //使用中国剩余定理变换
        delete[] ntt_ary3;
        if (ntt_ary4 != ntt_ary3)
        {
            delete[] ntt_ary4;
        }
    }
    template <typename T1, typename T2, typename T3>
    void trans_add(const T1 *in1, const T2 *in2, T3 *out, size_t len1, size_t len2, const hint::INT_64 base = 100) //可计算多项式的加法,默认为100进制
    {
        size_t result_len = std::max(len1, len2);
        hint::INT_64 tmp = 0;
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
            out[count] = static_cast<T3>(tmp % base);
            tmp /= base;
            count++;
        }
        if (tmp > 0)
        {
            out[count] = static_cast<T3>(tmp % base);
        }
    }
    template <typename T1, typename T2, typename T3>
    void trans_mul(const T1 *in1, const T2 *in2, T3 *out, size_t len1, size_t len2, const hint::INT_64 base = 100) //计算多项式的乘法
    {
        size_t out_len = len1 + len2;
        bool is_fft = false;
        if (out_len < 640)
        {
            hint::UINT_64 *con_result = new hint::UINT_64[out_len];
            ary_clr(con_result, out_len);
            normal_convolution(in1, in2, con_result, len1, len2);
            hint::UINT_64 tmp = 0;
            size_t pos = 0;
            while (pos < out_len)
            {
                tmp += con_result[pos];
                out[pos] = static_cast<T3>(tmp % base);
                tmp /= base;
                pos++;
            } //整理每一位
            delete[] con_result;
        }
        else if (out_len <= 2097152)
        {
            size_t fft_len = 1ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
            hint::Complex *fft_in1 = new hint::Complex[fft_len];
            hint::Complex *fft_in2 = new hint::Complex[fft_len];
            com_ary_copy(fft_in1, in1, len1);
            com_ary_copy(fft_in2, in2, len2);
            fft_convolution(fft_in1, fft_in2, fft_in1, fft_len);
            hint::UINT_64 tmp = 0;
            size_t pos = 0;
            while (pos < out_len)
            {
                tmp += static_cast<hint::UINT_64>(fft_in1[pos].real + 0.5);
                out[pos] = static_cast<T3>(tmp % base);
                tmp /= base;
                pos++;
            } //整理每一位
            delete[] fft_in1;
            delete[] fft_in2;
        }
        else if (out_len <= 1073741824)
        {
            size_t ntt_len = 1ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
            hint::UINT_64 *ntt_ary1 = new hint::UINT_64[ntt_len];
            hint::UINT_64 *ntt_ary2 = new hint::UINT_64[ntt_len];
            for (size_t i = 0; i < len1; i++)
            {
                ntt_ary1[i] = static_cast<hint::UINT_64>(in1[i]);
            }
            for (size_t i = 0; i < len2; i++)
            {
                ntt_ary2[i] = static_cast<hint::UINT_64>(in2[i]);
            }
            ary_clr(ntt_ary1 + len1, ntt_len - len1);
            ary_clr(ntt_ary2 + len2, ntt_len - len2);
            ntt_convolution(ntt_ary1, ntt_ary2, ntt_ary1, ntt_len);
            hint::UINT_64 tmp = 0;
            size_t pos = 0;
            while (pos < out_len)
            {
                tmp += static_cast<hint::UINT_64>(ntt_ary1[pos]);
                out[pos] = static_cast<T2>(tmp % base);
                tmp /= base;
                pos++;
            } //整理每一位
            delete[] ntt_ary1;
            delete[] ntt_ary2;
        }
        else
        {
            throw("Length error: too long");
        }
    }
    template <typename T1, typename T2>
    void trans_square(const T1 *in, T2 *out, size_t len, const hint::INT_64 base = 100) //平方
    {
        size_t out_len = hint::twice(len);
        if (out_len < 640)
        {
            hint::UINT_64 *con_result = new hint::UINT_64[out_len];
            ary_clr(con_result, out_len);
            normal_convolution(in, in, con_result, len, len);
            hint::UINT_64 tmp = 0;
            size_t pos = 0;
            while (pos < out_len)
            {
                tmp += con_result[pos];
                out[pos] = static_cast<T2>(tmp % base);
                tmp /= base;
                pos++;
            } //整理每一位
            delete[] con_result;
        }
        else if (out_len <= 2097152)
        {
            size_t fft_len = 1ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
            hint::Complex *fft_ary = new hint::Complex[fft_len];
            com_ary_copy(fft_ary, in, len);
            fft_convolution(fft_ary, fft_ary, fft_ary, fft_len);
            hint::UINT_64 tmp = 0;
            size_t pos = 0;
            while (pos < out_len)
            {
                tmp += static_cast<hint::UINT_64>(fft_ary[pos].real + 0.5);
                out[pos] = static_cast<T2>(tmp % base);
                tmp /= base;
                pos++;
            } //整理每一位
            delete[] fft_ary;
        }
        else if (out_len <= 1073741824)
        {
            size_t ntt_len = 1ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
            hint::UINT_64 *ntt_ary = new hint::UINT_64[ntt_len];
            for (size_t i = 0; i < len; i++)
            {
                ntt_ary[i] = static_cast<hint::UINT_64>(in[i]);
            }
            ary_clr(ntt_ary + len, ntt_len - len);
            ntt_convolution(ntt_ary, ntt_ary, ntt_ary, ntt_len);
            hint::UINT_64 tmp = 0;
            size_t pos = 0;
            while (pos < out_len)
            {
                tmp += static_cast<hint::UINT_64>(ntt_ary[pos]);
                out[pos] = static_cast<T2>(tmp % base);
                tmp /= base;
                pos++;
            } //整理每一位
            delete[] ntt_ary;
        }
        else
        {
            throw("Length error: too long");
        }
    }
    template <hint::UINT_32 BASE1, hint::UINT_32 BASE2, typename T, typename UNIT_T = hint::UINT_8>
    void base_conversion(T *data_ary, size_t &len) // 256进制转为100进制的数组方便打印和转字符串
    {
        if (len <= 1)
        {
            if (len == 0 || BASE1 == BASE2)
            {
                return;
            }
            hint::UINT_64 tmp = data_ary[0];
            data_ary[0] = tmp % BASE2;
            tmp /= BASE2;
            if (tmp > 0)
            {
                len++;
                data_ary[1] = tmp;
            }
            return;
        }
        size_t max_rank = 1ull << static_cast<hint::UINT_16>(ceil(log2(len)) - 1);         // unit_ary存储的base1的最高次幂
        const hint::UINT_64 base1to2_len = std::ceil(std::log2(BASE1) / std::log2(BASE2)); // base1到base2的数长度的比值
        size_t unit_ary_num = static_cast<hint::UINT_16>(log2(max_rank)) + 1;              // unit_ary存储的base1各个幂次的个数
        size_t result_len = static_cast<size_t>(base1to2_len * hint::twice(max_rank));     // 结果的长度
        hint::ary_clr(data_ary + len, result_len - len);                                   // 清零

        size_t unit_ary_len = (hint::twice(max_rank) - 1) * base1to2_len; // unit_ary的长度1+2+4+...max_rank
        UNIT_T *unit_ary = new UNIT_T[unit_ary_len];                      // 用一个数组存储(base1)^1,(base1)^2,(base1)^4...
        ary_clr(unit_ary, unit_ary_len);
        unit_ary[0] = BASE1 % BASE2;
        if (base1to2_len > 1)
        {
            hint::UINT_64 tmp = BASE1 / BASE2;
            size_t pos = 1;
            while (tmp > 0)
            {
                unit_ary[pos] = tmp % BASE2;
                tmp /= BASE2;
                pos++;
            }
            pos = len;
            while (pos > 0)
            {
                pos--;
                tmp = data_ary[pos];
                size_t trans_pos = pos * base1to2_len;
                for (size_t i = 0; i < base1to2_len; i++)
                {
                    data_ary[trans_pos + i] = tmp % BASE2;
                    tmp /= BASE2;
                }
            }
        }
        for (size_t i = 0, offset = 0; i < unit_ary_num - 1; i++)
        {
            size_t len = (1ull << i) * base1to2_len;
            hint::trans_square(unit_ary + offset, unit_ary + offset + len, len, BASE2);
            offset += len;
        }
        size_t pos = 0;
        hint::UINT_8 *tmp_product = new hint::UINT_8[base1to2_len * 2];
        for (size_t i = 0, offset = 0; i < unit_ary_num; i++)
        {
            size_t len = (1ull << i) * base1to2_len;
            size_t gap = hint::twice(len);
            pos = 0;
            tmp_product = hint::ary_realloc(tmp_product, gap);
            // hint::ary_clr(tmp_product, gap);
            while (pos < result_len)
            {
                hint::trans_mul(unit_ary + offset, data_ary + pos + len, tmp_product, len, len, BASE2);
                hint::trans_add(tmp_product, data_ary + pos, data_ary + pos, gap, len, BASE2);
                // hint::ary_copy(data_ary + pos, tmp_product, gap);
                pos += gap;
            }
            offset += len;
        }
        while (data_ary[result_len - 1] == 0)
        {
            result_len--;
        }
        len = result_len;
        delete[] tmp_product;
    }
    std::string ui64to_string(hint::UINT_64 input, hint::UINT_8 digits)
    {
        std::string result(digits, '0');
        for (hint::UINT_8 i = 0; i < digits; i++)
        {
            result[digits - i - 1] = static_cast<char>(input % 10 + '0');
            input /= 10;
        }
        return result;
    }
    hint::UINT_64 stoui64(const std::string &str, const hint::UINT_32 base = 10)
    {
        hint::UINT_64 result = 0;
        size_t len = str.size();
        for (size_t i = 0; i < len && i < 19; i++)
        {
            result *= 10;
            char c = str[i];
            if (c >= '0' && c <= '9')
            {
                result += static_cast<hint::UINT_64>(c - '0');
            }
        }
        return result;
    }
}
class HyperInt
{
private:
    hint::h_int data;
    constexpr void change_length(size_t new_length) //设置新的长度
    {
        if (new_length > data.size)
        {
            new_length = data.size;
        }
        if (new_length == 0)
        {
            neg_sign(false);
        }
        data.neg_n_len = (data.neg_n_len & hint::HINT_INT64_0X80) | new_length;
    }
    constexpr size_t generate_size(size_t new_size_input) const //生成1.5倍数组空间
    {
        if (new_size_input <= 2)
        {
            return 2;
        }
        size_t size1 = 1ull << static_cast<hint::UINT_16>(ceil(log2(new_size_input)));
        size_t size2 = hint::half(size1);
        size2 = size2 + hint::half(size2);
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
        neg_sign(false);
        hint::ary_clr(data.array, data.size);
    }
    void fill_element(hint::UINT_32 element) //填充相同的元素
    {
        std::fill(data.array, data.array + length(), element);
    }
    void quick_self_twice() //快速左移一位,只能操作大小为2^n的数
    {
        size_t len_tmp = length();
        if (len_tmp > 0)
        {
            hint::UINT_32 tmp = data.array[len_tmp - 1];
            hint::self_twice(tmp);
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
            hint::UINT_32 tmp = data.array[len_tmp - 1];
            hint::self_half(tmp);
            if (tmp == 0)
            {
                len_tmp--;
                if (len_tmp == 0)
                {
                    return;
                }
                change_length(len_tmp);
                data.array[len_tmp - 1] = hint::HINT_INT32_0X80;
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
    static HyperInt hint_mul(const HyperInt &input1, const HyperInt &input2)
    {
        return input1 * input2;
    }
    static HyperInt hint_square(const HyperInt &input)
    {
        return input.square();
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
        hint::UINT_64 tmp = 0, sum = 0;
        for (size_t pos1 = 0; pos1 < len1; pos1++)
        {
            for (size_t pos2 = 0; pos2 < len2; pos2++)
            {
                tmp = static_cast<hint::UINT_64>(data.array[pos1]) * input.data.array[pos2];
                for (size_t pos3 = pos1 + pos2; pos3 < result_len; pos3++)
                {
                    sum = tmp + result.data.array[pos3];
                    result.data.array[pos3] = static_cast<hint::UINT_32>(sum & hint::HINT_INT32_0XFF);
                    if ((sum >> hint::HINT_INT_BIT) > 0)
                    {
                        tmp = sum >> hint::HINT_INT_BIT;
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
        size_t result_len = hint::twice(len);
        result.reset_size(result_len);
        result.clear();
        result.change_length(result_len);
        hint::UINT_64 tmp = 0, sum = 0;
        for (size_t pos1 = 0; pos1 < len; pos1++)
        {
            for (size_t pos2 = 0; pos2 < len; pos2++)
            {
                tmp = static_cast<hint::UINT_64>(data.array[pos1]) * data.array[pos2];
                for (size_t pos3 = pos1 + pos2; pos3 < result_len; pos3++)
                {
                    sum = tmp + result.data.array[pos3];
                    result.data.array[pos3] = static_cast<hint::UINT_32>(sum & hint::HINT_INT32_0XFF);
                    if ((sum >> hint::HINT_INT_BIT) > 0)
                    {
                        tmp = sum >> hint::HINT_INT_BIT;
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
        size_t fft_len = 2ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
        bool is_qual_div = out_len > hint::HINT_FFT_MAX;
        if (is_qual_div)
        {
            hint::self_twice(fft_len);
        }
        hint::Complex *fft_ary1 = new hint::Complex[fft_len * 2];
        hint::Complex *fft_ary2 = fft_ary1 + fft_len; // new hint::Complex[fft_len];

        //每一位分解为短整数后存入复数组
        if (is_qual_div)
        {
            hint::UINT_8 *ary1_8 = reinterpret_cast<hint::UINT_8 *>(data.array);
            hint::UINT_8 *ary2_8 = reinterpret_cast<hint::UINT_8 *>(input.data.array);
            size_t data_len1 = len1 * 4;
            size_t data_len2 = len2 * 4;

            hint::com_ary_copy(fft_ary1, ary1_8, data_len1);
            hint::com_ary_copy(fft_ary2, ary2_8, data_len2);
        }
        else
        {
            hint::UINT_16 *ary1_16 = reinterpret_cast<hint::UINT_16 *>(data.array);
            hint::UINT_16 *ary2_16 = reinterpret_cast<hint::UINT_16 *>(input.data.array);
            size_t data_len1 = len1 * 2;
            size_t data_len2 = len2 * 2;

            hint::com_ary_copy(fft_ary1, ary1_16, data_len1);
            hint::com_ary_copy(fft_ary2, ary2_16, data_len2);
        }

        hint::fft_convolution(fft_ary1, fft_ary2, fft_ary1, fft_len);

        hint::UINT_64 tmp = 0;
        if (is_qual_div)
        {
            hint::UINT_8 *ary_8 = reinterpret_cast<hint::UINT_8 *>(result.data.array);
            size_t data_len = out_len * 4;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<hint::UINT_64>(fft_ary1[i].real + 0.5);
                ary_8[i] = static_cast<hint::UINT_8>(tmp & hint::HINT_INT8_0XFF);
                tmp >>= hint::HINT_CHAR_BIT;
            }
        }
        else
        {
            hint::UINT_16 *ary_16 = reinterpret_cast<hint::UINT_16 *>(result.data.array);
            size_t data_len = out_len * 2;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<hint::UINT_64>(fft_ary1[i].real + 0.5);
                ary_16[i] = static_cast<hint::UINT_16>(tmp & hint::HINT_INT16_0XFF);
                tmp >>= hint::HINT_SHORT_BIT;
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
        size_t out_len = hint::twice(len);

        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);
        size_t fft_len = 2ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));

        bool is_qual_div = out_len > hint::HINT_FFT_MAX;
        if (is_qual_div)
        {
            hint::self_twice(fft_len);
        }
        hint::Complex *fft_ary = new hint::Complex[fft_len];
        //每一位分解为短整数后存入复数组
        if (is_qual_div)
        {
            hint::UINT_8 *ary_8 = reinterpret_cast<hint::UINT_8 *>(data.array);
            size_t data_len = len * 4;

            hint::com_ary_copy(fft_ary, ary_8, data_len);
        }
        else
        {
            hint::UINT_16 *ary_16 = reinterpret_cast<hint::UINT_16 *>(data.array);
            size_t data_len = len * 2;

            hint::com_ary_copy(fft_ary, ary_16, data_len);
        }

        hint::fft_convolution(fft_ary, fft_ary, fft_ary, fft_len); //卷积

        hint::UINT_64 tmp = 0;
        if (is_qual_div)
        {
            hint::UINT_8 *ary_8 = reinterpret_cast<hint::UINT_8 *>(result.data.array);
            size_t data_len = out_len * 4;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<hint::UINT_64>(fft_ary[i].real + 0.5);
                ary_8[i] = static_cast<hint::UINT_8>(tmp & hint::HINT_INT8_0XFF);
                tmp >>= hint::HINT_CHAR_BIT;
            }
        }
        else
        {
            hint::UINT_16 *ary_16 = reinterpret_cast<hint::UINT_16 *>(result.data.array);
            size_t data_len = out_len * 2;
            for (size_t i = 0; i < data_len; i++)
            {
                tmp += static_cast<hint::UINT_64>(fft_ary[i].real + 0.5);
                ary_16[i] = static_cast<hint::UINT_16>(tmp & hint::HINT_INT16_0XFF);
                tmp >>= hint::HINT_SHORT_BIT;
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

        size_t ntt_len = 2ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
        hint::UINT_64 *ntt_ary1 = new hint::UINT_64[ntt_len * 2];
        hint::UINT_64 *ntt_ary2 = ntt_ary1 + ntt_len; // new hint::UINT_64[ntt_len];

        hint::ary_clr(ntt_ary1, ntt_len * 2);
        // hint::ary_clr(ntt_ary2, ntt_len);

        hint::UINT_16 *ary1_16 = reinterpret_cast<hint::UINT_16 *>(data.array); //每一位分解为短整数后存入数组
        hint::UINT_16 *ary2_16 = reinterpret_cast<hint::UINT_16 *>(input.data.array);
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

        hint::UINT_64 tmp = 0;
        hint::UINT_16 *ary_16 = reinterpret_cast<hint::UINT_16 *>(result.data.array);
        size_t data_len = out_len * 2;
        for (size_t i = 0; i < data_len; i++)
        {
            tmp += ntt_ary1[i];
            ary_16[i] = static_cast<hint::UINT_16>(tmp & hint::HINT_INT16_0XFF);
            tmp >>= hint::HINT_SHORT_BIT;
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
        size_t out_len = hint::twice(len);
        result.reset_size(out_len);
        result.clear();
        result.change_length(out_len);

        size_t ntt_len = 2ull << static_cast<hint::UINT_16>(ceil(log2(out_len)));
        hint::UINT_64 *ntt_ary = new hint::UINT_64[ntt_len];
        hint::ary_clr(ntt_ary, ntt_len);

        //每一位分解为短整数后存入数组
        hint::UINT_16 *ary_16 = reinterpret_cast<hint::UINT_16 *>(data.array);
        size_t data_len = len * 2;
        for (size_t i = 0; i < data_len; i++)
        {
            ntt_ary[i] = ary_16[i];
        }

        hint::ntt_convolution(ntt_ary, ntt_ary, ntt_ary, ntt_len);

        hint::UINT_64 tmp = 0;
        ary_16 = reinterpret_cast<hint::UINT_16 *>(result.data.array);
        data_len = out_len * 2;
        for (size_t i = 0; i < data_len; i++)
        {
            tmp += ntt_ary[i];
            ary_16[i] = static_cast<hint::UINT_16>(tmp & hint::HINT_INT16_0XFF);
            tmp >>= hint::HINT_SHORT_BIT;
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
            sub_len1 = hint::half(len1);
            sub_len_public = len1 - sub_len1;
            sub_len2 = 0;
            if (len2 > sub_len_public)
            {
                sub_len2 = len2 - sub_len_public;
            }
        }
        else
        {
            sub_len2 = hint::half(len2);
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
            std::future<HyperInt> sub_e_th = std::async(hint_mul, sub_a, sub_c);
            std::future<HyperInt> sub_f_th = std::async(hint_mul, sub_b, sub_d);
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
        hint::INT_64 tmp = 0;
        while (count < len_pub)
        {
            if (pos_1 < len_e)
            {
                tmp += sub_e.data.array[pos_1];
                pos_1++;
            }
            result.data.array[count] = static_cast<hint::UINT_32>(tmp + hint::HINT_INT32_0X10);
            tmp >>= hint::HINT_INT_BIT;
            count++;
        }
        hint::self_twice(len_pub);
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
            result.data.array[count] = static_cast<hint::UINT_32>(tmp + hint::HINT_INT32_0X10);
            tmp >>= hint::HINT_INT_BIT;
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
            result.data.array[count] = static_cast<hint::UINT_32>(tmp + hint::HINT_INT32_0X10);
            tmp >>= hint::HINT_INT_BIT;
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
        size_t result_len = hint::twice(len);
        if (result_len < 4)
        {
            return normal_square();
        }
        result.reset_size(result_len);
        result.clear();
        result.change_length(result_len);
        size_t sub_len, sub_len_public;
        sub_len = hint::half(len);
        sub_len_public = len - sub_len;

        HyperInt sub_a = split(0, sub_len_public);
        HyperInt sub_b = split(sub_len_public, sub_len);

        HyperInt sub_ab;
#ifdef MULTITHREAD
        if (hint::hint_threads > 1)
        {
            std::future<HyperInt> sub_a_th = std::async(hint_square, sub_a);
            std::future<HyperInt> sub_b_th = std::async(hint_square, sub_b);
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
        hint::INT_64 tmp = 0;
        while (count < len_pub)
        {
            if (pos_1 < len_aa)
            {
                tmp += sub_a.data.array[pos_1];
                pos_1++;
            }
            result.data.array[count] = static_cast<hint::UINT_32>(tmp & hint::HINT_INT32_0XFF);
            tmp >>= hint::HINT_INT_BIT;
            count++;
        }
        hint::self_twice(len_pub);
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
            result.data.array[count] = static_cast<hint::UINT_32>(tmp & hint::HINT_INT32_0XFF);
            tmp >>= hint::HINT_INT_BIT;
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
            result.data.array[count] = static_cast<hint::UINT_32>(tmp & hint::HINT_INT32_0XFF);
            tmp >>= hint::HINT_INT_BIT;
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
        const size_t precise_len = hint::twice(len1);

        HyperInt input_inv = HyperInt(1); //除数的倒数
        size_t input_inv_len = len2;      //小数部分长度
        constexpr size_t times = 640;
        for (size_t i = 0; i < times; i++)
        {
            HyperInt input_inv_tmp = input_inv;
            HyperInt input_inv_sq = input_inv.square(); //除数倒数的平方,小数部分为两倍input_inv_len
            input_inv_sq = input * input_inv_sq;
            size_t input_inv_sq_len = input_inv_len * 2;
            input_inv = input_inv.l_shift((input_inv_sq_len - input_inv_len) * hint::HINT_INT_BIT + 1);
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
            result.div_mod(input.first_int32());
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
        hint::UINT_64 tmp = 0, first_num2 = input.first_int64();
        hint::UINT_64 try_num = 0;
        while (!dividend.abs_smaller(input))
        {
            shift = dividend.length() - len2;
            hint::UINT_64 first_num1 = dividend.first_int64();
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
                double digit_3 = static_cast<double>(first_num1) * static_cast<double>(hint::HINT_INT32_0X10);
                tmp = static_cast<hint::UINT_64>(digit_3 / first_num2);
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
            result.data.array[shift] = static_cast<hint::UINT_32>(tmp & hint::HINT_INT32_0XFF);
        }
        result.neg_sign(static_cast<bool>(is_neg() ^ input.is_neg()));
        result.set_true_len();
        return result;
    }
    HyperInt newton_sqrt() const
    {
        size_t len = length();
        HyperInt shift_in = l_shift(128);
        HyperInt result((len + 5) / 2, 0XFFFFFFFF);
        HyperInt left, right;
        constexpr size_t times = 320;
        for (size_t i = 0; i < times; i++)
        {
            HyperInt tmp = shift_in / result + result;
            tmp.self_half();
            if (tmp.abs_compare(result) >= 0)
            {
                left = result.r_shift(64);
                right = tmp.r_shift(64);
                break;
            }
            result = std::move(tmp);
        }
        while (left.abs_smaller(right))
        {
            HyperInt mid = (left + right).half();
            if (left.abs_equal(mid))
            {
                return left;
            }
            hint::INT_32 cmp = abs_compare(mid.square());
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
    HyperInt normal_sqrt() const
    {
        size_t len = length();
        HyperInt left, right;
        left.reset_size(hint::half(len) + 1);
        left.clear();
        left.change_length(hint::half(len) + 1);
        left.data.array[hint::half(len) - 1] = 1;
        left.set_true_len();
        right.reset_size(hint::half(len) + 1);
        right.clear();
        right.change_length(hint::half(len) + 1);
        right.data.array[hint::half(len) - 1] = 2;
        right.set_true_len();
        while (!abs_smaller(right.square()))
        {
            right.quick_self_twice();
            left.quick_self_twice();
        }
        while (left.abs_smaller(right))
        {
            HyperInt mid = (left + right).half();
            if (left.abs_equal(mid))
            {
                return left;
            }
            hint::INT_32 cmp = abs_compare(mid.square());
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

public:
    HyperInt norsqrt()
    {
        return normal_sqrt();
    }
    HyperInt newsqrt()
    {
        return newton_sqrt();
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
    HyperInt(size_t new_length, hint::UINT_32 num) //填充new_length个num进行构造
    {
        data.size = generate_size(new_length);
        change_length(new_length);
        neg_sign(false);
        if (data.array != nullptr)
        {
            delete[] data.array;
        }
        data.array = new hint::UINT_32[data.size];
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
    HyperInt(const std::string &input) // string 参数构造
    {
        string_in(input);
    }
    HyperInt(char input[]) //字符串构造
    {
        string_in(input);
    }
    HyperInt(const char input[]) //字符串构造
    {
        string_in(input);
    }
    template <typename T>
    HyperInt(T input) //  通用构造
    {
        bool neg = hint::is_neg(input);
        hint::UINT_64 tmp = 0;
        if (neg)
        {
            tmp = static_cast<hint::UINT_64>(std::abs(static_cast<hint::INT_64>(input)));
        }
        else
        {
            tmp = static_cast<hint::UINT_64>(input);
        }
        data.size = 2;
        reset_size(2);
        change_length(2);
        data.array[0] = tmp & hint::HINT_INT32_0XFF;
        data.array[1] = tmp >> hint::HINT_INT_BIT;
        neg_sign(neg);
        set_true_len();
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
    HyperInt &operator=(char input[]) //字符串赋值
    {
        string_in(input);
        return *this;
    }
    template <typename T>
    HyperInt &operator=(T input) // hint::UINT_64 赋值
    {
        bool neg = hint::is_neg(input);
        hint::UINT_64 tmp = 0;
        if (neg)
        {
            tmp = static_cast<hint::UINT_64>(std::abs(static_cast<hint::INT_64>(input)));
        }
        else
        {
            tmp = static_cast<hint::UINT_64>(input);
        }
        data.size = 2;
        reset_size(2);
        change_length(2);
        data.array[0] = tmp & hint::HINT_INT32_0XFF;
        data.array[1] = tmp >> hint::HINT_INT_BIT;
        neg_sign(neg);
        set_true_len();
        return *this;
    }

    //基本操作
    void set_true_len();                            //去除前导0
    void neg_sign(bool neg);                        //设置符号是否为负
    hint::INT_64 div_mod(hint::UINT_32 divisor);    //自身除以divisor的同时返回余数
    hint::INT_64 mod(hint::UINT_32 divisor) const;   //返回对divisor的余数
    HyperInt power(hint::UINT_64 n) const;          //快速幂
    HyperInt square() const;                        //求自身的平方
    HyperInt square_root() const;                   //求自身的平方根
    bool is_neg() const;                            //返回符号是否为为负号
    static bool is_neg(const HyperInt &input);      //返回符号是否为为负号
    size_t length() const;                          //返回长度
    size_t size() const;                            //返回分配的数组空间
    HyperInt split(size_t begin, size_t len) const; //返回从下标begi开始长度为len的子数组
    hint::INT_64 to_int64() const;                  //转hint::INT_64
    hint::UINT_64 to_uint64() const;                //转hint::UINT_64
    std::string to_string() const;                  //转string,用10进制表示的字符串
    void string_in(const std::string &str);         //由十进制string字符串输入
    void string_in(const char str[]);               //由十进制C风格字符串输入
    void normal_string_in(const std::string &str);  //输入十进制字符串，普通算法
    void quick_string_in(const std::string &str);   //输入十进制字符串，快速算法
    void console_in();                              //从控制台读入十进制值
    void print_dec() const;                         //向控制台打印十进制值
    void print_hex() const;                         //向控制台打印十六进制值

    HyperInt add_sub(const HyperInt &input, bool is_add) const;                       //基础加减法a=b.add_sub(c,ture)->a=b+c;a=b.add_sub(c,fasle)->a=b-c,(b>c);
    void add_sub_inplace(const HyperInt &input, bool is_add, const size_t shift = 0); //就地加减 a+=b;a-=b,a加/减去左移位后的b，默认不移位
    void sub_inplace(const HyperInt &input);                                          //由减数调用,a.sub_inplace(b)->a=b-a;

    HyperInt &self_half();                  //自身右移一位
    HyperInt &self_twice();                 //自身左移一位
    HyperInt half() const;                  //返回右移一位后的值
    HyperInt twice() const;                 //返回左移一位后的值
    HyperInt r_shift(size_t n) const;       //右移n位
    HyperInt l_shift(size_t n) const;       //左移n位
    void reset_size(size_t new_size_input); //重新设定长度不小于new_size,1.5倍长度算法,在change_len()之前调用
    HyperInt abs() const;                   //返回绝对值
    hint::UINT_32 first_int32() const;      //返回开头一个元素转32位整数的结果
    hint::UINT_64 first_int64() const;      //返回开头两个元素转32位整数的结果

    hint::INT_32 abs_compare(const HyperInt &input, hint::INT_64 shift = 0) const; //自身和input移位shift比较，大于返回1，小于返回-1，等于返回0
    bool abs_larger(const HyperInt &input) const;                                  //绝对值是否大于input
    bool abs_smaller(const HyperInt &input) const;                                 //绝对值是否小于input
    bool abs_equal(const HyperInt &input) const;                                   //绝对值是否等于input
    bool equal_to_z() const;                                                       //判定是否为零
    bool is_even() const;                                                          //判断是否为偶数
    bool is_odd() const;                                                           //判断是否为奇数

    //逻辑运算
    bool operator==(const HyperInt &input) const;
    template <typename T>
    bool operator==(const T &input) const;

    bool operator!=(const HyperInt &input) const;
    template <typename T>
    bool operator!=(const T &input) const;

    bool operator>(const HyperInt &input) const;
    template <typename T>
    bool operator>(const T &input) const;

    bool operator>=(const HyperInt &input) const;
    template <typename T>
    bool operator>=(const T &input) const;

    bool operator<(const HyperInt &input) const;
    template <typename T>
    bool operator<(const T &input) const;

    bool operator<=(const HyperInt &input) const;
    template <typename T>
    bool operator<=(const T &input) const;

    //友元函数
    friend HyperInt abs(const HyperInt &input);
    friend void print(const HyperInt &input);

    friend bool operator==(const hint::INT_64 &input1, const HyperInt &input2);
    friend bool operator!=(const hint::INT_64 &input1, const HyperInt &input2);
    friend bool operator>(const hint::INT_64 &input1, const HyperInt &input2);
    friend bool operator>=(const hint::INT_64 &input1, const HyperInt &input2);
    friend bool operator<(const hint::INT_64 &input1, const HyperInt &input2);
    friend bool operator<=(const hint::INT_64 &input1, const HyperInt &input2);

    friend HyperInt operator+(const hint::INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator-(const hint::INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator*(const hint::INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator/(const hint::INT_64 &input1, const HyperInt &input2);
    friend HyperInt operator%(const hint::INT_64 &input1, const HyperInt &input2);

    friend hint::INT_64 &operator+=(hint::INT_64 &input1, const HyperInt &input2);
    friend hint::INT_64 &operator-=(hint::INT_64 &input1, const HyperInt &input2);
    friend hint::INT_64 &operator*=(hint::INT_64 &input1, const HyperInt &input2);
    friend hint::INT_64 &operator/=(hint::INT_64 &input1, const HyperInt &input2);
    friend hint::INT_64 &operator%=(hint::INT_64 &input1, const HyperInt &input2);

    friend std::string to_string(const HyperInt &input);
    friend std::ostream &operator<<(std::ostream &output, const HyperInt &input);
    friend std::istream &operator>>(std::istream &input, HyperInt &output);

    friend HyperInt randHyperInt(size_t len); //生成长度为len的随机数
    //算术运算

    HyperInt operator+(const HyperInt &input) const;
    template <typename T>
    HyperInt operator+(T input) const;
    HyperInt operator+() const;

    HyperInt operator-(const HyperInt &input) const;
    template <typename T>
    HyperInt operator-(T input) const;
    HyperInt operator-() const;

    HyperInt operator*(const HyperInt &input) const;
    template <typename T>
    HyperInt operator*(T input) const;

    HyperInt operator/(const HyperInt &input) const;
    template <typename T>
    HyperInt operator/(T input) const;

    HyperInt operator%(const HyperInt &input) const;
    template <typename T>
    HyperInt operator%(T input) const;

    HyperInt &operator+=(const HyperInt &input);
    template <typename T>
    HyperInt &operator+=(T input);

    HyperInt &operator-=(const HyperInt &input);
    template <typename T>
    HyperInt &operator-=(T input);

    HyperInt &operator*=(const HyperInt &input);
    template <typename T>
    HyperInt &operator*=(T input);

    HyperInt &operator/=(const HyperInt &input);
    template <typename T>
    HyperInt &operator/=(T input);

    HyperInt &operator%=(const HyperInt &input);
    template <typename T>
    HyperInt &operator%=(T input);

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

inline void HyperInt::set_true_len() //去除前导0
{
    size_t t_len = length();
    while (t_len > 0 && data.array[t_len - 1] == 0)
    {
        t_len--;
    }
    change_length(t_len);
}
inline void HyperInt::neg_sign(bool neg) //设置符号是否为负
{
    if (neg)
    {
        data.neg_n_len = data.neg_n_len | hint::HINT_INT64_0X80;
    }
    else
    {
        data.neg_n_len = data.neg_n_len & hint::HINT_INT64_0X7F;
    }
}
hint::INT_64 HyperInt::div_mod(hint::UINT_32 divisor) //自身除以divisor的同时返回余数
{
    lldiv_t div_tmp;
    if (divisor == 0)
    {
        return to_int64();
    }
    hint::UINT_64 last_rem = 0, tmp = 0, rem_num = 0;
    bool result_neg = is_neg() ^ hint::is_neg(divisor);
    neg_sign(result_neg);
    size_t pos = length();
    while (pos > 1)
    {
        pos--;
        tmp = (last_rem << hint::HINT_INT_BIT) + data.array[pos];
        div_tmp = lldiv(tmp, divisor);                              //一次性得到商和余数
        data.array[pos] = static_cast<hint::UINT_32>(div_tmp.quot); //当前数为变商
        last_rem = div_tmp.rem;                                     //得到余数
    }
    tmp = (last_rem << hint::HINT_INT_BIT) + data.array[0];
    div_tmp = lldiv(tmp, divisor);
    data.array[0] = static_cast<hint::UINT_32>(div_tmp.quot);
    rem_num = div_tmp.rem;
    set_true_len();
    return rem_num;
}
hint::INT_64 HyperInt::mod(hint::UINT_32 divisor) const //返回对divisor的余数
{
    //  lldiv_t div_tmp;
    if (divisor == 0)
    {
        return to_int64();
    }
    hint::INT_64 last_rem = 0, tmp = 0, rem_num = 0;
    size_t pos = length();
    while (pos > 1)
    {
        pos--;
        tmp = (last_rem << hint::HINT_INT_BIT) + data.array[pos];
        auto div_tmp = lldiv(tmp, divisor); //一次性得到商和余数
        last_rem = div_tmp.rem;             //得到余数
    }
    tmp = (last_rem << hint::HINT_INT_BIT) + data.array[0];
    auto div_tmp = lldiv(tmp, divisor);
    rem_num = div_tmp.rem;
    if (is_neg())
    {
        rem_num = -rem_num;
    }
    return rem_num;
}
inline HyperInt HyperInt::power(hint::UINT_64 n) const //快速幂
{
    HyperInt tmp(*this), result = HyperInt(1);
    if (!hint::is_odd(n))
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
        hint::self_half(n);
    }
    return result;
}
inline HyperInt HyperInt::square() const //求自身的平方
{
    size_t len = length();
    if (len <= 48)
    {
        return normal_square();
    }
    else if (len <= hint::HINT_QUAL_FFT_MAX / 2)
    {
        return fft_square();
    }
    else if (len < hint::HINT_NTT_MAX / 16)
    {
        return karatsuba_square();
    }
    else if (len <= hint::HINT_NTT_MAX / 2)
    {
        return ntt_square();
    }
    else
    {
        return karatsuba_square();
    }
}
inline HyperInt HyperInt::square_root() const
{
    size_t len = length();
    if (len < 2)
    {
        return HyperInt(static_cast<hint::UINT_32>(sqrt(first_int32())));
    }
    else if (len <= 4)
    {
        return normal_sqrt();
    }
    else
    {
        return newton_sqrt();
    }
}
inline bool HyperInt::is_neg() const //返回符号是否为为负号
{
    return data.neg_n_len & hint::HINT_INT64_0X80;
}
inline bool HyperInt::is_neg(const HyperInt &input)
{
    return input.is_neg();
}
inline size_t HyperInt::length() const //返回长度
{
    return data.neg_n_len & hint::HINT_INT64_0X7F;
}
inline size_t HyperInt::size() const //返回分配的数组空间
{
    return data.size;
}
inline HyperInt HyperInt::split(size_t begin, size_t len) const //返回从下标begi开始长度为len的子数组
{
    if (len == 0)
    {
        return HyperInt(0);
    }
    size_t data_len = length();
    if (begin >= data_len)
    {
        hint::INT_64 num = data.array[data_len - 1];
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
inline hint::INT_64 HyperInt::to_int64() const //转hint::INT_64
{
    if (to_uint64() == hint::HINT_INT64_0X80)
    {
        return static_cast<hint::INT_64>(hint::HINT_INT64_0X80);
    }
    hint::INT_64 out = 0;
    out = (data.array[1] & hint::HINT_INT32_0X7F);
    out <<= hint::HINT_INT_BIT;
    out += data.array[0];
    if (is_neg())
    {
        out = -out;
    }
    return out;
}
inline hint::UINT_64 HyperInt::to_uint64() const //转hint::UINT_64
{
    hint::UINT_64 out = 0;
    out = data.array[1];
    out <<= hint::HINT_INT_BIT;
    out += data.array[0];
    return out;
}
std::string HyperInt::to_string() const //转string,用10进制表示的字符串
{
    if (equal_to_z())
    {
        return std::string("0");
    }
    size_t in_len = length();
    if (in_len <= 2)
    {
        if (is_neg())
        {
            return std::string("-") + std::to_string(first_int64());
        }
        return std::to_string(first_int64());
    }
    std::string result_str;
    if (is_neg())
    {
        result_str += '-';
    }
    if (in_len <= 1536)
    {
        HyperInt input(*this);
        hint::UINT_64 rem_n = 0;
        constexpr hint::UINT_32 factor = 1e9;
        std::stack<std::string> str_stack;
        while (input.abs_larger(HyperInt(0)))
        {
            rem_n = input.div_mod(factor);
            str_stack.emplace(hint::ui64to_string(rem_n, 9));
        }
        if (!str_stack.empty())
        {
            str_stack.top() = std::to_string(rem_n);
        }
        while (!str_stack.empty())
        {
            result_str += str_stack.top();
            str_stack.pop();
        }
    }
    else
    {
        size_t result_len = in_len * 4;
        size_t ary_len = 2ull << static_cast<hint::UINT_16>(std::ceil(std::log2(result_len)));
        hint::UINT_8 *result = new hint::UINT_8[ary_len]; //结果数组
        memcpy(result, data.array, in_len * sizeof(*data.array));
        while (result[result_len - 1] == 0)
        {
            result_len--;
        }

        hint::base_conversion<hint::HINT_INT8_0X10, 100>(result, result_len); //进制转换

        result_str += std::to_string(static_cast<hint::UINT_32>(result[result_len - 1]));
        result_len--;
        while (result_len > 0) //整理每一位将100进制转为十进制
        {
            result_len--;
            result_str += hint::ui64to_string(static_cast<hint::UINT_64>(result[result_len]), 2);
        }
        delete[] result;
    }
    return result_str;
}
inline void HyperInt::string_in(const std::string &str)
{
    size_t len = str.size();
    if (len > 200000 && len <= 268435456)
    {
        quick_string_in(str);
    }
    else
    {
        normal_string_in(str);
    }
}
inline void HyperInt::string_in(const char str[])
{
    string_in(std::string(str));
}
inline void HyperInt::normal_string_in(const std::string &str) //输入十进制字符串
{
    clear();
    constexpr uint64_t factor = 1e19;
    size_t str_len = str.size(), pos = 0;
    if (str_len == 0)
    {
        clear();
        set_true_len();
        return;
    }
    if (str[0] == '-')
    {
        pos++;
    }
    while (str_len - pos >= 19)
    {
        *this *= factor;
        *this += hint::stoui64(std::string(str.begin() + pos, str.begin() + pos + 19));
        pos += 19;
    }
    if (str_len - pos > 0)
    {
        *this *= hint::qpow(10ull, str_len - pos);
        *this += hint::stoui64(std::string(str.begin() + pos, str.end()));
    }
    if (str[0] == '-')
    {
        neg_sign(true);
    }
}
void HyperInt::quick_string_in(const std::string &str)
{
    size_t in_len = str.size();
    if (in_len == 0)
    {
        clear();
        set_true_len();
        return;
    }
    size_t trans_len = in_len / 2 + 1;
    size_t ary_len = 1ull << static_cast<hint::UINT_16>(ceil(log2(trans_len)));
    hint::UINT_8 *trans_ary = new hint::UINT_8[ary_len];
    hint::ary_clr(trans_ary, ary_len);
    bool neg = (str[0] == '-');
    for (size_t pos = 0; pos < trans_len; pos++)
    {
        hint::UINT_16 tmp = 0, rank = 1;
        for (size_t i = 0; i < 2 && in_len > 0; i++)
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
    while (trans_ary[trans_len - 1] == 0)
    {
        trans_len--;
    }
    hint::base_conversion<100, hint::HINT_INT8_0X10>(trans_ary, trans_len);
    reset_size(1 + trans_len / 4);
    clear();
    change_length(1 + trans_len / 4);
    memcpy(data.array, trans_ary, trans_len * sizeof(*trans_ary));
    set_true_len();
    if (!equal_to_z())
    {
        neg_sign(neg);
    }
    delete[] trans_ary;
}
inline void HyperInt::console_in() //从控制台读入十进制值
{
    clear();
    char tmp = '0';
    bool head = true, neg = false;
    while (tmp != '\n' && '0' <= tmp && tmp <= '9')
    {
        if (tmp != '\n')
        {
            *this *= static_cast<hint::UINT_64>(10);
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
void HyperInt::print_dec() const //向控制台打印十进制值
{
    if (equal_to_z())
    {
        printf("0");
    }
    if (is_neg())
    {
        printf("-");
    }
    size_t in_len = length();
    if (in_len <= 2)
    {
        std::cout << first_int64();
    }
    else if (in_len <= 1536)
    {
        HyperInt input(*this);
        hint::UINT_64 rem_n = 0;
        constexpr hint::UINT_32 factor = 1e9;
        std::stack<hint::UINT_32> num_stack;
        while (input.abs_larger(HyperInt(0)))
        {
            rem_n = input.div_mod(factor);
            num_stack.push(rem_n);
        }
        if (!num_stack.empty())
        {
            printf("%u", num_stack.top());
            num_stack.pop();
        }
        while (!num_stack.empty())
        {
            printf("%09u", num_stack.top());
            num_stack.pop();
        }
    }
    else
    {
        size_t result_len = in_len * 4;
        size_t ary_len = 2ull << static_cast<hint::UINT_16>(std::ceil(std::log2(result_len)));
        hint::UINT_8 *result = new hint::UINT_8[ary_len]; //结果数组
        memcpy(result, data.array, in_len * sizeof(*data.array));
        while (result[result_len - 1] == 0)
        {
            result_len--;
        }

        hint::base_conversion<hint::HINT_INT8_0X10, 100>(result, result_len);
        if (result_len > 0)
        {
            result_len--;
            printf("%d", result[result_len]);
        }
        while (result_len > 0)
        {
            result_len--;
            printf("%02d", result[result_len]);
        }

        delete[] result;
    }
    printf("\n");
}
inline void HyperInt::print_hex() const //向控制台打印十六进制值
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
    printf("\n");
}

inline HyperInt HyperInt::add_sub(const HyperInt &input, bool is_add) const //基础加减法a=b.add_sub(c,ture)->a=b+c;a=b.add_sub(c,fasle)->a=b-c,(b>c);
{
    HyperInt result;
    size_t len1 = length(), len2 = input.length();
    if (is_add)
    {
        size_t result_len = std::max(len1, len2) + 1;
        result.reset_size(result_len);
        result.change_length(result_len);
        hint::INT_64 tmp = 0;
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
            result.data.array[count] = static_cast<hint::UINT_32>(tmp & hint::HINT_INT32_0XFF);
            tmp >>= hint::HINT_INT_BIT;
            count++;
        }
    }
    else
    {
        size_t result_len = std::max(len1, len2);
        result.reset_size(result_len);
        result.change_length(result_len);
        hint::INT_64 tmp = 0;
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
            result.data.array[count] = static_cast<hint::UINT_32>(tmp + hint::HINT_INT32_0X10);
            tmp >>= hint::HINT_INT_BIT;
            count++;
        }
    }
    result.set_true_len();
    return result;
}
inline void HyperInt::add_sub_inplace(const HyperInt &input, bool is_add, const size_t shift) //就地加减 a+=b;a-=b,a加/减去左移位后的b，默认不移位
{
    size_t len1 = length(), len2 = input.length();
    if (is_add)
    {
        size_t result_len = std::max(len1, len2) + 1;
        reset_size(result_len);
        hint::ary_clr(data.array + len1, result_len - len1);
        change_length(result_len);
        hint::INT_64 tmp = 0;
        size_t count = shift;
        while (count < result_len)
        {
            hint::UINT_32 first_num = 0;
            if (count < len1)
            {
                first_num = data.array[count];
                tmp += first_num;
            }
            if (count - shift < len2)
            {
                tmp += input.data.array[count - shift];
            }
            else if (tmp == static_cast<hint::INT_64>(first_num))
            {
                break;
            }
            data.array[count] = static_cast<hint::UINT_32>(tmp & hint::HINT_INT32_0XFF);
            tmp >>= hint::HINT_INT_BIT;
            count++;
        }
    }
    else
    {

        size_t result_len = std::max(len1, len2);
        reset_size(result_len);
        change_length(result_len);
        hint::INT_64 tmp = 0;
        size_t count = shift;
        while (count < result_len)
        {
            hint::UINT_32 first_num = 0;
            if (count < len1)
            {
                first_num = data.array[count];
                tmp += first_num;
            }
            if (count - shift < len2)
            {
                tmp -= input.data.array[count - shift];
            }
            else if (tmp == static_cast<hint::INT_64>(first_num))
            {
                break;
            }
            data.array[count] = static_cast<hint::UINT_32>(tmp + hint::HINT_INT32_0X10);
            tmp >>= hint::HINT_INT_BIT;
            count++;
        }
    }
    set_true_len();
}
inline void HyperInt::sub_inplace(const HyperInt &input) //由减数调用,a.sub_inplace(b)->a=b-a;
{
    size_t len1 = length(), len2 = input.length();
    size_t result_len = std::max(len1, len2);
    reset_size(result_len);
    change_length(result_len);
    hint::INT_64 tmp = 0;
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
        data.array[count] = static_cast<hint::UINT_32>(tmp + hint::HINT_INT32_0X10);
        tmp >>= hint::HINT_INT_BIT;
        count++;
    }
    set_true_len();
}
inline HyperInt &HyperInt::self_half() //自身右移一位
{
    hint::UINT_32 tmp1, tmp2 = 0;
    size_t pos = length();
    while (pos)
    {
        pos--;
        tmp1 = tmp2;
        tmp2 = data.array[pos];
        tmp1 = (tmp1 << (hint::HINT_INT_BIT - 1)) + (tmp2 >> 1);
        data.array[pos] = tmp1;
    }
    set_true_len();
    return *this;
}
inline HyperInt &HyperInt::self_twice() //自身左移一位
{
    hint::UINT_64 tmp = 0;
    size_t len = length();
    for (size_t pos = 0; pos < len; pos++)
    {
        tmp += (static_cast<hint::UINT_64>(data.array[pos]) << 1);
        data.array[pos] = static_cast<hint::UINT_32>(tmp);
        tmp >>= hint::HINT_INT_BIT;
    }
    if (tmp > 0)
    {
        len++;
        reset_size(len);
        change_length(len);
        data.array[len - 1] = hint::HINT_INT32_0X01;
    }
    return *this;
}
inline HyperInt HyperInt::half() const //返回右移后的值
{
    return r_shift(1);
}
inline HyperInt HyperInt::twice() const //返回左移后的值
{
    return l_shift(1);
}
inline HyperInt HyperInt::r_shift(size_t n) const //右移n位
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
    hint::UINT_32 tmp1, tmp2 = 0;
    result.reset_size(len);
    result.change_length(len);
    while (len > 0)
    {
        len--;
        tmp1 = tmp2;
        tmp2 = data.array[len + shift];
        tmp1 = (offset ? (tmp1 << (hint::HINT_INT_BIT - offset)) : 0) + (tmp2 >> offset);
        result.data.array[len] = tmp1;
    }
    result.set_true_len();
    return result;
}
inline HyperInt HyperInt::l_shift(size_t n) const //左移n位
{
    if (n == 0)
    {
        return *this;
    }
    HyperInt result;
    size_t shift = n / 32;
    size_t offset = n % 32;
    size_t len = length();
    hint::UINT_32 tmp1 = 0, tmp2 = 0;
    result.reset_size(len + shift + 1);
    result.clear();
    result.change_length(len + shift + 1);
    for (size_t pos = 0; pos < length(); pos++)
    {
        tmp1 = tmp2;
        tmp2 = data.array[pos];
        tmp1 = (offset ? (tmp1 >> (hint::HINT_INT_BIT - offset)) : 0) + (tmp2 << offset);
        result.data.array[pos + shift] = tmp1;
    }
    if (tmp2 >> (hint::HINT_INT_BIT - offset) & offset)
    {
        result.data.array[len + shift] = tmp2 >> (hint::HINT_INT_BIT - offset);
    }
    else
    {
        result.change_length(result.length() - 1);
    }
    return result;
}
inline void HyperInt::reset_size(size_t new_size_input) //重新设定长度不小于new_size,1.5倍长度算法,在change_len()之前调用
{
    size_t size_tmp = generate_size(new_size_input);
    if (data.array == nullptr)
    {
        data.size = size_tmp;
        data.array = new hint::UINT_32[data.size];
        change_length(0);
    }
    else if (data.size != size_tmp)
    {
        data.size = size_tmp;
        data.array = hint::ary_realloc(data.array, data.size);
        change_length(std::min(length(), data.size));
    }
}
inline HyperInt HyperInt::abs() const //返回绝对值
{
    HyperInt result(*this);
    result.neg_sign(false);
    return result;
}
inline hint::UINT_32 HyperInt::first_int32() const
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
inline hint::UINT_64 HyperInt::first_int64() const
{
    if (length() >= 2)
    {
        hint::UINT_64 result = data.array[length() - 1];
        result <<= hint::HINT_INT_BIT;
        result += data.array[length() - 2];
        return result;
    }
    else
    {
        return static_cast<hint::UINT_64>(first_int32());
    }
}
inline hint::INT_32 HyperInt::abs_compare(const HyperInt &input, hint::INT_64 shift) const //自身和input移位shift比较，大于返回1，小于返回-1，等于返回0
{
#if SIZE_T_BITS == 64
    hint::INT_64 len1 = static_cast<hint::INT_64>(length());
    hint::INT_64 len2 = static_cast<hint::INT_64>(input.length()) + shift;
#elif SIZE_T_BITS == 32
    hint::INT_32 len1 = static_cast<hint::INT_32>(length());
    hint::INT_32 len2 = static_cast<hint::INT_32>(input.length()) + shift;
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
    hint::UINT_32 num1 = 0, num2 = 0;
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
inline bool HyperInt::abs_larger(const HyperInt &input) const //绝对值是否大于input
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

    hint::UINT_32 num1 = 0, num2 = 0;
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
inline bool HyperInt::abs_smaller(const HyperInt &input) const //绝对值是否小于input
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
    hint::UINT_32 num1 = 0, num2 = 0;
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
inline bool HyperInt::abs_equal(const HyperInt &input) const //绝对值是否等于input
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
inline bool HyperInt::equal_to_z() const //判定是否为零
{
    return (0 == length());
}
inline bool HyperInt::is_even() const //判断是否为偶数
{
    if (length() == 0)
    {
        return true;
    }
    return !static_cast<bool>(data.array[0] & 1);
}
inline bool HyperInt::is_odd() const //判断是否为奇数
{
    if (length() == 0)
    {
        return false;
    }
    return (bool)data.array[0] & 1;
}

//逻辑运算
inline bool HyperInt::operator==(const HyperInt &input) const
{
    if (is_neg() != input.is_neg())
    {
        return false;
    }
    else
    {
        return abs_equal(input);
    }
}
template <typename T>
inline bool HyperInt::operator==(const T &input) const
{
    if (is_neg() != hint::is_neg(input))
    {
        return false;
    }
    else
    {
        return abs_equal(HyperInt(input));
    }
}
inline bool HyperInt::operator!=(const HyperInt &input) const
{
    if (is_neg() != input.is_neg())
    {
        return true;
    }
    else
    {
        return !abs_equal(input);
    }
}
template <typename T>
inline bool HyperInt::operator!=(const T &input) const
{
    if (is_neg() != hint::is_neg(input))
    {
        return true;
    }
    else
    {
        return !abs_equal(HyperInt(input));
    }
}
inline bool HyperInt::operator>(const HyperInt &input) const
{
    if (is_neg() != input.is_neg())
    {
        return input.is_neg();
    }
    else
    {
        return is_neg() != abs_larger(input);
    }
}
template <typename T>
inline bool HyperInt::operator>(const T &input) const
{
    if (is_neg() != hint::is_neg(input))
    {
        return hint::is_neg(input);
    }
    else
    {
        return is_neg() != abs_larger(HyperInt(input));
    }
}
inline bool HyperInt::operator>=(const HyperInt &input) const
{
    return !(*this < input);
}
template <typename T>
inline bool HyperInt::operator>=(const T &input) const
{
    return !(*this < input);
}
inline bool HyperInt::operator<(const HyperInt &input) const
{
    if (is_neg() != input.is_neg())
    {
        return is_neg();
    }
    else
    {
        return is_neg() != abs_smaller(input);
    }
}
template <typename T>
inline bool HyperInt::operator<(const T &input) const
{
    if (is_neg() != hint::is_neg(input))
    {
        return is_neg();
    }
    else
    {
        return is_neg() != abs_smaller(HyperInt(input));
    }
}
inline bool HyperInt::operator<=(const HyperInt &input) const
{
    return !(*this > input);
}
template <typename T>
inline bool HyperInt::operator<=(const T &input) const
{
    return !(*this > input);
}

//算术运算
inline HyperInt HyperInt::operator+(const HyperInt &input) const
{
    HyperInt result;
    if (!is_neg() != input.is_neg()) //是否同号
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
template <typename T>
inline HyperInt HyperInt::operator+(T input) const
{
    return *this + HyperInt(input);
}
inline HyperInt HyperInt::operator+() const
{
    return *this;
}

inline HyperInt HyperInt::operator-(const HyperInt &input) const
{
    HyperInt result;
    if (is_neg() != input.is_neg()) //是否异号
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
template <typename T>
inline HyperInt HyperInt::operator-(T input) const
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
        size_t min_len = std::min(length(), input.length());
        size_t sum_len = length() + input.length();
        size_t factor = 4;
        if (sum_len <= hint::HINT_FFT_MAX)
        {
            factor = 4;
        }
        else if (sum_len <= hint::HINT_QUAL_FFT_MAX)
        {
            factor = 9;
        }
        else
        {
            factor = 43;
        }
        if (sum_len <= 120 || (min_len <= factor * std::ceil(std::log2(sum_len))))
        {
            return normal_multiply(input);
        }
        else if (sum_len <= hint::HINT_QUAL_FFT_MAX)
        {
            return fft_multiply(input);
        }
        else if (sum_len < hint::HINT_NTT_MAX / 8)
        {
            return karatsuba_multiply(input);
        }
        else if (sum_len <= hint::HINT_NTT_MAX)
        {
            return ntt_multiply(input);
        }
        else
        {
            return karatsuba_multiply(input);
        }
    }
}
template <typename T>
HyperInt HyperInt::operator*(T input) const
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
    bool neg = hint::is_neg(input);
    size_t len = length();
    size_t result_len = len + 2;
    result.reset_size(result_len);
    result.change_length(result_len);
    result.neg_sign(is_neg() != neg);
    hint::UINT_32 tmp = 0, tmp1 = 0, tmp2 = 0;
    hint::UINT_64 abs_input = 0, sum1 = 0, sum2 = 0;
    if (neg)
    {
        abs_input = static_cast<hint::UINT_64>(std::abs(static_cast<hint::INT_64>(input)));
    }
    else
    {
        abs_input = static_cast<hint::UINT_64>(input);
    }
    const hint::UINT_64 input_num1 = abs_input & hint::HINT_INT32_0XFF;
    const hint::UINT_64 input_num2 = abs_input >> hint::HINT_INT_BIT;

    tmp2 = data.array[0];
    sum1 = tmp2 * input_num1;
    result.data.array[0] = sum1 & hint::HINT_INT32_0XFF;
    sum1 >>= hint::HINT_INT_BIT;
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
        tmp = sum1 & hint::HINT_INT32_0XFF;
        sum2 += input_num1 * tmp2 + tmp;
        tmp = sum2 & hint::HINT_INT32_0XFF;
        result.data.array[pos] = tmp;
        sum1 >>= hint::HINT_INT_BIT;
        sum2 >>= hint::HINT_INT_BIT;
    }
    result.set_true_len();
    return result;
}

inline HyperInt HyperInt::operator/(const HyperInt &input) const
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
template <typename T>
inline HyperInt HyperInt::operator/(T input) const
{
    return *this / HyperInt(input);
}

inline HyperInt HyperInt::operator%(const HyperInt &input) const
{
    HyperInt tmp = *this / input;
    tmp *= input;
    return *this - tmp;
}
template <typename T>
inline HyperInt HyperInt::operator%(T input) const
{
    return *this % HyperInt(input);
}

inline HyperInt &HyperInt::operator+=(const HyperInt &input)
{
    if (!is_neg() != input.is_neg()) //是否同号
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
template <typename T>
inline HyperInt &HyperInt::operator+=(T input)
{
    return *this += HyperInt(input);
}

inline HyperInt &HyperInt::operator-=(const HyperInt &input)
{
    if (is_neg() != input.is_neg()) //是否异号
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
template <typename T>
inline HyperInt &HyperInt::operator-=(T input)
{
    return *this *= HyperInt(input);
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
template <typename T>
inline HyperInt &HyperInt::operator*=(T input)
{
    if (input == 1)
    {
        return *this;
    }
    if (equal_to_z() || input == 0)
    {
        *this = HyperInt();
        return *this;
    }
    bool neg = hint::is_neg(input);
    size_t len = length();
    size_t result_len = len + 2;
    reset_size(result_len);
    change_length(result_len);
    neg_sign(is_neg() != neg);
    hint::UINT_32 tmp = 0, tmp1 = 0, tmp2 = 0;
    hint::UINT_64 abs_input = 0, sum1 = 0, sum2 = 0;

    if (neg)
    {
        abs_input = static_cast<hint::UINT_64>(std::abs(static_cast<hint::INT_64>(input)));
    }
    else
    {
        abs_input = static_cast<hint::UINT_64>(input);
    }
    const hint::UINT_64 input_num1 = abs_input & hint::HINT_INT32_0XFF;
    const hint::UINT_64 input_num2 = abs_input >> hint::HINT_INT_BIT;

    tmp2 = data.array[0];
    sum1 = tmp2 * input_num1;
    data.array[0] = sum1 & hint::HINT_INT32_0XFF;
    sum1 >>= hint::HINT_INT_BIT;
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
        tmp = sum1 & hint::HINT_INT32_0XFF;
        sum2 += input_num1 * tmp2 + tmp;
        tmp = sum2 & hint::HINT_INT32_0XFF;
        data.array[pos] = tmp;
        sum1 >>= hint::HINT_INT_BIT;
        sum2 >>= hint::HINT_INT_BIT;
    }
    set_true_len();
    return *this;
}

HyperInt &HyperInt::operator/=(const HyperInt &input)
{
    *this = *this / input;
    return *this;
}
template <typename T>
HyperInt &HyperInt::operator/=(T input)
{
    return *this /= HyperInt(input);
}

HyperInt &HyperInt::operator%=(const HyperInt &input)
{
    HyperInt tmp = *this / input;
    *this = *this - (tmp * input);
    return *this;
}
template <typename T>
HyperInt &HyperInt::operator%=(T input)
{
    return *this %= HyperInt(input);
}

HyperInt HyperInt::operator++(int)
{
    HyperInt tmp(*this);
    *this += HyperInt(1);
    return tmp;
}
HyperInt &HyperInt::operator++()
{
    *this += HyperInt(1);
    return *this;
}
HyperInt HyperInt::operator--(int)
{
    HyperInt tmp(*this);
    *this -= HyperInt(1);
    return tmp;
}
HyperInt &HyperInt::operator--()
{
    *this -= HyperInt(1);
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
    size_t result_len = std::max(len1, len2);
    size_t min_len = std::min(len1, len2);
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
    size_t result_len = std::min(len1, len2);
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
    size_t result_len = std::max(len1, len2);
    size_t min_len = std::min(len1, len2);
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
    size_t result_len = std::max(len1, len2);
    size_t min_len = std::min(len1, len2);
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
    size_t result_len = std::min(len1, len2);
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
    size_t result_len = std::max(len1, len2);
    size_t min_len = std::min(len1, len2);
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
bool operator!=(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 != input1;
}
bool operator==(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 == input1;
}
bool operator>(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 < input1;
}
bool operator>=(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 <= input1;
}
bool operator<(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 > input1;
}
bool operator<=(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 >= input1;
}
HyperInt operator+(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 + input1;
}
HyperInt operator-(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 - input1;
}
HyperInt operator*(const hint::INT_64 &input1, const HyperInt &input2)
{
    return input2 * input1;
}
HyperInt operator/(const hint::INT_64 &input1, const HyperInt &input2)
{
    return HyperInt(input1) / input2;
}
HyperInt operator%(const hint::INT_64 &input1, const HyperInt &input2)
{
    return HyperInt(input1) % input2;
}
hint::INT_64 &operator+=(hint::INT_64 &input1, const HyperInt &input2)
{
    return input1 += input2.to_int64();
}
hint::INT_64 &operator-=(hint::INT_64 &input1, const HyperInt &input2)
{
    return input1 -= input2.to_int64();
}
hint::INT_64 &operator*=(hint::INT_64 &input1, const HyperInt &input2)
{
    return input1 *= input2.to_int64();
}
hint::INT_64 &operator/=(hint::INT_64 &input1, const HyperInt &input2)
{
    return input1 /= input2.to_int64();
}
hint::INT_64 &operator%=(hint::INT_64 &input1, const HyperInt &input2)
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

HyperInt classic_factorial(hint::UINT_64 end, hint::UINT_64 start = 1) //累乘阶乘,可计算排列数A(n,m) n!/(n-m)!
{
    HyperInt result = 1;
    if (end < start)
    {
        return result;
    }
    for (hint::UINT_64 i = start; i <= end; i++)
    {
        result *= i;
    }
    result.set_true_len();
    return result;
}
HyperInt factorial(hint::UINT_64 end, hint::UINT_64 start = 1, const hint::UINT_32 rec_level = 0) //递归拆分,可计算排列数A(n,m) n!/(n-m)!
{
    if (end < start)
    {
        return HyperInt(1);
    }
    hint::UINT_64 len = end - start;
    if (len < 120)
    {
        HyperInt result = 1;
        for (hint::UINT_64 i = start; i <= end; i++)
        {
            result *= i;
        }
        return result;
    }
    hint::UINT_64 mid = start + (len * 3 / 5);
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
HyperInt combination(hint::UINT_64 n, hint::UINT_64 m) // return C(n,m)组合数公式n!/((n-m)!m!)
{
    if (m > hint::half(n))
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
HyperInt randHyperInt(size_t len)
{
    HyperInt result;
    if (len == 0)
    {
        return result;
    }
    result.reset_size(len);
    result.change_length(len);
    std::random_device seed;
    std::default_random_engine rand_num(seed());
    std::uniform_int_distribution<hint::UINT_32> uni(1, UINT32_MAX);
    result.data.array[len - 1] = uni(rand_num);
    if (len == 1)
    {
        result.set_true_len();
        return result;
    }
    uni = std::uniform_int_distribution<hint::UINT_32>(0, UINT32_MAX);
    for (size_t i = 0; i < len; i++)
    {
        result.data.array[i] = uni(rand_num);
    }
    result.set_true_len();
    return result;
}
#endif