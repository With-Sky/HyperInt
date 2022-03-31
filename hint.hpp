#ifndef HINT_HPP
#define HINT_HPP

#include<iostream>
#include<future>
#include<cstring>

#define _MAX_(x,y) (x>y?x:y)
#define _MIN_(x,y) (x<y?x:y)
#define _NEG_(x) (x<0)
#define _ABS_(x) (x<0?(-x):x)
#define _TWICE_(x) x<<=1
#define _HALF_(x) x>>=1
#define _ARYCOPY_(target,source,len) memcpy(target,source,len*sizeof(*target))

#define INT_BIT 32
#define INT32_0X7F 0x7FFFFFFF
#define INT32_0X01 0x00000001
#define INT32_0X80 0X80000000

#define UINT_32 unsigned int
#define UINT_64 unsigned long long
#define INT_64 long long
#define REG register

using namespace std;

struct h_int
{
    bool neg = false;
    UINT_32* array = NULL;
    size_t len, size;
};

class HyperInt
{
private:
    h_int data;
    short shorten_count, size_shorten_max = 64;//设定array长度不缩短次数阈值
    size_t true_len() const
    {
        REG size_t t_len = _MAX_(2, data.len);
        while (t_len > 2 && data.array[t_len - 1] == 0)
        {
            t_len--;
        }
        return t_len;
    }
    void quick_self_twice()
    {
        REG UINT_32 tmp;
        if (data.array[0])
        {
            tmp = data.array[0];
            data.array[0] = tmp << 1;
            data.array[1] = (tmp & INT32_0X80) >> (INT_BIT - 1);
        }
        else
        {
            tmp = data.array[data.len - 1];
            data.array[data.len - 1] = tmp << 1;
            if (data.size > data.len && (tmp << 1) == 0)
            {
                data.array[data.len] = INT32_0X01;
                data.len++;
            }
        }
    }
    void quick_self_half()
    {

        if (data.array[0])
        {
            data.array[0] = data.array[0] >> 1;
        }
        else
        {
            REG UINT_32 tmp = data.array[data.len - 1];
            data.array[data.len - 1] = tmp >> 1;
            data.array[data.len - 2] = (tmp & INT32_0X01) << (INT_BIT - 1);
        }
    }
    HyperInt quick_twice()
    {
        HyperInt result(*this);
        result.quick_self_twice();
        result.data.len = result.true_len();
        return result;
    }
    HyperInt quick_half()
    {
        HyperInt result(*this);
        result.quick_self_half();
        result.data.len = result.true_len();
        return result;
    }
    HyperInt add(h_int in1, h_int in2) const
    {
        if (!in1.len)
            return HyperInt(in2);
        else if (!in2.len)
            return HyperInt(in1);
        REG HyperInt result;
        result.reset_size(_MAX_(in1.len, in2.len) + 1);
        result.data.len = _MAX_(in1.len, in2.len) + 1;
        REG INT_64 tmp = 0;
        REG size_t pos1 = 0, pos2 = 0, count = 0;
        while (count < result.data.len)
        {
            if (pos1 < in1.len)
            {
                tmp += in1.array[pos1];
                pos1++;
            }
            if (pos2 < in2.len)
            {
                tmp += in2.array[pos2];
                pos2++;
            }
            result.data.array[count] = tmp;
            tmp >>= INT_BIT;
            count++;
        }
        result.data.len = result.true_len();
        return result;
    }
    HyperInt multiply(const HyperInt& input) const
    {
        REG HyperInt result;
        result.reset_size(data.len + input.data.len);
        result.clear();
        result.data.len = data.len + input.data.len;
        result.data.neg = data.neg ^ input.data.neg;
        REG UINT_64 tmp, sum;
        for (REG size_t pos1 = 0; pos1 < data.len; pos1++)
        {
            for (REG size_t pos2 = 0; pos2 < input.data.len; pos2++)
            {
                tmp = (UINT_64)data.array[pos1] * input.data.array[pos2];
                for (REG size_t pos3 = pos1 + pos2; pos3 < result.data.len; pos3++)
                {
                    sum = tmp + result.data.array[pos3];
                    result.data.array[pos3] = sum;
                    if ((sum >> 32) > 0)
                    {
                        tmp = sum >> 32;
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
    }
    HyperInt quick_multiply(h_int in1, h_int in2) const
    {
        REG HyperInt result;
        if (!in1.len || !in2.len)
            return result;
        if ((!HyperInt(in1)) || (!HyperInt(in2)))
            return result;
        result.reset_size(in1.len + in2.len);
        result.data.len = in1.len + in2.len;
        result.data.neg = in1.neg ^ in2.neg;
        if (_MIN_(in1.len, in2.len) < 3)
        {
            //cout<<"basic *"<<endl;
            REG UINT_64 tmp, sum;
            for (REG size_t pos1 = 0; pos1 < in1.len; pos1++)
            {
                for (REG size_t pos2 = 0; pos2 < in2.len; pos2++)
                {
                    tmp = (UINT_64)in1.array[pos1] * in2.array[pos2];
                    for (REG size_t pos3 = pos1 + pos2; pos3 < result.data.len; pos3++)
                    {
                        sum = tmp + result.data.array[pos3];
                        result.data.array[pos3] = sum;
                        if ((sum >> 32) > 0)
                        {
                            tmp = sum >> 32;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
            // REG UINT_64 tmp,sum;
            // tmp=(UINT_64)in1.array[0]*in2.array[0];
            // result.data.array[0]=tmp;
            // tmp>>=32;
            // tmp=tmp+(UINT_64)in1.array[0]*in2.array[1]+(UINT_64)in1.array[1]*in2.array[0];
            // result.data.array[1]=tmp;
            // tmp>>=32;
            // tmp+=(UINT_64)in1.array[1]*in2.array[1];
            // result.data.array[2]=tmp;
            // tmp>>=32;
            // result.data.array[3]=tmp;
            // tmp>>=32;
            // result.data.array[4]=tmp;
        }
        else
        {
            REG INT_64 tmp = 0;
            static UINT_32 threads_count;
            REG h_int h_sub_a, h_sub_b, h_sub_c, h_sub_d;
            REG size_t dig_len = _MAX_(in1.len, in2.len) >> 1, pos1, pos2, pos3, pos4, count;
            h_sub_a.array = h_sub_b.array = in1.array;
            h_sub_c.array = h_sub_d.array = in2.array;

            h_sub_a.len = _MIN_(in1.len, dig_len);
            h_sub_b.len = (in1.len > dig_len ? in1.len - dig_len : 0);
            h_sub_b.array += (h_sub_b.len ? dig_len : 0);

            h_sub_c.len = _MIN_(in2.len, dig_len);
            h_sub_d.len = (in2.len > dig_len ? in2.len - dig_len : 0);
            h_sub_d.array += (h_sub_d.len ? dig_len : 0);

            // future<HyperInt> th1=async(quick_multiply,h_sub_a,h_sub_b);
            // future<HyperInt> th2=async(quick_multiply,h_sub_a,h_sub_b);
            // future<HyperInt> th3=async(quick_multiply,h_sub_a,h_sub_b);

            REG HyperInt tmp_e = quick_multiply(h_sub_a, h_sub_c);
            REG HyperInt tmp_f = quick_multiply(h_sub_b, h_sub_d);
            REG HyperInt tmp_g = add(h_sub_a, h_sub_b) * add(h_sub_c, h_sub_d);//递归运算

            _ARYCOPY_(result.data.array, tmp_e.data.array, dig_len);
            pos1 = pos2 = pos3 = 0;
            pos4 = count = dig_len;
            //相加
            while (count < dig_len << 1)
            {
                if (pos1 < tmp_g.data.len)
                {
                    tmp += tmp_g.data.array[pos1];
                    pos1++;
                }
                if (pos2 < tmp_e.data.len)
                {
                    tmp -= tmp_e.data.array[pos2];
                    pos2++;
                }
                if (pos3 < tmp_f.data.len)
                {
                    tmp -= tmp_f.data.array[pos3];
                    pos3++;
                }
                if (pos4 < tmp_f.data.len)
                {
                    tmp += tmp_f.data.array[pos4];
                    pos4++;
                }
                result.data.array[count] = tmp;
                tmp >>= INT_BIT;
                count++;
            }
            pos1 = 0;
            pos2 = pos3 = pos4 = dig_len;
            count = dig_len << 1;
            while (count < result.data.len)
            {
                if (pos1 < tmp_f.data.len)
                {
                    tmp += tmp_f.data.array[pos1];
                    pos1++;
                }
                if (pos2 < tmp_g.data.len)
                {
                    tmp += tmp_g.data.array[pos2];
                    pos2++;
                }
                if (pos3 < tmp_e.data.len)
                {
                    tmp -= tmp_e.data.array[pos3];
                    pos3++;
                }
                if (pos4 < tmp_f.data.len)
                {
                    tmp -= tmp_f.data.array[pos4];
                    pos4++;
                }
                result.data.array[count] = tmp;
                tmp >>= INT_BIT;
                count++;
            }
        }
        result.data.len = result.true_len();
        return result;
    }
public:
    // HyperInt test(HyperInt in1, HyperInt in2)
    // {
    //     h_int a, b;
    //     a = in1.data;
    //     b = in2.data;
    //     return multiply(a);
    // }
    ~HyperInt()
    {
        if (data.array != NULL)
        {
            delete[] data.array;
            data.array = NULL;
        }
        //cout<<"_deleted_"<<endl;
    }//析构函数
    HyperInt()
    {
        data.len = data.size = 2;
        data.neg = false;
        shorten_count = 0;
        data.array = new UINT_32[data.size];
        data.array[1] = data.array[0] = 0;
        //cout<<"_1_"<<endl;
    }//无参数构造
    HyperInt(const HyperInt& input)
    {
        if (this != &input)
        {
            data.len = input.true_len();
            data.size = 2;
            while (data.size < data.len)
            {
                _TWICE_(data.size);
            }
            data.neg = input.data.neg;
            shorten_count = 0;
            data.array = new UINT_32[data.size];
            _ARYCOPY_(data.array, input.data.array, data.len);
        }
        //cout<<"_2_"<<endl;
    }//HyperInt拷贝构造
    HyperInt(HyperInt&& input)
    {
        if (this != &input)
        {
            data.len = input.true_len();
            data.size = input.data.size;
            shorten_count = 0;
            data.neg = input.data.neg;
            data.array = input.data.array;
            input.data.array = NULL;
        }
        //cout<<"_3_"<<endl;
    }//HyperInt移动构造
    HyperInt(const h_int& input)
    {
        data.len = _MAX_(input.len, 2);
        data.size = 2;
        while (data.size < data.len)
        {
            _TWICE_(data.size);
        }
        data.neg = input.neg;
        shorten_count = 0;
        data.array = new UINT_32[data.size];
        _ARYCOPY_(data.array, input.array, data.len);
        data.len = true_len();
        //cout<<"_4_"<<endl;
    }//h_int参数拷贝构造
    HyperInt(const INT_64 input)
    {
        data.len = data.size = 2;
        data.neg = (input < 0);
        shorten_count = 0;
        data.array = new UINT_32[data.size];
        data.array[0] = _ABS_(input);
        data.array[1] = (_ABS_(input) >> INT_BIT) & INT32_0X7F;
        //cout<<"_5_"<<endl;
    }//INT_64参数构造
    HyperInt(const string input)
    {
        string_in(input);
        //cout<<"_6_"<<endl;
    }//string参数构造
    HyperInt operator=(const HyperInt& input)
    {
        //cout<<"_=_1_"<<endl;
        if (this != &input)
        {
            data.len = input.true_len();
            if (data.size < input.data.len || shorten_count >= size_shorten_max)
            {
                data.size = 2;
                while (data.size < input.data.len)
                {
                    _TWICE_(data.size);
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
            data.neg = input.data.neg;
            _ARYCOPY_(data.array, input.data.array, data.len);
        }
        return *this;
    }//HyperInt拷贝赋值
    HyperInt operator=(HyperInt&& input)
    {
        //cout<<"_=_2_"<<endl;
        if (this != &input)
        {
            data.len = input.true_len();
            data.size = input.data.size;
            data.neg = input.data.neg;
            delete[] data.array;
            data.array = input.data.array;
            input.data.array = NULL;
        }
        return *this;
    }//HyperInt移动赋值
    HyperInt operator=(const h_int& input)
    {
        //cout<<"_=_3_"<<endl;
        data.len = _MAX_(input.len, 2);
        if (data.size < input.len || shorten_count >= size_shorten_max)
        {
            data.size = 2;
            while (data.size < input.len)
            {
                _TWICE_(data.size);
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
        data.len = true_len();
        return *this;
    }//h_int拷贝赋值
    HyperInt operator=(h_int&& input)
    {
        //cout<<"_=_3_"<<endl;
        data.len = _MAX_(input.len, 2);
        data.size = input.size;
        data.neg = input.neg;
        data.array = input.array;
        input.array = NULL;
        data.len = true_len();
        return *this;
    }//h_int移动赋值
    HyperInt operator=(const INT_64 input)
    {
        //cout<<"_=_3_"<<endl;
        data.len = 2;
        data.neg = _NEG_(input);
        data.array[0] = _ABS_(input);
        data.array[1] = (_ABS_(input) >> INT_BIT) & INT32_0X7F;
        return *this;
    }//longlong赋值
    HyperInt operator=(const string input)
    {
        //cout<<"_=_4_"<<endl;
        string_in(input);
        return *this;
    }//string赋值
    //基本操作
    void reset_size(const size_t& new_size)
    {
        data.size = 2;
        while (data.size < new_size)
        {
            _TWICE_(data.size);
        }
        data.len = _MIN_(data.len, data.size);
        shorten_count = 0;
        UINT_32* new_array = new UINT_32[data.size];
        _ARYCOPY_(new_array, data.array, data.len);
        delete[] data.array;
        data.array = new_array;
    }//设定长度为2^n且不小于new_size
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
    void console_in()
    {
        clear();
        char tmp = '0';
        while (tmp != '\n' && '0' <= tmp <= '9')
        {
            if (tmp != '\n')
            {
                *this *= 10;
                *this += ((long long)(tmp - '0'));
            }
            tmp = getchar();
        }
    }
    void string_in(string s)
    {
        clear();
        char tmp = '0';
        for (string::iterator i = s.begin(); i != s.end(); i++)
        {
            if ('0' <= *i <= '9')
            {
                *this *= 10;
                *this += ((long long)(*i - '0'));
            }
        }
    }
    void console_out_hex() const
    {
        if (data.neg)
        {
            printf("-");
        }
        for (REG size_t i = data.len - 1; i > 0; i--)
        {
            printf("%X ", data.array[i]);
        }
        printf("%X", data.array[0]);
    }
    void console_out_dec() const
    {
        REG HyperInt tmp1(*this), tmp2;
        REG size_t pos = 0;
        char* result = new char[10 * data.len];
        while (tmp1 > 0)
        {
            tmp2 = tmp1 / 10;
            result[pos] = (tmp1 - tmp2 * 10).to_int64() + '0';
            tmp1 = move(tmp2);
            pos++;
        }
        result[pos] = 0;
        while (pos > 0)
        {
            putchar(result[pos - 1]);
            pos--;
        }
    }
    HyperInt& self_half()
    {
        REG UINT_32 tmp1, tmp2 = data.array[0];
        for (REG size_t i = 1; i < data.len; i++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[i];
            _HALF_(tmp1);
            tmp1 += ((tmp2 & INT32_0X01) << INT_BIT - 1);
            data.array[i - 1] = tmp1;
        }
        data.array[data.len - 1] = tmp2 >> 1;
        return *this;
    }// /2
    HyperInt& self_twice()
    {
        REG UINT_32 tmp1 = data.array[0], tmp2;
        data.array[0] = tmp1 << 1;
        for (REG size_t i = 1; i < data.len; i++)
        {
            tmp2 = data.array[i];
            tmp1 = (tmp2 << 1) + ((tmp1 & INT32_0X80) >> INT_BIT - 1);
            data.array[i] = tmp1;
            tmp1 = tmp2;
        }
        if (data.size > data.len)
        {
            data.array[data.len] = (tmp1 & INT32_0X80) >> (INT_BIT - 1);
            data.len++;
        }
        return *this;
    }// *2
    HyperInt half() const
    {
        HyperInt result;
        result.reset_size(data.len);
        result.data.len = data.len;
        REG UINT_32 tmp1, tmp2 = data.array[0];
        for (REG size_t i = 1; i < result.data.len; i++)
        {
            tmp1 = tmp2;
            tmp2 = data.array[i];
            _HALF_(tmp1);
            tmp1 += ((tmp2 & INT32_0X01) << INT_BIT - 1);
            result.data.array[i - 1] = tmp1;
        }
        result.data.array[data.len - 1] = tmp2 >> 1;
        return result;
    }// /2
    HyperInt twice() const
    {
        HyperInt result;
        result.reset_size(data.len);
        result.data.len = data.len;
        REG UINT_32 tmp1 = data.array[0], tmp2;
        result.data.array[0] = tmp1 << 1;
        for (REG size_t i = 1; i < result.data.len; i++)
        {
            tmp2 = data.array[i];
            tmp1 = (tmp2 << 1) + ((tmp1 & INT32_0X80) >> INT_BIT - 1);
            result.data.array[i] = tmp1;
            tmp1 = tmp2;
        }
        if (result.data.size > result.data.len)
        {
            result.data.array[result.data.len] = (tmp1 & INT32_0X80) >> (INT_BIT - 1);
            result.data.len++;
        }
        return result;
    }// *2
    HyperInt abs() const
    {
        HyperInt result(*this);
        result.data.neg = false;
        return result;
    }
    bool abs_larger(const HyperInt& input) const
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
    bool abs_smaller(const HyperInt& input) const
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
    bool abs_equal(const HyperInt& input) const
    {
        REG size_t t_len1 = true_len(), t_len2 = input.true_len();
        if (t_len1 != t_len2)
            return false;
        else
        {
            REG UINT_32 num1, num2;
            while (t_len1 > 0)
            {
                t_len1--;
                if (data.array[t_len1] != input.data.array[t_len1])
                    return false;
            }
            return true;
        }
    }//比较是否为相反数
    void clear()
    {
        memset(data.array, 0, data.size * sizeof(*data.array));
    }

    //逻辑运算
    bool operator==(const HyperInt& input) const;
    bool operator==(const INT_64& input) const;
    bool operator>(const HyperInt& input) const;
    bool operator>(const INT_64& input) const;
    bool operator>=(const HyperInt& input) const;
    bool operator>=(const INT_64& input) const;
    bool operator<(const HyperInt& input) const;
    bool operator<(const INT_64& input) const;
    bool operator<=(const HyperInt& input) const;
    bool operator<=(const INT_64& input) const;
    // bool operator!();
    bool operator!() const;
    //operator bool();

    //算术运算
    HyperInt operator+(const HyperInt& input) const;
    HyperInt operator+() const;
    HyperInt operator-(const HyperInt& input) const;
    HyperInt operator-() const;
    HyperInt operator*(const HyperInt& input) const;
    HyperInt operator/(const HyperInt& input) const;
    HyperInt operator%(const HyperInt& input) const;
    HyperInt operator+=(const HyperInt& input);
    HyperInt operator-=(const HyperInt& input);
    HyperInt operator*=(const HyperInt& input);
    HyperInt operator/=(const HyperInt& input);
    HyperInt operator%=(const HyperInt& input);
    HyperInt operator ++ (int);
    HyperInt& operator ++ ();
    HyperInt operator -- (int);
    HyperInt& operator -- ();
    HyperInt power(UINT_32 n) const;

    //友元函数
    friend bool operator==(const long long& input1, const HyperInt& input2);
    friend bool operator>(const long long& input1, const HyperInt& input2);
    friend bool operator>=(const long long& input1, const HyperInt& input2);
    friend bool operator<(const long long& input1, const HyperInt& input2);
    friend bool operator<=(const long long& input1, const HyperInt& input2);
    friend HyperInt operator+(const long long& input1, const HyperInt& input2);
    friend HyperInt operator-(const long long& input1, const HyperInt& input2);
    friend HyperInt operator*(const long long& input1, const HyperInt& input2);
    friend HyperInt operator/(const long long& input1, const HyperInt& input2);
    friend HyperInt operator%(const long long& input1, const HyperInt& input2);
    friend HyperInt operator+=(long long& input1, const HyperInt& input2);
    friend HyperInt operator-=(long long& input1, const HyperInt& input2);
    friend HyperInt operator*=(long long& input1, const HyperInt& input2);
    friend HyperInt operator/=(long long& input1, const HyperInt& input2);
    friend HyperInt operator%=(long long& input1, const HyperInt& input2);
    friend HyperInt abs(HyperInt input);
    friend void print(HyperInt& input);
};
//逻辑运算
inline bool HyperInt::operator==(const HyperInt& input) const
{
    if (data.neg ^ input.data.neg)
        return false;
    else
        return abs_equal(input);
}
inline bool HyperInt::operator==(const INT_64& input) const
{
    if (data.neg ^ _NEG_(input))
        return false;
    else
        return abs_equal(input);
}
bool HyperInt::operator>(const HyperInt& input) const
{
    if (data.neg ^ input.data.neg)
        return input.data.neg;
    else
        return data.neg ^ abs_larger(input);
}
inline bool HyperInt::operator>(const INT_64& input) const
{
    return *this > HyperInt(input);
}
inline bool HyperInt::operator>=(const HyperInt& input) const
{
    return !(*this < input);
}
inline bool HyperInt::operator>=(const INT_64& input) const
{
    return !(*this < HyperInt(input));
}
bool HyperInt::operator<(const HyperInt& input) const
{
    if (data.neg ^ input.data.neg)
        return data.neg;
    else
        return data.neg ^ abs_smaller(input);
}
inline bool HyperInt::operator<(const INT_64& input) const
{
    return *this < HyperInt(input);
}
inline bool HyperInt::operator<=(const HyperInt& input) const
{
    return !(*this > input);
}
inline bool HyperInt::operator<=(const INT_64& input) const
{
    return !(*this > HyperInt(input));
}
inline bool HyperInt::operator!() const
{
    if (data.len > 2)
        return false;
    return(data.array[0] == 0 && data.array[1] == 0);
}
// inline HyperInt::operator bool()
// {
//     if(data.len>2)
//         return true;
//     return(data.array[0]!=0||data.array[1]!=0);
// }

//算术运算
HyperInt HyperInt::operator+(const HyperInt& input) const//算术运算
{
    REG HyperInt result;
    result.reset_size(_MAX_(data.len, input.data.len) + 1);
    result.data.len = _MAX_(data.len, input.data.len) + 1;
    REG INT_64 tmp = 0;
    REG size_t pos1 = 0, pos2 = 0, count = 0;
    if (!data.neg ^ input.data.neg)
    {
        result.data.neg = data.neg;
        while (count < result.data.len)
        {
            if (pos1 < data.len)
            {
                tmp += data.array[pos1];
                pos1++;
            }
            if (pos2 < input.data.len)
            {
                tmp += input.data.array[pos2];
                pos2++;
            }
            result.data.array[count] = tmp;
            tmp >>= INT_BIT;
            count++;
        }
    }
    else
    {
        if (abs_equal(input))
        {
            result.data.len = result.true_len();
            return result;
        }
        else if (abs_larger(input))
        {
            result.data.neg = data.neg;
            while (count < result.data.len)
            {
                if (pos1 < data.len)
                {
                    tmp += data.array[pos1];
                    pos1++;
                }
                if (pos2 < input.data.len)
                {
                    tmp -= input.data.array[pos2];
                    pos2++;
                }
                result.data.array[count] = tmp;
                tmp >>= INT_BIT;
                count++;
            }
        }
        else
        {
            result.data.neg = !data.neg;
            while (count < result.data.len)
            {
                if (pos2 < input.data.len)
                {
                    tmp += data.array[pos2];
                    pos2++;
                }
                if (pos1 < data.len)
                {
                    tmp -= input.data.array[pos1];
                    pos1++;
                }
                result.data.array[count] = tmp;
                tmp >>= INT_BIT;
                count++;
            }
        }
    }
    result.data.len = result.true_len();
    return result;
}
inline HyperInt HyperInt::operator+() const
{
    //cout<<"_+_"<<endl;
    return *this;
}
HyperInt HyperInt::operator-(const HyperInt& input) const
{
    REG HyperInt result;
    result.reset_size(_MAX_(data.len, input.data.len) + 1);
    result.data.len = _MAX_(data.len, input.data.len) + 1;
    REG INT_64 tmp = 0;
    REG size_t pos1 = 0, pos2 = 0, count = 0;
    if (data.neg ^ input.data.neg)
    {
        result.data.neg = data.neg;
        while (count < result.data.len)
        {
            if (pos1 < data.len)
            {
                tmp += data.array[pos1];
                pos1++;
            }
            if (pos2 < input.data.len)
            {
                tmp += input.data.array[pos2];
                pos2++;
            }
            result.data.array[count] = tmp;
            tmp >>= INT_BIT;
            count++;
        }
    }
    else
    {
        if (abs_equal(input))
        {
            result.data.len = result.true_len();
            return result;
        }
        if (abs_larger(input))
        {
            result.data.neg = data.neg;
            while (count < result.data.len)
            {
                if (pos1 < data.len)
                {
                    tmp += data.array[pos1];
                    pos1++;
                }
                if (pos2 < input.data.len)
                {
                    tmp -= input.data.array[pos2];
                    pos2++;
                }
                result.data.array[count] = tmp;
                tmp >>= INT_BIT;
                count++;
            }
        }
        else
        {
            result.data.neg = !data.neg;
            while (count < result.data.len)
            {
                if (pos2 < input.data.len)
                {
                    tmp += data.array[pos2];
                    pos2++;
                }
                if (pos1 < data.len)
                {
                    tmp -= input.data.array[pos1];
                    pos1++;
                }
                result.data.array[count] = tmp;
                tmp >>= INT_BIT;
                count++;
            }
        }
    }
    result.data.len = result.true_len();
    return result;
}
inline HyperInt HyperInt::operator-() const
{
    //cout<<"_const-_"<<endl;
    HyperInt out(*this);
    out.data.neg = !data.neg;
    return out;
}
HyperInt HyperInt::operator*(const HyperInt& input) const
{
    h_int in1 = data, in2 = input.data;
    //HyperInt result=quick_multiply(in1,in2);
    //result.data.neg=data.neg^input.data.neg;
    return multiply(input);
}
// HyperInt HyperInt::operator*(const HyperInt &input)
// {
//     return multiply(input);
// }
HyperInt HyperInt::operator/(const HyperInt& input) const
{
    if (!input)
        throw "Can't divided by zero";
    REG HyperInt result, mid, tmp = 1;
    if (abs_smaller(input))
        return result;
    else
    {
        tmp.reset_size(data.len - (input.data.len < 3 ? 1 : input.data.len) + 1);
        result.reset_size(data.len - (input.data.len < 3 ? 1 : input.data.len) + 1);
        result.data.neg = data.neg ^ input.data.neg;
        while (!abs_smaller(tmp * input))
        {
            tmp.quick_self_twice();
        }
        result = tmp.quick_half();
        mid = (result + tmp).half();
        while (!result.abs_equal(mid))
        {
            if (abs_smaller(mid * input))
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
HyperInt HyperInt::operator%(const HyperInt& input) const
{
    HyperInt tmp = *this / input;
    tmp *= input;
    return *this - tmp;
}
HyperInt HyperInt::operator+=(const HyperInt& input)
{
    *this = *this + input;
    return *this;
}
HyperInt HyperInt::operator-=(const HyperInt& input)
{
    *this = *this - input;
    return *this;
}
HyperInt HyperInt::operator*=(const HyperInt& input)
{
    *this = *this * input;
    return *this;
}
HyperInt HyperInt::operator/=(const HyperInt& input)
{
    *this = *this / input;
    return *this;
}
HyperInt HyperInt::operator%=(const HyperInt& input)
{
    *this = *this % input;
    return *this;
}
HyperInt HyperInt::operator ++ (int)
{
    HyperInt tmp(*this);
    *this = *this + 1;
    return tmp;
}
HyperInt& HyperInt::operator ++ ()
{
    *this = *this + 1;
    return *this;
}
HyperInt HyperInt::operator -- (int)
{
    HyperInt tmp(*this);
    *this = *this - 1;
    return tmp;
}
HyperInt& HyperInt::operator -- ()
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

bool operator==(const long long& input1, const HyperInt& input2)
{
    return input2 == input1;
}
bool operator>(const long long& input1, const HyperInt& input2)
{
    return input2 > input1;
}
bool operator>=(const long long& input1, const HyperInt& input2)
{
    return input2 >= input1;
}
bool operator<(const long long& input1, const HyperInt& input2)
{
    return input2 < input1;
}
bool operator<=(const long long& input1, const HyperInt& input2)
{
    return input2 <= input1;
}
HyperInt operator+(const long long& input1, const HyperInt& input2)
{
    return input2 + input1;
}
HyperInt operator-(const long long& input1, const HyperInt& input2)
{
    return input2 - input1;
}
HyperInt operator*(const long long& input1, const HyperInt& input2)
{
    return input2 * input1;
}
HyperInt operator/(const long long& input1, const HyperInt& input2)
{
    return input2 / input1;
}
HyperInt operator%(const long long& input1, const HyperInt& input2)
{
    return input2 % input1;
}
HyperInt operator+=(long long& input1, const HyperInt& input2)
{
    return input1 += input2.to_int64();
}
HyperInt operator-=(long long& input1, const HyperInt& input2)
{
    return input1 -= input2.to_int64();
}
HyperInt operator*=(long long& input1, const HyperInt& input2)
{
    return input1 *= input2.to_int64();
}
HyperInt operator/=(long long& input1, const HyperInt& input2)
{
    return input1 /= input2.to_int64();
}
HyperInt operator%=(long long& input1, const HyperInt& input2)
{
    return input1 %= input2.to_int64();
}
HyperInt abs(HyperInt input)
{
    HyperInt result(input);
    result.data.neg = false;
    return result;
}
void print(HyperInt& input)
{
    input.console_out_dec();
}
#endif