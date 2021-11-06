#ifndef UINT_H
#define UINT_H
#include<iostream>
#include<cstring>
using namespace std;


class UltraInt
{
private:
    char *data;
    size_t len=0;
public:
    ~UltraInt()
    {
        if(len)
        {
            delete[] data;
            data=NULL;
        }
        //cout<<"_deleted_"<<endl;
    }
    UltraInt()
    {
        len=3;
        data=new char[len];
        data[0]='+';
        data[1]='0';
        data[2]='\0';
        //cout<<"_0_"<<endl;
    }
    UltraInt(const UltraInt &input)
    {
        len=input.len;
        data=new char[len];
        memcpy(data,input.data,len*sizeof(char));
        //cout<<"_1_"<<endl;
    }
    UltraInt(const long long &input)
    {
        len=2;
        if(input>0)
            for(long long n=1;n<=input;n*=10)
                len++;
        else if(input <0)
            for(long long n=-1;n>=input;n*=10)
                len++;
        else
            len++;
        data=new char[len];
        if(input>=0)
            *data='+';
        else
            *data='-';
        data[len-1]='\0';
        if(input>0)
            for(long long n=1,m=len-2;n<=input;n*=10,m--)
                data[m]=(input/n)%10+'0';
        else if(input <0)
            for(long long n=-1,m=len-2;n>=input;n*=10,m--)
                data[m]=(input/n)%10+'0';
        else
            data[1]='0';
        //cout<<"_2_"<<endl;
    }
    UltraInt(const char *input)
    {
        size_t i=1;
        len=strlen(input)+2;
        if(len>3)
        {
            if(*input=='+'||*input=='-'||*input=='0')
            {
                while(input[i]=='0')
                    i++;
                len-=i;
                data=new char[len];
                if(*input=='0')
                    *data='+';
                else
                *data=*input;
                memcpy(data+1,input+i,(len-1)*sizeof(char));
            }
            else
            {
                data=new char[len];
                *data='+';
                memcpy(data+1,input,(len-1)*sizeof(char));
            }
        }
        else
        {
            len=3;
            data=new char[len];
            data[0]='+';
            data[1]='0';
            data[2]='\0';
        }
            
        //cout<<"_3_"<<endl;
    }
    UltraInt operator=(const UltraInt &input)
    {
        delete[] data;
        len=input.len;
        data=new char[len];
        memcpy(data,input.data,len*sizeof(char));
        //cout<<"_4_"<<endl;
        return *this;
    }
    UltraInt operator=(const long long &input)
    {
        *this=UltraInt(input);
        //cout<<"_5_"<<endl;
        return *this;
    }
    UltraInt operator=(const char *input)
    {
        *this=UltraInt(input);
        //cout<<"_6_"<<endl;
        return *this;
    }
    friend ostream &operator<<(ostream &output,UltraInt &input);
    friend istream &operator>>(istream &input,UltraInt &output);
    friend UltraInt abs(UltraInt input);
    friend bool operator==(const long long &input1,const UltraInt &input2);
    friend bool operator>(const long long &input1,const UltraInt &input2);
    friend bool operator>=(const long long &input1,const UltraInt &input2);
    friend bool operator<(const long long &input1,const UltraInt &input2);
    friend bool operator<=(const long long &input1,const UltraInt &input2);
    friend UltraInt operator+(const long long &input1,const UltraInt &input2);
    friend UltraInt operator-(const long long &input1,const UltraInt &input2);
    friend UltraInt operator*(const long long &input1,const UltraInt &input2);
    friend UltraInt operator/(const long long &input1,const UltraInt &input2);
    friend UltraInt operator%(const long long &input1,const UltraInt &input2);
    friend UltraInt half(const UltraInt &input);
    bool operator==(const UltraInt &input)
    {
        if(len!=input.len)
            return false;
        else
            for(size_t i=0;i<len-1;i++)
                if(data[i]!=input.data[i])
                    return false;
        return true;
    }
    bool operator==(const long long &input)
    {
        return *this==UltraInt(input);
    }
    bool operator>(const UltraInt &input)
    {
        if(*data!=*input.data)
        {
            if(*data=='+')
                return true;
            else
                return false;
        }
        else if(*data=='+')
        {
            if(len>input.len)
                return true;
            else if(len<input.len)
                return false;
            else
                for(size_t i=0;i<len-1;i++)
                {
                    if(data[i]>input.data[i])
                        return true;
                    else if(data[i]<input.data[i])
                        return false;
                }
        }
        else
        {
            if(len<input.len)
                return true;
            else if(len>input.len)
                return false;
            else
                for(size_t i=0;i<len-1;i++)
                {
                    if(data[i]<input.data[i])
                        return true;
                    else if(data[i]>input.data[i])
                        return false;
                }
        }
        return false;
    }
    bool operator>(const long long &input)
    {
        return *this>UltraInt(input);
    }
    bool operator>=(const UltraInt &input)
    {
        if(*this<input)
            return false;
        else
            return true;
    }
    bool operator>=(const long long &input)
    {
        if(*this<UltraInt(input))
            return false;
        else
            return true;
    }
    bool operator<(const UltraInt &input)
    {
        if(*data!=*input.data)
        {
            if(*data=='-')
                return true;
            else
                return false;
        }
        else if(*data=='-')
        {
            if(len>input.len)
                return true;
            else if(len<input.len)
                return false;
            else
                for(size_t i=0;i<len-1;i++)
                {
                    if(data[i]>input.data[i])
                        return true;
                    else if(data[i]<input.data[i])
                        return false;
                }
        }
        else
        {
            if(len<input.len)
                return true;
            else if(len>input.len)
                return false;
            else
                for(size_t i=0;i<len-1;i++)
                {
                    if(data[i]<input.data[i])
                        return true;
                    else if(data[i]>input.data[i])
                        return false;
                }
        }
        return false;
    }
    bool operator<(const long long &input)
    {
        return *this<UltraInt(input);
    }
    bool operator<=(const UltraInt &input)
    {
        if(*this>input)
            return false;
        else
            return true;
    }
    bool operator<=(const long long &input)
    {
        if(*this>UltraInt(input))
            return false;
        else
            return true;
    }
    UltraInt operator+(const UltraInt &input)
    {
        size_t len2,out_len=1+(len>input.len?len:input.len);
        char *out;
        short tmp;
        out=new char[out_len];
        if(*data==*input.data)
        {
            memset(out+1,'0',(out_len-len)*sizeof(char));//初始化
            memcpy(out+out_len-len+1,data+1,(len-1)*sizeof(char));
            *out=*data;
            len2=input.len;
            while(out_len>2)//逐位计算
            {
                tmp=out[out_len-2]-'0';
                if(len2>2)
                {
                    tmp+=input.data[len2-2]-'0';
                    len2--;
                }
                if(tmp<10)
                {
                    out[out_len-2]=tmp+'0';
                    if(len2<=2)
                        break;
                }
                else
                {
                    out[out_len-2]=tmp-10+'0';
                    out[out_len-3]++;
                }
                out_len--;
            }
        }
        else
        {
            if(abs(*this)==abs(input))
                return (long long)0;
            else if(abs(*this)>abs(input))
            {
                memset(out+1,'0',(out_len-len)*sizeof(char));//初始化
                memcpy(out+out_len-len+1,data+1,(len-1)*sizeof(char));
                *out=*data;
                len2=input.len;
                while(out_len>2)//逐位计算
                {
                    tmp=out[out_len-2]-'0';
                    if(len2>2)
                    {
                        tmp-=(input.data[len2-2]-'0');
                        len2--;
                    }
                    if(tmp>=0)
                    {
                        out[out_len-2]=tmp+'0';
                        if(len2<=2)
                            break;
                    }
                    else
                    {
                        out[out_len-2]=tmp+10+'0';
                        out[out_len-3]--;
                    }
                    out_len--;
                }
            }
            else
            {
                memset(out+1,'0',(out_len-input.len)*sizeof(char));//初始化
                memcpy(out+out_len-input.len+1,input.data+1,(input.len-1)*sizeof(char));
                *out=*input.data;
                len2=len;
                while(out_len>2)//逐位计算
                {
                    tmp=out[out_len-2]-'0';
                    if(len2>2)
                    {
                        tmp-=(data[len2-2]-'0');
                        len2--;
                    }
                    if(tmp>=0)
                    {
                        out[out_len-2]=tmp+'0';
                        if(len2<=2)
                            break;
                    }
                    else
                    {
                        out[out_len-2]=tmp+10+'0';
                        out[out_len-3]--;
                    }
                    out_len--;
                }
            }
        }
        return out;
    }
    UltraInt operator+(const long long &input)
    {
        return *this+UltraInt(input);
    }
    UltraInt operator-(const UltraInt &input)
    {
        size_t len2,out_len=1+(len>input.len?len:input.len);
        char *out;
        short tmp;
        out=new char[out_len];
        if(*data!=*input.data)
        {
            memset(out+1,'0',(out_len-len)*sizeof(char));//初始化
            memcpy(out+out_len-len+1,data+1,(len-1)*sizeof(char));
            *out=*data;
            len2=input.len;
            while(out_len>2)//逐位计算
            {
                tmp=out[out_len-2]-'0';
                if(len2>2)
                {
                    tmp+=input.data[len2-2]-'0';
                    len2--;
                }
                if(tmp<10)
                {
                    out[out_len-2]=tmp+'0';
                    if(len2<=2)
                        break;
                }
                else
                {
                    out[out_len-2]=tmp-10+'0';
                    out[out_len-3]++;
                }
                out_len--;
            }
        }
        else
        {
            if(abs(*this)==abs(input))
                return (long long)0;
            else if(abs(*this)>abs(input))
            {
                memset(out+1,'0',(out_len-len)*sizeof(char));//初始化
                memcpy(out+out_len-len+1,data+1,(len-1)*sizeof(char));
                *out=*data;
                len2=input.len;
                while(out_len>2)//逐位计算
                {
                    tmp=out[out_len-2]-'0';
                    if(len2>2)
                    {
                        tmp-=(input.data[len2-2]-'0');
                        len2--;
                    }
                    if(tmp>=0)
                    {
                        out[out_len-2]=tmp+'0';
                        if(len2<=2)
                            break;
                    }
                    else
                    {
                        out[out_len-2]=tmp+10+'0';
                        out[out_len-3]--;
                    }
                    out_len--;
                }
            }
            else
            {
                memset(out+1,'0',(out_len-input.len)*sizeof(char));//初始化
                memcpy(out+out_len-input.len+1,input.data+1,(input.len-1)*sizeof(char));
                *out='X'-*data;
                len2=len;
                while(out_len>2)//逐位计算
                {
                    tmp=out[out_len-2]-'0';
                    if(len2>2)
                    {
                        tmp-=(data[len2-2]-'0');
                        len2--;
                    }
                    if(tmp>=0)
                    {
                        out[out_len-2]=tmp+'0';
                        if(len2<=2)
                            break;
                    }
                    else
                    {
                        out[out_len-2]=tmp+10+'0';
                        out[out_len-3]--;
                    }
                    out_len--;
                }
            }
        }
        return out;
    }
    UltraInt operator-(const long long &input)
    {
        return *this-UltraInt(input);
    }
    UltraInt operator*(const UltraInt &input)
    {
        unsigned char tmp,sum;
        size_t out_len=len+input.len-2;
        char *out;
        out=new char[out_len];
        memset(out,'0',(out_len-1)*sizeof(char));//初始化
        out[out_len-1]='\0';
        if(*data==*input.data)
            *out='+';
        else
            *out='-';
        for(size_t i=len-2;i>0;i--)//逐位计算
            for(size_t j=input.len-2;j>0;j--)
            {
                tmp=(data[i]-'0')*(input.data[j]-'0');
                for(size_t k=i+j;k>0;k--)
                {
                    sum=tmp+out[k]-'0';
                    if(sum<10)
                    {
                        out[k]=sum+'0';
                        break;
                    }
                    else
                    {
                        out[k]=sum%10+'0';
                        tmp=sum/10;
                    }
                }
            }
        //cout<<"_*_"<<endl;
        return out;
    }
    UltraInt operator*(const long long &input)
    {
        return *this*UltraInt(input);
    }    
    UltraInt operator/(const UltraInt &input)
    {
        UltraInt in1=abs(*this),in2=abs(input),tmp="0",mid,out=1;
        if(input.len==3&&input.data[1]=='0')
            return tmp;
        if(in1<in2)
            return tmp;
        else
        {
            while(in2*out*2<=in1)
                out=out*2;
            tmp=out*2;
            mid=half(tmp+out);
            while(!(mid==out))
            {
                if(mid*in2>in1)
                    tmp=mid;
                else
                    out=mid;
                mid=half(tmp+out);
            }
        }
        if(*data!=*input.data)
            *out.data='-';
        return out;
    }
    UltraInt operator/(const long long &input)
    {
        return *this/UltraInt(input);
    }
    UltraInt operator%(const UltraInt &input)
    {
        if(abs(*this)<abs(input))
            return *this;
        return *this-((*this/input)*input);
    }
    UltraInt operator%(const long long &input)
    {
        if(abs(*this)<abs(input))
            return *this;
        return *this-((*this/UltraInt(input))*input);
    }
    UltraInt operator ++ (int)
    {
        UltraInt tmp=*this;
        *this=*this+1;
        return tmp;
    }
    UltraInt &operator ++ ()
    {
        *this=*this+1;
        return *this;
    }
    UltraInt operator -- (int)
    {
        UltraInt tmp=*this;
        *this=*this-1;
        return tmp;
    }
    UltraInt &operator -- ()
    {
        *this=*this-1;
        return *this;
    }
    void print()
    {
        cout << data << endl;
    }
};
ostream &operator<<(ostream &output,UltraInt &input)
{
    output<<input.data;
    return output;
}
istream &operator>>(istream &input,UltraInt &output)
{
    size_t i=1;
    char tmp[2000000];
    input>>tmp;
    output.len=strlen(tmp)+2;
    if(*tmp=='+'||*tmp=='-'||*tmp=='0')
    {
        while(tmp[i]=='0')
            i++;
        output.len-=i;
        output.data=new char[output.len];
        if(*tmp=='0')
                *output.data='+';
        else
            *output.data=*tmp;
        memcpy(output.data+1,tmp+i,(output.len-1)*sizeof(char));
    }
    else
    {
        output.data=new char[output.len];
        *output.data='+';
        memcpy(output.data+1,tmp,(output.len-1)*sizeof(char));
    }
    return input;
}
UltraInt abs(UltraInt input)
{
    *input.data='+';
    return input;
}
bool operator==(const long long &input1,const UltraInt &input2)
{
    return UltraInt(input1)==input2;
}
bool operator>(const long long &input1,const UltraInt &input2)
{
    return UltraInt(input1)>input2;
}
bool operator>=(const long long &input1,const UltraInt &input2)
{
    return UltraInt(input1)>=input2;
}
bool operator<(const long long &input1,const UltraInt &input2)
{
    return UltraInt(input1)<input2;
}
bool operator<=(const long long &input1,const UltraInt &input2)
{
    return UltraInt(input1)<=input2;
}
UltraInt operator+(const long long &input1,UltraInt &input2)
{
    return input2+input1;
}
UltraInt operator-(const long long &input1,UltraInt &input2)
{
    return input2-input1;
}
UltraInt operator*(const long long &input1,UltraInt &input2)
{
    return input2*input1;
}
UltraInt operator/(const long long &input1,UltraInt &input2)
{
    return UltraInt(input1)/input2;
}
UltraInt operator%(const long long &input1,UltraInt &input2)
{
    return UltraInt(input1)%input2;
}
UltraInt half(const UltraInt &input)
{
    size_t out_len=input.len;
    char *out=new char[out_len];
    short tmp=0;
    memset(out+1,'0',(out_len-2)*sizeof(char));
    *out=*input.data;
    out[out_len-1]='\0';
    for(size_t i=1;i<out_len-1;i++)
    {
        tmp=(tmp+input.data[i]-'0')/2;
        out[i]=tmp+'0';
        tmp=((input.data[i]-'0')%2)*10;
    }
    return out;
}
#endif
