#include <iostream>
#include <vector>
#include <sstream>
#include<string>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <complex>

const double PI = acos(-1);
typedef std::complex<double> base;





class LongInt {
public:
    class divide_by_zero : public std::exception {  };
    
    
    std::vector<int>digits;
    bool is_negative;
    static const int BASE = 10;
    
    void remove_leading_zeros() {
        while (this->digits.size() > 1 && this->digits.back() == 0) {
            this->digits.pop_back();
        }
        if (this->digits.size() == 1 && this->digits[0] == 0) this->is_negative = false;
    }
    
    
    LongInt(std::string str) {
        if (str.length() == 0) {
            this->is_negative = false;
        }
        else {
            if (str[0] == '-') {
                str = str.substr(1);
                this->is_negative = true;
            }
            else {
                this->is_negative = false;
            }
            for (long long i = str.length(); i > 0; i -= 1) {
                if (i < 1)
                    this->digits.push_back(atoi(str.substr(0, i).c_str()));
                else
                    this->digits.push_back(atoi(str.substr(i - 1, 1).c_str()));
            }
            this->remove_leading_zeros();
        }
    }


    LongInt(signed long long l) {
        if (l < 0) { this->is_negative = true; l = -l; }
        else this->is_negative = false;
        do {
            this->digits.push_back(l % BASE);
            l /= BASE;
        } while (l != 0);
    }

    LongInt (unsigned long long l) {
        this->is_negative = false;
        do {
            this->digits.push_back(l % BASE);
            l /= BASE;
        } while (l != 0);
    }


    LongInt(const LongInt& l) {
        this->is_negative = l.is_negative;
        long long k = l.digits.size();
        this->digits.resize(k);
        for (int i = 0; i < k; i++) {
            this->digits[i] = l.digits[i];
        }
    }


    const LongInt operator +() const {
        return LongInt(*this);
    }

    const LongInt operator -() const {
        LongInt copy(*this);
        copy.is_negative = !copy.is_negative;
        return copy;
    }

    



    
    

   
    LongInt& operator+=(const LongInt& value);

    LongInt& operator-=(const LongInt& value);

    const LongInt operator++();

    const LongInt operator++(int);

    const LongInt operator--();

    const LongInt operator--(int);

    LongInt& operator*=(const LongInt& value);

    void shift_right();

    LongInt& operator/=(const LongInt& value);

    bool odd() const;

    bool even() const;

    const LongInt pow(LongInt n) const;

    void normalize();

    int get_base();

    void set_base(int base);

    LongInt Karatsuba(const LongInt& left, const LongInt& right);
    
    

   

};


std::ostream& operator <<(std::ostream& os, const LongInt& bi) {
    if (bi.digits.empty()) os << 0;
    else {
        if (bi.is_negative) os << '-';
        os << bi.digits.back();
        char old_fill = os.fill('0');
        for (long long i = static_cast<long long>(bi.digits.size()) - 2; i >= 0; --i) {
            os << std::setw(1) << bi.digits[i];
        }

        os.fill(old_fill);
    }

    return os;
}






bool operator ==(const LongInt& left, const LongInt& right) {
    if (left.is_negative != right.is_negative) return false;
    if (left.digits.empty()) {
        if (right.digits.empty() || (right.digits.size() == 1 && right.digits[0] == 0)) return true;
        else return false;
    }

    if (right.digits.empty()) {
        if (left.digits.size() == 1 && left.digits[0] == 0) return true;
        else return false;
    }
    if (left.digits.size() != right.digits.size()) return false;
    for (size_t i = 0; i < left.digits.size(); ++i) if (left.digits[i] != right.digits[i]) return false;

    return true;
}


bool operator <(const LongInt& left, const LongInt& right) {
    if (left == right) return false;
    if (left.is_negative) {
        if (right.is_negative) return ((-right) < (-left));
        else return true;
    }
    else if (right.is_negative) return false;
    else {
        if (left.digits.size() != right.digits.size()) {
            return left.digits.size() < right.digits.size();
        }
        else {
            for (long long i = left.digits.size() - 1; i >= 0; --i) {
                if (left.digits[i] != right.digits[i]) return left.digits[i] < right.digits[i];
            }

            return false;
        }
    }
}


bool operator !=(const LongInt& left, const LongInt& right) {
    return !(left == right);
}

bool operator <=(const LongInt& left, const LongInt& right) {
    return (left < right || left == right);
}

bool operator >(const LongInt& left, const LongInt& right) {
    return !(left <= right);
}

bool operator >=(const LongInt& left, const LongInt& right) {
    return !(left < right);
}


const LongInt operator -(LongInt left, const LongInt& right);


const LongInt operator +(LongInt left, const LongInt& right) {
    LongInt res = LongInt(left);
    if (res.is_negative) {
        if (right.is_negative) return -(-res + (-right));
        else return right - (-res);
    }
    else if (right.is_negative) return res - (-right);
    int carry = 0;
    for (size_t i = 0; i < std::max(res.digits.size(), right.digits.size()) || carry != 0; ++i) {
        if (i == res.digits.size()) res.digits.push_back(0);
        res.digits[i] += carry + (i < right.digits.size() ? right.digits[i] : 0);
        carry = res.digits[i] >= LongInt::BASE;
        if (carry != 0) res.digits[i] -= LongInt::BASE;
    }
    return res;
}

const LongInt operator -(LongInt left, const LongInt& right) {
    LongInt res = LongInt(left);
    if (right.is_negative) return res + (-right);
    else if (res.is_negative) return -(-res + right);
    else if (res < right) return -(right - res);
    int carry = 0;
    for (size_t i = 0; i < right.digits.size() || carry != 0; ++i) {
        res.digits[i] -= carry + (i < right.digits.size() ? right.digits[i] : 0);
        carry = res.digits[i] < 0;
        if (carry != 0) res.digits[i] += LongInt::BASE;
    }
    res.remove_leading_zeros();
    return res;
}


LongInt& LongInt::operator +=(const LongInt& value) {
    return *this = (*this + value);
}

LongInt& LongInt::operator -=(const LongInt& value) {
    return *this = (*this - value);
}


const LongInt LongInt::operator++() {
    return (*this += LongInt(unsigned long long(1)));
}

const LongInt LongInt::operator ++(int) {
    *this += LongInt(long long(1));
    return *this - LongInt(unsigned long long(1));
}

const LongInt LongInt::operator --() {
    return *this -= LongInt(unsigned long long(1));
}

const LongInt LongInt::operator --(int) {
    *this -= LongInt(unsigned long long(1));
    return *this + LongInt(unsigned long long(1));
}


const LongInt operator *(const LongInt& left, const LongInt& right) {
    LongInt result("0");
    result.digits.resize(left.digits.size() + right.digits.size());
    for (size_t i = 0; i < left.digits.size(); ++i) {
        int carry = 0;
        for (size_t j = 0; j < right.digits.size() || carry != 0; ++j) {
            long long cur = result.digits[i + j] +
                left.digits[i] * 1LL * (j < right.digits.size() ? right.digits[j] : 0) + carry;
            result.digits[i + j] = static_cast<int>(cur % LongInt::BASE);
            carry = static_cast<int>(cur / LongInt::BASE);
        }
    }
    result.is_negative = left.is_negative != right.is_negative;
    result.remove_leading_zeros();
    return result;
}


LongInt& LongInt::operator *=(const LongInt& value) {
    return *this = (*this * value);
}

void LongInt::shift_right() {
    if (this->digits.size() == 0) {
        this->digits.push_back(0);
        return;
    }
    this->digits.push_back(this->digits[this->digits.size() - 1]);
    for (size_t i = this->digits.size() - 2; i > 0; --i) this->digits[i] = this->digits[i - 1];
    this->digits[0] = 0;
}


const LongInt operator /(const LongInt& left, const LongInt& right) {
    if (right == LongInt(0LL)) throw LongInt::divide_by_zero();
    LongInt b = right;
    b.is_negative = false;
    LongInt result("0"), current("0");
    result.digits.resize(left.digits.size());
    for (long long i = static_cast<long long>(left.digits.size()) - 1; i >= 0; --i) {
        current.shift_right();
        current.digits[0] = left.digits[i];
        current.remove_leading_zeros();
        int x = 0, l = 0, r = LongInt::BASE;
        while (l <= r) {
            int m = (l + r) / 2;
            LongInt t("0");
            t = b * LongInt(long long(m));
            if (t <= current) {
                x = m;
                l = m + 1;
            }
            else r = m - 1;
        }

        result.digits[i] = x;
        current = current - b * LongInt(long long(x));
    }

    result.is_negative = left.is_negative != right.is_negative;
    result.remove_leading_zeros();
    return result;
}


LongInt& LongInt::operator /=(const LongInt& value) {
    return *this = (*this / value);
}


const LongInt operator %(const LongInt& left, const LongInt& right) {
    LongInt result = left - (left / right) * right;
    if (result.is_negative) result += right;
    return result;
}


bool LongInt::odd() const {
    if (this->digits.size() == 0) return false;
    return this->digits[0] & 1;
}

bool LongInt::even() const {
    return !this->odd();
}


const LongInt LongInt::pow(LongInt n) const {
    LongInt a(*this), result(LongInt(1LL));
    while (n != LongInt(0LL)) {
        if (n.odd()) result *= a;
        a *= a;
        n /= LongInt(2LL);
    }

    return result;
}


void LongInt::normalize() {
    while (digits.size() > 1 && digits.back() == 0)
        digits.pop_back();

    if (digits.empty())
        digits.push_back(0);
}


void split_at_begin (LongInt& b, int m, LongInt& b1, LongInt& b2) {
    b2.digits.assign(b.digits.end() - m, b.digits.end());
    b1.digits.assign(b.digits.begin(), b.digits.end() - m);
}


void split_at_end(LongInt& b, int m, LongInt& b1, LongInt& b2) {
    b2.digits.assign(b.digits.begin(), b.digits.begin() + m);
    b1.digits.assign(b.digits.begin() + m, b.digits.end());
}


LongInt& long_pow_10(LongInt& b, int deg){
    b.normalize();

    int n = b.digits.size();


    for (int i = n; i < n + deg; ++i)
        b.digits.push_back(0);


    for (int i = n + deg - 1; i >= deg; --i) {
        b.digits[i] = b.digits[i - deg];
    }

    for (int i = 0; i < deg; ++i) {
        b.digits[i] = 0;
    }

    return b;
}





LongInt Karatsuba(LongInt& b1, LongInt& b2) {
    if (b1.digits.size() == 1 || b2.digits.size() == 1) {
        return b1*b2;
    }

    int m = std::min(b1.digits.size(), b2.digits.size()) / 2;
    int floor_m = (int)m;
    int ceil_m = (int)(m + 0.5);

    LongInt high1("0"), low1("0");
    LongInt high2("0"), low2("0");

    split_at_end(b1, m, high1, low1);
    split_at_end(b2, m, high2, low2);

    LongInt sum1 = low1 + high1;
    LongInt sum2 = low2 + high2;

    LongInt z0 = Karatsuba(low1, low2);
    LongInt z1 = Karatsuba(sum1, sum2);
    LongInt z2 = Karatsuba(high1, high2);

    LongInt temp1 = z1 - z2;
    LongInt temp2 = temp1 - z0;
    LongInt s1 = LongInt(long_pow_10(z2, m * 2));
    LongInt s2 = LongInt(long_pow_10(temp2, m));
    LongInt s3 = LongInt(s1 + s2);
    LongInt s4 = LongInt(s3 + z0);

    return s4;
}


LongInt Toom_Cook(const LongInt& left, const LongInt& right) {
if (right.digits.size() < left.digits.size()) return Toom_Cook(right, left);
LongInt res("0"), x0("0"), x1("0"), x2("0"), y0("0"), y1("0"), y2("0"), y2byx2("0"), y1byx1("0"), y0byx0("0"), z10("0"), z20("0"), z21("0"), y2_y1("0"), y2_y0("0");
LongInt y1_y0("0"), x1_x0("0"), x2_x0("0"), x2_x1("0"), res_pow_m1("0"), res_pow_m2("0"), res_pow_0("0"), res_pow_2m1("0"), res_pow_2m2("0"), res_pow_m2m1("0");
long long size = std::min(left.digits.size(), right.digits.size());
long long actual_size = size;
if (size < 4) { return left * right; }
size = size - size % 3;
x0.digits.resize(size / 3);
y0.digits.resize(size / 3);
x1.digits.resize(size / 3);
y1.digits.resize(size / 3);
x2.digits.resize(size / 3);
y2.digits.resize(right.digits.size()- actual_size);
for (long long  i = 0; i < size / 3; i++)
{
    x0.digits.push_back(left.digits[i]);
    y0.digits.push_back(right.digits[i]);
}
for (long long  i = size / 3; i < 2 * size / 3; i++)
{
    x1.digits.push_back(left.digits[i]);
    y1.digits.push_back(right.digits[i]);
}
for (long long  i = 2 * size / 3; i < actual_size; i++)
{
    x2.digits.push_back(left.digits[i]);
    y2.digits.push_back(right.digits[i]);
}
for (long long  i = actual_size; i < right.digits.size(); i++)
{
    y2.digits.push_back(right.digits[i]);
}
y2byx2 = LongInt(Toom_Cook(x2,y2));
y1byx1 = LongInt(Toom_Cook(x1, y1));
y0byx0 = LongInt(Toom_Cook(x0, y0));

x2_x1 = LongInt(x2+ x1);
x2_x0 = LongInt(x2+ x0);
x1_x0 = LongInt(x1+ x0);
y2_y1 = LongInt(y2+y1);
y2_y0 = LongInt(y2+ y0);
y1_y0 = LongInt(y1+y0);

z10 = LongInt(Toom_Cook(y1_y0, x1_x0));
z10 = LongInt(z10-y1byx1);
z10 = LongInt(z10-y0byx0);

z20 = LongInt(Toom_Cook(y2_y0, x2_x0));
z20 = LongInt(z20-y2byx2);
z20 = LongInt(z20- y0byx0);

z21 = LongInt(Toom_Cook(y2_y1, x2_x1));
z21 = LongInt(z21- y2byx2);
z21 = LongInt(z21-y1byx1);
res_pow_m1.digits.resize(size / 3);
res_pow_2m1.digits.resize(2 * size / 3);
res_pow_m2.digits.resize(2 * size / 3);
res_pow_2m2.digits.resize(4 * size / 3);
res_pow_m2m1.digits.resize(3 * size / 3);



for (long long  i = 0; i < size / 3; i++)
{
    res_pow_m1.digits.push_back(0);
}
for (long long  i = 0; i < 2 * size / 3; i++)
{
    res_pow_2m1.digits.push_back(0);
}
for (long long  i = 0; i < 2 * size / 3; i++)
{
    res_pow_m2.digits.push_back(0);
}
for (long long i = 0; i < 4 * size / 3; i++)
{
    res_pow_2m2.digits.push_back(0);
}
for (long long i = 0; i < 3 * size / 3; i++)
{
    res_pow_m2m1.digits.push_back(0);
}
res_pow_m1.digits.insert(end(res_pow_m1.digits), begin(z10.digits), end(z10.digits));
res_pow_2m1.digits.insert(end(res_pow_2m1.digits), begin(y1byx1.digits), end(y1byx1.digits));
res_pow_m2.digits.insert(end(res_pow_m2.digits), begin(z20.digits), end(z20.digits));
res_pow_2m2.digits.insert(end(res_pow_2m2.digits), begin(y2byx2.digits), end(y2byx2.digits));
res_pow_m2m1.digits.insert(end(res_pow_m2m1.digits), begin(z21.digits), end(z21.digits));
res_pow_0 = LongInt(y0byx0);
res = LongInt(res_pow_2m2 + res_pow_m2);
res = LongInt(res + res_pow_2m1);
res = LongInt(res + res_pow_m1);
res = LongInt(res + res_pow_m2m1);
res = LongInt(res + res_pow_0);

return res;

}


int LongInt::get_base() {
    return LongInt::BASE;
}












void fft(std::vector<base>& a, bool invert) {
    int n = (int)a.size();
    if (n == 1) return;

    std::vector<base> a0(n / 2), a1(n / 2);
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = a[i];
        a1[j] = a[i + 1];
    }
    fft(a0, invert);
    fft(a1, invert);

    double ang = 2 * PI / n * (invert ? -1 : 1);
    base w(1), wn(cos(ang), sin(ang));
    for (int i = 0; i < n / 2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (invert)
            a[i] /= 2, a[i + n / 2] /= 2;
        w *= wn;
    }
}


LongInt Shonhage(LongInt& b1, LongInt& b2) {
    LongInt res("0");
    

    std::vector<base> f1(b1.digits.begin(), b1.digits.end()), f2(b2.digits.begin(), b2.digits.end());
    size_t n = 1;

    while (n < std::max(f1.size(), f2.size())) n <<= 1;
    n <<= 1;
    f1.resize(n), f2.resize(n);

    fft(f1, false), fft(f2, false);
    for (size_t i = 0; i < n; ++i)
        f1[i] *= f2[i];
    fft(f1, true);

    int carry = 0;
    int temp;
    for (size_t i = 0; i < n; ++i) {
        temp = int(f1[i].real() + 0.5) + carry;
        carry = temp / res.get_base();
        temp %= res.get_base();

        res.digits.push_back(temp);
    }

    res.digits.push_back(carry);
    res.normalize();
    res.digits.erase(res.digits.begin());
    return res;
}





















int main() {
  /*  std::cout << "Enter two long numbers"<<std::endl;
    std::string str1, str2;
    std::cin >> str1 >> str2;
    LongInt first(str1), second(str2);
    std::cout << "Enter multiplication type: " << std::endl;
    std::cout << "\t1 - Naive Multiplication" << std::endl;
    std::cout << "\t2 - Karatsuba Multiplication" << std::endl;
    int choise;
    std::cin >> choise;
    while (choise != 0) {
        switch (choise) {
        case 1:
            std::cout << first * second;
            break;
        case 2:
            std::cout << Karatsuba(first, second);
            break;
        default:
            std::cout << "Wrong input!!!\nOnly 1-2 requires. Type 0 to exit or try again";
        }
        std::cout << std::endl;
        std::cin >> choise;
    }
    return 0;
    
    */
    
    
    std::string b = "527572472462465";
    LongInt a("12536646456546456456456546");
    LongInt c(b);
    LongInt e("0");
    
   
    e =( a * c);
    
    std::cout << e<<std::endl;
    
    LongInt f("0");
    f = Shonhage(a, c);
    std::cout << f <<"Tom" << std::endl;
    return 0;
}
