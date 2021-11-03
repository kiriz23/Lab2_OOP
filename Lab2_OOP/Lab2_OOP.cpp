#include <iostream>
#include <vector>
#include <sstream>
#include<string>
#include <iomanip>
#include <cmath>


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
    if (left.is_negative) {
        if (right.is_negative) return -(-left + (-right));
        else return right - (-left);
    }
    else if (right.is_negative) return left - (-right);
    int carry = 0;
    for (size_t i = 0; i < std::max(left.digits.size(), right.digits.size()) || carry != 0; ++i) {
        if (i == left.digits.size()) left.digits.push_back(0);
        left.digits[i] += carry + (i < right.digits.size() ? right.digits[i] : 0);
        carry = left.digits[i] >= LongInt::BASE;
        if (carry != 0) left.digits[i] -= LongInt::BASE;
    }

    return left;
}

const LongInt operator -(LongInt left, const LongInt& right) {
    if (right.is_negative) return left + (-right);
    else if (left.is_negative) return -(-left + right);
    else if (left < right) return -(right - left);
    int carry = 0;
    for (size_t i = 0; i < right.digits.size() || carry != 0; ++i) {
        left.digits[i] -= carry + (i < right.digits.size() ? right.digits[i] : 0);
        carry = left.digits[i] < 0;
        if (carry != 0) left.digits[i] += LongInt::BASE;
    }

    left.remove_leading_zeros();
    return left;
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













LongInt Karatsuba(const LongInt& left, const LongInt& right) {
    if (left<(LongInt("10")) or right < (LongInt("10"))) return (left*right);
    
    LongInt res("0"),x0("0"),x1("0"),y0("0"),y1("0"),z0("0"),z1("0"),z2("0"),a1("0"),b1("0");
    res.digits.resize(left.digits.size() + right.digits.size());
    long long k = floor(std::min(left.digits.size(),right.digits.size()) / 2);
    x1.digits.resize(left.digits.size()-k);
    x0.digits.resize(k);
    y1.digits.resize(right.digits.size()-k);
    y0.digits.resize(k);
    for (int i = 0; i < k;i++) {
        x1.digits[i] = left.digits[i];
    }
    for (int i = k; i <left.digits.size(); i++) {
            x0.digits[i-k] = left.digits[i];
    }
    for (int i = 0; i < k; i++) {
        y1.digits[i] = right.digits[i];
    }
    for (int i = k; i < right.digits.size(); i++) {
        y0.digits[i - k] = right.digits[i];
    }
    res.digits.resize(left.digits.size() + right.digits.size());
    z0.digits.resize(x0.digits.size() + y0.digits.size());
    z1.digits.resize((x0+x1).digits.size() + (y0+y1).digits.size());
    z2.digits.resize(x1.digits.size() + y1.digits.size());
    z0 = Karatsuba(x1, y1);
    z1 = Karatsuba(x0 + x1, y0 + y1);
    z2 = Karatsuba(x0, y0);
    
    res = (z2 * LongInt(static_cast<long long>(pow(10, (k * 2))))) + ((z1 - z2 - z0) * LongInt(static_cast<long long>(pow(10, (k))))) + z0;
    
    return res;

}









int main() {
    long long b = 12334311;
    LongInt a("22333311");
    LongInt c(b);
    LongInt e("0");
    LongInt x("0");
    x = (Karatsuba(a, c));
    std::cout <<x;
    std::cout << "Kara" << std::endl;
    e = a * c;
    std::cout << e;
    
    return 0;
}
