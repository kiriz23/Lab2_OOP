#include <iostream>


class LongInt {
public:
    std::string digits;

    LongInt(std::string initialDigits = "") :digits(initialDigits) {
    }

    LongInt add(LongInt other) {
        LongInt res;
        int ost = 0;
        int k = std::max(digits.length(), other.digits.length());
        res.digits.insert(res.digits.begin(), k + 1, '0');

        for (int i = k; i >= 0; i--) {
            int t = 0;
            if (i < other.digits.length()) {
                t += other.digits[i] - '0';
            }
            if (i < digits.length()) {
                t += digits[i] - '0';
            }
            t += ost;
            if (t >= 10) {
                ost = 1;
                t -= 10;
            }
            else ost = 0;

            res.digits[i + 1] = t + '0';
        }
        if (ost != 0) {
            res.digits[0] += ost;
        }
        else {
            res.digits.erase(0, 1);
        }
        return res;
    }

    void print() {
        for (int i = 0; i < digits.length(); ++i) {
            std::cout << digits[i];
        }
        std::cout << std::endl;
    }


    LongInt operator+(LongInt other) {
        return add(other);
    }

};






int main() {
    LongInt a("999"), b("210"), c, o;

    c = a.add(b);
    o = a + b;
    o.print();

    c.print();

    return 0;
}