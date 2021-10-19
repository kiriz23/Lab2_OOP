#include <iostream>
#include <vector>


class Multiplication;


class LongInt {
public:
    std::string digits;
    
    
    LongInt add(LongInt other) {
        int temp = 0;
        bool next_order = false;
        LongInt res("00");
        int k = std::max(digits.length(), other.digits.length());
        for (int i = 0; i < k; i++) {
            temp = other.digits[i] + digits[i] ;

            if (temp >= 10) {
                res.digits[i] = temp % 10 + '0';
            }
            else {
                res.digits[i] = temp + '0';
                temp = 0;
            }
            if (i == 0 && temp / 10 != 0)
                next_order = true;
            temp /= 10;
        }
        if (next_order)
            res.digits.insert(0, 1, '1');
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

    Multiplication* multiplication;
   
    LongInt operator*(LongInt& other) { return multiplication->multiply(*this, other); };


};


class Multiplication {
public:
    virtual LongInt multiply(LongInt a, LongInt b) = 0;
};


class KaratsubaMultiplication : public Multiplication {
public:
    LongInt multiply(LongInt a, LongInt b);
};


class NaiveMultiplication : public Multiplication {
public:
    LongInt multiply(LongInt a, LongInt b) {};
};


LongInt LongInt::operator*(LongInt& other) {
    return multiplication->multiply(*this, other);
}


LongInt KaratsubaMultiplication::multiply(LongInt a, LongInt b) {
    LongInt res;
    std::cout << "Karatsuba ok";
    return res;
}


LongInt NaiveMultiplication::multiply(LongInt a, LongInt b) {
    LongInt res;
    std::cout << "Naive ok";
    return res;
}












int main() {
    LongInt a, b, c, o;
    a.multiplication = new KaratsubaMultiplication();
    c = a * b;
    /* c = a.add(b);
    o = a + b;
    o.print();

    c.print();*/
    a.multiplication = new NaiveMultiplication();
    return 0;
}