#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double PI = 3.141592653589793238460;

class complex{
    public:

    double real;
    double imaginary;

    complex(){
        real = 0;
        imaginary = 0;
    }

    complex(double real, double imaginary){
        this->real = real;
        this->imaginary = imaginary;
    }

    complex operator+(const complex& b){
        return complex(real + b.real, imaginary + b.imaginary);
    }

    complex operator-(const complex& b){
        return complex(real - b.real, imaginary - b.imaginary);
    }

    complex operator*(const complex& b){
        return complex(real*b.real - b.imaginary*imaginary, real*b.imaginary + b.real*imaginary);
    }
};

//apply Cooley-Tukey Recursively
void cooleyTukey(vector<complex>& values){
    int n = values.size();
    
    if (n <= 1){ //Base case -> no more edits to the array
        return;
    }

    vector<complex> even(n/2), odd(n/2); //this assumes that n is always a power of 2
    for (int i = 0 ; i < n/2; i ++){
        even[i] = values[2*i];
        odd[i] = values[2*i+1];
    }

    cooleyTukey(even); //In place recurse on the even and odd
    cooleyTukey(odd);
    
    for (int i = 0 ; i < n/2; i ++){ //combine both even and odd
        //twidle
        double theta =  -2 * PI * i / n;
        complex t = complex(cos(theta), sin(theta)) * odd[i];
        values[i] = even[i] + t;
        values[i + n/2] = even[i] - t;
    }
}

int main(){
    // Example input: 8 complex numbers
    vector<complex> input = {
        {0, 0}, {1, 0}, {2, 0}, {3, 0},
        {4, 0}, {5, 0}, {6, 0}, {7, 0}
    };

    cooleyTukey(input);

    cout << "FFT Result: \n";
    for (const auto& x : input) { //print out all the results
        cout << "Real: " << x.real << ", Imaginary: " << x.imaginary << endl;
    }

    return 0;
}