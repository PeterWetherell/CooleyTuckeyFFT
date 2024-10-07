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
        //cout << odd[i].real << " " << even[i].real << " " << cos(theta) << ", ";
        complex t = complex(cos(theta), sin(theta)) * odd[i];
        values[i] = even[i] + t;
        values[i + n/2] = even[i] - t;
    }
}

//apply Cooley-Tukey Recursively
void CooleyTukeyInPlaceRec(int start, int indexing, vector<complex>& values, vector<complex>& twiddle, vector<complex>& savedTwiddle){//convert the algorithm to work in place
    if (values.size() == indexing){ //Base case -> no more edits to the array
        return;
    }
    
    //In place recurse on the even and odd
    CooleyTukeyInPlaceRec(start,indexing*2, values, twiddle, savedTwiddle); //Even
    CooleyTukeyInPlaceRec(start + indexing,indexing*2, values, twiddle, savedTwiddle); //Odd

    int halfN = values.size()/(2*indexing);
    
    for (int i = 0 ; i < halfN; i ++){ //combine both even and odd
        //cout << values[start+indexing + i*2*indexing].real << " " << values[start + i*2*indexing].real << " " << twiddle[i*indexing].real << ", ";

        complex t = twiddle[i*indexing] * values[start+indexing + i*2*indexing];

        //Do the first half
        values[start + indexing*i] = values[start + i*2*indexing] + t; //even + t

        //Save the twiddle factor
        savedTwiddle[i] = complex(-2.0,0)*t;
    }
    for (int i = 0 ; i < halfN; i ++){ //combine both even and odd
        values[start + indexing*(i+halfN)] = values[start + indexing*i] + savedTwiddle[i]; //even - t = (even+t) - 2*t
    }
}

//Setup for the inplace
void cooleyTukeyInPlace(vector<complex>& values){ 
    int n = values.size();
    vector<complex> twiddle(n/2);
    for (int i = 0 ; i < n/2; i ++){ //calculate all the twidles.
        double theta =  -2 * PI * i / n;
        twiddle[i] = complex(cos(theta), sin(theta));
    }
    vector<complex> savedTwiddle(n/2);
    CooleyTukeyInPlaceRec(0, 1, values, twiddle, savedTwiddle);
}

int main(){
    // Example input: 8 complex numbers
    vector<complex> input = {
        {0, 0}, {1, 0}, {2, 0}, {3, 0},
        {4, 0}, {5, 0}, {6, 0}, {7, 0}
    };

    //cooleyTukey(input);
    cooleyTukeyInPlace(input);

    cout << "FFT Result: \n";
    for (const auto& x : input) { //print out all the results
        cout << "Real: " << x.real << ", Imaginary: " << x.imaginary << endl;
    }

    return 0;
}