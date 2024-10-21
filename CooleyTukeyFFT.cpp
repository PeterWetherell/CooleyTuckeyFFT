#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <string.h>

using namespace std;
using std::string;
using std::fstream;


typedef struct  WAV_HEADER{
    char                RIFF[4];        // RIFF Header      Magic header
    unsigned long       ChunkSize;      // RIFF Chunk Size  
    char                WAVE[4];        // WAVE Header      
    char                fmt[4];         // FMT header       
    unsigned long       Subchunk1Size;  // Size of the fmt chunk                                
    unsigned short      AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM 
    unsigned short      NumOfChan;      // Number of channels 1=Mono 2=Sterio                   
    unsigned long       SamplesPerSec;  // Sampling Frequency in Hz                             
    unsigned long       bytesPerSec;    // bytes per second 
    unsigned short      blockAlign;     // 2=16-bit mono, 4=16-bit stereo 
    unsigned short      bitsPerSample;  // Number of bits per sample      
    char                Subchunk2ID[4]; // "data"  string   
    unsigned long       Subchunk2Size;  // Sampled data length    

} wav_hdr; 

const double PI = 3.141592653589793238460;

class complex{
    public:

    double real;
    double imaginary;

    complex(){
        real = 0;
        imaginary = 0;
    }

    double getMagnitude(){
        return real*real + imaginary*imaginary;
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

int blockSize = 1024;

//apply Cooley-Tukey Recursively
void CooleyTukeyInPlaceRec(int startIndex, int start, int indexing, vector<complex>& values, complex* twiddle, complex* savedTwiddle){//convert the algorithm to work in place
    if (blockSize == indexing){ //Base case -> no more edits to the array
        return;
    }
    
    //In place recurse on the even and odd
    CooleyTukeyInPlaceRec(startIndex, start,indexing*2, values, twiddle, savedTwiddle); //Even
    CooleyTukeyInPlaceRec(startIndex, start + indexing,indexing*2, values, twiddle, savedTwiddle); //Odd

    int halfN = blockSize/(2*indexing);
    
    for (int i = 0 ; i < halfN; i ++){ //combine both even and odd
        //cout << values[start+indexing + i*2*indexing].real << " " << values[start + i*2*indexing].real << " " << twiddle[i*indexing].real << ", ";

        complex t = twiddle[i*indexing] * values[startIndex + start + indexing + i*2*indexing];

        //Do the first half
        values[startIndex + start + indexing*i] = values[startIndex + start + i*2*indexing] + t; //even + t

        //Save the twiddle factor
        savedTwiddle[i] = complex(-2.0,0)*t;
    }
    for (int i = 0 ; i < halfN; i ++){ //combine both even and odd
        values[startIndex + start + indexing*(i+halfN)] = values[startIndex + start + indexing*i] + savedTwiddle[i]; //even - t = (even+t) - 2*t
    }
}

//Setup for the inplace
void cooleyTukeyInPlace(vector<complex>& values){ 
    int size = values.size();
    complex* twiddle = new complex[size/2];
    for (int i = 0 ; i < size/2; i ++){ //calculate all the twidles.
        double theta =  -2 * PI * i / size;
        twiddle[i] = complex(cos(theta), sin(theta));
    }
    complex* savedTwiddle = new complex[size/2];
    for (int i = 0; i < size-blockSize; i += blockSize){
        CooleyTukeyInPlaceRec(i, 0, 1, values, twiddle, savedTwiddle);
    }
}

void inverseCooleyTukeyInPlace(vector<complex>& values){
    int size = values.size();

    for (int i = 0; i < size ; i ++){ //get the complex conjugate
        values[i].imaginary = -values[i].imaginary; // Conjugate
    }

    cooleyTukeyInPlace(values); //Do FFT on conjugate

    for (int i = 0; i < size ; i ++){ //get the complex conjugate again
        values[i].imaginary = -values[i].imaginary; // Conjugate
    }
    for (int i = 0; i < (size/blockSize)*blockSize; i++) { //Normalize everything
        values[i].real /= blockSize;
        values[i].imaginary /= blockSize;
    }

}

void removeWaves(double percentage, vector<complex>& values){
    int size = values.size();

    int numRemove = blockSize*percentage;
    for (int i = 0; i < size-blockSize; i += blockSize){
        for (int j = 0; j < numRemove; j ++){
            int minIndex = -1;
            double minMag = -1;
            for (int k = 1; k < blockSize; k ++){
                double a = values[i + k].real; //.getMagnitude();
                if (a != 0 && (minIndex == -1 || a < minMag)){
                    minMag = a;
                    minIndex = k;
                }
            }
            if (minIndex != -1){
                values[i+minIndex] = complex(0,0);
            }
        }
    }
}


int getFileSize(FILE *inFile){
    int fileSize = 0;
    fseek(inFile,0,SEEK_END);

    fileSize=ftell(inFile);

    fseek(inFile,0,SEEK_SET);
    return fileSize;
}

void writeWavFile(const string& outputFilePath, const wav_hdr& header, const vector<unsigned char>& pcmSamples) {
    FILE* outFile = fopen(outputFilePath.c_str(), "wb");
    if (!outFile) {
        cerr << "Could not open output file for writing." << endl;
        return;
    }

    // Write the WAV header
    fwrite(&header, sizeof(wav_hdr), 1, outFile);

    // Write PCM samples
    fwrite(pcmSamples.data(), 1, pcmSamples.size(), outFile);

    fclose(outFile);
}

void convertToPCM(vector<complex>* samples, int numChannels, vector<unsigned char>& pcmSamples, int byteDepth) {
    pcmSamples.clear();

    cout << "first 10 samples in all channels" << endl;

    for (int i = 0; i < samples[0].size(); i ++){
        for (int j = 0; j < numChannels; j ++){
            int val = samples[j][i].real;
            for (int k = 0; k < byteDepth; k ++){
                unsigned char byte = (val >> (k * 8)) & 0xFF;
                if (i < 10){
                    cout << (int)byte << " ";
                }
                pcmSamples.push_back(byte);
            }
        }
    }
}

void readData(vector<complex>*& data, FILE *wavFile,  int bytesPerSample, int numChannels){
    unsigned char sampleBytes[bytesPerSample];

    bool keepReading = true;

    cout << "first 10 samples in all channels" << endl;

    for (int i = 0; keepReading; i ++){
        for (int j = 0; j < numChannels; j ++){
            memset(sampleBytes, 0 , sizeof(sampleBytes));
            if (fread(sampleBytes,1,bytesPerSample,wavFile) != bytesPerSample){
                cout << endl << "Error or end of file reached after reading " << i << " samples." << endl;
                keepReading = false;
                break;
            }
            int sample = 0;
            for (int k = 0; k < bytesPerSample; k ++){ //(sampleBytes[2] << 16) | (sampleBytes[1] << 8) | sampleBytes[0];
                if (i < 10){
                    cout << (int)sampleBytes[k] << " ";
                }
                sample |= sampleBytes[k] << (k*8); //shift the previous contents and load value into the sample
            }
            if (sample & (1 << (bytesPerSample * 8 -1))) { //check if negative
                sample |= ~((1 << (bytesPerSample * 8)) - 1);  // Sign extend
            }
            data[j].push_back( complex( (double) sample, 0.0) );
        }
    }
}

int main(int argc, char *argv[]){
    wav_hdr wavHeader;
    FILE *wavFile;
    int headerSize = sizeof(wav_hdr),filelength = 0;
    string input;

    const char* filePath;

    if (argc > 1){
        input =  argv[1];
    }
    else {
        cout << "Enter path to wav file: ";
        cin >> input;
        cin.get();
    }

    double percentage = .5; //automatically assumes 50% reduction if they haven't entered a value when calling the exe
    if (argc > 2){
        percentage = stod(argv[2]); //conver the second input value to decimal
    }

    cout << "File : " << input << " " << percentage;

    cout << endl;

    filePath = input.c_str();

    wavFile = fopen( filePath , "rb");

    if(wavFile == NULL){
        printf("Can not able to open wave file\n");
        return -1;
    }

    fread(&wavHeader,headerSize,1,wavFile);

    int bytesPerSample = wavHeader.bitsPerSample/8;
    int numChannels = wavHeader.NumOfChan;
    
    vector<complex>* data = new vector<complex>[numChannels];
    readData(data, wavFile, bytesPerSample, numChannels);
    fclose(wavFile);

    cout << "Finished reading file" << endl;

    for (int i = 0; i < numChannels; i ++){
        cooleyTukeyInPlace(data[i]);

        cout << "Finished converting channel " << i << " to fourier" << endl;

        removeWaves(percentage, data[i]);

        cout << "Finished removing the most important parts" << endl;
        
        inverseCooleyTukeyInPlace(data[i]);

        cout << "Finished converting channel " << i << " back to origional"  << endl;
    }

    vector<unsigned char> pcmSamples;
    convertToPCM(data, numChannels, pcmSamples, bytesPerSample);

    string outputFilePath = "output.wav";
    writeWavFile(outputFilePath,wavHeader,pcmSamples);

    delete[] data;
    return 0;
}