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

void convertToPCM(vector<int>* samples, int numChannels, vector<unsigned char>& pcmSamples, int byteDepth) {
    pcmSamples.clear();

    cout << "first 10 samples in all channels" << endl;

    for (int i = 0; i < samples[0].size(); i ++){
        for (int j = 0; j < numChannels; j ++){
            int val = samples[j][i];
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

void readData(vector<int>*& data, FILE *wavFile,  int bytesPerSample, int numChannels){
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
            data[j].push_back(sample);
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
        cout << endl;
    }

    filePath = input.c_str();

    wavFile = fopen( filePath , "rb");

    if(wavFile == NULL){
        printf("Can not able to open wave file\n");
        return -1;
    }

    fread(&wavHeader,headerSize,1,wavFile);

    int bytesPerSample = wavHeader.bitsPerSample/8;
    int numChannels = wavHeader.NumOfChan;
    
    vector<int>* data = new vector<int>[numChannels];
    readData(data, wavFile, bytesPerSample, numChannels);
    fclose(wavFile);

    int numSamples = wavHeader.Subchunk2Size / (numChannels * bytesPerSample);  
    cout << data[0].size() << " " << data[1].size() << " " << numSamples << endl;

    vector<unsigned char> pcmSamples;
    convertToPCM(data, numChannels, pcmSamples, bytesPerSample);

    string outputFilePath = "output.wav";
    writeWavFile(outputFilePath,wavHeader,pcmSamples);

    delete[] data;
    return 0;
}