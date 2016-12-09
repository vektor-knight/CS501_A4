/******************************************************************************
 *
 *     Program:       testtone
 *
 *     Description:   This program generates a two-second 440 Hz test tone and
 *                    writes it to a .wav file.  The sound file has 16-bit
 *                    samples at a sample rate of 44.1 kHz, and is monophonic.
 *
 *     Author:        Leonard Manzara
 *
 *     Date:          November 21, 2009
 *
 ******************************************************************************/


/*  HEADER FILES  ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <ctime>
#include "WaveFile.h"
using namespace std;

#define DEBUG_MODE

/*  CONSTANTS  ***************************************************************/
#define PI					3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY			440.0

/*  Test tone duration in seconds  */
#define DURATION			2.0

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE			44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE		16

/*  Standard sample size in bytes  */
#define BYTES_PER_SAMPLE	(BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC			1
#define STEREOPHONIC		2

#define SWAP(a, b) tempr = (a); (a) = (b); (b) = tempr

/*  FUNCTION PROTOTYPES  *****************************************************/
bool checkFileExist(char* fileName);
void overlapAdd(float* inputData, int inputSize, float* impulseData, int impulseSize, float* outputData, int outputSize);
void padSignal(float output[], float signal[], int signalLen, int size);
void unpadSignal(float padded[], float unpadded[], int size);	
void fft(float data[], int nn, int isign);
void complexMultiplication(float compInput[],float compIR[],float compResult[], int size);
void scaleFFT(float result[], int size);

/******************************************************************************
 *
 *	function:	main
 *
 *	purpose:	Creates the test tone and writes it to the
 *               specified .wav file.
 *
 *   arguments:	argv[1]:  the filename for the output .wav file
 *
 ******************************************************************************/

int main (int argc, char *argv[])
{
    std::clock_t start;
    start = std::clock();
    if(argc < 4 || argc > 4)
    {
        cout << "Error using the program" << endl;
        cout << "please use ./app [InputAudio].wav [ImpulseResponse].wav [OutputAudio].wav" << endl;
        return 0;
    }


    /*  Create the sine wave test tone, using the specified
        frequency, duration, and number of channels, writing to
        a .wav file with the specified output filename  */
    if(!checkFileExist(argv[1]) || !checkFileExist(argv[2]))
    {
        cout << "Error: one or both of input and impulse files don't exist " << endl;
        return 0;
    }
    WaveFile *inputWave = new WaveFile();
    WaveFile *impulseWave = new WaveFile();

    float* inputData;
    int inputSize;
    inputData = inputWave->ReadInput(argv[1], inputData, &inputSize);
    
    float* impulseData;
    int impulseSize;
    impulseData = impulseWave->ReadInput(argv[2], impulseData, &impulseSize);
    
    int outputSize = inputSize + impulseSize - 1;
    float* outputData = new float[outputSize];
    
    overlapAdd(inputData, inputSize, impulseData, impulseSize, outputData, outputSize);
    inputWave->writeWaveFile(argv[3], outputSize, outputData);

    cout << "Time took (ms): " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000)<< endl;
    /*  End of program  */
    return 0;
}

bool checkFileExist(char* fileName)
{
    ifstream f(fileName);
    if (f.good()) {
        f.close();
        return true;
    }
    else {
        f.close();
        return false;
    }
}

void overlapAdd(float* inputData, int inputSize, float* impulseData, int impulseSize, float* outputData, int outputSize)
{
    int totalSize = 0;
    int paddedTotalSize = 1;
    totalSize = inputSize + impulseSize - 1;
    
    int i = 0;
    while (paddedTotalSize < totalSize)
    {
        paddedTotalSize <<= 1; //times to pad    
        i++; 
    }
    
    cout << "Padded size: " << paddedTotalSize << "Exponent: " << i << endl;
    cout << "Input size: " << inputSize << endl;
    cout << "Impulse size: " << impulseSize << endl;
    
    float* compInput = new float[2*paddedTotalSize];
    float* compImpulse = new float[2*paddedTotalSize];
    float* compResult = new float[2*paddedTotalSize];
    
    padSignal(compInput, inputData, inputSize, 2*paddedTotalSize);
    padSignal(compImpulse, impulseData, impulseSize, 2*paddedTotalSize);
    memset(compResult, 0, 2*paddedTotalSize); // padding zeroes
    
    fft(compInput-1, paddedTotalSize, 1);
    fft(compImpulse-1, paddedTotalSize, 1);
    
    cout << "Complex multiplication..." << endl;
    complexMultiplication(compInput, compImpulse, compResult, paddedTotalSize);
    
    fft(compResult-1, paddedTotalSize, -1);
    scaleFFT(compResult, paddedTotalSize);
    unpadSignal(compResult, outputData, outputSize);
}

void padSignal(float output[], float signal[], int signalLen, int size)
{
    int i, k;
    for(i = 0, k = 0; i<signalLen; i++, k+=2)
    {
        output[k] = signal[i];
        output[k+1] = 0;    
    }
    i = k;
    memset(output + k, 0, size -1); //adding zeroes
}

void unpadSignal(float padded[], float unpadded[], int size)
{
    int i, j;
    for(i = 0, j = 0; i <size; i++, j+=2)
    {
        unpadded[i] = padded[j];    
    }
}

void fft(float data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            SWAP(data[j], data[i]);
    		   SWAP(data[j+1], data[i+1]);
    	  }
        m = nn;
        while (m >= 2 && j > m)
    	  {
            j -= m;
    		   m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax)
    {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

void complexMultiplication(float compInput[],float compImpulse[],float compResult[], int size)
{
    int i = 0;
    int indexUsed;
    for(i = 0; i<size; i++)
    {
        indexUsed = i+i;
        compResult[indexUsed] = compInput[indexUsed] * compImpulse[indexUsed] - compInput[indexUsed+1] * compImpulse[indexUsed+1];
	     compResult[indexUsed+1] = compInput[indexUsed+1] * compImpulse[indexUsed] + compInput[indexUsed] * compImpulse[indexUsed+1];
    }
}

void scaleFFT(float result[], int size)
{
    int i = 0;
    for(i = 0; i < size; i++) 
    {
        result[i*2] /= (float)size;
        result[(i*2)+1] /= (float)size;
    }
}
