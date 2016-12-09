// Header files referenced from:
// https://raw.githubusercontent.com/rsbarhey/CPSC501-A4/master/V1.0/main.cpp
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#include "WaveFile.h"

//void convolve(float x[], int N, float h[], int M, float y[], int P);

// Helpers
void complexMul(float x[], float h[], float y[], int P);
void scaleFFT(float result[], int size);

// Overlap-add method requires padding and unpadding an input signal before applying FFT
void padSignal(float output[], float signal[], int signalLen, int size);
void unpad(float padded[], float unpadded[], int size);

// "Heavy lifters"
#define SWAP(a, b) tempr = (a); (a) = (b); (b) = tempr; // From handout (reference below)
void overlapAdd(float* inputData, int inputSize, float* impulseData, int impulseSize, float* outputData, int outputSize);
void overlapFFT(float data[], unsigned long nn, int isign); // From handout (reference below)

int main(int argc, char *argv[]) // To capture CL arguments
{
	// Main algorithm to capture WAVE file inputs, and output convolved WAVE
	// referenced from: https://raw.githubusercontent.com/rsbarhey/CPSC501-A4/master/V1.0/main.cpp
	// Added my own annotations for later FFT development
	std::clock_t start;
	start = std::clock();

	WaveFile *input = new WaveFile();
	WaveFile *impulse = new WaveFile();

	float* inputData;
	int inputSize;
	inputData = input->ReadInput(argv[1], inputData, &inputSize);

	float* impulseData;
	int impulseSize;
	impulseData = impulse->ReadInput(argv[2], impulseData, &impulseSize);

	// Output size should be N + M - 1 according to
	// Smith (p. 112-115).
	int outputSize = inputSize + impulseSize - 1;
	float* outputData = new float[outputSize];

	cout << "Applying FFT" << endl;
    	//convolve(inputData, inputSize, impulseData, impulseSize, outputData, outputSize);
    	overlapAdd(inputData, inputSize, impulseData, impulseSize, outputData, outputSize);
	input->writeWaveFile(argv[3], outputSize, outputData);

   	cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << endl;

    	return 0;
	}

/* Input-side convolution from Smith (p. 112-115).
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
	int n, m;

	for (n = 0; n < P; n++)
	{
		y[n] = 0.0;
	}
	for (n = 0; n < N; n++)
	{
		for (m = 0; m < M; m++)
		{
			y[n+m] += x[n] * h[m];
		}
	}

} */


// 
void complexMul(float x[], float h[], float y[], int P)
{
	// Tuning 2: Partial loop unrolling
	for (int i = 0; i < P; i+=4)
	{
		y[i*2] = x[i*2] * h[i*2] - x[i*2+1] * h[i*2+1]; // Complex subtraction
		y[i*2+1] = x[i*2+1] * h[i*2] + x[i*2] * h[i*2+1]; // Complex addition (accumulate)

		y[i*2+1] = x[i*2] * h[i*2] - x[i*2+2] * h[i*2+2]; // Complex subtraction
                y[i*2+2] = x[i*2+2] * h[i*2] + x[i*2] * h[i*2+2]; // Complex addition (accumulate)

  		y[i*2+2] = x[i*2] * h[i*2] - x[i*2+3] * h[i*2+3]; // Complex subtraction
                y[i*2+3] = x[i*2+3] * h[i*2] + x[i*3] * h[i*2+3]; // Complex addition (accumulate)

	 	y[i*2+3] = x[i*2] * h[i*2] - x[i*2+4] * h[i*2+4]; // Complex subtraction
                y[i*2+4] = x[i*2+4] * h[i*2] + x[i*2] * h[i*2+4]; // Complex addition (accumulate)
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

// Zero-pad the input signal x[]
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

// Unpad the signal of zeroes, given a padded vector of them
void unpad(float padded[], float unpadded[], int size)
{
	int i, k;
	for (i = 0, k = 0; i < size; i++, k += 2)
	{
		unpadded[i] = padded[k];
	}
}

// Reference: Class handout (12.2 Fast Fourier Transform (FFT))
// Originally cited from "Numerical Recipes in C"; the Danielson-Lanczos Formula
void overlapFFT(float data[], unsigned long nn, int isign)
{
	// Use double precision for the trig recurrence during Danielson-Lanczos
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;

	n = nn << 1;
	j = 1;
	for (int i = 1; i < n; i += 2)
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
		theta = isign * (6.28318530717959/mmax); // Candidate for tuning?
		wtemp = sin(0.5 * theta);
		// Tuning 1: -2 -> -2.0 (signed int to double)
		wpr = -2.0*wtemp*wtemp; // Candidate for tuning? (-2 -> -2.0)
		wpi = sin(theta);
		wr = 1;
		wi = 0;
		for (m = 1; m < mmax; m += 2)
		{
			for (i = m; i <= n; i += istep)
			{
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr*data[j+1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr; // Trig recurrence
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

// Overlap-add component, Smith p. 311-318
void overlapAdd(float* inputData, int inputSize, float* impulseData, int impulseSize, float* outputData, int outputSize)
{
	int size = 0;
	int paddedSize = 1;
	size = inputSize + impulseSize - 1; // candidate for tuning?

	int i = 0;
	while (paddedSize < size)
	{
		paddedSize <<= 1;
		i++;
	}

	// Setting up padding and complex multiplication
	// First code tuning: eliminate common sub-expression
		// eg. Assignment 2*paddedSize to a variable

	float* cIn = new float[2*paddedSize]; //
	float* cImpulse = new float[2*paddedSize];//
	float* cResult = new float[2*paddedSize];//

	padSignal(cIn, inputData, inputSize, 2*paddedSize);
	padSignal(cImpulse, impulseData, impulseSize, 2*paddedSize);
	memset(cResult, 0, 2*paddedSize);

	overlapFFT(cIn - 1, paddedSize, 1);
	overlapFFT(cImpulse - 1, paddedSize, 1);
	complexMul(cIn, cImpulse, cResult, paddedSize);

	overlapFFT(cResult - 1, paddedSize, -1);
	scaleFFT(cResult, paddedSize);
	unpad(cResult, outputData, outputSize);
}
