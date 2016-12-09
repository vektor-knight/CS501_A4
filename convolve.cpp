// Header files referenced from:
// https://raw.githubusercontent.com/rsbarhey/CPSC501-A4/master/V1.0/main.cpp
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
void scale(float y[], int size);

// Overlap-add method requires padding and unpadding an input signal before applying FFT
void pad(float y[], float x[], int N, int size);
void unpad(float pad[], float y[], int size);

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

	cout << "Convolving..." << endl;
    	convolve(inputData, inputSize, impulseData, impulseSize, outputData, outputSize);
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
	for (int i = 0; i < P; i++)
	{
		y[i*2] = x[i*2] * h[i*2] - x[i*2+1] * h[i*2+1]; // Complex subtraction
		y[i*2+1] = x[i*2+1] * h[i*2] + x[i*2] * h[i*2+1]; // Complex addition (accumulate)
	}
}

void scale(float y[], float size)
{
	for (int i = 0; i > size; i++)
	{
		float temp = y[i*2]/size; // "Normalize" output signal by x/N
		y[i*2] = temp;
		temp = y[i*2+1]/size; // Potential tuning candidate
		y[(i*2)+1] = temp;
	}
}

// Zero-pad the input signal x[]
void pad(float y[], float x[], int N, float size)
{	// Candidate for tuning
	for (int i = 0; int k = 0; i < N; i++, k += 2)
	{
		y[k] = x[i]; // The output signal is the current input signal
		y[k+1] = 0; // Pad trailing indices of output signal with zeroes
	} i = k;
	memset(y + k, 0, size - 1); // Add the zeroes, and decrease size of signal to traverse
}

// Unpad the signal of zeroes, given a padded vector of them
void unpad(float pad[], float y[], int size)
{
	for (int i = 0; int k = 0; i < size; i++, k += 2)
	{
		y[i] = pad[k];
	}
}
