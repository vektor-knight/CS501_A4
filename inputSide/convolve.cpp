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

// Prototypes
void convolve(float x[], int N, float h[], int M, float y[], int P);

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

// Input-side convolution from Smith (p. 112-115).
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

}

