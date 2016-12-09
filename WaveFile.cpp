// New WaveFile.cpp for FFT-specific use case.
// Referenced from https://raw.githubusercontent.com/rsbarhey/CPSC501-A4/master/V2.0/WaveFile.cpp


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "WaveFile.h"

float* WaveFile::ReadInput(char *fileName, float *signal, int *Thesize)
{
    ifstream inFile( fileName, ios::in | ios::binary);

	cout << endl << "Reading file: " << fileName << endl;

	inFile.seekg(ios::beg);

	chunkID = new char[4];
	inFile.read( chunkID, 4);
	cout << "Chunk ID: " << chunkID << endl;

	inFile.read( (char*) &chunkSize, 4);
	cout << "Chunk Size: " << chunkSize << endl;

	format = new char[4];
	inFile.read( format, 4);
	cout << "Format: " << format << endl;

	subChunk1ID = new char[4];
	inFile.read( subChunk1ID, 4);
	cout << "Subchunk1ID: " << subChunk1ID << endl;

	inFile.read( (char*) &subChunk1Size, 4);
	cout << "Subchunk1Size: " << subChunk1Size << endl;

	inFile.read( (char*) &audioFormat, 2);
	cout << "AudioForm: " << audioFormat << endl;

	inFile.read( (char*) &numChannels, 2);
	cout << "Num channels: " << numChannels << endl;

	inFile.read( (char*) &sampleRate, 4);
	cout << "sampleRate: " << sampleRate << endl;

	inFile.read( (char*) &byteRate, 4);
	cout << "Byte rate: " << byteRate << endl;

	inFile.read( (char*) &blockAlign, 2);
	cout << "Block align: " << blockAlign << endl;

	inFile.read( (char*) &bitsPerSample, 2);
	cout << "BitsPerSample: " << bitsPerSample << endl;

	if(subChunk1Size == 18)
	{
		char *garbage;
		garbage = new char[2];
		inFile.read( garbage, 2);
	}

	subChunk2ID = new char[4];
	inFile.read( subChunk2ID, 4);
	cout << "Subchunk2ID: " << subChunk2ID << endl;

	//DataSize
	inFile.read( (char*)&dataSize, 4);
	cout << "Data size: " << dataSize << " bytes" << endl;

	//GetData
	*Thesize = dataSize / 2;
	int size = dataSize / 2;
	fileData = new short[size];
	
	int j;
	for(j = 0 ; j < size-2; j+=3)
	{
		inFile.read((char*) &fileData[j], 2);
		inFile.read((char*) &fileData[j+1], 2);
		inFile.read((char*) &fileData[j+2], 2);
	}
	if(j == size-2)
	{
		inFile.read((char*) &fileData[size-1], 2);
		inFile.read((char*) &fileData[size-2], 2);
	}
	else if(j == size-1)
	{
		inFile.read((char*) &fileData[size-1], 2);			
	}
	
	printf("\nDone reading...now producing signal\n");

	//ProduceSignal
	short val;
	signal = new float[size];
	for(int i = 0; i < size; i++)
	{
		val = fileData[i];
		signal[i] = (val * 1.0) / (pow(2,15) - 1);
		if(signal[i] < -1.0)
			signal[i] = -1.0;

	}
	inFile.close();

	return signal;
}

void WaveFile::writeWaveFile(char* fileName, int numSamples, float* signal)
{
	ofstream outFile( fileName, ios::out | ios::binary);

	cout << endl << "Writing to file: " << fileName << endl;

	/*  Calculate the total number of bytes for the data chunk  */
	int chunkSize = numChannels * numSamples * (bitsPerSample / 8);

	chunkID = "RIFF";
	outFile.write( chunkID, 4);
	cout << "Chunk ID: " << chunkID << endl;

	outFile.write( (char*) &chunkSize, 4);
	cout << "Chunk Size: " << chunkSize << endl;

	format = "WAVE";
	outFile.write( format, 4);
	cout << "Format: " << format << endl;

	outFile.write( subChunk1ID, 4);
	cout << "Subchunk1ID: " << subChunk1ID << endl;

	subChunk1Size = 16;
	outFile.write( (char*) &subChunk1Size, 4);
	cout << "Subchunk1Size: " << subChunk1Size << endl;

	outFile.write( (char*) &audioFormat, 2);
	cout << "AudioForm: " << audioFormat << endl;

	outFile.write( (char*) &numChannels, 2);
	cout << "Num channels: " << numChannels << endl;

	outFile.write( (char*) &sampleRate, 4);
	cout << "sampleRate: " << sampleRate << endl;

	outFile.write( (char*) &byteRate, 4);
	cout << "Byte rate: " << byteRate << endl;

	outFile.write( (char*) &blockAlign, 2);
	cout << "Block align: " << blockAlign << endl;

	outFile.write( (char*) &bitsPerSample, 2);
	cout << "BitsPerSample: " << bitsPerSample << endl;

	outFile.write( subChunk2ID, 4);
	cout << "Subchunk2ID: " << subChunk2ID << endl;

	//Data size
	dataSize = numSamples * 2;
	outFile.write( (char*)&dataSize, 4);
	cout << "Data size: " << dataSize << " bytes" << endl;

	int16_t val;
	for(int i = 0; i < numSamples; i++)
	{
		val = (int16_t)(signal[i] * (pow(2,15) - 1));
		outFile.write((char*)&val, 2);
	}
	outFile.close();
}
