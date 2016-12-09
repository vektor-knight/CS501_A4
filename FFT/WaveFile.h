// References:
// Function prototypes for read/write: https://raw.githubusercontent.com/rsbarhey/CPSC501-A4/master/V1.0/WaveFile.h
// Canonical WAVE file format: http://soundfile.sapp.org/doc/WaveFormat/wav-sound-format.gif

#ifndef WAVEFILE
#define WAVEFILE

using namespace std;

class WaveFile
{
public:
        float* ReadInput(char *filename, float *signal, int *Thesize);
        void writeWaveFile(char* fileName, int numSamples, float *size);

	// Segment 1: "RIFF" chunk descriptor
        char *chunkID;
	int chunkSize;
	char *format;

	// Segment 2: "fmt" sub-chunk
        char *subChunk1ID;
	int subChunk1Size;
	int16_t audioFormat;
	int16_t numChannels;
	int sampleRate;
	int byteRate;
	int16_t blockAlign;
	int16_t bitsPerSample;

	// Segment 3: "data" sub-chunk
        char *subChunk2ID;
        short* fileData;
	int dataSize;
};

#endif
