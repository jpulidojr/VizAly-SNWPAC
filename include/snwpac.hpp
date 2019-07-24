#ifndef SNWPAC_H
#define SNWPAC_H


#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#include <stdlib.h>
#endif

#include <cstdint>
#include <string>
#include <sstream>

// interface includes
//#include <wavelet.h>

namespace snwpac
{

	class snwpac
	{
	public:

		EXPORT snwpac();
		EXPORT snwpac(int * inparams);
        EXPORT snwpac(int * inparams, bool verbose);
		EXPORT ~snwpac();

		// data = n-D data array
		// dataType = sizeof(dataType)
		// output = byte stream of compressed data
		EXPORT size_t compress(void * data, size_t dataType, void *&output);

		// data = compressed byte stream
		// output = n-D data array
		EXPORT size_t decompress(void *data, void *& output);

		// Prints the internal state parameters for the compressor
		EXPORT void printParams();

		// Prints the header information for a compressed data stream
		EXPORT void printHeader(void * data);

		// Returns header information on a compressed data stream
		EXPORT int * peek(void * data);

	protected:
		int pcnt, lvl,  nthreads;
		int * params;
		int mode; 
		
		// Debug/print related variables
		int verbose; 
		std::stringstream coutBuff;
		std::streambuf * old;

		// RLE analysis function
		template <typename T>
		void analyze(T * data);

		// Coefficient encoding with RLE or LZ4 with int quant
		template <typename T>
		size_t encode(T * in, void *& out);

		// Coefficient decoding with RLE or LZ4 with int quant
		template <typename T>
		size_t decode(void * in, T *& out);
	};


}

#endif