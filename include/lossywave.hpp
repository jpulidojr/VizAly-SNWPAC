#ifndef LOSSYWAVE_H
#define LOSSYWAVE_H


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

namespace lossywave
{

	class lossywave
	{
	public:

		EXPORT lossywave();
		EXPORT lossywave(int * inparams);
		EXPORT ~lossywave();

		// data = n-D data array
		// dataType = sizeof(dataType)
		// output = byte stream of compressed data
		EXPORT size_t compress(void * data, size_t dataType, void *&output);

		// data = compressed byte stream
		// output = n-D data array
		EXPORT size_t decompress(void *data, void *& output);

		// Prints the internal parameters passed during init
		EXPORT void printParams();

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