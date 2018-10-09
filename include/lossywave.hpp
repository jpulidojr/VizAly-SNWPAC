#ifndef LOSSYWAVE_H
#define LOSSYWAVE_H


#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#include <stdlib.h>
#endif

#include <cstdint>

// interface includes
//#include <wavelet.h>

namespace lossywave
{

	class lossywave
	{
	public:

		EXPORT lossywave();
		EXPORT lossywave(int * inparams);

		// data = n-D data array
		// dataType = sizeof(dataType)
		// output = byte stream of compressed data
		EXPORT size_t compress(void * data, size_t dataType, void *&output);

		// data = compressed byte stream
		// output = n-D data array
		EXPORT size_t decompress(void *data, void *& output);

	protected:
		int pcnt, lvl,  nthreads;
		int * params;
		int mode;

		template <typename T>
		void analyze(T * data);

		template <typename T>
		size_t encode(T * in, void *& out);

		template <typename T>
		size_t decode(void * in, T *& out);
	};


}

#endif