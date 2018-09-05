#include <lossywave.hpp>

#include <iostream>
#include <stdio.h>
#include <complex>
#include <time.h>
#include <omp.h>

// Wavelet functions
#include "wavelet.h"

// Sorting algorithms
#include <algorithm>
#include <gsl/gsl_sort.h>
#include "quick_sort.h"
//#include "bitonic_sort.h"

// Resolves Issue: https://stackoverflow.com/questions/30412951/unresolved-external-symbol-imp-fprintf-and-imp-iob-func-sdl2
// Only occurs when using pre-build GSL libraries from VS2013--
// Solution: Rebuild GSL with VS2015++
#if defined(_MSC_VER) && (_MSC_VER >= 1900)
FILE _iob[] = { *stdin, *stdout, *stderr };
extern "C" FILE * __cdecl __iob_func(void)
{
	return _iob;
}
#endif

namespace lossywave
{

	lossywave::lossywave()
	{
		params = NULL;
		pcnt = 10;
		lvl = 0;
		nthreads = 1;
	}

	lossywave::lossywave(int * inparams)
	{
		params = inparams;
		pcnt = 10;
		lvl = 0;
		nthreads = 1;
	}

	size_t lossywave::compress(void * data, size_t dataType, void *&output)
	{
		gsl_wavelet *w;
		gsl_wavelet_workspace *work;
		size_t type = 303; // Cubic B-splines
		w = gsl_wavelet_alloc(gsl_wavelet_bspline, type);
		gsl_wavelet_print(w);

		// wavelet transform the data

		// sort data and threshold by pcnt

		// requantize and encoded coefficients
		
		size_t comp_size;
		return comp_size;
	}


	size_t lossywave::decompress(void *data, void * output)
	{
		// decode coefficients and requantize

		// wavelet transform the data


		size_t out_size;
		return out_size;
	}


}