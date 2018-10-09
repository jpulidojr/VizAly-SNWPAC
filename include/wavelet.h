#ifndef _WAVELET_H_
#define _WAVELET_H_

#include <gsl/gsl_wavelet.h>
//#include <gsl/gsl_wavelet2d.h>

// N-to-M conversion helper functions
template <typename T> T *** arrayto3D(T * data, size_t x, size_t y, size_t z);
template <typename T> T ** arrayto2D(T * data, size_t x, size_t y);
template <typename T> T * arrayto1D(T *** data, size_t x, size_t y, size_t z);
template <typename T> T * arrayto1D(T ** data, size_t x, size_t y);

// Helper functions to compute size of wavelet hierarchy at a level
static int binary_logn (const size_t n);
int binary_upto_logn (const size_t n);
int binary_upto_value (const size_t n, const int log_n);

// Debug output of discrete wavelet signal
void gsl_wavelet_print (const gsl_wavelet * w);

// Native stepping functions that apply signal to data
template <typename T> static void dwt_step (const gsl_wavelet * w, T *a, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);
template <typename T> static void dwt_step_proto(const gsl_wavelet * w, T *a, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);

// 1D-support
template <typename T> int wavelet_transform(const gsl_wavelet * w, T *&data, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);

// 2D-support 
template <typename T> int wavelet2d_transform(const gsl_wavelet * w, T *data, size_t tda, size_t size1, size_t size2, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);
template <typename T> int wavelet2d_nstransform(const gsl_wavelet * w, T *data, size_t tda, size_t size1, size_t size2, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);

// 3D-support versions of (GSL) wavelet transform
template <typename T> int gsl_wavelet3d_transform (const gsl_wavelet * w, T *data, size_t tda, size_t size1, size_t size2,size_t size3, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);
template <typename T> int gsl_wavelet3d_nstransform (const gsl_wavelet * w, T *data, size_t tda, size_t size1, size_t size2, size_t size3, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);
template <typename T> int gsl_wavelet3d_nstransform (const gsl_wavelet * w, T *data, size_t tda, size_t size1, size_t size2, size_t size3, gsl_wavelet_direction dir, gsl_wavelet_workspace * work, int levels);
template <typename T> int gsl_wavelet3d_nstransform_proto(const gsl_wavelet * w, T *data, size_t tda, size_t size1, size_t size2, size_t size3, gsl_wavelet_direction dir, gsl_wavelet_workspace * work);

// Functions to compute dimensions for each coefficient block (quadrant) in a wavelet hierarchy
int *** getCoeffDims (int dimx, int dimy, int dimz);
int *** getCoeffDims (int dimx, int dimy);

// Set/Extract functions used in conjunction with above to extract coefficient subregions
template <typename T> T * getCutout3d(T * in, int * codims, int * dims);
template <typename T> T * getCutout2d(T * in, int * codims, int * dims);
template <typename T> void setCutout3d(T * in, T * cutout, int * codims, int * dims);
template <typename T> void setCutout2d(T * in, T * cutout, int * codims, int * dims);

template <typename T> int grow(int val, T *** data, size_t x, size_t y, size_t z);
template <typename T> int shrink(int val, T *** data, size_t x, size_t y, size_t z);

#endif

#include "wavelet.inl"
