#ifndef _READER_H_
#define _READER_H_

int binto2Darray(float ** in, char * filename,int xl,int xu, int yl,int yu, int zl, int zu );
int binto3Darray(double *** in, char * filename,int xl,int xu, int yl,int yu, int zl, int zu );
int binto3DarrayMPI(double *** in, char * filename,int xl,int xu, int yl,int yu, int zl, int zu, int field);

template <typename T> T * read_coefficients(const char * filename, size_t total);
template <typename T> T * read_coefficients_md(const char * filename, int * args);
int * read_md(const char * filename);
template <typename T> T * read_adaptive_coefficients(const char * filename, size_t total);
size_t * read_sort(const char * filename, size_t total);

int readConfig();
double findByteAddress(int x, int y, int z, int field);
int * findXYZAddress(double byteLocation, int field);

int findIndexXYZ(int x, int y, int z);
int * findXYZIndex(int index);

int dbto3Darray(double *** in,int xl,int xu, int yl,int yu, int zl, int zu, int field);

#endif

#include "reader.inl"