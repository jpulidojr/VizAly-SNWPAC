#ifndef _WRITER_H_
#define _WRITER_H_

#include<vector>
#include<string>

int readConfigw();
int arr_to_pvti_2d(float ** in, char * filename);
int arr_to_pvti_3d(float *** in, char * filename,int dx, int dy, int dz);
template<typename T> int arr_to_vti_3d_range_scale(T **** in, const char * filename, int xl, int xu, int yl, int yu, int zl, int zu, std::vector<std::string> names);

int arr_to_pvti(std::vector<int*>in,int fields, char * filename,int dx, int dy, int dz);
int arr_to_bin_3d_range_scale(double **** in, int num_fields,int header_sz, const char * filename,int xl,int xu, int yl,int yu, int zl, int zu);


int * coeff_args(unsigned short type, unsigned short lvl, unsigned char reg, unsigned short pad, int lx, int ly, int lz, int gx, int gy, int gz, char prec);
int free_args(int * args);

template <typename T> int save_coefficients(T * data, const char * filename, size_t total);
template <typename T> int save_coefficients_md(T * data, const char * filename, int * args);
int save_sort(size_t * sort, const char * filename, size_t total);

#endif

#include "writer.inl"
