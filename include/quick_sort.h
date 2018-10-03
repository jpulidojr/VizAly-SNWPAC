#ifndef _QUICK_H
#define _QUICK_H

static int CmpDbl(const void * a, const void * b);
static int CmpDblInd(const void *a, const void *b);

void merge(double * A, double * B, size_t m, size_t n, int dir);
void arraymerge(double * a, size_t size, size_t * index, int N, int dir);
void quick_sort(double * output, double * input, int dir, size_t total);

void quick_sort_index(size_t * index, double * input, int dir, size_t total);
void arraymerge_index(std::pair<double,size_t> * a, size_t size, size_t * index, int N, int dir);
void merge_index(std::pair<double,size_t> * A, std::pair<double,size_t> * B, size_t m, size_t n, int dir);

struct ind_val 
{
    size_t index;
    double value;
};

#endif