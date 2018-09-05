#ifndef _BITONIC_H
#define _BITONIC_H

void bitonic_sort(double * output, double * input, size_t stride, size_t total);
void bitonic_sort_seq(size_t start, size_t length, double * seq, int flag);
void bitonic_sort_par(size_t start, size_t length, double * seq, int flag);
void swap(double *a, double *b);

void bitonic_sort_index(size_t * output, double * input, size_t stride, size_t total);
void bitonic_sort_seq_index(size_t start, size_t length, double *seq, size_t * index, int flag);
void bitonic_sort_par_index(size_t start, size_t length, double *seq, size_t * index, int flag);
void swap(double *a, double *b, size_t *ai, size_t *bi);


#endif