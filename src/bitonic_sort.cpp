/* Code based off of:
   http://www.cs.rutgers.edu/~venugopa/parallel_summer2012/bitonic_openmp.html
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "bitonic_sort.h"

//#define MAX(A, B) (((A) > (B)) ? (A) : (B))
//#define MIN(A, B) (((A) > (B)) ? (B) : (A))
#define UP 0
#define DOWN 1

size_t m;

void bitonic_sort(double * output, double * input, size_t stride, size_t total)
{
    // copy into output
    int i, j;
    for( i=0; i<total; i++)
    {
        output[i]=input[i];
    }

    int flag;

    double startTime, elapsedTime; /* for checking/testing timing */
    double clockZero = 0.0;

    int numThreads;

    // start
    //startTime = walltime( &clockZero );

    numThreads =  omp_get_max_threads();

    // making sure input is okay
    if ( total < numThreads * 2 )
    {
        printf("The size of the sequence is less than 2 * the number of processes.\n");
        exit(0);
    }

    // the size of sub part
    m = total / numThreads;

    // make the sequence bitonic - part 1
    for (i = 2; i <= m; i = 2 * i)
    {
        #pragma omp parallel for shared(i, output) private(j, flag)
        for (j = 0; j < total; j += i)
        {
            if ((j / i) % 2 == 0)
                flag = UP;
            else
                flag = DOWN;
            bitonic_sort_seq(j, i, output, flag);
        }
    }

    // make the sequence bitonic - part 2
    for (i = 2; i <= numThreads; i = 2 * i)
    {
        for (j = 0; j < numThreads; j += i)
        {
            if ((j / i) % 2 == 0)
                flag = UP;
            else
                flag = DOWN;
            bitonic_sort_par(j*m, i*m, output, flag);
        }
        #pragma omp parallel for shared(j)
        for (j = 0; j < numThreads; j++)
        {
            if (j < i)
                flag = UP;
            else
                flag = DOWN;
            bitonic_sort_seq(j*m, m, output, flag);
        }
    }

    // bitonic sort
    //bitonic_sort_par(0, n, seq, UP);
    //bitonic_sort_seq(0, n, seq, UP);

    //end
    //elapsedTime = walltime( &startTime );

    /*
    // print a sequence
    for (i = 0; i < n; i++){
      printf("%d ", seq[i]);
    }
    printf("\n");
    */

    //printf("Elapsed time = %.2f sec.\n", elapsedTime);

    //free(seq);
}

void bitonic_sort_index(size_t * output, double * input, size_t stride, size_t total){
    // copy into output
    int i, j;
    double * tmp = new double[total];
    for( size_t i=0; i<total; i++)
    {
        output[i]=i; // THis array keeps track of the sorting indeces
        tmp[i]=input[i]; // This array will be changed and sorted
    }

    int flag;

    double startTime, elapsedTime; /* for checking/testing timing */
    double clockZero = 0.0;

    int numThreads;

    // start
    //startTime = walltime( &clockZero );

    numThreads =  omp_get_max_threads();

    // making sure input is okay
    if ( total < numThreads * 2 )
    {
        printf("The size of the sequence is less than 2 * the number of processes.\n");
        exit(0);
    }

    // the size of sub part
    m = total / numThreads;

    // make the sequence bitonic - part 1
    for (i = 2; i <= m; i = 2 * i)
    {
        #pragma omp parallel for shared(i, output,tmp) private(j, flag)
        for (j = 0; j < total; j += i)
        {
            if ((j / i) % 2 == 0)
                flag = UP;
            else
                flag = DOWN;
            bitonic_sort_seq_index(j, i, tmp, output, flag);
        }
    }
    
    // make the sequence bitonic - part 2
    for (i = 2; i <= numThreads; i = 2 * i)
    {
        for (j = 0; j < numThreads; j += i)
        {
            if ((j / i) % 2 == 0)
                flag = UP;
            else
                flag = DOWN;
            bitonic_sort_par_index(j*m, i*m, tmp, output, flag);
        }
        #pragma omp parallel for shared(j)
        for (j = 0; j < numThreads; j++)
        {
            if (j < i)
                flag = UP;
            else
                flag = DOWN;
            bitonic_sort_seq_index(j*m, m, tmp, output, flag);
        }
    }

    //bitonic_sort_par_index(0, total, tmp,output, UP);
    //bitonic_sort_seq_index(0, total, tmp,output, UP);

    delete[] tmp;

}
void bitonic_sort_seq(size_t start, size_t length, double *seq, int flag)
{
    size_t i;
    size_t split_length;

    if (length == 1)
        return;

    if (length % 2 !=0 )
    {
        printf("error\n");
        exit(0);
    }

    split_length = length / 2;

    // bitonic split
    for (i = start; i < start + split_length; i++)
    {
        if (flag == UP)
        {
            if (seq[i] > seq[i + split_length])
                swap(&seq[i], &seq[i + split_length]);
        }
        else
        {
            if (seq[i] < seq[i + split_length])
                swap(&seq[i], &seq[i + split_length]);
        }
    }

    bitonic_sort_seq(start, split_length, seq, flag);
    bitonic_sort_seq(start + split_length, split_length, seq, flag);
}

void bitonic_sort_seq_index(size_t start, size_t length, double *seq, size_t *index, int flag)
{
    size_t i;
    size_t split_length;

    if (length == 1)
        return;

    if (length % 2 !=0 )
    {
        printf("error\n");
        exit(0);
    }

    split_length = length / 2;

    // bitonic split
    for (i = start; i < start + split_length; i++)
    {
        if (flag == UP)
        {
            if (seq[i] > seq[i + split_length])
                swap(&seq[i], &seq[i + split_length], &index[i], &index[i + split_length]);
        }
        else
        {
            if (seq[i] < seq[i + split_length])
                swap(&seq[i], &seq[i + split_length], &index[i], &index[i + split_length]);
        }
    }

    bitonic_sort_seq_index(start, split_length, seq, index, flag);
    bitonic_sort_seq_index(start + split_length, split_length, seq, index, flag);
}

void bitonic_sort_par(size_t start, size_t length, double *seq, int flag)
{
    int i;
    size_t split_length;

    if (length == 1)
        return;

    if (length % 2 !=0 )
    {
        printf("The length of a (sub)sequence is not divided by 2.\n");
        exit(0);
    }

    split_length = length / 2;

    // bitonic split
    #pragma omp parallel for shared(seq, flag, start, split_length) private(i)
    for (i = start; i < start + split_length; i++)
    {
        if (flag == UP)
        {
            if (seq[i] > seq[i + split_length])
                swap(&seq[i], &seq[i + split_length]);
        }
        else
        {
            if (seq[i] < seq[i + split_length])
                swap(&seq[i], &seq[i + split_length]);
        }
    }

    if (split_length > m)
    {
        // m is the size of sub part-> n/numThreads
        bitonic_sort_par(start, split_length, seq, flag);
        bitonic_sort_par(start + split_length, split_length, seq, flag);
    }

    return;
}

void bitonic_sort_par_index(size_t start, size_t length, double *seq,size_t *index, int flag)
{
    int i;
    size_t split_length;

    if (length == 1)
        return;

    if (length % 2 !=0 )
    {
        printf("The length of a (sub)sequence is not divided by 2.\n");
        exit(0);
    }

    split_length = length / 2;

    // bitonic split
    #pragma omp parallel for shared(seq, index, flag, start, split_length) private(i)
    for (i = start; i < start + split_length; i++)
    {
        if (flag == UP)
        {
            if (seq[i] > seq[i + split_length])
                swap(&seq[i], &seq[i + split_length], &index[i], &index[i + split_length] );
        }
        else
        {
            if (seq[i] < seq[i + split_length])
                swap(&seq[i], &seq[i + split_length], &index[i], &index[i + split_length] );
        }
    }

    if (split_length > m)
    {
        // m is the size of sub part-> n/numThreads
        bitonic_sort_par_index(start, split_length, seq, index, flag);
        bitonic_sort_par_index(start + split_length, split_length, seq, index, flag);
    }

    return;
}

void swap(double *a, double *b)
{
    double t;
    t = *a;
    *a = *b;
    *b = t;
}

void swap(double *a, double *b, size_t *ai, size_t *bi)
{
    double t;
    t = *a;
    *a = *b;
    *b = t;

    size_t ti;
    ti = *ai;
    *ai = *bi;
    *bi = ti;
}