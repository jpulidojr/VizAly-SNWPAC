// Author: Jesus Pulido
// jpulido@ucdavis.edu

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <functional>
#include <algorithm>
#include <vector>

#include "quick_sort.h"

static int CmpDbl(const void *a, const void *b){
    double aa = *(double*)a;
    double bb = *(double*)b;
    if(aa<bb)
        return -1;
    else if (aa>bb)
        return 1;
    else
        return 0;
   /* if ( *(double*)a < *(double*)b ) return -1;
    if ( *(double*)a > *(double*)b ) return 1;
    if ( *(double*)a == *(double*)b ) return 0;*/
}

static int CmpDblInd(const void *a, const void *b){
    ind_val aa = *(ind_val*)a;
    ind_val bb = *(ind_val*)b;
    if(aa.value<bb.value)
        return -1;
    else if (aa.value>bb.value)
        return 1;
    else
        return 0;
}

void merge(double * A, double * B, size_t m, size_t n, int dir){
    int i=0, j=0, k=0;
    size_t size = m+n;

    double *C = new double[size];

    if(dir>=1) //Ascending
    {
        while ( i < m && j < n )
        {
            if(A[i] <= B[j]) 
                C[k] = A[i++];
            else 
                C[k] = B[j++];
            k++;
        }
    }
    else //Descending
    {
        while ( i < m && j < n )
        {
            if(A[i] >= B[j]) 
                C[k] = A[i++];
            else 
                C[k] = B[j++];
            k++;
        }
    }
    if (i < m) 
        for (int p = i; p < m; p++,k++) 
            C[k] = A[p];
    else 
        for (int p = j; p < n; p++,k++) 
            C[k] = B[p];

    for( i = 0; i < size; i++ ) 
        A[i] = C[i];
    
    delete[] C;
}

void arraymerge(double *a, size_t size, size_t *index, int N, int dir){
#ifdef _WIN32
    long long i;
#else
	size_t i;
#endif

    while( N > 1 )
    {
        for ( i = 0; i < N; i++ ) 
            index[i] = i*size/N;
        index[N]=size;

#pragma omp parallel for private(i)
        for( i=0; i<N; i+=2 )
        {
            //fprintf(stderr,"merging %d and %d, index %d and %d (up to %d)\n",i,i+1,index[i],index[i+1],index[i+2]);
            merge(a+index[i], a+index[i+1], index[i+1]-index[i], index[i+2]-index[i+1], dir);
            //for(int i=0; i<size; i++) fprintf(stderr,"after: %d %f\n",i,a[i]);
        }
        N /= 2;

    }
}

void quick_sort(double * output, double * input, int dir, size_t total){
    // copy into output
#ifdef _WIN32
    long long i;
#else
	size_t i;
#endif
    /*for( i=0; i<total; i++)
    {
        output[i]=input[i];
    }*/

    // Lets speed things up
    memcpy(output,input,total*sizeof(double));

    //omp_set_num_threads(1);

    int numThreads = omp_get_max_threads();
    size_t * index = new size_t[numThreads+1];
    for(i=0; i<numThreads; i++) 
        index[i] = i*total/numThreads; 
    index[numThreads] = total;

    double start = omp_get_wtime();
    // Sort each individual thread/process
#pragma omp parallel for private(i)
    for(i=0; i<numThreads; i++)
    {
        if (dir>=1) //ascending
            std::sort(output+index[i],output+index[i+1]);
        else  //descending
            std::sort(output+index[i],output+index[i+1], std::greater<double>());

        // Old algorithm (Slower)
        //qsort(output+index[i], index[i+1]-index[i], sizeof(double), CmpDbl);
    }

    double middle = omp_get_wtime();

    // Merge the sorted blocks
    if( numThreads > 1 )
        arraymerge(output, total, index, numThreads, dir);

    double end = omp_get_wtime();

    fprintf(stderr,"sort time = %g s, where %g s was spent on merging\n", end-start, end-middle);

    delete[] index;
}


void quick_sort_index(size_t * index, double * input, int dir, size_t total){

#ifdef _WIN32
    long long i;
#else
	size_t i;
#endif

    std::pair<double,size_t> * temp = new std::pair<double,size_t>[total];
    
    double startcopy = omp_get_wtime();
    for( i=0; i<total; i++)
    {
        temp[i].first=input[i];
        temp[i].second=i;
    }
    double endcopy = omp_get_wtime();
    fprintf(stderr,"\ncopy time = %g s\n", endcopy-startcopy);

    // Lets speed things up
    //memcpy(output,input,total*sizeof(double));

    int numThreads = omp_get_max_threads();
    size_t * indlist = new size_t[numThreads+1];
    for(i=0; i<numThreads; i++) 
        indlist[i] = i*total/numThreads; 
    indlist[numThreads] = total;

    double start = omp_get_wtime();
    // Sort each individual thread/process
#pragma omp parallel for private(i)
    for(i=0; i<numThreads; i++)
    {
        if (dir>=1) //ascending
            std::sort(temp+indlist[i],temp+indlist[i+1]);//, &CmpDblInd);
        else  //descending
            std::sort(temp+indlist[i],temp+indlist[i+1], std::greater<std::pair<double,size_t> >());

        // Old algorithm (Slower)
        //qsort(output+index[i], index[i+1]-index[i], sizeof(double), CmpDbl);
    }

    double middle = omp_get_wtime();

    // Merge the sorted blocks
    if( numThreads > 1 )
        arraymerge_index(temp, total, index, numThreads, dir);

    double end = omp_get_wtime();

    fprintf(stderr,"sort time = %g s, where %g s was spent on merging\n", end-start, end-middle);

    for( i=0; i<total; i++)
    {
        index[i] = temp[i].second;
    }
    //memcpy(index,temp,total*sizeof(size_t));
    double endfinal = omp_get_wtime();
    fprintf(stderr,"copy time = %g s\n", endfinal-end);

    delete[] indlist;
    delete[] temp;
}

void arraymerge_index(std::pair<double,size_t> * a, size_t size, size_t *index, int N, int dir)
{
#ifdef _WIN32
    long long i;
#else
	size_t i;
#endif

    while( N > 1 )
    {
        for ( i = 0; i < N; i++ ) 
            index[i] = i*size/N;
        index[N]=size;

#pragma omp parallel for private(i)
        for( i=0; i<N; i+=2 )
        {
            //fprintf(stderr,"merging %d and %d, index %d and %d (up to %d)\n",i,i+1,index[i],index[i+1],index[i+2]);
            merge_index(a+index[i], a+index[i+1], index[i+1]-index[i], index[i+2]-index[i+1], dir);
            //for(int i=0; i<size; i++) fprintf(stderr,"after: %d %f\n",i,a[i]);
        }
        N /= 2;

    }

}
void merge_index(std::pair<double,size_t> * A, std::pair<double,size_t> * B, size_t m, size_t n, int dir)
{
    int i=0, j=0, k=0;
    size_t size = m+n;

    std::pair<double,size_t> *C = new std::pair<double,size_t>[size];

    if(dir>=1) //Ascending
    {
        while ( i < m && j < n )
        {
            if(A[i].first <= B[j].first) 
                C[k] = A[i++];
            else 
                C[k] = B[j++];
            k++;
        }
    }
    else //Descending
    {
        while ( i < m && j < n )
        {
            if(A[i].first >= B[j].first) 
                C[k] = A[i++];
            else 
                C[k] = B[j++];
            k++;
        }
    }
    if (i < m) 
        for (int p = i; p < m; p++,k++) 
            C[k] = A[p];
    else 
        for (int p = j; p < n; p++,k++) 
            C[k] = B[p];

    for( i = 0; i < size; i++ ) 
        A[i] = C[i];
    
    delete[] C;


}
