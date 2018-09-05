
#include <time.h>
#include <math.h>
#include <omp.h>

#include <iostream>
#include <iomanip>

#include <gsl/gsl_wavelet.h>

#define ELEMENT(a,stride,i) ((a)[(stride)*(i)])
#define ELEMENT3(a,stride,i,j,k) ((a)[(i)+(stride)*(j)+(stride)*(stride)*(k)])
#define ELEMENT3D(a,i,j,k) ((a)[(i)+(j)+(k)])
#define ELEMENT3R(a,stride1,stride2,i,j,k) ((a)[(i)+(stride1)*(j)+(stride1)*(stride2)*(k)])

template <typename T>
T *** arrayto3D(T * data, size_t x, size_t y, size_t z){
    T *** output = new T **[x];
    for(size_t i=0; i < x; i++)
    {
        output[i] = new T* [y];
        for(size_t j=0; j<y; j++)
            output[i][j]= new T[z];
    }
    size_t count = 0;
    for(size_t k=0; k<z; k++)
        for(size_t j=0; j<y; j++)
            for(size_t i=0; i<x;i++)
                { output[i][j][k]=data[count]; count++;}
                
    return output;
}
template <typename T>
T ** arrayto2D(T * data, size_t x, size_t y){
    T ** output = new T *[x];
    for(size_t i=0; i < x; i++)
    {
        output[i] = new T [y];
    }
    size_t count = 0;
    for(size_t j=0; j<y; j++)
        for(size_t i=0; i<x; i++)
            { output[i][j]=data[count]; count++;}
                
    return output;
}

template <typename T>
T * arrayto1D(T *** data, size_t x, size_t y, size_t z){
    T * output = new T [x*y*z];
    size_t count = 0;
    for(size_t k=0; k<z; k++)
        for(size_t j=0; j<y; j++)
            for(size_t i=0; i<x;i++)
            { output[count]=data[i][j][k]; count++; } // Was i,j,k
    return output;
}

template <typename T>
T * arrayto1D(T ** data, size_t x, size_t y){
    T * output = new T [x*y];
    size_t count = 0;
    for(size_t j=0; j<y; j++)
        for(size_t i=0; i<x; i++)
        { output[count]=data[i][j]; count++; } // Was i,j
    return output;
}

// Structure: [n levels][8 cubes][xmin xmax ymin ymax zmin zmax]
inline int *** getCoeffDims (int dimx, int dimy, int dimz)
{
    // Based on dimx for now
    int levels=0;
    int curdim=dimx;
    while( curdim >= 2 )
    {
        levels+=1;
        curdim = curdim/2; // Revisit this in order to add non/2 support
    }

    curdim=dimx;
    int *** cdims = new int **[levels];
    for(int i=levels-1; i>=0;i--)
    {
        cdims[i] = new int *[8];
        for(int c=0; c<8 ;c++)
        {
            // c['LLL','HLL','LHL','HHL','LLH','HLH','LHH','HHH']
            cdims[i][c] = new int[6];
        }
        cdims[i][0][0] = 0; // xmin
        cdims[i][0][1] = (curdim/2)-1; // xmax
        cdims[i][0][2] = 0; // ymin
        cdims[i][0][3] = (curdim/2)-1; // ymax
        cdims[i][0][4] = 0; // zmin
        cdims[i][0][5] = (curdim/2)-1; // zmax

        cdims[i][1][0] = curdim/2; // xmin
        cdims[i][1][1] = curdim-1; // xmax
        cdims[i][1][2] = 0; // ymin
        cdims[i][1][3] = (curdim/2)-1; // ymax
        cdims[i][1][4] = 0; // zmin
        cdims[i][1][5] = (curdim/2)-1; // zmax

        cdims[i][2][0] = 0; // xmin
        cdims[i][2][1] = (curdim/2)-1; // xmax
        cdims[i][2][2] = curdim/2; // ymin
        cdims[i][2][3] = curdim-1; // ymax
        cdims[i][2][4] = 0; // zmin
        cdims[i][2][5] = (curdim/2)-1; // zmax

        cdims[i][3][0] = curdim/2; // xmin
        cdims[i][3][1] = curdim-1; // xmax
        cdims[i][3][2] = curdim/2; // ymin
        cdims[i][3][3] = curdim-1; // ymax
        cdims[i][3][4] = 0; // zmin
        cdims[i][3][5] = (curdim/2)-1; // zmax

        cdims[i][4][0] = 0; // xmin
        cdims[i][4][1] = (curdim/2)-1; // xmax
        cdims[i][4][2] = 0; // ymin
        cdims[i][4][3] = (curdim/2)-1; // ymax
        cdims[i][4][4] = curdim/2; // zmin
        cdims[i][4][5] = curdim-1; // zmax

        cdims[i][5][0] = curdim/2; // xmin
        cdims[i][5][1] = curdim-1; // xmax
        cdims[i][5][2] = 0; // ymin
        cdims[i][5][3] = (curdim/2)-1; // ymax
        cdims[i][5][4] = curdim/2; // zmin
        cdims[i][5][5] = curdim-1; // zmax

        cdims[i][6][0] = 0; // xmin
        cdims[i][6][1] = (curdim/2)-1; // xmax
        cdims[i][6][2] = curdim/2; // ymin
        cdims[i][6][3] = curdim-1; // ymax
        cdims[i][6][4] = curdim/2; // zmin
        cdims[i][6][5] = curdim-1; // zmax

        cdims[i][7][0] = curdim/2; // xmin
        cdims[i][7][1] = curdim-1; // xmax
        cdims[i][7][2] = curdim/2; // ymin
        cdims[i][7][3] = curdim-1; // ymax
        cdims[i][7][4] = curdim/2; // zmin
        cdims[i][7][5] = curdim-1; // zmax
        
        curdim=curdim/2;
    }

    // Debug: Print the coefficent dimensions for all levels
    for(int i=0; i<levels;i++)
    {
        printf("Level %d: \n",i);
        for(int c=0; c<8 ;c++)
        {
            printf("%d:\t %d %d \t %d %d \t %d %d\n",c,
                    cdims[i][c][0], cdims[i][c][1],
                    cdims[i][c][2], cdims[i][c][3],
                    cdims[i][c][4], cdims[i][c][5]);
        }
        printf("\n");
    }

    return cdims;
}

inline int *** getCoeffDims (int dimx, int dimy)
{
    // Based on dimx for now
    int levels=0;
    int curdim=dimx;
    while( curdim >= 2 )
    {
        levels+=1;
        curdim = curdim/2; // Revisit this in order to add non/2 support
    }

    curdim=dimx;
    int *** cdims = new int **[levels];
    for(int i=levels-1; i>=0;i--)
    {
        cdims[i] = new int *[4];
        for(int c=0; c<4 ;c++)
        {
            // c['LL','HL','LH','HH']
            cdims[i][c] = new int[6];
        }
        cdims[i][0][0] = 0; // xmin
        cdims[i][0][1] = (curdim/2)-1; // xmax
        cdims[i][0][2] = 0; // ymin
        cdims[i][0][3] = (curdim/2)-1; // ymax

        cdims[i][1][0] = curdim/2; // xmin
        cdims[i][1][1] = curdim-1; // xmax
        cdims[i][1][2] = 0; // ymin
        cdims[i][1][3] = (curdim/2)-1; // ymax

        cdims[i][2][0] = 0; // xmin
        cdims[i][2][1] = (curdim/2)-1; // xmax
        cdims[i][2][2] = curdim/2; // ymin
        cdims[i][2][3] = curdim-1; // ymax

        cdims[i][3][0] = curdim/2; // xmin
        cdims[i][3][1] = curdim-1; // xmax
        cdims[i][3][2] = curdim/2; // ymin
        cdims[i][3][3] = curdim-1; // ymax
        
        curdim=curdim/2;
    }

    // Debug: Print the coefficent dimensions for all levels
    for(int i=0; i<levels;i++)
    {
        printf("Level %d: \n",i);
        for(int c=0; c<4 ;c++)
        {
            printf("%d:\t %d %d \t %d %d \n",c,
                    cdims[i][c][0], cdims[i][c][1],
                    cdims[i][c][2], cdims[i][c][3]);
        }
        printf("\n");
    }

    return cdims;
}

static int binary_logn (const size_t n)
{
  size_t ntest;
  size_t logn = 0;
  size_t k = 1;

  while (k < n)
    {
      k *= 2;
      logn++;
    }

  ntest = (1 << logn);

  if (n != ntest)
    {
      return -1;                /* n is not a power of 2 */
    }

  return logn;
}

// This function only works for values/powers of 2
// This will recursively divide until no longer divisible by 2
// and return the count of times
inline int binary_upto_logn (const size_t n)
{
  size_t logn = 0;
  size_t k = n;

  while (k>1 && k%2==0)
  {
      k /= 2;
      logn++;
  }

  return logn;
}

inline int binary_upto_value (const size_t n, const int log_n)
{
  int ret_val=(int)n;
  for(int i=1; i<log_n; i++)
      ret_val/=2;

  return ret_val;
}


inline void gsl_wavelet_print (const gsl_wavelet * w)
{
  size_t n = w->nc;
  size_t i;

  printf ("Wavelet type: %s\n", w->type->name);
  printf("nc: %d offset: %d \n", w->nc, w->offset);
  /*double val = w->h1[0];
  long double val2 = w->h1[0];
  std::cout << "Double: " << std::setprecision(64) << val << std::endl;
  std::cout << "Long Double: " << std::setprecision(64) << val2 << std::endl;*/

  printf
    (" h1(%d):%12.8f   g1(%d):%12.8f       h2(%d):%12.8f   g2(%d):%12.8f\n",
     0, w->h1[0], 0, w->g1[0], 0, w->h2[0], 0, w->g2[0]);

  for (i = 1; i < (n < 10 ? n : 10); i++)
    {
      printf
        (" h1(%d):%12.8f   g1(%d):%12.8f       h2(%d):%12.8f   g2(%d):%12.8f\n",
         i, w->h1[i], i, w->g1[i], i, w->h2[i], i, w->g2[i]);
    }

  for (; i < n; i++)
    {
      printf
        ("h1(%d):%12.8f  g1(%d):%12.8f      h2(%d):%12.8f  g2(%d):%12.8f\n",
         i, w->h1[i], i, w->g1[i], i, w->h2[i], i, w->g2[i]);
    }

}

template <typename T> 
static void dwt_step (const gsl_wavelet * w, T *a, size_t stride, size_t n, gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
{
    double ai, ai1;
    size_t i, ii;
    size_t jf;
    size_t k;
    size_t n1, ni, nh, nmod;

    for (i = 0; i < work->n; i++)
    {
        work->scratch[i] = 0.0;
    }

    //for (i = 0; i < n; i++)
    //{
    //    // You should only see traversal in the respective axis only
    //    delete findXYZIndex(stride*i);
    //}
    nmod = w->nc * n; //std::cout << nmod << "(nmod) = " << w->nc << "(nc) * " << n << "(n)\n";
    nmod -= w->offset;            // center support 
    //std::cout << "offset: " << w->offset << std::endl;
    //std::cout << "stride: " << stride << std::endl;
    n1 = n - 1;
    nh = n >> 1; ////////printf("nh: %d N: %d\n", nh, n);

    //std::cout << "n1: " << n1 << " nh: " << nh << std::endl;
    if (dir == gsl_wavelet_forward)
    {
        for (ii = 0, i = 0; i < n; i += 2, ii++)
        {
            ///////std::cout << "ii: " << ii << " i: " << i << std::endl;
            double h = 0, g = 0;

            ni = i + nmod;
          
            for (k = 0; k < w->nc; k++)
            {
                //std::cout << jf << "(jf)= " << n1 << "(n1)&( " << ni << "(ni) + " << k << "(k) )"<< std::endl;
               
                jf = n1 & (ni + k);
                h += w->h1[k] * ELEMENT (a, stride, jf);
                g += w->g1[k] * ELEMENT (a, stride, jf);
                ////std::cout << "at: "<< jf<< " = h:" << h << " g:" <<g<< std::endl;

            }

            work->scratch[ii] += h;
            work->scratch[ii + nh] += g;
        }
        //system("PAUSE");
    }
    else
    {
        for (ii = 0, i = 0; i < n; i += 2, ii++)
        {
            ai = ELEMENT (a, stride, ii);
            ai1 = ELEMENT (a, stride, ii + nh);
            ni = i + nmod;
            for (k = 0; k < w->nc; k++)
            {
                jf = (n1 & (ni + k));
                work->scratch[jf] += (w->h2[k] * ai + w->g2[k] * ai1);
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        ELEMENT (a, stride, i) = work->scratch[i];
    }
}


static void
dwt_step_proto (const gsl_wavelet * w, double *a, size_t stride, size_t n,
          gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
    {
    double ai, ai1;
    size_t i, ii;
    size_t jf;
    size_t k;
    size_t n1, ni, nh, nmod;

    for (i = 0; i < work->n; i++)
    {
        work->scratch[i] = 0.0;
    }

    //for (i = 0; i < n; i++)
    //{
    //    // You should only see traversal in the respective axis only
    //    delete findXYZIndex(stride*i);
    //}
    nmod = w->nc * n; 
    nmod -= w->offset;            // center support 
	//std::cout << "(nmod) = " << nmod << " \t (w->nc)= " << w->nc << "\t (n)= " << n << std::endl;
    //std::cout << "offset: " << w->offset << std::endl;
    //std::cout << "stride: " << stride << std::endl;
	
    n1 = n - 1;
    nh = n >> 1; // Divide by 2 (even dataset)

    //std::cout << "n1: " << n1 << " nh: " << nh << std::endl;
    if (dir == gsl_wavelet_forward)
    {
        for (ii = 0, i = 0; i < n; i += 2, ii++)
        {
            ///////std::cout << "ii: " << ii << " i: " << i << std::endl;
            double h = 0, g = 0;

            ni = i + nmod;
          
            for (k = 0; k < w->nc; k++)
            {
                // (jf) computes the current index to apply the wavelet filter. Because
				// This assumes the data is periodic, it will loop back to zero once the
				// end boundary is reached. jf begin index advances by 2 every K Loop end
                //jf = n1 & (ni + k); // Bitwise AND operator (Special case for powers of 2)
				jf = (ni + k) % n; // Modulus (Slow, but allows odd datasets)

				//std::cout << jf << "(jf)= " << n1 << "(n1)&( " << ni << "(ni) + " << k << "(k) )"<< std::endl;
                h += w->h1[k] * ELEMENT (a, stride, jf);
                g += w->g1[k] * ELEMENT (a, stride, jf);
                ////std::cout << "at: "<< jf<< " = h:" << h << " g:" <<g<< std::endl;

            }

            work->scratch[ii] += h;
            work->scratch[ii + nh] += g;
        }
        //system("PAUSE");
    }
    else
    {
        for (ii = 0, i = 0; i < n; i += 2, ii++)
        {
            ai = ELEMENT (a, stride, ii);
            ai1 = ELEMENT (a, stride, ii + nh);
            ni = i + nmod;
            for (k = 0; k < w->nc; k++)
            {
				//jf = n1 & (ni + k); // Bitwise AND operator (Special case for powers of 2)
				jf = (ni + k) % n; // Modulus (Slow, but allows odd datasets)
                work->scratch[jf] += (w->h2[k] * ai + w->g2[k] * ai1);
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        ELEMENT (a, stride, i) = work->scratch[i];
    }
}

	// Standard version of the transform (unused)
	template <typename T>
int
wavelet2d_transform (const gsl_wavelet * w, 
                         T *data, size_t tda, size_t size1,
                         size_t size2, gsl_wavelet_direction dir,
                         gsl_wavelet_workspace * work)
{
  size_t i;

  if (size1 != size2)
    {
      GSL_ERROR ("2d dwt works only with square matrix", GSL_EINVAL);
    }

  if (work->n < size1)
    {
      GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

  if (binary_logn (size1) == -1)
    {
      GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    }

  if (size1 < 2)
    {
      return GSL_SUCCESS;
    }

  if (dir == gsl_wavelet_forward)
    {
      for (i = 0; i < size1; i++)       /* for every row j */
        {
          gsl_wavelet_transform (w, &ELEMENT(data, tda, i), 1, size1, dir, work);
        }
      for (i = 0; i < size2; i++)       /* for every column j */
        {
          gsl_wavelet_transform (w, &ELEMENT(data, 1, i), tda, size2, dir, work);
        }
    }
  else
    {
      for (i = 0; i < size2; i++)       /* for every column j */
        {
          gsl_wavelet_transform (w, &ELEMENT(data, 1, i), tda, size2, dir, work);
        }
      for (i = 0; i < size1; i++)       /* for every row j */
        {
          gsl_wavelet_transform (w, &ELEMENT(data, tda, i), 1, size1, dir, work);
        }
    }

  return GSL_SUCCESS;
}

// Non-standard 2D wavelet transform x,y
template <typename T>
int wavelet2d_nstransform (const gsl_wavelet * w, 
                           T *data, size_t tda, size_t size1,
                           size_t size2, gsl_wavelet_direction dir,
                           gsl_wavelet_workspace * work)
{
  size_t i, j;

  //if (size1 != size2)
  //  {
  //    GSL_ERROR ("2d dwt works only with square matrix", GSL_EINVAL);
  //  }

  if (work->n < size1 || work->n < size2 )
    {
      GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

  if (binary_logn (size1) == -1 || binary_logn(size2))
    {
       printf("Warning: n is not a power of 2! n= %d,%d\n",size1,size2);
       printf("       : operations will be up to= %d,%d\n",binary_upto_logn(size1),binary_upto_logn(size2));
      //GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    }

  if (size1 < 2)
    {
      return GSL_SUCCESS;
    }

  int limit_x=binary_upto_logn(size1);
  int limit_y=binary_upto_logn(size2);

  if (dir == gsl_wavelet_forward)
    {
      
      for (i = size1; i >= 2; i >>= 1)
        {
          for (j = 0; j < i; j++)       /* for every row j */
            {
              dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work);
            }
          for (j = 0; j < i; j++)       /* for every column j */
            {
              dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work);
            }
        }
    }
  else
    {
      for (i = 2; i <= size1; i <<= 1)
        {
          for (j = 0; j < i; j++)       /* for every column j */
            {
              dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work);
            }
          for (j = 0; j < i; j++)       /* for every row j */
            {
              dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work);
            }
        }
    }

  return GSL_SUCCESS;
}

// Standard transform (deprecated)
template <typename T>
int gsl_wavelet3d_transform (const gsl_wavelet * w, 
                         T *data, size_t tda, size_t size1,
                         size_t size2,size_t size3, gsl_wavelet_direction dir,
                         gsl_wavelet_workspace * work, double *result)
{
    size_t i,j,k;

    if (size1 != size2 || size3 != size2 || size1 != size3)
    {
        GSL_ERROR ("3d dwt works only with square datasets", GSL_EINVAL);
    }

    if (work->n < size1)
    {
        GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

    if (binary_logn (size1) == -1)
    {
        GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    }

    if (size1 < 2)
    {
        return GSL_SUCCESS;
    }

    T * tmp = new T[size1*size2];

    clock_t start = clock();
    if (dir == gsl_wavelet_forward)
    {
        for (i = 0; i < size1; i++)       // for every row j 
        {
            //cout << i << " ";
            size_t count=0;
            for(j=0; j<size2; j++)
            {
                for(k=0; k<size3; k++)
                {
                    //tmp[count]=ELEMENT3D(data, i,j*tda,k*tda*tda);
                    tmp[count]=ELEMENT3(data,tda,i,j,k);
                    count++;
                    //cout << i<< " "<< j<<" " << k<<" "<<i+j*tda+k*tda*tda << " " << size1*size2*size3<< "\n";
                    //cout << i*tda << " " << j*tda << " " << k*tda*tda<< "\n";
                    //gsl_wavelet_transform (w, &ELEMENT3(data, tda, i,j,k*tda), 1, size1, dir, work);
                    //gsl_wavelet_transform (w, &ELEMENT3D(data, i,j*tda,k*tda*tda), 1, size1, dir, work);
                }
            }

            //start = clock();
            //cout << "Forward " << i;
            //gsl_wavelet2d_transform(w, tmp, tda, size2, size3, dir, work);
			wavelet2d_transform(w, tmp, tda, size2, size3, dir, work);
            //cout <<  " of " << size1 << "\n";

            //secs = (clock() - start) / (double) 1000;
		    //cout << endl << "Forward Transform done. Time elapsed = " << secs << " seconds" << endl;


            count=0;
            for(j=0; j<size2; j++)
            {
                for(k=0; k<size3; k++)
                {
                    //ELEMENT3D(data, i,j*tda,k*tda*tda)=tmp[count];
                    ELEMENT3(data,tda,i,j,k)=tmp[count];
                    count++;
                }
            }

        }
        for (j = 0; j < size2; j++)       // for every column j 
        {
            //gsl_wavelet_transform (w, &ELEMENT(data, 1, i), tda, size2, dir, work);
            //cout << j << " ";
            size_t count=0;
            for(i=0; i<size1; i++)
            {
                for(k=0; k<size3; k++)
                {
                    //tmp[count]=ELEMENT3D(data, i,j*tda,k*tda*tda);
                    tmp[count]=ELEMENT3(data,tda,i,j,k);
                    count++;
                }
            }

            //gsl_wavelet2d_transform(w, tmp, tda, size1, size3, dir, work);
			wavelet2d_transform(w, tmp, tda, size1, size3, dir, work);
            count=0;
            for(i=0; i<size2; i++)
            {
                for(k=0; k<size3; k++)
                {
                    //ELEMENT3D(data, i,j*tda,k*tda*tda)=tmp[count];
                    ELEMENT3(data,tda,i,j,k)=tmp[count];
                    count++;
                }
            }

        }
        /*for (k = 0; k < size3; k++)       // for every depth j 
        {
            //gsl_wavelet_transform (w, &ELEMENT(data, 1, i), tda*tda, size3, dir, work);
            //cout << k << " ";
            size_t count=0;
            for(i=0; i<size1; i++)
            {
                for(j=0; j<size2; j++)
                {
                    tmp[count]=ELEMENT3D(data, i,j*tda,k*tda*tda);
                    count++;
                }
            }

            gsl_wavelet2d_transform(w, tmp, tda, size1, size2, dir, work);
            count=0;
            for(i=0; i<size1; i++)
            {
                for(j=0; j<size2; j++)
                {
                    ELEMENT3D(data, i,j*tda,k*tda*tda)=tmp[count];
                    count++;
                }
            }

        }*/
    }
    else
    {
        /*for (k = 0; k < size3; k++)       // for every depth j 
        {
            //gsl_wavelet_transform (w, &ELEMENT(data, 1, i), tda*tda, size3, dir, work);
            size_t count=0;
            for(i=0; i<size1; i++)
            {
                for(j=0; j<size2; j++)
                {
                    tmp[count]=ELEMENT3D(data, i,j*tda,k*tda*tda);
                    count++;
                }
            }

            gsl_wavelet2d_transform(w, tmp, tda, size1, size2, dir, work);
            count=0;
            for(i=0; i<size1; i++)
            {
                for(j=0; j<size2; j++)
                {
                    ELEMENT3D(data, i,j*tda,k*tda*tda)=tmp[count];
                    count++;
                }
            }

        }*/
        for (j = 0; j < size2; j++)       // for every column j 
        {
            //gsl_wavelet_transform (w, &ELEMENT(data, 1, i), tda, size2, dir, work);
            size_t count=0;
            for(i=0; i<size1; i++)
            {
                for(k=0; k<size3; k++)
                {
                    //tmp[count]=ELEMENT3D(data, i,j*tda,k*tda*tda);
                    tmp[count]=ELEMENT3(data,tda,i,j,k);
                    count++;
                }
            }

            //gsl_wavelet2d_transform(w, tmp, tda, size1, size3, dir, work);
			wavelet2d_transform(w, tmp, tda, size1, size3, dir, work);
            count=0;
            for(i=0; i<size2; i++)
            {
                for(k=0; k<size3; k++)
                {
                    //ELEMENT3D(data, i,j*tda,k*tda*tda)=tmp[count];
                    ELEMENT3(data,tda,i,j,k)=tmp[count];
                    count++;
                }
            }

        }
        for (i = 0; i < size1; i++)       // for every row j 
        {
            //gsl_wavelet2d_transform(w, data, tda, size2, size3, dir, work);
            size_t count=0;
            for(j=0; j<size2; j++)
            {
                for(k=0; k<size3; k++)
                {
                    //tmp[count]=ELEMENT3D(data, i,j*tda,k*tda*tda);
                    tmp[count]=ELEMENT3(data,tda,i,j,k);
                    count++;
                }
            }

            //gsl_wavelet2d_transform(w, tmp, tda, size2, size3, dir, work);
			wavelet2d_transform(w, tmp, tda, size2, size3, dir, work);
            count=0;
            for(j=0; j<size2; j++)
            {
                for(k=0; k<size3; k++)
                {
                    //ELEMENT3D(data, i,j*tda,k*tda*tda)=tmp[count];
                    ELEMENT3(data,tda,i,j,k)=tmp[count];
                    count++;
                }
            }
        }

    }
    delete tmp;

    return GSL_SUCCESS;
}


// Non-Standard transform (x/y/z per level step)
template <typename T> 
int gsl_wavelet3d_nstransform (const gsl_wavelet * w, 
                           T *data, size_t tda, size_t size1,
                           size_t size2,size_t size3, gsl_wavelet_direction dir,
                           gsl_wavelet_workspace * work)			   
{
    //size_t i, j,k;
    long long i,j,k;
    
	std::cout << "Loaded type of size: " << sizeof(T) << std::endl;
	std::cout << "Size of Data pointer: " << sizeof(data[0]) << std::endl;

    //omp_set_num_threads(1);
    

    if (size1 != size2 || size2 != size3)
    {
        //GSL_ERROR ("3d dwt works only with square datasets", GSL_EINVAL);
        printf("Warning: 3d dwt isnt a square dataset, uneven operations enabled\n");
    }

    if (work->n < size1 || work->n < size2 || work->n < size3)
    {
        GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

    if (binary_logn (size1) == -1 || binary_logn (size2) == -1 || binary_logn (size3) == -1)
    {
        //GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
        printf("Warning: n is not a power of 2! n= %d,%d,%d\n",size1,size2,size3);
        printf("       : operations will be up to= %d,%d,%d\n",binary_upto_logn(size1),binary_upto_logn(size2),binary_upto_logn(size3));
        GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    }

    if (size1 < 2)
    {
        return GSL_SUCCESS;
    }
    // Compute the multiplier amount for 2 multiplication
   /* double am=size1;
    for (i = size1; i >= 4; i >>= 1)
    {
        am/=2.0;
    }
    std::cout << am << std::endl <<std::endl;*/
    if (dir == gsl_wavelet_forward)
    {
        
        //double pm;
        //for(pm = size1; pm >=140; pm/=2)
        for (i = size1; i >= 2; i >>= 1)
        {     
            double start2;
            start2= omp_get_wtime();
            //i = ceil(pm);
            //std::cout << i << std::endl;
            /* for every x */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size1);
                    //dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work); // Single slice
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,0,j,k), 1, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            } // End pragma

            /* for every y */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size2);
                    //dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work); // Single slice
                    
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,0,k), tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            } // End pragma

            /* for every z */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size3);
                    //dwt_step (w, &ELEMENT(data, tda, j), tda*tda, i, dir, work); // Single slice
                    
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,k,0), tda*tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            } // End pragma

            /*if (i > 16) //output intermediate coefficients for visualization
            {
                char filenm[200];
                sprintf_s(filenm, "D:/rstrt/128_rstrt_twelve_bior4_4_coeffs_at_%lld.vti",i/2);
                double**** output3 = new double***[1];
                output3[0] = arrayto3D(data,size1, size2, size3);
                arr_to_vti_3d_range_scale(output3,filenm,0,size1-1, 0,size2-1,0, size3-1, 1);
                //arr_to_vti_3d_range_scale(output3,filenm,0,(int)i/2-1, 0,(int)i/2-1,0, (int)i/2-1, 1); // This saves only the coarse coeffs
                printf( "Done writing intermediate coeff file\n" );
                for(size_t i=0; i < size1; i++)
                {
                    for(size_t j=0; j<size2; j++)
                        delete output3[0][i][j];
                    delete output3[0][i];
                } delete output3[0];
            }*/

            start2 = omp_get_wtime()-start2;
            std::cout.precision(10);
        std::cout << "Forward Transform level " << i << " done. OMP Time elapsed = " << start2 << " secs" << std::endl;
        }
    }
    else
    {
        //double pm;
        for (i = 2; i <= size1; i <<= 1)
        //for(pm=140; pm<=size1; pm*=2)
        {
            //i=ceil(pm);
            //std::cout << i << "\n";
            /* for every z */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)      
                {
                    work1 = gsl_wavelet_workspace_alloc (size3);
                    //dwt_step (w, &ELEMENT(data, tda, j), tda, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,k,0), tda*tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            }

            /* for every y */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size2);
                    //dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,0,k), tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            }

            /* for every x */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size1);
                    //dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,0,j,k), 1, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            }
        }
    }

    return GSL_SUCCESS;
}


// Non-Standard transform (x/y/z per level step)
template <typename T> 
int gsl_wavelet3d_nstransform_proto (const gsl_wavelet * w, 
                           T *data, size_t tda, size_t size1,
                           size_t size2,size_t size3, gsl_wavelet_direction dir,
                           gsl_wavelet_workspace * work)			   
{
    //size_t i, j,k;
    long long i,j,k;
    
	std::cout << "Loaded type of size: " << sizeof(T) << std::endl;
	std::cout << "Size of Data pointer: " << sizeof(data[0]) << std::endl;

    //omp_set_num_threads(1);
    

    //if (size1 != size2 || size2 != size3)
    //{
        //GSL_ERROR ("3d dwt works only with square datasets", GSL_EINVAL);
    //    printf("Warning: 3d dwt isnt a square dataset, uneven operations enabled\n");
    //}


    //if (binary_logn (size1) == -1 || binary_logn (size2) == -1 || binary_logn (size3) == -1)
    //{
    //    //GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    //    printf("Warning: n is not a power of 2! n= %d,%d,%d\n",size1,size2,size3);
    //    printf("       : operations will be up to= %d,%d,%d\n",binary_upto_logn(size1),binary_upto_logn(size2),binary_upto_logn(size3));
    //    GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    //}

	//Compute sizes for three dimensions

	bool enableOdd=true;

	size_t totx=0,toty=0,totz=0;
	std::vector<int> sizesx;
	std::vector<int> sizesy;
	std::vector<int> sizesz;
	size_t divx=size1;
	size_t divy=size2;
	size_t divz=size3;
 
	//Predict output size	  
	while(dixv>1)
	{
		if(divx%2 == 0)
		{
			sizesx.insert(begin(sizesx),divx);
			divx/=2;
		}
		else
		{
			sizesx.insert(begin(sizesx),divx);
			divx+=1;
			divx/=2;
			totx+=1;
		}
	}
	while(diyv>1)
	{
		if(divy%2 == 0)
		{
			sizesy.insert(begin(sizesy),divy);
			divy/=2;
		}
		else
		{
			sizesy.insert(begin(sizesy),divy);
			divy+=1;
			divy/=2;
			toty+=1;
		}
	}
	while(divz>1)
	{
		if(divz%2 == 0)
		{
			sizesz.insert(begin(sizesz),divz);
			divz/=2;
		}
		else
		{
			sizesz.insert(begin(sizesz),divz);
			divz+=1;
			divz/=2;
			totz+=1;
		}
	}
	totx+=size1;
	toty+=size2;
	totz+=size3;



	if (work->n < totx || work->n < toty || work->n < totz)
    {
        GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

    if (size1 < 2)
    {
        return GSL_SUCCESS;
    }

    if (dir == gsl_wavelet_forward)
    {
		int xcnt = sizesx.size()-1;
		int ycnt = sizesy.size()-1;
		int zcnt = sizesz.size()-1;

		int xloc = xcnt;
		int yloc = ycnt;
		int zloc = zcnt;

        //for (i = size1; i >= 2; i >>= 1)
		while( xcnt >= 0 || ycnt >= 0 || zcnt >= 0)
        {     
            double start2;
            start2= omp_get_wtime();

			// Create new data size
		    T * datnew = new T[totx*toty*totz]();

		    //Copy data with shift (if there's data growth)------------------
		    std::copy(&data[0], &data[n], &datnew[tot-n]); 

		    size_t loc=tot-sizes[sizes.size()-1];

            /* for every x */
			if( xcnt >= 0 )
			{
				#pragma omp parallel shared(data)
				{
					gsl_wavelet_workspace * work1;
					#pragma omp for private(j,k,work1)
					for (j = 0; j < sizesy[yloc]; j++)
					{
						work1 = gsl_wavelet_workspace_alloc (sizesx[xloc]+1);
						//dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work); // Single slice
						for(k = 0; k < sizesz[zloc]; k++)
						{
							dwt_step_proto<T>(w, &ELEMENT3(data, tda,0,j,k), 1, i, dir, work1);
						}
						gsl_wavelet_workspace_free (work1);
					}
				} // End pragma
				xcnt--; xloc--;
				if(xloc<0)
					xloc=0;
			}

            /* for every y */
			if( ycnt >= 0 )
			{
				#pragma omp parallel shared(data)
				{
					gsl_wavelet_workspace * work1;
					#pragma omp for private(j,k,work1)
					for (j = 0; j < sizesx[xloc]; j++)
					{
						work1 = gsl_wavelet_workspace_alloc (sizesy[yloc]+1);
						//dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work); // Single slice
                    
						for(k = 0; k < sizesz[zloc]; k++)
						{
							dwt_step_proto<T>(w, &ELEMENT3(data, tda,j,0,k), tda, i, dir, work1);
						}
						gsl_wavelet_workspace_free (work1);
					}
				} // End pragma
				ycnt--; yloc--;
				if(yloc<0)
					yloc=0;
			}

            /* for every z */
			if( zcnt >= 0 )
			{
				#pragma omp parallel shared(data)
				{
					gsl_wavelet_workspace * work1;
					#pragma omp for private(j,k,work1)
					for (j = 0; j < sizesx[xloc]; j++)
					{
						work1 = gsl_wavelet_workspace_alloc (sizesz[zloc]+1);
						//dwt_step (w, &ELEMENT(data, tda, j), tda*tda, i, dir, work); // Single slice
                    
						for(k = 0; k < sizesy[yloc]; k++)
						{
							dwt_step_proto<T>(w, &ELEMENT3(data, tda,j,k,0), tda*tda, i, dir, work1);
						}
						gsl_wavelet_workspace_free (work1);
					}
				} // End pragma
				zcnt--; zloc--;
				if(zloc<0)
					zloc=0;
			}

            /*if (i > 16) //output intermediate coefficients for visualization
            {
                char filenm[200];
                sprintf_s(filenm, "D:/rstrt/128_rstrt_twelve_bior4_4_coeffs_at_%lld.vti",i/2);
                double**** output3 = new double***[1];
                output3[0] = arrayto3D(data,size1, size2, size3);
                arr_to_vti_3d_range_scale(output3,filenm,0,size1-1, 0,size2-1,0, size3-1, 1);
                //arr_to_vti_3d_range_scale(output3,filenm,0,(int)i/2-1, 0,(int)i/2-1,0, (int)i/2-1, 1); // This saves only the coarse coeffs
                printf( "Done writing intermediate coeff file\n" );
                for(size_t i=0; i < size1; i++)
                {
                    for(size_t j=0; j<size2; j++)
                        delete output3[0][i][j];
                    delete output3[0][i];
                } delete output3[0];
            }*/

            start2 = omp_get_wtime()-start2;
            std::cout.precision(10);
        std::cout << "Forward Transform level " << i << " done. OMP Time elapsed = " << start2 << " secs" << std::endl;
        }
    }
    else
    {
        //double pm;
        for (i = 2; i <= size1; i <<= 1)
        //for(pm=140; pm<=size1; pm*=2)
        {
            //i=ceil(pm);
            //std::cout << i << "\n";
            /* for every z */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)      
                {
                    work1 = gsl_wavelet_workspace_alloc (size3);
                    //dwt_step (w, &ELEMENT(data, tda, j), tda, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step_proto<T>(w, &ELEMENT3(data, tda,j,k,0), tda*tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            }

            /* for every y */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size2);
                    //dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step_proto<T>(w, &ELEMENT3(data, tda,j,0,k), tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            }

            /* for every x */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size1);
                    //dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step_proto<T>(w, &ELEMENT3(data, tda,0,j,k), 1, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            }
        }
    }

    return GSL_SUCCESS;
}

template <typename T> int gsl_wavelet3d_nstransform (const gsl_wavelet * w, 
                           T *data, size_t tda, size_t size1,
                           size_t size2,size_t size3, gsl_wavelet_direction dir,
                           gsl_wavelet_workspace * work, int levels)			   
{
    //size_t i, j,k;
    long long i,j,k;
   
	//omp_set_num_threads(1);
    

    if (size1 != size2 || size2 != size3)
    {
        GSL_ERROR ("3d dwt works only with square datasets", GSL_EINVAL);
        //printf("Warning: 3d dwt isnt a square dataset, uneven operations enabled\n");
		return -1;
    }

    if (work->n < size1 || work->n < size2 || work->n < size3)
    {
        GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

    if (binary_logn (size1) == -1 || binary_logn (size2) == -1 || binary_logn (size3) == -1)
    {
        //GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
        printf("Warning: n is not a power of 2! n= %d,%d,%d\n",size1,size2,size3);
        printf("       : operations will be up to= %d,%d,%d\n",binary_upto_logn(size1),binary_upto_logn(size2),binary_upto_logn(size3));
        GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
		return -1;
    }

    if (size1 < 2)
    {
        return GSL_SUCCESS;
    }

	int lvl=0;
    if (dir == gsl_wavelet_forward)
    {
        
        for (i = size1; i >= 2; i >>= 1)
        {     
            double start2;
            start2= omp_get_wtime();

			/* for every x */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size1);
                    //dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work); // Single slice
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,0,j,k), 1, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            } // End pragma

            /* for every y */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size2);
                    //dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work); // Single slice
                    
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,0,k), tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            } // End pragma

            /* for every z */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size3);
                    //dwt_step (w, &ELEMENT(data, tda, j), tda*tda, i, dir, work); // Single slice
                    
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,k,0), tda*tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            } // End pragma

            start2 = omp_get_wtime()-start2;
            std::cout.precision(10);
	        std::cout << "Forward Transform level " << i << " done. OMP Time elapsed = " << start2 << " secs" << std::endl;

			lvl++;
			if(lvl == levels)
				return GSL_SUCCESS;
        }
    }
    else
    {
		long long start=size1; 
		for(int it=0; it<levels;it++) //For the JHTDB Dataset, this is levels-1
		{
			start=start/2;
		}

        for (i = start; i <= size1; i <<= 1)
        {
            /* for every z */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)      
                {
                    work1 = gsl_wavelet_workspace_alloc (size3);
                    //dwt_step (w, &ELEMENT(data, tda, j), tda, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,k,0), tda*tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            }

            /* for every y */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size2);
                    //dwt_step (w, &ELEMENT(data, 1, j), tda, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,j,0,k), tda, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }
            }

            /* for every x */
            #pragma omp parallel shared(data)
            {
                gsl_wavelet_workspace * work1;
                #pragma omp for private(j,k,work1)
                for (j = 0; j < i; j++)
                {
                    work1 = gsl_wavelet_workspace_alloc (size1);
                    //dwt_step (w, &ELEMENT(data, tda, j), 1, i, dir, work); // 2D
                    for(k = 0; k < i; k++)
                    {
                        dwt_step<T>(w, &ELEMENT3(data, tda,0,j,k), 1, i, dir, work1);
                    }
                    gsl_wavelet_workspace_free (work1);
                }

            }
        }
    }

    return GSL_SUCCESS;
}

template <typename T>
int grow(int val, T *** data, size_t x, size_t y, size_t z)
{
    std::cout << "Growing array edges by " << val << std::endl;
    // Declare new array and copy over data to center of it
    // Create new array
    T *** _dat = new T**[x+(2*val)];
    for(int i=0; i<(x+(2*val)); i++) {
	    _dat[i] = new T*[y+(2*val)];
        for(int j=0; j<(y+(2*val)); j++) {	
		    _dat[i][j] = new T[z+(2*val)];
        }
    }
    int a,b,c;
    
    // Copy over data (assumes periodic by default)
    int edge=x+val;
    for(int i=0; i<(x+(2*val)); i++) {
        a=i;
        if(i<val) {a=(x+i);}
        if(i>=edge) {a=i-(x);}
        for(int j=0; j<(y+(2*val)); j++) {	
            b=j;
            if(j<val) {b=(y+j);}
            if(j>=edge) {b=j-(y);}
		    for(int k=0; k<(z+(2*val)); k++) {
                c=k;
                if(k<val) {c=(z+k);}
                if(k>=edge) {c=k-(z);}

                //cout << "To " << i << " " << j << " " << k << " from " << a-val << " " << b-val << " " << c-val << endl;
                _dat[i][j][k] = data[a-val][b-val][c-val];

            }
        }
    }


     // delete old
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            delete [] data[i][j]; }
        delete [] data[i];
    }
    delete [] data; data=NULL;

    // Copy over new
    data = _dat;

    return 1;
} 

template <typename T>
int shrink(int val, T *** data, size_t x, size_t y, size_t z)
{
    std::cout << "Shrinking array edges by " << val << std::endl;
    // Create new array
    T *** _dat = new T**[x-(2*val)];
    for(int i=0; i<(x-(2*val)); i++) {
	    _dat[i] = new T*[y-(2*val)];
        for(int j=0; j<(y-(2*val)); j++) {	
		    _dat[i][j] = new T[z-(2*val)];
        }
    }
    int a,b,c;
    // Copy over data

    for(int i=val; i<(x-(val)); i++) {
        for(int j=val; j<(y-(val)); j++) {	
		    for(int k=val; k<(z-(val)); k++) {
                _dat[i-val][j-val][k-val] = data[i][j][k];
            }
        }
    }

     // delete old
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            delete [] data[i][j]; }
        delete [] data[i];
    }
    delete [] data; data=NULL;

    // Update to new
    data = _dat;

    return 1;
}

template <typename T>
T * getCutout3d(T * in, int * codims, int * dims)
{
    T * out = NULL;

    if(!in)
        return NULL;

    int dx=codims[1]-codims[0]+1;
    int dy=codims[3]-codims[2]+1;
    int dz=codims[5]-codims[4]+1;

    int cnt=0;
    int ind=0;
    out = new T[dx*dy*dz];
    for(int z=codims[4];z<=codims[5];z++)
        for(int y=codims[2];y<=codims[3];y++)
            for(int x=codims[0];x<=codims[1];x++)
            {
                ind = (z*dims[0]*dims[1])+(y*dims[0])+x;
                out[cnt]=in[ind];
                cnt++;
            }

    return out;
}

template <typename T>
T * getCutout2d(T * in, int * codims, int * dims)
{
    T * out = NULL;

    if(!in)
        return NULL;

    int dx=codims[1]-codims[0]+1;
    int dy=codims[3]-codims[2]+1;

    int cnt=0;
    int ind=0;
    out = new T[dx*dy];
    for(int y=codims[2];y<=codims[3];y++)
        for(int x=codims[0];x<=codims[1];x++)
        {
            ind = (y*dims[0])+x;
            out[cnt]=in[ind];
            cnt++;
        }

    return out;
}

template <typename T>
void setCutout3d(T * in, T * cutout, int * codims, int * dims)
{
    if(!in)
        return;
    if(!cutout)
        return;

    //int dx=codims[1]-codims[0]+1;
    //int dy=codims[3]-codims[2]+1;
    //int dz=codims[5]-codims[4]+1;

    size_t cnt=0;
    size_t ind=0;
    for(int z=codims[4];z<=codims[5];z++)
        for(int y=codims[2];y<=codims[3];y++)
            for(int x=codims[0];x<=codims[1];x++)
            {
                ind = (z*dims[0]*dims[1])+(y*dims[0])+x;
                in[ind]=cutout[cnt];
                cnt++;
            }

    return;
}

template <typename T>
void setCutout2d(T * in, T * cutout, int * codims, int * dims)
{
    if(!in)
        return;
    if(!cutout)
        return;

    //int dx=codims[1]-codims[0]+1;
    //int dy=codims[3]-codims[2]+1;

    size_t cnt=0;
    size_t ind=0;
    for(int y=codims[2];y<=codims[3];y++)
        for(int x=codims[0];x<=codims[1];x++)
        {
            ind = (y*dims[0])+x;
            in[ind]=cutout[cnt];
            cnt++;
        }

    return;
}

template <typename T>
int wavelet_transform (const gsl_wavelet * w, 
                       T *&data, size_t stride, size_t n,
                       gsl_wavelet_direction dir, 
                       gsl_wavelet_workspace * work)
{
  size_t i;
  size_t tot=0;
  std::vector<int> sizes;

  bool enableOdd=false;
  if( binary_logn (n) == -1 )
  {  
	enableOdd=true;

	//Predict output size	  
	size_t div=n;
	//std::cout << "Hierarchy: ";
	while(div>1)
	{
		if(div%2 == 0)
		{
			//std::cout << div << " ";
			sizes.insert(begin(sizes),div);
			div/=2;
		}
		else
		{
			//std::cout << " (" << div << "+1) ";
			sizes.insert(begin(sizes),div);
			div+=1;
			div/=2;
			tot+=1;
		}
	}

	tot+=n;
	//std::cout << " SZ_out: " << tot << std::endl;
  }

  if (work->n < tot)
  {
      GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
  }

  if (n < 2)
  {
      return GSL_SUCCESS;
  }
  

  if(!enableOdd)
  {
	  //int limit_x=binary_upto_logn(n);
	  if (dir == gsl_wavelet_forward)
	  {
		  int steps=0;
		  //for (i = n; i >= 2 && steps < limit_x; i /= 2)
		  for (i = n; i >= 2 ; i /= 2)
		  {
			  //printf("Starting step: sz %d of %d\n",i,(int)n);
			  dwt_step (w, data, stride, i, dir, work);
			  steps++;
		  }

		  //printf("Completed %d steps\n", steps);
	  }
	  else
	  {
		  int st=2;
		  //if(binary_logn(n) == -1)
		  //	  st=binary_upto_value(n,limit_x);

		  for (i = st; i <= n; i *= 2)
		  {
			  //printf("Starting step: sz %d of %d\n",i,(int)n);
			  dwt_step (w, data, stride, i, dir, work);
		  }
	  }
  }
  else
  {

	  if (dir == gsl_wavelet_forward)
	  {
		  // Create new data size
		  T * datnew = new T[tot]();

		  //Copy data with shift (if there's data growth)
		  std::copy(&data[0], &data[n], &datnew[tot-n]); 

		  size_t loc=tot-sizes[sizes.size()-1];
		  // Declare local int i, otherwize size_t will underflow to MAX
		  for(int i=sizes.size()-1; i>=0; i--) //Rather then iterate by sizes, iterate by position of sizes array
		  {
			  //printf("Starting step:  %d of %d : Size %d\n",i+1,(int)sizes.size(),sizes[i]);
			  //std::cout << "loc: " << loc << std::endl;

			  // Detect padding and duplicate
			  if(sizes[i]%2!=0)
			  {
				//Simple shift
				loc=loc-1;

				//Duplicate coarse
				//Simple copy (odd -> even)
				datnew[loc] = datnew[loc+1]; 

				dwt_step_proto (w, &datnew[loc], stride, (sizes[i]+1), dir, work);
			  }
			  else
			  {
				// Use existing location for next size
				dwt_step_proto (w, &datnew[loc], stride, sizes[i], dir, work);
			  }
		  }

		  /*if(loc != 0)
		  {  std::cout << "Decomposition Error: Non-zero end location: " << loc << "\n";
		     return 0;
		  }
		  */
		  delete [] data;
		  data = datnew;

	  }
	  else //(gsl_wavelet_backward)
	  {
		  //gsl_wavelet_workspace *work;
		  size_t loc=0;
		  for(int i=0; i<sizes.size(); i++) //Rather then iterate by sizes, iterate by position of sizes array
		  {
			  //std::cout << std::endl;
			  //printf("Starting step:  %d of %d : Size %d\n",i+1,(int)sizes.size(),sizes[i]);
			  //std::cout << "loc: " << loc << std::endl;
			  //std::cout << std::endl;

			  if(sizes[i]%2!=0)
			  {
				//work = gsl_wavelet_workspace_alloc (sizes[i]+1);
				dwt_step_proto (w, &data[loc], stride, (sizes[i]+1), dir, work);
				//gsl_wavelet_workspace_free (work);
				data[loc]=0; // optional

				loc+=1;

			  }
			  else
			  {

				  //work = gsl_wavelet_workspace_alloc (sizes[i]);
				  dwt_step_proto (w, &data[loc], stride, sizes[i], dir, work);
				  //gsl_wavelet_workspace_free(work);
			  }


		  }

		  // Create new data size
		  T* datnew = new T[n]();

		  //Copy data removing shifts (if there's data growth)
		  std::copy(&data[tot-n], &data[tot], &datnew[0]);

		  delete [] data;
		  data = datnew;
	  }
  }

  return GSL_SUCCESS;
}