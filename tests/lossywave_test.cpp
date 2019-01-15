// lossywave_test.cpp
// Author: Jesus Pulido
// Simple example on how to use LossyWave
// 
// Usage: ./lossywave_test {pcnt}
//	  	  {pcnt} = Threshold percentage

#include <stdio.h>
#include <iostream> 
//#include <math.h>
#include <complex>
#include <time.h>
#include <omp.h>

#include <lossywave.hpp>

using namespace std; 

double argv_pcnt=50;
int argv_lvl=0;

int main (int argc, char **argv)
{
	
    if(argc == 2)
    {
        argv_pcnt = atof(argv[1]);
        cout << "Input percentage is " << argv_pcnt << endl;
    }
    if(argc == 3)
    {
        argv_pcnt = atof(argv[1]);
        argv_lvl = atoi(argv[2]);
        if(argv_pcnt == 100)
        {
            cout << "Input level is " << argv_lvl << endl;
        }else{
            cout << "Input percentage is " << argv_pcnt << endl;
        }
    }
	

	// Create dummy dataset
    int * dims = new int[3];
    dims[0]=32; dims[1]=32; dims[2]=32;
	size_t total = dims[0] * dims[1] * dims[2];
    float * input3d = new float[total];

	size_t cnt=0;
    for(int i=0; i< dims[0]; i++)
        for(int j=0; j< dims[1]; j++)
			for (int k = 0; k < dims[2]; k++)
			{
				input3d[cnt] = i + j + k;
				cnt++;
			}

    size_t ogSize = total * sizeof(input3d[0]);
    cout << "LW: Original size: " << ogSize << endl;

	// Set compression parameters
	int args[13] = { 404, 0, 128, 0, 
					dims[0], dims[1], dims[2], 
					dims[0], dims[1], dims[2], 
					sizeof(input3d[0]), argv_pcnt, 0 };
	// -------- Parameters ----------
	// { wave_type:404, chunk_level:0, region+compression_type , padding,
	//	local_dimx, local_dimy, local_dimz,
	//	global_dimx, global_dimy, global_dimz,
	//  value_size, pcnt_threshold, level_threshold }

	// Declare LW instance, enable debugging output
	lossywave::lossywave lw(args,true);
    lw.printParams();

	// Allocate memory for compression
	void * compressed;
	compressed = std::malloc(ogSize);

	// Compress
	size_t cmpSize = lw.compress(input3d, sizeof(input3d[0]), compressed);
	cout << "LW: Compressed size: " << cmpSize << endl;

	// Allocate memory for decompression
	void * decompressed;
	decompressed = std::malloc(ogSize);

	// Decompress
	size_t dcmpSize = lw.decompress(compressed, decompressed);
	cout << "LW: Decompressed size: " << dcmpSize << endl;
	
    // Check if data is valid
    if (ogSize == dcmpSize)
        std::cout << "Data Verified!" << std::endl;
    else
        std::cout << "Invalid Data! Compression Error." << std::endl;
    
	// Compare compressed vs original
	float * output3d = static_cast<float *>(decompressed);

    cout << "--Comparison--" << endl;
    double tot_en=0;
    double tot_diff=0;
    double mse=0;
    double max=-999;
    for(size_t i=0; i<total; i++)
    {
	    if(output3d[i]>max)
	        max=output3d[i];
	    if(input3d[i]>max)
	        max=input3d[i];
	                        
	    complex<double> value(output3d[i], 0.0);
	    tot_en += norm(value);
	    mse+=(pow((double)output3d[i]-input3d[i],(double)2.0));
	    tot_diff += abs(output3d[i]-input3d[i]);

    }
    mse /= (total);

    cout.precision(10);
    double PSNR = 10*log10(pow(max,2.0)/(mse));
    cout << "DIFF: " << tot_diff << endl;
	cout << "MSE: " << mse << endl;
	cout << "PSNR: " << PSNR << endl;
	cout << "Energy is " << tot_en << endl;

    return 0;
}
