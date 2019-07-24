// snwpac_test.cpp
// Author: Jesus Pulido
// Simple example on how to use SnwPac Test
// 
// Usage: ./snwpac_test {pcnt} {quant} {thr_lvl}
//	  	  {pcnt} = Threshold percentage
//	  	  {quant} = level of quantization
//	  	  {thr_lvl} = Threshold by coeff hierarchy level. {pcnt} must be set to 100.

#include <stdio.h>
#include <iostream> 
//#include <math.h>
#include <complex>
#include <time.h>
#include <omp.h>

#include <snwpac.hpp>

using namespace std; 

double argv_pcnt=50; // percentage of coefficients to threshold
int argv_quant = 0;
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
        argv_quant = atoi(argv[2]);
        cout << "Input percentage is " << argv_pcnt << endl;
        cout << "Input quant is " << argv_quant << endl;
    }
    if (argc == 4)
    {
        argv_pcnt = atof(argv[1]);
        argv_quant = atoi(argv[2]);
        argv_lvl = atoi(argv[3]);
        if (argv_pcnt == 100)
        {
            cout << "Input level is " << argv_lvl << endl;
        }
        else {
            cout << "Input percentage is " << argv_pcnt << endl;
        }
        cout << "Input quant is " << argv_quant << endl;
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
                input3d[cnt] = i + j + k;// + (rand() / RAND_MAX);
				cnt++;
			}

    size_t ogSize = total * sizeof(input3d[0]);
    cout << "SP: Original size: " << ogSize << endl;

	// Set compression parameters
	int args[13] = { 404, 0, 128+argv_quant, 0,
					dims[0], dims[1], dims[2], 
					dims[0], dims[1], dims[2], 
					sizeof(input3d[0]), argv_pcnt, 0 };
	// -------- Parameters ----------
	// { wave_type:404, chunk_level:0, region+(compression_type (+/-) quantization), padding,
	//	local_dimx, local_dimy, local_dimz,
	//	global_dimx, global_dimy, global_dimz,
	//  value_size, pcnt_threshold, level_threshold }

	// Declare SP instance, enable debugging output
	snwpac::snwpac sp(args,true);
    sp.printParams();

	// Allocate memory for compression
	void * compressed;
	compressed = std::malloc(ogSize);

	// Compress
	size_t cmpSize = sp.compress(input3d, sizeof(input3d[0]), compressed);
	cout << "SP: Compressed size: " << cmpSize << endl;

	// Verify compressed header
	sp.printHeader(compressed);

	// Allocate memory for decompression
	void * decompressed;
	decompressed = std::malloc(ogSize);

	// Decompress
	size_t dcmpSize = sp.decompress(compressed, decompressed);
	cout << "SP: Decompressed size: " << dcmpSize << endl;
	
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
	double max_rel_err = 0;
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
		double rel = abs(output3d[i] - input3d[i]);
		if (std::abs(input3d[i]) > 1)
		{
			rel /= input3d[i];
		}
		if (rel > max_rel_err)
			max_rel_err = rel;

        if (i < 10)
            std::cout << "id: " << i << " og: " << input3d[i] << " comp: " << output3d[i] << std::endl;
    }
    mse /= (total);

    cout.precision(10);
    double PSNR = 10*log10(pow(max,2.0)/(mse));
    cout << "DIFF: " << tot_diff << endl;
	cout << "MSE: " << mse << endl;
	cout << "PSNR: " << PSNR << endl;
	cout << "Max Rel Err: " << max_rel_err << endl;
	cout << "Energy is " << tot_en << endl;
	

    return 0;
}
