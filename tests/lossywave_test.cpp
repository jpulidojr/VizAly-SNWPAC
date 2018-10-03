#include <stdio.h>
#include <iostream> 
//#include <math.h>
#include <complex>
#include <time.h>
#include <omp.h>

#include <lossywave.hpp>

using namespace std; 

double argv_pcnt=100;
int argv_lvl=0;

int main (int argc, char **argv)
{
	/* // legacy args input 
    bool useLevel=false;
	bool save_coeff=true;
	bool save_bin=false;
	bool save_vti=false;

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
            useLevel=true;
            cout << "Input level is " << argv_lvl << endl;
        }else{
            cout << "Input percentage is " << argv_pcnt << endl;
        }
    }
	*/

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
				input3d[cnt] = 1;// i + j + k;
				cnt++;
			}

	// Set compression parameters
	int args[13] = { 303, 0, 128, 0, 
					dims[0], dims[1], dims[2], 
					dims[0], dims[1], dims[2], 
					sizeof(input3d[0]), 100, 0 };
	// -------- Parameters ----------
	// { wave_type:303, chunk_level:0, region+compression_type , padding,
	//	local_dimx, local_dimy, local_dimz,
	//	global_dimx, global_dimy, global_dimz,
	//  value_size, pcnt_threshold, level_threshold }

	// Declare LW instance
	lossywave::lossywave lw(args);

	// Allocate memory for compression
	void * compressed;
	compressed = std::malloc(total * sizeof(float));

	// Compress
	size_t cmpSize = lw.compress(input3d, sizeof(float), compressed);
	cout << "Compressed size: " << cmpSize << endl;

	// Allocate memory for decompression
	void * decompressed;
	decompressed = std::malloc(total * sizeof(float));

	// Decompress
	size_t dcmpSize = lw.decompress(compressed, decompressed);
	cout << "Decompressed size: " << dcmpSize << endl;
	
	// Compare compressed vs original
	float * output3d = static_cast<float *>(decompressed);

    cout << "Beginning comparison" << endl;
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
