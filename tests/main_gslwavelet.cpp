#include <stdio.h>
#include <iostream> 
//#include <math.h>
#include <complex>
#include <time.h>
#include <omp.h>

#include <gsl/gsl_sort.h>
#include <algorithm>

#include "reader.h"
#include "writer.h"
#include "wavelet.h"

#include "quick_sort.h"

using namespace std; 

double argv_pcnt=100;
int argv_lvl=0;

int main (int argc, char **argv)
{
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

    size_t i=0; //int nc = 1;
    int fields = 1;
    int * pDims = new int[3];
    pDims[0]=512; pDims[1]=512; pDims[2]=512;
    int stride = pDims[0];
    int total = pDims[0]*pDims[1]*pDims[2];
    double *data = NULL;// = new double[total];
    double *abscoeff = NULL;// = new double[total];
    size_t *p = NULL; // Used for sorting
    
    //FILE * f;
    gsl_wavelet *w;
    gsl_wavelet_workspace *work;
    
    size_t type = 303;
    //w = gsl_wavelet_alloc( NULL, 301 );
    w = gsl_wavelet_alloc (gsl_wavelet_bspline, type);
    //w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
    //w = gsl_wavelet_alloc (gsl_wavelet_haar, 2);
    gsl_wavelet_print ( w );


    cout << "Total threads in this machine: " << omp_get_max_threads() << endl;

    // Set to largest direction
    int val=pDims[0];
    if(val<pDims[1]) val=pDims[1];
    if(val<pDims[2]) val=pDims[2];
    
	work = gsl_wavelet_workspace_alloc( val ); //Normal during 3d Decomp
    //work = gsl_wavelet_workspace_alloc (total);//(total); //Normal during 1D decomp

    double secs=0;

    clock_t start = clock();
    cout << "Reading in file\n";
    // read in data (or generate)
    double **** input3d = new double***[fields];
    for(int fld=0; fld< fields; fld++)
    {
    //double *** input3d = new double**[pDims[0]];
        input3d[fld] = new double**[pDims[0]];
        for(int l=0; l<pDims[0]; l++)
        {
            input3d[fld][l] = new double*[pDims[1]];
            for(int j=0; j<pDims[1]; j++)
                input3d[fld][l][j] = new double[pDims[2]];
        }
        int res=0;
        //res = binto3DarrayMPI(input3d[fld], const_cast<char *>("D:/rstrt/rstrt.0009.bin"),0,pDims[0]-1,0,pDims[1]-1,0,pDims[2]-1, fld );
        
       
        //res = dbto3Darray(input3d[fld], 0,pDims[0]-1,0,pDims[1]-1,0,pDims[2]-1, fld );


        res = binto3DarrayMPI(input3d[fld], const_cast<char *>("D:/data/rstrt.twelve.bin"),0,pDims[0]-1,0,pDims[1]-1,0,pDims[2]-1, fld );
        //res = binto3DarrayMPI(input3d[fld], const_cast<char *>("D:/rstrt/8x8.bin"),0,pDims[0]-1,0,pDims[1]-1,0,pDims[2]-1, fld );
        if(!res)
            return 0;
    }
    secs = (clock() - start) / (double) 1000;
    cout << endl << "Reading file complete. Time elapsed = " << secs << " secs" << endl;

    /*cout << "Testing byte search functions\n";
    double add;
    for (int z=0;z<pDims[2];z++)
        for (int y=0;y<pDims[1];y++)
            for (int x=0;x<pDims[0];x++)
            {
                add = findByteAddress(x, y, z, 0);
                //cout << "Address for 128,128,128 is " << add << endl;
                int * loc = findXYZAddress(add,0);
                //cout << "Location for Address " << add << " is " << loc[0] << "," << loc[1] << "," <<loc[2] << endl;
                if (x != loc[0] || y != loc[1] || z != loc[2])
                    cout << "Test failed at location: " << x << "," << y << "," << z <<": Wrong-> " << loc[0] << "," << loc[1] << "," <<loc[2] << endl;
                delete loc;
            }
    cout << "Test successful!\n";
    return 0;*/

    // Overwrites read and creates dummy dataset
    /*for(int fld=0; fld< fields; fld++)
        for(int i=0; i<pDims[0]; i++)
            for(int j=0; j<pDims[1]; j++)
                for(int k=0; k<pDims[2];k++)
                   input3d[fld][i][j][k]=0;//fld+i+j+k;
	*/
    double start2=0;

	// Loops for 0->100% coefficient thresholds, if needed
    for(int nc=(int)argv_pcnt; nc <=argv_pcnt; nc++)
    {
        double **** output3d = new double***[fields];
        
        for(int fld=0; fld < fields; fld++)
        {
			// 3D to 1D data array conversion
            data = arrayto1D(input3d[fld],pDims[0],pDims[1],pDims[2]);

            //gsl_wavelet2d_transform_forward (w, data, 1, 16,16, work);

			
			// --- Testing odd wavelet method----
			//Make a copy of data array
			/*int wn = total;
			double * dat2= new double[total];
			for(int i=0; i<wn; i++)
				dat2[i]=data[i];*/

			cout << "Beginning Transform.. \n";
			secs=0;
			start = clock();
			start2= omp_get_wtime();
			gsl_wavelet3d_nstransform (w, data, stride, pDims[0], pDims[1], pDims[2], gsl_wavelet_forward, work);
			//gsl_wavelet_transform (w, dat2, 1, wn, gsl_wavelet_forward, work);
			start2 = omp_get_wtime()-start2;
			secs = (clock() - start) / (double) 1000;
			cout << endl << "Forward Transform field " << fld+1 << " done. Time elapsed = " << secs << " secs" << endl;
			cout << "Forward Transform field " << fld+1 << " done. OMP Time elapsed = " << start2 << " secs" << endl;

			// We just did the original 1-D transform, now lets test our own 1-D version with odd support
			/*delete [] dat2;
			dat2= new double[total];
			for(int i=0; i<wn; i++)
				dat2[i]=data[i];

			secs=0;
			start = clock();
			start2= omp_get_wtime();
			//gsl_wavelet3d_nstransform (w, data, stride, pDims[0], pDims[1], pDims[2], gsl_wavelet_forward, work);
			wavelet_transform (w, dat2, 1, wn, gsl_wavelet_forward, work);
			start2 = omp_get_wtime()-start2;
			secs = (clock() - start) / (double) 1000;
			cout << endl << "Forward Transform field " << fld+1 << " done. Time elapsed = " << secs << " secs" << endl;
			cout << "Forward Transform field " << fld+1 << " done. OMP Time elapsed = " << start2 << " secs" << endl;
			
			return 0;

			// ---- Testing odd wavelet method
			//Verify integrity
			/*wavelet_transform (w, dat2, 1, wn, gsl_wavelet_backward, work);

			double verify=0;
			for(int j=0;j<wn;j++)
			{
				verify += abs(dat2[j]-data[j]);
			}
			std::cout << std::endl << "Wn: " << wn << " Verification Status: " << verify << std::endl;

			delete [] dat2;*/

            /*wavelet_transform (w, data, 1, total, gsl_wavelet_backward, work);

            double err=0.0;
            cout.precision(32);
            for(int i=0; i<total; i++){
                err+=abs(dat2[i]-data[i]);
                cout << "Mismatch: " << i << " og: " << dat2[i] << " vs rec: " << data[i] << endl;
            }
            cout << "Accumilated error: " << err << endl;
            return 0;*/

			//This begins %-based coefficient thresholding
            if(!useLevel)
            {
                cout << "Beginning sort\n";
                abscoeff = new double[total];
                p = new size_t[total];

                for (i = 0; i < total; i++)
                {
                    abscoeff[i] = fabs (data[i]);
                }

                start = clock();
                
                //gsl_sort_index (p, abscoeff, 1, total); // Serial, old and slow
                quick_sort_index (p, abscoeff, 1, total); // Fast, MPI Based sort
                secs = (clock() - start) / (double) 1000;
                //cout << "Ending serial sort. Time elapsed = " << secs << " secs" << endl;

                /*size_t *p2;
                p2 = new size_t[total];
                start = clock();
                
                secs = (clock() - start) / (double) 1000;*/
                cout << "Ending parallel sort. Time elapsed = " << secs << " secs" << endl;

                delete [] abscoeff;

                /*for(size_t cnt = 0; cnt < total; cnt++)
                {
                    
                    if (p[cnt] != p2[cnt])
                    {
                        cout << "At: " << cnt << " we have p= " << p[cnt] << " as coeff " << data[p[cnt]] << endl;
                        cout << "At: " << cnt << " we have p2= " << p2[cnt] << " as coeff " << data[p2[cnt]] << endl;
                        cout << "Sort mismatch, exiting!" << endl;
                        
                        return 0;
                    }

                }
                cout << "Confirmed correct sorting" << endl;*/
            }
            
			/*
            // Save all coefficients to disk (no encoding)
            int * args = coeff_args( type, 0, 0, 0, pDims[0], pDims[1], pDims[2], pDims[0], pDims[1], pDims[2],sizeof(data[0]));

			// Save all coefficients to disk (RLE encoding) (reg=64+region(0,1,2....))
			// SEE SECTION BELOW: First threshold!!!

            //save_coefficients(data, "D:/Data/rstrt_280.0012.dens.b4.4.coeff",total); //Deprecated
            save_coefficients_md(data, "D:/Data/rstrt_test_0012.dens.b4.4.coeff",args);

            //free_args(args);
            // Save sort to disk
            //save_sort(p, "D:/rstrt/rstrt_256.0012.dens.b4.4.sort",total);
            

            // read coefficients from disk
            start = clock();
            //data = read_coefficients("D:/rstrt/rstrt.0012.dens.b3.1.coeffs",total);
            double * data2 = read_coefficients_md<double>("D:/Data/rstrt_test_0012.dens.b4.4.coeff",args);
            total = (size_t)args[4]*(size_t)args[5]*(size_t)args[6];
             secs = (clock() - start) / (double) 1000;
            cout << "Ending coef read. Time elapsed = " << secs << " secs" << endl;

            //free_args(args);

			for(size_t sz=0; sz<total; sz++)
                    if(data[sz] != data2[sz])
                        std::cout << "Compression mismatch: " << sz << ", orig: " << data[sz] << " comp " << data2[sz] << std::endl;
                std::cout << "COMPLETE! Sanity check passed\n";

			delete [] data2;
			return 0;

            // read sort from disk
            /*start = clock();
            p = read_sort("D:/rstrt/rstrt.0012.dens.b3.1.sort",total);
            secs = (clock() - start) / (double) 1000;
            cout << "Ending sort read. Time elapsed = " << secs << " secs" << endl;*/

            // Compute available levels (possibly used in functions below)
            int max_levels=0;
            int curdim = pDims[0];
            while(curdim >= 2)
            {
                max_levels+=1;
                curdim/=2;
            }

            // Write individual coefficient blocks to disk
            /*int *** cdims = getCoeffDims(pDims[0],pDims[1],pDims[2]);

            //cout << "Thresholding by level: " << argv_lvl << " of " << max_levels << "\n";
            int curdims[3] = {pDims[0],pDims[1],pDims[2]};
            for(int lvl=1; lvl < max_levels; lvl++)
            {
                int lvlsz=cdims[lvl][0][1]+1; // xmax+1 assuming 3D cube
                lvlsz = lvlsz*lvlsz*lvlsz;
                for(int i=0; i<1;i++)
                {
                    double * LLL = getCutout3d(data,cdims[lvl][i], curdims);
                    int * args = coeff_args( type, lvl, i, 0, cdims[lvl][0][1]+1, cdims[lvl][0][3]+1, cdims[lvl][0][5]+1, pDims[0], pDims[1], pDims[2],sizeof(double));
                    char filenm[200];

                    sprintf_s( filenm, "D:/rstrt/rstrt_twelve_bior4_4_lvl_%d_rg_%d.coeff",lvl,i);
                    save_coefficients_md(LLL,filenm,args);

                    free_args(args);
                    delete [] LLL;

                }
            }
            
            delete [] cdims;*/

			// Use %-based coefficient thresholding
            if(!useLevel)
            {
                double num_coeff = (total*((double)argv_pcnt/100));
                cout << "Thresholding " << argv_pcnt << "%, " << num_coeff << endl;
                start = clock();
                for (i = 0; (i + num_coeff) < total; i++)
                    data[p[i]] = 0;

                secs = (clock() - start) / (double) 1000;
                cout << "Ending Threshold. Time elapsed = " << secs << " secs" << endl;
                delete [] p;

                // Experimental ----------------
                // Count contiguous sectors of 0s for non-uniform data compression
				// Based on run-length encoding
				float sizeof_val = sizeof(data[0]);
                cout << "Full coefficient file size (double): " << total*sizeof_val << " bytes or " << (total*sizeof_val)/(1024.0*1024.0) << " MB\n";
                double nonzero = 0.0;
                double zero_contiguous = 0.0;
                int times_contiguous = 0;
				

                int cont_block_cnt=0; // DEBUG
                double largest_cont_block = 0; // DEBUG
                int cnt_exceed = 0;
                for (i = 0; i < total; i++)
                {
                    
                    if (data[i] == 0)
                    {    
                        cont_block_cnt=0; // DEBUG
                        times_contiguous++;
                        while (data[i]==0 && i < total)
                        {
                            zero_contiguous += sizeof_val;//8.0;
                            i++;
                            cont_block_cnt++; // DEBUG
                        }
                        if(cont_block_cnt>largest_cont_block) // DEBUG
                            largest_cont_block=cont_block_cnt; // DEBUG
                        while(cont_block_cnt>65535)
                        {
                            cnt_exceed++;
                            cont_block_cnt-=65535;
                        }
                    }

                    if( i < total)
                        if (data[i] != 0)
                            nonzero += sizeof_val; //8.0; //Double

                }
				
                cout << "Total number of zero bytes: " <<  zero_contiguous << " - Contiguous times: " << times_contiguous << "\n";
                cout << "Total non-zero bytes: " << nonzero << "\n";
                cout << "Total filesize non-compress (w header): " << nonzero + zero_contiguous + 32<< "\n";
                double comp_size = 32 + nonzero + 6.0*(double)(times_contiguous+cnt_exceed); // (4 bytes for signal + 2 bytes for value)
                cout << "Total filesize compressed (w header): " << comp_size << " bytes or " << comp_size/(1024.0*1024.0) << " MB\n";
                cout << "Compression ratio: " << (comp_size/(32+(total*sizeof_val)))*100 << "%  or " << (32+(total*sizeof_val))/comp_size << "\n";
                double overhead_mb = (comp_size/(1024.0*1024.0)) - ((32+(total*sizeof_val))/(1024.0*1024.0))*(argv_pcnt/100);
                cout << "Compression scheme overhead: " << overhead_mb << " MB or " << overhead_mb/(comp_size/(1024.0*1024.0))*100.0 << "%\n";
                cout << "DEBUG: Longest continuous block count: " << largest_cont_block << " - Times exceeded max unsigned short: " << cnt_exceed << "\n";
                //system("PAUSE");

				if(save_coeff)
				{
					char coeffnm[200];
					sprintf(coeffnm,"D:/Data/rstrt.0012_%03.0f_float_compressed_lz4.coeff",argv_pcnt);
					// Test saving coefficients //64 is RLE // 118-138 is lz4/fpr // 128 is only lz4
					// When range is 128, native lz4 is used. 
					// Less than 128, good for large dynamic range. More than 128, good for small dynamic ranges.

					
					//int * args = coeff_args( type, 0, 0+64, 0, pDims[0], pDims[1], pDims[2], pDims[0], pDims[1], pDims[2],sizeof(data[0]));
					int * args = coeff_args(type, 0, 0 + 133, 0, pDims[0], pDims[1], pDims[2], pDims[0], pDims[1], pDims[2], sizeof(data[0]));
					save_coefficients_md(data, coeffnm, args);
					free_args(args);
                
					system("PAUSE");
                
					double * data2 = read_coefficients_md<double>(coeffnm,args);
					for(size_t id=total-10; id<total; id++)
						std::cout << data2[id] << " ";


					/*for(size_t sz=0; sz<total; sz++)
						if(data[sz] != data2[sz])
							std::cout << "Compression mismatch: " << sz << ", orig: " << data[sz] << " comp " << data2[sz] << std::endl;*/
					std::cout << "COMPLETE! Sanity check passed\n";
					//delete [] data2;
					data = data2;
					//delete[] data;
					system("PAUSE");
				}
            }
            else
            {
                
                // Calculate dimensions per each level
                int *** cdims = getCoeffDims(pDims[0],pDims[1],pDims[2]);

                cout << "Thresholding by level: " << argv_lvl << " of " << max_levels << "\n";
                start = clock();
                int curdims[3] = {pDims[0],pDims[1],pDims[2]};
                for(int lvl=argv_lvl; lvl < max_levels; lvl++)
                {
                    int lvlsz=cdims[lvl][0][1]+1; // xmax+1 assuming 3D cube
                    lvlsz = lvlsz*lvlsz*lvlsz;
                    for(int i=1; i<8;i++)
                    {
                        double * LLL = getCutout3d(data,cdims[lvl][i], curdims);
                        for(int index=0;index<lvlsz;index++)
                            LLL[index] = 0.0;

                        setCutout3d(data,LLL,cdims[lvl][i], curdims);
                        delete [] LLL;
                    }
                }

                
                // Adaptive testing for a region of the dataset
                /*if(this->useAdaptive && argv_lvl < max_levels)
                {
                    cout << "Adaptively reconstruction corner of dataset\n";
                    for(int lvl=argv_lvl; lvl < max_levels; lvl++)
                    {
                        int width = this->cdims[lvl][0][1]+1;
                        int lvlsz=width*width*width; // assuming its a 3d cube cutout
                        for(int i=1; i<8;i++)
                        {
                            double * LLL = getCutout(coeffs,this->cdims[lvl][i], curdims);
                            int lpos[3] = {0,0,0};
                            for(int index=0;index<lvlsz;index++)
                            {
                                if(this->waveType == 1)
                                {
                                    if(lpos[0] > width/2 || lpos[1] > width/2 || lpos[2] > width/2 )
                                        LLL[index] = 0.0;
                                }
                                else
                                {
                                    if(lpos[0] > width/3 || lpos[2] > width/3 )//&& lpos[0] < width-(width/10))
                                        LLL[index] = 0.0;
                                }

                                lpos[0] += 1;
                                if(lpos[0] == width)
                                {
                                    lpos[0]=0;lpos[1]+=1;
                                    if(lpos[1] == width)
                                    {
                                        lpos[1]=0;lpos[2]+=1;
                                    }
                                }
                                //cout << "X: " << lpos[0] << " Y: " << lpos[1] << " Z: " << lpos[2] << "\n";
                            }

                            setCutout(reconst,LLL,this->cdims[lvl][i], curdims);
                            delete LLL;
                        }
                    }
                }*/

                // Free cdims
                for(int i=max_levels-1; i>=0;i--)
                {
                    for(int c=0; c<8 ;c++)
                    {
                        delete cdims[i][c];
                    }
                    delete cdims[i];
                }
                delete cdims;
                cdims = NULL;

                secs = (clock() - start) / (double) 1000;
                cout << "Ending Threshold. Time elapsed = " << secs << " secs" << endl;
            }

            start = clock();
            gsl_wavelet3d_nstransform<double> (w, data, stride, pDims[0], pDims[1], pDims[2], gsl_wavelet_backward, work);
            secs = (clock() - start) / (double) 1000;
            cout << endl << "Inverse Transform " << fld+1 << " done. Time elapsed = " << secs << " secs" << endl;
			std::cout << std::endl;

            output3d[fld] = arrayto3D(data,pDims[0], pDims[1], pDims[2]);
	        cout << "Done arrayTo3D conversion" << endl;
            delete [] data;
        }
        
        // Output data to file (or print)
        //for (i = 0; i < n; i++)
        // {
        //   printf ("%g\n", data[i]);
        // }
        cout << "Beginning comparison" << endl;
        for(int fld =0; fld < fields; fld++)
        {
            //compare
            double tot_en=0;
            double tot_diff=0;
            double mse=0;
            double max=-999;
            for(int i=0; i<pDims[0]; i++)
            {
                for(int j=0; j<pDims[1]; j++)
                {
	                for(int k=0; k<pDims[2];k++)
	                {
	                    if(output3d[fld][i][j][k]>max)
	                        max=output3d[fld][i][j][k];
	                    if(input3d[fld][i][j][k]>max)
	                        max=input3d[fld][i][j][k];
	                        
	                    complex<double> value(output3d[fld][i][j][k], 0.0);
	                    tot_en += norm(value);
	                    mse+=(pow((double)output3d[fld][i][j][k]-input3d[fld][i][j][k],(double)2.0));
	                    tot_diff += abs(output3d[fld][i][j][k]-input3d[fld][i][j][k]);
	                }
	            }
            }
            mse /= (total);
            //printf("For field %d of %d\n", fld+1, fields);
	        cout << "For field " << fld+1 << " of " << fields << endl;
            cout.precision(10);
            //printf("DIFF: %.10f\n", tot_diff);
            //printf("MSE: %.10f\n", mse);
            double PSNR = 10*log10(pow(max,2.0)/(mse));
            //printf("PSNR: %.10f\n", PSNR);
            //printf("Energy is %.10f\n", tot_en);
            cout << "DIFF: " << tot_diff << endl;
	        cout << "MSE: " << mse << endl;
	        cout << "PSNR: " << PSNR << endl;
	        cout << "Energy is " << tot_en << endl;

        }
        
        
        char filenm[200];
		if(save_vti)
		{
#ifdef _WIN32
			sprintf_s(
#else
		    snprintf(
#endif
		    filenm, "D:/Data/rstrt_twelve_bior4_4_%.2f.vti",argv_pcnt);
		
			std::cout << "Writing out file: " << filenm << std::endl;

			std::vector<std::string> field_names;
			for(int fld=0;fld<fields;fld++)
			{
				char tmpnm[20];
				sprintf_s(tmpnm,"Scalar %d",fld+1);
				field_names.push_back(tmpnm);
			}
			//output
			//arr_to_vti_3d_range_scale(output3d,filenm,0,pDims[0]-1, 0,pDims[1]-1,0, pDims[2]-1, field_names);

		}

		if(save_bin)
		{
#ifdef _WIN32
			sprintf_s(
#else
		    snprintf(
#endif
		    filenm, "D:/Data/rstrt_twelve_bior4_4_%.2f.bin",argv_pcnt);
			std::cout << "Writing out file: " << filenm << std::endl;
			// output
			//arr_to_bin_3d_range_scale(output3d,fields,0,filenm,0,pDims[0]-1, 0,pDims[1]-1,0, pDims[2]-1);
		}

	    cout << "Done. Freeing memory." << endl;
        for (int fld=0; fld<fields; fld++)
        {
            for (int i = 0; i < pDims[0]; i++)
            {
                for (int j = 0; j < pDims[1]; j++)
                {
                    delete [] output3d[fld][i][j];
                }
                delete [] output3d[fld][i];
            }
            delete [] output3d[fld];
        }
        delete [] output3d;
        output3d = NULL;
        
    }
    
    gsl_wavelet_free (w);
    gsl_wavelet_workspace_free (work);

    for (int fld=0; fld<fields; fld++)
    {
        for (int i = 0; i < pDims[0]; i++)
        {
            for (int j = 0; j < pDims[1]; j++)
            {
                delete [] input3d[fld][i][j];
            }
            delete [] input3d[fld][i];
        }
        delete [] input3d[fld];
    }
    delete [] input3d;
    input3d = NULL;

    return 0;
}
