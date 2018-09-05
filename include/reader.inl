#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <math.h>
#include <sstream>
#include <string>

#include "lz4.h"

template <typename T> 
T * read_coefficients(const char * filename, size_t total)
{
    T * data = NULL;

    std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file: " << filename <<std::endl;
		return data;
	}
    data = new T[total];

	//std::cout << "Opening "<< filename << std::endl;
    //std::cout.precision(20);

    T dbl_val = 0.0;
    for(size_t sz=0; sz < total; sz++)
    {
        f.read(reinterpret_cast<char*>( &dbl_val ), sizeof(T));
		data[sz]=dbl_val;
    }
    //std::cout<<"Coefficient's loaded\n";

    f.close();
    return data;
}


template <typename T> T * read_adaptive_coefficients(const char * filename, size_t total)
{
    T * data = NULL;

    std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file: " << filename <<std::endl;
		return data;
	}
    data = new T[total];

	//std::cout << "Opening "<< filename << std::endl;
    //std::cout.precision(20);

    T dbl_val = 0.0;
    // lets attempt to read the center cube

    for(int sz=0; sz < total; sz++)
    {


    f.read(reinterpret_cast<char*>( &dbl_val ), sizeof(T));
		data[sz]=dbl_val;
    }

    // choose half to test adaptiveness


    //std::cout<<"Coefficient's loaded\n";

    f.close();
    return data;
}


template <typename T>
T * read_coefficients_md(const char * filename, int * args)
{
    T * data = NULL;

    std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file: " << filename <<std::endl;
		return data;
	}

    //Read metadata(header)
    args = new int[11];

    unsigned short * tmp_s = new unsigned short[1]; unsigned char * tmp_c = new unsigned char[1];
    int * tmp_i = new int[1];
    // Write the header
    f.read(reinterpret_cast<char*>( tmp_s ), sizeof(short)); args[0]= tmp_s[0];
    f.read(reinterpret_cast<char*>( tmp_s ), sizeof(short)); args[1]= tmp_s[0];
    f.read(reinterpret_cast<char*>( tmp_c ), sizeof(char));  args[2]= tmp_c[0];
    f.read(reinterpret_cast<char*>( tmp_s ), sizeof(short)); args[3]= tmp_s[0];
    f.read(reinterpret_cast<char*>( tmp_i ), sizeof(int));   args[4]= tmp_i[0];
    f.read(reinterpret_cast<char*>( tmp_i ), sizeof(int));   args[5]= tmp_i[0];
    f.read(reinterpret_cast<char*>( tmp_i ), sizeof(int));   args[6]= tmp_i[0];
    f.read(reinterpret_cast<char*>( tmp_i ), sizeof(int));   args[7]= tmp_i[0];
    f.read(reinterpret_cast<char*>( tmp_i ), sizeof(int));   args[8]= tmp_i[0];
    f.read(reinterpret_cast<char*>( tmp_i ), sizeof(int));   args[9]= tmp_i[0];
    f.read(reinterpret_cast<char*>( tmp_c ), sizeof(char)); args[10]= tmp_c[0];
    delete [] tmp_s; delete [] tmp_c; delete [] tmp_i;
    
     // debug output
    std::cout<< "Header Metadata read from coefficients:"<< std::endl;
    std::cout << "Type: " << args[0] << " Level: " << args[1] << std::endl;
    std::cout << "Region: " << args[2] << " Padding: " << args[3] << std::endl;
    std::cout << "Local dims: " << args[4] << " " << args[5] << " " << args[6] << std::endl;
    std::cout << "Global dims: " << args[7] << " " << args[8] << " " << args[9] << std::endl;
    std::cout << "Data precision: " << args[10] << " bytes." << std::endl;
    if (args[2] >=64 )
        std::cout << "Compression: Enabled" << std::endl;
	if (args[2] >= 128)
		std::cout << "LZ4 Enabled" << std::endl;


    size_t total = (size_t)args[4]*(size_t)args[5]*(size_t)args[6];

    data = new T[total]();

	if(sizeof(T) != args[10])
		std::cout << "Warning: Header Data Precision " << args[10] << " does not match target " << sizeof(T) << std::endl;

	//std::cout << "Opening "<< filename << std::endl;
    //std::cout.precision(20);
	if (args[2] >= 118)
	{
		int cmpBytes = 0;
		//char strm1[4];
		// Begin reading the number of bytes
		f.read(reinterpret_cast<char*>(&cmpBytes), sizeof(cmpBytes));
		
		std::cout << "LZ4: Reading in " << cmpBytes << " bytes..." << std::endl;

		// Optional: Reverse integer to fp requantization:
		// No Requantization if args[2]=128
		// When it's 125, it's 128-125 = 3. We multiply by 1000
		// When it's 131, it's 128-131 = -3. divide by 1000
		double mod_bits = 1;
		size_t oval_sz = sizeof(T);

		if (args[2] > 128) // Divide (undo requant)
		{
			int cnt = args[2];
			while (cnt != 128)
			{
				mod_bits /= 10.0;
				cnt--;
			}
			oval_sz = sizeof(int);
			std::cout << "requantizing coefficients: 1*10^-" << args[2] - 128  << std::endl;
		}
		if (args[2] < 128) // Multiply (undo requant)
		{
			int cnt = args[2];
			while (cnt != 128)
			{
				mod_bits *= 10.0;
				cnt++;
			}
			oval_sz = sizeof(int);
			std::cout << "requantizing coefficients: 1*10^" << 128 - args[2] << std::endl;
		}

		// Define the temporary type (requant test)
		//T dbl_val;
		//int dbl_val; // INT SUPPORT

		char * cmpBuf = new char[cmpBytes];
		char * outBuf = new char[total*oval_sz];
		f.read(cmpBuf, cmpBytes);

		std::cout << "Read in bytes..Decompressing LZ4.." << std::endl;
		
		// Allocate temporary memory
		LZ4_decompress_safe(cmpBuf, outBuf, cmpBytes, total * oval_sz);
		delete[] cmpBuf;

		std::cout << "...Done!";

		// Write byte memory into temp data array
		size_t sk = 0;

		if (oval_sz == 8) // Type DOUBLE
		{
			char conc[8];
			T dbl_val;

			for (size_t sz = 0; sz < total; sz++)
			{
				for (size_t ind = 0; ind < sizeof(dbl_val); ind++)
				{
					conc[ind] = outBuf[(sz * sizeof(dbl_val)) + ind];
				}
				dbl_val = *reinterpret_cast<T*>(conc);
				
				data[sz] = dbl_val;
				

				//outBuf[sz * sizeof(T)];

				//dbl_val = *reinterpret_cast<T*>(outBuf[sz * sizeof(dbl_val)]);
				//std::cout << " v: " << dbl_val << " ";
				//data[sz] = dbl_val;
				//dbl_val = data[sz];
				//sk = sz * sizeof(T);
				//char * cval = reinterpret_cast<char*>(&dbl_val);

				//for (size_t id = 0; id < sizeof(dbl_val); id++)
				//{
				//	inBuf[sk + id] = cval[id];
				//}
				//delete[] cval;
				//inBuf[sk] = reinterpret_cast<char*>(&dbl_val), sizeof(dbl_val); //Insert T val into character buffer
				//memcpy?
			}
		}
		else if (oval_sz == 4) // Type FLOAT input, may be float/double output
		{
			char conc[4]; // INT SUPPORT
			if (args[10] == 8) // Requantize ints to floats
			{
				int dbl_val;
				for (size_t sz = 0; sz < total; sz++)
				{
					for (size_t ind = 0; ind < 4; ind++)
					{
						conc[ind] = outBuf[(sz * sizeof(dbl_val)) + ind];
					}
					dbl_val = *reinterpret_cast<int*>(conc); // INT SUPPORT
					data[sz] = (T)dbl_val * mod_bits; // INT SUPPORT
				}
			}
			else if (args[10] == 4)
			{
				T dbl_val;
				for (size_t sz = 0; sz < total; sz++)
				{
					for (size_t ind = 0; ind < 4; ind++)
					{
						conc[ind] = outBuf[(sz * sizeof(dbl_val)) + ind];
					}
					data[sz] = *reinterpret_cast<T*>(conc);
				}
			}
			else
			{


			}
			
		}
		else if (args[10] == 2)
		{

		}
		else
		{
			std::cout << " Type mismatch. Not implemented!" << std::endl;
			return 0;
		}
		
		delete[] outBuf;
		std::cout << "..LZ4 finished data copy!\n";

	}
	else if( args[2] >= 64)
    {

        int skip_cnt=0;
        size_t val_read=0;
        double rbytes=0;
        std::cout<<"Encoding enabled...";
        //unsigned char signal = 0x24; //$ sign UTF-8
        char signal = '$';
        char signal2 = '@';
        char signal3 = '#';
        char signal4 = 'h';

        char strm1[4];
        char strm2[4]; //sizeof(var)-1
        unsigned short skip=0;
        for(size_t sz=0; sz < total; sz++)
        {
            f.read(strm1,sizeof(signal)*4);
            rbytes+=sizeof(signal)*4; //debug

            if(strm1[0] == signal && strm1[1]== signal2 && strm1[2]== signal3 && strm1[3]== signal4) //$ is unsigned, but read as a signed char. Conflict?
            //if(*reinterpret_cast<int*>(strm1) == signal)
            {
                skip_cnt++;
                
                f.read(reinterpret_cast<char*>( &skip ),sizeof(skip));
                rbytes+=sizeof(skip); //debug
                //if (skip_cnt < 10)
                //cout << "\nFound skip flag ["<<strm1[0]<<"]["<<strm1[1]<<"] at "<< rbytes <<" for size: " << skip << endl;
                if (skip > total)
                {
                    std::cout << "ERROR: Data read error: Skip size > total size\n";
                    exit(0);
                }
                //cout << "Found skip flag ["<<*reinterpret_cast<int*>(strm1)<<"] at "<< rbytes <<" for size: " << skip << endl;

                while(skip>1 && sz < total)
                {
                    data[sz]=0.0;
                    sz++;
                    skip=skip-1;
                }
            }
            else
            {

				if(args[10]==8) // Type DOUBLE
				{
					// Read second half of the data
					f.read(strm2,4);
					rbytes+=4;

					char conc[8];
					//Combine two sections
					//memcpy(&conc[0], strm1, 2);
					//memcpy(&conc[3], strm2, 6);
					std::copy(strm1,strm1+4,conc);
					std::copy(strm2,strm2+4,conc+4);

					
					//double dbl_val = *reinterpret_cast<double*>(conc);
					//std::cout << "Read double: " << data[sz] << std::endl;
					//data[sz] = dbl_val;

					data[sz] = *reinterpret_cast<T*>(conc); 
					val_read+=1;
				}
				else if(args[10]==4) // Type FLOAT
				{
					// 4 bytes already in memory, just convert
					char conc[4];
					std::copy(strm1,strm1+4,conc);
					data[sz] =  *reinterpret_cast<T*>(conc);
					val_read+=1;
				}
				/*else if(args[10]==2) // Type SHORT
				{
					// 4 bytes already in memory, split into 2 values
					char conc[2];
					std::copy(strm1,strm1+2,conc);
					data[sz] =  *reinterpret_cast<T*>(conc);
					val_read+=1;
					sz++;

					//TODO: These two bytes may be the beginning of the signal.
					//		Disable control block until above is fixed
					std::copy(strm1+2,strm1+4,conc);
					data[sz] =  *reinterpret_cast<T*>(conc);
					val_read+=1;
				}*/
				else
				{
					std::cout << "ERROR: Data encoded in unsupported type. Bytes per value: " << args[10] << std::endl;
					delete [] data;
					return NULL;
				}
            }


            /*
            if (data[sz] == 0)
            {    
                int cont_block_cnt=0;
                while (data[sz]==0 && sz < total)
                {
                    sz++;
                    cont_block_cnt++; // DEBUG
                }
                // Encode the number of steps to skip
                while(cont_block_cnt>65535) //Exceeded u_short_int
                {
                    unsigned short smax = 65535;
                    out.write(reinterpret_cast<char*>( &signal ), sizeof(char));
                    out.write(reinterpret_cast<char*>( &smax ), sizeof(short));

                    cont_block_cnt-=65535;
                }
                out.write(reinterpret_cast<char*>( &signal ), sizeof(char));
                out.write(reinterpret_cast<char*>( &cont_block_cnt ), sizeof(short));

            }

            double dbl_val = data[sz];
            out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(double));*/
        }
        std::cout<<"...Coefficient's loaded\n";

        std::cout << "Bytes Read (w header): " << rbytes+32 << endl;
        std::cout << "Contiguous times: " << skip_cnt << std::endl;
        std::cout << "Values Read: " << val_read << std::endl;

    }
    else
    {
        T dbl_val = 0.0;
        for(size_t sz=0; sz < total; sz++)
        {
            f.read(reinterpret_cast<char*>( &dbl_val ), sizeof(T));   
		    data[sz]=dbl_val;
        }
        //std::cout<<"Coefficient's loaded\n";
    }

    f.close();
    return data;
}