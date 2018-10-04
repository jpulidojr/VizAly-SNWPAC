
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
//#include <math.h>
#include <sstream>

//Use existing lossless compressors
#include "lz4.h"

/*static int * pfs1;
static int dimx1, dimy1, dimz1;
static int header_size1, scalar_fields1;
static double field_len1;*/

// These functions typically used by LZ4
static size_t write_uint16(FILE* fp, uint16_t i)
{
	return fwrite(&i, sizeof(i), 1, fp);
}

static size_t write_bin(FILE* fp, const void* array, int arrayBytes)
{
	return fwrite(array, 1, arrayBytes, fp);
}


template <typename T>
int arr_to_vti_3d_range_scale(T **** in, const char * filename,int xl,int xu, int yl,int yu, int zl, int zu, std::vector<std::string> names)
{
	int sizeof_val = sizeof(in[0][0][0][0]);

    //if(!readConfigw())
    //{
    //    std::cout << "Error reading configuration file. Exiting..." << std::endl;
    //    return 0;
    //}
    //int *prefs = pfs1;
    //std::string * field_names = fn1;

    int dimx1 = (xu-xl)+1;//prefs[0];
    int dimy1 = (yu-yl)+1;//prefs[1];
    int dimz1 = (zu-zl)+1;//prefs[2];
    //header_size1 = prefs[3];
    int scalar_fields1 = names.size();//fld;//prefs[4];
    //int num_blocks = prefs[5];
    //std::cout << "Dimensions "<< dimx1 << " " << dimy1 << " " << dimz1 <<std::endl;
    //std::cout << "Header " << header_size1 << " Fields " << scalar_fields1<<std::endl;
    //std::cout << "Num Blocks " << num_blocks << std::endl;


    //field_len = (double)((length-header_size) / scalar_fields);
    size_t field_len1 = (double)dimx1*dimy1*dimz1*8;
    //std::cout << "Length per field: " << field_len1 << std::endl;

        // Skip the first bytes of the header
        //f.seekg (header_size, std::ios::cur);
    size_t len;
    len = strlen(filename);

    //std::string * part_filenames = new std::string[num_blocks+1];
	int num_blocks=1;

    // Define arrays for extent lengths
    int ** ext_strt = new int*[num_blocks];
    int ** ext_end = new int*[num_blocks];
    for (int m = 0; m < num_blocks; ++m)
    {
        ext_strt[m] = new int[3];
        ext_end[m] = new int[3];
    }
    // Define a size for 1 block of data
    int * block_size=new int[3];

    int blockID=0;
    // Set the extents for the entire dataset
    ext_strt[blockID][0]=xl;
    ext_strt[blockID][1]=yl;
    ext_strt[blockID][2]=zl;
    ext_end[blockID][0]=xu;
    ext_end[blockID][1]=yu;
    ext_end[blockID][2]=zu;

    // Set block size
    block_size[0]=(xu-xl)+1;block_size[1]=(yu-yl)+1;block_size[2]=(zu-zl)+1;

    // Start creating each individual block
    //for(int g=0; g<num_blocks; g++)
    //{
    int g = 0;
        std::ofstream out;
        // Create the output file
        printf("Creating %s\n", filename);
        out.open(filename, std::ifstream::binary);
        out.precision(20);

                // Write VTK format header
        out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl;
        out << "  <ImageData WholeExtent=\""<<ext_strt[g][0]<<" "<<ext_end[g][0]<<" "<<ext_strt[g][1];
        out << " "<<ext_end[g][1]<<" "<<ext_strt[g][2]<<" "<<ext_end[g][2]<<"\"";
        out << " Origin=\"0 0 0\" Spacing=\"1 1 1\">"<<std::endl;
            out << "    <Piece Extent=\""<<ext_strt[g][0]<<" "<<ext_end[g][0]<<" "<<ext_strt[g][1];
        out << " "<<ext_end[g][1]<<" "<<ext_strt[g][2]<<" "<<ext_end[g][2]<<"\">"<<std::endl;
            out << "      <PointData Scalars=\"density_scalars\">"<<std::endl;

        int bin_offset = 0;
        int * ghost = new int[3];
        // Write the headers for each of the scalar fields (assuming they're doubles)
        for(int i=0; i<scalar_fields1; i++)
        {
            // Factor in ghost data into the blocks
            if(ext_strt[g][0]!=0){ghost[0]=0;}else{ghost[0]=0;}
            if(ext_strt[g][1]!=0){ghost[1]=0;}else{ghost[1]=0;}
            if(ext_strt[g][2]!=0){ghost[2]=0;}else{ghost[2]=0;}
            //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
			if( sizeof_val == 8 )
				out << "        <DataArray type=\"Float64\" Name=\""<<names[i]<<"\" format=\"appended\" offset=\""<< (double)(i*((double)(block_size[0]+ghost[0])*(block_size[1]+ghost[1])*(block_size[2]+ghost[2])*8))+bin_offset <<"\"/>"<<std::endl;
			else if( sizeof_val == 4)
				out << "        <DataArray type=\"Float32\" Name=\""<<names[i]<<"\" format=\"appended\" offset=\""<< (double)(i*((double)(block_size[0]+ghost[0])*(block_size[1]+ghost[1])*(block_size[2]+ghost[2])*4))+bin_offset <<"\"/>"<<std::endl;
			else
			{
				std::cout << "ERROR: Data size of " << sizeof_val << " bytes not supported for vti!" << std::endl;
				return -1;
			}
            bin_offset += 4;
        }
        out << "      </PointData>"<<std::endl;
        out << "      <CellData>" << std::endl;
        out << "      </CellData>" << std::endl;
        out << "    </Piece>"<<std::endl;
        out << "  </ImageData>"<<std::endl;
        out << "  <AppendedData encoding=\"raw\">"<<std::endl;
        // Denote that you are about to write in binary
        out << "    _";

        //char * bin_value= new char[sizeof_val];//[8];
        //std::cout<< "[";
        for(int sf=0; sf < scalar_fields1; sf++)
        {
            // As per the API, first encode a 32-bit unsigned integer value specifying
            // the number of bytes in the block of data following it
            unsigned int block_byte_size = ((block_size[0]+ghost[0])*(block_size[1]+ghost[1])*(block_size[2]+ghost[2]))*8;
            if(((double)block_size[0]*block_size[1]*block_size[2]*8) > 4294967296.0)
            {
                printf("Error: Scalar field size is too large. Choose a larger\n");
                printf("number of blocks to segment by.\n");
                return 0;
            }
            out.write( reinterpret_cast<char*>( &block_byte_size ), 4 );
            int c=0;
            for(int z=ext_strt[g][2]; z<=ext_end[g][2]; z++)
            {
                int b=0;
                for(int y=ext_strt[g][1]; y<=ext_end[g][1]; y++)
                {
                    int a=0;
                    for(int x=ext_strt[g][0]; x<=ext_end[g][0]; x++)
                    {
                        T dbl_val = (T)in[sf][a][b][c];
                        out.write(reinterpret_cast<char*>( &dbl_val ), sizeof_val);//sizeof(double));
                        a++;
                    }
                    b++;
                }
                c++;
            }
            //std::cout<<"] - "<<field_names[sf]<<" done!\n[";
        }
        //printf("Complete!");
        //std::cout<<"]\n";
        //Finish off the xml file
        out << std::endl;
        out << "  </AppendedData>"<<std::endl;
        out << "</VTKFile>";

        //printf("%s Done!\n", part_filenames[g+1].c_str());
        out.close();
    //}
    //std::cout<<"Done!"<<std::endl;
    return 0;
}

template <typename T>
int save_coefficients(T * data, const char * filename, size_t total)
{

    std::ofstream out;
    // Create the output file
    printf("Creating %s\n", filename);
    out.open(filename, std::ifstream::binary);
    out.precision(20);

    for(int sz=0; sz < total; sz++)
    {
        T dbl_val = data[sz];
        out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(T));
    }
    std::cout<<"Coefficient's saved\n";

    out.close();
    return 1;
}

template <typename T>
int save_coefficients_md(T * data, const char * filename, int * args)
{
	
    std::ofstream out;
    // Create the output file
    printf("Creating %s\n", filename);
    out.open(filename, std::ifstream::binary);
    //out.precision(20);

    short * tmp_s = new short[1]; char * tmp_c = new char[1];
    int * tmp_i = new int[1];
    // Write the header
    tmp_s[0] = args[0]; out.write(reinterpret_cast<char*>( tmp_s ), sizeof(short));
    tmp_s[0] = args[1]; out.write(reinterpret_cast<char*>( tmp_s ), sizeof(short));
    tmp_c[0] = args[2]; out.write(reinterpret_cast<char*>( tmp_c ), sizeof(char));
    tmp_s[0] = args[3]; out.write(reinterpret_cast<char*>( tmp_s ), sizeof(short));
    tmp_i[0] = args[4]; out.write(reinterpret_cast<char*>( tmp_i ), sizeof(int));
    tmp_i[0] = args[5]; out.write(reinterpret_cast<char*>( tmp_i ), sizeof(int));
    tmp_i[0] = args[6]; out.write(reinterpret_cast<char*>( tmp_i ), sizeof(int));
    tmp_i[0] = args[7]; out.write(reinterpret_cast<char*>( tmp_i ), sizeof(int));
    tmp_i[0] = args[8]; out.write(reinterpret_cast<char*>( tmp_i ), sizeof(int));
    tmp_i[0] = args[9]; out.write(reinterpret_cast<char*>( tmp_i ), sizeof(int));
    tmp_c[0] = args[10];out.write(reinterpret_cast<char*>( tmp_c ), sizeof(char));
    delete [] tmp_s; delete [] tmp_c; delete [] tmp_i;

    size_t total = (size_t)args[4]*(size_t)args[5]*(size_t)args[6];
	// Check if we are using lz4 on coefficients
	if (args[2] >= 118)
	{
		std::cout << "Beginning LZ4 Routines...";
		// Optional: Perform floating point to integer requantization:
		// No Requantization if args[2]=128
		// When it's 125, it's 128-125 = 3. We divide by 1000
		// When it's 131, it's 128-131 = -3. Multibly by 1000
		double mod_bits = 1;
		size_t oval_sz = sizeof(T);

		if( args[2] > 128 ) // Multiply (for small dyn range data)
		{
			int cnt = args[2];
			while (cnt != 128)
			{
				mod_bits *= 10.0;
				cnt--;
			}
			oval_sz = sizeof(int);
			std::cout << "requantizing coefficients: 1*10^" << args[2] - 128 << std::endl;
		}
		if(args[2] < 128 ) // Divide (for high dyn range data)
		{ 
			int cnt = args[2];
			while (cnt != 128)
			{
				mod_bits /= 10.0;
				cnt++;
			}
			oval_sz = sizeof(int);
			std::cout << "requantizing coefficients: 1*10^-" << 128 - args[2] << std::endl;
		}

		
		// Target output variable (type)
		T dbl_val;
		//int dbl_val; // INT SUPPORT

		const size_t totBytes = total * oval_sz;
		size_t datBytes = totBytes;
		// WARNING: LZ4_MAX_INPUT_SIZE is about 2MB so suggest using a sliding buffer
		// Restriction inside of srcSize

		if (totBytes > LZ4_MAX_INPUT_SIZE)
		{
			std::cout << "Warning: Data to write is larger than supported by Lz4. Use rolling buffer!!\n";
			datBytes = LZ4_MAX_INPUT_SIZE;
		}

		LZ4_stream_t* lz4Stream = LZ4_createStream();
		const size_t cmpBufBytes = LZ4_COMPRESSBOUND(datBytes); //messageMaxBytes
		
		//Preallocate buffers
		char * inBuf = new char[datBytes];
		char * cmpBuf = new char[cmpBufBytes];

		// Use a sliding window approach and reinterpert type cast to char * array.
		// Use a ring buffer????
		
		size_t max_vals = datBytes / oval_sz;
		std::cout << " Encoding " << max_vals << " values.\n";
		
		size_t sk = 0;

		char * cval;

		if (oval_sz == 4)
		{
			int dbl_val; // int for requantization
			for (size_t sz = 0; sz < max_vals; sz++)
			{
				//dbl_val = data[sz]*100000; // INT SUPPORT
				dbl_val = data[sz]*mod_bits; // INT SUPPORT

				sk = sz * sizeof(dbl_val);
				cval = reinterpret_cast<char*>(&dbl_val);
				for (size_t id = 0; id < sizeof(dbl_val); id++)
				{
					inBuf[sk + id] = cval[id];
				}
			}
		}
		else if (oval_sz == 8)
		{
			T dbl_val; // No requant
			for (size_t sz = 0; sz < max_vals; sz++)
			{
				dbl_val = data[sz];
				sk = sz * sizeof(dbl_val);
				cval = reinterpret_cast<char*>(&dbl_val);
				for (size_t id = 0; id < sizeof(dbl_val); id++)
				{
					inBuf[sk + id] = cval[id];
				}
				//delete[] cval;
				//inBuf[sk] = reinterpret_cast<char*>(&dbl_val), sizeof(dbl_val); //Insert T val into character buffer
				//memcpy?
			}
		}
		else
		{
			std::cout << "ERROR: Invalid byte size!" << std::endl; return -1;
		}

		//const int cmpBytes = LZ4_compress_fast_continue(lz4Stream, inBuf, cmpBuf, datBytes, cmpBufBytes, 1); /rolling
		int cmpBytes = LZ4_compress_default(inBuf, cmpBuf, datBytes, cmpBufBytes);

		//Write out the buffer size first
		//write_uint16(FILE* fp, uint16_t i)
		out.write(reinterpret_cast<char*>(&cmpBytes), sizeof(cmpBytes));
		std::cout << "LZ4: Encoding " << cmpBytes << " bytes.." << std::endl;
		// Write out bytestream
		out.write(cmpBuf, cmpBytes);

		delete[] inBuf;
		delete[] cmpBuf;
		LZ4_freeStream(lz4Stream);

		std::cout << "..Done LZ4!\n";
	}
    // Check if we are encoding coefficients with RLE
    else if(args[2] >= 64)
    {
        std::cout<<"Encoding enabled...";
        //unsigned char signal = 0x24; //$ sign UTF-8
        char signal = '$';
        char signal2 = '@';
        char signal3 = '#';
        char signal4 = 'h';
        //std::cout << "Size of char=" << sizeof(signal) << " value " << signal << std::endl;
        unsigned short smax=0;

        float wbytes=0;
        int skips=0;
        for(size_t sz=0; sz < total; sz++)
        {
            
            if (data[sz] == 0)
            {    
                skips++;
                int cont_block_cnt=0;
                while (data[sz]==0 && sz < total)
                {
                    sz++;
                    cont_block_cnt++; // DEBUG
                }
                // Encode the number of steps to skip
                while(cont_block_cnt>65535) //Exceeded u_short_int
                {
                    smax = 65535;
                    out.write(reinterpret_cast<char*>( &signal ), sizeof(signal));
                    out.write(reinterpret_cast<char*>( &signal2 ), sizeof(signal2));
                    out.write(reinterpret_cast<char*>( &signal3 ), sizeof(signal3));
                    out.write(reinterpret_cast<char*>( &signal4 ), sizeof(signal4));
                    out.write(reinterpret_cast<char*>( &smax ), sizeof(smax));
                    skips++;
                    cont_block_cnt-=65535;
                    wbytes+=sizeof(signal)*4+sizeof(smax);
                }
                smax = cont_block_cnt;
                //std::cout << "\nSkipping at byte_loc: " << wbytes << " for " << smax << std::endl;
                out.write(reinterpret_cast<char*>( &signal ), sizeof(signal));
                out.write(reinterpret_cast<char*>( &signal2 ), sizeof(signal2));
                out.write(reinterpret_cast<char*>( &signal3 ), sizeof(signal3));
                out.write(reinterpret_cast<char*>( &signal4 ), sizeof(signal4));
                out.write(reinterpret_cast<char*>( &smax ), sizeof(smax));
                wbytes+=sizeof(signal)*4+sizeof(smax);
            }

            T dbl_val = data[sz];
            out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(dbl_val));
            wbytes+=sizeof(dbl_val);
        }
        std::cout<<"...Coefficient's saved\n";
        std::cout<<"Bytes written (w header): " << wbytes+32 << std::endl;
        std::cout << "Contiguous skips: " << skips << std::endl;
    }
	// Just write coefficients to storage with zeros intact
    else
    {
        for(size_t sz=0; sz < total; sz++)
        {
            T dbl_val = data[sz];
            out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(T));
        }
        std::cout<<"Coefficient's saved\n";
    }

    out.close();
    return 1;
}
