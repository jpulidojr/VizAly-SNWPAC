// Author: Jesus Pulido
// Contact: pulido@lanl.gov
/*
    This file will read a binary restart file into an internal 3d array.
    To change settings, see config.txt and modify accordingly.
	The contents of this file have been taken from my own implementation of bintopvti.cpp
*/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <math.h>
#include <sstream>
#include <string>
//#include "mpi.h"
#include "reader.h"

using namespace std;

#ifdef __cplusplus_cli
// .Net includes
#using <System.dll>
#using <System.Data.dll>
#include <vcclr.h>


using namespace System::Collections::Generic;
using namespace System::Text;
using namespace System::Data;
using namespace System::Data::SqlClient;
using namespace System::Collections;
using namespace System::Runtime::InteropServices;
#endif

// Global declarations
std::string * fn;
int * pfs;

int dimx, dimy, dimz;
int header_size, scalar_fields;
double field_len;

/* Simply reads custom parameters from file before conversion starts.
    See the file for details on what parameters are accessible
*/
int readConfig()
{
    std::ifstream config;
    config.open("config.txt");
    if (config.fail()) {
        std::cout << "Error: Unable to open config.txt" << std::endl;
		return 0;
	}

    pfs = new int[6];

    // Start reading file
    std::string line;
    while(config.good())
    {
        std::getline(config,line);
        //std::cout << line << std::endl;
        if(line.find("dimension_size")!=std::string::npos)
        {
            //std::cout<< "dim found" << std::endl;
            int delim1 = line.find(' ');
            int delim2 = line.find(' ',delim1+1);
            int delim3 = line.find(' ',delim2+1);
            int delim4 = line.length();
            pfs[0] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            pfs[1] = atoi(line.substr(delim2+1,delim3-delim2-1).c_str());
            pfs[2] = atoi(line.substr(delim3+1,delim4-delim3-1).c_str());
            continue;
        }
        if(line.find("header_size")!=std::string::npos)
        {
            //std::cout<< "header found" << std::endl;
            int delim1 = line.find(' ');
            int delim2 = line.length();
            pfs[3] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            continue;
        }
        if(line.find("scalar_fields")!=std::string::npos)
        {
            //std::cout<< "scalar found" << std::endl;
            int delim1 = line.find(' ');
            int delim2 = line.length();
            pfs[4] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            fn = new std::string[pfs[4]];
            continue;
        }
        if(line.find("field_names")!=std::string::npos)
        {
            //std::cout<< "fieldnames found" << std::endl;
            int prev_delim=-1;
            for(int i=0; i<pfs[4];i++)
            {
                if(i==pfs[4]-1) // end string
                {
                    int delim1=line.find(' ',prev_delim+1);
                    int delim2=line.length();
                    fn[i]=line.substr(delim1+1,delim2-delim1-1);
                }
                else // any other string
                {
                    int delim1=line.find(' ',prev_delim+1);
                    int delim2=line.find(' ',delim1+1);
                    fn[i]=line.substr(delim1+1,delim2-delim1-1);
                    prev_delim = delim1;
                }
            }
            continue;
        }
        if(line.find("num_blocks")!=std::string::npos)
        {
            //std::cout<< "number of blocks found" << std::endl;
            int delim1 = line.find(' ');
            int delim2 = line.length();
            pfs[5] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            continue;
        }        
    }

    //Debug output
    /*for(int i=0; i<6;i++)
    {
        std::cout<<pfs[i]<<std::endl;
    }
    for(int i=0; i<pfs[4];i++)
    {
        std::cout << fn[i] <<std::endl;
    }*/

    return 1;
}
/* Given x,y,z coordinates and current scalar field, this will find the 
    byte location of the data within a file
*/
double findByteAddress(int x, int y, int z, int field)
{
    int byte_size = 8; //size of double

    // add the header
    double byteLocation = (double)header_size;

    // first look at the scalar_field
    byteLocation += (double)field*dimx*dimy*dimz*byte_size;
    //byteLocation += field*field_len;

    // Look at z
    byteLocation += (double)z*dimx*dimy*byte_size;

    // look at y
    byteLocation += (double)y*dimx*byte_size;

    // look at x
    byteLocation += (double)x*byte_size;

    //printf("%f\n", byteLocation);
    return byteLocation;
}

/* Given byte address, return the x,y,z coordinates
*/
int * findXYZAddress(double byteLocation, int field)
{
    int byte_size = 8; //size of double
    int * pos = new int[3]; //x,y,z
    pos[0]=0;pos[1]=0;pos[2]=0;
    // remove the header
    byteLocation -= (double)header_size;


    // first look at the scalar_field
    byteLocation -= (double)field*dimx*dimy*dimz*byte_size;

    // Look at z
    while(byteLocation-((pos[2]+1)*dimx*dimy*byte_size) >= 0)
        pos[2]+=1;
    byteLocation -= (double)pos[2]*dimx*dimy*byte_size;
    
    // look at y
    while(byteLocation-((pos[1]+1)*dimx*byte_size) >= 0)
        pos[1]+=1;
    byteLocation -= (double)pos[1]*dimx*byte_size;

    // look at x
    pos[0]= byteLocation / byte_size;

    byteLocation -= pos[0]*byte_size;

    cout << "final address: " << byteLocation << ", X: " << pos[0] << " Y: " << pos[1] << " Z: " << pos[2] << endl;
    //printf("%f\n", byteLocation);
    return pos;
}

/* Given a 1D index, return the x,y,z coordinates
*/
int * findXYZIndex(int index)
{
    if(index < 0)
    {
        cout << "ERROR: Bad Index, " << index << endl;
        return 0;
    }

    int * pos = new int[3]; //x,y,z
    pos[0]=0;pos[1]=0;pos[2]=0;

    // Look at z
    while(index-((pos[2]+1)*dimx*dimy) >= 0)
        pos[2]+=1;
    index -= pos[2]*dimx*dimy;
    
    // look at y
    while(index-((pos[1]+1)*dimx) >= 0)
        pos[1]+=1;
    index -= pos[1]*dimx;

    // look at x
    pos[0]= index;

    index -= pos[0];

    cout << "final address: " << index << ", X: " << pos[0] << " Y: " << pos[1] << " Z: " << pos[2] << endl;
    //printf("%f\n", byteLocation);
    return pos;
}

int findIndexXYZ(int x, int y, int z)
{
    int location = 0;

    // first look at the scalar_field
    //location += (double)dimx*dimy*dimz*byte_size;

    // Look at z
    location += z*dimx*dimy;

    // look at y
    location += y*dimx;

    // look at x
    location += x;

    //printf("%f\n", byteLocation);
    return location;
}

int binto2Darray(float ** in, char * filename,int xl,int xu, int yl,int yu, int zl, int zu ){
	if( xl&&xu&&yl&&yu&&zl&&zu )
	{
		printf("Error: three dimensions specified for a 2D function call. Defaulting to X-Y axis\n");
		zu=0;
		zl=0;
	}
    if(!readConfig())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs;
    std::string * field_names = fn;

    dimx = prefs[0];
    dimy = prefs[1];
    dimz = prefs[2];
    header_size = prefs[3];
    scalar_fields = prefs[4];
    int num_blocks = prefs[5];
    //std::cout << "Dimensions "<< dimx << " " << dimy << " " << dimz <<std::endl;
    //std::cout << "Header " << header_size << " Fields " << scalar_fields<<std::endl;
    //std::cout << "Num Blocks " << num_blocks << std::endl;

	std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file"<<std::endl;
		return 0;
	}

	std::cout << "Opening "<< filename << std::endl;
    std::cout.precision(20);
	// Get the length of file:
	f.seekg (0, std::ios::end);
	double length = f.tellg();
	f.seekg (0, std::ios::beg);
	std::cout << "Length of file: " << length << std::endl;

    //field_len = (double)((length-header_size) / scalar_fields);
    field_len = (double)dimx*dimy*dimz*8;
    std::cout << "Length per field: " << field_len << std::endl;

	// Skip the first bytes of the header
	//f.seekg (header_size, std::ios::cur);
    int len;
    len = strlen(filename);

    std::string * part_filenames = new std::string[num_blocks+1];

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


    // Just generate the name for 1 file
    part_filenames[1].assign(filename);
    std::string tempfilename = part_filenames[1].substr(0,len-4);
    tempfilename.append(".vti");
    part_filenames[1] = tempfilename;

    int blockID=0;
    // Set the extents for the entire dataset
    ext_strt[blockID][0]=0;
    ext_strt[blockID][1]=0;
    ext_strt[blockID][2]=0;
    ext_end[blockID][0]=dimx-1;
    ext_end[blockID][1]=dimy-1;
    ext_end[blockID][2]=dimz-1;

    // Set block size
    block_size[0]=dimx;block_size[1]=dimy;block_size[2]=dimz;

    // Start creating each individual block
    for(int g=0; g<num_blocks; g++)
    {   
        int bin_offset = 0;
        int * ghost = new int[3];
        // Write the headers for each of the scalar fields (assuming they're doubles)
        for(int i=0; i<scalar_fields; i++)
        {
            // Factor in ghost data into the blocks
            if(ext_strt[g][0]!=0){ghost[0]=1;}else{ghost[0]=0;}
            if(ext_strt[g][1]!=0){ghost[1]=1;}else{ghost[1]=0;}
            if(ext_strt[g][2]!=0){ghost[2]=1;}else{ghost[2]=0;}
            //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
            bin_offset += 4;
        }
        char * bin_value = new char[8];
        std::cout<< "[";
        for(int sf=0; sf < scalar_fields; sf++)
        {  
            // As per the API, first encode a 32-bit unsigned integer value specifying
            // the number of bytes in the block of data following it
            for(int z=ext_strt[g][2]; z<=ext_end[g][2]; z++)
            {
                for(int y=ext_strt[g][1]; y<=ext_end[g][1]; y++)
                {
                    for(int x=ext_strt[g][0]; x<=ext_end[g][0]; x++)
                    {
                        std::streampos seekVal = findByteAddress(x,y,z,sf);
                        f.seekg (seekVal);
                        //std::cout << "seekVal: "<<f.tellg() <<std::endl;
                        //if (f.eof()) //This kills performance
                        //{
                        //    printf("Error: Byte Address is out of bounds at %f\n",
                        //                                 findByteAddress(x,y,z,sf));
                        //    printf("%i %i %i %i\n", x,y,z,sf);
                        //    printf("Exiting.....\n");
                        //    return 0;
                        //}
						double value;
                        f.read(reinterpret_cast<char*>( &value ), sizeof(double));
						in[x][y]=(float)value;
						//printf("Value: %f\n", value);
                    }
                }
            }
            std::cout<<"] - "<<field_names[sf]<<" done!\n[";
        }
        printf("Complete!");
        std::cout<<"]\n";
        //Finish off the xml file

        printf("%s Done!\n", part_filenames[g+1].c_str());
    }
    f.close();


    std::cout<<"Done!"<<std::endl;
	return 1;
}

int binto3Darray(float *** in, char * filename,int xl,int xu, int yl,int yu, int zl, int zu ){

    if(!readConfig())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs;
    std::string * field_names = fn;

    dimx = prefs[0];
    dimy = prefs[1];
    dimz = prefs[2];
    header_size = prefs[3];
    scalar_fields = prefs[4];
    int num_blocks = prefs[5];
    std::cout << "Dimensions "<< dimx << " " << dimy << " " << dimz <<std::endl;
    std::cout << "Header " << header_size << " Fields " << scalar_fields<<std::endl;
    std::cout << "Num Blocks " << num_blocks << std::endl;

	std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file"<<std::endl;
		return 0;
	}

	std::cout << "Opening "<< filename << std::endl;
    std::cout.precision(20);
	// Get the length of file:
	f.seekg (0, std::ios::end);
	double length = f.tellg();
	f.seekg (0, std::ios::beg);
	std::cout << "Length of file: " << length << std::endl;

    //field_len = (double)((length-header_size) / scalar_fields);
    field_len = (double)dimx*dimy*dimz*8;
    std::cout << "Length per field: " << field_len << std::endl;

	// Skip the first bytes of the header
	//f.seekg (header_size, std::ios::cur);
    int len;
    len = strlen(filename);

    std::string * part_filenames = new std::string[num_blocks+1];

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


    // Just generate the name for 1 file
    part_filenames[1].assign(filename);
    std::string tempfilename = part_filenames[1].substr(0,len-4);
    tempfilename.append(".vti");
    part_filenames[1] = tempfilename;

    int blockID=0;
    // Set the extents for the entire dataset
    ext_strt[blockID][0]=xl;
    ext_strt[blockID][1]=yl;
    ext_strt[blockID][2]=zl;
    ext_end[blockID][0]=xu;
    ext_end[blockID][1]=yu;
    ext_end[blockID][2]=zu;

    // Set block size
    block_size[0]=dimx;block_size[1]=dimy;block_size[2]=dimz;

    // Start creating each individual block
    for(int g=0; g<num_blocks; g++)
    {   
        int bin_offset = 0;
        int * ghost = new int[3];
        // Write the headers for each of the scalar fields (assuming they're doubles)
        for(int i=0; i<scalar_fields; i++)
        {
            // Factor in ghost data into the blocks
            if(ext_strt[g][0]!=0){ghost[0]=1;}else{ghost[0]=0;}
            if(ext_strt[g][1]!=0){ghost[1]=1;}else{ghost[1]=0;}
            if(ext_strt[g][2]!=0){ghost[2]=1;}else{ghost[2]=0;}
            //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
            bin_offset += 4;
        }
        char * bin_value = new char[8];
        std::cout<< "[";
        for(int sf=0; sf < scalar_fields; sf++)
        {  
            // As per the API, first encode a 32-bit unsigned integer value specifying
            // the number of bytes in the block of data following it
            for(int z=ext_strt[g][2]; z<=ext_end[g][2]; z++)
            {
                for(int y=ext_strt[g][1]; y<=ext_end[g][1]; y++)
                {
                    for(int x=ext_strt[g][0]; x<=ext_end[g][0]; x++)
                    {
                        std::streampos seekVal = findByteAddress(x,y,z,sf);
                        f.seekg (seekVal);
                        //std::cout << "seekVal: "<<f.tellg() <<std::endl;
                        //if (f.eof()) //This kills performance
                        //{
                        //    printf("Error: Byte Address is out of bounds at %f\n",
                        //                                 findByteAddress(x,y,z,sf));
                        //    printf("%i %i %i %i\n", x,y,z,sf);
                        //    printf("Exiting.....\n");
                        //    return 0;
                        //}
						double value;
                        f.read(reinterpret_cast<char*>( &value ), sizeof(double));
						in[x][y][z]=(float)value;
						//printf("Value: %f\n", value);
                    }
                }
            }
            std::cout<<"] - "<<field_names[sf]<<" done!\n[";
        }
        printf("Complete!");
        std::cout<<"]\n";
        //Finish off the xml file

        printf("%s Done!\n", part_filenames[g+1].c_str());
    }
    f.close();


    std::cout<<"Done!"<<std::endl;
	return 1;
}

int binto3DarrayMPI(double *** in, char * filename,int xl,int xu, int yl,int yu, int zl, int zu, int field){
    if(!readConfig())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs;
    std::string * field_names = fn;

    int mpirank=0;//  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    int mpisize=1;//  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    dimx = prefs[0];
    dimy = prefs[1];
    dimz = prefs[2];
    header_size = prefs[3];
    scalar_fields = prefs[4];
    int num_blocks = mpisize;
    //std::cout << "Dimensions "<< dimx << " " << dimy << " " << dimz <<std::endl;
    //std::cout << "Header " << header_size << " Fields " << scalar_fields<<std::endl;
    //std::cout << "Num Blocks " << num_blocks << std::endl;

	std::ifstream f;


    if(field >= scalar_fields)
    {
        std::cout <<"Error: Invalid scalar field chosen"<< std::endl;
        return 0;
    }
	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file"<<std::endl;
		return 0;
	}

//	std::cout << "Opening "<< filename << std::endl;
//    std::cout.precision(20);
	// Get the length of file:
	f.seekg (0, std::ios::end);
	double length = f.tellg();
	f.seekg (0, std::ios::beg);
//	std::cout << "Length of file: " << length << std::endl;

    //field_len = (double)((length-header_size) / scalar_fields);
    field_len = (double)dimx*dimy*dimz*8;
//    std::cout << "Length per field: " << field_len << std::endl;
/*
	// Skip the first bytes of the header
	//f.seekg (header_size, std::ios::cur);
    int len;
    len = strlen(filename);

    std::string * part_filenames = new std::string[num_blocks+1];

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


    // Just generate the name for 1 file
    part_filenames[1].assign(filename);
    std::string tempfilename = part_filenames[1].substr(0,len-4);
    tempfilename.append(".vti");
    part_filenames[1] = tempfilename;

    // Set the extents for each block of the sataset
    // Compute the extents for each piece file
    int blockID=0;
    int i, j, k, b;
    b = num_blocks;
    int * block_dim=new int[3];
    block_dim[0] = 1;
    block_dim[1] = 1;
    block_dim[2] = 1;
    block_size=new int[3];
    block_size[0] = dimx;
    block_size[1] = dimy;
    block_size[2] = dimz;
    
    while(1)
    {
        // find longest dimension to cut first
        int max = 0;
        for(i=1; i<3; i++) {
            if (block_size[i] > block_size[max]) {
                max = i;
            }
        }
        // smallest factor remaining gets assigned to this direction
        for(j=2; j<=b; j++)
        {
            if (b % j == 0)
            {
                block_dim[max] *= j;
                block_size[max] /= j;
                b /= j;
                break;
            }
        }
        if (b == 1)
          break;
        if (j > b)
        {
            printf("Error. Unable to partition the data into %i blocks. Try another size of blocks by ^2.\n", num_blocks);
            return 0;
        }     
    } 

    for(int z=0; z<block_dim[2]; z++)
    {
        for(int y=0; y<block_dim[1]; y++)
        {
            for(int x=0; x<block_dim[0]; x++)
            {
                // Provide 1 level of ghost data to fill in the gaps? (DISABED FOR READIN)
                
                ext_strt[blockID][0]=x*block_size[0];    
                ext_strt[blockID][1]=y*block_size[1];
                ext_strt[blockID][2]=z*block_size[2];
            
                ext_end[blockID][0]=(x+1)*block_size[0]-1;
                ext_end[blockID][1]=(y+1)*block_size[1]-1;
                ext_end[blockID][2]=(z+1)*block_size[2]-1;

                if(mpirank == blockID)
                    printf("Extents for block %i are %i %i %i %i %i %i \n",
                    blockID, ext_strt[blockID][0],ext_end[blockID][0],ext_strt[blockID][1],
                    ext_end[blockID][1],ext_strt[blockID][2],ext_end[blockID][2]);
                blockID++;
            }
        }
    }

    int g = mpirank;
    // Start creating each individual block
    //for(int g=0; g<num_blocks; g++)
    //{   
        int bin_offset = 0;
        int * ghost = new int[3];
        // Write the headers for each of the scalar fields (assuming they're doubles)
        for(int i=0; i<scalar_fields; i++)
        {
            // Factor in ghost data into the blocks
            if(ext_strt[g][0]!=0){ghost[0]=0;}else{ghost[0]=0;}
            if(ext_strt[g][1]!=0){ghost[1]=0;}else{ghost[1]=0;}
            if(ext_strt[g][2]!=0){ghost[2]=0;}else{ghost[2]=0;}
            //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
            bin_offset += 4;
        }*/
        char * bin_value = new char[8];
//        std::cout<< "[";
        
        //for(int sf=0; sf < scalar_fields; sf++)
        //{  
            // As per the API, first encode a 32-bit unsigned integer value specifying
            // the number of bytes in the block of data following it
            int c=0;
            for(int z=zl; z<=zu; z++)
            {

                int b=0;
                for(int y=yl; y<=yu; y++)
                {
                    int a=0;
                    std::streampos seekVal = findByteAddress(xl,y,z,field);
                    f.seekg (seekVal);
                    double value;
                    for(int x=xl; x<=xu; x++)
                    {
						
                        f.read(reinterpret_cast<char*>( &value ), sizeof(double));
                        in[a][b][c]=(double)value;

                        a++;
                    }
                    b++;
                }
                c++;
            }
            //std::cout<<"] - "<<field_names[sf]<<" done!\n[";
        //}
        //printf("Complete!");
        //std::cout<<"]\n";
        //Finish off the xml file

    //}
    f.close();


    //std::cout<<"Done Reading Block!"<<std::endl;
	return 1;
}



int * read_md(const char * filename)
{
    std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file: " << filename <<std::endl;
		return 0;
	}

    //Read metadata(header)
    int * args = new int[11];

    short * tmp_s = new short[1]; char * tmp_c = new char[1];
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
    //std::cout<< "Header Metadata read from coefficients:"<< std::endl;
    //std::cout << "Type: " << args[0] << " Level: " << args[1] << std::endl;
    //std::cout << "Region: " << args[2] << " Padding: " << args[3] << std::endl;
    //std::cout << "Local dims: " << args[4] << " " << args[5] << " " << args[6] << std::endl;
    //std::cout << "Global dims: " << args[7] << " " << args[8] << " " << args[9] << std::endl;
    //std::cout << "Data precision: " << args[10] << " bytes." << std::endl;
    //f.close()

    return args;
}

size_t * read_sort(const char * filename, size_t total)
{
    size_t * sort = new size_t[total];

    std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout<<"Error opening file"<<std::endl;
		return 0;
	}

	std::cout << "Opening "<< filename << std::endl;
    std::cout.precision(20);

    unsigned int uint_val=0;
    for(int sz=0; sz < total; sz++)
    {
        f.read(reinterpret_cast<char*>( &uint_val ), sizeof(unsigned int));
        sort[sz] = uint_val;

    }
    std::cout<<"Sorted Positions loaded\n";

    f.close();
    return sort;
}

int dbto3Darray(double *** in,int xl,int xu, int yl,int yu, int zl, int zu, int field)
{

#ifdef __cplusplus_cli
    //gcroot<String^> serverName = gcnew String("gwwn1");
	//gcroot<String^> codeDbName = gcnew String("mhddev");
    //String^ cString = String::Format("server={0};database={1};User ID=turbquery;Password=aa2465ways2k;", serverName, codeDbName);
    System::String ^ cString = "server=gwwn1;database=mhddev;User ID=turbquery;Password=aa2465ways2k;";
    //System::String ^ cString = "server=dsp048;database=turblib;User ID=jpulido;Password=pulido;";
    // For D.Livescu's data, the data is fast axis (x) wide into 8 partitions
    SqlConnection ^ sqlcon = gcnew SqlConnection(cString);
    sqlcon->Open();
	if (sqlcon->State == ConnectionState::Open)
		cerr << "JHTDBReader: connected to database" << endl;
	else
		cerr << "JHTDBReader: error connecting to database " << endl;

	System::String^ queryBox = System::String::Format("box[{0},{1},{2},{3},{4},{5}]", 
		xl, yl, zl, xu + 1, yu + 1, zu + 1);
	
	SqlCommand ^ command = gcnew SqlCommand();
	command->Connection = sqlcon;
	command->CommandType = CommandType::StoredProcedure;
	command->CommandText = "GetDataCutout";
    //command->CommandText = "GetSplineCutout";
	command->Parameters->Add(gcnew SqlParameter("@serverName",SqlDbType::VarChar));
	command->Parameters["@serverName"]->Value = "gwwn1";//this->serverName;
	command->Parameters->Add(gcnew SqlParameter("@dbname",SqlDbType::VarChar));
	command->Parameters["@dbname"]->Value = "turbdb101";//this->dbName; //mixingdb01 <-> 08
	command->Parameters->Add(gcnew SqlParameter("@codedb",SqlDbType::VarChar));
	command->Parameters["@codedb"]->Value = "mhddev";//codeDbName;
	command->Parameters->Add(gcnew SqlParameter("@dataset",SqlDbType::VarChar));
	command->Parameters["@dataset"]->Value = "vel";//this->dataset;
	command->Parameters->Add(gcnew SqlParameter("@blobDim",SqlDbType::Int));
	command->Parameters["@blobDim"]->Value = "8";//this->atomSize;
	command->Parameters->Add(gcnew SqlParameter("@timestep",SqlDbType::Int));
	command->Parameters["@timestep"]->Value = "0";//timeStep; //100 - 600 daniels data
	command->Parameters->Add(gcnew SqlParameter("@queryBox",SqlDbType::VarChar));
	command->Parameters["@queryBox"]->Value = queryBox;

	SqlDataReader ^ reader = command->ExecuteReader();

    cout << "Query executed\n";
    int components = 3;
    int dataPointSize = sizeof(float);
	int x_width = (xu - xl + 1);
	int y_width = (yu - yl + 1);
	int z_width = (zu - zl + 1);
	int data_size = x_width * y_width * z_width * components * dataPointSize;
	gcroot<array<unsigned char>^> rawdata = gcnew array<unsigned char>(data_size); //x/y/z

	__int64 bytesRead = 0;
	int bufferIndex = 0;
	while (reader->Read())
	{
		bufferIndex = (int)bytesRead;
		bytesRead += reader->GetBytes(0, 0, rawdata, bufferIndex, data_size - bufferIndex);
	}
	reader->Close();

    cout << "Query Successful, reformatting to output array\n";

    float * values = new float[components];
    int source_offset;
    for (int z = 0; z < z_width; ++z)
	{
		for (int y = 0; y < y_width; ++y)
		{
			for (int x = 0; x < x_width; ++x)
			{
				//int source_offset = z * CachedYWidth * CachedXWidth + y * CachedXWidth + x;
                source_offset = z * y_width * x_width + y * x_width + x;
				//int dest_offset = z * y_width * x_width + y * x_width + x;
				for (int c = 0; c < components; c++)
				{
					values[c] = System::BitConverter::ToSingle(rawdata, (source_offset * components + c) * dataPointSize);
				}
				//scalars->SetTuple(dest_offset, values);
				//float value = BitConverter::ToSingle(rawdata, source_offset * components * dataPointSize);
				//scalars->SetTuple1(dest_offset, value);
                in[x][y][z] = sqrt(pow(values[0],2)+pow(values[1],2)+pow(values[2],2));
			}
		}
	}
    cout << "Sucessful read\n";
#else
    cout << "ERROR: Program compiled without /CLR. JHUDB Reader disabled\n:";
    return 0;
#endif
    return 1;
}
