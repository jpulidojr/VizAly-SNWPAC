
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
//#include <math.h>
#include <sstream>
#include "writer.h"

int arr_to_vti_3d_range_scale(float **** in, const char * filename,int xl,int xu, int yl,int yu, int zl, int zu);

// Global declarations
std::string * fn1;
int * pfs1;

int dimx1, dimy1, dimz1;
int header_size1, scalar_fields1;
double field_len1;

int arr_to_pvti_2d(float ** in, char * filename)
{
    
    if(!readConfigw())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs1;
    std::string * field_names = fn1;

    dimx1 = prefs[0];
    dimy1 = prefs[1];
    dimz1 = prefs[2];
    header_size1 = prefs[3];
    scalar_fields1 = prefs[4];
    int num_blocks = prefs[5];
    std::cout << "Dimensions "<< dimx1 << " " << dimy1 << " " << dimz1 <<std::endl;
    std::cout << "Header " << header_size1 << " Fields " << scalar_fields1<<std::endl;
    std::cout << "Num Blocks " << num_blocks << std::endl;


    //field_len = (double)((length-header_size) / scalar_fields);
    field_len1 = (double)dimx1*dimy1*dimz1*8;
    std::cout << "Length per field: " << field_len1 << std::endl;

	// Skip the first bytes of the header
	//f.seekg (header_size, std::ios::cur);
    size_t len;
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

    // Start writing out files
    if(num_blocks > 1)
    {
        // Write out a header pvti file for the number of blocks
        std::ofstream outp;
        part_filenames[0].assign(filename);
        std::string tempfilename = part_filenames[0].substr(0,len-3);
        tempfilename.append("pvti");
        part_filenames[0] = tempfilename;

        // Begin formatting the head file header
        printf("Creating %s\n", part_filenames[0].c_str());
        outp.open(part_filenames[0].c_str(), std::ifstream::binary);
        outp.precision(20);
        outp << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        outp << "  <PImageData WholeExtent=\""<<0<<" "<<dimx1-1<<" "<<0<<" ";
        outp << dimy1-1<<" "<<0<<" "<<dimz1-1<<"\"";
        outp << " GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
        outp << "    <PPointData Scalars=\"density_scalars\">"<<std::endl;
        for(int i=0;i<scalar_fields1;i++)
        {
            outp << "      <PDataArray type=\"Float64\" Name=\""<<field_names[i]<<"\"/>"<<std::endl;
        }
        outp << "    </PPointData>" << std::endl;
        // Compute the extents for each piece file
        int i, j, b;
        b = num_blocks;
        int * block_dim=new int[3];
        block_dim[0] = 1;
        block_dim[1] = 1;
        block_dim[2] = 1;
        block_size=new int[3];
        block_size[0] = dimx1;
        block_size[1] = dimy1;
        block_size[2] = dimz1;
        
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
        printf("Each block will be- X: %i Y: %i Z: %i\n",block_size[0],block_size[1], block_size[2]);
        printf("# of each block - #X: %i #Y: %i #Z: %i\n", block_dim[0], block_dim[1], block_dim[2]);

        int blockID=0;
        for(int z=0; z<block_dim[2]; z++)
        {
            for(int y=0; y<block_dim[1]; y++)
            {
                for(int x=0; x<block_dim[0]; x++)
                {
                    // Provide 1 level of ghost data to fill in the gaps
                    if(x==0)
                        ext_strt[blockID][0]=x*block_size[0];
                    else
                        ext_strt[blockID][0]=x*block_size[0]-1;

                    if(y==0)
                        ext_strt[blockID][1]=y*block_size[1];
                    else
                        ext_strt[blockID][1]=y*block_size[1]-1;

                    if(z==0)
                        ext_strt[blockID][2]=z*block_size[2];
                    else
                        ext_strt[blockID][2]=z*block_size[2]-1;

                    ext_end[blockID][0]=(x+1)*block_size[0]-1;
                    ext_end[blockID][1]=(y+1)*block_size[1]-1;
                    ext_end[blockID][2]=(z+1)*block_size[2]-1;

                    printf("Extents for block %i are %i %i %i %i %i %i \n",
                        blockID, ext_strt[blockID][0],ext_end[blockID][0],ext_strt[blockID][1],
                        ext_end[blockID][1],ext_strt[blockID][2],ext_end[blockID][2]);
                    blockID++;
                }
            }
        }


        // After extents, begin writing each files metadata into the main head file
        for(int i=0;i<num_blocks;i++)
        {
            // Produce a filename for the piece data
            part_filenames[i+1].assign(filename);
            std::stringstream tempval;
            tempval << i;
            std::string tempfilename = part_filenames[i+1].substr(0,len-4);
            tempfilename.append("_");
            tempfilename.append(tempval.str());
            tempfilename.append(".vti");
            part_filenames[i+1] = tempfilename;

            outp << "    <Piece Extent=\""<<ext_strt[i][0]<<" "<<ext_end[i][0]<<" "<<ext_strt[i][1];
            outp << " "<<ext_end[i][1]<<" "<<ext_strt[i][2]<<" "<<ext_end[i][2];
            outp << "\" Source=\"" << part_filenames[i+1] << "\"/>" << std::endl;
        }
        outp << "  </PImageData>" << std::endl;
        outp << "</VTKFile>" << std::endl;
        outp.close();
    }
    else
    {
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
        ext_end[blockID][0]=dimx1-1;
        ext_end[blockID][1]=dimy1-1;
        ext_end[blockID][2]=dimz1-1;

        // Set block size
        block_size[0]=dimx1;block_size[1]=dimy1;block_size[2]=dimz1;
    }

    // Start creating each individual block
    for(int g=0; g<num_blocks; g++)
    {
      	std::ofstream out;
        // Create the output file
        printf("Creating %s\n", part_filenames[g+1].c_str());
        out.open(part_filenames[g+1].c_str(), std::ifstream::binary);
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
            if(ext_strt[g][0]!=0){ghost[0]=1;}else{ghost[0]=0;}
            if(ext_strt[g][1]!=0){ghost[1]=1;}else{ghost[1]=0;}
            if(ext_strt[g][2]!=0){ghost[2]=1;}else{ghost[2]=0;}
            //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
            
    	    out << "        <DataArray type=\"Float64\" Name=\""<<field_names[i]<<"\" format=\"appended\" offset=\""<< (double)(i*((double)(block_size[0]+ghost[0])*(block_size[1]+ghost[1])*(block_size[2]+ghost[2])*8))+bin_offset <<"\"/>"<<std::endl;
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

        char * bin_value = new char[8];
        std::cout<< "[";
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
            for(int z=ext_strt[g][2]; z<=ext_end[g][2]; z++)
            {
                for(int y=ext_strt[g][1]; y<=ext_end[g][1]; y++)
                {
                    for(int x=ext_strt[g][0]; x<=ext_end[g][0]; x++)
                    {
						double dbl_val = (double)in[x][y];
                        out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(double));
                    }
                }
            }
            std::cout<<"] - "<<field_names[sf]<<" done!\n[";
        }
        printf("Complete!");
        std::cout<<"]\n";
        //Finish off the xml file
        out << std::endl;
        out << "  </AppendedData>"<<std::endl;
        out << "</VTKFile>";

        printf("%s Done!\n", part_filenames[g+1].c_str());
        out.close();
    }
    std::cout<<"Done!"<<std::endl;
    return 0;
}

int readConfigw()
{
    std::ifstream config;
    config.open("config.txt");
    if (config.fail()) {
        std::cout << "Error: Unable to open config.txt" << std::endl;
		return 0;
	}

    pfs1 = new int[6];

    // Start reading file
    std::string line;
    while(config.good())
    {
        std::getline(config,line);
        //std::cout << line << std::endl;
        if(line.find("dimension_size")!=std::string::npos)
        {
            //std::cout<< "dim found" << std::endl;
            size_t delim1 = line.find(' ');
            size_t delim2 = line.find(' ',delim1+1);
            size_t delim3 = line.find(' ',delim2+1);
            size_t delim4 = line.length();
            pfs1[0] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            pfs1[1] = atoi(line.substr(delim2+1,delim3-delim2-1).c_str());
            pfs1[2] = atoi(line.substr(delim3+1,delim4-delim3-1).c_str());
            continue;
        }
        if(line.find("header_size")!=std::string::npos)
        {
            //std::cout<< "header found" << std::endl;
            size_t delim1 = line.find(' ');
            size_t delim2 = line.length();
            pfs1[3] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            continue;
        }
        if(line.find("scalar_fields")!=std::string::npos)
        {
            //std::cout<< "scalar found" << std::endl;
            size_t delim1 = line.find(' ');
            size_t delim2 = line.length();
            pfs1[4] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
            fn1 = new std::string[pfs1[4]];
            continue;
        }
        if(line.find("field_names")!=std::string::npos)
        {
            //std::cout<< "fieldnames found" << std::endl;
            int prev_delim=-1;
            for(int i=0; i<pfs1[4];i++)
            {
                if(i==pfs1[4]-1) // end string
                {
                    size_t delim1=line.find(' ',prev_delim+1);
                    size_t delim2=line.length();
                    fn1[i]=line.substr(delim1+1,delim2-delim1-1);
                }
                else // any other string
                {
                    size_t delim1=line.find(' ',prev_delim+1);
                    size_t delim2=line.find(' ',delim1+1);
                    fn1[i]=line.substr(delim1+1,delim2-delim1-1);
                    prev_delim = delim1;
                }
            }
            continue;
        }
        if(line.find("num_blocks")!=std::string::npos)
        {
            //std::cout<< "number of blocks found" << std::endl;
            size_t delim1 = line.find(' ');
            size_t delim2 = line.length();
            pfs1[5] = atoi(line.substr(delim1+1,delim2-delim1-1).c_str());
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

int arr_to_pvti_3d(float *** in, char * filename,int dx, int dy, int dz)
{
    
    if(!readConfigw())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs1;
    std::string * field_names = fn1;

    dimx1 = dx;//prefs[0];
    dimy1 = dy;//prefs[1];
    dimz1 = dz;//prefs[2];
    header_size1 = prefs[3];
    scalar_fields1 = prefs[4];
    int num_blocks = prefs[5];
    std::cout << "Dimensions "<< dimx1 << " " << dimy1 << " " << dimz1 <<std::endl;
    std::cout << "Header " << header_size1 << " Fields " << scalar_fields1<<std::endl;
    std::cout << "Num Blocks " << num_blocks << std::endl;


    //field_len = (double)((length-header_size) / scalar_fields);
    field_len1 = (double)dimx1*dimy1*dimz1*8;
    std::cout << "Length per field: " << field_len1 << std::endl;

	// Skip the first bytes of the header
	//f.seekg (header_size, std::ios::cur);
    size_t len;
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

    // Start writing out files
    if(num_blocks > 1)
    {
        // Write out a header pvti file for the number of blocks
        std::ofstream outp;
        part_filenames[0].assign(filename);
        std::string tempfilename = part_filenames[0].substr(0,len-3);
        tempfilename.append("pvti");
        part_filenames[0] = tempfilename;

        // Begin formatting the head file header
        printf("Creating %s\n", part_filenames[0].c_str());
        outp.open(part_filenames[0].c_str(), std::ifstream::binary);
        outp.precision(20);
        outp << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        outp << "  <PImageData WholeExtent=\""<<0<<" "<<dimx1-1<<" "<<0<<" ";
        outp << dimy1-1<<" "<<0<<" "<<dimz1-1<<"\"";
        outp << " GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
        outp << "    <PPointData Scalars=\"density_scalars\">"<<std::endl;
        for(int i=0;i<scalar_fields1;i++)
        {
            outp << "      <PDataArray type=\"Float64\" Name=\""<<field_names[i]<<"\"/>"<<std::endl;
        }
        outp << "    </PPointData>" << std::endl;
        // Compute the extents for each piece file
        int i, j, b;
        b = num_blocks;
        int * block_dim=new int[3];
        block_dim[0] = 1;
        block_dim[1] = 1;
        block_dim[2] = 1;
        block_size=new int[3];
        block_size[0] = dimx1;
        block_size[1] = dimy1;
        block_size[2] = dimz1;
        
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
        printf("Each block will be- X: %i Y: %i Z: %i\n",block_size[0],block_size[1], block_size[2]);
        printf("# of each block - #X: %i #Y: %i #Z: %i\n", block_dim[0], block_dim[1], block_dim[2]);

        int blockID=0;
        for(int z=0; z<block_dim[2]; z++)
        {
            for(int y=0; y<block_dim[1]; y++)
            {
                for(int x=0; x<block_dim[0]; x++)
                {
                    // Provide 1 level of ghost data to fill in the gaps
                    if(x==0)
                        ext_strt[blockID][0]=x*block_size[0];
                    else
                        ext_strt[blockID][0]=x*block_size[0]-1;

                    if(y==0)
                        ext_strt[blockID][1]=y*block_size[1];
                    else
                        ext_strt[blockID][1]=y*block_size[1]-1;

                    if(z==0)
                        ext_strt[blockID][2]=z*block_size[2];
                    else
                        ext_strt[blockID][2]=z*block_size[2]-1;

                    ext_end[blockID][0]=(x+1)*block_size[0]-1;
                    ext_end[blockID][1]=(y+1)*block_size[1]-1;
                    ext_end[blockID][2]=(z+1)*block_size[2]-1;

                    printf("Extents for block %i are %i %i %i %i %i %i \n",
                        blockID, ext_strt[blockID][0],ext_end[blockID][0],ext_strt[blockID][1],
                        ext_end[blockID][1],ext_strt[blockID][2],ext_end[blockID][2]);
                    blockID++;
                }
            }
        }


        // After extents, begin writing each files metadata into the main head file
        for(int i=0;i<num_blocks;i++)
        {
            // Produce a filename for the piece data
            part_filenames[i+1].assign(filename);
            std::stringstream tempval;
            tempval << i;
            std::string tempfilename = part_filenames[i+1].substr(0,len-4);
            tempfilename.append("_");
            tempfilename.append(tempval.str());
            tempfilename.append(".vti");
            part_filenames[i+1] = tempfilename;

            outp << "    <Piece Extent=\""<<ext_strt[i][0]<<" "<<ext_end[i][0]<<" "<<ext_strt[i][1];
            outp << " "<<ext_end[i][1]<<" "<<ext_strt[i][2]<<" "<<ext_end[i][2];
            outp << "\" Source=\"" << part_filenames[i+1] << "\"/>" << std::endl;
        }
        outp << "  </PImageData>" << std::endl;
        outp << "</VTKFile>" << std::endl;
        outp.close();
    }
    else
    {
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
        ext_end[blockID][0]=dimx1-1;
        ext_end[blockID][1]=dimy1-1;
        ext_end[blockID][2]=dimz1-1;

        // Set block size
        block_size[0]=dimx1;block_size[1]=dimy1;block_size[2]=dimz1;
    }

    // Start creating each individual block
    for(int g=0; g<num_blocks; g++)
    {
      	std::ofstream out;
        // Create the output file
        printf("Creating %s\n", part_filenames[g+1].c_str());
        out.open(part_filenames[g+1].c_str(), std::ifstream::binary);
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
            if(ext_strt[g][0]!=0){ghost[0]=1;}else{ghost[0]=0;}
            if(ext_strt[g][1]!=0){ghost[1]=1;}else{ghost[1]=0;}
            if(ext_strt[g][2]!=0){ghost[2]=1;}else{ghost[2]=0;}
            //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
            
    	    out << "        <DataArray type=\"Float64\" Name=\""<<field_names[i]<<"\" format=\"appended\" offset=\""<< (double)(i*((double)(block_size[0]+ghost[0])*(block_size[1]+ghost[1])*(block_size[2]+ghost[2])*8))+bin_offset <<"\"/>"<<std::endl;
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

        char * bin_value = new char[8];
        std::cout<< "[";
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
            for(int z=ext_strt[g][2]; z<=ext_end[g][2]; z++)
            {
                for(int y=ext_strt[g][1]; y<=ext_end[g][1]; y++)
                {
                    for(int x=ext_strt[g][0]; x<=ext_end[g][0]; x++)
                    {
						double dbl_val = (double)in[x][y][z];
                        out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(double));
                    }
                }
            }
            std::cout<<"] - "<<field_names[sf]<<" done!\n[";
        }
        printf("Complete!");
        std::cout<<"]\n";
        //Finish off the xml file
        out << std::endl;
        out << "  </AppendedData>"<<std::endl;
        out << "</VTKFile>";

        printf("%s Done!\n", part_filenames[g+1].c_str());
        out.close();
    }
    std::cout<<"Done!"<<std::endl;
    return 0;
}


int arr_to_pvti(std::vector<int*>in,int fields, char * filename,int dx, int dy, int dz)
{
    //Parse config file to obtain field names
    if(!readConfigw())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs1;
    std::string * field_names = fn1;
    scalar_fields1 = prefs[4];
    if(scalar_fields1 < fields)
    {
    	std::cout << "Error: More fields requested than available\n";
        return 0;
    }

	std::cout << "Vector size: " << in.size() << std::endl;
    int num_blocks = (int)in.size();

    size_t len = strlen(filename);
    std::string * part_filenames = new std::string[num_blocks+1];
    // Write out a header pvti file for the number of blocks
    std::ofstream outp;
    part_filenames[0].assign(filename);
    std::string tempfilename = part_filenames[0].substr(0,len-3);
    tempfilename.append("pvti");
    part_filenames[0] = tempfilename;

    // Begin formatting the head file header
    printf("Creating %s\n", part_filenames[0].c_str());
    outp.open(part_filenames[0].c_str(), std::ifstream::binary);
    outp.precision(20);
    outp << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    outp << "  <PImageData WholeExtent=\""<<0<<" "<<dx-1<<" "<<0<<" ";
    outp << dy-1<<" "<<0<<" "<<dz-1<<"\"";
    outp << " GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
    outp << "    <PPointData Scalars=\"density_scalars\">"<<std::endl;
    for(int i=0;i<fields;i++)
    {
        outp << "      <PDataArray type=\"Float64\" Name=\""<<field_names[i]<<"\"/>"<<std::endl;
        //outp << "      <PDataArray type=\"Float64\" Name=\""<<"density_scalars"<<"\"/>"<<std::endl;
    }
    outp << "    </PPointData>" << std::endl;

    // After extents, begin writing each files metadata into the main head file
    for(int i=0;i<num_blocks;i++)
    {
        // parse the name of the file
        std::string filenm;
        filenm.assign(filename);
        size_t found = filenm.find_last_of("/\\");
        // Produce a filename for the piece data
        part_filenames[i+1].assign(filenm.substr(found+1));
        //printf("Parsed filename: %s\n", part_filenames[i+1].c_str());
        len = strlen(part_filenames[i+1].c_str());
        std::stringstream tempval;
        tempval << i;
        std::string tempfilename = part_filenames[i+1].substr(0,len-4);
        tempfilename.append("_");
        tempfilename.append(tempval.str());
        tempfilename.append(".vti");
        part_filenames[i+1] = tempfilename;
        //printf("%s\n",tempfilename.c_str());

        outp << "    <Piece Extent=\""<<in[i][0]<<" "<<in[i][1]<<" "<<in[i][2];
        outp << " "<<in[i][3]<<" "<<in[i][4]<<" "<<in[i][5];
        outp << "\" Source=\"" << part_filenames[i+1] << "\"/>" << std::endl;
    }
    outp << "  </PImageData>" << std::endl;
    outp << "</VTKFile>" << std::endl;
    outp.close();

    return 1;
}

int arr_to_bin_3d_range_scale(double **** in, int num_fields, int header_sz, const char * filename,int xl,int xu, int yl,int yu, int zl, int zu)
{
    if(!readConfigw())
    {
        std::cout << "Error reading configuration file. Exiting..." << std::endl;
        return 0;
    }
    int *prefs = pfs1;
    std::string * field_names = fn1;

    dimx1 = (xu-xl)+1;//prefs[0];
    dimy1 = (yu-yl)+1;//prefs[1];
    dimz1 = (zu-zl)+1;//prefs[2];
    header_size1 = header_sz;//prefs[3];
    scalar_fields1 = num_fields;//prefs[4];
    int num_blocks = 1;//prefs[5]; // This function has removed multi-block output

    
    //std::cout << "Dimensions "<< dimx1 << " " << dimy1 << " " << dimz1 <<std::endl;
    //std::cout << "Header " << header_size1 << " Fields " << scalar_fields1<<std::endl;
    //std::cout << "Num Blocks " << num_blocks << std::endl;


    //field_len = (double)((length-header_size) / scalar_fields);
    field_len1 = (double)dimx1*dimy1*dimz1*8;
    //std::cout << "Length per field: " << field_len1 << std::endl;

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

    int g = 0;
    std::ofstream out;
    // Create the output file
    printf("Creating %s\n", filename);
    out.open(filename, std::ifstream::binary);
    out.precision(20);

    int * ghost = new int[3];

    char * header_padding = new char[1];
    // Write the header
    for(int i=0; i<header_sz; i+=1)
    {
        out.write(reinterpret_cast<char*>( header_padding ), sizeof(char));
    }

    for(int i=0; i<scalar_fields1; i++)
    {
        // Factor in ghost data into the blocks
        if(ext_strt[g][0]!=0){ghost[0]=0;}else{ghost[0]=0;}
        if(ext_strt[g][1]!=0){ghost[1]=0;}else{ghost[1]=0;}
        if(ext_strt[g][2]!=0){ghost[2]=0;}else{ghost[2]=0;}
        //printf("%i %i %i\n", ghost[0], ghost[1], ghost[2]);
    }

    char * bin_value = new char[8];
    //std::cout<< "[";
    for(int sf=0; sf < scalar_fields1; sf++)
    {
        int c=0;
        for(int z=ext_strt[g][2]; z<=ext_end[g][2]; z++)
        {
            int b=0;
            for(int y=ext_strt[g][1]; y<=ext_end[g][1]; y++)
            {
                int a=0;
                for(int x=ext_strt[g][0]; x<=ext_end[g][0]; x++)
                {
                    double dbl_val = in[sf][a][b][c];
                    out.write(reinterpret_cast<char*>( &dbl_val ), sizeof(double));
                    a++;
                }
                b++;
            }
            c++;
        }
    }
    //printf("Complete!");
    //std::cout<<"]\n";

    out.close();
    //std::cout<<"Done!"<<std::endl;
    return 0;
}

int * coeff_args(unsigned short type, unsigned short lvl, unsigned char reg, unsigned short pad, int lx, int ly, int lz, int gx, int gy, int gz, char prec)
{
    int * args = new int[11];
    args[0] = type; args[1] = lvl;
    args[2] = reg; args[3] = pad;
    args[4] = lx; args[5] = ly;
    args[6] = lz; args[7] = gx;
    args[8] = gy; args[9] = gz;
    args[10] = prec;

    // debug output
    std::cout<< "Converting parameters into array:"<< std::endl;
    std::cout << "Type: " << args[0] << " Level: " << args[1] << std::endl;
    std::cout << "Region: " << args[2] << " Padding: " << args[3] << std::endl;
    std::cout << "Local dims: " << args[4] << " " << args[5] << " " << args[6] << std::endl;
    std::cout << "Global dims: " << args[7] << " " << args[8] << " " << args[9] << std::endl;
    std::cout << "Data precision: " << args[10] << " bytes." << std::endl;

    return args;
}

int free_args(int * args)
{
	if(args)
		delete [] args;

    return 1;
}

int save_sort(size_t * sort, const char * filename, size_t total)
{
    std::ofstream out;
    // Create the output file
    printf("Creating %s\n", filename);
    out.open(filename, std::ifstream::binary);
    out.precision(20);

    for(int sz=0; sz < total; sz++)
    {
        unsigned int uint_val = sort[sz];
        out.write(reinterpret_cast<char*>( &uint_val ), sizeof(unsigned int));
    }
    std::cout<<"Sorted Positions saved\n";

    out.close();
    return 1;
}
