// LWDecomp.cpp
// Author: Jesus Pulido
// Simple LW decompressor with several 3D file formats: Bin, HDF5, vti.
// Note: HDF5 support must be enabled via cmake to work.
// 
// Usage: ./LWDecomp {input} {output}
//	  	  {input} = input header file (text)
//	  	  {output} = output filename (.bin | .h5 | .vti)

#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "lossywave.hpp"

#include <omp.h> //Only for timings

using namespace std;

#ifdef USE_HDF5
	#include "H5Cpp.h"
#endif

int main (int argc, char **argv)
{
	bool verbose = false;

	bool writeh5=false;
	bool writevti=false;
	bool writebin=false;

	string outfile;
	string infile;
	if(argc == 2)
    {
        cout << "Input file is " << argv[1] << endl;
		cout << "Output will be auto-generated" << endl;

		outfile = strcat(argv[1],".bin");
		writebin=true;
    }
    else if(argc == 3)
    {
		cout << "Input file is " << argv[1] << endl;
		cout << "Output file is " << argv[2] << endl;
		outfile = argv[2];

		if(outfile.find(".vti")!=std::string::npos)
		{
			writevti=true;
		}
		else if(outfile.find(".h5")!=std::string::npos)
		{
			writeh5=true;
		}
		else
			writebin=true;
    }
	else
	{
		cout << "Usage: LWDecomp <input_header> <output> " << endl;

		return 0;
	}
	infile = argv[1];

	std::cout << "--Write Mode: ";
	if(writeh5)
		cout << "hdf5";
	if(writevti)
		cout << "vti";
	if(writebin)
		cout << "binary";
	cout << endl;

	//argv[1]=infile will be a header file
	vector<string> field_names;
	vector<string> paths;

	ifstream input(infile.c_str());
	if(!input.is_open())
	{
		cout << "Error: Unable to open file " << infile << endl;
		return 0;
	}

	string sfields;
	std::getline(input, sfields);
	int fields = atoi(sfields.c_str()); 

	// Read in the header file and populate metadata
	cout << "Detected: " << fields << " fields." << endl;
	for(int fld=0; fld<fields; fld++)
	{
		string sname;
		std::getline(input, sname);
		field_names.push_back(sname);
		string spath;
		std::getline(input, spath);

		// Concatenate path with root dir of header file
		string directory;
		size_t last_slash_idx = infile.rfind('\\');
		// Try backslash (windows)
		if(string::npos != last_slash_idx)
			directory = infile.substr(0, last_slash_idx+1);
		else
		{  //Try forward slash (linux)
			last_slash_idx = infile.rfind('/');
			if(string::npos != last_slash_idx)
				directory = infile.substr(0, last_slash_idx+1);
		}
		string finalpath = directory + spath;
		// -----------

		paths.push_back(finalpath);
		cout << "Field: " << sname << endl;
	}
	cout << endl;
	//fields=1; //overwrite for now


#ifdef USE_HDF5
	H5::H5File file;
	H5::Group group;
#endif

	int x_w,y_w,z_w;
	float ** outdata = new float*[fields];
	for(int fld = 0; fld<fields; fld++)
	{
		std::ifstream f;
		// Try opening the file
		f.open(paths[fld].c_str(), std::ifstream::binary);
		if (f.fail()) {
			std::cout << "Error opening file: " << paths[fld].c_str() << std::endl;
			return 0;
		}

		// Read in data stream from file
		void * compressed;
		//float * coeffs = read_coefficients_md<float>(paths[fld].c_str(), args);

		// Process arguments
		int args[13];// = { 404, 0, 128 + argv_quant, 0,
					//dims[0], dims[1], dims[2],
					//dims[0], dims[1], dims[2],
					//sizeof(input3d[0]), argv_pcnt, 0 };

		lossywave::lossywave lw(args, true);
		lw.printParams();

		/*int * args = read_md(paths[fld].c_str());
		if(args==0)
		{
			std::cout<<"Error: Could not open file: " << paths[fld] << std::endl;
			return 0;
		}
		else
		{
			std::cout<< "Header Metadata read from coefficients:"<< std::endl;
			std::cout << "Type: " << args[0] << " Level: " << args[1] << std::endl;
			std::cout << "Region: " << args[2] << " Padding: " << args[3] << std::endl;
			std::cout << "Local dims: " << args[4] << " " << args[5] << " " << args[6] << std::endl;
			std::cout << "Global dims: " << args[7] << " " << args[8] << " " << args[9] << std::endl;
			std::cout << "Data precision: " << args[10] << " bytes." << std::endl;
		}*/
		x_w=args[7];
		y_w=args[8];
		z_w=args[9];
		size_t total = x_w*y_w*z_w;
		size_t ogSize = total * args[10];

		if( fld == 0 ) //&& Rank == 0
		{
			if(writeh5)
			{
#ifdef USE_HDF5
			file= H5::H5File( outfile.c_str(), H5F_ACC_TRUNC );

			int size_dims[3] = {x_w,y_w,z_w};
			double size_ph[3] = {10, 10, 10};
		

			// These should be extracted from original /universe/attr
			double omega_b = 0.047;
			double omega_m = 0.3;
			double omega_l = 0.7;
			double h = 0.685;
			double z = 4.2; //This may change per file

			if(x_w == 512)
				z = 4.2;
			if(x_w == 2048)
				z = 3.0;

			//std::string format_text = "nyx-lyaf";
			H5std_string format_text("nyx-lyaf");


			//hid_t str_type = H5Tcopy(H5T_C_S1);
			//hid_t scalar_space = H5Screate(H5S_SCALAR);
			H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
			H5::DataSpace scalar_space(H5S_SCALAR);

			// Fix type length for string format
			//H5Tset_size(str_type, strlen(format_text.c_str()));
			H5::Attribute attr = file.createAttribute("format",str_type,scalar_space);
			attr.write(str_type,format_text);
			attr.close();

			// Domain dataspace
			hsize_t dims_sz[1] = { 3 };
			H5::DataSpace dataspace(1,dims_sz);
			H5::DataSpace dataspace2(1,dims_sz);	

			// Domain metadata group
			H5::Group pregroup(file.createGroup("/domain"));
			attr = pregroup.createAttribute( "shape", H5::PredType::NATIVE_INT, dataspace);
			attr.write( H5::PredType::NATIVE_INT, size_dims);
			attr.close();

			attr = pregroup.createAttribute( "size", H5::PredType::NATIVE_DOUBLE, dataspace2);
			attr.write( H5::PredType::NATIVE_DOUBLE, size_ph);
			attr.close();
	
			pregroup.close();
			// Universe metadata group
			H5::Group unigroup(file.createGroup("/universe"));
			//group = H5Gcreate(file, "universe", H5P_DEFAULT, H5P_DEFAULT,
			//              H5P_DEFAULT);
		
			attr = unigroup.createAttribute( "omega_b", H5::PredType::NATIVE_DOUBLE, scalar_space);
			attr.write(  H5::PredType::NATIVE_DOUBLE, &omega_b);
			attr.close();
			attr = unigroup.createAttribute( "omega_m", H5::PredType::NATIVE_DOUBLE, scalar_space);
			attr.write(  H5::PredType::NATIVE_DOUBLE, &omega_m);
			attr.close();
			attr = unigroup.createAttribute( "omega_l", H5::PredType::NATIVE_DOUBLE, scalar_space);
			attr.write(  H5::PredType::NATIVE_DOUBLE, &omega_l);
			attr.close();
			attr = unigroup.createAttribute( "hubble", H5::PredType::NATIVE_DOUBLE, scalar_space);
			attr.write(  H5::PredType::NATIVE_DOUBLE, &h);
			attr.close();
			attr = unigroup.createAttribute( "redshift", H5::PredType::NATIVE_DOUBLE, scalar_space);
			attr.write(  H5::PredType::NATIVE_DOUBLE, &z);
			attr.close();
		
			unigroup.close();

			// Field groups
			group = H5::Group(file.createGroup("/native_fields"));
#endif
			}
		}


		std::cout << "Decompressing field: " << fld+1 << " of " << fields << std::endl;
		// Allocate memory for decompression
		void * decompressed;
		decompressed = std::malloc(ogSize);

		double startcopy = omp_get_wtime();

		size_t dcmpSize = lw.decompress(compressed, decompressed);

		double endcopy = omp_get_wtime();
		std::cout << "Decompression time = " << endcopy-startcopy << " secs." <<std::endl;

		// Check if data is valid
		if (ogSize == dcmpSize)
			std::cout << "Data Verified!" << std::endl;
		else
			std::cout << "Invalid Data! Compression Error." << std::endl;

		float * output3d = static_cast<float *>(decompressed);
		// Print a few to test output
		if(verbose)
		{
			for(size_t id=0; id<10;id++)
				std::cout << output3d[id] << " ";
			std::cout << std::endl;
		}
		
		//NOTE: For hdf5 output, do not reformat to 3D

		//outdata[fld] = arrayto3D(coeffs,x_w,y_w,z_w);
		//delete [] coeffs;
		outdata[fld] = output3d;

		if(writeh5)
		{
#ifdef USE_HDF5
			hsize_t dims[3];
			dims[0]=x_w;
			dims[1]=y_w;
			dims[2]=z_w;

			double startwrite = omp_get_wtime();
			// Copy data into HDF5 format
			cout << "Beginning to write data to h5.." << endl;
			H5::DataSpace dataspace(3,dims);
			H5::DataSet dataset = file.createDataSet("/native_fields/"+field_names[fld], H5::PredType::NATIVE_FLOAT, dataspace);
			dataset.write(outdata[fld], H5::PredType::NATIVE_FLOAT);

			dataset.close();
			cout << "Finished writing: " << field_names[fld] << endl;
			double endwrite = omp_get_wtime();
			std::cout << "H5 Write time = " << endwrite-startwrite << " secs." <<std::endl;

			// Keep data in memory if we want to write a vti file later
			if(!writevti)
				delete [] outdata[fld];
#endif
		}
	}

	
	if(writeh5)
	{
#ifdef USE_HDF5
		group.close();

		H5::Group group2(file.createGroup("/derived_fields"));
		group2.close();

		file.close();

		std::cout << "Finishing file: " << outfile.c_str() << std::endl;
	/*	H5::H5File file( outfile, H5F_ACC_TRUNC );
	H5::Group pregroup(file.createGroup("/domain"));


	hsize_t dims[3];
	dims[0]=512;//pDims[0];
	dims[1]=512;//pDims[1];
	dims[2]=512;//pDims[2];

	int size_dims[3] = {512,512,512};
	double size_ph[3] = {10, 10, 10};
	hsize_t dims_sz[1] = { 3 };

	H5::DataSpace dataspace(1,dims_sz);
	H5::DataSpace dataspace2(1,dims_sz);	
	
	H5::Attribute att1 = pregroup.createAttribute( "shape", H5::PredType::NATIVE_INT, dataspace);
	att1.write( H5::PredType::NATIVE_INT, size_dims);
	att1.close();

	H5::Attribute att2 = pregroup.createAttribute( "size", H5::PredType::NATIVE_DOUBLE, dataspace2);
        att2.write( H5::PredType::NATIVE_DOUBLE, size_ph);
	att2.close();
	
	pregroup.close();
	H5::Group group(file.createGroup("/native_fields"));

	for(int fld=0; fld<fields; fld++)
	{
		H5::DataSpace dataspace(3,dims);
		H5::DataSet dataset = file.createDataSet("/native_fields/"+field_names[fld], H5::PredType::NATIVE_FLOAT, dataspace);

		dataset.write(outdata[fld], H5::PredType::NATIVE_FLOAT); //Otherwise its H5::PredType::NATIVE_DOUBLE
		dataset.close();

	}

	group.close();
	file.close();

        std::cout << "Writing out file: " << outfile << std::endl;*/

#endif
	}


	// Copy data into VTI format
	if(writevti)
	{
		double startwrite = omp_get_wtime();
		//Convert to 3d;
		float **** output3d = new float***[fields];
		for(int fld=0;fld<fields;fld++)
		{
			//output3d[fld] = arrayto3D(outdata[fld],x_w,y_w,z_w);
			delete outdata[fld];
		}
		//Write 1 for now
		//arr_to_vti_3d_range_scale(output3d,outfile.c_str(),0,x_w-1, 0,y_w-1,0, z_w-1, field_names);
		double endwrite = omp_get_wtime();
		std::cout << "VTI Write time = " << endwrite-startwrite << " secs." <<std::endl;
	}


	return 0;
}
