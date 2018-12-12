function [] = arr3dtoqbin(filenm,array,x,y,z,header,field)
% ARR3DTOQBIN  Writes out data array to a binary file with quantization
%
%   filenm = File target location
%   array = 1,2,3 dimensional data array
%   x = size of dim1
%   y = size of dim2
%   z = size of dim3
%   header = header padding in # bytes
%   field = number of present fields (multi-variate data)
%
%   See also QBINTO3D
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

a = zeros(x,y,z);

file = fopen(filenm,'w');

bytesheader = int32(header);
%insert header
dummy = zeros(5);

for c = 0:header 
    fwrite(file,dummy(1),'double');
end

%Move pointer ahead of header
fseek(file,bytesheader,'bof');
%b = fread(file,1,'double')

%Move pointer to appropriate field
bytesheader = bytesheader + (field*x*y*z*8);

%for m = 1:x
%    for n = 1:y
%        for p = 1:z
%            %fseek(file,bytesheader,'bof');
%            %a(m,n,p) = fread(file,1,'double');
%            fwrite(file,array(m,n,p),'double');
%            bytesheader = bytesheader + 8;
%        end
%    end
%end

%Try binary write
a = reshape(array,x*y*z,1);
%fwrite(file,a,'*float');

tic
bits=100;
%range(a)
a=a*bits; %Controls number of bits for integer
%range(a)
y = a-rem(a,1);
y = int32(y);
toc

fwrite(file,y,'*int32');

fclose(file);
