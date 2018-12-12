function a = qbinto3d(filenm,x,y,z,header,field)
% QBINTO3D  Reads a binary file into a matlab data array
%
%   filenm = File target location
%   x = size of dim1
%   y = size of dim2
%   z = size of dim3
%   header = header padding in # bytes
%   field = number of present fields (multi-variate data)
%
%   See also ARR3DTOQBIN
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

a = zeros(x,y,z);
a = double(a);

file = fopen(filenm,'r');

fieldlocation = field*x*y*z*8;

bytesheader = int64(header+fieldlocation);

fseek(file,bytesheader,'bof');
%b = fread(file,1,'double')


%for m=1:x
%    for n = 1:y
%        for p = 1:z
%            fseek(file,bytesheader,'bof');
%            a(m,n,p) = fread(file,1,'double');
%            bytesheader = bytesheader + 8;
%        end
%    end
%end

%Try binary read
a = fread(file,(x*y*z),'*int');
a = reshape(a,x,y,z);

bits=100; %enhance by 'N' bits. i.e. 1000=> +3 bits.
tic
% Quantize by n digits
%tcoeffs2=a*bits; %Bit control
%icoeffs2 = tcoeffs2-rem(tcoeffs2,1);
%icoeffs2 = int32(icoeffs2);
dcoeffs2 = double(a);
a = dcoeffs2/bits; %Bit control
%clear tcoeffs2;
toc
        
fclose(file);
