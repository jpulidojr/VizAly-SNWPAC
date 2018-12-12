function [image] = LWdecompress(file7z)
% LWDECOMPRESS  Decompresses file
%
%   file7z = File target location in 7z format
%
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

%   Example:
%   file7z = 'D:\my.7z';
%   myImage = LWdecompress(file7z);

image = 0;

wave_type = 'bior4.4';

[dir7z,~,~] = fileparts( file7z );
dir7z = ['"' dir7z '"'];
file7z = ['"' file7z '"'];

% Read contents of archive and list
% -ba command requires 7z v15.xx or greater
% Download 7-zip Zstandard from here: https://github.com/mcmilk/7-Zip-zstd/releases
% This has support for lz4 and Bzip2
[status,cmdout] = system( ['"C:\Program Files\7-Zip-Zstandard\7z.exe" l -ba ',file7z] );

delim_cmd = strsplit(cmdout,'\n');
num_files = size(delim_cmd,2);

for fl = 1:num_files
    delim_file = strsplit(delim_cmd{fl});
    if size(delim_file,2) > 1
        file = delim_file{size(delim_file,2)}; % this should be the coefficient bin file
    end
end

if status ~= 0
    cmdout
    return
end

% Extract archive
fprintf('Decoding...');
tic
[status,cmdout] = system( ['"C:\Program Files\7-Zip-Zstandard\7z.exe" x -y -o',dir7z,' ',file7z] );
toc

if status ~= 0
    cmdout
    return
end

% Read in coefficient file and reconstruct
% Read Quantized input
if exist('deep32_sizes.mat', 'file') == 2
    %load('sizes.mat'); %NOTE: This is a temporary solution, will not need file in the future
    load('deep32_sizes.mat'); %NOTE: This is a temporary solution, will not need file in the future
    %%%num_coeff = sum(sizes(:,1).*sizes(:,2))
    %coeffs2 = qbinto3d(file,16418841,1,1,0,1); % (LSST) Number of coeffs will be dynamic later
    coeffs2 = qbinto3d(file,16910127,1,1,0,1); % (DLS Data)
    fprintf('..Reconstructing..');
    image = waverec2(coeffs2,sizes,wave_type);
else
    fprintf('Error: Sizes file missing!\n');
end

fprintf('...Done!\n');

end
