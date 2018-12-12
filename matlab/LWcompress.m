function [] = LWcompress(filenm,data, sizes)
% LWDECOMPRESS  Compresses file
%
%   filenm = File target location to be written in 7z format
%   data = data array
%   sizes = sizes array output from wavelet routine (temporary)
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

%   Example:
%   filenm = 'D:\my_quantized.bin';
%   [] = LWcompress(filenm);

% Quantize the image first
save('deep32_sizes.mat','sizes'); %NOTE: This is a temporary solution, will not need file in the future
fprintf('Quantizing...');
if exist(filenm, 'file') ~= 2
    arr3dtoqbin(filenm,data,numel(data),1,1,0,1)
end

% separates filename into path/file/extension
[dir7z,~,~] = fileparts( filenm );
dir7z = ['"' dir7z '"']; %get Path
fileout = ['"' strcat(filenm,'.lz4') '"'];

% Download 7-zip Zstandard from here: https://github.com/mcmilk/7-Zip-zstd/releases
% Create archive with lz4, may default to Zstandard
fprintf('Encoding...');
tic
[status,cmdout] = system( ['"C:\Program Files\7-Zip-Zstandard\7z.exe" a -t7z ',fileout,' ',filenm] );
toc

if status ~= 0
    cmdout
    return
end

% Remove temporary bin file
if exist(filenm, 'file')==2
    delete(filenm)
end

fprintf('...Done!\n');

end
