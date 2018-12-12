% TEST_COMPRESS tests compression routine on input file
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

close all;

%Load in Fits file
%FileName='D:\Data\LSST\data\lsst_e_898670970_f0_R02_S00_E000.fits';
FileName='data\Deep_32.fits';
f=dir(FileName);
fsizeMB = f.bytes/1024/1024; %Query original file size

if exist(strcat(FileName,'.gz'), 'file') ~= 2
    gzip(FileName); %GZIP Fits file if not already
end

f2=dir(strcat(FileName,'.gz'));
f2sizeMB = f2.bytes/1024/1024; %Query compressed GZIP file size


ShowImage=0; %Enable when testing with 1 image, if not disable
ori_im=fitsread(FileName); 

if ShowImage figure,imshow(histeq(sqrt(mat2gray(ori_im))),[]),title('original image'); end
 
% Perform noise removal on original image (grid-based) if its from LSST
% since it lacks stacking of images
do_noise_removal=0;
if do_noise_removal==1
    H = fspecial('gaussian',20,3.0);
    load('gauss_mask.mat');
    J = roifilt2(H,ori_im,gauss_mask);
    if ShowImage figure,imshow(histeq(sqrt(mat2gray(J))),[]),title('denoised image'); end
    ori_im = J;
end

wave_type = 'bior4.4';
max_levels = 8;

do_pcnt=1;
if do_pcnt==1
    %Percentage-based thresholding
    for pcnt=5:5:100

        % Wave operation
        tic
        [coeffs,sizes] = wavedec2(ori_im,max_levels,wave_type);
        toc

        % Threshold
        tic
        coeffs2 = thrdwt2d_all(coeffs,sizes,pcnt);
        toc
        
        str_num = sprintf('%03d',pcnt);
        out_name = strcat('deep32_Q32coeffs_',str_num,'.bin');
        
        % Write Compressed output
        LWcompress(out_name,coeffs2,sizes);
        
        continue;
        
    end
end

do_level=0;
if do_level==1
    %Level-based thresholding
    for levels=1:max_levels
        [coeffs,sizes] = wavedec2(ori_im,max_levels,wave_type);
        coeffs2 = thrdwt2d(coeffs,sizes,levels);

        %Quantize and write coefficients to disc
        str_num = sprintf('%03d',levels);
        out_name = strcat('deep32_Q32coeffs_lvl_',str_num,'.bin');
        
        % Quantized output
        if exist(out_name, 'file') ~= 2
            arr3dtoqbin(out_name,coeffs2,numel(coeffs2),1,1,0,1)
        end
        
        % Compress output
        
        
        
        continue;
    end
end