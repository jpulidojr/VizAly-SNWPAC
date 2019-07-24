% TEST_DECOMPRESS tests decompression routine on input file
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

close all;

%Load in Fits file
%FileName='data\Deep_32.fits';
FileName='data\lsst_e_898670970_f0_R02_S00_E000_stacked.fits';

f=dir(FileName);
fsizeMB = f.bytes/1024/1024; %Query original file size

if exist(strcat(FileName,'.gz'), 'file') ~= 2
    gzip(FileName); %GZIP file if not already, to compare sizes
end

f2=dir(strcat(FileName,'.gz'));
f2sizeMB = f2.bytes/1024/1024; %Query compressed GZIP file size

ShowImage=1; %Enable when testing with 1 image, if not disable
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
    %for pcnt=5:5:100
    for pcnt=100
        %load coeffs from disc
        str_num = sprintf('%03d',pcnt);
        out_name = strcat('.\\lsst_HQ32coeffs_',str_num,'.bin.lz4'); %Load the fast LZ4 version
        
        fprintf('Processing file %s \n',out_name);
        % Quantized output
        if exist(out_name, 'file') == 2
            imgout = SPdecompress(out_name);
            if (imgout == 0)
               fprintf('Error: Decompression failed!\n')
               continue;
            end
        else
            fprintf('Error: Compressed file not found at %d pcnt!\n', pcnt);
            continue;
        end
        
        out_im = histeq(sqrt(mat2gray(imgout)));
        out_imdiff = histeq(sqrt(mat2gray(imgout-ori_im)));


        if ShowImage figure,imshow(mat2gray(out_im),[]),title(strcat('lossy compressed at ',num2str(pcnt),' pcnt (INT)')); end
        %imwrite(out_im,strcat('C:/tmp/',wave_type,'_',str_num,'.png'))
        
        
        %if ShowImage figure,imshow(mat2gray(out_imdiff),[]),title('Orig v Comp Error'); end

        continue;
        
        %Test feature extraction code here with Zheng et. al. method
        % https://github.com/zhengcx789/Object-Detection-in-Astronomical-Images
        
        %[Final_detection_result,Final_Label,Final_statistics]=mainGlobalDetection(imgout, ShowImage,str_num);
        %fprintf('Objects Detected: %d\n', size(Final_statistics,2));

        %pause
    end
end
