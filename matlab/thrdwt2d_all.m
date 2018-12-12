function [wtPrime] = thrdwt2d_all(varargin)
%function [wtPrime] = thrdwt2d(wt,sizes, pcnt )
%   wt = wavelet coefficients
%   sizes = array of coefficient sizes
%   pcnt = percent of coefficients to threshold inside of thr_lvl
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu
 
    wt = varargin{1};
    sizes = varargin{2};
    pcnt = varargin{3};

    tot = numel(wt);
    
    if(pcnt > 0)

        % Sort Coefficients 
        [~,index] = sort(abs(wt),2,'descend');
        
        % Define threshold amount
        thramount = (tot*(pcnt/100.0));
        fprintf('Thresholding %12.8f of %12.8f at %f pcnt\n', thramount,tot,pcnt);

        % Threshold coefficients
        for i = 1:tot
            if i > thramount
                wt(:,index(i)) = 0;
            end
        end
    end

    %sum(sum(abs(wt)))
    wtPrime = wt;
    
end

