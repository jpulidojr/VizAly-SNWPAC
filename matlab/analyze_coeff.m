function [wtPrime] = analyze_coeff(varargin)
%   wt = wavelet coefficients
%   sizes = array of coefficient sizes
%   pcnt = percent of coefficients to threshold inside of thr_lvl
 
% Notes: Deriving output coeffs
% total = numel(coeffs) = 16910127;
% qprod = (prod(sizes,2));
% sum(3*qprod(1:9))-(2*576)

    wt = varargin{1};
    sizes = varargin{2};
    quant=0;
    if(length(varargin)==3)
        qhier = varargin{3};
        wtPrime = int32(wt);
        quant=1;
    else
        qhier = ones(size(sizes,1),1);
        wtPrime = wt;
    end

    qhier
    
    tot = numel(wt);
    
    qprod = prod(sizes,2);
    
    end_bound = tot;
    start_bound=1;
    
    %Note: Ignore highest level, that represents og dimension size
    for i=size(sizes,1)-1:-1:1
    
        % Extract coefficients at res level
        numel_quad = qprod(i);
        
        if i==1
            start_bound = 1; %were at LL
        else
            start_bound = end_bound-3*(numel_quad); %Skip LH,HL,HH
        end
        
        if(quant==1)
            a = wt(start_bound:end_bound)*qhier(i,1);
            y = a - rem(a,1);
            y = int32(y);
        
            wtPrime(start_bound:end_bound) = y;
        end
        
        %analyze coefficient ranges
        qmin = min(min(double(wtPrime(start_bound:end_bound))/qhier(i,1)));
        qmax = max(max(double(wtPrime(start_bound:end_bound))/qhier(i,1)));
        qmean = mean(mean(double(wtPrime(start_bound:end_bound))/qhier(i,1)));
        
        fprintf("At lvl %d sizes:(%d %d), min: %f  max: %f mean: %f\n",i,sizes(i,1),sizes(i,2),qmin,qmax,qmean);
        
        end_bound = start_bound;
        % Define threshold amount
        %thramount = (tot*(pcnt/100.0));
        %fprintf('Thresholding %12.8f of %12.8f at %f pcnt\n', thramount,tot,pcnt);

        % Threshold coefficients
        %for i = 1:tot
        %    if i > thramount
        %        wt(:,index(i)) = 0;
        %    end
        %end
    end

    %sum(sum(abs(wt)))
    
    
end

