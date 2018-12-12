function [wtPrime] = thrdwt2d(varargin)
%function [wtPrime] = thrdwt2d(wt,sizes, thr_lvl, [pcnt] )
%   wt = wavelet coefficients
%   sizes = array of coefficient sizes
%   thr_lvl = level to threshold coefficients (includes lower levels)
%   pcnt = percent of coefficients to threshold inside of thr_lvl
%
%   Author: Jesus Pulido
%           jpulido@ucdavis.edu

    if(length(varargin) <= 3)
        pcnt = 0; 
    else
        pcnt = varargin{4};
        if pcnt >=1
            pcnt = 0;
        end
    end
    wt = varargin{1};
    sizes = varargin{2};
    thr_lvl = varargin{3};

    s_min = (size(sizes,1)-2);
    
    if thr_lvl == 0
        wtPrime = zeros(numel(wt),1);
    end
    
    sprod = prod(sizes,2);
    
    cur_lvl=1;
    zero_start = sprod(cur_lvl);
    
    zero_end = numel(wt);
    thr_lvl=thr_lvl-1;
    while(thr_lvl>0)
        cur_lvl=cur_lvl+1;
        
        if(cur_lvl >= size(sizes,1))
            %if ~exists(lvl_stop)
            %   lvl_stop = cur_lvl
            %end
            thr_lvl=thr_lvl-1;
            continue;
        end
        
        zero_start=zero_start+(sprod(cur_lvl)*3);
        thr_lvl=thr_lvl-1;
    end
    
    if(pcnt > 0)
        % Find beginning and end of the current level
        if (cur_lvl>1)
            prev_start=zero_start-(sprod(cur_lvl)*3);
        else
            prev_start=1;
        end
        tot = numel(wt(:,prev_start:zero_start));
        
        % Sort Coefficients for this level
        [~,index] = sort(abs(wt(:,prev_start:zero_start)),2,'descend');
        
        % Define threshold amount
        thramount = (tot*pcnt);

        % Threshold coefficients
        for i = 1:tot
            if i > thramount
                pos = prev_start+index(i)-1;
                wt(:,pos) = 0;
            end
        end
    end

    if(zero_start < size(wt,2))
        wt(:,zero_start:zero_end)=0;
    end
    %sum(sum(abs(wt)))
    wtPrime = wt;
    
end

