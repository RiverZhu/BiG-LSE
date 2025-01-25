classdef LSEOpt
    
    properties
        
        adaptStep = 1; 
        step_default = 0.8;
        stepMin = 0.5;
        stepMax = 0.9;
        stepWindow = 1;
        stepIncr = 1.15;
        stepDecr = 0.80;
        failCount = 0;
        maxBadSteps = 10; % 5
        
        iter_max = 50;
        Iter_ref_max = 3;
        Iter_in_max = 20;
        
        ref_loop = 1;
        debug = 0;
        N_len = 512;
        L_len = 512;
        Bit = inf;
        sigma_w;
        y_min;
        alpha;
        zeta;
    end
    
    methods
        
        function init = LSEOpt(varargin)
            if nargin==0
                return
            else
                names = fieldnames(init);
                for i = 1:2: nargin-1
                    if any(strcmpi(varargin{i}, names))
                        propName = names(strcmpi(varargin{i}, names));
                        init.(propName{1}) = varargin{i+1};
                    end
                end
                return
            end
        end
        
    end
    
end