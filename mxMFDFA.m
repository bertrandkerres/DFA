function [ H, F ] = mxMFDFA( data, scale, varargin )
%DFA Detrended fluctuation analysis using compiled C++ mex code
%
%   Hq = DFA (data, scale) performs detrended fluctuation on data; data
%   must be a vector; scale is a row vector given in samples
%   Hq = DFA (data, scale, [q1 q2 ... qn]) calculates the Hurst
%       exponents Hq for different q; default: q = 2
%   Hq = DFA (data, scale, [q1 q2 ... qn], m) uses a different order of
%       detrending; default: m = 1
%   Hq = DFA (data, scale, [q1 q2 ... qn], integrate)
%       integrates the signal before running DFA; allowed = -1 (differentiate), 0, 1, 2;
%       default = 1
    
    ip = inputParser ();
    ip.addRequired ('data', @(x) isvector(x) && isnumeric(x));
    ip.addRequired ('scale', @(x) all(x>0) && size(x,1) == 1);
    ip.addOptional ('q', 2, @(x) isnumeric(x));
    ip.addOptional ('m', 1, @(x) isscalar(x) && x >= 0);
    ip.addOptional ('integrate', 1, @(x) isscalar(x) && (x==-1 || x==0 || x==1 ||x==2));
    ip.addParameter ('R2warn', 0, @(x) isscalar(x) && x>= 0 && x <= 1);
    ip.parse (data, scale, varargin{:});
    
    qvec = double(ip.Results.q(:));

    data = data - mean (data);
    if isrow (data)
        data = data';
    end

    H = NaN * zeros (length(qvec), 1);
    
    % Add another integration
    if ip.Results.integrate == 2
        data = cumsum (data);
        data = cumsum (data-mean(data));
    elseif ip.Results.integrate == 1
        data = cumsum (data);
    elseif ip.Results.integrate == -1
        data = [0; diff(data)];
        data = data - mean(data);
    end
    
    % Call mex function to calculate the structure function of the data
    F = mxDFA(data, uint32(scale), qvec, uint32(ip.Results.m));
       
    % Calculate H by fitting lines on a log-log plot
    for qc = 1 : length(qvec)
        logF = log(F(qc,:));
        [C, S] = polyfit (log(scale), logF, 1);
        H(qc) = C(1) - ip.Results.integrate + 1;
        SS_tot = sum( (logF - mean(logF)).^2 );
        R2 = 1 - S.normr^2 / SS_tot;
        if R2 <= ip.Results.R2warn
            warning ('q = %f: Linear fit of flucuation function has an R2 value of %f', qv(qc), R2);
        end
    end

end    