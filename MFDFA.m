function [ H, F ] = MFDFA( data, scale, varargin )
%MFDFA Multifractal Detrended Fluctuation Analysis for a dataset
%
%   Hq = MFDFA (data, scale) calculates the Hurst exponent
%       for data; scale is a row vector given in samples
%   Hq = MFDFA (data, scale, [q1 q2 ... qn]) calculates the Hurst
%       exponents Hq for different q; default: q = 2
%   Hq = MFDFA (data, scale, [q1 q2 ... qn], m) uses a different order of
%       detrending; default: m = 1
%   Hq = MFDFA (data, scale, [q1 q2 ... qn], m, integrate)
%       integrates the signal before running DFA; allowed = -1 (differentiate), 0, 1, 2;
%       default = 1
%   Hq = MFDFA (..., 'R2warn', R2) displays a warning if the fit of the
%       fluctuation function to scales has an R2 value below the given
%       threshold. Note that R2 is very sensitive to errors if the slope is
%       close to zero. Default: R2warn = 0, i.e. never warn.
%   [Hq, Fq] = MFDFA (...) returns also the fluctuation function as matrix, where 
%       columns correspond to the scale and rows to the order q
%
%   Written by Bertrand Kerres, kerres@kth.se, last update 2016-10-31

    
    ip = inputParser ();
    ip.addRequired ('data', @(x) isvector(x) && isnumeric(x));
    ip.addRequired ('scale', @(x) all(x>0) && size(x,1) == 1);
    ip.addOptional ('q', 2, @(x) isnumeric(x));
    ip.addOptional ('m', 1, @(x) isscalar(x) && x >= 0);
    ip.addOptional ('integrate', 1, @(x) isscalar(x) && (x==-1 || x==0 || x==1 ||x==2));
    ip.addParameter ('R2warn', 0, @(x) isscalar(x) && x>= 0 && x <= 1);
    ip.parse (data, scale, varargin{:});
    
    % Make q into row vector
    qv = ip.Results.q(:)';

    data = data - mean (data);
    if isrow (data)
        data = data';
    end
    
    m = round(ip.Results.m);    % Trend to be removed
    scale = ip.Results.scale;   % Scale vector

    % Initialize return matrix
    H = NaN (length(qv), 1);
    
    % Calculate profile function by integrating (default == 1) or
    % differentiating
    if ip.Results.integrate == 2
        data = cumsum (data);
        data = cumsum (data-mean(data));
    elseif ip.Results.integrate == 1
        data = cumsum (data);
    elseif ip.Results.integrate == -1
        data = [0; diff(data)];
        data = data - mean(data);
    end
        

    N = length (data);
    F = zeros (length(qv), length(scale));

    sc = 1;
    ws = warning('off','all');    % For the polyfit
    
    for s = scale
        Ns = floor (N/s);   % Number of segments
        F2_DFA = zeros (2*Ns, 1);
        
        % Calculate Vandermonde matrix for this scale, and its QR
        % decomposition (which is constant for one scale)
        
        % Using a scale like this avoids badly conditioned Vandermonde
        % matrices; must be col vector!
        s_idx = (linspace(-s/2, s/2, s) / (s/2))';
        
        % Vdm = [s_idx^m, s_idx^(m-1), s_idx^(m-2), ..., s_idx^0]
        Vdm = zeros (length(s_idx), m+1);
        Vdm(:,m+1) = ones (length(s_idx),1);
        for j = m:-1:1
           Vdm(:,j) = s_idx .* Vdm(:,j+1);
        end
        % QR decomposition
        [Q,R] = qr (Vdm, 0);        

        % Loop through data and calculate squared error after detrending
        for v = 1 : Ns
            % idx1, idx2, yk_idx1, yk_idx2 have length s
            idx1 = (v-1)*s+1 : 1 : v*s;         % Beginning to end
            idx2 = N+1-v*s : 1 : N-(v-1)*s;    % End to beginning
            Yk_idx1 = data(idx1);
            Yk_idx2 = data(idx2);

            % Detrendend fluctuation analysis
            % This is in effect the same as:
            %   p1 = polyfit (s_idx, Yk_idx1, m)
            %   SS1 = sum ((Yk_idx1 - polyval(p1, s_idx)).^2);
            % but much faster since:
            %   no need to check for Vdm condition estimate
            %   QR decomposition done only once per scale
            p1 = R \ (Q'*Yk_idx1);                  % Polyfit coefficients
            SS1 = sum ((Yk_idx1 - Vdm*p1).^2);      % Sum of squared error
            
            p2 = R \ (Q'*Yk_idx2);                  % Polyfit coefficients
            SS2 = sum ((Yk_idx2 - Vdm*p2).^2);      % Sum of squared error

            F2_DFA(v) = SS1 / s;
            F2_DFA(Ns+v) = SS2 / s;
        end
        
        
        % Calculate qth root mean of the window variances
        qc = 1;
        for q = qv
            if q ~= 0
                F(qc,sc) = (mean( F2_DFA.^(q/2))).^(1/q);
            else
                F(qc,sc) = exp(1/2 * mean( log(F2_DFA) ) );
            end
            qc = qc + 1;
        end
        sc = sc + 1;
    end
    
    warning(ws);
    
    % Fit linear in the log-log plot of time scale vs fluctuation function
    % --> Hurst exponent
    for qc = 1 : length(qv)
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