function [ sces, sxcc, xs ] = findshift2( x, s, shiftmax )
%SXCORR  Computes the scale cross-correlation between vectors X and S
%   [SCES SXCC] = FINDSTRETCH(X,S) computes the scale or stretch estimate 
%   between two signals as well as the the scale-invariant correlation 
%   coefficient. The estimates occur in two stages, first estimating the 
%   scale factor through the scale cross-corrleation function. The estimate
%   is then refined in the scale domain, assuming convexity. 
%
%   INPUTS:
%       X: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%       S: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%
%   OUTPUTS:
%    SCES: A K-by-1 vector of scale estimate(s) between X and S
%    SXCC: A K-by-1 vector of scale-invariant correlation coefficients 
%          between X and S
%
%   see also: sxcorr, fmt, fmincon
%

% -------------------------------------------------------------------------
% Copyright (C) 2012,2014  Joel B. Harley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at 
% your option) any later version. You should have received a copy of the 
% GNU General Public License along with this program. If not, see 
% <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------- 
% IF THIS CODE IS USED FOR A RESEARCH PUBLICATION, please cite:
%   J.B. Harley, J.M.F. Moura, "Scale transform signal processing for 
%   optimal ultrasonic temperature compensation," IEEE Transactions on  
%   Ultrasonics, Ferroelectrics and Frequency Control, vol. 59, no. 10, 
%   October 2012.
% -------------------------------------------------------------------------
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%


% MEASURE SAMPLE SIZE
N = size(s,1);          % Number of samples
M = size(s,2);          % Number of samples

% ESTABLISH SCALE-FREQUENCY AXIS
r = floor(N/2)+1; cn = ifftshift((((1:N)-r)/(N))).'; 
r = floor(M/2)+1; cm = ifftshift((((1:M)-r)/(M))).'; 

% FIND SCALE TRANSFORMS AND NORMALIZE
X = fft2(x); S = fft2(s);

% MULTIPLY IN THE SCALE DOMAIN
SX = conj(X).*S;

% FIND THE SCALE CROSS-CORRELATION
sx = M*ifft2(SX); 
[sxm, sxi] = max(sx);
[~  , sxx] = max(sxm);
sxy = sxi(sxx);


sxx = sxx-1;
sxy = sxy-1;
sxx(sxx>M/2) = sxx(sxx>M/2)-M;
sxy(sxx>N/2) = sxy(sxy>N/2)-N;

if abs(sxi) > shiftmax
    sces = 0;
else

    % SET OPTIMIZATION OPTIONS 
    options = optimset('DiffMaxChange',.0001, ...
                       'Algorithm', 'active-set', 'TolX', 1e-12, ...
                       'TolFun', 1e-12, ...
                       'Display', 'notify', 'MaxFunEvals', 10000); 

    % IMPROVE ESTIMATE (ASSUMING CONVEXITY AROUND THE MAXIMUM)
    sces = fmincon( ...
        @(sces) optfun(X, S, cn, cm, sces), ... 
        ([sxy sxx]), [], [], [], [], ([sxy sxx])-1, ([sxy sxx])+1, [], options );

end

% FIND THE NEW SCALE-INVARIANT CORRELATION COEFFICIENT
%sxcc = real(S'*(X.*exp(-1j*2*pi*sces*c)))/norm(X)/norm(S);
%sxcc = real(S'*(X.*exp(-1j*2*pi*sces*c)))/norm(X)^2;

sxcc = real(sum(sum(conj(S).*(X.*(exp(-1j*2*pi*(sces(1)*cn))*exp(-1j*2*pi*(sces(2)*cm)).') ))))/norm(S,'fro')/norm(X,'fro');

xs = real(ifft2(X.*(exp(-1j*2*pi*(sces(1)*cn))*exp(-1j*2*pi*(sces(2)*cm)).')));

%ss = real(ifft(S));
%xx = real(ifft(X.*exp(-1j*2*pi*sces*c)))/sxcc2;

end


function [ cx ] = optfun( X, S, cn, cm, shift )
    cx = 1-real(sum(sum(conj(S).*(X.*(exp(-1j*2*pi*(shift(1)*cn))*exp(-1j*2*pi*(shift(2)*cm)).') ))))/norm(S,'fro')/norm(X,'fro');
end