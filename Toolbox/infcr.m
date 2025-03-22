function [lagsopt, initres] = infcr(y, maxlag, minlag, crt, prt, x)
% PURPOSE: to determine VARX optimal lag length using information criteria
%---------------------------------------------------
% USAGE:  lagsopt = infcr(y,maxlag,minlag,crt,prt,x)
% where:    y    = an (nobs x nvar) matrix of y-vectors
%           maxlag = the maximum lag length
%           minlag = the minimum lag length
%           crt = the information criterion, 'aic' or 'bic'
%                  (default = 'aic')
%           prt = flag for printing
%                    0 = no, 1 = yes
%                    (default = 0)
%           x    = matrix of input variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%---------------------------------------------------
% RETURNS: lagsopt, the optimum number of lags
%          initres = an (maxlag x neqs) matrix of initial residuals
%                    corresponding to the estimated VARXs of order 1,2,...,
%                    maxlag.
%---------------------------------------------------
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

if nargin > 6
    error('wrong # of arguments to infcr');
elseif nargin == 6
    xflag = 1;
elseif nargin == 5
    xflag = 0;
elseif nargin == 4
    prt = 0;
    xflag = 0;
elseif nargin == 3
    prt = 0;
    xflag = 0;
    crt = 'aic';
end
if ~strcmp(crt, 'aic') && ~strcmp(crt, 'bic')
    error('wrong crt to infcr');
end

if maxlag < minlag
    error('maxlag < minlag in infcr');
end

[nobs, nvar] = size(y);
if (xflag ==1)
    [mx, nx] = size(x);
    if (mx ~= nobs)
       error('nobs nonequal to mx in infcr')
    end
else
    nx = 0;
end


initresa = zeros(maxlag+1, nvar);
dn = double(nobs-maxlag); %we condition on the first maxlag observations


% loop over lag lengths and compute information criteria
lagsopt = 0;
critm = 1.d10;
for lags = minlag:maxlag
    if (xflag == 1)
       resid = varx_res(y(maxlag-lags+1:end, :), lags, x(maxlag-lags+1:end, :));
    else
       resid = var_res(y(maxlag-lags+1:end, :), lags);
    end
    % first get var-cov matrices for residuals
    sigmar = cov(resid, 1);
    cra = log(det(sigmar));
    dnparbn = double(lags*nvar^2 + lags*nvar*nx + nvar*nx) / dn;
    if strcmp(crt, 'aic')
        pent = 2;
    else
        pent = log(dn);
    end
    crit = cra + dnparbn * pent;
    if prt == 1
        fprintf(1, 'nlag = %2d%s%s = %16.4f\n', lags, ' ', crt, crit)
    end
    if crit < critm
        critm = crit;
        lagsopt = lags;
    end
end
if (xflag == 1)
   resid1 = varx_res(y, 0, x);
else
   resid1 = var_res(y, 0);
end
initresa(1, :) = resid1(1, :);
initres = initresa(1:maxlag, :);
if nargout == 2
    for i = 1:maxlag
        if (xflag == 1)
           resid = varx_res(y, i, x);
        else
           resid = var_res(y, i);
        end
        initres(i, :) = resid(1, :);
    end
end