function [lagsopt, a, ferror] = lratiopppt(y, x, seas, maxlag, minlag, prt)
% PURPOSE: performs sequential likelihood ratio tests in varmax(p,p,p)
%          models to determine optimal p, starting from p=minlag and
%          proceeding up to p=maxlag. All models are estimated
%          using the same sample size: nobs - maxlag.
%          The innovations are computed as in the Hannan-Rissanen method.
%          That is, using a VARX approximation with length given by a
%          formula that depends on the sample size. 
%---------------------------------------------------
% USAGE:  [lagsopt,ferror] = lratiopppt(y,x,seas,maxlag,minlag,prt)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%           x    = matrix of input variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%         seas   = seasonality
%         maxlag = the maximum lag length. If empty on entry, it is
%                  calculated by the program as the order of the VARX
%                  approximation. In this case, or if maxlag is > 0, this
%                  maximum lag length it used to identify a VARMA(p,p,p) 
%                  model.
%         minlag = the minimum lag length
%            prt = flag for printing
%                    0 = no, 1 = yes
%                    (default = 0)
%---------------------------------------------------
% RETURNS: lagsopt, the optimum number of lags
%                      a = residuals of the VARX approximation (nobs x
%                      neqs)
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

lagsopt = [];
ferror = 0;

if nargin ~= 6
    disp('wrong # of arguments to lratiopqr');
    ferror = 1;
    return
end;

if maxlag < minlag
    disp('maxlag < minlag in lratiopqr');
    ferror = 2;
    return
end;

[nobs, s] = size(y);
if ~isempty(x)
    [mx, nx] = size(x);
    if (mx ~= nobs)
        disp('lratiopqr: nobs in x-matrix not the same as y-matrix');
        ferror = 3;
        return
    end
else
    nx = 0;
end

%First, determination of optimum VARX length
if (nx > 0)
    minp = 5;
else
    minp = 3;
end
if (seas > 1)
    aa = 2.;
else
    aa = 1.25;
end
pt = max(8, seas+minp);
minlags = 0;
maxlags = ceil(max(log(nobs)^aa, pt)); %VARX length is given by this formula

%%%%%%%%%%%%%%%%%%%%%%%
%compute VARX residuals using a VARX length given by a formula that depends 
%on the sample size

% compute initial residuals using VARXs with different lengths
initres = zeros(maxlags, s);
for i = minlags:maxlags - 1
    if nx == 0
        resid = var_res(y, i);
        initres(i+1, :) = resid(1, :);
    else
        resid = varx_res(y, i, x);
        initres(i+1, :) = resid(1, :);
    end
end
% compute the rest of residuals using a VARX(maxlags)
if nx == 0
    [res] = var_res(y, maxlags);
else
    [res] = varx_res(y, maxlags, x);
end

a = [initres(1:maxlags, :); res]; % use all possible VARX residuals

%If maxlag is empty, we compute it as the VARX length chose by the bic
%criterium
if isempty(maxlag) 
    if nx == 0
        % lagsopt = lratiocr(y, maxlags, minlags, prt);
        crt = 'bic';
        lagsopt = infcr(y, maxlags, minlags, crt, prt);
    else
        % lagsopt = lratiocrx(y, maxlags, minlags, prt, x);
        crt = 'bic';
        lagsopt = infcr(y, maxlags, minlags, crt, prt, x);
    end
    if prt == 1
        fprintf(1, 'estimated order in VARX = %2d\n\n', lagsopt);
    end  
    maxlag = lagsopt;
end
%
%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %alternative to fixed VARX order to compute VARX residuals
% if nx == 0
%     % [lagsopt, initres] = lratiocr(y, maxlags, minlags, prt);
%     crt = 'bic';
%     [lagsopt, initres] = infcr(y, maxlags, minlags, crt, prt);
%     res = var_res(y, lagsopt);
% else
%     % [lagsopt, initres] = lratiocrx(y, maxlags, minlags, prt, x);
%     crt = 'bic';
%     [lagsopt, initres] = infcr(y, maxlags, minlags, crt, prt, x);
%     res = varx_res(y, lagsopt, x);
% end
% a = [initres(1:lagsopt, :); res]; %residuals of the VARX approximation
% if prt == 1
%    fprintf(1, 'estimated order in VARX = %2d\n\n', lagsopt);
% end
% 
% if isempty(maxlag) || (maxlag > lagsopt)
%     maxlag = lagsopt;
% end
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%identify a varmax(p,p,p)
incr = 1; 
lagsopt = lratiocrax(y, maxlag, minlag, incr, prt, a, x);
if prt == 1
    fprintf(1, 'estimated order in VARMAX(p,p,p) = %2d\n\n', lagsopt);
end
