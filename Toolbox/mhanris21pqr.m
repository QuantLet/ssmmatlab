function str = mhanris21pqr(y, l, res, x, str)
% PURPOSE: performs the second step of the multivariate Hannan-Rissanen
% method for VARMAX models with restrictions and returns the estimated parameters
%---------------------------------------------------
% USAGE: str = mhanris2(y,res,x,str)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           l    = an integer to specify the l-th variable of y
%           res    = an (nobs x neqs) matrix of residuals estimated with a long VARX model
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%           str    = a structure containing the structure of the VARMAX
%---------------------------------------------------
% RETURNS: str = a structure containing the previous structure plus
%                the estimated parameters
%---------------------------------------------------
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sgpg.meh.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

if nargin ~= 5
    error('wrong # of arguments to mhanris21pqr');
end;

[nobs, neqs] = size(y);
if ~isempty(x)
    [nobs2, nx] = size(x);
    if (nobs2 ~= nobs)
        error('mhanris21pqr: nobs in x-matrix not the same as y-matrix');
    end
else
    nx = 0;
end
[nobs3, nr] = size(res);
if (nobs3 ~= nobs)
    error('mhanris21pqr: nobs in res-matrix not the same as y-matrix');
end;


s = str.s;
m = str.m;
kro = str.kro;
p = str.p;
q = str.q;
r = str.r;

if (s ~= neqs)
    error('mhanris21pqr: s nonequal dim(y_t)');
end
if (m ~= nx)
    error('mhanris21pqr: m nonequal dim(x_t)');
end


yy = y(:, l);

% adjust nobs to feed the lags
nlag = max(kro);
nobse = nobs - nlag;

ylag = glags(y, p);
if (p < nlag)
    ylag = ylag(nlag-p+1:end, :);
end
reslag = glags(res, q);
if (q < nlag)
    reslag = reslag(nlag-q+1:end, :);
end
if (nx > 0)
    xlag = ltflag(x, r);
    if (r < nlag)
        xlag = xlag(nlag-r+1:end, :);
    end
else
    xlag = [];
end
%we form a linear combination of the variables that has minimum (p, q, r)
%orders (scalar component). 
% vp1lm1 = y(nlag+1:end, 1:l-1);
% vplm1e = y(nlag+1:end, l+1:end); 
% vp = [vp1lm1, vplm1e];
vp = y(nlag+1:end, 1:l-1) - res(nlag+1:end, 1:l-1);



% form x-matrix, constant is included
X = [-vp, -ylag, reslag, xlag, ones(nobse, 1)];


yt = yy(nlag+1:end);
[beta, tv] = btval([], [yt, X]); % regression
str.beta = beta; % estimated parameters, included the constant at the end
str.tv = tv; % t-values

resid2 = yt - X * beta; % vectorized residuals
str.resid2 = resid2;
