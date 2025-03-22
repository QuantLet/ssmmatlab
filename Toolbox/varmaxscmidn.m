function [order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt)
% PURPOSE: identifies a VARMAX model based on scalar component models for
% the series y with inputs x. The series is assumed to follow an invertible
% but possibly nonstationary model.
% The scalar components are identified after a VARMAX(p,q,r) has been
% previously identified and estimated. The estimated  max(p,q,r) is an
% estimate of the maximum Kronecker index.
% If maxorder is empty, the maximum order for the VARMAX(p,q,r) model is
% set equal to the order of a VARX approximation.
% If maxorder is positive, a VARMAX(p,q,r) model is estimated using
% maxorder as maximum for p, q and r.
% The procedure identifies the s.c._i, i, i=1,2,...,s, equation by
% equation. En each equation, first, the past innovations are replaced with
% the estimated innovations of the VARMAX(p,q,r) model of the first step.
% Then, a sequence of LR tests allows for the estimation of the s.c. of
% that equation.
%---------------------------------------------------
% USAGE: [order,kro,scm] = varmaxscmidn(y,x,seas,maxorder,hr3,prt)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%        seas      = seasonality
%    maxorder      = empty, the programm will compute the maximum order for 
%                    the VARMAX(p,q,r) model to be identified.
%                    >0, this maximum order is used as the initial maximum
%                    order for the VARMAX(p,q,r) model to be identified.
%                    <0, the maximum order is fixed to -maxlag
%         hr3      = 1 perform only the first two stages of the HR method
%                    0 perform the three stages of the HR method
%         prt      = 1 print results of the VARX and VARMAX(p,p,p) tests
%---------------------------------------------------
% RETURNS: order = the system order, that is, the McMillan degree
%            kro = the Kronecker indices
%            scm = an array containg the scalar component models
%---------------------------------------------------
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

[ny, s] = size(y);
if ~isempty(x)
    [nx, m] = size(x);
    if (nx ~= ny)
        error('varmaxscmidn: nobs in x-matrix not the same as y-matrix');
    end
end
if nargin < 6
    prt = 0;
end
if nargin < 5
    hr3 = 1;
    prt = 0;
end

if nargin < 4
    maxorder = [];
    hr3 = 1;
    prt = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify a VARMAX(p,q,r) model
minorder = 0;
[lagsopt, ferror] = lratiopqr(y, x, seas, maxorder, minorder, prt);
p = max(lagsopt); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimate the VARMAX(p,q,r) model
kro = repmat(p, 1, s);
if hr3 == 0
    str = estvarmaxpqrPQR(y, x, seas, lagsopt, [0, 0, 0], hr3, 1);
else
    str = estvarmaxpqrPQR(y, x, seas, lagsopt, [0, 0, 0], hr3);
end
residv = str.residv; %first residuals obtained with VARX
residv = zeros(size(residv));   %first residuals equal zero
if hr3 == 0
    resid = str.resid3;
else
    resid = str.resid2;
end
nr = size(resid, 1);
if ny > nr
    residv = [residv(1:ny-nr, :); resid];
else
    residv = resid;
end
clear str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%loop over all variables to determine the s.c.
if (seas > 1)
    minlag = seas;
else
    minlag = 0;
end
maxlag = max(p, seas);
scm = zeros(s, 3);
for i = 1:s
    if prt == 1
        fprintf(1, 'variable = %2d\n', i);
    end
    [lagsopt] = lratiopqr1(y, i, residv, x, maxlag, minlag, prt);
    kro(i) = max(lagsopt);
    scm(i, :) = lagsopt;
end
order = sum(kro);
