function [Pi, alpha, betap, betaor, ferror] = evarma11rurimp(y, x, a, ir, nord)
% PURPOSE: performs VARMAX(1,1) estimation with rank imposed
%          and returns estimated Pi matrix in the form Pi=alpha*betap,
%          where Pi = -phi(1) in the error correction form.
%---------------------------------------------------
% USAGE: [beta,res] = evarma11r11Ru(yd,a,seas,y,x)
% where:    y     = an (nobs x neqs) matrix of observations
%           x     = matrix of input variables (nobs x nx)
%           a     = an (nobs x neqs) matrix of residuals
%           ir    = 0,1, corresponding to I(0) or I(1) case
%          nord   = rank of Pi
%                  (NOTE: constant vector automatically included)
%---------------------------------------------------
% RETURNS: Pi     = an (nobs x nobs) matrix, equal to -phi(1), and such
%                   that Pi = alpha*betap
%          alpha  = an (nobs x nord) matrix
%          betap  = an (nord x nobs) matrix
%          betaor = an (nobs x nord) matrix orthogonal to betap
%         ferror  = a flag for errors
%---------------------------------------------------
% Copyright (c) 21 July 2014 by Victor Gomez
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

ferror = 0;
Pi = [];
alpha = [];
betap = [];
if nargin ~= 5
    ferror = 1;
    disp('wrong # of arguments to evarma11rurimp');
    return
end;


[nobs, s] = size(y);
[nobsa, sa] = size(a);
if (nobsa ~= nobs)
    ferror = 3;
    disp('evarma11rurimp: nobs in a-matrix not the same as in y-matrix');
end;


[mx, m] = size(x);
if (mx ~= nobs) && (mx > 0)
    ferror = 3;
    disp('evarma11r11Rurimp: nobs in x-matrix not the same as in y-matrix');
end;

seas = 1;

%variables should be centered
med = mean(y)';
y = y - repmat(med', nobs, 1); %center variables
meda = mean(a)';
a = a - repmat(meda', nobs, 1); %center variables
if (mx > 0)
    medx = mean(x)';
    x = x - repmat(medx', nobs, 1); %center variables
end

% adjust nobs to feed the lags
if (seas <= 1), seas = 0;
end
nlag = seas + 1;
nobse = nobs - nlag;
ylag = glags(y, nlag);
alag = glags(a, nlag);
if mx > 0
    xlag = ltflag(x, nlag);
else
    xlag = [];
end

%obtain differenced variables for error correction form
twos = 2 * s;
threes = 3 * s;
mxf = m * (nlag + 1);
if (ir == 1)
    yd = diferm(y, 1);
    yd = yd(nlag:end, :);
    A = [xlag, alag, ylag(:, 1:s), yd];
    ind2 = mxf + s + 1:mxf + twos;
    ind3 = mxf + twos + 1:mxf + threes;
end

%reduced rank regression
%A=[xf xp yp yf];
[Q, R] = qr(A./sqrt(double(nobse)), 0);
L = R'; % [L_{11} 0 0; L_{21} L_{22} 0; L_{31} L_{32} L_{33}]
% ind2=mx*f+1:mx*(f+p)+my*p;
% ind3=mx*(f+p)+my*p+1:mx*(f+p)+my*(f+p);
L32 = L(ind3, ind2);
L22 = L(ind2, ind2);
L33 = L(ind3, ind3);
Sffbu = L32 * L32' + L33 * L33'; % \Sigma_{ff|u}
% Sppbu=L22*L22';          % \Sigma_{pp|u}
Sfpbu = L32 * L22'; % \Sigma_{fp|u}
%we use svd instead of Cholesky to obtain \Sigma_{ff|u}^{1/2} because
%\Sigma_{ff|u} can be of less than full rank. If \Sigma_{ff|u} is of full
%rank, the orthogonal matrices U and V coincide.
[Uf, SVf, Vf] = svd(Sffbu);
SSVf = diag(sqrt(diag(SVf)));
Rffbu = Uf * SSVf * Vf'; % \Sigma_{ff|u}^{1/2}
Lffbu = Vf * pinv(SSVf) * Uf'; % \Sigma_{ff|u}^{-1/2}
%again, we use svd instead of Cholesky to obtain \Sigma_{pp|u}^{1/2}
%because \Sigma_{pp|u} can be of less than full rank. If \Sigma_{ff|u} is
%of full rank, the orthogonal matrices U and V coincide.
[Up, SSVp, Vp] = svd(L22);
SSVpi = pinv(SSVp);
Lppbu = Vp * SSVpi * Up'; % \Sigma_{pp|u}^{-1/2}
AA = Lffbu * Sfpbu * Lppbu'; % \Sigma_{ff|u}^{-1/2}\Sigma_{fp|u}\Sigma_{pp|u}^{-1/2 '}
[U, SV, V] = svd(AA);


alpha = Rffbu' * U(:, 1:nord) * diag(sqrt(diag(SV(1:nord, 1:nord))));
betap = diag(sqrt(diag(SV(1:nord, 1:nord)))) * V(:, 1:nord)' * Lppbu;
% Pi=Rffbu'*U(:,1:nord)*SV(1:nord,1:nord)*V(:,1:nord)'*Lppbu;  %reduced rank estimator
Pi = alpha * betap;
betaor = Lppbu \ V(:, nord+1:end);
