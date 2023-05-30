function [result, str] = usmestim(y, str)
%
% This function estimates a univariate structural model using the exact
% maximum likelihood method.
%
%
% Inputs:
%         y: matrix containing the output series
%       str: a structure containing the initial model information. It
%            should be input as well as output because the concentrated
%            parameter can change. The concentrated parameter should be the
%            greatest variance. The program performs a preliminary
%            estimation to check it.
%  Outputs:
%          result: a structure with the following fields
%           .xvf : estimated parameters
%           .xf : vector of fixed parameters
%       .sigma2c: the estimated standard error of the parameter that has
%                 been concentrated out
%          .Pevf: prediction error variance
%         .SPevf: square root of Pevf
%          .tv  : t-values of the estimated parameters
%           .e  : vector of standardized residuals at the end of estimation
%                 (Q'_2*y)
%           .Ss : residual sum of squares (e'*e)
%           .F  : vector of nonlinear functions whose sum of squares is
%                 minimized at the end of estimation
%           .Ff : the product F'*F
%           .h  : vector of estimated regression estimates
%           .M  : matrix of mse of h
%           .A  : estimated state vector, x_{t|t-1}, obtained with the
%                 Kalman filter at the end of the sample
%           .P  : Mse of A
%          .tvr : t-values of the estimated regression parameters
%          .ser : standard errors of the estimated regression parameters
%       .ferror : flag for errors
%           str : the same structure as input. It should be present because
%                 the concentrated parameter may change
%
% Copyright (c) January 2010 by Victor Gomez
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

%the following three lines added 13/04/2012
if (nargout < 2)
    error('Structure str should also be output in usmestim')
end
%end of addition 13/04/2012
ferror = 0;


%initial parameters
xv = str.xv;
xf = str.xf;
pfix = str.pfix;
pvar = str.pvar;
X = str.X;
% this line added 12-2011
W = str.W;
% x0=str.x;
nr = length(pvar);
s = str.freq;


%Levenberg-Marquardt method
f = 'smfest';
if ((str.arp > 0) || (str.cycle > 0)), tr = 0;
else tr = 1;
end
mvx = 1;
chb = 0;
prt = 2;
if isfield(str,'tolf')
    tolf = str.tolf;
else
    tolf = 1e-4;
end
nu0 = .01;
jac = 1;

%preliminary estimation to know what is the biggest variance.
% disp('preliminary estimation')
maxit0 = 15;
prt0 = 1;
clear infm
infm = minfm(f, tr, mvx, tolf, maxit0, nu0, jac, prt0, [], []);
[xvp, J, ff, g, iter] = marqdt(infm, xv, y, s, pfix, pvar, xf, chb, str); %xvp,pause
maxabsxvidx = 0;
maxabsxv = -1.d10;
stord = str.stord;
for i = 1:nr
    if (stord(pvar(i)) <= 6)
        absxv = abs(xvp(i));
        if (absxv > 1.) && (absxv > maxabsxv)
            maxabsxv = absxv;
            maxabsxvidx = i;
        end
    end
end
if (maxabsxvidx > 0) %change of concentration parameter
    oldconc = str.conc;
    str.conc = stord(pvar(maxabsxvidx));
    stord(pvar(maxabsxvidx)) = oldconc;
    oldstordc = stord;
    nstord = sort(stord);
    str.stord = nstord;
    % Update of the initial parameter vectors x and xv
    cc = {'str.comp.irreg', 'str.comp.level', 'str.comp.slope', ...
        'str.comp.seas', 'str.comp.ar', 'str.comp.cycle'};
    oldc = eval(cc{oldconc});
    %    oldc,pause
    oldx = oldc(2);
    %   oldx, pause
    xv(maxabsxvidx) = oldx; %xv,pause
    oldpvar = pvar;
    oldpfix = pfix;
    npvar = zeros(size(pvar));
    npfix = zeros(size(pfix));
    nf = length(pfix);
    for i = 1:length(stord)
        for j = 1:nr
            if nstord(i) == oldstordc(oldpvar(j))
                npvar(j) = i;
            end
        end
        for j = 1:nf
            if nstord(i) == oldstordc(oldpfix(j))
                npfix(j) = i;
            end
        end
    end
    [npvar, ipvar] = sort(npvar);
    [npfix, ipfix] = sort(npfix);
    % a new variance has been concentrated out, but the initial conditions
    % are the same.
    xv = xv(ipvar);
    xf = xf(ipfix);
    pvar = npvar;
    pfix = npfix;
    str.pvar = npvar;
    str.pfix = npfix;
    str.xv = xv;
    str.xf = xf;
    str.x(npvar) = xv;
    str.x(npfix) = xf;
    [X, Z, G, W, T, H, ins, ii, ferror] = pr2usm(xv, xf, str);
    str.X = X;
    str.Z = Z;
    str.G = G;
    str.W = W;
    str.T = T;
    str.H = H;
    str.i = ii;
    %   str,pause
end
% disp('end of preliminary estimation')
%end of preliminary estimation

if (maxabsxvidx > 0) || ((maxabsxvidx == 0) && (iter == maxit0))
    maxit = 100;
    clear infm
    infm = minfm(f, tr, mvx, tolf, maxit, nu0, jac, prt, [], []);
    [xvf, J, ff, g, iter, conf] = marqdt(infm, xv, y, s, pfix, pvar, xf, chb, str);
    %this part added 21-2-2011
    if (abs(sum(sum(J))) <= 1.d-8)
        error('estimation failed in usmestim')
    end
    %end of addition
else
    xvf = xvp;
end

if ~isempty(xv)
    %t-values
    V = -J; % negative derivatives for the regression
    res2t = ff + V * xvf'; % regression
    [~, tv, se, Cv] = btval([], [res2t, V]); % tv contains the t-values
else % all parameters are fixed
    tv = [];
    se = [];
    Cv = [];
end

% this line changed 12-2011
if ~isempty(X) || ~isempty(W)
    chb = 1;
end
%
% get residuals and estimator
%
[F, e, h, M, Pevf, A, P, olsres] = smfun(xvf, y, s, pfix, pvar, xf, chb, str);

Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
% disp('concentrated parameter:')
sigma2c = Ss / ne; %estimated sigma square; %line changed on Nov 20, 2019
% compute prediction error variance (finite sample)
% disp('prediction error variance (finite sample)')
Pevf = Pevf * sigma2c;

result.xvf = xvf;
result.xf = xf;
result.e = e;
result.Ss = Ss;
result.Ff = Ff;
result.sigma2c = sigma2c;
result.Pevf = Pevf;
SPevf = sqrt(diag(Pevf));
result.SPevf = SPevf;
result.tv = tv';
result.se = se';
result.Cv = Cv;
result.F = F;
result.h = h;
result.M = M;
result.A = A;
result.P = P;
result.olsres = olsres;
% this line changed 12-2011
if ~isempty(X) || ~isempty(W)
    %t-values of regression parameters
    M = M * sigma2c;
    ser = sqrt(diag(M));
    tvr = h ./ ser;
    result.tvr = tvr;
    result.ser = ser;
else
    result.tvr = [];
    result.ser = [];
end
result.ferror = ferror;
