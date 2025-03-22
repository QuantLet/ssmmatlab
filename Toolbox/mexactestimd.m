function [xvf, str, ferror] = mexactestimd(y, x, str, Y)
%
% This function estimates a VARMAX model with unit roots parameterized in
% terms of the model for the ``differenced'' series using the exact maximum 
% likelihood method with a fast square root filter. Some of its parameters 
% can possibly be restricted to zero.
% 
% IIt is assumed that initial values
% obtained with the three stages of the Hannan-Rissanen,the conditional or
% the exact method are available in str.beta3. These initial estimates can 
% be obtained using the functions MHANRIS, ESTVARMAXPQRPQR or ESTVARAMXKRO
% for the Hannan-Rissanen method or function MCONESTIM for the conditional
% method or function MEXACTESTIM for the maximum likelihood method. 
% The state space echelon form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = Y_{t}*beta + H*alpha_{t} + D*x_{t}  + a_{t}
%
% The VARMAX echelon form is

%  phi(B)*(y_{t} - Y_{t}*beta) = gamma(B)*x_{t} + theta(B)a_{t},
%
% where B is the backshift operator, B*y_{t} = y_{t-1}.
%
% Inputs: y      = an (nobs x s) matrix of y-vectors
%         x      = matrix of input variables (nobs x m)
%       str: a structure containing a preliminary model estimation obtained
%            with functions mhanris, estvaramxpqrPQR or estvarmaxkro.
%         Y: (nobs*s x nbeta) matrix for regression variables. If empty,
%            the variables are centered and the mean is not estimated.
%  Output: xvf: estimated parameters
%          str: a structure containing the estimated VARMAX model.
%       ferror: flag for errors
% The fields of structure str on output are the same as those on input plus
% the following
%         sigma2c: concentrated parameter estimate
%      sigmarexct: estimated covariance matrix of innovations
%                  regression paramters
%               e: stack of white noise residuals
%        phisexct: an (s x s x nlag) array containing the estimated AR
%                  parameters
%      thetasexct: an (s x s x nlag) array containing the estimated MA
%                  parameters
%      gammasexct: an (s x m x nlag) array containing the estimated
%                  exogenous parameters
%       phitvexct: an (s x s x nlag) array containing the p-values of the
%                  estimated AR parameters
%     thetatvexct: an (s x s x nlag) array containing the p-values of the
%                  estimated MA parameters
%     gammatvexct: an (s x m x nlag) array containing the p-values of the
%                  estimated exogenous parameters
%         musexct: an (s x 1) vector containing the estimated constant
%        mutvexct: an (s x 1) vector containing the p-values of the
%                  constant
%       phistexct: same as phiscon but with coefficient matrices
%                  premultiplied by phis(:,:,1)^{-1} (VARMAX model not in
%                  echelon form)
%     thetastexct: same as thetascon but with coefficient matrices
%                  premultiplied by phis(:,:,1)^{-1} (VARMAX model not in
%                  echelon form)
%     gammastexct: same as gammascon but with coefficient matrices
%                  premultiplied by phis(:,:,1)^{-1} (VARMAX model not in
%                  echelon form)
%          Fsexct: an (n x n) matrix containing the estimated F, where n is
%                  the McMillan degree = sum of the Kronecker indices
%          Ksexct: an (n x s) matrix containing the estimated K
%          Bsexct: an (n x m) matrix containing the estimated B
%          Dsexct: an (s x m) matrix containing the estimated D
%          Hsexct: an (s x n) matrix containing the estimated H
%      ferror    = flag for errors
%
% Copyright (c) January 20150 by Victor Gomez
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


nr = str.nr;
DA = str.DA;

[ny, my] = size(y); %each row is an observation, each column is a variable
[nx, mx] = size(x);

%check stationarity and invertibility.
iar = chkstainv(str.Fs3); %if iar >1, the model is not stationary
if iar > 1
    fprintf(1, 'initial model nonstationary in mexactestimd, iar = %2d\n', iar);
    str = sta3r(str);
    [str, ferror] = aurirvarmapqPQ(str, nr, DA);
end
ima = chkstainv(str.Fs3-str.Ks3*str.Hs3); %if ima >1, the model is not invertible
if ima > 1
    fprintf(1, 'inital model noninvertible in mexactestimd, ima = %2d\n', ima);
    str = inv3r(str);
    [str, ferror] = aurirvarmapqPQ(str, nr, DA);
end


%variables are centered
if isempty(Y)
    y = y - repmat(mean(y), ny, 1);
    if mx > 0
        x = x - repmat(mean(x), ny, 1);
    end
end
[nY, nbeta] = size(Y);


%initial parameters
xv = str.beta3'; %polynomial parameter
%initial unit root parameters
xvd = str.xvd;
xv = [xvd xv];

if ~isfield(str,'sigma2c')
    sigmar = str.sigmar3; %covariance matrix parameters
    [R, p] = chol(sigmar);
    L = R';
    if (p > 0)
        disp('covariance matrix of residuals singular in mexactestimd')
        ferror = 1;
        return
    end
    L = L ./ L(1, 1); %sigma(1,1) is concentrated out of the
    %likelihood. The last my parameters in beta3 (xv) are the estimated mean
    Lh = vech(L);
    xv = [xv(1:end-my), Lh(2:end)'];
end
tol = 1d-6; %tolerance for CKMS recursions
% xv0=xv;                             %initial value for penalization (beta0)
% lambda=.88;
% xv0=lambda*xv0;


f = 'exactmedfvd';
tr = 0;
mvx = 1; %parameters for Levenberg-Marquardt method:
tolf = 1e-4; %f  :   a function to evaluate the vector ff of individual
maxit = 100;
nu0 = .01; %       functions such that ff'*ff is minimized
jac = 1;
prt = 2; %tr :   >0 x is passed from marqdt to f but not passed from f
clear infm
infm = minfm(f, tr, mvx, tolf, maxit, nu0, jac, prt, [], []);
% infm.F='conmedfjac';
%Levenberg-Marquardt method
chb = 0;
[xvf, J, ff, g, iter, conf] = marqdt(infm, xv, y, x, str, tol, Y, chb);
%this part added 21-2-2011
if (abs(sum(sum(J))) <= 1.d-8)
    error('estimation failed in mexactestimd')
end
%end of addition
beta = xvf';
chb = 1; %estimate regression variables
[fff, xvft, e, f, str, hb, Mb] = exactmedfvd(xvf, y, x, str, tol, Y, chb);
sigma2c = e' * e / double(ny*my); %concentrated parameter estimate
str.sigma2c = sigma2c;
bind = str.bind;
[nbind, mbind] = size(bind); %nbind includes the mean parameters
nparma = nbind - my; %number of parameters in ar and ma parts
L = zeros(my);
l = 0;
betam = [1., xvf(nparma+1:end)];
for i = 1:my
    cont = my - i + 1;
    ind = l + 1:l + cont;
    L(i:end, i) = betam(ind);
    l = l + cont;
end
Lh = vech(L);
sigmar = L * L' * sigma2c; %estimated exact covariance matrix of residuals
str.sigmarexct = sigmar;
str.e = e;
% sigma2c,sigmar

%standard errors and t-values
V = -J; % negative derivatives for the regression
res2t = ff + V * beta; % regression
[beta3, tv3] = btval([], [res2t, V]);
str.beta3 = beta;
str.tv3 = tv3;
% insert estimated parameters in vgam
vgam = str.vgam;
bind = str.bind;
vgams3 = vgam;
vgamtv3 = vgam;
nxvd = length(xvd);
nxv = length(xv);
betad = beta(1:nxvd);
beta = beta(nxvd + 1:nxv);
for i = 1:nparma
    vgams3(bind(i)) = beta(i);
    vgamtv3(bind(i)) = tv3(i);
end
str.vgams3 = vgams3;
str.vgamtv3 = vgamtv3;
lLh = length(Lh);
lbeta = length(beta);
if lLh > 1
    str.Lh = beta(lbeta-length(Lh)+2:lbeta);
    str.Lhtv3 = tv3(end-length(Lh)+2:end);
end

%use an auxiliary structure to compute the state space and VARMAX
%echelon forms corresponding to the exact m.l. estimation
str3 = str;
str3.vgams = str.vgams3;
str3.vgamtv = str.vgamtv3;
str3 = param2sse(str3);
%store results
str.phisexct = str3.phis;
str.thetasexct = str3.thetas;
str.gammasexct = str3.gammas;
str.phitvexct = str3.phitv;
str.thetatvexct = str3.thetatv;
str.gammatvexct = str3.gammatv;
str.musexct = hb;
str.mutvexct = hb ./ sqrt(diag(Mb));
str.phistexct = str3.phist;
str.thetastexct = str3.thetast;
str.gammastexct = str3.gammast;
str.Fsexct = str3.Fs;
str.Ksexct = str3.Ks;
str.Bsexct = str3.Bs;
str.Dsexct = str3.Ds;
str.Hsexct = str3.Hs;
str.xvd = betad;
str.xvdtv = tv3(1:nxvd);