function [xvf, str, ferror] = mconestim(y, x, str)
%
% This function estimates a VARMAX model in echelon form with some of its
% parameters possibly restricted to zero. Estimation is performed
% using the conditional method. It is assumed that initial values obtained
% with the three stages of the Hannan-Rissanen method are available
% in str.beta3. These initial estimates can be obtained using the functions
% MHANRIS, ESTVARMAXPQRPQR or ESTVARAMXKRO.
% The state space echelon form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = H*alpha_{t} + D*x_{t}  + a_{t}
%
% The VARMAX echelon form is

%  phi(B)*y_{t} = gamma(B)*x_{t} + theta(B)a_{t},
%
% where B is the backshift operator, B*y_{t} = y_{t-1}.
%
%
% Inputs: y      = an (nobs x s) matrix of y-vectors
%         x      = matrix of input variables (nobs x m)
%                   (NOTE: constant vector automatically included)
%       str: a structure containing a preliminary model estimation obtained
%            with functions mhanris, estvaramxpqrPQR or estvarmaxkro.
%  Output: xvf: estimated parameters
%          str: a structure containing the estimated VARMAX model.
%       ferror: flag for errors
% The fields of structure str on output are the same as those on input plus
% the following
%       residcon: an [(nobs-nlag) x s] matrix of standardized residuals
%                 after estimation, where nlag = max(kro)
%      sigmarcon: estimated covariance matrix of innovations
%        phiscon: an (s x s x nlag) array containing the estimated AR
%                 parameters
%      thetascon: an (s x s x nlag) array containing the estimated MA
%                 parameters
%      gammascon: an (s x m x nlag) array containing the estimated
%                 exogenous parameters
%       phitvcon: an (s x s x nlag) array containing the p-values of the
%                 estimated AR parameters
%     thetatvcon: an (s x s x nlag) array containing the p-values of the
%                 estimated MA parameters
%     gammatvcon: an (s x m x nlag) array containing the p-values of the
%                 estimated exogenous parameters
%         muscon: an (s x 1) vector containing the estimated constant
%        mutvcon: an (s x 1) vector containing the p-values of the constant
%       phistcon: same as phiscon but with coefficient matrices
%                 premultiplied by phis(:,:,1)^{-1} (VARMAX model not in
%                 echelon form)
%     thetastcon: same as thetascon but with coefficient matrices
%                 premultiplied by phis(:,:,1)^{-1} (VARMAX model not in
%                 echelon form)
%     gammastcon: same as gammascon but with coefficient matrices
%                 premultiplied by phis(:,:,1)^{-1} (VARMAX model not in
%                 echelon form)
%          Fscon: an (n x n) matrix containing the estimated F, where n is
%                 the McMillan degree = sum of the Kronecker indices
%          Kscon: an (n x s) matrix containing the estimated K
%          Bscon: an (n x m) matrix containing the estimated B
%          Dscon: an (s x m) matrix containing the estimated D
%          Hscon: an (s x n) matrix containing the estimated H
%      ferror   = flag for errors
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

ferror = 0;


%check stationarity and invertibility.
iar = chkstainv(str.Fs3); %if iar >1, the model is not stationary
if iar > 1
    fprintf(1, 'initial model nonstationary in mconestim, iar = %2d\n', iar);
    str = sta3r(str);
end
ima = chkstainv(str.Fs3-str.Ks3*str.Hs3); %if ima >1, the model is not invertible
if ima > 1
    fprintf(1, 'inital model noninvertible in mconestim, ima = %2d\n', ima);
    str = inv3r(str);
end


%initial parameters
xv = str.beta3'; 
% xv0=xv;                             %initial value for penalization
% xv0=.85*xv0;


f = 'conmedfv';
tr = 0;
mvx = 1; %parameters for Levenberg-Marquardt method:
tolf = 1e-4; %f : a function to evaluate the vector ff of
%    individual functions such that ff'*ff is
%    minimized
maxit = 100;
nu0 = .01;
jac = 1;
prt = 2; %tr: >0 x is passed from marqdt to f but not
%       passed from f to marqdt
clear infm
infm = minfm(f, tr, mvx, tolf, maxit, nu0, jac, prt, [], []);
infm.F = 'conmedfjac';
%Levenberg-Marquardt method
[xvf, J, ff, g, iter, conf] = marqdt(infm, xv, y, x, str);
%this part added 21-2-2011
if (abs(sum(sum(J))) <= 1.d-8)
    error('estimation failed in varmapqPQestim')
end
%end of addition
beta = xvf';


%standard errors and t-values
V = -J; % negative derivatives for the regression
res2t = ff + V * beta; % regression
[beta3, tv3] = btval([], [res2t, V]);
str.beta3 = beta;
str.tv3 = tv3;
% insert estimated parameters in vgam
vgam = str.vgam;
bind = str.bind;
[nbind, mbind] = size(bind);
vgams3 = vgam;
vgamtv3 = vgam;
for i = 1:nbind
    vgams3(bind(i)) = beta(i);
    vgamtv3(bind(i)) = tv3(i);
end
str.vgams3 = vgams3;
str.vgamtv3 = vgamtv3;


%use an auxiliary structure to compute the state space and VARMAX
%echelon forms corresponding to the conditional estimation
str3 = str;
str3.vgams = str.vgams3;
str3.vgamtv = str.vgamtv3;
str3 = param2sse(str3);
%store results
str.phiscon = str3.phis;
str.thetascon = str3.thetas;
str.gammascon = str3.gammas;
str.phitvcon = str3.phitv;
str.thetatvcon = str3.thetatv;
str.gammatvcon = str3.gammatv;
str.muscon = str3.mu;
str.mutvcon = str3.mutv;
str.phistcon = str3.phist;
str.thetastcon = str3.thetast;
str.gammastcon = str3.gammast;
str.Fscon = str3.Fs;
str.Kscon = str3.Ks;
str.Bscon = str3.Bs;
str.Dscon = str3.Ds;
str.Hscon = str3.Hs;


%estimated conditional residuals and their covariance matrix
nlag = max(str.kro);

resid = compresde0(y, x, str);
residcon = resid(nlag+1:end, :);
sigmarcon = cov(residcon, 1);
str.residcon = residcon;
str.sigmarcon = sigmarcon;
