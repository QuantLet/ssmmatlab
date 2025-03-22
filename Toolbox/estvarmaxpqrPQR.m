function [str, ferror] = estvarmaxpqrPQR(y, x, seas, ordersr, orderss, hr3, ...
    finv2, mstainv, nsig, tsig)
% PURPOSE: estimates a seasonal VARMAX model using the Hannan-Rissanen
%          method (function MHANRIS). It returns a structure containing the
%          estimated model. The estimated model is forced to be stationary
%          and invertible.
%          After having estimated a VARMAX model with ESTVARMAXPQRPQR, the
%          user can impose some zero restrictions in the model and
%          reestimate using the MHANRIS function.
%---------------------------------------------------
% USAGE: [str,ferror] = estvarmaxpqrPQR(y,x,seas,ordersr,orderss,hr3,finv2,mstainv,
%                              nsig,tsig)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%        seas      = seasonality
%        ordersr   = a (1 x 3) array containing the regular VARMAX orders
%        orderss   = a (1 x 3) array containing the seasonal VARMAX orders
%          hr3     = 1, perform only the first two stages of the HR method
%                    0, perform the three stages of the HR method, but only
%                    if the second stage model is invertible.
%        finv2     = 1, make model invertible after second stage of HR
%                    0, leave model as it is after second stage of HR
%      mstainv     = 1, use the DARE for enforcing stationarity or
%                    invertibility. This can only be used when there are no
%                    restrictions in the model.
%                  = 0, use multiplication by a small number.
%        nsig      = a (1 x 2) array. If nsig(i)=1, eliminate
%                    nonsignificant parameters after the i-th stage of the
%                    HR method, i=1,2. Default nsig=[0 0];
%        tsig      = a (1 x 2) array. If the t-value is less than tsig(i),
%                    the parameter is eliminated after the i-th stage of
%                    the HR method and the model is reestimated, i=1,2.
%                    Default tsig=[.75 1.].
% If the three stages of the Hannan-Rissanen method are performed,
% residuals based on the difference equation are also obtained using the
% third stage estimates .
%---------------------------------------------------
% RETURNS: str = a structure containing the estimated parameters with the
% following fields (see function MHANRIS)
%         s: number of outputs (neqs)
%         m: number of inputs  (nx)
%       kro: a (1 x s) vector containing the Kronecker indices
%       phi: an (s x s x maxkro) array with NaNs as parameters
%     theta: an (s x s x maxkro) array with NaNs as parameters
%     gamma: an (s x m x maxkro) array with NaNs as parameters
%     nparm: number of parameters
%      npar: an (s x s) array to define the Kronecker indices
%         F: an (n x n) matrix with NaNs as parameters, where n is the
%            McMillan degree = sum of the Kronecker indices
%         H: an (s x n) matrix with NaNs as parameters
%         K: an (n x s) matrix with NaNs as parameters
%         B: an (n x m) matrix with NaNs as parameters
%         D: an (s x m) matrix with NaNs as parameters
%    residv: an (nobs x s) matrix containing the residuals obtained in the
%            first stage
%   sigmarv: an (s x s) covariance matrix of residv
%      vgam: a {[(2*nlag + 1)*s^2 + (nlag + 1)*s*m + neqs] x 1} vector
%            containing the stacks of phi (except the first matrix), theta,
%            gamma, and s NaNs to account for the mean, where nlag =
%            max(kro).
%      bind: an [(nparm + s) x 1] index vector for the parameters in vgam.
%      beta: an [(nparm + s) x 1] vector containing the parameters
%            estimated in the second stage
%        tv: an [(nparm + s) x 1] vector containing the t-values of beta
%     vgams: a vector like vgam but with the NaNs replaced with the
%            parameters estimated in the second stage.
%    vgamtv: a vector like vgam but with the NaNs replaced with the
%            t-values corresponding to the second stage.
%   noninv2: = 1, if model is noninvertible after the second stage
%            = 0, if model is invertible after the second stage
%    nonst2: = 1, if model is nonstationary after the second stage
%            = 0, if model is stationary after the second stage
%    resid2: an [(nobs-nlag) x s] matrix containing the residuals of the
%            second stage regression
%   sigmar2: covariance matrix of resid2
%    musers: mean corresponding to the constant estimated in the second
%            stage
%      phis: same as phi but with the NaNs replaced with the parameters
%            estimated in the second stage
%     phitv: same as phi but with the Nans replaced with the t-values
%            corresponding to the second stage
%    thetas: same as theta but with the NaNs replaced with the parameters
%            estimated in the second stage
%   thetatv: same as theta but with the Nans replaced with the t-values
%            corresponding to the second stage
%    gammas: same as gamma but with the NaNs replaced with the parameters
%            estimated in the second stage
%   gammatv: same as gamma but with the Nans replaced with the t-values
%            corresponding to the second stage
%        mu: an (s x 1) vector containing the constant estimated in the
%            second stage.
%      mutv: an (s x 1) vector containing the t-values of the constant
%            estimated in the second stage.
%     phist: same as phis but with coefficient matrices premultiplied by
%            phis(:,:,1)^{-1} (VARMAX model not in echelon form)
%   thetast: same as thetas but with coefficient matrices premultiplied by
%            phis(:,:,1)^{-1} (VARMAX model not in echelon form)
%   gammast: same as gammas but with coefficient matrices premultiplied by
%            phis(:,:,1)^{-1} (VARMAX model not in echelon form)
%        Fs: same as F but with the NaNs replaced with the parameters
%            estimated in the second stage
%        Hs: same as H but with the NaNs replaced with the parameters
%            estimated in the second stage
%        Ks: same as K but with the NaNs replaced with the parameters
%            estimated in the second stage
%        Bs: same as B but with the NaNs replaced with the parameters
%            estimated in the second stage
%        Ds: same as D but with the NaNs replaced with the parameters
%            estimated in the second stage
%   noninv3: = 1, if model is noninvertible after the third stage
%            = 0, if model is invertible after the third stage
%    nonst3: = 1, if model is nonstationary after the third stage
%            = 0, if model is stationary after the third stage
%   resid23: an (nobs x s) matrix of residuals obtained before the third
%            stage using the VARMAX difference equation estimated in the
%            second stage starting with zeros.
%  sigmar23: covariance matrix of resid23
%     beta3: same as beta but containing the parameters
%            estimated in the third stage
%       tv3: same as tv but containing the t-values corresponding to the
%            third stage
%    vgams3: same as vgams but with the parameters estimated in the third
%            stage
%   vgamtv3: same as vgamtv but with the t-values corresponding to the
%            third stage
%      mus3: same as mu but with the constant estimated in the third stage
%     phis3: same as phis but containing the parameters estimated in the
%            third stage
%   thetas3: same as thetas but containing the parameters estimated in the
%            third stage
%   gammas3: same as gammas but containing the parameters estimated in the
%            third stage
%    phitv3: same as phitv but containing the t-values corresponding to the
%            third stage
%  thetatv3: same as thetatv but containing the t-values corresponding to
%            the third stage
%  gammatv3: same as gammatv but containing the t-values corresponding to
%            the third stage
%     mutv3: same as mutv but containing the t-values corresponding to the
%            third stage
%    phist3: same as phist but containing the parameters estimated in the
%            third stage
%  thetast3: same as thetast but containing the parameters estimated in the
%            third stage
%  gammast3: same as gammast but containing the parameters estimated in the
%            third stage
%       Fs3: same as Fs but containing the parameters estimated in the
%            third stage
%       Ks3: same as Ks but containing the parameters estimated in the
%            third stage
%       Bs3: same as Bs but containing the parameters estimated in the
%            third stage
%       Ds3: same as Ds but containing the parameters estimated in the
%            third stage
%       Hs3: same as Hs but containing the parameters estimated in the
%            third stage
%    resid3: an [(nobs-nlag) x s] matrix of residuals corresponding to the
%            third stage regression
%   sigmar3: covariance matrix of resid3
%   resid3m: an (nobs x s) matrix containing the residuals obtained using
%            the VARMAX difference equation estimated in the third stage
%            starting with zeros.
%  sigmar3m: covariance matrix of resid3m
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

ferror = 0;
str = [];

[ny, s] = size(y);
if ~isempty(x)
    [nx, m] = size(x);
    if (nx ~= ny)
        ferror = 1;
        disp('estvarmaxpqrPQR: nobs in x-matrix not the same as y-matrix');
        return
    end
else
    m = 0;
end

[nor, mor] = size(ordersr);
[nos, mos] = size(orderss);

if (seas == 1) && ((mor ~= 3) || (mos ~= 3))
    ferror = 1;
    disp('ordersr and orderss should be a 1 x 3 array');
    return
end

if nargin < 6
    hr3 = 1;
    finv2 = 0;
    mstainv = 0;
    nsig = [0, 0];
    tsig = [.75, 1.];
end
if nargin < 7
    finv2 = 0;
    mstainv = 0;
    nsig = [0, 0];
    tsig = [.75, 1.];
end
if nargin < 8
    mstainv = 0;
    nsig = [0, 0];
    tsig = [.75, 1.];
end
if nargin < 9
    nsig = [0, 0];
    tsig = [.75, 1.];
end
if nargin < 10
    tsig = [.75, 1.];
end


%transform multiplicative to non-multiplicative model and put restrictions
%in the model
str = restrcmodel(s, m, seas, ordersr, orderss);

%estimate restricted model
str = mhanris(y, x, seas, str, hr3, finv2, mstainv, nsig, tsig); 

s2 = s * s;
if (hr3 == 0) && (finv2 == 1) %third step of the HR method
    if str.nonst3 == 1
        [mphi, nphi, pp1] = size(str.phi);
        p = pp1 - 1;
        cont = 0;
        for i = 1:p
            cont = cont + sum(sum(isnan(str.phi(:, :, i+1))));
        end
        if (mstainv == 1) && (seas == 1) && (cont == s2 * p)
            str = sta3(str);
        else
            str = sta3r(str);
        end
    end
    if str.noninv3 == 1
        [mtheta, ntheta, qp1] = size(str.theta);
        q = qp1 - 1;
        cont = 0;
        for i = 1:q
            cont = cont + sum(sum(isnan(str.theta(:, :, i+1))));
        end
        if (mstainv == 1) && (seas == 1) && (cont == s2 * q)
            str = inv3(str);
        else
            str = inv3r(str);
        end
    end
    % third step residuals obtained using the difference equation
    str3 = str;
    nlag = sum(ordersr) + sum(orderss);
    if (nlag > 0)
        str3.thetast = str.thetast3;
        str3.phist = str.phist3;
        str3.gammast = str.gammast3;
        str3.mu = str.mus3;
        str3.phis = str3.phis3;
        [resid3, sigmar3] = compresde0(y, x, str3);
    else
        resid3 = str.resid3;
        sigmar3 = str.sigmar3;
    end
    str.resid3m = resid3;
    str.sigmar3m = sigmar3;
else %second step of the HR method
    if str.nonst2 == 1
        [mphi, nphi, pp1] = size(str.phi);
        p = pp1 - 1;
        cont = 0;
        for i = 1:p
            cont = cont + sum(sum(isnan(str.phi(:, :, i+1))));
        end
        if (mstainv == 1) && (seas == 1) && (cont == s2 * p)
            str = sta2(str);
        else
            str = sta2r(str);
        end
    end
    if str.noninv2 == 1
        [mtheta, ntheta, qp1] = size(str.theta);
        q = qp1 - 1;
        cont = 0;
        for i = 1:q
            cont = cont + sum(sum(isnan(str.theta(:, :, i+1))));
        end
        if (mstainv == 1) && (seas == 1) && (cont == s2 * q)
            str = inv2(str);
        else
            str = inv2r(str);
        end
    end
end
