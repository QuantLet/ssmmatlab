function str = mhanris(y, x, seas, str, hr3, finv2, mstainv, nsig, tsig)
% PURPOSE: applies the Hannan-Rissanen method for VARMAX models in echelon
% form with restrictions to the series y with inputs x. It returns a
% structure containing the estimated model. The state space echelon form
% is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = H*alpha_{t} + D*x_{t}  + a_{t}
%
% The VARMAX echelon form is

%  phi(B)*y_{t} = gamma(B)*x_{t} + theta(B)a_{t},
%
% where B is the backshift operator, B*y_{t} = y_{t-1}.
%---------------------------------------------------
% USAGE: str = mhanris(y,x,seas,str,hr3,finv2,mstainv,nsig,tsig)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%        seas      = seasonality, =1 no seasonality, >1 seasonality
%         str      = a structure containing model information as given by
%                    function MATECHELON on output (model in echelon form)
%          hr3      = 1, perform only the first two stages of the HR method
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
%---------------------------------------------------
% RETURNS: str = a structure containing the the following
%                fields
%         s: number of outputs (neqs)
%         m: number of inputs  (nx)
%       kro: a (1 x s) vector containing the Kronecker indices
%       phi: an (s x s x maxkro) array with NaNs as parameters
%     theta: an (s x s x maxkro) array with NaNs as parameters
%     gamma: an (s x m x maxkro) array with NaNs as parameters
%     nparm: number of parameters
%      npar: an (s x s) array to define the Kronecker indices
%         F: an (n x n) matrix with NaNs as parameters
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
        error('mhanris: nobs in x-matrix not the same as y-matrix');
    end
else
    m = 0;
end

if nargin < 6
    finv2 = 0;
    mstainv = 0;
end

if nargin < 7
    mstainv = 0;
end


if nargin < 8
    nsig2 = 0;
    nsig3 = 0;
else
    nsig2 = nsig(1);
    nsig3 = nsig(2);
end
if nargin < 9
    tsig2 = .75;
    tsig3 = 1.;
else
    tsig2 = tsig(1);
    tsig3 = tsig(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First step of the HR method:
% parameters for the long VARX
if (m > 0)
    minp = 5;
else
    minp = 3;
end
if (seas > 1)
    a = 2.;
else
    a = 1.25;
end
pt = max(8, seas+minp);
minlags = 0;
maxlags = ceil(max(log(ny)^a, pt)); %VARX length is given by this formula

%%%%%%%%%%%%%%%%%%%%%%%
%compute VAR residuals

% compute initial residuals using VARXs with different lengths
initres = zeros(maxlags, s);
for i = minlags:maxlags - 1
    if m == 0
        resid = var_res(y, i);
        initres(i+1, :) = resid(1, :);
    else
        resid = varx_res(y, i, x);
        initres(i+1, :) = resid(1, :);
    end
end
% compute the rest of residuals using a VARX(maxlags)
if m == 0
    [res] = var_res(y, maxlags);
else
    [res] = varx_res(y, maxlags, x);
end

residv = [initres(1:maxlags, :); res]; % use all possible VARX residuals

%  %alternative to fixed VARX order
%  prt=1;
%  if m == 0
%   [lagsopt,initres] = lratiocr(y,maxlags,minlags,prt);
%   res = var_res(y,lagsopt);
%  else
%   [lagsopt,initres] = lratiocrx(y,maxlags,minlags,prt,x);
%   res = varx_res(y,lagsopt,x);
%  end
%  fprintf(1,'p in VAR(p) = %2d\n',lagsopt);
%  residv=[initres(1:lagsopt,:); res]; % use all possible VARX residuals
%  %end of alternative

% % use residuals of shorter length than series
% residv=res;  [nres,mres]=size(res);
% y=y(end-nres+1:end,:);
% if m > 0
%  x=x(end-nres+1:end,:);
% end


% %computation of residuals using Whittle's algorithm. If inputs are present,
% %we should first regress the output on the inputs and then apply Whittle's
% %algorithm to the residuals of this regression.
% %
% str=mautcov(y,maxlags,0);c0=str.c0; cv=str.cv;
% strw=whittle(c0,cv);
% % we have to write this function: resid=reswhittle(y,strw);
%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second step of HR method:
str = mhanris2(y, residv, x, str); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsignificant parameter elimination in second step of HR method.
% Put nsig=1 and choose tsig for insignificant t-value
% nsig2=0; tsig2=.75;
%parameter elimination after the second step of HR method
if nsig2
    [ct2, str] = nse2(y, residv, x, tsig2, str); 
    str.ct2 = ct2;
end
if finv2 == 1
    noninv2 = str.noninv2;
    if (noninv2 > 0)
        if mstainv == 1
            str = inv2(str);
        else
            str = inv2r(str);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third step of HR method: only if model is invertible and hr3 = 0
if (str.noninv2 == 0) && (hr3 == 0)
    %  tic
    str = mhanris3(y, x, str);  
    %  toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nonsignificant parameter elimination in third step of HR method.
    % Put nsig=1 and choose tsig for insignificant t-value
    %parameter elimination after the third step of HR method
    if nsig3
        invert2 = 0;
        [ct3, str] = nse3(y, x, tsig3, invert2, str); 
        str.ct3 = ct3;  
    end
else
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
