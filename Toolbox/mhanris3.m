function str = mhanris3(y, x, str)
% PURPOSE: performs the third step of the multivariate Hannan-Rissanen
% method for VARMAX models with restrictions and returns the estimated parameters
%---------------------------------------------------
% USAGE: str = mhanris3(y,x,str)
% where:    y      = an (nobs x neqs) matrix of y-vectors
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
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

if nargin ~= 3
    error('wrong # of arguments to mhanris3');
end;

[nobs, neqs] = size(y); 
if ~isempty(x)
    [nobs2, nx] = size(x);
    if (nobs2 ~= nobs)
        error('mhanris3: nobs in x-matrix not the same as y-matrix');
    end
else
    nx = 0;
end

s = str.s;
m = str.m;
kro = str.kro;

if (s ~= neqs)
    error('mhanris3: s nonequal dim(y_t)');
end
if (m ~= nx)
    error('mhanris3: m nonequal dim(x_t)');
end

str.noninv3 = 0;
str.nonst3 = 0;
nlag = max(kro);

if nlag == 0 %the only parameter is the mean
    str.resid3 = str.resid2;
    str.sigmar3 = str.sigmar2;
    return
end

%check invertibility
if ~isempty(str.noninv2) && (str.noninv2 == 1)
    str.resid3 = str.resid2;
    str.sigmar3 = str.sigmar2;
    return %if model is noninvertible, we do not continue
end


% %check stationarity
% nonst2=0;
% if ~isempty(str.nonst2)
%  if str.nonst2==1, nonst2=1; end
% end

if ~isfield(str, 'resid23') || (isfield(str, 'resid23') && isempty(str.resid23))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new second step residuals. We compute them starting with zeros. Do not
    % use the Kalman filter at this stage.
    [resid2, sigmar2] = compresde0(y, x, str); % constant is passed in str.mu
    % plot(resid2)
    % pause
    str.resid23 = resid2;
    str.sigmar23 = sigmar2; %store covariance matrix of residuals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    resid2 = str.resid23;
    sigmar2 = str.sigmar23;
end

% adjust nobs to feed the lags
nobse = nobs - nlag;

ylag = glags(y, nlag);
% reslag = glags(resid2,nlag);
% vp=y(nlag+1:end,1:neqs)-resid2(nlag+1:end,1:neqs);
if (nx > 0)
    xlag = ltflag(x, nlag);
else
    xlag = [];
end

vgam = str.vgam; % the parameters for the constant are at the end of vgam
bind = str.bind;
[nbind, mbind] = size(bind);

%
% Third step. Only if there is a moving average part
%
str.beta3 = [];
str.tv3 = [];
str.vgams3 = [];
str.vgamtv3 = [];
str.mus3 = [];
str.phis3 = [];
str.thetas3 = [];
str.gammas3 = [];
str.phitv3 = [];
str.thetatv3 = [];
str.gammatv3 = [];
str.mutv3 = [];


resid = resid2;
sigmar = sigmar2;
% parameters passed from second step
phis = str.phis;
thetast = str.thetast;
beta = str.beta;


reslag = glags(resid, nlag);
vp = y(nlag+1:end, 1:neqs) - resid(nlag+1:end, 1:neqs);

% form x-matrix, constant is included
% xmat = kron([-ylag -vp reslag xlag ones(nobse,1)],sparse(eye(neqs)));
A = [-ylag, -vp, reslag, xlag, ones(nobse, 1)];
[na, ma] = size(A);
xmat = zeros(na*neqs, ma*neqs);
for i = 1:na
    for j = 1:ma
        aa = repmat(A(i, j), 1, neqs);
        dd = diag(aa);
        xmat((i - 1)*neqs+1:i*neqs, (j - 1)*neqs+1:j*neqs) = dd;
    end
end


% form design matrix for the regression.
zr = [];
for i = 1:length(vgam)
    if isnan(vgam(i))
        zr = [zr, xmat(:, i)];
    end
end


[R, p] = chol(sigmar);
L = R';
if (p > 0)
    error('covariance matrix of residuals2 singular in mhanris3')
end
thzrt = [];
[nzr, mzr] = size(zr);
Phi0 = phis(:, :, 1);
for i = 1:nobse
    thzrt = [thzrt; Phi0 \ zr((i - 1)*neqs+1:i*neqs, :)];
end
Th = [];
for i = 2:nlag + 1
    Th = [Th, thetast(:, :, i)];
end

V = zeros(nobs*neqs, mzr); % partial derivatives
for i = nlag + 1:nobs
    Vp = [];
    for j = i - 1:-1:i - nlag
        Vp = [Vp; V((j - 1)*neqs+1:j*neqs, :)];
    end
    V((i - 1)*neqs+1:i*neqs, :) = -Th * Vp - thzrt((i - nlag - 1)*neqs+1:(i - nlag)*neqs, :);
end
for i = nlag + 1:nobs
    V((i - 1)*neqs+1:i*neqs, :) = L \ V((i - 1)*neqs+1:i*neqs, :);
end

V = -V(nlag*neqs+1:end, :); % negative derivatives for the regression
lres2 = L \ resid(nlag+1:end, :)';
res2t = vec(lres2); % standardized residuals


res2t = res2t + V * beta; % regression
[beta3, tv3] = btval([], [res2t, V]);
str.beta3 = beta3;
str.tv3 = tv3;

% insert estimated parameters in vgam
vgams3 = vgam;
vgamtv3 = vgam;
for i = 1:nbind
    vgams3(bind(i)) = beta3(i);
    vgamtv3(bind(i)) = tv3(i);
end
str.vgams3 = vgams3;
str.vgamtv3 = vgamtv3;

%use an auxiliary structure to compute the state space and VARMAX
%echelon forms corresponding to the third step
str3 = str;
str3.vgams = str.vgams3;
str3.vgamtv = str.vgamtv3;
str3 = param2sse(str3);
%store results
str.phis3 = str3.phis;
str.thetas3 = str3.thetas;
str.gammas3 = str3.gammas;
str.phitv3 = str3.phitv;
str.thetatv3 = str3.thetatv;
str.gammatv3 = str3.gammatv;
str.mus3 = str3.mu;
str.mutv3 = str3.mutv;
str.phist3 = str3.phist;
str.thetast3 = str3.thetast;
str.gammast3 = str3.gammast;
str.Fs3 = str3.Fs;
str.Ks3 = str3.Ks;
str.Bs3 = str3.Bs;
str.Ds3 = str3.Ds;
str.Hs3 = str3.Hs;

% test for stationarity and invertibility after the third step
ima = chkstainv(str.Fs3-str.Ks3*str.Hs3); %if ima >1, the model is not invertible
%  report noninvertibility if any
str.noninv3 = 0;
if ima > 1
    str.noninv3 = 1;
end
iar = chkstainv(str.Fs3); %if iar >1, the model is not stationary
% report nonstationarity if any
str.nonst3 = 0;
if iar > 1
    str.nonst3 = 1;
end

% residuals of the third step regression
resid3v = res2t - V * beta3;
resid3 = zeros(nobse, neqs);
for i = 1:nobse
    resid3(i, :) = (L * resid3v((i - 1)*neqs+1:i*neqs))';
end
sigmar3 = cov(resid3, 1);
str.resid3 = resid3;
str.sigmar3 = sigmar3;
