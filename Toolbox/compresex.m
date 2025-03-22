function [resid, E, rSigmat] = compresex(y, x, str, tol, maxupdt, Y)
% PURPOSE: given a structure, it computes the model residuals and their
% covariance matrices using the square root CKMS recursions
%---------------------------------------------------
% USAGE: [resid2,sigmar2]=compresex(y,x,str)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%           str    = a structure containing the model information
%           Y      = an (nobs x (neqs x nbeta)) regression matrix
%---------------------------------------------------
% RETURNS: resid  = the residuals
%          E      = the augmented part of the residuals if Y is not empty
%         rSigmat = the Cholesky factors of the covariance matrix of
%                   residuals
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

[nobs, neqs] = size(y);
[nobs2, nx] = size(x);
kro = str.kro;
nlag = max(kro);
if (nargin < 6)
    Y = [];
end
if (nargin < 5)
    maxupdt = [];
end
if (nx > 0)
    r = nlag;
else
    r = 0;
end

if nlag == 0
    resid = y;
    E = [];
    rSigmat = [];
    return
else
    thetas = str.thetas;
    phis = str.phis;
    gammas = str.gammas;
    Sigma = str.sigmar2;
end

%if inputs are present, subtract filtered inputs from output
if r > 0
    %decouple the state space model. First, find input filter
    [phige, thetage, krog, ierror] = pecheform(phis, gammas);
    strg = matechelon(krog, neqs, nx);
    strg = vecparwr(strg);
    strg.phis = phige;
    strg.thetas = phige;
    strg.gammas = thetage;
    strg = armaxe2sse(strg);
    %Fs, Hs, Bs and Ds are the filter matrices for the inputs
    Fs = strg.Fs;
    Hs = strg.Hs;
    Bs = strg.Bs;
    Ds = strg.Ds;
    inc = 1;
    [z, rx1] = varmafil(x, Fs, Hs, Bs, Ds, krog, inc);
    %subtract filtered inputs from output
    yf = y - z; %first, without initial conditions
    Yx1 = [rx1, repmat(eye(neqs), nobs, 1)]; % add elements to estimate a mean in the regression
    %  Yx1=rx1;  %no mean
    yfvec = vec(yf'); % put the observations in a vector
    beta = Yx1 \ yfvec; 
    yfvec = yfvec - Yx1 * beta; %estimate intial conditions and subtract
    yx = reshape(yfvec, neqs, nobs);
    y = yx'; % transform from vector to matrix
end

%obtain filter for the innovations (decoupled if there are inputs)
if (nx > 0)
    [phie, thetae, kro, ierror] = pecheform(phis, thetas);
else
    phie = phis;
    thetae = thetas;
end
gammae = [];
m = 0;
stre = matechelon(kro, neqs, m);
stre = vecparwr(stre);
stre.phis = phie;
stre.thetas = thetae;
stre.gammas = gammae;
stre.sigmar2 = Sigma;
stre = armaxe2sse(stre); 
% Y=[];
[resid, E, rSigmat] = sqrt_ckms(y(1:end, :), Y, stre, maxupdt, tol);
