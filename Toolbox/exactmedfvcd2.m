function [ff, xv, e, f, str, stx] = exactmedfvcd2(xv, xf, y, x, str, pr, qr, r, Y, chb)
% PURPOSE: given a structure, it computes first the SSF including the unit
% roots. Then, we compute y_t = V_t + U_t, where V_t is the exogenous part 
% (without the initial state part). This last part will be estimated along
% with the other regression parameters and an added mean. 
% We store in a structure the model for U_t, the regression estimator and  
% its mse if Y is not empty, the state space at the end of filtering and 
% its mse, and the state space model corresponding to the exogenous part, 
% V_t.
% The state space echelon form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = Y_{t}*beta + H*alpha_{t} + D*x_{t}  + a_{t}
%
%
%---------------------------------------------------
% USAGE: [ff,xv,e,f,str,stx,recrs]=exactmedfvcd2(xv,y,x,str,Y,chb)
% where:    xv   = the parameter vector
%           y      = an (nobs x neqs) matrix of y-vectors
%           str    = a structure containing the model information
%           Y      = an (nobs x (neqs x nbeta)) regression matrix other
%                       than the mean
%                    (neqs x nbeta) if it is time invariant
%           chb    = 1 compute hb and Mb in Kalman filter
%                    0 do not compute hb and Mb in Kalman filter
%                    2 compute recursive residuals updating regression p.
%                    3 compute recursive residuals with fixed regression p.
%---------------------------------------------------
% RETURNS: ff   = a vector containing the individual functions at the
%                 solution
%          xv = the parameter vector, possibly modified
%          e    = a vector containing the standardized residuals
%          f    = a scalar containing the determinantal term
%          str  = the input structure str, possibly modified
%          stx  = a structure containing the following fields
%       .X,.Z,.G,.W,.T,.H,.ins,.i are the matrices and initial conditions
%       information corresponding to the state space model for the
%       endogenous part
%         .hb   = the beta estimator
%         .Mb   = the Mse of the beta estimator
%         .A    = the estimated augmented state vector at the end of
%                 filtering
%         .P    = the Mse of A at the end of filtering
%        .kro   = the kronecker indices for the model
%         recrs = standardized recursive residuals
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
ff = [];
e = [];
f = [];
A = [];
P = [];
Mb = [];
%obtain differenced series and parameters for the VARMA model followed by
%the differenced series.
[yd, xvv, xff, DA, Dr, Ds, ferror] = pr2varmapqPQd(y, xv, xf, str);

%compute first the VARMAX model with unit roots and then obtain the SSF
nr = str.nr;
phit = pmatmul(str.phisexct,Dr);  %phi polynomial contains unit roots
[~, m] = size(x);
[nz, s] = size(y);
if (nr > 0)
    prd = pr + 1;
else
    prd = pr;
end
maxgr = max([prd qr r]);
kro = [maxgr maxgr maxgr];
%generate structure for SSF with Kronecker indices equal to the maximum
%degree
[strr, ferror] = matechelon(kro, s, m); 
strr.nr = nr;
strr.phis = phit;
strr.thetas = str.thetasexct;
strr.gammas = str.gammasexct;
for i = prd + 2: maxgr + 1
    strr.phis(:,:,i)= zeros(s);
end
for i = qr + 2: maxgr + 1
    strr.thetas(:,:,i)= zeros(s);
end
for i = r + 2: maxgr + 1
    strr.gammas(:,:,i)= zeros(s,m);
end
strr = armaxe2sse(strr);  %transfor VARMAX to SSF
%add some information
strr.hb = str.musexct;
strr.Mb = str.mutvexct;
strr.sigma2c = str.sigma2c;
strr.beta3 = str.beta3;
strr.bind = str.bind;
strr.sigmarexct = str.sigmarexct;
strr.xvd = str.xvd;


%We use z_t = V_t + U_t, where V_t is the exogenous part. We estimate U
T = strr.Fs;
Z = strr.Hs;
if m > 0
    inc = 1;
    Bs = strr.Bs;
    Ds = strr.Ds;
    kro = strr.kro;
    [z, rx1] = varmafil(x, T, Z, Bs, Ds, kro, inc);
    %subtract filtered inputs from output
    yf = y - z;
    U = yf;
else
    U = y;
    Bs = [];
    Ds = [];
    rx1 = [];
end
%regression matrix for regression variables and initial state
Ya = [Y, rx1];


%compute covariance matrix of residuals using the last parameters in
%beta
xv = strr.beta3;
nxv = length(xv);
nxvd = length(strr.xvd);
xvv = xv(nxvd+1:nxv);    
bind = strr.bind;
[nbind, mbind] = size(bind); %nbind includes the mean parameters
nparma = nbind - s; %number of parameters in ar and ma parts
Lparm = xvv(nparma+1:end); %save parameters for the Cholesky factor L
L = zeros(s);
l = 0;
betam = [1., Lparm'];
for i = 1:s
    cont = s - i + 1;
    ind = l + 1:l + cont;
    L(i:end, i) = betam(ind);
    l = l + cont;
end
%set up the SSF matrices before calling the KF that will be applied to the
%series filtered with inputs
T = strr.Fs;
Z = strr.Hs;
HH = strr.Ks;
HH = HH * L;
GG = L;
W = [];
X = Ya; 
[ins, ii, ferror] = incossm(T, HH, nr); 


%We use y_t = V_t + U_t, where V_t is the exogenous part. 
if (chb == 0) || (chb == 1)
    %likelihood evaluation
    [e, f, hb, Mb, A, P] = scakflesqrt(U, X, Z, GG, W, T, HH, ins, ii, chb);
end

if m > 0
    %subtract the effect of the input initial state
    [~,nbeta] = size(Y);
    [~,nm1] = size(rx1);
    m1 = hb(nbeta + 1:nbeta + nm1);
    yfvec = vec(U'); %put the observations in a vector
    yfvec = yfvec - Ya(:,nbeta + 1:end) * hb(nbeta + 1:end); 
    yx = reshape(yfvec, s, nz);    
    U = yx';
    X = Y;
    %filter again without the input initial effect
    [e, f, hb, Mb, A, P] = scakflesqrt(U, X, Z, GG, W, T, HH, ins, ii, chb);
end

%Store in a
%structure the model for U_t, the regression estimator and its mse if Y is
%not empty, the state space at the end of filtering and its mse, and the
%state space model corresponding to the exogenous part, V_t.
stx.X = Y;
stx.Z = Z;
stx.G = GG;
stx.W = W;
stx.T = T;
stx.H = HH;
stx.ins = ins;
stx.i = ii;

stx.hb = hb;
stx.Mb = Mb;
stx.A = A;
stx.P = P;
stx.B = Bs;
stx.D = Ds;
if (m > 0)
    stx.m1 = m1;
else
    stx.m1 = [];
end
stx.kro = strr.kro;
nxvv = length(xvv);
nxv = length(xv);
if nxv > nxvv
    xv(nxv-nxvv+1:nxv) = xvv;
else
    xv = xvv;
end
