function [ff, beta, e, f, str, hb, Mb, resid, E] = exactmedfvd(beta, y, x, str, tol, Y, chb)
% PURPOSE: given a structure, it computes the functions in ff such that the
% expression ff'*ff is minimized in the Levenberg-Marquardt method. It uses
% the fast square root version of the Kalman filter.
% The state space form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = Y_{t}*beta + H*alpha_{t} + D*x_{t}  + a_{t}
%
%---------------------------------------------------
% USAGE: [ff,beta,e,f,str,hb,Mb]=exactmedfvd(beta,y,x,str,tol,Y,chb)
% where:    beta   = the parameter vector
%           y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%           str    = a structure containing the model information
%           tol    = tolerance for not updating in the square CKMS rec.
%           Y      = an (nobs x (neqs x nbeta)) regression matrix
%           chb    = 1   compute hb and Mb
%                    0 do not compute hb and Mb
%---------------------------------------------------
% RETURNS: ff   = a vector containing the individual functions at the
%                 solution
%          beta = the parameter vector, possibly modified
%          e    = a vector containing the standardized residuals
%          f    = a scalar containing the determinantal term
%         str   = the input structure str, possibly modified
%          hb   = the beta estimator
%          Mb   = the Mse of the beta estimator
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

%
% obtain differenced series and parameters for the VARMA model followed by
%the differenced series.
xf = [];
[yd, betav, xf, DA, Dr, Ds, ferror] = pr2varmapqPQd(y, beta, xf, str);
%betav contiene los parametros del modelo sin la parte de las raices
%unitarias
[ny, my] = size(y);
Yd = Y;
if ~isempty(Y)
    Xd = Y;
    [nxd, mxd] = size(Xd);
    if nxd == my
        %we make Y always time varying
        Xd = zeros(ny*my, mxd);
        for i = 1:ny
            ip = (i - 1) * my + 1:i * my;
            Xd(ip, :) = Y;
        end
    end
    if isfield(str, 'nr')
        nr = str.nr;
    else
        nr = 0;
    end
    if isfield(str, 'ns')
        ns = str.ns;
    else
        ns = 0;
    end
    if (ns > 0)
        Xdcold = [];
        for i = 1:mxd
            Xdcol = Xd(:, i);
            Xdcolt = reshape(Xdcol, ny, my);
            [XX, ferror] = dtimesy(Ds, Xdcolt);
            Xdcold = [Xdcold, vec(XX)];
        end
        Xd = Xdcold;
    end
    if (nr > 0)
        Xdcold = [];
        for i = 1:mxd
            Xdcol = Xd(:, i);
            Xdcolt = reshape(Xdcol, ny, my);
            [XX, ferror] = dtimesy(Dr, Xdcolt);
            Xdcold = [Xdcold, vec(XX)];
        end
        Xd = Xdcold;
    end
    Yd = Xd;
end
y = yd;
Y = Yd;
nr = str.nr;
[nx,~] = size(x);
if (nr > 0)
   x = x(2:nx,:);
end




[nobs, neqs] = size(y);
[nY, nbeta] = size(Y);



vgams = str.vgam;
bind = str.bind;
[nbind, mbind] = size(bind); %nbind includes the mean parameters
nparma = nbind - neqs; %number of parameters in ar and ma parts
Lparm = betav(nparma+1:end); %save parameters for the Cholesky factor L
for i = 1:nparma
    vgams(bind(i)) = betav(i);
end
str.vgams = vgams;
str = param2sse(str);

%check stationarity and invertibility. If necessary, change parameters
iar = chkstainv(str.Fs); %if iar >1, the model is not stationary
if iar > 1
    %  fprintf(1,'model nonstationary, iar = %2d\n',iar);
    %convert Phi(z) into Phi(lambda*z) for an appropriate lambda
    vgam = enfstab(str, 'phi  ');
    str.vgams = vgam;
    str = param2sse(str);
    for i = 1:nparma
        betav(i) = vgam(bind(i));
    end
    betav(nparma+1:end) = Lparm; %insert parameters for the Cholesky factor
end
ima = chkstainv(str.Fs-str.Ks*str.Hs); %if ima >1, the model is not invertible
if ima > 1
    %  fprintf(1,'model noninvertible, ima = %2d\n',ima);
    %convert Th(z) into Th(lambda*z) for an appropriate lambda
    vgam = enfstab(str, 'theta');
    str.vgams = vgam;
    str = param2sse(str);
    for i = 1:nparma
        betav(i) = vgam(bind(i));
    end
    betav(nparma+1:end) = Lparm; %insert parameters for the Cholesky factor
end


%compute covariance matrix of residuals using the last parameters in beta
L = zeros(neqs);
l = 0;
betam = [1, Lparm];
for i = 1:neqs
    cont = neqs - i + 1;
    ind = l + 1:l + cont;
    L(i:end, i) = betam(ind);
    l = l + cont;
end
sigmar = L * L';
str.sigmar2 = sigmar;


[resid, E, rSigmat] = compresex(y, x, str, tol, [], Y);
%compute ff
f = 1;
fc = 0;
e = zeros(nobs*neqs, 1);
SQT = [];
for ii = 1:nobs
    ind = (ii - 1) * neqs + 1:ii * neqs;
    V = rSigmat(ind, :);
    e(ind) = V \ resid(ii, :)';
    if nbeta > 0
        SQT = [SQT; V \ E((ii - 1)*nbeta+1:ii*nbeta, :)'];
    end
    f = f * abs(prod(diag(V)));
    [f, fc] = updatef(f, fc);
end


% qyy=[];
R = [];
SQT = [SQT, e];
if nbeta > 0
    [ns, ms] = size(SQT);
    [Q, R] = qr(SQT(:, 1:nbeta));
    qy = Q' * SQT(:, nbeta+1);
    %      qyy=qy(1:nbeta);
    e = qy(nbeta+1:ns);
    if chb == 1
        hb = R(1:nbeta, :) \ qy(1:nbeta);
        Mb = inv(R(1:nbeta, :));
        Mb = Mb * Mb';
    else
        hb = [];
        Mb = [];
    end
else
    e = SQT;
    hb = [];
    Mb = [];
end


nbsqs = nobs * neqs;
f = (f^(1 / (nbsqs))) * (2^(fc / (nbsqs)));
ff = e .* f;
