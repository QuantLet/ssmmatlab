%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 7.2.6 of Lütkepohl (2005), pp. 302-305 and p. 309
%
% sample: 1972Q2 -- 1998Q4
% West German data until 1990Q2, all of Germany data afterwards
% R - nominal long term interest rate (Umlaufsrendite)
%                    (source: Monatsberichte der Deutschen Bundesbank,
%                     quarterly values are values of last month of quarter)
% Dp - \Delta log gdp deflator (source: Deutsches Institut für Wirtschaftsforschung,
%                               Volkswirtschaftliche Gesamtrechnung)
% In Lütkepohl (2005), pp. 302-304, a VAR model with 3 lags and
% cointegration rank equal to one is estimated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y = load(fullfile('data', 'e6.dat'));
y = y(1:107, :);
y = fliplr(y);

[ny, my] = size(y);

freq = 4; % quarterly data
bg_year = 1972;
bg_per = 2;
datei = cal(bg_year, bg_per, freq);
Y = seasdmom(ny, datei);

lag = 36;
cw = 1.96;
ds = 0;
tname = {'Dp', 'R'};
for i = 1:my
    for dr = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all


%ys: preliminary seasonally adjusted series
betah = mulols(y, Y);
ys = y - Y * betah;

disp(' ')
disp('identify a VAR model for the preliminary seasonally adjusted series')
disp(' ')
disp('press any key to continue')
pause


%identify a VAR model for the preliminary seasonally adjusted series
maxlag = 6;
minlag = 1;
prt = 1;
lagsopt = varident(ys, maxlag, minlag, prt);
disp('press any key to continue')
pause

% %Matlab-Econ function for Johansen cointegration test
% lags=lagsopt;
% [h,pValue,stat,cValue,mles] = jcitest(ys,'lags',lags);

[D, nr, yd, DA, ferror] = mcrcregr(ys);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause


x = [];
seas = 1;
%number of unit roots in the model
nr = 1;
prt = 0;
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(ys, x, prt, nr);

disp(' ')
disp('identify a VAR model for the differenced series')
disp(' ')
disp('press any key to continue')
pause


%identify a VAR model for the differenced series
maxlag = 6;
minlag = 1;
prt = 1;
lagsopt = varident(yd, maxlag, minlag, prt);
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a VAR(4) model for the differenced series')
disp('')
disp('press any key to continue')
pause

%model for the differenced series: VAR(4).
finv2 = 1;
mstainv = 1;
hr3 = 0;
pr = 4;
qr = 0;
[str, ferror] = estvarmaxpqrPQR(yd, x, seas, [pr, qr, 0], [0, 0, 0], hr3, finv2, mstainv);
phi(:, :, 1) = str.phis3(:, :, 1);
for i = 2:pr + 1
    phi(:, :, i) = str.phis3(:, :, i);
end
Phi(:, :, 1) = phi(:, :, 1);
th(:, :, 1) = str.thetas3(:, :, 1);
Th(:, :, 1) = phi(:, :, 1);
Sigma = str.sigmar3;
clear str


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 7.2.6
% Estimation of the model parameterized in error correction form
%
%pass AR part to error correction form
phit = pmatmul(phi, Phi);
[Pi, Lambda, alpha, betap, ferror] = mid2mecf(phit, D, DA);
%re-specify Lambda so that the model coincides with that of Lütkepohl
%(2006)
Lambda = Lambda(:, :, 1:end-1);


%set up varma model in error correction form. Parameters are defined in
%terms of the error correction form.
[str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, Sigma, seas);

%estimate the model in error correction form
%
% Set up matrix X so that it can be used in the state space model with
% regression variables. The observation equation of this model should be of
% the form
%
% y_t = X_t*beta + u_t,                                              (1)
%
% where beta is the vector of regression coefficients and u_t follows a
% VAR(4) model. If we write model (1) in the form
%
% y_t = C*S_t + u_t,
%
% where S_t = (s_{1t},s_{2t},s_{3t})' is the vector of seasonal dummies, we
% can apply the vec operator so that vec(C*S_t) = kron(S'_t,I_k)*vec(C).
% This implies X_t = kron(S'_t,I_k) and beta = vec(C), where k is the
% dimension of y_t. Applying the differencing operator, Dr(B), to (1), it
% is obtained that
%
% Dr(B)*y_t = Dr(B)*X_t*beta + Dr(B)*u_t,                            (2)
%
% where Dr(B)*u_t follows a VAR(3) model. We can add a constant to model
% (2) so that the model to be estimated is
%
% Dr(B)*y_t = mu + Dr(B)*X_t*beta + Dr(B)*u_t.                        (3)
%
% We have to pass to function varmapqPQestime the matrix X containing the
% X_t and a parameter constant = 1,0 indicating whether we want a constant
% in the model or not, respectively. The differencing polynomial is passed
% in the structure str.
% Note that model (3) is different from the model considered in Lütkepohl
% 7.2.6 for this example with respect to the treatment of seasonal dummies.
% In our model, if Phi(B)u_t = A_t is the VAR(4) or VAR(5) model followed
% by u_t, it holds that
%
% Phi(B)y_t = Phi(B)C*S_t + A_t,
%
% whereas in Lütkepohl 7.2.6 the equation is
%
% Phi(B)y_t = C*S_t + A_t.
%
% However, it is shown in Lütkepohl pp. 334-335 that both models are
% equivalent.
%

disp(' ')
disp('estimation of the VAR(4) model in error correction form')
disp('press any key to continue')
pause

% matrix X is (n*p x nbeta)
[nb, mb] = size(betah);
X = ones(ny*my, nb*mb);
for i = 1:ny
    ip = (i - 1) * my + 1:i * my;
    X(ip, :) = kron(Y(i, :), eye(my));
end
constant = 1;
result = varmapqPQestime(y, str, X, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause


%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(y, xvrf, xrf, str);
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvv, xff, str);
[Lambdaf, alphaf, betapf, thf, Thf, Lf, ferror] = pr2ecf(xvv, xff, DA, str);
disp(' ');
disp('***** AR part in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2:4), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betapf);
alphaf = alphaf * betapf(:, 1:nb);
betapf = betapf(:, 1:nb) \ betapf;
disp('estimated beta transposed')
disp(betapf)
disp('estimated alpha')
disp(alphaf)
%constant is phif(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
in.fmt = char('%12.4f');
mprintar((sum(phif, 3) * result.h(1:my))', in, tit);
% C matrix
C = reshape(result.h(my+1:end), my, freq-1);
disp(' ')
disp('C matrix in model: y_t = C*S_t + u_t')
disp(C)

disp(' ')
disp('Insignificant parameters are fixed to zero, ')
disp('cointegration relation is fixed to [1 -4], ')
disp('and the model is reestimated')
disp('press any key to continue')
pause

%fix insignificant parameters to zero using the field Lambdan
clear str
Lambdaf(2, 1, 2) = 0.0;
Lambdaf(1, 1, 3) = 0.0;
Lambdaf(2, 1, 3) = 0.0;
Lambdaf(2, 1, 4) = 0.0;
DAf = DA;
DAf(2, 1) = .25;
betapf = m2mor(DAf);
[str, ferror] = suvarmapqPQe(Lambdaf, alphaf, betapf, thf, Thf, result.Sigmar, seas);
str.Lambdan(2, 1, 2) = 0.0;
str.Lambdan(1, 1, 3) = 0.0;
str.Lambdan(2, 1, 3) = 0.0;
str.Lambdan(2, 1, 4) = 0.0;
str.DAn(2, 1) = .25;
[str, ferror] = fixvarmapqPQe(str);
%estimate model again
resultf = varmapqPQestime(y, str, X, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(resultf)
disp('press any key to continue')
pause


%create estimated model
xvrf = resultf.xvf;
xrf = resultf.xf;
[ydf, xvv, xff, DAff, Dr, Ds, ferror] = pr2varmapqPQd(y, xvrf, xrf, str);
[phiff, thff, Phiff, Thff, Lff, ferror] = pr2varmapqPQ(xvv, xff, str);
[Lambdaff, alphaff, betapff, thff, Thff, Lff, ferror] = pr2ecf(xvv, xff, DAff, str);
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaff(:, :, 2:4), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betapff);
alphaff = alphaff * betapff(:, 1:nb);
betapff = betapff(:, 1:nb) \ betapff;
disp('estimated beta transposed')
disp(betapff)
disp('estimated alpha')
disp(alphaff)
%constant is phif(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(resultf.Sigmar, in, tit);
tit = 'Constant:';
in.fmt = char('%12.4f');
mprintar((sum(phiff, 3) * resultf.h(1:my))', in, tit);
% C matrix
disp(' ')
C = reshape(resultf.h(my+1:end), my, freq-1);
disp('C matrix in model: y_t = C*S_t + u_t')
disp(C)
%
% End of Example 7.2.6
% Estimation of the model parameterized in error correction form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 7.2.6
% Estimation of the model parameterized in terms of the model for the
% differenced series
%


%set up varma model for the differenced series. Parameters are defined in
%terms of the model for the differenced series, contained in phi, th, Phi,
%Th and Sigma, and the parameters of the differencing matrix polynomial,
%contained in DA.
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, seas);


%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DA);

% Set up matrix Y so that it can be used in the state space model with
% regression variables. The observation equation of this model should be of
% the form
%
% y_t = X_t*beta + u_t,                                              (1)
%
% where beta is the vector of regression coefficients and u_t follows a
% VAR(4) or VAR(5) model. If we write model (1) in the form
%
% y_t = C*S_t + u_t,
%
% where S_t = (s_{1t},s_{2t},s_{3t})' is the vector of seasonal dummies, we
% can apply the vec operator so that vec(C*S_t) = kron(S'_t,I_k)*vec(C).
% This implies X_t = kron(S'_t,I_k) and beta = vec(C), where k is the
% dimension of y_t. Applying the differencing operator, Dr(B), to (1), it
% is obtained that
%
% Dr(B)*y_t = Dr(B)*X_t*beta + Dr(B)*u_t,                            (2)
%
% where Dr(B)*u_t follows a VAR(3) model. We can add a constant to model
% (2) so that the model to be estimated is
%
% Dr(B)*y_t = mu + Dr(B)*X_t*beta + Dr(B)*u_t.                        (3)
%
% We have to pass to function varmapqPQestimd the matrix X containing the
% X_t and a parameter constant = 1,0 indicating whether we want a constant
% in the model or not, respectively. The differencing polynomial is passed
% in the structure str.
% Note that model (3) is different from the model considered in Lütkepohl
% 7.2.6 for this example with respect to the treatment of seasonal dummies.
% In our model, if Phi(B)u_t = A_t is the VAR(4) or VAR(5) model followed
% by u_t, it holds that
%
% Phi(B)y_t = Phi(B)C*S_t + A_t,
%
% whereas in Lütkepohl 7.2.6 the equation is
%
% Phi(B)y_t = C*S_t + A_t.
%
% However, it is shown in Lütkepohl pp. 334-335 that both models are
% equivalent.
%

disp(' ')
disp('estimation of the model in terms of the model for ')
disp('the differenced series. Note that the overall model')
disp('has changed. It is now VAR(5).                     ')
disp('press any key to continue')
pause


% matrix X is (n*p x nbeta)
[nb, mb] = size(betah);
X = ones(ny*my, nb*mb);
for i = 1:ny
    ip = (i - 1) * my + 1:i * my;
    X(ip, :) = kron(Y(i, :), eye(my));
end
constant = 1;
result = varmapqPQestimd(y, str, X, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(y, xvrf, xrf, str);
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvv, xff, str);


%error correction matrices
[Pi, Lambda, alpha, betafp, ferror] = mid2mecf(phif, Dr, DAf);
disp('estimated Lambda matrix')
disp(Lambda)
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambda(:, :, 2:5), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafp);
alpha = alpha * betafp(:, 1:nb);
betafp = betafp(:, 1:nb) \ betafp;
disp('estimated beta transposed')
disp(betafp)
disp('estimated alpha')
disp(alpha)
%constant is phif(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
in.fmt = char('%12.4f');
mprintar((sum(phif, 3) * result.h(1:my))', in, tit);
% C matrix
C = reshape(result.h(my+1:end), my, freq-1);
disp(' ')
disp('C matrix in model: y_t = C*S_t + u_t')
disp(C)


% Estimation of the model parameterized in terms of the model for the
% differenced series

%estimate error correction model in echelon form

disp(' ')
disp('estimate Kronecker indices for the original and differenced series')
disp('press any key to continue')
pause

%add information about the Kronecker indices
%estimate the Kronecker indices for the original series
maxorder = [];  
hr3 = 0;
prt = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);

disp('estimated Kronecker Indices for the original series: ')
disp(kro)
disp('press any key to continue')
pause
%estimate the Kronecker indices for the differenced and seasonally
%corrected series
maxorder = [];  
prt = 0;
[order, kro, scm] = varmaxscmidn(ydf, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the differenced series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

disp(' ')
disp('estimate model in echelon form for the differenced series.')
disp('Kronecker Indices are [1 3].')
disp('press any key to continue')
pause

%preliminary estimation
%estimate VARMA model in echelon form for the 'differenced and seasonally
%corrected series'
kro = [1, 3];
hr3 = 0;
finv2 = 1;
mstainv = 1;
strv = estvarmaxkro(yd, x, seas, kro, hr3, finv2, mstainv);
%add unit root information to the model structure
[strv, ferror] = aurirvarmapqPQ(strv, nr, DA);
%add freq to the model structure
strv.freq = seas;
%estimate error correction model in echelon form
[xvfk, strv, ferror] = mexactestimcd(y, strv, X, constant);

%create estimated model
xffk = [];
[ydd, xvvk, xffk, DAfk, Drk, Dsk, ferrork] = pr2varmapqPQd(y, xvfk, xffk, strv);
phifk = strv.phisexct;
thfk = strv.thetasexct;
[nphi, mphi, prk] = size(phifk);

%error correction matrices
[Pik, Lambdak, alphak, betafpk, ferror] = mid2mecf(phifk, Drk, DAfk);
disp(' ');
disp('***** Estimated model in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdak(:, :, 2:4), in, tit, strt);
tit = 'Theta';
strt = 1;
mprintar(strv.thetasexct(:, :, 2:4), in, tit, strt);
disp('t-values of Phi and Theta matrices:')
tit = 'tv-Phi';
strt = 1;
mprintar(strv.phitvexct(:, :, 2:4), in, tit, strt);
tit = 'tv-Theta';
strt = 1;
mprintar(strv.thetatvexct(:, :, 2:4), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafpk);
alphak = alphak * betafpk(:, 1:nb);
betafpk = betafpk(:, 1:nb) \ betafpk;
disp('estimated beta transposed')
disp(betafpk)
disp('estimated alpha')
disp(alphak)
%constant is phifk(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(strv.sigmarexct, in, tit);
tit = 'Constant:';
mprintar((sum(phifk, 3) * strv.musexct(1:my))', in, tit);
% C matrix
C = reshape(strv.musexct(my+1:end), my, freq-1);
disp(' ')
disp('C matrix in model: y_t = C*S_t + u_t')
disp(C)
%
% end of Example 7.2.6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
