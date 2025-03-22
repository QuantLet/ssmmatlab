%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Examples 7.4.3 and 14.6 of Lütkepohl (2005), pp. 312-314 and pp. 526-528
%
% Series are U.S. macroeconomic data: 1) real money stock M1, 2) GNP in
% billions of 1982 dollars, 3),the discount interest rate on new issues of
% 91-days treasury bills and 4) the yield on long term (20 years) treasury
% bonds. Logarithms of seasonally adjusted GNP and M1 data are used.
% The period is: 1954-1 through 1987-4.
% The number of observations is ny=136 for Example 7.4.3
% The data are truncated for the last four years for Example 14.6.
% The number of observations is ny=136-16=120 for Example 14.6.
% In Example 14.6, a model in echelon form for the 'differenced series' is
% estimated instead of an error correction model in reverse echelon form as
% in Lütkepohl pp. 526-528.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y = load(fullfile('data', 'e3.dat'));
y = y(1:136, :);

%take logs of M1 and GNP
y(:, 1) = log(y(:, 1));
y(:, 2) = log(y(:, 2));
[ny, my] = size(y);

lag = 32;
cw = 1.96;
freq = 4;
ds = 0;
tname = {'M1', 'GNP', 'discount', 'yield'};
for i = 1:my
    for dr = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all

disp(' ')
disp('estimation of a VAR(2) model for the original series')
disp('')
disp('press any key to continue')
pause


%First, we estimate a VAR of order 2
test = 0;
lags = 2;
res = var_est(y, lags, test);

disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);

disp(' ');
disp('***** Estimated t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ');
tit = 'Sigma:';
in.fmt = char('%12.7f');
mprintar(res.sigmar, in, tit);


% %Matlab-Econ function for Johansen cointegration test
% lags=1:2;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);

[D, nr, yd, DA, ferror] = mcrcregr(y);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause


x = [];
seas = 1;
%number of unit roots in the model
nr = 3;
prt = 0;
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(y, x, prt, nr);

disp(' ')
disp('identify a VAR model for the differenced series')
disp(' ')
disp('press any key to continue')
pause

%identify a VAR model for the differenced series
maxlag = 3;
minlag = 1;
prt = 1;
lagsopt = varident(yd, maxlag, minlag, prt);
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a VAR(2) model for the differenced series')
disp('')
disp('press any key to continue')
pause

%estimation of the model for the differenced series
%model for the differenced series: VAR(2).
hr3 = 0;
finv2 = 1;
mstainv = 1;
pr = 2;
qr = 0; %VAR(2)
[str, ferror] = estvarmaxpqrPQR(yd, x, seas, [pr, qr, 0], [0, 0, 0], hr3, finv2, mstainv);

phi(:, :, 1) = str.phis3(:, :, 1);
for i = 1:pr
    phi(:, :, i+1) = str.phis3(:, :, i+1);
end
Phi(:, :, 1) = phi(:, :, 1);
th(:, :, 1) = str.thetas3(:, :, 1);
for i = 1:qr
    th(:, :, i+1) = str.thetas3(:, :, i+1);
end
Th(:, :, 1) = phi(:, :, 1);
Sigma = str.sigmar3;
clear str

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 7.4.3
% Estimation of the VAR(2) model parameterized in error correction form
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

disp(' ')
disp('estimation of the VAR(2) model in error correction form')
disp('press any key to continue')
pause
%estimate the model in error correction form
%
Y = [];
constant = 1;
result = varmapqPQestime(y, str, Y, constant);

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
[Lambdaf, alphaf, betapf, thf, Thf, Lf, ferror] = pr2ecf(xvv, xff, DAf, str);
[phif, De, DAe, ferror] = mecf2mid(Lambdaf, alphaf, betapf);
phit = pmatmul(phif, Dr);

disp(' ');
disp('***** Estimated overall AR part  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(phit(:, :, 2:3), in, tit, strt);
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phif, 3) * result.h)', in, tit);
disp('press any key to continue')
pause


%error correction matrices
disp(' ');
disp('***** AR part in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betapf);
alphaf = alphaf * betapf(:, 1:nb);
betapf = betapf(:, 1:nb) \ betapf;
disp('estimated beta transposed')
disp(betapf)
disp('estimated alpha')
disp(alphaf)
disp('press any key to continue')
pause
clear str
%
% End of Example 7.4.3
% Estimation of the VAR(2) model parameterized in error correction form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 7.4.3
% Estimation of the VAR(2) model parameterized in terms of the model for
% the differenced series
%

%set up varma model for the differenced series. Parameters are defined in
%terms of the model for the differenced series, contained in phi, th, Phi,
%Th and Sigma, and the parameters of the differencing matrix polynomial,
%contained in DA.
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, seas);

%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DA);

disp(' ')
disp('estimation of the  VAR(2) model for the differenced series. ')
disp('Note that the overall model has changed. It is now VAR(3).  ')
disp('press any key to continue')
pause

%estimate model with cointegration relations imposed
Y = [];
constant = 1;
result = varmapqPQestimd(y, str, Y, constant);

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
phit = pmatmul(phif, Dr);

disp(' ');
disp('***** Estimated overall AR part  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(phit(:, :, 2:end), in, tit, strt);
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phif, 3) * result.h)', in, tit);
disp('press any key to continue')
pause


%error correction matrices
[Pi, Lambda, alpha, betafp, ferror] = mid2mecf(phif, Dr, DAf);
disp(' ');
disp('***** AR part in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambda(:, :, 2:3), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafp);
alpha = alpha * betafp(:, 1:nb);
betafp = betafp(:, 1:nb) \ betafp;
disp('estimated beta transposed')
disp(betafp)
disp('estimated alpha')
disp(alpha)
disp('press any key to continue')
pause

%
% end of Example 7.4.3
% Estimation of the VAR(2) model parameterized in terms of the model for
% the differenced series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 14.6
%
% In this example, we estimate a model in echelon form for the differenced
% series instead of a model in reverse echelon form as in Lütkepohl (2006),
% pp. 526-528
%

%reduce sample size for Example 14.6
y = y(1:120, :);

%compute new differencing matrix polynomial and new 'differenced series'
prt = 0;
[D, DA, yd, ferror] = mdfestim1r(y, x, prt, nr);


%add information about the Kronecker indices
%estimate the Kronecker indices for the original series
maxorder = -3;  %we fix the maximum order 
hr3 = 0;
prt = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the original series ')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause


%estimate the Kronecker indices for the differenced series
maxorder = -3;  %we fix the maximum order 
[order, kro, scm] = varmaxscmidn(yd, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the differenced series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause


disp(' ')
disp('estimate model in Echelon Form for the Differenced Series: ')
disp('Kronecker Indices are [1 0 0 1]                   ')
disp('press any key to continue')
pause

%estimate VARMA model in echelon form for the differenced series and
%aliminate some nonsignificant parameters
kro = [1, 0, 0, 1];
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 0];
tsig = [1., 0.];
strv = estvarmaxkro(yd, x, seas, kro, hr3, finv2, mstainv, nsig, tsig);


%add unit root information to the model structure
[strv, ferror] = aurirvarmapqPQ(strv, nr, DA);
%add freq to the model structure
strv.freq = seas;
%estimate model
s = size(y, 2);
Y = [];
constant = 1;
[xvfk, strv, ferror] = mexactestimcd(y, strv, Y, constant);

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
strt = 0;
mprintar(Lambdak(:, :, 1:2), in, tit, strt);
tit = 'Theta';
strt = 0;
mprintar(strv.thetasexct(:, :, 1:2), in, tit, strt);
disp('t-values of Phi and Theta matrices:')
tit = 'tv-Phi';
strt = 0;
mprintar(strv.phitvexct(:, :, 1:2), in, tit, strt);
tit = 'tv-Theta';
strt = 0;
mprintar(strv.thetatvexct(:, :, 1:2), in, tit, strt);
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
mprintar((sum(phifk, 3) * strv.musexct)', in, tit);

%
% End of Example 14.6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
