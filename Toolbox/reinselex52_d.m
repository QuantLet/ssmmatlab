%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 5.2 of Reinsel (1997), pp. 170-174
%
% Series are: 1) weekly production schedule figures and 2) billing figures.
% The number of observations is n=100.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y = load(fullfile('data', 'Weeklyproshed.dat'));

[ny, s] = size(y);

lag = 20;
cw = 1.96;
freq = 1;
ds = 0;
tname = {'weekly production schedule figures', 'billing figures'};
for i = 1:2
    for dr = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all

%sample correlation matrices
lag = 10;
ic = 1;
disp('******** Original Series:     ********');
str = mautcov(y, lag, ic);
disp('press any key to continue')
pause

%Preliminary VAR analysis
prt = 1;
minlag = 0;
maxlag = 6;
lagsopt = varident(y, maxlag, minlag, prt);
pause


%Estimate a VAR of order 5
disp('estimation of a VAR of order 5:')
disp('press any key to continue')
pause
test = 1;
lags = 5;
res = var_est(y, lags, test);

disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2:6), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);

disp(' ');
disp('*****t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2:6), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ');
tit = 'Sigma:';
mprintar(res.sigmar, in, tit);
disp(' ');
in.cnames = char(' Granger causality prob.:', ' ');
mprint(res.fprob, in);
disp(' ')
disp('press any key to continue')
pause


lag = 12;
ic = 1;
nr = s^2 * lags;
disp('******** VAR residuals:     ********');
str = mautcov(res.resid, lag, ic, nr);
disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all

%identify a VARMA(p,q) model for the series
maxlag = 5;
minlag = 0;
prt = 0;
x = [];
seas = 1;
[lagsopt, ferror] = lratiopqr(y, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause


disp(' ')
disp('estimate VARMA(4,1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause


%estimate the model with two variables
%First, estimate a VARMAX(4,1,0) model by the Hannan-Rissanen method.
x = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(y, x, freq, [4, 1, 0], [0, 0, 0], hr3, finv2);

%Then, impose restrictions and estimate using the HR method again
strv.phi(1, 1, 2:5) = zeros(1, 4);
strv.phi(1, 2, 2:5) = zeros(1, 4);
strv.phi(2, 1, 2:3) = zeros(1, 2);
strv.phi(2, 2, 4:5) = zeros(1, 2);
strv.theta(1, 2, 2) = 0;
strv.theta(2, 1, 2) = 0;
strv.nparm = strv.nparm - 14;
strv = mhanris(y, x, freq, strv, hr3, finv2);
disp(' ');
disp('***** Estimated VARMA(4,1) Model using the *****');
disp('***** HR method after fixing some parameters  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strv.phis3(:, :, 2:5), in, tit, strt);
disp(' ');
tit = 'th';
strt = 1;
mprintar(strv.thetas3(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause


%estimate using exact ML
%setup model
Phi = eye(2);
Th = eye(2);
phi = strv.phis3;
th = strv.thetas3(:, :, 1:2);
Sigma = strv.sigmar3;

%create regression variable for the mean
Y = eye(2);
%fix insignificant parameters to zero and set up model
Sigma(2, 1) = 0;
Sigma(1, 2) = 0;
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);
str.phin(1, 1, 2:5) = zeros(1, 4);
str.phin(1, 2, 2:5) = zeros(1, 4);
str.phin(2, 1, 2:3) = zeros(1, 2);
str.phin(2, 2, 4:5) = zeros(1, 2);
str.thn(1, 2, 2) = 0;
str.thn(2, 1, 2) = 0;
str.Lhn(2) = 0;
[str, ferror] = fixvarmapqPQ(str);

disp(' ')
disp('estimation using the exact method')
disp('press any key to continue')
pause

%estimate model using the exact method
result = varmapqPQestim(y, str, Y);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause


%estimated and fixed parameters
xvf = result.xvf;
xf = result.xf;
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr

%create estimated model
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvf, xf, str);
Sigmar = result.Sigmar;
%obtain t-values
[phift, thft, Phift, Thft, Lft, ferrort] = pr2varmapqPQ(result.tv, xf, str);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 0;
mprintar(phif(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(thf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Mean';
mprintar(result.h', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(Sigmar, in, tit);
disp(' ')
disp(' ')
disp('t-values: ')
tit = 'tv-phi';
strt = 1;
mprintar(phift(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(thft(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Mean';
mprintar(result.tvr', in, tit);
disp(' ')
tit = 'tv-L';
mprintar(Lft, in, tit);
disp('press any key to continue')
pause


%compute recursive residuals
[strf, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, Sigmar, freq);
%set up regression matrices
X = Y;
W = [];
%set up system matrices
T = strf.T;
Z = strf.Z;
G = strf.G;
H = strf.H;
%set up initial conditions
ndelta = 0; %number of unit roots
[ins, i, ferror] = incossm(T, H, ndelta);


[Xt, Pt, g, M, initf, recrs, recr] = scakff(y, X, Z, G, W, T, H, ins, i);
%plot recursive residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 12;
ic = 1;
nr = length(xvf) - s * (s + 1) / 2 + 1;
disp(' ')
disp('******** Recursive Residuals:     ********');
str = mautcov(recr, lag, ic, nr);
disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all

disp(' ')
disp('******** Second part: Alignment *************')
disp(' ')


%alignment
y1a = y(:, 1);
y1a(3:end) = y(1:98, 1);
ya = [y1a(3:end), y(3:end, 2)];

%Preliminary VAR analysis
prt = 1;
minlag = 0;
maxlag = 6;
lagsoptlr = varident(ya, maxlag, minlag, prt);
disp('Preliminary VAR analysis: lag lenght estimated with LR')
disp(lagsoptlr)
pause

disp('estimation of a VAR of order 3 as in Reinsel p.173:')
disp('press any key to continue')
pause

%Estimate a VAR of order 3 as in Reinsel p.173
test = 1;
lags = 3;
res = var_est(ya, lags, test);

disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);

disp(' ');
disp('*****t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ');
tit = 'Sigma:';
mprintar(res.sigmar, in, tit);
disp(' ');
% in.cnames = char(' Granger causality prob.:',' ');
% mprint(res.fprob,in);
tit = 'Granger causality prob.:';
mprintar(res.fprob, in, tit);
disp(' ')
disp('press any key to continue')
pause


lag = 12;
ic = 1;
nr = s^2 * lags;
disp('******** VAR residuals:     ********');
str = mautcov(res.resid, lag, ic, nr);
disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all


%identify a VARMA(p,q) model for the series
maxlag = 3;
minlag = 0;
prt = 0;
x = [];
seas = 1;
[lagsopt, ferror] = lratiopqr(ya, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause


%estimate a VARMA(2,1) by the Hannan-Rissanen method
hr3 = 0;
finv2 = 1;
x = [];
[strv, ferror] = estvarmaxpqrPQR(ya, x, freq, [2, 1, 0], [0, 0, 0], hr3, finv2);

disp('estimation of a VARMA(2,1) model using the conditional method')
disp('press any key to continue')
pause


%estimate the model using the conditional method
[xvfc, strc, ferrorc] = mconestim(ya, x, strv);
%
disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strc.muscon', in, tit);
tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 2:3), in, tit, strt);
tit = 'th';
strt = 1;
mprintar(strc.thetascon(:, :, 2), in, tit, strt);
disp(' ');
tit = 'Sigma:';
mprintar(strc.sigmarcon, in, tit);
disp(' ')
disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strc.mutvcon', in, tit);
tit = 'tv-phi';
strt = 1;
mprintar(strc.phitvcon(:, :, 2:3), in, tit, strt);
tit = 'tv-th';
strt = 1;
mprintar(strc.thetatvcon(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause


%fix all insignifican parameters to zero and estimate again
strv.phi(1, 1:2, 2) = zeros(1, 2);
strv.phi(1, 1:2, 3) = zeros(1, 2);
strv.theta(2, 1, 2) = 0;
strv.nparm = strv.nparm - 5;
strv = mhanris(ya, x, freq, strv, hr3, finv2);

disp('estimation using the conditional method again ')
disp('after fixing some parameters')
disp('press any key to continue')
pause


%estimate the model using the conditional method
[xvfc, strc, ferrorc] = mconestim(ya, x, strv);


disp(' ');
disp('***** Estimated VARMA(2,1) Model using the conditional *****');
disp('***** method after fixing some parameters              *****');
disp(' ');

clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strc.muscon', in, tit);
tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 2:3), in, tit, strt);
tit = 'th';
strt = 1;
mprintar(strc.thetascon(:, :, 2), in, tit, strt);
disp(' ');
tit = 'Sigma:';
mprintar(strc.sigmarcon, in, tit);
disp(' ')
disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strc.mutvcon', in, tit);
tit = 'tv-phi';
strt = 1;
mprintar(strc.phitvcon(:, :, 2:3), in, tit, strt);
tit = 'tv-th';
strt = 1;
mprintar(strc.thetatvcon(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause


lag = 12;
ic = 1;
nr = length(xvf) - s * (s + 1) / 2 + 1;
disp(' ')
disp('******** Conditional Residuals:     ********');
str = mautcov(strc.residcon, lag, ic, nr);
disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all
