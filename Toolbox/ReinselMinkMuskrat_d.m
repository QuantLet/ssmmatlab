%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reinsel (1997), pp. 96 and 164. Series of logarithms of the annual sales
% of mink and muskrat furs by the Hudson's Bay Company for the years
% 1850-1911, with T=62 annual observations. Using the LR sequential tests,
% an AR(3) model is selected. Then, ARMA(p,q) models for p=1,2,3 and
% q=0,1,2 were estimated. An ARMA(2,1) is selected.
% On p. 168, an echelon form is identified using Akaike's method. This form
% has Kronecker indices k1=2, k2=1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

yy = load(fullfile('data', 'lminkmus.dat'));
y = yy(:, 2:3);
x = [];

seas = 1;
[ny, s] = size(y);

lag = 36;
cw = 1.96;
freq = 1;
ds = 0;
dr = 0;
tname = {'Mink', 'Muskrat'};
for i = 1:s
    for ds = 0:0
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all

lag = 10;
ic = 1;
disp(' ')
disp('******** cross correlation matrices:     ********');
str = mautcov(y, lag, ic);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('press any key to continue')
pause

disp(' ')
disp('identify a VAR model for the series')
disp(' ')
%VAR order identification
prt = 1;
minlag = 0;
maxlag = 6;
lagsopt = varident(y, maxlag, minlag, prt);
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a VAR of order 3:')
disp('press any key to continue')
pause
test = 1;
lags = 3;
res = var_est(y, lags, test);

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
disp(' ')
tit = 'Sigma:';
mprintar(res.sigmar, in, tit);

disp(' ');
disp('***** Estimated t-values  *****');
disp(' ');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ')
tit = 'Granger caus. prob.:';
mprintar(res.fprob, in, tit);

disp(' ')
disp('press any key to continue')
pause


disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause


%identify a VARMA(p,q) model for the series
maxlag = -2; %we fix the maximum order to two
minlag = 0;
prt = 0;
[lagsopt, ferror] = lratiopqr(y, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause

disp(' ')
disp('identify Kronecker indices')
disp(' ')
disp('press any key to continue')
pause


%estimate the Kronecker indices for the original series
maxorder = -2; %we fix the maximum order to two
hr3 = 0;
prt = 1;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

disp(' ')
disp('estimate model in Echelon Form and eliminate')
disp('some insignificant parameters: ')
disp('Kronecker Indices are [2 1].                             ')
disp('press any key to continue')
pause


kro = [2, 1];
%estimate model using HR method (K.i. = [2 1]) and eliminate some
%insignificant parameters
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [1., 1.];
strv = estvarmaxkro(y, x, seas, kro, hr3, finv2, mstainv, nsig, tsig);


disp(' ');
disp('***** Estimated Model using the HR method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 0;
mprintar(strv.phis3(:, :, 1:3), in, tit, strt);
disp(' ')
tit = 'th';
strt = 0;
mprintar(strv.thetas3(:, :, 1:3), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 0;
mprintar(strv.phitv3(:, :, 1:3), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 0;
mprintar(strv.thetatv3(:, :, 1:3), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause


disp(' ')
disp('estimation using the conditional method')
disp('press any key to continue')
pause
%estimate using the conditional method
[xvfc, strc, ferrorc] = mconestim(y, x, strv);

disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strc.muscon', in, tit);

tit = 'phi';
strt = 0;
mprintar(strc.phiscon(:, :, 1:3), in, tit, strt);

tit = 'th';
strt = 0;
mprintar(strc.thetascon(:, :, 1:3), in, tit, strt);

disp(' ');
tit = 'Sigma:';
mprintar(strc.sigmarcon, in, tit);

disp(' ')
disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strc.mutvcon', in, tit);

tit = 'tv-phi';
strt = 0;
mprintar(strc.phitvcon(:, :, 1:3), in, tit, strt);

tit = 'tv-th';
strt = 0;
mprintar(strc.thetatvcon(:, :, 1:3), in, tit, strt);
disp('press any key to continue')
pause


disp(' ')
disp('estimation using the exact method')
disp('Kronecker indices are [2 1]')
disp('press any key to continue')
pause

%estimate using the exact method
Y = repmat(eye(s), ny, 1);
[xvf, strx, ferror] = mexactestim(y, x, strv, Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strx.musexct', in, tit);

tit = 'phi';
strt = 0;
mprintar(strx.phisexct(:, :, 1:3), in, tit, strt);

tit = 'th';
strt = 0;
mprintar(strx.thetasexct(:, :, 1:3), in, tit, strt);
disp(' ');
tit = 'Lh:';
mprintar(strx.Lh, in, tit);
disp(' ')
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);
disp(' ')


disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strx.mutvexct', in, tit);

tit = 'tv-phi';
strt = 0;
mprintar(strx.phitvexct(:, :, 1:3), in, tit, strt);

tit = 'tv-th';
strt = 0;
mprintar(strx.thetatvexct(:, :, 1:3), in, tit, strt);

tit = 'tv-Lh:';
mprintar(strx.Lhtv3, in, tit);
disp(' ')
disp('press any key to continue')
pause


disp(' ')
disp('******** Computation of Recursive Residuals   ********');
%compute recursive residuals
%set up regression matrices
X = eye(s);
W = [];
Sigmax = strx.sigmarexct;
[L, p] = chol(Sigmax, 'lower');
Lm1 = pinv(L);
%set up system matrices
T = strx.Fsexct;
Z = strx.Hsexct;
G = Lm1;
H = strx.Ksexct * Lm1;
%set up initial conditions
ndelta = 0; %number of unit roots
[ins, i, ferror] = incossm(T, H, ndelta);

% [Xt,Pt,g,M,initf,recrs]=scakff(y,X,Z,G,W,T,H,ins,i);
[Xt, Pt, g, M, initf, recrs, recr] = scakffsqrt(y, X, Z, G, W, T, H, ins, i);
%plot recursive residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 12;
ic = 1;
nr = strx.nparm;
disp(' ')
disp('******** Residuals:     ********');
str = mautcov(recr, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('Q statistics:')
disp(str.qstat)

disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all
