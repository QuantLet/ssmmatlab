%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.1 of Tsay (2014), pp. 335-341
%
% Monthly housing data of the United States from
% January 1963 to July 2012. The two series employed are
% 1. z1t: Logarithm of new homes sold in thousands of units (new residential sales)
% 2. z2t: Logarithm of the total new privately owned housing units started in
%             thousands of units (new residential construction)
% The aim of this exercise is to identify a model for the differenced
% series. We first identify the regular part. To this end, a seasonal (1,1)
% model is estimated. With the residuals of this model, a series of
% likelihood ratio tests is performed to identify the regular part. With
% the regular part fixed, we perform a series of likelihood ratio tests to
% identify the seasonal part using the residuals of the regular model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

da = load(fullfile('data', 'm-hsoldhst6312.dat'));
zt = log(da(:, 3:4));
tdx = da(:, 1) + da(:, 2) / 12;
[nx, mx] = size(zt);
freq = 12;

% subplot(2,1,1)
% plot(tdx,zt(:,1))
% xlabel('time'); ylabel('Hsold'); axis('tight');
% subplot(2,1,2)
% plot(tdx,zt(:,2))
% xlabel('time'); ylabel('Hstart'); axis('tight');
% disp('press any key to continue')
% pause
% close all


% lag=36; cw=1.96;  dr=0;
% tname={'US housing sold','US housing starts'};
% for i=1:2
%  for ds=0:1
%   c0=sacspacdif(zt(:,i),tname(i),dr,ds,freq,lag,cw);
%   pause
%  end
% end
% close all


yd = diferm(zt, freq); %seasonal differencing
yd = diferm(yd, 1); %regular differencing

%compute autocovariance and autocorrelation matrices
lag = 24;
ic = 1;
nr = 0;
disp(' ')
disp('******** Sample cross correlation matrices      ********');
disp('******** for the differenced series:            ********');
stre = mautcov(yd, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(stre.r0)


disp(' ')
disp('estimate a seasonal VARMA(1,1) model using the Hannan-Rissanen ')
disp('method to obtain filtered series')
disp('press any key to continue')
pause


%estimate a seasonal VARMAX(1,1) model by the Hannan-Rissanen method.
x = [];
hr3 = 0;
finv2 = 1;
mstainv = 1;
qr = 0;
ps = 1;
qs = 1;
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [0, qr, 0], [ps, qs, 0], hr3, finv2, mstainv);


%identify a VARMA(p,q) model (regular part) for the residual series
maxlag = 4;
minlag = 0;
prt = 1;
x = [];
seas = 1;
[lagsopt, ferror] = lratiopqr(strv.resid3, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause

disp(' ')
disp('estimate a regular VARMA(0,3) model using the Hannan-Rissanen ')
disp('method to obtain filtered series')
disp('press any key to continue')
pause

%estimate a regular VARMAX(0,3) model by the Hannan-Rissanen method.
x = [];
hr3 = 0;
finv2 = 1;
mstainv = 1;
qr = 3;
ps = 0;
qs = 0;
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [0, qr, 0], [ps, qs, 0], hr3, finv2, mstainv);

lag = 24;
ic = 1;
nr = 0;
disp(' ')
disp('******** Sample cross correlation matrices      ********');
disp('******** for the filtered series:               ********');
stre = mautcov(strv.resid3, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(stre.r0)

disp('Clearly, a seasonal model (0,1)_12 is specified')
disp('press any key to continue')
disp(' ')
pause


disp('Identify a model using sequential likelihood ratio tests')
disp('for the filtered series')

% identify a VARMA(P,Q)_s model (seasonal part) for the residual series.
% To obtain estimates of the innovations, a long VARX model is estimated
% first. Then, sequential likelihood ratio tests are performed to obtain an
% optimal VARMAX(i,i,i), starting with i=minlag until i=maxlag.
maxlag = 24;
minlag = 12;
prt = 1;
x = [];
seas = 12;
[lagsopt, ferror] = lratiopqr(strv.resid3, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause
