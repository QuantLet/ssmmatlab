%Example of estimation of the flour price series using a
%VARMA model.
%
% Monthly Flour Price Indices for Three U.S. cities, Buffalo, Minneapolis,
% and Kansas City, for the Period August 1972 Through November 1980. The
% three series are I(1) and it seems that there are no cointegration
% relationships. These series have been used by Tiao y Tsay (1989) and
% Lütkepohl and Poskitt (1996).

%load data
y = load(fullfile('data', 'flour-price.dat'));
x = [];
%logs are taken
y = log(y);
seas = 1;

lag = 20;
cw = 1.96;
ds = 0;
dr = 1;
tname = {'flour price 1', 'flour price 2', 'flour price 3'};
for i = 1:3
    c0 = sacspacdif(y(:, i), tname(i), dr, ds, seas, lag, cw);
    pause
end
close all

disp(' ')
disp('identify a VAR model for the series')
disp(' ')
disp('press any key to continue')
pause


%identify a VAR model for the series
maxlag = 6;
minlag = 1;
prt = 1;
lagsopt = varident(y, maxlag, minlag, prt);
disp('press any key to continue')
pause

disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause


%identify a VARMA(p,q) model for the series
maxlag = [];
minlag = 0;
prt = 0;
[lagsopt, ferror] = lratiopqr(y, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause


% %Matlab-Econ function for Johansen cointegration test
% lags=lagsopt;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);

[D, nr, yd, DA, ferror] = mcrcregr(y);
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause

disp(' ')
disp('identify Kronecker indices')
disp(' ')
disp('press any key to continue')
pause


%estimate the Kronecker indices for the series
maxorder = [];
hr3 = 0;
prt = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the original series ')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause
