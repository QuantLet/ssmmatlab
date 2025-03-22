%script file for example3.9 in Tsay (2014), p. 175.
%

y = load(fullfile('data', 'ushog.dat'));

disp(' ')
disp('identify a VAR model for the series')
disp(' ')
disp('press any key to continue')
pause

%VAR order identification
prt = 1;
minlag = 0;
maxlag = 9;
lagsopt = varident(y, maxlag, minlag, prt);
pause

disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause

%identify a VARMA(p,q,r) model for the series
maxlag = [];
minlag = 0;
prt = 1;
x = [];
seas = 1;
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
prt = 0;
maxorder = [];
hr3 = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the series ')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause
