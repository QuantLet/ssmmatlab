%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 12.3.5 of Lütkepohl (2005), pp. 477-479
%
% Series are:  German income and consumption, s.a.,
% 1960Q1-1982Q4; source: Deutsche Bundesbank
% The number of observations is ny=91-16=75.
% For the first differences in logs, a VARMA(2,2) in echelon form
% with k1=0, k2=2, is estimated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y = load(fullfile('data', 'vf_lutincon.dat'));
y = diferm(log(y), 1);
%sample is 1960Q1-1978Q4
y = y(1:75, :);
[ny, my] = size(y);
%subtract mean
y = y - repmat(mean(y), ny, 1);
x = [];

lag = 36;
cw = 1.96;
freq = 1;
ds = 0;
dr = 0;
tname = {'Income', 'Consumption'};
for i = 1:my
    for ds = 0:0
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all

seas = 1;

disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause


%identify a VARMA(p,q) model for the series
maxlag = [];   
minlag = 0;
prt = 1;
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


%estimate the Kronecker indices for the series
prt = 1;
maxorder = [];   
hr3 = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

disp(' ')
disp('estimation using the Hannan-Rissanen method')
disp('Kronecker indices are [0 1]')
disp('press any key to continue')
pause


%estimate model using HR method (K.i. = [0 1])
kro = [0, 1];
hr3 = 0;
finv2 = 1;
strv = estvarmaxkro(y, x, seas, kro, hr3, finv2);

disp(' ');
disp('***** Estimated  simplified VARMA Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 0;
mprintar(strv.phis3, in, tit, strt);
disp(' ')
tit = 'th';
strt = 0;
mprintar(strv.thetas3, in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ');
disp('*****t-values  *****');
disp(' ');
tit = 'tv-phi';
strt = 0;
mprintar(strv.phitv3, in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 0;
mprintar(strv.thetatv3, in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
in.fmt = char('%12.6f');
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause


disp(' ')
disp('estimation using the conditional method')
disp('Kronecker indices are [0 1]')
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
mprintar(strc.phiscon, in, tit, strt);
tit = 'th';
strt = 0;
mprintar(strc.thetascon, in, tit, strt);
disp(' ');
tit = 'Sigma:';
in.fmt = char('%12.6f');
mprintar(strc.sigmarcon, in, tit);
disp(' ')
disp('t-values: ')
disp(' ');
in.fmt = char('%12.4f');
tit = 'tv-mu:';
mprintar(strc.mutvcon', in, tit);
tit = 'tv-phi';
strt = 0;
mprintar(strc.phitvcon, in, tit, strt);
tit = 'tv-th';
strt = 0;
mprintar(strc.thetatvcon, in, tit, strt);
disp('press any key to continue')
pause


disp(' ')
disp('estimation using the exact method')
disp('Kronecker indices are [0 1]')
disp('press any key to continue')
pause

%estimate using the exact method
Y = [];
[xvf, str, ferror] = mexactestim(y, x, strv, Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 0;
mprintar(str.phisexct, in, tit, strt);
disp(' ')
tit = 'th';
strt = 0;
mprintar(str.thetasexct, in, tit, strt);
disp(' ')
tit = 'Sigma';
in.fmt = char('%12.6f');
mprintar(str.sigmarexct, in, tit);
disp(' ')
disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
tit = 'tv-phi';
strt = 0;
mprintar(str.phitvexct, in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 0;
mprintar(str.thetatvexct, in, tit, strt);
disp(' ')
