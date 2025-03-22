%Example of estimation of the simulated series used by Nsiri and Roy
%(1996), that follows a VARMA model in echelon form. The model is
%estimated by the conditional method.

clear
%load data
y = load(fullfile('data', 'nsiri.dat'));
x = [];
seas = 1;
[ny, s] = size(y);

freq = 1;
lag = 20;
cw = 1.96;
ds = 0;
dr = 0;
tname = {'Nsiri and Roy (1996) 1', 'Nsiri and Roy (1996) 2'};
for i = 1:2
    c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
    pause
end
close all

disp(' ')
disp('estimate Kronecker indices for the series')
disp('press any key to continue')
pause
%estimate the Kronecker indices for the original series
prt = 1;
maxorder = 2; 
hr3 = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the original series ')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

disp('estimate model in Echelon Form using the Hanna-Rissanen method: ')
disp('Kronecker Indices are [2 1].    ')
disp('press any key to continue')
pause


%estimate model using HR method (K.i. = [2 1]) and eliminate some
%nonsignificant parameters
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [1., 1.];
strv = estvarmaxkro(y, x, seas, [2, 1], hr3, finv2, mstainv, nsig, tsig);

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
strt = 1;
mprintar(strv.thetas3(:, :, 1:3), in, tit, strt);
disp(' ')
disp(' ');
tit = 'Sigma:';
mprintar(strv.sigmar3, in, tit);
disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
tit = 'tv-phi';
strt = 0;
mprintar(strv.phitv3(:, :, 1:3), in, tit, strt);
disp(' ');
tit = 'tv-th';
strt = 0;
mprintar(strv.thetatv3(:, :, 1:3), in, tit, strt);
disp('press any key to continue')
pause

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
% in.cnames = char(' phi(0):',' ',' phi(1):',' ',' phi(2):',' ');
% mprint(strc.phiscon,in);
tit = 'phi';
strt = 0;
mprintar(strc.phiscon, in, tit, strt);
disp(' ');
% in.cnames = char(' th(0):',' ',' th(1):',' ',' th(2):',' ');
% mprint(strc.thetascon,in);
tit = 'th';
strt = 0;
mprintar(strc.thetascon, in, tit, strt);
disp(' ');
% in.cnames = char(' Sigma:',' ');
% mprint(strc.sigmarcon,in);
tit = 'Sigma:';
mprintar(strc.sigmarcon, in, tit);
disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
% in.cnames = char(' tv-phi(0):',' ',' tv-phi(1):',' ',' tv-phi(2):',' ');
% mprint(strc.phitvcon,in);
tit = 'tv-phi';
strt = 0;
mprintar(strc.phitvcon, in, tit, strt);
disp(' ');
% in.cnames = char(' tv-th(0):',' ',' tv-th(1):',' ',' tv-th(2):',' ');
% mprint(strc.thetatvcon,in);
tit = 'tv-th';
strt = 0;
mprintar(strc.thetatvcon, in, tit, strt);
disp('press any key to continue')
pause

disp('estimation using the exact method')
disp('press any key to continue')
pause

%estimate using the exact method
Y = [];
[xvf, strx, ferror] = mexactestim(y, x, strv, Y);

disp(' ');
disp('***** Estimated Model using the exact method *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
% in.cnames = char(' phi(0):',' ',' phi(1):',' ',' phi(2):',' ');
% mprint(strx.phisexct,in);
tit = 'phi';
strt = 0;
mprintar(strx.phisexct, in, tit, strt);
disp(' ');
% in.cnames = char(' th(0):',' ',' th(1):',' ',' th(2):',' ');
% mprint(strx.thetasexct,in);
tit = 'th';
strt = 0;
mprintar(strx.thetasexct, in, tit, strt);
disp(' ');
% in.cnames = char(' Sigma:',' ');
% mprint(strx.sigmarexct,in);
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);
disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
disp(' ');
% in.cnames = char(' tv-phi(0):',' ',' tv-phi(1):',' ',' tv-phi(2):',' ');
% mprint(strx.phitvexct,in);
tit = 'tv-phi';
strt = 0;
mprintar(strx.phitvexct, in, tit, strt);
disp(' ');
% in.cnames = char(' tv-th(0):',' ',' tv-th(1):',' ',' tv-th(2):',' ');
% mprint(strx.thetatvexct,in);
tit = 'tv-th';
strt = 0;
mprintar(strx.thetatvexct, in, tit, strt);
disp(' ');
disp('press any key to continue')
pause


disp(' ')
disp('******** Computation of Recursive Residuals   ********');
%compute recursive residuals
%set up regression matrices
X = Y;
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
[Xt, Pt, g, M, initf, recrs, recr] = scakff(y, X, Z, G, W, T, H, ins, i);
%plot recursive residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 10;
ic = 1;
nr = strv.nparm;
disp(' ')
disp('******** Recursive Residuals:     ********');
str = mautcov(recr, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('Q statistics:')
disp(str.qstat)

[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all
