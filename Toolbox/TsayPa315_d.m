%script file for Paragraph 3.15 in Tsay (2014)
%

da1 = load(fullfile('data', 'm-pce.dat'));
da2 = load(fullfile('data', 'm-dspi.dat'));
y = [da1(:, 4), da2(:, 4)];
y = log(y);
zt = diferm(y, 1) * 100;
tdx = da1(:, 1) + da1(:, 2) / 12;
subplot(2, 1, 1)
plot(tdx(2:639), zt(:, 1))
legend('Pceg');
axis('tight');
subplot(2, 1, 2)
plot(tdx(2:639), zt(:, 2))
legend('Dspig');
axis('tight');
pause
close all


%VAR order identification
disp(' ');
disp('***** VAR order identification  *****');
disp(' ');
prt = 1;
minlag = 0;
maxlag = 13;
lagsopt = varident(zt, maxlag, minlag, prt);
disp('press any key to continue')
pause


%estimate VAR(3)
disp(' ')
disp('estimation of a VAR(3) model for the series')
disp('')
disp('press any key to continue')
pause
nlag = 3;
res = var_est(zt, nlag);


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


%eliminate some nonsignificant  parameters in the model using function
%estvarmaxpqrPQR.
disp(' ')
disp('Elimination of nonsignificant parameters in the VAR(3) model ')
disp('')
disp('press any key to continue')
pause
freq = 1;
xx = [];
hr3 = 1;
finv2 = 0;
nsig = [1, 0];
tsig = [1., 0.];
mstainv = 0;
[strvr, ferror] = estvarmaxpqrPQR(zt, xx, freq, [3, 0, 0], [0, 0, 0], hr3, finv2, ...
    mstainv, nsig, tsig);
disp(' ');
disp('***** Estimated simplified VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(strvr.phis(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvr.mu', in, tit);
disp(' ');
tit = 'Sigma';
mprintar(strvr.sigmar2, in, tit);

disp('*****t-values  *****');
disp(' ');
tit = 'tv-AR';
strt = 1;
mprintar(strvr.phitv(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvr.mutv', in, tit);
disp(' ')
disp('press any key to continue')
pause


%model checking
nvar = 2;
lag = 24;
ic = 1;
%in Tsay (2014), p. 182, the number of paramters is zero. Different results
%are obtained with nr=strvr.nparm, the number of parameters in the model.
nr = 0;
disp('******** VAR residuals:     ********');
str = mautcov(strvr.resid2, lag, ic, nr);
disp('Q statistics:')
disp(str.qstat)
disp('press any key to continue')
pause
disp('p-values of Q statistics:')
disp(str.pval)
disp('press any key to continue')
pause
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all


%impulse response functions
nc = 12;
[K, ierror] = ptransfer(strvr.phis, strvr.thetas, nc);
mk = -.25;
Mk = 1.1;
subplot(2, 2, 1)
plot(squeeze(K(1, 1, :)))
title('Orig. innovations');
axis([0, nc + .5, mk, Mk]);
xlabel('Lag');
ylabel('IRF');
subplot(2, 2, 2)
plot(squeeze(K(1, 2, :)))
title('Orig. innovations');
axis([0, nc + .5, mk, Mk]);
xlabel('Lag');
ylabel('IRF');
subplot(2, 2, 3)
plot(squeeze(K(2, 1, :)))
title('Orig. innovations');
axis([0, nc + .5, mk, Mk]);
xlabel('Lag');
ylabel('IRF');
subplot(2, 2, 4)
plot(squeeze(K(2, 2, :)))
title('Orig. innovations');
axis([0, nc + .5, mk, Mk]);
xlabel('Lag');
ylabel('IRF');
disp('press any key to continue')
pause
close all

disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause

%identify a VARMA(p,q) model for the series
maxlag = -3;  %we fix the maximum order
minlag = 0;
prt = 0;
x = [];
seas = 1;
[lagsopt, ferror] = lratiopqr(zt, x, seas, maxlag, minlag, prt);
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
maxorder = -3; %we fix the maximum order
hr3 = 0;
[order, kro, scm] = varmaxscmidn(zt, x, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the original series ')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

disp(' ')
disp('estimate VARMA(3,1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause

freq = 1;
xx = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(zt, xx, freq, [3, 1, 0], [0, 0, 0], hr3, finv2);

%eliminate some nonsignificant  parameters in the model
nsig = [0, 1];
tsig = [0., 1.];
mstainv = 0;
[strvr, ferror] = estvarmaxpqrPQR(zt, xx, freq, [3, 1, 0], [0, 0, 0], hr3, finv2, ...
    mstainv, nsig, tsig);

disp(' ');
disp('***** Estimated VARMA Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strvr.phis3(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strvr.thetas3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvr.mus3', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(strvr.phitv3(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strvr.thetatv3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvr.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strvr.sigmar3, in, tit);
disp('press any key to continue')
pause


%estimate the model using the conditional method
[xvfc, strc, ferrorc] = mconestim(zt, xx, strvr);


disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strc.thetascon(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Mean';
mprintar(strc.muscon', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strc.sigmarcon, in, tit);
disp(' ')
disp(' ')
disp('t-values: ')
tit = 'tv-phi';
strt = 1;
mprintar(strc.phitvcon(:, :, 2:4), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strc.thetatvcon(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Mean';
mprintar(strc.mutvcon', in, tit);
disp(' ')
disp('press any key to continue')
pause


lag = 24;
ic = 1;
nr = 4;
disp(' ')
disp('******** Conditional Residuals:     ********');
str = mautcov(strc.residcon, lag, ic, nr);
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
