%script file for Paragraph 3.15 in Tsay (2014), housing starts and mortgage
%rate
%

da = load(fullfile('data', 'm-hsmort7112.dat'));
zt = da(:, 3:4);
dzt = diferm(zt, 1);
tdx = da(:, 1) + da(:, 2) / 12;
dzt(:, 1) = dzt(:, 1) / 1000;


[nda, mda] = size(da);
subplot(2, 1, 1)
plot(tdx(2:nda), dzt(:, 1))
subplot(2, 1, 2)
plot(tdx(2:nda), dzt(:, 2))
pause
close all


disp(' ');
disp('***** VAR order identification  *****');
disp(' ');

%VAR order identification
prt = 1;
minlag = 0;
maxlag = 13;
lagsopt = varident(dzt, maxlag, minlag, prt);
pause


%estimate VAR(4)
nlag = 4;
res = var_est(dzt, nlag);

disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2:5), in, tit, strt);
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
mprintar(res.phitv(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);


nvar = 2;
lag = 24;
ic = 1;
nr = 0;
disp('******** VAR residuals:     ********');
str = mautcov(res.resid, lag, ic, nr);
disp('Q statistics:')
disp(str.qstat)
disp('press any key to continue')
pause
disp('p-values of Q statistics:')
disp(str.pval)
disp('press any key to continue')
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all


disp(' ')
disp('Elimination of nonsignificant parameters in the model ')
disp('')
disp('press any key to continue')
pause


%eliminate some nonsignificant  parameters in the model
nsig = [1, 0];
tsig = [1., 0.];
mstainv = 0;
hr3 = 0;
finv2 = 1;
xx = [];
freq = 1;
[strvr, ferror] = estvarmaxpqrPQR(dzt, xx, freq, [4, 0, 0], [0, 0, 0], hr3, finv2, ...
    mstainv, nsig, tsig);

disp(' ');
disp('***** Estimated simplified VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(strvr.phis3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvr.mus3', in, tit);
disp(' ');
tit = 'Sigma';
mprintar(strvr.sigmar3, in, tit);

disp('*****t-values  *****');
disp(' ');
tit = 'tv-AR';
strt = 1;
mprintar(strvr.phitv3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvr.mutv3', in, tit);
disp(' ')
disp('press any key to continue')
pause


disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause

%identify a VARMA(p,q) model for the series
maxlag = -2; %we fix the maximum order
minlag = 0;
prt = 0;
x = [];
seas = 1;
[lagsopt, ferror] = lratiopqr(dzt, x, seas, maxlag, minlag, prt);
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
prt = 0;
maxorder = [];
hr3 = 0;
[order, kro, scm] = varmaxscmidn(dzt, x, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the original series ')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause


disp(' ')
disp('estimate simplified VARMA(2,1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause

freq = 1;
xx = [];
hr3 = 1;
finv2 = 0;
nsig = [1, 0];
tsig = [.8, 0];
mstainv = 0;
[strvr, ferror] = estvarmaxpqrPQR(dzt, xx, freq, [2, 1, 0], [0, 0, 0], hr3, finv2, ...
    mstainv, nsig, tsig);

disp(' ');
disp('***** Estimated VARMA Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strvr.phis(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strvr.thetas(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvr.mu', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(strvr.phitv(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strvr.thetatv(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvr.mutv', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strvr.sigmar2, in, tit);
disp('press any key to continue')
pause


%estimate the model using the conditional method
hr3 = 0;
finv2 = 1;
seas = 1;
strvr = mhanris(dzt, xx, seas, strvr, hr3, finv2);
[xvfc, strc, ferrorc] = mconestim(dzt, xx, strvr);

disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 2:3), in, tit, strt);
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
mprintar(strc.phitvcon(:, :, 2:3), in, tit, strt);
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
%nr is the number of estimated parameters
nr = strc.nparm;
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

disp(' ')
disp('estimate simplified VARMA(2,1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause

freq = 1;
xx = [];
hr3 = 1;
finv2 = 0;
nsig = [1, 0];
tsig = [.8, 0];
mstainv = 0;
[strvr, ferror] = estvarmaxpqrPQR(dzt, xx, freq, [2, 1, 0], [0, 0, 0], hr3, finv2, ...
    mstainv, nsig, tsig);


disp(' ');
disp('***** Estimated VARMA Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strvr.phis(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strvr.thetas(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvr.mu', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(strvr.phitv(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strvr.thetatv(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvr.mutv', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strvr.sigmar2, in, tit);
disp('press any key to continue')
pause


%estimate the model using the conditional method
%first, compute initial conditions
hr3 = 0;
finv2 = 1;
seas = 1;
strvr = mhanris(dzt, xx, seas, strvr, hr3, finv2);
%then, apply contional method
[xvfc, strc, ferrorc] = mconestim(dzt, xx, strvr);

disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 2:3), in, tit, strt);
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
mprintar(strc.phitvcon(:, :, 2:3), in, tit, strt);
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
nr = 0;
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
