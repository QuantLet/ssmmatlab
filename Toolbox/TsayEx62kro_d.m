%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example  6.2  in Tsay (2014), pp. 346-352
%
% We identify Kronecker indices for a VARMAX model in echelon form
%
%Series are:  1) monthly U.S. regular conventional gas price z1t and
% 2)  heating oil price z2t of New York Harbor. Both series are measured in dollars
% per gallon. These prices depend on the crude oil and natural gas prices;  3) x1t
% the spot oil price of West Texas Intermediate, dollars per barrel, and  4) x2t the
% natural gas price of Henry Hub, LA, measured in dollars per million BTU. Thus,
% yt =(z1t, z2t), xt = (x1t, x2t), and k = m = 2. The sample period is from
% November 1993 to August 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

da = load(fullfile('data', 'm-gasoil.dat'));
yt = da(:, 3:6);
tdx = da(:, 1) + da(:, 2) / 12;


subplot(2, 2, 1)
plot(tdx, yt(:, 1))
xlabel('time');
ylabel('Greg');
axis('tight');
subplot(2, 2, 2)
plot(tdx, yt(:, 3))
xlabel('time');
ylabel('oilp');
axis('tight');
subplot(2, 2, 3)
plot(tdx, yt(:, 2))
xlabel('time');
ylabel('hoil');
axis('tight');
subplot(2, 2, 4)
plot(tdx, yt(:, 4))
xlabel('time');
ylabel('gasp');
axis('tight');
disp('press any key to continue')
pause
close all

zt = yt(:, 1:2);
xt = yt(:, 3:4);

disp(' ')
disp('estimate Kronecker indices for the series')
disp('press any key to continue')
pause
%estimate the Kronecker indices for the original series
prt = 1;
seas = 1;
maxorder = -2; %we fix the maximum order to two
hr3 = 0;
[order, kro, scm] = varmaxscmidn(zt, xt, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the original series ')
disp('using function "varmaxscmidn":')
disp(scm)
disp(kro)
disp('press any key to continue')
pause

disp('estimate model in Echelon Form using the Hanna-Rissanen method: ')
disp('press any key to continue')
pause


%estimate model using HR method (K.i. = [2 1]) and eliminate some
%nonsignificant parameters
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [1., 1.];
strv = estvarmaxkro(zt, xt, seas, kro, hr3, finv2, mstainv, nsig, tsig);
mlag = max(kro) + 1; %length of the matrix polynomials

disp(' ');
disp('***** Estimated VARMAX Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(strv.phis3(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'X part';
strt = 0;
mprintar(strv.gammas3(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'MA';
strt = 1;
mprintar(strv.thetas3(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ')
tit = 'tv-AR';
strt = 1;
mprintar(strv.phitv3(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-X part';
strt = 0;
mprintar(strv.gammatv3(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-MA';
strt = 1;
mprintar(strv.thetatv3(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause


%compute autocovariance and autocorrelation matrices of residuals
lag = 24;
ic = 1;
nr = 0; %nr=strx.nparm;
disp(' ')
disp('******** Residuals:     ********');
str = mautcov(strv.resid3, lag, ic, nr);
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

disp(' ')
disp('estimate using the conditional method')
disp('press any key to continue')
pause

%estimate using the conditional method
[xvfc, strc, ferrorc] = mconestim(zt, xt, strv);

disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strc.thetascon(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'X-part';
strt = 0;
mprintar(strc.gammascon(:, :, 1:mlag), in, tit, strt);
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
mprintar(strc.phitvcon(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strc.thetatvcon(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-X-part';
strt = 0;
mprintar(strc.gammatvcon(:, :, 1:mlag), in, tit, strt);
disp(' ');
tit = 'tv-Mean';
mprintar(strc.mutvcon', in, tit);
disp(' ')
disp('press any key to continue')
pause

disp(' ')
disp('estimate using the exact method')
disp('press any key to continue')
pause

%estimate model using the exact method
Y = eye(2);
[xvfx, strx, ferror] = mexactestimc(zt, xt, strc, Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strx.phisexct(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strx.thetasexct(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'X-part';
strt = 0;
mprintar(strx.gammasexct(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'Mean';
mprintar(strx.musexct', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strx.sigmarexct, in, tit);
disp(' ')
disp(' ')
disp('t-values: ')
tit = 'tv-phi';
strt = 1;
mprintar(strx.phitvexct(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strx.thetatvexct(:, :, 1:mlag), in, tit, strt);
disp(' ')
tit = 'tv-X-part';
strt = 0;
mprintar(strx.gammatvexct(:, :, 1:mlag), in, tit, strt);
disp(' ');
tit = 'tv-Mean';
mprintar(strx.mutvexct', in, tit);
