%script file example in Paragraph 4.7.2 in Tsay (2014)
%
clear
data = load(fullfile('data', 'q-gdp-ukcaus.dat'));
gdp = log(data(:, 3:5));
zt = diferm(gdp, 1);
[nx, mx] = size(zt);
tdx = [1:nx] / 4 + 1980; %time index

subplot(3, 1, 1)
plot(tdx, zt(:, 1))
xlabel('year');
ylabel('UK');
axis('tight');
subplot(3, 1, 2)
plot(tdx, zt(:, 2))
xlabel('year');
ylabel('CA');
axis('tight');
subplot(3, 1, 3)
plot(tdx, zt(:, 3))
xlabel('year');
ylabel('US');
axis('tight');
disp('press any key to continue')
pause
close all

disp(' ')
disp('identify Kronecker indices')
disp(' ')
disp('press any key to continue')
pause


%estimate the Kronecker indices for the series
maxorder = [];
hr3 = 0;
x = [];
seas = 1;
prt = 0;
[order, kro, scm] = varmaxscmidn(zt, x, seas, maxorder, hr3, prt);
disp('estimated Kronecker Indices for the differenced series ')
disp('using function "varmaxscmidnt":')
disp(kro)
disp('press any key to continue')
pause

%estimate model using HR method and eliminate
%nonsignificant parameters
disp('estimate model with kronecker indices [1 1 0] ')
disp('and eliminate some nonsignificant parameters')
kro = [1, 1, 0];
hr3 = 0;
finv2 = 1;
nsig = [1, 1];
tsig = [1., 1.];
mstainv = 0;
strv = estvarmaxkro(zt, x, seas, kro, hr3, finv2, mstainv, nsig, tsig);
disp(' ');
disp('***** Estimated  simplified VARMA Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');

tit = 'phi';
strt = 0;
mprintar(strv.phis3(:, :, 1:2), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strv.thetas3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);

disp(' ');
disp('*****t-values  *****');
disp(' ');
tit = 'tv-phi';
strt = 0;
mprintar(strv.phitv3(:, :, 1:2), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strv.thetatv3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma (x 10^4)';
mprintar(strv.sigmar3*1d5, in, tit);
disp('press any key to continue')
pause


%estimate using the exact method
disp(' ')
disp('estimation using the exact method')
disp('Kronecker indices are [1 1 0]')
disp('press any key to continue')
pause
[ny, s] = size(zt);
Y = eye(s);
xx = [];
[xvf, strx, ferrorc] = mexactestim(zt, xx, strv, Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 0;
mprintar(strx.phisexct(:, :, 1:2), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strx.thetasexct(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Mean';
mprintar(strx.musexct', in, tit);
disp(' ')
tit = 'Sigma (x 10^4)';
mprintar(strx.sigmarexct*1d5, in, tit);
disp(' ')
disp(' ')
disp('t-values: ')
tit = 'tv-phi';
strt = 0;
mprintar(strx.phitvexct(:, :, 1:2), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strx.thetatvexct(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Mean';
mprintar(strx.mutvexct', in, tit);
disp('press any key to continue')
pause


disp(' ')
disp('******** Computation of OLS Residuals   ********');
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
[Xt, Pt, recrs, recr] = scakfff(zt, X, Z, G, W, T, H, ins, i, strx.musexct); %regression parameters fixed
%plot  residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
plot(recr(:, 3)), legend('recr(:,3)'), pause
close all
%compute autocovariance and autocorrelation matrices ofresiduals
lag = 24;
ic = 1;
nr = 0;
disp(' ')
disp('******** OLS Residuals:     ********');
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
disp('press any key to continue')
pause
close all
