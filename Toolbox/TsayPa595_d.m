%script file for paragraph 5.9.5 in Tsay (2014)
%

clear
data = load(fullfile('data', 'm-bnd.dat'));
bnd = (data(:, 4:5));
[nx, mx] = size(bnd);
tdx = [1:nx] / 12 + 1954; %time index

subplot(2, 1, 1)
plot(tdx, bnd(:, 1))
xlabel('time');
ylabel('Aaa');
axis('tight');
subplot(2, 1, 2)
plot(tdx, bnd(:, 2))
xlabel('time');
ylabel('Baa');
axis('tight');
disp('press any key to continue')
pause
close all

disp(' ')
disp('identify a VAR model for the series')
disp(' ')


%identify a VAR model for the series
maxlag = 13;
minlag = 1;
prt = 1;
lagsopt = varident(bnd, maxlag, minlag, prt);
disp('press any key to continue')
pause

% %Matlab-Econ function for Johansen cointegration test
% lags=1:2;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);

[D, nr, yd, DA, ferror] = mcrcregr(bnd);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause

disp('number of unit roots : 1')
disp(' ')
x = [];
seas = 1;
%number of unit roots in the model
nr = 1;
prt = 0;
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(bnd, x, prt, nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Identification and estimation of a VAR model for the differenced series
% parameterized in error correction form
%
disp(' ')
disp('identify a VAR model for the differenced series')
disp(' ')
disp('press any key to continue')
pause

maxlag = 5;
minlag = 1;
prt = 1;
lagsopt = varident(yd, maxlag, minlag, prt);
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a VAR(2) model for the differenced series')
disp('using the Hannan-Rissanen method')
disp('')
disp('press any key to continue')
pause

%estimation of the model for the differenced series
%model for the differenced series: VAR(2).
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [.8, 1.];
pr = 2;
qr = 0; %VAR(2)
[str, ferror] = estvarmaxpqrPQR(yd, x, seas, [pr, qr, 0], [0, 0, 0], hr3, finv2, mstainv, nsig, tsig);

phi(:, :, 1) = str.phis3(:, :, 1);
for i = 1:pr
    phi(:, :, i+1) = str.phis3(:, :, i+1);
end
Phi(:, :, 1) = phi(:, :, 1);
th(:, :, 1) = str.thetas3(:, :, 1);
for i = 1:qr
    th(:, :, i+1) = str.thetas3(:, :, i+1);
end
Th(:, :, 1) = phi(:, :, 1);
Sigma = str.sigmar3;
clear str

phit = pmatmul(phi, Phi);
[Pi, Lambda, alpha, betap, ferror] = mid2mecf(phit, D, DA);

%set up varma model in error correction form. Parameters are defined in
%terms of the error correction form.
[str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, Sigma, seas);


disp(' ')
disp('estimation of the VAR(3) model in error correction form')
disp('press any key to continue')
pause
%estimate the model in error correction form
%
Y = [];
constant = 1;
result = varmapqPQestime(bnd, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(bnd, xvrf, xrf, str);
[Lambdaf, alphaf, betapf, thf, Thf, Lf, ferror] = pr2ecf(xvv, xff, DAf, str);
[phif, De, DAe, ferror] = mecf2mid(Lambdaf, alphaf, betapf);
phit = pmatmul(phif, Dr);

disp(' ');
disp('***** Estimated overall AR part  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(phit(:, :, 2:4), in, tit, strt);
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phif, 3) * result.h)', in, tit);
disp('press any key to continue')
pause

%error correction matrices
disp(' ');
disp('***** AR part in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2:end), in, tit, strt);
% reparameterize betap and alpha
[nb, mb] = size(betapf);
alphaf = alphaf * betapf(:, 1:nb);
betapf = betapf(:, 1:nb) \ betapf;
disp('estimated beta transposed')
disp(betapf)
disp('estimated alpha')
disp(alphaf)
disp('press any key to continue')
pause
clear str
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimation of the VAR(3) model parameterized in terms of the model for
% the differenced series
%
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [.8, 1.];
pr = 2;
qr = 0; %VAR(2)
[str, ferror] = estvarmaxpqrPQR(yd, x, seas, [pr, qr, 0], [0, 0, 0], hr3, finv2, mstainv, nsig, tsig);

phi(:, :, 1) = str.phis3(:, :, 1);
for i = 1:pr
    phi(:, :, i+1) = str.phis3(:, :, i+1);
end
Phi(:, :, 1) = phi(:, :, 1);
th(:, :, 1) = str.thetas3(:, :, 1);
for i = 1:qr
    th(:, :, i+1) = str.thetas3(:, :, i+1);
end
Th(:, :, 1) = phi(:, :, 1);
Sigma = str.sigmar3;
clear str


%set up varma model for the differenced series. Parameters are defined in
%terms of the model for the differenced series, contained in phi, th, Phi,
%Th and Sigma, and the parameters of the differencing matrix polynomial,
%contained in DA.
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, seas);

%eliminate insignificant parameters
for k = 2:pr + 1
    for i = 1:mx
        for j = 1:mx
            if phi(i, j, k) == 0
                str.phin(i, j, k) = 0;
            end
        end
    end
end
[str, ferror] = fixvarmapqPQ(str);


%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DA);

disp(' ')
disp('estimation of the  VAR(2) model for the differenced series. ')
disp('press any key to continue')
pause

%estimate model with cointegration relations imposed
Y = [];
constant = 1;
result = varmapqPQestimd(bnd, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(bnd, xvrf, xrf, str);
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvv, xff, str);
phit = pmatmul(phif, Dr);

%estimated VAR model with cointegration rank imposed
disp(' ');
disp('***** Estimated overall AR part  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(phit(:, :, 2:end), in, tit, strt);
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phif, 3) * result.h)', in, tit);
disp('press any key to continue')
pause

%error correction matrices
[Pi, Lambda, alpha, betafp, ferror] = mid2mecf(phif, Dr, DAf);
disp(' ');
disp('***** AR part in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2:3), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafp);
alpha = alpha * betafp(:, 1:nb);
betafp = betafp(:, 1:nb) \ betafp;
disp('estimated beta transposed')
disp(betafp)
disp('estimated alpha')
disp(alpha)
disp('press any key to continue')
pause

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Identification and estimation of a VARMA model in echelon form for the
% differenced series
%
x = [];
seas = 1;
%number of unit roots in the model
nr = 1;
prt = 0;
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(bnd, x, prt, nr);

%add information about the Kronecker indices
%estimate the Kronecker indices for the differenced series
maxorder = [];  %we fix the maximum order
prt = 0;
[order, kro, scm] = varmaxscmidn(yd, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the differenced series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause


disp(' ')
disp('estimate model in Echelon Form for the Differenced Series: ')
disp('Kronecker Indices are [1 1]  eliminating some insignificant ')
disp('paramters')
disp('press any key to continue')
pause

%estimate VARMA model in echelon form for the differenced series and
%aliminate some nonsignificant parameters
kro = [1, 1];
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [1., 1.];
strv = estvarmaxkro(yd, x, seas, kro, hr3, finv2, mstainv, nsig, tsig);


%add unit root information to the model structure
[strv, ferror] = aurirvarmapqPQ(strv, nr, DA);
%add freq to the model structure
strv.freq = seas;
%estimate model
s = size(bnd, 2);
Y = [];
constant = 1;
[xvfk, strx, ferror] = mexactestimcd(bnd, strv, Y, constant);


%create estimated model
xffk = [];
[ydd, xvvk, xffk, DAfk, Drk, Dsk, ferrork] = pr2varmapqPQd(bnd, xvfk, xffk, strx);
phifk = strx.phisexct;
thfk = strx.thetasexct;
[nphi, mphi, prk] = size(phifk);

%error correction matrices
[Pik, Lambdak, alphak, betafpk, ferror] = mid2mecf(phifk, Drk, DAfk);
disp(' ');
disp('***** Estimated model in error correction form  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdak(:, :, 2), in, tit, strt);
tit = 'Theta';
strt = 1;
mprintar(thfk(:, :, 2), in, tit, strt);
disp('t-values of Phi and Theta matrices:')
tit = 'tv-Phi';
strt = 1;
mprintar(strx.phitvexct(:, :, 2), in, tit, strt);
tit = 'tv-Theta';
strt = 1;
mprintar(strx.thetatvexct(:, :, 2), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafpk);
alphak = alphak * betafpk(:, 1:nb);
betafpk = betafpk(:, 1:nb) \ betafpk;
disp('estimated beta transposed')
disp(betafpk)
disp('estimated alpha')
disp(alphak)
%constant is phifk(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.7f');
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);
tit = 'Constant:';
mprintar((sum(phifk, 3) * strx.musexct)', in, tit);
% in.cnames = char(' Sigma:',' ',' ',' ',' Constant:');
% mprint([strx.sigmarexct sum(phifk,3)*strx.musexct],in);

disp(' ')
disp('******** Computation of ols Residuals   ********');
clear str

m = 0;
[strc, ferror] = matechelon(kro, s, m);
strc.phis = phifk;
strc.thetas = thfk;
strc.gammas = [];
strc = armaxe2sse(strc);
strc.sigmar2 = strx.sigmarexct;
% %residuals with fixed regression parameters and square root covariance
% %filter using the fast CKMS recursions
% Y=eye(s); tol=1.d-10; maxupdt=[];
% % mkro=max(kro); Sigma=strx.sigmarexct;
% %  [c,ierror]=macgf(phifk,thfk,Sigma,mkro+1);
%  [e,E,rSigmat]=sqrt_ckms(ydd,Y,strc,maxupdt,tol);
%  [ne,me]=size(e); recr0=zeros(ne,me); nbeta=s;
% for ii=1:ne
%  ind=(ii-1)*nbeta+1:ii*nbeta;     % V=rSigmat(ind,:);
%  recr0(ii,:)=e(ii,:) - (E(ind,:)'*strx.musexct)';
% end

%compute OLS residuals
%set up regression matrices
X = eye(s);
W = [];
Sigmax = strc.sigmar2;
[L, p] = chol(Sigmax, 'lower');
Lm1 = pinv(L);
%set up system matrices
T = strc.Fs;
Z = strc.Hs;
G = Lm1;
H = strc.Ks * Lm1;
%set up initial conditions
ndelta = 0; %number of unit roots
[ins, i, ferror] = incossm(T, H, ndelta);
[Xt, Pt, g, M, initf, recrs1, recr1] = scakffsqrt(ydd, X, Z, G, W, T, H, ins, i);
%residuals with fixed regression parameters and square root covariance
%filter
[KKP, PT, recrs, recr] = scakfffsqrt(ydd, X, Z, G, W, T, H, ins, i, strx.musexct);

%plot OLS residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
disp('press any key to continue')
pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 15;
ic = 1;
nr = 0; % nr=strx.nparm+3;
disp(' ')
disp('******** OLS Residuals:     ********');
stre = mautcov(recr, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(stre.r0)
disp('Q statistics:')
disp(stre.qstat)

disp('p-values of Q statistics:')
disp(stre.pval)
[m, n] = size(stre.pval);
t = 1:m;
plot(t, stre.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all
