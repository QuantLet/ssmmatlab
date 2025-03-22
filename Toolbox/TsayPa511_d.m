%script file for paragraph 5.11 in Tsay (2014)
%

clear
data = load(fullfile('data', 'q-4macro.dat'));
da = (data(:, 4:7));
zt = [log(da(:, 1)), da(:, 2), log(da(:, 4)), da(:, 3)];
[nx, mx] = size(zt)
tdx = [1:nx] / 4 + 1959; %time index

subplot(2, 2, 1)
plot(tdx, zt(:, 1))
xlabel('time');
ylabel('ln(gnp)');
axis('tight');
subplot(2, 2, 2)
plot(tdx, zt(:, 3))
xlabel('time');
ylabel('lnm1');
axis('tight');
subplot(2, 2, 3)
plot(tdx, zt(:, 2))
xlabel('time');
ylabel('tb3m');
axis('tight');
subplot(2, 2, 4)
plot(tdx, zt(:, 4))
xlabel('time');
ylabel('gs10');
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
lagsopt = varident(zt, maxlag, minlag, prt);
disp('press any key to continue')
pause


% %Matlab-Econ function for Johansen cointegration test
% lags=1:2;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);

[D, nr, yd, DA, ferror] = mcrcregr(zt);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause

disp('number of unit roots according to Tsay (2014): 2')
disp('this will be used in the following')
disp(' ')
x = [];
seas = 1;
%number of unit roots in the model
nr = 2;
prt = 0;
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(zt, x, prt, nr);

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

maxlag = 8;
minlag = 1;
prt = 1;
lagsopt = varident(yd, maxlag, minlag, prt);
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a VAR(5) model for the differenced series')
disp('using the Hannan-Rissanen method')
disp('')
disp('press any key to continue')
pause

%estimation of the model for the differenced series
%model for the differenced series: VAR(5).
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [.8, 1.];
pr = 5;
qr = 0; %VAR(5)
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
disp('estimation of the VAR(6) model in error correction form')
disp('press any key to continue')
pause
%estimate the model in error correction form
%
Y = [];
constant = 1;
result = varmapqPQestime(zt, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(zt, xvrf, xrf, str);
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
mprintar(phit(:, :, 2:end-1), in, tit, strt);
disp(' ')
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
% Estimation of the VAR(6) model parameterized in terms of the model for
% the differenced series
%
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [.8, 1.];
pr = 5;
qr = 0; %VAR(5)
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
disp('estimation of the  VAR(5) model for the differenced series. ')
disp('press any key to continue')
pause

%estimate model with cointegration relations imposed
Y = [];
constant = 1;
result = varmapqPQestimd(zt, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(zt, xvrf, xrf, str);
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvv, xff, str);
phit = pmatmul(phif, Dr);

%estimated VAR model with cointegration rank imposed
disp(' ');
disp('***** Estimated overall AR part  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
% in.cnames = char(' AR(1):',' ',' ',' ',' AR(2):',' ',' ',' ',' AR(3):',' ',' ',' ',...
%                                ' AR(4):',' ',' ',' ',' AR(5):',' ',' ',' ',' AR(6):',' ',' ',' ');
% mprint(phit(:,:,2:end),in);
tit = 'AR';
strt = 1;
mprintar(phit(:, :, 2:7), in, tit, strt);
disp(' ')
clear in
in.fid = 1;
in.fmt = char('%12.7f');
% in.cnames = char(' Sigma:',' ',' ',' ',' Constant:');
% mprint([result.Sigmar sum(phif,3)*result.h],in);
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
% in.cnames = char(' Lambda(1):',' ',' ',' ',' Lambda(2):',' ',' ',' ',...
%                                ' Lambda(3):',' ',' ',' ',' Lambda(4):',' ',' ',' ',...
%                               ' Lambda(5):',' ',' ',' ' );
% mprint(Lambda(:,:,2:6),in);
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2:6), in, tit, strt);
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


disp(' ')
disp('******** Computation of OLS Residuals   ********');
clear str

[str, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, result.Sigmar, seas);
s = mx;
kro = pr * ones(1, s);
m = 0;
[strc, ferror] = matechelon(kro, s, m);
thh(:, :, pr+1) = zeros(s);
thh(:, :, 1) = eye(s);
strc.phis = phif;
strc.thetas = thh;
strc.gammas = [];
strc = armaxe2sse(strc);
strc.sigmar2 = str.Sigma;
Y = eye(s);
tol = 1.d-10;
maxupdt = [];
[e, E, rSigmat] = sqrt_ckms(ydf, Y, strc, maxupdt, tol);
[ne, me] = size(e);
recr = zeros(ne, me);
nbeta = s;
for ii = 1:ne
    ind = (ii - 1) * nbeta + 1:ii * nbeta; % V=rSigmat(ind,:);
    recr(ii, :) = e(ii, :) - (E(ind, :)' * result.h)';
end

% %compute recursive residuals
% %set up regression matrices
% X=eye(s); W=[];
% Sigmax=strc.sigmar2;
% [L,p] = chol(Sigmax,'lower'); Lm1=pinv(L);
% %set up system matrices
% T=strc.Fs; Z=strc.Hs; G=Lm1; H=strc.Ks*Lm1;
% %set up initial conditions
% ndelta=0;                            %number of unit roots
% [ins,i,ferror]=incossm(T,H,ndelta);
% [Xt,Pt,g,M,initf,recrs,recr]=scakffsqrt(ydf,X,Z,G,W,T,H,ins,i);
%plot recursive residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
plot(recr(:, 3)), legend('recr(:,3)'), pause
plot(recr(:, 4)), legend('recr(:,4)'), pause
disp('press any key to continue')
pause
close all
%compute autocovariance and autocorrelation matrices of ols residuals
lag = 24;
ic = 1;
nr = 0;
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
pause
close all
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Identification and estimation of a VARMA model in echelon form for the
% differenced series
%
disp('the absolute values of the eigenvalues of ')
disp('the transition matrix in the former VAR')
disp('model are:')
abs(eig(str.T))
pause
disp('there is one eigenvalue close to one. This')
disp('suggests that there is only one cointegration')
disp('relationaship. In fact, ')
[D, nr, yd, DA, ferror] = mcrcregr(zt);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause

x = [];
seas = 1;
%number of unit roots in the model
nr = 3;
prt = 0;
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(zt, x, prt, nr);

%add information about the Kronecker indices
%estimate the Kronecker indices for the differenced series
maxorder = [];  
prt = 0;
[order, kro, scm] = varmaxscmidn(yd, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the differenced series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

%identify a VARMA(p,q) model for the differenced series
disp('identify a VARMA(p,q) model for the differenced series')
disp(' ')
maxlag = [];  
minlag = 0;
prt = 0;
[lagsopt, ferror] = lratiopqr(yd, x, seas, maxlag, minlag, prt);
disp('estimated orders in VARMA(p,q,r):')
disp(lagsopt)
disp('press any key to continue')
pause


disp(' ')
disp('estimate model in Echelon Form for the Differenced Series: ')
disp('Kronecker Indices are [1 1 1 0] ')
disp('press any key to continue')
pause

%estimate VARMA model in echelon form for the differenced series 
kro = [1, 1, 1, 0];
hr3 = 0;
finv2 = 1;
mstainv = 1;
strv = estvarmaxkro(yd, x, seas, kro, hr3, finv2, mstainv);


%add unit root information to the model structure
[strv, ferror] = aurirvarmapqPQ(strv, nr, DA);
%add freq to the model structure
strv.freq = seas;
%estimate model
s = size(zt, 2);
Y = [];
constant = 1;
[xvfk, strx, ferror] = mexactestimcd(zt, strv, Y, constant);


%create estimated model
xffk = [];
[ydd, xvvk, xffk, DAfk, Drk, Dsk, ferrork] = pr2varmapqPQd(zt, xvfk, xffk, strx);
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
% in.cnames = char(' Lambda(1):',' ',' ',' ',' Lambda(2):',' ',' ',' ');
% mprint(Lambdak(:,:,2:3),in);
tit = 'Lambda';
strt = 1;
mprintar(Lambdak(:, :, 2), in, tit, strt);
% in.cnames = char(' Theta(1):',' ',' ',' ',' Theta(2):',' ',' ',' ');
% mprint(thfk(:,:,2:3),in);
tit = 'Theta';
strt = 1;
mprintar(thfk(:, :, 2), in, tit, strt);
disp('t-values of Phi and Theta matrices:')
% in.cnames = char(' tv-Phi(1):',' ',' ',' ',' tv-Phi(2):',' ',' ',' ');
% mprint(strx.phitvexct(:,:,2:3),in);
tit = 'tv-Phi';
strt = 1;
mprintar(strx.phitvexct(:, :, 2), in, tit, strt);
% in.cnames = char(' tv-Theta(1):',' ',' ',' ',' tv-Theta(2):',' ',' ',' ');
% mprint(strx.thetatvexct(:,:,2:3),in);
tit = 'tv-theta';
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
% in.cnames = char(' Sigma:',' ',' ',' ',' Constant:');
% mprint([strx.sigmarexct sum(phifk,3)*strx.musexct],in);
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);
tit = 'Constant:';
mprintar((sum(phifk, 3) * strx.musexct)', in, tit);


disp(' ')
disp('******** Computation of OLS Residuals   ********');
clear str

m = 0;
[strc, ferror] = matechelon(kro, s, m);
strc.phis = phifk;
strc.thetas = thfk;
strc.gammas = [];
strc = armaxe2sse(strc);
strc.sigmar2 = strx.sigmarexct;
Y = eye(s);
tol = 1.d-10;
maxupdt = [];
[e, E, rSigmat] = sqrt_ckms(ydd, Y, strc, maxupdt, tol);
[ne, me] = size(e);
recrs = zeros(ne, me);
nbeta = s;
for ii = 1:ne
    ind = (ii - 1) * nbeta + 1:ii * nbeta; % V=rSigmat(ind,:);
    recrs(ii, :) = e(ii, :) - (E(ind, :)' * strx.musexct)';
end
%plot residuals
plot(recrs(:, 1)), legend('recrs(:,1)'), pause
plot(recrs(:, 2)), legend('recrs(:,2)'), pause
plot(recrs(:, 3)), legend('recrs(:,3)'), pause
plot(recrs(:, 4)), legend('recrs(:,4)'), pause
disp('press any key to continue')
pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 24;
ic = 1;
nr = 0; % nr=strx.nparm+3;
disp(' ')
disp('********  Residuals:     ********');
stre = mautcov(recrs, lag, ic, nr);
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
