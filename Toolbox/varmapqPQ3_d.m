%Example of estimation of the flour price series using a
%VARMA model in echelon form for the differenced series. 
%
% Monthly Flour Price Indices for Three U.S. cities, Buffalo, Minneapolis,
% and Kansas City, for the Period August 1972 Through November 1980. The 
% three series are I(1) and it seems that there are no cointegration
% relationships. These series have been used by Tiao y Tsay (1989) and
% Lütkepohl and Poskitt (1996).

%load data
y=load(fullfile('data','flour-price.dat')); x=[];
%logs are taken
y=log(y); 
seas=1;

% lag=20; cw=1.96; ds=0; dr=1;
% tname={'flour price 1','flour price 2','flour price 3'};
% for i=1:3      
%   c0=sacspacdif(y(:,i),tname(i),dr,ds,seas,lag,cw);
%   pause
% end
% close all

%identify a VAR model for the series
maxlag=6; minlag=1; prt=1;
lagsopt = varident(y,maxlag,minlag,prt);
disp('press any key to continue')
pause

%identify a VARMA(p,q) model for the series
maxlag=6; minlag=0; prt=0;
[lagsopt,ferror] = lratiopqr(y,x,seas,maxlag,minlag,prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause



% %Matlab-Econ function for Johansen cointegration test
% lags=lagsopt;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);

[D,nr,yd,DA,ferror]=mcrcregr(y);
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause


%series is differenced
dr=1;
yd=diferm(y,dr);


%estimate the Kronecker indices for the differenced series
maxorder = []; 
hr3 = 0; 
prt = 0;
[order,kro,scm] = varmaxscmidn(yd,x,seas,maxorder,hr3,prt);
disp('estimated Kronecker Indices for the differenced series ')
disp('using function "varmaxscmidnt":')
disp(kro)
disp('press any key to continue')
pause

disp('estimate model in Echelon Form for the differenced series ')
disp('using the Hannan-Rissanen method and eliminate some ')
disp('insignificant paramters')                         
disp('We change the Kronecker Indices to [1 0 0].               ')
disp('press any key to continue')
pause



%estimate model using HR method (K.i. = [1 0 0]) and eliminate some
%nonsignificant parameters
kro=[1 0 0];
hr3=0; finv2=1; mstainv=1;  nsig=[1 0]; tsig=[1. 0.];
strv = estvarmaxkro(yd,x,seas,kro,hr3,finv2,mstainv,nsig,tsig); 


disp(' ');
disp('***** Estimated Model using the HR method  *****');
disp(' ');
clear in
in.fid=1;
in.fmt = char('%12.4f'); 
tit = 'phi'; strt = 0;
mprintar(strv.phis3,in,tit,strt);
disp(' ')
tit = 'th'; strt = 0;
mprintar(strv.thetas3,in,tit,strt);
disp(' ')
disp(' ');
tit = 'Sigma:';  
mprintar(strv.sigmar3,in,tit);
disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f'); 
tit = 'tv-phi'; strt = 0;
mprintar(strv.phitv3,in,tit,strt);
disp(' ');
tit = 'tv-th'; strt = 0;
mprintar(strv.thetatv3,in,tit,strt);
disp('press any key to continue')
pause






%estimate using the conditional method
[xvfc,strc,ferrorc]=mconestim(yd,x,strv);


disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid=1;
in.fmt = char('%12.4f'); 
tit = 'phi';  strt = 0;
mprintar(strc.phiscon(:,:,1),in,tit,strt);
tit = 'th';  strt = 0;
mprintar(strc.thetascon(:,:,1:2),in,tit,strt);
disp(' ');
tit = 'Sigma:'; 
mprintar(strc.sigmarcon,in,tit);
disp(' ')
disp('t-values: ')
disp(' ');
tit = 'tv-mu:'; 
mprintar(strc.mutvcon',in,tit);
tit = 'tv-phi';  strt = 0;
mprintar(strc.phitvcon(:,:,1),in,tit,strt);
tit = 'tv-th';  strt = 0;
mprintar(strc.thetatvcon(:,:,1:2),in,tit,strt);
disp('press any key to continue')
pause


%compute autocovariance and autocorrelation matrices of residuals
lag=6; ic=1;
disp('******** Conditional Residuals:     ********');
str=mautcov(strc.residcon,lag,ic);  
disp('press any key to continue')
pause

disp('estimation using the exact method')
disp('press any key to continue')
pause

%estimate using the exact method
Y=[];                                         
[xvf,strx,ferror]=mexactestim(yd,x,strv,Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid=1;
in.fmt = char('%12.4f'); 
tit = 'phi'; strt = 0;
mprintar(strx.phisexct(:,:,1),in,tit,strt);
disp(' ')
tit = 'th'; strt = 0;
mprintar(strx.thetasexct,in,tit,strt);
disp(' ')
tit = 'Sigma'; 
mprintar(strx.sigmarexct,in,tit);
disp(' ')
tit = 'Lh'; 
mprintar(strx.Lh',in,tit);
disp(' ')
disp('t-values: ')
tit = 'tv-phi'; strt = 1;
mprintar(strx.phitvexct(:,:,1),in,tit,strt);
disp(' ')
tit = 'tv-th'; strt = 1;
mprintar(strx.thetatvexct,in,tit,strt);
disp(' ')
tit = 'tv-Lh'; 
mprintar(strx.Lhtv3',in,tit);
disp('press any key to continue')
pause


disp(' ')
disp('******** Computation of Recursive Residuals   ********');
%compute recursive residuals
%set up regression matrices
X=Y; W=[]; 
Sigmax=strx.sigmarexct;
[L,p] = chol(Sigmax,'lower'); Lm1=pinv(L);
%set up system matrices
T=strx.Fsexct; Z=strx.Hsexct; G=Lm1; H=strx.Ksexct*Lm1;
%set up initial conditions
ndelta=0;                            %number of unit roots
[ins,i,ferror]=incossm(T,H,ndelta);
[Xt,Pt,g,M,initf,recrs,recr]=scakff(yd,X,Z,G,W,T,H,ins,i);
%plot recursive residuals
plot(recr(:,1)), legend('recr(:,1)'),pause
plot(recr(:,2)), legend('recr(:,2)'),pause
plot(recr(:,3)), legend('recr(:,3)'),pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag=6; ic=1; nr=strv.nparm;
disp(' ')
disp('******** Recursive Residuals:     ********');
str=mautcov(recr,lag,ic,nr);  
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('Q statistics:')
disp(str.qstat)

[m,n]=size(str.pval); t=1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all


