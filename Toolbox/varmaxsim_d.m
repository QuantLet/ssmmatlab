%script file to simulate a VARMAX model
%
% Let the VARMAX model be
%
%   Phi(B) y_t = Gamma(B) x_t + Theta(B) a_t,
%
% where B is the backshift operator, By_t = y_{t-1}. Then, the series y_t can be
% simulated as the sum of two variables, u_t and v_t, namely
%
%   u_t = Phi(B)^{-1}Gamma(B)x_t,  v_t = Phi(B)^{-1}Theta(B)a_t,
%
% where x_t  is given or follows an arbitrary model, phi(B) x_t = theta(B) b_t.
%
% The VARMAX model can also be given as a decoupled model,
%
%   y_t = Alpha(B)^{-1}Gamma(B)x_t +  Phi(B)^{-1}Theta(B)a_t
%
% In this case, make the appropiate changes in the filter for the inputs.
%

clear
l = 50; %number of initial observations to be discarded
N = 300; %number of observations of the simulated series
m = 2; %number of inputs
s = 3; %number of outputs
seed = 20; %this is to always generate the same series

%polynomial matrices Phi and Theta
Phi(:, :, 1) = eye(s);
Theta(:, :, 1) = eye(s);
Phi(:, :, 2) = -[0.8, 0.4, 0.; -0.3, 0.6, 0.; 0., 0., -1.];

%polynomial matrix Gamma
Gamma(:, :, 1) = [.1, -.3; -.5, .8; .2, -.7];

%first, simulate v_t
%covariance matrix of the a_t innovartions
S = [2.0, 0.5, 0.; 0.5, 1.0, 0.; 0., 0., .3];

%simulate v_t
[v, ferror] = varmasim(l, N, Phi, Theta, S, seed);

%second, simulate x_t
%polynomial matrices phi and theta
phi(:, :, 1) = eye(m);
phi(:, :, 2) = -eye(m);
theta(:, :, 1) = eye(m);
theta(:, :, 2) = [-.2, -.3; .6, -1.1];
%covariance matrix of the b_t innovartions
sigma = .2;
Sx = eye(m) * sigma;

%simulate x_t
seed = seed + 2;
[x, ferror] = varmasim(l, N, phi, theta, Sx, seed);


%third, filter inputs by Phi(B)^{-1}Gamma(B)

% %filter inputs without using the input model
% [u,su] = varmafilp(x,Phi,Gamma);
%filter inputs using the input model
freq = 1;
phix = phi;
thx = theta;
Phix(:, :, 1) = eye(m);
Thx(:, :, 1) = eye(m);
[u, su] = varmafilp(x, Phi, Gamma, phix, thx, Phix, Thx, Sx, freq);

%simulate y_t
y = u + v;

%estimate model using HR method
seas = 1;
[strv, ferror] = estvarmaxpqrPQR(y, x, seas, [1, 0, 0], [0, 0, 0]);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'Phi';
strt = 1;
mprintar(strv.phis(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Gamma';
strt = 0;
mprintar(strv.gammas(:, :, 1), in, tit, strt);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar2, in, tit);
disp(' ')
disp('***** T-values  *****');
disp(' ');
tit = 'tv-Phi';
strt = 1;
mprintar(strv.phitv(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Gamma';
strt = 0;
mprintar(strv.gammatv(:, :, 1), in, tit, strt);

%estimate the Kronecker indices for the series
maxorder = 3;
prt = 0;
seas = 1;
hr3 = 0;
[order, kro, scm] = varmaxscmidn(y, x, seas, maxorder, hr3, prt);
disp(' ')
disp('estimated Kronecker Indices for the series')
disp('using function "varmaxscmidn":')
disp(kro)
disp('press any key to continue')
pause

disp('Estimate a simplified VARMAX model with these Kronecker indices')
disp(' ')
%estimate VARMA model in echelon form for the and
%aliminate some nonsignificant parameters
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [1., 1.];
strv = estvarmaxkro(y, x, seas, kro, hr3, finv2, mstainv, nsig, tsig);

maxkro = max(kro) + 1;
disp(' ');
disp('***** Estimated VARMAX Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 0;
mprintar(strv.phis3(:, :, 1:maxkro), in, tit, strt);
disp(' ')
tit = 'gamma';
strt = 0;
mprintar(strv.gammas3(:, :, 1:maxkro), in, tit, strt);
disp(' ')
tit = 'th';
strt = 0;
mprintar(strv.thetas3(:, :, 1:maxkro), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 0;
mprintar(strv.phitv3(:, :, 1:maxkro), in, tit, strt);
disp(' ')
tit = 'tv-gamma';
strt = 0;
mprintar(strv.gammatv3(:, :, 1:maxkro), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 0;
mprintar(strv.thetatv3(:, :, 1:maxkro), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar3, in, tit);
