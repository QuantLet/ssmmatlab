function z = varimass(l, N, H, F, K, A, Sigma, Xi, stda, seed)
%
%        This function generates a VARMA model
%        It uses the state space representation
%        stda: standard deviation of the innovations in percentage of the levels (default: 1)
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%


z = [];
[nh, mh] = size(H);
[nx, mx] = size(Xi);
if nargin < 8
    disp('number of arguments must be at least 9 in varimass')
    return
elseif nargin == 8
    stda = eye(nh);
elseif nargin == 10
    if isOctave
        randn('state', seed);
    else
        rng(seed); %use the same start for random number generator
    end
end

if isempty(A)
    d = 0;
else
    [na, d] = size(A);
end

if d > 0
    delta = randn(d, 1);
else
    delta = 0;
end
U = randn(mx, 1);
[r, junk] = size(F);
if ~isempty(Sigma)
    % Sigma=sparse(Sigma);
    [rs, ~] = size(Sigma);
    [LSigma, p] = chol(Sigma);
    % LSigma=full(LSigma);
else
    p = 0;
    LSigma = [];
end
if p > 0
    LSigmat = zeros(rs);
    LSigmat(1:p-1, 1:p-1) = LSigma';
    LSigma = LSigmat;
    % LSigma = [LSigma', zeros(rs, rs-p+1)];
else
    LSigma = LSigma';
end
% format long
% A
% delta
% p
% LSigma
% Sigma=full(Sigma)
% LSigma*LSigma'
% U
% Xi
% stda
% format short


%initial state vector
if (d == 0)
    xt = Xi * LSigma * U;
elseif (d == mh)
    xt = A * delta;
else
    xt = A * delta + Xi * LSigma * U;
end

% format long
% xt
% format short

zz = zeros(l+N, nh);
% z=zeros(N,nh);
%
% Run the state space model
%
H = sparse(H);
F = sparse(F);
for i = 1:l + N
    at = stda * randn(nh, 1);
    % general state space model
    zz(i, :) = (H * xt + at)';
    xt = F * xt + K * at;
end
z = zz(l+1:N+l, :);
% % make the series positive if it is nonstationary
% mz=min(z); mmz=min(mz);
% if (d > 0) & (mmz < 0)
%  z=z-mmz*1.1*ones(size(z));
% end
