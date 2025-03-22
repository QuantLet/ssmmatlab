function str = mautcov(y, lag, ic, nr)
%
%
%      This function computes the autocovariance matrices of y(t) up to lag lag.
% ic = 1: compute the autocorrelation matrices
%      0: do not compute the autocorrelation matrices
% nr = number of parameters for Hosking's portmanteau statistic
%
%      returns: a structure str containing
%               c0: covariance at lag 0
%               cv: three dimensional array containing the autocovariance
%                    matrices up to lag lag
%                r: the autocorrelation matrices up to lag lag
%              sgn: matrices containing the significance of the
%                   autocorrelations according to the 2/sqrt(n) limits
%             sgnt: matrix containing all sgn matrices together
%              r0: the autocorrelation matrix for lag zero
%  In addition, if nr is input,
%             qstat: the Q statistics up to lag
%             pval : the p-values of the Q statistics
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
[n, m] = size(y);
cv = zeros(m, m, lag);
% mean
mu = mean(y);
% for i=1:n
%  y(i,:)=y(i,:)-mu;
% end
y = y - repmat(mu, n, 1); % mean corrected series
% variance
c0 = y' * y / n;
% autocovariances
for i = 1:lag
    cv(:, :, i) = y(1:n-i, :)' * y(1+i:n, :) / n;
end
% autocorrelation matrices
if ic == 1
    r = zeros(size(cv));
    sgn = blanks(m);
    sgnt = [];
    sl = 2 / sqrt(n);
    stdv = sqrt(diag(c0));
    for i = 1:lag
        for j = 1:m
            for k = 1:m
                rr = cv(j, k, i) / (stdv(j) * stdv(k));
                if (rr > sl)
                    sgn(j, k, i) = '+';
                elseif (rr < -sl)
                    sgn(j, k, i) = '-';
                else
                    sgn(j, k, i) = '.';
                end
                r(j, k, i) = rr;
            end
        end
        sgnt = [sgnt, sgn(:, :, i)];
    end
    disp(' ');
    disp('******** Autocorrelation function ********');
    disp('         Signs        ACF')
    tit = [];
    for i = 1:m;
        tit = char(tit, ['ser #', int2str(i)]);
    end
    tit = tit(2:m+1, :);
    blk = abs(' ');
    tit = [tit, char(ones(m, 10-size(tit, 2))*blk)];
    tits = [tit(:, 1:10), char(ones(m, 2)*blk)];
    sep = char(ones(m, 3)*blk);
    sep(:, 2) = char(ones(m, 1)*abs('|'));
    for i = 1:lag
        a = char(ones(m, 7*m)*blk);
        for j = 1:m
            a(j, :) = sprintf('%6.2f ', r(j, :, i));
        end
        disp([' k = ', int2str(i)]);
        disp([tits, sgn(:, :, i), sep, a]);
        disp(' ');
    end
    disp('All signs of ACF together ');
    disp(sgnt);
    disp(' ');
else
    r = [];
    sgn = [];
    sgnt = [];
end
str.c0 = c0;
str.cv = cv;
str.r = r;
str.sgn = sgn;
str.sgnt = sgnt;
r0 = zeros(m);
if (ic == 1)
    for j = 1:m
        for k = 1:m
            r0(j, k) = c0(j, k) / (stdv(j) * stdv(k));
        end
    end
end
str.r0 = r0;


%Hosking's (1980) portmanteau statistics
if nargin == 4
    r0i = r0 \ eye(size(r0));
    V = kron(r0i, r0i);
    qstat = zeros(lag, 1);
    sr = zeros(lag, 1);
    df = zeros(lag, 1);
    pval = zeros(lag, 1);
    for i = 1:lag
        rv = vec(r(:, :, i)');
        sr(i) = rv' * V * rv;
        qstat(i) = sum(sr(1:i)./(n - 1:-1:n - i)') * n * n;
        df(i) = max(0, m^2*i-nr);
        if df(i) > 0
            pval(i) = 1 - gammp(df(i)*.5, qstat(i)*.5);
        else
            pval(i) = 1;
        end
    end
    str.qstat = qstat;
    str.pval = pval;
end
