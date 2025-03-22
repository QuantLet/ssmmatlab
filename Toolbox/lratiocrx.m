function [lagsopt, initres] = lratiocrx(y, maxlag, minlag, prt, x)
% PURPOSE: performs likelihood ratio tests for varx model to determine
% optimal lag length. All models are estimated  using the same sample size:
% nobs - maxlag.
%---------------------------------------------------
% USAGE:  [lagsopt,initres] = lratiocrx(y,maxlag,minlag,prt,x)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           maxlag = the maximum lag length
%           minlag = the minimum lag length
%           prt    = flag for printing
%                    0 = no, 1 = yes (default = 0)
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%---------------------------------------------------
% RETURNS: lagsopt = the optimum number of lags
%          initres = an (maxlag x neqs) matrix of initial residuals
%                    corresponding to the estimated VARXs of order 1,2,...,
%                    maxlag.
%---------------------------------------------------
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

if nargin > 5
    disp('wrong # of arguments to lratiocrx');
    return
end;

if maxlag < minlag
    disp('maxlag < minlag in lratiocrx');
    return
end;


if isempty(x)
    [lagsopt, initres] = lratiocr(y, maxlag, minlag, prt);
    %  lagsopt=[]; initres=[];
    %  disp('matrix x should not be empty en lratiocrx')
    %  disp('use function lratiocr instead')
    return
else
    [mx, nx] = size(x);
end
[nobs, neqs] = size(y);
initresa = zeros(maxlag+1, neqs);


% loop over lag lengths and do likelihood ratio tests
% the sample size is in all cases nobs - maxlag (Tsay, 2014, p. 62)
lagsopt = maxlag;
first = 1;
if minlag == 0, minlag = 1;
end
% for i=maxlag:-1:minlag
for i = minlag:maxlag
    % adjust nobs to feed the lags
    %  nobse = nobs-i;
    %  resid1 = varx_res(y,i,x); initresa(i+1,:)=resid1(1,:);
    %  resid2 = varx_res(y(2:end,:),i-1,x(2:end,:)); % restricted model
    resid1 = varx_res(y(maxlag-i+1:end, :), i, x(maxlag-i+1:end, :)); %initresa(i+1,:)=resid1(1,:);
    resid2 = varx_res(y(maxlag-i+2:end, :), i-1, x(maxlag-i+2:end, :)); % restricted model
    % compute likelihood ratio test
    % first get var-cov matrices for residuals
    epe1 = resid1' * resid1;
    epe2 = resid2' * resid2;
    %Reinsel's correction, p.107
    nbrofvar2 = neqs + nx;
    nbrofvar1 = nbrofvar2 * (i - 1) + nx;
    const = double(nobs-maxlag-nbrofvar1-(nbrofvar2 + neqs + 1)/2);
    lratio = const * (log(det(epe2)) - log(det(epe1)));
    % find marginal probability
    %     lprob = chis_prb(lratio,neqs*neqs);
    lprob = gammp(nbrofvar2*neqs*.5, lratio*.5); % the p-value is 1-lprob
    if prt == 1
        out = [i, i - 1, lratio, 1 - lprob];
        fprintf(1, 'nlag = %2d %2d, LR statistic = %16.4f, probability = %6.4g \n', out);
    end
    if (1 - lprob > 0.05) && (first == 1) % low p-values reject the null hypothesis that the i-th lag coefficient is zero
        lagsopt = i - 1;
        first = 0;
        %  if (1-lprob < 0.05) && (first == 1) % low p-values reject the null hypothesis
        %   % we reject the null hypothesis that the i-th lag coefficient is zero
        %   lagsopt=i; first=0;
    end
end
% if (first == 1)
%     lagsopt = minlag;
% end
resid1 = varx_res(y, 0, x);
initresa(1, :) = resid1(1, :);
initres = initresa(1:maxlag, :);
if nargout == 2
    for i = 1:maxlag
        resid = varx_res(y, i, x);
        initres(i, :) = resid(1, :);
    end
end
