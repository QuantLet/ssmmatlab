function [lg, initres] = varident(y, maxlag, minlag, prt, x)
% PURPOSE: performs likelihood ratio test, bic,
%                   and aic for var model to determine
%                   optimal lag length
%---------------------------------------------------
% USAGE:  [lagsopt,initres] = varident(y,maxlag,minlag,prt,x)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           maxlag = the maximum lag length
%           minlag = the minimum lag length
%           prt    = flag for printing
%                    0 = no, 1 = yes (default = 0)
%           x      = optional matrix of variables (nobs x nx)
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
    disp('wrong # of arguments to varident');
    return
elseif nargin == 5
    xflag = 1;
elseif nargin == 4
    xflag = 0;
elseif nargin == 3
    prt = 0;
    xflag = 0;
end;

if maxlag < minlag
    disp('maxlag < minlag in varident');
    return
end;


[nobs, neqs] = size(y);
if xflag == 0
    nx = 0;
elseif xflag == 1
    if isempty(x)
        nx = 0;
        xflag = 0;
    else
        [mx, nx] = size(x);
    end
end
initresa = zeros(maxlag+1, neqs);

dn = double(nobs-maxlag); %we condition on the first maxlag observations

% loop over lag lengths and do likelihood ratio tests
% the sample size is in all cases nobs - maxlag (Tsay, 2014, p. 62)
% intialization
if xflag == 1 % case of deterministic variables
    resid1 = var_res(y(maxlag+1:end, :), 0, x(maxlag+1:end, :)); %initresa(i+1,:)=resid1(1,:);
else % case of no deterministic variables
    resid1 = var_res(y(maxlag+1:end, :), 0); %initresa(i+1,:)=resid1(1,:);
end;
sigmar1 = resid1' * resid1 / dn;
cra = log(det(sigmar1));
critaicm = cra;
critbicm = cra;
lratio = 0.;
lprob = 1.;
if prt == 1
    fprintf(1, '                          AIC                 BIC                LR                p-val\n');
    out = [0, critaicm, critbicm, lratio, 1 - lprob];
    fprintf(1, 'nlag = %2d%16.4f%16.4f%16.4f%16.4f \n', out);
end
lagsoptaic = 0;
lagsoptbic = 0;
lagsopt = minlag;

first = 1;
if minlag == 0, minlag = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=maxlag:-1:minlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = minlag:maxlag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if xflag == 1 % case of deterministic variables
        resid1 = var_res(y(maxlag-i+1:end, :), i, x(maxlag-i+1:end, :)); %initresa(i+1,:)=resid1(1,:);
        resid2 = var_res(y(maxlag-i+2:end, :), i-1, x(maxlag-i+2:end, :)); % restricted model
    else % case of no deterministic variables
        resid1 = var_res(y(maxlag-i+1:end, :), i); %initresa(i+1,:)=resid1(1,:);
        resid2 = var_res(y(maxlag-i+2:end, :), i-1); % restricted model
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %compute information criteria
    sigmar1 = resid1' * resid1 / dn;
    cra = log(det(sigmar1));
    dnparbn = double(i*neqs^2) / dn;
    pentaic = 2;
    pentbic = log(dn);
    critaic = cra + dnparbn * pentaic;
    critbic = cra + dnparbn * pentbic;
    
    % compute likelihood ratio test
    % first get var-cov matrices for residuals
    epe1 = resid1' * resid1;
    epe2 = resid2' * resid2;
    %Reinsel's correction, pp.107,
    const = double(nobs-maxlag-1.5-neqs*i);
    lratio = const * (log(det(epe2)) - log(det(epe1)));
    % find marginal probability
    %     lprob = chis_prb(lratio,neqs*neqs);
    lprob = gammp(neqs*neqs*.5, lratio*.5); % the p-value is 1-lprob
    if prt == 1
        out = [i, critaic, critbic, lratio, 1 - lprob];
        fprintf(1, 'nlag = %2d%16.4f%16.4f%16.4f%16.4f \n', out);
    end
    if critaic < critaicm
        critaicm = critaic;
        lagsoptaic = i;
    end
    if critbic < critbicm
        critbicm = critbic;
        lagsoptbic = i;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (1 - lprob > 0.05) && (first == 1) % low p-values reject the null hypothesis that the i-th lag coefficient is zero
        lagsopt = i - 1;
        first = 0;
        %  if (1-lprob < 0.05) && (first == 1) % low p-values reject the null hypothesis
        %   % we reject the null hypothesis that the i-th lag coefficient is zero
        %   lagsopt=i; first=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
lg = [lagsoptaic, lagsoptbic, lagsopt];
if prt == 1
    fprintf(1, 'selected orders by AIC, BIC and LR  = %3d%3d%3d\n', lg);
end
% if xflag == 1 % case of deterministic variables
%   resid1 = var_res(y,0,x); initresa(1,:)=resid1(1,:);
% else % case of no deterministic variables
%   resid1 = var_res(y,0);   initresa(1,:)=resid1(1,:);
% end;
initres = initresa(1:maxlag, :);
if nargout == 2
    for i = 1:maxlag
        if xflag == 0
            resid = var_res(y, i);
            initres(i, :) = resid(1, :);
        else
            resid = varx_res(y, i, x);
            initres(i, :) = resid(1, :);
        end
    end
end
