function [lagsopt, initres] = lratiocrax(y, maxlag, minlag, incr, prt, a, x)
% PURPOSE: performs sequential likelihood ratio tests in varmax(p,p,p)
%          models to determine optimal p, starting from
%          varmax(minlag,minlag,minlag) and proceeding up to
%          varmax(maxlag,maxlag,maxlag). All models are estimated using the
%          same sample size: nobs - maxlag.To estimjate the models, it
%          requires as input the estimated innovations.
%---------------------------------------------------
% USAGE:  lrratio(y,maxlag,minlag,incr,prt,a,x,seas)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%           maxlag = the maximum lag length
%           minlag = the minimum lag length
%           incr = the increment for the LR test
%           prt = flag for printing
%                    0 = no, 1 = yes
%                    (default = 0)
%           x    = matrix of input variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%           a    = matrix of estimated innovations (nobs x neqs)
%---------------------------------------------------
% RETURNS: lagsopt, the optimum number of lags
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

if nargin ~= 7
    error('wrong # of arguments to lratiocrax');
end

if maxlag < minlag
    error('maxlag < minlag in lratiocrax');
end


[nobs, neqs] = size(y);
if isempty(x)
    nx = 0;
else
    [mx, nx] = size(x);
end
initresa = zeros(maxlag+1, neqs);
% sigmarv=cov(a,1);

% loop over lag lengths and do likelihood ratio tests
% the sample size is in all cases nobs - maxlag (Tsay, 2014, p. 62)
lagsoptlr = maxlag;
first = 1;
if  (minlag == 0)
    minlag = 1;
end

% BIC 
% start with zero lags
lagsoptbic = 0;
dn = double(nobs-maxlag); %we condition on the first maxlag observations
pentbic = log(dn);
resid1 = y(maxlag+1:end, :); 
sigmar1 = resid1' * resid1 / dn;
cra = log(det(sigmar1));
critbicm = cra;
lratio = 0.;
lprob = 1.;
if prt == 1
    fprintf(1, '                  BIC              LR              p-val\n');
    out = [0, critbicm, lratio, 1 - lprob];
    fprintf(1, 'nlag = %2d%16.4f%16.4f%16.4f \n', out);
end



% for i=maxlag:-incr:minlag
for i = minlag:incr:maxlag
    % adjust nobs to feed the lags
    %  nobse = nobs-i;
    %  resid1 = var_resax(y,i,a,x);
    %  resid2 = var_resax(y(2:end,:),i-1,a(2:end,:),x(2:end,:)); % restricted model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % i
    imincr = i - incr;
    if (imincr < 0)
        break
    end
    r1 = imincr;
    r2 = incr;   
    clear str;
    kro = repmat(i, 1, neqs);
    str = matechelon(kro, neqs, nx);
    %  str = mhanris2(y,a,x,str);
    str = mhanris2(y(maxlag-i+1:end, :), a(maxlag-i+1:end, :), ...
        x(maxlag-i+1:end, :), str);
    resid1 = str.resid2; %residuals of the unrestricted HR regression
    clear str;
    kro = repmat(imincr, 1, neqs);
    str = matechelon(kro, neqs, nx);
    %  str = mhanris2(y(incr+1:end,:),a(incr+1:end,:),x(incr+1:end,:),str);
    str = mhanris2(y(maxlag-i+incr+1:end, :), ...
        a(maxlag-i+incr+1:end, :), x(maxlag-i+incr+1:end, :), str);
    resid2 = str.resid2; %residuals of the restricted HR regression
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store initial residual
    %  initresa(i+1,:)=resid1(1,:);
    % compute likelihood ratio test
    % first get var-cov matrices for residuals
    epe1 = resid1' * resid1;
    epe2 = resid2' * resid2;
    %Reinsel's correction, p.107:
    %Y=X_1B_1 + X_2B_2 + e,
    %X_1 is (T x r_1), X_2 is (T x r_2). H_0: B_2=0
    %M=-[T - r_1 - (r_2 + k +1)/2]log(U).
    %  nbrofvar=neqs+neqs; nbrofvar1=nx+neqs+nbrofvar*r1+nx*imincr;
    %  nbrofvar2=nbrofvar*r2+nx*incr;
    %  const = double(nobse-nbrofvar1-(nbrofvar2+neqs+1)/2);
    nbrofvar = neqs + neqs + nx;
    nbrofvar1 = nx + nbrofvar * r1; 
    %the extra nx above is to account for the zero lag in inputs
    nbrofvar2 = nbrofvar * r2;
    const = double(nobs-maxlag-nbrofvar1-(nbrofvar2 + neqs + 1)/2);
    lratio = const * (log(det(epe2)) - log(det(epe1)));
    % find marginal probability
    %     lprob = chis_prb(lratio,neqs*neqs);
    degf = double(nbrofvar2*neqs);
    if degf > 0
        lprob = gammp(degf*.5, lratio*.5); % the p-value is 1-lprob
    else
        lprob = 1;
    end
    % if prt == 1
    %     out = [i, imincr, lratio, 1 - lprob];
    %     fprintf(1, 'nlag = %2d %2d, LR statistic = %16.4f, probability = %6.4g \n', out);
    % end
    if (1 - lprob > 0.05) && (first == 1) % low p-values reject the null hypothesis that the i-th lag coefficient is zero
        lagsoptlr = i - 1;
        first = 0; 
        %  if (1-lprob < 0.05) && (first == 1) % low p-values reject the null hypothesis
        %   % we reject the null hypothesis that the i-th lag coefficient is zero
        % %   if (noninv1 == 0) & (noninv2 == 0)
        %    lagsopt=i; first=0;
        %   end
    end


    %compute BIC information criteria
    sigmar1 = epe1 / dn;
    cra = log(det(sigmar1));
    dnparbn = double(i*(neqs^2 + neqs^2) + i*neqs*nx + neqs*nx) / dn;
    critbic = cra + dnparbn * pentbic;
    if critbic < critbicm
        critbicm = critbic;
        lagsoptbic = i;
    end
    if prt == 1
        out = [i, critbic, lratio, 1 - lprob];
        fprintf(1, 'nlag = %2d%16.4f%16.4f%16.4f \n', out);
    end

end
%use the bic to decide
lagsopt = lagsoptbic;
if (lagsopt == 0)
    lagsopt = 1;
end
%use lr to decide
% lagsopt = lagsoptlr;
% % if (first == 1)
% %     lagsopt = minlag;
% % end
if prt == 1
    out = [lagsoptbic, lagsoptlr];
    fprintf(1, 'selected orders by BIC and LR  = %3d%3d\n', out);
end
% resid1 = var_resax(y,0,a,x); initresa(1,:)=resid1(1,:);
initres = initresa(1:maxlag, :);
if nargout == 2
    for i = 1:maxlag
        clear str;
        kro = repmat(i, 1, neqs);
        str = matechelon(kro, neqs, nx);
        str = mhanris2(y, a, x, str);
        resid1 = str.resid2; %residuals of the unrestricted HR regression
        initres(i, :) = resid1(1, :);
    end
end
