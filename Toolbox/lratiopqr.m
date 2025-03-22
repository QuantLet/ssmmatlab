function [lagsopt, ferror] = lratiopqr(y, x, seas, maxlag, minlag, prt)
% PURPOSE: performs sequential likelihood ratio tests in varmax(p,q,r)
%          models to determine optimal p, q and r. To obtain estimates of
%          the innovations, a long VARX model is estimated first. Then,
%          sequential likelihood ratio tests are performed to obtain an
%          optimal VARMAX(i,i,i), starting with i=minlag until i=maxlag.
%          Finally, with the same estimated innovations of the long VARX
%          model, sequential lr tests or bic computations are performed to
%          obtain an optimal VARMAX(p,q,r), starting with p, q and r equal 
%          to maxlag or i, depending on the user instructions, down until 
%          p, q and r equal zero.
%          All models are estimated using the same sample size: nobs -
%          maxlag or nobs - i - 1.
%---------------------------------------------------
% USAGE:  [lagsopt,ferror] = lratiopqr(y,x,seas,maxlag,minlag,prt)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%           x    = matrix of input variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%         seas   = seasonality
%           maxlag = the maximum lag length. If empty or positive on entry, 
%                    it is calculated by the program as the order of a
%                    VARMA(p,p,p) approximation.
%                    If negative, it is fixed to -maxlag
%           minlag = the minimum lag length
%           prt = flag for printing
%                    0 = no, 1 = yes
%                    (default = 0)
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

lagsopt = [];
ferror = 0;
flagmxl = 0;

if nargin ~= 6
    disp('wrong # of arguments to lratiopqrt');
    ferror = 1;
    return
end;

if (maxlag < 0)
    maxlag = - maxlag;
    flagmxl = 1;
end
if (maxlag < minlag)
    disp('maxlag < minlag in lratiopqrt');
    ferror = 2;
    return
end;

[nobs, neqs] = size(y);
if ~isempty(x)
    [mx, nx] = size(x);
    if (mx ~= nobs)
        disp('lratiopqrt: nobs in x-matrix not the same as y-matrix');
        ferror = 3;
        return
    end
else
    nx = 0;
end

%first identify a varmax(p,p,p)
[lagsopt, a, ferror] = lratiopppt(y, x, seas, maxlag, minlag, prt);
%flagmxl = 1 implies maxlag is fixed by the user
if isempty(maxlag) && (flagmxl == 0)  
    maxlag = lagsopt;
end

if (maxlag > lagsopt) && (flagmxl == 0)
    maxlag = lagsopt;
end

%initial model orders
p = maxlag;
q = maxlag;
pflag = 1;
qflag = 1;
if (nx > 0) 
    r = maxlag;
    rflag = 1;
else
    r = 0;
    rflag = 0;
end


% loop over lag lengths and do likelihood ratio tests
% the sample size is in all cases nobs - maxlag (Tsay, 2014, p. 62)
flag = 1;
while (flag == 1)
    if (p == 0), pflag = 0;
    end
    if (q == 0), qflag = 0;
    end
    if (r == 0), rflag = 0;
    end
    if (pflag == 1)
        maxordr1 = max([p, q, r]);
        clear str; %kro=repmat(maxordr1,1,neqs); str = matechelon(kro,neqs,nx);
        %   str = restrcmodel(str,neqs,nx,maxordr1,p,q,r);
        str = restrcmodel(neqs, nx, 1, [p, q, r], [0, 0, 0]);
        %   str = mhanris2(y,a,x,str);
        str = mhanris2(y(maxlag-maxordr1+1:end, :), a(maxlag-maxordr1+1:end, :), ...
            x(maxlag-maxordr1+1:end, :), str);
        resid1 = str.resid2; %residuals of the first HR regression
        maxordr2 = max([p - 1, q, r]);
        clear str; %kro=repmat(maxordr2,1,neqs); str = matechelon(kro,neqs,nx);
        %   str = restrcmodel(str,neqs,nx,maxordr2,p-1,q,r);
        str = restrcmodel(neqs, nx, 1, [p - 1, q, r], [0, 0, 0]);
        %   dmdr=maxordr1-maxordr2;
        %   str = mhanris2(y(1+dmdr:end,:),a(1+dmdr:end,:),x(1+dmdr:end,:),str);
        str = mhanris2(y(maxlag-maxordr2+1:end, :), a(maxlag-maxordr2+1:end, :), ...
            x(maxlag-maxordr2+1:end, :), str);
        resid2 = str.resid2; %residuals of the second HR regression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute likelihood ratio test
        % first get var-cov matrices for residuals
        epe1 = resid1' * resid1;
        epe2 = resid2' * resid2;
        %Reinsel's correction, p.107
        %   nbrofvar1=neqs*(p-1)+neqs*q+nx*(r+1)+neqs;
        nbrofvar1 = neqs * (p - 1) + neqs * q + nx * (r + 1);
        nbrofvar2 = neqs;
        nobse = nobs - maxordr1;
        const = double(nobse-nbrofvar1-(nbrofvar2 + neqs + 1)/2);
        lratio = const * (log(det(epe2)) - log(det(epe1)));
        % find marginal probability
        %     lprob = chis_prb(lratio,neqs*neqs);
        degf = double(nbrofvar2*neqs);
        if degf > 0
            lprob = gammp(degf*.5, lratio*.5); % the p-value is 1-lprob
        else
            lprob = 1;
        end
        if prt == 1
            out = [p, p - 1, lratio, 1 - lprob];
            fprintf(1, ...
                'p = %2d %2d, LR statistic = %16.4f, probability = %6.4g \n', out);
        end
        if (1 - lprob >= 0.05) % low p-values reject the null
            %  hypothesis that the i-th lag coefficient is zero
            p = p - 1;
        else
            pflag = 0;
        end
    end
    if (qflag == 1)
        maxordr1 = max([p, q, r]);
        clear str; %kro=repmat(maxordr1,1,neqs); str = matechelon(kro,neqs,nx);
        %   str = restrcmodel(str,neqs,nx,maxordr1,p,q,r);
        str = restrcmodel(neqs, nx, 1, [p, q, r], [0, 0, 0]);
        %   str = mhanris2(y,a,x,str);
        str = mhanris2(y(maxlag-maxordr1+1:end, :), a(maxlag-maxordr1+1:end, :), ...
            x(maxlag-maxordr1+1:end, :), str);
        resid1 = str.resid2; %residuals of the first HR regression
        maxordr2 = max([p, q - 1, r]);
        clear str; %kro=repmat(maxordr2,1,neqs); str = matechelon(kro,neqs,nx);
        %   str = restrcmodel(str,neqs,nx,maxordr2,p,q-1,r);
        str = restrcmodel(neqs, nx, 1, [p, q - 1, r], [0, 0, 0]);
        %   dmdr=maxordr1-maxordr2;
        %   str = mhanris2(y(1+dmdr:end,:),a(1+dmdr:end,:),x(1+dmdr:end,:),str);
        str = mhanris2(y(maxlag-maxordr2+1:end, :), a(maxlag-maxordr2+1:end, :), ...
            x(maxlag-maxordr2+1:end, :), str);
        resid2 = str.resid2; %residuals of the second HR regression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute likelihood ratio test
        % first get var-cov matrices for residuals
        epe1 = resid1' * resid1;
        epe2 = resid2' * resid2;
        %Reinsel's correction, p.107
        %   nbrofvar1=neqs*p+neqs*(q-1)+nx*(r+1)+neqs;
        nbrofvar1 = neqs * p + neqs * (q - 1) + nx * (r + 1);
        nbrofvar2 = neqs;
        nobse = nobs - maxordr1;
        const = double(nobse-nbrofvar1-(nbrofvar2 + neqs + 1)/2);
        lratio = const * (log(det(epe2)) - log(det(epe1)));
        % find marginal probability
        degf = double(nbrofvar2*neqs);
        if degf > 0
            lprob = gammp(degf*.5, lratio*.5); % the p-value is 1-lprob
        else
            lprob = 1;
        end
        if prt == 1
            out = [q, q - 1, lratio, 1 - lprob];
            fprintf(1, ...
                'q = %2d %2d, LR statistic = %16.4f, probability = %6.4g \n', out);
        end
        if (1 - lprob >= 0.05) % low p-values reject the null
            %  hypothesis that the i-th lag coefficient is zero
            q = q - 1;
        else
            qflag = 0;
        end
    end
    if (rflag == 1)
        maxordr1 = max([p, q, r]);
        clear str; %kro=repmat(maxordr1,1,neqs); str = matechelon(kro,neqs,nx);
        %   str = restrcmodel(str,neqs,nx,maxordr1,p,q,r);
        str = restrcmodel(neqs, nx, 1, [p, q, r], [0, 0, 0]);
        %   str = mhanris2(y,a,x,str);
        str = mhanris2(y(maxlag-maxordr1+1:end, :), a(maxlag-maxordr1+1:end, :), ...
            x(maxlag-maxordr1+1:end, :), str);
        resid1 = str.resid2; %residuals of the first HR regression
        maxordr2 = max([p, q, r - 1]);
        clear str; %kro=repmat(maxordr2,1,neqs); str = matechelon(kro,neqs,nx);
        %   str = restrcmodel(str,neqs,nx,maxordr2,p,q,r-1);
        str = restrcmodel(neqs, nx, 1, [p, q, r - 1], [0, 0, 0]);
        %   dmdr=maxordr1-maxordr2;
        %   str = mhanris2(y(1+dmdr:end,:),a(1+dmdr:end,:),x(1+dmdr:end,:),str);
        str = mhanris2(y(maxlag-maxordr2+1:end, :), a(maxlag-maxordr2+1:end, :), ...
            x(maxlag-maxordr2+1:end, :), str);
        resid2 = str.resid2; %residuals of the second HR regression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute likelihood ratio test
        % first get var-cov matrices for residuals
        epe1 = resid1' * resid1;
        epe2 = resid2' * resid2;
        %Reinsel's correction, p.107
        %   nbrofvar1=neqs*p+neqs*q+nx*r+neqs;
        nbrofvar1 = neqs * p + neqs * q + nx * r;
        nbrofvar2 = nx;
        nobse = nobs - maxordr1;
        const = double(nobse-nbrofvar1-(nbrofvar2 + neqs + 1)/2);
        lratio = const * (log(det(epe2)) - log(det(epe1)));
        % find marginal probability
        %     lprob = chis_prb(lratio,neqs*neqs);
        degf = double(nbrofvar2*neqs);
        if degf > 0
            lprob = gammp(degf*.5, lratio*.5); % the p-value is 1-lprob
        else
            lprob = 1;
        end
        if prt == 1
            out = [r, r - 1, lratio, 1 - lprob];
            fprintf(1, ...
                'r = %2d %2d, LR statistic = %16.4f, probability = %6.4g \n', out);
        end
        if (1 - lprob >= 0.05) % low p-values reject the null
            %  hypothesis that the i-th lag coefficient is zero
            r = r - 1;
        else
            rflag = 0;
        end
    end
    if (pflag == 0) && (qflag == 0) && (rflag == 0)
        flag = 0;
    end
end
lagsopt = [p, q, r];
if prt == 1
    fprintf(1, 'Estimated orders in VARMAX(p,q,r): ');
    fprintf(1, 'p = %2d, q = %2d, r = %2d \n\n', lagsopt);
end
