function [x, fjac, ff, g, iter, conf] = marqdt(info, x, varargin)
%
% This function minimizes a non-linear sum of squares function using
% Levenberg-Marquard's method.
% Input arguments:
%   info structure containing function names and optimization options
%   .f  :   a function to evaluate the vector ff of individual functions
%           such that ff'*ff is minimized
%   .tr :   >0 x is passed from marqdt to f but not passed from f to marqdt
%           =0 x is passed from marqdt to f and passed from f to marqdt
%   .tolf:  a parameter used for stopping
%   .jac:   =1 evaluation of jacobian and gradient at the solution is performed
%           =0 no evaluation of jacobian and gradient at the solution is performed
% .maxit:   maximum number of iterations
%   .nu0:   initial value of the nu parameter
%   .prt:   =1 printing of results
%           =0 no printing of results
% x:        a vector containing the initial parameter values
% varargin: arguments to be passed to function f
%
% Output arguments:
% x        a vector of (untransformed) parameters containing the solution
% fjac     a matrix containing the jacobian at the solution
% ff       a vector containing the individual functions at the solution
% g        a vector containing the (1/2)gradient at the solution
% iter     a scalar whose value is the number of iterations
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


infoz.f = info.f;
if isfield(info, 'F')
    infoz.F = info.F;
    anjac = 1;
else
    anjac = 0;
end
% set defaults
infoz.tr = 0;
infoz.tol = 1e-4;
infoz.jac = 0;
infoz.maxit = 100;
infoz.nu0 = 0;
infoz.prt = 0;

if ~isempty(info)
    if ~isstruct(info)
        error('marqdt: options should be in a structure variable');
    end;
    % parse options
    fields = fieldnames(info);
    nf = length(fields);
    for i = 2:nf
        if strcmp(fields{i}, 'tr')
            infoz.tr = info.tr;
        elseif strcmp(fields{i}, 'tolf')
            infoz.tol = info.tolf;
        elseif strcmp(fields{i}, 'jac')
            infoz.jac = info.jac;
        elseif strcmp(fields{i}, 'maxit')
            infoz.maxit = info.maxit;
        elseif strcmp(fields{i}, 'nu0')
            infoz.nu0 = info.nu0;
        elseif strcmp(fields{i}, 'prt')
            infoz.prt = info.prt;
        end;
    end;
else
    % rely on default options
end;


nu = infoz.nu0;
conv = 0;
n = length(x);
iter = 0;
conf = 0;

fHandle = str2func(infoz.f);
[ff] = fHandle(x, varargin{:});

%  ff=feval(infoz.f,x,varargin{:});
if anjac == 0
    [fjac, g] = fdjac2(infoz, x, ff, varargin{:});
    %[fjac,ff,g] = fdjac(infoz,x,varargin{:});
else
    fHandleF = str2func(infoz.F);
    [fjac, g] = fHandleF(x, varargin{:});
    %   [fjac,g] = feval(infoz.F,x,varargin{:});
end

conf = conf + n + 1;
while conv == 0
    iter = iter + 1;
    rs0 = ff' * ff;
    flag = 0; 
    while flag == 0
        if nu > 0, wk = [fjac; sqrt(nu) * eye(n)];
        else wk = fjac;
        end
        [Q, RR] = qr(wk);
        R = RR(1:n, 1:n);
        if rank(R) == n, flag = 1;
        else
            if nu == 0, nu = 0.01;
            else nu = 4 * nu;
            end
        end
    end
    delta = R \ (R' \ (-g));
    xd = x + delta';
    %  xd
    if infoz.tr == 0
        fHandle = str2func(infoz.f);
        [ffd, xd] = fHandle(xd, varargin{:});
        %    [ffd,xd]=feval(infoz.f,xd,varargin{:});
        conf = conf + 1;
    else
        fHandle = str2func(infoz.f);
        ffd = fHandle(xd, varargin{:});
        %    ffd=feval(infoz.f,xd,varargin{:});
        conf = conf + 1;
    end
    rs = ffd' * ffd;
    rr = abs(rs0-rs) / (1 + abs(rs0));
    rp = max(abs(x-xd)./(ones(size(x)) + abs(x)));
    if infoz.prt == 2
        fprintf(1, ' it= %2i', iter);
        fprintf(1, ' rs= %12.8f', rs);
        fprintf(1, '    rr= %12.8f', rr);
        fprintf(1, '    rp= %12.8f', rp);
        fprintf(1, '    nu= %12.8f\n', nu);
    end
    if (rr <= infoz.tol) || (rp <= sqrt(infoz.tol)) || (iter >= infoz.maxit)
        conv = 1;
        x = xd;
        ff = ffd;
        if infoz.jac == 1
            %  [fjac,g] = fdjac2(infoz,x,ff,varargin{:});
            if anjac == 0
                [fjac, g] = fdjac2(infoz, x, ff, varargin{:});
            else
                fHandleF = str2func(infoz.F);
                [fjac, g] = fHandleF(x, varargin{:});
                %      [fjac,g] = feval(infoz.F,x,varargin{:});
            end
            conf = conf + n;
        end
    else
        s = (rs0 - rs) / (nu * (delta' * delta) - delta' * g);
        if s < 0.25
            if (nu == 0), nu = .001;
            else nu = 4 * nu;
            end
        elseif s > .75
            nu = nu / 2;
        end
        if s > 0
            x = xd;
            ff = ffd;
            %   [fjac,g] = fdjac2(infoz,x,ff,varargin{:});
            if anjac == 0
                [fjac, g] = fdjac2(infoz, x, ff, varargin{:});
            else
                fHandleF = str2func(infoz.F);
                [fjac, g] = fHandleF(x, varargin{:});
                %      [fjac,g] = feval(infoz.F,x,varargin{:});
            end
            conf = conf + n;
        end
    end
end
if infoz.prt > 0
    fprintf(1, '\n rs= %12.8f', rs);
    fprintf(1, '    rr= %12.8f', rr);
    fprintf(1, '    rp= %12.8f', rp);
    fprintf(1, '    nu= %12.8f\n', nu)
end
