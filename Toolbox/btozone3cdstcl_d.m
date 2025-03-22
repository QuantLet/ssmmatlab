
%
% Ozone series from Box and Tiao (1973)
%
% 1) An ARIMA model is estimated first (specification file btozone3.m in
%    subdirectory spec).
% 2) The canonical decomposition is obtained
% 3) A low-pass filter (HP filter) is applied to the canonical trend. A
%    smooth trend an a cycle are obtained whose sum is the canonical trend.
% 4) A band-pass filter is applied to the canonical trend. A smooth cycle
%    and a trend are obtained whose sum is the canonical trend.
%
% Steps 3) and 4) are based on Gómez, V. (2001). ''The use of butterworth
% filters for trend and cycle estimation in economic time series''.
% Journal of Business and Economic Statistics, 19, 365–373. The algorithms
% for many of the procedures are described in Gómez (2016), ''Multivariate
% Time Series With Linear State Space Structure''. Springer International
% Publishing AG: Switzerland
%

clear;

%arima estimation
out = arimaestos('btozone3');

%
%Signal extraction using the canonical decomposition
%
mY = size(out.model.Y, 2);
Ycomp = cell(1, mY);
for i = 1:mY
    Ycomp{i} = [];
end
for i = mY - 3:mY - 3
    Ycomp{i} = 'trend';
end
Ycomp{mY} = 'irreg'; %assign the outlier to the irregular
outa = arimasigex(out, Ycomp);

%plot spectra
if out.gft == 1
    plotspcd(outa)
end

thrc = outa.compcd.ptnum; % trend-cycle numerator
phir = outa.compcd.ptden; % trend-cycle denominator
%   outa.compcd.ptnur;    % number of nonstationary roots in phir
sigma2r = outa.compcd.ptvar; % variance of the trend-cycle innovations (*)
thsc = outa.compcd.stnum; % seasonal numerator
phis = outa.compcd.stden; % seasonal denominator
%   outa.compcd.stnur; % number of nonstationary roots in phis
sigma2s = outa.compcd.stvar; % variance of the seasonal innovations (*)
thtc = outa.compcd.rt; % transitory component (MA term)
sigma2t = outa.compcd.rtvar; % variance of the transitory component innovations (*)
sigma2i = outa.compcd.itvar; % variance of the irregular component (*)
phitc = outa.compcd.phi; % stationary AR trend polynomial
%(*) in units of the series model innovations

disp('trend-cycle numerator:')
disp(thrc)
disp('trend-cycle denominator:')
disp(phir)
disp('variance of the trend-cycle innovations (*)')
disp(sigma2r)
pause
disp('seasonal numerator:')
disp(thsc)
disp('seasonal denominator:')
disp(phis)
disp('variance of the seasonal innovations (*)')
disp(sigma2s)
pause
disp('stationary AR trend polynomial')
disp(phitc)
pause
disp('variance of the irregular component (*)')
disp(sigma2i)
if ~isempty(thtc)
    disp('transitory numerator:')
    disp(thtc)
    phit = 1.;
    disp('transitory denominator:')
    disp(phit)
    disp('variance of the transitory innovations (*)')
    disp(sigma2t)
end
disp('(*) in units of var(A)')

% ***********************************************************************
%                        LOW-PASS FILTERS
%
% design of a low-pass filter to obtain a smooth trend and a cycle. The
% filter is specified giving Lambda and Di. The sine But. filter is the
% Hodrick-Prescott filter. The filter will be applied to the canonical
% trend.
Lambda = 127864.652;
Di = 2;
% Parameter lambda for monthly data, see Gomez (2001)
%"The Use of Butterworth Filters for Trend and Cycle Estimation in Economic
%Time Series"
% select from the two following lines the tangent or the sine But. filter
% by uncommenting or commenting the appropriate lines.
% [compf,ferror]=dtanbut([],[],[],Di,[],Lambda);
[compbst, ~] = dsinbut([], [], [], Di, [], Lambda);
% plot gain function of both low-pass filters
figure
ggsintanbut([], [], [], compbst.Di, compbst.Thetac)
pause
close all

%estimation of smooth trend and cycle by applying a low-pass filter to the
%trend-cycle component of the canonical decomposition.
filter = 'lp';
outb = arimasigextc(outa, compbst, filter);

%************************************************************************
%                    BAND-PASS FILTER
%
% design of band-pass filter to obtain a well defined cycle and a
% relatively smooth trend. Frequencies are expressed divided by pi. The
% filter will be applied to the canonical trend.
D(1) = .12;
D(2) = .12;
xp1 = .02;
xp2 = .08;
xs = .15; %Di=3
% Tangent band-pass filter
[compbp, ferror] = dbptanbut(D, xp1, xp2, xs);
% plot gain function of the tangent band-pass filter
figure
ggbptanbut(D, xp1, xp2, xs, compbp.Di, compbp.Alph, compbp.Lambda)
pause
close all

%estimation of smooth cycle and trend by applying a band-pass filter to the
%trend-cycle component of the canonical decomposition.
filter = 'bp';
outc = arimasigextc(outa, compbp, filter);


%additional plots
datei = outa.datei;
% plot both the smooth and the canonical trend
trend = outa.StochCc(:, 1);
bsttrend = outb.StochCctc(:, 1);
vnames = char('USIPI trend', 'USIPI bsttrend (low-pass filter)');
figure
tsplot([trend, bsttrend], datei, vnames);
pause

%plot both the band-pass filter and the canonical trend
figure
bptrend = outc.StochCctc(:, 1);
vnames = char('USIPI trend', 'USIPI bptrend (band-pass filter)');
figure
tsplot([trend, bptrend], datei, vnames);
pause

%plot both the smooth and the band-pass filter trend
figure
vnames = char('USIPI bsttrend (low-pass filter)', 'USIPI bptrend (band-pass filter)');
figure
tsplot([bsttrend, bptrend], datei, vnames);
pause

% plot both the band-pass filter and the low-pass filter cycle
bstcycle = outb.StochCctc(:, 2);
bpcycle = outc.StochCctc(:, 2);
figure
vnames = char('USIPI btcycle', 'USIPI bpcycle');
tsplot([bstcycle, bpcycle], datei, vnames);
pause
close all
