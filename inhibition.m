function [sigs, T2]=inhibition(para)

% generate simulaiton data

% The configuration should contain:
% para.num                    :   number of alpha cycles 
% para.Fs                     :   sampling freqeucy 
% para.hf                     :   gamma frequency 
% para.shift                  :   create directionality either alpha leads gamma or alpha lags gamma 
                                  % 'lead', 'lag' or [];
% para.timdiff                :   alpha and gamma time lag in ms
% para.a(slope)/c(threshold)  :   Sigmoid function parameter sigmf(x, [a, c]) = 1./(1 + exp(-a*(x-c)))
% para.ampmean/ampstd         :   alpha amplitude modulation gassuain
%                                 distribution.ampmean: mean(default 2); ampstd; std(default 0.2).
% para.tmean/tstd             :   alpha cycles duration gassuain
%                                 distribution. tmean: mean(default 0.1s); tstd: std(defacult 0.02s) .

% para.n1                     :   Gaussian white noise level
% para.n2                     :   pinck noise level  

% output: 
% sigs                        :   four signals {alpha0;alphas;gamma;mix};we are looking CFC at mix channel sigs(4,:)
% T2                          :   time series index 

% copyright @Haiteng Jiang; When you use this code ,please cite this paper
% Jiang, H., Bahramisharif, A, van Gerven, M.A.J. and Jensen, O (2015). 
% Measuring Directionality between Neuronal Oscillations of Different Frequencies.
% NeuroImage, 118:359-367

 if isfield(para, 'num')     
     num = para.num; 
 end
 
if isfield(para, 'Fs')     
     Fs = para.Fs; 
end

if isfield(para, 'hf')     
     hf = para.hf; 
end


if isfield(para, 'shift')     
     shift  = para.shift; 
else 
      shift = [];
end

if isfield(para, 'ampmean')     
     ampmean  = para.ampmean; 
else 
     ampmean =  2;
end


if isfield(para, 'ampstd')     
     ampstd  = para.ampstd; 
else 
     ampstd =  0.2;
end

if isfield(para, 'tmean')     
     tmean  = para.tmean; 
else 
     tmean =  0.1;
end


if isfield(para, 'tstd')     
     tstd  = para.tstd; 
else 
     tstd =  0.02;
end

if isfield(para, 'timdiff')     
     timdiff = para.timdiff; 
end
 

if isfield(para, 'a')     
     a = para.a; 
end

if isfield(para, 'c')     
     c = para.c; 
end

if isfield(para, 'n1')     
     n1 = para.n1;     
end

if isfield(para, 'n2')     
     n2 = para.n2;     
end


dt         =  1/Fs;
Time       =  tmean + tstd.*randn(num,1);        % alpha range cycle length 80-120 ms
amps       =  ampmean+ampstd*randn(num,1);     % fluctuated alpha amplitude 
alpha      =  [];
alpha2     =  [];

% generate fluctuated amplitude alpha signal(phase jittered)
for i = 1:num
T    = Time(i);
t    = 0:dt:T-dt;
sig  = (1+sin(1/T* 2*pi*t+1.5*pi));
alpha2 = [alpha2 sig];
sig   = amps(i)*sig;       % fluctuated alpha amplitude  
alpha = [alpha sig];
end

alpha      = 5* alpha;     % enhanced amplitude to be more realistic
T1         = 0:dt:length(alpha)*dt-dt;
T2         = 0:dt:(length(alpha)-timdiff)*dt-dt;

%  sigmoid threshold gamma 
gammas     = (1-1./(1 + exp(-a*(alpha-c)))).*(sin(2*pi*hf*T1)+1);
% gammas     = 1./(1 + exp(-a*(alpha-c))).*(sin(2*pi*hf*T1)+1);

% shift to creat directionality not mean for CFC 
switch shift                         
case {'lead'}        % alpha leads gamma 
alphastemp = zeros(1,length(alpha));
for i = 1:length(alpha)-timdiff
alphastemp(i) = alpha(i+timdiff);
end
alphas = alphastemp(1:end-timdiff);
gamma  = gammas(1:end-timdiff);
   
case{'delay'}        % alpha  lags gamma 
alphastemp = zeros(1,length(alpha));
for i = timdiff+1:length(alpha)
alphastemp(i) = alpha(i-timdiff);
end
 alphas(1:length(T2))= alphastemp(timdiff+1:length(alpha));
 gamma  = gammas(timdiff+1:end);
    otherwise         % no directionality 
 alphas = alpha(1:end-timdiff); 
 gamma  = gammas(1:end-timdiff);
end

% alphas      = 5* alphas;     % enhanced amplitude to be more realistic
alpha0  = alpha2(1:length(T2));
mix    =  alphas+gamma+n1*randn([1,length(T2)])+n2*pinknoise(numel(T2)); % signal we analysis
sigs = [alpha0;alphas;gamma;mix];
 
function x = pinknoise(Nx)
% pink noise 
B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
A = [1 -2.494956002   2.017265875  -0.522189400];
nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
x = filter(B,A,v);    % Apply 1/F roll-off to PSD
x = x(nT60+1:end);    % Skip transient response
end


end
