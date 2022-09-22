%% test simulaiton in different conditons 

% copyright @Haiteng Jiang; When you use this code ,please cite this paper
% Jiang, H., Bahramisharif, A, van Gerven, M.A.J. and Jensen, O (2015). 
% Measuring Directionality between Neuronal Oscillations of Different Frequencies.
% NeuroImage, 118:359-367

clear all;
close all;
clc;

% generate simulation data 
% (3 channels): 
                          %chan1: alpha lags gamma; 
                          %chan2: alpha leads gamma;
                          %chan3:  no lags ;

nchanX =  3;        % nr of chann
M     =  500000;   % M points data
X  =  zeros(nchanX,M);

 for i = 1:nchanX

para               = [];
para.num           = 6000;                 % number of alpha cycles
para.Fs            = 1000;                 % sampling frequency 
para.hf            = 70;                   % gamma frequency 
if mod(i,3)  ==1
para.shift         = 'delay';              % for directionality only: 'lead' or 'delay'
elseif mod(i,3)  ==2
para.shift         = 'lead';
else
para.shift         = 'no lag';
end
para.timdiff       = 10;                   % time lag 
para.a             = 10;    
para.c             = 6;                    % a(slope)/c(threshold)  :  Sigmoid function parameter sigmf(x, [a, c]) = 1./(1 + exp(-a*(x-c)))
para.n1            = 0.8;                    % Gaussian white noise level 
para.n2            = 0.2;                 % pink noise level 

[sig,T]= inhibition(para);
X(i,:) = sig(4,1:M);

 end

 
 para          = [];
para.Fs       = 1000;             
para.v        = 40:5:150;
para.K        = 2; 
para.L        = 1;           
para.Ncycles   = 5;          
para.nFFT     = 2;       
para.detrend  = 'yes';
CFCdist_X     = cfc_dist_onetrial(X,para);


para          = [];
para.method    = 'coh'; 
CFCs_X = cfc_quantification(CFCdist_X,[],para);

para            =[];
para.freqVec_x   =  5:30;
para.fwidth      =  2 ;       
para.jack        = 'no';
CFDs_X = cfd_quantification(CFCdist_X,[],para);


CFD_X1=squeeze(CFDs_X.CFD(1,:,:));
M1   = max(abs(CFD_X1(:)));

CFD_X2=squeeze(CFDs_X.CFD(2,:,:));
M2   = max(abs(CFD_X2(:)));

CFD_X3=squeeze(CFDs_X.CFD(3,:,:));
M3   = max(abs(CFD_X3(:)));


M =max([M1 M2 M3]);


%%
figure();
subplot(2,3,1)
imagesc(CFCs_X.freqVec_X,CFCs_X.freqVec_Y,squeeze(CFCs_X.CFC(1,:,:)));axis xy;colorbar;
xlabel('phase frequency (Hz)  ','fontsize',10)
ylabel('amplitude frequency(Hz) ','fontsize',10)
xlim([5  30])
title('alpha lags gamma CFC')

subplot(2,3,2)
imagesc(CFCs_X.freqVec_X,CFCs_X.freqVec_Y,squeeze(CFCs_X.CFC(2,:,:)));axis xy;colorbar;
xlabel('phase frequency (Hz)','fontsize',10)
ylabel('amplitude frequency(Hz)','fontsize',10)
xlim([5  30])
title('alpha leads gamma CFC')

subplot(2,3,3)
imagesc(CFCs_X.freqVec_X,CFCs_X.freqVec_Y,squeeze(CFCs_X.CFC(3,:,:)));axis xy;colorbar;
xlabel('phase frequency (Hz)','fontsize',10)
ylabel('amplitude frequency(Hz)','fontsize',10)
xlim([5  30])
title('no lags CFC ')


subplot(2,3,4)
imagesc(CFDs_X.LF,CFDs_X.HF,CFD_X1,[-M M]);axis xy;colorbar
xlabel('phase frequency (Hz)','fontsize',10)
ylabel('amplitude frequency(Hz)','fontsize',10)
xlim([5  30])
title('alpha lags gamma CFD')

subplot(2,3,5)
imagesc(CFDs_X.LF,CFDs_X.HF,CFD_X2,[-M M]);axis xy;colorbar
xlabel('phase frequency (Hz)','fontsize',10)
ylabel('amplitude frequency(Hz)','fontsize',10)
xlim([5  30])
title('alpha leads gamma CFD')

subplot(2,3,6)
imagesc(CFDs_X.LF,CFDs_X.HF,CFD_X3,[-M M]);axis xy;colorbar
xlabel('phase frequency (Hz)','fontsize',10)
ylabel('amplitude frequency(Hz)','fontsize',10)
xlim([5  30])
title('no lags CFD')

