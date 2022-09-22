% function test_cfc_singletrial: test the cfc for a single trial


clear all;
close all;
clc;

% generate simulation data 
%%  region  X  (3 channels): 
                          %chan1: alpha lags gamma; 
                          %chan2: alpha leads gamma;
                          %chan3: alpha lags gamma;

nchanX =  3;        % nr of chann
M     =  200000;   % M points data
X  =  zeros(nchanX,M);

 for i = 1:nchanX

para               = [];
para.num           = 3000;                 % number of alpha cycles
para.Fs            = 1000;                 % sampling frequency 
para.hf            = 60;                   % gamma frequency 
if mod(i,2)  ==1
para.shift         = 'delay';              % for directionality only: 'lead' or 'delay'
else 
para.shift         = 'lead';
end
para.timdiff       = 10;                   % time lag 
para.a             = 10;    
para.c             = 6;                    % a(slope)/c(threshold)  :  Sigmoid function parameter sigmf(x, [a, c]) = 1./(1 + exp(-a*(x-c)))
para.n1            = 0.5;                    % Gaussian white noise level 
para.n2            = 0.1;                 % pink noise level 

[sig,T]= inhibition(para);
X(i,:) = sig(4,1:M);

end


para          = [];
para.Fs       = 1000;             
para.v        = 30:5:150;
para.xlim     = [1 60];
para.K        = 2; 
para.L        = 1;           
para.Ncycles   = 5;          
para.nFFT     = 4;       
para.detrend  = 'yes';
CFCdist_X = cfc_dist_onetrial(X,para);


%%  region  Y  (3 channels): 
                          %chan1: alpha lags gamma; 
                          %chan2: alpha leads gamma;
                           %chan3: alpha lags gamma;
nchanY =  3;        % nr of chann
M     =  200000;   % M points data
Y  =  zeros(nchanY,M);

 for i = 1:nchanX

para               = [];
para.num           = 3000;                 % number of alpha cycles
para.Fs            = 1000;                 % sampling frequency 
para.hf            = 60;                   % gamma frequency 
if mod(i,2)  ==1
para.shift         = 'delay';              % for directionality only: 'lead' or 'delay'
else 
para.shift         = 'lead';
end
para.timdiff       = 10;                   % time lag 
para.a             = 10;    
para.c             = 6;                    % a(slope)/c(threshold)  :  Sigmoid function parameter sigmf(x, [a, c]) = 1./(1 + exp(-a*(x-c)))
para.n1            = 0.5;                    % Gaussian white noise level 
para.n2            = 0.1;                 % pink noise level 

[sig,T]= inhibition(para);
Y(i,:) = sig(4,1:M);

end


para          = [];
para.Fs       = 1000;             
para.v        = 30:5:150;
para.xlim     = [1 60];
para.K        = 2; 
para.L        = 1;           
para.Ncycles   = 5;          
para.nFFT     = 4;       
para.detrend  = 'yes';
CFCdist_Y = cfc_dist_onetrial(Y,para);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % within-region CFC %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para          = [];
para.method    = 'coh'; 
CFCs_Y = cfc_quantification(CFCdist_Y,[],para);

figure(1);
subplot(2,1,1)
imagesc(CFCs_Y.freqVec_X,CFCs_Y.freqVec_Y,squeeze(CFCs_Y.CFC(1,:,:)));axis xy;colorbar
xlabel('phase frequency (Hz) from Region Y chan 1 ','fontsize',10)
ylabel('amplitude frequency(Hz) from Region Y chan 1','fontsize',10)
title(['within-region CFC: Region Y chan 1 ',CFCs_Y.method],'fontsize',10)
xlim([4  40])

subplot(2,1,2)
imagesc(CFCs_Y.freqVec_X,CFCs_Y.freqVec_Y,squeeze(CFCs_Y.CFC(2,:,:)));axis xy;colorbar
xlabel('phase frequency (Hz) from Region Y chan 2','fontsize',10)
ylabel('amplitude frequency(Hz) from Region Y chan 2','fontsize',10)
title(['within-region CFC: Region Y chan 2 ',CFCs_Y.method],'fontsize',10)
xlim([4  40])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % within-region CFD %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para            =[];
para.freqVec_x   =  4:40;
para.fwidth      =  2 ;       
para.jack        = 'no';

CFDs_Y = cfd_quantification(CFCdist_Y,[],para);


figure(2);
subplot(2,1,1)
CFD_Y1=squeeze(CFDs_Y.CFD(1,:,:));
M1   = max(abs(CFD_Y1(:)));
imagesc(CFDs_Y.LF,CFDs_Y.HF,CFD_Y1,[-M1 M1]);axis xy;colorbar
xlabel('phase frequency (Hz) from Region Y chan 1 ','fontsize',10)
ylabel('amplitude frequency(Hz) from Region Y chan 1','fontsize',10)
title('within-region CFD: Region Y chan 1','fontsize',10)
xlim([4  40])

subplot(2,1,2)
CFD_Y2=squeeze(CFDs_Y.CFD(2,:,:));
M2   = max(abs(CFD_Y2(:)));
imagesc(CFDs_Y.LF,CFDs_Y.HF,CFD_Y2,[-M2 M2]);axis xy;colorbar
xlabel('phase frequency (Hz) from Region Y chan 2 ','fontsize',10)
ylabel('amplitude frequency(Hz) from Region Y chan 2','fontsize',10)
title('within-region CFD: Region Y chan 2','fontsize',10)
xlim([4  40])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % cross-region CFC %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para          = [];
para.method    = 'coh'; 
CFCs_XY = cfc_quantification(CFCdist_X,CFCdist_Y,para);  % no cross-region CFC between X and Y

%%%%%%%%%%%%% plot cross-region CFC between X and Y: phase from X and amplitude from Y:CFCs_XY.CFC1 %%%%

figure(3);
M = 0.1;
for i= 1:nchanX  
subplot(nchanX,nchanY,(i-1)*nchanY+1)
imagesc(CFCs_XY.freqVec_X,CFCs_XY.freqVec_Y,squeeze(CFCs_XY.CFC1(i,1,:,:)),[0 M]);axis xy;colorbar
xlabel(['phase frequency from Region X chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 1')
xlim([4  40])

subplot(nchanX,nchanY,(i-1)*nchanY+2)
imagesc(CFCs_XY.freqVec_X,CFCs_XY.freqVec_Y,squeeze(CFCs_XY.CFC1(i,2,:,:)),[0 M]);axis xy;colorbar
xlabel(['phase frequency from Region X chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 2')
xlim([4  40])

end


%%%%%%%%%%%%%plot  cross-region CFC between Y and X: phase from Y and amplitude from X: CFCs_XY.CFC2 %%%%
figure(4);
for i= 1:nchanY  
subplot(nchanY,nchanX,(i-1)*nchanX+1)
imagesc(CFCs_XY.freqVec_X,CFCs_XY.freqVec_Y,squeeze(CFCs_XY.CFC2(i,1,:,:)),[0 M]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region X chan 1')
xlim([4  40])

subplot(nchanY,nchanX,(i-1)*nchanX+2)
imagesc(CFCs_XY.freqVec_X,CFCs_XY.freqVec_Y,squeeze(CFCs_XY.CFC2(i,2,:,:)),[0 M]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region X chan 2')
xlim([4  40])

subplot(nchanY,nchanX,(i-1)*nchanX+3)
imagesc(CFCs_XY.freqVec_X,CFCs_XY.freqVec_Y,squeeze(CFCs_XY.CFC2(i,3,:,:)),[0 M]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region X chan 3')
xlim([4  40])

end

para          = [];
para.method    = 'coh'; 
CFCs_YY = cfc_quantification(CFCdist_Y,CFCdist_Y,para);  % cross-region CFC: cross channels in region Y : 

%%%%%%%%%%%%%plot  cross-region CFC between Y and Y: phase from Y and amplitude from Y (cross channels) %%%%

figure(5);
for i= 1:nchanY  
subplot(nchanY,nchanY,(i-1)*nchanY+1)
imagesc(CFCs_YY.freqVec_X,CFCs_YY.freqVec_Y,squeeze(CFCs_YY.CFC1(i,1,:,:)));axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 1')
xlim([4  40])

subplot(nchanY,nchanY,(i-1)*nchanY+2)
imagesc(CFCs_YY.freqVec_X,CFCs_YY.freqVec_Y,squeeze(CFCs_YY.CFC1(i,2,:,:)));axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 2')
xlim([4  40])

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % cross-region CFD %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


para            =[];
para.freqVec_x   =  4:40;
para.fwidth      =  2 ;       
para.jack        = 'no';

CFDs_XY = cfd_quantification(CFCdist_X,CFCdist_Y,para); % cross-region CFD between Y and X


%%%%%%%%%%%%%plot cross-region CFD between Y and X: phase from X and amplitude from Y: CFDs_XY.CFD1 %%%%
figure(6)
for i= 1:nchanX  
subplot(nchanX,nchanY,(i-1)*nchanY+1)
Mi1 = squeeze(CFDs_XY.CFD1(i,1,:,:));
M1   = max(abs(Mi1(:)));
imagesc(CFDs_XY.LF,CFDs_XY.HF,Mi1,[-M1 M1]);axis xy;colorbar
xlabel(['phase frequency from Region X chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 1')
xlim([4  40])
subplot(nchanX,nchanY,(i-1)*nchanY+2)
Mi2 = squeeze(CFDs_XY.CFD1(i,2,:,:));
M2   = max(abs(Mi2(:)));
imagesc(CFDs_XY.LF,CFDs_XY.HF,Mi2,[-M2 M2]);axis xy;colorbar
xlabel(['phase frequency from Region X chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 2')
xlim([4  40])
end

%%%%%%%%%%%%%plot cross-region CFD between Y and X: phase from Y and amplitude from X: CFDs_XY.CFD2 %%%%

figure(7)
for i= 1:nchanY  
subplot(nchanY,nchanX,(i-1)*nchanX+1)
Mi1 = squeeze(CFDs_XY.CFD2(i,1,:,:));
M1   = max(abs(Mi1(:)));
imagesc(CFDs_XY.LF,CFDs_XY.HF,Mi1,[-M1 M1]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region X chan 1')
xlim([4  40])

subplot(nchanY,nchanX,(i-1)*nchanX+2)
Mi2 = squeeze(CFDs_XY.CFD2(i,2,:,:));
M2   = max(abs(Mi2(:)));
imagesc(CFDs_XY.LF,CFDs_XY.HF,Mi2,[-M2 M2]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region X chan 2')
xlim([4  40])

subplot(nchanY,nchanX,(i-1)*nchanX+3)
Mi3 = squeeze(CFDs_XY.CFD2(i,3,:,:));
M3   = max(abs(Mi3(:)));
imagesc(CFDs_XY.LF,CFDs_XY.HF,Mi3,[-M3 M3]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region X chan 2')
xlim([4  40])

end


% cross-region CFD in region Y: cross channels

para            =[];
para.freqVec_x   =  4:40;
para.fwidth      =  2 ;       
para.jack        = 'no';

CFDs_YY = cfd_quantification(CFCdist_Y,CFCdist_Y,para); 

%%%%%%%%%%%%%plot  cross-region CFD between Y and Y: phase from Y and amplitude from Y (cross channels) %%%%
figure(8)
for i= 1:nchanY 
subplot(nchanY,nchanY,(i-1)*nchanY+1)
Mi1 = squeeze(CFDs_YY.CFD1(i,1,:,:));
M1   = max(abs(Mi1(:)));
imagesc(CFDs_YY.LF,CFDs_YY.HF,Mi1,[-M1 M1]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 1')
xlim([4  40])
subplot(nchanY,nchanY,(i-1)*nchanY+2)
Mi2 = squeeze(CFDs_YY.CFD1(i,2,:,:));
M2   = max(abs(Mi2(:)));
imagesc(CFDs_YY.LF,CFDs_YY.HF,Mi2,[-M2 M2]);axis xy;colorbar
xlabel(['phase frequency from Region Y chan ',num2str(i)])
ylabel('amplitude frequency from Region Y chan 2')
xlim([4  40])
end



