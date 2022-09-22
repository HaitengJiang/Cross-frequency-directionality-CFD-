function CFC = cfc_dist_multitrial(x,para)

% compute  FFT of low frequency signal, high frequency envelope signal and
% cross-spectra for CFC/CFD : multitrial mutilchannel 
% 
% input 
% 
%  x :                   ntrial*nchan*npoints signal
% para.Fs:               sample frequency 
% para.v :               Vector of interested  high frequency
% para.xlim              range of interested low frequency range. e.g.[1
%                        60]; default is all. This tries to save space
% para.K:                length of window to estimate in second
% para.L:                sliding steps in second  
% para.Ncycles:           number of cycles to estimate high frequency power
% para.nFFT              nfft points windows of fft, which determines frequency
%                        resolution; in second 
% para.detrend           detrend or not 'yes' or 'no'        
% 
% 
% Output:
% 
% CFC.Xffts           FFT of low frequency signal
% CFC.Yffts           FFT of high frequency power 
% CFC.CSs             cross-spectra between low frequency signal and high frequency power 
% CFC.freqVec_X       phase frequency bin 
% CFC.freqVec_Y       amplitude frequency bin 
% CFC.Fs              sample frequency
% CFC.dfreq           frequency resoulation
% 
% copyright @Haiteng Jiang; When you use this code ,please cite this paper
% Jiang, H., Bahramisharif, A, van Gerven, M.A.J. and Jensen, O (2015). 
% Measuring Directionality between Neuronal Oscillations of Different Frequencies.
% NeuroImage, 118:359-367

CFC = [];
if isfield(para, 'Fs')     
     Fs = para.Fs; 
end
 
 if isfield(para, 'v')     
     v = para.v; 
 end
 
 if isfield(para, 'K')     
     K = floor(Fs*para.K);  
 end
 
 if isfield(para, 'L')     
     L = floor(Fs*para.L); 
 end 
 
 if isfield(para, 'Ncycles')     
     Ncycles = para.Ncycles; 
 end
 
 if isfield(para, 'nFFT')     
    nFFT = 2^nextpow2(Fs*para.nFFT);  % 2^nextpow2 to make FFT faster
 end
 

 
[ntrial,nchan, npoints] = size(x);

S = length(1:L:npoints-K+1) ;     %  S segments within one trial

if isfield(para, 'xlim')     
     
  f = para.xlim;          %  interested frequency range 
 
  Vec_X = (Fs/nFFT)*(0:(nFFT-1)/2);
 [~,fminid] = min(abs(Vec_X-f(1)));
 [~,fmaxid] = min(abs(Vec_X-f(2)));
                
  Xffts       = zeros(nchan,numel(v),fmaxid-fminid+1,S*ntrial);
  Yffts       = zeros(nchan,numel(v),fmaxid-fminid+1,S*ntrial);
  CSs         = zeros(nchan,numel(v),fmaxid-fminid+1,S*ntrial);

else
   
  Xffts       = zeros(nchan,numel(v),floor(nFFT/2)-1,S*ntrial);
  Yffts       = zeros(nchan,numel(v),floor(nFFT/2)-1,S*ntrial);
  CSs         = zeros(nchan,numel(v),floor(nFFT/2)-1,S*ntrial);

    
end
   

for k=1:ntrial
    if nchan==1        % to support one chan in the loop 
    sig =squeeze(x(k,:,:))';
    else
    sig = squeeze(x(k,:,:));
    end
    CFCtemp = cfc_dist_onetrial(sig,para);
    Xffts(:,:,:,(k-1)*S+1:k*S)        = CFCtemp.Xffts;
    Yffts(:,:,:,(k-1)*S+1:k*S)         = CFCtemp.Yffts;
    CSs(:,:,:,(k-1)*S+1:k*S)          = CFCtemp.CSs;
end

CFC.Xffts  =Xffts; 
CFC.Yffts = Yffts;
CFC.CSs   = CSs;
if isfield(para, 'xlim')   
CFC.freqVec_X = Vec_X(fminid:fmaxid);
else
CFC.freqVec_X = (Fs/nFFT)*(1:(nFFT-1)/2);
end
CFC.freqVec_Y = v;
CFC.Fs        = Fs;
CFC.dfreq     = Fs/nFFT;



