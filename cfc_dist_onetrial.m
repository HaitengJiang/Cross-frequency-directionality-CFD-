function CFC = cfc_dist_onetrial(x,para)

% compute  FFT of low frequency signal, high frequency envelope signal and
% cross-spectra for CFC/CFD : single trial mutilchannel 
% 
% input
% 
%  x :                   nchan*npoints signal
% para.Fs:               sample frequency 
% para.v :               Vector of interested  high frequency 
% para.xlim              range of interested low frequency range. e.g.[1
%                        60]; default is all. This tries to save space
% para.K:                length of window to estimate (in second)
% para.L:                sliding steps (in second)  
% para.Ncycles:           number of cycles to estimate high frequency power
% para.nFFT              nfft points windows of fft, which determines frequency
%                        resolution(in second)
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
% copyright @Haiteng Jiang; When you use this code ,please cite this paper
% Jiang, H., Bahramisharif, A, van Gerven, M.A.J. and Jensen, O (2015). 
% Measuring Directionality between Neuronal Oscillations of Different Frequencies.
% NeuroImage, 118:359-367

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
 

 
[nchan, npoints] = size(x);

% Define how segmented the trials into K time windows with the 'nshift'

nstart = 1:L:npoints-K+1;
nstop = nstart +K-1; 
S = length(nstart) ;     %  S segments



Xffts = zeros(nchan,length(v),nFFT,S);  
Yffts = zeros(nchan,length(v),nFFT,S);
CSs = zeros(nchan,length(v),nFFT,S);


for chani=1:nchan
    
sig = x(chani,:);

% The Hannings taper applied for low frequencies when calculating the phase
% consistency

HannLowFreq = hanning(K);

Xfft = zeros(length(v),nFFT,S);  
Yfft = zeros(length(v),nFFT,S);
CS = zeros(length(v),nFFT,S);

for vi=1:length(v)    
    
    % Calculate the temporal development of the high-frequnce power using a sliding 
    % time-window N points long after applying a Hanning taper. Do this
    % when looping over the frequencies on the y-axis. 
    
    f = v(vi);    
    N = floor(Ncycles*Fs/f);
    HannHighFreq = hanning(N)';
    
    % extract amplitude envelope of  high frequency  
    
    if isfield(para, 'detrend') && strcmp(para.detrend, 'yes')
        
    y = abs(conv(detrend(sig),HannHighFreq.*exp(i*2*pi*f.*(1:N)/Fs),'same'))/sumsqr(HannHighFreq);
    
    else
     
    y = abs(conv(sig,HannHighFreq.*exp(i*2*pi*f.*(1:N)/Fs),'same'))/sumsqr(HannHighFreq);
    
    end    
    
    
    % Loop over the segments of the trials. For each segment calculate the
    % phase difference in the complex domain, fft(sig) .* conj(fft(y)) and 
    % accumulate those values. 
    
    for s = 1:S
        
        if isfield(para, 'detrend') && strcmp(para.detrend, 'yes')
        
       sigs = detrend(sig(nstart(s):nstop(s)));
         ys  = detrend(y(nstart(s):nstop(s)));
         
        else            
        sigs = sig(nstart(s):nstop(s));
         ys   =  y(nstart(s):nstop(s));   
        end       

         
        Xfft(vi,:,s) =  fft(HannLowFreq'.*sigs,nFFT);
        Yfft(vi,:,s) =  fft(HannLowFreq'.*ys,nFFT);
        CS(vi,:,s)   = Xfft(vi,:,s).*conj(Yfft(vi,:,s));        
        
    end
  
   
end

Xffts(chani,:,:,:) = Xfft;  
Yffts(chani,:,:,:) = Yfft; 
CSs(chani,:,:,:) = CS; 

    
end


if isfield(para, 'xlim')     
     
  f = para.xlim;  %  interested frequency range 
 
 Vec_X = (Fs/nFFT)*(0:(nFFT-1)/2);
 [~,fminid] = min(abs(Vec_X-f(1)));
 [~,fmaxid] = min(abs(Vec_X-f(2)));

CFC.Xffts     = Xffts(:,:,fminid:fmaxid,:);  
CFC.Yffts     = Yffts(:,:,fminid:fmaxid,:); 
CFC.CSs       = CSs(:,:,fminid:fmaxid,:);
CFC.freqVec_X = Vec_X(fminid:fmaxid);

else
  
CFC.Xffts    = Xffts(:,:,2:nFFT/2,:);   % exclude DC
CFC.Yffts    = Yffts(:,:,2:nFFT/2,:); 
CFC.CSs      = CSs(:,:,2:nFFT/2,:);
CFC.freqVec_X = (Fs/nFFT)*(1:(nFFT-1)/2);

end



CFC.freqVec_Y = v;
CFC.Fs        = Fs;
CFC.dfreq     = Fs/nFFT;

end





