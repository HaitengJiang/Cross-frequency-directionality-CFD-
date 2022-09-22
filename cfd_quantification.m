function CFDs = cfd_quantification(dist1,dist2,para)

% calculate within region CFD (dist2 = []) and cross-region
% CFD; Input from  cfc_dist_onetrial/multitrial
% 
% input:
% 
% dist1   input from cfc_dist_onetrial/multitrial signal 1
% dist2   input from cfc_dist_onetrial/multitrial signal 2
% para.freqVec_x:      low frequency range to calculate phase slope index
% para.fwidth:         frequency bandwidth to estimate phase slope(in Hz) e.g., [f-fwidth/2,f+fwidth/2]
% para.jack:           'yes': Jackknife(time comsuming,more robust) 'no': no jackknife
% 
% output:
% 
% %%%% within region CFC: cfc_quantification(dist1,[],para)%%%%%%
% 
% CFDs.CFD          cross freuqency directionality map
% CFDs.method       method 'plv','coh','kl'
% CFDs.freqVec_X     phase frequency bin for CFD 
% CFDs.freqVec_Y     amplitude frequency bin for CFD
% 
% %%%% cross region CFD: cfc_quantification(dist1,dist2,para)%%%%%%
% 
% CFDs.CFD1          phase frequency from dist1 and amplitude frequency
%                     from dist2
% CFDs.CFD2          phase frequency from dist2 and amplitude frequency
%                     from dist1
% CFDs.LF            phase frequency bin for CFD 
% CFDs.HF            amplitude frequency bin for CFD
% 
% copyright @Haiteng Jiang; When you use this code ,please cite this paper
% Jiang, H., Bahramisharif, A, van Gerven, M.A.J. and Jensen, O (2015). 
% Measuring Directionality between Neuronal Oscillations of Different Frequencies.
% NeuroImage, 118:359-367


CFDs = [];

if isfield(para, 'freqVec_x')     
     LF = para.freqVec_x; 
end

if isfield(para, 'fwidth')     
    fwidth = para.fwidth; 
end

if isfield(para, 'jack') && strcmp(para.jack, 'yes')
 jack = 1;
else    
  jack = 0;        
end

HF     = dist1.freqVec_Y;
dealtf = dist1.dfreq;

if isempty(dist2)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % within-region CFD %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('within-region CFD') 
[nchan,nv,nfs,~] = size(dist1.CSs);
CFD =  zeros(nchan,numel(HF),numel(LF));


for nchani = 1:nchan
  CFDtempCS      = squeeze(dist1.CSs(nchani,:,:,:));
  CFDtempYffts   = squeeze(dist1.Yffts(nchani,:,:,:)); 
  CFDtempXffts   = squeeze(dist1.Xffts(nchani,:,:,:));
  CFD(nchani,:,:) = dist2ps(CFDtempXffts,CFDtempYffts,CFDtempCS,jack);
  
end


CFDs.CFD  = CFD;
CFDs.LF  = LF;
CFDs.HF  = HF;


else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % cross-region CFD %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
disp('cross-region CFD') 

[nchan1,~,~,~] = size(dist1.CSs);
[nchan2,~,~,~] = size(dist2.CSs);  
      
  % conditon1  phase freq from dist1 amplitude freq from dist2   
   CFD1 = zeros(nchan1,nchan2,numel(HF),numel(LF));   
  for nchan1i = 1:nchan1   % phase freq from dist1      
  for nchan2i = 1:nchan2   % amplitude freq from dist2
      
   CFDtempX1   = squeeze(dist1.Xffts(nchan1i,:,:,:));
   CFDtempY2   = squeeze(dist2.Yffts(nchan2i,:,:,:)); 
   CStemp12    = CFDtempX1.*conj(CFDtempY2);    
   CFD1(nchan1i,nchan2i,:,:) = dist2ps(CFDtempX1,CFDtempY2,CStemp12 ,jack);

  end
  end
  
  % conditon2   phase freq from dist2, amplitude freq from dist1
  if(isequaln(dist1,dist2)==0)     % dist1~=dist2
   CFD2 = zeros(nchan2,nchan1,numel(HF),numel(LF));   
  for nchan2i = 1:nchan2   % phase freq from dist2      
  for nchan1i = 1:nchan1   % amplitude freq from dist1
      
   CFDtempX2   = squeeze(dist2.Xffts(nchan2i,:,:,:));
   CFDtempY1   = squeeze(dist1.Yffts(nchan1i,:,:,:)); 
   CStemp21    = CFDtempX2.*conj(CFDtempY1);
   CFD2(nchan2i,nchan1i,:,:) = dist2ps(CFDtempX2,CFDtempY1,CStemp21,jack);
  
  
  end
  end  
  end

CFDs.CFD1  = CFD1;
if(isequaln(dist1,dist2)==0)     % dist1~=dist2
CFDs.CFD2  = CFD2;  
end
CFDs.LF  = LF;
CFDs.HF  = HF;
    
end

function CFD = dist2ps(CFDtempXffts,CFDtempYffts,CFDtempCS,jack)


nep = size(CFDtempCS,3);            % n epochs

% whether to use Jackknife to calculate phase slope deviation (time comsuming)
if jack ==0           %  no  jackknife 
 CS   = mean(CFDtempCS,3)./sqrt((mean(abs(CFDtempXffts).^2,3) .* mean(abs(CFDtempYffts).^2,3))); 
 CFD   = cs2ps(CS,LF,HF,fwidth,dealtf); 
else                   %  jackknife 
jackid = 1:nep;
CFDjack = zeros(numel(HF),numel(LF),nep);
for i = 1:nep
 jacktemp = jackid;
 jacktemp(i) = [];
 Xtemp   =  CFDtempXffts(:,:,jacktemp);
 Ytemp  = CFDtempYffts(:,:,jacktemp);
 Ctemp    = CFDtempCS(:,:,jacktemp);
 CSjack   = mean(Ctemp,3)./sqrt((mean(abs(Xtemp).^2,3) .* mean(abs(Ytemp).^2,3))); 
 CFDjack(:,:,i) = cs2ps(CSjack,LF,HF,fwidth,dealtf);   
end 
CFD  = squeeze(std(CFDjack,0,3))*sqrt(nep);

end


end



function CFD = cs2ps(CS,LF,HF,fwidth,dealtf)
% complex coherency to phase slope

CFD   = zeros(numel(HF),numel(LF));
 
for k = 1:numel(LF)
    maxfreqbin = floor(LF(k)/dealtf);         % fequency range to estimate phase slope
    step       = round(0.5*fwidth/dealtf);    % bandwidth to calculate phase slope    
    for j =1:numel(HF)        
    CStemp = squeeze(CS(j,:));
    CFDtemp = imag(conj(CStemp(1:end-1)).*CStemp(2:end));
    CFD(j,k)= mean(CFDtemp(maxfreqbin-step:maxfreqbin+step));   
    end
end
     
 
end    


end




