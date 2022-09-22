function CFCs = cfc_quantification(dist1,dist2,para)

% calculate within region CFC (dist2 = []) and cross-region
% CFC; supporting method : phase locking value (plv), cohernce (coh)
% imaginary part of coherence (icoh) and  KL divergence (kl). Input from  cfc_dist_onetrial/multitrial
% 
% input:
% 
% dist1   input from cfc_dist_onetrial/multitrial signal 1
% dist2   input from cfc_dist_onetrial/multitrial signal 2
% para.method:  plv:  phase locking value
%               coh:  coherence
%               ioch: imaginary part of coherence
%               kl:   KL divergence 
%               
% output:
% 
%  %%%%% within region CFC: cfc_quantification(dist1,[],para)%%%%%%
% 
% CFCs.CFC          cross freuqency coupling map
% CFCs.method       method 'plv','coh','icoh','kl'
% CFCs.freqVec_X     phase frequency bin for CFC 
% CFCs.freqVec_Y     amplitude frequency bin for CFD
% 
%  %%%%% cross region CFC: cfc_quantification(dist1,dist2,para)%%%%%%
% 
% CFCs.CFC1          phase frequency from dist1 and amplitude frequency
%                     from dist2
% CFCs.CFC2          phase frequency from dist2 and amplitude frequency
%                     from dist1
% CFCs.method       method 'plv','coh','icoh','kl'
% CFCs.freqVec_X     phase frequency bin for CFC 
% CFCs.freqVec_Y     amplitude frequency bin for CFC
% 
% copyright @Haiteng Jiang; When you use this code ,please cite this paper
% Jiang, H., Bahramisharif, A, van Gerven, M.A.J. and Jensen, O (2015). 
% Measuring Directionality between Neuronal Oscillations of Different Frequencies.
% NeuroImage, 118:359-367


CFCs    =[];

if isempty(dist2)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % within-region CFC %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('within-region CFC') 

switch para.method
 
  case 'plv'
      
  [nchan,nv,nfs,~] = size(dist1.CSs);
  CFC =  zeros(nchan,nv,nfs);
  
 for chani = 1:nchan
   CFCtempCS = squeeze(dist1.CSs(chani,:,:,:));
   CFC(chani,:,:) = abs(mean(CFCtempCS./ abs(CFCtempCS),3));
      
 end 
 
CFCs.CFC   =  CFC;
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y; 
  
    
    
 case 'coh'      
  
  [nchan,nv,nfs,~] = size(dist1.CSs);
  CFC =  zeros(nchan,nv,nfs);
  
 for chani = 1:nchan
   CFCtempCS      = squeeze(dist1.CSs(chani,:,:,:));
   CFCtempXffts   = squeeze(dist1.Xffts(chani,:,:,:));
   CFCtempYffts   = squeeze(dist1.Yffts(chani,:,:,:));   
   CFC(chani,:,:) = abs(mean(CFCtempCS,3)).^2 ./(mean(abs(CFCtempXffts).^2,3) .* mean(abs(CFCtempYffts).^2,3));

end 
 
CFCs.CFC = CFC;
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y;

case 'icoh'      
  
  [nchan,nv,nfs,~] = size(dist1.CSs);
  CFC =  zeros(nchan,nv,nfs);
  
 for chani = 1:nchan
   CFCtempCS      = squeeze(dist1.CSs(chani,:,:,:));
   CFCtempXffts   = squeeze(dist1.Xffts(chani,:,:,:));
   CFCtempYffts   = squeeze(dist1.Yffts(chani,:,:,:));   
   CFC(chani,:,:) = imag(mean(CFCtempCS,3)).^2 ./(mean(abs(CFCtempXffts).^2,3) .* mean(abs(CFCtempYffts).^2,3));

end 
 
CFCs.CFC = CFC;
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y;
       
    
case 'kl'
 
 [nchan,nv,nfs,~] = size(dist1.CSs);
  CFC =  zeros(nchan,nv,nfs);
  
 for chani = 1:nchan
   CFCtempCS      = squeeze(dist1.CSs(chani,:,:,:));
   CFCtempYffts   = squeeze(dist1.Yffts(chani,:,:,:));  
   
 % Calculate the KL divergence 
   CFC(chani,:,:) = kl_cal(CFCtempCS,CFCtempYffts) ;  

end 
 
CFCs.CFC = CFC;
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y;
      
end
   
    
    
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % cross-region CFC %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
disp('cross-region CFC') 

[nchan1,nv,nfs,~] = size(dist1.CSs);
[nchan2,nv,nfs,~] = size(dist2.CSs);


switch para.method  
  
 
  case 'plv' 
      
  % conditon1 
   CFC1 = zeros(nchan1,nchan2,nv,nfs);   
  for nchan1i = 1:nchan1   % phase freq from dist1      
  for nchan2i = 1:nchan2   % amplitude freq from dist2
      
   CFCtempX1   = squeeze(dist1.Xffts(nchan1i,:,:,:));
   CFCtempY2   = squeeze(dist2.Yffts(nchan2i,:,:,:)); 
   CStemp12    = CFCtempX1.*conj(CFCtempY2);
   CFC1(nchan1i,nchan2i,:,:) = abs(mean(CStemp12./ abs(CStemp12),3));
  
  end
  end
  
  % conditon2 
    if(isequaln(dist1,dist2)==0)    % dist1~=dist2
   CFC2 = zeros(nchan2,nchan1,nv,nfs);    
  for nchan2i = 1:nchan2   % phase freq from dist2      
  for nchan1i = 1:nchan1   % amplitude freq from dist1
      
   CFCtempX2   = squeeze(dist2.Xffts(nchan2i,:,:,:));
   CFCtempY1   = squeeze(dist1.Yffts(nchan1i,:,:,:)); 
   CStemp21    = CFCtempX2.*conj(CFCtempY1);
   CFC2(nchan2i,nchan1i,:,:) = abs(mean(CStemp21./ abs(CStemp21),3));
  
  end
  end
  end
 
CFCs.CFC1   =  CFC1;
  if(isequaln(dist1,dist2)==0)  % dist1~=dist2
CFCs.CFC2   =  CFC2;
 end
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y; 

case 'coh'   
  
 
      
  % conditon1  phase freq from dist1 amplitude freq from dist2   
   CFC1 = zeros(nchan1,nchan2,nv,nfs);   
  for nchan1i = 1:nchan1   % phase freq from dist1      
  for nchan2i = 1:nchan2   % amplitude freq from dist2
      
   CFCtempX1   = squeeze(dist1.Xffts(nchan1i,:,:,:));
   CFCtempY2   = squeeze(dist2.Yffts(nchan2i,:,:,:)); 
   CStemp12    = CFCtempX1.*conj(CFCtempY2);   
   CFC1(nchan1i,nchan2i,:,:) = abs(mean(CStemp12,3)).^2 ./(mean(abs(CFCtempX1).^2,3) .* mean(abs(CFCtempY2).^2,3));

  end
  end
  
  % conditon2   phase freq from dist2, amplitude freq from dist1
  
   if(isequaln(dist1,dist2)==0)  % dist1~=dist2
  CFC2 = zeros(nchan2,nchan1,nv,nfs);   
  for nchan2i = 1:nchan2   % phase freq from dist2      
  for nchan1i = 1:nchan1   % amplitude freq from dist1
      
   CFCtempX2   = squeeze(dist2.Xffts(nchan2i,:,:,:));
   CFCtempY1   = squeeze(dist1.Yffts(nchan1i,:,:,:)); 
   CStemp21    = CFCtempX2.*conj(CFCtempY1);
   CFC2(nchan2i,nchan1i,:,:) = abs(mean(CStemp21,3)).^2 ./(mean(abs(CFCtempX2).^2,3) .* mean(abs(CFCtempY1).^2,3));

  
  end
  end
  end
 
CFCs.CFC1   =  CFC1;
 if(isequaln(dist1,dist2)==0)  % dist1~=dist2
CFCs.CFC2   =  CFC2;
end
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y; 

case 'icoh'
    
  
  
      
  % conditon1  phase freq from dist1 amplitude freq from dist2   
  CFC1 = zeros(nchan1,nchan2,nv,nfs);    
  for nchan1i = 1:nchan1   % phase freq from dist1      
  for nchan2i = 1:nchan2   % amplitude freq from dist2
      
   CFCtempX1   = squeeze(dist1.Xffts(nchan1i,:,:,:));
   CFCtempY2   = squeeze(dist2.Yffts(nchan2i,:,:,:)); 
   CStemp12    = CFCtempX1.*conj(CFCtempY2);   
   CFC1(nchan1i,nchan2i,:,:) = imag(mean(CStemp12,3)).^2 ./(mean(abs(CFCtempX1).^2,3) .* mean(abs(CFCtempY2).^2,3));

  end
  end
  
  % conditon2   phase freq from dist2, amplitude freq from dist1
   if(isequaln(dist1,dist2)==0)  % dist1~=dist2
   CFC2 = zeros(nchan2,nchan1,nv,nfs);     
  for nchan2i = 1:nchan2   % phase freq from dist2      
  for nchan1i = 1:nchan1   % amplitude freq from dist1
      
   CFCtempX2   = squeeze(dist2.Xffts(nchan2i,:,:,:));
   CFCtempY1   = squeeze(dist1.Yffts(nchan1i,:,:,:)); 
   CStemp21    = CFCtempX2.*conj(CFCtempY1);
   CFC2(nchan2i,nchan1i,:,:) = imag(mean(CStemp21,3)).^2 ./(mean(abs(CFCtempX2).^2,3) .* mean(abs(CFCtempY1).^2,3));

  
  end
  end
  
  end
 
CFCs.CFC1   =  CFC1;
 if(isequaln(dist1,dist2)==0)  % dist1~=dist2
CFCs.CFC2   =  CFC2;
end
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y; 

case  'kl'    
    
   
  
      
  % conditon1  phase freq from dist1 amplitude freq from dist2   
   CFC1 = zeros(nchan1,nchan2,nv,nfs);  
  for nchan1i = 1:nchan1   % phase freq from dist1      
  for nchan2i = 1:nchan2   % amplitude freq from dist2
      
   CFCtempX1   = squeeze(dist1.Xffts(nchan1i,:,:,:));
   CFCtempY2   = squeeze(dist2.Yffts(nchan2i,:,:,:)); 
   CStemp12    = CFCtempX1.*conj(CFCtempY2); 
   
    % Calculate the KL divergence   
    
   CFC1(nchan1i,nchan2i,:,:) = kl_cal(CStemp12,CFCtempY2);

  end
  end
  
  % conditon2   phase freq from dist2, amplitude freq from dist1
   if(isequaln(dist1,dist2)==0)   % dist1~=dist2
  CFC2 = zeros(nchan2,nchan1,nv,nfs);   
  for nchan2i = 1:nchan2   % phase freq from dist2      
  for nchan1i = 1:nchan1   % amplitude freq from dist1
      
   CFCtempX2   = squeeze(dist2.Xffts(nchan2i,:,:,:));
   CFCtempY1   = squeeze(dist1.Yffts(nchan1i,:,:,:)); 
   CStemp21    = CFCtempX2.*conj(CFCtempY1);
   
      % Calculate the KL divergence 
CFC2(nchan2i,nchan1i,:,:) =  kl_cal(CStemp21,CFCtempY1);  
  
  end
  end
  end
 
CFCs.CFC1   =  CFC1;
 if(isequaln(dist1,dist2)==0)   % dist1~=dist2
CFCs.CFC2   =  CFC2;
end
CFCs.method = para.method;
CFCs.freqVec_X = dist1.freqVec_X;
CFCs.freqVec_Y = dist1.freqVec_Y;  
 
              
end

end

%%%%%%%%%%%%%% Calculate the KL divergence%%%%%%%%%% 


function KLtemp = kl_cal(CFCtempCS,CFCtempYffts)
 
 
nbins = 20; 
KLtemp = zeros(size(CFCtempCS,1),size(CFCtempCS,2));

for k=1:size(CFCtempCS,1)
    for l=1:size(CFCtempCS,2)
        Ang = pi+squeeze(angle(CFCtempCS(k,l,:)));    
        m = 1+(floor(nbins*(Ang/(2*pi))));
        
        
        P = abs(CFCtempYffts(k,l,:));
        Pdist = zeros(1,nbins); 
        for n=1:length(Ang)
            if m(n) > nbins
                m(n) = nbins;
            end
            if m(n) < 1
                m(n) = 1;
            end
            Pdist(m(n)) = Pdist(m(n)) + P(n); 
        end
        
        Pdist = Pdist/sum(Pdist);        
        Q = mean(Pdist)*ones(1,nbins);      
        

        for n=1:nbins
            if Pdist(n) ~= 0
                KLtemp(k,l) = KLtemp(k,l) + Pdist(n)*log(Pdist(n)/Q(n));
            end
        end
    end
end
   
 KLtemp= KLtemp/log(nbins);       

end


end





      
 








    
 


