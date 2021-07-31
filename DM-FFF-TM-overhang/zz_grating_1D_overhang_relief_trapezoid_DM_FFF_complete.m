% planar diffraction by a sinusoidal relif grating
clear all
clc
lam    = 1;
%permW = 1.0;
%permW = 1.72^2;
permW = 1.5^2;
%permW  = (1+5i)^2;   %NOTE that the permitivity is complex!!!
thI    = 0*(pi/180);%+1e-3*(pi/180);
epsB   = 1;
Lam    = 1.5;
d      = 1.5;
epsW   = permW;
epsS   = epsW;
k0     = 2*pi/lam;
q      = lam/Lam;
q0     = sin(thI);

Phi=0*(pi/180);
%Phi=30*(pi/180); 

%% Control parameters
StrucParam.CS='CC_CS'; 
                                      %'CC_CS': default,cos^2 and cos*sin: according to Popov, et.al. , and Habib Mohamad et. al.
                                                                          
StrucParam.filtering=0; %0: default, no filtering; 1: filtering; According to Nikolay M. Lyndin, et. al.                     
StrucParam.threshold=20;  %for eigenvalue
%Near field plot;
StrucParam.PlotField='NO';  %'NO':default,not plot field;
                                        %'YES',plot near field, quite slow.
StrucParam.h_upper=1;      %h_upper=StrucParam.h_upper*h_grating;
StrucParam.h_lower=1;      %h_lower=StrucParam.h_lower*h_grating;
StrucParam.Number_of_Period=2;
StrucParam.Resolution_x=200;
StrucParam.Resolution_z=200;                                        
%% Functions of sin and cos tangent to the profile

K1=2*pi/(Lam*cos(Phi));

%% Control parameters
StrucParam.CS='CC_CS'; 
                                      %'CC_CS': default,cos^2 and cos*sin: according to Popov, et.al. , and Habib Mohamad et. al.
                                                                          
%Near field plot;
StrucParam.PlotField='NO';  %'NO':default,not plot field;
                                        %'YES',plot near field, quite slow.
StrucParam.h_upper=1;      %h_upper=StrucParam.h_upper*h_grating;
StrucParam.h_lower=1;      %h_lower=StrucParam.h_lower*h_grating;
StrucParam.Number_of_Period=2;
StrucParam.Resolution_x=200;
%StrucParam.Resolution_z=200;                                        
%% Functions of sin and cos tangent to the profile

K1=2*pi/(Lam*cos(Phi));

%%{
%NVM
%simple profile which can be discribed by a single-value function
%exaple: simple sinusodial profile.
%fsin=@(x) 1./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce cos(phi(x))
%fcos=@(x) (((pi*d)/(Lam)).*(sin(K1.*x)))./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce sin(phi(x))

%comple profile described by parametrized function, like circle, overhang sinusodial function etc.
%Parametrized derivatives.
%fsin=(dx/dt)/norm;
%fcos=(dy/dt)/norm=sqrt(1-fsin^2);
    
%NOTE that here cos and sin have been exchanged relative to the definition above
%fcos=C^2; fsin=C*S
fcos=@(t) 1./...
              (1+ (0*(t>=0 && t<Lam/4)+8/3*(t>=Lam/4 && t<5/8*Lam)+0*(t>=5/8*Lam && t<3/4*Lam)+8*(t>=3/4*Lam && t<=Lam)).^2);
fsin=@(t) (0*(t>=0 && t<Lam/4)+8/3*(t>=Lam/4 && t<5/8*Lam)+0*(t>=5/8*Lam && t<3/4*Lam)+8*(t>=3/4*Lam && t<=Lam))./...
             (1+ (0*(t>=0 && t<Lam/4)+8/3*(t>=Lam/4 && t<5/8*Lam)+0*(t>=5/8*Lam && t<3/4*Lam)+8*(t>=3/4*Lam && t<=Lam)).^2);

%}

%{  
%NOTE that here cos and sin have been exchanged relative to the definition above
%fcos=C^2; fsin=C*S    
fcos=@(x) 0.*x;
fsin=@(x) 0.*x;
%}

%% Preallocate fields for reflection amplitudes and errors
errorS_R=zeros(1,1);
errorP_R=zeros(1,1);
errorS_T=zeros(1,1);
errorP_T=zeros(1,1);

%% Reference values for sinusiodal grating, computed by C-Method

%%n=1.5,theta=0[deg],Phi=14.03 [deg]; C method,
RP_ref=[0.00524116141559641,0.00319977467116945,0.000577897490217867];
TP_ref=[0.0340806792370740,0.746840505991768,0.0475129752805627,0.150213890203357,0.0123931488599989];
%}

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =10;         % lowest number of modes
nMax_u    =10;         % highest number of modes
nMax_step = 1;          % length of the step - modes

%% %%%%%%%% Better to choose some constant value, and adjust it according to instability %%%
Nly_l    =10;           % lower number of layers
Nly_u    = 10;           % upper number of layers
Nly_step = 1;             % length of the step - layers
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Npt_l    = 10;           % points
Npt_u    = 10;           % pints
Npt_step =1;             % length of the step - layers

tic                   % start time count
RPvec=[];
TPvec=[];

RSvec=[];
TSvec=[];

for nMax =  nMax_l:nMax_step:nMax_u
    
    for N = Nly_l:Nly_step:Nly_u
        for Npoints=Npt_l:Npt_step:Npt_u
        tic
        nMax
        Npoints
    %% Relief profile
      %d1 = d/N;
      %%zV = linspace(d1/2,d-d1/2,N);
      %zV=(0:N)*d1;   %In fact it is boundary position.
      Ntot=N*Npoints+1;
      zV=linspace(0,d,Ntot);

      t=zeros(2,Ntot);
     xt=zeros(2,Ntot);
     xt_original=zeros(2,Ntot);
     wV=zeros(1,Ntot);
     epst = zeros(2,Ntot);
     
     cosx=zeros(2,Ntot);
     sinx=zeros(2,Ntot);
     x_cs=zeros(1,Ntot);

     epst(1,:)= epsB;
     epst(2,:) = epsW;

    for k=1:Ntot
    %t(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    %t(2,k)=Lam*cos(Phi)-t(1,k);
    t(1,k)=3/8*zV(k)+Lam/4;
    t(2,k)=1/8*zV(k)+3/4*Lam;
    
    cosx(1,k)=fcos(t(1,k));
    cosx(2,k)=fcos(t(2,k));
    
    sinx(1,k)=fsin(t(1,k));
    sinx(2,k)=fsin(t(2,k));
    
    xt(1,k)=sec(Phi)*t(1,k)+tan(Phi)*zV(k);
    xt(2,k)=sec(Phi)*t(2,k)+tan(Phi)*zV(k);
    
    xt_original(1,k)=xt(1,k);
    xt_original(2,k)=xt(2,k);
    
    xt(1,k)=xt(1,k)/Lam;
    xt(2,k)=xt(2,k)/Lam;
    if xt(1,k)<1 && xt(2,k)>1
        x0=xt(1,k);
        xt(1,k)=mod(xt(2,k),1);
        xt(2,k)=x0;
        epst(1,k)= epsW;
        epst(2,k) = epsB;
        
        cos0=cosx(1,k);
        cosx(1,k)=cosx(2,k);
        cosx(2,k)=cos0;
       
       sin0=sinx(1,k);
       sinx(1,k)=sinx(2,k);
       sinx(2,k)=sin0;
    
    elseif xt(1,k)>1
        xt(1,k)=mod(xt(1,k),1);
        xt(2,k)=mod(xt(2,k),1);
    end
   end
    %}
     x_cs=(xt(1,:)+xt(2,:))/2;
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(StrucParam.PlotField,'NO')
        [RP,s0V,TP,sSubV] = computeScatMatNVM_DM(lam,thI,epsB,Lam,epsS,sinx,cosx,xt,zV,x_cs,epst,nMax,N,Npoints);
        else
        [RP,s0V,TP,sSubV,H,HE_m,Z_p1_mat,F,T_mat,zV_tot] =computeScatMatNVM_DM_field(lam,thI,epsB,Lam,epsS,sinx,cosx,xt,xt_original,d,zV,x_cs,epst,nMax,N,Npoints,StrucParam) ;
        end
 
        % RESULTS
        index_R=find(imag(s0V)==0);
        RPvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)  = (abs((RP(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        RP=RPvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:) ;
        RP=RP(:);
        %errorP_R((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step) =1/numel(RP_ref)*sum(sum(sum((RPvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)-RP_ref).^2./RP_ref.^2)));      
         
        %%{
        index_T=find(imag(sSubV)==0);
        TPvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)  =epsB/epsS* (abs((TP(index_T)').^2).*sSubV(1,index_T))./s0V(1,nMax+1);
        TP= TPvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:) ;
        TP=TP(:);
        
        RT=[RP;TP];
        RT_ref=[RP_ref(:);TP_ref(:)];
        %errorP_T((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step) =1/numel(TP_ref)*sum(sum(sum((TPvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)-TP_ref).^2./TP_ref.^2)));      
       
        error((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step)=max(abs(RT-RT_ref)./RT_ref);
        %}
        
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step) = toc;
 
end
end
end
       tot_Run_time=sum(sum(c_time));
%       error_study=error(:,:,7);
%% Save results to a file
%filename = 'dielectric_overhang_relief_trapezoid_TM_DM_1percent.mat';
filename = 'dielectric_overhang_relief_trapezoid_TM_DM_1percent_N=21_L=10_Npt=10.mat';
save(filename)
RPvec
RP_ref
TPvec
TP_ref
tot_Run_time