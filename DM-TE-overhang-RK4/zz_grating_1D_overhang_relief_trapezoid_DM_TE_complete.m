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

%% Preallocate fields for reflection amplitudes and errors
errorS_R=zeros(1,1);
errorS_T=zeros(1,1);

%% Reference values for sinusiodal grating, computed by C-Method

%% Reference values for sinusiodal grating, computed by C-Method
%%{
%%n=1.5,theta=0[deg],Phi=14.03 [deg]; C method,
RS_ref=[0.0116110339181890	0.00397254128892009	0.00269456222006578];
TS_ref=[0.0984302860119821	0.476629669959978	0.0994663070054246	0.268814629169694	0.0384145463940903];
%}


%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =20;         % lowest number of modes
nMax_u    =20;         % highest number of modes
nMax_step = 1;          % length of the step - modes

Nly_l    = 25;           % lower number of layers
Nly_u    = 25;           % upper number of layers
Nly_step = 1;             % length of the step - layers

Npt_l    = 40;           % lower number of layers
Npt_u    = 40;           % upper number of layers
Npt_step = 1;             % length of the step - layers

tic                   % start time count
RSvec=[];
TSvec=[];

RSvec=[];
TSvec=[];

for nMax =  nMax_l:nMax_step:nMax_u
    
    for N = Nly_l:Nly_step:Nly_u
        for Npoints=Npt_l:Npt_step:Npt_u
        tic
        
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
     
     
     zVm=(zV(1:end-1)+zV(2:end))/2;

      tm=zeros(2,Ntot-1);
     xtm=zeros(2,Ntot-1);
     xt_originalm=zeros(2,Ntot-1);
     wVm=zeros(1,Ntot-1);
     epstm = zeros(2,Ntot-1);
     
     cosxm=zeros(2,Ntot-1);
     sinxm=zeros(2,Ntot-1);
     x_csm=zeros(1,Ntot-1);

     epstm(1,:)= epsB;
     epstm(2,:) = epsW;
     

    for k=1:Ntot
     %t(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    %t(2,k)=Lam*cos(Phi)-t(1,k);
    t(1,k)=3/8*zV(k)+Lam/4;
    t(2,k)=1/8*zV(k)+3/4*Lam;
        
    xt(1,k)=sec(Phi)*t(1,k)+tan(Phi)*zV(k);
    xt(2,k)=sec(Phi)*t(2,k)+tan(Phi)*zV(k);
       
    xt(1,k)=xt(1,k)/Lam;
    xt(2,k)=xt(2,k)/Lam;
    if xt(1,k)<1 && xt(2,k)>1
        x0=xt(1,k);
        xt(1,k)=mod(xt(2,k),1);
        xt(2,k)=x0;
        epst(1,k)= epsW;
        epst(2,k) = epsB;
       
    elseif xt(1,k)>1
        xt(1,k)=mod(xt(1,k),1);
        xt(2,k)=mod(xt(2,k),1);
    end
    end
    %}
     x_cs=(xt(1,:)+xt(2,:))/2;
     
     for k=1:Ntot-1
     %tm(1,k)=acos(-(2*(zVm(k))/d-1))/(2*pi/Lam/cos(Phi));
    %tm(2,k)=Lam*cos(Phi)-tm(1,k);
    tm(1,k)=3/8*zVm(k)+Lam/4;
    tm(2,k)=1/8*zVm(k)+3/4*Lam;
        
    xtm(1,k)=sec(Phi)*tm(1,k)+tan(Phi)*zVm(k);
    xtm(2,k)=sec(Phi)*tm(2,k)+tan(Phi)*zVm(k);
       
    xtm(1,k)=xtm(1,k)/Lam;
    xtm(2,k)=xtm(2,k)/Lam;
    if xtm(1,k)<1 && xtm(2,k)>1
        x0=xtm(1,k);
        xtm(1,k)=mod(xtm(2,k),1);
        xtm(2,k)=x0;
        epstm(1,k)= epsW;
        epstm(2,k) = epsB;
       
    elseif xtm(1,k)>1
        xtm(1,k)=mod(xtm(1,k),1);
        xtm(2,k)=mod(xtm(2,k),1);
    end
    end
    %}
     x_csm=(xtm(1,:)+xtm(2,:))/2;
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(StrucParam.PlotField,'NO')
        [RS0,s0V,TS0,sSubV,T_dd,R_ud] = computeScatMatNVM_DM(lam,thI,epsB,Lam,epsS,xt,zV,epst,xtm,zVm,epstm,nMax,N,Npoints);
        else
        [RS0,s0V,TS0,sSubV,E]  =computeScatMatNVM_DM_Field(lam,thI,epsB,Lam,d,epsS,xt,xt_original,zV,epst,nMax,N,Npoints,StrucParam);
        end
 
        % RESULTS
        index_R=find(imag(s0V)==0);
        RSvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)  = (abs((RS0(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        RS=RSvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:) ;
        RS=RS(:);
        %errorS_R((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step) =1/numel(RS_ref)*sum(sum(sum((RSvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)-RS_ref).^2./RS_ref.^2)));      
         
        %%{
        index_T=find(imag(sSubV)==0);
        TSvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)  = (abs((TS0(index_T)').^2).*sSubV(1,index_T))./s0V(1,nMax+1);
        TS= TSvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:) ;
        TS=TS(:);
        
        RT=[RS;TS];
        RT_ref=[RS_ref(:);TS_ref(:)];
        %errorS_T((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step) =1/numel(TS_ref)*sum(sum(sum((TSvec((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step,:)-TS_ref).^2./TS_ref.^2)));      
       
        error((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step)=max(abs(RT-RT_ref)./RT_ref);
        %}
        
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(Npoints-Npt_l+Npt_step)/Npt_step) = toc;
end
end
end
       tot_Run_time=sum(sum(c_time))
%% Save results to a file
filename = 'dielectric_overhang_relief_trapezoid_TE_DM_1percent.mat';
%filename = 'dielectric_overhang_relief_trapezoid_TE_DM_1percent_N=21_L=10_Npt=10.mat';
save(filename)
%{
RSvec
RS_ref
TSvec
TS_ref
tot_Run_time
%}