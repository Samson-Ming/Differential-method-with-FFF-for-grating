% planar diffraction by a sinusoidal relif grating
clear all
clc
lam    = 1;
%permW = 1.0;
%permW = 1.72^2;
%permW = 3^2;
%permW  = (1+5i)^2;   %NOTE that the permitivity is complex!!!
permW  = (0.1+10i)^2;   %NOTE that the permitivity is complex!!!
%permW  = (1+10i)^2;   %NOTE that the permitivity is complex!!!

thI    = 15*(pi/180);%+1e-3*(pi/180);
epsB   = 1;
Lam    = 425e-3/525e-3;
d      = 360e-3/525e-3;
epsW   = permW;
epsS   = epsW;
k0     = 2*pi/lam;
q      = lam/Lam;
q0     = sin(thI);

%Phi=0*(pi/180);
Phi=30*(pi/180); 

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
fcos=@(t) (sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2./((sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2+(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2); %predpis funkce cos(phi(x))
fsin=@(t) (sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).*(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t)))./((sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2+(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2); %predpis funkce cos(phi(x))

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

%{
%%n=1.72,theta=15[deg],Phi=30 [deg]; C method,WRONG,sort beta, need rotation
RP_ref=[0.004048697696454   0.002273649038408];
RS_ref=[0.020351069199529   0.002595456005422];
TP_ref=[0.204563103140346   0.715328171657040   0.073786421337080];
TS_ref=[0.465951079541316   0.429741994958952   0.081360436389336];
%}

%{
%%n=1.72,theta=15[deg],Phi=30 [deg]; LPEM
RP_ref=[0.013847619286624   0.001321278113111];
RS_ref=[0.033479360118408   0.003225124752801];
TP_ref=[0.276865218496503   0.698838843903608   0.009143687581078];
TS_ref=[0.429215218389913   0.503996524935058   0.030083766291384];
%}

%%{
%%n=1+5i,theta=15[deg],Phi=30 [deg]; 
%FDTD %PLEASE DONOT DO THIS TEST.
%RP_ref=[];
%RS_ref=[];
%LPEM
RP_ref=[0.690696929434143   0.044679728260437];
RS_ref=[0.271923709974736   0.502654339036197];
%FEM
RP_ref=[0.6969397576849351	0.042593465794829075];
RS_ref=[0.27262164737017386	0.5021319693110445];
%Stardard FFM:RCWA-1D+S method:N=401,M=200;
RP_ref=[0.685220257896775 0.0451030975250511];
RS_ref=[0.271998824578871 0.502563618626655];
%Stardard FFM:Reticolo+modified T method:N=401,M=200;
RP_ref=[0.685220257895786 0.045103097525231];
RS_ref=[0.271998824578443 0.502563618627200];
%FFM-FFF:N=401,M=200;
RP_ref=[0.685220257896644 0.045103097524921];  %Staircase_CS

RP_ref=[0.685220257896636  0.045103097524920];  %Staircase_CCCS

RP_ref=[0.679735717417338 0.043012742570048];  %FFF_CS,layer by layer

RP_ref=[0.680169300567980 0.043733149380730];  %FFF4_CCCS,layer by layer


%FFM-FFF-new:N=61,M=100;m=12;
PR_ref=[0.593521420734058 0.033795070317326];  %Standard RCWA-1D+S
RP_ref=[0.593521420734065 0.033795070317335];  %Staircase

RP_ref=[0.583566274577521 0.032372514988748];  %FFF,once

RP_ref=[0.667510806093971, 0.034303475841890];  %FFF_CS,layer by layer

RP_ref=[0.668736568932016 0.041550849121362];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=81,M=200;m=12;
RP_ref=[0.611794654180626 0.036396764482898];  %Standard RCWA-1D+S
RP_ref=[0.611794654180615 0.036396764482898];  %Staircase_CS
RP_ref=[0.611794654180613 0.036396764482897];  %Staircase_CCCS

RP_ref=[0.608290185707760 0.035151844834948];  %FFF,once
RP_ref=[0.682107033399329, 0.038013455352558];  %FFF_CS,layer by layer
RP_ref=[0.681599347907159, 0.042963578909171];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=101,M=300;m=12;
RP_ref=[0.626801354406051 0.037930934014082];  %Standard RCWA-1D+S,W-S
RP_ref=[0.626801354406098 0.037930934014083];  %Staircase,C-S,W-S


RP_ref=[0.624271398463843 0.037365482267794];  %FFF,once

RP_ref=[0.684959865382469 0.039917326428615];  %FFF_CS,layer by layer
RP_ref=[0.685056825574264 0.043589805451875];  %FFF4_CCCS,layer by layer

%Differential method
RP_ref=[0.690420241599464 0.044574616761142];  %N=35,N_layer=20,N_int_pt=50;

%C method
RP_ref=[0.691131155045779   0.044826802871665];
RS_ref=[0.271937696305567   0.502697416029626];
%}

RP_ref=[0.956070372866905 0.036126429767618];

%RP_ref=[0.897562186983547 0.028743965469486];

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =50;         % lowest number of modes
nMax_u    =50;         % highest number of modes
nMax_step = 1;          % length of the step - modes

Nly_l    =45;           % lower number of layers
Nly_u    = 45;           % upper number of layers
Nly_step = 1;             % length of the step - layers

Npt_l    = 60;           % lower number of layers
Npt_u    = 60;           % upper number of layers
Npt_step = 1;             % length of the step - layers

tic                   % start time count
RPvec=[];
TPvec=[];

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

    for k=1:Ntot
    t(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    t(2,k)=Lam*cos(Phi)-t(1,k);
    
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
        RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-Nly_l+Nly_step)/Nly_step,:)  = (abs((RP(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        errorP_R((nMax-nMax_l+nMax_step)/nMax_step,(N-Nly_l+Nly_step)/Nly_step) =1/numel(RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-Nly_l+Nly_step)/Nly_step,:))*sum(sum(sum((RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-Nly_l+Nly_step)/Nly_step,:)-RP_ref).^2./RP_ref.^2)));      
        %RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((RS(index_R)').^2).*s0V(1,index_R)./s0V(1,nMax+1);
        %errorS_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));
        
        %{
        index_T=find(imag(sSubV)==0);
        TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  =epsB/epsS* (abs((TP(index_T,nMax+1)').^2).*sSubV(1,index_T))./s0V(1,nMax+1);
        TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((TS(index_T,nMax+1)').^2).*sSubV(1,index_T)./s0V(1,nMax+1);
        errorS_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TS_ref).^2./TS_ref.^2)));
        errorP_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TP_ref).^2./TP_ref.^2)));
        %}
        
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-Nly_l+Nly_step)/Nly_step) = toc;
    end
end
end
%end
       tot_Run_time=sum(sum(c_time));
%% Save results to a file
%filename = 'test.mat';
%save(filename)
%exit
%save lossy_overhang_sinusodial_FMM_FFF_vs_LPEM4.mat
RPvec
RP_ref