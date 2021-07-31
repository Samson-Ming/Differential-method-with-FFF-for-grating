function [RP,s0V,TP,sSubV,H] =computeScatMatNVM_DM_Field(lam,thI,epsB,Lam,d,epsS,sinx,cosx,xt,xt_original,zV,x_cs,epst,nMax,N,Npoints,StrucParam) 
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate


nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sqrt(epsB)*sin(thI);
I=eye(nDim);
zero=zeros(nDim);
nV=-nMax:nMax;

%% Prepare Factorization Matrices
qV  = q0+nV*q;
kX   = diag(qV);  
    
s0V = sqrt(epsB-qV.^2);
sSubV = sqrt(epsS-qV.^2);

%% S0
    
    R_ud=zeros(nDim);
    R_du=zeros(nDim);
    T_uu=I;
    T_dd=I;
 %---------------------------------------------------------  
 t=zeros(2*nDim,2*nDim);

 %for the first layer
    zpt=linspace(zV(1),zV(2),Npoints+1);   % in fact bounary(1)-  and boundary(1)+
    h=zpt(2)-zpt(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %For the substrate
    epsMn = epsMatrix(epst(1,1),epst(2,1),xt(:,1),nMax);
    etaMn = epsMatrix(1/epst(1,1),1/epst(2,1),xt(:,1),nMax);
    [cosM,sinM] = generateSinCosMat(sinx(:,1),cosx(:,1),x_cs(1),nMax); 
   
    [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
    CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    M=1i*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    
   Mj=
    for ipt=1:Npoints
        

 
 %---------------------------------------------------------
        d_p1=zeros(nDim,1);
        d_p1(nMax+1,1)=1;
        
        u_0=zeros(nDim,1);
%---------------------------------------------------------

        RP=R_ud(:,:,end)*d_p1;
        TP=T_dd(:,:,end)*d_p1; 

%H=0;
%%{
%--------------------------------------------------------- 
d_3=StrucParam.h_upper*d;      %h_upper=StrucParam.h_upper*h_grating;
d_0=StrucParam.h_lower*d;      %h_lower=StrucParam.h_lower*h_grating;
d_tot=d_0+d+ d_3;
a_0=[diag(exp(-1i*k0*sSubV*d_0))*u_0;TP];
a_p1=[RP;diag(exp(-1i*k0*s0V*d_3))*d_p1];  %

a_2=zeros(2*nDim,N);

%{
di=T_dd(:,:,1)\TP;
a_2(:,1)=[R_ud(:,:,1)*di;di];
flag=1
if N>2
for iL=1:N-1
%     di=T_dd(:,:,iL)\TP;
%     a_2(:,iL)=[R_ud(:,:,iL)*di;di];

     T_matrix=[W_2_11(:,:,iL+1),W_2_12(:,:,iL+1);W_2_21(:,:,iL+1),W_2_22(:,:,iL+1)]\[W_2_11(:,:,iL),W_2_12(:,:,iL);W_2_21(:,:,iL),W_2_22(:,:,iL)]*[Phi_P(:,:,iL),0*I;0*I,Phi_N(:,:,iL)];
    
     a_2(1:nDim,iL+1)=T_matrix(1:nDim,1:nDim)*a_2(1:nDim,iL)+T_matrix(1:nDim,nDim+1:2*nDim)*a_2(nDim+1:2*nDim,iL);
     a_2(nDim+1:2*nDim,iL+1)=T_matrix(nDim+1:2*nDim,1:nDim)*a_2(1:nDim,iL)+T_matrix(nDim+1:2*nDim,nDim+1:2*nDim)*a_2(nDim+1:2*nDim,iL);
     flag=iL
end
end
%}

%%{
a_2(1:nDim,1)=s_0(1:nDim,nDim+1:2*nDim)*(s_0(nDim+1:2*nDim,nDim+1:2*nDim)\TP);

% Pomocí translace volným prostorem provedu transformaci u_1_0 na u_1_1
%u_1_1=X1(:,:,1)*c_m(:,:,N);

% Výpočet c_p(:,:,1) --- z=0

%a_2(1:nDim,N)=s_Np1(1:nDim,1:nDim)\(RP-s_Np1(1:nDim,1+nDim:2*nDim)*d_p1);
u_n_n=s_Np1(1:nDim,1:nDim)\(RP-s_Np1(1:nDim,1+nDim:2*nDim)*d_p1);

a_2(nDim+1:2*nDim,N)=s_Np1(nDim+1:2*nDim,1:nDim)*u_n_n+s_Np1(nDim+1:2*nDim,nDim+1:2*nDim)*d_p1;

d_n_plus_1_n(:,N)=Phi_N_inv(:,:,N)*a_2(nDim+1:2*nDim,N); 
%nepřehledný for-cyklus, nejlépe kouknou do poznámek
if N>=2    
    for i=1:N-1
        a_2(1:nDim,N-i+1)=R_ud(:,:,N-i+1)*d_n_plus_1_n(:,N-i+1); %využívám S-matice --- koeficient R_ud           
        
        inv_s2=inv(s_2(1:nDim,1:nDim,N-i));
        
       u_n_n= inv_s2*(a_2(1:nDim,N-i+1)-s_2(1:nDim,1+nDim:2*nDim,N-i)*a_2(1+nDim:2*nDim,N-i+1));
        
        a_2(nDim+1:2*nDim,N-i)= s_2(1+nDim:2*nDim,1:nDim,N-i)*u_n_n+ s_2(nDim+1:2*nDim,nDim+1:2*nDim,N-i)*d_n_plus_1_n(:,N-i+1);
       
        d_n_plus_1_n(:,N-i)=Phi_N_inv(:,:,N-i)*a_2(nDim+1:2*nDim,N-i); 
    end
end
%}


%Plot field
x=(linspace(0,StrucParam.Number_of_Period*Lam,StrucParam.Resolution_x+1))';

exp_x=(exp(1i*k0*x*qV)).';

% boundary of layers
%{
 boundary=zeros(1,N+1);
 for i=1:1:N
      boundary(i+1)=sum(layer_thickness(1:i))*1e6;    
 end
%}
 boundary=(0:N)*d1;
 
z= (linspace(-d_0,d_tot-d_0,StrucParam.Resolution_z+1))';

d1 = d/N;
boundary=(0:N)*d1;
boundary=[-d_0,boundary,d_tot-d_0];
a_2=[a_0,a_2,a_p1];

H=zeros(StrucParam.Resolution_z+1,StrucParam.Resolution_x+1);

Hymq_u=zeros(nDim,nDim,N+2);
Hymq_u(:,:,1)=I;
Hymq_u(:,:,N+2)=I;
Hymq_u(:,:,2:N+1)=W_2_11(:,:,1:N);

Hymq_d=zeros(nDim,nDim,N+2);
Hymq_d(:,:,1)=I;
Hymq_d(:,:,N+2)=I;
Hymq_d(:,:,2:N+1)=W_2_12(:,:,1:N);
Qu=k0*[sSubV.',sVP,s0V.'];
Qd=k0*[-sSubV.',sVN,-s0V.'];
 
 for z_index=1:StrucParam.Resolution_z+1
        for iL=1:N+2
            if z(z_index)>=boundary(iL) && z(z_index)<boundary(iL+1)
            Hym=Hymq_u(:,:,iL)*(a_2(1:nDim,iL).*exp(1i*Qu(:,iL)*(z(z_index)-boundary(iL))))+Hymq_d(:,:,iL)*(a_2(nDim+1:2*nDim,iL).*exp(1i*Qd(:,iL)*(z(z_index)-boundary(iL+1))));
            H(z_index,:)=(Hym.')*exp_x;
            end
        end
 end
 
 figure;
 imagesc(x,z,abs(H));
 axis xy
 colorbar
 colormap jet

 hold on
 
 for ix=1:StrucParam.Number_of_Period
     mid=(xt_original(1,end)+xt_original(2,end))/2;
     xt_original1=[[0;Lam],xt_original,[mid;mid]];
     zV1=[0,zV,d];
      plot(xt_original1-Lam+(ix-1)*Lam,zV1,'k','linewidth',1.5);
      plot(xt_original1+(ix-1)*Lam,zV1,'k','linewidth',1.5);
 end
 hold off
%}            