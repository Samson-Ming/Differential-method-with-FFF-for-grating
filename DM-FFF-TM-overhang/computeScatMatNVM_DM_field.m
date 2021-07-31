function [RP,s0V,TP,sSubV,H,HE_m,Z_p1_mat,F,T_mat,zV_tot] =computeScatMatNVM_DM_field(lam,thI,epsB,Lam,epsS,sinx,cosx,xt,xt_original,d,zV,x_cs,epst,nMax,N,Npoints,StrucParam)      
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate
nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sqrt(epsB)*sin(thI);
I=eye(nDim);
II=eye(2*nDim);
zero=eye(nDim);

%% Prepare Factorization Matrices
nV=-nMax:nMax;
qV  = q0+nV*q;
kX   = diag(qV);  
    
s0V = sqrt(epsB-qV.^2);
sSubV = sqrt(epsS-qV.^2);
 
%{
   epsMn = epsB*I;
    etaMn = 1/epsB*I;
    %[cosM,sinM] = generateSinCosMat([0,0],[0,0],x_cs(1),nMax); 
   cosM=zero;
   sinM=zero;
    [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    %CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    CP=-[B*D*kX,(B*D*B-A-etaMn\I);(kX*D*kX-I),kX*D*B];
    M0=1i*k0*CP;
    %}
%% W0 and W_Np1 
    W_0_11=I;                             % W_1_11=zero,%W_1_21=I
    W_0_12=I;
    W_0_21=diag(sSubV/epsS);
    W_0_22=-W_0_21;
    W_0=[W_0_11,W_0_12;W_0_21,W_0_22];
    
    W_Np1_11=I;
    W_Np1_12=I;
    W_Np1_21=diag(s0V/epsB);
    W_Np1_22=-W_Np1_21;
    W_Np1=[W_Np1_11,W_Np1_12;W_Np1_21,W_Np1_22];

%% S0    
    R_ud0=zeros(nDim);
    R_du0=zeros(nDim);
    T_uu0=I;
    T_dd0=I;
    
   R_ud=zeros(nDim,nDim,N);
   T_dd=zeros(nDim,nDim,N);
 %---------------------------------------------------------  
    N_tot=N*Npoints+1;
    F=zeros(2*nDim,2*nDim,N_tot);  
   
    %F(:,:,1)=W_0;
    %for layer 1, bourndary 1, i.e. z=0;
    % in fact bounary(1)-  and boundary(1)+
    h=zV(2)-zV(1);
    z_index0=(1-1)*Npoints+1;
    %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index0),epst(2,z_index0),xt(:,z_index0),nMax);
    etaMn = epsMatrix(1/epst(1,z_index0),1/epst(2,z_index0),xt(:,z_index0),nMax);
    [cosM,sinM] = generateSinCosMat(sinx(:,z_index0),cosx(:,z_index0),x_cs(z_index0),nMax); 
   
    [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    %CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    CP=-[B*D*kX,(B*D*B-A-etaMn\I);(kX*D*kX-I),kX*D*B];
    M_j=1i*k0*CP;
    
    
    T_mat=zeros(2*nDim,2*nDim,N);
    Z_p1_mat=zeros(nDim,nDim,N);
    
    T=W_0;
    F(:,:,z_index0)=T;
    for z_index=z_index0+1:z_index0+Npoints
        %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index),epst(2,z_index),xt(:,z_index),nMax);
    etaMn = epsMatrix(1/epst(1,z_index),1/epst(2,z_index),xt(:,z_index),nMax);
    [cosM,sinM] = generateSinCosMat(sinx(:,z_index),cosx(:,z_index),x_cs(z_index),nMax); 
   
    [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    %CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    CP=-[B*D*kX,(B*D*B-A-etaMn\I);(kX*D*kX-I),kX*D*B];
    
    M_j1=1i*k0*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    MM=(II-1/2*h*M_j1)\(II+1/2*h*M_j);
    T=MM*T;
    F(:,:,z_index)=T;
    M_j=M_j1;
    end
    T=W_Np1\T;

    Z_p1=I/(T(nDim+1:2*nDim,1:nDim)*R_ud0+T(nDim+1:2*nDim,nDim+1:2*nDim));
    
    
    R_ud(:,:,1)=(T(1:nDim,nDim+1:2*nDim)+T(1:nDim,1:nDim)*R_ud0)*Z_p1;                   %S+-, (7.153)
    T_dd(:,:,1)=T_dd0*Z_p1;                     %S--, (7.156)   
    
    T_mat(:,:,1)=T;
    Z_p1_mat(:,:,1)=Z_p1;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %for layer 1:N
    for iLy=2:N
         T=W_Np1;
         z_index0=(iLy-1)*Npoints+1;
         F(:,:,z_index0)=T;
    for z_index=z_index0+1:z_index0+Npoints
        %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index),epst(2,z_index),xt(:,z_index),nMax);
    etaMn = epsMatrix(1/epst(1,z_index),1/epst(2,z_index),xt(:,z_index),nMax);
    [cosM,sinM] = generateSinCosMat(sinx(:,z_index),cosx(:,z_index),x_cs(z_index),nMax); 
   
    [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    CP=-[B*D*kX,(B*D*B-A-etaMn\I);(kX*D*kX-I),kX*D*B];
    M_j1=1i*k0*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    MM=(II-1/2*h*M_j1)\(II+1/2*h*M_j);
    T=MM*T;
    F(:,:,z_index)=T;
    M_j=M_j1;
    end
    T=W_Np1\T;

    Z_p1=I/(T(nDim+1:2*nDim,1:nDim)*R_ud(:,:,iLy-1)+T(nDim+1:2*nDim,nDim+1:2*nDim));
    
    
    R_ud(:,:,iLy)=(T(1:nDim,nDim+1:2*nDim)+T(1:nDim,1:nDim)*R_ud(:,:,iLy-1))*Z_p1;                   %S+-, (7.153)
    T_dd(:,:,iLy)=T_dd(:,:,iLy-1)*Z_p1;                     %S--, (7.156)
    T_mat(:,:,iLy)=T;
     Z_p1_mat(:,:,iLy)=Z_p1;
    end
    
    
    %{
    %for the superstrate if dispersive
    T=II;
    Z_p1=(T(nDim+1:2*nDim,1:nDim)*R_ud+T(nDim+1:2*nDim,nDim+1:2*nDim))\I;
    
    R_ud=(T(1:nDim,nDim+1:2*nDim)+T(1:nDim,1:nDim)*R_ud)*Z_p1;                   %S+-, (7.153)
    T_dd=T_dd*Z_p1;                     %S--, (7.156) 
    %}  
        

 %---------------------------------------------------------
   % R and T
 %---------------------------------------------------------
        d_p1=zeros(nDim,1);
        d_p1(nMax+1,1)=1;
%---------------------------------------------------------

        RP=R_ud(:,:,end)*d_p1;
        TP=T_dd(:,:,end)*d_p1;

%---------------------------------------------------------   
%Near field
        H=zeros(N_tot,StrucParam.Resolution_x+1);        
        
        A_m=zeros(2*nDim,N_tot);
        HE_m=zeros(2*nDim,N_tot);
                
       A_m(:,1)=[d_p1*0;TP];
       A_m(:,end)=[RP;d_p1];
       
       %for Layer N
       iLy=N;
       z_index0=(iLy-1)*Npoints+1;
       A_m(nDim+1:2*nDim,z_index0)=Z_p1_mat(:,:,iLy)*d_p1;%-Z_p1_mat(:,:,N)*T_mat(Dim+1:2*nDim,1:nDim,iLy)*(dp1*0)
       A_m(1:nDim,z_index0)=R_ud(:,:,iLy-1)*A_m(nDim+1:2*nDim,z_index0);
       
       F_0=A_m(:,z_index0);
       
    for z_index=z_index0:z_index0+Npoints-1
         HE_m(:,z_index)= F(:,:,z_index)*F_0;
    end  
    
    %for layer 2:N-1
    
    for iLy=N-1:-1:2
        z_index0=(iLy-1)*Npoints+1;
        A_m(nDim+1:2*nDim,z_index0)=Z_p1_mat(:,:,iLy)*F_0(nDim+1:2*nDim);%-Z_p1_mat(:,:,N)*T_mat(Dim+1:2*nDim,1:nDim,iLy)*(dp1*0)
        A_m(1:nDim,z_index0)=R_ud(:,:,iLy-1)*A_m(nDim+1:2*nDim,z_index0);
       
       F_0=A_m(:,z_index0);
       
   for z_index=z_index0:z_index0+Npoints-1
         HE_m(:,z_index)= F(:,:,z_index)*F_0;
    end
    end
    
    %for layer 1
    
    iLy=1;
    z_index0=(iLy-1)*Npoints+1;
    A_m(nDim+1:2*nDim,z_index0)=TP;%-Z_p1_mat(:,:,N)*T_mat(Dim+1:2*nDim,1:nDim,iLy)*(dp1*0)
   A_m(1:nDim,z_index0)=0*TP;
   F_0=A_m(:,z_index0);
   for z_index=z_index0:z_index0+Npoints-1
         HE_m(:,z_index)= F(:,:,z_index)*F_0;
    end
   
   
   %----------------------------------------------------------------------------
    
       H_m=HE_m(1:nDim,:);
        
        x=(linspace(0,StrucParam.Number_of_Period*Lam,StrucParam.Resolution_x+1))';

        exp_x=(exp(1i*k0*x*qV)).';
        
        for iL=1:N_tot         
       H(iL,:)=H_m(:,iL).' * exp_x;
        end
        
        d_3=StrucParam.h_upper*d;      %h_upper=StrucParam.h_upper*h_grating;
        N_3=round(d_3/h)+1;
        zV_3=d+(0:N_3-1)*h;
        d_0=StrucParam.h_lower*d;      %h_lower=StrucParam.h_lower*h_grating;
        N_0=round(d_0/h)+1;
        zV_0=(1-N_0:0)*h;
        
        H_3=zeros(N_3,StrucParam.Resolution_x+1);
        d_p1_u=diag(exp(-1i*k0*s0V*(zV_3(end)-zV_3(1))))*d_p1;
        for iL=1:N_3
            H_3_m=RP.*exp(1i*k0*s0V.'*(zV_3(iL)-zV_3(1)))  + d_p1_u.*exp(-1i*k0*s0V.'*(zV_3(iL)-zV_3(end)));           
            H_3(iL,:)=H_3_m.' * exp_x;
        end
        
        H_0=zeros(N_0,StrucParam.Resolution_x+1);
        for iL=1:N_0
            H_0_m=TP.*exp(-1i*k0*sSubV.'*(zV_0(iL)-zV_0(end)));
            H_0(iL,:)=H_0_m.' * exp_x;
        end
        
        H_tot=[H_0;H(2:end-1,:);H_3];
        zV_tot=[zV_0,zV(2:end-1),zV_3];        
        
        
figure;
% imagesc(x,zV_tot,abs(H_tot));
pcolor(x,zV_tot,abs(H_tot));
shading interp;
 %axis equal
 axis xy
 colorbar
 colormap jet

hold on
 
 for ix=1:StrucParam.Number_of_Period         
      plot(xt_original-Lam+(ix-1)*Lam,zV,'w','linewidth',1.5);
      plot(xt_original+(ix-1)*Lam,zV,'w','linewidth',1.5);
 end
 hold off
        