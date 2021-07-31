function [RS,s0V,TS,sSubV,T_dd,R_ud] =computeScatMatNVM_DM(lam,thI,epsB,Lam,epsS,xt,zV,epst,xtm,zVm,epstm,nMax,N,Npoints)    
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
 
%% W0 and W_Np1 
    % %{
    W_0_11=I;                             % W_1_11=zero,%W_1_21=I
    W_0_12=I;
    W_0_21=-diag(sSubV);
    W_0_22=-W_0_21;
    W_0=[W_0_11,W_0_12;W_0_21,W_0_22];
    
    W_Np1_11=I;
    W_Np1_12=I;
    W_Np1_21=-diag(s0V);
    W_Np1_22=-W_Np1_21;
    W_Np1=[W_Np1_11,W_Np1_12;W_Np1_21,W_Np1_22];
%}

%{
    W_0_11=-diag(1./sSubV);                             % W_1_11=zero,%W_1_21=I
    W_0_12=-diag(1./sSubV);  
    W_0_21=I;
    W_0_22=-W_0_21;
    W_0=[W_0_11,W_0_12;W_0_21,W_0_22];
    
    W_Np1_11=-diag(1./s0V);
    W_Np1_12=-diag(1./s0V);
    W_Np1_21=I;
    W_Np1_22=-W_Np1_21;
    W_Np1=[W_Np1_11,W_Np1_12;W_Np1_21,W_Np1_22];
%}

%% S0    
    R_ud=zeros(nDim);
    R_du=zeros(nDim);
    T_uu=I;
    T_dd=I;
 %---------------------------------------------------------  
    
    %for layer 1, bourndary 1, i.e. z=0;
    % in fact bounary(1)-  and boundary(1)+
    h=zV(2)-zV(1);
    z_index0=(1-1)*Npoints+1;
    %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index0),epst(2,z_index0),xt(:,z_index0),nMax);
   
    CP=[zero,-1i*I;1i*(kX^2-epsMn),zero];
    M_j=k0*CP;   
    
    T=W_0;
    for z_index=z_index0+1:z_index0+Npoints
        %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index),epst(2,z_index),xt(:,z_index),nMax);
    
    CP=[zero,-1i*I;1i*(kX^2-epsMn),zero];
    M_j1=k0*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMnm = epsMatrix(epstm(1,z_index-1),epstm(2,z_index-1),xt(:,z_index-1),nMax);
    
    CPm=[zero,-1i*I;1i*(kX^2-epsMnm),zero];
    M_jm=k0*CPm; 
    %%%%%%%%%%%%%%%%%%%%%%%
    %T=(II-1/2*h*M_j1)\(II+1/2*h*M_j)*T;
    k1=h*M_j*T;
    k2=h*M_jm*(T+1/2*k1);
    k3=h*M_jm*(T+1/2*k2);
    k4=h*M_jm*(T+k3);
    T=T+1/6*k1+1/3*k2+1/3*k3+1/6*k4;  
    %%%%%%%%%%%%%%%%%%%%%%%
    M_j=M_j1;
    end
    T=W_Np1\T;
    
    Z_p1=(T(nDim+1:2*nDim,1:nDim)*R_ud+T(nDim+1:2*nDim,nDim+1:2*nDim))\I;
    
    R_ud=(T(1:nDim,nDim+1:2*nDim)+T(1:nDim,1:nDim)*R_ud)*Z_p1;                   %S+-, (7.153)
    T_dd=T_dd*Z_p1;                     %S--, (7.156)   
    
    %%%%%%%%%%%%%%%%%%%%%%%
   
    %for layer 1:N
    for iLy=2:N
         T=W_Np1;
         z_index0=(iLy-1)*Npoints+1;
    for z_index=z_index0+1:z_index0+Npoints
        %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index),epst(2,z_index),xt(:,z_index),nMax);
  
    CP=[zero,-1i*I;1i*(kX^2-epsMn),zero];
    M_j1=k0*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMnm = epsMatrix(epstm(1,z_index-1),epstm(2,z_index-1),xt(:,z_index-1),nMax);
    
    CPm=[zero,-1i*I;1i*(kX^2-epsMnm),zero];
    M_jm=k0*CPm; 
    %%%%%%%%%%%%%%%%%%%%%%%
    %T=(II-1/2*h*M_j1)\(II+1/2*h*M_j)*T;
    k1=h*M_j*T;
    k2=h*M_jm*(T+1/2*k1);
    k3=h*M_jm*(T+1/2*k2);
    k4=h*M_jm*(T+k3);
    T=T+1/6*k1+1/3*k2+1/3*k3+1/6*k4;  
    %%%%%%%%%%%%%%%%%%%%%%%
    M_j=M_j1;
    end
    T=W_Np1\T;
    
    Z_p1=I/(T(nDim+1:2*nDim,1:nDim)*R_ud+T(nDim+1:2*nDim,nDim+1:2*nDim));
    
    R_ud=(T(1:nDim,nDim+1:2*nDim)+T(1:nDim,1:nDim)*R_ud)*Z_p1;                   %S+-, (7.153)
    T_dd=T_dd*Z_p1;                     %S--, (7.156)   
    end
    
    
    %{
    %for the superstrate
    T=II;
    Z_p1=(T(nDim+1:2*nDim,1:nDim)*R_ud+T(nDim+1:2*nDim,nDim+1:2*nDim))\I;
    
    R_ud=(T(1:nDim,nDim+1:2*nDim)+T(1:nDim,1:nDim)*R_ud)*Z_p1;                   %S+-, (7.153)
    T_dd=T_dd*Z_p1;                     %S--, (7.156) 
    %}  
        

 
 %---------------------------------------------------------
    % R and T
 %---------------------------------------------------------
   % R and T
 %---------------------------------------------------------
        d_p1=zeros(nDim,1);
        d_p1(nMax+1,1)=1;
%---------------------------------------------------------

        RS=R_ud*d_p1;
        TS=T_dd*d_p1;  
        