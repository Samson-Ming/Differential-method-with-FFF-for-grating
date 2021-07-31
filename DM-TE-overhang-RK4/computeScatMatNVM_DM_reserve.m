function [RS,s0V,TS,sSubV,T_dd,R_ud] =computeScatMatNVM_DM_reserve(lam,thI,epsB,Lam,epsS,xt,zV,epst,nMax,N,Npoints)      
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

xt=fliplr(xt);
zV=fliplr(zV);
epst=fliplr(epst);
%% W0 and W_Np1 
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

%% S0    
    R_ud=zeros(nDim);
    R_du=zeros(nDim);
    T_uu=I;
    T_dd=I;
 %---------------------------------------------------------  
    
    %for layer 1, bourndary 1, i.e. z=0;
    % in fact bounary(1)-  and boundary(1)+
    h=abs(zV(2)-zV(1));
    z_index0=(1-1)*Npoints+1;
    %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index0),epst(2,z_index0),xt(:,z_index0),nMax);
    %etaMn = epsMatrix(1/epst(1,z_index0),1/epst(2,z_index0),xt(:,z_index0),nMax);
    %[cosM,sinM] = generateSinCosMat(sinx(:,z_index0),cosx(:,z_index0),x_cs(z_index0),nMax); 
   
    %[A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    %CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    CP=[zero,-1i*I;1i*(kX^2-epsMn),zero];
    M_j=k0*CP;
    
    T=W_Np1;
    for z_index=z_index0+1:z_index0+Npoints
        %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index),epst(2,z_index),xt(:,z_index),nMax);
    %etaMn = epsMatrix(1/epst(1,z_index),1/epst(2,z_index),xt(:,z_index),nMax);
    %[cosM,sinM] = generateSinCosMat(sinx(:,z_index),cosx(:,z_index),x_cs(z_index),nMax); 
   
    %[A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    %CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    %CP=-[B*D*kX,(B*D*B-A-etaMn\I);(kX*D*kX-I),kX*D*B];
    CP=[zero,-1i*I;1i*(kX^2-epsMn),zero];
    M_j1=k0*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    T=(II-1/2*h*M_j1)\(II+1/2*h*M_j)*T;
    M_j=M_j1;
    end
    T=W_0\T;
    
    Z_p1=(T(1:nDim,1:nDim)-R_ud*T(nDim+1:2*nDim,1:nDim))\I;
    
    R_ud=-Z_p1*(T(1:nDim,nDim+1:2*nDim)-R_ud*T(1:nDim,1:nDim));                   %S+-, (7.166)
    T_dd=T_dd*(T(nDim+1:2*nDim,nDim+1:2*nDim)+T(nDim+1:2*nDim,1:nDim)*R_ud);                     %S--, (7.171)   
    
    %%%%%%%%%%%%%%%%%%%%%%%
   
    %for layer 1:N
    for iLy=2:N
         T=W_Np1;
         z_index0=(iLy-1)*Npoints+1;
    for z_index=z_index0+1:z_index0+Npoints
        %%%%%%%%%%%%%%%%%%%%%%%
    %M matrix at each z;
    epsMn = epsMatrix(epst(1,z_index),epst(2,z_index),xt(:,z_index),nMax);
    %etaMn = epsMatrix(1/epst(1,z_index),1/epst(2,z_index),xt(:,z_index),nMax);
    %[cosM,sinM] = generateSinCosMat(sinx(:,z_index),cosx(:,z_index),x_cs(z_index),nMax); 
   
    %[A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    %CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];
    %CP=-[B*D*kX,(B*D*B-A-etaMn\I);(kX*D*kX-I),kX*D*B];
    CP=[zero,-1i*I;1i*(kX^2-epsMn),zero];
    
    M_j1=k0*CP;   
    %%%%%%%%%%%%%%%%%%%%%%%
    T=(II-1/2*h*M_j1)\(II+1/2*h*M_j)*T;
    M_j=M_j1;
    end
    T=W_0\T;
    
    Z_p1=(T(1:nDim,1:nDim)-R_ud*T(nDim+1:2*nDim,1:nDim))\I;
    
    R_ud=-Z_p1*(T(1:nDim,nDim+1:2*nDim)-R_ud*T(1:nDim,1:nDim));                   %S+-, (7.166)
    T_dd=T_dd*(T(nDim+1:2*nDim,nDim+1:2*nDim)+T(nDim+1:2*nDim,1:nDim)*R_ud);                     %S--, (7.171)   
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
        