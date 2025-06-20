clc
clear

tic
maxiter=100;

r = 3;
T = 100;
N = 100;
T1 = T/2;
Break0=T1;
tau1=0.3;
tau2=0.7;

parfor  niter=1:maxiter
SSR_QML=[]
F=randn(T,r);           
L_1=normrnd(0,1,N,r)'; 
L_2=normrnd(0,1,N,r)';
E=randn(T,N);     

X=zeros(T,N);
X(1:T1,:)=F(1:T1,:)*L_1+E(1:T1,:);
X(T1+1:T,:)=F(T1+1:T,:)*L_2+E(T1+1:T,:);

T_min=T*tau1;
T_max=T*tau2;
[F_hat,L_hat,VNT]=panelFactorNew(X,r);

for k=T_min:T_max  
    
Sigma1=F_hat(1:k,:)'*F_hat(1:k,:)./k;
Sigma2=F_hat((1+k):T,:)'*F_hat((1+k):T,:)./(T-k);
SSR_QML(k-T_min+1)= k*log(det(Sigma1))+(T-k)*log(det(Sigma2));

end

Break_QML(niter,:)=find(SSR_QML==min(SSR_QML))-1+T_min;

end



hist(Break_QML)
title('QML')



toc