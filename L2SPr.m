function [result]= L2SPr( K,y,alpha,beta,mu)
%initialization
n=length(unique(y));
maxIter=300;
nn=length(y);
Y1=zeros(nn);
Y2=zeros(nn);
Y3=zeros(nn);
Z=eye(nn);
%Z=rand(nn);
H=Z;

%main function
for i=1:maxIter
    Zold=Z;
    J=(K+mu*eye(nn))\(mu*Z+K-Y1);
     J(find(J<0))=0;
W=(2*alpha*K*H*H'*K'+mu*eye(nn))\(mu*Z-Y2+2*alpha*K*H*K');
 W(find(W<0))=0;
% J=J-diag(diag(J));

H=(2*alpha*K'*W*W'*K+mu*eye(nn))\(mu*Z-Y3+2*alpha*K'*W*K);
% W=W-diag(diag(W));
 H(find(H<0))=0;  
 
D=(W+Y1/mu+J+Y2/mu+H+Y3/mu)/3;
[U S V] = svd(D, 'econ');
       
        diagS = diag(S);
        svp = length(find(diagS > beta/(3*mu)));
        diagS = max(0,diagS - beta/(3*mu));
        
        if svp < 0.5 %svp = 0
            svp = 1;
        end
        Z = U(:,1:svp)*diag(diagS(1:svp))*V(:,1:svp)';
Z(find(Z<0))=0;

Y1=Y1-mu*(Z-J);
Y2=Y2-mu*(Z-W);
Y3=Y3-mu*(Z-H);
 mu=mu*1.1;

if((i>5)&(norm(Z-Zold,'fro') < norm(Zold,'fro') * 1e-5))  
        break
    end
end

L=(Z+Z')/2;
actual_ids = spectral_clustering(L, n);

result=ClusteringMeasure(actual_ids ,y);