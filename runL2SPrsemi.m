clear all
load YALE_165n_1024d_15c_uni;
load('data/YALE_165n_1024d_15c_zscore_uni_kernel/YALE_165n_1024d_15c_zscore_uni_kernel_gaussian_0.1_post_Sample-Scale.mat');
warning off

para1=[ 1e-5 1e-3 1e-2  .1 10  ];
para2=[1e-3 1e-2  .1 1];

for ij=1:length(para1)
lambda=para1(ij);
for iji=1:length(para2)
mu=para2(iji);
%for ji=1:length(para3)
%mu=para3(ji);

% [res] =localmulticluster(X,K,y,lgamma,1,mu(j))
fprintf('params%12.6f%12.6f%\n',lambda,mu)
%[result]=preserveclustering(K,D,y,gamma,mu)
%[result]=slkes(K,y,alpha,beta,mu)
[result]= L2SPrsemi('/Users/apple/Desktop/PreserveClusteringExperiment/result/slkersemi/yale/polyplus2.txt', K,y,alpha,beta,mu);
%dlmwrite('D:\kernelpreserve\result\slkes-ba\gaussian100.txt',[alpha,beta,mu,result],'-append','delimiter','\t','newline','pc');
    
end
end