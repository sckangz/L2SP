clear all
load YALE_165n_1024d_15c_uni;
load('data/YALE_165n_1024d_15c_zscore_uni_kernel/YALE_165n_1024d_15c_zscore_uni_kernel_gaussian_0.1_post_Sample-Scale.mat');
warning off

para1=[ 1e-5 1e-3 1e-2  .1 10  ];
para2=[1e-3 1e-2  .1 1];
para3=[1e-3 1e-2  .1 1];


for ij=1:length(para1)
alpha=para1(ij);
for iji=1:length(para2)
beta=para2(iji);
for ji=1:length(para3)
mu=para3(ji);


fprintf('params%12.6f%12.6f%\n',alpha,beta,mu)

[result]=L2SPs(K,y,alpha,beta,mu)
dlmwrite('',[alpha,beta,mu,result],'-append','delimiter','\t','newline','pc');
end
end
end