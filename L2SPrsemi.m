function [meanacc]= L2SPrsemi(filename, K,y,alpha,beta,mu)
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
Z=max(abs(D)-beta/(mu*3),0).*sign(D);

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

%loop
for r = 0.1:0.2:0.5
 dlmwrite(filename,r,'-append','delimiter','\t','newline','pc');
    for it = 1:20

        [m,n]=size(K);
        c=length(unique(y)); % number of class
        numperc=floor(n/c); % number of data per class
        labelperc = floor(r*numperc); % number of labeled data per class
        labelindperc = sort(randperm(numperc,labelperc)); % index of labeled data selected
        labelind = []; % labelind: index of known label
        for i = 1:c
            labelind = [labelind labelindperc+(i-1)*numperc];
        end
        ulabel = setdiff(1:n,labelind);%index of unlabeled data
        %sync notation
        mm = length(labelind);
        labels = zeros(mm, c);
        for i = 1:mm
            labels(i,y(labelind(i))) = 1;
        end
        known_nodes = labelind;
        nodes_to_predict = ulabel;
        options.alpha = 0.99;
        %label
        if ~isfield(options, 'alpha')
            options.alpha = 0.99;
        end
        if ~isfield(options, 'precision')
            options.precision = 1e-5;
        end
        if ~isfield(options, 'maxiter')
            options.maxiter = 100;
        end

        [n, m] = size(L);

        if n ~= m
            error('label_diffusion:A_square', 'Adjacency matrix must be square');
        end
        if options.alpha <= 0 || options.alpha >= 1
            error('label_diffusion:alpha_range', 'alpha must belong to ]0,1[');
        end

        % Normalization of the adjacency matrix
        D = diag(sum(L,2).^(-0.5));
        diffusion_matrix = D*L*D;

        k = size(labels, 2);
        ini_scores = zeros(n, k);
        ini_scores(known_nodes, :) = labels;
        predictions_score = ini_scores;

        % Propagation
        for i=1:options.maxiter
            last_score = predictions_score;
            predictions_score = options.alpha*diffusion_matrix*predictions_score + (1-options.alpha)*ini_scores;
            if max(max(abs(last_score-predictions_score))) < options.precision
                break;
            end
        end

        % keep only predictions for nodes specified by 'nodes_to_predict'
        predictions_score = predictions_score(nodes_to_predict, :);

        % calculate accuracy
        [ur,uc]=size(ulabel);
        [max_value,max_ind] = max(predictions_score,[],2);
        cnt = 0;
        for i = 1:uc
            if max_ind(i) == y(ulabel(i))
                cnt = cnt+1;
            end
        end
        result(it) = cnt/uc;
        
        dlmwrite(filename,[result(it)],'-append','delimiter','\t','newline','pc');
    end
    meanacc=mean(result);
    stdacc=std(result);
    fprintf('%12.6f%12.6f\n',meanacc,stdacc)
    dlmwrite(filename,[alpha,beta,mu,meanacc,stdacc],'-append','delimiter','\t','newline','pc');
   
end

