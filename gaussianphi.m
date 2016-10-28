A=eye(10)
for i=1:10
    for j=1:10
        if i~=j
            A(i,j)=0.9;
        end
    end
end
    A;
    covsigma=blkdiag(A,A,A,A,A,A,A,A,A,A)
    size(covsigma);

    for i=1:100
    for j=1:100
        if(i==j)
                     mu(i,j)=0;
        end
       
    end
    end
size(mu)
    R2=mvnrnd(mu,covsigma,100);
    size(R2)
[N,M] = size(R2);

cohMatrix = zeros(M,M);
for i = 1 : M
    for j = i+1 : M
        cohMatrix(i,j) = R2(:,i)'*R2(:,j)/(norm(R2(:,i)) * norm(R2(:,j)));
    end
end

minCoh = min(min(cohMatrix));

cohmatrix = cohMatrix + cohMatrix' + eye(M);
maxCoh = max(max(cohMatrix))

