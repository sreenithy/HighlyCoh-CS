clear;
Fe2=[0.01 0.05 0.1 0.2 0.4 0.6 0.8 0 1 2];
    
sigmae=0.01;
for i=1:20
    Z(i)=randn(10,1);
    
end

for i=1:20
 for j=1:10
    X(i,j)= randn(10,1)*sqrt(Se2);
 end
end
     for i=1:20
        for j=1:10
            a(i,j)=Z(i)+X(i,j)
        end
    end
    a
    size(a)
   a'
    A=a';
    size(A)
    [N,M] = size(A);

cohMatrix = zeros(M,M);
for i = 1 : M
    for j = i+1 : M
        cohMatrix(i,j) = A(:,i)'*A(:,j)/(norm(A(:,i)) * norm(A(:,j)));
    end
end

minCoh = min(min(cohMatrix));
maxCoh = max(max(cohMatrix))
cohmatrix = cohMatrix + cohMatrix' + eye(M);

    

