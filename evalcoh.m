function coh=evalcoh(D)
[N,M] = size(D);

cohMatrix = zeros(M,M);
for i = 1 : M
    for j = i+1 : M
        cohMatrix(i,j) = D(:,i)'*D(:,j)/(norm(D(:,i)) * norm(D(:,j)));
    end
end

minCoh = min(min(cohMatrix));

cohmatrix = cohMatrix + cohMatrix' + eye(M);
maxCoh = max(max(cohMatrix))
return 