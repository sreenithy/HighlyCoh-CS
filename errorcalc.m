function p=errorcalc(xest,indice,M)
indice
indexx=zeros(M,1);
indexx(indice,:)=1;
indexx
[a,b]=sort(abs(xest),'descend');
ind_xest=b(1:4);
ind_xest
indexxest=zeros(M,1);
indexxest(ind_xest)=1;
indexxest

p=(sum(abs(indexx-indexxest))>0);
p
end

