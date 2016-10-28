clear;
sigmae=0.01;
for i=1:20
    Z(i)=randn(1)
end

for i=1:20
 for j=1:10
    X(i,j)= randn(1)*sigmae
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
% clear;
% sigmae=0.01;
% for i=1:20
%     Z(i)=randn(0,1)
% end
% 
% for i=1:20
%  for j=1:10
%     X(i,j)= mvnrnd(0,(sigmae*i))
%  end
% end
%      for i=1:20
%         for j=1:10
%             a(i,j)=Z(i)+X(i,j)
%         end
%     end
%     a
%     size(a)
%    a'
%     A=a';
% 
%     
% 

    

