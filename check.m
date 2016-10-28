load Adct
Pi=A;
% Pi=Mat(1:30,1:300);
Pi = Pi*diag(1./sqrt(diag(Pi'*Pi)));

A=Pi'*Pi;
            %# Create a colored plot of the matrix values


h = imagesc(A)

impixelregion(h)
