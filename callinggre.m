X = load('dataXforAR1.mat');
  opt = struct('K',32, 'samet','mexomp','saopt',struct('tnz',4));
  D=4;
  M=100;
  nonzeroW = sign(randn(D,1)).* ( rand(D,1)*0.5 + 0.5 );      % nonzero Rows
        ind = randperm(M);                      % select active sources at random locations
        indice = ind(1:D);
        Wgen = zeros(M,1);
        Wgen(indice,:) = nonzeroW;
        X=Wgen;
        
  Ds = dictlearn_mb('X',X, opt);    
  figure(1);clf;plot(Ds.ptab(:,1),Ds.ptab(:,2),'b-');xlabel(Ds.ptc{1});ylabel(Ds.ptc{2}); 
