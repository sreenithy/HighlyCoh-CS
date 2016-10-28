
path(path, './functions');

path(path, './ApproximateMessagePassing');

clear all;
  % problem dimension
D = 4; % diversity (number of active sources)
D
SNR = 20;               % SNR
Fe=[0.1 0.05 0.01 0.005 0.001];
for i=1
    
    i
    Se2=Fe(i);
    for it = 1 : 1
        it
        a1=randn(10,5);
        a2=randn(10,5);
        a3=randn(10,5);
        a4=randn(10,5);
        b1=randn(10,1);
        
       a11=randn(10,1)*sqrt(Se2);
       a12=randn(10,1)*sqrt(Se2);
       a13=randn(10,1)*sqrt(Se2);
       a14=randn(10,1)*sqrt(Se2);
       a15=randn(10,1)*sqrt(Se2);
       A1=[a11 a12 a13 a14 a15];
%        Phi=[A1 A2 A3 A4];
%        
%        
%        y21=randn(10,1)*sqrt(Se2);
%        y22=randn(10,1)*sqrt(Se2);
%        y23=randn(10,1)*sqrt(Se2);
%        y24=randn(10,1)*sqrt(Se2);
%        y25=randn(10,1)*sqrt(Se2);
%        y31=randn(10,1)*sqrt(Se2);
%        y32=randn(10,1)*sqrt(Se2);
%        y33=randn(10,1)*sqrt(Se2);
%        y34=randn(10,1)*sqrt(Se2);
%        y35=randn(10,1)*sqrt(Se2);
%        y41=randn(10,1)*sqrt(Se2);
%        y42=randn(10,1)*sqrt(Se2);
%        y43=randn(10,1)*sqrt(Se2);
%        y44=randn(10,1)*sqrt(Se2);
%        y45=randn(10,1)*sqrt(Se2);
%        x11=z1+y11;
%        x12=z1+y12;
%        x13=z1+y13;
%        x14=z1+y14;
%        x15=z1+y15;
%        
%        x21=z2+y21;
%        x22=z2+y22;
%        x23=z2+y23;
%        x24=z2+y24;
%        x25=z2+y25;
%        
%        x31=z2+y31;
%        x32=z2+y32;
%        x33=z2+y33;
%        x34=z2+y34;
%        x35=z2+y35;
%        
%         x41=z2+y41;
%        x42=z2+y42;
%        x43=z2+y43;
%        x44=z2+y44;
%        x45=z2+y45;
%        
%       pi= [x11 x12 x13 x14 x15 x21 x22 x23 x24 x25 x31 x32 x33 x34 x35 x41 x42 x43 x44 x45];
%         pi = pi*diag(1./sqrt(diag(pi'*pi)));  
%         Phi=pi;
%         [N M]=size(Phi);
%         nonzeroW = sign(randn(D,1)).* ( rand(D,1)*0.5 + 0.5 );      % nonzero Rows
%         %ind = randperm(M);                      % select active sources at random locations
%         indice = [1 2 3 4];
%         Wgen = zeros(M,1);
%         Wgen(indice,:) = nonzeroW;
%         Wgen
%         signal = Phi * Wgen;                    % noiseless signal
%         stdnoise = std(signal)*10^(-SNR/20);    % observation noise
%         noise = randn(N,1).*(ones(N,1)*stdnoise);
%         T = signal + noise; % noisy signal
%         %============================TMSBL================================
%         X_tsbl = TMSBL(Phi, T, 'noise','mild');
%         [X_tsbl, Wgen];
% %                         ind_x=find(abs(Wgen)>0);
% %         index_x=zeros(20,1);
% %         index_x(ind_x)=1;
% %         
% %         [a,b]=sort(abs(X_tsbl),'descend');
% %         ind_est=b(1:4);
% %         index_est=zeros(20,1);
% %         index_est(ind_est)=1;
% %         
% %         error(it)=1*(sum(abs(index_x-index_est))>0)
% %      %  mse_TMSBL(it) = (norm(Wgen - X_tsbl,'fro')/norm(Wgen,'fro'))^2;
% %         
% %                % ============================== ExCoV ========================
% %                 X_excov = ExCoVapp(Phi,T,'Visibility',0);
% %         
% %         ind_x=find(abs(Wgen)>0);
% %         index_x=zeros(20,1);
% %         index_x(ind_x)=1;
% %         
% %         [a,b]=sort(abs(X_excov),'descend');
% %         ind_est=b(1:4);
% %         index_est=zeros(20,1);
% %         index_est(ind_est)=1;
% %         
% %         error(it)=1*(sum(abs(index_x-index_est))>0)
% %        % mse_EXCOV(it) = (norm(Wgen - X_excov,'fro')/norm(Wgen,'fro'))^2;
% %               
%         %%============================== OMP ================
% %                 [x,r,normR,residHist, errHist] = OMP( Phi, T, D,[] )
% %                 X_omp=x;
% %         ind_x=find(abs(Wgen)>0);
% %         index_x=zeros(20,1);
% %         index_x(ind_x)=1;
% %         
% %         [a,b]=sort(abs(X_omp),'descend');
% %         ind_est=b(1:4);
% %         index_est=zeros(20,1);
% %         index_est(ind_est)=1;
% %         
% %         error(it)=1*(sum(abs(index_x-index_est))>0)
% %          mse_omp(it) = (norm(Wgen - X_omp,'fro')/norm(Wgen,'fro'))^2;
% %         
% %          
% %               %============================== Bayesian Compressive Sensing ========================
% %         initsigma2 = std(T)^2/1e2;
% %         [X_coef,used6] = BCS_fast_rvm(Phi,T,initsigma2,1e-8);
% %         X_bcs = zeros(M,1);
% %         X_bcs(used6) = X_coef;
% %         
% %         ind_x=find(abs(Wgen)>0);
% %         index_x=zeros(20,1);
% %         index_x(ind_x)=1;
% %         
% %         [a,b]=sort(abs(X_bcs),'descend');
% %         ind_est=b(1:4);
% %         index_est=zeros(20,1);
% %         index_est(ind_est)=1;
% %         
% %         error(it)=1*(sum(abs(index_x-index_est))>0)
%       %  mse_bcs(it) = (norm(Wgen - X_bcs,'fro')/norm(Wgen,'fro'))^2;
%        
%         % ========================= L1_LS ==============
% %                 [x,status,history] = l1_ls(Phi,T,0.01);
% %                 X_l1ls=x;
% %         ind_x=find(abs(Wgen)>0);
% %         index_x=zeros(20,1);
% %         index_x(ind_x)=1;
% %         
% %         [a,b]=sort(abs(X_l1ls),'descend');
% %         ind_est=b(1:4);
% %         index_est=zeros(20,1);
% %         index_est(ind_est)=1;
% %         
% %         error(it)=1*(sum(abs(index_x-index_est))>0);
%       % mse_l1ls(it) = (norm(Wgen - X_l1ls,'fro')/norm(Wgen,'fro'))^2;
%         
%     end
%     p(i)=sum(error)/1000;
%     
%    end
% 
%  
%  p
% % Fe
% % figure(1);
% % plot(Fe,p,'b--o');
% % xlabel('Variance');
% % ylabel('P(e)');
% % title('Performance of L1LS');
% 
