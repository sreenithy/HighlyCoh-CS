
path(path, './functions');

path(path, './ApproximateMessagePassing');

clear all;
load GeneratedDic;               % load the dictionary matrix
Dtn = A;
[N,M] = size(Dtn)    % problem dimension
D = 4; % diversity (number of active sources)
D
SNR = 20;               % SNR
iterNum = 1000;            % number of experiments
Phi = Dtn./(ones(N,1)*sqrt(sum(Dtn.^2)));    % dictionary matrix with columns draw uniformly from the surface of a unit hypersphere

Fe=[0.01 0.03 0.05 0.1 0.2 0.4 0.6 0.8 1 2];
for i=1:10
    i
    Se2=Fe(i);
    for it = 1 : 1000
        it
        z1=randn(10,1);
        z2=randn(10,1);
        z3=randn(10,1);
        z4=randn(10,1);
       y11=randn(10,1)*sqrt(Se2);
       y12=randn(10,1)*sqrt(Se2);
       y13=randn(10,1)*sqrt(Se2);
       y14=randn(10,1)*sqrt(Se2);
       y15=randn(10,1)*sqrt(Se2);
       y21=randn(10,1)*sqrt(Se2);
       y22=randn(10,1)*sqrt(Se2);
       y23=randn(10,1)*sqrt(Se2);
       y24=randn(10,1)*sqrt(Se2);
       y25=randn(10,1)*sqrt(Se2);
       y31=randn(10,1)*sqrt(Se2);
       y32=randn(10,1)*sqrt(Se2);
       y33=randn(10,1)*sqrt(Se2);
       y34=randn(10,1)*sqrt(Se2);
       y35=randn(10,1)*sqrt(Se2);
       y41=randn(10,1)*sqrt(Se2);
       y42=randn(10,1)*sqrt(Se2);
       y43=randn(10,1)*sqrt(Se2);
       y44=randn(10,1)*sqrt(Se2);
       y45=randn(10,1)*sqrt(Se2);
       x11=z1+y11;
       x12=z1+y12;
       x13=z1+y13;
       x14=z1+y14;
       x15=z1+y15;
       
       x21=z2+y21;
       x22=z2+y22;
       x23=z2+y23;
       x24=z2+y24;
       x25=z2+y25;
       
       x31=z2+y31;
       x32=z2+y32;
       x33=z2+y33;
       x34=z2+y34;
       x35=z2+y35;
       
        x41=z2+y41;
       x42=z2+y42;
       x43=z2+y43;
       x44=z2+y44;
       x45=z2+y45;
       
      pi= [x11 x12 x13 x14 x15 x21 x22 x23 x24 x25 x31 x32 x33 x34 x35 x41 x42 x43 x44 x45];
        pi = pi*diag(1./sqrt(diag(pi'*pi)));  
        Phi=pi;
        [N M]=size(Phi);
        nonzeroW = sign(randn(D,1)).* ( rand(D,1)*0.5 + 0.5 );      % nonzero Rows
        ind = randperm(M);                      % select active sources at random locations
        indice = ind(1:D);
        Wgen = zeros(M,1);
        Wgen(indice,:) = nonzeroW;
        signal = Phi * Wgen;                    % noiseless signal
        stdnoise = std(signal)*10^(-SNR/20);    % observation noise
        noise = randn(N,1).*(ones(N,1)*stdnoise);
        T = signal + noise; % noisy signal
        %============================TMSBL================================
%         X_tsbl = TMSBL(Phi, T, 'noise','mild');
%         [X_tsbl, Wgen];
%         ind_x=find(abs(Wgen)>0);
%         index_x=zeros(20,1);
%         index_x(ind_x)=1;
%         
%         [a,b]=sort(abs(X_tsbl),'descend');
%         ind_est=b(1:4);
%         index_est=zeros(20,1);
%         index_est(ind_est)=1;
%         
%         error(it)=1*(sum(abs(index_x-index_est))>0)
%         mse_TMSBL(it) = (norm(Wgen - X_tsbl,'fro')/norm(Wgen,'fro'))^2;
        
        %        % ============================== ExCoV ========================
%                 X_excov = ExCoVapp(Phi,T,'Visibility',0);
%         
%         ind_x=find(abs(Wgen)>0);
%         index_x=zeros(20,1);
%         index_x(ind_x)=1;
%         
%         [a,b]=sort(abs(X_excov),'descend');
%         ind_est=b(1:4);
%         index_est=zeros(20,1);
%         index_est(ind_est)=1;
%         
%         error(it)=1*(sum(abs(index_x-index_est))>0)
%         mse_EXCOV(it) = (norm(Wgen - X_excov,'fro')/norm(Wgen,'fro'))^2;
        %             %============================== CoSaMP ========================
        %               X_cosamp = cosamp(Phi,T,D,1e-5);
        %               F3 = perfSupp(X_cosamp,indice,'firstlargest', D);
        %               fail_cosamp(it) = (F3~=1);
        %               mse_cosamp(it) = (norm(Wgen - X_cosamp,'fro')/norm(Wgen,'fro'))^2;
        %%============================== OMP ================
                [x,r,normR,residHist, errHist] = OMP( Phi, T, D,[] )
                X_omp=x;
        ind_x=find(abs(Wgen)>0);
        index_x=zeros(20,1);
        index_x(ind_x)=1;
        
        [a,b]=sort(abs(X_omp),'descend');
        ind_est=b(1:4);
        index_est=zeros(20,1);
        index_est(ind_est)=1;
        
        error(it)=1*(sum(abs(index_x-index_est))>0)
         mse_omp(it) = (norm(Wgen - X_omp,'fro')/norm(Wgen,'fro'))^2;
        
        %============================== Subspace Pursuit ========================
        %               Xsub = CSRec_SP(D,Phi,T);
        %               F4 = perfSupp(Xsub.x_hat,indice,'firstlargest', D);
        %               fail_subspace(it) = (F4~=1);
        %               mse_subspace(it) = (norm(Wgen - Xsub.x_hat,'fro')/norm(Wgen,'fro'))^2;
        %============================== Approximate Message Passing ========================
        %         Xamp = reconstructAmp(Phi,T, 10);
        %         F5 = perfSupp(Xamp,indice,'firstlargest', D);
        %         fail_amp(it) = (F5~=1);
        %         mse_amp(it) = (norm(Wgen - Xamp,'fro')/norm(Wgen,'fro'))^2;
        %         %============================== Bayesian Compressive Sensing ========================
%         initsigma2 = std(T)^2/1e2;
%         [X_coef,used6] = BCS_fast_rvm(Phi,T,initsigma2,1e-8);
%         X_bcs = zeros(M,1);
%         X_bcs(used6) = X_coef;
%         
%         ind_x=find(abs(Wgen)>0);
%         index_x=zeros(20,1);
%         index_x(ind_x)=1;
%         
%         [a,b]=sort(abs(X_bcs),'descend');
%         ind_est=b(1:4);
%         index_est=zeros(20,1);
%         index_est(ind_est)=1;
%         
%         error(it)=1*(sum(abs(index_x-index_est))>0)
%         mse_bcs(it) = (norm(Wgen - X_bcs,'fro')/norm(Wgen,'fro'))^2;
        %         %          %============================== EM-SBL ========================
        %                         ini_lambda = std(T)^2/1e2;
        %                         X_emsbl = EMSBL(Phi, T, ini_lambda, 1);
        %                         F7 = perfSupp(X_emsbl,indice,'firstlargest', D);
        %                         fail_emsbl(it) = (F7~=1);
        %                         mse_emsbl(it) = (norm(Wgen - X_emsbl,'fro')/norm(Wgen,'fro'))^2;
        %============================== Magic L1 ========================
        %         path(path, './l1magic-1.1');
        %         x0 = Phi'*T;
        %         epsilon =  stdnoise*sqrt(D)*sqrt(1 + 2*sqrt(2)/sqrt(D));
        %         X_l1 = l1qc_logbarrier(x0, Phi, [], T, epsilon, 1e-3);
        %         X_l1=l1magic(Phi, T, Wgen);
        % %         F8 = perfSupp(X_l1,indice,'firstlargest', D);
        % %         fail_l1(it) = (F8~=1);
        %         mse_l1(it) = (norm(Wgen - X_l1,'fro')/norm(Wgen,'fro'))^2;
        %============================== Hard Thresholding Pursuit algorithm ========================
        %                   % implemented as a noiseless case (have removed the noise)
        %                   [X_htp] = HTP(signal,Phi,D);
        %                   F10 = perfSupp(X_htp,indice,'firstlargest', D);
        %                   fail_htp(it) = (F10~=1);
        %                   mse_htp(it) = (norm(Wgen - X_htp,'fro')/norm(Wgen,'fro'))^2;
        
        %         %========================= Fast Bayesian Matching Pursuit (know SNR) ========================
        %
        %                         p1 = D/M;
        %                                 sig2s = [0; 0.01];
        %                         mus = [0;0];
        %                         variance=10^(-SNR(i)/10);
        %                         Phi = Phi*diag(1./sqrt(diag(Phi'*Phi)));
        %
        %                         %              [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max] = ...
        %                         %                     fbmpc_fxn_reduced(T, Phi, p1, variance, sig2s, mus, D);
        %                         %
        %                     [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max, e_tot] = ...
        %                             fbmpr_gem_refine_fxn(T, Phi, p1, variance, sig2s, mus);
        %                         F11 = perfSupp(xmmse,indice,'firstlargest', D);
        %                         fail_fbmp(it) = (F11~=1);
        %                         mse_fbmp(it) = (norm(Wgen - xmmse,'fro')/norm(Wgen,'fro'))^2;
        %         %========================= FOCUSS (know SNR) ========================
        %                         variance=10^(-SNR(i)/10);
        %                         [X_focuss] = MFOCUSS(Phi, T,variance);
        %                         %F12 = perfSupp(X_focuss,indice,'firstlargest', D);
        %                         %fail_focuss(it) = (F12~=1);
        %                         mse_focuss(it) = (norm(Wgen - X_focuss,'fro')/norm(Wgen,'fro'))^2;
        % ========================= SL0 (removed noise) ==============
        %                    sigma_min = 0.001;
        %             sigma_decrease_factor = 0.5;
        %             X_sl0 = SL0(Phi, signal, sigma_min, sigma_decrease_factor);
        %            F13 = perfSupp(X_sl0,indice,'firstlargest', D);
        %            fail_sl0(it) = (F13~=1);
        %          mse_sl0(it) = (norm(Wgen - X_sl0,'fro')/norm(Wgen,'fro'))^2;
        
        
        % ========================= L1_LS ==============
%                 [x,status,history] = l1_ls(Phi,T,0.01);
%                 X_l1ls=x;
%         ind_x=find(abs(Wgen)>0);
%         index_x=zeros(20,1);
%         index_x(ind_x)=1;
%         
%         [a,b]=sort(abs(X_l1ls),'descend');
%         ind_est=b(1:4);
%         index_est=zeros(20,1);
%         index_est(ind_est)=1;
%         
%         error(it)=1*(sum(abs(index_x-index_est))>0)
%         mse_l1ls(it) = (norm(Wgen - X_l1ls,'fro')/norm(Wgen,'fro'))^2;
    end
    p(i)=sum(error)/1000
    %              mean_p_TMSBL(i)=sum(p);
    %           mean_fail_TMSBL(i)=mean(fail_TMSBL);
    %  mean_MSE_TMSBL(i)=mean(mse_TMSBL);
    %               mean_fail_EXCOV(i)=mean(fail_EXCOV);
    %               mean_MSE_EXCOV(i)=mean(mse_EXCOV);
    %       mean_fail_cosamp(i)=mean(fail_cosamp);
    %           mean_MSE_cosamp(i)=mean(mse_cosamp);
    %       mean_fail_subspace(i)=mean(fail_subspace);
    %      mean_mse_subspace(i)=mean(mse_subspace);
    %     mean_fail_amp(i)=mean(fail_amp);
    %     mean_mse_amp(i)=mean(mse_amp);
    %   mean_fail_bcs(i)=mean(fail_bcs);
    %   mean_mse_bcs(i)=mean(mse_bcs);
    %   mean_fail_emsbl(i)=mean(fail_emsbl);
    %   mean_mse_emsbl(i)=mean(mse_emsbl);
    %  mean_fail_l1(i)=mean(fail_l1);
    %  mean_mse_l1(i)=mean(mse_l1);
    %  mean_fail_htp(i)=mean(fail_htp);
    % mean_mse_htp(i)=mean(mse_htp);
    %  mean_fail_fbmp(i)=mean(fail_fbmp);
    %  mean_mse_fbmp(i)=mean(mse_fbmp);
    % mean_fail_focuss(i)=mean(fail_focuss);
    % mean_mse_focuss(i)=mean(mse_focuss);
    % mean_fail_sl0(i)=mean(fail_sl0);
    %  mean_mse_sl0(i)=mean(mse_sl0);
    %  mean_fail_omp(i)=mean(fail_omp);
    %  mean_mse_omp(i)=mean(mse_omp);
    %  mean_fail_l1ls(i)=mean(fail_l1ls);
    %  mean_mse_l1ls(i)=mean(mse_l1ls);
end



%mean_fail_TMSBL
%mean_MSE_TMSBL
% mean_fail_EXCOV
% mean_MSE_EXCOV
%mean_fail_cosamp
% mean_MSE_cosamp
% mean_fail_subspace
% mean_mse_subspace
% mean_fail_amp
% mean_mse_amp
% mean_fail_htp
% mean_mse_htp
% mean_fail_l1
%mean_mse_l1
% mean_fail_emsbl
% mean_mse_emsbl
% mean_fail_fbmp
% mean_mse_fbmp
% mean_fail_focuss
% mean_mse_focuss
% mean_fail_sl0
% mean_mse_sl0
% mean_fail_bcs
% mean_mse_bcs
% mean_fail_l1ls
% mean_mse_l1ls
%

p
Fe
figure(1);
plot(Fe,p,'b--o');
xlabel('Variance');
ylabel('P(e)');
title('Performance of OMP');

