
path(path, './functions');

path(path, './ApproximateMessagePassing');

clear all;
load DicOversamp;               % load the dictionary matrix
Dtn = D;
[N,M] = size(Dtn);      % problem dimension
D = 4; % diversity (number of active sources)
D
SNR = [20 10 0];               % SNR
iterNum = 800;            % number of experiments
Phi = Dtn./(ones(N,1)*sqrt(sum(Dtn.^2)));    % dictionary matrix with columns draw uniformly from the surface of a unit hypersphere


for i=1:1
    i
    for it = 1 : iterNum
        nonzeroW = sign(randn(D,1)).* ( rand(D,1)*0.5 + 0.5 );      % nonzero Rows
        ind = randperm(M);                      % select active sources at random locations
        indice = ind(1:D);
        Wgen = zeros(M,1);
        Wgen(indice,:) = nonzeroW;
        signal = Phi * Wgen;                    % noiseless signal
        stdnoise = std(signal)*10^(-SNR(i)/20);    % observation noise
        noise = randn(N,1).*(ones(N,1)*stdnoise);
        T = signal + noise; % noisy signal
        %============================TMSBL================================
        
        %         X_tsbl = TMSBL(Phi, T, 'noise','mild');
        %
        %
        %      %   [F1] = calc(X_tsbl,indice,Wgen,'firstlargest',D);
        %
        % %         p(it)=errorcalc(X_tsbl,indice,M);
        % %        [Wgen, X_tsbl]
        % %                         %fail_TMSBL(it) = (F1~=1);
        %         mse_TMSBL(it) = (norm(Wgen - X_tsbl,'fro')/norm(Wgen,'fro'))^2;
        %
        %============================== ExCoV ========================
        %        X_excov = ExCoVapp(Phi,T,'Visibility',0);
        %         [F2] = perfSupp(X_excov,indice,'firstlargest', D);
        %         fail_EXCOV(it) = (F2~=1);
        %         mse_EXCOV(it) = (norm(Wgen - X_excov,'fro')/norm(Wgen,'fro'))^2;
        %     %============================== CoSaMP ========================
        %       X_cosamp = cosamp(Phi,T,D,1e-5);
        %       F3 = perfSupp(X_cosamp,indice,'firstlargest', D);
        %       fail_cosamp(it) = (F3~=1);
        %       mse_cosamp(it) = (norm(Wgen - X_cosamp,'fro')/norm(Wgen,'fro'))^2;
        %
        %============================== Subspace Pursuit ========================
        %       Xsub = CSRec_SP(D,Phi,T);
        %       F4 = perfSupp(Xsub.x_hat,indice,'firstlargest', D);
        %       fail_subspace(it) = (F4~=1);
        %       mse_subspace(it) = (norm(Wgen - Xsub.x_hat,'fro')/norm(Wgen,'fro'))^2;
        %============================== Approximate Message Passing ========================
                Xamp = reconstructAmp(Phi,T, 10);
                F5 = perfSupp(Xamp,indice,'firstlargest', D);
                fail_amp(it) = (F5~=1);
                mse_amp(it) = (norm(Wgen - Xamp,'fro')/norm(Wgen,'fro'))^2;
        %============================== Bayesian Compressive Sensing ========================
        %         initsigma2 = std(T)^2/1e2;
        %         [X_coef,used6] = BCS_fast_rvm(Phi,T,initsigma2,1e-8);
        %         X_bcs = zeros(M,1);
        %         X_bcs(used6) = X_coef;
        %         F6 = perfSupp(X_bcs,indice,'firstlargest', D);
        %         fail_bcs(it) = (F6~=1);
        %         mse_bcs(it) = (norm(Wgen - X_bcs,'fro')/norm(Wgen,'fro'))^2;
        %          %============================== EM-SBL ========================
        %                 ini_lambda = std(T)^2/1e2;
        %                 X_emsbl = EMSBL(Phi, T, ini_lambda, 1);
        %                 F7 = perfSupp(X_emsbl,indice,'firstlargest', D);
        %                 fail_emsbl(it) = (F7~=1);
        %                 mse_emsbl(it) = (norm(Wgen - X_emsbl,'fro')/norm(Wgen,'fro'))^2;
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
        
        %========================= Fast Bayesian Matching Pursuit (know SNR) ========================
        %
        %                 p1 = D/M;
        %                         sig2s = [0; 0.01];
        %                 mus = [0;0];
        %                 LP = 30;
        %                 variance=10^(-SNR(i)/10);
        %                 stop=0;
        %                 E=5;
        %                 r_stop=0.02;
        %                 %              [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max] = ...
        %                 %                     fbmpc_fxn_reduced(T, Phi, p1, variance, sig2s, mus, D);
        %                 %
        %             [xmmse, xmmse_star, psy_star, nu_star, T_star, d_tot, d_max, e_tot] = ...
        %                     fbmpr_gem_refine_fxn(T, Phi, p1, variance, sig2s, mus);
        %                 F11 = perfSupp(xmmse,indice,'firstlargest', D);
        %                 fail_fbmp(it) = (F11~=1);
        %                 mse_fbmp(it) = (norm(Wgen - xmmse,'fro')/norm(Wgen,'fro'))^2;
        %========================= FOCUSS (know SNR) ========================
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
        
    end
    
    %       mean_p_TMSBL(i)=sum(p);
    %   %  mean_fail_TMSBL(i)=mean(fail_TMSBL);
    %     mean_MSE_TMSBL(i)=mean(mse_TMSBL);
    %       mean_fail_EXCOV(i)=mean(fail_EXCOV);
    %       mean_MSE_EXCOV(i)=mean(mse_EXCOV);
    %       mean_fail_cosamp(i)=mean(fail_cosamp);
    %       mean_MSE_cosamp(i)=mean(mse_cosamp);
    %       mean_fail_subspace(i)=mean(fail_subspace);
    %      mean_mse_subspace(i)=mean(mse_subspace);
              mean_fail_amp(i)=mean(fail_amp);
              mean_mse_amp(i)=mean(mse_amp);
    %       mean_fail_bcs(i)=mean(fail_bcs);
    %       mean_mse_bcs(i)=mean(mse_bcs);
    %       mean_fail_emsbl(i)=mean(fail_emsbl);
    %mean_mse_emsbl(i)=mean(mse_emsbl);
    %                mean_fail_l1(i)=mean(fail_l1);
    %                mean_mse_l1(i)=mean(mse_l1);
    %           mean_fail_htp(i)=mean(fail_htp);
    %           mean_mse_htp(i)=mean(mse_htp);
    %           mean_fail_fbmp(i)=mean(fail_fbmp);
    %           mean_mse_fbmp(i)=mean(mse_fbmp);
    % mean_fail_focuss(i)=mean(fail_focuss);
    % mean_mse_focuss(i)=mean(mse_focuss);
    %           mean_fail_sl0(i)=mean(fail_sl0);
    %           mean_mse_sl0(i)=mean(mse_sl0);
end

% mean_p_TMSBL
% %mean_fail_TMSBL
% mean_MSE_TMSBL
% mean_fail_EXCOV
% mean_MSE_EXCOV
% mean_fail_cosamp
% mean_MSE_cosamp
% mean_fail_subspace
% mean_mse_subspace
mean_fail_amp
mean_mse_amp
% mean_fail_htp
% mean_mse_htp
% mean_fail_l1
% mean_mse_l1
% mean_fail_emsbl
% mean_mse_emsbl
% mean_fail_fbmp
% mean_mse_fbmp
% mean_fail_focuss
% mean_mse_focuss
% mean_fail_sl0
% mean_mse_sl0





