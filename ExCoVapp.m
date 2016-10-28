function [s_hat,sigma2_hat,A_index,Count]=ExCoVapp(H,y,varargin)
%Approx. ExCoV routine: 
%Coded by Kun Qiu (kqiu@iastate.edu)
%updated Feb 14, 2009

%Function usage
%=======================================================
%INPUT(compulsory):
%H:                      the sensing matrix
%y:                       the measurement column vector
%
%INPUT(optional):
%'Noi_Cov':        the N by N noise covariance structure matrix
%                         (default=eye(N))
%'FIR_len':         the size of moving average window
%                         (default=10)
%'Init_mode':      type of initial sparsity support (valued in {0,1})
%                         0: initialize with empty set 
%                         1: initialize with a crude high level signal components support estimate
%                         (default=1)
%'Visibility':         Option to see visually the reconstrution process (valued in {0,1})
%                         0: work silently
%                         1: work openly
%                         (default=1)
%========================================================
%OUTPUT:
%s_hat:              the signal esitmate
%sigma2_hat:    the noise variance estimate
%A_index:         the estimated set of high level signal components
%Count:             Count of number of iterations
%========================================================

if (nargin-length(varargin))~=2
    error('Missing required inputs!');
end

[N,m]=size(H);

%Setting default values for the optional inputs
C=eye(N);
FIR_len=10;
Init_mode=1;
Visibility=1;

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('Noi_Cov')
                C=varargin{i+1};
            case upper('FIR_len')
                FIR_len=varargin{i+1};
            case upper('Init_mode')
                Init_mode=varargin{i+1};
            case upper('Visibility')
                Visibility=varargin{i+1};     
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end

Inv_HH=inv(H*H');
C_inv=inv(C);

scale=sqrt(N/(y'*C_inv*y));
y=scale*y;

%Initialization
if Init_mode
    ma=floor(0.5*N/log(m/N));
else
    ma=1;
end
mb=m-ma;
u=H'*Inv_HH*y;
u_sort=sort(abs(u),'descend');
A_index=find(abs(u)>=u_sort(ma))';
B_index=1:m;
B_index(A_index)=[];
HA=H(:,A_index);
HB=H(:,B_index);
s_A_init=u(A_index);
sigma2=(y-HA*s_A_init)'*C_inv*(y-HA*s_A_init)/N;
delta2_a=10*sigma2*ones(ma,1)./diag(HA'*HA);

expand=1;
compress=0;
exit_flag=0;
Count=0;
p=0;
A_index_record_pre=[];
GML_star=-inf;
GML_record=-inf;

if Visibility
    figure(1)
    set(gcf,'color','w');
    subplot(2,1,1);
    plot(-inf,-inf);
    axis([0,N,0,500]);
    xlabel('Sparsity','Fontsize',12);
    ylabel('GL','Fontsize',12);
    title('Now expanding...','Fontsize',16);
    box off
    hold on
    subplot(2,1,2);
    axis([1,m,-1.5,+1.5]);
    title('ExCoV is estimating the signal...','Fontsize',16);
    xlabel('signal index','Fontsize',12);
    ylabel('signal value','Fontsize',12);
    box off
end

%Main iteration
while ~exit_flag
    Count=Count+1;
    p=p+1;
    
    %One EM step
    s_A=inv(HA'*HA+diag(sigma2./delta2_a))*HA'*y;
    upsilon=(y-HA*s_A)/sigma2;
    s_B=HB'*upsilon;
    sigma2=norm(y-HA*s_A)^2/N;
    HAtHAdiag=diag(HA'*HA);    
    delta2_a=max(s_A.^2,sigma2/10./HAtHAdiag);
    
    GML(p)=real(0.5*(-log(0.5*(N-ma))-(N-ma-2)*log(sigma2)-y'*(y-HA*s_A)/sigma2-sum(log((HAtHAdiag.^2)/2./(sigma2+HAtHAdiag.*(delta2_a))))));    

    if p>=FIR_len
        MvAvgGML(p)=mean(GML(p-FIR_len+1:p));
    else
        MvAvgGML(p)=-inf;
    end

    if GML(p)>GML_star
        GML_star=GML(p);
        A_index_star=A_index;
        B_index_star=B_index;
        s_A_star=s_A;
        sigma2_star=sigma2;
        ma_star=ma;
        mb_star=mb;
    end
    
    if Visibility
        figure(1)
        subplot(2,1,2);
        stem(A_index,s_A/scale,'Marker','none','linewidth',2,'color','k');
        axis([1,m,-1.5,+1.5]);
        xlabel('signal index','Fontsize',12);
        ylabel('signal value','Fontsize',12);
        title('ExCoV is estimating the signal...','Fontsize',16);
        box off
    end

    GML_plot=GML(p);
    if expand
        if Visibility
            figure(1)
            subplot(2,1,1);
            plot_index=ma;
            plot(plot_index,GML_plot,'ko');
            axis([0,N,0,500]);
            title('Expanding...','Fontsize',16);
        end
        if p<=FIR_len|MvAvgGML(p)>=MvAvgGML(p-1)|GML(p)>=MvAvgGML(p-1)
            %Expansion
            k=find(abs(s_B)==max(abs(s_B)),1,'first');
            A_index=[A_index,B_index(k)];
            B_index(k)=[];
            delta2_a=[delta2_a;sigma2/(norm(HB(:,k))^2)];
            HA=[HA,HB(:,k)];
            HB(:,k)=[];
            ma=ma+1;
            mb=mb-1;
        else
            expand=0;
            clear GML
            clear MvAvgGML
            p=0;
            compress=1;
            expand=0;
        end
    end

    if compress
        if Visibility
            figure(1)
            subplot(2,1,1);
            plot_index=ma;
            plot(plot_index,GML_plot,'o','color',[0.66,0.66,0.66]);
            axis([0,N,0,500]);
            title('Compressing...','Fontsize',16);
        end
        if (p<=FIR_len|MvAvgGML(p)>=MvAvgGML(p-1)|GML(p)>=MvAvgGML(p-1))&ma>2            
            %compress
            k=find(delta2_a==min(delta2_a),1,'first');
            B_index=[B_index,A_index(k)];
            A_index(k)=[];
            delta2_a(k)=[];
            HB=[HB,HA(:,k)];
            HA(:,k)=[];
            ma=ma-1;
            mb=mb+1;
        else
            compress=0;
            A_index=A_index_star;
            B_index=B_index_star;
            s_A=s_A_star;
            ma=ma_star;
            mb=mb_star;            
            sigma2=sigma2_star;
            clear GML
            clear MvAvgGML
            GML_plot=GML_star;
        end
    end%if compress

    %Reset step
    if ~(expand|compress)
        if GML_record<GML_star
            GML_record=GML_star;
            s_A_record=s_A;
            ma_record=ma;
            A_index_record=A_index;
            B_index_record=B_index;
            sigma2_record=sigma2;
            GML_star=-inf;
            
            s_init=zeros(m,1);
            s_init(A_index_record)=s_A_record;
            u=s_init+H'*Inv_HH*(y-H*s_init);
            u_sort=sort(abs(u),'descend');
            A_index=find(abs(u)>=u_sort(ma))';
            B_index=1:m;
            B_index(A_index)=[];
            HA=H(:,A_index);
            HB=H(:,B_index);
            
            s_A_init=u(A_index);
            sigma2=(y-HA*s_A_init)'*C_inv*(y-HA*s_A_init)/N;
            delta2_a=10*sigma2*ones(ma,1)./diag(HA'*HA);
            
            expand=1;
            compress=0;
            p=0;
        end
        if length(union(A_index_record_pre,A_index_record))==length(intersect(A_index_record_pre,A_index_record))
            exit_flag=1;
        end
        A_index_record_pre=A_index_record;
    end
end

A_index=A_index_record;
s_hat=zeros(m,1);
s_hat(A_index_record)=s_A_record;
sigma2_hat=sigma2_record;

s_hat=s_hat/scale;
sigma2_hat=sigma2_hat/(scale^2);

if Visibility
    GML_plot=GML_record;
    plot_index=ma_record;
    figure(1)
    subplot(2,1,1);
    plot(plot_index,GML_plot,'rp-','markersize',20);    
    title('Done','Fontsize',16);
    hold off
    subplot(2,1,2);
    stem(A_index_record,s_A_record/scale,'Marker','none','linewidth',2,'color','k');
    axis([1,m,-1.5,+1.5]);
    xlabel('signal index','Fontsize',12);
    ylabel('signal value','Fontsize',12);
    title('ExCoV estimated signal','Fontsize',16);
    figure(1)
    close(1);
end
