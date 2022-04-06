%% Clear

clear
clc
close all

%% Add dirs into path

currentFolder = pwd;
addpath(genpath(currentFolder));

%% Configurations

% for simulation
N_T     = 1e6;          % Number of realization for averaging (in order to obtain a smooth curve)
L       = 1e-4;         % Parameter of Possion distribution for each individual
M       = 100;          % Number of measurements
N       = 1000;         % Signal length
K       = ceil(L*N);    % Estimated sparsity level
S       = 2;            % number of indices selected in each iteration of MOLS
SNR     = (0:5:50);     % SNR in db
P       = 1/50;         % Probability for Bernoulli matrix A
X_mu    = log10(10^6);  % lnX~N(log(1e6),log(1e3)/3).
X_sigma = log10(10^3)/3;
Y_thres = 500;          % Thresholding for RT-PCR

%% algorithms

fields  = {'name','type','func','color','linestyle','marker','markerSize'};

methods = {
    'CGT, $c=1$', 1, @(y,Psi)MOLS_cK(y,Psi,1,K,eps,S), 'k', '-', 's', 5;...
    'CGT-Bin, $c=1$', 2, @(y_bin,Psi)MOLS_cK(y_bin,Psi,1,K,eps,S), 'k', '-', 'none', 5;...
    'CGT, $c=2$', 1, @(y,Psi)MOLS_cK(y,Psi,2,K,eps,S), 'b', '-', 'o', 6;...
    'CGT-Bin, $c=2$', 2, @(y_bin,Psi)MOLS_cK(y_bin,Psi,2,K,eps,S), 'b', '--', 'none', 5;...
    'CGT, $c=4$', 1, @(y,Psi)MOLS_cK(y,Psi,4,K,eps,S), 'r', '-', '+', 5;...
    'CGT-Bin, $c=4$', 2, @(y_bin,Psi)MOLS_cK(y_bin,Psi,4,K,eps,S), 'r', ':', 'none', 5;...
    'CGT, $c=8$', 1, @(y,Psi)MOLS_cK(y,Psi,8,K,eps,S), [0,0.6,0.6], '-', '^', 4;...
    'CGT-Bin, $c=8$', 2, @(y_bin,Psi)MOLS_cK(y_bin,Psi,8,K,eps,S), [0,0.6,0.6], '-.', 'none', 5;...
    };

methods     = cell2struct(methods,fields,2);
num_method  = numel(methods);

%% Test

% loop for SNR
num_snr     = length(SNR);
precision   = zeros(num_method,num_snr);
recall      = zeros(num_method,num_snr);

for idx_snr = 1:num_snr
    
    tic1    = tic;
    
    tp          = zeros(num_method,1);
    fp          = zeros(num_method,1);
    fn          = zeros(num_method,1);
    
    % compute the std of Gaussian noise (without dilute)
    snr         = 10 ^ (SNR(idx_snr) / 10);
    snr_sqrt    = sqrt(snr);

    % loop for averaging
    for idx_itr = 1:N_T
        
        [x,~,T] = GenSparseVec_COVID19(L,X_mu,X_sigma,N);
        if isempty(T)
           continue
        end
        A       = generate_A(M,N,P);
        dilute  = sum(A > eps,2);

        % synthesis data for CGT
        z       = A * x;
        z_pos   = z > eps;
        W_sigma = min(z(z_pos)) / (sqrt(sum(z_pos)) *snr_sqrt);
        z(z_pos)= z(z_pos) + randn(sum(z_pos),1) .* W_sigma;
        z_bin   = double(z > Y_thres);

        % Pretreatments
        % we omit scaling here, since it will not affect the algorithm
        D       = diag(dilute);
        Phi     = D * A;
        u       = D * z;
        u_bin   = D * z_bin;

        % subtract mean
        Psi     = Phi - mean(Phi);
        y       = u - mean(u);
        y_bin   = u_bin - mean(u_bin);
      
        % loop for methods
        for idx_method = 1:num_method
            tmp_method = methods(idx_method);
            
            switch tmp_method.type
                case 1
                    [~,T_hat,~] = tmp_method.func(y,Psi);
                case 2
                    [~,T_hat,~] = tmp_method.func(y_bin,Psi);
            end
            
            tp(idx_method)  = tp(idx_method) + length( intersect(T, T_hat) );
            fp(idx_method)  = fp(idx_method) + length( setdiff(T_hat, T) );
            fn(idx_method)	= fn(idx_method) + length( setdiff(T, T_hat) );
        end

    end
    
    precision(:,idx_snr)  = tp ./ (tp + fp + eps);
    recall(:,idx_snr)     = tp ./ (tp + fn + eps);
    
    disp(['***' num2str(idx_snr) '-th process took ' num2str(toc(tic1)) ' (sec.).']);

end

%% Plot

figure(1);
save('fig_SNR.mat','precision','recall')

subplot(1,2,1)
hold on
for idx_method = 1:num_method
    tmp_method = methods(idx_method);
    smooth = spcrv([[SNR(1),SNR,SNR(end)];...
    [precision(idx_method,1),precision(idx_method,:),precision(idx_method,end)]],4,50);
    plot(smooth(1,:),smooth(2,:),...
        'Color',tmp_method.color,'LineStyle',tmp_method.linestyle,'Marker',tmp_method.marker,...
        'linewidth',1,'MarkerSize',tmp_method.markerSize,'MarkerIndices',1:5:size(smooth,2))
end
xlabel('SNR (dB)','Interpreter','latex')
ylabel('Precision','Interpreter','latex')
ylim([0,0.6])
% legend(methods.name,'Location','southeast','Interpreter','latex')
set(gca,'Fontname','times new Roman');
axis square

subplot(1,2,2)
hold on
for idx_method = 1:num_method
    tmp_method = methods(idx_method);
    smooth = spcrv([[SNR(1),SNR,SNR(end)];...
    [recall(idx_method,1),recall(idx_method,:),recall(idx_method,end)]],4,50);
    plot(smooth(1,:),smooth(2,:),...
        'Color',tmp_method.color,'LineStyle',tmp_method.linestyle,'Marker',tmp_method.marker,...
        'linewidth',1,'MarkerSize',tmp_method.markerSize,'MarkerIndices',1:5:size(smooth,2))
end
xlabel('SNR (dB)','Interpreter','latex')
ylabel('Recall','Interpreter','latex')
ylim([0.8,1])
legend(methods.name,'Location','southeast','Interpreter','latex')
set(gca,'Fontname','times new Roman');
axis square
