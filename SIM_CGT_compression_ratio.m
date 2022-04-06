%% Clear

clear
clc
close all

%% Add dirs into path

currentFolder = pwd;
addpath(genpath(currentFolder));

%% Configurations

% for simulation
N_T     = 1e6;          % Number of realization for averaging
L       = 1e-4;         % Parameter of Possion distribution for each individual
M       = 100;          % Number of measurements
N       = 500:500:4500; % Signal length
K       = ceil(L*N);    % Estimated sparsity level
S       = 2;            % number of indices selected in each iteration of MOLS
P       = 1/50;         % Probability for Bernoulli matrix A
X_mu    = log10(10^6);  % lnX~N(log(1e6),log(1e3)/3).
X_sigma = log10(10^3)/3;

%% algorithms

fields  = {'name','type','func','color','linestyle','marker'};

methods = {
    'CGT, $c=1$', 1, @(y,psi,k)MOLS_cK(y,psi,1,k,eps,S), 'k', '-', '+';...
    'CGT-Bin, $c=1$', 2, @(y_bin,psi,k)MOLS_cK(y_bin,psi,1,k,eps,S), 'k', '-', 'o';...
    'CGT, $c=2$', 1, @(y,psi,k)MOLS_cK(y,psi,2,k,eps,S), 'b', '-', 'x';...
    'CGT-Bin, $c=2$', 2, @(y_bin,psi,k)MOLS_cK(y_bin,psi,2,k,eps,S), 'b', '--', 'd';...
    'CGT, $c=4$', 1, @(y,psi,k)MOLS_cK(y,psi,4,k,eps,S), 'r', '-', 's';...
    'CGT-Bin, $c=4$', 2, @(y_bin,psi,k)MOLS_cK(y_bin,psi,4,k,eps,S), 'r', ':', '^';...
    'CGT, $c=8$)', 1, @(y,Psi,k)MOLS_cK(y,Psi,8,k,eps,S), [0,0.6,0.6], '-', 'h';...
    'CGT-Bin, $c=8$', 2, @(y_bin,Psi,k)MOLS_cK(y_bin,Psi,8,k,eps,S), [0,0.6,0.6], '-.', '>';...
    };

methods     = cell2struct(methods,fields,2);
num_method  = numel(methods);

%% Test

% loop for SNR
num_N     = length(N);
precision   = zeros(num_method,num_N);
recall      = zeros(num_method,num_N);

for idx_N = 1:num_N
    
    tic1    = tic;
    
    tp      = zeros(num_method,1);
    fp      = zeros(num_method,1);
    fn      = zeros(num_method,1);
    
    N_tmp   = N(idx_N);
    K_tmp   = K(idx_N);
    
    % loop for averaging
    for idx_itr = 1:N_T
        
        [x,~,T] = GenSparseVec_COVID19(L,X_mu,X_sigma,N_tmp);
        if isempty(T)
           continue
        end
        A       = generate_A(M,N_tmp,P);
        dilute  = sum(A > eps,2);

        % synthesis data for CGT
        z       = A * x;
        z_bin   = double(z > eps);

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
                    [~,T_hat,~] = tmp_method.func(y,Psi,K_tmp);
                case 2
                    [~,T_hat,~] = tmp_method.func(y_bin,Psi,K_tmp);
            end
            
            tp(idx_method)  = tp(idx_method) + length( intersect(T, T_hat) );
            fp(idx_method)  = fp(idx_method) + length( setdiff(T_hat, T) );
            fn(idx_method)	= fn(idx_method) + length( setdiff(T, T_hat) );
        end

    end
    
    precision(:,idx_N)  = tp ./ (tp + fp + eps);
    recall(:,idx_N)     = tp ./ (tp + fn + eps);
    
    disp(['***' num2str(idx_N) '-th process took ' num2str(toc(tic1)) ' (sec.).']);

end

%% Plot

figure(1);
save('fig_compression.mat','precision','recall')
subplot(1,2,1)
hold on
for idx_method = 1:num_method
    tmp_method = methods(idx_method);
    smooth = spcrv([[N(1),N,N(end)]./M; ...
    [precision(idx_method,1),precision(idx_method,:),precision(idx_method,end)]],4,50);
    plot(smooth(1,:),smooth(2,:),...
        'Color',tmp_method.color,'LineStyle',tmp_method.linestyle,'Marker',tmp_method.marker,...
        'linewidth',1,'MarkerSize',5,'MarkerIndices',1:5:size(smooth,2))
end
xlabel('Compression ratio $n/m$','Interpreter','latex')
ylabel('Precision','Interpreter','latex')
axis([5,45,0,0.6])
% legend(methods.name,'Location','southeast','Interpreter','latex')
set(gca,'Fontname','times new Roman');
axis square

subplot(1,2,2)
hold on
for idx_method = 1:num_method
    tmp_method = methods(idx_method);
    smooth = spcrv([[N(1),N,N(end)]./M; ...
    [recall(idx_method,1),recall(idx_method,:),recall(idx_method,end)]],4,50);
    plot(smooth(1,:),smooth(2,:),...
        'Color',tmp_method.color,'LineStyle',tmp_method.linestyle,'Marker',tmp_method.marker,...
        'linewidth',1,'MarkerSize',5,'MarkerIndices',1:5:size(smooth,2))
end
xlabel('Compression ratio $n/m$','Interpreter','latex')
ylabel('Recall','Interpreter','latex')
axis([5,45,0.7,1])
legend(methods.name,'Location','southeast','Interpreter','latex')
set(gca,'Fontname','times new Roman');
axis square
