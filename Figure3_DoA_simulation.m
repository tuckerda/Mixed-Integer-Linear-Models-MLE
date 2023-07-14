% Simulate planar array output to compare DoA estimation performance of  maximum likelihood
% phase unwrapping estimation to grid search-based maximum likelihood estimation with complex data
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
% 
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree. 
%
%%
clearvars; close all;
n = 2;
L = 5; % Number of sensors
m = L*(L-1)/2; % Number of sensor pairs
T = 10; % Number of snapshots
num_trials = 1e4;
SNR_db = (5:1:14).'; % Signal-to-noise ratio (dB)
theta = (1/4)*2*pi*ones([1,num_trials]); % Azimuth angle of arrival (radians)
phi = (1/2)*pi/2*ones([1,num_trials]); % Elevation angle of arrival (radians)
omega = ang2dir([theta; phi],n);

% Array element positions in units of one-half wavelength
positions = [0,0; -4,0; 0,4; 0,-5; 5,0]; % Holdsworth 2005 cross-shaped array
pair_idx = nchoosek(1:L,2);
Ac = pi*(positions(pair_idx(:,2),:) - positions(pair_idx(:,1),:));
Mc = 2*pi*eye(m);
V = form_Lambda_basis(Ac,Mc);

%% Simulate array output 
sp = 1; % Signal power
SNR_lin = 10.^(SNR_db/10);

% Form grids for grid search-based MLE
grid_res = [1e-2,1e-3];
omega_scan = cell([length(grid_res),1]);
for k=1:length(grid_res)
    if k==1
        omega_bounds = [-ones([n,1]), ones([n,1])];
    else
        omega_bounds = [-grid_res(k-1)*ones([n,1]), grid_res(k-1)*ones([n,1])];
    end
    omega_grid = cell(n,1);
    for i = 1:n
        omega_grid(i) = { linspace( omega_bounds(i,1), omega_bounds(i,2), ceil( ( omega_bounds(i,2) - omega_bounds(i,1))/grid_res(k) + 1 ) ) };
    end
    tmp = cell(size(omega_grid));
    [tmp{:}] = ndgrid(omega_grid{:});
    tmp = reshape(cat(n+1,tmp{:}),[],n).';
    tmp(:, vecnorm(tmp,2,1) > 1 ) = [];
    omega_scan{k} = tmp;
end

omega_iqml = zeros([n,numel(SNR_db),num_trials]);
omega_pml = zeros([n,numel(SNR_db),num_trials]);
CRB = zeros([numel(SNR_db),1]);
pml_rmse = zeros([numel(SNR_db),1]);
iqml_rmse = zeros([numel(SNR_db),1]);
rmse_fun = @(omega_hat) sqrt( 1/num_trials * sum( (omega_hat(:) - omega(:)).^2 ) );
for ii=1:length(SNR_db)
    n_var = sp / SNR_lin(ii); % Noise power
    disp("SNR value "+ii+" of "+length(SNR_db)+": SNR = "+SNR_db(ii)+" dB")
    Ry = zeros([L,L,num_trials]);
    phi_hat = zeros([m,num_trials]);

    Sigma = awgn_phase_covariance(pair_idx,SNR_lin(ii),T);
    Prec = inv(Sigma);
    for k=1:num_trials
        s = sqrt(sp).*exp(1i*2*pi*rand([1,T])); %Signal
        y = exp(1i*pi*positions*omega(:,k))*s + sqrt(n_var/2)*(randn([L,T])+1i*randn([L,T])); % Array output
        Ry(:,:,k) = 1/T*(y*y'); % Sample covariance
        phi_hat(:,k) = tril_vec(angle(Ry(:,:,k)));
    end

    % DoA estimation through maximum likelihood phase unwrapping
    precompute = milm_mle_precompute(Ac,Mc,V,Prec);
    precompute = doa_pue_precompute(precompute);
    for k=1:num_trials
        omega_pml(:,ii,k) = doa_pue(Ac,phi_hat(:,k),Mc,V,Prec,precompute);
    end

    % Grid search-based maximum likelihood DoA estimation using complex data
    parfor k=1:num_trials
        omega_iqml(:,ii,k) =  doa_mle(Ry(:,:,k),positions,omega_scan);
    end

    iqml_rmse(ii) = rmse_fun( squeeze(omega_iqml(:,ii,:)) );
    pml_rmse(ii) = rmse_fun( squeeze(omega_pml(:,ii,:)) );
    [~,tmp] = crb_c(positions,T,[.1;.1],s,n_var);
    CRB(ii) = sqrt(trace(tmp));
end

%% Display results
figure(1); clf;
ax = axes;
set(ax,'YScale','log');
xlabel('SNR (dB)')
ylabel('RMSE')
axis(ax,[min(SNR_db),max(SNR_db),2.4e-3,0.25])
hold on;grid on;

line_crb = plot(SNR_db, CRB,'LineStyle','-','Marker','none','MarkerSize',6);
line_iqml = plot(SNR_db, iqml_rmse,'LineStyle','-.','Marker','.','MarkerSize',13);
line_pml = plot(SNR_db, pml_rmse,'LineStyle','none','Marker','o','MarkerSize',5);
legend_text = {'MLE','MLPUE','$\sqrt{ \rm{tr} \left( \rm{CRB} \right) }$'};
legend(ax,[line_iqml,line_pml,line_crb],legend_text, 'Location', 'Northeast','Interpreter','latex');
