% Simulate PC-MRI velocity estimation
% Compare runtimes for MLPUE and MLE from IQ data
%
% Copyright (c) 2023, David Tucker, Shen Zhao, Lee C. Potter
% All rights reserved.
% 
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree. 

%%
clearvars;close all

L=4;%four-point encoding
m = L*(L-1)/2; % Number of phase difference pairs
num_trials=1e2;% (=1e4) Number of random draws for avg timing
num_MLgrid=100;% (=100 pts for reporting, although Carrillo suggests 1000)
    % points per edge of polytope
T=8;%number of coils
s=10;%signal to noise ratio
% Comment: about 10 minutes on MacAir 2GHz Intel i5 for 10000 trials at 100
% points per edge of 3D polytope

% For array & A matrix: balanced four point encoding
gammarho = pi/100;
Moment = [-1,-1,-1;1,1,-1;1,-1,1;-1,1,1];
[Nencode,Ndim]=size(Moment);

%% construct operators
pair_idx = nchoosek(1:L,2);
Ac=Moment(pair_idx(:,2),:) - Moment(pair_idx(:,1),:);
Mc=2*pi/gammarho*eye(m);
P91=[1,0,0,0,0,-1; 1,0,0,0,0,1; 0,0,1,1,0,0];%Pelc 1991 preprocessing
Sigma = awgn_phase_covariance(pair_idx,s,T);
Prec = pinv(Sigma);
V = Lambda_basis(Ac,Mc);
Unimodular = [0,1,-1;0,0,-1;1,0,2];
V = V*Unimodular; %change the lattice basis to have 2 orthogonal columns
pinvV=pinv(V);
pre=[];

% sampling grid for IQ MLE
grid=linspace(-0.50,0.50,num_MLgrid);
[gridx, gridy,gridz] = ndgrid(grid,grid,grid);
vgrid = V*[gridx(:),gridy(:),gridz(:)]';

TrialsStart=tic;
%% Draw points across polytope and simulate noisy data
TimeML=zeros(num_trials,1);
TimeIQ=zeros(num_trials,1);
for ntrial = 1:num_trials
    x = V*rand(3,1);%uniform across polytope
    ynoiseless=s*exp(1i*gammarho*Moment*x);
    % unit variance iid complex gaussian noise
    Y = repmat(ynoiseless,1,T) + ...
        (1/sqrt(2))*(randn(L,T)+1i*randn(L,T));
    R = Y*Y';%data covariance
    y = angle(tril_vec(R));%lower triangular phases
    % MLPUE from phase differences
    tic
    [VhatML,~,pre] = milm_mle(Ac,y/gammarho,Mc,V,Prec,1,pre);
    %put onto centered interval
    coordinate = pinvV*VhatML;
    coordinate(coordinate >0.5) = coordinate(coordinate >0.5)-1;
    VhatML = V*coordinate;
    TimeML(ntrial)=toc;
    %IQ-MLE from IQ data
    a =exp(1j*gammarho*Moment*vgrid);%precompute
    tic;
    R2=sqrtm(R);
    [~,indx]=max( sum(abs(R2*a).^2, 1));%matlab friendly
    VhatIQ=vgrid(:,indx);
    TimeIQ(ntrial)=toc;
    if(ntrial==5)%view the estimators
        truth = pinvV*x;%place on centered interval
        truth(truth >0.5) = truth(truth >0.5)-1;
        truth = V*truth;
        fprintf('Example trial:\n')
        fprintf('     Truth     MLPUE     IQMLE\n')
        disp([truth VhatML VhatIQ]);
    end
end
%omit first call--don't include precompute time
ComputeRatio = mean(TimeIQ(2:end)) / mean(TimeML(2:end));
display(ComputeRatio)
TotalTime = toc(TrialsStart) 
