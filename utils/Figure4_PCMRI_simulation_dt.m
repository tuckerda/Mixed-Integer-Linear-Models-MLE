% Simulate PC-MRI velocity estimation to compare to the Pelc 1991 estimator
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
num_trials=10;% (=1e3 for figure)
num_grid=50;% (=50 for figure) points-1 per edge of polytope
T=8;%number of coils
s=10;%signal to noise ratio
% Comment: about 60 minutes on MacAir 2GHz Intel i5 for 1000 trials at 51
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
[Pc_lb,Pc_ub]=Pc_bounds(Ac,Mc,Prec);
Unimodular = [0,1,-1;0,0,-1;1,0,2];
V = V*Unimodular; %change the lattice basis to have 2 orthogonal columns
pre=[];
fprintf('Approximate lower and upper bounds on wrap error:\n %g, %g\n',[1-Pc_lb,1-Pc_ub]);
fprintf('\nLattice basis V for velocity\n');disp(V)

%% Draw points across polytope and simulate noisy data
grid=linspace(-0.50,0.50,num_grid+1);
[gridx, gridy,gridz] = ndgrid(grid,grid,grid);
v = V*[gridx(:),gridy(:),gridz(:)]';
RMSEPelc=zeros(1,size(v,2));
RMSEML=zeros(1,size(v,2));
pinvV=pinv(V);
for n=1:size(v,2)%for each candidate true velocity
    ynoiseless=s*exp(1i*gammarho*Moment*v(:,n));
    sqerrPelc=zeros(num_trials,1);%clear
    sqerrMLPUE=zeros(num_trials,1);%clear
    for ntrial = 1:num_trials %for each noise trial
        Y = repmat(ynoiseless,1,T) + ...
            1/sqrt(2)*(randn(L,T)+1i*randn(L,T));
        %unit variance iid complex gaussian noise
        R = Y*Y';%data covariance 
        y = angle(tril_vec(R));%lower triangular phases
        VhatPelc = P91*y/(4*gammarho);%note scaling factor not in A, P91
        % MLE from phase differences
        [VhatML,~,pre] = milm_mle(Ac,y/gammarho,Mc,V,Prec,1,pre);  
        %put onto centered interval
        coordinate = pinvV*VhatML;
        coordinate(coordinate >0.5) = coordinate(coordinate >0.5)-1;
        VhatML = V*coordinate;
        % squared error
        sqerrPelc(ntrial) = sum( (v(:,n)-VhatPelc).^2);
        sqerrMLPUE(ntrial) = sum( (v(:,n)-VhatML).^2);
    end
    %RMSE
    RMSEPelc(n) = sqrt(mean(sqerrPelc));
    RMSEML(n) = sqrt(mean(sqerrMLPUE));
end

%% Construct Figures
fsize=20;viewAz=10;viewEl=30;plotsize=30;maxc=100;%110,10
opaq=0.4;

C = int2bit(0:2^3-1,3);
Vertices = V*C;
Faces = [1 5 7 3;
    2 6 8 4;
    1 2 6 5;
    3 4 8 7;
    1 3 4 2;
    5 7 8 6];
Colors = [0.2157    0.4941    0.7216];%linspecer(1);

figure;
hold on;
scatter3(v(1,:), v(2,:), v(3,:),plotsize,RMSEPelc);
alpha(opaq)
axis image;
view(viewAz,viewEl);%view(10,25);%view(10,70)
colorbar
caxis([0,maxc])
patch('Vertices',Vertices'-sum(V,2)'/2,'Faces',Faces,...
    'FaceColor',Colors(1,:),'FaceAlpha',0,'EdgeAlpha',1);
set(gca,'FontSize',fsize,'FontName','Times New Roman')
%exportgraphics(gca,'Figure4_Pelc.eps');

figure;
hold on;
scatter3(v(1,:), v(2,:), v(3,:),plotsize,RMSEML);
alpha(opaq)
axis image
view(viewAz,viewEl);%view(0,25)
colorbar
caxis([0,maxc])
patch('Vertices',Vertices'-sum(V,2)'/2,'Faces',Faces,...
    'FaceColor',Colors(1,:),'FaceAlpha',0,'EdgeAlpha',1);
set(gca,'FontSize',fsize,'FontName','Times New Roman')
%exportgraphics(gca,'Figure4_MLPUE.eps');


%% error variances & wrapping
SigmaML=  pinv(gammarho*Ac'*Prec*Ac*gammarho);
SigmaPelc = P91*Sigma*P91'/(4*gammarho)^2;
ErrorEllipseVolRatio=sqrt(det(SigmaPelc))/sqrt(det(SigmaML));
display(ErrorEllipseVolRatio)
PctUnaliasedNotWrapped = [length(find(RMSEPelc < 3))/length(RMSEPelc) ...
length(find(RMSEML< 3))/length(RMSEML)]*100;
display(PctUnaliasedNotWrapped)
