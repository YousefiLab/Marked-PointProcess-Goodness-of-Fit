clear
close all
clc

%% choose the time and mark resolutions according to the Lambda function and desired intervals
dt=0.001;
t2=0:dt:8;  % this is the simulation interval 
M1_l = 5:0.06:11;
M2_l = 15:0.03:18;
M3_l = 24:0.1:34;
M_l = {M1_l M2_l M3_l};
%% compute ground intensity
GndIntensity = CalGndIntensity(t2,M_l);   % compute ground intensity

%% generate spikes for simulation, does not needed in real cases
Spikes = GenSpike(t2,GndIntensity);        % get the spike info from user, this is provided by user, too

%% IRCM
M1_h = 5:0.006:11;  % about 250 points seems good
M2_h = 15:0.003:18;  % more than 500 points ok
M3_h = 24:0.01:34;  % more than 500 points ok
M_h = {M1_h M2_h M3_h};
[u_1,v_1] = IRCM(GndIntensity,Spikes,M_h,M_l,dt);

%% Uniformity test
% Multivariate KS test:
X = [u_1',v_1];
C_a = 1.12; % it is be calculated by Monte Carlo method (Table 2 in paper)
MKS_1 = Gof_MKS(X,C_a);
if (MKS_1 == 1)
    disp('MKS:The null hypothesis is passed under MDCI')
else
    disp('MKS:The null hypothesis is rejected under MDCI')
end
DF = 26; % Degree of freedom
P = 0.05; % p-value 
chi_squ_1 = Gof_Pearson(X,DF,P);
if (chi_squ_1 == 1)
    disp('ChiS:The null hypothesis is passed under IRCM')
else
    disp('ChiS:The null hypothesis is rejected under IRCM')
end

%% Visualising
figure;%% 2-D rescaled time samples via v_i^(1)
plot(u_1,v_1(:,1),'.');
xlabel('u_i')
ylabel('v_i^{(1)}')
title('Rescaled time and mark samples')

figure;% 3-hypercube rescaled mark samples v_i
plot3(v_1(:,1),v_1(:,2),v_1(:,3),'.','Markersize',12)
grid on
xlabel('v_i^{(1)}');ylabel('v_i^{(2)}');zlabel('v_i^{(3)}')
title('Rescaled mark samples')

%% compute Mark intensity
MarkIntensity = CalMarkIntensity(t2,M_l,dt);

%% MDCI
[u_2,v_2] = MDCI(t2,MarkIntensity,Spikes,M_h,M_l,dt);

%% Uniformity test
% Multivariate KS test:
X = [u_2',v_2];
C_a = 1.12; % it should be calculated by Monte Carlo method
MKS_2 = Gof_MKS(X,C_a);
if (MKS_2 == 1)
    disp('MKS:The null hypothesis is passed under MDCI')
else
    disp('MKS:The null hypothesis is rejected under MDCI')
end
DF = 26;
P = 0.05;
chi_squ_2 = Gof_Pearson(X,DF,P);
if (chi_squ_2 == 1)
    disp('ChiS:The null hypothesis is passed under MDCI')
else
    disp('ChiS:The null hypothesis is rejected under MDCI')
end

%% Visualising
figure;%% 2-D rescaled time samples via v_i^(1)
plot(u_2,v_2(:,1),'.');
xlabel('u_i')
ylabel('v_i^{(1)}')
title('Rescaled time and mark samples')

figure;% 3-hypercube rescaled mark samples v_i
plot3(v_2(:,1),v_2(:,2),v_2(:,3),'.','Markersize',12)
grid on
xlabel('v_i^{(1)}');ylabel('v_i^{(2)}');zlabel('v_i^{(3)}')
title('Rescaled mark samples')
