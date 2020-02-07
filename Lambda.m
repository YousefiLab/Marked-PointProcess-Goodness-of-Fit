% Behzad Nazari
% 2019, Nov 1

% this function is modified by user to generate the desired Lambda

% input arguement: 
% Time: a vector represnting the real time samples that lambda should be samples
% M1, M2, M3: a one dimensional array, representing where the Mark space
% should be in each dimension

% output arguement:
% Lambda1: a 3 dimensional array representing the Lambda for the provided
% time and mark space

function [Lambda]=Lambda(Time,M1,M2,M3)
a1   = log(0.3);% Neuron1 fixed rate parameter  %0.15
miu  = -2;% The gaussian mean of lambda_neuron1
s_x  = sqrt(.5);% The gaussian variance of lambda
% xt   = 0.8*sin(8*Time*2*pi/8);
K = length(Time);
rng(1)
xt = zeros(1,K);
for k=2:K
    xt(k) = 0.98*xt(k-1)+0.3*randn;
end
GroundIntensity = exp(a1-(xt-miu).^2/2/s_x/s_x );  
DataSize = length(Time)*length(M1)*length(M2)*length(M3)*8/1e6;
if(DataSize>1000)
    s45 = sprintf('%6d mega bytes needed, like to proceed?\ny: to continue\nany other key: terminate m file\n',DataSize)
    key1 = input(s45,'s');
    if(key1 ~= 'y')
        error('too much memory needed') ;
    end
end

Lambda = zeros(length(Time),length(M1),length(M2),length(M3));
for t=1:length(Time)
%     tt1 = rem(Time(t),4);
%     tt2 = 2-abs(tt1-2);
    mu  = [8   16   28];
    sigma = [0.3 .001 0.002 ; 0.001 0.04 0.003 ; 0.002 0.003 0.4];

    [X1,X2,X3] = ndgrid(M1,M2,M3);
    X  = [X1(:) X2(:) X3(:)];
    M  = mvnpdf(X,mu,sigma);
    MM = reshape(M,length(M1),length(M2),length(M3));
    Lambda(t,:,:,:) = MM*GroundIntensity(t);
end
