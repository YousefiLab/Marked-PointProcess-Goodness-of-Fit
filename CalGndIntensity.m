% Behzad Nazari
% 2019, Nov 1

% this function is modified by user to generate the desired Lambda
% if GroundIntensity is available in any other way , rewrite this code to
% provide ground intensity without extra calculations
% If ground intensity is not available , this function fetches Lambda in
% the given time and mark points and integrate to calculate ground
% intensity

% input arguement: 
% Time: a vector represnting the real time samples that lambda should be samples
% M1, M2, M3: a one dimensional array, representing where the Mark space
% should be in each dimension
% Memory: memoty limit to be used in MegaByte
% warning: high memory usage may freeze the computer

% output arguement:
% a vector representing ground intensity in the given Time samples

function [GndIntensity]=CalGndIntensity(Time,M)
sprintf('calculating ground intensity')
GndIntensity=zeros(1,length(Time));
M1 = M{1,1};
M2 = M{1,2};
M3 = M{1,3};
dms=(M1(2)-M1(1))*(M2(2)-M2(1))*(M3(2)-M3(1));

BlockNo= ceil(8*length(Time)*length(M1)*length(M2)*length(M3)/(500*1e6));
BlockWidth=round(length(Time)/BlockNo);

for k=0:BlockNo-1
    sprintf(' processing block %d out of %d blocks',k,BlockNo)
    if(k<BlockNo-1)
        tt=Time(k*BlockWidth+1:k*BlockWidth+BlockWidth);
        [Lambda1]=Lambda(tt,M1,M2,M3);    
    else 
        tt=Time(k*BlockWidth+1:length(Time));
        [Lambda1]=Lambda(tt,M1,M2,M3);    
    end
    for kk=1:length(tt)
        GndIntensity(k*BlockWidth+kk)=dms*sum(sum(sum(squeeze(Lambda1(kk,:,:,:)))));
    end
end

