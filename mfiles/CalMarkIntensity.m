% Yalda Amidi
% 2020, Jan 12

% input arguement: 
% Time: a vector represnting the real time samples that lambda should be samples
% M1, M2, M3: a one dimensional array, representing where the Mark space
% should be in each dimension
% dt: time resolution
% Memory: memoty limit to be used in MegaByte
% warning: high memory usage may freeze the computer

% output arguement:
% a array representing mark intensity in the given Time samples

function [MarkIntensity]=CalMarkIntensity(Time,M,dt)
sprintf('calculating mark intensity')
M1 = M{1,1};
M2 = M{1,2};
M3 = M{1,3};
MarkIntensity=zeros(length(M1),length(M2),length(M3));


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
    
    MarkIntensity = MarkIntensity+dt*squeeze(sum(Lambda1,1));
    
end

