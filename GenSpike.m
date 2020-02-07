
function [Spikes]=GenSpike(Time,GroundIntensity)
s  = (GroundIntensity>rand(1,length(Time)));

event_mark       = zeros(length(find(s)),3);
event_mark(:,1)  = Time(s==1);% Spike times, go to column 1

[m,~]=size(event_mark);

for k=1:m
    mu=[8   16   28];
    sigma=[0.3 0.001 0.002 ; 0.001 0.04 .003 ; 0.002 0.003 0.4];
    rand1=mvnrnd(mu,sigma);
    event_mark(k,2:4)=rand1;
end
Spikes=event_mark;