function [u,v]=IRCM(GndIntensity,Spikes,M_h,M_l,dt)

[m,~]=size(Spikes);
u=zeros(1,m);
s2=round(Spikes(1,1)/dt+1);
s=sum(GndIntensity(1:s2-1));
u(1)=1-exp(-s);
for k=2:m
    s1=round(Spikes(k-1,1)/dt+1);
    s2=round(Spikes(k,1)/dt+1);
    s=sum(GndIntensity(s1:s2-1));
    u(k)=1-exp(-s);
end
% high resolution mark intervals
M1_h = M_h{1,1};
M2_h = M_h{1,2};
M3_h = M_h{1,3};
% low resolution mark intervals
% M1_l = M_l{1,1};
M2_l = M_l{1,2};
M3_l = M_l{1,3};
v1=zeros(1,m);
for i=1:m
    sprintf('IRCM: v1 %d/%d',i,m)
    Landa1=Lambda(Spikes(i,1),M1_h,M2_l,M3_l);
    fm=zeros(1,length(M1_h));
    for m1=1:length(M1_h)
        tt=squeeze(Landa1(1,m1,:,:));
        fm(m1)=sum(sum(tt));
    end
    fm=fm/sum(fm);
    m2=(Spikes(i,2)-M1_h(1))/(M1_h(2)-M1_h(1))+1;
    m2=round(m2);
    v1(i)=sum(fm(1:m2));
end
v2=zeros(1,m);
for i=1:m
    sprintf('IRCM: v2 %d/%d',i,m)
    Landa1=Lambda(Spikes(i,1),Spikes(i,2),M2_h,M3_l);
    fm=zeros(1,length(M2_h));
    for m2=1:length(M2_h)
        tt=squeeze(Landa1(1,1,m2,:));
        fm(m2)=sum(tt);
    end
    fm=fm/sum(fm);
    m2=(Spikes(i,3)-M2_h(1))/(M2_h(2)-M2_h(1))+1;
    m2=round(m2);
    v2(i)=sum(fm(1:m2));
end
v3=zeros(1,m);
for i=1:m
    sprintf('IRCM: v3 %d/%d',i,m)
    Landa1=Lambda(Spikes(i,1),Spikes(i,2),Spikes(i,3),M3_h);
    fm=squeeze(Landa1(1,1,1,:));
    fm=fm/sum(fm);
    m2=(Spikes(i,4)-M3_h(1))/(M3_h(2)-M3_h(1))+1;
    m2=round(m2);
    v3(i)=sum(fm(1:m2));
end
v = [v1' v2' v3'];