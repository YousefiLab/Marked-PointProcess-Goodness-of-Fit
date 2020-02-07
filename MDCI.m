function [u,v]=MDCI(Time,MarkIntensity,Spikes,M_h,M_l,dt)

M1_l = M_l{1,1};
M2_l = M_l{1,2};
M3_l = M_l{1,3};

M1_h = M_h{1,1};
M2_h = M_h{1,2};
M3_h = M_h{1,3};
[m,~]=size(Spikes);
u=zeros(1,m);
for i=1:m
    m1=(Spikes(i,2)-M1_l(1))/(M1_l(2)-M1_l(1))+1;
    m2=(Spikes(i,3)-M2_l(1))/(M2_l(2)-M2_l(1))+1;
    m3=(Spikes(i,4)-M3_l(1))/(M3_l(2)-M3_l(1))+1;
    si=round(Spikes(i,1)/dt)+1;
    Landa=Lambda(Time(Time<=Spikes(i,1)),Spikes(i,2),Spikes(i,3),Spikes(i,4));
    m1=round(m1);
    m2=round(m2);
    m3=round(m3);
    Gamma = MarkIntensity(m1,m2,m3);
    ft_m = Landa/Gamma;
    u(i)= min(1,sum(ft_m(1:si)*dt));
end
v1=zeros(1,m);
SMI = sum(MarkIntensity,'all');% calculate the summation of mark intensity over all dimensions
fm=(sum(sum(MarkIntensity,2),3))/SMI;
for i=1:m
    sprintf('MDCI: v1 %d/%d',i,m)
    m2=(Spikes(i,2)-M1_h(1))/(M1_h(2)-M1_h(1))+1;
    m2=round(m2);
    v1(i)=sum(fm(1:m2));
end
fm  = squeeze(sum(MarkIntensity,3))/SMI;
fm_1 =(sum(sum(MarkIntensity,2),3))/SMI;
v2=zeros(1,m);
for i=1:m
    sprintf('MDCI: v2 %d/%d',i,m)
    m1=(Spikes(i,2)-M1_h(1))/(M1_h(2)-M1_h(1))+1;
    m1=round(m1);
    m2=(Spikes(i,3)-M2_h(1))/(M2_h(2)-M2_h(1))+1;
    m2=round(m2);
    v2(i)=sum(fm(m1,1:m2),2)/(fm_1(m1));
end
fm_12  = squeeze(sum(MarkIntensity,3))/SMI;
fm = MarkIntensity/SMI;
v3=zeros(1,m);
for i=1:m
    sprintf('MDCI: v3 %d/%d',i,m)
    m1=(Spikes(i,2)-M1_h(1))/(M1_h(2)-M1_h(1))+1;
    m1=round(m1);
    m2=(Spikes(i,3)-M2_h(1))/(M2_h(2)-M2_h(1))+1;
    m2=round(m2);
    m3=(Spikes(i,4)-M3_h(1))/(M3_h(2)-M3_h(1))+1;
    m3=round(m3);
    v3(i)=sum(fm(m1,m2,1:m3),3)/fm_12(m1,m2);
end
v = [v1' v2' v3'];