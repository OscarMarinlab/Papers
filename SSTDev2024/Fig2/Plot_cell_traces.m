clear
%load Fallall
[FileName3,PathName3,FilterIndex3]=uigetfile('*.mat');
addpath(PathName3);
load(FileName3,'Fnew');
load(FileName3,'dF_traces');

%load Selectionall Clusters
[FileName4,PathName4,FilterIndex4]=uigetfile('*.mat');
addpath(PathName4);
load(FileName4,'idd2');

other=1:15:150;

for i=1:length(other)
    plot(spike_rates(other(i),:)+i);
    hold on;
    
end

SST=idd2(1:5:size(idd2,2));
for i=1:length(SST)
    plot(spike_rates(SST(i),:)+length(other)+i);
    hold on;
end

other=1:15:150;

for i=1:length(other)
    plot(dF_traces(other(i),:)+i);
    hold on;
    
end

SST=idd2;
for i=1:length(SST)
    plot(dF_traces(SST(i),:)+length(other)+i);
    hold on;
end