clear all

%load red cells (cluster.mat)
[FileName5,PathName5,FilterInde5]=uigetfile('*.mat');
addpath(PathName5);
load(FileName5,'idd2');

%Load spikes and df/f (Fall.mat)
[FileName6,PathName6,FilterIndex6]=uigetfile('*.mat');
addpath(PathName6);
load(FileName6,'Fnew');
load(FileName6,'spike_rates');
load(FileName6,'locs');
load(FileName6,'ab');
load(FileName6,'ab2');
load(FileName6,'thresh');
load(FileName6,'iscell');


a1=iscell(:,1)==1;
a1=find(iscell(:,1)==1);
Fnew=Fnew(a1,:);

%Load H-events(Fallallparticipation.mat)
[FileName7,PathName7,FilterIndex7]=uigetfile('*.mat');
addpath(PathName7);
load(FileName7,'long2');
load(FileName7,'seq');

%% Mean Cell spikes (probability of spiking per frame)RED CELLS
RedCellSpike=nansum(spike_rates(idd2,:)')/(size(spike_rates,2)-64);
MeanRedCellSpike=mean(RedCellSpike)

%% Mean Cell df/f (average of calcium trasce) RED CELLS
RedCelldf=nansum(Fnew(idd2,:)')/(size(Fnew,2));
MeanRedCelldf=mean(RedCelldf)

%% Mean Cell spikes (probability of spiking per frame) Other
Other=1:size(spike_rates,1);
Other2=setdiff(Other, idd2);
OtherCellSpike=nansum(spike_rates(Other2,:)')/(size(spike_rates,2)-64);
MeanOtherCellSpike=mean(OtherCellSpike)

%% Mean Cell df/f (average calcium trace) for only green cells (Other)
Other=1:size(spike_rates,1);
Other2=setdiff(Other, idd2);
OtherCelldf=nansum(Fnew(Other2,:)')/(size(Fnew,2));
MeanOtherCelldf=mean(OtherCelldf)

%% Spike_probability in events RED cells (probability of spiking per frame)

mean_tempRed=[];
for i=1:length(locs);
    temp=spike_rates(idd2,locs(i)-40:locs(i)+40);
    mean_tempRed(i)=mean(temp(:));
end

grandMeanRed=mean(mean_tempRed);

%% Spike_probability in events All cells (probability of spiking per frame)
Other=1:size(spike_rates,1);
Other2=setdiff(Other, idd2);

mean_temp=[];
for i=1:length(locs);
    temp=spike_rates(Other2,locs(i)-40:locs(i)+40);
    mean_temp(i)=mean(temp(:));
end

grandMean=mean(mean_temp);

%% Median Cell spikes (probability of spiking per frame)RED CELLS
MedianRedCellSpike=median(spike_rates(idd2,33:end-32)');
GrandMedianRedCellSpike=mean(MedianRedCellSpike)

%% Median Cell df/f (median of calcium trasces) RED CELLS
RedCelldf=median(Fnew(idd2,:)');
MedianRedCelldf=mean(RedCelldf)

%% Median Cell spikes (probability of spiking per frame)Other
MedianOtherCellSpike=median(spike_rates(Other2,33:end-32)');
GrandMedianOtherCellSpike=mean(MedianOtherCellSpike)

%% Median Cell df/f (median of calcium traces) Other
OtherCelldf=median(Fnew(Other2,:)');
MedianOtherCelldf=mean(OtherCelldf)

%% SD Cell spikes (probability of spiking per frame)RED CELLS
sdRedCellSpike=std(spike_rates(idd2,33:end-32)');
GrandsdRedCellSpike=mean(sdRedCellSpike)

%% SD Cell df/f (probability of spiking per frame)RED CELLS
sdRedCelldf=std(Fnew(idd2,:)');
GrandsdRedCelldf=mean(sdRedCelldf)

%% SD Cell spikes (probability of spiking per frame)Other
sdOtherCellSpike=std(spike_rates(Other2,33:end-32)');
GrandsdOtherCellSpike=mean(sdOtherCellSpike)

%% SD Cell df/f (median of calcium trasces) Other
OtherCelldf=std(Fnew(Other2,:)');
sdOtherCelldf=mean(OtherCelldf)

%% Calculate spike probability of Red cells and Other in H-events

mean_tempHevents=[];
for i=1:length(long2);
    tempHevents=spike_rates(Other2,seq{long2(i)});
    mean_tempHevents(i)=mean(tempHevents(:));
end

grandMeanHevents=mean(mean_tempHevents);

mean_tempHeventsRED=[];
for i=1:length(long2);
    tempHeventsRED=spike_rates(idd2,seq{long2(i)});
    mean_tempHeventsRED(i)=mean(tempHeventsRED(:));
end

grandMeanHeventsRED=mean(mean_tempHeventsRED);
