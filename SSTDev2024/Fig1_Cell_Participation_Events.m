Ag=9;  %Age of the Animal/change when necessary
Dep='Sup';  %either 'Sup' or 'Deep'

BasePath='F:\Calcium Imaging\tac1\Spontaneous_taken_from _veh';
a=Age==Ag;
tf1 = strcmp('Sup',Dep);
if tf1==1
b=Depth<=450;
else
b=Depth>=200;
end

cases=a&b;
cases=find(cases==1);

sumCasesCellParticipation={};
Mean_cell_participation_Events=[];
Median_cell_participation_Events=[];
std_cell_participation_Events=[];

for kk=1:length(cases)
    
  a=fileList(cases(kk)).folder;
    cd(a)      
load('Fallallparticipation.mat','participation','len_seq');
mean_cell_participation=zeros(1,length(len_seq));
participation2=participation;
participation2(participation2>0)=1;
participation2=sum(participation2);
participation2=participation2/size(participation,1);
Mean_cell_participation_Events(kk)=mean(participation2);
Median_cell_participation_Events(kk)=median(participation2);
std_cell_participation_Events(kk)=std(participation2);
sumCasesCellParticipation{kk}=participation2;

clearvars -except sumCasesCellParticipation participation2 kk cases fileList BasePath Mean_cell_participation_Events Median_cell_participation_Events std_cell_participation_Events
end

cd(BasePath)
SumCasesCellParticipation2_InjectionV_P8=cell2mat(sumCasesCellParticipation);%Change according to variable
histogram(SumCasesCellParticipation2_InjectionV_P8, 0:0.1:1,'Normalization','probability')%Change according to variable

Mean_cell_participation_InjectionV_P8=Mean_cell_participation_Events; %Change according to variable
Median_cell_participation_InjectionV_P8=Median_cell_participation_Events;%Change according to variable
std_cell_participation_InjectionV_P8=std_cell_participation_Events;%Change according to variable
sumCasesCellParticipation{kk}=participation2;

save('MeanMedianST_InjectionCNO_P11','Mean_cell_participation_InjectionCNO_P11','Median_cell_participation_InjectionCNO_P11','std_cell_participation_InjectionCNO_P11');

save('SumCasesCellParticipation2_InjectionCNO_P11')

%%histogram with pooled ages
SumCasesCellParticipation2_CNO_P11_12= cat(2,SumCasesCellParticipation2_InjectionCNO_P12, SumCasesCellParticipation2_InjectionCNO_P11);%Change according to variable
histogram(SumCasesCellParticipation2_CNO_P11_12, 0:0.02:1,'Normalization','probability')%Change according to variable

%%histogram with single cases (i.e CNO or Vehicle, not pooled ages)
histogram(SumCasesCellParticipation2_InjectionV_P8, 0:0.02:1,'Normalization','probability')%Change according to variable

save('SumCasesCellParticipation2_InjectionV_P8') % change name accordingly

[FileName3,PathName3,FilterIndex3]=uigetfile('*.mat');
addpath(PathName3);

SS = FileName3(1:end-4);
SSS='Histogrm_Cell_participationp_Baseline.fig';
St=strcat(SS,SSS);
savefig(St)
close all

%KS test between the different ages
[h1,p1] = kstest2(SumCasesCellParticipation2_InjectionCNO,SumCasesCellParticipation2_InjectionV_P8,'Alpha',0.05)
[h2,p2] = kstest2(SumCasesCellParticipation2_p7_8,SumCasesCellParticipation2_p11_12,'Alpha',0.05)
[h3,p3] = kstest2(SumCasesCellParticipation2_p9_10,SumCasesCellParticipation2_p11_12,'Alpha',0.05)

save('KSstats_ComparisonCNO_VS_vehicle','h1','p1','h2','p2','h3','p3')

%Find peaks in the distributions
KSdensity_InjectionCNO=ksdensity(SumCasesCellParticipation2_InjectionCNO);
KSdensity_p11_12=ksdensity(SumCasesCellParticipation2_p11_12_V);
KSdensity_p11_12=ksdensity(SumCasesCellParticipation2_p11_12);

plot(KSdensity_InjectionCNO)

[pks1,locs1]=findpeaks(KSdensity_InjectionCNO)
[pks2,locs2]=findpeaks(KSdensity_p9_10)
[pks3,locs3]=findpeaks(KSdensity_p11_12)

save('Peaks_InjectionV','pks1','locs1','pks2','locs2','pks3','locs3')

TF1 = islocalmin(KSdensity_InjectionV);
find(TF1==1)

TF2 = islocalmin(KSdensity_p9_10);
find(TF2==1)

TF3 = islocalmin(KSdensity_p11_12);
find(TF3==1)

%Minimum values in the data
save('Minimum','TF1','TF2','TF3')



