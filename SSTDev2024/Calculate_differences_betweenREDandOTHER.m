clear
%load Fallall
[FileName3,PathName3,FilterIndex3]=uigetfile('*.mat');
addpath(PathName3);
load(FileName3,'spikenums');
load(FileName3,'thresh');

%load Selectionall Clusters
[FileName4,PathName4,FilterIndex4]=uigetfile('*.mat');
addpath(PathName4);
load(FileName4,'idd2');

for i=1:length(spikenums)
y(i)=sum(spikenums(:,i));
end
y = medfilt1(y,10);
up=find(y>thresh);
A=up;
A(end+1)=2;   % adds new endpoint to very end of A so code picks up end of last group of consecutive values
I_1=find(diff(A)~=1);  % finds where sequences of consecutive numbers end
[m,n]=size(I_1);   % finds dimensions of I_1 i.e. how many sequences of consecutive numbers you have
startpoint=1;    % sets start index at the first value in your array
seq=cell(1,n);  % had to preallocate because without, it only saved last iteration of the for loop below

% used n because this array is a row vector
for i=1:n
End_Idx=I_1(i);   %set end index
seq{i}=A(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
startpoint=End_Idx+1;   %update start index for the next consecutive sequence
end

for i=1:length(seq)
len_seq(i)=length(seq{i});
end

IEI=600/length(len_seq)
length(len_seq)% # of events

participation=zeros(length(seq),length(spikenums(:,1)));
for i=1:length(seq)
aaa=[];
bbb=[];
aaa=seq{i};
bbb=spikenums(:,aaa);
if length(aaa)>1
participation(i,:)=sum(bbb')/length(aaa);
else
participation(i,:)=bbb;
end
end

participation=participation';
firing_red=participation(idd2,:);
nn=1:length(spikenums(:,1));
novl=setdiff(nn,idd2);

firing_nonred=participation(novl,:);
large=median(len_seq)+(2*iqr(len_seq));
long=find(len_seq>=large);
short=find(len_seq<large);

firing_red_long=participation(idd2,long);
firing_red_short=participation(idd2,short);
firing_nonred_long=participation(novl,long);
firing_nonred_short=participation(novl,short);
red_l=reshape(firing_red_long,[],1);
red_s1=reshape(firing_red_short,[],1);
nonred_l=reshape(firing_nonred_long,[],1);
nonred_s=reshape(firing_nonred_short,[],1);
x = [red_l; red_s1; nonred_l; nonred_s];
g = [zeros(length(red_l'), 1); ones(length(red_s1'), 1); 2*ones(length(nonred_l'), 1) ; 3*ones(length(nonred_s'), 1)];
boxplot(x, g,'Labels',{'red_long','red_short','nonred_long','nonred_short'})
SS = FileName3(1:end-4);
SSS='participationbylength.fig';
St=strcat(SS,SSS);
savefig(St)
close all

maxpeakcells=zeros(length(seq),1);
for i=1:length(seq)
aaa=[];
bbb=[];
ccc=[];
aaa=seq{i};
bbb=spikenums(:,aaa);
ccc=sum(bbb');
ccc=ccc';
ccc=find(ccc>0);
ccc=length(ccc)/length(spikenums(:,1));
maxpeakcells(i)=ccc;
end

large2=0.8;
long2=find(maxpeakcells>=large2);
short2=find(maxpeakcells<large2);

firing_red_long2=participation(idd2,long2);
firing_red_short2=participation(idd2,short2);
firing_nonred_long2=participation(novl,long2);
firing_nonred_short2=participation(novl,short2);

red_l2=reshape(firing_red_long2,[],1);
red_s2=reshape(firing_red_short2,[],1);

nonred_l2=reshape(firing_nonred_long2,[],1);
nonred_s2=reshape(firing_nonred_short2,[],1);

x = [red_l2; red_s2; nonred_l2; nonred_s2];
g = [zeros(length(red_l2'), 1); ones(length(red_s2'), 1); 2*ones(length(nonred_l2'), 1) ; 3*ones(length(nonred_s2'), 1)];
boxplot(x, g,'Labels',{'red_high','red_low','nonred_high','nonred_low'})
SS = FileName3(1:end-4);
SSS='participationby80percent.fig';
St=strcat(SS,SSS);
savefig(St)
close all

tms=[];
for i=1:length(seq)
a=seq{i};
tms(i)=mean(a);
end
tms=tms*(1/30);
df_tms=diff(tms);
std_tms=std(df_tms);

figure;
histogram(df_tms,1:1:50,'Normalization','probability');%% percentage of events is normalized to the probability
SS = FileName3(1:end-4);
SSS='histogram_iei.fig';
St=strcat(SS,SSS);
savefig(St)

IEIlong=mean(diff(tms(long2)))
IEIshort=mean(diff(tms(short2)))
close all

%Fraction of long and short events:
f_long_events=((length(long2)/length(seq))*100)
f_short_events=((length(short2)/length(seq))*100)

%%participation in events
for i=1:length(spikenums(:,1))
a=find(participation(i,:)>0);
partinevents(i)=length(a)/length(seq);
end
part_events_red=partinevents(idd2);
part_events_nonred=partinevents(novl);
x = [part_events_red'; part_events_nonred'];
g = [zeros(length(part_events_red'), 1); ones(length(part_events_nonred'), 1)];
boxplot(x, g,'Labels',{'sst_participation','all_participation'})
SS = FileName3(1:end-4);
SSS='partcipationinevents.fig';
St=strcat(SS,SSS);
savefig(St)

sst_participation = mean(part_events_red)
all_participation = mean(part_events_nonred)
%sst_participation_std2 = std(part_events_red)
%all_participation_std2 = std(part_events_nonred)
save('participation.mat')
%%'redmeanhigh','redsthigh','redmeanlow','redstlow','allmeanhigh','allsthigh','allmeanlow','allstlow');
close all

partinevents2=participation;
partinevents2(partinevents2>0)=1;
for i=1:length(spikenums(:,1))
a=partinevents2(i,long2);
b=partinevents2(i,short2);
a=find(a>0);
b=find(b>0);
parthigh(i)=length(a)/(length(a)+length(b));
partlow(i)=length(b)/(length(a)+length(b));
end

firing_red_high=parthigh(idd2);
firing_red_low=partlow(idd2);
firing_nred_high=parthigh(novl);
firing_nred_low=partlow(novl);
x = [firing_red_high'; firing_red_low'; firing_nred_high'; firing_nred_low'];
g = [zeros(length(firing_red_high'), 1); ones(length(firing_red_low'), 1); 2*ones(length(firing_nred_high'), 1) ; 3*ones(length(firing_nred_low'), 1)];
boxplot(x, g,'Labels',{'red_high','red_low','nonred_high','nonred_low'})
SS = FileName3(1:end-4);
SSS='percentageinevents.fig';
St=strcat(SS,SSS);
savefig(St)

redmeanhigh = mean(firing_red_high)
%redsthigh = std2(firing_red_high)
redmeanlow = mean(firing_red_low)
%redstlow = std2(firing_red_low)
allmeanhigh = mean(firing_nred_high)
%allsthigh = std2(firing_nred_high)
allmeanlow = mean(firing_nred_low)
%allstlow = std2(firing_nred_low)
save('participation.mat','redmeanhigh','redmeanlow','allmeanhigh','allmeanlow','IEI', 'long2','short2', 'df_tms');

length(len_seq)
close all

%%corrrelation length vs participation
scatter(maxpeakcells,len_seq)
xlabel('maxpeak')
ylabel('length')
SS = FileName3(1:end-4);
SSS='peakvslength.fig';
St=strcat(SS,SSS);
savefig(St)
close all
SS = FileName3(1:end-4);
SSS='participation.mat';
St=strcat(SS,SSS);
save(St,'-v7.3');
clear all


