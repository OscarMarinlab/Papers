clear all

ii=1;
ana=0;

CellParticipation_injection={};
CellParticipation_Baseline={};
EventsIEI_injection={};
EventsIEI_Baseline={};
EventsLength_injection={};
EventsLength_Baseline={};

while ii<2
    
   clearvars -except CellParticipation_injection CellParticipation_Baseline EventsIEI_injection EventsIEI_Baseline  EventsLength_injection EventsLength_Baseline ii FileName3

prompt = 'Add case? ';
ii = input(prompt);
if ii<2
[FileName3,PathName3,FilterIndex3]=uigetfile('*.mat');
addpath(PathName3);


load(FileName3,'spikenums');
load(FileName3,'thresh');
spikenumsB=spikenums;
threshB=thresh;

[FileName3,PathName3,FilterIndex3]=uigetfile('*.mat');
addpath(PathName3);

load(FileName3,'spikenums');
load(FileName3,'thresh');

spikenumsI=spikenums;
threshI=thresh;


for i=1:length(spikenumsB)
    y(i)=sum(spikenumsB(:,i));
end

y = medfilt1(y,10);


up=find(y>threshB);

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

len_seq=[];

for i=1:length(seq)
    len_seq(i)=length(seq{i});
end

%IEI=600/length(len_seq)% #events

%%Frequency distribution
tms=[];
for i=1:length(seq)
    a=seq{i};
    tms(i)=mean(a);
end

tms=tms*(1/30);

df_tms_baseline=diff(tms);
std_tms_baseline=std(df_tms_baseline);


maxpeakcells_baseline=zeros(length(seq),1);
for i=1:length(seq)
    aaa=[];
    bbb=[];
    ccc=[];
    aaa=seq{i};
    bbb=spikenumsB(:,aaa);
    ccc=sum(bbb');
    ccc=ccc';
    ccc=find(ccc>0);
    ccc=length(ccc)/length(spikenumsB(:,1));
    maxpeakcells_baseline(i)=ccc;
end

%%Distribution of event Length
len_seqB=len_seq;

y=[];
for i=1:length(spikenumsI)
    y(i)=sum(spikenumsI(:,i));
end

y = medfilt1(y,10);


up=find(y>threshI);
seq={};
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

len_seq=[];

for i=1:length(seq)
    len_seq(i)=length(seq{i});
end


%%Frequency distribution
tms=[];
for i=1:length(seq)
    a=seq{i};
    tms(i)=mean(a);
end

tms=tms*(1/30);

df_tms_injection=diff(tms);
std_tms_injection=std(df_tms_injection);


maxpeakcells_injection=zeros(length(seq),1);
for i=1:length(seq)
    aaa2=[];
    bbb2=[];
    ccc2=[];
    aaa2=seq{i};
    bbb2=spikenumsI(:,aaa2);
    ccc2=sum(bbb2');
    ccc2=ccc2';
    ccc2=find(ccc2>0);
    ccc2=length(ccc2)/length(spikenumsI(:,1));
    maxpeakcells_injection(i)=ccc2;
end

%%Distribution of event Length
len_seqI=len_seq;


CellParticipation_injection{end+1}=maxpeakcells_injection;
CellParticipation_Baseline{end+1}=maxpeakcells_baseline;

EventsIEI_injection{end+1}=df_tms_injection;
EventsIEI_Baseline{end+1}=df_tms_baseline;

EventsLength_injection{end+1}=len_seqI;
EventsLength_Baseline{end+1}=len_seqB;


end
end

CellParticipation_injection2=cell2mat(CellParticipation_injection');
CellParticipation_Baseline2=cell2mat(CellParticipation_Baseline');
EventsIEI_injection2=cell2mat(EventsIEI_injection);
EventsIEI_Baseline2=cell2mat(EventsIEI_Baseline);
EventsLength_injection2=cell2mat(EventsLength_injection);
EventsLength_Baseline2=cell2mat(EventsLength_Baseline);

hist1= histogram(EventsIEI_Baseline2,1:1:50,'Normalization','probability','EdgeColor', 'red', 'FaceColor',  'red');
hold on
hist2 = histogram(EventsIEI_injection2,1:1:50,'Normalization','probability', 'EdgeColor', 'green', 'FaceColor',  'green', 'FaceAlpha', 0.2);
legend('Baseline','Injection')


SS = FileName3(1:end-4);
SSS='histogram_IEI.fig';
St=strcat(SS,SSS);
savefig(St);
close all

hist1= histogram(CellParticipation_Baseline2,20,'Normalization','probability','EdgeColor', 'red', 'FaceColor',  'red');
hold on
hist2 = histogram(CellParticipation_injection2,20,'Normalization','probability', 'EdgeColor', 'green', 'FaceColor',  'green', 'FaceAlpha', 0.2);
legend('Baseline','Injection')


SS = FileName3(1:end-4);
SSS='histogram_participation.fig';
St=strcat(SS,SSS);
savefig(St);
close all

hist1= histogram(EventsLength_Baseline2,1:5:250,'Normalization','probability','EdgeColor', 'red', 'FaceColor',  'red');
hold on
hist2 = histogram(EventsLength_injection2,1:5:250,'Normalization','probability', 'EdgeColor', 'green', 'FaceColor',  'green', 'FaceAlpha', 0.2);
legend('Baseline','Injection')


SS = FileName3(1:end-4);
SSS='histogram_lenght.fig';
St=strcat(SS,SSS);
savefig(St);
close all


%% Compare distributions (baseline vs. injection)
[h,p,ks2stat]=kstest2(EventsLength_Baseline2,EventsLength_injection2)

[h,p,ks2stat]=kstest2(CellParticipation_Baseline2,CellParticipation_injection2)

[h,p,ks2stat]=kstest2(EventsIEI_Baseline2,EventsIEI_injection2)

save('Comparison_events.mat');
