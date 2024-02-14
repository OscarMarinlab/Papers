
%% Load data
clear
[FileName,PathName,FilterIndex]=uigetfile('*.mat');
addpath(PathName);
load 'Fall.mat';
cd (PathName);


%% Make partition of longitudinal data
% Make partition (change according to the number of days)
day{1}=F(:,1:end);
day{2}=F(:,18001:36000);
day{3}=F(:,36001:54000);
day{4}=F(:,54001:72000);
day{5}=F(:,72001:90000);

dayneu{1}=Fneu(:,1:end);
dayneu{2}=Fneu(:,18001:36000);
dayneu{3}=Fneu(:,36001:54000);
dayneu{4}=Fneu(:,54001:72000);
dayneu{5}=Fneu(:,72001:90000);

for ii=1:length(day)
F1=day{ii};
Fneu1=dayneu{ii};

Fnew=zeros(size(F1));
for i=1:length(iscell(:,1))
    Fnew(i,:)=(F1(i,:)-0.7*Fneu1(i,:))/median(Fneu1(i,:));
end

Fnew(Fnew==Inf)=0;
Fnew(Fnew==-Inf)=0;

dF_traces=Fnew;  % Rename for cascade
save 'D:/matlab.mat' dF_traces
[status] = system('python C:/Users/User/Downloads/Cascade/Demo_predict.py')

load('D:/full_prediction_matlab.mat.mat');

%To use discrete spike inference uncomment:
%[status] = system('python C:/Users/User/Downloads/Cascade/Demo_discrete_spikes.py','-echo')

%load('D:/discrete_spikes_matlab.mat');

%att=spike_time_estimates;
%savename=strcat(PathName,'spikes',num2str(ii),'.mat');
%save(savename, 'spike_rates','discrete_approximation','spike_time_estimates') 

%To use thresholded method uncomment:

savename=strcat(PathName,'spikes',num2str(ii),'.mat');
save(savename, 'spike_rates')  
clearvars spike_rates
end

%% make spikes
%set up threshold for each cell( change if necessary)
spikenums=spike_rates;

%spikenums(spikenums>0.3)=1;
%spikenums(spikenums<1)=0;

%%Spikes with 95:
for i=1:length(iscell(:,1))
    tr=prctile(spike_rates(i,:),95);
    sp=spike_rates(i,:);
    sp(sp<tr)=0;
    sp(sp>=tr)=1;
    sp(isnan(sp))=0;
    spikenums(i,:)=sp;
end

%set up threshold for all( change if necessary)
%tr=prctile(spike_rates(:),95);
%spikenums=spike_rates;;
%spikenums(spikenums<tr)=0;
%spikenums(spikenums>=tr)=1;
%spikenums(isnan(spikenums))=0;
for i=1:length(Fneu(:,1)) 
    aa=[];
    aa(1,:)=stat{1,i}.xpix;
    aa(2,:)=stat{1,i}.ypix;
    Coor{i}=aa;
end

a1=iscell(:,1)==1;
a1=find(iscell(:,1)==1);
spikenums=spikenums(a1,:);
spike_rates=spike_rates(a1,:);
Coor=Coor(a1);
%CellCl=CellCl(:,a1);


%% make surrogate distribution
spikenums22=spikenums(:,1:end);
for j=1:500
thresh=zeros(length(spike_rates(:,1)),length(spikenums22(1,:)));
parfor i=1:length(spike_rates(:,1))
thresh(i,:)=circshift(spikenums22(i,:), randi([1, length(spikenums22(1,:))]));
end
summ=sum(thresh,1);
summ2{j}=summ;
end
tress = cell2mat(summ2);
thresh = prctile(tress,99);

for i=1:length(spikenums)
    y(i)=sum(spikenums(:,i));
end


y = medfilt1(y,10);
[pks,locs] = findpeaks(y,'MinPeakDistance',40,'MinPeakHeight',thresh);

figure;
imagesc(spikenums)
hold on;
plot(y,'color','r')

findpeaks(y,'MinPeakDistance',40,'MinPeakHeight',thresh);
SS = FileName(1:end-4);
SSS='Events.fig';
St=strcat(SS,SSS);
savefig(St)

%% postprocessing
lb=40; %frames before
hb=40; %frames after vector of concatenated vector should be center at 0
if (locs(end)+hb)>length(spikenums(1,:))
    locs=locs(1:end-1);
    pks=pks(1:end-1);
end

if (locs(1)-lb)<1
    locs=locs(2:end);
    pks=pks(2:end);
end

att=length(find(iscell(:,1)==1));

t1=NaN(length(locs),att);
t2=zeros(length(locs),att);


for i=1:length(pks)
    for j=1:att
        axxx=spikenums(j,locs(i)-lb:locs(i)+hb);
        ax=min(find(axxx==1));
        if ax>=1
            t1(i,j)=ax;
        end
    end
end

tttt=round(length(-lb:hb)/2);
t11=t1-tttt;
t12=t11;
t12(~isnan(t11))=1;
t12(isnan(t12))=0;

%% ab and ab2 variables
t1=NaN(length(locs),length(spikenums(:,1)));
t2=zeros(length(locs),length(spikenums(:,1)));


for i=1:length(pks)
    for j=1:length(spikenums(:,1))
        axxx=spikenums(j,locs(i)-9:locs(i)+9);
        ax=min(find(axxx==1));
        if ax>=1
            t1(i,j)=ax;
        end
    end
end

t11=t1-9;


for i=1:length(t11(1,:))
    ab(i)=((length(t11(:,i))-sum(isnan(t11(:,i))))/length(locs))*100;
    ab2(i)=nanmean(t11(:,i));
end

clearvars t1 t2 tttt axxx ax lb hb summ summ2 tress

%%
%start cluster algo, run time should probably be improved
Race = t12';
[NCell,NRace] = size(Race);
[IDX2,sCl,M,S] = kmeansopt(Race,100,'var');
% M = CovarM(Race);
% IDX2 = kmedoids(M,NCl);
NCl = max(IDX2);

[~,x2] = sort(IDX2);
MSort = M(x2,x2);

%Race clusters
R = cell(0);
CellScore = zeros(NCell,NCl);
CellScoreN = zeros(NCell,NCl);
for i = 1:NCl
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end

%Assign cells to cluster with which it most likely spikes
[~,CellCl] = max(CellScoreN,[],2);
%Remove cells with less than 2 spikes in a given cluster
CellCl(max(CellScore,[],2)<2) = 0;
[X1,x1] = sort(CellCl);

figure
subplot(1,2,1)
imagesc(MSort)
colormap jet
axis image
xlabel('RACE #')
ylabel('RACE #')

subplot(1,2,2)
imagesc(Race(x1,x2),[-1 1.2])
axis image
xlabel('RACE #')
ylabel('Cell #')



SS = FileName(1:end-4);
SSS='Clusters.fig';
St=strcat(SS,SSS);

savefig(St)

%% Remove cluster non-statistically significant

sClrnd = zeros(1,1000);
for i = 1:1000
    sClrnd(i) = kmeansoptrnd(Race,10,NCl);
end
% TODO: replace max by percentile (done)
NClOK = sum(sCl>prctile(sClrnd,95));
sClOK = sCl(1:NClOK)';

%save('NClustersOK.mat','NClOK')

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);

%% Statistical definition of cell assemblies

NShuf = 5000;
%Count number of participation to each cluster
CellP = zeros(NCell,NCl); CellR = zeros(NCell,NCl);
for i = 1:NCl
    CellP(:,i) = sum(Race(:,IDX2 == i),2);
    CellR(:,i) = CellP(:,i)/sum(IDX2 == i);
end

%Test for statistical significance
CellCl = zeros(NCl,NCell); %Binary matrix of cell associated to clusters
for j = 1:NCell
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(j,:));
    parfor l = 1:NShuf
        Random = randperm(NRace);
        Random = Random(1:Nrnd);
        Racer = zeros(1,NRace);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(:,IDX2 == i),2);
        end
    end
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl))); 
    for i = 1:NCl
        CellCl(i,j) = double(CellP(j,i)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
end

A0 = find(sum(CellCl) == 0); %Cells not in any cluster
A1 = find(sum(CellCl) == 1); %Cells in one cluster
A2 = find(sum(CellCl) >= 2); %Cells in several clusters

%Keep cluster where they participate the most
for i = A2
    [~,idx] = max(CellR(i,:));
    CellCl(:,i) = 0;
    CellCl(idx,i) = 1;
end
C0 = cell(0);
k = 0;
for i = 1:NCl
    if length(find(CellCl(i,:)))>5
        k = k+1;
        C0{k} = find(CellCl(i,:));
    end
end

SS = FileName(1:end-4);
SSS='CellList.mat';
St=strcat(SS,SSS);
save(St,'C0');

%Participation rate to its own cluster
SS = FileName(1:end-4);
SSS='ClusPart.mat';
St=strcat(SS,SSS);
CellRCl = max(CellR([A1 A2],:),[],2);
save(St,'CellRCl')
%% Assign RACE to groups of cells

NCl = length(C0);
[NCell,NRace] = size(Race);
%Cell count in each cluster
RCl = zeros(NCl,NRace);
PCl = zeros(NCl,NRace);
for i = 1:NCl
    RCl(i,:) = sum(Race(C0{i},:));
end

try
RCln = zeros(NCl,NRace);
for j = 1:NRace
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(:,j));
    parfor l = 1:NShuf
        Random = randperm(NCell);
        Random = Random(1:Nrnd);
        Racer = zeros(NCell,1);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(C0{i}));
        end
    end
    %         ThMin = mean(Random) - 2*std(Random);
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));
    for i = 1:NCl
        PCl(i,j) = double(RCl(i,j)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
    %Normalize (probability)
    RCln(:,j) = RCl(:,j)/sum(Race(:,j));
end
catch
    PCl=0;
end

SS = FileName(1:end-4);
SSS='CellCluster.mat';
St=strcat(SS,SSS);

save(St,'PCl');

%% Show Sorted Rasterplot

%load('RaceCl2')
%load('CellClList2.mat')
%load('Race.mat')
[NCell,NRace] = size(Race);
NCl = length(C0);

%Recreate CellCl (equivalent of RCl for the cells)
CellCl = zeros(NCl,NCell);
for i = 1:NCl
    CellCl(i,C0{i}) = 1;
end

NCl = length(C0);

Cl0 = find(sum(PCl,1) == 0);
Cl1 = find(sum(PCl,1) == 1);
Cl2 = find(sum(PCl,1) == 2);
Cl3 = find(sum(PCl,1) == 3);
Cl4 = find(sum(PCl,1) == 4);

Bin = 2.^(0:NCl-1);

%Sort Cl1
[~,x01] = sort(Bin*PCl(:,Cl1));
Cl1 = Cl1(x01);

%Sort Cl2
[~,x02] = sort(Bin*PCl(:,Cl2));
Cl2 = Cl2(x02);

%Sort Cl3
[~,x03] = sort(Bin*PCl(:,Cl3));
Cl3 = Cl3(x03);

RList = [Cl0 Cl1 Cl2 Cl3 Cl4];
%x1 from DetectRace

[X1,x1] = sort(Bin*CellCl);

figure
imagesc(Race(x1,RList))
colormap hot
axis image

SS = FileName(1:end-4);
SSS='CellCluster.fig';
St=strcat(SS,SSS);
savefig(St)


SS = FileName(1:end-4);
SSS='all.mat';
St=strcat(SS,SSS);

save(St,'-v7.3');

% Convert to binary
spks2=single(imbinarize(spks));
imagesc(spks2)
savefig('RasterAll');