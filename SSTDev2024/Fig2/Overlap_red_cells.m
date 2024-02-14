%%start manual selection
manual_selectionrrrrr

clear all;

%[FileName,PathName,FilterIndex]=uigetfile('*.tif');
%addpath(PathName);
%tic; RS = read_file(FileName); toc;
%aa2=size(RS);
aa2=[1024,1024];

%load selecion from red image
[FileName,PathName,FilterIndex]=uigetfile('*.mat');
addpath(PathName);
load(FileName);

%load all file somas only
[FileName2,PathName2,FilterIndex2]=uigetfile('*.mat');
addpath(PathName2);
load(FileName2,'ops');
load(FileName2,'Coor');
load(FileName2,'iscell');
load(FileName2,'ab2');
load(FileName2,'ab');
load(FileName2,'spikenums');
load(FileName2,'locs');
load(FileName2,'CellCl');
load(FileName2,'PCl');
load(FileName2,'t11');
load(FileName2,'spike_rates');
load(FileName2,'thresh');

Cn=ops.meanImg;
aa1=size(Cn);
iks = mldivide(aa1(2),aa1(1));

d1=aa1(1);
d2=aa1(2);

for i=1:length(sel)
tmp2=sel(i,:);
tempI=floor(tmp2(2)*iks);
tempJ=floor(tmp2(1)*iks);
tempI(tempI+4>aa1(1))=aa1(1)-4;
tempJ(tempJ+4>aa1(1))=aa1(1)-4;
tempI(tempI-4<0)=4;
tempJ(tempJ-4<0)=4;
BW1=zeros(d1,d2);
for k=0:9
for kk=0:9
BW1(tempI+(k-4),tempJ+(kk-4))=1;
end
end
dwns{i}=BW1;
end

figure;
imagesc(ops.meanImg)
hold on;
colormap('gray');

for i=1:length(dwns)
BW2 = dwns{i};
visboundaries(BW2,'LineWidth',0.2,'Color','m')
end

for i=1:length(Coor)
cont = medfilt1(Coor{i}')';
for h=1:size(cont,2)
BW1(cont(2,h),cont(1,h))=1;
end
visboundaries(BW1,'LineWidth',0.2,'Color','g')
end

SS = FileName(1:end-4);
SSS='Overlap_SST.fig';
St=strcat(SS,SSS);
savefig(St);
d1=aa1(1);
d2=aa1(2);
BW4=[];

for i=1:length(Coor)
i
BW1=zeros(d1,d2);
cont = medfilt1(Coor{i}')';
cont(cont==0)=1;
for h=1:size(cont,2)
BW1(cont(2,h),cont(1,h))=1;
end
for j=1:length(dwns)
BW2=[];
BW2 = dwns{j};
BW3=BW2.*BW1;
BW4(i,j)=sum(sum(BW3));
end
end

a=[];

for i=1:length(BW4(:,2));
a(i)=sum(BW4(i,:));
if a(i)>35
idd(i)=1;
b555=find(BW4(i,:)==max(BW4(i,:)));
b1(i)=b555(1);
else
idd(i)=0;
b1(i)=0;
end
end

mn1=[];
mn2=[];
mmn1=[];
mmn2=[];
figure;
hold on;
for i = 1:length(Coor)
if idd(i)==1
scatter(ab2(i),ab(i),'r')
mn1(end+1)=ab2(i);
mmn1(end+1)=ab(i);
end
if idd(i)==0
scatter(ab2(i),ab(i),'b')
mn2(end+1)=ab2(i);
mmn2(end+1)=ab(i);
end
end

pd = fitdist(mn2','Normal');
x_pdf = [-20:0.5:20];
y = pdf(pd,x_pdf);
pd2 = fitdist(mn1','Normal');
x_pdf2 = [-20:0.5:20];
y2 = pdf(pd2,x_pdf2);
hold on
scale = max(mmn2);
scale2=max(mmn1);
h=area((x_pdf),(y.*scale))
h(1).FaceColor = [0 0 1];
alpha(h(1),0.2)
h=area((x_pdf2),(y2.*scale2));
h(1).FaceColor = [1 0 0];
alpha(h(1),0.2)

SS = 'red';
SSS='OverlapScatter_SST.fig';
St=strcat(SS,SSS);
savefig(St);

mn1=mn1(~isnan(mn1));
mn2=mn2(~isnan(mn2));

[h,p,ci,stats]=ttest2(mn1,mn2)
ordercells=mean(mn2)-mean(mn1)

SS = FileName(1:end-4);
SSS='stats_SST.mat';
St=strcat(SS,SSS);
save(St,'h','p','ci','stats','ordercells');

%%sparseness
temp1=spikenums;
for i=1:length(locs)
temp1(:,locs(i)-40:locs(i)+40)=0;
end
temp1sum=sum(sum(temp1));
sparsenessFactor=(temp1sum/((length(spikenums(1,:))-length(locs)*7)*length(spikenums(:,1))))*100

SS = FileName(1:end-4);
SSS='sparsegeneral.mat';
St=strcat(SS,SSS);
save(St,'sparsenessFactor')

%%redcells in clusters
for i=1:length(CellCl(:,1))
cluu{i}=find(CellCl(i,:)==1);
end

idd2=find(idd==1);

for i=1:length(CellCl(:,1))
ovlclus{i}=ismember(cluu{i},idd2);
end

try
for i=1:length(CellCl(:,1))
cluu2{i}=find(PCl(i,:)==1);
end
catch
cluu2=0;
end

try
for i=1:length(CellCl(:,1))
tnew=[];
aaa=cluu2{i};
bbb=cluu{i};
tnew=t11(aaa,:);
tnew=tnew(:,bbb);
tclus{i}=tnew;
end
catch
tclus=0;
end

try
for j=1:length(CellCl(:,1))
t13=tclus{j};
abb=[];
abb2=[];
for i=1:length(t13(1,:))
abb(i)=((length(t13(:,i))-sum(isnan(t13(:,i))))/length(cluu2{j}))*100;
abb2(i)=nanmean(t13(:,i));
end
mab{j}=abb;
mab2{j}=abb2;
end
catch
end

clsnormalab=[];
clsnormalab2=[];
clsredab=[];
clsredab2=[];
lp=ceil(sqrt(length(CellCl(:,1))));

for i=1:length(CellCl(:,1))
subplot(lp,lp,i)
hold on
atmp=ovlclus{i};
abtmp=mab{i};
ab2tmp=mab2{i};
for j=1:length(atmp)
if atmp(j)==0
scatter(ab2tmp(j),abtmp(j),'b')
clsnormalab2(end+1)=ab2tmp(j);
clsnormalab(end+1)=abtmp(j);
else
scatter(ab2tmp(j),abtmp(j),'r')
clsredab2(end+1)=ab2tmp(j);
clsredab(end+1)=abtmp(j);
end
end
end

redinclus=[];
allinclus=0;

for i=1:length(ovlclus)
redinclus=[redinclus find(ovlclus{1,i}==1)];
allinclus=allinclus+length(ovlclus{1,i});
end

length(redinclus)
allinclus=allinclus-length(redinclus);
allinclus

SS = FileName(1:end-4);
SSS='allClusters.fig';
St=strcat(SS,SSS);
savefig(St)
SS = FileName(1:end-4);
SSS='allClusters.mat';
St=strcat(SS,SSS);
save(St,'clsnormalab2','clsnormalab','clsredab2','clsredab','cluu','cluu2','idd2','tclus','redinclus','allinclus','-mat','-v7.3');
%save(St,'idd2','-mat','-v7.3');

%%redcycle
ovls=find(idd==1);
for i=1:length(ovls)
rcls(i,:)=spikenums(ovls(i),:);
end

%find events
for i=1:length(rcls)
y(i)=sum(rcls(:,i));
end

%%make surrogate distribution
for j=1:500
thresh2=zeros(length(rcls(:,1)),length(rcls(1,:)));
parfor i=1:length(rcls(:,1))
thresh2(i,:)=circshift(rcls(i,:), randi([1, length(rcls(1,:))]));
end
summ=sum(thresh2,1);
summ2{j}=summ;
end

tress = cell2mat(summ2);
thresh2 = prctile(tress,99);
[pks2,locs2] = findpeaks(y,'MinPeakDistance',40,'MinPeakHeight',thresh2);
redmeanCycle = mean(diff(locs2))*(1/30)
allCycle = mean(diff(locs))*(1/30)
figure;
imagesc(rcls)
hold on;
plot(y,'color','r')
SS = FileName(1:end-4);
SSS='redcycle.mat';
St=strcat(SS,SSS);
save(St,'redmeanCycle','allCycle','locs','locs2','-mat','-v7.3');
close all

%%corr vs distance
xcen=[];
ycen=[];

for i=1:length(spike_rates(:,1))
try
cnttemp=Coor{i};
xcen(i)=mean(cnttemp(2,:));
ycen(i)=mean(cnttemp(1,:));
catch
xcen(i)=250;
ycen(i)=250;
end
end

nn=1:length(spike_rates(:,1));
novl=setdiff(nn,ovls);
dxcen=xcen(ovls);
dycen=ycen(ovls);
ndxcen=xcen(novl);
ndycen=ycen(novl);
dspike_rates=spike_rates(ovls,:);
ndspike_rates=spike_rates(novl,:);
dDistMap=zeros(length(dxcen),length(dxcen));

for i=1:length(dxcen)
for j=1:length(dxcen)
X = [dxcen(i),dycen(i);dxcen(j),dycen(j)];
d = pdist(X,'euclidean');
dDistMap(i,j)=d;
end
end

ndDistMap=zeros(length(ndxcen),length(ndxcen));
for i=1:length(ndxcen)
parfor j=1:length(ndxcen)
X = [ndxcen(i),ndycen(i);ndxcen(j),ndycen(j)];
d = pdist(X,'euclidean');
ndDistMap(i,j)=d;
end
end

dCorrMap=zeros(length(dxcen),length(dxcen));
for i=1:length(dxcen)
for j=1:length(dxcen)
X = dspike_rates(i,:);
Y = dspike_rates(j,:);
d = nancorr(X',Y');
dCorrMap(i,j)=d;
end
end

ndCorrMap=zeros(length(ndxcen),length(ndxcen));
for i=1:length(ndxcen)
for j=1:length(ndxcen)
X = ndspike_rates(i,:);
Y = ndspike_rates(j,:);
d = nancorr(X',Y');
ndCorrMap(i,j)=d;
end
end

dDistMap=triu(dDistMap);
ndDistMap=triu(ndDistMap);
dCorrMap=triu(dCorrMap);
ndCorrMap=triu(ndCorrMap);
dCorrMap = dCorrMap - diag(diag(dCorrMap));
ndCorrMap = ndCorrMap - diag(diag(ndCorrMap));
mask = triu(true(size(dCorrMap)),1);
dCM = dCorrMap(mask);
dDM=dDistMap(mask);
mask = triu(true(size(ndCorrMap)),1);
ndCM = ndCorrMap(mask);
ndDM=ndDistMap(mask);
scatter(ndDM,ndCM,'b');
hold on;
scatter(dDM,dCM,'r');
P1 = polyfit(ndDM,ndCM,1);
yfit1 = P1(1)*ndDM+P1(2);
hold on;
plot(ndDM,yfit1,'black-.');
P2 = polyfit(dDM,dCM,1);
yfit2 = P2(1)*dDM+P2(2);
hold on;
plot(dDM,yfit2,'r-.');
slope_cells=P1(1)
slope_redcells=P2(1)

SS = FileName(1:end-4);
SSS='corr_vd_dist.mat';
St=strcat(SS,SSS);
save(St,'ndCM','ndDM','dCM','dDM','P1','P2','-mat','-v7.3');
SS = FileName(1:end-4);
SSS='corr_vs_dist.fig';
St=strcat(SS,SSS);
savefig(St);
close all

%% Calculate mean correlation
CorrMap=zeros(size(spike_rates,1),size(spike_rates,1));

for i=1:size(spike_rates,1)
    for j=i:size(spike_rates,1)
    X = spike_rates(i,:);
    Y = spike_rates(j,:);
    d = nancorr(X',Y');
    CorrMap(i,j)=d;
    end
end

CorrMap = CorrMap - diag(diag(CorrMap));

mask = triu(true(size(CorrMap)),1);
CM = CorrMap(mask);
MeanCM=mean(CM) % mean correlation between all cells

SS = FileName(1:end-4);
SSS='corr.mat';
St=strcat(SS,SSS);
save(St,'CM','MeanCM','-mat','-v7.3');

%%
for i=1:length(spikenums)
y(i)=sum(spikenums(:,i));
end

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

pd3=fitdist(len_seq','Exponential');
pd_ci = paramci(pd3);
aaa=pd_ci(2)-pd3.mu;
large=pd3.mu+(3*aaa);
histfit(len_seq,25,'exponential')

SS = FileName(1:end-4);
SSS='Events_histo.fig';
St=strcat(SS,SSS);
savefig(St);