close all, clear all; 
filename = '2017_09_07_0012.atf'; %write filename 
[~,NAME,EXT] = fileparts(filename); 
%%
A = importdata(filename, '\t', 11); %import data
current = A.data;
cellcode =NAME; 
cellspot =31; % closest spot where cell soma is
maptype = 5.; %laser intensity
layer = 'L2_'; %what layer is the pyramid in? 

%%
%loop thru

for c=1:42
    data(c, :)  = sCRACM_v2(current, c); 
end
%%

%%%MAP EVENTS INTO GRID NOW%%%
MapAVG(:,1) = data.spot; 
MapAVG(:,2) = data.mean_Amp; 
%recreates map grid pattern
grid = [1,11,7,3,35,9,40;13,23,17,21,15,25,19;27,33,39,31,37,29,5;6,2,12,8,4,42,10;20,26,14,24,18,22,16;34,30,41,28,36,32,38];
%%
%fills a new variable, map, with amplitudes from average map in correct spatial location
%spatial location

map=zeros(6,7);
[~,lingrid] = ismember(MapAVG(:,1),grid);
MapAVG(:,3)=lingrid;
gridmap=sortrows(MapAVG,3);
imap=42;
while imap>0
map(imap)=gridmap(imap,2);
imap=imap-1;
end;
%%
%%Plot raw map
climsraw=[0 round(max(map(:)),-2)];
imagesc(map, climsraw);colormap(hot);axis image; axis off;   
[rowcell,colcell]=find(grid==cellspot);hold on;
plot(colcell,rowcell,'b.','MarkerSize',40);
colorbar 
%text(0.5,1,cellcode,'Color','w','FontSize',14);
%text(0.5,2,layer,'Color','w','FontSize',14);
%text(0.5,18,'pA','Color','w','FontSize',14);
scale=num2str(climsraw); text(7,18,scale,'Color','w','FontSize',14);
%export pA map png file and matrix
%export_fig(strcat(layer,cellcode,'_pA'));
saveas(gcf,strcat(layer,cellcode,'_pA','.png')); 
dlmwrite(strcat(cellcode,'_pA','.csv'),map);
hold off

%%Calculate chessboard distance between each spot and cell body
distgrid=zeros(6,7);
distgrid(rowcell,colcell)=1;
dist=70*bwdist(distgrid,'chessboard');
distfreqs(:,1)=unique(dist);
distfreqs(:,2)=histc(dist(:),distfreqs(:,1));

dplot(:,1)=MapAVG(:,1);
dplot(:,2)=MapAVG(:,2);
Y=1;
X=[rowcell,colcell];

while Y<43
[T,U]=find(grid==dplot(Y,1));
X(2,1)=T;;X(2,2)=U;
dplot(Y,4)=pdist(X,'chebychev')*70;
Y=Y+1;
clear T; clear U;
end;


