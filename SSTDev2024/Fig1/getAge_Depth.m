clear 
fileList = dir('**/Fallall.mat');

Depth=[];
Age=[];

for i=1:124 %124 || Change according to sample size
    startIndex=[];
    a=fileList(i).folder;
    a = lower(a);
    startIndex = regexp(a,'um');
    startIndex=startIndex(2);
    a1=startIndex-3;
    a2=startIndex-1;
    Depth(i)=str2num(a(a1:a2));
    
    startIndex={};
    startIndex{1} = regexp(a,'p7');
    startIndex{2} = regexp(a,'p8');
    startIndex{3} = regexp(a,'p9');
    startIndex{4} = regexp(a,'p10');
    startIndex{5} = regexp(a,'p11');
    startIndex{6} = regexp(a,'p12');
    
    tf = cellfun('isempty',startIndex);
    startIndex(tf) = {0};
    start=cell2mat(startIndex);
    b=max(start);
    c=find(start==b);
    c=c+6;
    Age(i)=c;
end

Age=Age';
Depth=Depth';

clearvars -except fileList Depth Age