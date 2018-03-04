%Mafb fl/fl analyses with SST:RCE

function mafb_e16(filename, brainno, age, geno, in, px) 
%%
[~,imagename, ~] = fileparts(filename); 
Brainno = brainno; %define the root of your file here
Age = age;
Genotype = geno; 
px = px^2; 

%name the files you are interested in analyzing
s = struct('inputfile', {strcat(imagename, '.tif')}, 'outputfile', strcat(imagename, 'mafb.csv')); 
%se1 = strel('line', 10, 0);
%se1b = strel('line', 10, 40);
se2 = strel('disk', 1, 8);
se5 = strel(ones(1,1));
se3 = strel('disk', 2, 8);
se4 = strel ('disk',4);

% read images and background subtract 
Ich1 = imread(s(1).inputfile, in); %channel reporter staining
bkg1 = mean(mode ((double(Ich1))));
Ich1_s = Ich1 - bkg1;

%defining MZ index
[r , c] = size (Ich1);
r = int16(r);
Imz = round(r*0.2);
Isvz = round(r*0.3);
%%
%ch1 image segmentation, reporter channel
Ich1_ed = imerode (Ich1_s, se3); 
Ich1_ed = imdilate (Ich1_ed, se3); 
Ich1_ed_d = imreconstruct (Ich1_ed, Ich1_s);
c1_fgm = imregionalmax(Ich1_ed_d); 
c1_fgm = imfill(c1_fgm, 'holes'); 
c1_fgm2 = bwareaopen (c1_fgm, 400);
c1_fgm = c1_fgm - c1_fgm2; 
[~, bw_ch1] = thresh_tool(Ich1); 
bw_ch1 = bwareaopen (bw_ch1, 20); 
ch1_total = bw_ch1 & c1_fgm;
ch1_total = imdilate(ch1_total, se2); 


%identify centroids
stats = regionprops(ch1_total, 'Centroid');
centroids = cat(1, stats.Centroid);
cr = round (centroids); 
x1 = cr (:, 1);
y1 = cr (:, 2); 
num_objects = length (stats); 
indcent = sub2ind((size(ch1_total)), y1, x1); 
centroidbw2 = false (size (ch1_total)); 
centroidbw2(indcent) = true; 
ch1_c = centroidbw2 & ch1_total; 
ch1_cd = imdilate (ch1_c, se4); 
%%
%defining GFP+ cells in MZ
MZ = ch1_c;
MZ(Imz:r, :) = 0;
SVZ = ch1_c - MZ; 
MZ = imdilate (MZ, se4);

%%

%defining GFP+ cells in IZ/SVZ
IZSVZ = ch1_c;
IZSVZ(1:Isvz, :) = 0;
IZSVZ = imdilate (IZSVZ, se4);
%%

%subplots of all 4 figures
hfig1 = figure(1); 
set(hfig1, 'Position', [123    18   570   687]); 
subplot_tight(2,2,1);
imshow (imadjust (Ich1_s, [0 0.05]));
title ('reporter postive cells'); 

subplot_tight(2,2,2);
imshow(ch1_total);
title ('all reporter+');

subplot_tight(2,2,3);
imshow(MZ);
title('MZ')

subplot_tight(2,2,4);
imshow(IZSVZ);
title('IZ-SVZ');

%calculating area
r = double (r);
c = double (c);
Imz = double (Imz); 
Ta = r*c/px/10^6; 
MZa= Imz*c/px/10^6;
SVZa= (r-Isvz)*c/px/10^6;

% object labeling MZ % extract features (Intensity this case) for each object 
[~, num_m] = bwlabel(MZ); 
[~, num_svz] = bwlabel(IZSVZ); 
[~, totalcell] = bwlabel(ch1_c); 
data_mz = {'MZ', num_m, MZa};
data_svz = {'SVZ', num_svz, SVZa}; 
data_total = {'All', totalcell, Ta}; 
data_all = vertcat(data_mz, data_svz, data_total); 

% consolidating data
data_all(:,4) = {Brainno};
data_all (:,5) = {imagename}; 
data_all (:,6) = {Age};
data_all (:,7) = {Genotype};
A = cell2dataset(data_all, 'VarNames', {'Region', 'GFP+', 'Areamm2', 'Brainno', 'Fieldname', 'Age', 'Genotype'});
export (A, 'File', s(1).outputfile,'Delimiter',','); 
end
