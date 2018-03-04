%MafB analyses

close all, clear all;
imagename = 'b2_E18_s008'; %define the root of your file here
Brainno = 'B2';
Age = 'E18.5';
%name the files you are interested in analyzing
s = struct('inputfile', {strcat(imagename, '.tif')}, 'outputfile', strcat(imagename, '_sstMafb.csv')); 
se1 = strel('line', 10, 0);
se1b = strel('line', 10, 40);
se2 = strel('disk', 1, 8);
se5 = strel(ones(1,1));
se3 = strel('disk', 2, 8);
se4 = strel ('disk',4);

% read images and background subtract 
Ich1 = imread(s(1).inputfile, 3); %channel (tdt) reporter staining
Ich2 = imread(s(1).inputfile, 1); %channel Mafb (green)
Idapi = imread(s(1).inputfile, 2); %channel DAPI
bkg1 = mean(mode ((double(Ich1))));
bkg2 = mean(mode ((double(Ich2))));
bkg3 = mean (mode ((double(Idapi))));
Ich1_s = Ich1 - bkg1;
Ich2_s = Ich2 - bkg2; 
Idapi_s = Idapi - bkg3; 
Ich2_a = imadjust (Ich2_s, [0 0.01]); 

%ch1 image segmentation, reporter channel

Ich1_ed = imerode (Ich1_s, se3); 
Ich1_ed = imdilate (Ich1_ed, se3); 
Ich1_ed_d = imreconstruct (Ich1_ed, Ich1_s);
c1_fgm = imregionalmax(Ich1_ed_d); 
c1_fgm = imfill(c1_fgm, 'holes'); 
c1_fgm2 = bwareaopen (c1_fgm, 800);
c1_fgm = c1_fgm - c1_fgm2; 
Level1 = graythresh(Ich1_s); 
bw_ch1 = im2bw(Ich1_s, Level1); 
bw_ch1 = bwareaopen (bw_ch1, 10); 
ch1_total = bw_ch1 & c1_fgm;
ch1_total = imdilate(ch1_total, se2); 
Ich1_a = (imadjust (Ich1_s, [0 0.01])); 
figure, imshow (Ich1_a), title ('reporter postive cells'); 
figure, imshow(ch1_total), title ('all reporter+ ');


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

%calculating how big the Dapi image is
[r , c] = size (Idapi);
r = int16(r);

%plot profile of dapi
avg = mean (Idapi_s');
avg_s= smooth (avg, 2, 'rlowess');
max_s = max (avg_s);
df_s = diff(avg_s);
df_s2 = smooth(df_s, 5, 'rlowess')*0.1* max_s;
avg_df= mean(df_s2);
figure, hold on
plot(avg_s, 'r');
plot (df_s2, 'g');
hold off


[Cmax, Imz] = max(df_s2); %1st max = boundary to MZ
Imz = round(1.2*Imz); 
[Cmin, Ivz] = min(df_s2); %1st min = boundry to IZ/ VZ

%defining MZ
MZ = ch1_cd;
MZ(Imz:r, :) = 0;
figure, imshow (MZ), title ('MZ boundary');
CP = ch1_cd; 
CP (Ivz:r, :) = 0;
CP = CP - MZ; 
figure, imshow (CP), title ('CP boundary');
SVZ = ch1_cd - MZ - CP; 
figure, imshow (SVZ),title ('IZ/SVZ boundary');

%visualizing double positive
MZperim = bwperim(MZ);
CPperim = bwperim(CP);
SVZperim = bwperim(SVZ); 
overlay_ch2 = imoverlay(Ich2_a, MZperim, [1 0 0]);
overlay_ch2 = imoverlay(overlay_ch2, CPperim, [0 1 0]);
overlay_ch2 = imoverlay(overlay_ch2, SVZperim, [0 0 1]);
figure, imshow(overlay_ch2, []), title ('tdt positive area in Maf');

Ich1_a = imadjust(Ich1, [0 0.04]);
overlay_ch1 = imoverlay(Ich1_a, MZperim, [1 0 0]);
overlay_ch1 = imoverlay(overlay_ch1, SVZperim, [0 0 1]);
overlay_ch1 = imoverlay(overlay_ch1, CPperim, [0 1 0]);
figure, imshow(overlay_ch1, []), title ('tdt positive area in tdt');


% object labeling MZ % extract features (Intensity this case) for each object 
[M1, num_m] = bwlabel(MZ); 
stats = regionprops(M1, 'PixelIdxList');  
for i=1:num_m 
    Imz_int_ch1(i)= mean(Ich1(stats(i).PixelIdxList));
    Imz_int_ch2(i)= mean(Ich2(stats(i).PixelIdxList)); 
end


% CP
[CP1, num_cp] = bwlabel(CP); 
stats = regionprops(CP1, 'PixelIdxList');  
for i=1:num_cp 
    Icp_int_ch1(i)= mean(Ich1(stats(i).PixelIdxList));
    Icp_int_ch2(i)= mean(Ich2(stats(i).PixelIdxList)); 
end

% svz
[svz1, num_svz] = bwlabel(SVZ); 
stats = regionprops(svz1, 'PixelIdxList');  
for i=1:num_svz 
    Isvz_int_ch1(i)= mean(Ich1(stats(i).PixelIdxList));
    Isvz_int_ch2(i)= mean(Ich2(stats(i).PixelIdxList)); 
end

% labeling MZ data
label_mz(1:num_m,1) = {'MZ'};
data_mz1 = Imz_int_ch1';
data_mz1 = num2cell(data_mz1);
label_mz1(1:num_m, 1) = {'tdt_int'}; 
data_mz2 = Imz_int_ch2';
data_mz2 = num2cell(data_mz2);
label_mz2(1:num_m, 1) = {'Mafb_int'};  
data_mz1 = horzcat(data_mz1, label_mz1, label_mz); 
data_mz2 = horzcat(data_mz2, label_mz2, label_mz); 

% labeling CP data
label_cp(1:num_cp,1) = {'CP'};
data_cp1 = Icp_int_ch1';
data_cp1 = num2cell(data_cp1);
label_cp1(1:num_cp, 1) = {'tdt_int'}; 
data_cp2 = Icp_int_ch2';
data_cp2 = num2cell(data_cp2);
label_cp2(1:num_cp, 1) = {'Mafb_int'};  
data_cp1 = horzcat(data_cp1, label_cp1, label_cp); 
data_cp2 = horzcat(data_cp2, label_cp2, label_cp); 

% labeling SVZ data
label_svz(1:num_svz,1) = {'SVZ'};
data_svz1 = Isvz_int_ch1';
data_svz1 = num2cell(data_svz1);
label_svz1(1:num_svz, 1) = {'tdt_int'}; 
data_svz2 = Isvz_int_ch2';
data_svz2 = num2cell(data_svz2);
label_svz2(1:num_svz, 1) = {'Mafb_int'};  
data_svz1 = horzcat(data_svz1, label_svz1, label_svz); 
data_svz2 = horzcat(data_svz2, label_svz2, label_svz);

%combine all area data & experorting
data_all = vertcat(data_mz1, data_mz2, data_cp1, data_cp2, data_svz1, data_svz2); 
data_all(:,4) = {Brainno};
data_all (:,5) = {imagename}; 
data_all (:, 6) = {Age};
A = cell2dataset(data_all, 'VarNames', {'Intensity', 'Channel', 'Region', 'Brainno', 'Fieldname', 'Age'});
export (A, 'File', s(1).outputfile,'Delimiter',','); 
