close all, clear all;
imagename = 'b4_E18_s027'; %define the root of your file here
Brainno = 'B1';
Age = 'E18.5';
s = struct('inputfile', {strcat(imagename, '.tif')}, 'outputfile', strcat(imagename, '_Mafb.csv')); 

% read images and background subtract 
Idapi = imread(s(1).inputfile, 2); %channel DAPI
bkg3 = mean (mode ((double(Idapi))));
Idapi_s = Idapi - bkg3; 

%calculating how big the Dapi image is
[r , c] = size (Idapi);
r = int16(r);
d = r/4;

%plot profile of dapi
avg = mean (Idapi_s');
avg_s= smooth (avg, 10, 'rlowess');
max_s = max (avg_s);
df_s = diff(avg_s);
df_s2 = smooth(df_s, 20, 'rlowess')*0.01* max_s;
avg_df= mean(df_s2);
figure, hold on
plot (df_s2, 'b');
hold off


[Cmax, Imz] = max(df_s2); %1st max = boundary to MZ
Imz = 1.2*Imz; 
[Cmin, Ivz] = min(df_s2); %1st min = boundry to IZ/ VZ
figure, imshow (imadjust (Idapi_s, [0 0.05]));

%checking ROI defination

MZ = Idapi;
MZ(Imz:r, :) = 0;
figure, imshow (imadjust (MZ, [0 0.05])), title ('MZ boundary');
CP = Idapi; 
CP (Ivz:r, :) = 0;
CP = CP - MZ; 
figure, imshow (imadjust (CP, [0 0.05])), title ('CP boundary');
SVZ = Idapi - MZ - CP; 
figure, imshow (imadjust (SVZ, [0 0.05])),title ('IZ/SVZ boundary');
