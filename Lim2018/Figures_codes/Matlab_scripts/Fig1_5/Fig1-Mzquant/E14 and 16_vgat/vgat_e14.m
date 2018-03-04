% May 2016- for vgat-cre X RCE - using confocal to optical section. this
% reduces the number of cells "merged" together thus allowing for more accurate
% quantification.


function data = vgat_e14(filename, A)


%%
I_gfp = imread (filename, A.ch);
bkg1 = mean (mode ((double(I_gfp))));
I_gfp = I_gfp-bkg1; %background subtraction of gfp
% opening-closing by reconstruction of GFP image
se2 = strel('disk', 1); %structral element 
se2B = strel ('disk', 1);
I_gfp_e = imerode(I_gfp, se2); %eroding center of mass of cell
%figure, imshow(I_gfp_e, []), title('Image eroded');
I_gfp_obr = imreconstruct(I_gfp_e, I_gfp); %reconstruction after eroding
%figure, imshow(I_gfp_obr, []), title('GFPImage reconstructed');
I_gfp_obrd = imdilate(I_gfp_obr, se2); %measuring local max
%figure, imshow(I_gfp_obrd, []), title('I_gfp_obrd'); %showing local max
I_gfp_obrcbr = imreconstruct(imcomplement(I_gfp_obrd), imcomplement(I_gfp_obr)); 
%reconstruct image with eroded and dilate image
I_gfp_obrcbr = imcomplement(I_gfp_obrcbr);
%inverse of eroded and dilated image, orginal image eroded & dilate
%figure, imshow(I_gfp_obrcbr, []), title('I_gfp_obrcbr')
%%
%gating GFP image by locate regional maxima
fgm = imregionalmax(I_gfp_obrd); %don't change
%figure, imshow(fgm, []), title ('fgm');
fgm2 = bwareaopen(fgm,5); %gating for size
figure, imshow(fgm2);
%%
[~, bw_gfp] = thresh_tool(imadjust(I_gfp)); %image to black and white - binary operation
bw_total = bw_gfp & fgm2;
bw_total = bwareaopen (bw_total, 5); 
%figure, imshow(bw_total); 
bw_gfp = bwareaopen(bw_gfp,10); % morphological opening to gate for the miniminal size of objects 
bwf_total = fgm2 & bw_gfp; % logical "AND" to only consider objects that overlap between fgm6 and bw3

%%
[r , col] = size (I_gfp); 
IndMZ = round(0.1*r); 
MZ_gfp= bwf_total; 
MZ_gfp(IndMZ:r, :) = 0;
svz_gfp= bwf_total; 
svz_gfp(1:IndMZ, :) = 0;
r = double (r);
col = double (col); 
Ta = r*col/(A.px^2)/10^6;
A_MZ = 0.1*Ta; 
A_SVZ = 0.90*Ta; 
%%
% object labeling 
[~, MZ] = bwlabel(MZ_gfp); 
[~, svz] = bwlabel(svz_gfp); 
%%
%subplot of figures
%subplots
hfig1 = figure(1); 
set(hfig1, 'Position', [400 200 1200 600]); 

%first subplot
subplot_tight(1,3,1); imshow(imadjust(I_gfp)); title ('All reporter'); 

%second subplot
subplot_tight(1,3,2); imshow (MZ_gfp); title ('MZ+ cells'); 

%3rd subplot 
subplot_tight(1,3,3); imshow(svz_gfp); title ('SVZ+ cells'); 
%

% write the data in excel spreasheet
%%
d(:, 1)= [MZ, svz];
d (:, 2)  = [A_MZ , A_SVZ]; 
data1= num2cell(d);
%%
[~,imagename, ~] = fileparts(filename);
data3 =  cell(2, 4); 
data3(:,4) = {imagename};
data3(:, 2) = {A.Brainno}; 
data3(:, 3) = {A.Age}; 
data3(:, 1)  = {'MZ', 'SVZ'}';
%%
C1 = horzcat (data1, data3); 
data = cell2dataset(C1, 'VarNames', {'GFP', 'Areamm2', 'Region', 'Brainno', 'Age', 'imagename'}); 
end



 
