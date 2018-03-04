%%
MZperim = bwperim(MZ_gfp); 
allperim = bwperim(bwf_total); 
cmap = colormap ('gray(200)');
res = grs2rgb(imadjust(I_gfp),cmap);
MZbr = imoverlay(res, MZperim, [1 0 0]);
%figure, imshow (MZbr)
allbr = imoverlay(res, allperim, [0 1 0]);

hfig1 = figure(2); 
set(hfig1, 'Position', [400 200 900 800]); 

%first subplot
subplot_tight(1,2,1); imshow(allbr); title ('All reporter'); 

%second subplot
subplot_tight(1,2,2); imshow (MZbr); title ('MZ+ cells'); 
