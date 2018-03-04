%subplots scripts
%get (hfig1, 'Position')
hfig1 = figure(1); 
set(hfig1, 'Position', [1145 275 841 755]); 
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
imshow(SVZ);
title('SVZ')
