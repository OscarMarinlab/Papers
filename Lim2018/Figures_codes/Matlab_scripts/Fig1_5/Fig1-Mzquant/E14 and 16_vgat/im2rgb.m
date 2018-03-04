  function rgb = im2rgb(im,full_map); %nested
    %coerce intensities into gray range [0,1]
    gray = imadjust(im,[],[0 1]);
    %generate indexed image
    num_colors = size(full_map,1);
    ind = gray2ind(gray,num_colors);
    %convert indexed image to RGB
    rgb = ind2rgb(ind,full_map);
  end %im2rgb