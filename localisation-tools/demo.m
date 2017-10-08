%% demo.m

% get some simulated localisation data
[xy, ref_image] = simulate_data(.01, 1000);

% pad the ref image to be the correct size for image alignment
[ref_imageX] = rot90(padarray(ref_image, [2,20]));

% make an image structure and display a downsampled image
localisation_image(xy);
[im, binDat] = localisation_image(xy);
diffraction_limited_image(binDat);

% make a perfect user input?
x = [];
y = [];
for i=1:100
    t = (2.*pi) * ((100-i)/99.);
    x = cat(1,x, 5.*(128.+100.*cos(t)));
    y = cat(1,y, 5.*(128.+100.*sin(t)));
end
[binDat] = autopick_pores(binDat, x, y);

% pick the structures using a manual starting input
% [binDat] = autopick_pores(binDat);

% align them
[binDat] = align_structures(binDat, ref_imageX);

% output a montage
plot_montage(binDat);