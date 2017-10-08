%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inspect_structure
%
% Plot details of an individual structure, including the original image,
% Localisations, and transformations applied.
%
% INPUT:
%   molecules data from impy in the form of a binDat
%
% OUTPUT:
%   None
%
% Lowe, A.R. 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = inspect_structure(binDat, idx)

opts = get_options;

% make some images of the structures
[raw_xy] = binDat.structures(idx).raw;
[rot_xy] = binDat.structures(idx).rotated;
[raw_image] = localisation_image(crop_localisations(raw_xy), opts.alignment_bin_size);
[rot_image] = localisation_image(crop_localisations(rot_xy), opts.alignment_bin_size);


figure
subplot(2,5,[1,2])
imagesc(binDat.image, [0 opts.display_max])
hold on
plot(binDat.structures(idx).centroid(1),binDat.structures(idx).centroid(2),'gs'); 
colormap(hot)
axis image
title('Total image')

env_hist = plot_envelope(binDat, binDat.poly_hist, opts);

subplot(2,5,[6,7])
hold on
plot(raw_xy(:,1), raw_xy(:,2), 'kx')
patch(  env_hist(:,1).*opts.bin_size, env_hist(:,2).*opts.bin_size, ...
        [.2 .2 .2],'EdgeColor',[0. 0. 0.],'FaceColor',[0.,0.,0.], ...
        'FaceAlpha',0.1); 
hold off
axis image

xlim([binDat.structures(idx).centroid(1)*opts.bin_size-max(opts.window_size), ...
    binDat.structures(idx).centroid(1)*opts.bin_size+max(opts.window_size)]); 

ylim([binDat.structures(idx).centroid(2)*opts.bin_size-max(opts.window_size), ...
    binDat.structures(idx).centroid(2)*opts.bin_size+max(opts.window_size)]); 


subplot(2,5,[3])
imagesc(raw_image)
axis image

subplot(2,5,[4])
imagesc(rot90(rot_image))
axis image

return

function [crop_xy] = crop_localisations(xy)
crop_xy = xy - repmat(min(xy),size(xy,1),1);
return