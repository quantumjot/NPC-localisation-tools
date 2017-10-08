%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_montage
%
% Plot the localisation microscopy image with the user defined nuclear
% envelope, envelope histogram and putative pore structures overlaid
%
% INPUT:
%   molecules data from impy in the form of a binDat
%   envelope vector
%   options
%
% OUTPUT:
%   None
%
% Lowe, A.R. 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_montage(binDat)

% go through each set of localisations, rotate and translate them and make
% an image for the montage

all_xy = [];
for i=1:length(binDat.structures)

    [xy] = binDat.structures(i).rotated;
    [rot] = binDat.structures(i).alignment_rotation;
    [trans] = binDat.structures(i).alignment_translation;
    
    [rotated_xy] = (rot*xy')';
    [translated_xy] = rotated_xy;
    translated_xy(:,1) = translated_xy(:,1)-trans(1);  
    translated_xy(:,2) = translated_xy(:,2)-trans(2);
    
    if i<6
        [im] = localisation_image(subtract_min2d(translated_xy), 0.01);
        imwrite(rot90(uint8(im)), hot, strcat('NPC_',num2str(i),'.tif'));
    end
    % store it
    all_xy = cat(1,all_xy, translated_xy);
end

% all_xy(:,1) = all_xy(:,1) - min(all_xy(:,1));
% all_xy(:,2) = all_xy(:,2) - min(all_xy(:,2));

all_xy = subtract_min2d(all_xy);

% make a localisation image
[im] = localisation_image(all_xy, 0.01);

% make the projections
y_bins = linspace(0,max(all_xy(:,2)), size(im,2));
x_bins = linspace(0,max(all_xy(:,1)), size(im,1));
[y_projection, ~] = hist(all_xy(:,2), y_bins);
[x_projection, ~] = hist(all_xy(:,1), x_bins);

figure
h1 = subplot(4,4,[1,2,3,5,6,7,9,10,11]);
imagesc(rot90(im))
colormap hot
axis image
h2 = subplot(4,4,[4,8,12]);
plot(y_projection,y_bins,'k-');
h3 = subplot(4,4,[13,14,15]);
plot(x_bins,x_projection,'k-');

return

function [xy] = subtract_min2d(xy)
min_xy = min(xy);
xy = xy - repmat(min_xy, size(xy,1), 1);
return
