%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_envelope
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
% NOTE:
%   At some point in the future it would be nice to extend the envelope
%   interpolation to utilise a spline defined by the surface normals. At
%   the moment this is a linear interpolation between the vertices of the
%   user input.
%
% Lowe, A.R. 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = plot_envelope(binDat, poly_hist, options)

if nargout > 0
    display_envelope = 0;
else
    display_envelope = 1;
end





if display_envelope
    % first plot the STORM image and the user input envelope vector
    figure
    imagesc(binDat.image,[0 options.display_max]), axis image, colormap(hot);
    hold on
    for i = 1:length(binDat.envelope)-1
        plot(binDat.envelope([i,i+1],2),binDat.envelope([i,i+1],1),'bo-');
    end
end


% now plot the envelope histogram as a symmetrical histogram rotated about
% the normal vector and aligned with the envelope
max_counts = max(poly_hist(:,5))/25.;
inner = [];
outer = [];
for i = 1:length(poly_hist)
    dx = (poly_hist(i,4)*poly_hist(i,5)/max_counts);
    dy = (poly_hist(i,3)*poly_hist(i,5)/max_counts);
    inner = cat(1,inner,[poly_hist(i,2)-dx,poly_hist(i,1)-dy]);
    outer = cat(1,outer,[poly_hist(i,2)+dx,poly_hist(i,1)+dy]);
    % plot a vector arrow denoting the direction of the surface normal
    if (mod(i,20) == 0 && display_envelope)
        edge_length = sqrt(dx*dx + dy*dy);
        quiver(poly_hist(i,2),poly_hist(i,1),20.*dx./edge_length,20.*dy./edge_length,'w');
    end
end

% plot the actual envelope histogram as a patch
xp = cat(1,inner(:,1),flipud(outer(:,1)));
yp = cat(1,inner(:,2),flipud(outer(:,2)));

if display_envelope == 0
    varargout{1} = [xp, yp];
end

if display_envelope
    patch(xp,yp,[.2 .2 .2],'EdgeColor',[1. 1. 1.],'FaceColor',[1.,1.,1.],'FaceAlpha',0.1); 
    % plot the location of found pores
    for i=1:length(binDat.structures)
        [xy] = binDat.structures(i).centroid;
        plot(xy(2), xy(1), 'g.', 'MarkerSize',8)
    end
    title('Interpolated NE with found NPCs')
    hold off
end

return