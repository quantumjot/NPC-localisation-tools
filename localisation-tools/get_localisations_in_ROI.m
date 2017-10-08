%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract localisations from an ROI of the binned image
%
% Uses the hash table for fast look-up. Optionally, will rotate the
% extracted ROI based on the applied vector
%
% INPUT:
%   molecules binDat
%   ROI centroid
%   Vector desribing ROI orientation
%   Width and 
%   Height of the box/ROI
%
%
% OUTPUT:
%   counts
%   raw localisations in ROI
%   rotated and centred normalisations
%
% Lowe, A.R. 2010-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = get_localisations_in_ROI(binDat, centroid, vector, width, height)
 
    % get the normal vector and normalise it to unit
    n_length = sqrt(vector(1)^2+vector(2)^2);
    n1 = [-vector(2),vector(1)]./n_length;
    n2 = [vector(1),vector(2)]./n_length;
        
    % calculate the vertices of polygon
    Ax = centroid(1) - width*n2(1)-height*n1(1);
    Ay = centroid(2) - width*n2(2)-height*n1(2);
    Bx = centroid(1) - width*n2(1)+height*n1(1);
    By = centroid(2) - width*n2(2)+height*n1(2);
    Cx = centroid(1) + width*n2(1)+height*n1(1);
    Cy = centroid(2) + width*n2(2)+height*n1(2);
    Dx = centroid(1) + width*n2(1)-height*n1(1);
    Dy = centroid(2) + width*n2(2)-height*n1(2);
    
    window_x = [Ax Bx Cx Dx Ax];
    window_y = [Ay By Cy Dy Ay];

    % now crop out the localisations
    [crop_region] = round([ min(window_x)-1, max(window_x)+1,...
                            min(window_y)-1, max(window_y)+1 ]);
                        
    crop_localisations = [];
    for y = crop_region(3):crop_region(4)
        for x = crop_region(1):crop_region(2)
            crop_localisations = cat(1,crop_localisations, get_localisations(binDat, x,y));
        end
    end
    
   
    if isempty(crop_localisations)
        varargout{1} = 0;
        return
    end
    
    % pick out points from original dataset
    [in] = inpolygon(   crop_localisations(:,1),...
                        crop_localisations(:,2),...
                        window_x.*binDat.bin_size,...
                        window_y.*binDat.bin_size );
    
    if nargout > 0
        % number of localisations in window
        varargout{1} = sum(in); 
    end
    
    if nargout > 1
        % the actual localisations data
        localisations_in_ROI = crop_localisations(in,:);
        varargout{2} = localisations_in_ROI;
    end
    
    if nargout > 2
        % if we want a rotated version, do it here
        rotated_ROI = localisations_in_ROI;
        
        % calcuate the angle relative to the origin
        origin = [0.0 1.0];
        theta = acos(dot(origin,n1));
        if n1(1)<0
            theta = -theta;
        end
        theta = theta + pi/2;
        
        % centre the localisations
        rotated_ROI(:,1) = rotated_ROI(:,1) - centroid(1)*binDat.bin_size;
        rotated_ROI(:,2) = rotated_ROI(:,2) - centroid(2)*binDat.bin_size;
        
        
        % rotate the localisations using this matrix
        rotation_matrix = [ cos(theta) -sin(theta); 
                            sin(theta)  cos(theta)   ];
                        
        rotated_ROI = rotation_matrix*rotated_ROI';
        
        % the rotated version of the data
        varargout{3} = rotated_ROI';
    end

return