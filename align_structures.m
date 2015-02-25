%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align_structures
%
% Align structures identified from localisation microscopy using a
% reference object and rigid body transformations.  Should work with both
% image based methods (normalised cross-correlation) and point cloud based
% methods such as Iterative Closest Point (ICP).
%
% INPUT:
%   molecules data from impy in the form of a binDat
%   reference structure
%
% OUTPUT:
%   transformation matrices added to binDat
%
% NOTE:
%   Creates an affine transformation matrix (rigid body, i.e. rotation and
%   translation only) which can be applied both to 'images' and to the
%   original point clouds for alignment.
%
%   Can use Wilm ICP algorithm:
%   http://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point
%
%
%   Transformation matrices are of the following format:
%
%   Rotation =      |   cos(theta) -sin(theta)  0 |
%                   |   sin(theta)  cos(theta)  0 |
%                   |   0           0           1 |
%
%
%   Translation =   | t_x |
%                   | t_y |
%                   |  0  |
%
%
%   Also note that we could combine both of these into a single affine
%   transformation matrix, but I prefer to keep them separate for
%   interface uniformity with other algorithms.
%
%
% Lowe, A.R. 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [binDat] = align_structures(binDat, ref_structure)

options = get_options();

num_structures = length(binDat.structures);

if options.use_rotations
    % set up the rotations to be +/- a small angular range
    [rotations] = linspace(-options.alignment_angle_range,options.alignment_angle_range,21);
else
    [rotations] = 0;
end

% go through each structure and compare it with the reference. Determine
% the rigid body transformation required to align it

for i = 1:num_structures

    % get the pre-rotated localisations
    [localisations] = binDat.structures(i).rotated;
    
    % since the localisations are centred at (0,0), return to a positive
    % coordinate space
    localisations(:,1) = localisations(:,1) + options.window_size(2);
    localisations(:,2) = localisations(:,2) + options.window_size(1);
    
    if strcmp(options.alignment_method,'ICP')
        % use the ICP algorithm (here, corr is the RMS error)
        [rot, trans, corr] = icp(ref_structure, cmp, 15);
        
    elseif strcmp(options.alignment_method, 'xcorr')
        % use normalised cross-correlation for the alignment
        [cmp] = uint16(zeros(   ceil(2.*options.window_size(1)/options.alignment_bin_size),...
                                ceil(2.*options.window_size(2)/options.alignment_bin_size)  ));
        
        % do some horrible cropping to make sure everything is correct size
        [cmp_raw] = rot90(localisation_image(localisations, options.alignment_bin_size));
        cmp(1:size(cmp_raw,1),1:size(cmp_raw,2)) = cmp(1:size(cmp_raw,1),1:size(cmp_raw,2)) + imcrop(cmp_raw,[1,1,size(cmp,2),size(cmp,1)]);
                     
        % use the normalised cross correlation method
        [rot, trans, corr] = align_image(ref_structure, cmp, rotations, options.alignment_bin_size);
    end
    
    
    
    % store those values for later
    binDat.structures(i).alignment_method = options.alignment_method;
    binDat.structures(i).alignment_rotation = rot;
    binDat.structures(i).alignment_translation = trans;
    binDat.structures(i).alignment_error = 1.-corr;     % note, the type of ERROR is different for each of the algorithms
        
end



return




function [rot, trans, corr] = align_image(cmp_structure, ref_structure, rotations, bin_size)

% start correlation poor
best_correlation = -1e99;
theta = 0;

for r=1:length(rotations)

    % make a rotated version if necessary
    [cmp] = imrotate(cmp_structure, rotations(r), 'bilinear', 'crop');

    
    % normalised cross correlation
    [c] = normxcorr2(ref_structure, cmp);
    [max_c, imax] = max(abs(c(:)));

    % if we're doing better, store the info
    if max_c > best_correlation
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        offset = [ (xpeak-size(ref_structure,2)) (ypeak-size(ref_structure,1)) ];
        best_correlation = max_c;
        theta = rotations(r);
    end

end

% now make the transforms 
trans = [ -offset(1)*bin_size; -offset(2)*bin_size;];

rot =   [ cosd(theta) -sind(theta)
          sind(theta)  cosd(theta) ];
      
corr = best_correlation;

return

