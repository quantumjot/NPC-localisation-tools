%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopick_pores
%
% Function to isolate putative NPC structures from the nuclear envelope of
% cells in localisation microscopy data
%
% INPUT:
%   molecules data from impy in the form of a binDat
%   envelope vector (optional)
%
% OUTPUT:
%   envelope_X, envelope_Y - the original user input
%   pores - positions of structures in the envelope and vectors for align
%
% NOTE:
%   At some point in the future it would be nice to extend the envelope
%   interpolation to utilise a spline defined by the surface normals. At
%   the moment this is a linear interpolation between the vertices of the
%   user input.
%
% Lowe, A.R. 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [binDat] = autopick_pores(binDat, varargin)

% grab the data and options
options = get_options();


% if we're given a predifined envelope vector, use it, else ask the user
% for input
if (nargin > 2)
    envelope_X = varargin{1};
    envelope_Y = varargin{2};
else
    % get the user to define the nuclear envelope
    disp('-----');
    disp('Please select the outline of the nuclear envelope, using the');
    disp('mouse. Note that the outline should be desribed in a clockwise');
    disp('manner to preserve the orientation of the inside/outside vector');
    disp('-----');
    figure, imagesc(binDat.image,[0 options.display_max]), axis image, colormap(gray)
   
    h = warndlg('Please pick the nuclear rim outline...');
    waitfor(h);
    
    [~, envelope_Y, envelope_X] = roipoly();
end


% refine the envelope by calculating the histogram along the envelope
exempt = []; % can exempt sections of envelope here
[poly_hist, pore_locations] = polygon_histogram(binDat, [envelope_X, envelope_Y], exempt);

% now go back and extract the pores as necessary
fprintf('Found %d candidate structures. Extracting... \n',length(pore_locations));

pores(length(pore_locations)) = struct();
for i = 1:length(pore_locations)
    candidate = poly_hist(pore_locations(i),:);

    pores(i).centroid = candidate(1:2);
    pores(i).vector = candidate(3:4);
    
    % now get the actual data
    [counts, raw, rotated] = get_localisations_in_ROI(  binDat,...
                                                        pores(i).centroid,...
                                                        pores(i).vector,...
                                                        options.window_size(1)/binDat.bin_size,...
                                                        options.window_size(2)/binDat.bin_size );
    pores(i).counts = counts;
    pores(i).raw = raw;
    pores(i).rotated = rotated;
    pores(i).window_size = options.window_size;
end

binDat.envelope = [envelope_X,envelope_Y];
binDat.structures = pores;
binDat.poly_hist = poly_hist;

% plot some output
plot_envelope(binDat, poly_hist, options)

return






% calculate the center of mass of the polygon
function [cmass] = center_of_mass(polygon)

area = polyarea(polygon(:,1),polygon(:,2));
cx=0;
cy=0;
for i = 1:length(polygon)
    if (i < length(polygon))
        xy_a = [polygon(i,1)   polygon(i,2)  ];
        xy_b = [polygon(i+1,1) polygon(i+1,2)];
    else
        xy_a = [polygon(i,1) polygon(i,2)];
        xy_b = [polygon(1,1) polygon(1,2)];
    end
    cx = cx + (xy_a(1)+xy_b(1))*((xy_a(1)*xy_b(2))-(xy_b(1)*xy_a(2)));
    cy = cy + (xy_a(2)+xy_b(2))*((xy_a(1)*xy_b(2))-(xy_b(1)*xy_a(2)));
end

cmass = [(1 / (6*area))*cx, (1 / (6*area))*cy];
return






% pass a scanning window around a polygon and calculate number of localisations
% within that window, then parse it for potential locations
function [poly_hist, pore_locations] = polygon_histogram(binDat, polygon, exempt)

poly_hist = [];
options = get_options();

% calculate the normals at each vertex of the envelope
[normals] = calculate_surface_normals(polygon);

new_envelope = [];
disp('Generating envelope histogram...');

for i = 1:length(polygon)-1
    
    % exempt sections of the envelope if necessary
    if any(exempt==i)
        continue
    end
    
    % make sure the vector is a closed loop
    if (i < length(polygon))
        xy_a = [polygon(i,1)   polygon(i,2)  ];
        xy_b = [polygon(i+1,1) polygon(i+1,2)];
    else
        xy_a = [polygon(i,1) polygon(i,2)];
        xy_b = [polygon(1,1) polygon(1,2)];
    end
    
    % interpolate along the line
    edge_length = 1. * sqrt((xy_a(1) - xy_b(1))^2 + (xy_a(2) - xy_b(2))^2);
    dx = (xy_b(1) - xy_a(1)) / edge_length;
    dy = (xy_b(2) - xy_a(2)) / edge_length;
    
    % start at the first vertex
    x = xy_a(1);
    y = xy_a(2);
    
    % normal of the first vertex of the edge
    nx = normals(i,1);
    ny = normals(i,2);
    
    % normal of the second vertex of the edge
    if (i < length(polygon)-1)
        nx_b = normals(i+1,:);
    else
        nx_b = normals(1,:);
    end
    
    % calculate the rate of change of the normal
    dnx = (nx_b(1) - nx) / edge_length;
    dny = (nx_b(2) - ny) / edge_length;
       
    % interpolate along the envelope segment
    d = 0;
    while d<=edge_length;

        % TODO: make window width/length image bin size independent
        
        % count the number of localisations in the window
        [counts] = get_localisations_in_ROI(binDat, [x,y], [nx,ny], 10., 1.);
                
        % store all of the data here:
        % x,y, normal_x, normal_y, number_of_localisations, polygon_segment
        poly_hist = cat(1,poly_hist, [x,y,nx,ny,counts,i]);

        % move the centroid of the window    
        x=x+dx;
        y=y+dy;
        
        % interpolate the surface normal vector
        nx=nx+dnx;
        ny=ny+dny;
        
        % increment the segment counter
        d = d + 1;
    end
    
% give the user some feedback as to the progress   
disp(sprintf(' - Completed %d of %d segments of envelope...',i,length(polygon)-1));
            
end


% now perform a savitzky-golay filtering of the envelope
% new_envelope(:,5:6) = sgolayfilt(new_envelope(:,5:6),1,21);

% now find peaks within the histogram (these are the putative pore
% positions)

% minimum separation between peaks in histogram (based on image bin size)
min_separation_pixels = ceil((options.min_separation/binDat.pixels_2_nm)/binDat.bin_size);
[pks, locs] = findpeaks(poly_hist(:,5),'MINPEAKHEIGHT',options.min_localisations,'MINPEAKDISTANCE',min_separation_pixels);



% finally, return the locations of the putative pores
pore_locations = locs;

return
