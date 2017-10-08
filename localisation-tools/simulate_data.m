%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_envelope
%
% Generates a simulated localisation data set for oriented structures on a
% circular envelope
%
% INPUT:
%   Noise level (std of gaussian noise added to localisations)
%
% OUTPUT:
%   Simulated localisation data
%
% NOTE:
%   We add normally distributed noise to generate the localisation data
%
% Lowe, A.R. 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [simulated_data, ref_image] = simulate_data(noise, num_mols)

% get some options
options = get_options();

% make a simulated arrow
[arrow_xy] = simulated_arrow(num_mols, 1);

% now make a circular envelope with arrows rotated and pointed toward the centre
nucleus_radius = 100.;
num_pores = 200;
CCD_size = 256.;

simulated_data = [CCD_size, CCD_size]; % add one localisation at the extreme for visual purposes

for p = 1:num_pores;
    
    theta = ((2.*pi*p) / num_pores);
    
    centre = [  CCD_size/2. + nucleus_radius*cos(theta),... 
                CCD_size/2. + nucleus_radius*sin(theta) ];
    
    % get the simulated localisations
    [new_arrow] = arrow_xy;
    
    % make sure the arrow is normal to the envelope and pointed inward
    theta = theta - pi/2.;
    
    % rotate
    [rotation_matrix] = [ cos(theta) -sin(theta); 
                          sin(theta)  cos(theta) ];
                        
    [new_arrow] = (rotation_matrix*new_arrow')';
    
    % add some normally distributed noise
    [new_arrow] = new_arrow + normrnd(0., noise, [size(new_arrow,1), size(new_arrow,2)]);
    
    % translate
    new_arrow(:,1) = new_arrow(:,1) + centre(1);
    new_arrow(:,2) = new_arrow(:,2) + centre(2);
    
    simulated_data = cat(1,simulated_data, new_arrow);
end

% now make a reference image for alignment
[ref_arrow] = arrow_xy;
ref_arrow(:,1) = ref_arrow(:,1) + .5;
ref_arrow(:,2) = ref_arrow(:,2) + 1.;
[ref_image] = localisation_image(ref_arrow, options.alignment_bin_size);
                   
return

% generate coordinates for a simulated arrow
function [xy] = simulated_arrow(num_localisations, size)

num_localisations_body = ceil(num_localisations/2);
num_localisations_arms = ceil(num_localisations/4);

% make the body of the arrow
y = linspace(-size,size, num_localisations_body)';
x = zeros(length(y),1);
xy = [x,y];

% now make the arms
dy = size / num_localisations_arms;
dx = size / (2*num_localisations_arms);
for i = 1:num_localisations_arms
    xy = cat(1,xy, [-i*dx, -size+i*dy]);
    xy = cat(1,xy, [i*dx, -size+i*dy]);
end



return
