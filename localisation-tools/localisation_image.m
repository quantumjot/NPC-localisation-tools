%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binned image creation for Single-Molecule Localisation Microscopy
%
% Generates a simple binned image from localisation microscopy data
%
% INPUT:
%   molecules data from impy (columns 2,3 are X,Y)
%   bin size
%
% OUTPUT:
%   binDat structure containing:
%   - uncorrected, binned image
%   - hash table for lookup
%
% NOTE:
%   Constantly converting between (display) image space and coordinate 
%   space is a pain since the arrays are stored and displayed differently 
%   (i,j) vs (j,i).
%
% Lowe, A.R. 2010-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [varargout] = localisation_image(varargin)

binDat = [];
options = get_options();

% check the type of the input object
if isa(varargin{1},'char')
    path_to_file = varargin{1};
    try
        [data] = load(path_to_file);
        binDat.original_file = path_to_file;
        molecules = data.xy;
        
        % TODO: use other sources of data here
    catch
        disp('Cannot load specified file...');
    end
else
    molecules = varargin{1};
end


if (nargin > 1)
    options.bin_size = varargin{2};
end

if options.verbose
    fprintf('Calculating localisation image (bin size: %2.2f)...\n', options.bin_size);
end

% calculate the bins in a fast way
binx = 1+floor((1.0/options.bin_size)*(molecules(:,1)-0.5*options.bin_size));
biny = 1+floor((1.0/options.bin_size)*(molecules(:,2)-0.5*options.bin_size));
bins = [binx, biny];

% sanity check that we haven't exceeded bounds
min_bin = 1;
max_bin = max(bins);

bins_to_exclude = (bins(:,1)<min_bin | bins(:,2)<min_bin | bins(:,1)>max_bin(1) | bins(:,2)>max_bin(2));
bins(bins_to_exclude,:) = [];
molecules(bins_to_exclude,:) = [];

% now make the 2D histogram of the data, using the sparse matrix trick
[output_image] = uint16(full(sparse(bins(:,1), bins(:,2), 1)));


% if we've got an output, calculate a hash table for later use
do_display = 1;

if nargout > 0
    varargout{1} = output_image;
    do_display = 0;
end

if nargout > 1
    
    % use this opportunity to create a hash table for fast lookup later on
    fprintf('Calculating hash table...\n');
    hash_indices = sortrows([uint32(bins(:,1)+(bins(:,2)-1)*max_bin(2)), (1:size(bins,1))'],1);
    unique_hash_indices = unique(hash_indices(:,1));

    % make a table of ranges so that we don't have to do lots of lookups
    % while creating the hash table
    boundaries = [1];
    d = diff(hash_indices(:,1))>0;
    boundaries = cat(1, boundaries, find(d==1)+1);
    end_boundaries = cat(1,boundaries(2:end,1)-1,size(hash_indices,1));
    boundaries = [boundaries, end_boundaries];
    
    % now use a MATLAB container class to hold the hash table
    % basically a wrapper around a Java Map class
    keys = num2cell(unique_hash_indices);
    vals = cell(length(unique_hash_indices),1);
    for i = 1:length(unique_hash_indices)
        rows = boundaries(i,1):boundaries(i,2);
        vals{i} = hash_indices(rows,2);
    end
         
    hash_table = containers.Map(keys,vals);
    
    % output
    binDat.image = output_image;
    binDat.bin_size = options.bin_size;
    binDat.pixels_2_nm = options.pixels_2_nm;
    binDat.hashtable = hash_table;
    binDat.localisations = molecules;
    binDat.title = 'dSTORM image';
    binDat.date = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
    
    varargout{2} = binDat;
    
    % set up a future display of the image and the binDat structure
    do_display = 0;
end


%%
% display the image

if do_display

    % add a scale bar
    display_image = output_image;

    scale_bar_width = (options.scale_bar*1000.) / (options.pixels_2_nm * options.bin_size);
    scale_bar_height = (options.scale_bar*50.) / (options.pixels_2_nm * options.bin_size);
    scale_bar_offset = (10./options.bin_size); % 10 CCD pixels from edge

    scale_bar = [ size(output_image,1) - round(scale_bar_height) - scale_bar_offset,...
                  size(output_image,1) - scale_bar_offset,...
                  size(output_image,2) - round(scale_bar_width) - scale_bar_offset,...
                  size(output_image,2) - scale_bar_offset];

    display_image(display_image > options.display_max) = options.display_max;

    % set the pixel values of the image in the region of the scale bar
    display_image(scale_bar(1)-1:scale_bar(2)+1,scale_bar(3)-1:scale_bar(4)+1) = 0;
    display_image(scale_bar(1):scale_bar(2),scale_bar(3):scale_bar(4)) = options.display_max;


    figure
    imagesc(display_image,[0, options.display_max]);
    if exist('binDat.original_file')
        title(binDat.original_file, 'Interpreter','None');
    end
    colormap(hot);
    colorbar();
    axis image;
end



return