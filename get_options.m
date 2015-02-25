function [options] = get_options

% general features of the system
options.bin_size = 0.2;                 % fraction of a CCD pixel for rendering
options.pixels_2_nm = 100.;             % conversion between CCD pixels and nm
options.scale_bar = 5;                  % scale bar in microns
options.display_max = 50;               % maximum pixel intensity for display purposes
options.wavelength = 670.;              % wavelength of fluorophore emission
options.numerical_aperture = 1.49;      % NA of objective lens
options.verbose = 0;                    % set this option for verbose output

% options for NPC finding
options.scanning_window = [20.,1.];     % scanning window length (nm)
options.min_localisations = 100.;       % minimum number of localisations to be considered
options.min_separation = 50.;           % separation between peaks in envelope histogram (nm)

% options for structure alignment
options.window_size = [2., .6];         % crop window size in CCD pixels
options.alignment_bin_size = 0.05;      % alignment bin size
options.alignment_method = 'xcorr';     % use image alignment as a default
options.use_rotations = 0;              % allow angular adjustment during alignment
options.alignment_angle_range = 8.;     % range of angles to vary

return