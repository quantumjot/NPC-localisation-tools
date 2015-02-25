function [varargout] = diffraction_limited_image(binDat)

options = get_options();

% convolve the image with a 'PSF'
PSF_FWHM = (.61*options.wavelength)/options.numerical_aperture;
PSF_sigma = PSF_FWHM / (2.355 * options.pixels_2_nm * options.bin_size);
PSF = fspecial('gaussian', ceil(9.*PSF_sigma), PSF_sigma);
convolved_image = conv2(double(binDat.image), PSF);

% downsample back to the camera pixel size
diffraction_limited = imresize(convolved_image, binDat.bin_size);

if nargout > 0
    varargout{1} = diffraction_limited;
else
    figure
    imagesc(diffraction_limited)
    colormap gray
    axis image
    title('Simulated diffraction-limited image');
end

end