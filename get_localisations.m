%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract localisations from a pixel of the binned image
%
% INPUT:
%   molecules binDat
%   x,y coordinate of pixel in binDat binned image space
%
% OUTPUT:
%   raw localisations using hash table
%
% Lowe, A.R. 2010-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [localisations] = get_localisations(binDat, x, y)

% use the hash to return the localisations data
hash = uint32(x+(y-1)*size(binDat.image,2));


% check that the key exists
if isKey(binDat.hashtable, hash)
    rows = binDat.hashtable(hash);
    [localisations] = binDat.localisations(rows,:);
    return
else
    % disp('Warning... no localisations found in bin');
    localisations = [];
end

return