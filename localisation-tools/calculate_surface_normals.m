%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_surface_normals
%
% Function to calculate the surface normals for the closed path envelope
% vector
%
% INPUT:
%   a polygon representing a closed loop of the envelope vector
%
% OUTPUT:
%   vector of surface normals for the closed envelope vector
%
% Lowe, A.R. 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [normals] = calculate_surface_normals(polygon)
% remove the last point?
polygon(length(polygon),:) = [];
normals = [];
for i = 0:length(polygon)-1
    
    pre = mod(i-1,length(polygon))+1;
    post = mod(i+1,length(polygon))+1;
    
    currP = polygon(i+1,:);
    preP  = polygon(pre,:);
    postP = polygon(post,:);
    
    vec1  = currP - preP;
    vec2  = postP - currP;
    
    vecSum = (vec1 + vec2);
    normal = [-vecSum(2), vecSum(1)] ./ sqrt((vecSum(2)^2 + vecSum(1)^2));
    
    % store the normal of this vertex
    normals = cat(1,normals,normal);
end

return