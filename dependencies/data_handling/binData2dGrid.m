function [zmeans, counts, zs, xidx, yidx] = ...
    binData2dGrid(uvz, uminmax, vminmax, nU, nV, uniqueOnly)
% binData2dGrid(uvz, uminmax, vminmax, nU, nV)
%   Bin 3d data into 2d grid and take means of z values. 
%   Maps values into bins around bin centers (umin, ..., umax) and 
%   (vmin, ..., vmax)
%
% Parameters
% ----------
% uvz : N x 3 numeric array
%   the (uvz) data to bin in (u,v)
% zs : 
% uu : Q x 1 numeric array
%   sorted 
% vv : R x 1 numeric array
%   sorted indices of bin
% nU : int
%   number of bins in x/u direction
% nV : int
%   number of bins in y/v direction
%   
%
% Returns
% -------
% zmeans : nU * nV numeric array 
%   mean values in each bin. Note that the output x dimension is the u
%   dimension (ie the first dimension), so to show the result in an axis
%   use imshow(zmeans') -- ie, transpose the result
% counts : nU * nV int array
%   how many scattered values appeared in each bin
% zs : nU * nV cell array 
%   lists of the z values that fall into each unique x/y combination
% xidx : nU x 1 int array
%   the "column" (x bin) at which each sample value resides
% yidx : nV x 1 int array
%   the "row" (y bin) at which each sample value resides
%
% NPMitchell 2020 adaptation of Walter Robinson's suggestion from ref: 
% https://www.mathworks.com/matlabcentral/answers/
%       322113-binning-data-with-2d-coordinates

if nargin < 6
    uniqueOnly = false ;
end

% Round uvz(:, [1,2]) onto u_vals, v_vals
umin = uminmax(1) ;
umax = uminmax(2) ;
vmin = vminmax(1) ;
vmax = vminmax(2) ;
urounded = round((uvz(:, 1) - umin) * (nU-1) / (umax - umin)) ;
vrounded = round((uvz(:, 2) - vmin) * (nV-1) / (vmax - vmin)) ;
urounded(urounded < 0 ) = 0 ;
urounded(urounded > nU-1) = nU-1 ;
vrounded(vrounded < 0 ) = 0 ;
vrounded(vrounded > nV-1) = nV-1 ;

% Check it
% scatter(urounded, vrounded, 5, uvz(:, 3))

% Now find indices of the rounded data to accumulate
if uniqueOnly
    [~, ~, xidx] = unique(urounded);
    [~, ~, yidx] = unique(vrounded);

    % count the number of points at each unique x/y combination
    counts = accumarray([xidx(:), yidx(:)], 1);  
    % average the z that fall into each unique x/y combination
    sums = accumarray([xidx(:), yidx(:)], uvz(:,3));
    zmeans = sums ./ counts ;

else
    xidx = discretize(urounded, -0.5:nU+0.5) ;
    yidx = discretize(vrounded, -0.5:nV+0.5) ;
    
    % count the number of points at each unique x/y combination
    counts = accumarray([xidx(:), yidx(:)], 1);  
    % average the z that fall into each unique x/y combination
    sums = accumarray([xidx(:), yidx(:)], uvz(:,3));
    zmeans = sums ./ counts ;
end

% Check it
% scatter(xidx, yidx, 5, uvz(:, 3))

% Nans at locations in which no measurements were made
zmeans(counts == 0) = NaN ;

% Ensure that output is nU * nV
if any(size(zmeans) ~= [nU, nV])
    % First check that xidx (yidx) runs within [1, nU] ([1,nV])
    assert(all(xidx >= 1) & all(xidx <= nU))
    assert(all(yidx >= 1) & all(yidx <= nV))
    
    % Pad zmeans with nans 
    zmeans_new = nan(nU, nV) ;
    counts_new = zeros(nU, nV) ;
    if uniqueOnly
        zmeans_new(min(xidx):max(xidx), min(yidx):max(yidx)) = zmeans ;
        counts_new(min(xidx):max(xidx), min(yidx):max(yidx)) = counts ;
    else
        zmeans_new(1:max(xidx), 1:max(yidx)) = zmeans ;    
        counts_new(1:max(xidx), 1:max(yidx)) = counts ;    
    end
    zmeans = zmeans_new ;
    counts = counts_new;
    zmeans(counts == 0) = NaN ;
end

% Check it
% scatter(xidx, yidx, 1+counts(:))
% imagesc(isnan(zmeans))
% imagesc(zmeans) ;
% hold on;
% uu = (uvz(:, 1) - umin) * nU / (umax-umin) ; 
% vv = (uvz(:, 2) - vmin) * nV / (vmax-vmin) ; 
% scatter(uu, vv, 5, uvz(:, 3))

if nargout > 2
    %create a list of the z that fall into each unique x/y combination
    if uniqueOnly
        zs = accumarray([xidx(:), yidx(:)], uvz(:,3), [], @(V) {V}, {});
    else
        zs_full = cell(nU, nV) ;
        zs = accumarray([xidx(:), yidx(:)], uvz(:,3), [], @(V) {V}, {});
        zs_full(1:size(zs, 1), 1:size(zs, 2)) = zs ;
        zs = zs_full ;
    end
end

