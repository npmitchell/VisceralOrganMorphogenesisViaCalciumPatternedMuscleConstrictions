function [phi0s, coscoeffs] = phiOffsetsFromPrevMeshWithDilation(TF, TV2D, TV3Drs, uspace, ...
    vspace, prev3d_sphi, visualize, options)
%PHIOFFSETSFROMPREVMESH(TF, TV2D, TV3Drs, nU, vpsace, prev3d_sphi) 
%   Find the offset in phi (the y dimension of the 2d pullback) that
%   minimizes the difference in 3D of the positions of each DV hoop from
%   that of the previous timepoint.
%
% Parameters
% ----------
% TF : nU*nV x 3 int array
%   The mesh connectivity list, indexing into the vertex arrays TV2D and
%   TV3Drs
% TV2D : nU*nV x 2 float array
%   The mesh vertex locations in 2d
% TV3Drs : nU*nV x 3 float array
%   The mesh vertex locations in 3d
% uspace : nU float array
%   The values of u for each line of constant v in pullback space
% vspace : nV float array
%   The values of v for each line of constant u in pullback space
% prev3d_sphi : nU x nV x 3 float array
%   The 3D coordinates of the embedding for the reference timepoint
%   (previous timepoint, for ex) at the 2D locations given by uspace and 
%   vspace. Note that uspace is not used explicitly, only nU = length(uspace) 
% 
% Returns
% -------
% phi0s : nV x 1 float array
%   the offsets in V dimension that minimize distances of each DV hoop from
%   analogous hoop in previous xyz (from previous timepoint, for ex)
%
%
% NPMitchell 2019

% Interpret vargin as boolean for visualization 
if nargin > 6
    if visualize
        fig = figure('visible', 'on') ;
    end
else
    visualize = false ;
end

if nargin < 8
    % options = optimset('PlotFcns','optimplotfval','TolX',1e-7);
    options = optimset() ; 
end

% Consider each value of u in turn
% Fit for phi0 such that v = phi - phi0
nU = length(uspace) ;
nV = length(vspace) ;
phi0s = zeros(nU, 1) ;
coscoeffs = zeros(nU, 1) ;
% The V values here are (0...1)
vqq = vspace ;

prog = repmat('.', [1 floor(nU/10)]) ;
for qq = 1:nU
    tic
    % curve = curves3d(qq, :) ;
    % The previous 3d embedding values are stored 
    prev3dvals = squeeze(prev3d_sphi(qq, :, :)) ;
    
    % Check it
    % plot3(prev3dvals(:, 1), prev3dvals(:, 2), prev3dvals(:, 3), '.')
    % hold on;
    % tmpx = prev3d_sphi(:, :, 1) ;
    % tmpy = prev3d_sphi(:, :, 2) ;
    % tmpz = prev3d_sphi(:, :, 3) ;
    % scatter3(tmpx(:), tmpy(:), tmpz(:), 2, 'MarkerFaceAlpha', 0.1,...
    %     'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'c')
        
    % Note: interpolate2Dpts_3Dmesh(cutMeshrs.f, cutMeshrs.u, cutMeshrs.v, uv) 
    
    % 
    phiopt = fminsearch(@(phi0)...
        sum(vecnorm(...
        interpolate2Dpts_3Dmesh(TF, TV2D, ...
            TV3Drs, [uspace(qq) * ones(nV, 1),...
            mod(vqq + phi0(1) + phi0(2) * cos(vqq * 2 * pi), 1)]) ...
            - prev3dvals, 2, 2).^ 2), [0., 0.], options);
    phi0s(qq) = phiopt(1) ;
    coscoeffs(qq) = phiopt(2) ;
    
    % Old version
    % phi0s(qq) = fminbnd(@(phi0)...
    % sum(vecnorm(...
    % interpolate2Dpts_3Dmesh(TF, TV2D, ...
    %     TV3Drs, [uspace(qq) * ones(nV, 1), mod(vqq + phi0(1), 1)]) ...
    %     - prev3dvals, 2, 2) .^ 2), lowerbound, upperbound, options);
        
    % Visualize the minimization output values
    if visualize
        plot(phi0s)
        title(['computed $phi_0$' num2str(qq) ' / ' num2str(nU)], ...
            'Interpreter', 'Latex')
        pause(0.000000001)
    end  
    runtimeIter = toc ;
    % Display progress bar
    if mod(qq, 10) == 1
        prog(min(length(prog), max(1, floor(qq/10)))) = '*' ;
        fprintf([prog '(' num2str(runtimeIter) 's per u value)\n'])
    end
end
if visualize
    close all
end
end

