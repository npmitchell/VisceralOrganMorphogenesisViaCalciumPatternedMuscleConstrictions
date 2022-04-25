function [length_lobes, area_lobes, volume_lobes] = ...
    aux_compute_lobe_dynamics(folds, ssfold, ssmax, lobeDir, timePoints, ...
    t0, timeInterval, timeUnits, spaceUnits, ...
    spcutMeshBase, nV, nU, rot, trans, resolution, flipy, xyzlim, colors, save_ims, overwrite_lobeims)
%AUX_COMPUTE_LOBE_DYNAMICS auxiliary function for Generate_Axisymmetric_Pullbacks_Orbifold.m
%   Compute the lobe dynamics for all timepoints
% 
%   todo: generalize to #lobes = other than 4
%
% Parameters
% ----------
% folds : 
% ssfold :
% ssmax : 
% lobeDir :  str
%   the path to where lobe information/data/images are stored
% timePoints : #timepoints x 1 float array
%   the time in minutes of each frame
% spcutMeshBase : str
%   full file path to the (s,phi) coord sys cylinder cut mesh
% nV : int
%   sampling number along DV axis
% nU : int 
%   sampling number along AP axis
% rot : 3x3 float array
%   rotation matrix transforming data into APDV coord sys
% trans : 1x3 float array
%   translation to perform after rotation to move APDV coord sys to origin
% resolution : float 
%   um / pixel conversion for data 
% colors : (>=#lobes)x3 float array
%   colors to use to color each lobe
% save_ims : bool
%   save an image of the mesh in 3d, with each lobe colored differently
% 
%
%
% NPMitchell 2019

disp('Computing length_lobes...')
% Length is given by ssfold.
ss_lobes = cat(2, ssfold, ssmax) ;
length_lobes = ss_lobes ;
for qq=2:4
    length_lobes(:, qq) = ss_lobes(:, qq) - ss_lobes(:, qq-1) ;
end

disp('Computing area_lobes and volume_lobes...')
% Now compute Surface area and Volume in one pass.
lobeImDir = fullfile(lobeDir, 'images_lobes3D') ;
if ~exist(lobeImDir, 'dir')
    mkdir(lobeImDir)
end

% Preallocate arrays
area_lobes = zeros(length(timePoints), 4) ;
volume_lobes = zeros(length(timePoints), 4) ;

disp('Computing surface area and volume of each lobe...')
for kk = 1:length(timePoints)
    % Translate to which timestamp
    t = timePoints(kk) ;
    timestr = sprintf('%04d', t) ;
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;

    % Load the centerline too
    % fn = sprintf(clineDVhoopBase, t) ;
    % load(fn, 'avgpts')
    avgpts = spcutMesh.avgpts ;

    % preallocate
    vol_kk = zeros(size(folds, 2) + 1, 1) ;
    area_kk = zeros(size(folds, 2) + 1, 1) ;

    % Create figure for plotting lobes
    lobeimfn = fullfile(lobeImDir, ['lobes_' timestr '.png']) ;
    redo_lobeims = save_ims && (~exist(lobeimfn, 'file') || overwrite_lobeims) ;
    if redo_lobeims
        close all
        fig = figure('visible', 'off') ;
    end

    % Announce which timestamp we consider
    if mod(t, 40) == 0
        disp(['SA and Vol: t = ', num2str(t)])
        if redo_lobeims
            disp([' ... and saving each figure to ' lobeimfn])
        end
    end

    % Define rotated, translated mesh
    vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
    % use QS to allow flipy
    if flipy
        vrs(:, 2) = -vrs(:, 2) ;
    end
    
    for lobe = 1:size(folds, 2) + 1
        % Note that the vertices are ordered in AP strips, with
        % nU elements for each value of nV. 

        % This is for transposed ordering
        % if lobe == 1
        %     rmvtx = ((f1+1)*nV):size(vrs, 1) ;
        %     rear = ((f1-1) * nV + 1):(((f1 + 1) * nV) - 1) ;
        %     front = 1:nV ;
        % elseif lobe == 2
        %     rmvtx = [1:(f1*nV), ((f2 + 1) * nV):size(vrs, 1) ] ;
        %     rear = ((f2-1) * nV + 1):(((f2 + 1) * nV) - 1) ;
        %     front = ((f1-1) * nV + 1):(((f1 + 1) * nV) - 1) ;
        % elseif lobe == 3
        %     rmvtx = [1:(f2*nV), ((f3 + 1) * nV):size(vrs, 1) ] ;
        %     rear = ((f3-1) * nV + 1):(((f3 + 1) * nV) - 1) ;
        %     front = ((f2-1) * nV + 1):(((f2 + 1) * nV) - 1) ;
        % elseif lobe == 4
        %     rmvtx = 1:(f3*nV) ;
        %     rear = ((f4-1) * nV + 1):(((f4 + 1) * nV) - 1) ;
        %     front = ((f3-1) * nV + 1):(((f3 + 1) * nV) - 1) ;
        % end

        % Create matrix of indices, each row identical (y location)
        allvs = (0:(nV-1)) * nU ;
        if lobe == 1
            % rename the fold indices (in U)
            front_id = 1 ;
            rear_id = folds(kk, 1) ;
            urm = (rear_id+1):nU ;
        elseif lobe == size(folds, 2) + 1
            front_id = folds(kk, lobe-1) ;
            rear_id = nU ;
            urm = 1:(front_id-1) ;
        else
            front_id = folds(kk, lobe-1) ;
            rear_id = folds(kk, lobe) ;
            urm = [1:(front_id-1), (rear_id+1):nU] ;
        end
        rear_id = double(rear_id) ;
        front_id = double(front_id) ;
        urm = double(urm) ;

        % Find all indices on the rear hoop and front hoop
        rear = (rear_id - front_id + 1) + ...
            (0:(nV-1)) * (rear_id - front_id + 1) ;
        front = 1 + (0:(nV-1)) * (rear_id - front_id + 1) ;
        % Create brick of indices, each column identical (x location)
        rmvtx = (urm' .* ones(length(urm), nV))' ;
        % Add y location and unravel
        rmvtx = rmvtx + allvs' .* ones(nV, length(urm)) ;
        rmvtx = rmvtx(:) ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute area of mesh up to first fold
        % extract faces of first lobe
        % remove the faces that do not include 1:fold1
        [ newF, newV, oldVertexIDx ] = remove_vertex_from_mesh(spcutMesh.f, ...
            vrs, rmvtx) ;
        % Compute area without closing the volume
        area_kk(lobe) = meshArea(newV, newF) ;

        % Close the lobe volume with additional triangles
        % first close the rear
        addpt = avgpts(rear_id, :) ;
        newV = cat(1, newV, addpt) ;
        addptidx = length(newV);
        backf = zeros(length(rear), 3) ;
        for qq = 1:length(rear)
            nid = mod(qq + 1, length(rear)) ;
            if nid == 0
                nid = length(rear) ;
            end
            backf(qq, :) = [rear(qq), rear(nid), addptidx];
        end
        newF = cat(1, newF, backf) ;

        % Close the front face
        addpt = avgpts(front_id, :) ;
        newV = cat(1, newV, addpt) ;
        addptidx = length(newV);
        frontf = zeros(length(front), 3) ;
        for qq = 1:length(front)
            nid = mod(qq + 1, length(front)) ;
            if nid == 0
                nid = length(front) ;
            end
            frontf(qq, :) = [front(qq), addptidx, front(nid)];
        end
        newF = cat(1, newF, frontf) ;

        % Finally, compute volume
        [vol_kk(lobe), ~] = meshVolumeArea(newV, newF) ;

        % Plot lobes in 3D
        if redo_lobeims
            trisurf(newF, newV(:, 1), newV(:, 2), newV(:, 3), ...
                'EdgeColor', colors(lobe, :) * 0.5, 'FaceColor', colors(lobe, :));
            hold on;
        end
    end

    % Add area_kk and vol_kk to lists
    area_lobes(kk, :) = area_kk ;
    volume_lobes(kk, :) = abs(vol_kk) ;   

    % Save plot of lobes, each a different color
    if redo_lobeims
        axis equal
        title(['Lobes, t = ' sprintf('%03d', (t-t0)*timeInterval) ' ' timeUnits])
        xlabel(['x [' spaceUnits ']'])
        ylabel(['y [' spaceUnits ']'])
        zlabel(['z [' spaceUnits ']'])
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        disp(['Saving figure: ' lobeimfn])
        saveas(fig, lobeimfn)
        close all
    end
end
