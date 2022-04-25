% Analyze the GCAMP pulses for width of anterior fold using Mef2GAL4
% NPMitchell 2021

clear all
preview = false ;
overwrite = true ;
resDirFn = 'anteriorMutantResults' ;
pausetime = 3 ;

gutMatlabDir = '/mnt/data/code/gut_matlab/' ;
% gutMatlabDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/gut_matlab/' ;

addpath(fullfile(gutMatlabDir, 'addpath_recurse')) ;
addpath_recurse(gutMatlabDir) 

% datdir = '/mnt/data/confocal_data/gut/Mef2GAL4klarUASGCAMP6sIII/analysis_mef2G4kGCaMP6sIII/' ;
% datdir = '/Users/npmitchell/Desktop/gut/GCaMP6sIII_analysis/Mef2GAL4klarGCaMP6sIII_analysis' ;
datdir='/mnt/data/confocal_data/gut/48YGAL4GCaMP6sIIAntpNSRC3/analysis_AntpNSRC3_48YGAL4GCaMP6sII/';

expts = {...
    '202107142134_AntpNSRC3_e2',...
    '202107151435_AntpNSRC3_e3',...
    '202107172300_AntpNSRC3_e3', ...
    '202107181225_AntpNSRC3_e2',...
    '202107182209_AntpNSRC3_e2', ...
    } ;

% poster frames: 202106211253_e3, 202106211440_e1, 202106211440_e2

clipY0s = {[14,14+265], ...             % 202107142134 e2 
    [69,69+270], ...   % 202107151435 e3 --> if not rot, then [56,250];
    [63,63+260], ...             % 202107172300_AntpNSRC3_e3
    [77,232], ...             % 202107181225_AntpNSRC3_e2
    [44,44+265], ...       %   202107182209_AntpNSRC3_e2
    }; 
clipYAdjs = {[0, 0], [50,-50]}; 
clipXs = {[1, Inf], [1, Inf], [1, Inf], ...
    [1, Inf], [1, Inf]} ;
% Timestep in MINUTES
dts = 1/60 * [91,91,93,93,93] ;
% Conversion in micron / XY pixels
pix2um = [0.3022826,...           % 202107142134_AntpNSRC3 e2
    0.33215031289, ...      % 202107151435 e3
    0.3305353,...           % 202107172300_AntpNSRC3_e3
    0.33215, ...            % 202107181225_AntpNSRC3_e2
    0.329999, ...           % 202107182209_AntpNSRC3_e2
    ] ;

pcPower = [12, 12, 15, 15, 15] ;
dates = [15, 14, 17, 18, 18, 18] ;
dz = [2.5, 2.5, 3, 3, 3 ] ;
% DEFINE FOLD TIME AS time of 10 micron indentation on ventral side at
% 16*2.5um = 40 um into the tissue -- this is basically the saggital plane
% fold X position at frame foldT, in pixels
anteriorXs = [79, ...             % 202107142134 e2
    118, ...            % 202107151435 e3
    147, ...             % 202107172300_AntpNSRC3_e3
    132, ...             % 202107181225_AntpNSRC3_e2
    173, ...            %   202107182209_AntpNSRC3_e2
    ] ;          
% fold frame # at which position is measured
posteriorTs = [14, ...              % 202107142134 e2
    10, ...             % 202107151435 e3
    32, ...             % 202107172300_AntpNSRC3_e3
    54, ...             % 202107181225_AntpNSRC3_e2
    64, ...             %   202107182209_AntpNSRC3_e2
    ]; 
% onset of POSTERIOR folding in minutes
t0 = dts .* posteriorTs ;

% Fixed xlimits in microns
xlimFix = [0, 100] ;
ylimFix = [] ;
ylimFixAll = [-15, 90] ;
fixXticks = [0, 25, 50, 75, 100] ;
cLimits = [0, 10] ;


% Filtering of drift deltaX options
capDeltaX = [0, 0, 0, ...
    0, 0] ;
medFiltDeltaXW = [0, 0, 0, ...
    0, 0] ;

xfixed = linspace(0, 100, 100/0.33) ;

save(fullfile(datdir, resDirFn, 'AntpNSRC3HeterozSettings.mat'), ...
    'clipY0s', 'dts', 'pix2um', 'anteriorXs', 'posteriorTs', ...
    'pcPower', 'dz', 'dates', 'expts', 'xfixed') ;

% Antp domain is about 40 microns, broadens on dorsal side.
close all

% If drift in X, define it here as [frame#, offsetX] 
for ee = 1:length(expts)
    fns1 = dir(fullfile(datdir, expts{ee}, 'images', '*c001.png')) ;
    fns2 = dir(fullfile(datdir, expts{ee}, 'images', '*c002.png')) ;
    fns3 = dir(fullfile(datdir, expts{ee}, 'images', '*c003.png')) ;
    
    % For each yrange
    for clipyPairIdx = 1:length(clipYAdjs)
        outfn = sprintf([expts{ee} '_results_Yrange%d.mat'], clipyPairIdx) ;
        if ~exist(fullfile(datdir, resDirFn, outfn), 'file') || overwrite
            if ~exist(fullfile(datdir, resDirFn, outfn), 'file')  
                disp(['results not on disk: ee = ' num2str(ee)])
            else
                disp(['overwriting: ee = ' num2str(ee)])
            end
            % this y range is defined
            clipY = clipY0s{ee} + clipYAdjs{clipyPairIdx} ;
            clipX = clipXs{ee} ;
            foldX = anteriorXs(ee) ;
            foldT = posteriorTs(ee) ;

            % Look at differences along y for each timepoint
            ntps = length(fns1) ;
            cmap = viridis(ntps) ;
            timestamps = dts(ee) * (0:ntps-1) ;
            timestamps = timestamps - t0(ee) ;
            
            % Hack
            % if clipyPairIdx == 1
            %     deltaX = zeros(ntps, 1) ;
            % end
            
            clf
            for ii = 1:ntps
                im1 = double(imread(fullfile(fns1(ii).folder, fns1(ii).name))) ;
                im2 = double(imread(fullfile(fns3(ii).folder, fns2(ii).name))) ;
                im3 = double(imread(fullfile(fns3(ii).folder, fns3(ii).name))) ;

                if isfinite(clipX(2))
                    c1 = im1(clipY(1):clipY(2), clipX(1):clipX(2)) ;
                    c2 = im2(clipY(1):clipY(2), clipX(1):clipX(2)) ;
                    c3 = im3(clipY(1):clipY(2), clipX(1):clipX(2)) ;
                else
                    c1 = im1(clipY(1):clipY(2), clipX(1):end) ;
                    c2 = im2(clipY(1):clipY(2), clipX(1):end) ;
                    c3 = im3(clipY(1):clipY(2), clipX(1):end) ;
                end
                
                % mean image
                curr = (c1 + c2 + c3) / 3 ; 
                
                if ii == 50
                    pause(1)
                end
                
                % top-hat TRANSFORM 
                % First do top-hat FILTERING, then subtract
                % -- make an elliptical strel
                % d12 = abs(c1 - c2) ; 
                % d13 = abs(c1 - c3) ;
                % d23 = abs(c2 - c3) ;

                c1fg = imgaussfilt(c1, 1 / pix2um(ee)) ;
                c2fg = imgaussfilt(c2, 1 / pix2um(ee)) ;
                c3fg = imgaussfilt(c3, 1 / pix2um(ee)) ;
                d12 = abs(c1fg - c2fg) ; 
                d13 = abs(c1fg - c3fg) ;
                d23 = abs(c2fg - c3fg) ;

                
                % KEEP BIGGER
                % lOpening = 2 ;  % DV extent of opening
                % wOpening = 2 ;  % AP extent of opening
                % se_nlong = round(lOpening / pix2um(ee)) ;
                % se_nwide = round(wOpening / pix2um(ee)) ;
                % % se0 = strel('disk', se_nlong) ;
                % % se_nwide = round(wObjects / pix2um(ee)) ; % 2um diameter / pix2um
                % % se_sampl = round(size(se0.Neighborhood, 2) / se_nwide) ;
                % % seN = se0.Neighborhood(:, 1:se_sampl:end) ;
                % se = strel('arbitrary', ones([se_nlong, se_nwide])) ;
                % hatt1 = imtophat(d12, se) ;
                % d12new = d12 - hatt1 ;
                % figure(1); 
                % subplot(2, 2, 1); imagesc(d12); axis equal; colorbar
                % subplot(2, 2, 2); imagesc(hatt1); axis equal; colorbar
                % subplot(2, 2, 3); imagesc(d12-hatt1); axis equal; colorbar
                % subplot(2, 2, 4); imagesc(hatt1-d12); axis equal; colorbar
                
                % BLUR to KEEP BIGGER
                d12new = imgaussfilt(d12, 0.5 / pix2um(ee)) ;
                d13new = imgaussfilt(d13, 0.5 / pix2um(ee)) ;
                d23new = imgaussfilt(d23, 0.5 / pix2um(ee)) ;
                % d12blur = (d12blur - mean(d12blur) > 0) .* d12blur ;
                
                % KEEP SMALLER
                lOpening = 25 ;  % DV extent of opening
                wOpening = 25 ;  % AP extent of opening
                se_nlong = round(lOpening / pix2um(ee)) ;
                se_nwide = round(wOpening / pix2um(ee)) ;
                % se0 = strel('disk', se_nlong) ;
                % se_nwide = round(wObjects / pix2um(ee)) ; % 2um diameter / pix2um
                % se_sampl = round(size(se0.Neighborhood, 2) / se_nwide) ;
                % seN = se0.Neighborhood(:, 1:se_sampl:end) ;
                se = strel('arbitrary', ones([se_nlong, se_nwide])) ;
                hatt12 = imtophat(d12new, se) ;
                hatt23 = imtophat(d23new, se) ;
                hatt13 = imtophat(d13new, se) ;
                
                
                if preview && mod(ii, 25) == 0
                    figure(1); 
                    subplot(2, 2, 1); imagesc(d12new); axis equal; colorbar
                    subplot(2, 2, 2); imagesc(hatt12); axis equal; colorbar
                    subplot(2, 2, 3); imagesc(d12new-hatt12); axis equal; colorbar
                    subplot(2, 2, 4); imagesc(hatt12-d12new); axis equal; colorbar
                    sgtitle('d12 and hatt12')
                    pause(pausetime)
                    subplot(2, 2, 1); imagesc(d23new); axis equal; colorbar
                    subplot(2, 2, 2); imagesc(hatt23); axis equal; colorbar
                    subplot(2, 2, 3); imagesc(d23new-hatt23); axis equal; colorbar
                    subplot(2, 2, 4); imagesc(hatt23-d23new); axis equal; colorbar
                    sgtitle('d23 and hatt23')
                    pause(pausetime)
                    subplot(2, 2, 1); imagesc(d13new); axis equal; colorbar
                    subplot(2, 2, 2); imagesc(hatt13); axis equal; colorbar
                    subplot(2, 2, 3); imagesc(d13new-hatt13); axis equal; colorbar
                    subplot(2, 2, 4); imagesc(hatt13-d13new); axis equal; colorbar
                    sgtitle('d13 and hatt13')
                    pause(pausetime)
                    clf
                    plot(sum(hatt12, 1)); hold on; 
                    plot(sum(hatt23, 1)); hold on; 
                    plot(sum(hatt13, 1)); hold on; 
                    sgtitle('hatt12, hatt23, hatt13')
                    pause(pausetime)
                end
                
                if ii == 1 
                    if isfinite(clipX(2))
                        wX = clipX(2)-clipX(1) + 1 ;
                    else
                        wX = size(c1, 2)-clipX(1) + 1 ;
                    end
                    dd = ones(length(fns1), wX) ;
                    ss = ones(length(fns1), wX) ;
                    tophats = ones(length(fns1), wX) ;
                    xshifted = zeros(length(fns1), wX) ;
                    xx = pix2um(ee) * (0:wX-1) ;
                    dinterp = zeros(ntps, size(xfixed, 2)) ; 
                    sinterp = zeros(ntps, size(xfixed, 2)) ;
                    deltaX = zeros(ntps, 1) ;
                elseif clipyPairIdx < 3
                    % Find offset in x wrt previous image
                    assert(any(any(abs(curr-prev))))
                    try
                        [deltaX(ii), ~] = ExtPhaseCorrelation(curr, prev) ;
                        LX = clipXs{ee}(2) - clipXs{ee}(1) ;
                        if abs(deltaX(ii)) > capDeltaX(ee)
                            deltaX(ii) = 0 ;
                        end
                        if deltaX(ii) > LX * 0.5 
                            deltaX(ii) = mod(deltaX(ii), LX) ;
                        elseif deltaX(ii) < -LX* 0.5 
                            deltaX(ii) = -mod(abs(deltaX(ii)), LX) ;                    
                        end
                    catch
                        disp('WARNING: could not compute phase corr!')
                        deltaX(ii) = 0 ;
                    end
                else
                    deltaXfn = sprintf([expts{ee} '_results_Yrange%d.mat'], 2) ;
                    load(fullfile(datdir, resDirFn, deltaXfn), 'deltaX')
                end

                % s1 = sum(c1, 1) ;
                % s2 = sum(c2, 1) ;
                % s3 = sum(c3, 1) ;                 
                % d12 = abs(s1 - s2); 
                % d13 = abs(s1 - s3) ;
                % d23 = abs(s2 - s3) ;
                % d12 = abs(sum(c1 - c2, 1)) ; 
                % d13 = abs(sum(c1 - c3, 1)) ;
                % d23 = abs(sum(c2 - c3, 1)) ;

                % Save as image
                diffDir = sprintf(fullfile(fns1(ii).folder, 'diffs_clipY%02d'), clipyPairIdx) ;
                if ~exist(diffDir, 'dir')
                    mkdir(diffDir)
                end
                fn12 = fullfile(diffDir, [fns1(ii).name(1:end-4) '_d12.png']) ;
                fn23 = fullfile(diffDir, [fns1(ii).name(1:end-4) '_d23.png']) ;
                fn13 = fullfile(diffDir, [fns1(ii).name(1:end-4) '_d13.png']) ;
                assert(max(hatt12(:)) < 255)
                % imwrite(uint16(hatt12 * 2^8), fn12) ;
                % imwrite(uint16(hatt23 * 2^8), fn23) ;
                % imwrite(uint16(hatt13 * 2^8), fn13) ;
                
                % Save transient signal in matrix
                sh12 = sum(hatt12, 1) ;
                sh23 = sum(hatt23, 1) ; 
                sh13 = sum(hatt13, 1) ;
                ddii = mean([sh12; sh23; sh13], 1) ;
                tophats(ii, :) = ddii ;
                
                % Background estimation
                sii = min([sh12; sh23; sh13], [], 1) ;
                ss(ii, :) = sii ;
                dd(ii, :) = ddii - sii ;

                % Save image for phase correlation with next image
                prev = curr ;
            end

            % Median filter deltaX
            if medFiltDeltaXW(ee) > 0
                deltaX = movmedian(deltaX, medFiltDeltaXW(ee)) ;
            end
            
            % Register the position so that (1) stabilized and (2) 0um is anterior fold
            for ii = 1:ntps
                % recall dd and ss for this tp
                ddii = dd(ii, :) ;
                sii = ss(ii, :) ;

                % Save shifted ap positions in matrix
                xshiftii = xx+sum(deltaX(1:ii)) ;
                xshifted(ii, :) = xshiftii ;
            end
            % Now overall offset
            xshiftPix = xshifted(foldT, foldX) ;
            xshifted = xshifted - xshiftPix ;
            
            for ii = 1:ntps
                % recall dd and ss for this tp
                ddii = dd(ii, :) ;
                sii = ss(ii, :) ;

                % Save shifted ap positions in matrix
                xshiftii = xshifted(ii, :) ;
                
                % Interpolate shifted results
                dinterp(ii, :) = interp1(xshiftii, ddii, xfixed) ;
                sinterp(ii, :) = interp1(xshiftii, sii, xfixed) ;

                % Factor in translation of the embryo over the timecourse
                plot(xx+sum(deltaX(1:ii)), dd(ii, :), 'color', cmap(ii, :)) ; hold on;

            end

            % First plot -- raw curves
            set(gcf, 'Visible', 'off')
            xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('GCAMP transient signal [a.u.]', 'interpreter', 'latex')    
            title('Transient GCAMP activity')
            cb = colorbar ;
            caxis([min(timestamps), max(timestamps)])
            ylabel(cb, 'time [min]')
            outfigfn = sprintf(fullfile(datdir, expts{ee}, '00_difference_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn) 

            % Kymograph -- pre shift
            clf
            imagesc(xx, timestamps, dd)
            xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('time [min]', 'interpreter', 'latex')
            title('Transient GCAMP activity')
            cb = colorbar ;
            ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
            outfigfn = sprintf(fullfile(datdir, expts{ee}, '01_difference_kymo_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn)

            % Shifted kymo
            clf
            for row = 1:ntps
                imagesc(xshifted(row, :), timestamps(row), dd(row, :));
                hold on;
            end
            xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('time [min]', 'interpreter', 'latex')
            title('Transient GCAMP signal', ...
                'interpreter', 'latex')
            xlim(xlimFix)
            if ~isempty(ylimFix)
                ylim(ylimFix)
            else
                ylim([min(timestamps), max(timestamps)])
            end
            cb = colorbar ;
            ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
            outfigfn = sprintf(fullfile(datdir, expts{ee}, '02_difference_kymoShift_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn)
            
            
            % PLot deltaX
            clf
            subplot(2, 1, 1)
            plot(timestamps, deltaX, '.-') ;
            ylabel('$\delta x$ [$\mu$m]', 'interpreter', 'latex')
            xlabel('time [min]', 'interpreter', 'latex')
            subplot(2, 1, 2)
            plot(timestamps, cumsum(deltaX), '.-') ;
            ylabel('$\int \delta x$ [$\mu$m]', 'interpreter', 'latex')
            xlabel('time [min]', 'interpreter', 'latex')
            outfigfn = sprintf(fullfile(datdir, expts{ee}, '02b_shift_deltaX_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn)

            % Plot interpolated results
            clf
            timeGrid = timestamps' .* ones(size(dinterp)) ;
            apGrid = xfixed .* ones(size(dinterp)) ;
            imagesc(xfixed, timestamps, dinterp) ;
            xlim(xlimFix)
            if ~isempty(ylimFix)
                ylim(ylimFix)
            else
                ylim([min(timestamps), max(timestamps)])
            end
            xlabel('stabilized, re-gridded ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('time [min]', 'interpreter', 'latex')
            title('Transient GCAMP signal [registered frame]', ...
                'interpreter', 'latex')
            cb = colorbar ;
            ylabel(cb, '$ \delta I $ [a.u.]', 'interpreter', 'latex')
            outfigfn = sprintf(fullfile(datdir, expts{ee}, ...
                '03_difference_kymoShiftInterp_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn)
            xticks(fixXticks)
            ylim(ylimFixAll)
            saveas(gcf, ...
                fullfile(datdir, resDirFn, sprintf('expt%02d_kymo_clipY%d.png', ...
                ee, clipyPairIdx)))

            % Average over time
            % Get t=0 idx
            [~, tidx0] = min(abs(timestamps)) ;
            end5 = find(timestamps - 5 > 0, 1);
            end10 = find(timestamps - 10 > 0, 1);
            end15 = find(timestamps - 15 > 0, 1);
            end20 = find(timestamps - 20 > 0, 1);
            end25 = find(timestamps - 25 > 0, 1);
            end30 = find(timestamps - 30 > 0, 1);
            end45 = find(timestamps - 45 > 0, 1);
            end60 = find(timestamps - 60 > 0, 1);
            activity00 = mean(dinterp(1:tidx0, :))' ;
            activity05 = mean(dinterp(tidx0:end5, :))' ;
            activity10 = mean(dinterp(tidx0:end10, :))' ;
            activity15 = mean(dinterp(tidx0:end15, :))' ;
            activity20 = mean(dinterp(tidx0:end20, :))' ;
            activity25 = mean(dinterp(tidx0:end25, :))' ;
            activity30 = mean(dinterp(tidx0:end30, :))' ;
            activity45 = mean(dinterp(tidx0:end45, :))' ;
            activity60 = mean(dinterp(tidx0:end60, :))' ;
            legends = {'$\langle t<$0 minutes$\rangle$', ...
                '$\langle 0<t<$5 minutes$\rangle$', ...
                '$\langle 0<t<$10 minutes$\rangle$', ...
                '$\langle 0<t<$15 minutes$\rangle$', ...
                '$\langle 0<t<$20 minutes$\rangle$', ...
                '$\langle 0<t<$25 minutes$\rangle$', ...
                '$\langle 0<t<$30 minutes$\rangle$'} ;
            colors = flipud(viridis(7)) ;

            clf
            plot(xfixed, activity00, 'color', colors(1, :)); hold on;
            plot(xfixed, activity05, 'color', colors(2, :))
            plot(xfixed, activity10, 'color', colors(3, :))
            plot(xfixed, activity15, 'color', colors(4, :))
            plot(xfixed, activity20, 'color', colors(5, :))
            plot(xfixed, activity25, 'color', colors(6, :))
            plot(xfixed, activity30, 'color', colors(7, :))
            legend(legends, 'interpreter', 'latex')
            xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('GCAMP transient activity [a.u.]', 'interpreter', 'latex')    
            title('Transient GCAMP activity', ...
                'interpreter', 'latex')
            xlim(xlimFix)
            outfigfn = sprintf(fullfile(datdir, expts{ee}, '07_difference_meanBgSub_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn) 

            % Filter in space
            close all
            w1um = round(1 / pix2um(ee)) ;
            plot(xfixed, movmedian(activity00, w1um), 'color', colors(1, :)); 
            hold on;
            plot(xfixed, movmedian(activity05, w1um), 'color', colors(2, :)); 
            plot(xfixed, movmedian(activity10, w1um), 'color', colors(3, :)); 
            plot(xfixed, movmedian(activity15, w1um), 'color', colors(4, :)); 
            plot(xfixed, movmedian(activity20, w1um), 'color', colors(5, :)); 
            plot(xfixed, movmedian(activity25, w1um), 'color', colors(6, :)); 
            plot(xfixed, movmedian(activity30, w1um), 'color', colors(7, :)); 
            legend(legends, 'interpreter', 'latex')
            xlabel('stabilized ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('GCAMP transient activity [a.u.]', 'interpreter', 'latex')    
            title('Transient GCAMP activity, filtered in space 1 $\mu$m', ...
                'interpreter', 'latex')
            xlim(xlimFix)
            outfigfn = sprintf(fullfile(datdir, ...
                expts{ee}, ...
                '08_difference_meanBgSubSm_Yrange%d.png'), clipyPairIdx) ;
            saveas(gcf, outfigfn) 
            saveas(gcf, fullfile(datdir, resDirFn, ...
                sprintf('expt%02d_means_clipY%d_%s.png', ...
                ee, clipyPairIdx, expts{ee})))

            % Save result
            save(fullfile(datdir, resDirFn, outfn), 'xx', 'xshifted', 'dd', 'deltaX', ...
                'xfixed', 'dinterp', 'sinterp', 'timeGrid', 'apGrid', ...
                'activity00', 'activity05', 'activity10', ...
                'activity15', 'activity20', ...
                'activity25', 'activity30', ...
                'activity45', 'activity60')
        else
            load(fullfile(datdir, resDirFn, outfn), 'xx', 'xshifted', 'dd', 'deltaX', ...
                'xfixed', 'dinterp', 'sinterp', 'timeGrid', 'apGrid')
        end
    end
    
end


%% Average all experiments together
expts2include = 1:length(expts);
for clipyPairIdx = 1:2
    colors = define_colors ;
    fixTimeStamps = -20:1.5:90.5 ;  % minutes
    kymoM = zeros(length(xfixed), length(fixTimeStamps)) ;
    nsamples = kymoM ;
    
    % For different averaging
    avgMin = [15, 20, 25] ;
    allres = [] ;
    for avgID = 1:length(avgMin)
        edmy = 1 ;
        disp(['Performing ' num2str(avgMin(avgID)) ' min average'])
        close all
        statAll = zeros(length(xfixed), length(expts2include)) ;
        for ee = expts2include
            disp(['loading ee = ' num2str(ee)])
            outfn = sprintf([expts{ee} '_results_Yrange%d.mat'], clipyPairIdx) ;
            resfn = fullfile(datdir, resDirFn, outfn) ;
            load(resfn, 'xfixed', 'activity15', ...
                'activity20', 'activity25', 'activity30', ...
                'timeGrid', 'dinterp')
            timestamps = unique(timeGrid) ;
            acts = {activity15, activity20, activity25, activity30} ;

            w1um = round(1 / pix2um(ee)) ;
            activity = movmedian(acts{avgID}, w1um) ;
            padd = round(5/pix2um(ee)) ;
            
            % no need for this anymore -- all x are xfixed already
            % keep = xfixed > xlimFix(1) & xfixed < xlimFix(2) ;
            % idx2keep = find(keep) ;
            % normIdx = (idx2keep(1) : idx2keep(1)+padd) ;
            % normIdx = [normIdx, idx2keep(end)-padd:idx2keep] ;
            
            normVal = nanmedian(activity) ;
            bgVal = mean(mink([activity(1:padd); activity(end-padd:end)], 10)) ;
            maxVal = maxk(activity(1+padd:end-padd), 10) ;
            maxVal = mean(maxVal) ;
            anorm = (activity - bgVal) / (maxVal - bgVal); 

            % AVERAGE KYMO
            if avgID == 1
                dnorm = (dinterp - normVal) / (maxVal - normVal) ;
                intrp = griddedInterpolant({timestamps, xfixed'}, dnorm, 'nearest', 'none') ;
                [tt, xx] = ndgrid(fixTimeStamps, xfixed) ;
                newKymo = intrp(tt, xx) ;
                nsamples = nsamples + ~isnan(newKymo)' ;
                newKymo(isnan(newKymo)) = 0 ;
                kymoM = kymoM + newKymo' ;
            end
            
            if avgID == 4
                pause(0.001)
            end
            hold on; 
            plot(xfixed, anorm, 'color', colors(edmy, :)); 
            statAll(:, edmy) = anorm ;
            hold on;

            if ee == expts2include(1)
                xall = xfixed ;
            end
            % min(xfixed(keep))
            % max(xfixed(keep))

            % Collate all results into table
            % keep4all = ismember(xall, xfixed(keep)) ;
            allres(:, edmy) = anorm ;
            edmy = edmy + 1;
        end
        xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
        ylabel('normalized transient GCAMP activity [a.u.]', ...
            'interpreter', 'latex')
        title(sprintf('transient calcium activity, %d min average', avgMin(avgID)), ...
            'interpreter', 'latex')
        resfn = fullfile(datdir, resDirFn, 'gcamp_activity_results') ;
        saveas(gcf, [resfn sprintf('_clipY%d_avg%d.png', clipyPairIdx, avgMin(avgID))])
        saveas(gcf, [resfn sprintf('_clipY%d_avg%d.pdf', clipyPairIdx, avgMin(avgID))])

        
        % Take stats
        close all
        fig = figure('units','centimeters','position',[0,0,6,10]);
        
        % kymo panel
        subplot(2, 1, 1)
        kymoM(~isfinite(kymoM)) = NaN;
        kymoMean = (kymoM - min(kymoM(:))) ./ nsamples  ;
        imagesc(xfixed, fixTimeStamps, kymoMean' ); 
        caxis(cLimits)
        % xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
        ylabel('time from fold onset [min]', 'interpreter', 'latex')
        
        hold on;
        redcolor = colors(2, :) ;
        plot(xlimFix, 0*xlimFix + 5, '--', 'color', redcolor)
        plot(xlimFix, 0*xlimFix - 5, '--', 'color', redcolor)
        plot(xlimFix, 0*xlimFix, '-', 'color', redcolor)
        xticks(fixXticks)
        set(gca, 'xticklabels', [])
        ylim(ylimFixAll)
        colormap(viridis)
            
        % curve panel
        subplot(2, 1, 2)
        avgact = mean(statAll, 2, 'omitnan') ;
        stdact = std(statAll, 0, 2, 'omitnan') ;
        lineProps = {'-','color', colors(1, :)} ;
        factor = 1.0 / max(avgact(:)) ;
        avgact = movmean(avgact, 3) ;
        stdact = movmedian(stdact, 15) ;
        h1=shadedErrorBar(xfixed, avgact*factor, stdact*factor, 'lineProps', lineProps) ;
        % ylim([min(avgact*factor - stdact*factor), max(avgact*factor + stdact*factor)])
        ylim([-0.25, 1.25])
        xticks(fixXticks)
        xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
        ylabel('GCaMP activity [a.u.]', ...
            'interpreter', 'latex')
        title(sprintf('%d min average', avgMin(avgID)), ...
            'interpreter', 'latex')
        resfn = fullfile(datdir, resDirFn, 'gcamp_mean_results') ;
        saveas(gcf, [resfn sprintf('_clipY%d_avg%d.png', clipyPairIdx, avgMin(avgID))])
        saveas(gcf, [resfn sprintf('_clipY%d_avg%d.pdf', clipyPairIdx, avgMin(avgID))])
        
        % Plot mean kymo for first pass
        if avgID == 1
            close all; 
            chelixMap = cubehelix(128,0.43,-0.68,1.3,0.4,[0,1.0],[0,1.0]) ;
            kymoMean = (kymoM - min(kymoM(:))) ./ nsamples  ;
            imagesc(xfixed, fixTimeStamps, kymoMean' ); 
            cb = colorbar;
            caxis(cLimits)
            xlabel('ap position from anterior fold [$\mu$m]', 'interpreter', 'latex')
            ylabel('time from fold onset [min]', 'interpreter', 'latex')
            ylabel(cb,'transient GCAMP activity [a.u.]')
            xticks(fixXticks)
            resfn = fullfile(datdir, resDirFn, 'gcamp_kymo_results') ;
            colormap(viridis)
            saveas(gcf, [resfn sprintf('_clipY%d.png', clipyPairIdx)])
            saveas(gcf, [resfn sprintf('_clipY%d.pdf', clipyPairIdx)])
            colormap(chelixMap)
            saveas(gcf, [resfn sprintf('_clipY%d_cubeHelix.png', clipyPairIdx)])
            saveas(gcf, [resfn sprintf('_clipY%d_cubeHelix.pdf', clipyPairIdx)])
            
        end
    end
end


disp('done')
