%% Segment out a Monge-form "slab" [z0(x,y), z1(x,y)] from confocal data 
% NPMitchell 2020
%
% NOTE: THIS DATASET IS SINGLE COLOR
%
% This is a pipeline to segment confocal data & take MIPs

% temporary path def
% cd /mnt/data/confocal_data/gut/2020/mef2GAL4klarUASsqhGFPUASCAAXmCh/202010191649_mef2GAL4klarUASsqhGFPUASCAAXmCh_40x1p6x_240spf_0p4um_5to10pc3p5to7pc_oilImmersol/

addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
addpath_recurse('/mnt/data/code/gut_matlab/')
addpath('/mnt/data/code/TexturePatch_for_git/TexturePatch/')

dataDir = '/mnt/data/optogenetics_confocal/' ;
dataDir = [dataDir 'WTcontrol_48YG4kCAAXmCh/'] ;

% WT
dataDir = [dataDir '202106061730_48YCAAXmCh_0p75um_1p25x40x_2p5t5pc_lav3_615ns/'] ;

% noKlar
% dataDir = [dataDir '202103181730_antpG4OCRLGap43mCh_40x1p6x_5mpf_4pc3pc_to_12pc9pc_600ns_lav3_DC/'] ;

% klar
% dataDir = fullfile(dataDir, '202105252111_AntpG4kOCRLgap43_0p75um_1p25x40x_lav3_3t6pc3t6pc_5mpf_480ns_LED4') ;
% dataDir = fullfile(dataDir, '202105311838_AntpG4kOCRL_0p75um_1p5x40_3t6pc_lav3_86s_5mpf') ;
cd(dataDir)
gutDir = '/mnt/data/code/gut_matlab/' ;
addpath(fullfile(gutDir, 'addpath_recurse'))
addpath_recurse(gutDir)


% We start by clearing the memory and closing all figures
clear; close all; clc;

% Configuration metadata for this dataset
tp0 = 2 ;  % which timepoint (index) is the onset of folding?
embryoView = 'RD' ; % 'RD'/RV/LD/LV --> 
                     % which way is out of page (decreasing z) (R/L) 
                     % and up (increasing Y in Fiji) (D/V)

%% INITIALIZE ImSAnE PROJECT ==============================================
%
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.

dataDir    =  cd; 
projectDir = dataDir ;

% A filename base template - to be used throughout this script
% the 32 bit fn
fn = '' ;
% the 16 bit fn
%file16name = 'antpOCRLgap43_T%03d' ;     
file16name = 'wt_T%03d' ;     

% 1.5x40x with 800 pixels
pix2um = 0.242753316645 ;   % 06061730
% pix2um = 0.2431478 ;      % 05311838
% 1.5x40x with 1024 pixels
% pix2um = 0.18990684474 ;
% 1.25x40x with 1024 pixels
% pix2um = 0.22669667644183772 ;
um2pix = 1 / pix2um ;
resolution = [pix2um, pix2um, 0.75] ;
timepoints = 0:16 ;

%% Join data into stacks
data_is_split = false ;
if data_is_split
    fns0 = dir('./splitChannels/*C0*.tif') ;
    fns1 = dir('./splitChannels/*C1*.tif') ;
    for tidx = 1:length(fns0)
        fn0 = fullfile(fns0(tidx).folder, fns0(tidx).name) ;
        fn1 = fullfile(fns1(tidx).folder, fns1(tidx).name) ;
        fn2 = sprintf([file16name '.tif'], tidx) ;

        if ~exist(fn2, 'file')
            disp(['stacking file for tidx = ' num2str(tidx)])

            im0 = readTiff4D(fn0, 1) ;
            im1 = readTiff4D(fn1, 1) ;

            % Reorient the image
            if strcmpi(embryoView, 'rv')
                im0 = fliplr(im0) ;
                im0 = flipud(im0) ;
                im1 = fliplr(im1) ;
                im1 = flipud(im1) ;
            elseif strcmpi(embryoView, 'rd')
                disp('no rotation or flip necessary')
            else
                error('handle here')
            end

            im = cat(4, im0, im1) ;

            writeTiff5D(permute(im, [1, 2, 4, 3]), fn2)

            % % check it
            % for qq = 1:size(im0, 3)
            %   imagesc(squeeze(im0(:, :, qq))) ;
            %   pause(0.01)
            % end
        else
            disp('file exists')
        end
    end
    fn = file16name ;
else
    fn = file16name ;
end

%% CREATE EXPERIMENT
% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = project.Experiment(projectDir, dataDir);

% Set file and experiment meta data
% Set required additional information on the files.
% We assume on individual image stack for each time point, labeled by time.
%  To be able to load the stack, we need to tell the project wehre the data
%  is, what convention is assumed for the file names, available time
%  points, and the stack resolution.  Options for modules in ImSAnE are
%  organized in MATLAB structures, i.e a pair of field names and values are
%  provided for each option.
%
% The following file metadata information is requireL:
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of times available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
%
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th
              

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = 2;
fileMeta.timePoints         = timepoints;
fileMeta.stackResolution    = resolution ; % the px resolution (found in the .lif; 4 dec places)
fileMeta.swapZT             = 0;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is requireL:
%
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB = 123
% * 'dynamicSurface'    , Not implemented yet, future plan: boolean, false: static surface
% * 'detectorType'      , name of detector class, e.g. radielEdgeDetector
%                         ,(user threshholded), fastCylinderDetector
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
%
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = [1];
expMeta.channelColor        = 1;
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 0;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
% expMeta.detectorType = 'surfaceDetection.planarEdgeDetector';
% expMeta.fitterType = 'surfaceFitting.tpsFitter';
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';
expMeta.fitterType          = 'surfaceFitting.tpsFitter';

% %% 32 to 16 BIT CONVERSION
% % Check that data is 16 bit. If not, convert to 16bit
% pc2use = 99.9;
% scalemethod = 'user-defined' ; % use either 'user-defined' or 'prctile'
% scalemax = 0.05 ;
% 
% for channel_check = 1: fileMeta.timePoints(end) % last tp to convert to 16 bit
%     fullFileName = [sprintf(fn, channel_check) '.tiff'] ;
%     info = imfinfo(fullFileName) ;
%     full16fn = [sprintf(file16name, channel_check) '.tiff'] ;
%     bitDepth = info.BitDepth ;
% 
%     if (bitDepth == 32) && ~isfile(full16fn)
% 
%         disp([fullFileName ' is not 16bit, converting...'])
% 
%         % Note that imread only loads a single frame
%         % A = imread(fullFileName) ;
%         % scalemin = double(min(A(:))) ;
%         % scalemax = double(max(A(:))) ;
%         disp('')
%         disp('Reading 32 bit file to convert...')
%         A = readSingleTiff(fullFileName) ;
%         tmpA = A(:) ;
%         disp('')
%         disp('Computing scalemin, scalemax')
% 
%         % Optional step here to figure out what the cutoff
%         % intensity should be
%         % tmpA_no_ouliers = tmpA(tmpA < pcntile(tmpA, 99)) ;
%         % thisstd = std(tmpA_no_ouliers) ;
%         % check it using histogram(tmpA)
%         thismedian = median(tmpA) ;
% 
%         %goodmedian = 2559.00;
%         %worstmedian = 420.00;
%         %range2correct = goodmedian - worstmedian ;
%         %normal_pc2use = 99.9999 ;
%         %worstcase_pc2use = 99.99 ;
%         %diffpc = normal_pc2use - worstcase_pc2use ;
%         %pc2use = normal_pc2use + diffpc * (thismedian - goodmedian) / range2correct ;
%         %pc2use = max(worstcase_pc2use, pc2use) ;
%         %pc2use = min(normal_pc2use, pc2use) ;
%         chanpositionstart = strfind(fullFileName,'Ch');
%         chanposition = fullFileName(chanpositionstart+2);
%         chanposition = str2num(chanposition);
%         disp('determining prctile')
%         if strcmp('scalemethod', 'prctile')
%             scalemax = double(prctile(tmpA, pc2use)) ;
%         else
%             disp('using user supplied scalemax')
%         end
%             
%         scalemin = double(min(tmpA(tmpA > 0))) ;
%         disp(['scalemax = ', num2str(scalemax)])
%         disp(['scalemin = ', num2str(scalemin)])
% 
%         disp('Showing slice of scaled image')
%         % Note to self to check scale:
%         close all
%         imagesc(squeeze(A(:, 300, :)))
%         title('Checking scale of image')
%         waitfor(gcf)
%         histogram(A(:))
%         title('Histogram of values')
%         waitfor(gcf)
% 
%         % data = readSingleTiff(fullFileName);
%         im2 = mat2gray(A,[scalemin scalemax]);
%         im2 = uint16(2^16*im2);
%         imSize = size(im2);
%         
%         % Check scale:
%         imagesc(squeeze(im2(:, 300, :)))
%         title('Checking scale of image. Close image to continue')
%         axis equal
%         colorbar
%         waitfor(gcf)
% 
%         % Save the 16 bit image
%         disp(['Saving 16bit volume to ' full16fn])
%         imwrite(im2(:,:,1),full16fn,'tiff','Compression','none');
%         for z = 2 : imSize(3)
%             imwrite(im2(:,:,z),full16fn,'tiff','Compression','none','WriteMode','append');
%         end
%         disp('done saving 16bit volume')
% 
%         fn = file16name ;
% 
%     elseif isfile(full16fn)
%         % the data is already 16 bit, so we're good
%         fullFileName = [sprintf(fn, channel_check) '.tiff'] ;
%         disp([fullFileName ' has been converted.'])
% 
%         fn = file16name ;
%     else
%         disp('File is 16bit.')
%     end
% end

%% INSTANTIATE EXPERIMENT CLASS
% Now set the meta data in the experiment.
fileMeta.filenameFormat = [ file16name '.tif' ] ;
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%%%%%%%%%%%%%%%%%%%%%%%%%
% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%
mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
msls_axis_order = 'xyzc';
% Mesh marching options
normal_step = 10;

% Define the surface detection parameters
channel = 2;
foreGroundChannel = 2;
ssfactor = 4;
niter = 25 ;
niter0 = 115 ;
ofn_smoothply = 'mesh_' ;
ofn_ply = 'mesh_ms_' ; 
ofn_ls = 'msls_' ;
meshlabCodeDir = '/mnt/data/code/meshlab_codes/' ;
ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper' ;
lambda1 = 1 ;
lambda2 = 1 ;
exit_thres = 0.00001 ;
smoothing = 0.1 ;
nu = 0.0 ;
pre_nu = -5 ;
pre_smoothing = 1 ;
post_nu = 2;
post_smoothing = 4 ;
radius_guess = 10 ;
center_guess = 'empty_string' ;

% Name the output mesh directory ------------------------------------------
if projectDir(end) ~= filesep
    projectDir = [projectDir filesep];
end
mslsDir = fullfile(projectDir, 'msls_output');

%% LOAD THE FIRST TIME POINT ==============================================
%xp.loadTime(xp.fileMeta.timePoints(first_tp));
%xp.rescaleStackToUnitAspect();

%% DETECT THE SURFACE =====================================================
% Surface detection parameters --------------------------------------------
detectOptions = struct('channel', 1, ...
            'ssfactor', ssfactor,... % subsampling factor: downsampling of raw data
            'niter', 40, ... % how many iterations before exit if no convergence
            'niter0', 40, ... % how many iterations before exit if no convergence for first timepoint
            'lambda1', 1, ...  % lambda1/lambda2 decides weight of inclusion/exclusion of interior/exterior
            'lambda2', 1, ...  % lambda1/lambda2 decides weight of inclusion/exclusion of interior/exterior
            'nu', nu, ... % float: how many pressure (dilation/erosion) steps per iteration
            'smoothing', smoothing,... % float: how many smoothing steps per iteration (can be <1)
            'post_nu', post_nu, ... % how many iterations to dilate (if positive) or erode (if negative) after convergence
            'post_smoothing', post_smoothing,... % how many iterations of smoothing after convergence
            'exit_thres', 1e-6, ... % convergence thresholL: maximum difference between subsequent level sets upon which to exit algorithm ('close enough')
            'foreGroundChannel',foreGroundChannel, ... % the index of the first dimension of the 4d input data (if 4d)
            'fileName', sprintf( fn, xp.currentTime ), ... % the filename of h5 to train on
            'mslsDir', mslsDir, ...  % the directory for all output data/images
            'ofn_ls', ofn_ls, ...  % the output filename for level sets
            'ofn_ply', ofn_ply, ... % the output filename for PLY files
            'ms_scriptDir', ms_scriptDir, ... % the directory containing run_morphsnakes.py
            'timepoint', 0, ... % which timepoint in the data to consider
            'zdim', 3, ... % Which dimension is the z dimension
            'pre_nu', pre_nu, ... % number of dilation/erosion passes for positive/negative values
            'pre_smoothing', pre_smoothing, ... % number of smoothing passes before running MS
            'ofn_smoothply', ofn_smoothply,... % the output file name (not including path directory)
            'mlxprogram', fullfile(meshlabCodeDir, mlxprogram), ... % the name of the mlx program to use to smooth the results. Note that if mesh_from_pointcloud==true, should take obj as input and mesh as output.
            'init_ls_fn', 'mesh_initguess', ... % the name of the initial level set to load, if any
            'run_full_dataset', false, ... % run MS on a time series, not just one file
            'radius_guess', radius_guess, ... % radius of the initial guess sphere
            'dset_name', 'exported_data', ... % the name of the dataset to load from h5
            'save', true, ... % whether to save intermediate results
            'center_guess', 'empty_string', ... % xyz of the initial guess sphere ;
            'plot_mesh3d', false, ...  % if save is true, plot intermediate results in 3d 
            'dtype', 'h5', ... % h5 or npy: use hdf5 or numpy file format for input and output ls
            'mask', 'none', ... % filename for mask to apply before running MS
            'mesh_from_pointcloud', false, ... % use a pointcloud from the marching cubes algorithm rather than a mesh to create smoothed mesh
            'prob_searchstr', '_Probabilities.h5', ... % if dataset mode, what string to seek for loading all probabilities in data directory (glob datadir/*searchstr)
            'physicalaxisorder', 'yxzc', ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
            'preilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'ilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'include_boundary_faces', true, ... % keep faces along the boundaries of the data volume if true
            'smooth_with_matlab', -1, ... % if <0, use meshlab. If >0, smooth the mesh after marching cubes mesh creation using matlab instead of mlxprogram, with diffusion parameter lambda = this value. If =0, no smoothing.
            'pythonVersion', '') ; 
        
% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

% clear msls_exten imwriteOptions saveDir
% clear channel foreGroundChannel
% clear niter niter0 lambda1 lambda2
% clear exit_thres smoothing nu
% clear post_nu post_smoothing

%% Adjust LUT to be approx constant over time & Isotropic resolution
convert2IsotropicCombined = false ;
if convert2IsotropicCombined
    % Use Fiji BioFormats Exporter to export each channel and each timepoint to
    % a separate TIFF
    adjustlow_pctile = 0;
    adjusthigh_pctile = 99.5 ;
    for tt = xp.fileMeta.timePoints
        name_out = sprintf('./combined_T%02d.tif', tt) ;
        if ~exist(name_out, 'file')
            xp.loadTime(tt)
            xp.rescaleStackToUnitAspect();
            IV = xp.stack.image.apply() ;
            % Save a TIFF of the 2color data
            for ii = 1:length(IV)
                IVii = IV{ii} ;
                vlo = double(prctile( IVii(:) , adjustlow_pctile )) / double(max(IVii(:))) ;
                vhi = double(prctile( IVii(:) , adjusthigh_pctile )) / double(max(IVii(:))) ;
                disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                    ' for ', num2str(adjustlow), '/', num2str(adjusthigh)])
                IV{ii} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
            end
            IV1 = IV{1} ;
            IV2 = IV{2} ;
            im = uint16(zeros(size(IV1, 1), size(IV1,2), 2, size(IV1,3), 1)) ;
            im(:, :, 1, :, 1) = reshape(uint16(IV1), [size(IV1, 1), size(IV1, 2), 1, size(IV1, 3), 1]) ;
            im(:, :, 2, :, 1) = reshape(uint16(IV2), [size(IV1, 1), size(IV1, 2), 1, size(IV1, 3), 1]) ;
            dimensionOrder = 'XYZCT' ;   % must be in {'XYZCT', 'XYZTC', 'XYCTZ', 'XYCZT', 'XYTCZ', 'XYTZC'}           
            % bfsave(im, name_out, 'dimensionOrder', dimensionOrder, 'BigTiff', false)

            writeTiff5D(im, name_out)
        end
    end

    %% INSTANTIATE EXPERIMENT CLASS AGAIN WITH COMBINED STACKS
    % Now set the meta data in the experiment.
    fileMeta.filenameFormat = 'combined_T%02d.tif' ;
    fileMeta.stackResolution = min(fileMeta.stackResolution) * [1,1,1] ;
    xp.setFileMeta(fileMeta);
    xp.setExpMeta(expMeta);
    xp.initNew();
end


%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done
for t = xp.fileMeta.timePoints(1:end)
    if ~exist(fullfile(projectDir, [(sprintf(file16name, t)) '.h5']), 'file')
        disp(['Did not find file: ', fullfile(projectDir, [sprintf(file16name, t) '.h5'])])
        
        % Only load and rescale if multiple timepoints/channels
        % if length(xp.fileMeta.timePoints) > 1 && t > xp.fileMeta.timePoints(1)
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();
        % end
        
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.fileName = sprintf( file16name, xp.currentTime ) ;
        xp.setDetectOptions( detectOpts2 );
        xp.detector.prepareIlastik(xp.stack);
        disp(['done outputting downsampled data h5: tp=' num2str(t) ' for surface detection'])
    else
        disp(['h5 ' num2str(t) ' was already output, skipping...'])
    end
end
disp('Open with ilastik if not already done')

%% TRAIN DATA IN ILASTIK TO IDENTIFY SURFACE ==========================
% % open ilastik, train until probabilities and uncertainty are satisfactory
% 

%% Create MorphSnakesLevelSet from the Probabilities from ilastik ========
% Skip if already done
% Now detect all surfaces
run_full_dataset_ms = false ;
detectOptions.run_full_dataset = 'none' ;  % projectDir ; % 'none' ;  % override here
if strcmp(detectOptions.run_full_dataset, projectDir)
    % assert(run_full_dataset_ms)
    disp('Running dataset mode')
    xp.setTime(xp.fileMeta.timePoints(1));
    detectOpts2 = detectOptions ;
    detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
    detectOpts2.nu = 4 ;
    detectOpts2.niter0 = 5 ;
    xp.setDetectOptions( detectOpts2 );
    xp.detectSurface();
else
    assert(~run_full_dataset_ms)
    assert(strcmp(detectOptions.run_full_dataset, 'none'))
    % Morphosnakes for all remaining timepoints INDIVIDUALLY ==============
    for tp = xp.fileMeta.timePoints
        % try
            xp.setTime(tp);
            % xp.loadTime(tp) ;
            % xp.rescaleStackToUnitAspect();

            % make a copy of the detectOptions and change the fileName
            detectOpts2 = detectOptions ;
            % detectOpts2.post_smoothing = 1 ;
            detectOpts2.timepoint = xp.currentTime ;
            detectOpts2.fileName = sprintf( fn, xp.currentTime );
            % detectOpts2.mlxprogram = fullfile(meshlabCodeDir, ...
            %      'surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4') ;
            %  _octree12.mlx') ;
            % detectOpts2.mlxprogram = fullfile(meshlabCodeDir, ...
            %     'laplace_surface_rm_resample25k_reconstruct_LS3_1p2pc_ssfactor4_vip8test.mlx') ;
            detectOpts2.mlxprogram = fullfile(meshlabCodeDir, ...
                 'laplace_surface_rm_resample25k_reconstruct_LS3_wu13_ssfactor4.mlx') ;
            xp.setDetectOptions( detectOpts2 );
            xp.detectSurface();
            % For next time, use the output mesh as an initial mesh
            detectOpts2.init_ls_fn = 'none' ;
        % catch
        %    error('here')
        %     disp('Could not create mesh -- skipping for now')
        %     % On next timepoint, use the tp previous to current time
        %     detectOptions.init_ls_fn = [detectOptions.ofn_ls, ...
        %             num2str(tp - 1, '%06d' ) '.' detectOptions.dtype] ;
        % end
    end
end


%%
xp.loadTime(timepoints(1)) ;

%% Load mesh, subsample, triangulate, Texturepatch
IV = xp.stack.image.apply() ;
sz1 = size(IV{1}, 1) ;
sz2 = size(IV{1}, 2) ;
sz3 = size(IV{1}, 3) ;
nU = round(150 * sz1 / sz2) ;
nV = 150 ;
uminmax = [0, sz1] ;
vminmax = [0, sz2] ;
%%
layerWidth = 5 ; % 1 for RFP
subsample = 10 ;
zOffset = 5 ; % 20 ;
gaussSigma = 0.5 ; % 2;
axisOrder = 'zyx' ;

axisOrder = erase(axisOrder, 'c') ;
zdim = find(axisOrder== 'z') ;


%% Take median z(t) for each (x,y) and remake stacks
preview = true ;
overwrite = true ;
filteredMeshFn = './msls_output/medianFilteredMeshes.mat' ;
if ~exist(filteredMeshFn, 'file') || overwrite 
    prevMesh0 = [] ;
    prevMesh1 = [] ;
    xp.loadTime(timepoints(1));
    xp.rescaleStackToUnitAspect();
    IV = xp.stack.image.apply() ;
    sz1 = size(IV{1}, 1) ;
    sz2 = size(IV{1}, 2) ;
    sz3 = size(IV{1}, 3) ;
    if contains(lower(axisOrder), 'xyz')
        szX = sz1 ; szY = sz2 ; szZ = sz3 ;
        zdim = 3 ;
    elseif contains(lower(axisOrder), 'zyx') 
        szX = sz3 ; szY = sz2 ; szZ = sz1 ;
        zdim = 1 ;
    else
        error('handle here')
    end
    %% Consider each mesh, resample. OPTIONAL: filter with adjacent timepoints 
    for tidx = 1:length(timepoints)
        tp = timepoints(tidx) ;
        disp(['t = ' num2str(tp)])
        meshfn = fullfile(mslsDir, sprintf('mesh_ms_%06d.ply', tp)) ;
        mesh = read_ply_mod(meshfn) ;
        % scatter3(mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), mesh.v(:, 1))

        % Floor is max(x)
        plane = find(mesh.v(:, zdim) > max(mesh.v(:, zdim))-0.1) ;
        % surfId = setdiff(1:length(mesh.v), plane) ;

        % % Sample some of the planar points
        % p4bnd = plane(1:10:end) ;
        % % Keep only boundary of planar points
        % p4bndV = mesh.v(p4bnd, :) ;
        % p4bndF = delaunay(p4bndV(:, 2), p4bndV(:, 3)) ;
        % tri = triangulation(p4bndF, p4bndV)  ;
        % bnd = freeBoundary(tri) ;
        % bndV = p4bndV(bnd(:, 1), :) ;

        % check it
        if preview
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
            axis equal
            xlabel('x'); ylabel('y'); zlabel('z') ;
            pause(0.1);
        end
        
        % Subsample the surface 
        [face, vertex] = remove_vertex_from_mesh(mesh.f, mesh.v, plane) ;
        [ surfF, surfV, oldVertexIDx ] = ...
             remove_isolated_mesh_components( face, vertex ) ;

        % check it
        if preview
            trisurf(triangulation(surfF, surfV), surfV(:, 3), 'edgecolor', 'none')
            axis equal;
            view(2)
            xlabel('x'); ylabel('y'); zlabel('z') ;
            pause(0.1);
        end

        % surfV = mesh.v(surfId, :) ;
        if contains(lower(axisOrder), 'zyx')
            uvz = [surfV(:, 3), surfV(:, 2), surfV(:, 1)] ;
        elseif contains(lower(axisOrder), 'xyz')
            if strcmpi(embryoView, 'rd')        
                uvz = [surfV(:, 1), surfV(:, 2), surfV(:, 3)] ;
            else
                error('handle here')
            end

            if preview
                scatter3(uvz(:, 1), uvz(:, 2), uvz(:, 3), 5, uvz(:, 3))
                axis equal
                xlabel('x'); ylabel('y'); zlabel('z') ;
            end
        end

        % subsample by binning and finding medians
        [zmeans, counts, zs, xidx, yidx] = ...
            binData2dGrid(uvz, uminmax, vminmax, nU, nV, false) ;
        medz = zeros(nU, nV) ;
        for pp = 1:size(zs, 1)
            for qq = 1:size(zs, 2)
                medz(pp, qq) = median(zs{pp, qq}) - zOffset ;
            end
        end
        Z0 = sz3; 
        out_of_frame = find(isnan(medz)) ;
        medz(out_of_frame) = Z0 ;
        % local minimum filter (lower z is further from midsagittal plane)
        medz = imerode(medz, true(5)) ;
        % medz = imgaussfilt(medz, gaussSigma) ;
        [yy, xx] = meshgrid(linspace(1, sz1, nU), linspace(1, sz2, nV)) ;
        xx = xx';
        yy = yy';

        % NOTE: This should look correct
        if preview
            imagesc( linspace(1, sz2, nV), linspace(1, sz1, nU), medz)
            axis equal ;
            title('rescaled resampled surface')
            pause(1)
        end    

        % Create ring of zeros around shape
        % se = strel('disk', 2);
        % keep = find(imdilate(~isnan(zmeans), se)) ;
        % assert(all(size(medz) == size(zmeans)))

        xr = xx(:) ;
        yr = yy(:) ;
        zr = medz(:) ;
        scatter3(xr, yr, zr, 4, zr); view(2); axis equal

        % xr = xr(keep) ;
        % yr = yr(keep) ;
        % zr = zr(keep) ;
        % xr = [xr; 0; sz1; sz1; 0] ;
        % yr = [yr; 0; 0; sz2; sz2] ;
        % zr = [zr; Z0; Z0; Z0; Z0] ;

        % Inspect cross-section
        if preview && false
            % Make this number larger to sample more of the nearby mesh
            width = 4 ;
            leaves = 1:50:szY ;
            % Show the cross-section
            for leaf = leaves
                inds = find(abs(yy - leaf) < width) ;
                clf
                if strcmpi(erase(axisOrder, 'c'), 'xyz')
                    im = squeeze(IV{2}(:, leaf, :))' ;
                elseif strcmpi(erase(axisOrder, 'c'), 'xyz')
                    im = squeeze(IV{2}(leaf, :, :))' ;
                end
                imshow(mat2gray(im, [0, double(rms1d(im(:)))]))
                hold on; 
                if any(inds)
                    hold on;
                    plot(xx(inds), medz(inds), 'c.')
                end
                pause(0.1)
            end
        end

        % Triangulate the result
        faces = delaunay(xr, yr) ;
        mesh = struct() ;
        mesh.f = faces ;
        mesh.v = [xr, yr, zr] ;

        % SMOOTH IN TIME (SKIPPING THIS SINCE MESHES ARE CLEAN)
        % if tidx == 1
        %     % Do not define mesh until next timepoint
        % elseif tidx == 2
        %     % First mesh is average of first two unsmoothed meshes
        %     meshFilt{tidx-1} = mesh.v ;
        %     assert(all(mesh.v(:, 1) == prevMesh1.v(:, 1))) ;
        %     assert(all(mesh.v(:, 2) == prevMesh1.v(:, 2))) ;
        %     meshFilt{tidx-1}(:, 3) = mean( [mesh.v(:, 3), prevMesh1.v(:, 3)], 2) ;
        % else
        %     % Take median of previous two and current meshes
        %     meshFilt{tidx-1} = mesh.v ;
        %     assert(all(mesh.v(:, 1) == prevMesh1.v(:, 1))) ;
        %     assert(all(mesh.v(:, 2) == prevMesh1.v(:, 2))) ;
        %     meshFilt{tidx-1}(:, 3) = mean( [mesh.v(:, 3),  ...
        %         prevMesh0.v(:, 3), prevMesh1.v(:, 3)], 2) ;
        % end
        
        meshFilt{tidx} = mesh.v ;
        out_of_frame_idx{tidx} = out_of_frame ;

        % Prepare for next timepoint
        prevMesh0 = prevMesh1 ;
        prevMesh1 = mesh ;

        if preview && tidx > 1
            clf
            tri = triangulation(mesh.f, meshFilt{tidx-1}) ;
            trisurf(tri, 'edgecolor', 'none')
            view(2); title('triangulation of rescaled resampling')
            axis equal
            pause(1) ;
        end
    end
    save(filteredMeshFn, 'meshFilt', 'faces', 'out_of_frame_idx')
else
    load(filteredMeshFn, 'meshFilt', 'faces', 'out_of_frame_idx')    
end


%% Onion stacks
nPos = 25 ;
nNeg = 10 ;
dataAxisOrder = 'xyz' ;

layerWidth = 4 ;
midLayerOffset = 7 ;
lam = 0.01 ;  % 1 ; 

preview = false ;
preview_final = true ;
overwrite = true ;

timepoints = 0:16 ;
tidx2do = 1:5:length(timepoints) ;
tidx2do = [tidx2do, setdiff(1:length(timepoints), tidx2do)] ;

% Load all filtered meshes
load('./msls_output/medianFilteredMeshes.mat', 'meshFilt', 'faces', 'out_of_frame_idx')
for tidx = tidx2do
    tp = timepoints(tidx) ;
    imfn = sprintf('./texturePatches/slice_T%03d_c%01d.tif', tp, 2) ;
    if ~exist(imfn, 'file') || overwrite
        xp.loadTime(tp);
        xp.rescaleStackToUnitAspect();
        IV = xp.stack.image.apply() ;
        sz1 = size(IV{1}, 1) ;
        sz2 = size(IV{1}, 2) ;
        sz3 = size(IV{1}, 3) ;
        if contains(lower(dataAxisOrder), 'xyz')
            szX = sz1 ; szY = sz2 ; szZ = sz3 ;
            zdim = 3 ;
        elseif contains(lower(dataAxisOrder), 'xyz')
            szX = sz3 ; szY = sz2 ; szZ = sz1 ;
            zdim = 1 ;
        else
            error('handle here')
        end
        
        % Could load original mesh
        % meshfn = fullfile(mslsDir, sprintf('mesh_ms_%06d.ply', tp)) ;
        % mesh = read_ply_mod(meshfn) ;
        
        % Instead use filtered mesh
        mesh = struct() ;
        mesh.f = faces ;
        mesh.v = meshFilt{tidx} ;
        out_of_frame = out_of_frame_idx{tidx} ;         
        
        % Fix vertices that are off mesh or near the mesh boundary
        bndId = neighbors(tri, out_of_frame); % find(zr > max(medz(:))-5)) ;
        bndId = unique(bndId(~isnan(bndId))) ;
        bndId = tri.ConnectivityList(bndId, :) ;
        bndId = unique(bndId(:)) ;
        vsm0 = mesh.v ;
        mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', bndId, lam, ...
            'implicit', mesh.v, 1000);

        if preview
            clf
            tri = triangulation(mesh.f, mesh.v) ;
            trisurf(tri, mesh.v(:, 3), 'edgecolor', 'none')
            view(2); title('smoothed triangulation of rescaled resampling')
            axis equal
            xlabel('x'); ylabel('y'); zlabel('z') ;
            pause(0.1)
        end
        
        m2d = mesh ;
        m2d.v = [xr, yr] ;
        mesh.vn = zeros(size(mesh.v)) ;
        mesh.vn(:, zdim) = 1 ;
        % mesh.vn(out_of_frame, 3) = 0 ;
        mesh.v = mesh.v ;
        if preview
            trisurf(triangulation(mesh.f, mesh.v), mesh.vn(:, 3), 'edgecolor', 'none')
            title('normal vectors')
            axis equal
            caxis([-1, 1])
            xlabel('x'); ylabel('y'); zlabel('z') ;
            pause(0.1)
        end

        % TexturePatch
        Opts = struct() ;
        if strcmpi(axisOrder, 'xyz')
            Opts.imSize = [sz2, sz1] ;
        elseif strcmpi(axisOrder, 'zyx')
            Opts.imSize = [sz1, sz2] ;
        else
            error('handle here')
        end
        Opts.numLayers = [nPos, nNeg] ;
        Opts.layerSpacing = 1 ;
        Opts.vertexNormal = mesh.vn ;
        Opts.extrapolationMethod = 'nearest' ; % 'nearest' ;'none'; 
        if contains(lower(axisOrder), 'zyx')
            [ patchIm, imref, zeroID, MIP, SIP ] = ...
                texture_patch_to_image(mesh.f, m2d.v, mesh.f, ...
                mesh.v(:, [2,1,3]), IV, Opts) ;
        elseif contains(lower(axisOrder), 'xyz')
            IVtmp = IV ;
            IVtmp{1} = permute(IVtmp{1}, [2,1,3]) ;
            IVtmp{2} = permute(IVtmp{2}, [2,1,3]) ;
            [ patchIm, imref, zeroID, MIP, SIP ] = ...
                texture_patch_to_image(mesh.f, m2d.v(:, [2,1]), mesh.f, ...
                mesh.v(:, [1,2,3]), IVtmp, Opts) ;            
        else    
            error('handle here')
        end
        if preview
            imshow(permute(squeeze(max(patchIm, [], 3)), [2, 1, 3]))
        end

        % texture_patch_3d(mesh.f, mesh.v, mesh.f, mesh.v(:, [2,1,3]), IV{2})

        if preview_final
            figure(1); clf
            for ii = 1:size(patchIm, 3)
                imshow(mat2gray(squeeze(patchIm(:, :, ii)), [0, 0.25])')
                title(['ii= ' num2str(ii)])
                pause(0.01)
            end
        end
        
        % Look at figure
        if preview
            figure(2)
            trisurf(triangulation(mesh.f, mesh.v), mesh.v(:, 3), 'edgecolor', 'none')
            view(2); axis equal
            pause(0.01)
        end
        
        % imshow(patchIm(:,:,2, zeroID))
        % caxis([0, 0.2])

        if ~exist('./texturePatches/', 'dir')
            mkdir('./texturePatches/')
        end
        ch1 = patchIm ;
        
        % % Preview cross-section of stack
        % for ii = 100:1:450 % size(ch2, 2)
        %     imshow(mat2gray(squeeze(ch2(ii, :, :)), [0, 0.1 * max(ch2(:))])')
        %     title(['ii= ' num2str(ii)])
        %     pause(0.1)
        % end

        ch1 = uint8(mat2gray(ch1, [min(ch1(:)), max(ch1(:))]) * 2^8) ;
        % ch1 = permute(ch1, [2, 1, 3]) ;
        chs = {ch1} ;

        % Save image stack
        for ch = 1
            outfn = sprintf('./texturePatches/slice_T%03d_c%01d.tif', tp, ch) ;
            writeTiff5D(chs{ch}, outfn, 8) ;
        end

        % Save mesh
        save(sprintf('./texturePatches/mesh_T%03d.mat', tp), 'mesh', 'm2d', 'Opts') ;
    else
        
        disp(['already on disk: t=' num2str(tp)])
    end
    
    outfn1 = sprintf('./texturePatches/layer0_T%03d_c%01d.png', tp, 1) ;
    if ~exist(outfn1, 'file') || overwrite
        imfn1 = sprintf('./texturePatches/slice_T%03d_c%01d.tif', tp, 1) ;
        ims1 = loadtiff(imfn1) ;
        minLayer = max(nNeg-layerWidth+midLayerOffset, 1) ;
        maxLayer = min(nNeg+layerWidth+midLayerOffset, size(ims1, 3)) ;
        im01 = squeeze(max(ims1(:, :, minLayer:maxLayer), [], 3)) ;
        im01 = uint8(255 * mat2gray(im01, double([min(im01(:)), max(im01(:))]))) ;
        imwrite(im01, outfn1)
    end
end

%% Open images in iLastik, train on membrane for segmentation


%% Segment the cells by selecting watershed seed points 
overwrite = false ;
preThres = 0.3 ;
cellSize = 20;
strelRadius = 0;
gaussKernel = 0 ;
heightMiminum = 1;
useFilteredMesh = true ;

mkdir('./cellSegmentation/')
tidx2do = 1:27 ;
for tidx = tidx2do
    tp = timepoints(tidx) ;
    disp(['tp = ', num2str(tp)])
    segfn = sprintf('./cellSegmentation/automask_T%03d.png', tp) ;
    % segfn2 = sprintf('./cellSegmentation/automask2_T%03d.png', tp) ;
    % segfn3 = sprintf('./cellSegmentation/automask3_T%03d.png', tp) ;
    
    if ~exist(segfn, 'file') || overwrite
        %% Load image ilastik probabilities
        imfn = sprintf('./texturePatches/layer0_T%03d_c1_Probabilities.h5', tp) ;
        prob = h5read(imfn, '/exported_data/') ;
        im2segment = squeeze(prob(1, :, :)) ;
        
        %% Use raw
        % imfn = sprintf('./texturePatches/rgb/layer0_T%03d_rgb.png', tp) ;
        % im = imread(imfn) ;
        % im2segment = squeeze(im(:, :, 2)) ;

        % NEW watershed segmentation
        adaphisteqClip = 0.5 ;
        strelRadius = 1 ;
        [skel, DL] = segmentImageWatershed(im2segment, adaphisteqClip, strelRadius) ;
        
        % skel = segmentImageWatershedSimple(img, )
        
        % Mask the segmentation by depth of surface
        if useFilteredMesh
            surface = meshFilt{tidx} ;
        else
            load(sprintf('./texturePatches/mesh_T%03d.mat', tp), 'mesh') ;
            surface = mesh.v ;
        end
        % resample surface onto image grid
        [xx,yy] = meshgrid(1:size(DL, 1), 1:size(DL, 2)) ;
        resurf = interpolate2Dpts_3Dmesh(mesh.f, surface(:, 1:2), surface, [xx(:), yy(:)]) ;
        maxZ = max(resurf(:, 3)) ;
        keep = reshape(resurf(:, 3) < (maxZ - 1), [size(DL, 2), size(DL, 1)]) ;
        
        mask1 = skel' .* keep ;
        imshow(mask1) ; pause(0.001)
        
        % Save the result
        imwrite(mask1, segfn)
        
        disp(['Done with seg ' num2str(tp)])
    end
end



%% Load in GIMP and export binary masks

%% Read binary masks
thres = 160 ; % 160
timeInterval = 5 ;
timeUnits = 'min' ;
overwrite = false ;
colors = define_colors ;
imDir = './cellSegmentation/images/' ;
imDir_seg = fullfile(imDir, 'segmentation') ;
imDir_bnd = fullfile(imDir, 'boundaries2d') ;
imDir_bnd3d = fullfile(imDir, 'boundaries3d') ;
segDir = './cellSegmentation/' ;
textureDir = './texturePatches/'; 
if ~exist(imDir, 'dir')
    mkdir(imDir)
end
if ~exist(imDir_seg, 'dir')
    mkdir(imDir_seg)
end
if ~exist(imDir_bnd, 'dir')
    mkdir(imDir_bnd)
end
if ~exist(imDir_bnd3d, 'dir')
    mkdir(imDir_bnd3d)
end

tidx2do = 1:18 ;

for tidx = tidx2do
    tp = timepoints(tidx) ;
    outfn = fullfile(segDir, sprintf('T%03d_polygons3d.mat', tp)) ; 
    
    if ~exist(outfn, 'file') || overwrite
        disp(['Processing t=' num2str(tp) ': ' outfn])
        % try
            % segfn = fullfile(segDir, sprintf('T%03dmask.png', tp)) ;
            segfn = fullfile(segDir, sprintf('automask1_T%03d.png', tp)) ;
            seg = imread(segfn) ;
            load(fullfile(textureDir, sprintf('mesh_T%03d.mat', tp)), ...
                'mesh', 'm2d', 'Opts') ;
            imfn = fullfile(textureDir, sprintf('slice_T%03d_c%01d.tif', tp, 2)) ;
            im = loadtiff(imfn) ;

            % Threshold for initial segmentation
            bw = seg > thres ;

            % Convert segmentation into polygon ids
            se1 = strel('disk', 1) ;
            
            % If linewidth was 1 or 2, use this:
            % se2 = strel('disk', 5) ;
            % bw2 = imdilate(bw, se2) ;
            
            bw2 = bwskel(bw) ; 
            
            % Check it 
            clf
            bwcheck = imoverlay(0.5*bw, bw2, 'r') ;
            imshow(bwcheck)
            pause(1)

            % Get all boundaries    
            CC = bwconncomp(~bw2, 4);
            segIm = labelmatrix(CC) ;
            c2d = cell(max(segIm(:)), 1) ;
            props = regionprops(segIm, 'centroid') ;
            
            % Load single page image
            singlePage = imread( ...
                fullfile(textureDir, sprintf('layer0_T%03d_c2.png', tp))) ;

            % cmpQ = brewermap(max(segIm(:)), 'Paired') ;
            nCells = double(max(segIm(:))) ;
            modulation = repmat([-1,-1,-1; 1,1,1], ceil(nCells * 0.5), 1) ;
            cmpQ = phasemap(double(max(segIm(:)))) + 0.1 * modulation(1:nCells, :) ;
            cmpQ(cmpQ > 1) = 1;
            cmpQ(cmpQ < 0) = 0;
            % cmpQ = hsv(double(max(segIm(:)))) ;
            coloredLabelsImage = label2rgb(segIm, cmpQ, 'k', 'shuffle');
            % Display the pseudo-colored image.
            coloredLabelsImage = (coloredLabelsImage) .* uint8(segIm > 1) ;
            outim = uint8(0.5 * coloredLabelsImage + singlePage) ;
            imfn = fullfile(imDir_seg, sprintf('T%03d_segcolor.png', tp)) ;
            imwrite(outim, imfn)
            clf
            imshow(outim) ;
            centroids2d = zeros(max(segIm(:)), 2) ;
            for pId = 2:max(segIm(:))    
                centroids2d(pId, :) = props(pId).Centroid ;
                pgon = false(size(bw)) ;
                pgon(segIm == pId) = true ;
                pgon = imdilate(pgon, se1) ;
                tmp = bwboundaries(pgon, 8);

                hold on;
                plot(tmp{1}(:, 2), tmp{1}(:, 1), '.-', 'color', colors(2, :)) 
                plot(centroids2d(:, 1), centroids2d(:, 2), '.', 'color', colors(3, :))
                % pause(0.05) ;

                bnd = [tmp{1}(:, 2), tmp{1}(:, 1)] ;

                c2d{pId} = bnd ;
            end 
            imfn = fullfile(imDir_bnd, sprintf('T%03d_boundaries.png', tp)) ;
            title(['t = ' num2str((tp - tp0) * timeInterval) ' ' timeUnits])
            saveas(gcf, imfn)

            % Extract polygon shapes
            % Get z position for each cell polygon
            [c3d, centroids3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d, fieldfaces] = ...
                polygon3dMeasurements(m2d.f, mesh.v, m2d.v, c2d, centroids2d) ;

            % % Check this 
            % Color by aspect ratio
            clf
            aratioImage = double(0*segIm) ;
            qstrength = sqrt(moment2 ./ moment1) -1 ;
            for cid = 1:length(c2d)
                if ~isnan(moment2(cid))
                    aratioImage(segIm == cid) = qstrength(cid) * cos(2*ang1(cid)) ;
                end
            end        
            aratioImage(segIm < 2) = NaN ;
            [xL, yL] = size(aratioImage) ;
            imagesc(1:yL*pix2um, 1:xL*pix2um, aratioImage); 
            axis equal
            grid off
            cb = colorbar ;
            xlabel('ap position, [$\mu$m]', 'interpreter', 'latex')
            ylabel('dv position, [$\mu$m]', 'interpreter', 'latex')
            caxis(max(abs(aratioImage(:))) * [-1,1])
            colormap(bwr)
            ylabel(cb, '$Q_{xx}$', 'interpreter', 'latex')
            title(['t = ' num2str((tp - tp0) * timeInterval) ' ' timeUnits], ...
                'interpreter', 'latex')
            imfn = fullfile(imDir_bnd3d, sprintf('T%03d_Qxx3d.pdf', tp)) ;
            saveas(gcf, imfn)
            imfn = fullfile(imDir_bnd3d, sprintf('T%03d_Qxx3d.png', tp)) ;
            saveas(gcf, imfn)

            
            % Check it against 2d
            % nPolygons = length(c2d) ;
            % areas2d = zeros(nPolygons, 1) ;
            % perim2d = zeros(nPolygons, 1) ;
            % moment1_2d = zeros(nPolygons, 1) ;
            % moment2_2d = zeros(nPolygons, 1) ;
            % ang1_2d = zeros(nPolygons, 1) ;
            % ang2_2d = zeros(nPolygons, 1) ;
            % moinertia_2d = zeros(nPolygons, 2, 2) ;
            % for cid = 1:nPolygons
            %     pgon = c2d{cid} ;
            %     if ~isempty(pgon)
            %         [ geom, iner, cpmo ] = polygeom( pgon(:, 1), ...
            %             pgon(:, 2) ) ;
            %         areas2d(cid) = geom(1) ;
            %         perim2d(cid) = geom(4) ;
            %         moment1_2d(cid) = cpmo(1) ;
            %         ang1_2d(cid) = cpmo(2) ;
            %         moment2_2d(cid) = cpmo(3) ;
            %         moinertia_2d(cid, :, :) = [iner(4) -iner(6); -iner(6) iner(5)] ;
            %     end
            % end
            % qstrength2d = sqrt(moment2_2d ./ moment1_2d) - 1 ;
            % subplot(1, 2, 1)
            % plot(qstrength2d, (qstrength- qstrength2d(:)) ./ qstrength(:), '.')
            % subplot(1, 2, 2)
            % plot(areas2d, (areas- areas2d(:)) ./ areas2d, '.')
            
            % View 3d polygons
            clf
            minzum = min(mesh.v (:, 3) * pix2um) ;
            trisurf(triangulation(mesh.f, ...
                mesh.v * pix2um - [0, 0, minzum]), ...
                'edgecolor', 'none')
            for pId = 2:max(segIm(:))
                hold on;
                plot3(c3d{pId}(:, 1) * pix2um, ...
                    c3d{pId}(:, 2) * pix2um, ...
                    c3d{pId}(:, 3) * pix2um - minzum, '-')
            end
            view(2)
            axis equal
            colormap gray
            grid off
            xlabel('ap position, [$\mu$m]', 'interpreter', 'latex')
            ylabel('dv position, [$\mu$m]', 'interpreter', 'latex')
            cb = colorbar ;
            ylabel(cb, '$z$ depth, [$\mu$m]', 'interpreter', 'latex')
            imfn = fullfile(imDir_bnd3d, sprintf('T%03d_boundaries3d.pdf', tp)) ;
            title(['t = ' num2str((tp - tp0) * timeInterval) ' ' timeUnits], ...
                'interpreter', 'latex')
            saveas(gcf, imfn)
            imfn = fullfile(imDir_bnd3d, sprintf('T%03d_boundaries3d.png', tp)) ;
            saveas(gcf, imfn)

            disp(['Saving measurement to ', outfn])
            readme = struct() ;
            readme.centroids2d = '2d evaluation coordinates for tangent plane on embedding -- could/should be the 2d polygon centroids' ;
            readme.centroids3d = 'centroids2d in push-forward / embedding space' ;
            readme.areas = 'area of each polygon in tangent plane';
            readme.perim = 'perimeter of each polygon in tangent plane';
            readme.moment1 = 'maximum moment of inertia (extremal #1)' ;
            readme.ang1 = 'angle of eigenvector of moment of inertia tensor of polygon in 3d in a frame rotated to the 2d pullback coordinate frame' ;
            readme.moment2 = 'minimum moment of inertia (extremal #2)';
            readme.moinertia = 'embeddingspace moment of inertia of pushed forward polygon in basis aligned with pullback' ;
            readme.cellQ2d = 'polygon as quasi-2d' ;
            readme.fieldfaces = 'indices of faces on 3d mesh of evaluation points (centroids)';
            save(outfn, 'c3d', 'centroids2d', 'centroids3d', ...
                'areas', 'perim', 'moment1', 'ang1', 'moment2', 'ang2', ...
                'moinertia', 'cellQ2d', 'fieldfaces', 'readme') ;
        % catch
        %     disp(['Skipping tp= ' num2str(tp)])
        % end
            results_in_memory = false ;
    else
        results_in_memory = false ;
    end    
end

%% Identify lobe2-3 divide via surfaces
% for now just draw it
tidx2do = 1:length(timepoints) ;
% fold2 = 665 * ones(length(timepoints), 1) ;
fold2 = 575 * ones(length(timepoints), 1) ;
time_offset = 0 * timepoints ;
time_offset(2) = 3 ;
time_offset(3:end) = 6 ;
t0 = tp0 * timeInterval - time_offset(tp0) ;

%% plot Q as function of ap position
clf
exciteIdx = 2 ;
% cmap0 = cividis() ;
% cmap0 = cmap0(round(linspace(1, length(cmap0), length(tidx2do))), :) ;
cmap0 = viridis(length(timepoints)) ;
% Load an im for size
imfn = fullfile(textureDir, sprintf('slice_T%03d_c%01d.tif', timepoints(1), 2)) ;
im = loadtiff(imfn) ;
xe = linspace(0, size(im, 2), 20) ;

Qct_a = nan(length(tidx2do), 2) ;
Qct_s = nan(length(tidx2do), 2) ;
Qct_e = nan(length(tidx2do), 2) ;
Qst_a = nan(length(tidx2do), 2) ;
Qst_s = nan(length(tidx2do), 2) ;
Qst_e = nan(length(tidx2do), 2) ; 
for tidx = tidx2do
    tp = timepoints(tidx) ;
    outfn = fullfile(segDir, sprintf('T%03d_polygons3d.mat', tp)) ; 
    
    % What is the division between lobes12 and 34 for this timepoint?
    lobe12 = [0, fold2(tidx), size(im, 2)];
    
    if exist(outfn, 'file') 
        disp(['loading t=' num2str(tp) ': ' outfn])
        tmp = load(outfn) ;
        cx = tmp.centroids2d(:, 1) ;
        weights = tmp.areas ;
        m1 = abs(tmp.moment1) ;
        m2 = abs(tmp.moment2) ;
        a1 = tmp.ang1 ;
        qcos2 = (sqrt(m2./m1)-1) .* cos(2 * a1) ;
        % qcos2 = -(sqrt(m2./m1)-1) ; % .* cos(2 * a1) ;
        qsin2 = (sqrt(m2./m1)-1) .* sin(2 * a1) ;
        [midx, meanC, stdC, nC] = binDataMeanStdWeighted(cx, qcos2, xe, weights) ;
        [~, meanS, stdS, nS] = binDataMeanStdWeighted(cx, qsin2, xe, weights) ;
        
        [~, meanC_L12, stdC_L12, nC_L12] = binDataMeanStdWeighted(cx, qcos2, lobe12, weights) ;
        [~, meanS_L12, stdS_L12, nS_L12] = binDataMeanStdWeighted(cx, qsin2, lobe12, weights) ;
        
        assert(all(nC == nS)) ;
        steC = stdC ./ max(1, sqrt(nC-1)) ;
        steS = stdS ./ max(1, sqrt(nS-1)) ;
        steC_L12 = stdC_L12 ./ max(1, sqrt(nC_L12-1)) ;
        steS_L12 = stdS_L12 ./ max(1, sqrt(nS_L12-1)) ;
        
        lineProps = {'.', 'color', cmap0(tidx, :)} ;
        subplot(1, 2, 1)
        h1=shadedErrorBar(midx, meanC, steC, 'lineProps', lineProps) ;
        subplot(1, 2, 2)
        h2=shadedErrorBar(midx, meanS, steS, 'lineProps', lineProps) ;
        hold on;
        
        statStruct = struct() ;
        % Raw data sub-structure
        statStruct.raw = struct() ;
        statStruct.raw.centroids2d = cx ;
        statStruct.raw.areas = weights ;
        statStruct.raw.moment1 = m1 ;
        statStruct.raw.moment2 = m2 ;
        statStruct.raw.angle = tmp.ang1 ;
        % Now ap statistics and lobe12 specifically
        statStruct.ap = struct() ;
        % all cells
        statStruct.ap.midx = midx ;
        statStruct.ap.Qcos2t_avg = meanC ;
        statStruct.ap.Qcos2t_std = stdC ;
        statStruct.ap.Qcos2t_ste = steC ;
        statStruct.ap.nC = nC ;
        statStruct.ap.Qsin2t_avg = meanS ;
        statStruct.ap.Qsin2t_std = stdS ;
        statStruct.ap.Qsin2t_ste = steS ;
        statStruct.ap.nS = nS ;
        % cells in lobes 1 and 2
        statStruct.lobe12 = struct() ;
        statStruct.lobe12.Qcos2t_avg = meanC_L12(1) ;
        statStruct.lobe12.Qcos2t_std = stdC_L12(1) ;
        statStruct.lobe12.Qcos2t_ste = steC_L12(1) ;
        statStruct.lobe12.nC = nC_L12(1) ;
        statStruct.lobe12.Qsin2t_avg = meanS_L12(1) ;
        statStruct.lobe12.Qsin2t_std = stdS_L12(1) ;
        statStruct.lobe12.Qcos2t_ste = steC_L12(1) ;
        statStruct.lobe12.nS = nS_L12(1) ;
        
        stats{tidx} = statStruct ;
        
        Qct_a(tidx, :) = meanC_L12 ;
        Qct_s(tidx, :) = stdC_L12 ;
        Qct_e(tidx, :) = steC_L12 ;
        Qst_a(tidx, :) = meanS_L12 ;
        Qst_s(tidx, :) = stdS_L12 ;
        Qst_e(tidx, :) = steS_L12 ;
    end
    
end
% Save figure of Qx(x,t) and Qy(x,t)
subplot(1, 2, 1)
xlabel('ap position [pixels]', 'interpreter', 'latex')
ylabel('$\langle Q_{xx} \rangle$', 'interpreter', 'latex')
subplot(1, 2, 2)
xlabel('ap position [pixels]', 'interpreter', 'latex')
ylabel('$\langle Q_{yy} \rangle$', 'interpreter', 'latex')
imfn = fullfile(segDir, 'stats') ;
disp(['Saving statistics (x,t): ' imfn])
saveas(gcf, [imfn '.pdf']) ;
saveas(gcf, [imfn '.png']) ;

% Plot over time
keep = ~isnan(Qct_a(:, 1)) ;
Qct_ak = Qct_a(keep, 1) ;
Qct_sk = Qct_s(keep, 1) ;
Qct_ek = Qct_e(keep, 1) ;
Qst_ak = Qst_a(keep, 1) ;
Qst_sk = Qst_s(keep, 1) ;
Qst_ek = Qst_e(keep, 1) ;
timestamps = (timepoints(keep) * timeInterval) - t0 - time_offset(keep) ;
clf
hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
lineProps = {'.-', 'color', colors(2, :)} ;
h2=shadedErrorBar(timestamps, Qst_ak, Qst_sk, 'lineProps', lineProps) ;

plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(1, :)) ;
plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(1, :)) ;
plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(2, :)) ;
plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(2, :)) ;

legend([h1.patch, h2.patch], {'$Q_{xx}$', '$Q_{yy}$'}, 'interpreter', 'latex')

xlabel('time [min]', 'interpreter', 'latex')
ylabel('$Q_{xx}$, $Q_{yy}$', 'interpreter', 'latex')
imfn = fullfile(segDir, 'stats_dynamics') ;
disp(['Saving statistics (x,t): ' imfn])
saveas(gcf, [imfn '.pdf']) ;
saveas(gcf, [imfn '.png']) ;

% Save the data
disp('saving the data')
save('./cellSegmentation/results.mat', ...
    'keep', ...
    'Qct_ak', 'Qct_sk', 'Qst_ak', 'Qst_sk', ...
    'Qct_ek', 'Qst_ek', ...
    'stats', ...
    'timestamps', 'exciteIdx')

%% Compare to WT
imfn = fullfile(segDir, 'stats_compareWT') ;
% WTfn =  '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary.mat' ;
WTfn =  '/mnt/data/48Ygal4UASCAAXmCherry/201902072000_excellent/Time6views_60sec_1p4um_25x_obis1p5_2/data/deconvolved_16bit/msls_output/cellSegmentation/seg3d_corrected/stats_summary_L12.mat' ;
tmp = load(WTfn) ;

% Qxx
close all
subplot(1, 2, 1)
hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1=shadedErrorBar(timestamps, Qct_ak, Qct_sk, 'lineProps', lineProps) ;
plot(timestamps, Qct_ak-Qct_ek, '--', 'color', colors(1, :)) ;
plot(timestamps, Qct_ak+Qct_ek, '--', 'color', colors(1, :)) ;

lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tmp.timeStamps, tmp.meanxsL12, tmp.stdxs, 'lineProps', lineProps) ;
plot(tmp.timeStamps, tmp.meanxsL12 + tmp.stdmeanxs, '--', 'color', colors(2, :)) ;
plot(tmp.timeStamps, tmp.meanxsL12 - tmp.stdmeanxs, '--', 'color', colors(2, :)) ;
xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
ylabel('$Q_{xx}$', 'interpreter', 'latex')

% Qyy
subplot(1, 2, 2)
cla ;
hold on;
lineProps = {'.-', 'color', colors(1, :)} ;
h1=shadedErrorBar(timestamps, Qst_ak, Qst_sk, 'lineProps', lineProps) ;
plot(timestamps, Qst_ak-Qst_ek, '--', 'color', colors(1, :)) ;
plot(timestamps, Qst_ak+Qst_ek, '--', 'color', colors(1, :)) ;

lineProps = {'.-', 'color', colors(2, :)} ;
h2 = shadedErrorBar(tmp.timeStamps, tmp.meanysL12, tmp.stdys, 'lineProps', lineProps) ;
plot(tmp.timeStamps, tmp.meanysL12 + tmp.stdmeanys, '--', 'color', colors(2, :)) ;
plot(tmp.timeStamps, tmp.meanysL12 - tmp.stdmeanys, '--', 'color', colors(2, :)) ;
xlabel(['time [' timeUnits ']'], 'interpreter', 'latex')
ylabel('$Q_{yy}$', 'interpreter', 'latex')
saveas(gcf, [imfn '.png'])
saveas(gcf, [imfn '.pdf'])


