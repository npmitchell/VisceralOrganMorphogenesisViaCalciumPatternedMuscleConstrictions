%% GENERATE_AXISYMMETRIC_PULLBACK =========================================
% Pipeline for creating 'axisymmetric' pullbacks of the growing Drosophila 
% midgut about a centerline (previously computed). Uses orbifold method to
% guarantee conformality except at x=0 and x=L
% Run from anywhere, just specify dataDir below.
%
% Uses dorsal points for constructing curve loaded from: 
%    dpFile = fullfile( cylCutDir, 'ap_boundary_dorsalpts.h5' )
% which are saved in slice_mesh_endcaps.m
%
% Requirements
% ------------
% imsane
%
% Inputs
% ------
% meshes as topological cylinders, in meshDir/cylindercut/ 
%       --> run slice_mesh_endcaps.m
% meshDir/translation_APDV.txt
% meshDir/rotation_APDV.txt% meshDir/xyzlim_APDV_um.txt
% 
% Saves
% -------
% meshDir/meshStack.mat
% PullbackImages_###step/ (images)
% spcutMeshSmRS: struct with fields
%       f:  2*(nU-1)*(nV-1) x 3 int array
%               connectivity list (faces indexing vertices)
%       v:  nU*nV x 3 float array
%               3d coordinates of the mesh vertices
%       u:  nU*nV x 2 float array
%               2d coordinates in units of (avg centerline pathlength of this DVhoop, 1/2pi radians of hoop angle)
%       vn: nU*nV x 3 float array
%               smoothed vertex normals
%
% By Dillon Cislo and Noah Mitchell
%==========================================================================

clear; close all; clc;

%% ------------------------------------------------------------------------
% Run setup.m from imsane
% -------------------------------------------------------------------------

%% Parameters
% Overwriting options
overwrite_cleancenterlines = false ;
overwrite_pullbacks = false ;
overwrite_cleanCylMesh = false ;
overwrite_cutMesh = false ;
overwrite_spcutMesh = false ;
overwrite_writhe = false ;
overwrite_SmRSIms = false ;
overwrite_spcutMeshSm = false ;
overwrite_folds = false ;
overwrite_lobedynamics = false ;
overwrite_foldims = false ;
overwrite_lobeims = false ;
overwrite_spcutMesh_smoothradii = false ;
overwrite_piv = false ;
% Other options for what to do
phi_method = 'texture' ; % options are 'texture' and 'curves3d'
generate_sphi_coord = true ;
generate_uphi_coord = false ;
save_ims = true ;
nV = 100 ;
nU = 100 ;
dvexten = sprintf('_nU%04d_nV%04d', nU, nV) ;
nCurves_yjitter = 100 ;
nCurves_sphicoord = 1000 ;
% Plotting params
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;
% Parameters for cutMesh creation
nsegs4path = 5 ;
a_fixed = 2 ;
% Parameters for spcutMesh creation
normal_shift = 10 ;
maxJitter = 100 ;
maxTwChange = 0.15 ;
% for phi0 calculation via texture matching
lowerboundy = -350 ;
upperboundy = 350 ;
step_phi0tile = 25 ;
width_phi0tile = 150 ;
potential_sigmay = 350 ;


%% Add paths
% Add some necessary code to the path (ImSAnE should also be setup!) ------
if strcmp(computer, 'MACI64')
    rootdir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/'; 
    gutpath = fullfile(rootdir, 'gut_matlab/') ;
    addpath(genpath(fullfile(gutpath, 'euclidean_orbifolds')));
    addpath(genpath(fullfile(gutpath, 'gptoolbox')));
    addpath_recurse(fullfile(rootdir, 'imsane/external/')) ;
    
    % Setup a working directory for the project, where extracted surfaces,
    % metadata and debugging output will be stored.  Also specifiy the
    % directory containing the data.
    dataDir = ['/Volumes/Pal/48YGal4UASCAAXmCh/201902072000_excellent/', ...
        'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit/'] ;
    cd(dataDir)
else
    codepath = '/mnt/data/code/' ;
    gutpath = '/mnt/data/code/gut_matlab/' ;
    addpath(genpath('/mnt/crunch/djcislo/MATLAB/euclidean_orbifolds'));
    addpath(genpath('/mnt/data/code/gptoolbox'));
    addpath_recurse('/mnt/data/code/imsaneV1.2.3/external/') ;
    
    % Setup a working directory for the project, where extracted surfaces,
    % metadata and debugging output will be stored.  Also specifiy the
    % directory containing the data.
    dataDir = [ '/mnt/crunch/48Ygal4UASCAAXmCherry/201902072000_excellent/', ...
        'Time6views_60sec_1.4um_25x_obis1.5_2/data/deconvolved_16bit_smallbox/' ];
    cd(dataDir)
end

addpath(genpath(fullfile(gutpath, 'TexturePatch')));
addpath(gutpath) ;
addpath_recurse(fullfile(gutpath, ['basics' filesep])) ;
addpath_recurse(fullfile(gutpath, ['axisymmetric_pullbacks' filesep])) ;
addpath_recurse(fullfile(gutpath, ['plotting' filesep])) ;
addpath_recurse(fullfile(gutpath, ['mesh_handling' filesep])) ;
addpath_recurse(fullfile(gutpath, ['h5_handling' filesep])) ;
addpath_recurse(fullfile(gutpath, ['geometry' filesep])) ;
addpath_recurse(fullfile(gutpath, ['curve_functions' filesep])) ;
addpath(fullfile(gutpath, ['ExtPhaseCorrelation' filesep])) ;
addpath(fullfile(gutpath, 'savgol')) ;
addpath(fullfile(codepath, 'DEC')) ;
% addpath(genpath('/mnt/crunch/djcislo/MATLAB/TexturePatch'));

%% Define some colors
[colors, color_names] = define_colors() ;
blue = [0 0.4470 0.7410] ;
orange = [0.8500 0.3250 0.0980] ;
yellow = [0.9290, 0.6940, 0.1250] ;
purple = [0.4940, 0.1840, 0.5560] ;
green = [0.4660, 0.6740, 0.1880] ;
sky = [0.3010, 0.7450, 0.9330] ;
red = [0.6350, 0.0780, 0.1840] ;
brick = [0.800000 0.250000 0.330000] ;
light_green =[0.560000 0.930000 0.560000] ;
light_gray = [0.830000 0.830000 0.830000] ;
bwr = diverging_cmap([0:0.01:1], 1, 2) ;

%% Initialize ImSAnE Project ==============================================
projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd( projectDir );

% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = project.Experiment(projectDir, dataDir);

% Set file and experiment meta data
%
% Set required additional information on the files.
%
% We assume on individual image stack for each time point, labeled by time.
%  To be able to load the stack, we need to tell the project wehre the data
%  is, what convention is assumed for the file names, available time
%  points, and the stack resolution.  Options for modules in ImSAnE are
%  organized in MATLAB structures, i.e a pair of field names and values are
%  provided for each option.
%
% The following file metadata information is required:
%
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of itmes available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
%
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th

% A filename base template - to be used throughout this script
fileNameBase = 'Time_%06d_c1_stab';

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fileNameBase, '.tif'];
fileMeta.nChannels          = 1;
fileMeta.timePoints         = 110:263;
fileMeta.stackResolution    = [.2619 .2619 .2619];
fileMeta.swapZT             = 1;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is required:
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

expMeta                     = struct();
expMeta.channelsUsed        = 1;
expMeta.channelColor        = 1;
expMeta.description         = 'Apical membrane in Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% Now set the meta data in the experiment.
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

clear fileMeta expMeta

%% Initialize Some Directory Definitions ==================================
% The top level data directory
meshDir = [ dataDir ...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1' ];

% The dataset name base for the AD points
ADBase = '/mesh_apical_stab_%06d/adorsal';

% The dataset name base for the PD points
PDBase = '/mesh_apical_stab_%06d/pdorsal';

% The folder where the pullback images will be saved
nshift = strrep(sprintf('%03d', normal_shift), '-', 'n') ;
shiftstr = ['_' nshift 'step'] ;
imFolder = fullfile(meshDir, ['PullbackImages' shiftstr] ) ;
imFolder_e = [imFolder '_extended'] ;
imFolder_r = [imFolder '_relaxed'] ;
imFolder_re = [imFolder '_relaxed_extended'] ;
pivDir = fullfile(meshDir, 'piv') ;
meshBase = fullfile(fullfile(meshDir, 'meshes'), 'mesh_apical_stab_%06d.ply') ;
alignedMeshBase = fullfile(fullfile(meshDir, 'aligned_meshes'), 'mesh_apical_stab_%06d_APDV_um.ply') ;
cutFolder = fullfile(meshDir, 'cutMesh') ;
cutMeshBase = fullfile(cutFolder, [fileNameBase, '_cutMesh.mat']) ;
cutMeshImagesDir = fullfile(cutFolder, 'images') ;
cylCutDir = fullfile(meshDir, 'cylindercut') ;
cylCutMeshOutDir = fullfile(cylCutDir, 'cleaned') ;
cylCutMeshOutImDir = fullfile(cylCutMeshOutDir, 'images') ;
centerlineDir = fullfile(meshDir, 'centerline') ;
centerlineBase = fullfile(centerlineDir, 'mesh_apical_stab_%06d_centerline_exp1p0_res*.txt') ;
cntrsFileName = fullfile(centerlineDir, 'mesh_apical_stab_%06d_centerline_scaled_exp1p0_res*.txt') ;

% centerline from DV hoops
clineDVhoopDir = fullfile(centerlineDir, ['centerline_from_DVhoops' shiftstr dvexten]) ;
clineDVhoopBase = fullfile(clineDVhoopDir, 'centerline_from_DVhoops_%06d.mat');
clineDVhoopImDir = fullfile(clineDVhoopDir, 'images') ;
clineDVhoopFigBase = fullfile(clineDVhoopImDir, 'clineDVhoop_%06d.png') ;
radiusDir = fullfile(meshDir, ['radiusDVhoops' shiftstr dvexten]) ;
radiusImDir = fullfile(radiusDir, 'images') ;

% Writhe dir for DV hoops
writheDir = fullfile(clineDVhoopDir, 'writhe') ;

% The file name base for the cylinder meshes
cylinderMeshBase = fullfile( cylCutDir, ...
    'mesh_apical_stab_%06d_cylindercut.ply' );
cylinderMeshCleanBase = fullfile( cylCutMeshOutDir, ...
    'mesh_apical_stab_%06d_cylindercut_clean.ply' );
cylinderMeshCleanFigBase = fullfile( cylCutMeshOutImDir, ...
    'mesh_apical_stab_%06d_cylindercut_clean.png' );

% Lobe identification paths
lobeDir = fullfile(meshDir, 'lobes') ;
foldHoopImDir = fullfile(lobeDir, 'constriction_hoops/') ;

% Define cutMesh directories
sphiDir = fullfile(meshDir, ['sphi_cutMesh' shiftstr dvexten]) ;
spcutMeshBase = fullfile(sphiDir, 'mesh_apical_stab_%06d_spcutMesh.mat') ;
sphiSmDir = fullfile(sphiDir, 'smoothed') ;
sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
spcutMeshSmBase = fullfile(sphiSmDir, '%06d_spcutMeshSm.mat') ;
spcutMeshSmRSBase = fullfile(sphiSmRSDir, '%06d_spcutMeshSmRS.mat') ;
spcutMeshSmRSCBase = fullfile(sphiSmRSCDir, '%06d_spcMSmRSC.mat') ;
imFolder_sp = [imFolder '_sphi' dvexten] ;
imFolder_sp_e = [imFolder '_sphi' dvexten '_extended'] ;
imFolder_up = [imFolder '_uphi' dvexten] ;
imFolder_up_e = [imFolder '_uphi' dvexten '_extended'] ;
% time-averaged meshes
imFolder_spsm = fullfile(imFolder_sp, 'smoothed') ;
imFolder_spsm_e = fullfile(imFolder_sp_e, 'smoothed') ;
imFolder_spsm_e2 = fullfile(imFolder_sp_e, 'LUT_smoothed') ;  % raw LUT, no histeq
imFolder_rsm = fullfile(imFolder_r, 'smoothed') ;
imFolder_rsm_e = fullfile(imFolder_re, 'smoothed') ;
% For streadying the images via phi0
phi0fitBase = fullfile(sphiDir, 'phi0s_%06d_%02d.png') ; 

% The file containg the AD/PD points
dpFile = fullfile( cylCutDir, 'ap_boundary_dorsalpts.h5' );

% Folder for curvature measurements
KHSmDir = fullfile(sphiSmRSCDir, 'curvature') ;
KSmDir = fullfile(KHSmDir, 'gauss') ;
HSmDir = fullfile(KHSmDir, 'mean') ;

tomake = {imFolder, imFolder_e, imFolder_r, imFolder_re,...
    pivDir, cutFolder, cutMeshImagesDir, cylCutMeshOutDir,...
    cylCutMeshOutImDir, clineDVhoopDir, clineDVhoopImDir, writheDir, ...
    sphiDir, imFolder_sp, imFolder_sp_e, imFolder_sp, imFolder_sp_e,  ...
    lobeDir, radiusDir, radiusImDir, ...
    sphiSmDir, sphiSmRSImDir, sphiSmRSDir, sphiSmRSCDir, ...
    sphiSmRSPhiImDir, ...
    imFolder_spsm, imFolder_rsm, imFolder_spsm_e, ...
    foldHoopImDir, KHSmDir, KSmDir, HSmDir, ...
    fullfile(KSmDir, 'dorsal'), fullfile(KSmDir, 'ventral'), ...
    fullfile(KSmDir, 'latleft'), fullfile(KSmDir, 'latleft'), ...
    fullfile(HSmDir, 'dorsal'), fullfile(HSmDir, 'ventral'), ...
    fullfile(HSmDir, 'latleft'), fullfile(HSmDir, 'latleft')} ;
for i = 1:length(tomake)
    dir2make = tomake{i} ;
    if ~exist( dir2make, 'dir' )
        mkdir(dir2make);
    end
end
clearvars tomake dir2make i nshift 

%% Load rotation, translation, resolution
rot = dlmread(fullfile(meshDir, 'rotation_APDV.txt')) ;
buff = 20 ;
xyzlim_raw = dlmread(fullfile(meshDir, 'xyzlim.txt'), ',', 1, 0) ;
xyzlim = dlmread(fullfile(meshDir, 'xyzlim_APDV_um.txt'), ',', 1, 0) ;
xyzlim(:, 1) = xyzlim(:, 1) - buff ;
xyzlim(:, 2) = xyzlim(:, 2) + buff ;
xyzlim_APDV = xyzlim ;
trans = dlmread(fullfile(meshDir, 'translation_APDV.txt'));
resolution = dlmread(fullfile(meshDir, 'resolution.txt'), ',', 1, 0) ;
if resolution(:) == resolution(1)
    resolution = resolution(1) ;
else
    error('Have not handled case for anisotropic resolution')
end
clearvars buff

%% Identify anomalies in centerline data
ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold
disp('Searching for anomalous centerlines')
timePoints = xp.fileMeta.timePoints ;

cleanclinefn = fullfile(centerlineDir, 'centerlines_anomalies_fixed.mat') ;
if exist(cleanclinefn, 'file') || ~overwrite_cleancenterlines
    load(cleanclinefn, 'cntrlines')
else
    ssrs = zeros(length(timePoints) - 1, 1) ;
    cntrlines = cell(length(timePoints), 1) ;
    cntrlines_rs = cell(length(timePoints), 1) ;
    anomalous = false(length(timePoints), 1) ;
    for kk = 1:length(timePoints)
        t = timePoints(kk) ;

        % load the centerline for this timepoint
        clfns = dir(sprintf(centerlineBase, t)) ;
        cline = dlmread(fullfile(clfns(1).folder, clfns(1).name)) ;
        cntrlines{t} = cline ; 
        tmp = ((rot * cline')' + trans) * resolution ; 
        cntrlines_rs{t} = cat(2, ss_from_xyz(tmp), tmp) ;

        % compare to previous
        if t > timePoints(1)
            idx = pointMatch(cline, prevcline) ;
            ssr = sum(vecnorm(cline - prevcline(idx, :), 2, 2)) / length(idx) ;
            anomal = ssr > ssr_thres ;
            % Declare that a bad centerline was found
            if anomal
                disp(['Found t= ' num2str(t) ' to be anomalous'])
            end
        else
            ssr = 0 ;
            anomal = false ;
        end

        % rename as previous
        if ~anomal
            prevcline = cline ;
        end

        % Store for a recap
        ssrs(kk) = ssr * resolution ;
        anomalous(kk) = anomal ;
    end
    clearvars kk cline

    % Now adjust the anomalous ones to fix them up
    ssrfix = zeros(length(anomalous), 1) ;
    for qq = find(anomalous)
        % Find the last timepoint that was NOT anomalous and precedes the
        % current one
        previd = find(~anomalous & (timePoints' < timePoints(qq)), 1, 'last') ;
        prevtp = timePoints(previd) ;
        nextid = find(~anomalous & (timePoints' > timePoints(qq)), 1, 'first') ;
        cl1 = cntrlines{timePoints(previd)} ;
        cl2 = cntrlines{timePoints(nextid)} ;
        id = pointMatch(cl1, cl2) ;
        cline = (cl2(id, :) + cl1) * 0.5 ;
        cntrlines{timePoints(qq)} = cline ;

        % Do the same for rs centerline
        cl1 = cntrlines_rs{timePoints(previd)} ;
        cl2 = cntrlines_rs{timePoints(nextid)} ;
        id = pointMatch(cl1(:, 2:4), cl2(:, 2:4)) ;
        cline = (cl2(id, 2:4) + cl1(:, 2:4)) * 0.5 ;
        cntrlines_rs{timePoints(qq)} = cat(2, ss_from_xyz(cline), cline) ;

        % Prepare to plot the fixed SSR    
        % load the centerline for prev timepoint
        clfns = dir(sprintf(centerlineBase, prevtp)) ;
        prevcline = dlmread(fullfile(clfns(1).folder, clfns(1).name)) ;
        % Find the SSR
        idx = pointMatch(cline, prevcline) ;
        ssrfix(qq) = sum(vecnorm(cline - prevcline(idx, :), 2, 2)) / length(idx) ;
    end
    clearvars cl1 cl2 qq

    % Plot the anomalous determination
    close all; figure('visible', 'off')
    plot(timePoints, ssrs, '-')
    hold on;
    plot(timePoints(anomalous), ssrs(anomalous), 'o')
    plot(timePoints(anomalous), ssrfix(anomalous), 's')
    legend({'S.S.R.', 'anomalous', 'adjusted'})
    xlabel('time [min]')
    ylabel('Sum of squared residuals [\mum]')
    title('Identification of anomalous centerlines')
    saveas(gcf, fullfile(centerlineDir, 'centerlines_anomalies_fixed.png'))
    close all

    % Save the corrected centerlines as a cell array
    save(cleanclinefn, 'cntrlines') 

    % View the corrected centerlines
    for qq = timePoints
        first = max(min(timePoints), qq - 1) ;
        last = min(max(timePoints), qq + 1) ;
        for kk = first:last
            cl = cntrlines{kk} ;
            plot3(cl(:, 1), cl(:, 2), cl(:, 3))
            hold on;
        end
        hold off
        xlim(xyzlim_raw(1, :))
        ylim(xyzlim_raw(2, :))
        zlim(xyzlim_raw(3, :))
        title('Visualizing the corrected centerlines')
        pause(0.001)
    end
    clearvars cl

    % View the corrected centerlines
    for qq = timePoints
        first = max(min(timePoints), qq - 1) ;
        last = min(max(timePoints), qq + 1) ;
        for kk = first:last
            cl = cntrlines_rs{kk} ;
            plot3(cl(:, 2), cl(:, 3), cl(:, 4))
            hold on;
        end
        hold off
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        title('Visualizing the corrected RS centerlines')
        pause(0.001)
    end
    clearvars cl
    close all
    disp('done')
    clearvars ssrs ssr anomal anomalous
end


%% Iterate Through Time Points to Create Pullbacks ========================
% Check if cutmeshes already saved
mstckfn = fullfile(meshDir, 'meshStack_orbifold.mat') ;
spmstckfn = fullfile(meshDir, 'spmeshStack_orbifold.mat') ;
outcutfn = fullfile(cutFolder, 'cutPaths_%06d.txt') ;
outadIDxfn = fullfile(cylCutMeshOutDir, 'adIDx.h5') ;
outpdIDxfn = fullfile(cylCutMeshOutDir, 'pdIDx.h5') ;

% if exist(mstckfn, 'file') && ~overwrite_meshStack
%     % The cutPaths.h5 is loading here
%     disp(['Not loading meshStack since we can load those one by one'])
%     % disp(['Loading meshStack: ' mstckfn])
%     % load(mstckfn)
%     disp(['Not loading spmeshStack since we load those one by one'])
%     
%     if resave_imswhile
%         for t = xp.fileMeta.timePoints
%             aux_resave_cutpath_figs
%         end
%     end
% else
%     disp('meshStack is not on disk or is to be overwritten, compute...')
%     meshStack = cell( length(xp.fileMeta.timePoints), 1 );
%     spmeshStack = cell( length(xp.fileMeta.timePoints), 1 );

for t = xp.fileMeta.timePoints(154:end)
    disp(['NOW PROCESSING TIME POINT ', num2str(t)]);
    tidx = xp.tIdx(t);
    
    % Load the data for the current time point ------------------------
    xp.setTime(t) ;
    % Load or compute clean cylindrical mesh
    mesh3dfn =  sprintf( cylinderMeshCleanBase, t ) ;
    if ~exist(mesh3dfn, 'file') || overwrite_cleanCylMesh
        disp('Overwriting/Computing clean cylinderMesh')
        % Load the cylinder mesh
        cylmeshfn = sprintf( cylinderMeshBase, t ) ;
        mesh = read_ply_mod( cylmeshfn );
        % Clean the mesh, including averaging normals onto vertices using
        % angle-weighted scheme
        mesh = cleanCylMesh(mesh) ;   
        [adIDx, pdIDx] = aux_adjust_dIDx(mesh, t, dpFile, ADBase, PDBase, cylinderMeshCleanBase, outadIDxfn, outpdIDxfn, xp) ;

        %% Save the 3d cut mesh with new indices
        % This is saving the cylinder meshes with no ears. Also adIDx
        % is saved in h5file.
        plywrite_with_normals(mesh3dfn, mesh.f, mesh.v, mesh.vn)
        % Save adIDx with new indices
        save_to_h5(outadIDxfn, ['/' sprintf('%06d', t) ], adIDx, ['adIDx for t=' num2str(t) ' already exists'])
        % Save pdIDx with new indices
        save_to_h5(outpdIDxfn, ['/' sprintf('%06d', t) ], pdIDx, ['pdIDx for t=' num2str(t) ' already exists'])
        disp('done with cylindermesh cleaning')

         % View results --------------------------------------------------------
        mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
        if save_ims
            aux_plot_cleanCylMesh
        end
    else
        mesh = read_ply_mod(mesh3dfn) ;
        adIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t)]) ;
        pdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t)]) ;

        % View results --------------------------------------------------------
        mesh3dfigfn = sprintf( cylinderMeshCleanFigBase, t ) ;
        if (~exist(mesh3dfigfn, 'file') || overwrite_cleanCylMesh) && save_ims
            aux_plot_cleanCylMesh
        end
    end

    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(cutMeshBase, t) ;
    if ~exist(cutMeshfn, 'file') || overwrite_cutMesh
        if overwrite_cutMesh
            fprintf('Generating to overwrite cutMesh...')
        else
            fprintf('cutMesh not saved. Generating cutMesh... ');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % try geodesic if first timepoint
        if t == xp.fileMeta.timePoints(1)
            cutOptions.method = 'fastest' ;
            disp(['Cutting mesh using method ' cutOptions.method])
            cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
            cutP = cutMesh.pathPairs(:, 1) ;
            adIDx = cutP(1) ;
            pdIDx = cutP(end) ;
            prevTw = twist(mesh.v(cutP, :), cntrlines{t}) ;
            compute_pullback = true ;
        else 
            % If a previous Twist is not held in RAM, compute it
            % if ~exist('prevTw', 'var')
            % Load previous mesh and previous cutP
            prevcylmeshfn = sprintf( cylinderMeshCleanBase, t-1) ;
            prevmesh = read_ply_mod( prevcylmeshfn ); 
            prevcutP = dlmread(sprintf(outcutfn, t-1), ',', 1, 0) ;
            previousP = prevmesh.v(prevcutP, :) ;
            % Load previous centerline in raw units
            prevcntrfn = sprintf(centerlineBase, t-1) ;
            prevcline = cntrlines{t-1} ;
            % Compute Twist for this previous timepoint
            prevTw = twist(previousP, prevcline) ;
            
            % [edgelen, annulusv2d] = DiscreteRicciFlow.EuclideanRicciFlow(mesh.f, mesh.v, ...
            %     'BoundaryType', 'fixed', 'BoundaryShape', 'Circles', ...
            %     'MaxIter', 25, 'MaxCircIter', 21, ...
            %     'Tolerance', 1e-6, 'CircTolerance', 1e-4, 'PCGTolerance', 1e-4) ;
            % findAnnularPathZeroWindingNumber(mesh.f, annulusv2d, adIDx, pdIDx)
            
            % Which path to match this one to: choose previous timepoint
            % Load previous mesh and previous cutP
            prevcylmeshfn = sprintf( cylinderMeshCleanBase, t-1) ;
            prevmesh = read_ply_mod( prevcylmeshfn ); 
            prevcutP = dlmread(sprintf(outcutfn, t-1), ',', 1, 0) ;
            previousP = prevmesh.v(prevcutP, :) ;
            
            [cutMesh, adIDx, pdIDx, cutP, prevTw] = ...
                generateCutMeshFixedTwist(mesh, adIDx, pdIDx, ...
                cntrlines{t},...  % supply the current corrected centerline
                nsegs4path, prevTw, previousP, ...
                'MaxTwChange', maxTwChange, 'MaxJitter', maxJitter, ...
                'PrevCntrline', prevcline) ;
            compute_pullback = true ;  
        end

        % Store this path for the next one to be nearby
        % Limit the number of segments to nsegs4path
        % previousP = cutMesh.v(cutP, :) ;
        % pstep = round(length(cutP) / nsegs4path ) ;
        % previousP = previousP(1:pstep:end, :) ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Done with generating initial 3D CutMesh with cutPath\n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save the cutPath to txt file
        header = 'cutP (path of cut), indexing into vertices' ;
        write_txt_with_header(sprintf(outcutfn, t), cutP, header)  
        clearvars header

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate pullback to rectangular domain ---------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The surface parameterization algorithm (optionally) takes four vertex IDs
        % as input to specify the corners of the square parameterization domain.
        % Maddeningly, the order in which these points are specified to not seem to
        % effect the output. For consistency, we perform a post-hoc correction so
        % that the final output has the following geometric ordering
        %
        %   (AD1)-------(PD1)
        %     |           |
        %     |           |
        %     |           |
        %     |           |
        %   (AD2)-------(PD2)
        %
        % Note that the pathPairs variable has the following columns:
        %   [ ( AD1 -> PD1 ), ( AD2 -> PD2 ) ]
        %--------------------------------------------------------------------------

        % View results --------------------------------------------------------
        % P = cutMesh.pathPairs(:,1);
        % 
        % trisurf( triangulation( mesh.f, mesh.v ) );
        % 
        % hold on
        %
        % line( mesh.v(P,1), mesh.v(P,2), mesh.v(P,3), ...
        %     'Color', 'c', 'LineWidth',2);
        % 
        % scatter3( mesh.v(adIDx,1), mesh.v(adIDx,2), mesh.v(adIDx,3), ...
        %     'filled', 'r' );
        % scatter3( mesh.v(pdIDx,1), mesh.v(pdIDx,2), mesh.v(pdIDx,3), ...
        %     'filled', 'm' );
        % 
        % hold off
        % 
        % axis equal
        % 
        % clear P

        %----------------------------------------------------------------------
        % Generate Pullback to Annular Orbifold Domain
        %----------------------------------------------------------------------
        fprintf('Relaxing network via Affine transformation... ');
        cutMesh = flattenAnnulus( cutMesh );

        % Find lateral scaling that minimizes spring network energy
        ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
        % Assign scaling based on options: either a0 or a_fixed
        % if tidx == 1 && ~a_fixed
        %     a_fixed = ar ;
        % end      
        % a = a_fixed ;

        % Scale the x axis by a or ar
        uvtx = cutMesh.u ;
        cutMesh.u = [ a_fixed .* uvtx(:,1), uvtx(:,2) ];
        cutMesh.ar = ar ;
        cutMesh.umax = a_fixed ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Done flattening cutMesh. Now saving.\n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Save cutMesh
        save(cutMeshfn, 'cutMesh', 'adIDx', 'pdIDx', 'cutP')
        disp('done with plotting & saving cut')
    else
        fprintf('Loading Cut Mesh from disk... ');
        load(cutMeshfn, 'cutMesh') 
        cutP = dlmread(sprintf(outcutfn, t), ',', 1, 0) ;
        compute_pullback = ~isempty(cutP) ;
    end

    % Plot the cutPath (cutP) in 3D
    if save_ims && overwrite_cutMesh 
        disp('Saving cutP image')
        aux_plot_cutP
    end

    % Compute the pullback if the cutMesh is ok
    if compute_pullback
        disp(['Evolving mesh along normal shift for pullback images: shift=' num2str(normal_shift)])
        % Displace normally ---------------------------------------------------
        cutMesh.v = cutMesh.v + cutMesh.vn * normal_shift ;

        % todo !!!!
        % if 
        % else
        %     % Load the cutMesh
        %     load(sprintf(cutMeshBase, t), 'cutMesh')
        % end

        % View results --------------------------------------------------------
        % 
        % patch( 'Faces', cutMesh.f, 'Vertices', cutMesh.u, ...
        %     'FaceVertexCData', cutMesh.v(:,3), 'FaceColor', 'interp', ...
        %     'EdgeColor', 'k' );
        % 
        % hold on
        % 
        % cornerColors = [ 1 0 0; 1 0 1; 0 1 1; 0 1 0 ];
        % corners = [ adIDx cutMesh.pathPairs(1,2), ...
        %     cutMesh.pathPairs(end,1) pdIDx ];
        % 
        % scatter( cutMesh.u( corners, 1 ), cutMesh.u( corners, 2 ), [], ...
        %     cornerColors, 'filled' );
        % 
        % hold off
        % 
        % axis equal
        % 
        % clear cornerColors corners

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate s,phi coord system for rotated,scaled mesh (rs)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Establishing s,phi coord system\n');

        if ~exist(sprintf(spcutMeshBase, t), 'file') || overwrite_spcutMesh
            if overwrite_spcutMesh
                disp('Overwriting spcutMesh...')
            else
                disp('spcutMesh not on disk. Generating ...')
            end

            % Transform from u,v coordinates to s, phi coordinates
            % [scoords, phicoords] = generateSPhiFromUV();

            %----------------------------------------------------------------------
            % Generate tiled orbifold triangulation
            %----------------------------------------------------------------------
            tileCount = [1 1];  % how many above, how many below
            cutMeshrs = cutMesh;
            % Rotate and translate TV3D
            cutMeshrs.v = ((rot * cutMesh.v')' + trans) * resolution ;
            cutMeshrs.vn = (rot * cutMesh.vn')' ;
            [ ~, ~, TV3D, TVN3D ] = tileAnnularCutMesh( cutMesh, tileCount );
            [ TF, TV2D, TV3Drs ] = tileAnnularCutMesh( cutMeshrs, tileCount );

            %----------------------------------------------------------------------
            % Calculate abbreviated centerline from cutMesh boundaries
            %----------------------------------------------------------------------
            % Load centerline in raw units
            % cntrfn = sprintf(cntrsFileName, t) ;
            % cline = dlmread(cntrfn, ',') ;
            cline = cntrlines_rs{t} ; 
            ss = cline(:, 1) ;
            cline = cline(:, 2:end) ;
            disp('Finding relevant segment of centerline')
            [cseg, acID, pcID, bdLeft, bdRight] = centerlineSegmentFromCutMesh(cline, TF, TV2D, TV3Drs) ;

            %----------------------------------------------------------------------
            % Generate surface curves of constant s
            %----------------------------------------------------------------------
            % For lines of constant phi
            disp('Creating crude uv curves with du=const to define uspace by ds(u)')
            % Make grid
            eps = 1e-14 ;
            uspace0 = linspace( eps, cutMesh.umax - eps, nU )' ;
            vspace = linspace( eps, 1-eps, nV )' ;

            disp('Casting crude (equal dU) points into 3D...')
            crude_ringpath_ss = ringpathsGridSampling(uspace0, vspace, TF, TV2D, TV3Drs) ;

            % Resample crude_ringpath_ds made from uspace0 (equal du, not equal ds_3D in u direction)
            [uspace, eq_ringpath_ss] = equidistantSampling1D(linspace(0, 1, nU)', crude_ringpath_ss, nU, 'linear') ;
            % ensure that uspace is nU x 1, not 1 x nU
            uspace = reshape(uspace, [nU, 1]) ; 
            % hedge the first and last point to avoid NaNs
            eps = 1e-13 ;
            uspace(1) = uspace(1) + eps ;
            uspace(end) = uspace(end) - eps ;
            clearvars dsuphi curves3d uspace0

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Casting resampled points into 3D (approx equal ds_3D in u dir, but variable ds_3D in v dir)...')
            % NOTE: first dimension indexes u, second indexes v
            curves3d = zeros(nU, nV, 3) ;  % in units of um
            for kk = 1:nU
                if mod(kk, 50) == 0
                    disp(['u = ' num2str(kk / nU)])
                end
                uv = [cutMesh.umax * uspace(kk) * ones(size(vspace)), vspace] ;
                curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv) ;
            end 

            % Check the 3d curves 
            if preview
                figure ; hold on;
                for kk = 1:nU
                    plot3(curves3d(kk, :, 1), curves3d(kk, :, 2), curves3d(kk, :, 3), '.') 
                end
                axis equal
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Compute s(u) and radius(u) for "uniform"--> evenly sample each DV hoop (0,1) so ds_3D=const \n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Resample at evenly spaced dphi in embedding space (rs, in um)
            fprintf('Resampling curves...\n')
            c3d_dsv = zeros(size(curves3d)) ;  % in units of um
            for i=1:nU
                % Note: no need to add the first point to the curve
                % since the endpoints already match exactly in 3d and
                % curvspace gives a curve with points on either
                % endpoint (corresponding to the same 3d location).
                c3d_dsv(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d(i, :, :)), nV, true) ;
                if vecnorm(squeeze(c3d_dsv(i, 1, :)) - squeeze(c3d_dsv(i, end, :))) > 1e-7
                    error('endpoints do not join! Exiting')
                end

                % Visualization for Troubleshooting:
                % triplot(TF, TV2D(:, 1), TV2D(:, 2))
                % hold on;
                % plot(uv(:, 1), uv(:, 2), '.')
            end

            % Check the 3d curves 
            if preview
                figure ; hold on;
                for kk = 1:nU
                    plot3(c3d_dsv(kk, :, 1), c3d_dsv(kk, :, 2), c3d_dsv(kk, :, 3), '.') 
                end
                axis equal
            end

            fprintf('Finding s(u) and r(u) of resampled "uniform" c3ds [uniform ds in V dir]...\n')
            % mcline is the resampled centerline, with mss
            % avgpts is the raw Nx3 averaged hoops, with avgpts_ss
            [mss, mcline, radii_from_mean_uniform_rs, avgpts_ss, avgpts] = srFromDVCurves(c3d_dsv) ;

            % Used to find radius using original centerline
            % [ssv, radii, avgpts, cids] = srFromDVCurvesGivenCenterline(ss, cline, c3ds) ;
            % Could operate just on the centerline segment
            cseg_ss = ss(acID:pcID) ;
            % [ssv, radii, avgpts, cids] = srFromDVCurves(cseg_ss, cseg, c3ds) ;
            % 
            % Adjust the centerline indices to index into the full
            % centerline. Note that cseg_ss already does this for ss.
            % cids = cids + acID ;

            % Plot new centerline
            aux_plot_clineDVhoop(avgpts, avgpts_ss, cseg, cline, cseg_ss, curves3d, xyzlim, clineDVhoopFigBase, t)

            % Optional: clean curve with polynomial and point match
            % avgpts onto cleaned curve. Skipping for later.

            % Compute ringpath_ss, the mean distance traveled from one
            % line of constant u to the next
            disp('Computing ringpath_ss in "uniform" resampling (equal ds along DV)...')
            % The distance from one hoop to another is the
            % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
            dsuphi = reshape(vecnorm(diff(c3d_dsv), 2, 3), [nU-1, nV]) ;
            ringpath_ds = nanmean(dsuphi, 2) ;
            ringpath_ss = cumsum([0; ringpath_ds]) ;
            clearvars dsuphi ringpath_ds

            % Save new centerline in rotated translated units
            fn = sprintf(clineDVhoopBase, t) ;
            disp(['Saving new centerline to ' fn])
            save(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')

            % Note: radii_from_mean_uniform_rs is the radius of 
            % interpolated hoops, not the actual points

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Done making new centerline using uniformly sampled hoops\n') ;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Create new3d, the regridded pts at UV, moved to sphi')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            onesUV = ones(nU, nV) ;
            uu = uspace * cutMesh.umax .* onesUV ;
            vv = (vspace .* onesUV')' ;
            uv = [uu(:), vv(:)] ;
            % Note: here interpolate uv in the TV2D coord system, then
            % use uphi as the actual 2D coordinates for these vertices
            % NOTE: unlike curves3d, new3d is NOT rotated/translated/scaled
            new3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uv) ;
            
            IVloaded = false ;
            if t == xp.fileMeta.timePoints(1)
                % Store for next timepoint
                phiv = (vspace .* ones(nU, nV))' ;
                phi0s = zeros(size(uspace)) ;
                phi0_fit = phi0s ;
            else
                % Load previous sphi vertices in 3d 
                plotfn = sprintf(phi0fitBase, t, 0);
                if strcmp(phi_method, '3dcurves')
                    % Load the previous spcutMesh and call it prev3d_sphi
                    % Also note the previous spcutMesh pullback image's fn
                    tmp = load(sprintf(spcutMeshBase, ...
                        xp.fileMeta.timePoints(tidx-1)), 'spcutMesh') ;
                    prevf = tmp.spcutMesh.f ;
                    prev3d_sphi = reshape(tmp.spcutMesh.v, [nU, nV, 3]) ; 
                    imfn_sp_prev = sprintf( ...
                        fullfile([imFolder_sp, '/', fileNameBase, '.tif']), ...
                        xp.fileMeta.timePoints(tidx-1) ) ;

                    % fit the shifts in the y direction
                    dmyk = 0 ;
                    phi0_fit = zeros(size(uspace)) ;
                    phi0s = zeros(size(uspace)) ;
                    phi0_fit_kk = 1 ; % for first pass                
                    phiv_kk = (vspace .* ones(nU, nV))' ;
                    ensureDir([sphiDir, '/phi0_correction/'])
                    while any(phi0_fit_kk > 0.002) && dmyk < 6
                        disp(['Iteration ' num2str(dmyk)])
                        plotfn = sprintf(phi0fitBase, t, dmyk);
                        
                        % Will we save check pullbacks to preview the algo?
                        if save_phi0patch
                            patchImFn = sprintf( ...
                                fullfile(sphiDir, 'phi0_correction', [fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                                xp.fileMeta.timePoints(tidx-1) )  ;
                            geomImFn = sprintf( ...
                                fullfile(sphiDir, 'phi0_correction', ['3d' fileNameBase '_prephi0_' num2str(dmyk) '.tif']), ...
                                xp.fileMeta.timePoints(tidx-1) )  ;

                            % Load the intensity data for this timepoint
                            if ~IVloaded
                                % (3D data for coloring mesh pullback)
                                xp.loadTime(t);
                                xp.rescaleStackToUnitAspect();

                                % Raw stack data
                                IV = xp.stack.image.apply();
                                IV = imadjustn(IV{1});         
                                IVloaded = true ;
                            end
                            
                            % Texture patch options
                            Options.PSize = 5;
                            Options.EdgeColor = 'none';
                            % Texture image options
                            Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
                            Options.yLim = [0 1];
                            
                            % Roll options into a struct
                            patchOpts.patchImFn = patchImFn ;
                            patchOpts.imfn_sp_prev = imfn_sp_prev ;
                            patchOpts.IV = IV ;
                            patchOpts.ringpath_ss = ringpath_ss ;
                            patchOpts.Options = Options ;
                            patchOpts.v3d = new3d ;
                        else
                            patchOpts = [] ;
                        end
                        
                        % Minimize difference in DV hoop positions wrt
                        % previous pullback mesh                        
                        [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
                            uspace * cutMesh.umax, phiv_kk, prev3d_sphi, -0.45, 0.45, ...
                            save_ims, plotfn, save_phi0patch, preview, patchOpts) ;
                                              
                        % Update the result
                        dmyk = dmyk + 1;
                        phi0_fit = phi0_fit + phi0_fit_kk ;
                        phi0s = phi0s + phi0s_kk ;
                        phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
                        
                        
                        % plot mesh colored by the phase phi 
                        % previous timepoint
                        xtmp = prev3d_sphi(:, :, 1) ;
                        ytmp = prev3d_sphi(:, :, 2) ;
                        ztmp = prev3d_sphi(:, :, 3) ;
                        phitmp = (vspace .* ones(nU, nV))' ;
                        colormap parula ;
                        % cmap = parula ;
                        % colors = cmap(max(1, uint8(colortmp(:) * length(parula))), :) ;
                        trisurf(prevf, xtmp(:), ytmp(:), ztmp(:), phitmp(:), ...
                            'FaceColor', 'interp',...
                            'EdgeColor', 'none', 'FaceAlpha', 0.25)
                        axis equal
                        % freezeColors
                        
                        % before fitting
                        hold on;
                        pe0 = find(phitmp(:) < 1e-4 | phitmp(:) > 0.99) ;
                        plot3(new3d(pe0, 1), new3d(pe0, 2), ...
                            new3d(pe0, 3), '.')
                        % colormap copper
                        % trimesh(prevf, new3d(inds, 1), new3d(:, 2), new3d(:, 3),...
                        %     phitmp(:), 'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                        % freezeColors

                        % after fitting
                        hold on;                   
                        pekk = find(mod(phiv_kk(:), 1) < 1e-4 | mod(phiv_kk(:), 1) > 0.99) ;
                        plot3(new3d(pekk, 1), new3d(pekk, 2), ...
                            new3d(pekk, 3), '^')                        
                        % colormap summer
                        % trimesh(prevf, new3d(:, 1), new3d(:, 2), new3d(:, 3),...
                        %     mod(phiv_kk(:), 1), 'FaceColor', 'interp', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                        % freezeColors
                        xlabel('x [\mum]')
                        ylabel('y [\mum]')
                        zlabel('z [\mum]')
                        view(2)
                        saveas(gcf, geomImFn)
                    end
                    
                elseif strcmp(phi_method, 'texture')
                    imfn_sp_prev = sprintf( ...
                        fullfile([imFolder_sp, '/', fileNameBase, '.tif']), ...
                        xp.fileMeta.timePoints(tidx-1) ) ;
                    
                    % Load the intensity data            
                    % Load 3D data for coloring mesh pullback
                    xp.loadTime(t);
                    xp.rescaleStackToUnitAspect();

                    % Raw stack data
                    IV = xp.stack.image.apply();
                    IV = imadjustn(IV{1});         
                    IVloaded = true ;
            
                    % Texture patch options
                    Options.PSize = 5;
                    Options.EdgeColor = 'none';
                    % Texture image options
                    Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
                    Options.yLim = [0 1];

                    % fit the shifts in the y direction
                    % todo: save uncorrected patchIms,
                    % could try tiling twice...
                    dmyk = 0 ;
                    phi0_fit = zeros(size(uspace)) ;
                    phi0s = zeros(size(uspace)) ;
                    phi0_fit_kk = 1 ; % for first pass                
                    phiv_kk = (vspace .* ones(nU, nV))' ;
                    ensureDir([sphiDir, '/phi0_correction/'])
                    while any(phi0_fit_kk > 0.002) && dmyk < 6
                        disp(['Iteration ' num2str(dmyk)])
                        plotfn = sprintf(phi0fitBase, t, dmyk);
                        patchImFn = sprintf( ...
                            fullfile([sphiDir, '/phi0_correction/', fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                            xp.fileMeta.timePoints(tidx-1) )  ;
                        [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevPullback(IV, ...
                            new3d, cutMesh.umax, uspace * cutMesh.umax, phiv_kk, ...
                            ringpath_ss, imfn_sp_prev, lowerboundy, upperboundy, ...
                            save_ims, plotfn, Options, ...
                            step_phi0tile, width_phi0tile, potential_sigmay, 'integer', ...
                            patchImFn) ;
                        
                        % Update the result
                        dmyk = dmyk + 1;
                        phi0_fit = phi0_fit + phi0_fit_kk ;
                        phi0s = phi0s + phi0s_kk ;
                        phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
                    end
                else
                    error("Could not recognize phi_method: must be 'texture' or '3dcurves'")
                end
                close all

                % Store to save at this timepoint
                phiv = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
            end

            % NOTE: We have coordinates u,phiv that we associate with
            % the 3d coordinates already mapped to uv
            uphi = [uu(:), phiv(:)] ;
            
            % plot(uphi(:, 1), uphi(:, 2), '.')
            % xlabel('u')
            % ylabel('\phi')
            % waitfor(gcf)
            
            % Recompute radii_from_mean_uniform_rs as radii_from_avgpts 
            % NOTE: all radius calculations done in microns, not pixels
            sphi3d_rs = ((rot * new3d')' + trans) * resolution ;
            radii_from_avgpts = zeros(size(sphi3d_rs, 1), size(sphi3d_rs, 2)) ;
            for jj = 1:nU
                % Consider this hoop
                hoop = squeeze(sphi3d_rs(jj, :, :)) ;
                radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
            end

            % Triangulate the sphigrid and store as its own cutMesh
            % sphiv = zeros(nU, nV, 2) ;
            % sphiv(:, :, 1) = sv ;
            % sphiv(:, :, 2) = phiv ;
            sv = ringpath_ss .* onesUV ;
            % % Triangulate the mesh (topology is already known):
            % tmptri = defineFacesRectilinearGrid(sp, nU, nV) ;
            % % Old version did not assume topology as given:
            % tmptri = delaunay(sv(:), phiv(:)) ;
            % disp('orienting faces of delaunay triangulation (s,phi)')
            % tmptri = bfs_orient( tmptri );

            % Define path pairs for tiling the (s,phi) cut mesh
            spcutP1 = 1:nU;
            spcutP2 = nU*nV - fliplr(0:(nU-1)) ;
            spcutMesh.pathPairs = [ spcutP1', spcutP2' ];

            % Check to see if any members of pathPairs connect to
            % non-Nearest Neighbors. Not necessary now that we assume
            % known gridded mesh topology
            % cleantri = cleanBoundaryPath2D(tmptri, [sv(:), phiv(:)], spcutMesh.pathPairs(:), true) ;

            spcutMesh.f = defineFacesRectilinearGrid(uv, nU, nV) ;
            spcutMesh.nU = nU ;
            spcutMesh.nV = nV ;
            % First resampling
            spcutMesh.v0 = new3d ;
            % spcutMesh.vrs0 = ((rot * new3d')' + trans) * resolution ;
            % Define normals based on the original mesh normals
            spvn03d = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uphi) ;
            spvn03d = spvn03d ./ vecnorm(spvn03d, 2, 2) ;
            spcutMesh.vn0 = spvn03d ;
            spcutMesh.sphi0 = [sv(:), phiv(:)] ;
            spcutMesh.uphi0 = uphi ;
            % Note: uv has no direct relation with cutMesh, just a grid
            % for utility and reference, but it does have unequal 
            % spacing in u in anticipation of building sphi0 as a 
            % near perfect grid.
            spcutMesh.uv = uv ;  

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SECOND RESAMPLING
            % Make a new grid
            slin = linspace(0, max(spcutMesh.sphi0(:, 1)), nU) ;
            plin = linspace(0, 1, nV) ;
            [ss, pp] = meshgrid(slin, plin) ;
            % Push the endpoints on each boundary in by epsilon to
            % avoid NaNs
            eps = 1e-14 ;
            ss(:, 1) = eps ;
            ss(:, end) = ss(:, end) - eps ;
            % Transpose so that x increases with increasing index first
            ss = ss' ;
            pp = pp' ;
            sp = [ss(:), pp(:)] ;

            % Tile the spcutMesh
            tileCount = [2, 2] ;
            spcutMesh.u = spcutMesh.sphi0 ;
            spcutMesh.v = spcutMesh.v0 ;
            spcutMesh.vn = spcutMesh.vn0 ;
            [ faces, v2d, v3d, vn3d ] = tileAnnularCutMesh( spcutMesh, tileCount );
            spcutMesh = rmfield(spcutMesh, 'u') ;
            spcutMesh = rmfield(spcutMesh, 'v') ;
            spcutMesh = rmfield(spcutMesh, 'vn') ;
            spv3d = interpolate2Dpts_3Dmesh(faces, v2d, v3d, sp) ;
            % check the pts
            % plot3(spv3d(:, 1), spv3d(:, 2), spv3d(:, 3))  

            % also interpolate the normals
            spvn3d = interpolate2Dpts_3Dmesh(faces, v2d, vn3d, sp) ;
            spvn3d = spvn3d ./ vecnorm(spvn3d, 2, 2) ;

            % Define new faces for second rectilinear resampling
            % NOTE: not necessary since we already defined the topology
            % from the guess [sv(:), phiv(:)] stored as spcutMesh.sphi0
            % spcutMesh.f = defineFacesRectilinearGrid(sp, nU, nV) ;
            spcutMesh.sphi = sp ;
            spcutMesh.v = spv3d ;
            spcutMesh.vn = spvn3d ;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            spcutMesh.ringpath_ss = ringpath_ss ;
            spcutMesh.radii_from_mean_uniform_rs = radii_from_mean_uniform_rs ;  % from uniform DV sampling
            spcutMesh.radii_from_avgpts = radii_from_avgpts ;
            spcutMesh.mss = mss ;       % from uniform DV sampling, also stored in centerline
            spcutMesh.mcline = mcline ; % from uniform DV sampling, also stored in centerline
            spcutMesh.avgpts = avgpts ; % from uniform DV sampling, also stored in centerline
            spcutMesh.avgpts_ss = avgpts_ss ; % from uniform sampling, also stored in centerline

            % Define optimal isoareal Affine dilation factor in s
            % tmp = spcutMesh.sphi ;
            % tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
            % arsp = minimizeIsoarealAffineEnergy( spcutMesh.f, spcutMesh.v, tmp );
            % clearvars tmp
            spcutMesh.ar = cutMesh.ar ;

            % todo: check that u coords have not shifted upon
            % redefinition of sphi0 -> sphi

            % Save s,phi and their 3D embedding
            spcutMesh.phi0s = phi0s ;
            spcutMesh.phi0_fit = phi0_fit ;
            save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
        else
            disp('Loading spcutMesh from disk...')
            load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
            IVloaded = false ;

            % Load new centerline
            fn = sprintf(clineDVhoopBase, t) ;
            disp(['Loading new centerline from ' fn])
            load(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
        end
        fprintf('Done with generating S,Phi coords \n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Create pullback using S,Phi coords \n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %--------------------------------------------------------------
        % Generate Output Image Files
        %--------------------------------------------------------------
        imfn = sprintf( fullfile([imFolder, '/', fileNameBase, '.tif']), t ); 
        imfn_r = sprintf( fullfile([imFolder_r, '/', fileNameBase, '.tif']), t ) ;
        imfn_sp = sprintf( fullfile([imFolder_sp, '/', fileNameBase, '.tif']), t ) ;
        imfn_up = sprintf( fullfile([imFolder_up, '/', fileNameBase, '.tif']), t ) ;
        pullbacks_exist1 = exist(imfn, 'file') && exist(imfn_r, 'file') ;
        pullbacks_exist2 = exist(imfn_sp, 'file') && (exist(imfn_up, 'file') || ~generate_uphi_coord) ;
        if (~pullbacks_exist1 || ~pullbacks_exist2 || overwrite_pullbacks) && ~IVloaded
            % Load 3D data for coloring mesh pullback
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();

            % Raw stack data
            IV = xp.stack.image.apply();
            IV = imadjustn(IV{1});
        end

        if ~exist(imfn_sp, 'file') || overwrite_pullbacks
            fprintf(['Generating SP output image: ' imfn_sp]);
            % Assigning field spcutMesh.u to be [s, phi] (ringpath
            % and azimuthal angle)
            spcutMesh.u = spcutMesh.sphi ;
            aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp)
            spcutMesh = rmfield(spcutMesh, 'u') ;
        end

        if (~exist(imfn_up, 'file') || overwrite_pullbacks) && generate_uphi_coord
            fprintf(['Generating uphi output image: ' imfn_up]);
            % Assigning field spcutMesh.u to be [s, phi] (ringpath
            % and azimuthal angle)
            spcutMesh.u = spcutMesh.uphi ;
            aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up)
            spcutMesh = rmfield(spcutMesh, 'u') ;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate Output Image File -- regular UV coordinates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist(imfn, 'file') || overwrite_pullbacks
            % Generate output image in uv
            fprintf(['Generating output image: ' imfn]);
            aux_generate_orbifold(cutMesh, a_fixed, IV, imfn)
        else
            disp('Skipping pullback image generation since exists')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save relaxed image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if ~exist(imfn_r, 'file')
            disp('Generating relaxed image for sphi coords...')
            spcutMesh.u = spcutMesh.sphi ;
            aux_generate_orbifold(spcutMesh, spcutMesh.ar, IV, imfn_r)
            spcutMesh = rmfield(spcutMesh, 'u') ;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save submesh array. Each cell element contains all the 
        % submeshes for that TP, which in this case is just one.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % meshStack{tidx} = cutMesh ;
        % if generate_sphi_coord
        %     spmeshStack{tidx} = spcutMesh ;
        % end
        fprintf('Done\n');
    else
        disp('Skipping computation of pullback')
    end
    clear Options IV

%     %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
%     disp(['Saving meshStack to disk: ' mstckfn])
%     save(mstckfn, 'meshStack') ;
%     
%     %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
%     disp(['Saving spmeshStack to disk: ' spmstckfn])
%     if generate_sphi_coord
%         save(spmstckfn, 'spmeshStack') ;
%     end
end
clearvars t cutP spcutMesh spvn3d ss pp uv tileCount slin plin plotfn
clearvars IVloaded IV uphi sphi
disp('Done with generating spcutMeshes and cutMeshes')

%% Preview results ========================================================
check = false ;
if check
    aux_preview_results
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE IMAGES IN Y AND RESAVE ============================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDirs = {imFolder, imFolder_sp, imFolder_up} ;
imDirs_e = {imFolder_e, imFolder_sp_e, imFolder_up_e} ;
ntiles = 50 ;   %   The number of bins in each dimension for histogram equilization for a
                %   square original image. That is, the extended image will have (a_fixed *
                %   ntiles, 2 * ntiles) bins in (x,y).
for qq = 1:2
    options.histeq = true ;
    options.a_fixed = a_fixed ;
    options.ntiles = ntiles ;
    options.overwrite = overwrite_pullbacks;
    extendImages(imDirs{qq}, imDirs_e{qq}, fileNameBase, options)
    disp(['done ensuring extended tiffs for ' imDirs{qq} ' in ' imDirs_e{qq}])
end
disp('done')

%% FIND THE FOLDS SEPARATING COMPARTMENTS =================================
% First compute using the avgpts (DVhoop means)
guess123 = [0.2, 0.5, 0.8] ;
max_wander = 20 ; % max amount that DV hoop id can wander per timepoint
disp('Identifying lobes...')
foldfn = fullfile(lobeDir, ['fold_locations_sphi' dvexten '_avgpts.mat']) ;
if exist(foldfn, 'file') && ~overwrite_folds
    disp('Loading lobes')
    % Save the fold locations as a mat file
    load(foldfn, 'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset', ...
        'rssfold', 'rssfold_frac', 'rssmax', 'rmax')
    
else
    [folds, ssfold, ssfold_frac, ssmax, rmax, fold_onset] = ...
        identifyLobes(xp.fileMeta.timePoints,...
            spcutMeshBase, guess123, max_wander, preview, 'avgpts') ;
    
    % Compute ringpath pathlength for results found using centerline
    disp('Converting folds to ringpath_ss locations...')
    [rssfold, rssfold_frac, rssmax] = rssFromFoldID(folds, xp.fileMeta.timePoints, spcutMeshBase) ;

    % Save the fold locations as a mat file
    save(foldfn, 'rssfold', 'rssfold_frac', 'rssmax', 'rmax', ...
        'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset')
end
clearvars guess123 maxwander

% Plot results as both avgpts and ringpath distances
fold_ofn = dir(fullfile(lobeDir, ['radii_folds' dvexten '_avgpts*.png'])) ;
if (length(fold_ofn) == length(timePoints)) || overwrite_foldims
    if save_ims
        disp('Plotting ss folds...')
        aux_plot_folds(folds, ssfold, ssfold_frac, ssmax, rmax, nU, ...
            xp.fileMeta.timePoints, lobeDir, dvexten, spcutMeshBase, ...
            'avgpts', overwrite_foldims)
    end
end
fold_ofn = dir(fullfile(lobeDir, ['radii_folds' dvexten '_avgpts*.png'])) ;
if (length(fold_ofn) == length(timePoints)) || overwrite_foldims
    if save_ims
        disp('Plotting rss folds...')
        aux_plot_folds(folds, rssfold, rssfold_frac, rssmax, rmax, nU, ...
            xp.fileMeta.timePoints, lobeDir, dvexten, spcutMeshBase, ...
            'ringpath', overwrite_foldims)
    end
end
disp('done')

%% RECOMPUTE WRITHE OF MEANCURVE CENTERLINES ==============================
% First compute using the avgpts (DVhoop means)
disp('Computing/Loading writhe...')
omit_endpts = 4 ;
wrfn = fullfile(lobeDir, ['writhe_sphi' dvexten '_avgpts.mat']) ;
if ~exist(wrfn, 'file') || overwrite_writhe
    [Wr, Wr_density, dWr, Length_t, clines_resampled] = ...
        aux_compute_writhe(clineDVhoopBase, xp.fileMeta.timePoints, true, omit_endpts, preview) ;
    
    % Save the fold locations as a mat file
    save(wrfn, 'Wr', 'Wr_density', 'dWr', 'Length_t', 'clines_resampled')
else
    load(wrfn, 'Wr', 'Wr_density', 'dWr', 'Length_t', 'clines_resampled')
end
Wr_style = 'Levitt' ;
tmpfn = fullfile(writheDir, ['writhe_' Wr_style '_vs_time_comparison_DVhoop.png']) ;
if ~exist(tmpfn, 'file') || overwrite_writhe
    % Compute ringpath pathlength for results found using centerline
    area_volume_fn = fullfile(meshDir, 'surfacearea_volume_stab.mat') ;
    aux_plot_writhe(xp.fileMeta.timePoints, clines_resampled, ...
        Wr, Wr_density, dWr, Length_t, writheDir, area_volume_fn, ...
        fold_onset, Wr_style, xyzlim, clineDVhoopBase, ...
        cylinderMeshCleanBase, rot, trans, resolution, omit_endpts, false)
    
end
disp('done')

%% Compute surface area and volume for each compartment ===================
lobe_dynamics_fn = fullfile(lobeDir, ['lobe_dynamics' dvexten '.mat']) ;
if exist(lobe_dynamics_fn, 'file') && ~overwrite_lobedynamics
    % Load length, surface area, and volume dynamics for each lobe
    disp('Loading lobe length, area, and volume...')
    load(lobe_dynamics_fn, 'length_lobes', 'area_lobes', 'volume_lobes')
    tp = xp.fileMeta.timePoints - min(fold_onset) ;
else
    [length_lobes, area_lobes, volume_lobes] = aux_compute_lobe_dynamics(folds, ssfold, ssmax, lobeDir, xp.fileMeta.timePoints, ...
        spcutMeshBase, nV, nU, rot, trans, resolution, xyzlim, colors, save_ims, overwrite_lobeims) ;
    % Save surface area and volume dynamics for each lobe
    save(fullfile(lobeDir, ['lobe_dynamics' dvexten '.mat']), ...
        'length_lobes', 'area_lobes', 'volume_lobes')
end
clearvars fold_ofn

%% plot length, area, and volume for each lobe ============================
lobe_dynamics_figfn = fullfile(lobeDir, ['lobe_dynamics' dvexten '.png']) ;
fig1exist = exist(lobe_dynamics_figfn, 'file') ;
% scaled version of same plot
lobe_dynamics_figfn = fullfile(lobeDir, ['lobe_dynamics' dvexten '_scaled.png']) ;
fig2exist = exist(lobe_dynamics_figfn, 'file') ;
if save_ims && (~fig2exist || ~fig2exist || overwrite_lobeims)
    disp('Plotting lobe dynamics...')
    aux_plot_lobe_dynamics(length_lobes, area_lobes, volume_lobes, ...
            timePoints, fold_onset, colors, lobe_dynamics_figfn)
else
    disp('Skipping lobe dynamics plot since it exists...')
end

%% Plot motion of avgpts at folds in yz plane over time ===================
aux_plot_avgptcline_lobes(folds, fold_onset, lobeDir, dvexten, save_ims, ...
    overwrite_lobeims, tp, timePoints, spcutMeshBase, clineDVhoopBase)
disp('done')

%% Plot motion of DVhoop at folds in yz plane over time ===================
aux_plot_constriction_DVhoops(folds, fold_onset, foldHoopImDir, dvexten, save_ims, ...
    overwrite_lobeims, tp, timePoints, spcutMeshBase, alignedMeshBase, ...
    normal_shift, rot, trans, resolution, colors, xyzlim)
disp('done')

%% SMOOTH MEAN CENTERLINE RADIUS ==========================================
% todo: rework this section so that it takes place after smoothing meshes
aux_smooth_avgptcline_radius_before_mesh_smoothing(overwrite_spcutMesh_smoothradii, ...
    timePoints, spcutMeshBase, clineDVhoopBase, radiusImDir, ...
    rot, trans, resolution, xyzlim, nU, nV)

%% Smooth the sphi grid meshes in time ====================================
% Check if all already exist
redo_meshsmooth = overwrite_spcutMeshSm ;
qq = 1 ;
disp('Checking if smoothed meshes already exist on file')
while ~redo_meshsmooth && qq < length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(qq) ;
    % Check that all meshes are saved--if not, declare that we must compute
    smfn = sprintf(spcutMeshSmBase, t) ;
    smrsfn = sprintf(spcutMeshSmRSBase, t) ;
    smrscfn = sprintf(spcutMeshSmRSCBase, t) ;
    if ~exist(smfn, 'file') || ~exist(smrsfn, 'file') || ~exist(smrscfn, 'file')
        redo_meshsmooth = true ;
    end
    qq = qq + 1 ;
end

if redo_meshsmooth
    disp('Mesh smoothing does not exist on file, computing')
    % Load all spcutMesh objects 
    timePoints = xp.fileMeta.timePoints ;
    vM = zeros(length(timePoints), nU*nV, 3);
    nsmM = zeros(length(timePoints), nU*nV, 3) ;
    for i = 1:length(timePoints)
        t = timePoints(i) ;
        % Load the spcutMesh for this timepoint
        disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
        load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
        vM(i, :, :) = spcutMesh.v ;
    end
    disp('built v3d matrix')
    % Filter in time axis
    disp('Building tripulse filter equivalent to tripuls(-0.5:0.1:0.5)')
    tripulse = 0:0.2:1 ;
    tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
    tripulse = tripulse ./ sum(tripulse(:)) ;
    tripulse = reshape(tripulse, [length(tripulse), 1]) ;
    % linfilt = 0.1 * ones(10, 1, 1) ;
    % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
    v3dsmM = imfilter(vM, tripulse, 'replicate') ;
    % vsmM = permute(vsmM, [2,1,3]) ;
    % nsmM = permute(nsmM, [2,1,3]) ;

    % Alternative is to use Gaussian filter but seems no padding option here:
    % smoothdata(vM, 1, 'gaussian', 10)
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the smoothed meshes, then smoothed/rotated/scaled meshes 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for qq = 1:length(xp.fileMeta.timePoints)
        t = xp.fileMeta.timePoints(qq) ;

        % First build normals from smoothed vertices as facenormals
        vqq = squeeze(v3dsmM(qq, :, :)) ;
        nqq = per_vertex_normals(vqq, spcutMesh.f, 'Weighting', 'angle') ;
        % pack it into a struct mesh
        mesh.v = vqq ;
        mesh.vn = nqq ;
        mesh.f = spcutMesh.f ;

        % Average normals with neighboring normals
        disp('Averaging normals with neighboring normals')
        nqq = average_normals_with_neighboring_vertices(mesh, 0.5) ;
        nsmM(qq, :, :) = nqq ;

        % Next save the mesh
        smfn = sprintf(spcutMeshSmBase, t) ;
        smrsfn = sprintf(spcutMeshSmRSBase, t) ;
        smrscfn = sprintf(spcutMeshSmRSCBase, t) ;
        if ~exist(smfn, 'file') || ~exist(smrsfn, 'file') || ~exist(smrscfn, 'file') || overwrite_spcutMeshSm
            vqq = squeeze(v3dsmM(qq, :, :)) ;
            nqq = squeeze(nsmM(qq, :, :)) ;
            nsmM(qq, :, :) = nqq ./ vecnorm(nqq, 2, 2) ;

            % rotate and scale
            nqqrs = (rot * nqq')' ;
            vqqrs = ((rot * vqq')' + trans) * resolution;

            spcutMeshSm.f = spcutMesh.f ;
            spcutMeshSm.v = vqq ;
            spcutMeshSm.vn = nqq ;
            spcutMeshSm.u = spcutMesh.sphi ;
            spcutMeshSm.nU = spcutMesh.nU ;
            spcutMeshSm.nV = spcutMesh.nV ;
            spcutMeshSm.pathPairs = spcutMesh.pathPairs ;

            % Resave s,phi and their 3D embedding
            disp(['Saving ' sprintf(spcutMeshSmBase, t)])
            save(sprintf(spcutMeshSmBase, t), 'spcutMeshSm') ;

            % Also save rotated and scaled (RS) copy of the time-smoothed mesh
            spcutMeshSmRS = spcutMeshSm ;
            spcutMeshSmRS.v = vqqrs ;
            spcutMeshSmRS.vn = nqqrs ;
            % Resave s,phi and their 3D embedding
            save(sprintf(spcutMeshSmRSBase, t), 'spcutMeshSmRS') ;
            clearvars vqq vqqrs

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CLOSED & GLUED SMOOTHED MESHES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To close the mesh, do the following:
            tmp = spcutMeshSmRS ;
            tmp.u = spcutMesh.sphi ;
            spcutMeshSmRSC = glueCylinderCutMeshSeam(tmp) ;
            save(sprintf(spcutMeshSmRSCBase, t), 'spcutMeshSmRSC') ;

            % check it
            % fig = figure ;
            % % triplot(spcutMeshSmRSC.f, spcutMeshSmRSC.u(:, 1), spcutMeshSmRSC.u(:, 2))
            % trisurf(spcutMeshSmRSC.f, spcutMeshSmRSC.v(:, 1),...
            %     spcutMeshSmRSC.v(:, 2), spcutMeshSmRSC.v(:, 3))
            % To fully close the mesh use:
            % anewpt = mean(spcutMeshSmRS.v(1:nV:end, :), 1)
            % pnewpt = mean(spcutMeshSmRS.v(nU:nV:end, :), 1)
            % waitfor(fig)
        end
    end
else
    disp('Mesh smoothing already exists on file, loading...')
    timePoints = xp.fileMeta.timePoints ;
    v3dsmM = zeros(length(timePoints), nU*nV, 3);
    nsmM = zeros(length(timePoints), nU*nV, 3) ;
    
    % Load each mesh into v3dsmM and nsmM    
    for qq = 1:length(xp.fileMeta.timePoints)
        t = xp.fileMeta.timePoints(qq) ;
        load(sprintf(spcutMeshSmBase, t), 'spcutMeshSm') ;
        v3dsmM(qq, :, :) = spcutMeshSm.v ;
        nsmM(qq, :, :) = spcutMeshSm.vn ;
    end
end
disp('done smoothing meshes in time')

%% Plot the time-smoothed meshes
pdir = ensureDir(fullfile(sphiSmRSPhiImDir, 'perspective')) ;
ddir = ensureDir(fullfile(sphiSmRSPhiImDir, 'dorsal')) ;
vdir = ensureDir(fullfile(sphiSmRSPhiImDir, 'ventral')) ;
ldir = ensureDir(fullfile(sphiSmRSPhiImDir, 'latL')) ;
rdir = ensureDir(fullfile(sphiSmRSPhiImDir, 'latR')) ;
for qq = 1:length(timePoints)
    t = timePoints(qq) ;
    disp(['checking figures for time-smoothed meshes: t = ' num2str(t)])
    % prep directories
    fig1fn = fullfile(sphiSmRSImDir, [sprintf('%04d', t ) '.png']) ;
    fp0fn = fullfile(pdir, [sprintf('%04d', t ) '.png']) ;
    fp1fn = fullfile(ddir, ['dorsal_' sprintf('%04d', t ) '.png']) ;
    fp2fn = fullfile(vdir, ['ventral_' sprintf('%04d', t ) '.png']) ;
    fp3fn = fullfile(ldir, ['latL_' sprintf('%04d', t ) '.png']) ;
    fp4fn = fullfile(rdir, ['latR_' sprintf('%04d', t ) '.png']) ;
    e0 = ~exist(fig1fn, 'file') ;
    e1 = ~exist(fp0fn, 'file') ;
    e2 = ~exist(fp1fn, 'file') || ~exist(fp2fn, 'file') ;
    e3 = ~exist(fp3fn, 'file') || ~exist(fp4fn, 'file') ;
    if e0 || e1 || e2 || e3 || overwrite_SmRSIms
        disp(['Consider smoothed sphi gridmesh for time ' num2str(t)])
        % Load the spcutMesh for this timepoint
        vqq = squeeze(v3dsmM(qq, :, :)) ;
        nqq = squeeze(nsmM(qq, :, :)) ;

        % rotate and scale
        nqqrs = (rot * nqq')' ;
        vqq = ((rot * vqq')' + trans) * resolution;
    end
    
    % Plot embedding RS colored by normal vector in y direction
    if e0 || overwrite_SmRSIms
        fig = figure('visible', 'off') ;
        % Color figure with normals taken via faceNormals  
        trisurf(spcutMeshSm.f, vqq(:, 1), vqq(:, 2), vqq(:, 3), ...
            nqqrs(:, 2), 'EdgeColor', 'none')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        title(['Smoothed mesh embedding, t=' num2str(t)])
        disp(['Saving figure: ' fig1fn])
        saveas(fig, fig1fn)
        close all
    end
    
    % Plot embedding colored by phi
    if e1 || e2 || e3 || overwrite_SmRSIms 
        fig = figure('visible', 'off') ;
        % Plot embedding in 3d color coding (color=phi)
        pc = (1:nV) .* ones(nU, nV) / nV ;
        trisurf(spcutMeshSm.f, vqq(:, 1), vqq(:, 2), vqq(:, 3), ...
            pc(:), 'EdgeColor', 'none')
        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        zlabel('z [\mum]')
        c = colorbar() ;
        c.Label.String = '\phi / 2\pi' ;
        title(['Smoothed mesh embedding, t=' num2str(t)])
        disp(['Saving figure: ' fp0fn])
        saveas(fig, fp0fn)
        view(0, 90)
        saveas(fig, fp1fn)
        view(0, 270)
        saveas(fig, fp2fn)
        view(0, 0)
        saveas(fig, fp3fn)
        view(0, 180)
        saveas(fig, fp4fn)
        close all
    end
end
disp('done')
clearvars fig vqq nqqrs e0 e1 e2 e3 pdir ddir vdir rdir ldir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE MEAN AND GAUSSIAN CURVATURES OF SMOOTHED MESHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getderivatives=0;
disp('Computing/Loading Curvatures...')
for t = xp.fileMeta.timePoints
    Kfns = {fullfile(KSmDir, 'latleft', sprintf('gausscurv_latleft_%06d.png', t)), ...
        fullfile(KSmDir, 'dorsal', sprintf('gausscurv_dorsal_%06d.png', t)), ...
        fullfile(KSmDir, 'latright', sprintf('gausscurv_latright_%06d.png', t)), ...
        fullfile(KSmDir, 'ventral', sprintf('gausscurv_ventral_%06d.png', t)) };
    Hfns = {fullfile(HSmDir, 'latleft', sprintf('meancurv_latleft_%06d.png', t)), ...
        fullfile(HSmDir, 'dorsal', sprintf('meancurv_dorsal_%06d.png', t)), ...
        fullfile(HSmDir, 'latright', sprintf('meancurv_latright_%06d.png', t)),...
        fullfile(HSmDir, 'ventral', sprintf('meancurv_ventral_%06d.png', t))} ;
    KHfn = fullfile(KHSmDir, sprintf('gauss_mean_curvature_%06d.mat', t)) ;
    if ~exist(KHfn, 'file') || overwrite_curvatures
        % load glued / closed mesh for curvature computation
        load(sprintf(spcutMeshSmRSCBase, t), 'spcutMeshSmRSC') ;
        fv.faces = spcutMeshSmRSC.f ;
        fv.vertices = spcutMeshSmRSC.v ;
        PrincipalCurvatures = GetCurvatures( fv, getderivatives);
        gaussCurv = (PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:))';
        meanCurv = -0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';

        opts.outfn = Kfns ;
        opts.label = 'Gaussian curvature, $K$' ;
        opts.sscale = 0.01 ;
        opts.cbarPosition = [.85 .333 .02 .333] ;
        opts.xlim = xyzlim(1, :) ;
        opts.ylim = xyzlim(2, :) ;
        opts.zlim = xyzlim(3, :) ;
        scalarFieldOnSurface(fv.faces, fv.vertices, gaussCurv, opts)
        opts.outfn = Hfns ;
        opts.sscale = 0.3 ;
        opts.label = 'mean curvature, $H$' ;
        scalarFieldOnSurface(fv.faces, fv.vertices, meanCurv, opts)
        % Save the curvatures as a mat file
        save(KHfn, 'meanCurv', 'gaussCurv')
    else
        disp('Curvatures already computed & saved on disk')
    end
end
error('here after curvature calc')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redo Pullbacks with time-smoothed meshes ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Create pullback using S,Phi coords with time-averaged Meshes \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for qq = 1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(qq) ;
    disp(['checking pullbacks for time-smoothed meshes: t = ' num2str(t)])
    
    % Load time-smoothed mesh
    load(sprintf(spcutMeshSmBase, t), 'spcutMeshSm') ;
    
    %--------------------------------------------------------------
    % Generate Output Image File
    %--------------------------------------------------------------
    imfn_spsm = sprintf( fullfile( imFolder_spsm, [fileNameBase, '.tif']), t ) ;
    imfn_rsm  = sprintf( fullfile( imFolder_rsm, [fileNameBase, '.tif']), t ) ;
    pullbacks_exist1 = exist(imfn_spsm, 'file') ;
    pullbacks_exist2 = exist(imfn_rsm, 'file') ;
    if ~pullbacks_exist1 || ~pullbacks_exist2 || overwrite_pullbacks
        % Load 3D data for coloring mesh pullback
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();

        % Raw stack data
        IV = xp.stack.image.apply();
        IV = imadjustn(IV{1});
    end

    if ~exist(imfn_spsm, 'file') || overwrite_pullbacks
        fprintf(['Generating SP output image for sm mesh: ' imfn_spsm]);
        % Assigning field spcutMesh.u to be [s, phi] (ringpath
        % and azimuthal angle)
        aux_generate_orbifold( spcutMeshSm, a_fixed, IV, imfn_spsm)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save relaxed image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if ~exist(imfn_rsm, 'file') || overwrite_pullbacks
        disp('Generating relaxed image for sm mesh...')
        tmp = spcutMesh.sphi ;
        tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
        arspsm = minimizeIsoarealAffineEnergy( spcutMeshSm.f, spcutMeshSm.v, tmp );
        aux_generate_orbifold(spcutMeshSm, arspsm, IV, imfn_rsm)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TILE SMOOTHED IMAGES IN Y AND RESAVE ===================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDirs = {imFolder_spsm, imFolder_rsm, } ;
imDirs_e = {imFolder_spsm_e, imFolder_rsm_e, } ;
ntiles = 50 ;   %   The number of bins in each dimension for histogram equilization for a
                %   square original image. That is, the extended image will have (a_fixed *
                %   ntiles, 2 * ntiles) bins in (x,y).
for qq = 1:length(imDirs)
    options.histeq = true ;
    options.a_fixed = a_fixed ;
    options.ntiles = ntiles ;
    options.overwrite = overwrite_pullbacks;
    extendImages(imDirs{qq}, imDirs_e{qq}, fileNameBase, options)
    disp(['done ensuring extended tiffs for ' imDirs{qq} ' in ' imDirs_e{qq}])
end
clearvars qq ntiles imDirs_e imDirs
disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEASURE THICKNESS ======================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure thickness from images of imFolder_spsm_e2 (extendedLUT_smoothed)
fprintf('Create pullback nstack using S,Phi coords with time-averaged Meshes \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer_spacing = 1 / resolution ;
n_outward = 10 ;
n_inward = 8 ;
for qq = 1:length(xp.fileMeta.timePoints)
    t = xp.fileMeta.timePoints(qq) ;
    disp(['t = ' num2str(t)])
    
    % Load time-smoothed mesh
    load(sprintf(spcutMeshSmBase, t), 'spcutMeshSm') ;
    
    %--------------------------------------------------------------
    % Generate Output Image File
    %--------------------------------------------------------------
    spacingstr = strrep(sprintf('%0.2fum', layer_spacing * resolution), '.', 'p') ;
    imfn_spsm = sprintf( fullfile( imFolder_spsm_e2, ...
        [fileNameBase, '_%02d_%02d_' spacingstr '.tif']), t, n_outward, n_inward ) ;
    
    if ~exist(imfn_spsm, 'file') || overwrite_pullbacks
        % Load 3D data for coloring mesh pullback
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();

        % Raw stack data
        IV = xp.stack.image.apply();
        IV = imadjustn(IV{1});
        
        fprintf(['Generating SP output image for sm mesh: ' imfn_spsm]);
        % Assigning field spcutMesh.u to be [s, phi] (ringpath
        % and azimuthal angle)
        Options.numLayers = [n_outward, n_inward] ;
        Options.layerSpacing = layer_spacing ;
        Options.smoothIter = 1 ;
        Options.yLim = [-0.5, 1.5] ;
        % Note that we pass a_fixed * 0.5 since the image is extended by a
        % factor of two
        aux_generate_orbifold( spcutMeshSm, a_fixed * 0.5, IV, imfn_spsm, Options)
    end
    clear Options
end

%% TRAIN IN ILASTIK ON MIDGUT TISSUE TO GET THICKNESS
% Read thickness training output


%% DUMP OR LOAD HERE [break]
clearvars fold_ofn fig1exist fig2exist id idx mcline nsmM prevcline tmp
clearvars fig1exist fig1fn fig2exist fp1fn fp2fn fp3fn fp4fn IV 
clearvars n_inward n_outward TV2D TV3D 
clearvars layer_spacing 
dumpfn = fullfile(meshDir, 'orbifold_dump_before_piv.mat') ;
save(dumpfn)
load(dumpfn)
clearvars dumpfn

%% PERFORM PIV ============================================================
% Select all frames in PullbackImages_extended_shifted/ 
% Select Sequencing style 1-2, 2-3, ... 
% Image Preprocessing > Select All
% PIV settings: 128 (32 step), 64 (16 step) for two passes
% Post-processing: vector validation: Select velocity limits, 5 stdev
% filter, local median filter: thres=2, eps=0.1
% with default values.
% File > Save > Export all frames as pivresults.mat
disp('Loading pivresults_sphi_pass0.mat ...')
piv = load(fullfile(pivDir, 'pivresults_sphi_pass0.mat')) ;
disp('done')

%% Get timestamps for the images with pullbacks
fns = dir(fullfile(imFolder_e, '*.tif')) ;
npiv = length(fns) ;
time = zeros(npiv, 1);
meshidx = zeros(npiv, 1) ;
for i=1:npiv
    tmp = split(fns(i).name, '.tif') ;
    tmp = split(tmp{1}, '_') ;
    time(i) = str2double(tmp{2}) ;
    meshidx(i) = time(i) - time(1) + 1 ;
end
dt = diff(time) ;
disp('done building dt, meshidx') 
save(fullfile(pivDir, 'timestamps_orbifold.mat'), 'time', 'meshidx')

%% Subtract off the mean flow in y for each frame ========================= 
% meanv = zeros(length(v_filtered), 1) ;
% for i=1:length(v_filtered)
%     tmp = v_filtered{i} ;
%     meanv(i) = mean(tmp(:)) ;
% end
% dy = round(meanv) ;
% shifty = [0; -cumsum(dy)] ;
% disp('computed shifty')
% save(fullfile(pivDir, 'shifty.mat'), 'shifty')

%% Shift each frame by shifty and resave ==================================
% fns = dir(strrep(fullfile([imFolder_sp_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
% imFolder_es = [imFolder_sp_e '_shifted' filesep] ;
% if ~exist(imFolder_es, 'dir')
%     mkdir(imFolder_es) ;
% end
% 
% % Write the images to disk and get their sizes while being written 
% disp('Resaving shifted images if necessary...')
% imsizes = zeros(length(fns), 2) ;
% for i=1:length(fns)
%     outfn = fullfile(imFolder_sp_es, fns(i).name) ;
%     if ~exist(outfn, 'file')
%         disp(['Reading ' fns(i).name])
%         fileName = split(fns(i).name, '.tif') ;
%         fileName = fileName{1} ;
%         im = imread(fullfile(fns(i).folder, fns(i).name)) ;
%         % Note adaptive histogram equilization
%         im = adapthisteq(circshift(im, shifty(i), 1)) ;
%         imsizes(i, :) = size(im) ;
%         imwrite( im, outfn, 'TIFF' );
%     else
%         im = imread(outfn) ;
%         imsizes(i, :) = size(im) ;
%     end
% end
% disp('done with shifted images')

%% PERFORM PIV ON SHIFTED FRAMES ==========================================
% % Select all frames in PullbackImages_extended_shifted/ 
% % Select Sequencing style 1-2, 2-3, ... 
% % Load settings: piv_set_pass1.mat
% % Image Preprocessing > Select All
% % PIV settings: 128 (32 step), 64 (16 step), 32 (16 step) for three passes
% disp('Loading pivresults_orbifold_pass1.mat ...')
% load(fullfile(pivDir, 'pivresults_orbifold_pass1.mat'))


%% MAKE MAP FROM PIXEL TO XYZ =============================================
disp('Making map from pixel to xyz to compute velocities in 3d for smoothed meshes...')
% Get position of embedding points associated with velocities 
pivOutDir = fullfile(pivDir, 'images') ;
if ~exist(pivOutDir, 'dir')
    mkdir(pivOutDir)
end
use_shifty = false ;
show_v3d_on_data = false ;

% Compute size of domain (0,1), (0, 1)
if use_shifty
    fns = dir(strrep(fullfile([imFolder_spsm_es, '/', fileNameBase, '.tif']), '%06d', '*')) ;
else
    fns = dir(strrep(fullfile([imFolder_spsm_e, '/', fileNameBase, '.tif']), '%06d', '*')) ;
end

im = imread(fullfile(fns(i).folder, fns(i).name)) ;
% size of extended image
esize = size(im) ;
% map extended image size to (0, 1), (-0.5, 1.5)
xesz = esize(2) ;
yesz = esize(1) ;
% map from pixel y to network y (sphi)
Ypix2y = @(ypix) 2.0 * ypix / yesz - 0.5 ;
Xpix2x = @(xpix) 2.0 * xpix / yesz ;
% map from network xy to pixel xy
% Note that we need to flip the image (Yscale - stuff) since saved ims had
% natural ydirection.
% Assume here a coord system xy ranging from (0, xscale) and (0, 1) 
% maps to a coord system XY ranging from (0, Yscale * 0.5) and (0, Yscale)
x2Xpix = @(x, Yscale, xscale) (Yscale * 0.5) * a_fixed * x / xscale ;
dx2dX = @ (y, Yscale, xscale) (Yscale * 0.5) * a_fixed * x / xscale ;
dy2dY = @ (y, Yscale) (Yscale*0.5)*y ;
if use_shifty
    % y2Ypix = @(y, h, Yscale) Yscale - (Yscale*0.5)*(y+0.5) + h ;
    y2Ypix = @(y, h, Yscale) (Yscale*0.5)*(y+0.5) + h ;
else
    % y2Ypix = @(y, Yscale) Yscale - (Yscale*0.5)*(y+0.5) ;
    y2Ypix = @(y, Yscale) (Yscale*0.5)*(y+0.5) ;
    shifty = [] ;
end

% piv3d name base for loading each timepoint separately
piv3dfn = fullfile(fullfile(pivDir, 'piv3d'), 'piv3d_%04d.mat') ;
ensureDir(fullfile(pivDir, 'piv3d')) ;

% Check if all timepoints' piv3d exist already
redo_piv3d = overwrite_piv ; 
i = 1 ;
while ~redo_piv3d && i < length(time)
    redo_piv3d = ~exist(sprintf(piv3dfn, time(i)), 'file') ;
    i = i + 1 ;
end
if ~redo_piv3d
    disp(['Loading piv3d from ' piv3dfn])
    piv3d = cell(length(fns)) ;
    for i=1:(length(fns)-1)
        if mod(time(i), 10) == 0
            disp(['loading piv3d for time t=' num2str(time(i))])
        end
        load(sprintf(piv3dfn, time(i)), 'piv3dstruct')
        piv3d{i} = piv3dstruct ;
    end
else
    disp('Overwriting/computing piv3d data...')
    piv3d = cell(length(fns), 1) ;
    
    % Iterate over all images with flow fields ----------------------------
    for i=1:(length(fns) - 1)
        timestr = sprintf('%03d', time(i)) ;
        disp(['t = ' timestr])
        t = time(i) ;

        % Get scale of image
        Xsc0 = xesz ;
        Ysc0 = yesz ;  % imsizes(i, 2) ;
        Ysc1 = yesz ;  % imsizes(i+1, 2) ;

        % Load spcutMesh
        mesh0 = load(sprintf(spcutMeshSmBase, time(i)), 'spcutMeshSm') ;
        mesh0 = mesh0.spcutMeshSm ;
        mesh1 = load(sprintf(spcutMeshSmBase, time(i + 1)), 'spcutMeshSm') ;
        mesh1 = mesh1.spcutMeshSm ;
        assert(time(i) + dt(i) == time(i+ 1))
        xsc0 = max(mesh0.u(:, 1)) ;
        xsc1 = max(mesh1.u(:, 1)) ;
        mesh0x = x2Xpix(mesh0.u(:, 1), Ysc0, xsc0) ;
        mesh0y = y2Ypix(mesh0.u(:, 2), Ysc0) ;
        mesh1x = x2Xpix(mesh1.u(:, 1), Ysc1, xsc1) ;
        mesh1y = y2Ypix(mesh1.u(:, 2), Ysc1) ;

        % Load the positions of the velocity vectors in pixels
        x0 = piv.x{i} ;
        y0 = piv.y{i} ;
        uu = piv.u_filtered{i} ;
        vv = piv.v_filtered{i} ; 
        % Get position in next timepoint in pixels (in xy plane)
        % Clip the x position to the size of the image, and wrap the y position
        x1 = x0 + uu ;
        eps = 1e-8 ;
        x1 = max(x1, eps) ;
        x1 = min(x1, xesz-eps) ;
        y1 = mod(y0 + vv, Ysc1) ;  

        % Two methods of obtaining the 3d vel evaluation pts. 
        % Barycenteric version:
        % [pt0, pt1, tr0, tr0_orig, tr1, tr1_orig, tria] = aux_barycentricInterpTiledMesh(mesh0, mesh1,...
        %     Ysc0, xsc0, Ysc1, xsc1, use_shifty, shifty, i, x0, y0, x1, y1, x2Xpix, y2Ypix, dx2dX, dy2dY) ;
        % mf0c = tr0.ConnectivityList ;
        % m0xy = tr0.Points ;
        % mesh0x = tr0_orig.Points(:, 1) ;
        % mesh0y = tr0_orig.Points(:, 2) ;
        % mesh1x = tr1_orig.Points(:, 1) ;
        % mesh1y = tr1_orig.Points(:, 2) ;

        % % Instead, do interpolation:
        disp('Interpolating 3d vertices for tiled mesh 0')
        tileCount = [2, 2] ;
        if any(isnan(mesh0.v))
            % fill missing x,y,z values in mesh0.v
            for dim = 1:3
                dgrid = reshape(mesh0.v(:, dim), [mesh0.nU, mesh0.nV]) ;
                dgrid = fillmissing(dgrid, 'linear') ;
                mesh0.v(:, dim) = dgrid(:) ;
            end
            clearvars dgrid
        end
        [ tm0f, tm0v2d, tm0v3d, ~ ] = tileAnnularCutMesh( mesh0, tileCount );
        tm0X = x2Xpix(tm0v2d(:, 1), Ysc0, xsc0) ;
        tm0Y = y2Ypix(tm0v2d(:, 2), Ysc0) ;
        tm0XY = [tm0X, tm0Y] ;
        [pt0, fieldfaces, tr0] = interpolate2Dpts_3Dmesh(tm0f, tm0XY, tm0v3d, [x0(:), y0(:)]) ;
        
        % note: for interpolation of velocity onto face centers, use
        m0X = x2Xpix(mesh0.u(:, 1), Ysc0, xsc0) ;
        m0Y = y2Ypix(mesh0.u(:, 2), Ysc0) ;
        m0XY = [m0X, m0Y] ;
        mesh_for_interp.pathPairs = mesh0.pathPairs ;
        mesh_for_interp.u = m0XY ;
        mesh_for_interp.f = mesh0.f ;
        
        % If any fieldfaces are NaN, jitter the points a bit
        if any(isnan(fieldfaces))
            % Fill in missing data. Will this work?
            facegrid = reshape(fieldfaces, size(x0)) ;
            facegrid = fillmissing(facegrid, 'linear') ;
            fieldfaces = uint8(facegrid(:)) ;
            
            if preview
                % which are the bad ones?
                badID = find(isnan(fieldfaces)) ;

                % check it in a plot
                triplot(tr0.ConnectivityList, tr0.Points(:, 1), tr0.Points(:, 2), 0*tr0.Points(:, 2), 'Edgecolor', 'none')
                hold on;
                bd = freeBoundary(tr0) ;
                plot3(tr0.Points(bd(:, 1), 1), tr0.Points(bd(:, 1), 2), 0*tr0.Points(bd(:,1), 1), '.')
                hold on;
                plot3(x0(badID), y0(badID), 0*x0(badID), 'o')
                xlim([min(x0(badID)) - 1, min(x0(badID)) + 1])
                title('NaNs in fieldfaces')
                waitfor(fig)
            end
        end
                
        % Interpolate for next timepoint
        disp('Interpolating 3d vertices for tiled mesh 1')
        tileCount = [2, 2] ;
        [ tm1f, tm1v2d, tm1v3d, ~] = tileAnnularCutMesh( mesh1, tileCount );
        tm1XY = [x2Xpix(tm1v2d(:, 1), Ysc1, xsc1), y2Ypix(tm1v2d(:, 2), Ysc1)] ;
        pt1 = interpolate2Dpts_3Dmesh(tm1f, tm1XY, tm1v3d, [x1(:), y1(:)]) ;
        % Old version
        % Xbi = scatteredInterpolant(tm1X, tm1Y, tm1v3d(:, 1)) ;
        % Ybi = scatteredInterpolant(tm1X, tm1Y, tm1v3d(:, 2)) ;
        % Zbi = scatteredInterpolant(tm1X, tm1Y, tm1v3d(:, 3)) ;
        % pt1 = [Xbi(x1(:), y1(:)), Ybi(x1(:), y1(:)), Zbi(x1(:), y1(:))] ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (preview || save_ims) && show_v3d_on_data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot the advected mesh
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load 3D data for coloring mesh pullback
            disp('loading timepoint data')
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();
            % Raw stack data
            IV0 = xp.stack.image.apply();
            IV0 = imadjustn(IV0{1});
            
            % also plot the next timepoint
            t1 = t + 1 ;
            disp('loading subsequent timepoint data')
            xp.loadTime(t1);
            xp.rescaleStackToUnitAspect();
            % Raw stack data
            IV1 = xp.stack.image.apply();
            IV1 = imadjustn(IV1{1});
        end
        
        if preview
            disp('previewing')
            % evaluate flow field at (x, Ysz - y). Add (u,v) to (x, Ysz - y).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot the advected mesh
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot the advected tissue on top of the next frame
            % interpolate velocities (x0, y0) onto mesh (mesh0x, mesh0y)
            uuinterp = griddedInterpolant(x0', y0', uu') ; % in inverted Y
            vvinterp = griddedInterpolant(x0', y0', vv') ; 
            umeshpix = uuinterp(mesh0x, mesh0y) ; % interpolate in correct Y
            vmeshpix = vvinterp(mesh0x, mesh0y) ;
            addx = umeshpix ;
            addy = vmeshpix ;
            mesh0adv_pix = [ mesh0x + addx, mesh0y + addy] ;

            % Texture image options
            Options.imSize = ceil( Xsc0 .* [ a_fixed 1 ] );  % row, col pix
            Options.yLim = [0 Ysc0];
            Options.xLim = [0 Xsc0];
            % original image RED
            im0 = texture_patch_to_image(tm0f, [tm0X, tm0Y], ...
                tm0f, tm0v3d(:, [ 2 1 3]), IV0, Options );
            % im0 = imread(fullfile(fns(i).folder, fns(i).name)) ;
            % advected image GREEN
            
            patchIma = texture_patch_to_image(mesh0.f, mesh0adv_pix, ...
                mesh0.f, mesh0.v(:, [ 2 1 3]), IV0, Options );
            % patchIma = adapthisteq(patchIma, 'NumTiles', [round(a_fixed * ntiles), round(2 * ntiles)]) ;
            % Next timepoint BLUE
            im1 = texture_patch_to_image(tm1f, [tm1X, tm1Y], ...
                tm1f, tm1v3d(:, [ 2 1 3]), IV1, Options );
            % Load next timepoint BLUE
            % im1 = imread(fullfile(fns(i+1).folder, fns(i+1).name)) ;
            % patchImRGB = cat(3, im0, uint8(patchIma * 255), im1) ;
            patchImRGB = cat(3, im0, patchIma, im1) ;

            % Plot the image and advected image, in inverted Y space but
            % with increasing Y going down
            close all; 
            fig1 = figure(1) ; 
            h = imshow( patchImRGB );
            hold on;  
            quiver(mesh0x, mesh0y, addx, addy, 0, 'color', yellow, 'linewidth', 2)
            % plot(mesh0x, mesh0y, 'o')
            % plot(mesh0adv_pix(:, 1), mesh0adv_pix(:, 2), 's')
            triplot(mesh0.f, mesh0x, mesh0y, 'color', red, 'linewidth', 2)
            triplot(mesh0.f, mesh0x + addx, mesh0y + addy, 'color', green, 'linewidth', 2)
            % axis equal
            
            % Check the displacement by toggling between frames
            % fig = figure(2) ;
            % pressed_enter = false ; 
            % kk = 0 ;
            % while ~pressed_enter
            %     if kk == 0
            %         imshow(flipud(im0))
            %         %set( gca, 'YDir', 'Normal' );
            %         titlestr = '0: <Enter> to exit, <-> to toggle' ;
            %     else
            %         imshow(flipud(im1))
            %         %set( gca, 'YDir', 'Normal' );
            %         titlestr = '1: <Enter> to exit, <-> to toggle' ;
            %     end
            %     title(titlestr)
            %     was_a_key = waitforbuttonpress;
            %     left = strcmp(get(fig, 'CurrentKey'), 'leftarrow');
            %     rght = strcmp(get(fig, 'CurrentKey'), 'rightarrow') ;
            %     if was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
            %         pressed_enter = true ;
            %     elseif was_a_key && left || rght 
            %         kk = mod(kk + 1, 2) ;
            %     end
            % end

            % Check differences
            d0 = mat2gray(im1, [0, double(max(im1(:)))]) - mat2gray(im0, [0, double(max(im0(:)))]);
            d0pos = 0*d0 ;
            d0pos(d0 > 0) = d0(d0 > 0) ;
            d0neg = 0*d0 ;
            d0neg(d0 < 0) = abs(d0(d0 < 0)) ;
            pd0 = cat(3, d0pos, 0*d0, d0neg) ; 

            % diff between next and advected 
            da = mat2gray(im1, [0, double(max(im1(:)))]) - patchIma ;
            dapos = 0*da ;
            dapos(da > 0) = da(da > 0) ;
            daneg = 0*da ;
            daneg(da < 0) = abs(da(da < 0)) ;
            pda = cat(3, dapos, 0*da, daneg) ;
            % plot both
            fig1 = figure(1) ;
            imshow(pd0) 
            title(['\langle|t1 - t0|\rangle = ' num2str(100 * mean(pd0(:))) '%'])
            waitfor(fig1)

            fig2 = figure(2) ;
            imshow(pda) 
            title(['\langle|t1 - advected t0|\rangle = ' num2str(100 * mean(pda(:))) '%'])
            waitfor(fig2)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % if save_ims
        %     disp('visualize 2d triangulation')
        %     % Check out the triangulation    
        %     close all
        %     fig = figure('Visible', 'Off');
        % 
        %     % second option
        %     % tr0 = triangulation(tm0f, tm0v2d) ;
        %     tr0_orig = triangulation(mesh0.f, [mesh0x, mesh0y]) ;
        % 
        %     hc = triplot(tr0, 'Color', blue) ;
        %     % hold on
        %     ho = triplot(tr0_orig, 'Color', orange) ;
        %     axis equal
        %     title('Triangulation in pullback image space')
        %     ylim([0 yesz])
        %     xlabel('x [pix]')
        %     ylabel('y [pix]')
        %     legend({'tiled', 'original'})
        %     saveas(fig, fullfile(pivOutDir, ['tri_' timestr '.png']))
        %     if preview
        %         set(fig, 'Visible', 'On')
        %         waitfor(fig)
        %     end
        %     close all
        % end
        
        %% Take the difference to get the velocity field ------------------
        v0 = (pt1 - pt0) / dt(i) ;
        
        % Resolve tangential and normal velocities, obtain jacobian =======
        [v0n, v0t, v0t2d, jac, facenormals, g_ab, dilation] = ...
            resolveTangentNormalVelocities(tm0f, tm0v3d, v0, fieldfaces, tm0XY) ;
        
        % I have checked that det(jac * jac') = det(jjac * jjac') for a
        % triangle. 
                
        % todo: check that the dilation = ratio of areas of triangles2d /
        % triangles3d
        
        if preview
            % Checking the dilation
            plot(tm0X, tm0Y, '.')
            hold on;
            scatter(x0(:), y0(:), 10, dilation)
            waitfor(gcf)

            % Now independently measure the areas of triangles in 3d
            % Evaluate for every face in tm0
            fa3d = doublearea(tm0v3d, tm0f) * 0.5 ;
            % also do 2d areas
            fa2d = doublearea(tm0XY, tm0f) * 0.5 ;
            arearatio = fa2d(fieldfaces) ./ fa3d(fieldfaces) ;
            figure;
            plot(dilation, arearatio, '.') ;
            hold on; 
            plot([0, 5], [0, 5], '--')
            xlabel('dilation from jacobian')
            ylabel('area ratio of triangles')
            title('Checking dilation from 3d->2d: should be y=x')
            waitfor(gcf)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign one velocity per face via interpolation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Obtain circumcenters for the current timepoint's 2D mesh
        % tr0 = triangulation(tm0f, tm0XY) ;  <-- already done in interpolate2dpts_3Dmesh
        % bc = circumcenter(tr0) ;

        % Calculate centroid of faces -- here evaluate only on untiled mesh
        bc = cat( 3, mesh_for_interp.u(mesh_for_interp.f(:,1), :), ...
            mesh_for_interp.u(mesh_for_interp.f(:,2), :),...
            mesh_for_interp.u(mesh_for_interp.f(:,3), :) );
        bc = mean( bc, 3 );
                
        % cast as NxMx3, then as N*M x 1 arrays for vx,vy,vz separately
        v3dgrid = reshape(v0, [size(piv.x{1}, 1), size(piv.x{1}, 2), 3]) ;
        xvel = squeeze(v3dgrid(:, :, 1)) ;
        yvel = squeeze(v3dgrid(:, :, 2)) ;
        zvel = squeeze(v3dgrid(:, :, 3)) ;
        % Find the velocity at the advected position using only 2d coords
        Fx = griddedInterpolant(x0', y0', xvel') ;
        Fy = griddedInterpolant(x0', y0', yvel') ;
        Fz = griddedInterpolant(x0', y0', zvel') ;
        v3dfaces = [Fx(bc(:, 1), bc(:, 2)), ...
            Fy(bc(:, 1), bc(:, 2)), Fz(bc(:, 1), bc(:, 2))] ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign one tangential velocity per face via interpolation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v3dvertices = [Fx(mesh_for_interp.u(:, 1), mesh_for_interp.u(:, 2)), ...
            Fy(mesh_for_interp.u(:, 1), mesh_for_interp.u(:, 2)), ...
            Fz(mesh_for_interp.u(:, 1), mesh_for_interp.u(:, 2))] ;
        % %% Convert pivSimAvg to evaluated at vertices of simple smoothed meshes
        % piv3dfn = fullfile(fullfile(pivDir, 'piv3d'), 'piv3d_%04d.mat') ;
        % % load piv3d as a cell array
        % disp(['Loading piv3d from ' piv3dfn])
        % piv3d = cell(length(fns)) ;
        % for i=1:(length(fns)-1)
        %     if mod(time(i), 10) == 0
        %         disp(['loading piv3d for time t=' num2str(time(i))])
        %     end
        %     load(sprintf(piv3dfn, time(i)), 'piv3dstruct')
        %     piv3d{i} = piv3dstruct ;
        % end
        
        %% Save the results in datstruct ----------------------------------
        % v0, v0n, v0t are in units of um/min,  
        % while v0t2d, g_ab, and jacobian are in pullback pixels
        datstruct.m0XY = m0XY ;
        datstruct.m0f = mesh0.f ;
        datstruct.m0v3d = mesh0.v ;
        datstruct.x0 = x0 ;
        datstruct.y0 = y0 ;
        datstruct.pt0 = pt0 ;
        datstruct.pt1 = pt1 ;
        datstruct.v0 = v0 / dt(i) ;
        datstruct.v0n = v0n / dt(i) ;
        datstruct.v0t = v0t / dt(i) ;
        datstruct.facenormals = facenormals ;
        
        % Face-centered velocities
        datstruct.v3dvertices = v3dvertices ;
        datstruct.v3dfaces = v3dfaces ;
        datstruct.v3dfaces_rs = (rot * v3dfaces')' * resolution / dt(i) ; 
        assert(all(size(v3dfaces) == size(datstruct.m0f)))
        
        % rotated and scaled velocities
        datstruct.v0_rs = (rot * v0')' * resolution / dt(i) ;
        datstruct.v0n_rs = v0n * resolution / dt(i) ;
        datstruct.v0t_rs = (rot * v0t')' * resolution / dt(i) ;
        datstruct.normals_rs = (rot * facenormals')' ;
        
        % 2d pullback velocities
        datstruct.v0t2d = v0t2d ;           % in pullback pix / min
        datstruct.g_ab = g_ab ;             % pullback metric
        datstruct.dilation = dilation ;     % dilation of face from 3d to 2d
        datstruct.dilation_rs = dilation / resolution ;     % dilation of face from 3d to 2d
        datstruct.jacobian = jac ;          % jacobian of 3d->2d transformation
        datstruct.fieldfaces = fieldfaces ; % faces where v defined from PIV
        
        %        
        readme = struct ;
        readme.m0XY = "Px2 float: 2d mesh vertices in pullback image pixel space" ;
        readme.m0f = "#facesx3 float: mesh connectivity list" ;
        readme.m0v3d = "Px3 float: 3d mesh coordinates in embedding pixel space" ;
        readme.x0 = "QxR float: x value of 2d velocity evaluation coordinates in pullback image pixel space" ;
        readme.y0 = "QxR float: y value of 2d velocity evaluation coordinates in pullback image pixel space" ;
        readme.pt0 = "Nx3 float: 3d location of evaluation points [mesh pix]" ;
        readme.pt1 = "Nx3 float: 3d location of advected point in next mesh [mesh pix]" ;
        readme.v0 = "Nx3 float: 3d velocity [mesh pix / min]" ;
        readme.v0n = "Nx1 float: normal velocity [mesh pix / min]" ;
        readme.v0t = "Nx3 float: tangential velocity in 3d [mesh pix / min]" ;
        readme.facenormals = "#faces x 3 float: face normals for all faces" ;
        readme.v3dfaces = "#faces x 3 float: face-centered velocities in embedding space" ;
        readme.v3dvertices = "#mesh vertices x 3 float: vertex-centered velocities in embedding space" ;
        readme.v0_rs = "Nx3 float: rotated/scaled 3d velocities [um/min]" ;
        readme.v0n_rs = "Nx3 float: scaled normal velocities [um/min]" ;
        readme.v0t_rs = "Nx3 float: rotated/scaled tangential 3d velocities [um/min]" ;
        readme.normals_rs = "Nx3 float: normal vector of field faces in 3d [unitless]" ;
        readme.v0t2d = "Nx2 float: in-plane velocity in pullback pix / min" ;
        readme.g_ab = "Nx2x2 float: pullback metric" ;
        readme.dilation = "Nx1 float: dilation of face from 3d to 2d (A_2d / A_3d)" ;
        readme.jac = "#faces x 1 cell: jacobian of 3d->2d transformation" ;
        readme.fieldfaces = "Nx1 int: face indices into spcutMesh where v defined" ;
        datstruct.readme = readme ;
        piv3d{i} = datstruct ;

        %% Visualize the flow in 3d --------------------------------------- 
        piv3dimagefn = fullfile(pivOutDir, ['piv3d_' timestr '.png']) ;
        if save_ims && (overwrite_piv || ~exist(piv3dimagefn, 'file'))
            disp('visualize flow in 3d')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            close all
            fig = figure('Visible', 'Off') ;
            hold on            
            % for checking purposes, grab the first few indices
            m0 = mesh0 ;
            % rmIDx = 600:length(m0.v) ;
            % [m0.f, m0.v, oldvIdx] = remove_vertex_from_mesh(m0.f, m0.v, rmIDx) ;
            % m0.vn = m0.vn(oldvIdx, :) ;
            m1 = mesh1 ;
            
            if show_v3d_on_data
                % Instead of using texture_patch_3d
                % create griddedInterpolant 
                IV0interp = griddedInterpolant(single(32768 + 0.5 * IV0), 'cubic') ;
                v0i = IV0interp(m0.v(:, 2), m0.v(:, 1), m0.v(:, 3)) ;
                patch('Faces', m0.f, 'Vertices', m0.v, ...
                    'FaceVertexCData', v0i, 'FaceColor', 'interp', ...
                    'EdgeColor', 'none') ; 
                hold on;
                IV1i = griddedInterpolant(single(32768 - 0.5 * IV1)) ;
                v1i = IV1i(m1.v(:, 2), m1.v(:, 1), m1.v(:, 3)) ;
                patch('Faces', m1.f, 'Vertices', m1.v, ...
                    'FaceVertexCData', v1i, 'FaceColor', 'interp', ...
                    'EdgeColor', 'none') ;
                axis equal
                colormap(bwr)
                % Option 2:
                % least squares vertex baking : intead of evaluating vertices
                % in image volume and interpolating values, baking will
                % minimize an energy functional that makes the intensities
                % closer to the proper values from texture mapping
                % 

                % Option 3 (expensive): texture mapping
                % clearvars Options
                % Options.PSize = 2;
                % Options.EdgeColor = 'none';
                % % Options.Rotation = rot ;
                % % Options.Translation = trans ;
                % texture_patch_3d( m0.f, m0.v, ...
                %     m0.f, m0.v(:, [2 1 3]), 32768 + 0.5 * IV0, Options );
                % hold on            
                % texture_patch_3d( m1.f, m1.v, ...
                %     m1.f, m1.v(:, [2 1 3]), 32768 - 0.5 * IV1, Options );
                % hold on            
                % axis equal
                % nearby = false(size(pt0(:, 1), 1), 1) ;
                % for qq = 1:length(pt0(:, 1))
                %     if min(vecnorm(pt0(qq, :) - m0.v, 2, 2)) < 10
                %         nearby(qq) = true ;
                %     end
                % end
                % nearby = find(nearby) ;
                % quiver3(pt0(nearby, 1), pt0(nearby, 2), pt0(nearby, 3), ...
                %     v0(nearby, 1), v0(nearby, 2), v0(nearby, 3), 0, 'color', green) 
            else
                m0vr = ((rot * m0.v')' + trans) * resolution ;
                m1vr = ((rot * m1.v')' + trans) * resolution ;
                trisurf(m0.f, m0vr(:, 1), m0vr(:, 2), m0vr(:, 3), ...
                    'FaceColor', blue, 'EdgeColor', 'none', 'FaceAlpha', 0.1)
                hold on;
                trisurf(m1.f, m1vr(:, 1), m1vr(:, 2), m1vr(:, 3), ...
                    'FaceColor', red, 'EdgeColor', 'none', 'FaceAlpha', 0.1)
            end
            hold on      
            pt0r = ((rot * pt0')' + trans) * resolution ;
            quiver3(pt0r(:, 1), pt0r(:, 2), pt0r(:, 3), ...
                datstruct.v0_rs(:, 1), ...
                datstruct.v0_rs(:, 2), ...
                datstruct.v0_rs(:, 3), 0, 'color', green)
            
            % scatter3(pt0(:, 1), pt0(:, 2), pt0(:, 3))
            % plot3(pt0(:, 1), pt0(:, 2), pt0(:, 3), 'o', 'color', yellow)
            % plot3(pt1(:, 1), pt1(:, 2), pt1(:, 3), 's', 'color', yellow)
            axis equal
            title(['t=' timestr]) 
            xlabel('x')
            ylabel('y')
            zlabel('z')
            xlim(xyzlim_APDV(1, :))
            ylim(xyzlim_APDV(2, :))
            zlim(xyzlim_APDV(3, :))
            % plot3(m1.v(:, 1), m1.v(:, 2), m1.v(:, 3), '.')
            
            % saveas(gcf, fullfile('/Users/npmitchell/Desktop/tmp/', ['piv3d_' timestr '.png']))
            saveas(gcf, piv3dimagefn)
            close all
        end
        
        %% Draw 2D flows
        vtdir2d = fullfile(pivDir, 'vt2d') ;
        vndir2d = fullfile(pivDir, 'vn2d') ;
        if ~exist(vtdir2d, 'dir')
            mkdir(vtdir2d)
        end
        if ~exist(vndir2d, 'dir')
            mkdir(vndir2d)
        end
        
        % check if figure output exists already
        fileName = split(fns(i).name, '.tif') ;
        fileName = fileName{1} ;
        outimfn = fullfile(vndir2d, [fileName '.png']) ;
        if save_ims && (~exist(outimfn, 'file') || overwrite_piv)

            % Load the image to put flow on top
            im = imread(fullfile(fns(i).folder, fns(i).name)) ;
            im = cat(3, im, im, im) ;  % convert to rgb for no cmap change
            im = im * washout2d + max(im) * (1-washout2d) ;
            % Control for dilation
            v0t2dsc = v0t2d ./ dilation ;
            options.label = '$|v_t|$ [$\mu$m / min]' ;
            options.xlim = [0, size(im, 2)] ;
            options.ylim = [0, size(im, 1)] ;
            options.outfn = fullfile(vtdir2d, [fileName '.png']) ;
            % plot the 2d tangential flow, adjusted for dilation
            vectorFieldHeatPhaseOnImage(im, x0(1,:), y0(:,1)', ...
                v0t2dsc(:,1), v0t2dsc(:, 2), 10, options) ;
            
            % xlims = xlim ;
            % ylims = ylim ;
            % hold on
            % quiver(x0(:), y0(:), v0t2dsc(:, 1), v0t2dsc(:, 2), 0) ;
            % axis equal
            % xlim(xlims) ;
            % ylim(ylims) ;
            % % Extract image from figure axes
            % patchIm = getframe(gca);
            % outimfn = fullfile(vtdir2d, [fileName '.png']) ;
            % imwrite( patchIm.cdata, outimfn );
            
            %% Now draw normal flow as heatmap
            close all; clear alpha 
            alphaVal = 0.5 ;
            % image = im * washout2d + max(im) * (1-washout2d) ;
            % xlims = xlim ;
            % ylims = ylim ;
            % hold on
            % % normalflow = pcolor(x0, y0, reshape(v0n, size(x0)));
            % normalflow = scatter(x0(:), y0(:), 10, v0n, 'filled') ;
            % axis equal
            % colormap( normalflow, bwr);
            % colorbar( normalflow);
            % set( normalflow, 'AlphaData', alpha );   
            options.caxis = [-10 10] ;
            options.cmap = bwr ;
            options.alpha = 0.5 ;
            v0ngrid = reshape(v0n, size(x0)) ;
            % Create the figure
            figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off')
            
            % Option 1: heatmap on alpha image
            % --------------------------------
            % xfield = x0(1, :)' ;
            % yfield = y0(:, 1) ;
            % image = mat2gray(im, [0, 256]) ;
            % c_handle = heatmap_on_alphaimage(image, yfield, xfield, v0grid, options) ;
            % % Extract image from figure axes
            % patchIm = getframe(gca);
            % % Write figure to file
            % outimfn = fullfile(vndir2d, [fileName '.png']) ;
            % % print('-dpng','-r300', outimfn)
            % imwrite( patchIm.cdata, outimfn );
            
            % Option 2: rgb + overlay
            % --------------------------------
            imshow(im) ; hold on;
            pcolor(x0, y0, v0ngrid)
            shading interp
            colormap(bwr)
            alpha(alphaVal)
            set(gca, 'clim', [-5 5])
            ylim([500 1500])
            c = colorbar();
            % Manually flush the event queue and force MATLAB to render the colorbar
            % necessary on some versions
            drawnow
            % Get the color data of the object that correponds to the colorbar
            cdata = c.Face.Texture.CData;
            % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
            cdata(end,:) = uint8(alphaVal * cdata(end,:));
            % Ensure that the display respects the alpha channel
            c.Face.Texture.ColorType = 'truecoloralpha';
            % Update the color data with the new transparency information
            c.Face.Texture.CData = cdata;
            c.Label.String = 'v_n' ;
            title(['t = ' num2str(t)])
            % Write figure to file
            saveas(gcf, outimfn)
            
        end
        
        if preview
            aux_check_barycentric_triangle
        end
        
        clear x0 y0 uu vv x1 y1
        clear pt0 v0 mesh0x mesh0y 
        clear pt1 v1 mesh1x mesh1y
        clear meshxy meshabove meshbelow 
        clear m0xy mf0 mf0c tr0 fieldfaces baryc0
        clear m1xy mf1 mf1c tr1 t1_contain baryc1
        clear vxa vya vza x123a y123a z123a
        clear vxb vyb vzb x123b y123b z123b
        clear datstruct
        
        % Save this timepoint's piv3d
        disp(['Saving to file: ' sprintf(piv3dfn, time(i))])
        piv3dstruct = piv3d{i} ;
        save(sprintf(piv3dfn, time(i)), 'piv3dstruct')
        
    end 
end
disp('done with piv3dstruct')


%% First do very simpleminded averaging of velocities
do_simpleavg = true ;
plot_vxyz = true ;
if do_simpleavg
    disp('Performing/Loading simple averaging')
    % Create directories
    pivSimAvgDir = fullfile(pivDir, 'simpleAvg') ;
    pivSimAvgImXDir = fullfile(pivSimAvgDir, 'vx') ;
    pivSimAvgImYDir = fullfile(pivSimAvgDir, 'vy') ;
    pivSimAvgImZDir = fullfile(pivSimAvgDir, 'vz') ;
    pivSimAvgImTDir = fullfile(pivSimAvgDir, 'vtH') ;  % Heatmap
    pivSimAvgImGDir = fullfile(pivSimAvgDir, 'vtG ') ;  % Gaussian smoothed in space
    pivSimAvgImSDir = fullfile(pivSimAvgDir, 'vmag ') ;  % speed |v_3D|
    pivSimAvgImQDir = fullfile(pivSimAvgDir, 'vtQ') ;  % Quiverplot
    pivSimAvgImNDir = fullfile(pivSimAvgDir, 'vn') ;
    pivSimAvgImDvgDir = fullfile(pivSimAvgDir, 'dvg') ;
    pivSimAvgImCurlDir = fullfile(pivSimAvgDir, 'curl') ;
    pivSimAvgImShearDir = fullfile(pivSimAvgDir, 'shear_dvphi_ds') ;
    dilDir = fullfile(pivDir, 'dilation') ;
    vxyorigDir = fullfile(pivDir, 'vxyorig') ;
    ensureDir(pivSimAvgDir)
    if plot_vxyz
        ensureDir(pivSimAvgImXDir)
        ensureDir(pivSimAvgImYDir)
        ensureDir(pivSimAvgImZDir)
    end
    dirs2make = {pivSimAvgImTDir, pivSimAvgImGDir,...
        pivSimAvgImNDir, dilDir, vxyorigDir, pivSimAvgImSDir, ...
        pivSimAvgImShearDir, pivSimAvgImDvgDir, pivSimAvgImCurlDir} ;
    for pp = 1:length(dirs2make)
        ensureDir(dirs2make{pp}) ;
    end
    
    % Declare plotting options for limits
    vtscale = 5 ;      % um / min
    vnscale = 2 ;       % um / min
    vscale = 2 ;        % um / min
    alphaVal = 0.7 ;    % alpha for normal velocity heatmap
    qsubsample = 10 ;   % quiver subsampling in pullback space 
    
    % Check if the time smoothed velocities exist already
    v2dsmMumfn = fullfile(pivSimAvgDir, 'v2dMum_simpletimeavg.mat') ;
    v2dsmMfn = fullfile(pivSimAvgDir, 'v2dM_simpletimeavg.mat') ;
    vnsmMfn = fullfile(pivSimAvgDir, 'vnM_simpletimeavg.mat') ;
    vsmMfn = fullfile(pivSimAvgDir, 'vM_simpletimeavg.mat') ;
    % vertex-based velocities
    vvsmMfn = fullfile(pivSimAvgDir, 'vvM_simpletimeavg.mat') ;
    % face-based velocities
    vfsmMfn = fullfile(pivSimAvgDir, 'vfM_simpletimeavg.mat') ;
    if exist(v2dsmMumfn, 'file') && exist(v2dsmMfn, 'file') && ...
            exist(vnsmMfn, 'file') && exist(vsmMfn, 'file') && ...
            exist(vfsmMfn, 'file') && exist(vvsmMfn, 'file') && ...
            ~overwrite_piv
        % load the velocities smoothed in time, (vM v2dM and vnM)
        disp('Loading the time-smoothed velocities...')
        load(v2dsmMfn, 'v2dsmM') ;
        load(v2dsmMumfn, 'v2dsmMum') ;
        load(vnsmMfn, 'vnsmM') ;
        load(vsmMfn, 'vsmM') ;
        load(vvsmMfn, 'vvsmM') ;
        load(vfsmMfn, 'vfsmM') ;
        disp('loaded!')
    else
        disp('Could not load time-smoothed velocities from disk')
        disp('Computing them...')
        first = true ;
        for i = 1:length(piv3d)
            if ~isempty(piv3d{i})
                if first 
                    vM = zeros(length(piv3d), size(piv3d{i}.v0_rs, 1), size(piv3d{i}.v0_rs, 2));
                    vfM = zeros(length(piv3d), size(piv3d{i}.v3dfaces, 1), size(piv3d{i}.v3dfaces, 2)); 
                    vvM = zeros(length(piv3d), ...
                        size(piv3d{i}.v3dvertices, 1), ...
                        size(piv3d{i}.v3dvertices, 2)); 
                    vnM = zeros(length(piv3d), size(piv3d{i}.v0n_rs, 1), size(piv3d{i}.v0n_rs, 2));
                    v2dM = zeros(length(piv3d), size(piv3d{i}.v0t2d, 1), size(piv3d{i}.v0t2d, 2));
                    v2dMum = zeros(length(piv3d), size(piv3d{i}.v0t2d, 1), size(piv3d{i}.v0t2d, 2));
                    first = false ;
                end
                vM(i, :, :) = piv3d{i}.v0_rs ;          % in um/min rs
                % try
                vfM(i, :, :) = piv3d{i}.v3dfaces_rs ;   % in um/min rs
                % catch
                %     v3dfaces = piv3d{i}.v3dfaces ;
                %     vfM(i, :, :) = (rot * v3dfaces')' * resolution / dt(i) ; 
                % end
                vnM(i, :, :) = piv3d{i}.v0n_rs ;        % in rs coords, unit length
                vvM(i, :, :) = (rot * piv3d{i}.v3dvertices')' * resolution / dt(i) ;
                v2dM(i, :, :) = piv3d{i}.v0t2d ;        % in pixels/ min
                v2dMum(i, :, 1) = piv3d{i}.v0t2d(:, 1) ./ piv3d{i}.dilation ; % in pix/min, scaled as um/min
                v2dMum(i, :, 2) = piv3d{i}.v0t2d(:, 2) ./ piv3d{i}.dilation ;
                v2dMum(i, :, 2) = piv3d{i}.v0t2d(:, 2) ./ piv3d{i}.dilation ;
            end
        end
        clearvars first 
        disp('built v0 matrix')
        % Filter in time axis
        % linfilt = 0.1 * ones(10, 1, 1) ;
        % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
        disp('Building tripulse filter equivalent to tripuls()')
        tripulse3 = [ 0.3333; 0.6666; 1; 0.6666; 0.3333];
        tripulse3 = tripulse3 ./ sum(tripulse3(:)) ;
        tripulse3 = reshape(tripulse3, [length(tripulse3), 1]) ;

        vsmM = imfilter(vM, tripulse3, 'replicate');         % in um/min, rs
        vvsmM = imfilter(vvM, tripulse3, 'replicate');       % in um/min, rs
        vnsmM = imfilter(vnM, tripulse3, 'replicate');       % in um/min
        v2dsmM = imfilter(v2dM, tripulse3, 'replicate');     % in pix/min
        v2dsmMum = imfilter(v2dMum, tripulse3, 'replicate'); % in scaled pix/min, proportional to um/min  
        vfsmM = imfilter(vfM, tripulse3, 'replicate') ;      % in um/min, rs

        % Save the simpleminded averaging
        disp('Saving the time-smoothed velocities to disk')
        save(fullfile(pivSimAvgDir, 'v2dMum_simpletimeavg.mat'), 'v2dsmMum') ;  % in scaled pix/min, proportional to um/min 
        save(fullfile(pivSimAvgDir, 'v2dM_simpletimeavg.mat'), 'v2dsmM') ;      % in pix/min
        save(fullfile(pivSimAvgDir, 'vnM_simpletimeavg.mat'), 'vnsmM') ;        % in um/min
        save(fullfile(pivSimAvgDir, 'vM_simpletimeavg.mat'), 'vsmM') ;          % in um/min
        save(fullfile(pivSimAvgDir, 'vvM_simpletimeavg.mat'), 'vvsmM') ;          % in um/min, rs
        save(fullfile(pivSimAvgDir, 'vfM_simpletimeavg.mat'), 'vfsmM') ;        % in um/min, rs
        
        disp('done')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUTOCORRELATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    overwrite_autocorrelations = false ;
    do_acorr = false ;
    redo_acorr = ~exist(fullfile(pivSimAvgDir, 'autocorr_velocities.png'), 'file') ;
    if (redo_acorr || overwrite_autocorrelations) && do_acorr
        aux_autocorrelations
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VELOCITY PLOTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    overwrite_vsm_plots = false ;
    aux_vsm_plots
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot texture map in 3d with velocities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:size(vM, 1)
        t = time(i) ;
        vtqfn0 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_perspective.png']) ;
        vtqfn1 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_dorsal.png']) ;
        vtqfn2 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_ventral.png']) ;
        vtqfn3 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_lateralL.png']) ;
        vtqfn4 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_lateralR.png']) ;
        
        % Plot the tangential velocity as quiver on top of the embedding,
        % colored by normal velocity (all in one in embedding!)
        if ~exist(vtqfn0, 'file') || overwrite_vsm_plots
            disp(['Creating ' vtqfn])
            tic
            
            v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
            vnsm_ii = squeeze(vnsmM(i, :, :)) ;

            % Load the data for the current time point ------------------------
            xp.setT=ime(t) ;
            % Load 3D data for coloring mesh pullback
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();

            fig = figure('units', 'normalized', ...
                    'outerposition', [0 0 1 1], 'visible', 'off') ;
            
            % Psize is the linear dimension of the grid drawn on each triangular face
            Options.PSize = 5;
            Options.EdgeColor = 'none';
            Options.Rotation = rot ;
            Options.FaceScalarField = vn_interp ;

            % Raw stack data
            IV = xp.stack.image.apply();
            IV = imadjustn(IV{1});

            % First args are physical vertices, then texture faces (same number as 
            % physical faces, but connectivity can be different), then texture
            % vertices, which can also be different. The vertices are in row, column, 
            % page format, so x and y get flipped. IV is the texture volume.
            % Options.PSize 
            [ patchIm, imref, zeroID, MIP, SIP ] = ...
                texture_patch_3d( mesh.f, mesh.v, ...
                mesh.f, mesh.v(:, [2 1 3]), IV, Options );
            
            % Add quiver to image
            % Downsample by a factor of qsubsample
            xx = piv.x{i}(1, :)' ;
            yy = piv.y{i}(:, 1) ;
            ww = length(xx) ;
            hh = length(yy) ;
            vx = reshape(v2dsmum_ii(:, 1), [ww, hh]) ;
            vy = reshape(v2dsmum_ii(:, 2), [ww, hh]) ;
            QX = imresize(vx, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
            QY = imresize(vy, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
            xq = 1:qsubsample:ww ;
            yq = 1:qsubsample:hh ;
            [xg, yg] = meshgrid(xq, yq) ;
            h2 = quiver(xg(:), yg(:), QX(:), QY(:), 0) ;
            
            % Format the figure and axis
            xlim(xyzmin(1, :))
            ylim(xyzmin(2, :))
            zlim(xyzmin(3, :))
            title(['$t = $' num2str(t) ' min'], 'Interpreter', 'Latex', 'Color', 'white') 
            xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
            ylabel('DV position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
            zlabel('lateral position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

            % Make background black
            set(gca,'Color','k')
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf, 'Color', 'k')

            % Make tick labels white
            labeltype = {'XTickLabel', 'YTickLabel', 'ZTickLabel'} ;
            for q = 1:2
                % get the current tick labeks
                ticklabels = get(gca, labeltype{q});
                % prepend a color for each tick label
                ticklabels_new = cell(size(ticklabels));
                for i = 1:length(ticklabels)
                    ticklabels_new{i} = ['\color{white} ' ticklabels{i}];
                end
                % set the tick labels
                set(gca, labeltype{q}, ticklabels_new);
            end

            % Capture all five views (perspective + DVLR
            disp(['saving figure...' num2str(t, '%06d')])
            saveas(fig, vtqfn0)
            view(0, 90)
            % DORSAL
            saveas(fig, vtqfn1)
            view(0, 270)
            % VENTRAL
            saveas(fig, vtqfn2)

            % Lateral views
            view(0, 0)
            % make white z tick labels
            q = 3;
            % get the current tick labeks
            ticklabels = get(gca, labeltype{q});
            % prepend a color for each tick label
            ticklabels_new = cell(size(ticklabels));
            for i = 1:length(ticklabels)
                ticklabels_new{i} = ['\color{white} ' ticklabels{i}];
            end
            % set the tick labels
            set(gca, labeltype{q}, ticklabels_new);

            % LEFT
            saveas(fig, vtqfn3)
            view(0, 180)
            % RIGHT
            saveas(fig, vtqfn4)
            close all
            toc
            clear Options IV
        end
        
    end

% Check the orientation of the phasebar
% imshow(im)
% hold on;
clf
xsize = 2000 ;
ysize = 2000 ;
step = 100 ;
imsz = [ xsize ysize ] ;
[xx, yy] = meshgrid(1:step:xsize, 1:step:ysize) ;
ucheck = xx - xsize * 0.5 ;
vcheck = yy - ysize * 0.5 ;
vangle = mod(atan2(vcheck, -ucheck), 2* pi) ;
imshow(im * washout2d + max(im) * (1-washout2d)) ;
hold on ;
h2 = imagesc(xx(:), yy(:), vangle) ;
hold on ;
quiver(xx, yy, ucheck, vcheck, 5) ;
colormap phasemap
phasebar
caxis([0, 2*pi])
axis on
set(gcf, 'visible', 'on')
waitfor(gcf)

end



%% Divergence and Curl (Helmholtz-Hodge)
qsubU = 5 ; 
qsubV = 10 ;
niter_smoothing = [4, 10] ;
plot_dec_pullback = true ;
plot_dec_texturepatch = false ;


% define the output dirs
decDir = fullfile(pivSimAvgDir, 'dec') ;
dvgDir2d = fullfile(pivSimAvgDir, 'dec_div2D') ;
dvgDir3d = fullfile(pivSimAvgDir, 'dec_div3D') ;
dvgDir3dt = fullfile(pivSimAvgDir, 'dec_div3Dt') ;
curlDir2d = fullfile(pivSimAvgDir, 'dec_curl2D') ;
curlDir3d = fullfile(pivSimAvgDir, 'dec_curl3D') ;
curlDir3dt = fullfile(pivSimAvgDir, 'dec_curl3Dt') ;
harmDir2d = fullfile(pivSimAvgDir, 'dec_harm2D') ;
harmDir3d = fullfile(pivSimAvgDir, 'dec_harm3D') ;
% create the output dirs
dirs2do = {decDir, dvgDir2d, dvgDir3d, dvgDir3dt, ...
    curlDir2d, curlDir3d, curlDir3dt, ...
    harmDir2d, harmDir3d} ;
for i = 1:length(dirs2do)
    ensureDir(dirs2do{i}) ;
end

% Consider each timepoint and plot the div and curl
for i = 1:length(piv3d)
    if ~isempty(piv3d{i})
        t = xp.fileMeta.timePoints(i) ;
        % Prepare filenames
        div2dfn = fullfile(dvgDir2d, [sprintf(fileNameBase, t) '_div2d.png']) ; 
        div3dfn = fullfile(dvgDir3d, [sprintf(fileNameBase, t) '_div3d.png']) ;
        div3dtfn = fullfile(dvgDir3dt, [sprintf(fileNameBase, t) '_divt3d.png']) ;
        curl2dfn = fullfile(curlDir2d, [sprintf(fileNameBase, t) '_curl2d.png']) ;
        curl3dfn = fullfile(curlDir3d, [sprintf(fileNameBase, t) '_curl3d.png']) ;
        curl3dtfn = fullfile(curlDir3dt, [sprintf(fileNameBase, t) '_curlt3d.png']) ;
        harm2dfn = fullfile(harmDir2d, [sprintf(fileNameBase, t) '_harm2d.png']) ; 
        harm3dfn = fullfile(harmDir3d, [sprintf(fileNameBase, t) '_harm3d.png']) ;
        
        % Obtain smoothed velocities on all faces
        vfsm = squeeze(vfsmM(i, :, :)) ;
        v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
        
        % Use current time's tiled smoothed mesh
        % Note: vfsmM is in um/min rs
        FF = piv3d{i}.m0f ;   % #facesx3 float: mesh connectivity list
        V2D = piv3d{i}.m0XY ; % Px2 float: 2d mesh vertices in pullback image pixel space
        v3drs = ((rot * piv3d{i}.m0v3d')' + trans) * resolution ;
        cutM.f = FF ;
        cutM.u = V2D ;
        cutM.v = v3drs ;
        cutM.nU = nU ;
        cutM.nV = nV ;
        disp('Decomposing flow into div/curl...')
        [divs, rots, harms, glueMesh] = ...
            helmHodgeDECRectGridPullback(cutM, vfsm, ...
            'niterSmoothing', niter_smoothing, ...
            'clipDiv', [-5, 5], 'clipRot', [-0.5, 0.5], ...
            'preview', preview, 'method', 'both') ;
        
        % save divs, rots, and harms
        save(fullfile(decDir, [sprintf(fileNameBase, t) '_dec.mat']), ...
            'divs', 'rots', 'harms')
        
        % Plot results
        disp('Plotting div/curl...')
        if plot_dec_pullback
            im = imread(fullfile(fns(i).folder, fns(i).name)) ;
            im = cat(3, im, im, im) ;  % convert to rgb for no cmap change
            addTitleStr = ['$t=$', num2str(t - min(fold_onset))] ;
            Options.addTitleStr = addTitleStr ;
            Options.div2dfn = div2dfn ;
            Options.div3dfn = div3dfn ;
            Options.rot2dfn = curl2dfn ;
            Options.rot3dfn = curl3dfn ;
            Options.qsubU = 5 ; 
            Options.qsubV = 5 ;
            Options.sscaleDiv = 0.5 ;
            Options.sscaleRot = 0.2 ;
            Options.qscaleDiv = 50 ;
            Options.qscaleRot = 50 ;
            Options.xyzlim = xyzlim ;
            opts2d.xlim = [0, size(im, 1)] ;
            opts2d.ylim = [0.25 * size(im, 2), 0.75 * size(im, 2) ] ;
            xy = {piv3d{i}.x0, piv3d{i}.y0} ;
            plotHelmHodgeDECPullback(im, cutM, vfsm, xy, v2dsmum_ii, ...
                divs, rots, Options, opts2d)
        end
        if plot_dec_texturepatch
            % load current timepoint
            % (3D data for coloring mesh pullback)
            xp.loadTime(t);
            xp.rescaleStackToUnitAspect();

            % Raw stack data
            IV = xp.stack.image.apply();
            IV = imadjustn(IV{1});     
            IV = max(IV(:)) - IV ;

            addTitleStr = ['$t=$', num2str(t)] ;
            Options.addTitleStr = addTitleStr ;
            Options.div3dfn = div3dtfn ;
            Options.rot3dfn = curl3dtfn ;
            Options.sscaleDiv = 0.5 ;
            Options.sscaleRot = 0.2 ;
            Options.xyzlim = xyzlim ;
            
            % Load the cutmesh vertices and normals
            cutM.v = piv3d{i}.m0v3d ;
            cutM.v3drs = ((rot * cutM.v')' + trans) * resolution ;
            plotHelmHodgeDECTexture3d(IV, cutM, divs, rots, rot, trans, Options)
        end
    end
end


%% Simpleminded streamlines from velocity scaled by dilation
% Build 3d grid of positions and velocities
% Assume that x0, y0 are fixed for all time
ntimes = length(piv.x) ;
x0 = piv.x{i} ;
y0 = piv.y{i} ;
xlen = size(v2dsmM, 1) ;
ylen = size(v2dsmM, 2) ;
vdat = v2dsmM(:) ;
vxdat = reshape(vdat(:, 1), [ntimes, xlen, ylen]) ;
vydat = reshape(vdat(:, 2), [ntimes, xlen, ylen]) ;
vzdat = ones(size(vydat)) ; % we march through time at 1 index / timestep
% Define positions we track through the streamlines
startx = x0(1:200:end) ;
starty = y0(1:200:end) ;
startz = zeros(size(starty)) ;
streamline(x0, y0, z0, vxdat, vydat, vzdat, startx, starty, startz)
view(2)
save(gcf, fullfile(pivDir, 'streamlines.png')) 
error('break')


%% Smooth velocities in time ==============================================
% Make correspondences between faces in t_i, t_{i-1}, and t_{i+1} 
% Interpolate velocities on grid to smooth them in time
piv3dfn = fullfile(pivDir, 'piv3d.mat') ;
if exist(piv3dfn, 'file')
    load(piv3dfn)
else
    piv3d = cell(length(fns), 1) ;
    for i=1:length(fns) - 1
        % Average in time by advecting v3d along correlation vector 
        x0 = piv.x{i} ;
        y0 = piv.y{i} ;
        uu = piv.u_filtered{i} ;
        v0 = piv.v_filtered{i} ; 
        v3d0 = piv3d{i}.v0 ;
        v3d0grid = reshape(v3d0, [size(piv.x{1}, 1), size(piv.x{1}, 2), 3]) ;
        x0vel = squeeze(v3d0grid(:, :, 1)) ;
        y0vel = squeeze(v3d0grid(:, :, 2)) ;
        z0vel = squeeze(v3d0grid(:, :, 3)) ;
                
        % Advect 0 -> 1
        x0ad1 = x0 + u0 ;
        y0ad1 = y0 + v0 ;
        
        % Advect 0 -> -1
        x0ad1 = x0 - u0 ;
        y0ad1 = y0 - v0 ;
        
        % next timept
        tpid = i + 1 ;
        if tpid < length(fns)
            x1 = piv.x{tpid} ;
            y1 = piv.y{tpid} ;
            u1 = piv.u_filtered{tpid} ;
            v1 = piv.v_filtered{tpid} ; 
            v3d1 = piv3d{tpid}.v0 ;
            % cast as NxMx3, then as N*M x 1 arrays for vx,vy,vz separately
            v3d1grid = reshape(v3d1, [size(piv.x{1}, 1), size(piv.x{1}, 2), 3]) ;
            x1vel = squeeze(v3d1grid(:, :, 1)) ;
            y1vel = squeeze(v3d1grid(:, :, 2)) ;
            z1vel = squeeze(v3d1grid(:, :, 3)) ;
            % Find the velocity at the advected position using only 2d coords
            Fx = griddedInterpolant(x1', y1', x1vel') ;
            Fy = griddedInterpolant(x1', y1', y1vel') ;
            Fz = griddedInterpolant(x1', y1', z1vel') ;
            v3d1 = [Fx(x0ad1(:), y0ad1(:)), ...
                Fy(x0ad1(:), y0ad1(:)), Fz(x0ad1(:), y0ad1(:))] ;
        else
            v3d1 = v3d0 ;
        end
        
        % previous timept
        tpid = i - 1; 
        if tpid > 0
            x2 = piv.x{tpid} ;
            y2 = piv.y{tpid} ;
            u2 = piv.u_filtered{tpid} ;
            v2 = piv.v_filtered{tpid} ; 
            v3d2 = piv3d{tpid}.v0 ;
            % cast as NxMx3, then as N*M x 1 arrays for vx,vy,vz separately
            v3d2grid = reshape(v3d2, [size(piv.x{1}, 1), size(piv.x{1}, 2), 3]) ;
            x2vel = squeeze(v3d2grid(:, :, 1)) ;
            y2vel = squeeze(v3d2grid(:, :, 2)) ;
            z2vel = squeeze(v3d2grid(:, :, 3)) ;
            % Find the velocity at the advected position using only 2d coords
            Fx = griddedInterpolant(x0', y0', x2vel') ;
            Fy = griddedInterpolant(x0', y0', y2vel') ;
            Fz = griddedInterpolant(x0', y0', z2vel') ;
            v3d2 = [Fx(x0ad2(:), y0ad2(:)), ...
                Fy(x0ad2(:), y0ad2(:)), Fz(x0ad2(:), y0ad2(:))] ;
        else
            v3d2 = v3d0 ;
        end
        
        % Perform a boxcar average & store the result
        v0avg = (1/3) * (v3d0 + v3d1 + v3d2) ;
        piv3d{i}.v0avg = v0avg ;
        
        % Decompose the result into tangential and normal components
        
        % Plot the tangential velocity as heatmap on top of the image
        fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'off') ;
        imshow(im * washout2d + max(im) * (1-washout2d)) ;
        hold on;
        h2 = imagesc(piv.x{i}(1, :), piv.y{i}(:, 1), vangle) ;
        colormap phasemap
        % phasebar
        set(h2, 'AlphaData', speed)
        plot([foldx; foldx], [0, 0, 0; yesz, yesz, yesz], 'k--')
        saveas(fig, fullfile(pivSimAvgImTDir, [sprintf('%04d', time(i)) '.png'])) ;    
        close all
    end
end



%% Check phi0s
close all
for ii = 1:length(xp.fileMeta.timePoints)
    t = time(ii) ;
    % Load the spcutMesh for this timepoint
    disp(['Loading spcutMesh from disk... [t = ' num2str(t) ']'])
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;

    % plot cut paths
    vv = spcutMesh.v ;
    inds = spcutMesh.pathPairs(:, 1) ;
    clf;
    trisurf(spcutMesh.f, spcutMesh.v(:, 1), spcutMesh.v(:, 2), ...
        spcutMesh.v(:, 3), spcutMesh.v(:, 1), 'edgecolor', 'none', 'Facealpha', 0.2); 
    hold on
    plot3(vv(inds, 1), vv(inds, 2), vv(inds, 3), '.-')
    pause(0.1)
    
    
    % % plot phi0s
    % plot(spcutMesh.phi0s, '.')
    % hold on;
    % plot(spcutMesh.phi0_fit, '-')
    % pause(0.1)
    % clf
end


