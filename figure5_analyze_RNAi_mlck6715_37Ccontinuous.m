% analyze_RNAi_results
% Here, 
%   13 is bands forming, late 13 is bands formed
%   14a is closing
%   14b is mostly closed but not fully
%   15a is after closure
%   15b is after/during first fold

% 0: tub67;tub15
% 0.1: 48YGAL4 (no UAS)
% 0.2: 48YGAL4 klar x UAS-CAAX mCh
% 0.3: mef2Gal4 (no UAS)
% 0.4: UAS-SERCA (no GAL4)
% 1: tub67;tub15 x UAS-MLCK RNAi TRiP#4
% 2: mef2G4klar x UAS-SERCA.R751Q

% 0 --> WT
% 1 --> disrupted but folds
% 2 --> incomplete Fold
% 3 --> missing one Fold
% 4 --> missing two folds
% 5 --> missing all folds
% 6 --> no gut closure
% 7 --> dead
close all
clear
outdir = '/mnt/data/analysis/mlckRNAi';
x = 100 ;
scoreLowerLimit = 1 ;
scoreUpperLimit = 6 ;
genoM = 1 ;
genoC = 0 ;

%% continuous 37C heatshock
% Date, EmbryoID, Genotype, Stage, Result
date = 1 ;
embryoID = 2 ; 
genotype = 3 ;
stage = 4;
score = 5; 
minStage = 13.1 ;
maxStage = 15.2 ;


stage_result_mlck = [...
... % 202107201643_6715MLCKRNAi4_37C_4mpf_15um_7pc488BioII_10x2p66x_bioII CONTINUOUS
202107201643, 1, 1, 14.1, 2; ...
202107201643, 2, 1, 13, 4; ...
202107201643, 3, 1, 15.1, 0; ...
202107201643, 4, 1,  11, 6; ...
202107201643, 5, 1,  13, 3; ...
202107201643, 6, 1, 12, 5; ...
202107201643, 7, 1,  11,  6; ...
202107201643, 8, 1, 16.1, 0; ...
202107201643, 9, 1, 16.1, 0; ...
202107201643, 10, 1, 15.1 , 4; ...
202107201643, 11, 1, 14.2, 5; ...
202107201643, 12, 1, 16.1, 0; ...
202107201643, 13, 1, 13, 3 ; ...  %<-- no posterior fold, anterior is incomplete
202107201643, 14, 1, 13, 6; ...
202107201643, 15, 1, 13, 5 ; ...
202107201643, 16, 1, 13, 7; ...
202107201643, 17, 1, 15.1, 4 ; ...
202107201643, 18, 1,  11, 6; ...
202107201643, 19, 1,  14.1, 6; ...
202107201643, 20, 1, 11, 6; ...
202107201643, 21, 1,  13, 5; ...
202107201643, 22, 1,  15.1, 5; ...
202107201643, 23, 1,  10, 6; ...
202107201643, 24, 1,  13, 4; ...
202107201643, 25, 1,  15.2, 1; ...
202107201643, 26, 1,  16.1, 7; ...
202107201643, 27, 1,  x, x ; ...
202107201643, 28, 1, 14.2, 4; ...   % e28
202107201643, 29, 1,  15.1, 0; ...
202107201643, 30, 1, 14.2, 1; ...
202107201643, 31, 1, 13, 6; ...
202107201643, 32, 1, 10, x; ...
202107201643, 33, 1, 11, x; ...
202107201643, 34, 1, 16.1, 0; ...
202107201643, 35, 1, x, x; ...
202107201643, 36, 1, 13, x; ...
202107201643, 37, 1, 13, x; ...
202107201643, 38, 1, 15.1, 5; ...
202107201643, 39, 1, 15.2, 4; ...  %<-- no further folding once HS begins
202107201643, 40, 1, 5.2, 3; ...  % <-- no anterior fold
202107201643, 41, 1, x, 7; ...
202107201643, 42, 1, 15.2, 0; ...
202107201643, 43, 1, 11, 6; ...
202107201643, 44, 1, 15.1, 1; ...
202107201643, 45, 1, 14.1, 4; ...
202107201643, 46, 1, 15.1, 1; ...
202107201643, 47, 1, 13, 5; ...
202107201643, 48, 1, 16.1, 0; ...
 ];

% load controls

    
stage_result_control_continuous = [...
... % 20210803 bioII mef2 table: date, embryoID, initial HS stage, all normal development (0)
    202108031752, 1, 0.3, 13, 0 ; ... 
    202108031752, 2, 0.3, 14, 0 ; ... 
    202108031752, 3, 0.3, 16.1, 0 ; ... 
    202108031752, 4, 0.3, 12, 0 ; ... 
    202108031752, 5, 0.3, 16.1, 0 ; ... 
    202108031752, 6, 0.3, 16.1, 0 ; ... 
    202108031752, 7, 0.3, 13, 0 ; ... 
    202108031752, 8, 0.3, 16.1, 0 ; ... 
    202108031752, 9, 0.3, 15.1, 0 ; ... 
... % 20210803 bioII tub67;15 table:
    202108031752, 10, 0, 15.1, 0 ; ... 
    202108031752, 11, 0, 15.1, 0 ; ... 
    202108031752, 12, 0, 15.1, 0 ; ... 
    202108031752, 13, 0, 15.1, 0 ; ... 
    202108031752, 14, 0, 11, x ;... <---- not done yet at time of viewing (21:32)
    202108031752, 15, 0, 15.1, 0 ; ... 
    202108031752, 16, 0, 15.2, 0 ; ... 
    202108031752, 17, 0, 13, 0 ; ... 
    202108031752, 18, 0, 15.2, 0 ; ... 
    202108031752, 19, 0, 15.1, 0 ; ... 
    202108031752, 20, 0, 13, 0 ; ... 
    202108031752, 21, 0, 15.1, 0 ; ... 
...% 202108101608_mef2G4k_6715cntrl_37Ccontinuous_boiII_4mpf_20um_7pc488
202108101608, 23, 0.3, 15.1, 1 ;
202108101608, 24, 0.3, 13, 4 ;
202108101608, 25, 0.3, 14.2, 2;
202108101608, 26, 0.3, 13, 1;
202108101608, 27, 0.3, 15.1, 1;
202108101608, 28, 0.3, 15.1 , 1;
202108101608, 29, 0.3, 13.2, 1;
202108101608, 30, 0.3, 15.1, 0;
202108101608, 31, 0.3, 13, 0;
...% tub67;15 table:
...% 202108101608_mef2G4k_6715cntrl_37Ccontinuous_boiII_4mpf_20um_7pc488
202108101608, 1, 0, 16.1, 1;
202108101608, 2, 0, 15.1, 2 ; % incomplete posterior
202108101608, 3, 0, 15.1, 2;
202108101608, 4, 0, 16.1, 1;
202108101608, 5, 0, 14.2, 3;  % difficult to see, dark, but definitely not normal. At the very least incomplete posterior
202108101608, 6, 0, 15.1, 1 ;
202108101608, 7, 0, 13, 6;
202108101608, 8, 0, 14.1, 1;
202108101608, 9, 0, 13.2,  5;
202108101608, 10, 0, 15.1, 2;
202108101608, 11, 0, 15.2, 1;
202108101608, 12, 0, 15.1, 0; 
202108101608, 13, 0, 14.1, 0; 
202108101608, 14, 0, 15.1, 0; 
202108101608, 15, 0, 15.1, 0;
202108101608, 16, 0, 15.1, 0
202108101608, 17, 0, 13, 0;
202108101608, 18, 0, 14.1, 0;
202108101608, 19, 0, 15.1, 0;
202108101608, 20, 0, 15.1, 0 ;
202108101608, 21, 0, 14.2, 0;
202108101608, 22, 0, 15.1, 0;
];

% %%  202107241336 CONTROL 55min --  date, genotype, embryoID, starting stage, result
% stage_result_control_55min = [...
% 202107241336, 1, 0,  13, 0; ...
% 202107241336, 2, 0,  15.1,0 ; ...
% 202107241336, 3, 0,  13, 4; ...
% 202107241336, 4, 0,  14.2, 0; ...
% 202107241336, 5, 0,  15.1, 0; ...
% 202107241336, 6, 0,  14.1, 2; ...
% 202107241336, 7, 0,  10, 5 ; ...
% 202107241336, 8, 0,  15.1, 0; ...
% 202107241336, 9, 0,  11, 5; ...
% 202107241336, 10, 0,  14.2, 0 ; ...
% 202107241336, 11, 0,  16.1, 0; ...
% 202107241336, 12, 0,  13, 3; ...
% 202107241336, 13, 0,  x, x; ...
% 202107241336, 14, 0,  15.1, 0; ...
% 202107241336, 15, 0,  14.2, 3; ...
% 202107241336, 16, 0,  11, 0; ...
% 202107241336, 17, 0,  11, 0; ...
% 202107241336, 18, 0,  14.2, 3; ...
% 202107241336, 19, 0,  12 , 5; ... % <--- initially upside down within vitelline
% 202107241336, 20, 0,  15.1, 0; ...
% 202107241336, 21, 0,  15.2, 0; ...
% 202107241336, 22, 0,  15.2, 2; ... % <-- incomplete anterior forld
% 202107241336, 23, 0,  15.1, 0; ...
% 202107241336, 24, 0,  12,   3; ... % <-- irregular folds, as if gut is sideways but embryo is not
% 202107241336, 25, 0,  15.1, 0; ...
% 202107241336, 26, 0,  13,   5; ...
% 202107241336, 27, 0,  15.1, 0; ...
% 202107241336, 28, 0,  x ,   x; ...
% 202107241336, 29, 0,  14.1, 3; ...
% 202107241336, 30, 0,  16.1, 3; ...
% 202107241336, 31, 0,  15.1, 0; ...
% 202107241336, 32, 0,  15.1, 0; ...
% 202107241336, 33, 0,  15.1, 0 ; ...
% 202107241336, 34, 0,  12,   2;  ...% <-- folds out of order, anterior is very late as to be incomplete
% 202107241336, 35, 0,  15.1, 0 ; ...
% 202107241336, 36, 0,  15.1, 0; ...
% 202107241336, 37, 0,  15.1, 0 ; ...
% 202107241336, 38, 0,  15.1,1 ; ...
% 202107241336, 39, 0,  12, 0; ...
% 202107241336, 40, 0,  13, 0; ...
% 202107241336, 41, 0,  15.1, 0 ; ...
% 202107241336, 42, 0,  15.1, 0; ...
% 202107241336, 43, 0,  11, 3 ;...
% 202107241336, 44, 0,  10, 0 ;...
% 202107241336, 45, 0,  15.1, 0 ;...
% 202107241336, 46, 0,  13, 7 ;...
% 202107241336, 47, 0,  15.2, 0;...
% ];


results = [stage_result_mlck; stage_result_control_continuous];
    
for ii = 1:2
    
    inBin = find(results(:, stage) > minStage & ...
        results(:, stage) < maxStage) ;
    mutantIdx = intersect(inBin, find(results(:, genotype) == genoM)) ;
    
    if ii == 1
        exten = sprintf('_restricted_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) == genoC)) ;
    else
        exten = sprintf('_allControls_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) < 1)) ;
    end
    exten = strrep(exten, '.', 'p') ;
    

    clf; 
    close all;

    % Clean out those that are x or 7
    keepM = find(results(mutantIdx, score) < scoreUpperLimit) ;
    keepC = find(results(controlIdx, score) < scoreUpperLimit) ;
    mutantIdx = mutantIdx(keepM) ;
    controlIdx = controlIdx(keepC) ;

    % N mutant and N control
    nM = length(mutantIdx) ;
    nC = length(controlIdx) ;

    misMutant = intersect(mutantIdx, find(results(:, score) > scoreLowerLimit)) ;
    foldMutant = intersect(mutantIdx, find(results(:, score) <= scoreLowerLimit)) ;
    misControl = intersect(controlIdx, find(results(:, score) > scoreLowerLimit)) ;
    foldControl = intersect(controlIdx, find(results(:, score) <= scoreLowerLimit)) ;
    fracBadFolds_m = length(misMutant) / nM ;
    fracBadFolds_c = length(misControl) / nC ;

    assert(nM == length(misMutant) + length(foldMutant))
    assert(nC == length(misControl) + length(foldControl))

    nfoldM = length(foldMutant) ;
    nmisM = length(misMutant) ;
    nfoldC = length(foldControl) ;
    nmisC = length(misControl) ;

    % Significance
    xx = table([nfoldM; nmisM],[nfoldC;nmisC],...
        'VariableNames',{'RNAi','control'},'RowNames',{'folded','misfolded'}) ;

    [h,pval,stats] = fishertest(xx) ;

    success = [nfoldM / nM, nfoldC/nC] ;
    bar(success)
    hold on;
    yerr0 = sqrt((success .* (1-success)) ./ [nM, nC]) ;
    yneg = min(yerr0, success) ;
    ypos = min(yerr0, 1-success) ;
    errorbar([1,2], success, yneg, ypos , 'LineStyle','none')
    set(gca,'xticklabel',{'UAS-MLCK RNAi', 'control'});
    ylabel('probablity of forming three folds')

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C' exten '.png'])) ;

    %% superbar plot
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    clf;
    Y = success;
    E = cat(3, yneg, ypos) ;

    Colors = [
        0.90    0.55    0.55
        0.62    0.76    0.84
        0.89    0.10    0.11
        0.12    0.47    0.70
        ];
    Colors = reshape(Colors, [2 2 3]);

    P = [pval, pval; pval, pval];
    % Make P symmetric, by copying the upper triangle onto the lower triangle
    % PT = P';
    % lidx = tril(true(size(P)), -1);
    % P(lidx) = PT(lidx);

    superbar(Y, 'E', E, 'P', P, 'BarFaceColor', Colors, 'Orientation', 'v', ...
        'ErrorbarStyle', 'I', 'PLineOffset', 0.1, 'PStarShowGT', false);

    xlim([0.5 2.5]);
    % ylim([0 1]);
    set(gca, 'YTick', [0, 0.5, 1])
    set(gca, 'XTick', [1, 2]);
    ylims = ylim;
    expVal = sprintf('%e', pval) ;
    keepDigits = expVal(end-2:end) ;
    text(1.5, ylims(2)-0.03, ...
        ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pval)], 'interpreter', 'latex', 'horizontalalignment', 'center')

    set(gca,'xticklabel',{'UAS-MLCK RNAi', 'control'});
    ylabel('probablity of forming three folds')
    title('tub>67;tub>15 x UAS-MLCK RNAi')

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_superbar' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_superbar' exten '.png'])) ;
    

    %% Print contingency table
    xx

    fn = fullfile(outdir, ['rnai_continuous37C_contingency' exten '.txt']) ;
    writetable( xx, fn, 'WriteRowNames', true)


end



%% Quantized stats
for ii = 1:2
    
    close all
    inBin = find(results(:, stage) > minStage & ...
        results(:, stage) < maxStage) ;
    mutantIdx = intersect(inBin, find(results(:, genotype) == genoM)) ;
    
    if ii == 1
        exten = sprintf('_restricted_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) == genoC)) ;
    else
        exten = sprintf('_allControls_min%0.1f_max%0.1f', minStage, maxStage) ;
        controlIdx = intersect(inBin, find(results(:, genotype) < 1)) ;
    end
    exten = strrep(exten, '.', 'p') ;
    
    % Clean out those that are x or 7
    keepM = find(results(mutantIdx, score) < scoreUpperLimit) ;
    keepC = find(results(controlIdx, score) < scoreUpperLimit) ;
    mutantIdx = mutantIdx(keepM) ;
    controlIdx = controlIdx(keepC) ;

    % N mutant and N control
    nM = length(mutantIdx) ;
    nC = length(controlIdx) ;
    nT = nM + nC ;

    mutant0folds = intersect(mutantIdx, find(results(:, score) == 5)) ;
    mutant1folds = intersect(mutantIdx, find(results(:, score) == 4)) ;
    mutant2folds = intersect(mutantIdx, find(results(:, score) == 3 | results(:, score) == 2)) ;
    mutant3folds = intersect(mutantIdx, find(results(:, score) <= 1)) ;
    control0folds = intersect(controlIdx, find(results(:, score) == 5)) ;
    control1folds = intersect(controlIdx, find(results(:, score) == 4)) ;
    control2folds = intersect(controlIdx, find(results(:, score) == 3 | results(:, score) == 2)) ;
    control3folds = intersect(controlIdx, find(results(:, score) <= 1)) ;

    nM0 = numel(mutant0folds) ;
    nM1 = numel(mutant1folds) ;
    nM2 = numel(mutant2folds) ;
    nM3 = numel(mutant3folds) ;
    nC0 = numel(control0folds) ;
    nC1 = numel(control1folds) ;
    nC2 = numel(control2folds) ;
    nC3 = numel(control3folds) ;
    hisMutant = [nM0, nM1, nM2, nM3] ;
    hisControl = [nC0, nC1, nC2, nC3] ;
    
    nfolds = [0, 1, 2, 3] ;
    colors = define_colors ;
    colors = [
        0.90    0.55    0.55 ;
        0.62    0.76    0.84 ] ;
    bb = barh(nfolds, [ hisMutant', hisControl'], 'histc') ;
    bb(1).FaceColor = colors(1, :);
    bb(2).FaceColor = colors(2, :);
    xlabel('count', 'interpreter', 'latex')
    ylabel('number of folds', 'interpreter', 'latex')
    legend({'MLCK RNAi', 'control'})
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_barh' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_barh' exten '.png'])) ;
    
    boxMutant = [0*ones(nM0, 1); ...
                 1*ones(nM1, 1); ...
                 2*ones(nM2, 1); ...
                 3*ones(nM3, 1)] ;
    boxControl = [0*ones(nC0, 1); ...
                 1*ones(nC1, 1); ...
                 2*ones(nC2, 1); ...
                 3*ones(nC3, 1)] ;
    g1 = repmat({'MLCK RNAi'},length(boxMutant),1);
    g2 = repmat({'control'},length(boxControl),1);
    gg = [g1; g2];
    boxplot([boxMutant; boxControl], gg)
    ylabel('number of folds')
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_boxp' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_boxp' exten '.png'])) ;
   
    % attempt 1--a bit awkward, uses barh and subtightplot
    %     clf
    %     marg_h = 0.1 ;
    %     marg_w = 0.1 ;
    %     ax1 = subtightplot(1,2, 1,[0,0], marg_h, marg_w)
    %     vMutant = hisMutant ./ sum(hisMutant) ;
    %     bm = barh(nfolds,0.5*[vMutant;-vMutant;-vMutant]', 'stacked')
    %     xlim([-0.5, 0.5])
    %     bm(1).FaceColor = colors(1, :);
    %     bm(2).FaceColor = colors(1, :);
    %     bm(3).FaceColor = colors(1, :);
    %     bm(1).EdgeColor = 'none';
    %     bm(2).EdgeColor = 'none';
    %     bm(3).EdgeColor = 'none';
    %     ylabel('number of folds')
    %     xlabel('frequency')
    %     xticks([])
    %     
    %     ax2 = subtightplot(1,2, 2,[0,0], marg_h, marg_w)
    %     vControl = hisControl ./ sum(hisControl) ;
    %     bc = barh(nfolds,0.5*[vControl;-vControl;-vControl]', 'stacked') ;
    %     xlim([-0.5, 0.5])
    %     bc(1).FaceColor = colors(2, :);
    %     bc(2).FaceColor = colors(2, :);
    %     bc(3).FaceColor = colors(2, :);
    %     bc(1).EdgeColor = 'none';
    %     bc(2).EdgeColor = 'none';
    %     bc(3).EdgeColor = 'none';
    %     xticks([])
    %     
    %     linkaxes([ax1, ax2], 'y');
    
    
    %%
    close all
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    
    nMutant = sum(hisMutant) ;
    nControl = sum(hisControl) ;
    vMutant = hisMutant ./ sum(hisMutant) ;
    vControl = hisControl ./ sum(hisControl) ;
    gg = [0,0,0,0, 1, 1, 1, 1];
    % boxyViolinPlot(gg, [nfolds, nfolds], ...
    %     [vMutant, vControl], 1, 'center') ;
    boxyViolinPlot([0,0,0,0], nfolds, ...
        vMutant, 1, 'center', colors(1, :), 'none') ;
    boxyViolinPlot([1,1,1,1], nfolds, ...
        vControl, 1, 'center', colors(2, :), 'none') ;
    
    % Significance
    num = mean(boxMutant) -  mean(boxControl)  ;
    denom = sqrt(nanvar(boxControl) / nControl + ...
        nanvar(boxMutant) / nMutant) ;
    zscoreNumFolds = num / denom ;
    pvalueNumFolds = normcdf(zscoreNumFolds) ;
    
    success = [ mean(boxMutant) ,  mean(boxControl) ] ;
    yerr0 = [std(boxMutant) / sqrt(nMutant), std(boxControl) / sqrt(nControl) ];
    yneg = min(yerr0, success) ;
    ypos = min(yerr0, 3-success) ;
    
    
    % finish plot with labels
    ylim([-0.5, 4])
    set(gca, 'YTick', [0, 1, 2, 3])
    set(gca, 'XTick', [0, 1]);
    expVal = sprintf('%e', pvalueNumFolds) ;
    keepDigits = expVal(end-1:end) ;
    text(0.5, 3.8, ...
        ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pvalueNumFolds)], 'interpreter', 'latex', 'horizontalalignment', 'center')

    set(gca,'xticklabel',{'MLCK RNAi', 'control'});
    ylabel('number of folds')
    
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_violin' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_violin' exten '.png'])) ;
    
    
    %% 
    hold on;
    errorbar([0, 1], [mean(boxMutant), mean(boxControl)], ...
        yneg, ypos, 'LineWidth', 2, 'LineStyle','none', 'color', [0.5, 0.5,0.5])
   
    xlim([-0.5 1.5]);
    ylim([-0.5, 4])
    set(gca, 'YTick', [0, 1, 2, 3])
    set(gca, 'XTick', [0, 1]);
    
    set(gca,'xticklabel',{'MLCK RNAi', 'control'});
    ylabel('number of folds')

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_violinErr' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_violinErr' exten '.png'])) ;

    %% Flip the histogram sideways
    
    close all
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    pos = get(gca, 'Position') ;
    set(gca, 'Position', pos + [0, 0.15, 0, -0.1])
    boxyViolinPlot([0,.2,.4,.6], 0, ...
        0.2, vMutant, 'left', colors(1, :), 'none') ;
    boxyViolinPlot(1 + [0,.2,.4,.6], 0, ...
        0.2, vControl, 'left', colors(2, :), 'none') ;
    
    text(0.9, 0.7, ...
        ['$p=$' sprintf(['%0.' num2str(keepDigits) 'f'],pvalueNumFolds)], 'interpreter', 'latex', 'horizontalalignment', 'center')

    xticks([0.1, 0.3, 0.5, 0.7, 1.1, 1.3, 1.5, 1.7])
    xticklabels([0,1,2,3, 0,1,2,3])
    yticks([0, 0.5, 1])
    xlabel('number of folds')
    ylabel('frequency')
    xlim([-0.2, 2])
    ylim([0, 1])
    text(0.4, 0.6, ['N=', num2str(nMutant)], 'HorizontalAlignment', 'center')
    text(1.4, 0.6, ['N=', num2str(nControl)], 'HorizontalAlignment', 'center')
    text(0.4, -0.2, 'MLCK RNAi', 'HorizontalAlignment', 'center')
    text(1.4, -0.2, 'control', 'HorizontalAlignment', 'center')
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_clusterHist' exten '.pdf'])) ;
    saveas(gcf, fullfile(outdir, ['rnai_continuous37C_clusterHist' exten '.png'])) ;

end
