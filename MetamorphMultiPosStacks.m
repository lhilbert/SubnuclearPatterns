clear all

%% Choose input directory and compile file list

plotSegmentation = false;

sourceDirectory = ...
    '/Users/hilbert/academia/VastenhouwLab/Flavopiridol_Sphere/FP_Inhibition_Sphere_Nov2016/';

filepath_cell = {...
    'Flavopiridol_EU_RNA_17Nov2016/CTRL_0001.nd',...
    'Flavopiridol_EU_RNA_17Nov2016/CTRL_0002.nd',...
    'Flavopiridol_EU_RNA_17Nov2016/FP_0002.nd',...
    'Flavopiridol_EU_RNA_17Nov2016/FP_0004.nd'};

seriesExcludeCells = {[],[],[],[],[],[],[],[],[],[],[]};

cond_vec = [0,0,1,1]; % [Flavopiridol] in uM

useSets = [1:4];
% useSets = [4];

filepath_cell = filepath_cell(useSets);
seriesExcludeCells = seriesExcludeCells(useSets);
cond_vec = cond_vec(useSets);

numSets = numel(filepath_cell);

fullPaths = cell(1,numSets);

for kk = 1:numSets
    
    fullPaths{kk} = fullfile(sourceDirectory,filepath_cell{kk});
    
end

%% --- analysis parameters

segChannel = 3;
histStepSize = 1;

% % unit: square micrometers
minArea = 40; % minimum area for nucleus to be recognized
maxArea = 210;% maximum area for nucleus to be recognized

%% --- read and segment 3D stacks

% error tracking
errorflag = false(1,numSets);
errorMessages = cell(1,numSets);

% Store different images of nucleus
maxContrastCell = cell(1,numSets);
segCell = cell(1,numSets);
xNucSectionCell = cell(1,numSets);
shellImageCell = cell(1,numSets);

% Quantification outcomes
nucVol = zeros(1,numSets);
nucCent = cell(1,numSets);

nucArea_cell = cell(1,numSets);
nucInt_cell = cell(1,numSets);
cytoInt = cell(1,numSets);

voxelSizes = zeros(1,numSets);

spectra = cell(1,numSets);

% hist_ints = cell(1,numSets);
% hist_counts = cell(1,numSets);
% int_var = zeros(1,numSets);
% fft_cell = cell(1,numSets);
neighbor_matching_cell = cell(1,numSets);
AClength_cell = cell(1,numSets);
nucImage_cell = cell(1,numSets);

for ff = 1:numSets
    
    fprintf('Processing file %d of %d\n', ...
        ff,numSets)
    
    filepath = fullPaths{ff};
    
    % --- Make a reader instance
    reader = bfGetReader(filepath);
    
    % --- extract stack and microscope info from meta data
    omeMeta = reader.getMetadataStore();
    
    numChannels = omeMeta.getChannelCount(0);
    numImages = omeMeta.getImageCount();
    
    % --- get the voxel edge sizes
    
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0);
    voxelSizeX = voxelSizeX.value(ome.units.UNITS.MICROM);
    voxelSizeX = voxelSizeX.doubleValue();
    voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0);
    voxelSizeY = voxelSizeY.value(ome.units.UNITS.MICROM);
    voxelSizeY = voxelSizeY.doubleValue();

    voxelSizes(ff) = voxelSizeX;
    
    % ---
    rawStackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    rawStackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    rawStackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % image height, pixels

    pixelArea = voxelSizeX.*voxelSizeY;

    nucArea_cell{ff} = zeros(numChannels,numImages);
    nucInt_cell{ff} = zeros(numChannels,numImages);
    neighbor_matching_cell{ff} = zeros(numChannels,numImages);

    nucImage_cell{ff} = cell(numChannels,numImages);
    
    parfor_progress(numImages);
    
    for gg = 1:numImages
        
        reader.setSeries(gg-1);
        
        % Read in the raw stack
        
        rawStack = cell(1,numChannels);
        
        for cc = 1:numChannels
            
            rawStack{cc} = ...
                zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ);
            
            for ll = 1:rawStackSizeZ
                
                % Direct assigment
                planeInd = reader.getIndex(ll-1,cc-1,0)+1;
                rawStack{cc}(:,:,ll) = bfGetPlane(reader,planeInd);
                
            end
            
        end
                            
        [nuc_int,nuc_area,neighbor_matching,AClength,nucImage] = ...
            analyzeSingleStack(rawStack,segChannel,pixelArea,...
            minArea,maxArea);
        
        nucArea_cell{ff}(gg) = nuc_area;
        nucInt_cell{ff}(:,gg) = nuc_int;
        neighbor_matching_cell{ff}(:,gg) = neighbor_matching;
        AClength_cell{ff}(:,gg) = AClength;
        if ~isempty(nucImage)
            nucImage_cell{ff}(:,gg) = nucImage;
        end
        
        parfor_progress;
        
    end
        
    reader.close();

    parfor_progress(0);

end


%% -- plot distributions

figure(1)

clf

txn_division = 1500;
plotString = {'ko','r+'};

pooledMatching = {[],[]};
pooledpolII_Int = {[],[]};
pooledRNA_Int = {[],[]};
pooledAClength = {[],[]};

allMatching = [];
allPolII_Int = [];
allRNA_Int = [];
allAClength = [];

for kk = 1:numSets
    
    [~,numImages] = size(neighbor_matching_cell{kk});

    allMatching = [allMatching,neighbor_matching_cell{kk}(3,:)];
    allPolII_Int = [allPolII_Int,nucInt_cell{kk}(2,:)];
    allRNA_Int = [allRNA_Int,nucInt_cell{kk}(3,:)];
    allAClength = [allAClength,AClength_cell{kk}(3,:)];
    
    if cond_vec(kk) == 0
        thisPlotStyle = plotString{1};
        pooledpolII_Int{1} = [pooledpolII_Int{1},nucInt_cell{kk}(2,:)];
        pooledRNA_Int{1} = [pooledRNA_Int{1},nucInt_cell{kk}(1,:)];
        pooledMatching{1} = ...
            [pooledMatching{1},neighbor_matching_cell{kk}(3,:)];
        pooledAClength{1} = ...
            [pooledAClength{1},AClength_cell{kk}(3,:)];
    else
        thisPlotStyle = plotString{2};
        pooledpolII_Int{2} = [pooledpolII_Int{2},nucInt_cell{kk}(2,:)];
        pooledRNA_Int{2} = [pooledRNA_Int{2},nucInt_cell{kk}(1,:)];
        pooledMatching{2} = ...
            [pooledMatching{2},neighbor_matching_cell{kk}(3,:)];
        pooledAClength{2} = ...
            [pooledAClength{2},AClength_cell{kk}(3,:)];
    end
    
end

keepMask = ~isnan(allMatching) & ~isnan(allPolII_Int) ...
    & ~isnan(allRNA_Int) & ~isnan(allAClength);

allMatching = allMatching(keepMask);
allPolII_Int = allPolII_Int(keepMask);
allRNA_Int = allRNA_Int(keepMask);
allAClength = allAClength(keepMask);

CTRL_PolII_int = pooledpolII_Int{1};
CTRL_RNA_int = pooledRNA_Int{1};
CTRL_matching = pooledMatching{1};
CTRL_AC_length = pooledAClength{1};

CTRL_PolII_int = ...
    CTRL_PolII_int(~isnan(CTRL_PolII_int));
CTRL_RNA_int = ...
    CTRL_RNA_int(~isnan(CTRL_RNA_int));
CTRL_matching = ...
    CTRL_matching(~isnan(CTRL_matching));
CTRL_AC_length = ...
    CTRL_AC_length(~isnan(CTRL_AC_length));


includeInds = CTRL_PolII_int>txn_division;

CTRL_active_PolII_int = ...
    CTRL_PolII_int(includeInds);
CTRL_active_RNA_int = ...
    CTRL_RNA_int(includeInds);
CTRL_active_matching = ...
    CTRL_matching(includeInds);
CTRL_active_AC_length = ...
    CTRL_AC_length(includeInds);



FP_PolII_int = pooledpolII_Int{2};
FP_RNA_int = pooledRNA_Int{2};

FP_AC_length = pooledAClength{2};

FP_PolII_int = ...
    FP_PolII_int(~isnan(FP_PolII_int));
FP_RNA_int = ...
    FP_RNA_int(~isnan(FP_RNA_int));
FP_matching = ...
    FP_matching(~isnan(FP_matching));
FP_AC_length = ...
    FP_AC_length(~isnan(FP_AC_length));


subplot(1,2,1)

plot(CTRL_active_PolII_int,CTRL_active_RNA_int,...
    'ko','MarkerFaceColor',[0,0,0])
hold on
plot(FP_PolII_int,FP_RNA_int,'r+')
hold off

xlabel('Intensity Ser2Phos (a.u.)')
ylabel('Intensity RNA (a.u.)')



subplot(1,2,2)


plot(CTRL_active_PolII_int,...
    CTRL_active_matching,'ko','MarkerFaceColor',[0,0,0])
hold on
plot(FP_PolII_int,FP_matching,'r+')
hold off

xlabel('Intensity RNA (a.u.)')
ylabel('L_{corr} DNA (a.u.)')




%% hkjhlhbli.ub

CTRL_noTxn_int = pooledpolII_Int{1}(pooledpolII_Int{1}<txn_division);
CTRL_withTxn_int = pooledpolII_Int{1}(pooledpolII_Int{1}>=txn_division);

FP_PolII_int = FP_PolII_int(~isnan(FP_PolII_int));
CTRL_noTxn_int = CTRL_noTxn_int(~isnan(CTRL_noTxn_int));
CTRL_withTxn_int = ...
    CTRL_withTxn_int(~isnan(CTRL_withTxn_int));



FP_matching = pooledMatching{2};
CTRL_noTxn_matching = pooledMatching{1}(pooledpolII_Int{1}<txn_division);
CTRL_withTxn_matching = pooledMatching{1}(pooledpolII_Int{1}>=txn_division);

FP_matching = FP_matching(~isnan(FP_matching));
CTRL_noTxn_matching = CTRL_noTxn_matching(~isnan(CTRL_noTxn_matching));
CTRL_withTxn_matching = ...
    CTRL_withTxn_matching(~isnan(CTRL_withTxn_matching));

FP_AC = pooledAClength{2}(pooledpolII_Int{2}<txn_division);
CTRL_noTxn_AC = pooledAClength{1}(pooledpolII_Int{1}<txn_division);
CTRL_withTxn_AC = pooledAClength{1}(pooledpolII_Int{1}>=txn_division);

FP_AC = FP_AC(~isnan(FP_AC));
CTRL_noTxn_AC = CTRL_noTxn_AC(~isnan(CTRL_noTxn_AC));
CTRL_withTxn_AC = ...
    CTRL_withTxn_AC(~isnan(CTRL_withTxn_AC));



% bootstrap confidence intervals

FP_CI_matching = bootci(5000,@(vals)mean(vals),FP_matching);
CTRL_noTxn_CI_matching = ...
    bootci(5000,@(vals)mean(vals),CTRL_noTxn_matching);
CTRL_withTxn_CI_matching = ...
    bootci(5000,@(vals)mean(vals),CTRL_withTxn_matching);



densitySupport = linspace(0,1,500);

[FP_density,~] = ksdensity(FP_matching,densitySupport,'BandWidth',0.1);
[CTRL_noTxn_density,~] = ...
    ksdensity(CTRL_noTxn_matching,densitySupport,'BandWidth',0.1);
[CTRL_withTxn_density,~] = ...
    ksdensity(CTRL_withTxn_matching,densitySupport,'BandWidth',0.1);

binEdges = linspace(0,1,30);
binCenters = binEdges(1:end-1) + 0.5.*(binEdges(2)-binEdges(1));


subplot(3,5,5)

[CTRL_withTxn_count,~] = ...
    histc(CTRL_withTxn_matching,binEdges);
CTRL_withTxn_count = CTRL_withTxn_count(1:end-1);


[YLimVals] = [0,1.5.*max(CTRL_withTxn_count)];

patch([CTRL_withTxn_CI_matching(1),CTRL_withTxn_CI_matching(2),...
    CTRL_withTxn_CI_matching(2),CTRL_withTxn_CI_matching(1)],...
    [YLimVals(1),YLimVals(1),YLimVals(2),YLimVals(2)],...
    [0.65,0.65,0.65],'EdgeColor','none')

hold on

bar(binCenters,CTRL_withTxn_count,1.0,...
    'EdgeColor','none','FaceColor',[0,0,0])

hold off

title('Control, interphase')

ylabel('Frequency')

set(gca,'Box','on')

subplot(3,5,10)

[CTRL_noTxn_count,~] = ...
    histc(CTRL_noTxn_matching,binEdges);
CTRL_noTxn_count = CTRL_noTxn_count(1:end-1);

bar(binCenters,CTRL_noTxn_count,1.0,...
    'EdgeColor','none','FaceColor',[0,0,0])

[YLimVals] = [0,1.5.*max(CTRL_noTxn_count)];

patch([CTRL_noTxn_CI_matching(1),CTRL_noTxn_CI_matching(2),...
    CTRL_noTxn_CI_matching(2),CTRL_noTxn_CI_matching(1)],...
    [YLimVals(1),YLimVals(1),YLimVals(2),YLimVals(2)],...
    [0.65,0.65,0.65],'EdgeColor','none')

hold on

bar(binCenters,CTRL_noTxn_count,1.0,...
    'EdgeColor','none','FaceColor',[0,0,0])

hold off

title('Control, not interphase')

ylabel('Frequency')

set(gca,'Box','on')


subplot(3,5,15)

[FP_count,~] = histc(FP_matching,binEdges);
FP_count = FP_count(1:end-1);

bar(binCenters,FP_count,1.0,...
    'EdgeColor','none','FaceColor',[0,0,0])


[YLimVals] = [0,1.5.*max(FP_count)];

patch([FP_CI_matching(1),FP_CI_matching(2),...
    FP_CI_matching(2),FP_CI_matching(1)],...
    [YLimVals(1),YLimVals(1),YLimVals(2),YLimVals(2)],...
    [0.65,0.65,0.65],'EdgeColor','none')

hold on

bar(binCenters,FP_count,1.0,...
    'EdgeColor','none','FaceColor',[0,0,0])

hold off


title('1 \muM Flavopiridol')

xlabel('DNA demixing')
ylabel('Frequency')

set(gca,'Box','on')



meanFP = mean(FP_matching);
stdFP = std(FP_matching);

meanCTRL_withTxn = mean(CTRL_withTxn_matching);
stdCTRL_withTxn = std(CTRL_withTxn_matching);

meanCTRL_noTxn = mean(CTRL_noTxn_matching);
stdCTRL_noTxn = std(CTRL_noTxn_matching);


% -- linear fit to get slopes

FP_fittedModel = fitlm(FP_PolII_int,FP_matching);
FP_Rsquared = FP_fittedModel.Rsquared;
FP_slope= table2array(FP_fittedModel.Coefficients(2,1));
FP_slope_pValue = table2array(FP_fittedModel.Coefficients(2,4));

CTRL_fittedModel = fitlm(CTRL_withTxn_int,CTRL_withTxn_matching);
CTRL_Rsquared = CTRL_fittedModel.Rsquared;
CTRL_slope = table2array(CTRL_fittedModel.Coefficients(2,1));
CTRL_slope_pValue = table2array(CTRL_fittedModel.Coefficients(2,4));


subplot(1,5,1:3)

plot([1,1].*txn_division,[0,1],'k--')

hold on

FPfit_handle = plot(FP_PolII_int(FP_PolII_int<txn_division),...
    FP_fittedModel.feval(FP_PolII_int(FP_PolII_int<txn_division)),'r-',...
    'LineWidth',1.4);
CTRLfit_handle = plot(CTRL_withTxn_int,...
    CTRL_fittedModel.feval(CTRL_withTxn_int),'k-',...
    'LineWidth',1.4);

FP_handle = plot(FP_PolII_int,FP_matching,'r+','MarkerSize',3);
CTRL_noTxn_handle = ...
    plot(CTRL_noTxn_int,CTRL_noTxn_matching,'ks',...
    'Color',[0.5,0.5,0.5],'MarkerSize',3);
CTRL_withTxn_handle = ...
    plot(CTRL_withTxn_int,CTRL_withTxn_matching,'ko','MarkerSize',3);


hold off

xlabel('Pol II Ser2Phos intensity (a.u.)')
ylabel('DNA demixing')

legend([CTRL_withTxn_handle;CTRL_noTxn_handle;FP_handle;...
    FPfit_handle;CTRLfit_handle],...
    {'Control','Control, no transcription',...
    '1 \muM Flavopiridol',...
    sprintf('Slope: %0.2g (p=%0.2g)',FP_slope,FP_slope_pValue),...
    sprintf('Slope: %0.2g (p=%0.2g)',CTRL_slope,CTRL_slope_pValue)},...
    'EdgeColor','none')


subplot(1,5,4)


mean_vals = [meanCTRL_withTxn,meanFP,meanCTRL_noTxn];
upper_CI = ...
    [CTRL_withTxn_CI_matching(2),...
    FP_CI_matching(2),CTRL_noTxn_CI_matching(2)]-mean_vals;
lower_CI = ...
    -[CTRL_withTxn_CI_matching(1),...
    FP_CI_matching(1),CTRL_noTxn_CI_matching(1)]+mean_vals;

errorbar(1:3,mean_vals,upper_CI,lower_CI,'ko')

set(gca,'XLim',[0.5,3.5],...
    'XTickLabel',{'Ctrl','FP','Ctrl, no txn'})
ylabel('DNA demixing')


%% --- Putting out example images

figure(2)

clf

channel2Max = 5500;

plotSpan = 17;

ExampleCalls = {[11,32],[11,30],[8,15],...
    [10,7],[10,35],[10,9]};
ExampleShifts = {[4.5,5.0],[4.5,4],[5,5]...
    [5.0,3.5],[1.5,0.5],[5,4]};


numExamples = numel(ExampleCalls);

for kk = 1:numExamples
    
    call_kk = ExampleCalls{kk}(1);
    call_ll = ExampleCalls{kk}(2);
    
    plotScale = plotSpan.*voxelSizes(call_kk);
    
    if ~isempty(nucImage_cell{call_kk}{1,call_ll})
        
        subplot(2,numExamples,kk)
        thisImage = nucImage_cell{call_kk}{1,call_ll};
        
        [totalY,totalX] = size(thisImage);
        
        imagesc([0,totalX.*voxelSizes(call_kk)],...
            [0,totalX.*voxelSizes(call_kk)],...
            thisImage./max(thisImage(:)))
        axis tight
        axis equal
        
        set(gca,'XLim',[0,plotSpan]+ExampleShifts{kk}(2),...
            'YLim',[0,plotSpan]+ExampleShifts{kk}(1))
        
        xlabel('')
        ylabel('')
        
        title(sprintf('%2.2f',...
            neighbor_matching_cell{call_kk}(1,call_ll)))
        
        colormap(gray)
        
        set(gca,'YTickLabels',[],'XTickLabels',[])

        subplot(2,numExamples,numExamples+kk)
        thisImage = nucImage_cell{call_kk}{2,call_ll};
        imagesc([0,totalX.*voxelSizes(call_kk)],...
            [0,totalX.*voxelSizes(call_kk)],...
            thisImage./channel2Max,[0,1])
        axis tight
        axis equal
        
        xlabel('')
        ylabel('')
        
        set(gca,'XLim',[0,plotSpan]+ExampleShifts{kk}(2),...
            'YLim',[0,plotSpan]+ExampleShifts{kk}(1))
        
        
        title(sprintf('%0.0f',...
            nucInt_cell{call_kk}(2,call_ll)))
        
        colormap(gray)
        
        set(gca,'YTickLabels',[],'XTickLabels',[])
        
    else
        subplot(2,numExamples,kk)
        cla
        subplot(2,numExamples,numExamples+kk)
        cla
    end
    
end

subplot(2,numExamples,1)
hold on
plot([1.5,6.5]+ExampleShifts{1}(1),...
    [14.5,14.5]+ExampleShifts{1}(2),'w-','LineWidth',5)
hold off