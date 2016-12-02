clear all

%% Choose input directory and compile file list

plotSegmentation = false;

dataset = 2;

if dataset == 1
    
    sourceDirectory = ...
        '/Users/hilbert/academia/VastenhouwLab/Flavopiridol_Sphere/FP_Inhibition_Sphere_Nov2016/';
    
    filepath_cell = {...
        'Flavopiridol_EU_RNA_17Nov2016/CTRL_0001.nd',...
        'Flavopiridol_EU_RNA_17Nov2016/CTRL_0002.nd',...
        'Flavopiridol_EU_RNA_17Nov2016/FP_0002.nd',...
        'Flavopiridol_EU_RNA_17Nov2016/FP_0004.nd'};
    
elseif dataset == 2
    
    
    sourceDirectory = ...
        '/Users/hilbert/academia/VastenhouwLab/Flavopiridol_Sphere/FP_Inhibition_Sphere_Nov2016/';
    
    filepath_cell = {...
        'FP_EURNAStaining_HighNumber_25Nov2016/CTRL_0001.nd',...
        'FP_EURNAStaining_HighNumber_25Nov2016/CTRL_0002.nd',...
        'FP_EURNAStaining_HighNumber_25Nov2016/CTRL_0003.nd',...
        'FP_EURNAStaining_HighNumber_25Nov2016/CTRL_0004.nd',...
        'FP_EURNAStaining_HighNumber_25Nov2016/FP_0002.nd'};
    
    
end

seriesExcludeCells = {[],[],[],[],[],[],[],[],[],[],[]};

cond_vec = [0,0,0,0,1]; % [Flavopiridol] in uM

useSets = [1:numel(filepath_cell)];

% useSets = [1,3,4];

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
cytoInt_cell = cell(1,numSets);
CoV_cell = cell(1,numSets);

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
    cytoInt_cell{ff} = zeros(numChannels,numImages);
    CoV_cell{ff} = zeros(numChannels,numImages);
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
                            
        [nuc_int,cyto_int,nuc_area,neighbor_matching,...
            AClength,CoV,nucImage] = ...
            analyzeSingleStack(rawStack,segChannel,pixelArea,...
            minArea,maxArea);
        
        nucArea_cell{ff}(gg) = nuc_area;
        nucInt_cell{ff}(:,gg) = nuc_int;
        cytoInt_cell{ff}(:,gg) = cyto_int;
        neighbor_matching_cell{ff}(:,gg) = neighbor_matching;
        AClength_cell{ff}(:,gg) = AClength;
        CoV_cell{ff}(:,gg) = CoV;
        if ~isempty(nucImage)
            nucImage_cell{ff}(:,gg) = nucImage;
        end
        
        parfor_progress;
        
    end
        
    reader.close();

    parfor_progress(0);

end


%% -- plot distributions

if dataset == 1

mitotic_cell_inds = {[8,50,118],...
    [8,35,88,102],[],[]};

exclude_prophase_inds = {[19,25,86,131],[],[],[]};

elseif dataset == 2;
    
    mitotic_cell_inds = {...
        [4,10,87,145,151,158,174,190,219],...
        [9,14,56,67], ...
        [9,39,58,80,122,169,185],...
        [7,9,16,21,26,41,54,68,86,88,95,99,113,122,123,246,249,273],[]};
   
    exclude_prophase_inds = {...
        [47,118,119,120],...
        [],...
        [35,79,156,233],...
        [131,133,169,170,182,187],[]};
    
end

figure(1)

clf

txn_division = 800;
plotString = {'ko','r+'};

pooledMatching = {[],[]};
pooledDNA_Int = {[],[]};
pooledpolII_Int = {[],[]};
pooledpolII_cytoInt = {[],[]};
pooledRNA_Int = {[],[]};
pooledRNA_cytoInt = {[],[]};
pooledAClength = {[],[]};
pooledCoV = {[],[]};

mitoticMatching = [];
mitoticDNA_Int = [];
mitoticpolII_Int = [];
mitoticpolII_cytoInt = [];
mitoticRNA_Int = [];
mitoticRNA_cytoInt = [];
mitoticAClength = [];
mitoticCoV = [];

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
    
    [~,numCells] = size(nucImage_cell{kk});
    interphase_mask = true(1,numCells);
    interphase_mask(mitotic_cell_inds{kk}) = false;
    interphase_mask(exclude_prophase_inds{kk}) = false;

    if cond_vec(kk) == 0
        thisPlotStyle = plotString{1};
        pooledDNA_Int{1} = [pooledDNA_Int{1},nucInt_cell{kk}(3,...
            interphase_mask)];
        pooledpolII_Int{1} = [pooledpolII_Int{1},nucInt_cell{kk}(2,...
            interphase_mask)];
        pooledRNA_Int{1} = [pooledRNA_Int{1},nucInt_cell{kk}(1,...
            interphase_mask)];
        pooledpolII_cytoInt{1} = ...
            [pooledpolII_cytoInt{1},cytoInt_cell{kk}(2,...
            interphase_mask)];
        pooledRNA_cytoInt{1} = ...
            [pooledRNA_cytoInt{1},cytoInt_cell{kk}(1,...
            interphase_mask)];
        pooledMatching{1} = ...
            [pooledMatching{1},neighbor_matching_cell{kk}(3,...
            interphase_mask)];
        pooledAClength{1} = ...
            [pooledAClength{1},AClength_cell{kk}(3,...
            interphase_mask)];
        pooledCoV{1} = [pooledCoV{1},CoV_cell{kk}(3,...
            interphase_mask)];
    else
        thisPlotStyle = plotString{2};
        pooledDNA_Int{2} = [pooledDNA_Int{2},nucInt_cell{kk}(3,...
            interphase_mask)];
        pooledpolII_Int{2} = [pooledpolII_Int{2},nucInt_cell{kk}(2,...
            interphase_mask)];
        pooledRNA_Int{2} = [pooledRNA_Int{2},nucInt_cell{kk}(1,...
            interphase_mask)];
        pooledpolII_cytoInt{2} = ...
            [pooledpolII_cytoInt{2},cytoInt_cell{kk}(2,...
            interphase_mask)];
        pooledRNA_cytoInt{2} = ...
            [pooledRNA_cytoInt{2},cytoInt_cell{kk}(1,...
            interphase_mask)];
        pooledMatching{2} = ...
            [pooledMatching{2},neighbor_matching_cell{kk}(3,...
            interphase_mask)];
        pooledAClength{2} = ...
            [pooledAClength{2},AClength_cell{kk}(3,...
            interphase_mask)];
        pooledCoV{2} = [pooledCoV{2},CoV_cell{kk}(3,...
            interphase_mask)];
    end
    
    mitoticDNA_Int = [mitoticDNA_Int,nucInt_cell{kk}(3,...
        mitotic_cell_inds{kk})];
    mitoticpolII_Int = [mitoticpolII_Int,nucInt_cell{kk}(2,...
        mitotic_cell_inds{kk})];
    mitoticRNA_Int = [mitoticRNA_Int,nucInt_cell{kk}(1,...
        mitotic_cell_inds{kk})];
    mitoticpolII_cytoInt = ...
        [mitoticpolII_cytoInt,cytoInt_cell{kk}(2,...
        mitotic_cell_inds{kk})];
    mitoticRNA_cytoInt = ...
        [mitoticRNA_cytoInt,cytoInt_cell{kk}(1,...
        mitotic_cell_inds{kk})];
    mitoticMatching = ...
        [mitoticMatching,neighbor_matching_cell{kk}(3,...
        mitotic_cell_inds{kk})];
    mitoticAClength = ...
        [mitoticAClength,AClength_cell{kk}(3,...
        mitotic_cell_inds{kk})];
    mitoticCoV = [mitoticCoV,CoV_cell{kk}(3,...
        mitotic_cell_inds{kk})];
    
end

keepMask = ~isnan(allMatching) & ~isnan(allPolII_Int) ...
    & ~isnan(allRNA_Int) & ~isnan(allAClength);

allMatching = allMatching(keepMask);
allPolII_Int = allPolII_Int(keepMask);
allRNA_Int = allRNA_Int(keepMask);
allAClength = allAClength(keepMask);

CTRL_DNA_int = pooledDNA_Int{1};
CTRL_PolII_int = pooledpolII_Int{1}-1.*pooledpolII_cytoInt{1};
CTRL_RNA_int = pooledRNA_Int{1}-1.*pooledRNA_cytoInt{1};
CTRL_matching = pooledMatching{1};
CTRL_AC_length = pooledAClength{1};
CTRL_CoV = pooledCoV{1};

CTRL_DNA_int = ...
    CTRL_DNA_int(~isnan(CTRL_DNA_int));
CTRL_PolII_int = ...
    CTRL_PolII_int(~isnan(CTRL_PolII_int));
CTRL_RNA_int = ...
    CTRL_RNA_int(~isnan(CTRL_RNA_int));
CTRL_matching = ...
    CTRL_matching(~isnan(CTRL_matching));
CTRL_AC_length = ...
    CTRL_AC_length(~isnan(CTRL_AC_length));
CTRL_CoV = ...
    CTRL_CoV(~isnan(CTRL_CoV));

includeInds = CTRL_PolII_int>txn_division;

CTRL_active_DNA_int = ...
    CTRL_DNA_int(includeInds);
CTRL_active_PolII_int = ...
    CTRL_PolII_int(includeInds);
CTRL_active_RNA_int = ...
    CTRL_RNA_int(includeInds);
CTRL_active_matching = ...
    CTRL_matching(includeInds);
CTRL_active_AC_length = ...
    CTRL_AC_length(includeInds);
CTRL_active_CoV = ...
    CTRL_CoV(includeInds);

includeInds = CTRL_PolII_int<=txn_division;

CTRL_noTxn_DNA_int = ...
    CTRL_DNA_int(includeInds);
CTRL_noTxn_PolII_int = ...
    CTRL_PolII_int(includeInds);
CTRL_noTxn_RNA_int = ...
    CTRL_RNA_int(includeInds);
CTRL_noTxn_matching = ...
    CTRL_matching(includeInds);
CTRL_noTxn_AC_length = ...
    CTRL_AC_length(includeInds);
CTRL_noTxn_CoV = ...
    CTRL_CoV(includeInds);

FP_DNA_int = pooledDNA_Int{2};
FP_PolII_int = pooledpolII_Int{2}-1.*pooledpolII_cytoInt{2};
FP_RNA_int = pooledRNA_Int{2}-1.*pooledRNA_cytoInt{2};
FP_matching = pooledMatching{2};
FP_AC_length = pooledAClength{2};
FP_CoV = pooledCoV{2};

FP_DNA_int = ...
    FP_DNA_int(~isnan(FP_DNA_int));
FP_PolII_int = ...
    FP_PolII_int(~isnan(FP_PolII_int));
FP_RNA_int = ...
    FP_RNA_int(~isnan(FP_RNA_int));
FP_matching = ...
    FP_matching(~isnan(FP_matching));
FP_AC_length = ...
    FP_AC_length(~isnan(FP_AC_length));
FP_CoV = FP_CoV(~isnan(FP_CoV));



subplot(2,2,1)

cla

pdfParams = struct();
pdfParams.h = [200,50];

pp = gkde2([CTRL_PolII_int.',CTRL_RNA_int.'],...
    pdfParams);

contourf(pp.x,pp.y,pp.pdf,'LineColor','none')
colormap([linspace(0.5,0,32).',...
    linspace(0.5,0,32).',linspace(0.5,0,32).'])
modColorMap = flipud(gray);
modColorMap = modColorMap(1:20,:);
colormap(modColorMap)
hold on

interphase_handle = plot(CTRL_PolII_int,CTRL_RNA_int,...
    'k.','MarkerFaceColor','none',...
    'MarkerEdgeColor',[0,0,0],'MarkerSize',2);

hold on

mitotic_handle = plot(mitoticpolII_Int-mitoticpolII_cytoInt,...
    mitoticRNA_Int-mitoticRNA_cytoInt,...
    'bo','MarkerFaceColor','none',...
    'MarkerEdgeColor',[1,0,1],'MarkerSize',3);

legend([interphase_handle;mitotic_handle],{'Interphase','Mitotic'},...
    'EdgeColor','none','Color','none','Location','NorthWest')

xlabel('\DeltaI_{nucleus} Pol II Ser2Phos (a.u.)')
ylabel('\DeltaI_{nucleus} RNA (a.u.)')

title('Control')

carryOverXLim = get(gca,'XLim');
carryOverYLim = get(gca,'YLim');

subplot(2,2,2)

pp = gkde2([FP_PolII_int.',FP_RNA_int.'],...
    pdfParams);

contourf(pp.x,pp.y,pp.pdf,'LineColor','none')

hold on

plot(FP_PolII_int,FP_RNA_int,'r.',...
    'MarkerEdgeColor',[1,0,0],'MarkerSize',2)

xlabel('\DeltaI_{nucleus} Pol II Ser2Phos (a.u.)')
ylabel('\DeltaI_{nucleus} RNA (a.u.)')

title('1 \muM Flavopiridol')

set(gca,'XLim',carryOverXLim,'YLim',carryOverYLim)

subplot(2,2,3)


if dataset == 1

    rna_min_val = 800;
    rna_max_val = 1800;

elseif dataset == 2
    
    rna_min_val = 250;
    rna_max_val = 800;
    
end

selection_area_handle = plot(...
    [rna_min_val,rna_min_val,rna_max_val,rna_max_val,rna_min_val],...
    [0.11,0.2,0.2,0.11,0.11],'k-',...
    'Color',[0.65,0.65,0.65]);

hold on

% --- do running window percentiles

prctl_vals = [25,75];

min_count = 20;
window_width = 400;
numWindows = 40;
window_edges_left = linspace( ...
    min(min(CTRL_RNA_int),min(FP_RNA_int)), ...
    max(max(CTRL_RNA_int),max(FP_RNA_int))-window_width, ...
    numWindows);
window_edges_right = window_edges_left+window_width;
window_centers = (window_edges_left+window_edges_right)./2;

FP_mean = zeros(1,numWindows);
FP_SEM = zeros(1,numWindows);
FP_count = zeros(1,numWindows);
CTRL_mean = zeros(1,numWindows);
CTRL_SEM = zeros(1,numWindows);
CTRL_count = zeros(1,numWindows);

num_prctl_vals = numel(prctl_vals);
CTRL_prctls = zeros(num_prctl_vals,numWindows);
FP_prctls = zeros(num_prctl_vals,numWindows);

for ww = 1:numWindows
   
    includeMask = CTRL_RNA_int>=window_edges_left(ww) ...
        & CTRL_RNA_int<window_edges_right(ww);

    CTRL_count(ww) = sum(includeMask);
    
    CTRL_mean(ww) = mean(CTRL_CoV(includeMask));
    CTRL_SEM(ww) = std(CTRL_CoV(includeMask))./sqrt(CTRL_count(ww));
    CTRL_prctls(:,ww) = prctile(CTRL_CoV(includeMask),prctl_vals);
    
    includeMask = FP_RNA_int>=window_edges_left(ww) ...
        & FP_RNA_int<window_edges_right(ww);

    FP_count(ww) = sum(includeMask);
    
    FP_mean(ww) = mean(FP_CoV(includeMask));
    FP_SEM(ww) = std(FP_CoV(includeMask))./sqrt(FP_count(ww));
    FP_prctls(:,ww) = prctile(FP_CoV(includeMask),prctl_vals);
end

CTRL_handle = errorbar(window_centers(CTRL_count>=min_count),...
    CTRL_mean(CTRL_count>=min_count),...
    CTRL_SEM(CTRL_count>=min_count),'k-','LineWidth',1.0);

hold on

FP_handle = errorbar(window_centers(FP_count>min_count),...
    FP_mean(FP_count>min_count),...
    FP_SEM(FP_count>min_count),...
    'r-','LineWidth',1.0);

mit_handle = errorbar(mean(mitoticRNA_Int-mitoticRNA_cytoInt),...
    mean(mitoticCoV),std(mitoticCoV)./sqrt(numel(mitoticCoV)),'b^');

hold off

legend([CTRL_handle,FP_handle,mit_handle],...
    {'Control','1 \mum FP','Mitotic'},...
    'EdgeColor','none','Color','none')

xlabel('\Delta I_{nucleus} RNA (a.u.)')
ylabel('CoV DNA intensity')


% legend([CTRL_handle;FP_handle;...
%     FPfit_handle;CTRLfit_handle],...
%     {'Control','1 \muM Flavopiridol',...
%     sprintf('Slope: %0.2g (p=%0.2g)',FP_slope,FP_slope_pValue),...
%     sprintf('Slope: %0.2g (p=%0.2g)',CTRL_slope,CTRL_slope_pValue)},...
%     'EdgeColor','none')


%%%%%%%%%%%%%%%%%%%%
% --- dependence of correlation length on Pol II activity


CTRL_PolII_in = CTRL_PolII_int(...
    CTRL_RNA_int>rna_min_val & CTRL_RNA_int<rna_max_val);
CTRL_AC_Length_in = CTRL_AC_length(...
    CTRL_RNA_int>rna_min_val & CTRL_RNA_int<rna_max_val);
CTRL_CoV_in = CTRL_CoV(...
    CTRL_RNA_int>rna_min_val & CTRL_RNA_int<rna_max_val);

FP_PolII_in = FP_PolII_int(...
    FP_RNA_int>rna_min_val & FP_RNA_int<rna_max_val);
FP_AC_Length_in = FP_AC_length(...
    FP_RNA_int>rna_min_val & FP_RNA_int<rna_max_val);
FP_CoV_in = FP_CoV( ...
    FP_RNA_int>rna_min_val & FP_RNA_int<rna_max_val);




subplot(2,2,4)

% --- do running window percentiles

prctl_vals = [25,75];

min_count = 10;
window_width = 500;
numWindows = 80;
window_edges_left = linspace( ...
    min(min(CTRL_PolII_int),min(FP_PolII_int)), ...
    max(max(CTRL_PolII_int),max(FP_PolII_int))-window_width, ...
    numWindows);
window_edges_right = window_edges_left+window_width;
window_centers = (window_edges_left+window_edges_right)./2;

FP_mean = zeros(1,numWindows);
FP_SEM = zeros(1,numWindows);
FP_count = zeros(1,numWindows);
CTRL_mean = zeros(1,numWindows);
CTRL_SEM = zeros(1,numWindows);
CTRL_count = zeros(1,numWindows);

num_prctl_vals = numel(prctl_vals);
CTRL_prctls = zeros(num_prctl_vals,numWindows);
FP_prctls = zeros(num_prctl_vals,numWindows);

for ww = 1:numWindows
   
    includeMask = CTRL_PolII_in>=window_edges_left(ww) ...
        & CTRL_PolII_in<window_edges_right(ww);

    CTRL_count(ww) = sum(includeMask);
    
    CTRL_mean(ww) = mean(CTRL_CoV_in(includeMask));
    CTRL_SEM(ww) = std(CTRL_CoV_in(includeMask))./sqrt(CTRL_count(ww));
    CTRL_prctls(:,ww) = prctile(CTRL_CoV_in(includeMask),prctl_vals);
    
    includeMask = FP_PolII_in>=window_edges_left(ww) ...
        & FP_PolII_in<window_edges_right(ww);

    FP_count(ww) = sum(includeMask);
    
    FP_mean(ww) = mean(FP_CoV_in(includeMask));
    FP_SEM(ww) = std(FP_CoV_in(includeMask))./sqrt(FP_count(ww));
    FP_prctls(:,ww) = prctile(FP_CoV_in(includeMask),prctl_vals);
end

errorbar(window_centers(CTRL_count>=min_count),...
    CTRL_mean(CTRL_count>=min_count),...
    CTRL_SEM(CTRL_count>=min_count),'k-','LineWidth',1.0)

hold on

errorbar(window_centers(FP_count>min_count),...
    FP_mean(FP_count>min_count),...
    FP_SEM(FP_count>min_count),...
    'r-','LineWidth',1.0)

hold off

legend({'Control','1 \mum FP'},...
    'EdgeColor','none','Color','none')

xlabel('\DeltaI_{nucleus} Pol II Ser2Phos (a.u.)')
ylabel('CoV DNA intensity')


% legend([CTRL_handle;FP_handle;...
%     FPfit_handle;CTRLfit_handle],...
%     {'Control','1 \muM Flavopiridol',...
%     sprintf('Slope: %0.2g (p=%0.2g)',FP_slope,FP_slope_pValue),...
%     sprintf('Slope: %0.2g (p=%0.2g)',CTRL_slope,CTRL_slope_pValue)},...
%     'EdgeColor','none')



%% -- plot example images

figure(2)

clf

image_box = [180,180];

example_cells = {[1,175],[1,7],[3,6],[3,61],[1,8]};
example_offsets = {[40,50],[55,110],[51,56],[13,35],[36,60]};
%4,58,81,92,93,99,146,175
% [8,35,88,102]
% 51, 60, 61
num_examples = numel(example_cells);

for kk = 1:num_examples
    
      
    subplot(numChannels,num_examples,...
        kk)
    
    showImg = nucImage_cell{example_cells{kk}(1)}{1,...
        example_cells{kk}(2)}(...
        example_offsets{kk}(1)+(1:image_box(1)),...
        example_offsets{kk}(2)+(1:image_box(2)));
    imagesc([0,image_box(2).*voxelSizeX],[0,image_box(1).*voxelSizeY],...
        showImg,[1500,6000])
%     title(...
%         sprintf('\\DeltaI_{nucleus} %0.0f a.u.',...
%         nucInt_cell{example_cells{kk}(1)}(1,example_cells{kk}(2))))
    colormap(gray)
    
    xlabel('')
    set(gca,'XTickLabel',[])

    if kk == 1
        ylabel('y [\mum]')
    else
        ylabel('')
        set(gca,'YTickLabel',[])
    end
    axis equal
    axis tight
    
    
    subplot(numChannels,num_examples,...
        kk+num_examples)
    
    showImg = nucImage_cell{example_cells{kk}(1)}{2,...
        example_cells{kk}(2)}(...
        example_offsets{kk}(1)+(1:image_box(1)),...
        example_offsets{kk}(2)+(1:image_box(2)));
    imagesc([0,image_box(2).*voxelSizeX],[0,image_box(1).*voxelSizeY],...
        showImg,[500,7000])
%     title(...
%         sprintf('\\DeltaI_{nucleus} %0.0f a.u.',...
%         nucInt_cell{example_cells{kk}(1)}(2,example_cells{kk}(2))))
    colormap(gray)
    
    xlabel('')
    set(gca,'XTickLabel',[])

    if kk == 1
        ylabel('y [\mum]')
    else
        ylabel('')
        set(gca,'YTickLabel',[])
    end
    axis equal
    axis tight
    
    
    subplot(numChannels,num_examples,...
        kk+2.*num_examples)
    
    showImg = nucImage_cell{example_cells{kk}(1)}{3,...
        example_cells{kk}(2)}(...
        example_offsets{kk}(1)+(1:image_box(1)),...
        example_offsets{kk}(2)+(1:image_box(2)));
    imagesc([0,image_box(2).*voxelSizeX],[0,image_box(1).*voxelSizeY],...
        showImg)
    title(...
        sprintf('CoV DNA %3.3f',...
        CoV_cell{example_cells{kk}(1)}(3,example_cells{kk}(2))))
    colormap(gray)
    
    xlabel('x [\mum]')
    if kk == 1
        ylabel('y [\mum]')
    else
        ylabel('')
        set(gca,'YTickLabel',[])
    end
    
    axis equal
    axis tight
    
    
    
    
    
end







%% -- Investigate individual nucleus analysis

figure(3)

clf

cond = 4;
nucleus = 1;

for kk = 1:numel(nucImage_cell{cond})

    for cc = 1:numChannels
        
        subplot(1,numChannels,cc)
        
        showImg = nucImage_cell{cond}{cc,kk};
        imagesc(showImg)
        title(...
            sprintf('CoV DNA %3.3f',CoV_cell{cond}(cc,kk)))
        colormap(gray)
        
        axis equal
        axis tight
        
    end

    fprintf('Condition %d, cell %d\n',cond,kk)
    
    waitforbuttonpress

end




