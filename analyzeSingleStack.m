function [nuc_int,nuc_area,neighbor_matching,AC_length,nucImage] = ...
    analyzeSingleStack(rawStack,segChannel,pixelArea,minArea)

plotSegmentation = false;

% --- find center slice by maximum contrast

[rawStackSizeY,rawStackSizeX,rawStackSizeZ] = ...
    size(rawStack{segChannel});

numChannels = numel(rawStack);

contrastVec = zeros(1,rawStackSizeZ);

% use only center 20% square for contrast calculation

yLow = round(rawStackSizeY.*0.4);
yHigh = round(rawStackSizeY.*0.6);

xLow = round(rawStackSizeX.*0.4);
xHigh = round(rawStackSizeX.*0.6);

for ll = 1:rawStackSizeZ
    
    section = squeeze(rawStack{segChannel}(yLow:yHigh,xLow:xHigh,ll));
    
    diffMatrix = (section(1:end-1,1:end-1)...
        - section(2:end,2:end)).^2;
    
    contrastVec(ll) = sqrt(mean(diffMatrix(:)));
    
end

[~,maxContrastInd] = max(contrastVec);

maxContrastImage = ...
    cellfun(@(elmt)squeeze(elmt(:,:,maxContrastInd)),...
    rawStack,'UniformOutput',false);

% --- segmentation

medianFiltImage = ...
    medfilt2(maxContrastImage{segChannel},[5 5]);

maxInt = max(medianFiltImage(:));
binEdges = 0:1:maxInt;

counts = histc(medianFiltImage(:),binEdges);

counts = counts(1:end-1);
binEdges = binEdges(1:end-1);

binEdges = (double(binEdges(counts>0))).';
counts = double(counts(counts>0));

counts = counts./max(counts);

[threshVal,otsuMetric] = otsuLimit(binEdges,counts,[0,Inf]);

segImage = medianFiltImage>=threshVal;

% --- erode/dilate treatment with hole-filling

dilateCycles = 6;
dilateMeasure = 6;
erodeCycles = 6;

% Make neighborhood mask for operations
neighborhood = ...
    [0,1,0;...
    1,1,1;...
    0,1,0];

for kk = 1:dilateCycles
    segImage = imdilate(segImage,neighborhood);
end

% Fill holes in segmented image
segImage = imfill(segImage,'holes');

for kk = 1:erodeCycles
    segImage = imerode(segImage,neighborhood);
end

% Fill holes in segmented image
segImage = imfill(segImage,'holes');

%--- erode/dilate treatment over



if plotSegmentation
    
    figure(1)
    
    subplot(3,2,1)
    
    plot(binEdges,counts,'k-')
    hold on
    plot([1,1].*threshVal,[0,1],'k--')
    hold off
    
    xlabel('Intensity')
    ylabel('Frequency')
    
    subplot(3,2,3)
    
    plot(1:rawStackSizeZ,contrastVec,'k-')
    
    hold on
    
    hold off
    
    xlabel('z section')
    ylabel('Contrast')
    
    
    
    subplot(3,2,2)
    
    imagesc([0,rawStackSizeX],[0,rawStackSizeY],...
        maxContrastImage{segChannel})
    
    axis equal
    
    title('Segmentation')
    
    set(gca,'XLim',[0,rawStackSizeX],'YLim',[0,rawStackSizeY],...
        'YDir','normal')
    
    colormap(gray)
    
    
    
    subplot(3,2,4)
    
    imagesc([0,rawStackSizeX],[0,rawStackSizeY],...
        maxContrastImage{segChannel})
    
    axis equal
    
    title('Quantification')
    
    set(gca,'XLim',[0,rawStackSizeX],'YLim',[0,rawStackSizeY],...
        'YDir','normal')
    
    colormap(gray)
    
    
    
    
    subplot(3,2,6)
    
    imagesc([0,rawStackSizeX],[0,rawStackSizeY],...
        maxContrastImage{segChannel}.*segImage)
    
    axis equal
    
    set(gca,'XLim',[0,rawStackSizeX],'YLim',[0,rawStackSizeY],...
        'YDir','normal')
    
    title(sprintf('%3.3f um',sum(segImage(:))))
    
    colormap(jet)
    
end


% --- Segmentation and quantification

segIntImg = maxContrastImage;
segBinImg = segImage;

% --- segment out nucleus closest to center of slice

regions = bwconncomp(segBinImg);
regionProps = regionprops(regions,'Area');
regionAreas = [regionProps(:).Area];


% Restrict to objects in area range
validVolInds = regionAreas.*pixelArea>=minArea;
nucleiRegions = struct;
nucleiRegions.Connectivity = regions.Connectivity;
nucleiRegions.ImageSize = regions.ImageSize;
nucleiRegions.NumObjects = sum(validVolInds);
nucleiRegions.PixelIdxList = regions.PixelIdxList(validVolInds);

% Extract further properties only for the valid regions

nucleiProps = regionprops(nucleiRegions,segIntImg{segChannel},...
    'Area','WeightedCentroid','MeanIntensity','BoundingBox','Image');

% Get object closest to center of the image
nucleiCentroid = {nucleiProps(:).WeightedCentroid};
[~,minInd] = min(cellfun(@(elmt)...
    sum((elmt-[rawStackSizeX./2,rawStackSizeY./2]).^2),nucleiCentroid));

if numel(nucleiCentroid)<1
    
    
    subplot(1,2,1)
    imagesc(maxContrastImage{segChannel})
    axis equal
    
    subplot(1,2,2)
    imagesc(segBinImg)
    axis equal
    
    drawnow
    
    neighbor_matching = zeros(1,numChannels);
    neighbor_matching(:) = NaN;
    
    AC_length = zeros(1,numChannels);
    AC_length(:) = NaN;

    nuc_area = NaN;
    nuc_int = zeros(1,numChannels);
    nuc_int(:) = NaN;
    
    nucImage = {};
    
    return;
    
end

maxCentroid = nucleiCentroid{minInd};

maxArea = nucleiProps(minInd).Area;
maxInt = nucleiProps(minInd).MeanIntensity;
maxBBox = nucleiProps(minInd).BoundingBox;
maxImage = nucleiProps(minInd).Image;

nuc_area = maxArea;
nucCent = maxCentroid;

segInds = nucleiRegions.PixelIdxList(minInd);
segInds = segInds{1};

segMask = false(size(segBinImg));
segMask(segInds) = true;

%     % --- intensity histogram of segmented section
%
%     rawIntVals = maxContrastCell{ff}{quantChannel}(segMask);
%     rawIntVals = rawIntVals(:);
%
%     int_var = var(rawIntVals);
%
%     hist_ints = 0:max(rawIntVals);
%     hist_counts = zeros(size(hist_ints));
%
%     for kk = 1:numel(hist_ints)
%
%         hist_counts(kk) = sum(rawIntVals==kk-1);
%
%     end
%
%     [quantThreshVal,quantOtsuMetric] = ...
%         otsuLimit(hist_ints,hist_counts,[0,Inf]);
%
% --- quantification of N/C itensity ratios

nucBBox = maxBBox;

% Make mask to determine cytoplasmic intensity

totalDil = dilateCycles+erodeCycles;

% Determine small bounding box (without dilate)
smallMinY = nucBBox(2)+0.5;
smallMinX = nucBBox(1)+0.5;
smallMaxY = nucBBox(2)+nucBBox(4)-0.5;
smallMaxX = nucBBox(1)+nucBBox(3)-0.5;


% Determine extended bounding box (after dilate)
fullExtMinY = smallMinY-totalDil;
fullExtMinX = smallMinX-totalDil;
fullExtMaxY = smallMaxY+totalDil;
fullExtMaxX = smallMaxX+totalDil;

% Limit extended bounding box to within image limits
extMinY = max(1,fullExtMinY);
yLoDiff = extMinY - fullExtMinY;
extMinX = max(1,fullExtMinX);
xLoDiff = extMinX - fullExtMinX;
extMaxY = min(rawStackSizeY,fullExtMaxY);
yHiDiff = fullExtMaxY - extMaxY;
extMaxX = min(rawStackSizeX,fullExtMaxX);
xHiDiff = fullExtMaxX - extMaxX;

% Extended bounding box size
extSizeY = extMaxY - extMinY + 1;
extSizeX = extMaxX - extMinX + 1;

% Inclusion mask
inclMask = zeros(extSizeY,extSizeX);
inclMask((1+totalDil-yLoDiff):(end-totalDil+yHiDiff),...
    (1+totalDil-xLoDiff):(end-totalDil+xHiDiff))...
    = maxImage;

% Exclusion mask
exclMask = segMask(extMinY:extMaxY,extMinX:extMaxX);

dilateNeighborhood = zeros(3,3);
dilateNeighborhood(2,2) = 1;
dilateNeighborhood(1,2) = 1;
dilateNeighborhood(3,2) = 1;
dilateNeighborhood(2,1) = 1;
dilateNeighborhood(2,3) = 1;

for dd = 1:dilateCycles
    inclMask = imdilate(inclMask,dilateNeighborhood);
    exclMask = imdilate(exclMask,dilateNeighborhood);
end

for dd = 1:dilateMeasure
    inclMask = imdilate(inclMask,dilateNeighborhood);
end

measureMask = inclMask & ~exclMask;

cyto_int = zeros(1,numChannels);
nuc_int = zeros(1,numChannels);

for cc = 1:numChannels
    
    % -- Determine nuclear to cytoplasmic ratios
    
    wholeImage = maxContrastImage{cc};
    
    cutoutImage = wholeImage(extMinY:extMaxY,...
        extMinX:extMaxX);
    
    cyto_int(cc) = mean(cutoutImage(measureMask));
    nuc_int(cc) = mean(wholeImage(segMask));
    
end

addMask = zeros(size(segBinImg));
addMask(extMinY:extMaxY,...
    extMinX:extMaxX) = measureMask;

% shellImage = segBinImg | addMask;




neighbor_matching = zeros(1,numChannels);
AC_length = zeros(1,numChannels);

nucImage = cell(1,numChannels);

for cc = 1:numChannels
    
    % --- quantification of neighbor differences

    quantImg = maxContrastImage{cc};
    
    nucImage{cc} = quantImg.*segMask;
    
    shift = 1;
    
    YMask = segMask(1:end-shift,:) & segMask(1+shift:end,:);
    XMask = segMask(:,1:end-shift) & segMask(:,1+shift:end);
    
    YImg = quantImg(1:end-shift,:);
    XImg = quantImg(:,1:end-shift);
    
    YShiftImg = quantImg(1+shift:end,:);
    XShiftImg = quantImg(:,1+shift:end);
    
    YVals = YImg(YMask);
    XVals = XImg(XMask);
    
    YShiftVals = YShiftImg(YMask);
    XShiftVals = XShiftImg(XMask);
    
    YDiff = mean(abs(YVals-YShiftVals));
    XDiff = mean(abs(XVals-XShiftVals));
    
    YRandDiff = mean(abs(YVals-YShiftVals(randperm(numel(YShiftVals)))));
    XRandDiff = mean(abs(XVals-XShiftVals(randperm(numel(XShiftVals)))));
    
    neighbor_matching(cc) = ...
        ((XRandDiff-XDiff)./XRandDiff ...
        + (YRandDiff-YDiff)./YRandDiff)./2;
    
    % --- quantification of spatial correlation
    
    shiftVec = 0:10;
    
    numShifts = numel(shiftVec);
    XCorr = zeros(size(shiftVec));
    YCorr = zeros(size(shiftVec));
    
    for kk = 1:numShifts
        
        shift = shiftVec(kk);
        
        YMask = segMask(1:end-shift,:) & segMask(1+shift:end,:);
        XMask = segMask(:,1:end-shift) & segMask(:,1+shift:end);
        
        YImg = quantImg(1:end-shift,:);
        XImg = quantImg(:,1:end-shift);
        
        YShiftImg = quantImg(1+shift:end,:);
        XShiftImg = quantImg(:,1+shift:end);
        
        YVals = YImg(YMask);
        XVals = XImg(XMask);
        
        YShiftVals = YShiftImg(YMask);
        XShiftVals = XShiftImg(XMask);
        
        YMean = mean(YVals);
        YVar = var(YVals);
        
        XMean = mean(XVals);
        XVar = var(XVals);
        
        YCorr(kk) = mean((YVals-YMean).*(YShiftVals-YMean))./YVar;
        XCorr(kk) = mean((XVals-XMean).*(XShiftVals-XMean))./XVar;
        
    end
    
    meanACfun = (XCorr+YCorr)./2;
    
    %     firstBelow = find(meanACfun<0,1,'first');
    %     AUC = sum(meanACfun(2:(firstBelow-1)));
    %
    AC_length(cc) = sum(meanACfun(2:end)).*sqrt(pixelArea);
    
end

%
%     if plotSegmentation
%
%         subplot(3,2,5)
%
%                 plot([0,shiftVec(end)],[0,0],'k--','Color',[0.6,0.6,0.6])
%                 hold on
%                 plot(shiftVec,(XCorr{ff}+YCorr{ff})./2)
%         %         area(shiftVec(2:firstBelow-1),meanACfun(2:firstBelow-1))
%                 hold off
%
%                 xlabel('Image shift [pixels]')
%                 ylabel('Autocorrelation function')
%
%                 title(sprintf('Sum till crossing: %3.3f',AC_length(ff)))
%
%                 set(gca,'YLim',[-1,1])
%
%                 %         plot(hist_ints{ff},hist_counts{ff},'k-')
%                 %
%                 %         xlabel('Intensity [AU]')
%                 %         ylabel('Frequency')
%                 %
%                 %         title(sprintf('Variance %3.3f AU',int_var(ff)))
%
% %                 fft_diag = diag(fft_cell{ff});
% %
% %                 area(smooth(abs((fft_diag)),1),1e3)
% %
% %                 set(gca,'YScale','log','YLim',[1e2,1e8])
% %
% %                 disp(sum(log10(abs(fft_diag(10:25)))))

end


%% --- plot FFT results
%
% subplot(1,2,1)
%
% for ff = 1:numSeries
%
%    plot(smooth(abs((fft_cell{ff}(2:rawStackSizeY./2+1,1))),1),'k-')
%
%    set(gca,'YScale','log','XScale','linear')
%
%    hold on
%
%    xlabel('f [pixel^{-1}]')
%    ylabel('A [AU]')
%
% end
%
% hold off
%
% subplot(1,2,2)
%
% plot(nucInt{txnChannel}(validSeries),...
%     cellfun(@(elmt)sum(abs(elmt(10:40,1))),fft_cell(validSeries)),'ko')