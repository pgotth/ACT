%% Punkt 1: Ziel: Zeigen, dass ACT-Schnipsel  ungefähr gleichgroß - Correlation Function Approach
MicronsPerPixel_TMAB = 0.2495;
addpath(genpath('C:\Users\PabiG\Documents\MATLAB'))
addpath(genpath('C:\Users\PabiG\Desktop\MasterThesis\MATLAB'))
addpath(genpath(('C:\Users\PabiG\Desktop\Bandscheiben - Biochemie\Code Jürgen')))
% cd('C:\Users\PabiG\Desktop\ACT Paper Auswertung\TMA_B_MM#3_SingleCores\MM3b_5')
TMA_B_Dir_3 = 'C:\Users\PabiG\Desktop\ACT Paper Auswertung\TMA_B_MM#3_SingleCores\Alle';

All_Images_3 = strcat(TMA_B_Dir_3, "\*.png");
All_Images_3 = dir(All_Images_3);
FirstPic_3 = ((imread(strcat(TMA_B_Dir_3,"\",All_Images_3(1).name))));

TileIndices = round(0.2*size(FirstPic_3, 1)):round(0.8*size(FirstPic_3, 1));


EffectiveLengthScales_3 = zeros(1, length(All_Images_3));
EffectiveLengthScales_Sigma_3 = zeros(1, length(All_Images_3));
SpatialCorrelation = cell(1, length(All_Images_3));


for PicInd = 1:length(All_Images_3)
    ACT_Core = ((imread(strcat(TMA_B_Dir_3,"\",All_Images_3(PicInd).name))));
    ACT_Core = double(ACT_Core)./255;
    ACT_Core = ACT_Core(TileIndices, TileIndices, :);
    Background = zeros(size(ACT_Core, 1), size(ACT_Core, 2));
    Background(mean(ACT_Core, 3)>0.9) = 1;
    Background = bwareaopen(Background, round((MicronsPerPixel_TMAB)^(-1)*4^2*pi));
    Schnipsel = ~imdilate(Background, strel('disk', 2));
    
    
    autoi=Autocorr2D_SG(Schnipsel);
    Autocorr_Schnipsel= UnwrapAutocorr_JL( autoi,1:round(5*(1/MicronsPerPixel_TMAB)));
    P = polyfit(Autocorr_Schnipsel(:,1),Autocorr_Schnipsel(:,2),1);
    estimated_cluster_radius=-Autocorr_Schnipsel(1,2)/P(1)*MicronsPerPixel_TMAB;
    EffectiveLengthScales_3(PicInd) = estimated_cluster_radius;
    EffectiveLengthScales_3_Normed(PicInd) = estimated_cluster_radius./(length(find(Schnipsel>0))/length(find(Schnipsel==0)));
    
    
    
    %calculate spatial "variance" of lengthscales
    ScalesForVariance = zeros(1, 4);
    ImageLength = size(Schnipsel, 1);
    Goto = round(5*(1/MicronsPerPixel_TMAB));
    for i = 1:4
        switch i
            case 1
                ScalesForVariance(i) = calcEffectiveLength(Schnipsel(1:floor(ImageLength/2), 1:floor(ImageLength/2)), Goto, MicronsPerPixel_TMAB);
            case 2
                ScalesForVariance(i) = calcEffectiveLength(Schnipsel(1:floor(ImageLength/2), floor(ImageLength/2):2*floor(ImageLength/2)), Goto, MicronsPerPixel_TMAB);
            case 3
                ScalesForVariance(i) = calcEffectiveLength(Schnipsel(floor(ImageLength/2):2*floor(ImageLength/2), 1:floor(ImageLength/2)), Goto, MicronsPerPixel_TMAB);
            case 4
                ScalesForVariance(i) = calcEffectiveLength(Schnipsel(floor(ImageLength/2):ImageLength, floor(ImageLength/2):ImageLength), Goto, MicronsPerPixel_TMAB);
        end
    end
    EffectiveLengthScales_Sigma_3(PicInd) = std(ScalesForVariance);
    RadiallyAveragedCorrelation = calcSpatialCorrelation(Schnipsel, round(200*(1/MicronsPerPixel_TMAB)), 1000);
    SpatialCorrelation{PicInd} = RadiallyAveragedCorrelation;
    % figure, imshow(Schnipsel)
    % title(num2str(round(estimated_cluster_radius)))
    % ax = gca;
    % ax.FontSize = 16;
end
%
%
% LengthScaleVariance_MM3 = EffectiveLengthScales_Sigma_3;
% LengthScale_MM3 = EffectiveLengthScales_3;
% clearvars -except LengthScaleVariance_MM2 LengthScale_MM2 LengthScaleVariance_MM3 LengthScale_MM3
% plot the length scales

%plot probability of length scales
figure, histogram(LengthScale_MM2, 'Normalization', 'probability')
hold on
histogram(LengthScale_MM3, 'Normalization', 'probability')
legend(strcat('MM2 ($$N_{Core} = 27$$): $$\bar{\xi}=$$ ', num2str(round(mean(LengthScale_MM2))), '$$\mu m$$'), strcat('MM3 ($$N_{Core} = 29$$): $$\bar{\xi} = $$ ', num2str(round(mean(LengthScale_MM3))), '$$\mu m$$'), 'Interpreter','latex')
xlabel('Length Scale $$\xi [ \mu m ] $$', 'Interpreter','latex')
ylabel('Probability')
ax = gca;
ax.FontSize = 16;
%
% figure, histogram(LengthScaleVariance_MM2, 'Normalization', 'probability')
% hold on
% histogram(LengthScaleVariance_MM3, 'Normalization', 'probability')
% legend(strcat('MM2 ($$N_{Core} = 27$$): $$\bar{\sigma}_{\xi}^{Core}=$$ ', num2str(round(mean(LengthScaleVariance_MM2))), '$$\mu m$$'), strcat('MM3 ($$N_{Core} = 29$$): $$ \bar{\sigma}_{\xi}^{Core}= $$ ', num2str(round(mean(LengthScaleVariance_MM3))), '$$\mu m$$'), 'Interpreter','latex')
% xlabel('Intra-Core Length Scale Variance $$\sigma_{\xi}^{Core} [ \mu m ]$$', 'Interpreter', 'latex')
% ylabel('Probability')
% ax = gca;
% ax.FontSize = 16;

% plot CDF
%
% figure, histogram(LengthScale_MM2, 'Normalization', 'cdf')
% hold on
% histogram(LengthScale_MM3, 'Normalization', 'cdf')
% legend(strcat('MM2 ($$N_{Core} = 27$$)'), strcat('MM3 ($$N_{Core} = 29$$)'), 'Interpreter','latex')
% xlabel('Length Scale $$\xi [ \mu m ] $$', 'Interpreter','latex')
% ylabel('CDF')
% ax = gca;
% ax.FontSize = 16;


%% PUNKT 1: Wir behalten die morphologischen Klassen, gezeigt anhand von H&E Bildern von 'ACT und TMA

%% zunächst wird die korrelationsfunktion zwischen TMA und ACT verglichen
addpath(genpath('C:\Users\PabiG\Documents\MATLAB'))
addpath(genpath('C:\Users\PabiG\Desktop\MasterThesis\MATLAB'))
addpath('C:\Users\PabiG\Desktop\Improvement Algo\StainEstimation\Optimize')
Single_TMA_Dir = 'C:\Users\PabiG\Desktop\ACT Paper Auswertung\Punkt1\TMA_A\MM#2_A_HE_2_QP_H_Channel'; %'E:\MM#3_A_HE3_Cores'; 
MicronsPerPixel = 0.3275;
All_Images = strcat(Single_TMA_Dir, "\*.tif");
Single_TMA_Pics = dir(All_Images);
intensity_edges = linspace(0, 1, 51); %60 bins
cluster_edges = [0:2:5000];
Radius = 1000;
ClusterObjectNumber = 4000;
ClusterSizeScale = [0:2:5000];
NumberAngles = 100;
TMA_Correllation = zeros(Radius+1, 1);
ACT_Correllation = zeros(Radius+1, 1);
Numeric_Correllation = zeros(Radius+1, 1);

TMA_ColorHist = zeros(1,length(intensity_edges)-1);
ACT_ColorHist = zeros(1,length(intensity_edges)-1);

TMA_ClusterSizeHist = zeros(1,length(cluster_edges)-1);
ACT_ClusterSizeHist = zeros(1,length(cluster_edges)-1);

TMA_CumulativeClusterSizeSum  = [];
ACT_CumulativeClusterSizeSum  = [];
TMA_SortedClusterSize = [];
ACT_SortedClusterSize = [];


TMA_NoObjects = 0;
ACT_NoObjects = 0;


TMA_count = 0;
ACT_count = 0;

AreaToNeglect = 12*ceil((1/MicronsPerPixel)^2); %12 µm^2 ist gut!!

f = waitbar(0, 'Compare TMA, ACT, Numerical Randomized');



for Pic_Ind = 1:length(Single_TMA_Pics)
    Current_Pic = ((imread(strcat(Single_TMA_Dir,"\",Single_TMA_Pics(Pic_Ind).name))));
    %     if false % -> do colordeconvolution
    %     [Channel_1, Channel_2] = OptimizeColorDeconv(Current_Pic, 4, true);
    
    %     imwrite(Channel_1, strcat('E:\MM#3_A_HE3_Cores\ch1_opt_color_dec\',Single_TMA_Pics(Pic_Ind).name ));
    %     imwrite(Channel_2, strcat('E:\MM#3_A_HE3_Cores\ch2_opt_color_dec\',Single_TMA_Pics(Pic_Ind).name ))
    Current_Pic = double(Current_Pic)./255;
    Current_Pic = nanmean(Current_Pic, 3);
    Current_Pic = Current_Pic./max(Current_Pic(:));
    Current_Pic = adpmedian(Current_Pic, 2*floor(10/2)+1);
    Current_Pic = Current_Pic-imgaussfilt(Current_Pic, 100);
    Current_Pic = Current_Pic./max(Current_Pic(:));


    ThreshedPic = Current_Pic;
    ThreshedPic(Current_Pic < 0.5*graythresh(Current_Pic(Current_Pic > 0.05))) = 0;
    ThreshedPic(~logical(bwareaopen(ThreshedPic, round((1/MicronsPerPixel)^2*pi)))) = 0;
    
    %shuffle pixels and create numerically randomized image
    
    CoreIndices = find(logical(imdilate(ThreshedPic, strel('disk', round(5*(1/MicronsPerPixel))))));
    RandomIndices = randperm(numel(CoreIndices));
    
    NumericRandomized = zeros(size(ThreshedPic));
    NumericRandomized(CoreIndices) = ThreshedPic(CoreIndices(RandomIndices));
    
    
    % pixel correllation
    
%     RadiallyAveragedCorrelation = calcSpatialCorrelation(ThreshedPic, Radius, NumberAngles);
%     RadiallyAveragedCorrelation_Numeric = calcSpatialCorrelation(NumericRandomized, Radius, NumberAngles);
    
    % color histogram
    h = histogram(ThreshedPic(ThreshedPic > 0), 'BinEdges',intensity_edges, 'Normalization', 'probability', 'HandleVisibility', 'off', 'Visible', 'off');
    ColorHist = h.Values;
    close all;
    % cluster sizes
    ClusterPic = Current_Pic;
    ClusterPic(Current_Pic < 0.5*graythresh(Current_Pic(Current_Pic > 0.05))) = 0;
    ClusterPic(~logical(bwareaopen(ClusterPic, AreaToNeglect))) = 0;
    CC = bwconncomp(ClusterPic);
    %% Cumulative Cluster Sum
    ClusterSize=cellfun(@length,CC.PixelIdxList).*(MicronsPerPixel)^2;
    SortedClusterSize=sort(ClusterSize);
    % bin data for
    [Binned_SortedClusterSize_Index] = discretize(SortedClusterSize, ClusterSizeScale);
    SortedClusterSize = ClusterSizeScale(Binned_SortedClusterSize_Index(~isnan(Binned_SortedClusterSize_Index)));
    
    
    
    CurrentClusterNumber = length(SortedClusterSize);
    
    CumulativeClusterSum = [];
    if length(SortedClusterSize) > ClusterObjectNumber
        for i=1:ClusterObjectNumber
            CumulativeClusterSum(i)=sum(SortedClusterSize(1:i))/sum(SortedClusterSize(1:ClusterObjectNumber));
        end
        SortedClusterSize = SortedClusterSize(1:ClusterObjectNumber);
    else
        for i=1:length(SortedClusterSize)
            CumulativeClusterSum(i)=sum(SortedClusterSize(1:i))/sum(SortedClusterSize);
        end
        CumulativeClusterSum = [CumulativeClusterSum, nan(1, ClusterObjectNumber-CurrentClusterNumber)];
        SortedClusterSize = [SortedClusterSize, nan(1, ClusterObjectNumber-CurrentClusterNumber)];
    end
    %%
    
    ClusterSize = histogram(cellfun(@length,CC.PixelIdxList).*(MicronsPerPixel)^2, 'BinEdges',cluster_edges, 'Normalization', 'probability', 'HandleVisibility', 'off');
    ClusterSizeHistogram = ClusterSize.Values;
    NumberOfObjects = CC.NumObjects;
    
    if Pic_Ind > 1 && Pic_Ind < length(Single_TMA_Pics)
%         TMA_Correllation = TMA_Correllation + RadiallyAveragedCorrelation;
%         Numeric_Correllation = Numeric_Correllation + RadiallyAveragedCorrelation_Numeric;
        TMA_ColorHist = TMA_ColorHist + ColorHist;
        TMA_ClusterSizeHist = TMA_ClusterSizeHist + ClusterSizeHistogram;
        TMA_NoObjects = TMA_NoObjects + NumberOfObjects;
        TMA_CumulativeClusterSizeSum  = cat(1, TMA_CumulativeClusterSizeSum, CumulativeClusterSum);
        TMA_SortedClusterSize = cat(1,TMA_SortedClusterSize , SortedClusterSize);
        
        
        TMA_count = TMA_count + 1;
    else
%         ACT_Correllation = ACT_Correllation + RadiallyAveragedCorrelation;
        ACT_ColorHist = ACT_ColorHist + ColorHist;
        ACT_ClusterSizeHist = ACT_ClusterSizeHist + ClusterSizeHistogram;
        ACT_NoObjects = ACT_NoObjects + NumberOfObjects;
        ACT_SortedClusterSize = cat(1,ACT_SortedClusterSize , SortedClusterSize);
        ACT_CumulativeClusterSizeSum  = cat(1, ACT_CumulativeClusterSizeSum, CumulativeClusterSum);
        ACT_count = ACT_count + 1;
        ACT_Pic = ThreshedPic;
    end
    waitbar(Pic_Ind/length(Single_TMA_Pics), f)
end
% Numeric_Correllation = Numeric_Correllation./TMA_count;
% TMA_Correllation = TMA_Correllation./TMA_count;
% ACT_Correllation = ACT_Correllation./ACT_count;
TMA_ColorHist = TMA_ColorHist./TMA_count;
ACT_ColorHist = ACT_ColorHist./ACT_count;
ACT_ClusterSizeHist = ACT_ClusterSizeHist./ACT_count;
TMA_ClusterSizeHist = TMA_ClusterSizeHist./TMA_count;
ACT_NoObjects = ACT_NoObjects./ACT_count;
TMA_NoObjects = TMA_NoObjects./TMA_count;

%% plot correllation
subplot(2,3,1)
plot([1:length(TMA_Correllation)].*MicronsPerPixel, (TMA_Correllation), 'LineWidth', 3)
ylim([0 1])
legend(strcat(" \langle TMA \rangle_{N=", num2str(TMA_count), "}^{H}"))
title('TMA')
ylabel('Spatial Correlation')
xlabel('r [\mum]')
subplot(2,3,2)
plot([1:length(TMA_Correllation)].*MicronsPerPixel, (ACT_Correllation), 'LineWidth', 3)
ylim([0 1])
legend(strcat(" \langle ACT \rangle_{N=", num2str(ACT_count), "}^{H}"))
title('ACT')
ylabel('Spatial Correlation')
xlabel('r [\mum]')
subplot(2,3,3)
plot([1:length(TMA_Correllation)].*MicronsPerPixel, (Numeric_Correllation), 'LineWidth', 3)
ylim([0 1])
title('NR')
ylabel('Spatial Correlation')
xlabel('r [\mum]')
legend(strcat(" \langle NR \rangle_{N=", num2str(TMA_count), "}^{H}"))
subplot(2,3,4)
TMA_Pic = ((imread(strcat(Single_TMA_Dir,"\",Single_TMA_Pics(length(Single_TMA_Pics)-6).name))));
TMA_Pic = mean(double(TMA_Pic)./255,3);
imagesc(TMA_Pic)
axis off;
subplot(2,3,5)
imagesc(ACT_Pic)
axis off;
subplot(2,3,6)
TMA_Pic(TMA_Pic < 0.3*graythresh(TMA_Pic(TMA_Pic > 0.05))) = 0;
TMA_Pic(~logical(bwareaopen(TMA_Pic, ceil((1/MicronsPerPixel)^2*pi)))) = 0;
CoreIndices = find(logical(imdilate(TMA_Pic, strel('disk', round(5*(1/MicronsPerPixel))))));
RandomIndices = randperm(numel(CoreIndices));
NumericRandomized = zeros(size(TMA_Pic));
NumericRandomized(CoreIndices) = ThreshedPic(RandomIndices);
imagesc(NumericRandomized)
axis off;
% cluster radius estimation
% LinearFitUntil = round(5*1/MicronsPerPixel);
% P = polyfit([1:LinearFitUntil].*MicronsPerPixel, (ACT_Correllation(1:LinearFitUntil))',1);
% estimated_cluster_radius=-ACT_Correllation(1)/P(1)*MicronsPerPixel;
% EffectiveLengthScales_3(PicInd) = estimated_cluster_radius;


% cluster sizes
Goto = round(0.5*100);
figure, plot(cluster_edges(1:Goto), (TMA_ClusterSizeHist(1:Goto)), 'LineWidth', 1)
hold on, plot(cluster_edges(1:Goto), (ACT_ClusterSizeHist(1:Goto)), 'LineWidth', 1)
legend(strcat(" \langle TMA \rangle_{N=", num2str(TMA_count), "}^{H}"),strcat(" \langle ACT \rangle_{N=", num2str(ACT_count), "}^{H}"))
xlabel('Cluster Size [\mum^2]')
ylabel('Probability')
title(strcat("$$\int \limits_{0}^{", num2str(Goto*2), "} C_{ACT} = ", num2str(sum(ACT_ClusterSizeHist(1:Goto))), "$$", ";   ","$$\int \limits_{0}^{", num2str(Goto*2), "} C_{TMA} = ", num2str(sum(TMA_ClusterSizeHist(1:Goto))), "$$"), 'interpreter','latex')
ax = gca;
ax.FontSize = 16;
% plot P*<C>
figure, plot(cluster_edges(1:Goto), smooth(TMA_ClusterSizeHist(1:Goto).*cluster_edges(1:Goto)), 'LineWidth', 2)
hold on, plot(cluster_edges(1:Goto), smooth(ACT_ClusterSizeHist(1:Goto).*cluster_edges(1:Goto)), 'LineWidth', 2)
legend(strcat(" \langle TMA \rangle_{N=", num2str(TMA_count), "}^{H}"),strcat(" \langle ACT \rangle_{N=", num2str(ACT_count), "}^{H}"))
xlabel('Cluster Size [\mum^2]')
ylabel('Probability \cdot \langle C \rangle')
ax = gca;
ax.FontSize = 16;
% plot cumulative distribution function
TMA_ClusterSize_CDF = [];
ACT_ClusterSize_CDF = [];
NumberOfBins = 2000;
for i = 1:NumberOfBins
    TMA_ClusterSize_CDF(i) = sum(TMA_ClusterSizeHist(1:i))/sum(TMA_ClusterSizeHist(:));
    ACT_ClusterSize_CDF(i) = sum(ACT_ClusterSizeHist(1:i))/sum(ACT_ClusterSizeHist(:));
end
figure, plot(cluster_edges(1:NumberOfBins), TMA_ClusterSize_CDF)
hold on, plot(cluster_edges(1:NumberOfBins), ACT_ClusterSize_CDF)
% plot cumulative cluster sizes sum

% core-independent cumulative cluster size plots
Linear_TMA_SortedClusterSize = sort(TMA_SortedClusterSize(:));
Linear_TMA_SortedClusterSize(isnan(Linear_TMA_SortedClusterSize)) = [];
Linear_ACT_SortedClusterSize = sort(ACT_SortedClusterSize(:));
Linear_ACT_SortedClusterSize(isnan(Linear_ACT_SortedClusterSize)) = [];
TMA_CumulativeClusterSum = zeros(length(Linear_TMA_SortedClusterSize), 1);
for i=1:length(Linear_TMA_SortedClusterSize)
    TMA_CumulativeClusterSum(i)=sum(Linear_TMA_SortedClusterSize(1:i))/sum(Linear_TMA_SortedClusterSize);
end
ACT_CumulativeClusterSum = zeros(length(Linear_ACT_SortedClusterSize),1);
for i=1:length(Linear_ACT_SortedClusterSize)
    ACT_CumulativeClusterSum(i)=sum(Linear_ACT_SortedClusterSize(1:i))/sum(Linear_ACT_SortedClusterSize);
end
figure('Renderer', 'painters', 'Position', [50 50 800 700])
semilogx( (Linear_TMA_SortedClusterSize), (TMA_CumulativeClusterSum), 'b*')
hold on
semilogx((Linear_ACT_SortedClusterSize), (ACT_CumulativeClusterSum), 'r*')
xlim([0 15000])
xlabel('Cluster Size C [$$\mu m^2$$]', 'Interpreter', 'latex')
ylabel('$$F_C$$', 'Interpreter', 'latex')
legend('$$\langle F_C \rangle_{TMA}^{N=91}$$', '$$\langle F_C \rangle_{ACT}^{N=2}$$', 'Interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;





%% plot cumulative cluster sizes for same scales --> paper fig
figure('Renderer', 'painters', 'Position', [50 50 800 700])
for i = 1:91
    TMA_Handle = plot( (TMA_SortedClusterSize(i, :)),  TMA_CumulativeClusterSizeSum(i, :), 'b--', 'LineWidth', 1);
    xlim([0 5000])

    hold on
    
end
for i = 1:2
    ACT_Handle = plot( (ACT_SortedClusterSize(i, :)),  ACT_CumulativeClusterSizeSum(i, :), 'r-', 'LineWidth', 3.5);
    xlim([0 5000])

    hold on
end
set(gca, 'XScale', 'log')
set(gca,'XTick',([0 10 100 1000 5000]))

legend( [ACT_Handle,TMA_Handle], {'ACT (N=2)', strcat('TMA (N=', num2str(TMA_count), ')')}, 'Interpreter', 'latex')
xlabel('Cluster Size [$$\mu m^2$$]', 'Interpreter', 'latex')
ylabel('Cumulative Sum $$F_C$$', 'Interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;

%% plot cumulative sums of ACT from #2 and #3
figure('Renderer', 'painters', 'Position', [50 50 800 700])
for i = 1:2
    ACT_Handle_MM3 = plot( (ACT_SortedClusterSize(i, :)),  ACT_CumulativeClusterSizeSum(i, :), '-', 'LineWidth', 3.5, 'Color', [217/255 83/255 25/255]);
    xlim([0 5000])
    hold on
end
for i = 1:2
    ACT_Handle_MM2 = plot( (ACT_SortedClusterSize_MM2(i, :)),  ACT_CumulativeClusterSizeSum_MM2(i, :), '-', 'LineWidth', 3.5, 'Color', [0 114/255 189/255]);
    xlim([0 5000])
    hold on
end
set(gca, 'XScale', 'log')
set(gca,'XTick',([0 10 100 1000 5000]))
legend( [ACT_Handle_MM2,ACT_Handle_MM3], {'MM#2 ACT (N=2)', 'MM#3 ACT (N=2)'}, 'Interpreter', 'latex')
xlabel('Cluster Size [$$\mu m^2$$]', 'Interpreter', 'latex')
ylabel('Cumulative Sum $$F_C$$', 'Interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;

%% plot color histograms
figure('Renderer', 'painters', 'Position', [50 50 800 700])
plot((1:50)./50, TMA_ColorHist,'LineWidth', 5, 'Color', [0 0 1])
hold on, plot((1:50)./50, ACT_ColorHist,'LineWidth', 5, 'Color',[1 0 0])
xlabel('Intensity $$I_H$$', 'Interpreter', 'latex')
ylabel('$$P(I_H)$$', 'Interpreter', 'latex')
legend(strcat(" $$\langle I_H \rangle_{N=", num2str(TMA_count), "}^{TMA}$$"),strcat("$$ \langle I_H \rangle_{N=", num2str(ACT_count), "}^{ACT}$$"), 'Interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;




%% PUNKT 2
% 2.	In einem weiteren Schritt sollte dann gezeigt werden, dass die Herstel-
% lung einer solchen Grundmasse reproduzierbar ist, dass also bei Grundmasse
% 2 sechs verschiedene Stanzen aus dieser Grundmasse wiederum in Bezug auf
% Größe und Form einen hohen Ähnlichkeitsgrad aufweisen.
% (TMA B aus Grundmasse #3 = MM#3_TMA B)

%% initialize data from raw images
Mum_per_PXL = 0.249;
Names3 = {'AE3_1', 'CK7', 'E-Cad', 'ER', 'Her2', 'Ki67', 'PR'};
Names2 = {'CK7', 'E-Cad', 'ER', 'Ki67', 'PR'};
Grundmasse = 3;
TresholdToThresh = 0.05;
ThresholdFactor = 0.5;
ColorsFromHowManyCores = 6;
[Store_ACT_TMAB, Store_ACT_Histograms_TMAB]  = InitializeTMAB(ColorsFromHowManyCores,Grundmasse , TresholdToThresh, ThresholdFactor);

%% plot histograms
AllHists = cat(1,Store_ACT_Histograms_TMAB2, Store_ACT_Histograms_TMAB);
figure
for i = 1:size(AllHists, 1)
    subplot(6,2,i)
    for j = 1:6
        helper = AllHists{i,j,1};
        if helper(end) < 0.2
            hold on
            plot([1:length(helper)]./length(helper),helper)
        end
    end
    title(Names{i})
end
% saveas(gcf, 'TMAB_Histograms.png')
% saveas(gcf, 'TMAB_Histograms.fig')

%% Plot properties
LogScale = 0;
figure,
Point2_Plot_twoVariabels((Store_ACT_TMAB), Names2, 8,9, LogScale)

%% Calculate Equivalence
addpath('C:\Users\PabiG\Desktop\MasterThesis\MATLAB')
addpath('C:\Users\PabiG\Desktop\ACT Paper Auswertung\MATLAB')

Names3 = {'AE31', 'CK7', 'E-Cad', 'ER', 'Her2', 'Ki67', 'PR'};
Names2 = {'CK7', 'E-Cad', 'ER', 'Ki67', 'PR'};
%init parameters for equivalence test

Color{1} = [0 0 0]; Color{2} = [1 0 0];Color{3} = [0 1 0];Color{4} = [0 0 1];Color{5} = [1 128/255 0];Color{6} = [0 128/255 1]; Color{7} = [128/255 128/255 1];

% init parameters for Core-Permutations
[M, ~] = permn([1:6], 2);
Permutations = M;
Permutations(Permutations(:,1) == Permutations(:,2), :) = []; %keine gleichen stanzen vergleichen
Permutations = unique(Permutations, 'rows');

Range =  0:0.01:2;
SampleRange = [100]; %for generating normal distributed means

format long
DifferenceIntensityAbsolute_Permutation = zeros(length(Permutations), size(Store_ACT_TMAB, 1), length(Range));
AllPs = cell(length(SampleRange),1);
AllAbsDiff = cell(length(SampleRange));
SampleInd = 0;

f = waitbar(0, 'Calculate Equivalence');
NumberSteps = length(SampleRange)*size(Store_ACT_TMAB, 1)*length(Range)*length(Permutations);
WaitbarIndex = 0;
for SampleSize = SampleRange
    SampleInd = SampleInd + 1;
    Helper_p = [];
    Helper_x = [];
    for StainInd = 1:size(Store_ACT_TMAB, 1)
        EffectSizeInd = 0;
        for min_effectsize_intensity = Range%init
            if min_effectsize_intensity == Range(1)
                PermMean_p = [];
                DiffIntensity = [];
            end
            EffectSizeInd = EffectSizeInd +1;
            pValues_Intensities_Permutation = [];
            DifferenceIntensityAbsolute_Permutation = [];
            for PermIndCores = 1:length(Permutations)
                WaitbarIndex = WaitbarIndex + 1;
                PermInd1 = Permutations(PermIndCores, 1);
                PermInd2 = Permutations(PermIndCores, 2);
                Intensity1 = (Store_ACT_Histograms_TMAB{StainInd,PermInd1,3} );
                Intensity2 = (Store_ACT_Histograms_TMAB{StainInd,PermInd2,3});
                if Intensity1(end)<0.2 && Intensity2(end)<0.2 %kick out false contrast images
                    [Intensity1, ~, ~] = SampleMeans(Intensity1, SampleSize);
                    [Intensity2, ~, ~] = SampleMeans(Intensity2, SampleSize);
                    [~, pValues, pooledSD] = TOST(Intensity1, Intensity2, min_effectsize_intensity);
                    pValues_Intensities_Permutation = [pValues_Intensities_Permutation, max(pValues)];
                    DifferenceIntensityAbsolute_Permutation =[DifferenceIntensityAbsolute_Permutation, pooledSD*min_effectsize_intensity];
                else
                    DifferenceIntensityAbsolute_Permutation =[DifferenceIntensityAbsolute_Permutation, NaN];
                    pValues_Intensities_Permutation = [pValues_Intensities_Permutation, NaN];
                end
                waitbar(WaitbarIndex/NumberSteps, f)
            end
            PermutationMean_p_OneStain = nanmean(pValues_Intensities_Permutation(:));
            DiffInt_p_OneStain = nanmean(DifferenceIntensityAbsolute_Permutation(:));
            PermMean_p  = cat(1,PermMean_p,PermutationMean_p_OneStain);
            DiffIntensity = cat(1,DiffIntensity,DiffInt_p_OneStain);
        end
        Helper_p = cat(2,Helper_p, PermMean_p);
        Helper_x= cat(2,Helper_x, DiffIntensity);
    end
    AllPs{SampleInd}  = Helper_p;
    AllAbsDiff{SampleInd} = Helper_x;
end
if Grundmasse == 2
    P_Values_Different_MM2 =  AllPs{1}; %first dimension: intensity diff, permutation mean, second: Stains, third: sub-sample size
    Absolute_Intensity_Diff_MM2 = AllAbsDiff{1};
elseif Grundmasse == 3
    P_Values_Different_MM3 = cat(3, AllPs{1}); %first dimension: intensity diff, permutation mean, second: Stains, third: sub-sample size
    Absolute_Intensity_Diff_MM3 = cat(3, AllAbsDiff{1});
end
%% generate data for both grundmassen
clearvars -except P_Values_Different_MM2 P_Values_Different_MM3 Absolute_Intensity_Diff_MM2 Absolute_Intensity_Diff_MM3 Names3 Names2
P_Values_Different = cat(2, P_Values_Different_MM2, P_Values_Different_MM3);
Absolute_Intensity_Diff=cat(2, Absolute_Intensity_Diff_MM2, Absolute_Intensity_Diff_MM3);
Names = {'MM2 CK7', 'MM2 E-Cad', 'MM2 ER', 'MM2 Ki67', 'MM2 PR', 'MM3 AE31', 'MM3 CK7', 'MM3 E-Cad', 'MM3 ER', 'MM3 Her2', 'MM3 Ki67', 'MM3 PR'};

%% check with sample means give best results --> the higher the better --> here take 100
% CheckSignificance = abs(P_Values_Different - 0.05);
% [~, I] = min(CheckSignificance, [], 1);
% first = Absolute_Intensity_Diff(:,:,1);
% second = Absolute_Intensity_Diff(:,:,2);
% third = Absolute_Intensity_Diff(:,:,3);
% [~, bestSample] = min([mean(first(I(:,:,1))), mean(second(I(:,:,2))), mean(third(I(:,:,3)))]);
%% Plot p-Values over effectsize and total intensity differences
figure
plot(Absolute_Intensity_Diff(:,:,1), P_Values_Different(:,:,1), 'LineWidth', 3)
legend(Names)
hold on
line( [linspace(0,0.13)],[linspace(0.05,0.05)], 'LineStyle' ,'--','LineWidth', 2, 'Color', 'black' )
ylabel('p (TMA B Permutation Mean)')
xlabel('Absolute Intensity Difference |\mu_1 - \mu_2|')
%% make table
SubSampleIndex = 1;
pValues = P_Values_Different(:,:,SubSampleIndex);
IntDiff = Absolute_Intensity_Diff(:,:,SubSampleIndex);
[~, SignificantEqual_Ind] = min(abs(pValues - 0.05));
Stain = Names';

EqualFrom = zeros(12,1);
for i = 1:12
    EqualFrom(i) = IntDiff(SignificantEqual_Ind(i), i);
end
T = table(Stain, EqualFrom)

%% PUNKT 3

%% PUNKT 4
% 4.	In einem weiteren Ansatz geht es dann darum, die beiden 100er TMAs
% (MM#2_TMA A und MM#3_TMA A) zu analysieren. Der Hintergedanke hinter diesen
% 100er TMAs ist der, dass die Blockpräparate, von denen jeweils ein Schnitt
% zur Grundmasse des entsprechenden ACTs herangezogen wurde, in diesem TMA
% repräsentiert sind. Diese TMAs sind mit insgesamt 10 Antikörpern (ER, PR,
% Her2neu, Ki67, AE1/3, CK7, CD45, Aktin, CK5/6 und Vimentin) gefärbt worden
% (eine entsprechende Liste mit diesen Infos hatte Ihnen bereits Frau Rosin erstellt).
All_TMA_Dir = 'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_A';
ColorDeconvolutionPrecalculated = true;

[Store_ACT, Store_TMA] = calculate_TMA_A_Data(All_TMA_Dir, ColorDeconvolutionPrecalculated);

%%

Names = {'AE13', 'Aktin', 'CD31', 'CD45', 'CK514', 'D240', 'E-Cadherin', 'ER', 'Her2', 'Ki67', 'p63', 'PR', 'Vimentin'};
clearvars -except Store_ACT Store_ACT_Histograms Store_TMA Store_TMA_Histograms Names

Store_ACT_MM2 = Store_ACT;
Store_ACT_Histograms_MM2 = Store_ACT_Histograms;
Store_TMA_MM2 = Store_TMA;
Store_TMA_Histograms_MM2 = Store_TMA_Histograms;
clearvars -except Store_ACT_MM2 Store_ACT_Histograms_MM2 Store_TMA_MM2 Store_TMA_Histograms_MM2 Names

%% ACT Stanzen eines TMAs vergleichen -> mittlere relative Abweichung;
%--> hierfür als referenz: Varianz der TMA Stanzen untereinander  (Varianz der ACTs untereinander sollte deutlcih niedriger sein)
Average_Relative_Deviation = cell(4,1);
for i = 1:4
    Average_Relative_Deviation{i} = .5.*abs(Store_ACT(:,1,i) - Store_ACT(:,2,i))./Store_ACT(:,1,i) + .5.*abs(Store_ACT(:,1,i) - Store_ACT(:,2,i))./Store_ACT(:,2,i);
end
figure(10)
T = table(Average_Relative_Deviation{1},'RowNames',Names');
T.Properties.VariableNames{1} = 'Average Relative Deviation';
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized',  'Position',[0, 0, 1, 1]);

%% plot mittlere cluszerintensität für alle ACT und TMA: 16.12: noch machen
%--> plot mit konfidenzintervall
figure(11)
for i = 1:13
    format long
    CI = Calculate_CI_For_UnknownDistribution(Store_TMA(i,:,3), 0.05);
    line([nanmean(Store_ACT(i,:,3), 2) nanmean(Store_ACT(i,:,3), 2)], [CI(1) CI(2)], 'LineWidth', 2, 'Color', 'black')
    hold on
end
hold on
for i = 1:13
    hold on
    scatter(nanmean(Store_ACT(i,:,3), 2), nanmean(Store_TMA(i,:,3), 2), 50 ,'filled', 'MarkerFaceColor',[1 0 0])
    text(nanmean(Store_ACT(i,:,3), 2), nanmean(Store_TMA(i,:,3), 2), Names{i})
end
xlabel('\langle I_C \rangle_{ACT}^{n=2}')
ylabel('\langle I_C \rangle_{TMA}^{n=96}')
legend('Confidence Interval')
ax = gca;
ax.FontSize = 16;
%add line with slope 1
% hold on
% line([5 9].*10^(-3), [5 9].*10^(-3), 'LineWidth', 2, 'Color', 'blue')




%% plot mittlere Klustergröße
Mum_per_Pixel = 0.3277;
figure
for i = 1:13
    format long
    CI = Calculate_CI_For_UnknownDistribution(Store_TMA(i,:,1).*(Mum_per_Pixel)^2, 0.05);
    %     CI_Plot(CI(1),CI(2),nanmean(Store_ACT(i,:,3), 2));
    
    line([nanmean(Store_ACT(i,:,1), 2).*(Mum_per_Pixel)^2 nanmean(Store_ACT(i,:,1), 2).*(Mum_per_Pixel)^2], [CI(1) CI(2)], 'LineWidth', 2, 'Color', 'black')
    hold on
    %     errorbar(nanmean(Store_ACT(i,:,3), 2), nanmean(Store_TMA(i,:,3)), 2, 'o', 'MarkerSize', 20)
end
hold on
for i = 1:13
    hold on
    scatter(nanmean(Store_ACT(i,:,1), 2).*(Mum_per_Pixel)^2, nanmean(Store_TMA(i,:,1)).*(Mum_per_Pixel)^2, 50 ,'filled','MarkerFaceColor',[1 0 0])
    text(nanmean(Store_ACT(i,:,1), 2).*(Mum_per_Pixel)^2, nanmean(Store_TMA(i,:,1)).*(Mum_per_Pixel)^2, Names{i})
end
xlabel('\langle C \rangle_{ACT}^{n=2} [\mum^2]')
ylabel('\langle C \rangle_{TMA}^{n=96} [\mum^2]')
legend('Confidence Interval')
ax = gca;
ax.FontSize = 16;
% hold on
% line([10 30], [10 30], 'LineWidth', 2, 'Color', 'blue')
% set(gcf, 'Position', [100, 100, 1000, 700]);



%% plot mittlere Klusterverianz
Mum_per_Pixel = 0.3277;
figure
for i = 1:13
    format long
    CI = Calculate_CI_For_UnknownDistribution(Store_TMA(i,:,2).*(Mum_per_Pixel)^2, 0.05);
    %     CI_Plot(CI(1),CI(2),nanmean(Store_ACT(i,:,3), 2));
    
    line([nannanmean(Store_ACT(i,:,2), 2).*(Mum_per_Pixel)^2 nannanmean(Store_ACT(i,:,2), 2).*(Mum_per_Pixel)^2], [CI(1) CI(2)], 'LineWidth', 2, 'Color', 'black')
    hold on
    %     errorbar(nannanmean(Store_ACT(i,:,3), 2), nannanmean(Store_TMA(i,:,3)), 2, 'o', 'MarkerSize', 20)
end
hold on
for i = 1:13
    hold on
    scatter(nannanmean(Store_ACT(i,:,2), 2).*(Mum_per_Pixel)^2, nannanmean(Store_TMA(i,:,2)).*(Mum_per_Pixel)^2, 50 ,'filled', 'MarkerFaceColor',[0 1 0])
    text(nannanmean(Store_ACT(i,:,2), 2).*(Mum_per_Pixel)^2, nannanmean(Store_TMA(i,:,2)).*(Mum_per_Pixel)^2, Names{i})
end
xlabel('\langle \sigma_C \rangle_{ACT}^{n=2} [\mum^2]')
ylabel('\langle \sigma_C \rangle_{TMA}^{n=96} [\mum^2]')
legend('Confidence Interval')
ax = gca;
ax.FontSize = 16;

legend('Confidence Intervall')
% set(gcf, 'Position', [100, 100, 1000, 700]);
hold on
% line([10 120], [10 400], 'LineWidth', 2, 'Color', 'blue')

%% plot VARIANZ DER MITTLEREN KLUSTERINTENSITÄTEN
Mum_per_Pixel = 0.3277;
figure
for i = 1:13
    format long
    CI = Calculate_CI_For_UnknownDistribution(Store_TMA(i,:,4), 0.05);
    %     CI_Plot(CI(1),CI(2),nannanmean(Store_ACT(i,:,3), 2));
    
    line([nannanmean(Store_ACT(i,:,4), 2) nannanmean(Store_ACT(i,:,4), 2)], [CI(1) CI(2)], 'LineWidth', 2, 'Color', 'black')
    hold on
    %     errorbar(nannanmean(Store_ACT(i,:,3), 2), nannanmean(Store_TMA(i,:,3)), 2, 'o', 'MarkerSize', 20)
end
hold on
for i = 1:13
    hold on
    scatter(nanmean(Store_ACT(i,:,4), 2), nanmean(Store_TMA(i,:,4)), 50 ,'filled', 'MarkerFaceColor',[1 0 0])
    text(nanmean(Store_ACT(i,:,4), 2), nanmean(Store_TMA(i,:,4)), Names{i})
end
xlabel('\langle \sigma_{I_C} \rangle_{ACT}^{n=2}')
ylabel('\langle \sigma_{I_C} \rangle_{TMA}^{n=96}')
legend('Confidence Interval')
ax = gca;
ax.FontSize = 16;
hold on
% line([6 12].*10^(-3), [6 12].*10^(-3), 'LineWidth', 2, 'Color', 'blue')
% set(gcf, 'Position', [100, 100, 1000, 700]);

%variabellabel = 1 => nanmean(ClusterSizes);
%             2 => nanstd(ClusterSizes)
%             3 => nanmean(Intensities);
%             4 => nanstd(Intensities);
%             5 => max(ClusterSizes);
%             6 => median(ClusterSizes);
%               7 = mad(ClusterSizes);
%               8) = median(Intensities);
%               9) = mad(Intensities);
TMA_VariabelStrings = {'\langle C \rangle_{TMA}^{n=96}';...
    ' \langle \sigma_C \rangle_{TMA}^{n=96}'; ...
    '\langle I_C \rangle_{TMA}^{n=96}'; ...
    '\langle \sigma_{I_C} \rangle_{TMA}^{n=96}';
    'max\{ C \}_{TMA}^{n=96}'};
ACT_VariabelStrings = {'\langle C \rangle_{ACT}^{n=2}';...
    ' \langle \sigma_C \rangle_{ACT}^{n=2}'; ...
    '\langle I_C \rangle_{ACT}^{n=2}'; ...
    '\langle \sigma_{I_C} \rangle_{ACT}^{n=2}';
    'max\{ C \}_{ACT}^{n=2}'};
Single_TMA_Core_VariabelStrings = {'\langle C \rangle_{TMA}^{n=1}';...
    ' \langle \sigma_C \rangle_{TMA}^{n=1}'; ...
    '\langle I_C \rangle_{TMA}^{n=1}'; ...
    '\langle \sigma_{I_C} \rangle_{TMA}^{n=1}';
    'max\{ C \}_{TMA}^{n=1}'};
Store_SingleTMACores_ACT_Strings = {'SingleTMACores vs ACT_nanmeanClustersize';...
    'SingleTMACores vs ACT_VarianceClustersize';
    'SingleTMACores vs ACT_nanmeanClusterintensity';
    'SingleTMACores vs ACT_VarianceClusterintensity';
    'SingleTMACores vs ACT_MaxClustersize'};

%% plot TMA and ACT (MM2 --> green, MM3 --> red)


VariabelLabel = 8;

% figure('Renderer', 'painters', 'Position', [50 50 800 700])

PlotCIandTMAvsACT(( Store_ACT_MM3), (Store_TMA_MM3), NewNames, VariabelLabel, [1 0 0])
xlabel('$$ \langle \overline{\overline{I}}_p  \rangle_{ACT}^{n=2} $$', 'interpreter', 'latex')
ylabel('$$ \langle \overline{\overline{I}}_p  \rangle_{TMA}^{n=96} $$', 'interpreter', 'latex')
legend('CI', 'interpreter', 'latex')
% title("$$ X = \overline{\overline{I}}_p $$", 'Interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;

ax.XTick = [0.1:0.05:0.2];
ax.YTick = [0.1:0.05:0.2];
axis('equal')
% hold on
% line([0.05 0.2], [0.05 0.2], 'LineWidth', 2, 'Color', 'blue')


%%


VariabelLabel = 4;

% figure('Renderer', 'painters', 'Position', [50 50 800 700])

PlotCIandTMAvsACT(( Store_ACT_MM3), (Store_TMA_MM3), NewNames, VariabelLabel, [1 0 0])
xlabel('$$ \langle \sigma_{I_p}  \rangle_{ACT}^{n=2} $$', 'interpreter', 'latex')
ylabel('$$ \langle \sigma_{I_p}  \rangle_{TMA}^{n=96} $$', 'interpreter', 'latex')
legend('CI', 'interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;
ax.XTick = [0:0.03:0.1];
ax.YTick = [0:0.03:0.1];
axis('equal')
% hold on
% line([0.0 0.1], [0.0 0.1], 'LineWidth', 2, 'Color', 'blue')


%% compare ACT1 vs ACT2
Mum_per_Pixel = 0.3277;
VariabelLabel = 6;
PlotCIandTMAvsACT_withoutCI(Store_ACT, Store_ACT_MM2, Names, VariabelLabel)
xlabel('  \langle X \rangle_{ACT_{MM#3}}^{n=2}')
ylabel(' \langle X \rangle_{ACT_{MM#2}}^{n=2}')
title(strcat("X = ", X{VariabelLabel}))
axis('equal')
ax = gca;
ax.FontSize = 16;
% hold on
% line([20 40], [20 40], 'LineWidth', 2, 'Color', 'blue')
%% compare TMA 1 vs TMA2
X = {'nanmean Clustersize', 'nanstd Clustersize', 'nanmean Intensity', 'nanstd Intensity', 'max Clustersize', 'median Clustersize', 'mad Clustersize', 'median Intensity', 'mad Intensity'};
figure('Renderer', 'painters', 'Position', [50 50 800 700])

VariabelLabel = 8;
PlotCIandTMAvsACT_doubleCI(Store_TMA, Store_TMA_MM2, Names, VariabelLabel)
xlabel('\langle X \rangle_{TMA_{MM#3}}^{n=96}')
ylabel('\langle X \rangle_{TMA_{MM#2}}^{n=96}')
title(strcat("X = ", X{VariabelLabel}))
ax = gca;
ax.FontSize = 16;
% hold on
% line([0.0 170],[0.0 170], 'LineWidth', 2, 'Color', 'blue')


%% plot single tma cores against act nanmean (n=2)
VariabelLabel = 5;
TMAvsACTPlot(Store_ACT, Store_TMA, VariabelLabel, ACT_VariabelStrings, Single_TMA_Core_VariabelStrings)
saveas(gcf,Store_SingleTMACores_ACT_Strings{VariabelLabel}, 'png')
saveas(gcf,Store_SingleTMACores_ACT_Strings{VariabelLabel}, 'fig')
close all


%% compare act-normalized TMA 1 vs TMA 2
X = {'mean Clustersize', 'std Clustersize', 'mean Intensity', 'std Intensity', 'max Clustersize', 'median Clustersize', 'mad Clustersize', 'median Intensity', 'mad Intensity'};

figure('Renderer', 'painters', 'Position', [50 50 800 700])
inds = 1:13;
inds([8 12 10]) = [];
Names = NewNames; Names([8 12 10]) = [];

VariabelLabel = 8;
PlotCIandNormalizedTMAs_doubleCI(Store_TMA_MM3(inds,:,:), Store_TMA_MM2(inds,:,:), Store_ACT_MM3(inds,:,:), Store_ACT_MM2(inds,:,:), Names, VariabelLabel, [1 0 0]);%Store_TMA_MM3(:,randperm(96,15),:), Store_TMA_MM2(:,randperm(96,15),:) , Names, VariabelLabel, [1 0 1])
xlabel('$$\langle \overline{\overline{I}}_p^{stand} \rangle_{TMA_{MM\#3}}^{n=96}$$', 'interpreter', 'latex')
ylabel('$$\langle \overline{\overline{I}}_p^{stand} \rangle_{TMA_{MM\#2}}^{n=96}$$', 'interpreter', 'latex')
% title(strcat("X = ", X{VariabelLabel}), 'interpreter', 'latex')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 40;
axis('equal')



%% Paper Fig: --> plot TMA As with more than just 2 ACT Cores
%% get additioanl ACT Data
%order: #2ER PR Ki67, #3ER PR Ki67
precalculatedColors = 0;
All_TMA_B_Dir = {'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#2_TMA_B\aGanze TMA Bs\ER_Ganze_TMABs', ...
    'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#2_TMA_B\aGanze TMA Bs\PR_Ganze_TMABs', ...
    'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#2_TMA_B\aGanze TMA Bs\Ki67_Ganze_TMABs', ...
    'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_B\ER_Ganzer_TMAB', ...
    'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_B\PR_Ganzer_TMAB', ...
    'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_B\Ki67_Ganzer_TMAB'};


All_StoreAdditionalACT = cell(1, length(All_TMA_B_Dir));

for i = 1:length(All_TMA_B_Dir)
    TMA_B_Dir = All_TMA_B_Dir{i};
    StoreAdditionalACT = getAdditionalACTCores(TMA_B_Dir,precalculatedColors );
    All_StoreAdditionalACT{i} = StoreAdditionalACT;
end


%% now load mm2 and mm3 data and add aditional acts

ER_PR_Ki67_Indices = [8 12 10];
%fill #2
AllACTs_MM2_cell = [];
ACT_count = 0;
AllACTs_MM2_double = [];

for i = 1:size(Store_ACT, 1) 
    if ismember(i, ER_PR_Ki67_Indices)
        ACT_count = ACT_count + 1;
        AdditionalACTStore = reshape(All_StoreAdditionalACT{ACT_count}, 36, 9);
    else % do not fill
        AdditionalACTStore = nan(36,9);
    end
    AllACTs_MM2_cell{i} = cat(1,squeeze(Store_ACT_MM2(i,:,:)), AdditionalACTStore);
    AllACTs_MM2_double = cat(3, AllACTs_MM2_double, AllACTs_MM2_cell{i});
end
%get in right order
AllACTs_MM2_double = permute(AllACTs_MM2_double, [3 1 2]);
%fill #3
AllACTs_MM3_cell = [];
ACT_count = 0;
AllACTs_MM3_double = [];

for i = 1:size(Store_ACT, 1) 
    if ismember(i, ER_PR_Ki67_Indices)
        ACT_count = ACT_count + 1;
        AdditionalACTStore = reshape(All_StoreAdditionalACT{ACT_count}, 36, 9);
    else % do not fill
        AdditionalACTStore = nan(36,9);
    end
    AllACTs_MM3_cell{i} = cat(1,squeeze(Store_ACT_MM3(i,:,:)), AdditionalACTStore);
    AllACTs_MM3_double = cat(3, AllACTs_MM3_double, AllACTs_MM3_cell{i});
end
AllACTs_MM3_double = permute(AllACTs_MM3_double, [3 1 2]);




%%
VariabelLabel = 8;
PlotCIandTMAvsACT(( AllACTs_MM3_double(:,:,:)), (Store_TMA_MM3), NewNames, VariabelLabel, [0 1 1])
xlabel('\langle X \rangle_{ACT}^{n=2}')
ylabel('\langle X \rangle_{TMA}^{n=96}')
legend('Confidence Interval')
% title(strcat("X = ", X{VariabelLabel}))
ax = gca;
ax.FontSize = 16;
axis('equal')

hold on
line([0.05 0.2], [0.05 0.2], 'LineWidth', 2, 'Color', 'blue')



%% comapre act and tma histogram
PrintIt =0;
CompareACT_TMA_Histograms(2,Names,PrintIt , Store_ACT_Histograms, Store_TMA_Histograms , 'ACT_{MM#3}', 'TMA_{MM#3}')



%% plot function
function PlotCIandTMAvsACT(Store_ACT, Store_TMA, Names, VariabelLabel, Color)

Mum_per_Pixel = 0.3277;
hold on
for i = 1:size(Store_ACT, 1)
    format long
    CI = Calculate_CI_For_UnknownDistribution(Store_TMA(i,:,VariabelLabel), 0.05);
    line([nanmean(Store_ACT(i,:,VariabelLabel), 2) nanmean(Store_ACT(i,:,VariabelLabel), 2)], [CI(1) CI(2)], 'LineWidth', 3, 'Color', 'black')
    hold on
end
hold on
for i = 1:size(Store_ACT, 1)
    hold on
    scatter(nanmean(Store_ACT(i,:,VariabelLabel), 2), nanmean(Store_TMA(i,:,VariabelLabel)), 400 ,'filled', 'MarkerFaceColor',Color)
    text(nanmean(Store_ACT(i,:,VariabelLabel), 2), nanmean(Store_TMA(i,:,VariabelLabel)), Names{i}, 'FontSize', 20)
end

end

%% plot function with double CI
function PlotCIandTMAvsACT_doubleCI(Store_ACT, Store_TMA, Names, VariabelLabel)

hold on
for i = 1:13
    format long
    %vertical CI
    CI = Calculate_CI_For_UnknownDistribution(Store_TMA(i,:,VariabelLabel), 0.05);
    line([nanmean(Store_ACT(i,:,VariabelLabel), 2) nanmean(Store_ACT(i,:,VariabelLabel), 2)], [CI(1) CI(2)], 'LineWidth', 2, 'Color', 'black')
    hold on
    %horizontal CI
    CI = Calculate_CI_For_UnknownDistribution(Store_ACT(i,:,VariabelLabel), 0.05);
    line( [CI(1) CI(2)], [nanmean(Store_TMA(i,:,VariabelLabel), 2) nanmean(Store_TMA(i,:,VariabelLabel), 2)], 'LineWidth', 2, 'Color', 'black')
    
end
hold on
for i = 1:13
    hold on
    scatter(nanmean(Store_ACT(i,:,VariabelLabel), 2), nanmean(Store_TMA(i,:,VariabelLabel)), 50 ,'filled', 'MarkerFaceColor',[0 0 0])
    text(nanmean(Store_ACT(i,:,VariabelLabel), 2), nanmean(Store_TMA(i,:,VariabelLabel)), Names{i})
end

end

%% plot-function without CI


function PlotCIandTMAvsACT_withoutCI(Store_ACT, Store_TMA, Names, VariabelLabel)
figure,
for i = 1:13
    hold on
    scatter(nanmean(Store_ACT(i,:,VariabelLabel), 2), nanmean(Store_TMA(i,:,VariabelLabel)), 50 ,'filled', 'MarkerFaceColor',[0 0 0])
    text(nanmean(Store_ACT(i,:,VariabelLabel), 2), nanmean(Store_TMA(i,:,VariabelLabel)), Names{i})
end

end


%% plot single tma cores against act

function TMAvsACTPlot(Store_ACT, Store_TMA, VariabelLabel, ACT_VariabelStrings, Single_TMA_Core_VariabelStrings)
figure,
for i = 1:13
    for j = 1:size(Store_TMA,2)
        hold on
        scatter(nanmean(Store_ACT(i,:,VariabelLabel), 2), Store_TMA(i,j,VariabelLabel), 50 ,'filled', 'MarkerFaceColor',[0 0 0])
    end
end

xlabel(ACT_VariabelStrings{VariabelLabel})
ylabel(Single_TMA_Core_VariabelStrings{VariabelLabel})
ax = gca;
ax.FontSize = 16;
end

%% plot act_normalized tmas against each other
function PlotCIandNormalizedTMAs_doubleCI(Store_TMA1, Store_TMA2, Store_ACT1, Store_ACT2, Names, VariabelLabel, Color)
hold on
for i = 1:size(Store_TMA1, 1)
    format long
    Rescaled_TMA1 = (Store_TMA1(i,:,VariabelLabel)-nanmean(Store_ACT1(i,:,VariabelLabel)))./nanmean(Store_ACT1(i,:,VariabelLabel+1));
    Rescaled_TMA2 = (Store_TMA2(i,:,VariabelLabel)-nanmean(Store_ACT2(i,:,VariabelLabel)))./nanmean(Store_ACT2(i,:,VariabelLabel+1));
    %vertical CI
    CI = Calculate_CI_For_UnknownDistribution(Rescaled_TMA1, 0.05);
    line([nanmean(Rescaled_TMA2, 2) nanmean(Rescaled_TMA2, 2)], [CI(1) CI(2)], 'LineWidth', 2, 'Color', Color)
    hold on
    %horizontal CI
    CI = Calculate_CI_For_UnknownDistribution(Rescaled_TMA2, 0.05);
    line( [CI(1) CI(2)], [nanmean(Rescaled_TMA1, 2) nanmean(Rescaled_TMA1, 2)], 'LineWidth', 2, 'Color', Color)
    
end
hold on
for i = 1:size(Store_TMA1, 1)
    Rescaled_TMA1 = (Store_TMA1(i,:,VariabelLabel)-nanmean(Store_ACT1(i,:,VariabelLabel)))./nanmean(Store_ACT1(i,:,VariabelLabel+1));
    Rescaled_TMA2 = (Store_TMA2(i,:,VariabelLabel)-nanmean(Store_ACT2(i,:,VariabelLabel)))./nanmean(Store_ACT2(i,:,VariabelLabel+1));
    hold on
    scatter(nanmean(Rescaled_TMA2, 2), nanmean(Rescaled_TMA1, 2), 150 ,'filled', 'MarkerFaceColor',Color)
    text(nanmean(Rescaled_TMA2, 2), nanmean(Rescaled_TMA1, 2), Names{i}, 'FontSize', 20)
end

end

%% plot histograms

function CompareACT_TMA_Histograms( ThreshLevel, Names, PrintIt,Store_ACT_Histograms, Store_TMA_Histograms, legend1, legend2)

for StainNumber = 1:13
    ACT = Store_ACT_Histograms(StainNumber, :, ThreshLevel);
    TMA = Store_TMA_Histograms(StainNumber, :, ThreshLevel);
    edges = [0:0.02:1];
    if StainNumber < 13
        subplot(7,2,StainNumber)
        bar(edges(1:end-1),nanmean(vertcat(ACT{:})), 'FaceColor', [1 0 0], 'FaceAlpha', 0.5)
        hold on,
        bar(edges(1:end-1),nanmean(vertcat(TMA{:})), 'FaceColor', [0 1 0], 'FaceAlpha', 0.5)
        title(Names{StainNumber})
        ax = gca;
        ax.FontSize = 13;
    else
        subplot(7,2,[13,14])
        bar(edges(1:end-1),nanmean(vertcat(ACT{:})), 'FaceColor', [1 0 0], 'FaceAlpha', 0.5)
        hold on,
        bar(edges(1:end-1),nanmean(vertcat(TMA{:})), 'FaceColor', [0 1 0], 'FaceAlpha', 0.5)
        legend('ACT', 'TMA')
        title(Names{StainNumber})
        xlabel('Intensity')
        ylabel('Probability')
        legend(legend1, legend2)
        ax = gca;
        ax.FontSize = 13;
    end
    
end
if PrintIt == 1
    print(fig, strcat("C:\Users\PabiG\Desktop\ACT Paper Auswertung\Punkt4\Plots\", Names{StainNumber}), '-dpng')
    close all
end
end

%% initialize Punkt 2 Data
function [Store_ACT_TMAB, Store_ACT_Histograms_TMAB] =  InitializeTMAB(Colors, Grundmasse, TresholdToThresh, ThresholdFactor) %Grundmasse is either 2 or 3
if Grundmasse == 2
    All_TMA_Dir = 'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#2_TMA_B';
    
elseif Grundmasse == 3
    if Colors == 6
        All_TMA_Dir = 'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_B\Deconvolution_ColorsFromSix';
    elseif Colors == 1
        All_TMA_Dir =   'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_B\Deconvolution_SingleCoreColors';
    end
    
else
    disp('FALSCHE GRUNDMASSE ANGEGEBEN')
end
AllFolders = dir(All_TMA_Dir);
NumberOfSubFolder = 1;
ACT_Cores_To_Save = 6;
Store_ACT_TMAB = NaN(length(AllFolders)-2, ACT_Cores_To_Save, 9);
Store_ACT_Histograms_TMAB = cell(length(AllFolders)-2, ACT_Cores_To_Save, 3);
%iterate all TMAs
for TMA_Ind = 3:length(AllFolders)
    
    ACTSubFolder = strcat(All_TMA_Dir,"\" ,AllFolders(TMA_Ind).name);
    AllSubFolders = dir(ACTSubFolder);
    All_Images = strcat(ACTSubFolder,"\",AllSubFolders(2+NumberOfSubFolder).name ,"\*.tif");
    Single_TMA_Pics = dir(All_Images);
    %     mkdir(strcat(All_TMA_Dir, "\", AllFolders(TMA_Ind).name, "_Thresh"))
    tic
    
    
    for Pic_Ind = 1:length(Single_TMA_Pics)
        Current_Pic = ((imread(strcat(Single_TMA_Pics(Pic_Ind).folder,"\",Single_TMA_Pics(Pic_Ind).name))));Current_Pic = double(Current_Pic)./double(max(Current_Pic(:)));
        Current_Pic = nanmean(Current_Pic, 3);
        %threshold threshing inly valued pixels
        T = graythresh(Current_Pic(Current_Pic > TresholdToThresh)); %default 0.05
        %store normlized histograms
        edges = [0:0.02:1];
        [h1,~] = histcounts(Current_Pic(Current_Pic > ThresholdFactor*T),edges ,'Normalization', 'probability');
        Store_ACT_Histograms_TMAB{TMA_Ind-2, Pic_Ind, 3} = Current_Pic(Current_Pic > ThresholdFactor*T);
        Store_ACT_Histograms_TMAB{TMA_Ind-2, Pic_Ind, 1} = h1;
        %set background to zero
        Current_Pic(Current_Pic < ThresholdFactor*T ) = 0;
        CC = bwconncomp(Current_Pic);
        PXLList = CC.PixelIdxList;
        %delete clusters smaller than ~0.5 mum^2
        PXL_TO_DELETE = PXLList(cellfun(@length, PXLList) < 2);
        PXL_TO_DELETE = cat(1, PXL_TO_DELETE{:});
        Current_Pic(PXL_TO_DELETE) = 0;
        CC = bwconncomp(Current_Pic);
        ClusterSizes = cellfun(@length, CC.PixelIdxList);
        Store_ACT_Histograms_TMAB{TMA_Ind-2, Pic_Ind, 2} = ClusterSizes;
        Intensities = cellfun(@(Index)  CalcMeanClusterIntensity(Current_Pic, Index), CC.PixelIdxList);
        %store statistical properties
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 1) = nanmean(ClusterSizes);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 2) = nanstd(ClusterSizes);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 3) = nanmean(Intensities);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 4) = nanstd(Intensities);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 5) = max(ClusterSizes,[],'omitnan');
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 6) = nanmedian(ClusterSizes);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 7) = mad(ClusterSizes);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 8) = nanmedian(Intensities);
        Store_ACT_TMAB(TMA_Ind-2, Pic_Ind, 9) = mad(Intensities);
    end
    toc
end
end

%% function plot every act core in mean and std space

function  Point2_Plot_twoVariabels(ACT, Names, VariabelLabel1, VariabelLabel2, LogScale)
X = {'mean Clustersize', 'std Clustersize', 'mean Intensity', 'std Intensity', 'max Clustersize', 'median Clustersize', 'mad Clustersize', 'median Intensity', 'mad Intensity'};
Color{1} = [0 0 0]; Color{2} = [1 0 0];Color{3} = [0 1 0];Color{4} = [0 0 1];Color{5} = [1 128/255 0];Color{6} = [0 128/255 1]; Color{7} = [128/255 128/255 1];
for StainInd = 1:max(size(ACT(:,1,1)))
    for PicInd =  1:max(size(ACT(1,:,1)))
        hold on
        scatter(ACT(StainInd, PicInd, VariabelLabel1), ACT(StainInd, PicInd, VariabelLabel2), 'MarkerFaceColor', Color{StainInd}, 'MarkerEdgeColor', Color{StainInd}, 'LineWidth', 3)
        text(ACT(StainInd, PicInd, VariabelLabel1), ACT(StainInd, PicInd, VariabelLabel2), Names{StainInd})
    end
end
if LogScale == 1
    xlabel(strcat("log[ ", X{VariabelLabel1}, " ]"))
    ylabel(strcat("log[ ", X{VariabelLabel2}, " ]"))
else
    xlabel(X{VariabelLabel1})
    ylabel(X{VariabelLabel2})
end
ax = gca;
ax.FontSize = 16;
end

% calculate spatial correlation
function RadiallyAveragedCorrelation = calcSpatialCorrelation(ThreshedPic, Radius, NumberAngles)
addpath(genpath('C:\Users\PabiG\Documents\MATLAB'))
addpath(genpath('C:\Users\PabiG\Desktop\MasterThesis\MATLAB'))

result = Autocorr2D_SG(ThreshedPic);
[C_H,I] = max(result(:));
[y,x] = ind2sub(size(result), I);
UnwIm = Unwrapper(result, x, y, Radius, NumberAngles);
RadiallyAveragedCorrelation = UnwIm;% (mean(UnwIm, 2))./max(max(mean(UnwIm, 2)));

end

% estimate cluster radius
function EffectiveLengthScale = calcEffectiveLength(Schnipsel, Goto, MicronsPerPixel_TMAB)

autoi=Autocorr2D_SG(Schnipsel);
Autocorr_Schnipsel= UnwrapAutocorr_JL( autoi,1:Goto);
P = polyfit(Autocorr_Schnipsel(:,1),Autocorr_Schnipsel(:,2),1);
EffectiveLengthScale=-Autocorr_Schnipsel(1,2)/P(1)*MicronsPerPixel_TMAB;

end


%get TMA A Data

function [Store_ACT, Store_TMA] = calculate_TMA_A_Data(All_TMA_Dir, ColorSeparationPrecalculated)
addpath(genpath('C:\Users\PabiG\Documents\MATLAB'))
addpath(genpath('C:\Users\PabiG\Desktop\MasterThesis\MATLAB'))
addpath(genpath(('C:\Users\PabiG\Desktop\Bandscheiben - Biochemie\Code Jürgen')))
addpath('C:\Users\PabiG\Desktop\Improvement Algo\StainEstimation')
addpath('C:\Users\PabiG\Desktop\ACT Paper Auswertung\MATLAB')
addpath('C:\Users\PabiG\Desktop\Steffen Paper\KernShapes')

AllFolders = dir(All_TMA_Dir);
Store_TMA = NaN(length(AllFolders)-2, 96, 9);
Store_ACT = NaN(length(AllFolders)-2, 2, 9);

Store_TMA_Histograms = cell(length(AllFolders)-2, 96, 3);
Store_ACT_Histograms = cell(length(AllFolders)-2, 2, 3);
%iterate all TMAs
for TMA_Ind = 3:length(AllFolders)
    Single_TMA_Dir = strcat(All_TMA_Dir,"\" ,AllFolders(TMA_Ind).name);
    All_Images = strcat(Single_TMA_Dir, "\*.tif");
    Single_TMA_Pics = dir(All_Images);
    %     mkdir(strcat(All_TMA_Dir, "\", AllFolders(TMA_Ind).name, "_Thresh"))
    tic
    for Pic_Ind = 1:length(Single_TMA_Pics)
        if ColorSeparationPrecalculated == true
            Current_Pic = ((imread(strcat(Single_TMA_Dir,"\",Single_TMA_Pics(Pic_Ind).name))));Current_Pic = double(Current_Pic)./255;
            
            %for Punkt1 TMA_A S
            Current_Pic = nanmean(Current_Pic, 3);
            Current_Pic = Current_Pic./max(Current_Pic(:));
            Current_Pic = adpmedian(Current_Pic, 2*floor(10/2)+1);
            Current_Pic = Current_Pic-imgaussfilt(Current_Pic, 100);
            Current_Pic = Current_Pic./max(Current_Pic(:));
%             figure, imshow(Current_Pic)
            
            
        elseif ColorSeparationPrecalculated == false
            Current_Pic = ((imread(strcat(Single_TMA_Dir,"\",Single_TMA_Pics(Pic_Ind).name))));
            disp(num2str(TMA_Ind))
            disp(num2str(Pic_Ind))
            [ch1, ch2, ~,~,~] = ColorDeconvolution.calcColorTransform(Current_Pic, 0.99, 1, false, 0.05);
            Current_Pic = ch1./max(ch1(:));
        else
            disp('Color Deconvolution Precalculated or not???')
        end
        
        %calculate nanmean clustersize, variance clustersize, nanmean
        %clusterintensity, variance clusterintensity, Ratio
        %Kluster/background
        %         ~strcmp(Single_TMA_Pics(Pic_Ind).name,"A-8.tif") || ~strcmp(Single_TMA_Pics(Pic_Ind).name,"N-1.tif")
        if  Pic_Ind > 2
            
            T = graythresh(Current_Pic(Current_Pic > 0.05));
            
            %             PicToPrint = Current_Pic;
            %             PicToPrint(PicToPrint < 0.4*T) = 0;
            %             figure('Visible', 'off')
            %             imshow(PicToPrint)
            %             print(strcat(All_TMA_Dir, "\", AllFolders(TMA_Ind).name, "_Thresh", "\", Single_TMA_Pics(Pic_Ind).name(1:end-4)), '-dpng')
            %             close all
            
            intensity_edges = [0:0.02:1];
            [h1,~] = histcounts(Current_Pic(Current_Pic > 0.5*T),intensity_edges ,'Normalization', 'probability');
            [h2,~] = histcounts(Current_Pic(Current_Pic > 0.4*T),intensity_edges ,'Normalization', 'probability');
            [h3,~] = histcounts(Current_Pic(Current_Pic > 0.3*T),intensity_edges ,'Normalization', 'probability');
            
            Store_TMA_Histograms{TMA_Ind-2, Pic_Ind-2, 1} = h1;
            Store_TMA_Histograms{TMA_Ind-2, Pic_Ind-2, 2} = h2;
            Store_TMA_Histograms{TMA_Ind-2, Pic_Ind-2, 3} = h3;
            
            
            
            Current_Pic(Current_Pic < 0.4*T ) = 0;
            CC = bwconncomp(Current_Pic);
            PXLList = CC.PixelIdxList;
            %delete clusters smaller than ~0.5 mum^2
            PXL_TO_DELETE = PXLList(cellfun(@length, PXLList) < 4);
            PXL_TO_DELETE = cat(1, PXL_TO_DELETE{:});
            Current_Pic(PXL_TO_DELETE) = 0;
            
            CC = bwconncomp(Current_Pic);
            ClusterSizes = cellfun(@length, CC.PixelIdxList);
            Intensities = cellfun(@(Index)  CalcMeanClusterIntensity(Current_Pic, Index), CC.PixelIdxList);
            
            
            
            
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 1) = nanmean(ClusterSizes);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 2) = nanstd(ClusterSizes);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 3) = nanmean(Intensities);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 4) = nanstd(Intensities);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 5) = max(ClusterSizes,[],'omitnan');
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 6) = nanmedian(ClusterSizes);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 7) = mad(ClusterSizes);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 8) = nanmedian(Intensities);
            Store_TMA(TMA_Ind-2, Pic_Ind-2, 9) = mad(Intensities);
        else
            intensity_edges = [0:0.02:1];
            T = graythresh(Current_Pic(Current_Pic > 0.05));
            [h1,~] = histcounts(Current_Pic(Current_Pic > 0.5*T),intensity_edges ,'Normalization', 'probability');
            [h2,~] = histcounts(Current_Pic(Current_Pic > 0.4*T),intensity_edges ,'Normalization', 'probability');
            [h3,~] = histcounts(Current_Pic(Current_Pic > 0.3*T),intensity_edges ,'Normalization', 'probability');
            Store_ACT_Histograms{TMA_Ind-2, Pic_Ind, 1} = h1;
            Store_ACT_Histograms{TMA_Ind-2, Pic_Ind, 2} = h2;
            Store_ACT_Histograms{TMA_Ind-2, Pic_Ind, 3} = h3;
            
            
            Current_Pic(Current_Pic < 0.4*T) = 0;
            CC = bwconncomp(Current_Pic);
            PXLList = CC.PixelIdxList;
            %delete clusters smaller than ~0.5 mum^2
            PXL_TO_DELETE = PXLList(cellfun(@length, PXLList) < 4);
            PXL_TO_DELETE = cat(1, PXL_TO_DELETE{:});
            Current_Pic(PXL_TO_DELETE) = 0;
            CC = bwconncomp(Current_Pic);
            ClusterSizes = cellfun(@length, CC.PixelIdxList);
            Intensities = cellfun(@(Index)  CalcMeanClusterIntensity(Current_Pic, Index), CC.PixelIdxList);
            
            
            
            Store_ACT(TMA_Ind-2, Pic_Ind, 1) = nanmean(ClusterSizes);
            Store_ACT(TMA_Ind-2, Pic_Ind, 2) = nanstd(ClusterSizes);
            Store_ACT(TMA_Ind-2, Pic_Ind, 3) = nanmean(Intensities);
            Store_ACT(TMA_Ind-2, Pic_Ind, 4) = nanstd(Intensities);
            Store_ACT(TMA_Ind-2, Pic_Ind, 5) = max(ClusterSizes, [], 'omitnan');
            Store_ACT(TMA_Ind-2, Pic_Ind, 6) = nanmedian(ClusterSizes);
            Store_ACT(TMA_Ind-2, Pic_Ind, 7) = mad(ClusterSizes);
            Store_ACT(TMA_Ind-2, Pic_Ind, 8) = nanmedian(Intensities);
            Store_ACT(TMA_Ind-2, Pic_Ind, 9) = mad(Intensities);
            
        end
    end
    toc
end
toc
end

%% get additional ACT Data for TMA A s

function StoreAdditionalACT = getAdditionalACTCores(TMA_B_Dir, precalculatedColors)

addpath(genpath('C:\Users\PabiG\Documents\MATLAB'))
addpath(genpath('C:\Users\PabiG\Desktop\MasterThesis\MATLAB'))
addpath(genpath(('C:\Users\PabiG\Desktop\Bandscheiben - Biochemie\Code Jürgen')))
addpath('C:\Users\PabiG\Desktop\Improvement Algo\StainEstimation')
addpath('C:\Users\PabiG\Desktop\ACT Paper Auswertung\MATLAB')
% Goto TMA Bs
% TMA_B_Dir = 'C:\Users\PabiG\Desktop\ACT Paper Auswertung\MM#3_TMA_BColorDec';
All_Images_3 = strcat(TMA_B_Dir, "\*.tif");
All_Images_3 = dir(All_Images_3);

%constant variables
MicronsPerPixel_TMAB  = 0.2495.*2; % times 2 cause of compression for whole TMA Bs


% Measure
StoreAdditionalACT = NaN(length(All_Images_3), 6, 9);

% figure
for TMAInd = 1%:length(All_Images_3)
    
    
    im_uint8 = ((imread(strcat(TMA_B_Dir,"\",All_Images_3(TMAInd).name))));
    if ~precalculatedColors %has just 1 channel
        im_uint8 = im_uint8(:,:,1:3);
    end
    im_double = double(im_uint8)./max(max(double(im_uint8)));
    
    
    Cores = imdilate(~(mean(im_double, 3) > 0.94), strel('disk', round(30*1/MicronsPerPixel_TMAB)));
    CC_Cores = bwconncomp(Cores);
    [Size, CoreIndex] = sort(cellfun(@length, CC_Cores.PixelIdxList), 'descend');
    try
    Size = Size(1:6);
    CoreIndex = CoreIndex(1:6);
    catch
    end
    
    TooSmall = Size.*MicronsPerPixel_TMAB.^2 < 300^2*pi;
    Size(TooSmall) = [];
    CoreIndex(TooSmall) = [];
    Cores = zeros(size(Cores));
    Cores(cat(1, CC_Cores.PixelIdxList{CoreIndex})) = 1;
%      figure, imshow(Cores)
    
     
     
     if ~precalculatedColors
         if isempty(gcp('nocreate'))
            parpool(8) 
         end
         
%          [Channel_1, Channel_2, Residuum, st1, st2] = useOptimizedColorDeconvolution(im_uint8, 0, 0.2, 3, 60);
          [Channel_1, Channel_2, Residuum, st1, st2] = ColorDeconvolution.calcColorTransform(im_uint8, 1, 1, 0, 0.05);
         StainPositiveChannel = Channel_1;
     else
         StainPositiveChannel = im_double;
     end
     StainPositiveChannel(~logical(Cores)) = 0;
     
     
    
   %like for TMA As
    TresholdToThresh = 0.05;
    ThresholdFactor = 0.4;
    
    
    ACT_Cores_To_Save = length(CoreIndex);
    
    %iterate all TMAs
    for CoreInd = 1:ACT_Cores_To_Save
        Current_Pic = zeros(size(StainPositiveChannel));
        Current_Pic(CC_Cores.PixelIdxList{CoreIndex(CoreInd)}) = StainPositiveChannel(CC_Cores.PixelIdxList{CoreIndex(CoreInd)});
        Current_Pic = (Current_Pic)./(max(Current_Pic(:)));
        
        
%         figure, imshow(Current_Pic)
        %% jürgen enhace for TMA A standardization
        Current_Pic = Current_Pic./max(Current_Pic(:));
        Current_Pic = adpmedian(Current_Pic, 2*floor(MicronsPerPixel_TMAB^(-1)*3/2)+1);
        Current_Pic = Current_Pic-imgaussfilt(Current_Pic, round(MicronsPerPixel_TMAB^(-1)*30)); %gaussian filter 30 mum
        Current_Pic = Current_Pic./max(Current_Pic(:));
        
        %threshold threshing inly valued pixels
        T = graythresh(Current_Pic(Current_Pic > TresholdToThresh)); %default 0.05
        %store normlized histograms
        Current_Pic(Current_Pic < ThresholdFactor*T)  = 0;
        CC = bwconncomp(Current_Pic);
        PXLList = CC.PixelIdxList;
        
        
        %% measure like TMA A
        PXL_TO_DELETE = PXLList(cellfun(@length, PXLList) < 4);
        PXL_TO_DELETE = cat(1, PXL_TO_DELETE{:});
        Current_Pic(PXL_TO_DELETE) = 0;
        CC = bwconncomp(Current_Pic);
        ClusterSizes = cellfun(@length, CC.PixelIdxList);
        Intensities = cellfun(@(Index)  CalcMeanClusterIntensity(Current_Pic, Index), CC.PixelIdxList);
        StoreAdditionalACT(TMAInd, CoreInd, 1) = nanmean(ClusterSizes);
        StoreAdditionalACT(TMAInd, CoreInd, 2) = nanstd(ClusterSizes);
        StoreAdditionalACT(TMAInd, CoreInd, 3) = nanmean(Intensities);
        StoreAdditionalACT(TMAInd, CoreInd, 4) = nanstd(Intensities);
        StoreAdditionalACT(TMAInd, CoreInd, 5) = max(ClusterSizes, [], 'omitnan');
        StoreAdditionalACT(TMAInd, CoreInd, 6) = nanmedian(ClusterSizes);
        StoreAdditionalACT(TMAInd, CoreInd, 7) = mad(ClusterSizes);
        StoreAdditionalACT(TMAInd, CoreInd, 8) = nanmedian(Intensities);
        StoreAdditionalACT(TMAInd, CoreInd, 9) = mad(Intensities);
    end
end
end
