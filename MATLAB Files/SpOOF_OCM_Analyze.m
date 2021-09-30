%%
clear all;
OCMIntMap = OCMMap_viridis(256);
OCMPhaseMap = circshift(OCMMap_twilight(128),64);
OCMDivMap = SPoOFMap_colorcet('D5');
%%
fbase = 'Dish1_Glut_18MinLater_125.00usExp_4000fps_0';
% bbase = 'Back_200.00usExp_4000fps_24';
OutName1 = 'Processed_P1';
OutName2 = 'Processed_P2';

numF = 10900;

numx = 400;
numy = 400;

FID1 = fopen([pwd '\' fbase '\' OutName1 '.bin'],'rb');
fseek(FID1,numx*numy*2*1200*4,'bof');
Out_P1 = fread(FID1,numx*numy*2,'float');
Out_P1 = reshape( Out_P1(1:2:end-1) + 1i.*Out_P1(2:2:end),[numx numy]);
fclose(FID1);

FID2 = fopen([pwd '\' fbase '\' OutName2 '.bin'],'rb');
fseek(FID2,numx*numy*2*1200*4,'bof');
Out_P2 = fread(FID2,numx*numy*2,'float');
Out_P2 = reshape( Out_P2(1:2:end-1) + 1i.*Out_P2(2:2:end),[numx numy]);
fclose(FID2);

%%
OCMPhaseMap = circshift( SPoOFMap_colorcet('C3'), 0);

RaiMap = [linspace(0,0.4,64)' linspace(0,0,64)' linspace(0,0.5,64)';...
    linspace(0.4,1,64)' linspace(0,0,64)' linspace(0.5,0.5,64)';...
    linspace(1,1,64)' linspace(0,1,64)' linspace(0.5,0.5,64)';...
    linspace(1,1,64)' linspace(1,1,64)' linspace(0.5,1,64)';...
    ];

Sigma = 0.015;
[X, Y] = meshgrid( linspace(-1,1,numx),linspace(-1,1,numy));
RHO = sqrt(X.^2  + Y.^2);
G2D = exp( - ( RHO.*RHO)/(2.*Sigma.^2));
G2D = G2D / max(G2D(:));
G2D = 1 - G2D;
Win = (RHO.^2 <=1);

y1 = 0.5;
y2 = -15;
x1 = 1.5;
x2 = -12;
xyc = -6.5;
ycen = 0;
xcen = 0;
off = -0.5;
[x,y] = meshgrid( linspace(-1,1,numx),linspace(-1,1,numy));
PhaseMask1 = y1.*y + x1.*x + y2.*(y-ycen).^2 + x2.*(x-xcen).^2 + xyc.*x.*y + off;

y1 = 0;
y2 = 15;
x1 = -3;
x2 = 11.5;
xyc = -6.5;
ycen = 0;
xcen = 0;
off = -2.5;
[x,y] = meshgrid( linspace(-1,1,numx),linspace(-1,1,numy));
PhaseMask2 = y1.*y + x1.*x + y2.*(y-ycen).^2 + x2.*(x-xcen).^2 + xyc.*x.*y + off;

figure(1);
clf;
set(gcf,'Color',[1 1 1]);
subplot(321);
hold all;
imagesc( abs( Out_P1'/3).^0.8 .* Win,'AlphaData',Win);
colormap(gca,OCMIntMap);
axis image;
% colorbar;
caxis([0 1]*1500)
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('|E_{P1}|','Color',[1 1 1]*0);
colorbar('Location','EastOutside','LineWidth',1.5);


subplot(322);
hold all;
cla;
imagesc( abs( Out_P2'/2).^0.8 .* Win,'AlphaData',Win);
colormap(gca,OCMIntMap);
caxis([0 1]*1500)
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('|E_{P2}|','Color',[1 1 1]*0);
axis image;
colorbar('Location','EastOutside','LineWidth',1.5);


subplot(323);
imagesc( Win.*-angle(exp(-1i*PhaseMask1).*Out_P1)','AlphaData', Win);
colormap(gca,OCMPhaseMap);
axis image;
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('\angleE_{P1}','Color',[1 1 1]*0);
caxis([-pi pi]);
c = colorbar('Location','EastOutside','LineWidth',1.5,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-?','-?/2','0','?/2','?'});

subplot(324);
imagesc( Win.* angle(exp(-1i*PhaseMask2).*Out_P2)','AlphaData', Win);
colormap(gca,OCMPhaseMap);
axis image;
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('\angleE_{P2}','Color',[1 1 1]*0);
caxis([-pi pi]);
c = colorbar('Location','EastOutside','LineWidth',1.5,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-?','-?/2','0','?/2','?'});


subplot(325);
imagesc( Win.*angle(exp(-1i*PhaseMask1).*Out_P1.*exp(-1i*PhaseMask2).*Out_P2)','AlphaData', Win);
colormap(gca,OCMDivMap);
axis image;
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('\angleE_{P1} - \angleE_{P2}','Color',[1 1 1]*0);
c = colorbar('Location','EastOutside','LineWidth',1.5,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-?','-?/2','0','?/2','?'});
caxis([-pi pi])


subplot(326);
imagesc( Win.*(angle(exp(-1i*PhaseMask1).*Out_P1.*exp(1i*PhaseMask2).*conj(Out_P2))'),'AlphaData', Win);
colormap(gca,OCMPhaseMap);
axis image;
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('\angleE_{P1} + \angleE_{P2}','Color',[1 1 1]*0);
caxis([-pi pi]);
c = colorbar('Location','EastOutside','LineWidth',1.5,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-?','-?/2','0','?/2','?'});
savefig(gcf,[pwd '\' fbase '\JustBeautifulImages_20MAY2021.fig']);


%%
DFWin = 4;
Out_P1_DF = fftshift(fft2(Out_P1)).*G2D;
Out_P1_DF(end/2-DFWin+1:end/2+DFWin,end/2-DFWin+1:end/2+DFWin) = 0;
Out_P1_DF = (ifft2(ifftshift(Out_P1_DF)));

Out_P2_DF = fftshift(fft2(Out_P2)).* G2D;
Out_P2_DF(end/2-DFWin+1:end/2+DFWin,end/2-DFWin+1:end/2+DFWin) = 0;
Out_P2_DF = (ifft2(ifftshift(Out_P2_DF)));

AbsWin = 4;

ElecLoc = [215 236];
Orig = [135 150];

[LocListX, LocListY] = meshgrid((1:5:150)+Orig(1),(1:5:150)+Orig(2));
LocList = [LocListX(:) LocListY(:)];

DistToElec = sqrt( (LocList(:,1) - ElecLoc(1)).^2 + (LocList(:,2) - ElecLoc(2)).^2);
[~,I] = sort(DistToElec);

XL = LocList(I,1)-AbsWin/2;
YL = LocList(I,2)-AbsWin/2;

OCMPhaseMap2 = circshift(SPoOFMap_colorcet('C3'),128);
ROIColors = SPoOFMap_colorcet('I1','N',length(XL)).^2;
if length(XL) > 256
   ColorETemp = ROIColors;
   ROIColors = repmat(ROIColors,[ceil(length(XL)/256) 1]);
   for ccccidx = 1:ceil(length(XL)/256)
        ROIColors(ccccidx:ceil(length(XL)/256):end,:) = ColorETemp;
   end
end
    
figure(2);
set(gcf,'Color',[1 1 1]);
clf;
subplot(221);
hold all;
imagesc( abs( Out_P1_DF'/3).^0.8.*Win,'AlphaData', Win*0.4);
colormap(gca,OCMIntMap);
axis image;
caxis([0.1 1]*1500)
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('|E_{P1}|','Color',[1 1 1]*0);
for lidx = 1:length(XL)
    plot([YL(lidx)+AbsWin YL(lidx)+AbsWin],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
    plot([YL(lidx) YL(lidx)],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
    plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx) XL(lidx)],'Color',ROIColors(lidx,:),'LineWidth',1.5);
    plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx)+AbsWin XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
%     text(YL(lidx)-8, XL(lidx)-8, num2str(lidx),'Color',ROIColors(lidx,:),'FontSize',12,'FontName','Segoe UI');
end
colorbar('Location','EastOutside','LineWidth',1.5);

subplot(222);
hold all;
cla;
imagesc( abs( Out_P2_DF'/2).^0.8.*Win,'AlphaData', Win);
colormap(gca,OCMIntMap);
caxis([0.1 1]*1500)
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('|E_{P2}|','Color',[1 1 -2]*0);
axis image;
% for lidx = 1:length(XL)
%     plot([YL(lidx)+AbsWin YL(lidx)+AbsWin],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
%     plot([YL(lidx) YL(lidx)],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
%     plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx) XL(lidx)],'Color',ROIColors(lidx,:),'LineWidth',1.5);
%     plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx)+AbsWin XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
%     text(YL(lidx)-5, XL(lidx)-5, num2str(lidx),'Color',ROIColors(lidx,:),'FontSize',12,'FontName','Segoe UI');
% end
colorbar('Location','EastOutside','LineWidth',1.5);

OCMDivMap = SPoOFMap_colorcet('D5');

subplot(223);
hold all; cla;
imagesc( -angle(exp(-1i*PhaseMask1).*Out_P1)','AlphaData', Win*0.4);
colormap(gca,OCMPhaseMap);
axis image;
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI','Color',[1 1 1]);
title('\angleE_{P1}','Color',[1 1 1]*0);
for lidx = 1:length(XL)
    plot([YL(lidx)+AbsWin YL(lidx)+AbsWin],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',2);
    plot([YL(lidx) YL(lidx)],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',2);
    plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx) XL(lidx)],'Color',ROIColors(lidx,:),'LineWidth',2);
    plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx)+AbsWin XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',2);
%     if rem(lidx,10) == 1
%         text(YL(lidx)-8, XL(lidx)+14, num2str(lidx),'Color',ROIColors(lidx,:),'FontSize',15,'FontName','Segoe UI');
%     end
end
c = colorbar('Location','EastOutside','LineWidth',1.5,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-?','-?/2','0','?/2','?'});

subplot(224);
imagesc( angle(exp(-1i*PhaseMask2).*Out_P2)','AlphaData', Win);
colormap(gca,OCMPhaseMap);
axis image;
set(gca,'XTick',[],'YTick',[],'FontSize',15,'YDir','reverse','XColor',[1 1 1],'YColor',[1 1 1],'FontName','Segoe UI');
title('\angleE_{P2}','Color',[1 1 1]*0);
c = colorbar('Location','EastOutside','LineWidth',1.5,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-?','-?/2','0','?/2','?'});

savefig(gcf,[pwd '\' fbase '\' fbase '_ROIPlots_20MAY2021.fig']);

%%
FID1 = fopen([pwd '\' fbase '\' OutName1 '.bin'],'rb');
FID2 = fopen([pwd '\' fbase '\' OutName2 '.bin'],'rb');

numFrames = numF;
ResP1 = zeros([AbsWin+1 AbsWin+1 length(XL) numFrames]);
ResP2 = zeros([AbsWin+1 AbsWin+1 length(XL) numFrames]);
wBar = waitbar(0,'Starting...');
tic;
for frameIdx = 1:numFrames
    Out_P1 = fread(FID1,numx*numy*2,'float');
    Out_P1 = reshape( Out_P1(1:2:end-1) + 1i.*Out_P1(2:2:end),[numx numy]);

    Out_P2 = fread(FID2,numx*numy*2,'float');
    Out_P2 = reshape( Out_P2(1:2:end-1) + 1i.*Out_P2(2:2:end),[numx numy]);

    Out_P1 = exp(-1i*PhaseMask1).*Out_P1;
    Out_P2 = exp(-1i*PhaseMask2).*Out_P2;
    for lidx = 1:length(XL)
        ResP1(:,:,lidx,frameIdx) = Out_P1(YL(lidx):YL(lidx)+AbsWin,XL(lidx):XL(lidx)+AbsWin);
        ResP2(:,:,lidx,frameIdx) = Out_P2(YL(lidx):YL(lidx)+AbsWin,XL(lidx):XL(lidx)+AbsWin);
    end
    waitbar(frameIdx/numFrames,wBar,['Loaded ' num2str(frameIdx) ' frames in ' num2str(toc,'%1.0f') ' s...']);
end
close(wBar);
fclose(FID2);
fclose(FID1);
%% Butterworth filtered P1 and P2 phase Plots

FileDetails.fbase = fbase;
FileDetails.fDEEET = 'Spontaneous'; 
FileDetails.FileExt = 'TodaysDate';
N = length(XL);

Stim = 0;
fs = 4000;

Filt.ShowStim = 0;
Filt.ShowEPhys = 1;
Filt.LocEPhys = 0.5;
Filt.EPhysCShift = -1000;
Filt.SaveFig = 0;

Filt.fop = 4;
Filt.PlotOption = 2;

Filt.Gap = 25;
Filt.windowSizeMean = 5; 
Filt.AbsWin = AbsWin;
Filt.PhaseOff = 1;
Filt.xOff = 0;
Filt.yOff = 0;

[T3_N, taxis, StringOpt, FiltOpt, CorrVals_Plots] = SPoOFOCM_PlotResponse(ResP1, ResP2, N, fs, numFrames, FileDetails, Filt, Stim);
xlim([1600 2400]-400);
set(gca,'YTick',[],'YColor',[1 1 1]);
plot([402 402],[30 50],'k','LineWidth',3);
text(398,40,'20 nm','Color','k','FontName','Segoe UI','FontSize',15,'HorizontalAlignment','right');
ylim([-20 1500])
% ylim([-20 (N+3).*Filt.Gap]);

%%
T3Chart = (T3_N');
T3Chart(2:end,:) = (T3_N(:,2:end)')- mean(T3_N(:,2:end),2,'omitnan')';
figure(4);
clf; set(gcf,'Color',[1 1 1]);
imagesc((taxis-0.4)*1000,0:N,T3Chart)
colormap(SPoOFMap_colorcet('D4'));
% caxis([1 3])
caxis([-10 10])
box on;
xlabel('Time (ms)');
title({['NE-4C cells | ' FileDetails.fDEEET];[StringOpt ' | ' FiltOpt]});
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal');
xlim([1 1700]);

%%

T3Chart_Norm = zeros(size(T3Chart));
for lidx = 2:N
    NoiseVal = std( T3Chart(lidx,2000:8000));
    T3Chart_Norm(lidx,:) = T3Chart(lidx,:)./NoiseVal ;
end

%%
FID1 = fopen([pwd '\' fbase '\' OutName1 '.bin'],'rb');


% v = VideoWriter('Dish1_18minLater_PlotOpt2_StandardDev.mp4','MPEG-4');
% v.Quality = 99;
% v.FrameRate = 10;
% open(v);

tStart = 550*4;
tEnd = 2250*4;
tSpace = 100;
tArray = tStart:tSpace:tEnd;
numTBins = length(tArray) - 1;

fseek(FID1,numx*numy*2*tStart*4,'bof');
OutP_Show_Prev = zeros([numx numy]);

figure(20); set(gcf,'Color',[1 1 1]);clf; 
set(gcf,'InvertHardCopy','off');

for tidx = 1:numTBins
    clf;
    Out_PShow = fread(FID1,numx*numy*2,'float');
    fseek(FID1,numx*numy*2*(tSpace-1)*4,'cof');
    Out_PShow = reshape( Out_PShow(1:2:end-1) + 1i.*Out_PShow(2:2:end),[numx numy]);   
    
    T3ToPlot = std( (T3Chart_Norm(2:end,tStart + (tidx-1)*tSpace:tStart + (tidx)*tSpace)).^2,...
        1,2,'omitnan');

    
    YLL = YL-min(YL);
    XLL = XL-min(XL);
    
%     subplot(3,1,2);
%     imagesc( reshape(T3ToPlot, [30 30])');
%     colormap(gca,SPoOFMap_colorcet('CBL1'));
%     axis image;
%     set(gca,'YTick',[],'XTick',[ ],'FontName','Segoe UI','FontSize',15,'YDir','reverse');
%     caxis([0 4]);
    
    T3ToPlot(T3ToPlot <= 2.65) = 0;
    T3ToPlot(T3ToPlot >= 4) = 4;
    T3ToPlot(isnan(T3ToPlot)) = 0;
    XLL(T3ToPlot == 0) = [];
    YLL(T3ToPlot == 0) = [];
    T3ToPlot(T3ToPlot == 0) = [];
    
    subplot(3,1,[1 2]); hold all;
    set(gca,'Color','k');
    imagesc( angle(exp(-1i*PhaseMask1(min(YL):max(YL),min(XL):max(XL)) + 1i).*...
        Out_PShow(min(YL):max(YL),min(XL):max(XL)))','AlphaData',0.8);
    scatter(YLL,XLL,(T3ToPlot(:).^2)*5,T3ToPlot(:),'s',...
        'MarkerEdgeColor','none','MarkerFaceColor','y');%,'MarkerFaceAlpha',T3ToPlot(:));   
    caxis([-pi pi]);
    colormap(gca,SPoOFMap_colorcet('C5'));
    axis image; box on;
    set(gca,'YTick',[],'XTick',[ ],'FontName','Segoe UI','FontSize',15,'YDir','reverse','LineWidth',1.5);
    text(5,5,[num2str(taxis((tidx-1)*tSpace+1)*1000) ' to ' ...
        num2str(taxis((tidx)*tSpace+1)*1000) ' ms'],'Color','w','FontName','Segoe UI','FontSize',15);
    text(32,109,'?_1','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
    text(82,83,'?_2','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
    text(111,101,'?_3','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
    text(59,34,'?_4','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
%     plot([115 135],[5 5],'LineWidth',6,'Color','w');
%     text(116,9,'10 µm','FontName','Segoe UI','FontSize',15,'Color','w');
    
    ylim([min(XL) max(XL)]-min(XL));
    xlim([min(YL) max(YL)]-min(YL));
    
    
        
%     for lidx = 2:50%N
%         plot( [taxis(tStart + (tidx-1)*tSpace+1) taxis(tStart + (tidx)*tSpace+1)]*1000,...
%             [T3ToPlot(lidx) T3ToPlot(lidx)]+lidx*10,'Linewidth',1,'Color',ROIColors(lidx,:));
%     end
%     set(gca,'Ytick',[],'FontName','Segoe UI','FontSize',15,'YDir','normal','LineWidth',1);
%     xlim([tStart tEnd]/4);

 
    subplot(3,1,3);
    cla;
    plot( taxis(1:1 + (tidx)*tSpace)*1000, T3Chart(1,tStart:tStart + (tidx)*tSpace)/5,'k','LineWidth',1);
    xlim(([tStart tEnd] - tStart)/4);
    box on;
    set(gca,'Ytick',[],'FontName','Segoe UI','FontSize',15,'YDir','normal','LineWidth',1);
    OutP_Show_Prev = Out_PShow;
    ylim([-6 8]);
    xlabel('Time (ms)');
    ylabel('EPhys recording (a.u.)');
    
    pause(0.1);
    
    print(gcf,[pwd '\' fbase '\PicsSeries\StandardDev_PlotOpt2_timeBin12p5ms_' num2str(tidx) '.jpg'],'-djpeg','-r600');
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end
% close(v);
fclose(FID1);

%%

tStart = 550*4;
tEnd = 2250*4;
tSpace = 100;
tArray = tStart:tSpace:tEnd;
numTBins = length(tArray) - 1;


figure(20); set(gcf,'Color',[1 1 1]);clf; 
set(gcf,'InvertHardCopy','off');
clf;
subplot(3,1,[1 2]); hold all;
set(gca,'Color','k');
imagesc( angle(exp(-1i*PhaseMask1(min(YL):max(YL),min(XL):max(XL)) + 1i).*...
    Out_PShow(min(YL):max(YL),min(XL):max(XL)))','AlphaData',0.6);
TimeColors = SPoOFMap_colorcet('I1', 'N',numTBins).^1.1;

for tidx = 1:numTBins
    
    T3ToPlot = std( (T3Chart_Norm(2:end,tStart + (tidx-1)*tSpace:tStart + (tidx)*tSpace)).^2,...
        1,2,'omitnan');

    
    YLL = YL-min(YL);
    XLL = XL-min(XL);
    
    T3ToPlot(T3ToPlot < 1) = 0;
    T3ToPlot(T3ToPlot >= 2) = 2;
    T3ToPlot(isnan(T3ToPlot)) = 0;
    
    XLL(T3ToPlot == 0) = [];
    YLL(T3ToPlot == 0) = [];
    T3ToPlot(T3ToPlot == 0) = [];
    
    subplot(3,1,[1 2]); hold all;
    scatter(YLL,XLL,(T3ToPlot(:).^2)*12,T3ToPlot(:),'s',...
        'MarkerEdgeColor','none','MarkerFaceColor',TimeColors(tidx,:),'MarkerFaceAlpha',0.4);   
    caxis([-pi pi]);
    colormap(gca,SPoOFMap_colorcet('C5'));
    axis image; box on;
    set(gca,'YTick',[],'XTick',[ ],'FontName','Segoe UI','FontSize',15,'YDir','reverse','LineWidth',1.5);
%     text(5,5,[num2str(taxis((tidx-1)*tSpace+1)*1000) ' to ' ...
%         num2str(taxis((tidx)*tSpace+1)*1000) ' ms'],'Color','w','FontName','Segoe UI','FontSize',15);
%     text(32,109,'?_1','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
%     text(82,83,'?_2','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
%     text(111,101,'?_3','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
%     text(59,34,'?_4','FontName','Segoe UI','FontSize',24,'Color',([232 74 39]/256).^0.5);
%     plot([115 135],[5 5],'LineWidth',6,'Color','w');
%     text(116,9,'10 µm','FontName','Segoe UI','FontSize',15,'Color','w');
    
    ylim([min(XL)-5 max(XL)+5]-min(XL));
    xlim([min(YL)-5 max(YL)+5]-min(YL));
    
    
        
%     for lidx = 2:50%N
%         plot( [taxis(tStart + (tidx-1)*tSpace+1) taxis(tStart + (tidx)*tSpace+1)]*1000,...
%             [T3ToPlot(lidx) T3ToPlot(lidx)]+lidx*10,'Linewidth',1,'Color',ROIColors(lidx,:));
%     end
%     set(gca,'Ytick',[],'FontName','Segoe UI','FontSize',15,'YDir','normal','LineWidth',1);
%     xlim([tStart tEnd]/4);

 
    subplot(3,1,3); hold all;
    plot( taxis((tidx-1)*tSpace + 1:1 + (tidx)*tSpace)*1000, ...
        T3Chart(1,tStart + (tidx-1)*tSpace:tStart + (tidx)*tSpace)/5,'Color',TimeColors(tidx,:).^2,'LineWidth',1);
    xlim(([tStart tEnd] - tStart)/4);
    box on;
    set(gca,'Ytick',[],'FontName','Segoe UI','FontSize',15,'YDir','normal','LineWidth',1);
    OutP_Show_Prev = Out_PShow;
    ylim([-6 8]);
    xlabel('Time (ms)');
    ylabel('EPhys recording (a.u.)');
    
    pause(0.01);
end



%%
figure(14);
Filt.Gap = 7;
% ColorPlotMap = SPoOFMap_colorcet('I1','N',N).^2;
ColorPlotMap = [0 0 0; ROIColors];
clf;
hold all;
set(gcf,'Color',[1 1 1]);
NoiseVal = std( T3Chart(1,2000:8000));
plot(taxis*1000 - 400,T3Chart(1,:)/5 + Filt.Gap,'Color',ColorPlotMap(1,:),...
        'LineWidth',1);
for lidx = 2:50%N+1
    NoiseVal = std( T3Chart(lidx,2000:8000));
    T3SigNorm = T3Chart(lidx,:)./NoiseVal ;
    plot(taxis*1000 - 410, T3SigNorm + Filt.Gap*lidx,'Color',ColorPlotMap(lidx,:),...
        'LineWidth',1);
    if rem(lidx-2,10) == 0
       text(1700,mean(T3SigNorm(2000:end-2000))+ Filt.Gap*lidx,...
           ['ROI ' num2str(lidx-1)],'FontSize',15,'FontName','Segoe UI');
    end
end
xlim([1 1700]);
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal','YTick',[],'YColor',[1 1 1]);
ylim([0 375]);
plot([1 1],[50 60],'Color','k','LineWidth',4);
text(-85,55,'10 nm','FontSize',15,'FontName','Segoe UI');
%%
%
[bMean,aMean] = butter(6,[25]/(fs/2),'low');
Sig1 = median( T3Chart(1,4800:8400),1,'omitnan')/7;
% Sig1 = Sig1 - filter(bMean, aMean, Sig1);
% Sig1 = (Sig1 - min(Sig1))./(max(Sig1) - min(Sig1));
Sig2 = mean( T3Chart(2:26,4800:8400),1,'omitnan');
% Sig2 = circshift( (Sig2 - min(Sig2))./(max(Sig2) - min(Sig2)),0);

[bStop,aStop] = butter(4,[115 125]/(fs/2),'stop');
[Acor, LagE] = xcorr( Sig1, Sig2, ...
    'coeff');
Acor = filter(bStop,aStop,Acor);

[~,I] = max((Acor( round(length(Acor)/2-400):round(length(Acor)/2+400))));
lagDiff = LagE(I + round(length(Acor)/2 - 400))/4;
disp([num2str(lagDiff,'%1.3f'),' ms']);
figure(6); clf; box on;
subplot(211);
plot(LagE/4,Acor)
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal');

subplot(212);
plot(linspace(-2000,2000,length(LagE)),log10( abs( fftshift(fft(Acor)))))
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal');

figure(5);clf; set(gcf,'Color',[1 1 1]);
subplot(2,1,1);hold all; box off;
plot(taxis(4800:8400)*1000 - 410, Sig2,'Color',[232 74 39]/256,'LineWidth',1.0); 
plot(taxis(4800:8400)*1000 - 400, Sig1 - 6,'k','LineWidth',1.0); 
% plot(taxis*1000, mean( T3Chart([87],:),1,'omitnan')-6,'Color',[0.2 0.4 1],'LineWidth',1.0); 
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal','YTick',[],'YColor',[1 1 1]);
xlim([1000 1700]);
ylim([-9 3]);

% CorrVals_Plots (isnan(CorrVals_Plots)) = 0;
subplot(2,3,4); 
scatter (CorrVals_Plots(:,1),CorrVals_Plots(:,2),25,1:size(CorrVals_Plots,1),'filled'); 
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal');
colormap(gca,ROIColors);
box on;

subplot(2,3,5);
plot(CorrVals_Plots(:,1),'LineWidth',1.0); 
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal');
CorrVals_Plots(isnan(CorrVals_Plots(:,1)),2) = 0;
CorrVals_Plots((CorrVals_Plots(:,1) < -0.5),2) = NaN;


subplot(2,3,6);
plot(1:N,CorrVals_Plots(:,2),'LineWidth',1.0); 
set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI','YDir','normal');
% ylim([-100 20])
% xlim([1000 1700]);
%% cross correlation on filtered birefringence

NTwin = 25;
LenWin = floor((numFrames - 899)/NTwin);
CrossCorrMat = zeros([N N NTwin]);
TimeDiff = zeros([N N NTwin]);
StdDisp = zeros([N N NTwin]);
DisThresh = 0.4;
tshort = taxis(1:LenWin);
for tidx = 2:NTwin
    for lidx1 = 1:N        
        for lidx2 = 1:(lidx1-1)
            Diss1 = squeeze(T3_N((tidx- 1)*LenWin +1:(tidx + 0)*LenWin,lidx1));
            Diss2 = squeeze(T3_N((tidx- 1)*LenWin +1:(tidx + 0)*LenWin,lidx2));
            
            [CCorr] = corrcoef( Diss1, Diss2);
            [Acor, LagE] = xcorr( Diss1, Diss2,50*4,'unbiased');
%             [Acor, LagE] = xcorr( Diss1, Diss2);
             [~,I] = max(abs(Acor));
            lagDiff = LagE(I);
    
            StdDisp(lidx1,lidx2, tidx) = (( std(Diss1) - DisThresh) > 0)*(( std(Diss2) - DisThresh) > 0);        
            CrossCorrMat(lidx1,lidx2,tidx) = CCorr(2);
            TimeDiff(lidx1,lidx2,tidx) = 1000*lagDiff/fs;
        end
    end
end
%
figure(10);
set(gcf,'Color',[1 1 1]);
clf;
for tidx = 2:NTwin
    subplot(4,(NTwin-1)/4,tidx-1);
%     subplot(2,(NTwin)/2,tidx-1);
%     imagesc(CrossCorrMat(:,:,tidx))%,'AlphaData',StdDisp(:,:,tidx)>0);
%     colormap(gca,SPoOFMap_colorcet('L17'));
%     caxis([0.35 0.65]);
    
    imagesc( TimeDiff(:,:,tidx),'AlphaData',(abs(CrossCorrMat(:,:,tidx))-0.2)*5);   
    colormap(gca,SPoOFMap_colorcet('D5'));
    caxis([-1 1]*10);

    if tidx > 25
        colorbar('Location','southoutside');
    end
    text(12,6,[num2str(taxis((tidx-1)*LenWin +1)*1000,'%1.0f') ...
        ' - ' num2str(taxis((tidx)*LenWin)*1000,'%1.0f') ' ms'],'Color','k','FontSize',15,'FontName','Segoe UI');
    
    axis image;
    set(gca,'FontSize',15,'Color',[1 1 1],'XColor',[1 1 1]-1,'YColor',[1 1 1]-1,'FontName','Segoe UI');
    set(gca,'LineWidth',1.5,'YDir','reverse','XAxisLocation','bottom');
    box off;
    pause(0.05);
end

% savefig(gcf,[pwd '\' FileDetails.fbase '\SumOfAngles_' FiltOpt '_' FileDetails.FileExt '_CorrCoeff_0p25To0p75.fig']);
% savefig(gcf,[pwd '\' FileDetails.fbase '\SumOfAngles_' FiltOpt '_' FileDetails.FileExt '_Lag_m2To2.fig']);
%% SHOW CONN MAP TIME DIFF


IMShoE = -Win.*(angle(Out_P1.*exp(2i))');
Lcol = [0.1 0.1 0.1];

figure(10);
clf; set(gcf,'Color',[1 1 1]);
CCRLim = 0.7;
TLagLim = 10;
TLagFact = 10;

TlagMap = SPoOFMap_colorcet('D5','N',TLagLim*2 + 1);
for tidx = 2:NTwin
    subplot(4,(NTwin-1)/4,tidx-1);
    
    hold all;
    imagesc( IMShoE,'AlphaData',0.8*Win);  
    axis image;
    xlim([0 400]);
    ylim([0 400]);
    box on;
    set(gca,'FontSize',15,'LineWidth',1.5,'Color',[1 1 1],'XTick',[],'YTick',[],'YDir','reverse');
    for lidx1 = 1:N
        for lidx2 = 1:(lidx1-1)
%             StdDisp(lidx1,lidx2,tidx)
            if( CrossCorrMat(lidx1,lidx2,tidx) >= CCRLim)% && StdDisp(lidx1,lidx2,tidx) > 0)
                TvalM = round((TimeDiff(lidx1,lidx2,tidx)*TLagFact + TLagLim));
                if TvalM < 1
                    TvalM = 1;
                elseif TvalM > TLagLim*2+1
                    TvalM = TLagLim*2+1;
                end
                plot( [YL(lidx1)+AbsWin/2 YL(lidx2)+AbsWin/2],[XL(lidx1)+AbsWin/2 XL(lidx2)+AbsWin/2],'Color',TlagMap(TvalM,:),'LineWidth',2);
            end
        end
    end
    for lidx = 1:length(XL)
    plot([YL(lidx)+AbsWin YL(lidx)+AbsWin],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
    plot([YL(lidx) YL(lidx)],[XL(lidx) XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
    plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx) XL(lidx)],'Color',ROIColors(lidx,:),'LineWidth',1.5);
    plot([YL(lidx) YL(lidx)+AbsWin],[XL(lidx)+AbsWin XL(lidx)+AbsWin],'Color',ROIColors(lidx,:),'LineWidth',1.5);
%     text(YL(lidx)+4, XL(lidx)-15, num2str(lidx),'Color',Lcol,'FontSize',10);
    end
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1]);
    colormap(gca,SPoOFMap_colorcet('C5'));
%     title(['Conn Map b/n t = ' num2str(taxis((tidx-1)*LenWin +1),'%1.2f') ...
%         ' s and ' num2str(taxis((tidx)*LenWin),'%1.2f') ' s']);
    text(200,-15,[num2str(taxis((tidx-1)*LenWin +1)*1000,'%1.0f') ...
      ' - ' num2str(taxis((tidx)*LenWin)*1000,'%1.0f') ' ms'],...
      'Color','k','FontSize',14,'FontName','Segoe UI','HorizontalAlignment','center');
  pause(0.1);
end