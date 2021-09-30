function [T3_N, taxis, StringOpt, FiltOpt, CorrVals] = SPoOFOCM_PlotResponse(ResP1, ResP2, N, fs, numFrames, FileDetails, Filt, Stim)
    figure(7);
    clf; hold all;
    set(gcf,'Color',[1 1 1]);
    BlueMap = SPoOFMap_colorcet('L13','N',N);
    RedMap = SPoOFMap_colorcet('L14','N',N);
    ColorEMap = SPoOFMap_colorcet('I1','N',N).^2;
    if N > 256
       ColorETemp = ColorEMap;
       ColorEMap = repmat(ColorEMap,[ceil(N/256) 1]);
       for ccccidx = 1:ceil(N/256)
            ColorEMap(ccccidx:ceil(N/256):end,:) = ColorETemp;
       end
    end
    CorrVals = zeros([N 2]);
%%  
    taxis = (0:(numFrames-1))/fs;
    
    if Filt.ShowStim == 1
        Stim.ONT = 15;
        Stim.OFFT = 25;
        Stim.NPulses1 = 50;
        Stim.NGap = 0;
        Stim.NPulses2 = 0;
        Stim.Delay = 450;
        Stim.fs = 4000;

        Stim.Stim = [ linspace(0,0,Delay/1000*fs) ...
            repmat([linspace(1,1,ONT/1000*fs) linspace(0,0,OFFT/1000*fs)],[1 NPulses1])...
            repmat([linspace(0,0,ONT/1000*fs) linspace(0,0,OFFT/1000*fs)],[1 NGap])...
            repmat([linspace(1,1,ONT/1000*fs) linspace(0,0,OFFT/1000*fs)],[1 NPulses2])]; 
    end
%%    

    if Filt.fop == 1
        fc = [50 70];
        [b,a] = butter(3,fc/(fs/2),'stop');
        FiltOpt = '60HzF';
    elseif Filt.fop == 2
        FiltOpt = '1000HzLPF';
        [b,a] = butter(6,[1000]/(fs/2),'low');
    elseif Filt.fop == 3
        FiltOpt = '70HzHPF';
        [b,a] = butter(6,[70]/(fs/2),'high');  
    elseif Filt.fop == 4
        FiltOpt = 'Gauss';
        windowSize = 8; 
        b = gausswin(windowSize,2);%(1/windowSize)*ones(1,windowSize);
        b = b/sum(b);
        a = 1;
    elseif Filt.fop == 5
        FiltOpt = 'Rect';
        windowSize = 10; 
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
    end
    
%     bMean = (1/Filt.windowSizeMean)*ones(1,Filt.windowSizeMean);
%     aMean = 1;
    [bMean,aMean] = butter(6,[Filt.windowSizeMean]/(fs/2),'low');
%         bMean = gausswin(Filt.windowSizeMean,Filt.windowSizeMean/2);%(1/windowSize)*ones(1,windowSize);
%         aMean = 1;

    T1Med = zeros([N numFrames]);
    T2Med = zeros([N numFrames]);
    for lidx = 1:N
        RP1 = circshift(circshift( ResP1, Filt.xOff,1), Filt.yOff,2);
        RP2 = circshift(circshift( ResP2, Filt.xOff,1), Filt.yOff,2);
        T1 = reshape( RP1(Filt.PhaseOff+1:end-Filt.PhaseOff,Filt.PhaseOff+1:end-Filt.PhaseOff,lidx,:),...
            [((Filt.AbsWin-2*Filt.PhaseOff)+1).^2 numFrames]);
        T2 = reshape( circshift(circshift( RP2(Filt.PhaseOff+1:end-Filt.PhaseOff,Filt.PhaseOff+1:end-Filt.PhaseOff,lidx,:), Filt.xOff,1), Filt.yOff,2),...
            [((Filt.AbsWin-2*Filt.PhaseOff)+1).^2 numFrames]);
        T1Med(lidx,:) = mean([T1],1);
        T2Med(lidx,:) = mean([T2],1);
    end
    %
    T1Med = mean(T1Med,1);
    T2Med = mean(T2Med,1);
    
    %%
    if Filt.ShowEPhys == 1   
        FIDE = fopen([pwd '\' FileDetails.fbase '\EPhysRec.bin'],'rb');
        EPhysRec = fread(FIDE,109000,'double');
        fclose(FIDE);

        Divfactor = 1;
        taxis2 = linspace(0,109000/40000,109000/Divfactor);
        EPhysRec = reshape( EPhysRec, [Divfactor 109000/Divfactor]);
        EPhysRec = EPhysRec(1,:);

        fstop = 400;
        fc = [fstop-10 fstop+10];
        [bE,aE] = butter(2,[1000]/((40000/Divfactor)/2),'low');
        EPhysFilt = filtfilt(bE,aE,EPhysRec);
        EPhysFilt = mean( reshape( EPhysFilt,[10 10900]),1);
        EPhysFilt = circshift((EPhysFilt(1:numFrames) - mean(EPhysFilt(1:numFrames))),Filt.EPhysCShift);%-972);%,180);

        
%         EPhysFilt = filtfilt(bStop,aStop,EPhysFilt);
%         plot(taxis2*1000,EPhysRec*Filt.Gap,'LineWidth',0.5,'Color',[0.25 0.25 0.25]); 
        plot(taxis*1000,(EPhysFilt)*Filt.Gap/2 + (Filt.LocEPhys-1)*Filt.Gap,'LineWidth',0.5,'Color',[0.25 0.25 0.25]);    
    end
    %%
    T3_N = zeros([numFrames N+1]);

    for lidx = 1:N%[1:N]
        RP1 = circshift(circshift( ResP1, Filt.xOff,1), Filt.yOff,2);
        RP2 = circshift(circshift( ResP2, Filt.xOff,1), Filt.yOff,2);
        T1 = reshape( RP1(Filt.PhaseOff+1:end-Filt.PhaseOff,Filt.PhaseOff+1:end-Filt.PhaseOff,lidx,:),...
            [((Filt.AbsWin-2*Filt.PhaseOff)+1).^2 numFrames]);
        T2 = reshape( circshift(circshift( RP2(Filt.PhaseOff+1:end-Filt.PhaseOff,Filt.PhaseOff+1:end-Filt.PhaseOff,lidx,:), Filt.xOff,1), Filt.yOff,2),...
            [((Filt.AbsWin-2*Filt.PhaseOff)+1).^2 numFrames]);

        if Filt.PlotOption == 1

            uwOff = 0;
            TMed = 0*(sqrt(T1Med.^2 + T2Med.^2));
            T = (sqrt(T1.^2 + T2.^2));
%             T = T./abs(T);
            Tfilt = filtfilt(b,a,angle( mean( T,1).*exp(1i.*uwOff) .*exp(-1i*angle(TMed))));

            count1 = 1;
            while( (max( Tfilt(2000:end-2000)) - min( Tfilt(2000:end-2000))) > 1 && count1 <= 200)
                uwOff = uwOff - 0.2;
                Tfilt = filtfilt(b,a,angle( mean(T,1).*exp(uwOff .*1i).*exp(-1i.*angle(TMed))));
                count1 = count1 + 1;
            end

            taxis2 = taxis;
            TFiltShow = abs(Tfilt - filtfilt( bMean, aMean, Tfilt))*850/(4*pi*1.36);
            NoiseLev = 1;%std(TFiltShow(2000:end-1000));
            XFFF = NoiseLev.*ones( size( TFiltShow));
            if max(TFiltShow(2000:end-2000)) - min(TFiltShow(2000:end-2000)) > 200
                TFiltShow = 0*TFiltShow;
                XFFF = sin(2.*pi.*(1:length(XFFF))/100);
                XFFF( XFFF >= 0) = 1;
                XFFF( XFFF < 0) = NaN;
            end
%             
            NoiseLev = std(TFiltShow(2000:end-1000));
            Pulses = abs(TFiltShow) > NoiseLev*3;

            PulseLoc = taxis(Pulses);

%             scatter( PulseLoc*1000, ones([length(PulseLoc) 1]) + Filt.Gap*(lidx + 0.9), 8,[0.5 0.5 0.5],'filled','MarkerEdgeColor','none');
            plot(taxis2*1000, Filt.Gap/2 + TFiltShow./XFFF + Filt.Gap*lidx,'Color',ColorEMap(lidx,:),'LineWidth',1);
            
            if rem(lidx-1,10) == 0
                text(2200,mean(TFiltShow(500:end-500))+ Filt.Gap*lidx + Filt.Gap/2,['ROI ' num2str(lidx)],'FontSize',12,'FontName','Segoe UI');
            end

            T3_N(:,lidx) = TFiltShow./XFFF;%TFiltShow/ NoiseLev;
            StringOpt = 'Filtered \angle \surd( E_{P1}^2 + E_{P2}^2) (nm)';
            
            [acor,lag] = xcorr( EPhysFilt(4400:numFrames-2400)/3,TFiltShow(4400:end-2400)/NoiseLev,400,'unbiased');

            [CorrVals(lidx,1),I] = max(abs(acor));
            CorrVals(lidx,2) = lag(I)/4;
            disp([num2str(lidx) '|' num2str(CorrVals(lidx,1),'%1.3f') ' | ' num2str(CorrVals(lidx,2),'%1.3f') ' ms']);

        elseif Filt.PlotOption == 2
   
            uwOff = 0;
            T1filt = filtfilt(b,a,angle( mean(T1,1).*exp(1i.*uwOff) .*exp(-1i*angle(0*T1Med))));

            count1 = 1;
            while( (max( T1filt(1000:end-1000)) - min( T1filt(1000:end-1000))) > 2 && count1 <= 200)
                uwOff = uwOff - 0.2;
                T1filt = filtfilt(b,a,angle( mean(T1,1).*exp(uwOff .*1i).*exp(-1i.*angle(0*T1Med))));
                count1 = count1 + 1;
            end

            T2filt = filtfilt(b,a,angle( mean(T2,1).*exp(1i.*uwOff) .*exp(-1i*angle(0*T2Med))));

            count1 = 1;
            while( (max( T2filt(1000:end-1000)) - min( T2filt(1000:end-1000))) > 2 && count1 <= 200)
                uwOff = uwOff - 0.2;
                T2filt = filtfilt(b,a,angle( mean(T2,1).*exp(uwOff .*1i).*exp(-1i.*angle(0*T2Med))));
                count1 = count1 + 1;
            end

            taxis2 = taxis;
            TFiltShow = -( abs(T1filt - filtfilt( bMean, aMean, T1filt)) + ...
                abs(T2filt - filtfilt( bMean, aMean, T2filt)))*850/(4*pi*1.36);
%             TFiltShow = abs(TFiltShow - mean( TFiltShow(5000:8800)));
            
%             TFiltShow = filtfilt(bStop,aStop,TFiltShow);
            if max(TFiltShow(1000:end-1000)) - min(TFiltShow(1000:end-1000)) > 100
                TFiltShow = NaN*TFiltShow;
            end
             
            NoiseLev = 1;%std(T2FiltShow(1000:end-1000));
            Pulses = abs(TFiltShow) > NoiseLev*5;

            plot(taxis2*1000, Filt.Gap/2 + TFiltShow./NoiseLev + Filt.Gap*lidx,'Color',ColorEMap(lidx,:),'LineWidth',1);

            if rem(lidx-1,10) == 0
                text(2500,mean(TFiltShow(500:end-500))+ Filt.Gap*lidx + Filt.Gap/2,['ROI ' num2str(lidx)],'FontSize',12,'FontName','Segoe UI');
            end
            
            T3SS = zeros(size(TFiltShow));
            T3_N(:,lidx) = TFiltShow./NoiseLev;
            StringOpt = 'Filtered phase plots (\angleE_{P1} + \angleE_{P2}) in (nm)';
            
            Sig1 = EPhysFilt(4800:8400);
            Sig2 = TFiltShow(4800:8400);
            [acor,lag] = xcorr(  (Sig1 - min(Sig1))./(max(Sig1) - min(Sig1)), ...
                (Sig2 - min(Sig2))./(max(Sig2) - min(Sig2)),...
                400);

            [CorrVals(lidx,1),I] = max(abs(acor));
            CorrVals(lidx,2) = lag(I)/4;
            disp([num2str(lidx) '|' num2str(CorrVals(lidx,1),'%1.3f') ' | ' num2str(CorrVals(lidx,2),'%1.3f') ' ms']);

        elseif Filt.PlotOption == 3
   
            uwOff = 0;
            T1filt = (filtfilt(b,a,unwrap(angle( mean(T1,1).*exp(uwOff .*1i).*exp(-1i.*angle(T1Med))))));

            count1 = 1;
%             while( std( abs(T1filt(300:end-300))) > 20 && count < 20)
            while( (max( T1filt(1000:end-1000)) - min( T1filt(1000:end-1000))) > 4 && count1 <= 20)
                uwOff = uwOff - 0.5;
                T1filt = (filtfilt(b,a,unwrap(angle( mean(T1,1).*exp(uwOff .*1i).*exp(-1i.*angle(T1Med))))));
                count1 = count1 + 1;
            end
            
            uwOff = 0;
            count2 = 1;
            T2filt = (filtfilt(b,a,unwrap(angle( mean(T2,1).*exp(uwOff .*1i).*exp(-1i.*angle(T2Med))))));
%             while( std( abs(T2filt(300:end-300))) > 20 && count < 20)
            while( (max( T2filt(1000:end-1000)) - min( T2filt(1000:end-1000))) > 4 && count2 <= 20)
               uwOff = uwOff - 0.5;
               T2filt = (filtfilt(b,a,unwrap(angle( mean(T2,1).*exp(uwOff .*1i).*exp(-1i.*angle(T2Med))))));
               count2 = count2 + 1;
            end
            
            if count1 >= 20 || count2 >= 20
                T1filt = T1filt*0;
                T2filt = T2filt*0;
            end
          
            T1filt = cumsum(angle( exp(1i.*T1filt).*exp(-1i.*circshift(T1filt,-1))));
%             T1LFP = circshift(filtfilt( bMean, aMean, T1filt),-0*Filt.windowSizeMean/2);
            T2filt = cumsum(angle( exp(1i.*T2filt).*exp(-1i.*circshift(T2filt,-1))));
%             T2LFP = circshift(filtfilt( bMean, aMean, T2filt),-0*Filt.windowSizeMean/2);
%             T2FiltShow = (abs(T2filt - T2LFP) + abs(T1filt - T1LFP))*850/(4*pi*1.36);
            TFiltShow = (angle( exp(1i.*T2filt).*exp(1i.*T1filt)))*850/(4*pi*1.36);
            TFiltShow_LFP = circshift(filtfilt( bMean, aMean, TFiltShow),0);
            TFiltShow = abs( TFiltShow - TFiltShow_LFP);
            NoiseLev = std(TFiltShow(1000:end-1000));
            
            if max(TFiltShow(1000:end-1000)) - min(TFiltShow(1000:end-1000)) > 20
                TFiltShow = 0*TFiltShow;
            end
                       
            plot(taxis*1000, Filt.Gap/2 + (TFiltShow - mean(TFiltShow(1000:end-1000))+ Filt.Gap*lidx),'Color',ColorEMap(lidx,:),'LineWidth',1);

            if rem(lidx-1,10) == 0
                text(2500,mean(TFiltShow(1000:end-1000))+ Filt.Gap*lidx + Filt.Gap/2,['ROI ' num2str(lidx)],'FontSize',13,'FontName','Segoe UI');
            end
            
%             ylabel('High-pass filtered displacement from Pol. 2 (nm)','Color','k');
            T3_N(:,lidx) = (TFiltShow - mean(TFiltShow(1000:end-1000)))./NoiseLev;
            StringOpt = 'Filtered phase plots (\angleE_{P1} + \angleE_{P2}) in (nm)';
            [acor,lag] = xcorr( EPhysFilt(6400:numFrames-1000)/3,TFiltShow(6400:numFrames-1000)/mean( abs(TFiltShow(6400:end-1000))),'unbiased');

            [CorrVals(lidx,1),I] = max(abs(acor));
            CorrVals(lidx,2) = lag(I)/4;

        elseif Filt.PlotOption == 4
            MF = [1 1];
            T1 = medfilt2( real(T1),MF) + 1i.*medfilt2( imag(T1),MF);
            T1D = (real(mean(T1,1)))/100;
            T1D = (T1D - circshift(T1D,1))/1;
            T1filt = cumsum(filtfilt(b,a,T1D));

            T2 = medfilt2( real(T2),MF) + 1i.*medfilt2( imag(T2),MF);
            T2D = ((mean(T2,1)))/100;
            T2D = cumsum(T2D - circshift(T2D,1))/1;
            T2filt = real(filtfilt(b,a,T2D))./abs(filtfilt(b,a,T2D));
            
            plot(taxis*1000, Filt.Gap/2 + (T2filt - mean(T2filt(500:end-500))+ Filt.Gap*lidx),'Color',ColorEMap(lidx,:),'LineWidth',1);

            if rem(lidx-1,10) == 0
                text(2300,mean(T1filt(500:end-500))+ Filt.Gap*lidx + Filt.Gap/2,['ROI ' num2str(lidx)],'FontSize',12,'FontName','Segoe UI');
            end
            
            ylabel('High-pass filtered real differential from Pol. 2 (nm)','Color','k');
            T3_N(:,lidx) = T1filt - mean(T1filt(200:end-200));
            StringOpt = 'Real differential';

            [acor,lag] = xcorr( EPhysFilt(1:numFrames),T2filt/mean( abs(T2filt(:))) );

            [MM,I] = max(abs(acor));
            lagDiff = lag(I);
            disp([num2str(lidx) '|' num2str(MM,'%1.3f') ' | ' num2str(lagDiff/4000,'%1.3f') ' ms']);
        elseif Filt.PlotOption == 5
            T1filt = (mean(T1,1));
            T2filt = (mean(T2,1));

    %         T3 = 10*(abs( T1filt) - abs(T2filt))./(abs(T1filt)+abs(T2filt));
            T3 = 10*( abs( T1filt + T2filt)) / 1E4;
            T3LFP = circshift(filtfilt( bMean, aMean, T3),-windowSizeMean/2);

            T3 = abs(T3 - T3LFP );

            plot(taxis*1000, (T3 - T3(100)  + Filt.Gap*lidx),'Color',BlueMap(lidx,:),'LineWidth',0.5);   

            text(10,10+(T3(400) - T3(end) + Filt.Gap*lidx),['ROI ' num2str(lidx)],'FontSize',15);
            ylabel('Difference between the two magnitudes (a.u.)');

            StringOpt = 'Magnitude difference';
            T3_N(:,lidx) = T3 - mean(T3(500:end));

        elseif Filt.PlotOption == 6
            T1filt = mean(filtfilt(b,a,abs(T1)),1);
            T2filt = mean(filtfilt(b,a,abs(T2)),1);

            T3 = 120*(abs( T1filt) + abs(T2filt))/20000;%./(abs(T1filt)+abs(T2filt));

            plot(taxis*1000, (T3 - T3(2000)  + Filt.Gap*lidx),'Color',BlueMap(lidx,:),'LineWidth',1);   

            text(10,Filt.Gap/4+(T3(2500) - T3(2000) + Filt.Gap*lidx),['ROI ' num2str(lidx)],'FontSize',15);
            ylabel('Filtered Normalized Difference between the two magnitudes (a.u.)');

            StringOpt = 'Filtered Magnitude difference';
            T3_N(:,lidx) = T3 - mean(T3);
        elseif Filt.PlotOption == 7
            MF = [1 1];
            uwOff = -1;
            T1 = medfilt2( real(T1),MF) + 1i.*medfilt2( imag(T1),MF);
            T1D = (angle(mean(T1,1).*exp(uwOff .*1i)) - angle(T1Med) + 0.2)*850/(4*pi*1.36);
            count = 1;
            while( mean( T1D(300:end) > 50)&& count < 30)
                uwOff = uwOff - 1;
                T1D = (angle(mean(T1,1).*exp(uwOff .*1i)) - angle(T1Med) + 0.2)*850/(4*pi*1.36);
                count = count + 1;
            end
            T1D = T1D - circshift(T1D,1);
            T1filt = cumsum(filtfilt(b,a,T1D));
            T1LFP = circshift(filtfilt( bMean, aMean, T1filt),-Filt.windowSizeMean/2);
            T1filt = T1filt - T1LFP;

            uwOff = -1;
            count = 1;
            T2 = medfilt2( real(T2),MF) + 1i.*medfilt2( imag(T2),MF);
            T2D = (angle(mean(T2,1).*exp(uwOff .*1i)) - angle(T2Med) + 0.2)*850/(4*pi*1.36);
            while( mean( T2D(300:end) > 50)&& count < 30)
                uwOff = uwOff - 1;
                T2D = (angle(mean(T2,1).*exp(uwOff .*1i)) - angle(T2Med) + 0.2)*850/(4*pi*1.36);
                count = count + 1;
            end
            T2D = T2D - circshift(T2D,1);
            T2filt = cumsum(filtfilt(b,a,T2D));
            T2LFP = circshift(filtfilt( bMean, aMean, T2filt),-Filt.windowSizeMean/2);
            T2filt = T2filt - T2LFP;

            plot(taxis*1000, (T1filt - T1filt(100)  + Filt.Gap*lidx),'Color',BlueMap(lidx,:),'LineWidth',0.5);   
            plot(taxis*1000, Filt.Gap/2 + (T2filt - T2filt(100)+ Filt.Gap*lidx),'Color',RedMap(lidx,:),'LineWidth',0.5);

            text(10,6+(T2filt(2000) - T2filt(end) + Filt.Gap*lidx),['ROI ' num2str(lidx)],'FontSize',15);
            ylabel('Filtered Displacement from phase of P1 (green) or P2 (pink) images (nm)');  
            T3_N(:,lidx) = T1filt - mean(T1filt(500:end-500));
            StringOpt = 'Phase-diff Filtered Phase plots';    
        elseif Filt.PlotOption == 8
            MF = [2 2];
            uwOff = -1;
            count = 1;
%             T1 = medfilt2( real(T1),MF) + 1i.*medfilt2( imag(T1),MF);
            T1D = (angle(mean(T1,1).*exp(uwOff .*1i)) - angle(0*T1Med))*850/(4*pi*1.36);
            while( (mean( abs(T1D(1600:6400))) > 20)&& count < 20)
                uwOff = uwOff - 1;
                T1D = (angle(mean(T1,1).*exp(uwOff .*1i)) - angle(0*T1Med))*850/(4*pi*1.36);
                count = count + 1;
            end
            
            uwOff = -1;
%             T2 = medfilt2( real(T2),MF) + 1i.*medfilt2( imag(T2),MF);
            T2D = (angle(mean(T2,1).*exp(uwOff .*1i)) - angle(0*T2Med))*850/(4*pi*1.36);
            count = 1;
            while( (mean( abs(T2D(1600:6400))) > 20) && count < 20)
                uwOff = uwOff - 1;
                T2D = (angle(mean(T2,1).*exp(uwOff .*1i)) - angle(0*T2Med))*850/(4*pi*1.36);
                count = count + 1;
            end
            if (std( abs(T2D(1600:6400))) > 15) || (std( abs(T1D(1600:6400))) > 15)
                T2D = 0*T2D;
                T1D = 0*T1D;
            end
            T1D = (T1D - circshift(T1D,1))/1;
            T1filt = (filtfilt(b,a,T1D));
            T2D = (T2D - circshift(T2D,1))/1;
            T2filt = (filtfilt(b,a,T2D));

%             T3 = (cumsum(angle(exp(1i.*T2filt).*exp(-1i.*T1filt))));
%             T3 = cumsum(0*T2filt - T1filt);
            T3 = (cumsum(T2filt + T1filt));% - T1filt;
            T3LFP = circshift(filtfilt( bMean, aMean, T3), 0);
            T3 = abs(T3 - T3LFP);
            
            NoiseLev = std(T3(1600:6400));
            Pulses = abs(T3) > NoiseLev*4;
            PulseLoc = taxis(Pulses);
            
%             scatter( PulseLoc*1000, ones([length(PulseLoc) 1]) + + Filt.Gap*(lidx + 0.7), 5,[0.5 0.5 0.5],'filled','MarkerEdgeColor','none')
            plot(taxis*1000, (T3   + Filt.Gap*lidx),'Color',ColorEMap(lidx,:),'LineWidth',0.5);      
            
            if rem(lidx,10) == 0
                text(2600,mean(T3(1600:6400))+ Filt.Gap*lidx,['ROI ' num2str(lidx)],'FontSize',0.75*Filt.Gap,'FontName','Segoe UI');
            end
                ylabel('Difference between the two phases (nm)');     

            T3_N(:,lidx) = T3;
            StringOpt = 'Phase-diff Filtered Phase diff plots';      
        
        elseif Filt.PlotOption == 9
%             MF = [1 1];
%             T1 = medfilt2( real(T1),MF) + 1i.*medfilt2( imag(T1),MF);
            T1D = T1;
%             T1D = (T1D - circshift(T1D,1,2))/1;
            T1filt = (mean((filtfilt(b,a,T1D,[],2)),1));

%             T2 = medfilt2( real(T2),MF) + 1i.*medfilt2( imag(T2),MF);
%             T2D = (abs(mean(T2,1)) - 0*abs(T2Med))/100;
%             T2D = (T2D - circshift(T2D,1))/1;
%             T2filt = abs(cumsum(filtfilt(b,a,T2D)));
            T2D = T2;
%             T2D = (T2D - circshift(T2D,1,2))/1;
            T2filt = (mean(filtfilt(b,a,T2D,[],2),1));

%             T3 = abs(T2filt -  T1filt)./abs(T2filt + T1filt);
%             T3 = abs(sqrt(T1filt.^2 + T2filt.^2))/1000;
            T3 = abs(T1filt) - abs(T2filt);

            T3LFP = filtfilt( bMean, aMean, T3);

            T3 = abs(T3 - T3LFP)/250 + 0*T3LFP; % abs(T3 - median(T3(500:end)));

            NoiseLev = std(T3(500:end-500));
            Pulses = abs(T3) > NoiseLev*4;

            PulseLoc = taxis(Pulses);

%            scatter( PulseLoc*1000, ones([length(PulseLoc) 1]) + Filt.Gap*(lidx + 0.75), 8,[0.5 0.5 0.5],'filled','MarkerEdgeColor','none')
            plot(taxis*1000, (T3  + Filt.Gap*lidx),'Color',ColorEMap(lidx,:),'LineWidth',1);
            if rem(lidx-1,10) == 0
                text(2500,mean(T3(500:end-500))+ Filt.Gap*lidx,['ROI ' num2str(lidx)],'FontSize',13,'FontName','Segoe UI');
            end
            ylabel('Difference between the two Abs (nm)');

            T3_N(:,lidx) = T3;
            T3_N( (abs(T3) < NoiseLev),lidx) = 0; 
            StringOpt = 'Difference between the two Abs (nm)';   
            
%             [acor,lag] = xcorr( EPhysFilt(400:numFrames-400)/3,T3(400:end-400)/mean( abs(T3(400:end-400))) );
% 
%             [MM,I] = max(abs(acor));
%             lagDiff = lag(I);
%             disp([num2str(lidx) '|' num2str(MM,'%1.3f') ' | ' num2str(lagDiff/4,'%1.3f') ' ms']);
        end
        pause(0.01);
    end
    xlim([30 2680]);

    set(gca,'LineWidth',1.5,'FontSize',15,'FontName','Segoe UI');
    box on;
    xlabel('Time (ms)');
    set(gca,'Position',get(gca,'Position').*[1 1 1.1 1.0])
    
    if Filt.ShowEPhys == 1
        T3_N(:,2:end) = T3_N(:,1:end-1);
        T3_N(:,1) = EPhysFilt(1:numFrames)*Filt.Gap/2;
    end
    title({['NE-4C cells | ' FileDetails.fDEEET] ;[StringOpt ' | ' FiltOpt]},'Interpreter','tex');
    set(gca,'Position',get(gca,'Position')-[0.025 0 0 0]);
    xlim([100 2600]);
%     ylim([0 1600]);
    set(gca,'YTick',[],'YColor',[1 1 1]);
%     plot([902 902],[20 30],'k','LineWidth',3);
    ylabel(StringOpt, 'Color','k','Interpreter','tex','FontName','Segoe UI');
    if Filt.SaveFig == 1
        savefig(gcf,[pwd '\' FileDetails.fbase '\' StringOpt '_' FiltOpt '_' FileDetails.FileExt '.fig']);
    end

end