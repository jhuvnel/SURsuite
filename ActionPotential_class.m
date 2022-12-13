classdef ActionPotential_class <handle
    properties
        DAT
    end
    
    methods
        % Constructor takes a SpikeDATA object.
        function obj = ActionPotential_class(SpikeDATA_object)
            obj.DAT = SpikeDATA_object;
        end
        
        function STR = GetFindAPs(obj, DAT, THRESH)
            % Subtract out baseline offset.
            DAT = DAT - mean(DAT);
            
            % SMOOTH the raw head-stage data first...
            SMOO = smooth(DAT,21,'sgolay');
            
            % THEN "differentiate" using the shifted subtraction method.
            SS = obj.DoShiftSubtract(SMOO);
            
            % SQUARE the signal, but PRESERVE the sign.
            SS = (SS.*SS) .* sign(SS);
            
            % This is HARD CODED to look for NEGATIVE PEAKS!!
            % NEED TO make this auto-detect the AP polarity.
            
            % If a specific threshold is not given then use a Moving
            % Maximum to set the "threshold" for detecting peaks.
            if nargin < 3
                MM = -movmax(-SS, 10000);
                THRESH = 0.3*MM;
            end
            
            % "Threshold" based on moving maximum.
            HIT = SS < THRESH;
            HIT_idx = find(HIT);
            
            % Find the "clumps" of data points that are above the
            % threshold. Each Clump may be an action potential.
            CLUMP_idx = diff(HIT_idx) > 100;
            CLUMP = HIT_idx(CLUMP_idx);
            
            % Loop through the clumps, and save an index for the Action
            % Potential, that is at the max() value of the smoothed data
            % (i.e., the "peak" of the Action Potential).
            AP = zeros(length(CLUMP), 1);
            for i=1:length(CLUMP)
                vals = SMOO((1:200)+CLUMP(i));
                apidx = find(vals == max(vals));
                apidx = apidx(1);
                AP(i) = apidx + CLUMP(i);
            end
            
            % Get rid of any low-amplitude outliers, whose peak is well
            % below the mean of all the peaks.
            AP(SMOO(AP) < 0.35*mean(SMOO(AP))) = [];
            
            % Compute the AP Template, as the average of all the APs. Plot
            % it.
            IDX = -100:1000;
            AVG = zeros(length(IDX), 1);
            N = length(AP);
            for i=1:N
                if AP(i) + IDX(1) < 1; continue; end
                if AP(i) + IDX(end) > length(SMOO); break; end
                SIG = SMOO(AP(i)+IDX);
                %line(IDX, SIG);
                AVG = AVG + SIG;
            end
            
            AVG = AVG ./ N;
            
            % Return a bunch of internal values in a struct, for debugging
            % or plotting.
            STR = MakeStruct(DAT, SS, CLUMP, AP, SMOO, AVG);
        end
        
        % Find/mark Action Potentials in analog data.
        %
        % Optionally pass a fixed Threshold. Otherwise, a moving maximum
        % threshold is computed automatically.
        %
        % Return a struct with the various internal computed values.
        function [APT_unzeroed, APT] = FindAPs(obj, DAT, THRESH)
            % Subtract out baseline offset.
            DAT_unzeroed = DAT;
            DAT = DAT - mean(DAT);

            % SMOOTH the raw head-stage data first...
                        SMOO_unzeroed = smooth(DAT_unzeroed,21,'sgolay');
                        SMOO = smooth(DAT,21,'sgolay');
                        [pospks,lposocs,posw,posp] = findpeaks(SMOO,1:length(SMOO),'minpeakheight',.05,'minpeakdistance',900);
            %             [negpks,neglocs,negw,negp] = findpeaks(-SMOO,1:length(SMOO),'minpeakheight',.05,'minpeakdistance',900);
            %
            %             if obj.DAT.APorientation == 0
            %                 if mean(pospks) > mean(negpks)
            %                     obj.DAT.APorientation = 1;
            %                 else
            %                     obj.DAT.APorientation = -1;
            %                 end
            %             end
            % THEN "differentiate" using the shifted subtraction method.
            SS_unzeroed = obj.DoShiftSubtract(SMOO_unzeroed);
            SS = obj.DoShiftSubtract(SMOO);

            % SQUARE the signal, but PRESERVE the sign.
            SS_unzeroed = (SS_unzeroed.*SS_unzeroed) .* sign(SS_unzeroed);
            SS = (SS.*SS) .* sign(SS);


            % If a specific threshold is not given then use a Moving
            % Maximum to set the "threshold" for detecting peaks.

            if nargin < 3
                MM_unzeroed = -movmax(-SS_unzeroed, 10000);
                THRESH_unzeroed = 0.3*MM_unzeroed;

                MM = -movmax(-SS, 10000);
                THRESH = 0.3*MM;
            end

            % "Threshold" based on moving maximum.
            HIT_unzeroed = SS_unzeroed < THRESH_unzeroed;
            HIT_idx_unzeroed = find(HIT_unzeroed);

            HIT = SS < THRESH;
            HIT_idx = find(HIT);



            % Find the "clumps" of data points that are above the
            % threshold. Each Clump may be an action potential.
            CLUMP_idx_unzeroed = diff([1; HIT_idx_unzeroed]) > 100;
            CLUMP_unzeroed = HIT_idx_unzeroed(CLUMP_idx_unzeroed);

            CLUMP_idx = diff([1; HIT_idx]) > 100;
            CLUMP = HIT_idx(CLUMP_idx);

            % Loop through the clumps, and save an index for the Action
            % Potential, that is at the max() value of the smoothed data
            % (i.e., the "peak" of the Action Potential).
            AP_unzeroed = zeros(length(CLUMP_unzeroed), 1);
            AP = zeros(length(CLUMP), 1);

            for i=1:length(CLUMP)
                bd = (1:200)+CLUMP(i);
                if bd(end) > length(SMOO)
                    bd = 1+CLUMP(i):length(SMOO);
                end
                vals = SMOO(bd);
                apidx = find(vals == max(vals));
                apidx = apidx(1);
                AP(i) = apidx + CLUMP(i);
            end

            for i=1:length(CLUMP_unzeroed)
                bdU = (1:200)+CLUMP_unzeroed(i);
                if bdU(end) > length(SMOO)
                    bdU = 1+CLUMP_unzeroed(i):length(SMOO);
                end
                vals_unzeroed = SMOO_unzeroed(bdU);
                apidx_unzeroed = find(vals_unzeroed == max(vals_unzeroed));
                apidx_unzeroed = apidx_unzeroed(1);
                AP_unzeroed(i) = apidx_unzeroed + CLUMP_unzeroed(i);
            end

            % Get rid of any low-amplitude outliers, whose peak is well
            % below the mean of all the peaks.
            AP_unzeroed(SMOO_unzeroed(AP_unzeroed) < 0.35*mean(SMOO_unzeroed(AP_unzeroed))) = [];
            AP_unzeroed(SMOO_unzeroed(AP_unzeroed) > 3*mean(SMOO_unzeroed(AP_unzeroed))) = [];
            AP(SMOO(AP) < 0.35*mean(SMOO(AP))) = [];
            AP(SMOO(AP) > 3*mean(SMOO(AP))) = [];

            % Compute the AP Template, as the average of all the APs. Plot
            % it.
            IDX = -100:1000;
            AVG = zeros(length(IDX), 1);
            N = length(AP);
            AVG_RAW = zeros(length(IDX), 1);

            for i=1:N
                if AP(i) + IDX(1) < 1; continue; end
                if AP(i) + IDX(end) > length(SMOO); break; end
                SIG = SMOO(AP(i)+IDX);
                SIG_RAW = DAT(AP(i)+IDX);
                %line(IDX, SIG);
                AVG = AVG + SIG;
                AVG_RAW = AVG_RAW + SIG_RAW;
            end

            AVG_unzeroed = zeros(length(IDX), 1);
            AVG_RAW_unzeroed = zeros(length(IDX), 1);
            N_unzeroed = length(AP_unzeroed);
            for i=1:N_unzeroed
                if AP_unzeroed(i) + IDX(1) < 1; continue; end
                if AP_unzeroed(i) + IDX(end) > length(SMOO_unzeroed); break; end
                SIG_unzeroed = SMOO_unzeroed(AP_unzeroed(i)+IDX);
                SIG_RAW_unzeroed = DAT_unzeroed(AP_unzeroed(i)+IDX);
                %line(IDX, SIG);
                AVG_unzeroed = AVG_unzeroed + SIG_unzeroed;
                AVG_RAW_unzeroed = AVG_RAW_unzeroed + SIG_RAW_unzeroed;
            end

            if max(AVG ./ N) < 0.05
                N = length(pospks)*2;
            end
            AVG = AVG ./ N;
            AVG_RAW = AVG_RAW ./ N;
            AVG_unzeroed = AVG_unzeroed ./ N_unzeroed;
            AVG_RAW_unzeroed = AVG_RAW_unzeroed ./ N_unzeroed;

            % Return a bunch of internal values in a struct, for debugging
            % or plotting.
            APT_unzeroed = MakeStruct(SS_unzeroed, CLUMP_unzeroed, AP_unzeroed, SMOO_unzeroed, AVG_unzeroed, AVG_RAW_unzeroed);
            APT = MakeStruct(SS, CLUMP, AP, SMOO, AVG, AVG_RAW);
        end
        
        function [samps, time, startPt, endPt, peakPolarization1, peakPolarization2, maxHyperpol] = findAPWidth(obj, apTemp, chan, debugPlot, debugPlots)
            if isprop(obj.DAT,'APRefProperites')
                if obj.DAT.APRefProperites.(chan).peakDePol2 == 0
                    chanUsed = obj.DAT.GetChan(chan);
                    ap = apTemp.AVG_RAW;
                    apfilt = smooth(ap,15,'sgolay');
                    if debugPlot
                        startPt=obj.DAT.APRefProperites.(chan).startPt;

                        pt = obj.DAT.APRefProperites.(chan).peakDePol1;
                        [~,ml]=max(ap(pt-5:pt+5));
                        peakPolarization1 = ml+pt-6;
                        peakPolarization2 = 0;

                        mt = obj.DAT.APRefProperites.(chan).maxHyperpol;
                        [~,ml]=min(apfilt(mt-6:mt+6));
                        maxHyperpol = ml+mt-7;

                        endPt = find(apfilt(maxHyperpol:end)>mean(apfilt)/5,1,'first')+maxHyperpol;
                        if endPt - maxHyperpol < 5
                            [aa,aab]=min(abs(apfilt(maxHyperpol+5:400)-mean(apfilt(1:600))));
                            endPt = aab+maxHyperpol+5-1;
                        elseif endPt-maxHyperpol > 30
                            endPt = maxHyperpol+30;
                        elseif isempty(endPt)
                            endPt = maxHyperpol+30;
                        end

                        samps = endPt-startPt;
                        time = samps/chanUsed.SampleRate;

                        if ~isempty(debugPlots)
                            plot(debugPlots.apPoints,ap)
                            hold(debugPlots.apPoints,'on')
                            plot(debugPlots.apPoints,apfilt)
                            plot(debugPlots.apPoints,startPt,apfilt(startPt),'r*')
                            plot(debugPlots.apPoints,peakPolarization1,apfilt(peakPolarization1),'g*')
                            plot(debugPlots.apPoints,maxHyperpol,apfilt(maxHyperpol),'k*')
                            plot(debugPlots.apPoints,endPt,apfilt(endPt),'m*')

                        else
                            figure
                            plot(ap)
                            hold on
                            plot(apfilt)
                        end
                    end
                else
                end
            else
                %% If AP Properties were not picked manually 
                apTemp = apTemp.AVG_RAW;
                load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\APStyle1.mat');
                load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\APStyle2.mat');
                load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\APStyle3.mat');
                chanUsed = obj.DAT.GetChan(chan);
                ap = double(apTemp);
                [~,~,dt1]=findsignal(APStyle1,ap,'Normalization','zscore');
                [~,~,dt2]=findsignal(APStyle2,ap,'Normalization','zscore');
                [~,~,dt3]=findsignal(APStyle3,ap,'Normalization','zscore');
                %             tf = 3;
                %             apfilt = filtfilt(ones(1,3)/3,1,ap);
                %             snrv = snr(apfilt);
                %             while abs(snrv) > 4
                %                 tf = tf + 2;
                %                 apfilt = filtfilt(ones(1,tf)/tf,1,ap);
                %                 snrv = snr(apfilt);
                %             end
                [~, closestN] = min([dt1 dt2 dt3]);

                switch closestN
                    case 1
                        apfilt = smooth(ap,15,'sgolay');
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,ap)
                                hold(debugPlots.apPoints,'on')
                                plot(debugPlots.apPoints,apfilt)
                            else
                                figure
                                plot(ap)
                                hold on
                                plot(apfilt)
                            end
                        end
                        [a,b]=islocalmax(apfilt);
                        [B,I] = max(b);
                        maxpt1 = I;
                        maxpt2 = [];
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,maxpt1,apfilt(maxpt1),'g*')
                            else
                                plot(maxpt1,apfilt(maxpt1),'g*'); %max pt of AP
                            end
                        end
                        [t, p] = islocalmin(apfilt(1:maxpt1));
                        [firstminV, firstminL] = max(find(t));
                        sigmean = mean(ap)*10;
                        if apfilt(firstminV) > sigmean
                            realLim = find(apfilt(1:maxpt1)<sigmean,1,'last');
                            [t, p] = islocalmin(apfilt(1:realLim));
                            [firstminV, firstminL] = max(find(t));
                        end
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,firstminV,apfilt(firstminV),'r*')
                            else
                                plot(firstminV,apfilt(firstminV),'r*'); %First corner pt (start of AP)
                            end
                        end

                        [t1, p1] = islocalmin(apfilt);
                        [secondminV, secondminL] = max(p1);
                        if secondminL > 600
                            [t1, p1] = islocalmin(apfilt(1:600));
                            [secondminV, secondminL] = max(p1);
                        end
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,secondminL,apfilt(secondminL),'k*')
                            else
                                plot(secondminL,apfilt(secondminL),'k*'); %Min after peak of AP (Hyperpolarization)
                            end
                        end

                        endP = find(apfilt(secondminL:end)>mean(apfilt)/5,1,'first')+secondminL;
                        if endP - secondminL < 5
                            endP = find(apfilt(secondminL:end)>apfilt(1)+mean(apfilt),1,'first')+secondminL;
                        end
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,endP,apfilt(endP),'m*')
                                hold(debugPlots.apPoints,'off')
                            else
                                plot(endP,apfilt(endP),'m*'); % end point of recovery period of AP
                            end
                        end

                        samps = secondminL - firstminV;
                        time = samps/chanUsed.SampleRate;
                        startPt = firstminV;
                        endPt = endP;
                        maxHyperpol = secondminL;%minPt = secondminL;
                        peakPolarization1 = maxpt1;
                        peakPolarization2 = [];

                    case 2
                        apfilt = smooth(ap,15,'sgolay');

                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,ap)
                                hold(debugPlots.apPoints,'on')
                                plot(debugPlots.apPoints,apfilt)
                            else
                                figure
                                plot(ap)
                                hold on
                                plot(apfilt)
                            end
                        end

                        [a,b]=islocalmax(apfilt);
                        [B,I] = maxk(b,2);
                        if I(1) < I(2)
                            maxpt1 = I(1);
                            maxpt2 = I(2);
                        else
                            maxpt1 = I(2);
                            maxpt2 = I(1);
                        end

                        %             [vm, lm] = max(apfilt);
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,maxpt1,apfilt(maxpt1),'g*')
                                plot(debugPlots.apPoints,maxpt2,apfilt(maxpt2),'b*')
                            else
                                plot(maxpt1,apfilt(maxpt1),'g*');
                                plot(maxpt2,apfilt(maxpt2),'b*');
                            end
                        end

                        [t, p] = islocalmin(apfilt(1:maxpt1));
                        [firstminV, firstminL] = max(find(t));

                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,firstminV,apfilt(firstminV),'r*')
                            else
                                plot(firstminV,apfilt(firstminV),'r*');
                            end
                        end

                        [t1, p1] = islocalmin(apfilt);
                        [secondminV, secondminL] = max(p1);
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,secondminL,apfilt(secondminL),'k*')
                            else
                                plot(secondminL,apfilt(secondminL),'k*');
                            end
                        end

                        difap = [diff(apfilt);0];
                        thirdminL = find(difap(maxpt2:end)>0,1)+maxpt2;
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,thirdminL,apfilt(thirdminL),'y*')
                            else
                                plot(thirdminL,apfilt(thirdminL),'y*');
                            end
                        end

                        endP = find(apfilt(thirdminL:end)>mean(apfilt)/5,1,'first')+thirdminL;

                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,endP,apfilt(endP),'m*')
                                hold(debugPlots.apPoints,'off')
                            else
                                plot(endP,apfilt(endP),'m*');
                            end
                        end

                        samps = thirdminL - firstminV;
                        time = samps/chanUsed.SampleRate;
                        startPt = firstminV;
                        endPt = endP;
                        maxHyperpol = thirdminL;%minPt = thirdminL;
                        peakPolarization1 = maxpt1;
                        peakPolarization2 = maxpt2;
                    case 3
                        apfilt = smooth(ap,15,'sgolay');
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,ap)
                                hold(debugPlots.apPoints,'on')
                                plot(debugPlots.apPoints,apfilt)
                            else
                                figure
                                plot(ap)
                                hold on
                                plot(apfilt)
                            end
                        end
                        [a,b]=islocalmax(apfilt);
                        [B,I] = max(b);
                        maxpt1 = I;
                        maxpt2 = [];
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,maxpt1,apfilt(maxpt1),'g*')
                            else
                                plot(maxpt1,apfilt(maxpt1),'g*'); %max pt of AP
                            end
                        end
                        [t, p] = islocalmin(apfilt(1:maxpt1));
                        [firstminV, firstminL] = max(find(t));
                        sigmean = mean(ap)*10;
                        if (maxpt1-firstminV) < 15
                            t(firstminV) = 0;
                            [firstminV, firstminL] = max(find(t));
                        end
                        if apfilt(firstminV) > sigmean
                            realLim = find(apfilt(1:maxpt1)<sigmean,1,'last');
                            if ~isempty(realLim)
                                [t, p] = islocalmin(apfilt(1:realLim));
                                [firstminV, firstminL] = max(find(t));
                            end
                        end
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,firstminV,apfilt(firstminV),'r*')
                            else
                                plot(firstminV,apfilt(firstminV),'r*'); %First corner pt (start of AP)
                            end
                        end

                        [t1, p1] = islocalmin(apfilt(maxpt1:end));
                        [secondminV, secondminL] = max(p1);
                        secondminL = secondminL+maxpt1-1;
                        if secondminL > 300
                            [t1, p1] = islocalmin(apfilt(maxpt1:300));
                            [secondminV, secondminL] = max(p1);
                            secondminL = secondminL+maxpt1-1;
                        end
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,secondminL,apfilt(secondminL),'k*')
                            else
                                plot(secondminL,apfilt(secondminL),'k*'); %Min after peak of AP (Hyperpolarization)
                            end
                        end

                        endP = find(apfilt(secondminL:end)>mean(apfilt)/5,1,'first')+secondminL;
                        if endP - secondminL < 5
                            [aa,aab]=min(abs(apfilt(secondminL+5:400)-mean(apfilt(1:600))));
                            endP = aab+secondminL+5-1;
                        elseif isempty(endP)
                            endP = secondminL+30;
                        end
                        if debugPlot
                            if ~isempty(debugPlots)
                                plot(debugPlots.apPoints,endP,apfilt(endP),'m*')
                                hold(debugPlots.apPoints,'off')
                            else
                                plot(endP,apfilt(endP),'m*'); % end point of recovery period of AP
                            end
                        end

                        samps = secondminL - firstminV;
                        time = samps/chanUsed.SampleRate;
                        startPt = firstminV;
                        endPt = endP;
                        maxHyperpol = secondminL;%minPt = secondminL;
                        peakPolarization1 = maxpt1;
                        peakPolarization2 = [];
                end
            end

        end
        
        function findAPSAfterSubtraction(obj, idxss, chan, app, debugFlg, findFlg)
            warning('off','signal:findpeaks:largeMinPeakHeight')
            if isempty(app)
                D = obj.DAT;
            else
                D = app.DAT;
                app.ProgressBar.Title.String = 'Marking Action Potentials';
                app.PBarObj.Position(3) = 0*1000;
                app.PBarTxt.String = [num2str(round(0*100)),'%'];
                aa = app.FigPlot;
                app.APLocPlot = aa;
                drawnow
            end
            %%ss = smooth(DAT,501,'sgolay',101);
            for idxs = idxss
                sz = size(D.StimBlocks(idxs).([chan,'_ALLSTIM']));
                APLoc = zeros(sz);
                APCount = 0;
                uu = D.StimBlocks(idxs).([chan,'_ZAPPED']);
                apSamps = D.StimBlocks(idxs).([chan,'_APTemplate']).apSamps;
                startPt = D.StimBlocks(idxs).([chan,'_APTemplate']).startPt;
                endPt = D.StimBlocks(idxs).([chan,'_APTemplate']).endPt;
                maxPt1 = D.StimBlocks(idxs).([chan,'_APTemplate']).peakDePol1;
                maxPt2 = D.StimBlocks(idxs).([chan,'_APTemplate']).peakDePol2;
                minPt = D.StimBlocks(idxs).([chan,'_APTemplate']).maxHyperpol;
                ds = [];
                fullap = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt+10:minPt) - mean(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW);
                fAP = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt)';
                
                zappedDataY = D.StimBlocks(idxs).([chan,'_ZAPPED']);
                zappedDataY = reshape(zappedDataY',1,198000);
                [istart,istop,mindist] = findsignal(zappedDataY,fullap,'Normalization','center','NormalizationLength',100);
                
                if isempty(app)
                    
                else
                    app.ProgressBar.Title.String = ['Marking: ','Stim ',num2str(D.StimBlocks(idxs).Stim_E),' Ref ',num2str(D.StimBlocks(idxs).Ref_E),' Current ',num2str(D.StimBlocks(idxs).Current_uA),' Block Number ',num2str(idxs)];
                    drawnow
                end
                
                if findFlg
                    figure
                    findsignal(zappedDataY,fullap,'Normalization','center','NormalizationLength',100)
                    aaga = gcf;
                    aaga.Position = [1          41        1920         963];
                end

                uur = reshape(uu',1,198000);
                shiftV = mean(uur);
                uur = uur-shiftV;
                uu = reshape(uur,1000,198)';
                mm = movmax(uur,9000)*.65;
                mmr = reshape(mm,1000,198)';
                sz = size(uu);
                if fAP(1)>0.1
                    fAP = fAP-fAP(1);
                end
                mAP = max(fAP);
                if mAP > mean(movmax(uur,9000))*1.5
                    mAP = mAP*.7;
                end
                skipNext = 0;
                for ii = 1:sz(1)
                    if ~skipNext
                        if any(uu(ii,:)>mmr(ii,:))
                            ids2u = find(uu(ii,:)>mmr(ii,:));
                            [mVal,mLoc]=max(uu(ii,ids2u));
                            if mVal>mAP*.6
                                if ids2u(mLoc)>950
                                    if ii+1<sz(1)
                                        un = uu(ii+1,1:150);
                                        un2 = [uu(ii,:) un];
                                        [mVal,mLoc]=max(un2(950:end));
                                        mLoc = mLoc+950-1;
                                        if mLoc <= 1000
                                            APLoc(ii,mLoc) = 1;
                                            APCount = APCount + 1;
                                            skipNext = 1;
                                        end
                                    end
                                elseif ids2u(mLoc)<50
                                    if length(ids2u) >10
                                        id = ids2u(mLoc);
                                        APLoc(ii,id) = 1;
                                        APCount = APCount + 1;
                                    end
                                else
                                    id = ids2u(mLoc);
                                    if mVal < 0.1
                                        y=lowpass(uu(ii,:),200,200000);
                                        if y(id)>mAP*.6
                                            APLoc(ii,id) = 1;
                                            APCount = APCount + 1;
                                        end
                                    else
                                        APLoc(ii,id) = 1;
                                        APCount = APCount + 1;
                                    end
                                end
                                
                            end
                        end
                    else
                        skipNext = 0;
                    end
                end

                if debugFlg
                    figure
                    ALLSTIMPlot = reshape(D.StimBlocks(idxs).([chan,'_ALLSTIM'])',1,198000);
                    ZO = reshape(D.StimBlocks(idxs).([chan,'_ZAPPED'])',1,198000);
                    apl = reshape(D.StimBlocks(idxs).([chan,'_APLoc'])',1,198000);
                    plot(ALLSTIMPlot)
                    hold on
                    plot(ZO)
                    plot(find(apl),ZO(find(apl)),'g*')
                    abc = gcf;
                    abc.Position = [1          41        1920         963];
                    hold off
                    name = ['Stim ',num2str(D.StimBlocks(idxs).Stim_E),' Ref ',num2str(D.StimBlocks(idxs).Ref_E),' Current ',num2str(D.StimBlocks(idxs).Current_uA),' Block Number ',num2str(idxs)];
                    title(name);
                    uiwait
                    
                end
                if isempty(app)
                    D.StimBlocks(idxs).([chan,'_APLoc']) = APLoc;
                    D.StimBlocks(idxs).([chan,'_APCount']) = APCount;
                    obj.DAT = D;
                else
                    app.tempAPLoc = APLoc;
                    app.tempAPCount = APCount;
                    app.CurrID = idxs;
                    app.CurrChannel = chan;
                    ALLSTIMPlot = reshape(D.StimBlocks(idxs).([chan,'_ALLSTIM'])',1,198000);
                    ZO = reshape(D.StimBlocks(idxs).([chan,'_ZAPPED'])',1,198000);
                    apl = reshape(app.tempAPLoc',1,198000);
%                     aa = axes(app.AnalysisPlotsTab);
                    plot(aa,ALLSTIMPlot)
                    hold(aa,'on')
                    plot(aa,ZO)
                    plot(aa,find(apl),ZO(find(apl)),'g*')
                    hold(aa,'off')
                    aa.XLim = [-5000 50000]; %200000
                    aa.YLim = [mean(ZO(find(apl)))*-.8 mean(ZO(find(apl)))*1.3];
                    aa.Title.String = [num2str(APCount),' Total APs'];
%                     name = ['Stim ',num2str(D.StimBlocks(idxs).Stim_E),' Ref ',num2str(D.StimBlocks(idxs).Ref_E),' Current ',num2str(D.StimBlocks(idxs).Current_uA),' Block Number ',num2str(idxs)];
%                     aa.Title.String(name);
                    app.figure1.WindowButtonDownFcn = @(src,event)app.APPosCheck;
                    app.figure1.WindowKeyPressFcn = @(src,event)app.KeyDownFcn;
                    switch app.PlotRawDataSwitch.Value
                        case 'On'
                            aa.Children(3).Visible = 'on';
                        case 'Off'
                            aa.Children(3).Visible = 'off';
                    end
                    uiwait
                    if app.stopButton.Value
                        break
                    end
                    app.KeyPressed = 0;
                    D.StimBlocks(idxs).([chan,'_APLoc']) = app.tempAPLoc;
                    D.StimBlocks(idxs).([chan,'_APCount']) = app.tempAPCount;
                    obj.DAT = D;
                    app.figure1.WindowButtonDownFcn = [];
                    app.figure1.WindowKeyPressFcn = [];
                    p = find(idxs==idxss)/length(idxss);
                    app.PBarObj.Position(3) = p*1000;
                    app.PBarTxt.String = [num2str(round(p*100)),'%'];
                    cla(app.FigPlot)
                    drawnow
                end
                
            end
            if isempty(app)
                
            else
                app.figure1.WindowButtonDownFcn = [];
                cla(app.FigPlot)
                app.ProgressBar.Title.String = 'Finished';
                app.PBarObj.Position(3) = 1*1000;
                app.PBarTxt.String = ['Waiting'];
                drawnow
            end
            
        end
        %                 zdy = zappedDataY;
%                 
%                 isVs = [];
%                 istVs = [];
%                 inds = [];
%                 dbl = 0;
%                 strat = 0;
%                 smallTs = 0;
%                 resid = [];
%                 while (mindist < 0.15) && (~dbl)
%                     indd = fix(istart/1000);
%                     indd2 = fix(istop/1000);
%                     opt = 0;
%                     indd2Use = indd;
%                     if indd == 3
%                         db = 1;
%                     end
%                     if indd ~= indd2
%                         if any(ismember(inds,indd+1)) || any(ismember(inds,indd2+1))
%                             indd2Use = inds(find(ismember(inds,[indd+1 indd2+1])))-1;
%                         else
%                             [maxv,maxl]=max(zdy(istart:istop));
%                             maxll = maxl+istart-(indd*1000)-1;
%                             if maxll > 1000
%                                 opt = 2;
%                                 indd2Use = indd2;
%                             else
%                                 opt = 1;
%                                 indd2Use = indd;
%                             end
%                         end
%                     end
%                     subtt = mean(zappedDataY((indd*1000)+1:((indd*1000)+1000)));
%                     tf=11;
%                     u = zdy(istart:istop)-subtt;
%                     filt1=filtfilt(ones(1,tf)/tf,1,u);
%                     [istart2,istop2,dist2] = findsignal(filt1,fullap);
%                     [maxv2,maxl2]=max(filt1);
%                     [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05);%,'NPeaks',1);
%                     wT = max(wT);
%                     if any(ismember(inds,indd2Use+1))
%                         dbl = 1;
%                     elseif smallTs == 3
%                         dbl = 1;
%                     elseif ((max(zdy(istart:istop))-subtt)/max(fullap) < 0.65) || (maxv2/max(fullap) < 0.65)
%                         switch opt
%                             case 0
%                                 bds = ((indd*1000)+1):(((indd+1)*1000));
%                                 zdy(bds) = zdy(bds(1));
%                             case 1
%                                 bds = ((indd*1000)+1):(((indd2+1)*1000));
%                                 zdy(bds) = zdy(bds(1));
%                             case 2
%                                 bds = ((indd*1000)+1):(((indd2+1)*1000));
%                                 zdy(bds) = zdy(bds(1));
%                         end
%                         smallTs = smallTs +1;
%                     else
%                         smallTs = 0;
%                         switch opt
%                             case 0
%                                 inds = [inds indd+1];
%                                 isVs = [isVs istart-(indd*1000)];
%                                 istVs = [isVs istop-(indd*1000)];
%                                 [maxv,maxl]=max(zdy(istart:istop));
%                                 maxll = maxl+istart-(indd*1000)-1;
%                                 APLoc(indd+1,maxll) = 1;
%                                 APCount = APCount + 1;
%                                 bds = ((indd*1000)+1):(((indd+1)*1000));
%                                 msu = mean(zdy(bds));
%                                 resid = [resid; zdy(bds(1:85))-msu];
%                                 zdy(bds) = zdy(bds(1));
%                                 
%                             case 1
%                                 inds = [inds indd+1];
%                                 isVs = [isVs istart-(indd*1000)];
%                                 istVs = [isVs istop-(indd*1000)];
%                                 [maxv,maxl]=max(zdy(istart:istop));
%                                 maxll = maxl+istart-(indd*1000)-1;
%                                 APLoc(indd+1,maxll) = 1;
%                                 APCount = APCount + 1;
%                                 bds = ((indd*1000)+1):(((indd2+1)*1000));
%                                 resid = [resid; zdy(bds(1:85))];
%                                 zdy(bds) = zdy(bds(1));
%                             case 2
%                                 inds = [inds indd2+1];
%                                 isVs = [isVs istart-(indd*1000)];
%                                 istVs = [isVs istop-(indd*1000)];
%                                 [maxv,maxl]=max(zdy(istart:istop));
%                                 maxll = maxl+istart-(indd*1000)-1000;
%                                 APLoc(indd2+1,maxll) = 1;
%                                 APCount = APCount + 1;
%                                 bds = ((indd*1000)+1):(((indd2+1)*1000));
%                                 resid = [resid; zdy(bds(1:85))];
%                                 zdy(bds) = zdy(bds(1));
%                         end
%                     end
%                     if findFlg
%                         close(aaga)
%                     end
%                     [istart,istop,mindist] = findsignal(zdy,fullap,'Normalization','center','NormalizationLength',100);
%                     if findFlg
%                         figure
%                         findsignal(zdy,fullap,'Normalization','center','NormalizationLength',100)
%                         aaga = gcf;
%                         aaga.Position = [1          41        1920         963];
%                     end
%                 end
%                 if findFlg
%                     close(aaga)
%                 end
%                 bestMatch = zappedDataY(istart:istop);
%                 
%                 
%                 if (D.StimBlocks(idxs-1).([chan,'_APCount']) < 198) || (D.StimBlocks(idxs).Current_uA == 20)
%                     for pos = 1:sz(1)
%                         if pos == 29
%                             tttth = 1;
%                         end
%                         ynflg = 0;
%                         if ~any(APLoc(pos,:))
%                             u = uu(pos,:);
%                             
%                             [istart,istop,dist] = findsignal(u,fullap);
% %                             [v, w] = unique( u, 'stable' );
% %                             duplicate_indices = setdiff( 1:numel(u), w );
% %                             if any(duplicate_indices < 60)
% %                                 toTest = duplicate_indices(duplicate_indices < 60);
% %                                 uvs = unique(u(toTest),'stable');
% %                                 minTab = [];
% %                                 maxTab = [];
% %                                 for ijk = 1:length(uvs)
% %                                     if length(find(u==uvs(ijk)))>10
% %                                         uthg = find(u==uvs(ijk));
% %                                         minTab = [minTab min(uthg)];
% %                                         maxTab = [maxTab max(uthg)];
% %                                     end
% %                                 end
% %                                 testTab = 1;
% %                             else
% %                                 testTab = 0;
% %                             end
%                             szs = size(resid);
%                             cDs = [];
%                             iss = [];
%                             isst = [];
%                             for jj = 1:szs(1)
%                                 [istartDi,istopDi,dist] = findsignal(u(1:150)-u(1),resid(jj,:)-resid(jj,1));
%                                 cDs = [cDs dist];
%                                 iss = [iss istartDi];
%                                 isst = [isst istopDi];
%                             end
%                             [residDis, residL] = min(cDs);
%                             if ~isempty(residL)
%                                 if max(resid(residL,:)-resid(residL,1)) > 0.5
%                                     u(iss(residL):isst(residL)) = u(iss(residL):isst(residL))-resid(residL,:)-resid(residL,1);
%                                 end
%                             end
%                             mean(u);
%                             mean(u(200:250));
%                             if abs(mean(u)) < abs(mean(u(200:250)))
%                                 subF = mean(u);
%                             elseif any(ismember(istart:istop,200:250))
%                                 subF = mean(u);
%                             else
%                                 subF = mean(u(200:250));
%                             end
%                             u = u-subF;
%                             %tf=11;
%                             %u=filtfilt(ones(1,tf)/tf,1,u);
%                             [pksT,locsT,wT,pT] = findpeaks(u,'MinPeakHeight',.05,'NPeaks',1);
%                             [pks,locs,w,p] = findpeaks(fullap,'MinPeakHeight',.05,'NPeaks',1);
%                             [istart,istop,dist] = findsignal(u,fullap);
%                             if istart>150
%                                 [pksT,locsT,wT,pT] = findpeaks(u(150:end),'MinPeakHeight',.05,'NPeaks',1);
%                             end
%                             [maxv,maxl]=max(u);
%                             if findFlg
%                                 figure
%                                 findsignal(u,fullap);
%                                 abc = gcf;
%                                 abc.Position = [1          41        1920         963];
%                                 
%                                 
%                                 hold on
%                                 plot(maxl,u(maxl),'r*')
%                                 hold off
%                                 s = 100;
%                                 pos
%                             end
%                             if pos == 182
%                                 sttttt= 1;
%                             end
%                             if (maxl > istart) && (maxl < istop)
%                                 if maxv/max(fullap) > 0.7 %dist<0.05
%                                     if maxl > 1000-length(fullap)
%                                         if pos < 198
%                                             u2 = uu(pos+1,1:150);
%                                             u2 = u2-subF;%mean(u2);
%                                             %u2=filtfilt(ones(1,tf)/tf,1,u2);
%                                             newU = [u u2];
%                                             [istart,istop,dist] = findsignal(newU,fullap);
%                                             [maxv,maxl]=max(newU);
%                                         else
%                                             [Vals_out, Stim_ticks_out, Zero_tick] = D.DataBetweenTicks(D.StimBlocks(idxs).End_tick, D.StimBlocks(idxs).End_tick+D.TicksPerSample*1000+1000);
%                                             Vals_out = Vals_out'-Vals_out(1);
%                                             newU = [u Vals_out(1001:1150)];
%                                             [istart,istop,dist] = findsignal(newU,fullap);
%                                             [maxv,maxl]=max(newU);
%                                         end
%                                         if findFlg
%                                             close(abc)
%                                             figure
%                                             findsignal(newU,fullap);
%                                             abc = gcf;
%                                             abc.Position = [1          41        1920         963];
%                                             
%                                             hold on
%                                             plot(maxl,newU(maxl),'r*')
%                                             hold off
%                                         end
%                                         if (maxl <= 1000) && (dist <0.1)
%                                             APLoc(pos,maxl) = 1;
%                                             APCount = APCount + 1;
%                                             ds = [ds dist];
%                                             ynflg = 1;
%                                         end
%                                     elseif (maxl < length(fullap)) && (pos >1)
%                                         u2 = uu(pos-1,851:1000);
%                                         u2 = u2-subF;%mean(u2);
%                                         %u2=filtfilt(ones(1,tf)/tf,1,u2);
%                                         newU = [u2 u];
%                                         
%                                         [istart,istop,dist] = findsignal(newU,fullap);
%                                         [maxv,maxl]=max(newU);
%                                         dnu = [diff(newU(istart:istop)) mean(newU(istart:istop))/2];
%                                         warning('off','signal:findpeaks:largeMinPeakHeight')
%                                         [pksV,stpks,~,p] = findpeaks(-dnu,'MinPeakHeight',0.05);
%                                         [pksV,endpks,~,p] = findpeaks(dnu,'MinPeakHeight',0.05);
%                                         if ~isempty(stpks) && ~isempty(endpks)%(length(flatMin)>1) || (length(flatMax)>1)
%                                             if length(stpks) == length(endpks)
%                                                 bbnds = [];
%                                                 for ij = 1:length(stpks)
%                                                     bbnds = [bbnds (stpks(ij)+1):endpks(ij)];
%                                                 end
%                                                 time = round(0:obj.DAT.SampleInterval_sec:(length(newU)-1)*obj.DAT.SampleInterval_sec,7); %normalized time of the artifact template
%                                                 timeTrim = time;
%                                                 timeTrim(bbnds+istart-1) = []; %Blanking time data
%                                                 newU(bbnds+istart-1) = []; %Blanking artifact template
%                                                 [curve, goodness, output] = fit(timeTrim',newU','linearinterp'); %Fitting a spline to the smoothed and blanked artifact template
%                                                 x = time;
%                                                 y = curve(time)';
%                                                 newU = y;
%                                                 [istart,istop,dist] = findsignal(newU,fullap);
%                                                 [maxv,maxl]=max(newU);
%                                             end
%                                         end
%                                         %                                 endpks = find(dnu>.06);
%                                         %                                 stpks = find(dnu<-.06);
%                                         %                                 [maxWIv,maxWIl]=max(newU(istart:istop));
%                                         %                                 [minWIv,minWIl]=min(newU(istart:istop));
%                                         %                                 flatMax = find(maxWIv==newU(istart:istop));
%                                         %                                 flatMin = find(minWIv==newU(istart:istop));
%                                         if findFlg
%                                             close(abc)
%                                             figure
%                                             findsignal(newU,fullap);
%                                             abc = gcf;
%                                             abc.Position = [1          41        1920         963];
%                                             
%                                             
%                                             hold on
%                                             plot(maxl,newU(maxl),'r*')
%                                             hold off
%                                         end
%                                         if (maxl > 150) && (dist < 0.5)
%                                             if dist < 0.075
%                                                 ds = [ds dist];
%                                                 APLoc(pos,maxl-150) = 1;
%                                                 APCount = APCount + 1;
%                                                 ynflg = 1;
%                                             else
%                                                 tf=11;
%                                                 filt1=filtfilt(ones(1,tf)/tf,1,newU);%%%%%%%%%%
%                                                 [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05,'NPeaks',1);
%                                                 [pks,locs,w,p] = findpeaks(fullap,'MinPeakHeight',.05,'NPeaks',1);
%                                                 if ~isempty(locsT)
%                                                     if locsT < istart || locsT > istop
%                                                         [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05);
%                                                         [wT,locsTT] = max(wT);
%                                                         locsT = locsT(locsTT);
%                                                     end
%                                                 end
%                                                 if wT < 10
%                                                     a = find([diff(abs(flip(filt1(1:locsT)))) 0]>0,1,'first');
%                                                     [dd,ddf]=min(abs(filt1(locsT+20:locsT+90)));
%                                                     wT = ddf+a;
%                                                 end
%                                                 if wT/w > .85
%                                                     if istart-20 < 0
%                                                         minvs = min(newU(1:istop+20));
%                                                     else
%                                                         minvs = min(newU(istart-20:istop+20));
%                                                     end
%                                                     if abs(minvs/max(fullap)) < 0.5 && ~ismember({chan},'NerveCon')
% %                                                         if testTab
% %                                                             if (minTab(1)<=maxl-150) && (maxl-150<=maxTab(1))
% %                                                                 ynflg = 0;
% %                                                             else
% %                                                                 ds = [ds dist];
% %                                                             APLoc(pos,maxl-150) = 1;
% %                                                             APCount = APCount + 1;
% %                                                             ynflg = 1;
% %                                                             end
% %                                                         else
%                                                             ds = [ds dist];
%                                                             APLoc(pos,maxl-150) = 1;
%                                                             APCount = APCount + 1;
%                                                             ynflg = 1;
% %                                                         end
%                                                     else
%                                                         ds = [ds dist];
%                                                         APLoc(pos,maxl-150) = 1;
%                                                         APCount = APCount + 1;
%                                                         ynflg = 1;
%                                                     end
%                                                 end
%                                             end
%                                         elseif ~isempty(stpks) && ~isempty(endpks)%(length(flatMin)>1) || (length(flatMax)>1)
%                                             %if length(flatMin>1)
%                                             if length(stpks) == length(endpks)
%                                                 bbnds = [];
%                                                 for ij = 1:length(stpks)
%                                                     bbnds = [bbnds (stpks(ij)+1):endpks(ij)];
%                                                 end
%                                                 time = round(0:obj.DAT.SampleInterval_sec:(length(newU)-1)*obj.DAT.SampleInterval_sec,7); %normalized time of the artifact template
%                                                 timeTrim = time;
%                                                 timeTrim(bbnds+istart-1) = []; %Blanking time data
%                                                 newU(bbnds+istart-1) = []; %Blanking artifact template
%                                                 
%                                                 [curve, goodness, output] = fit(timeTrim',newU','linearinterp'); %Fitting a spline to the smoothed and blanked artifact template
%                                                 x = time;
%                                                 y = curve(time)';
%                                                 %                                         tf = 5;
%                                                 %                                         filt1=filtfilt(ones(1,tf)/tf,1,newU);
%                                                 [is,is,dist] = findsignal(y,fullap);
%                                             end
%                                             if (maxl > 150) && (dist < 0.075)
%                                                 ds = [ds dist];
%                                                 APLoc(pos,maxl-150) = 1;
%                                                 APCount = APCount + 1;
%                                                 ynflg = 1;
%                                             end
%                                             
%                                             %end
%                                         end
%                                     elseif dist < 0.35
%                                         if (maxl < 50) && (wT/w < .95)
%                                             ynflg = 0;
%                                         else
%                                             APLoc(pos,maxl) = 1;
%                                             APCount = APCount + 1;
%                                             ds = [ds dist];
%                                             ynflg = 1;
%                                         end
%                                     end
%                                     
%                                 elseif dist < 0.05
%                                     APLoc(pos,maxl) = 1;
%                                     APCount = APCount + 1;
%                                     ds = [ds dist];
%                                     ynflg = 1;
%                                 elseif abs(u(istop)-fullap(end)) > 0.02
%                                     %                                     u = u-(u(istop)-fullap(end));
%                                     u = u-(mean(u(1:50)));
%                                     [istart,istop,dist] = findsignal(u,fullap);
%                                     [maxv,maxl]=max(u);
%                                     if dist < 0.1
%                                         APLoc(pos,maxl) = 1;
%                                         APCount = APCount + 1;
%                                         ds = [ds dist];
%                                         ynflg = 1;
%                                     end
%                                 elseif maxv/max(fullap) > 0.6
%                                     uuu = u;
%                                     uuu(istart:istop) = u(istart);
%                                     tf=11;
%                                     filt1=filtfilt(ones(1,tf)/tf,1,uuu);
%                                     if max(filt1)/max(fullap) < 0.2
% %                                         if testTab
% %                                             if (minTab(1)<=maxl) && (maxl<=maxTab(1))
% %                                                 ynflg = 0;
% %                                             else
% %                                                 ds = [ds dist];
% %                                                 APLoc(pos,maxl) = 1;
% %                                                 APCount = APCount + 1;
% %                                                 ynflg = 1;
% %                                             end
% %                                         else
%                                             ds = [ds dist];
%                                             APLoc(pos,maxl) = 1;
%                                             APCount = APCount + 1;
%                                             ynflg = 1;
% %                                         end
%                                     end
%                                 elseif (maxv-mean(u(1:50)))/max(fullap) > 0.65
%                                     uuu = u;
%                                     uuu(istart:istop) = u(istart);
%                                     if max(uuu)/max(fullap) < 0.2
%                                         APLoc(pos,maxl) = 1;
%                                         APCount = APCount + 1;
%                                         ds = [ds dist];
%                                         ynflg = 1;
%                                     end
%                                 end
%                             else
%                                 
%                                 [maxv2,maxl2]=max(u(istart:istop));
%                                 if (maxv2/max(fullap) > 0.7) && (dist < 0.5) && (maxl2>10)
%                                     [pks,locs,w,p] = findpeaks(fullap,'MinPeakHeight',.05,'NPeaks',1);
%                                     if istart-10 < 1
%                                         [pksT,locsT,wT,pT] = findpeaks(u(istart:end),'MinPeakHeight',.05,'NPeaks',1);
%                                     else
%                                         [pksT,locsT,wT,pT] = findpeaks(u(istart-10:end),'MinPeakHeight',.05,'NPeaks',1);
%                                         
%                                     end
%                                     if wT < 5
%                                         tf=11;
%                                         filt1=filtfilt(ones(1,tf)/tf,1,u);
%                                         [istart2,istop2,dist2] = findsignal(filt1,fullap);
%                                         [maxv2,maxl2]=max(filt1);
%                                         [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05);%,'NPeaks',1);
%                                         wT = max(wT);
%                                     end
%                                     if (round(wT/w,3) >=.7)
%                                         if (dist > 0.25) && (istart < 85) && (residDis < .15)
%                                             ynflg = 0;
%                                         else
%                                             flatMax = find(maxv==u(1:100));
%                                             if flatMax > 1
%                                                 if istart-20 < 1
%                                                     ttu1 = 1;
%                                                     dd3 = 0;
%                                                 else
%                                                     ttu1 = istart-20;
%                                                     dd3 = istart-21;
%                                                 end
%                                                 [maxvn,maxln]=max(u(ttu1:istop));
% %                                                 if testTab
% %                                                     if (minTab(1)<=maxln+dd3) && (maxln+dd3<=maxTab(1))
% %                                                         ynflg = -1;
% %                                                     else
% %                                                         APLoc(pos,maxln+dd3) = 1;
% %                                                     end
% %                                                 else
%                                                     APLoc(pos,maxln+dd3) = 1;
% %                                                 end
%                                             else
% %                                                 if testTab
% %                                                     if (minTab(1)<=maxl) && (maxl<=maxTab(1))
% %                                                         ynflg = -1;
% %                                                     else
% %                                                         APLoc(pos,maxl) = 1;
% %                                                     end
% %                                                 else
%                                                     APLoc(pos,maxl) = 1;
% %                                                 end
%                                             end
%                                             if ynflg < 0
%                                                 ynflg =0;
%                                             else
%                                                 APCount = APCount + 1;
%                                                 ds = [ds dist];
%                                                 ynflg = 1;
%                                             end
%                                         end
%                                     elseif istart > (mean(isVs)-std(isVs)-2)
%                                         APLoc(pos,maxl2+istart-1) = 1;
%                                         APCount = APCount + 1;
%                                         ds = [ds dist];
%                                         ynflg = 1;
%                                     elseif (maxl > 1000-length(fullap)) && (maxv/max(fullap) > 0.7)
%                                         u2 = uu(pos+1,1:150);
%                                         u2 = u2-subF;
%                                         newU = [u u2];
%                                         [istart3,istop3,dist3] = findsignal(newU,fullap);
%                                         [maxv3,maxl3]=max(newU);
%                                         if (maxl3 <= 1000) && (dist3 < 0.075)
%                                             APLoc(pos,maxl3) = 1;
%                                             APCount = APCount + 1;
%                                             ds = [ds dist3];
%                                             ynflg = 1;
%                                         elseif maxl3 <= 1000
%                                             if (maxv3/max(fullap) > 0.7) && (dist3 < 0.5)
%                                                 ds = [ds dist3];
%                                                 APLoc(pos,maxl3) = 1;
%                                                 APCount = APCount + 1;
%                                                 ynflg = 1;
%                                             end
%                                         end
%                                     end
%                                 elseif (maxl > 1000-length(fullap)) && (maxv/max(fullap) > 0.7)
%                                     if pos == 198
%                                         vss = D.ValsBetweenTimes(D.StimBlocks(idxs).Start_sec,D.StimBlocks(idxs).End_sec+1,chan);
%                                         u2 = vss(198001:198150)';
%                                     else
%                                         u2 = uu(pos+1,1:150);
%                                     end
%                                     u2 = u2-subF;%mean(u2);
%                                     %u2=filtfilt(ones(1,tf)/tf,1,u2);
%                                     newU = [u u2];
%                                     [istart2,istop2,dist2] = findsignal(newU,fullap);
%                                     [maxv2,maxl2]=max(newU);
%                                     
%                                     if findFlg
%                                         close(abc)
%                                         figure
%                                         findsignal(newU,fullap);
%                                         abc = gcf;
%                                         abc.Position = [1          41        1920         963];
%                                         
%                                         hold on
%                                         plot(maxl2,newU(maxl2),'r*')
%                                         hold off
%                                     end
%                                     if (maxl2 <= 1000) && (dist2 < 0.075)
%                                         APLoc(pos,maxl2) = 1;
%                                         APCount = APCount + 1;
%                                         ds = [ds dist2];
%                                         ynflg = 1;
%                                     elseif maxl2 <= 1000
%                                         tf=11;
%                                         filt1=newU;%filtfilt(ones(1,tf)/tf,1,newU);
%                                         [istart3,istop3,dist3] = findsignal(filt1,fullap);
%                                         [maxv3,maxl3]=max(filt1);
%                                         if (maxv3/max(fullap) > 0.7) && (dist3 < 0.5)
%                                             ds = [ds dist3];
%                                             APLoc(pos,maxl2) = 1;
%                                             APCount = APCount + 1;
%                                             ynflg = 1;
%                                         end
%                                     end
%                                 elseif (maxl < length(fullap)) && (maxv/max(fullap) > 0.7) && (pos>1)
%                                     u2 = uu(pos-1,851:1000);
%                                     u2 = u2-subF;%mean(uu(pos-1,:));
%                                     %                                     tf = 11;
%                                     %                                     u2=filtfilt(ones(1,tf)/tf,1,u2);
%                                     newU = [u2 u];
% %                                     if testTab
% %                                         twe = 1:length(newU);
% %                                         ff = fit(twe',newU','smoothingspline','exclude',[minTab(1):maxTab(1)+1]+150);
% %                                         newU = ff(twe);
% %                                     end
%                                     [istart2,istop2,dist2] = findsignal(newU,fullap);
%                                     [maxv2,maxl2]=max(newU);
%                                     if findFlg
%                                         close(abc)
%                                         figure
%                                         findsignal(newU,fullap);
%                                         abc = gcf;
%                                         abc.Position = [1          41        1920         963];
%                                         
%                                         
%                                         hold on
%                                         plot(maxl2,newU(maxl2),'r*')
%                                         hold off
%                                     end
%                                     if (maxl2 > 150) && (dist2 < 0.075)
%                                         ds = [ds dist2];
%                                         APLoc(pos,maxl2-150) = 1;
%                                         APCount = APCount + 1;
%                                         ynflg = 1;
%                                     elseif maxl2 > 150
%                                         filt1 = u(1:150);
%                                         [istart2,istop2,dist2] = findsignal(filt1,fullap);
%                                         [maxv2,maxl2]=max(filt1);
%                                         [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05,'NPeaks',1);
%                                         [pks,locs,w,p] = findpeaks(fullap,'MinPeakHeight',.05,'NPeaks',1);
%                                         if wT <7
%                                             tf=11;
%                                             filt1=filtfilt(ones(1,tf)/tf,1,u(1:150));
%                                             [istart2,istop2,dist2] = findsignal(filt1,fullap);
%                                             [maxv2,maxl2]=max(filt1);
%                                             [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05);%,'NPeaks',1);
%                                             wT = max(wT);
%                                         end
%                                         if findFlg
%                                             close(abc)
%                                             figure
%                                             findsignal(filt1,fullap);
%                                             abc = gcf;
%                                             abc.Position = [1          41        1920         963];
%                                             
%                                             
%                                             hold on
%                                             plot(maxl2,filt1(maxl2),'r*')
%                                             hold off
%                                         end
%                                         if isempty(dist2)
%                                             ynflg = 0;
%                                         else
%                                             if (maxv2/max(fullap) > 0.7) && (dist2 < 0.5) && (wT/w > 0.6)
%                                                 if (dist > 0.25) && (istart < 85) && (residDis < .15)
%                                                     ynflg = 0;
%                                                 else
%                                                     ds = [ds dist2];
%                                                     APLoc(pos,maxl2) = 1;
%                                                     APCount = APCount + 1;
%                                                     ynflg = 1;
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 elseif (maxv/max(fullap) > 0.7)
%                                     tf=11;
%                                     filt1=filtfilt(ones(1,tf)/tf,1,u);
%                                     [istart2,istop2,dist2] = findsignal(filt1,fullap);
%                                     [maxv2,maxl2]=max(filt1);
%                                     [pksT,locsT,wT,pT] = findpeaks(filt1,'MinPeakHeight',.05,'NPeaks',1);
%                                     [pks,locs,w,p] = findpeaks(fullap,'MinPeakHeight',.05,'NPeaks',1);
%                                     
%                                     if (dist < 0.2)
%                                         ds = [ds dist];
%                                         APLoc(pos,maxl) = 1;
%                                         APCount = APCount + 1;
%                                         ynflg = 1;
%                                     elseif isempty(wT)
%                                         ynflg = 0;
%                                     elseif (dist2 < 0.5) && (wT/w > 0.6)
%                                         ds = [ds dist];
%                                         APLoc(pos,maxl) = 1;
%                                         APCount = APCount + 1;
%                                         ynflg = 1;
%                                     end
%                                 end
%                                 %                         maxl;
%                                 %                         maxv/max(fullap);
%                                 %                         if (maxl < 85) && (maxv/max(fullap) < 0.7)
%                                 %                             %noresponse Maybe add saving this stim
%                                 %
%                                 %                         elseif maxl > 150
%                                 %                             %noresponse
%                                 %                         else
%                                 %                             %response = [response; ALLSTIMPlot(bound)];
%                                 %                             %                             figure
%                                 %                             %                             findsignal(u(1:250),fullap);
%                                 %                             %                             abc = gcf;
%                                 %                             %                             abc.Position = [1          41        1920         963];
%                                 %                             %                             close(abc)
%                                 %                             [istart,istop,dist] = findsignal(u(1:250),fullap);
%                                 %                             u2 = ALLSTIMPlot2(bound);
%                                 %                             fAP = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt)';
%                                 %                             u2(istart:(istart+length(fAP)-1)) = u2(istart:(istart+length(fAP)-1)) - fAP;
%                                 %                             fornextTemp = [fornextTemp; u2];
%                                 %                         end
%                             end
%                             if findFlg
%                                 if ynflg == 1
%                                     'yes'
%                                 else
%                                     'no'
%                                 end
%                                 close(abc)
%                             end
%                         end
%                     end
%                 else
%                     zappedDataY = D.StimBlocks(idxs).([chan,'_ZAPPED']);
%                     DAT = reshape(zappedDataY',1,198000);
%                     tf = 501;
%                     ss = filtfilt(ones(1,tf)/tf,1,DAT);
%                     DAT = DAT-ss;
%                     % SMOOTH the raw head-stage data first...
%                     SMOO = smooth(DAT,21,'sgolay');
%                     
%                     % THEN "differentiate" using the shifted subtraction method.
%                     SS = obj.DoShiftSubtract(SMOO);
%                     
%                     % SQUARE the signal, but PRESERVE the sign.
%                     SS = (SS.*SS) .* sign(SS);
%                     
%                     % This is HARD CODED to look for NEGATIVE PEAKS!!
%                     % NEED TO make this auto-detect the AP polarity.
%                     
%                     % If a specific threshold is not given then use a Moving
%                     % Maximum to set the "threshold" for detecting peaks.
%                     [pks,locs,w,p] = findpeaks(SMOO,1:length(SMOO),'minpeakheight',.05,'minpeakdistance',900);
%                     if length(pks) == 198
%                         sz = size(D.StimBlocks(idxs).([chan,'_ALLSTIM']));
%                         APLoc = zeros(sz);
%                         APCount = 198;
%                         for pos = 1:sz(1)
%                             ls = locs(pos)-1000*(pos-1);
%                             APLoc(pos,ls) = 1;
%                         end
%                     end
%                 end
%                 
%                 
%                 mean(ds);
%                 max(ds);

        function [i1, i2, iPk, d] = findAPs_Brian(obj, zappedDataX, zappedDataY, apTemplate, unzappedY, debugPlotFlag)
            go = 1;
            SMOOc = zappedDataY;
            i1 = [];
            i2 = [];
            iPk = [];
            i1Temp = [];
            i2Temp = [];
            iPkTemp = [];
            d = [];
            corD = [];
            apDiffA = [];
            hDiff = [];
            negAP = 0;
            
            % Put the offset back in!!
            %apTemplate = apTemplate + mean(zappedDataY);
            AVG = apTemplate;
            dblChek = zeros(1,length(zappedDataY));
            while go
                III = length(i2);
                %                 if III == 190
                %                     III
                %
                %                 end
                %fprintf('%d ', III);
                %if mod(III,20) == 0; fprintf('\n'); end
                % %
                if debugPlotFlag
                    findsignal(SMOOc,AVG);%'Normalization','power','Annotate','all')
                    drawnow
                    abc = gcf;
                    hold(abc.Children(2), 'on')
                    bb = plot(abc.Children(2),1:length(SMOOc),zappedDataY,'g');
                    uistack(bb,'bottom')
                    cc = plot(abc.Children(2),1:length(SMOOc),unzappedY,'k');
                    uistack(cc,'bottom')
                    drawnow
                    hold(abc.Children(2),'off')
                end
                
                %[istart,istop,dist]=findsignal(SMOOc,apTemplate);%,'Normalization','power');
                [istart,istop,dist]=findsignal(SMOOc,apTemplate);%,'Normalization','power');
                if debugPlotFlag
                    abc.Children(2).YLim = [min(SMOOc(istart:istop))-0.06 max(SMOOc(istart:istop))+0.08];
                    abc.Children(2).XLim = [istart-2600 istop+2600];
                end
                [~,w] = max(SMOOc(istart:istop));
                w = w+istart-1;
                apDiff = SMOOc(w)-SMOOc(w+20);
                [~,wTemp]=max(apTemplate);
                tempDiff = apTemplate(wTemp)-apTemplate(wTemp+20);
                tempDiff*.5;
                zappedDataY(w)-max(apTemplate);
                [iPkTemp PkO] = sort(iPk);
                i1Temp = i1(PkO);
                i2Temp = i2(PkO);
                startWin = (istart > i1Temp & istart < i2Temp);
                stopWin = (istop > i1Temp & istop < i2Temp);
                [c,lags] = xcorr([0;diff(SMOOc(istart:istop))],[0;diff(apTemplate)],0,'coeff');
                if isempty(i1)
                    i1 = [i1 istart];
                    i2 = [i2 istop];
                    d = [d dist];
                    corD = [corD c];
                    [~,w] = max(SMOOc(istart:istop));
                    w = w+istart-1;
                    iPk = [iPk w];
                    btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                    aD = SMOOc(istart):btwP:SMOOc(istop);
                    SMOOc(istart:istop) = aD;
                    apDiffA = [apDiffA apDiff];
                    hDiff = [hDiff (zappedDataY(w)-max(apTemplate))];
                    
                elseif length(i1)>200
                    go = 0;
                elseif negAP > 3
                    if debugPlotFlag
                        keyboard
                    end
                    mPk = mean(zappedDataY(iPkTemp))*(3/4);
                    if any(dblChek)
                        SMOOc(logical(dblChek)) = zappedDataY(logical(dblChek));
                    end
                    chunks = SMOOc>mPk;
                    stInds = find([0;diff(chunks)]>0);
                    eInds = find([0;diff(chunks)]<0);
                    if ~isempty(stInds)
                        if stInds(end)==198100
                            %                         chunks = 0;
                            stInds(end) = [];
                        end
                    end
                    if ~isempty(stInds)
                        if length(stInds) > length(eInds)
                            if stInds(end)>198050
                                stInds(end) = [];
                            end
                        else
                            if stInds(1) > eInds(1)
                                eInds(1) = [];
                            else
                            end
                        end
                    end
                    
                    if length(stInds) ~= length(eInds)
                        keyboard
                    end
                    if any(chunks)
                        if length(stInds) == length(eInds)
                            for qw = 1:length(stInds)
                                tSMOOc = SMOOc;
                                if stInds(qw)-214<0
                                    b1 = 1;
                                else
                                    tSMOOc(1:stInds(qw)-215) = tSMOOc(stInds(qw)-215);
                                    b1 = stInds(qw)-214;
                                end
                                if eInds(qw)+214>length(SMOOc)
                                    b2 = length(SMOOc);
                                else
                                    tSMOOc(eInds(qw)+215:end) = tSMOOc(eInds(qw)+215);
                                    b2 = eInds(qw)+214;
                                end
                                if debugPlotFlag
                                    findsignal(tSMOOc,AVG);%'Normalization','power','Annotate','all')
                                    drawnow
                                    abc = gcf;
                                    hold(abc.Children(2), 'on')
                                    bb = plot(abc.Children(2),1:length(tSMOOc),zappedDataY,'g');
                                    uistack(bb,'bottom')
                                    cc = plot(abc.Children(2),1:length(SMOOc),unzappedY,'k');
                                    uistack(cc,'bottom')
                                    drawnow
                                    t = 1:length(tSMOOc);
                                    plot(abc.Children(2),t(iPk),zappedDataY(iPk),'r*','MarkerSize',10);
                                    drawnow
                                    hold(abc.Children(2),'off')
                                end
                                [istart,istop,dist]=findsignal(tSMOOc,apTemplate);%,'Normalization','power');
                                if debugPlotFlag
                                    abc.Children(2).YLim = [min(SMOOc(istart:istop))-0.06 max(SMOOc(istart:istop))+0.08];
                                    abc.Children(2).XLim = [istart-2600 istop+2600];
                                    keyboard
                                end
                                
                                [~,w] = max(tSMOOc(istart:istop));
                                w = w+istart-1;
                                apDiff = tSMOOc(w)-tSMOOc(w+20);
                                [~,wTemp]=max(apTemplate);
                                tempDiff = apTemplate(wTemp)-apTemplate(wTemp+20);
                                zappedDataY(w)-max(apTemplate);
                                [iPkTemp PkO] = sort(iPk);
                                i1Temp = i1(PkO);
                                i2Temp = i2(PkO);
                                startWin = (istart > i1Temp & istart < i2Temp);
                                stopWin = (istop > i1Temp & istop < i2Temp);
                                if (istart < b1 | istart > b2)
                                elseif (istop < b1 | istop > b2)
                                elseif (tSMOOc(w)-max(apTemplate))<-max(apTemplate)/5
                                elseif (stInds(qw) < istart | stInds(qw) > istop)
                                elseif (eInds(qw) < istart | eInds(qw) > istop)
                                elseif obj.D.SamplesPerSec/min(abs(iPkTemp-w))>650
                                elseif (tSMOOc(w)-tSMOOc(istart))>(tSMOOc(w)-tSMOOc(istop))
                                elseif dist/d(end) > 3.5
                                    if (istart<2000) && (SMOOc(w)>max(apTemplate))
                                        i1 = [i1 istart];
                                        i2 = [i2 istop];
                                        d = [d dist];
                                        [~,w] = max(SMOOc(istart:istop));
                                        w = w+istart-1;
                                        iPk = [iPk w];
                                        btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                                        aD = SMOOc(istart):btwP:SMOOc(istop);
                                        SMOOc(istart:istop) = aD;
                                        apDiffA = [apDiffA apDiff];
                                        hDiff = [hDiff (zappedDataY(w)-max(apTemplate))];
                                    end
                                else
                                    i1 = [i1 istart];
                                    i2 = [i2 istop];
                                    d = [d dist];
                                    [~,w] = max(SMOOc(istart:istop));
                                    w = w+istart-1;
                                    iPk = [iPk w];
                                    btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                                    aD = SMOOc(istart):btwP:SMOOc(istop);
                                    SMOOc(istart:istop) = aD;
                                    apDiffA = [apDiffA apDiff];
                                    hDiff = [hDiff (zappedDataY(w)-max(apTemplate))];
                                end
                            end
                        end
                    end
                    go = 0;
                elseif dist/d(end) > 2.2
                    if debugPlotFlag
                        keyboard
                    end
                    if (istart<2000) && (SMOOc(w)>max(apTemplate))
                        saveFlg = 1;
                    elseif ((zappedDataY(w)-max(apTemplate))>-max(apTemplate)/6) && ~any(startWin) && ~any(stopWin) && ((mean(corD)/c)<2)
                        
                        saveFlg = 1;
                    else
                        negAP = negAP + 1;
                        btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                        aD = SMOOc(istart):btwP:SMOOc(istop);
                        SMOOc(istart:istop) = aD;
                        saveFlg = 0;
                    end
                    if saveFlg
                        i1 = [i1 istart];
                        i2 = [i2 istop];
                        d = [d dist];
                        corD = [corD c];
                        [~,w] = max(SMOOc(istart:istop));
                        w = w+istart-1;
                        iPk = [iPk w];
                        btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                        aD = SMOOc(istart):btwP:SMOOc(istop);
                        SMOOc(istart:istop) = aD;
                        apDiffA = [apDiffA apDiff];
                        hDiff = [hDiff (zappedDataY(w)-max(apTemplate))];
                        negAP = 0;
                    end
                elseif (zappedDataY(w)-max(apTemplate))<-max(apTemplate)/5.7
                    if debugPlotFlag
                        keyboard
                    end
                    if (zappedDataY(w)-max(apTemplate))>-max(apTemplate)/5
                        dblChek(istart:istop) = 1;
                    end
                    negAP = negAP + 1;
                    btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                    aD = SMOOc(istart):btwP:SMOOc(istop);
                    SMOOc(istart:istop) = aD;
                    
                    %                 elseif apDiff < tempDiff*.5
                    %                     go = 0;
                    %         elseif (dist/d(end))>1.5
                    %             go = 0;
                    % 				elseif any([abs(zappedDataX(iPk)-zappedDataX(w))]<0.0025)
                    % 					btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                    % 					aD = SMOOc(istart):btwP:SMOOc(istop);
                    % 					SMOOc(istart:istop) = aD;
                elseif any(startWin)
                    if debugPlotFlag
                        keyboard
                    end
                    negAP = negAP + 1;
                    btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                    aD = SMOOc(istart):btwP:SMOOc(istop);
                    SMOOc(istart:istop) = aD;
                elseif any(stopWin)
                    if debugPlotFlag
                        keyboard
                    end
                    negAP = negAP + 1;
                    btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                    aD = SMOOc(istart):btwP:SMOOc(istop);
                    SMOOc(istart:istop) = aD;
                else
                    i1 = [i1 istart];
                    i2 = [i2 istop];
                    d = [d dist];
                    corD = [corD c];
                    [~,w] = max(SMOOc(istart:istop));
                    w = w+istart-1;
                    iPk = [iPk w];
                    btwP = (SMOOc(istop)-SMOOc(istart))/(length(istart:istop)-1);
                    aD = SMOOc(istart):btwP:SMOOc(istop);
                    SMOOc(istart:istop) = aD;
                    apDiffA = [apDiffA apDiff];
                    hDiff = [hDiff (zappedDataY(w)-max(apTemplate))];
                    negAP = 0;
                end
                
                % 				         close(gcf)
            end
        end
        
        % Find the APs before and during a stimulation block. If the neuron
        % is "responding" during the stimulation block (far fewer APs left
        % after initial zapping), re-zap the block with a "corrected"
        % zapping template (with the AP removed from the template).
        function [P,S] = FindAPsBlock(obj, blockidx, debugPlotFlag)
            [DPRE,DSTIM,StimIdx] = obj.GetBlockData(blockidx);
            P = obj.FindAPs(DPRE);
            S = obj.FindAPs(DSTIM, -0.01);
            
            % Are there a LOT fewer APs within the stim block, than before
            % the stim block? If so, then we assume our APs also got
            % zapped, and try to subtract the AP template out of the
            % zapping template, and re-zap the data.
            NPRE = length(P.CLUMP);
            NST = length(S.CLUMP);
            
            %if NST < NPRE*0.6   % Not enough APs detected in stim block?
            if true
                fprintf('"Fixing" stim block, adding in AP template...\n');
                
                
                ZT_in = obj.AZT(:, blockidx);
                cu = obj.D.StimBlocks(blockidx).Current_uA;
                % 				[preMid, postMid, timeTrim, sigTrim, xTemplate, yTemplate] = obj.SplineZapTemplate(ZT_in, obj.D, cu);
                ZT_fixed = ZT_in;
                %ZT_fixed = obj.FixZapTemplate(ZT_in, P.AVG);
                
                % Get the unzapped data, and re-zap it.
                [a, b, c] = obj.GetBlockData(blockidx, obj.D.Values);
                [~,DSTIM] = obj.GetBlockData(blockidx, obj.D.Values);
                DSTIMO = DSTIM;
                ApTemp = P.AVG(1:250);
                [maxPV,maxP] = max(ApTemp);
                [minPV,minP] = min(ApTemp);
                if abs(maxPV)>abs(minPV)
                    direct = 1;
                    dAP = [0;diff(P.AVG(1:maxP))];
                    fp = max(find(dAP<0));
                    ApTemp = P.AVG(fp-5:250);
                    zt2 = ZT_fixed;
                    if debugPlotFlag
                        ff = figure;
                        plot(zt2)
                        hold on
                        keyboard
                    end
                    [mzv,mz] = max(ZT_fixed(1:200));
                    chkMax = ZT_fixed(1:200);
                    chkMax(chkMax<mzv-.05) = 0;
                    TF = find(islocalmax(chkMax,'FlatSelection','all'));
                    if length(TF) > 1
                        mz = max(TF);
                        mzv = ZT_fixed(mz);
                        if length(TF) > 3
                            withInFlg = 1;
                        else
                            withInFlg = 0;
                        end
                    else
                        withInFlg = 0;
                    end
                    if ~withInFlg
                        zt3 = zt2;
                        zt3(1:mz) = mzv;
                        if debugPlotFlag
                            plot(zt3)
                        end
                        tf=5;
                        filt1=filtfilt(ones(1,tf)/tf,1,zt3(mz+1:end)); %Smooth part of the orignial artifact template
                        zt3(mz+1:end) = filt1;
                        if debugPlotFlag
                            plot(zt3)
                        end
                        dZT3 = [0;diff(zt3)];
                        dZT3(1:mz) = -mzv;
                        pl = min(find(dZT3(1:200)>-0.0005));
                        if isempty(pl)
                            [~,pl] = max(dZT3(1:150));
                        end
                        if debugPlotFlag
                            plot(pl,zt3(pl),'b*')
                        end
                        zt3(1:pl(1)) = -2;
                        if debugPlotFlag
                            plot(zt3)
                        end
                        [locMaxV, locMax] = max(zt3(1:200));
                        if debugPlotFlag
                            plot(locMax,zt3(locMax),'g*')
                        end
                        zt3(1:locMax) = 2;
                        if debugPlotFlag
                            plot(zt3)
                        end
                        [locMin2V, locMin2] = min(zt3(1:350));
                        if debugPlotFlag
                            plot(locMin2,zt3(locMin2),'r*')
                        end
                        ApTemp2 = P.AVG(fp-5:end-300);
                        [mpv,mp] = max(ApTemp2);
                    else
                        ApTemp2 = P.AVG(fp:end-300);
                        [mpv,mp] = max(ApTemp2);
                        
                        locMaxV = 500;
                        locMax = mz;
                        toShift = min(TF);
                        
                    end
                    if debugPlotFlag
                        hold off
                    end
                else
                    %                     direct = -1;
                    %                     dAP = [0;diff(P.AVG(1:minP))];
                    %                     fp = max(find(dAP>0));
                    %                     ApTemp = P.AVG(fp-5:250);
                    %                     zt2 = ZT_fixed;
                    %                     [mzv,mz] = max(ZT_fixed);
                    %                     chkMax = ZT_fixed;
                    %                     chkMax(chkMax<mzv-.05) = 0;
                    %                     TF = find(islocalmax(chkMax,'FlatSelection','all'));
                    %                     if length(TF) > 1
                    %                         mz = max(TF);
                    %                         mzv = ZT_fixed(mz);
                    %                     end
                    %                     zt3 = zt2;
                    %                     zt3(1:mz) = mzv;
                    %                     tf=5;
                    %                     filt1=filtfilt(ones(1,tf)/tf,1,zt3(mz+1:end)); %Smooth part of the orignial artifact template
                    %                     zt3(mz+1:end) = filt1;
                    %                     dZT3 = [0;diff(zt3)];
                    %                     pl = min(find(dZT3>0));
                    %                     zt3(1:pl(1)) = -2;
                    %                     [locMaxV, locMax] = max(zt3(1:200));
                    %                     zt3(1:locMax) = 2;
                    %                     [locMin2V, locMin2] = min(zt3(1:350));
                    %                     ApTemp2 = P.AVG(fp-5:end-300);
                    %                     [mpv,mp] = max(ApTemp2);
                end
                
                
                
                if debugPlotFlag
                    keyboard
                end
                if locMaxV/mpv>500
                    tp = (1:length(ApTemp2))+toShift+20;
                    zt2(tp) = zt2(tp)-ApTemp2;
                    ZT_fixed = zt2;
                    partialFlag = 0;
                    fprintf('Edited template...\n');
                    obj.AZT_fixed = [obj.AZT_fixed 2];
                elseif locMaxV/mpv>0.4
                    tp = locMax-mp:locMax-mp+length(ApTemp2)-1;
                    zt2(tp) = zt2(tp)-ApTemp2;
                    ZT_fixed = zt2;
                    partialFlag = 0;
                    fprintf('Edited template...\n');
                    obj.AZT_fixed = [obj.AZT_fixed 2];
                elseif locMaxV/mpv>0.2
                    partialFlag = 1;
                    ZT_fixedO = ZT_fixed;
                    tp = locMax-mp:locMax-mp+length(ApTemp2)-1;
                    zt2(tp) = zt2(tp)-ApTemp2;
                    ZT_fixedSub = zt2;
                    obj.AZT_fixed = [obj.AZT_fixed 1];
                else
                    partialFlag = 0;
                    ZT_fixed = ZT_fixed;
                    obj.AZT_fixed = [obj.AZT_fixed 0];
                end
                for i=1:length(StimIdx)
                    %fprintf('%d  idx %d\n', i, StimIdx(i));
                    IDX = round( (1:length(ZT_fixed)) + StimIdx(i));
                    if IDX(end) > length(DSTIM)
                        break;
                    end
                    
                    if ~partialFlag
                        DSTIM(IDX) = DSTIM(IDX) - ZT_fixed;
                    else
                        dd = DSTIM(IDX);
                        
                        zr = mean(dd(end-300:end-100));
                        dd(1:pl) = -2;
                        [locMaxVTest, locMaxTest] = max(dd(1:450)-zr);
                        %                         f = figure;
                        %                         a = axes(f);
                        %                         plot(DSTIM(IDX))
                        %                         a.XLim = [5.7604 220.0461];
                        %                         a.YLim = [-0.5831 1.2362];
                        %                         uiwait(f)
                        if locMaxVTest/mpv>0.3
                            DSTIM(IDX) = DSTIM(IDX) - ZT_fixedSub;
                            ZT_fixed = ZT_fixedSub;
                        else
                            DSTIM(IDX) = DSTIM(IDX) - ZT_fixedO;
                            ZT_fixed = ZT_fixedO;
                        end
                    end
                end
                % Re-find the APs in the newly zapped data.
                S.DAT = DSTIM;
                %S = obj.FindAPs(DSTIM, -0.015);
                %S = obj.FindAPs(DSTIM, -0.007);
                
                XIDX = (1:length(DSTIM)) / 200000.0;
                [i1, i2, iPk, d] = obj.findAPs_Brian(XIDX, DSTIM, P.AVG(30:250), DSTIMO, debugPlotFlag);
                S.AP = sort(iPk);
                
                S.ZT_fixed = ZT_fixed;
            end
            S.StimIdx = StimIdx;  % For convenience, save the stim indices.
            S.BlockIdx = blockidx;
            
            % Save the stim index for each AP.
            S.APStimIdx = obj.MatchAPToStim(S.AP, S.StimIdx);
            
            % Save the "modulo 5mSec" AP indices. Subtract out the stim
            % index from each AP index.
            S.APModulo = S.AP(:) - S.StimIdx(S.APStimIdx);
            %             f = figure('units','normalized','outerposition',[0 0 1 1]);
            %             axt = axes(f);
            %             t = 1:length(DSTIMO);
            %             grid(axt,'on')
            %             plot(t,DSTIMO,'b')
            %             hold(axt,'on')
            %             plot(t,DSTIM,'k')
            %             plot(t(S.StimIdx),DSTIMO(S.StimIdx),'r*')
            %             plot(t(S.StimIdx(S.APStimIdx)),DSTIMO(S.StimIdx(S.APStimIdx)),'g*')
            %             plot(t(S.AP),DSTIM(S.AP),'y*')
            %             hold(axt,'off')
            %             uiwait(f)
            if debugPlotFlag
                close(ff)
            end
            
        end
        
        % For initial testing/debug, just mark pre-block APs.
        %function [DPRE,SS,CLUMP,AP, SMOO, AVG] = MarkPre(obj, blockidx)
        function varargout = MarkPre(obj, blockidx)
            DPRE = obj.GetBlockData(blockidx);
            [varargout{1:nargout}] = obj.FindAPs(DPRE);
        end
        
        function varargout = MarkPost(obj, blockidx)
            [~,DPOST] = obj.GetBlockData(blockidx);
            
            % For now, use arbitrary threshold of 0.01 on the ShiftSubtract
            % function output. In future, use knowlege of the pre-stim
            % block AP amplitude.
            [varargout{1:nargout}] = obj.FindAPs(DPOST, -0.01);
        end
        
        % For APs that are within a stim block, match each AP with the stim
        % pulse that comes before it, so we can subtract out the stim time
        % to create the "modulo 5mSec" histogram plots.
        %
        % Pass in arra of AP sample indices, and stim pulse sample indices.
        function APStimIdx = MatchAPToStim(obj, APs, Stims)
            APStimIdx = zeros(length(APs),1);
            
            % We add a "fake" stim pulse at the end that is guaranteed to
            % be larger than the last AP. Makes the loop logic simpler.
            Stims(end+1) = APs(end)+100;
            
            % We step through the Stim indices, and make sure we always use
            % the one that is just before the current AP.
            Stim_cur = 1;
            Stim_next = 2;
            for i = 1:length(APs)
                % Find stim that comes just before this AP.
                while APs(i) >= Stims(Stim_next)
                    Stim_cur = Stim_cur+1;
                    Stim_next = Stim_next+1;
                end
                APStimIdx(i) = Stim_cur;
            end
        end
        
        % Perform the actual shifted subtraction that gives us the peaks at the
        % action potentials.
        function [APPEAK,APPEAK2] = DoShiftSubtract(obj, DAT)
            
            % By eye, for a particular file (Nancy 3/4/2020 Data4) found the timing to
            % be about 0uSec, 280uSec, 680uSec, and 1500uSec. 200K samples/sec, 0.2
            % samples/usec.
            NSAMP1 =  280 * 0.2;
            NSAMP2 =  680 * 0.2;
            NSAMP3 = 1500 * 0.2;
            
            % Shift the data to line up the 4 main "event points" of the action
            % potential (start, first peak, second peak, end), but maintain the
            % overall data array length.
            SHIFT1 = [DAT(NSAMP1:end); repelem(DAT(end), NSAMP1-1)'];
            SHIFT2 = [DAT(NSAMP2:end); repelem(DAT(end), NSAMP2-1)'];
            SHIFT3 = [DAT(NSAMP3:end); repelem(DAT(end), NSAMP3-1)'];
            
            % Using 4 points is a bit better than using just the two peaks. Not by
            % much, but probably worth it to get the 10-20% increase in the peaks. Note
            % that, assuming that the unshifted DAT and final point SHIFT3 are roughly
            % the baseline, then they cancel, and the whole thing simplifies to the
            % second version using just SHIFT2-SHIFT1, with a factor of 2 since they
            % appear twice here.
            APPEAK = (DAT-SHIFT1) + (SHIFT2-SHIFT1) + (SHIFT2-SHIFT3);
            
            % Almost as good as the 3-term thing above.
            APPEAK2 = 2*(SHIFT2-SHIFT1);
            
        end
        
    end
    
end