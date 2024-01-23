classdef Template_class <handle
    properties
		Artifact_Template
        AP_Template
        DAT
    end
    
    methods
		% Constructor takes a SpikeDATA object.
		function obj = Template_class(SpikeDATA_object)
			obj.DAT = SpikeDATA_object;
        end
        
        function [AllTEMPLATE_unzeroed, AllTEMPLATE] = GetAllArtifactTemplate(obj, chan, handles)
            temp = obj.GetArtifactTemplate(1, chan, []);
            NSAMPS = length(temp);
            NBLOCKS = length(obj.DAT.StimBlocks);
            
            AllTEMPLATE_unzeroed = zeros(NSAMPS, NBLOCKS);
            AllTEMPLATE = zeros(NSAMPS, NBLOCKS);
            if ~isempty(handles)
                handles.ProgressBar.Title.String = ['Getting Initial Artifact Templates'];
                drawnow
            end
            for i=1:NBLOCKS
                [AllTEMPLATE_unzeroed(:,i), AllTEMPLATE(:,i)] = obj.GetArtifactTemplate(i, chan, []);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_unzeroed_original']) = AllTEMPLATE_unzeroed(:,i)';
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_original']) = AllTEMPLATE(:,i)';
                sz = size(AllTEMPLATE_unzeroed(:,i)');
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_unzeroed_used']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_used']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_unzeroed_All']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_All']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_unzeroed_t']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_t']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_ArtifactTemplate_confirmed']) = 0;
                obj.DAT.StimBlocks(i).([chan,'_APLoc']) = zeros(sz);
                obj.DAT.StimBlocks(i).([chan,'_APCount']) = 0;
                if ~isempty(handles)
                    handles.PBarObj.Position(3) = (i)/NBLOCKS*1000;
                    handles.PBarTxt.String = [num2str(round((i)/NBLOCKS*100)),'%'];
                    drawnow
                end

            end
            if ~isempty(handles)
                handles.ProgressBar.Title.String = ['Finished'];
                handles.PBarObj.Position(3) = 0;
                handles.PBarTxt.String = 'Waiting';
                drawnow
            end
            id = obj.DAT.GetChan(chan);
            obj.DAT.CHANS{id.Num}.AllTEMPLATE_unzeroed = AllTEMPLATE_unzeroed;
            obj.DAT.CHANS{id.Num}.AllTEMPLATE = AllTEMPLATE;
        end
        % Pass in the index of the block for which we want to find the Zap
        % template.
        function [TEMPLATE_unzeroed,TEMPLATE] = GetArtifactTemplate(obj, Block_idx, chan, ALLSTIM)
            if ~isempty(chan)
                DAT = obj.DAT; % For convenience, copy over.
                ch = DAT.GetChan(chan);
                BLOCK = DAT.StimBlocks(Block_idx);
                STIMS = find((DAT.Stim_ticks >= [BLOCK.Start_tick]) ...
                    & (DAT.Stim_ticks <= [BLOCK.End_tick]));
                
                % ALWAYS zap across THE WHOLE STIM PULSE INTERVAL!!!
                % NEVER SHORTER!!
                % We rely on the full "wrap around" of the template, when
                % taking the mean/zero below as the last few samples of the
                % template.
                ZAP_DURATION_sec = 0.005;
                NSAMPS = round(ZAP_DURATION_sec * ch.SampleRate);
                
                % ALLSTIM is a matrix, one row per stim, with each row containing
                % the analog samples of the stim artifact.
                ALLSTIM = zeros(length(STIMS), NSAMPS);
                
                % Loop over the stim pulses, and gather into ALLSTIM matrix.
                for i=1:length(STIMS)
                    tick = DAT.Stim_ticks(STIMS(i));
                    Vals = obj.DAT.ValsAtTick(tick, NSAMPS, chan);
                    if isempty(Vals); continue; end
                    % Now plunk in the array of data for each stim.
                    ALLSTIM(i,:) = Vals;
                end
            end
            % Average across ALL STIM PULSE data to create the template.
            TEMPLATE_unzeroed = trimmean(ALLSTIM, 20);
            
            % NOW we have a single averaged "trace", and we need to "zero"
            % it. Take the mean of the last samples towards the end of the
            % template as the "zero" for the template.
            N_AVG_ZERO = 50;
            TEMPLATE_ZERO = mean(TEMPLATE_unzeroed(end-N_AVG_ZERO:end));
            TEMPLATE = TEMPLATE_unzeroed - TEMPLATE_ZERO;
        end
        
        function [tss, tss2] = InternalCompare(obj, idxs, debugflag, chan)
            D = obj.DAT;
            tss = struct();
            tss2 = struct();
            origStim = D.StimBlocks(idxs).([chan,'_ALLSTIM']);
            APMax = max(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG);
            currVal = D.StimBlocks(idxs).Current_uA;
            tss.t1 = origStim(1,:);
            [~,Loc] = max(origStim(1,90:end));
            tss2.t1Loc = Loc+89;
            passed100 = 0;
            for i = 2:198
                fs = fields(tss);
                ds = [];
                stdVals = [];
                for j = 1:length(fs)
                    %         findsignal(origStim(i,65:565),mean(tss.(fs{j})(:,65:565),1))
                    [~,~,dist] = findsignal(origStim(i,65:565),mean(tss.(fs{j})(:,65:565),1));
                    ds = [ds dist];
                    stdVals = [stdVals std(tss2.([fs{j},'Loc']))];
                end
                [mindV, mindL] = min(ds);
                if mindV < 1
                    name = ['t',num2str(mindL)];
                    tss.(name) = [tss.(name);origStim(i,:)];
                    [~,Loc] = max(origStim(i,90:end));
                    tss2.([name,'Loc']) = [tss2.([name,'Loc']);Loc+89];
                elseif (~passed100) && (max(origStim(i,100:end))/APMax>0.6)
                    name = ['t',num2str(length(fs)+1)];
                    tss.(name) = origStim(i,:);
                    [~,Loc] = max(origStim(i,90:end));
                    tss2.([name,'Loc']) = [];
                    tss2.([name,'Loc']) = [tss2.([name,'Loc']);Loc+89];
                    passed100 = 1;
                elseif (passed100) && (max(origStim(i,100:end))/APMax>0.6)
                    name = ['t',num2str(mindL)];
                    tss.(name) = [tss.(name);origStim(i,:)];
                    [~,Loc] = max(origStim(i,90:end));
                    tss2.([name,'Loc']) = [tss2.([name,'Loc']);Loc+89];
                else
                    name = ['t',num2str(length(fs)+1)];
                    tss.(name) = origStim(i,:);
                    [~,Loc] = max(origStim(i,90:end));
                    tss2.([name,'Loc']) = [];
                    tss2.([name,'Loc']) = [tss2.([name,'Loc']);Loc+89];
                end
            end
        end
        
        function [TemplateStyle, Positivepks, cornerpt] = GetArtTempDets(obj, idxs,  chan, artTemp, apLength)
            D = obj.DAT;
            origStim = D.StimBlocks(idxs).([chan,'_ALLSTIM']);
            APMax = max(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG);
            currVal = D.StimBlocks(idxs).Current_uA;
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle1.mat');
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle2.mat');
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle3.mat');
            [~,~,dt1]=findsignal(ArtifactStyle1,artTemp,'Normalization','zscore');
            [~,~,dt2]=findsignal(ArtifactStyle2,artTemp,'Normalization','zscore');
            [~,~,dt3]=findsignal(ArtifactStyle3,artTemp,'Normalization','zscore');
            [~,~,dt1n]=findsignal(-ArtifactStyle1,artTemp,'Normalization','zscore');
            [~,~,dt2n]=findsignal(-ArtifactStyle2,artTemp,'Normalization','zscore');
            [~,~,dt3n]=findsignal(-ArtifactStyle3,artTemp,'Normalization','zscore');
            [~, closestN] = min([dt1 dt2 dt3 dt1n dt2n dt3n]);
            switch closestN
                case 1
                    [mv, ml] = max(artTemp);
                    p = find(artTemp==mv,1,'last')+5;
                    TemplateStyle = 1;
                    Positivepks = p-5;
                case 2
                    smTemp = smooth(artTemp,7,'sgolay');
                    [TF, P] = islocalmax(smTemp); %artTemp
                    %[pkV,pkLoc]=maxk(P,2);
                    %p2u = max(pkLoc);
                    pnz = (P~=0);
                    ptemp = find(P>mean(P(pnz)));
                    if length(ptemp) > 2
                        mv = max(artTemp);
                        compare2max = artTemp(ptemp)/mv;
                        ptemp(compare2max<0.1) = [];
                        [~,minL]=min(smTemp);
                        if length(ptemp) > 2
                            p2v = ptemp(find(ptemp>minL,1,'first'));
                            [~,ll]=min(abs(ptemp-find(artTemp==mv)));
                            if ptemp(ll) == p2v
                                p1v = ptemp(find(ptemp<minL,1,'last'));
                            else
                                p1v = ptemp(ll);
                            end

                            ptemp(~ismember(ptemp,[p1v p2v])) = [];
                            %                             if length(ptemp) == 3
                            %                                 ptemp(2) = [];
                            %                             else
                            %                                 ptemp(P(ptemp)<.1) = [];
                            %                             end

                        end
                        flatM = find(artTemp==artTemp(ptemp(1)));
                        pmax1 = flatM(end);
                        if ptemp == pmax1
                            p = ptemp+5;
                            Positivepks = ptemp;
                        else
                            sp = find(ptemp>pmax1,1);
                            flatM2 = find(artTemp==artTemp(ptemp(sp)));
                            pmax2 = flatM2(end);
                            p = pmax2+5;
                            Positivepks = [pmax1 pmax2];
                        end
                    else
                        pint = find(ptemp<200,1,'last');
                        p2u = ptemp(pint);
                        p = find(artTemp==artTemp(p2u),1,'last')+5;
                        Positivepks = ptemp;
                    end
                    TemplateStyle = 2;
                case 3
                    smTemp = smooth(artTemp,7,'sgolay');
                    [TF, P] = islocalmax(smTemp);
                    pnz = (P~=0);
                    ptemp = find(P>mean(P(pnz)));
                    ptempt = ptemp;
                    if length(ptemp) > 2
                        mv = max(artTemp);
                        compare2max = artTemp(ptemp)/mv;
                        ptemp(compare2max<0.1) = [];
                        if length(ptemp) > 2
                            if length(ptemp) == 3
                                ptemp(2) = [];
                            else
                                ptemp(P(ptemp)<.1) = [];
                            end

                        end
                        if isempty(ptemp)
                            [a,d]=maxk(P(ptempt),2);
                            ptemp = ptempt(d);
                        end
                        flatM = find(artTemp==artTemp(ptemp(1)));
                        pmax1 = flatM(end);
                        if length(ptemp)>1
                            sp = find(ptemp>pmax1,1);
                        else
                            sp=1;
                        end
                        flatM2 = find(artTemp==artTemp(ptemp(sp)));
                        pmax2 = flatM2(end);
                        p = pmax2+5;
                        Positivepks = [pmax1 pmax2];
                    else
                        [~, mms] = max(artTemp(ptemp));
                        p2u = ptemp(mms);
                        p = find(artTemp==artTemp(p2u),1,'last')+5;
                        Positivepks = ptemp;
                    end
                    TemplateStyle = 3;
                case 4
                    [mv, ml] = max(artTemp);
                    p = find(artTemp==mv,1,'last')+5;
                    TemplateStyle = 4;
                    Positivepks = p-5;
                case 5
                    smTemp = smooth(artTemp,7,'sgolay');
                    [TF, P] = islocalmax(smTemp); %artTemp
                    %[pkV,pkLoc]=maxk(P,2);
                    %p2u = max(pkLoc);
                    pnz = (P~=0);
                    ptemp = find(P>mean(P(pnz)));
                    if length(ptemp) > 2
                        mv = max(artTemp);
                        compare2max = artTemp(ptemp)/mv;
                        ptemp(compare2max<0.1) = [];
                        [~,minL]=min(smTemp);
                        if length(ptemp) > 2
                            p2v = ptemp(find(ptemp>minL,1,'first'));
                            [~,ll]=min(abs(ptemp-find(artTemp==mv)));
                            if ptemp(ll) == p2v
                                p1v = ptemp(find(ptemp<minL,1,'last'));
                            else
                                p1v = ptemp(ll);
                            end

                            ptemp(~ismember(ptemp,[p1v p2v])) = [];
                            %                             if length(ptemp) == 3
                            %                                 ptemp(2) = [];
                            %                             else
                            %                                 ptemp(P(ptemp)<.1) = [];
                            %                             end

                        end
                        flatM = find(artTemp==artTemp(ptemp(1)));
                        pmax1 = flatM(end);
                        if ptemp == pmax1
                            p = ptemp+5;
                            Positivepks = ptemp;
                        else
                            sp = find(ptemp>pmax1,1);
                            flatM2 = find(artTemp==artTemp(ptemp(sp)));
                            pmax2 = flatM2(end);
                            p = pmax2+5;
                            Positivepks = [pmax1 pmax2];
                        end
                    else
                        pint = find(ptemp<200,1,'last');
                        p2u = ptemp(pint);
                        p = find(artTemp==artTemp(p2u),1,'last')+5;
                        Positivepks = ptemp;
                    end
                    TemplateStyle = 5;
                case 6
                    smTemp = smooth(artTemp,7,'sgolay');
                    [TF, P] = islocalmax(smTemp);
                    pnz = (P~=0);
                    ptemp = find(P>mean(P(pnz)));
                    ptempt = ptemp;
                    if length(ptemp) > 2
                        mv = max(artTemp);
                        compare2max = artTemp(ptemp)/mv;
                        ptemp(compare2max<0.1) = [];
                        if length(ptemp) > 2
                            if length(ptemp) == 3
                                ptemp(2) = [];
                            else
                                ptemp(P(ptemp)<.1) = [];
                            end

                        end
                        if isempty(ptemp)
                            [a,d]=maxk(P(ptempt),2);
                            ptemp = ptempt(d);
                        end
                        flatM = find(artTemp==artTemp(ptemp(1)));
                        pmax1 = flatM(end);
                        if length(ptemp)>1
                            sp = find(ptemp>pmax1,1);
                            if isempty(sp)
                                sp = 1;
                            end
                        else
                            sp=1;
                        end
                        flatM2 = find(artTemp==artTemp(ptemp(sp)));
                        pmax2 = flatM2(end);
                        p = pmax2+5;
                        Positivepks = [pmax1 pmax2];
                    else
                        [~, mms] = max(artTemp(ptemp));
                        p2u = ptemp(mms);
                        p = find(artTemp==artTemp(p2u),1,'last')+5;
                        Positivepks = ptemp;
                    end
                    TemplateStyle = 6;
            end
            
            artTemp2 = artTemp;
            artTemp2(1:p) = artTemp(p);
            
            artTemp2D = [0 diff(artTemp2)];

            smTemp = smooth(artTemp,7,'sgolay');
            [we,we2] = islocalmin(smTemp(p+1:300));
            locMinLT = find(we2>0,1,'first');
            locMinL = find(artTemp2D(p+1:300)>-0.005,1,'first');
            lmlt = locMinLT+p;
            lml = locMinL+p;
            
            if lmlt <135
                locMinL = locMinLT;
            end
            
            if isempty(locMinL)
                [locMinV,locMinL] = max(artTemp2D(70:115));
                locMinL = locMinL+70;
                pt2Start = locMinL-1;
            else
                if currVal < 200
                    pt2Start = locMinL+p-1;
                else
                    pt2Start = locMinL+p;
                end
            end
            
            if ismember({chan},'NerveCon')
                davg = mean(diff(artTemp(pt2Start - 17:pt2Start - 13)));
                subtpt = 15;
                if obj.DAT.StimBlocks(idxs).Current_uA ~= 200
                    while davg < -0.01
                        subtpt = subtpt -1;
                        davg = mean(diff(artTemp(pt2Start - (subtpt+2):pt2Start - (subtpt-2))));
                    end
                end
                pt2Start = pt2Start - subtpt;
               
            end
            pt2End = pt2Start+apLength;%110;
            cornerpt = pt2Start;
        end
        
        function [timeTrim, sigTrim, Template, TemplateStyle, Positivepks, cornerpt] = SplineArtifactTemplate(obj, artTemp, apLength, currVal, idxs, debugflag, chan, debugPlots)
            % 1. find the last max point
            % 2. spline from 10 data points after max on
            % 3. find local min after max
            % 4. find point between consistent steep slope and local min
            % 5. blank from that mid point to 100th data point
            % 6.
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle1.mat');
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle2.mat');
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle3.mat');
            load('R:\Morris, Brian\MATLAB\SUR Data Analysis\QuickZap\ArtifactStyle4.mat');
            [~,~,dt1]=findsignal(ArtifactStyle1,artTemp,'Normalization','zscore');
            [~,~,dt2]=findsignal(ArtifactStyle2,artTemp,'Normalization','zscore');
            [~,~,dt3]=findsignal(ArtifactStyle3,artTemp,'Normalization','zscore');
            [~,~,dt4]=findsignal(ArtifactStyle4,artTemp,'Normalization','zscore');
            ato = artTemp;
            [~, closestN] = min([dt1 dt2 dt3 dt4]);
            if (obj.DAT.StimBlocks(idxs).Ref_E == 10) && (ismember({chan},'ZapOut'))
                if closestN < 4
                    closestN = 1;
                end
            end
            chanInfo = obj.DAT.GetChan(chan);
            switch closestN
                case 1
                    %% If it matches ArtifactStyle1
                    sign2U = 1;
                    [mv, ml] = max(artTemp);
                    p = find(artTemp==mv,1,'last')+5;
                    TemplateStyle = 1;
                    Positivepks = p-5;
                case 2
                    %% If it matches ArtifactStyle2
                    sign2U = 1;
                    smTemp = smooth(artTemp,7,'sgolay');
                    [TF, P] = islocalmax(smTemp); %artTemp
                    pnz = (P~=0);
                    ptemp = find(P>mean(P(pnz)));
                    if length(ptemp) > 2
                        mv = max(artTemp);
                        compare2max = artTemp(ptemp)/mv;
                        ptemp(compare2max<0.1) = [];
                        [~,minL]=min(smTemp);
                        if length(ptemp) > 2
                            p2v = ptemp(find(ptemp>minL,1,'first'));
                            [~,ll]=min(abs(ptemp-find(artTemp==mv)));
                            if ptemp(ll) == p2v
                                p1v = ptemp(find(ptemp<minL,1,'last'));
                            else
                                p1v = ptemp(ll);
                            end

                            ptemp(~ismember(ptemp,[p1v' p2v])) = [];
                            %                             if length(ptemp) == 3
                            %                                 ptemp(2) = [];
                            %                             else
                            %                                 ptemp(P(ptemp)<.1) = [];
                            %                             end

                        end
                        flatM = find(artTemp==artTemp(ptemp(1)));
                        pmax1 = flatM(end);

                        if ptemp == pmax1
                            p = ptemp+5;
                            Positivepks = ptemp;
                        else
                            if all(ptemp<pmax1)
                                pmax2 = pmax1;
                                pmax1 = max(ptemp);
                                p = pmax2+5;
                                Positivepks = [pmax1 pmax2];
                            else
                                sp = find(ptemp>pmax1,1);
                                if isempty(sp)
                                    sp = length(ptemp);
                                end
                                flatM2 = find(artTemp==artTemp(ptemp(sp)));
                                pmax2 = flatM2(end);
                                p = pmax2+5;
                                Positivepks = [pmax1 pmax2];
                            end
                            
                        end
                    else
                        pint = find(ptemp<200,1,'last');
                        p2u = ptemp(pint);
                        p = find(artTemp==artTemp(p2u),1,'last')+5;
                        Positivepks = ptemp;
                    end
                    TemplateStyle = 2;
                case 3
                    %% If it matches ArtifactStyle3
                    sign2U = 1;
                    smTemp = smooth(artTemp,7,'sgolay');
                    [TF, P] = islocalmax(smTemp);
                    pnz = (P~=0);
                    ptemp = find(P>mean(P(pnz)));
                    if length(ptemp) > 2
                        mv = max(artTemp);
                        compare2max = artTemp(ptemp)/mv;
                        ptemp(compare2max<0.1) = [];
                        ptempt = ptemp;
                        if length(ptemp) > 2
                            if length(ptemp) == 3
                                ptemp(2) = [];
                            else
                                ptemp(P(ptemp)<.1) = [];
                            end

                        end
                        if isempty(ptemp)
                            [a,d]=maxk(P(ptempt),2);
                            ptemp = ptempt(d);
                        end
                        flatM = find(artTemp==artTemp(ptemp(1)));
                        pmax1 = flatM(end);
                        if length(ptemp)>1
                            sp = find(ptemp>pmax1,1);
                            if isempty(sp)
                                sp = length(ptemp);
                            end
                        else
                            sp=1;
                        end
                        flatM2 = find(artTemp==artTemp(ptemp(sp)));
                        pmax2 = flatM2(end);
                        p = pmax2+5;
                        Positivepks = [pmax1 pmax2];
                    else
                        [~, mms] = max(artTemp(ptemp));
                        p2u = ptemp(mms);
                        p = find(artTemp==artTemp(p2u),1,'last')+5;
                        Positivepks = ptemp;
                    end
                    TemplateStyle = 3;
                case 4
                    %% If it matches the negative of ArtifactStyle4
                    sign2U = -1;
                    artTemp = artTemp*-1;
                    [mv, ml] = max(artTemp);
                    p = find(artTemp==mv,1,'last')+5;
                    TemplateStyle = 4;
                    Positivepks = p-5;
            end

            if debugflag
                if ~isempty(debugPlots)
                    plot(debugPlots.artPoints,artTemp*sign2U)
                    hold(debugPlots.artPoints,'on')
                    plot(debugPlots.artPoints,p,artTemp(p)*sign2U,'r*')
                else
                    figure
                    plot(artTemp*sign2U)
                    hold on
                    plot(p,artTemp(p)*sign2U,'r*')
                end
            end
            %tf = 5;
            %splTemp = filtfilt(ones(1,tf)/tf,1,artTemp(p:end));
            artTemp2 = artTemp;
            artTemp2(1:p) = artTemp(p);
            
            if debugflag
                if ~isempty(debugPlots)
                    plot(debugPlots.artPoints,artTemp2*sign2U)
                else
                    plot(artTemp2)
                end
            end
            
            artTemp2D = [0 diff(artTemp2)];
            if debugflag
                if ~isempty(debugPlots)
                    plot(debugPlots.artPoints,artTemp2D*sign2U)
                else
                    plot(artTemp2D)
                end
            end
            
            %locMinL = find(artTemp2D(p+1:100)>0,1);
            smTemp = smooth(artTemp,7,'sgolay');
            [we,we2] = islocalmin(smTemp(p+1:300));
            locMinLT = find(we2>0,1,'first');
            locMinL = find(artTemp2D(p+1:300)>-0.005,1,'first');
            lmlt = locMinLT+p;
            lml = locMinL+p;
            
            if (lmlt <135) && (TemplateStyle ~= 4)
                locMinL = locMinLT;
            end
            
            if isempty(locMinL)
                [locMinV,locMinL] = max(artTemp2D(70:115));
                locMinL = locMinL+70;
                pt2Start = locMinL-1;
            else
                if debugflag
                    if ~isempty(debugPlots)
                        plot(debugPlots.artPoints,locMinL+p, artTemp(locMinL+p)*sign2U, 'g*')
                    else
                        plot(locMinL+p, artTemp(locMinL+p), 'g*')
                    end
                end
                %[locMinV,locMinL]=min(abs(artTemp2D(p+1:70))); may need to add
                %an if statement
                %locMinL = 91-p-1;
                if currVal < 200
                    pt2Start = locMinL+p-1;
                else
                    pt2Start = locMinL+p;
                end
            end
            
%             if TemplateStyle == 3
%                 pt2Start=pt2Start+7;
%             end
            
            if ismember({chan},'NerveCon')
                if (pt2Start-17) < 0
                    davg = mean(diff(artTemp(1:pt2Start)));
                else
                    davg = mean(diff(artTemp(pt2Start - 17:pt2Start - 13)));
                end
                subtpt = 15;
                if obj.DAT.StimBlocks(idxs).Current_uA ~= 200
                    if obj.DAT.StimBlocks(idxs).Current_uA == 150
                        mAP = max(obj.DAT.StimBlocks(idxs).([chan,'_APTemplate']).AVG);
                        time = round(0:chanInfo.SampleInterval_sec:(length(artTemp)-1)*chanInfo.SampleInterval_sec,7); %normalized time of the artifact template
                        artTempT = smooth(artTemp,7,'sgolay');
                        bbd = [pt2Start-subtpt:pt2Start+apLength-subtpt];
                        [curve, goodness, output] = fit(time',artTempT,'smoothingspline','Exclude',bbd); %Fitting a spline to the smoothed and blanked artifact template
                        yTemplate = curve(time(bbd))';
                        mvmm = max(artTempT(bbd)'-yTemplate);
                        limm = -0.015;
                        if davg < limm
                            while (davg < limm) && (mvmm < mAP*.9)
                                subtpt = subtpt -1;
                                davg = mean(diff(artTemp(pt2Start - (subtpt+2):pt2Start - (subtpt-2))));
                                bbd = [pt2Start-subtpt:pt2Start+apLength-subtpt];
                                [curve, goodness, output] = fit(time',artTempT,'smoothingspline','Exclude',bbd); %Fitting a spline to the smoothed and blanked artifact template
                                yTemplate = curve(time(bbd))';
                                mvmm = max(artTempT(bbd)'-yTemplate);
                            end
                        else
                            while (davg > limm) || (mvmm <= mAP)
                                subtpt = subtpt +1;
                                if (pt2Start - (subtpt+2))== 0
                                    break
                                else
                                    davg = mean(diff(artTemp(pt2Start - (subtpt+2):pt2Start - (subtpt-2))));
                                    bbd = [pt2Start-subtpt:pt2Start+apLength-subtpt];
                                    [curve, goodness, output] = fit(time',artTempT,'smoothingspline','Exclude',bbd); %Fitting a spline to the smoothed and blanked artifact template
                                    yTemplate = curve(time(bbd))';
                                    mvmm = max(artTempT(bbd)'-yTemplate);
                                end
                                
                            end
                        end
                        
                    else
                        limm = -0.01;
                        while davg < limm
                            subtpt = subtpt -1;
                            davg = mean(diff(artTemp(pt2Start - (subtpt+2):pt2Start - (subtpt-2))));
                        end
                    end
                else
                    mAP = max(obj.DAT.StimBlocks(idxs).([chan,'_APTemplate']).AVG);
                    time = round(0:chanInfo.SampleInterval_sec:(length(artTemp)-1)*chanInfo.SampleInterval_sec,7); %normalized time of the artifact template
                    artTempT = smooth(artTemp,7,'sgolay');
                    bbd = [pt2Start-subtpt:pt2Start+apLength-subtpt];
                    [curve, goodness, output] = fit(time',artTempT,'smoothingspline','Exclude',bbd); %Fitting a spline to the smoothed and blanked artifact template
                    yTemplate = curve(time(bbd))';
                    mvmm = max(artTempT(bbd)'-yTemplate);
                    if mvmm > mAP*2.2
                        while mvmm > mAP*2.2
                            subtpt = subtpt -1;
                            davg = mean(diff(artTemp(pt2Start - (subtpt+2):pt2Start - (subtpt-2))));
                            bbd = [pt2Start-subtpt:pt2Start+apLength-subtpt];
                            [curve, goodness, output] = fit(time',artTempT,'smoothingspline','Exclude',bbd); %Fitting a spline to the smoothed and blanked artifact template
                            yTemplate = curve(time(bbd))';
                            mvmm = max(artTempT(bbd)'-yTemplate);
                        end
                    elseif mvmm < mAP
                        while mvmm < mAP
                            subtpt = subtpt +1;
                            davg = mean(diff(artTemp(pt2Start - (subtpt+2):pt2Start - (subtpt-2))));
                            bbd = [pt2Start-subtpt:pt2Start+apLength-subtpt];
                            [curve, goodness, output] = fit(time',artTempT,'smoothingspline','Exclude',bbd); %Fitting a spline to the smoothed and blanked artifact template
                            yTemplate = curve(time(bbd))';
                            mvmm = max(artTempT(bbd)'-yTemplate);
                        end
                    end
                end
                
                if pt2Start - subtpt < p
                    subtpt = 0;
                end
                
                if obj.DAT.StimBlocks(idxs).Current_uA == 150
                    if length(find(max(artTemp)==artTemp)) > 5
                        subtpt = subtpt-3;
                    end
                end
                pt2Start = pt2Start - subtpt;
                
               
            end
            pt2End = pt2Start+apLength;%110;
            cornerpt = pt2Start;
            
            if debugflag
                if ~isempty(debugPlots)
                        plot(debugPlots.artPoints,pt2Start, artTemp(pt2Start)*sign2U,'b*')
                else
                    plot(pt2Start, artTemp(pt2Start)*sign2U,'b*')
                end
            end
            
            
            
            if debugflag
                if ~isempty(debugPlots)
                        plot(debugPlots.artPoints, pt2End, artTemp(pt2End)*sign2U,'k*')
                else
                    plot(pt2End, artTemp(pt2End),'k*')
                end
            end
            %             [TF,P] = islocalmax(artTemp(pt2Start:pt2End));
            %             if any(P>0.03)
            %                 lml = find(P>0.03,1);
            %                 localmax = lml + pt2Start;
            %                 [locMinV,locMinL]=min(abs(artTemp2D(localmax+1:100)));
            %                 pt2Start = locMinL+localmax-3;
            %             end
            
            time = round(0:chanInfo.SampleInterval_sec:(length(artTemp)-1)*chanInfo.SampleInterval_sec,7); %normalized time of the artifact template
            % 			[~,preMid]=min(abs(time-(0.0007-(b4mid)))); % Start index of blank section where a confounded AP would be
            % 			[~,postMid]=min(abs(time-(0.0007+(aftermid)))); % End index of blank section where a confounded AP would be
            tt = 1:length(time);
            tt = [tt(1:pt2Start) tt(pt2End:end)];
            artTemp(p:end) = smooth(artTemp2(p:end),7,'sgolay');
            timeTrim = [time(1:pt2Start) time(pt2End:end)]; %Blanking time data
            sigTrim = [artTemp(1:pt2Start) artTemp(pt2End:end)]; %Blanking artifact template
            
            if debugflag
                if ~isempty(debugPlots)
                    plot(debugPlots.artPoints, tt, sigTrim*sign2U)
                else
                    plot(tt, sigTrim*sign2U)
                end
            end
            
            [curve, goodness, output] = fit(timeTrim',sigTrim','smoothingspline'); %Fitting a spline to the smoothed and blanked artifact template
            xTemplate = time;
            yTemplate = curve(time)';
            
            if debugflag
                if ~isempty(debugPlots)
                    plot(debugPlots.artPoints, yTemplate*sign2U)
                else
                    plot(yTemplate*sign2U)
                end
            end
            yTemplate= yTemplate';
            finalT = artTemp;
            finalT(pt2Start:pt2End) = yTemplate(pt2Start:pt2End);
            Template = finalT*sign2U;
            if debugflag
                if ~isempty(debugPlots)
                    plot(debugPlots.artPoints, Template)
                else
                    plot(Template)
                end
            end
            
            if TemplateStyle == 3
                if closestN == 4
                    artTemp = -ato;
                else
                    artTemp = ato;
                end
                if apLength > 90
                    apLength = 90;
                end
                pt2End = pt2Start+apLength;%110;
                while mean(artTemp(pt2Start:pt2End)-yTemplate(pt2Start:pt2End)') > 0.025
                    if closestN == 4
                        artTemp = -ato;
                    else
                        artTemp = ato;
                    end
                    artTemp(p:end) = smooth(artTemp2(p:end),7,'sgolay');
                    pt2Start = pt2Start+1;
                    time = round(0:chanInfo.SampleInterval_sec:(length(artTemp)-1)*chanInfo.SampleInterval_sec,7); %normalized time of the artifact template
                    % 			[~,preMid]=min(abs(time-(0.0007-(b4mid)))); % Start index of blank section where a confounded AP would be
                    % 			[~,postMid]=min(abs(time-(0.0007+(aftermid)))); % End index of blank section where a confounded AP would be
                    tt = 1:length(time);
                    tt = [tt(1:pt2Start) tt(pt2End:end)];
                    timeTrim = [time(1:pt2Start) time(pt2End:end)]; %Blanking time data
                    sigTrim = [artTemp(1:pt2Start) artTemp(pt2End:end)]; %Blanking artifact template
                    [curve, goodness, output] = fit(timeTrim',sigTrim','smoothingspline'); %Fitting a spline to the smoothed and blanked artifact template
                    xTemplate = time;
                    yTemplate = curve(time);
                end
                redoflg = 1;
            else
                [mm3,mm4] = min(yTemplate(pt2Start:pt2End));
                pt = mm4+pt2Start-1;
                artTemp(pt)-mm3;
                redoflg = 0;
                ff = find(we2>0);
                if (artTemp(pt)-mm3 > 0.25) && (~ismember({chan},'NerveCon'))%|| (mean(yTemplate(pt2Start:pt2End)) > mean(artTemp(pt2Start:pt2End)))
                    locMinLT = find(we2>0);
                    if we2(locMinLT(1)) <1
                        locMinL = locMinLT(2)+p;
                        pt2Start = locMinL;
                        cornerpt = pt2Start;
                        redoflg = 1;
                    end
                elseif (mean(yTemplate(pt2Start:pt2End)) > mean(artTemp(pt2Start:pt2End)))
                    pt2Start = pt2Start-1;
                    cornerpt = pt2Start;
                    redoflg = 1;
                elseif length(ff)>1
                    if (we2(ff(2))> we2(ff(1))) && (~ismember({chan},'NerveCon'))
                        locMinLT = find(we2>0);
                        if locMinLT(2)-locMinLT(1) > 22
                            locMinL = locMinLT(2)+p;
                            pt2Start = locMinL;
                            cornerpt = pt2Start;
                            redoflg = 1;
                        end
                    end
                end
                pt2End = pt2Start+apLength;%110;
            end
            
            
            
            
            
            if redoflg
                
                time = round(0:chanInfo.SampleInterval_sec:(length(artTemp)-1)*chanInfo.SampleInterval_sec,7); %normalized time of the artifact template
                % 			[~,preMid]=min(abs(time-(0.0007-(b4mid)))); % Start index of blank section where a confounded AP would be
                % 			[~,postMid]=min(abs(time-(0.0007+(aftermid)))); % End index of blank section where a confounded AP would be
                tt = 1:length(time);
                tt = [tt(1:pt2Start) tt(pt2End:end)];
                timeTrim = [time(1:pt2Start) time(pt2End:end)]; %Blanking time data
                sigTrim = [artTemp(1:pt2Start) artTemp(pt2End:end)]; %Blanking artifact template
                
                if debugflag
                    if ~isempty(debugPlots)
                        plot(debugPlots.artPoints, tt, sigTrim*sign2U)
                    else
                        plot(tt, sigTrim)
                    end
                end
                
                [curve, goodness, output] = fit(timeTrim',sigTrim','smoothingspline'); %Fitting a spline to the smoothed and blanked artifact template
                xTemplate = time;
                yTemplate = curve(time)';
                
                if debugflag
                    if ~isempty(debugPlots)
                        plot(debugPlots.artPoints, yTemplate*sign2U)
                    else
                        plot(yTemplate)
                    end
                end
                yTemplate= yTemplate';
                finalT = artTemp;
                finalT(pt2Start:pt2End) = yTemplate(pt2Start:pt2End);
                Template = finalT*sign2U;
                if debugflag
                    if ~isempty(debugPlots)
                        plot(debugPlots.artPoints, Template)
                        hold(debugPlots.artPoints,'off')
                    else
                        plot(Template)
                    end
                end
            end
            obj.DAT.StimBlocks(idxs).ArtifactTempPKS = Positivepks;
            obj.DAT.StimBlocks(idxs).ArtifactTempCorner = cornerpt;
            %             [locminval, locminloc] = min(Template(pt2Start:pt2End));
            %             locminloc = locminloc+pt2Start-1;
            %
            %             while (locminval-artTemp2D(locminloc))>-0.011
            %                 pt2End = pt2End+1;
            %
            %                 time = round(0:obj.DAT.SampleInterval_sec:(length(artTemp)-1)*obj.DAT.SampleInterval_sec,7); %normalized time of the artifact template
            %                 % 			[~,preMid]=min(abs(time-(0.0007-(b4mid)))); % Start index of blank section where a confounded AP would be
            %                 % 			[~,postMid]=min(abs(time-(0.0007+(aftermid)))); % End index of blank section where a confounded AP would be
            %                 tt = 1:length(time);
            %                 tt = [tt(1:pt2Start) tt(pt2End:end)];
            %                 timeTrim = [time(1:pt2Start) time(pt2End:end)]; %Blanking time data
            %                 sigTrim = [artTemp(1:pt2Start)' artTemp(pt2End:end)']; %Blanking artifact template
            %                 [curve, goodness, output] = fit(timeTrim',sigTrim','smoothingspline'); %Fitting a spline to the smoothed and blanked artifact template
            %                 xTemplate = time;
            %                 yTemplate = curve(time)';
            %
            %                 Template = yTemplate';
            %
            %                 [locminval, locminloc] = min(Template(pt2Start:pt2End));
            %                 locminloc = locminloc+pt2Start-1;
            %             end
            
        end
        
        function [ALLAPTEMPLATE_unzeroed,ALLAPTEMPLATE] = GetAllAPTemplate(obj, chan, handles)
            APC = ActionPotential_class(obj.DAT);
            NBLOCKS = length(obj.DAT.StimBlocks);
            
            ALLAPTEMPLATE_unzeroed = struct();
            ALLAPTEMPLATE = struct();
            if ~isempty(handles)
                handles.ProgressBar.Title.String = ['Getting Action Potential Templates'];
                drawnow
            end
            for i=1:NBLOCKS
                [ALLAPTEMPLATE_unzeroedt, ALLAPTEMPLATEt] = obj.GetAPTemplate(APC, i, chan);
                if i == 1
                    fsz = fields(ALLAPTEMPLATE_unzeroedt);
                    fs = fields(ALLAPTEMPLATEt);
                    for jj = 1:length(fs)
                        ALLAPTEMPLATE_unzeroed.(fsz{jj}) = [];
                        ALLAPTEMPLATE.(fs{jj}) = [];
                    end
                end
                ALLAPTEMPLATE_unzeroed(i) = ALLAPTEMPLATE_unzeroedt;
                ALLAPTEMPLATE(i) = ALLAPTEMPLATEt;
                
                obj.DAT.StimBlocks(i).([chan,'_APTemplate_unzeroed']) = ALLAPTEMPLATE_unzeroedt;
                obj.DAT.StimBlocks(i).([chan,'_APTemplate']) = ALLAPTEMPLATEt;
                
                if (i > 1) && (obj.DAT.StimBlocks(i).NStim > 100) 
                    mc = max(ALLAPTEMPLATEt.AVG);
                    mp = max(obj.DAT.StimBlocks(i-1).([chan,'_APTemplate']).AVG);
                    if mp/mc > 5
                        obj.DAT.StimBlocks(i).([chan,'_APTemplate_unzeroed']) = obj.DAT.StimBlocks(i-1).([chan,'_APTemplate_unzeroed']);
                        obj.DAT.StimBlocks(i).([chan,'_APTemplate']) = obj.DAT.StimBlocks(i-1).([chan,'_APTemplate']);
                    end
                end
                if ~isempty(handles)
                    handles.PBarObj.Position(3) = (i)/NBLOCKS*1000;
                    handles.PBarTxt.String = [num2str(round((i)/NBLOCKS*100)),'%'];
                    drawnow
                end
                
            end
            if ~isempty(handles)
                handles.ProgressBar.Title.String = ['Finished'];
                handles.PBarObj.Position(3) = 0;
                handles.PBarTxt.String = 'Waiting';
                drawnow
            end
            id = obj.DAT.GetChan(chan);
            obj.DAT.CHANS{id.Num}.ALLAPTEMPLATE = ALLAPTEMPLATE;
            obj.DAT.CHANS{id.Num}.ALLAPTEMPLATE_unzeroed = ALLAPTEMPLATE_unzeroed;
        end

        function [APTEMPLATE_unzeroed,APTEMPLATE] = GetAPTemplate(obj, APC, Block_idx, chan)
            DAT = obj.DAT; % For convenience, copy over.
            [DPRE,DSTIM,StimIdx] = DAT.GetBlockData(Block_idx, chan);
            % [APTEMPLATE_unzeroed,APTEMPLATE] = APC.FindAPs(DPRE);
            [APTEMPLATE_unzeroed,APTEMPLATE] = APC.SimpleFindAPs(DPRE);
        end

    end
end