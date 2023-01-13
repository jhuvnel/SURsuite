classdef ZapQuick_class <handle
	
	properties
		DAT
		Values_zapped
        APLoc
               
	end
	
	methods
		% Constructor takes a SpikeDATA object.
		function obj=ZapQuick_class(SpikeDATA_object)
			obj.DAT = SpikeDATA_object;
			if ~isprop(obj.DAT, 'SamplesPerSec')
				obj.DAT.addprop('SamplesPerSec');
			end
			obj.DAT.SamplesPerSec = obj.DAT.TicksPerSec / obj.DAT.TicksPerSample;
			obj.Values_zapped = obj.DAT.Values;
            obj.APLoc = zeros(length(obj.DAT.Values),1);
        end
		
        function ClearArtifact(obj, chan, idxss, app)
            if isempty(app)
                D = obj.DAT;
            else
                D = app.DAT;
            end
            sz = size(D.StimBlocks(1).([chan,'_ArtifactTemplate_original']));
            for idxs = idxss
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_unzeroed_used']) = zeros(sz);
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = zeros(sz);
                
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_unzeroed_All']) = zeros(sz);
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = zeros(sz);
                
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_unzeroed_t']) = zeros(sz);
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_t']) = zeros(sz);
                
                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 0;

            end
        end
        
        function SubtractArtifact(obj, chan, idxss, APC, TC, app, debugFlg, showFindFlg)
            if isempty(app)
                D = obj.DAT;
                ffWait = waitbar(0,'Please wait...');
                debugPlots = [];
            else
                D = app.DAT;
                debugPlots = struct();
                app.ProgressBar.Title.String = 'Please Wait...';
                debugPlots.apPoints = subtightplot(3,5,[1:2],[0.025 0.000000000005],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.apPoints.Title.String = 'AP Points';
                debugPlots.artPoints = subtightplot(3,5,[6],[0.025 0.000000000005],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.artPoints.Title.String = 'Spline Artifact Points';
                debugPlots.ReartPoints = subtightplot(3,5,[7],[0.025 0.000000000005],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.ReartPoints.Title.String = 'ReSize Artifact Points';
                debugPlots.NoResponsePlot = subtightplot(3,5,11,[0.025 0.000000000005],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.NoResponsePlot.Title.String = 'No Response';
                debugPlots.ResponsePlot = subtightplot(3,5,12,[0.025 0.000000000005],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.ResponsePlot.Title.String = 'Response';
                debugPlots.first = subtightplot(3,5,[3:5],[0.025 0.05],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.first.Title.String = 'Spline Template';
                debugPlots.second = subtightplot(3,5,[8:10],[0.025 0.05],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.second.Title.String = 'Binned Template';
                debugPlots.third = subtightplot(3,5,[13:15],[0.025 0.05],[0.10 0.075],[.05 .05],'Parent', app.AnalysisPlotsTab);
                debugPlots.third.Title.String = 'Resized Template';
                drawnow
            end

            
            li = length(idxss);
            nns = 1;
            for idxs = idxss
                Name = ['Stim ',num2str(D.StimBlocks(idxs).Stim_E),' Ref ',num2str(D.StimBlocks(idxs).Ref_E), ' Current ', num2str(D.StimBlocks(idxs).Current_uA)];
                if isempty(app)
                    waitbar((nns-1)/li,ffWait,['Subtracting Artifact From: ',Name])
                else
                    app.ProgressBar.Title.String = ['Subtracting Artifact From: ',Name];
                    app.idxs = idxs;
                    
                end
                if D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed'])
                    
                else
                    %D.StimBlocks(idxs).AllArtTemps = [];
                    %[D.StimBlocks(idxs).Stim_E D.StimBlocks(idxs).Ref_E D.StimBlocks(idxs).Current_uA]
                    ALLSTIM = obj.GetAllStim(idxs, chan);
                    ALLSTIMPlot = reshape(ALLSTIM',1,198000);
                    mvavg = movmean(ALLSTIMPlot,1000);
                    mvavg2 = mvavg-(mvavg(1)-ALLSTIMPlot(1))*.7;
                    ASP = ALLSTIMPlot;
                    ASP(1:197400) = ASP(1:197400)-mvavg2(1:197400); %Subtracting DC Offset
                    ASP(197401:end) = ASP(197401:end)+(ASP(197400)-ASP(197401));
                    ALLSTIM = reshape(ASP,1000,198)';
                    D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM;
                    app.ProgressBar.Title.String = ['Subtracting Artifact From: ',Name, ' - Getting AP Info'];
                    [apSamps, apTime, startPt, endPt, peakDePol1, peakDePol2, maxHyperpol] = APC.findAPWidth(D.StimBlocks(idxs).([chan,'_APTemplate']), chan, 1, debugPlots);
                    if 1
                        if ~isempty(debugPlots)
                            debugPlots.apPoints.Title.String = 'AP Points';
                            drawnow
                        else
                            bbb = gcf;
                        end
                    end
                    D.StimBlocks(idxs).([chan,'_APTemplate']).apSamps = apSamps;
                    D.StimBlocks(idxs).([chan,'_APTemplate']).apTime = apTime;
                    D.StimBlocks(idxs).([chan,'_APTemplate']).startPt = startPt;
                    D.StimBlocks(idxs).([chan,'_APTemplate']).endPt = endPt;
                    D.StimBlocks(idxs).([chan,'_APTemplate']).peakDePol1 = peakDePol1;
                    D.StimBlocks(idxs).([chan,'_APTemplate']).peakDePol2 = peakDePol2;
                    D.StimBlocks(idxs).([chan,'_APTemplate']).maxHyperpol = maxHyperpol;
                    app.ProgressBar.Title.String = ['Subtracting Artifact From: ',Name, ' - Getting Spline Template'];
                    [a,b,c, TemplateStyle, Positivepks, cornerpt] = TC.SplineArtifactTemplate(D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']), apSamps, D.StimBlocks(idxs).Current_uA, idxs, 1, chan, debugPlots);
                    
                    if 1
                        if ~isempty(debugPlots)
                            debugPlots.artPoints.Title.String = 'Spline Artifact Points';
                            drawnow
                        else
                            aa = gcf;
                        end
                    end
                    temp2Use = c;
                    
                    if ismember({chan},'NerveCon')
                        [tss, tss2] = TC.InternalCompare(idxs, debugFlg, chan);
                        fss = fields(tss2);
                        maxvF = [];
                        szvF = [];
                        for i = 1:length(fss)
                            maxvF = [maxvF std(tss2.(fss{i}))];
                            szv = size(tss2.(fss{i}));
                            szvF = [szvF szv(1)];
                        end
                        sttemp = find(maxvF>200);
                        sizChk = find((szvF(sttemp)>20) & (szvF(sttemp)<190));
                        if isempty(sizChk)
                        elseif length(sizChk) > 1
                            %[stdVal,toUsenum] = max(maxvF);
                            [ol,oi]=max(szvF(sttemp(sizChk)));
                            toUsenum = sttemp(oi);
                            toUse = tss.(['t',num2str(toUsenum)]);
                            artTemp = mean(toUse,1);
                            [TemplateStyle, Positivepks, cornerpt] = TC.GetArtTempDets(idxs, chan, artTemp, apSamps);
                            temp2Use = artTemp;
                        else
                            toUsenum = sttemp(sizChk);
                            toUse = tss.(['t',num2str(toUsenum)]);
                            artTemp = mean(toUse,1);
                            [TemplateStyle, Positivepks, cornerpt] = TC.GetArtTempDets(idxs, chan, artTemp, apSamps);
                            temp2Use = artTemp;
                        end
                    end
                    
                    %temp2Use = D.DAGAN2templ(:,idxs(1));
                    [ALLSTIM, ZAPPED] = obj.ZapBlock(idxs,temp2Use,chan, Positivepks, cornerpt);
                    D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED;
                    D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM;
                    ALLSTIMPlot = reshape(ALLSTIM',1,198000);
                    ZAPPEDPlot = reshape(ZAPPED',1,198000);
                    if 1
                        if ~isempty(debugPlots)
                            plot(debugPlots.first,ALLSTIMPlot)
                            hold(debugPlots.first, 'on')
                            plot(debugPlots.first,ZAPPEDPlot)
                            hold(debugPlots.first, 'off')
                            debugPlots.first.Title.String = 'Spline Template';
                            debugPlots.first.YLim = [-0.5 0.5];
                            drawnow
                            
                        else
                            a = figure;
                            plot(ALLSTIMPlot)
                            hold on
                            plot(ZAPPEDPlot)
                            hold off
                            a.Position = [1          41        1920         963];
                        end
                    end
                    obj.DAT = D;
                    app.ProgressBar.Title.String = ['Subtracting Artifact From: ',Name, ' - Getting Binned Template'];
                    drawnow
                    [response, noresponse] = obj.BinnedTemplate(idxs, chan, ZAPPEDPlot, ALLSTIM, debugPlots, showFindFlg);
                    
                    
                    [ALLSTIM2, ZAPPED2] = obj.ZapBlock(idxs,mean(noresponse)',chan); % 6a.
                    ALLSTIMPlot2 = reshape(ALLSTIM2',1,198000);
                    ZAPPEDPlot2 = reshape(ZAPPED2',1,198000);
                    if 1
                        if ~isempty(debugPlots)
                            plot(debugPlots.second,ALLSTIMPlot2)
                            hold(debugPlots.second,'on')
                            plot(debugPlots.second,ZAPPEDPlot2)
                            hold(debugPlots.second,'off')
                            debugPlots.second.Title.String = 'Binned Template';
                            debugPlots.second.YLim = [-0.5 0.5];
                            drawnow
%                             uiwait
                        else
                            c=figure;
                            plot(ALLSTIMPlot2)
                            hold on
                            plot(ZAPPEDPlot2)
                            hold off
                            c.Position = [1          41        1920         963];
%                             uiwait
%                             close(a)
%                             close(b)
%                             close(aa)
%                             close(bbb)
                        end
                    end
                    

                    if (D.StimBlocks(idxs-1).([chan,'_ArtifactTemplate_confirmed']) && (D.StimBlocks(idxs).Current_uA ~= 20))
                        app.ProgressBar.Title.String = ['Subtracting Artifact From: ',Name, ' - Getting Resized Template'];
                        o = D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']);
                        if TemplateStyle == 1
                            
                            [mv, ml] = max(o);
                            [minv, minl] = min(o);
                            pe = find(o==mv,1,'last')+1;
                            ps = find(o==mv,1,'first')-1;
                            pme = find(o==minv,1,'last')+1;
                            pms = find(o==minv,1,'first')-1;
                            [wt,w2t] = islocalmin(o(pe:300));
                            ll2t = find(w2t>0,1,'first');
                            %[ll,ll2]=max(w2);
                            do = [diff(o) 0];
                            ll2 = find(do(pe:300)>-0.01,1,'first');
                            if ll2t-ll2<5
                                ll2 = ll2t;
                            end
                        elseif TemplateStyle == 2
                            [TF, P] = islocalmax(o);
                            pnz = (P~=0);
                            ptemp = find(P>mean(P(pnz)));
                            
                            if length(Positivepks)<2
                                pe = Positivepks(1)+1;
                            else
                                pe = Positivepks(2)+1;
                            end
                            [wt,w2t] = islocalmin(o(pe:300));
                            ll2t = find(w2t>0,1,'first');
                            %[ll,ll2]=max(w2);
                            do = [diff(o) 0];
                            ll2 = find(do(pe:300)>-0.01,1,'first');
                            if ll2t-ll2<5
                                ll2 = ll2t;
                            end
                        else
                            pe = Positivepks(1)+1;
                            [wt,w2t] = islocalmin(o(pe:300));
                            ll2 = find(w2t>0,1,'first');
                            
                        end
                        if all(D.StimBlocks(idxs).([chan,'_ArtifactTemplate_t'])==0)
                            temp2use2 = D.StimBlocks(idxs-1).([chan,'_ArtifactTemplate_used']);
                        else
                            temp2use2 = D.StimBlocks(idxs).([chan,'_ArtifactTemplate_t']);
                        end
                        
                        v1 = max(o);
                        v2 = max(temp2use2);
                        mfactor = 1;%v1/v2;
                        t = temp2use2*mfactor;
                        
                        
                        %                 [wa,w2a] = islocalmin(temp2use(pe:100));
                        %                 [lla,ll2a]=max(w2a);
                        mm1 = o(ll2+pe-1);
                        mm2 = temp2use2(ll2+pe-1);
                        if TemplateStyle == 3
                            tt2 = temp2use2*((mm1/mm2));
                            tt2(1:ll2+pe-3) = o(1:ll2+pe-3);
                        elseif ismember({chan},'NerveCon')
                            ptt=ll2+pe-1;
                            if ptt-50 < 1
                                mm1 = o(1:ptt);
                                mm2=temp2use2(1:ptt);
                            else
                                mm1 = o(ptt-50:ptt);
                                mm2=temp2use2(ptt-50:ptt);
                            end
                            mdm = mm1./mm2;
                            %                                 mdm(mdm>4) = [];
                            mdm(mdm<1) = [];
                            mdm(mdm>mean(mdm)) = [];
                            if isempty(mdm)
                                tt2=temp2use2;
                            else
                                tt2=temp2use2*(mean(mdm));
                            end
                            [xi,yi] = polyxpoly(1:1000,tt2,1:1000,o);
                            xi = round(xi);
                            cnrPt2use = xi(find(xi>pe,1,'first'));
                            if cnrPt2use > ptt
                                tt2(1:ptt-15) = o(1:ptt-15);
                            else
                                tt2(1:cnrPt2use) = o(1:cnrPt2use);
                            end
                        else
                            tt2 = temp2use2*((mm1/mm2)-0.1);
                            tt2(1:ll2+pe-3) = o(1:ll2+pe-3);
                        end
                        %             mindif = 1;
                        %             tt2(1:mindif) = [];
                        %             tt2 = [tt2 repmat(tt2(end),1,mindif)];
                        %                             tt2(1:ll2+pe-3) = o(1:ll2+pe-3);
                        
                        t = tt2;
                        if 1
                            if ~isempty(debugPlots)
                                plot(debugPlots.ReartPoints,D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']))
                                hold(debugPlots.ReartPoints,'on')
                                plot(debugPlots.ReartPoints, temp2use2)
                                plot(debugPlots.ReartPoints, tt2)
                                plot(debugPlots.ReartPoints, 1:200,abs(o(1:200)-tt2(1:200)))
                                plot(debugPlots.ReartPoints, ll2+pe-1,o(ll2+pe-1),'b*')
                                hold(debugPlots.ReartPoints,'off')
                                debugPlots.ReartPoints.Title.String = 'ReSized Template Points';
                                drawnow
                                
                            else
                                dd = figure;
                                plot(D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']))
                                hold on
                                plot(temp2use2)
                                plot(tt2)
                                plot(1:200,abs(o(1:200)-tt2(1:200)))
                                plot(ll2+pe-1,o(ll2+pe-1),'b*')
                                hold off
                                dd.Position = [1          41        1920         963];
                            end
                        end
                        
                        [ALLSTIM3, ZAPPED3] = obj.ZapBlock(idxs,t',chan);
                        ALLSTIMPlot3 = reshape(ALLSTIM3',1,198000);
                        ZAPPEDPlot3 = reshape(ZAPPED3',1,198000);
                        if 1
                            if ~isempty(debugPlots)
                                plot(debugPlots.third,ALLSTIMPlot3)
                                hold(debugPlots.third, 'on')
                                plot(debugPlots.third,ZAPPEDPlot3)
                                hold(debugPlots.third, 'off')
                                debugPlots.third.Title.String = 'Resized Template';
                                debugPlots.third.YLim = [-0.5 0.5];
                                drawnow
                            else
                                e = figure;
                                plot(ALLSTIMPlot3)
                                hold on
                                plot(ZAPPEDPlot3)
                                hold off
                                e.Position = [1          41        1920         963];
                                
                                uiwait
                                close(a)
                                close(b)
                                close(dd)
                                close(aa)
                                close(bbb)
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED3;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM3;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All'])  = t;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = t;
                            end
                        end                
                    end
                    app.figure1.WindowKeyPressFcn = @(src,event)app.KeyDownFcn;
                    app.figure1.WindowButtonDownFcn = @(src,event)app.ArtifactSubMouseDown;
                    app.figure1.WindowButtonUpFcn = @(src,event)app.ArtifactSubMouseUp;
                    app.debugPlots = debugPlots;
                    app.TemplateStyle = TemplateStyle;
                    app.Positivepks = Positivepks;
                    app.ProgressBar.Title.String = ['Subtracting Artifact From: ',Name, ' - PLEASE CHOOSE TEMPLATE TO USE'];
                    uiwait
                    app.KeyPressed = 0;
                    app.figure1.WindowKeyPressFcn = [];
                    app.figure1.WindowButtonDownFcn = [];
                    app.figure1.WindowButtonUpFcn = [];
                    if app.stopButton.Value
                        break
                    end
                    if ~isempty(debugPlots)
                        if app.SplineTempButton.Value
                            if isempty(app.ALLSTIM)
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = temp2Use;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = temp2Use;
                            else
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = app.ZAPPED;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = app.ALLSTIM;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = app.temp2Use;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = app.temp2Use;
                            end
                            D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
                            
                            app.SplineTempButton.Value = 0;
                        elseif app.BinTempButton.Value
                            if isempty(app.ALLSTIM2)
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED2;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM2;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = noresponse;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = mean(noresponse);
                            else
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = app.ZAPPED2;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = app.ALLSTIM2;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = app.noresponse;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = mean(app.noresponse);
                            end
                            D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
                            
                            app.BinTempButton.Value = 0;
                        elseif app.ReSizeTempButton.Value
                            if isempty(app.ALLSTIM3)
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED3;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM3;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All'])  = t;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = t;
                            else
                                D.StimBlocks(idxs).([chan,'_ZAPPED']) = app.ZAPPED3;
                                D.StimBlocks(idxs).([chan,'_ALLSTIM']) = app.ALLSTIM3;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All'])  = app.temp2Use2;
                                D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = app.temp2Use2;
                            end
                            app.ReSizeTempButton.Value = 0;
                            %D.DAGAN2templ_used(:,idxs) = t;
                            D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
                            
                            fullap = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:maxHyperpol+10) - mean(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW);
                            if D.StimBlocks(idxs).Current_uA ~= 200
                                fornextTemp = [];
                                BLOCK = D.StimBlocks(idxs);
                                STIMS = find((D.Stim_ticks >= [BLOCK.Start_tick]) ...
                                    & (D.Stim_ticks <= [BLOCK.End_tick]));
                                ticks = D.Stim_ticks(STIMS)./2;
                                ticks = ticks - (ticks(1)-1);
                                for pos = 1:length(ticks)
                                    if pos == length(ticks)
                                        bound = ticks(pos):length(ALLSTIMPlot);
                                    else
                                        bound = ticks(pos):(ticks(pos+1)-1);
                                    end
                                    if ~any(isnan(ZAPPEDPlot2))
                                        u = ZAPPEDPlot2(bound);% - mean(ZAPPEDPlot2);
                                    else
                                        u = ZAPPEDPlot3(bound);
                                    end
                                    u = u-mean(u);
                                    %                                             figure
                                    %                                             findsignal(u,fullap);
                                    %                                             abc = gcf;
                                    %                                             abc.Position = [1          41        1920         963];
                                    %                                             close(abc)
                                    [istart,istop,dist] = findsignal(u,fullap);
                                    [maxv,maxl]=max(u(40:end));
                                    maxl = maxl +39;
                                    %                                     maxls = [maxls maxl];
                                    if (maxl > istart) && (maxl < istop)
                                        if maxl < 300
                                            u2 = ALLSTIMPlot2(bound);
                                            fAP = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt)';
                                            u2(istart:(istart+length(fAP)-1)) = u2(istart:(istart+length(fAP)-1)) - fAP;
                                            fornextTemp = [fornextTemp; u2];
                                        end
                                    else
                                        maxl;
                                        maxv/max(fullap);
                                        if (maxl < 85) && (maxv/max(fullap) < 0.7)
                                            %noresponse Maybe add saving this stim
                                            
                                        elseif maxl > 150
                                            %noresponse
                                        else
                                            %response = [response; ALLSTIMPlot(bound)];
                                            %                             figure
                                            %                             findsignal(u(1:250),fullap);
                                            %                             abc = gcf;
                                            %                             abc.Position = [1          41        1920         963];
                                            %                             close(abc)
                                            [istart,istop,dist] = findsignal(u(1:250),fullap);
                                            u2 = ALLSTIMPlot2(bound);
                                            fAP = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt)';
                                            u2(istart:(istart+length(fAP)-1)) = u2(istart:(istart+length(fAP)-1)) - fAP;
                                            fornextTemp = [fornextTemp; u2];
                                        end
                                    end
                                end
                                D.StimBlocks(idxs+1).([chan,'_ArtifactTemplate_t']) = mean(fornextTemp);
                            end
                        end
                        app.ALLSTIM = [];
                        app.ZAPPED = [];
                        app.ALLSTIM2 = [];
                        app.ZAPPED2 = [];
                        app.ALLSTIM3 = [];
                        app.ZAPPED3 = [];
                        app.temp2Use = [];
                        app.temp2Use2 = [];
                        app.noresponse = [];

                        cla(debugPlots.apPoints)
                        cla(debugPlots.artPoints)
                        cla(debugPlots.ReartPoints)
                        cla(debugPlots.NoResponsePlot)
                        cla(debugPlots.ResponsePlot)
                        cla(debugPlots.first)
                        cla(debugPlots.second)
                        cla(debugPlots.third)
                        
                    else
                        D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED2;
                        D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM2;
                        D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = noresponse;
                        D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = mean(noresponse);
                        D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
                    end
                    
                    
%                     if (sz(1) >= 20) && (sz(1) ~= 198) % 5a.
%                                                 
% %                         [ALLSTIM2, ZAPPED2] = obj.ZapBlock(idxs,mean(noresponse)',chan); % 6a.
% %                         ALLSTIMPlot2 = reshape(ALLSTIM2',1,198000);
% %                         ZAPPEDPlot2 = reshape(ZAPPED2',1,198000);
% %                         if 1
% %                             if ~isempty(debugPlots)
% %                                 plot(debugPlots.second,ALLSTIMPlot2)
% %                                 hold(debugPlots.second,'on')
% %                                 plot(debugPlots.second,ZAPPEDPlot2)
% %                                 hold(debugPlots.second,'off')
% %                                 debugPlots.second.Title.String = 'Binned Template';
% %                                 drawnow
% %                                 uiwait
% %                             else
% %                                 c=figure;
% %                                 plot(ALLSTIMPlot2)
% %                                 hold on
% %                                 plot(ZAPPEDPlot2)
% %                                 hold off
% %                                 c.Position = [1          41        1920         963];
% %                                 uiwait
% %                                 close(a)
% %                                 close(b)
% %                                 close(aa)
% %                                 close(bbb)
% %                             end
% %                         end
% %                         
% %                         if ~isempty(debugPlots)
% %                             if app.SplineTempButton.Value
% %                                 D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED;
% %                                 D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM;
% %                                 D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = temp2Use;
% %                                 D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = temp2Use;
% %                                 D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
% %                                 
% %                                 app.SplineTempButton.Value = 0;
% %                             elseif app.BinTempButton.Value
% %                                 D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED2;
% %                                 D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM2;
% %                                 D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = noresponse;
% %                                 D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = mean(noresponse);
% %                                 D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
% %                                 
% %                                 app.BinTempButton.Value = 0;
% %                             end
% %                             cla(debugPlots.apPoints)
% %                             cla(debugPlots.artPoints)
% %                             cla(debugPlots.ReartPoints)
% %                             cla(debugPlots.NoResponsePlot)
% %                             cla(debugPlots.ResponsePlot)
% %                             cla(debugPlots.first)
% %                             cla(debugPlots.second)
% %                             cla(debugPlots.third)
% %                             
% %                         else
% %                             D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED2;
% %                             D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM2;
% %                             D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All']) = noresponse;
% %                             D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = mean(noresponse);
% %                             D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
% %                         end
%                         
%                         
%                     else
%                         if D.StimBlocks(idxs-1).([chan,'_ArtifactTemplate_confirmed'])
%                             o = D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']);
%                             if TemplateStyle == 1
%                                 
%                                 [mv, ml] = max(o);
%                                 [minv, minl] = min(o);
%                                 pe = find(o==mv,1,'last')+1;
%                                 ps = find(o==mv,1,'first')-1;
%                                 pme = find(o==minv,1,'last')+1;
%                                 pms = find(o==minv,1,'first')-1;
%                                 [wt,w2t] = islocalmin(o(pe:300));
%                                 ll2t = find(w2t>0,1,'first');
%                                 %[ll,ll2]=max(w2);
%                                 do = [diff(o) 0];
%                                 ll2 = find(do(pe:300)>-0.01,1,'first');
%                                 if ll2t-ll2<5
%                                     ll2 = ll2t;
%                                 end
%                             elseif TemplateStyle == 2
%                                 [TF, P] = islocalmax(o);
%                                 pnz = (P~=0);
%                                 ptemp = find(P>mean(P(pnz)));
%                                 
%                                 %                     if length(ptemp) > 2
%                                 %                         mv = max(o);
%                                 %                         compare2max = o(ptemp)/mv;
%                                 %                         ptemp(compare2max<0.5) = [];
%                                 %                         flatM = find(o==o(ptemp(1)));
%                                 %                         pmax1 = flatM(end);
%                                 %
%                                 %                         sp = find(ptemp>pmax1,1);
%                                 %                         flatM2 = find(o==o(ptemp(sp)));
%                                 %                         pmax2 = flatM2(end);
%                                 %                         pe = pmax2+1;
%                                 %                     else
%                                 %                         pmax1 = ptemp(1);
%                                 %                         pmax2 = ptemp(2);
%                                 %                         pe = pmax2+1;
%                                 %                     end
%                                 pe = Positivepks(2)+1;
%                                 [wt,w2t] = islocalmin(o(pe:300));
%                                 ll2t = find(w2t>0,1,'first');
%                                 %[ll,ll2]=max(w2);
%                                 do = [diff(o) 0];
%                                 ll2 = find(do(pe:300)>-0.01,1,'first');
%                                 if ll2t-ll2<5
%                                     ll2 = ll2t;
%                                 end
%                             else
%                                 pe = Positivepks(1)+1;
%                                 [wt,w2t] = islocalmin(o(pe:300));
%                                 ll2 = find(w2t>0,1,'first');
%                                 
%                             end
%                             if all(D.StimBlocks(idxs).([chan,'_ArtifactTemplate_t'])==0)
%                                 temp2use = D.StimBlocks(idxs-1).([chan,'_ArtifactTemplate_used']);
%                             else
%                                 temp2use = D.StimBlocks(idxs).([chan,'_ArtifactTemplate_t']);
%                             end
%                             
%                             v1 = max(o);
%                             v2 = max(temp2use);
%                             mfactor = 1;%v1/v2;
%                             t = temp2use*mfactor;
%                             
%                             
%                             %                 [wa,w2a] = islocalmin(temp2use(pe:100));
%                             %                 [lla,ll2a]=max(w2a);
%                             mm1 = o(ll2+pe-1);
%                             mm2 = temp2use(ll2+pe-1);
%                             if TemplateStyle == 3
%                                 tt2 = temp2use*((mm1/mm2));
%                                 tt2(1:ll2+pe-3) = o(1:ll2+pe-3);
%                             elseif ismember({chan},'NerveCon')
%                                 ptt=ll2+pe-1;
%                                 if ptt-50 < 1
%                                     mm1 = o(1:ptt);
%                                     mm2=temp2use(1:ptt);
%                                 else
%                                     mm1 = o(ptt-50:ptt);
%                                     mm2=temp2use(ptt-50:ptt);
%                                 end
%                                 mdm = mm1./mm2;
% %                                 mdm(mdm>4) = [];
%                                 mdm(mdm<1) = [];
%                                 mdm(mdm>mean(mdm)) = [];
%                                 tt2=temp2use*(mean(mdm));
%                                 [xi,yi] = polyxpoly(1:1000,tt2,1:1000,o);
%                                 xi = round(xi);
%                                 cnrPt2use = xi(find(xi>pe,1,'first'));
%                                 if cnrPt2use > ptt
%                                     tt2(1:ptt-15) = o(1:ptt-15);
%                                 else
%                                 tt2(1:cnrPt2use) = o(1:cnrPt2use);
%                                 end
%                             else
%                                 tt2 = temp2use*((mm1/mm2)-0.1);
%                                 tt2(1:ll2+pe-3) = o(1:ll2+pe-3);
%                             end
%                             %             mindif = 1;
%                             %             tt2(1:mindif) = [];
%                             %             tt2 = [tt2 repmat(tt2(end),1,mindif)];
% %                             tt2(1:ll2+pe-3) = o(1:ll2+pe-3);
% 
%                             t = tt2;
%                             if 1
%                                 if ~isempty(debugPlots)
%                                     plot(debugPlots.ReartPlots,D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']))
%                                     hold(debugPlots.ReartPlots,'on')
%                                     plot(debugPlots.ReartPlots, temp2use)
%                                     plot(debugPlots.ReartPlots, tt2)
%                                     plot(debugPlots.ReartPlots, 1:200,abs(o(1:200)-tt2(1:200)))
%                                     plot(debugPlots.ReartPlots, ll2+pe-1,o(ll2+pe-1),'b*')
%                                     hold(debugPlots.ReartPlots,'off')
%                                     debugPlots.ReartPlots.Title.String = 'ReSized Template Points';
%                                     drawnow
%                                     
%                                 else
%                                     dd = figure;
%                                     plot(D.StimBlocks(idxs).([chan,'_ArtifactTemplate_original']))
%                                     hold on
%                                     plot(temp2use)
%                                     plot(tt2)
%                                     plot(1:200,abs(o(1:200)-tt2(1:200)))
%                                     plot(ll2+pe-1,o(ll2+pe-1),'b*')
%                                     hold off
%                                     dd.Position = [1          41        1920         963];
%                                 end
%                             end
%                             
%                             [ALLSTIM2, ZAPPED2] = obj.ZapBlock(idxs,t',chan);
%                             ALLSTIMPlot2 = reshape(ALLSTIM2',1,198000);
%                             ZAPPEDPlot2 = reshape(ZAPPED2',1,198000);
%                             if 1
%                                 if ~isempty(debugPlots)
%                                     plot(debugPlots.third,ALLSTIMPlots2)
%                                     hold(debugPlots.third, 'on')
%                                     plot(debugPlots.third,ZAPPEDPlots2)
%                                     hold(debugPlots.third, 'off')
%                                     debugPlots.third.Title.String = 'Resized Template';
%                                     drawnow
%                                     uiwait
%                                 else
%                                     e = figure;
%                                     plot(ALLSTIMPlot2)
%                                     hold on
%                                     plot(ZAPPEDPlot2)
%                                     hold off
%                                     e.Position = [1          41        1920         963];
%                                     
%                                     uiwait
%                                     close(a)
%                                     close(b)
%                                     close(dd)
%                                     close(aa)
%                                     close(bbb)
%                                     D.StimBlocks(idxs).([chan,'_ZAPPED']) = ZAPPED2;
%                                     D.StimBlocks(idxs).([chan,'_ALLSTIM']) = ALLSTIM2;
%                                     D.StimBlocks(idxs).([chan,'_ArtifactTemplate_All'])  = t;
%                                     D.StimBlocks(idxs).([chan,'_ArtifactTemplate_used']) = t;
%                                 end
%                             end
%                             
%                             
%                             %D.DAGAN2templ_used(:,idxs) = t;
%                             D.StimBlocks(idxs).([chan,'_ArtifactTemplate_confirmed']) = 1;
%                             fullap = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:maxHyperpol+10) - mean(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW);
%                             if D.StimBlocks(idxs).Current_uA ~= 200
%                                 fornextTemp = [];
%                                 for pos = 1:length(ticks)
%                                     if pos == length(ticks)
%                                         bound = ticks(pos):length(ALLSTIMPlot);
%                                     else
%                                         bound = ticks(pos):(ticks(pos+1)-1);
%                                     end
%                                     u = ZAPPEDPlot2(bound);% - mean(ZAPPEDPlot2);
%                                     u = u-mean(u);
%                                     %                                             figure
%                                     %                                             findsignal(u,fullap);
%                                     %                                             abc = gcf;
%                                     %                                             abc.Position = [1          41        1920         963];
%                                     %                                             close(abc)
%                                     [istart,istop,dist] = findsignal(u,fullap);
%                                     [maxv,maxl]=max(u(40:end));
%                                     maxl = maxl +39;
% %                                     maxls = [maxls maxl];
%                                     if (maxl > istart) && (maxl < istop)
%                                         if maxl < 300
%                                             u2 = ALLSTIMPlot2(bound);
%                                             fAP = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt)';
%                                             u2(istart:(istart+length(fAP)-1)) = u2(istart:(istart+length(fAP)-1)) - fAP;
%                                             fornextTemp = [fornextTemp; u2];
%                                         end
%                                     else
%                                         maxl;
%                                         maxv/max(fullap);
%                                         if (maxl < 85) && (maxv/max(fullap) < 0.7)
%                                             %noresponse Maybe add saving this stim
%                                             
%                                         elseif maxl > 150
%                                             %noresponse
%                                         else
%                                             %response = [response; ALLSTIMPlot(bound)];
%                                             %                             figure
%                                             %                             findsignal(u(1:250),fullap);
%                                             %                             abc = gcf;
%                                             %                             abc.Position = [1          41        1920         963];
%                                             %                             close(abc)
%                                             [istart,istop,dist] = findsignal(u(1:250),fullap);
%                                             u2 = ALLSTIMPlot2(bound);
%                                             fAP = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt)';
%                                             u2(istart:(istart+length(fAP)-1)) = u2(istart:(istart+length(fAP)-1)) - fAP;
%                                             fornextTemp = [fornextTemp; u2];
%                                         end
%                                     end
%                                 end
%                                 D.StimBlocks(idxs+1).([chan,'_ArtifactTemplate_t']) = mean(fornextTemp);
%                             end
%                         end
%                     end
                end

                if isempty(app)
                    waitbar((nns)/li,ffWait,['Subtracting Artifact From: ',Name])
                else
                    app.PBarObj.Position(3) = (nns)/li*1000;
                    app.PBarTxt.String = [num2str(round((nns)/li*100)),'%'];
                    drawnow
                end
                obj.DAT = D;
                nns = nns +1;
            end
            
            if isempty(app)
                waitbar((nns)/li,ffWait,['Finished'])
                close(ffWait)
            else
                app.ProgressBar.Title.String = 'Finished';
                app.PBarObj.Position(3) = (nns)/li*1000;
                app.PBarTxt.String = [num2str(round((nns)/li*100)),'%'];
                drawnow
            end
            obj.DAT = D;
        end
		
        function [response, noresponse] = BinnedTemplate(obj, idxs, chan, ZAPPEDPlot, ALLSTIM, debugPlots, showFindFlg)
            D = obj.DAT; % Just for convenience.
            peakDePol1 = D.StimBlocks(idxs).([chan,'_APTemplate']).peakDePol1;
            startPt = D.StimBlocks(idxs).([chan,'_APTemplate']).startPt;
            endPt = D.StimBlocks(idxs).([chan,'_APTemplate']).endPt;
            ALLSTIMPlot = reshape(ALLSTIM',1,198000);

            trimAPTemp = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(peakDePol1-10:peakDePol1+35);
            trimAPTemp = trimAPTemp - mean(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW);
            if ismember({chan},'NerveCon')
                trimAPTemp = D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW(startPt:endPt);
                trimAPTemp = trimAPTemp - mean(D.StimBlocks(idxs).([chan,'_APTemplate']).AVG_RAW);
            end
            BLOCK = D.StimBlocks(idxs);
            STIMS = find((D.Stim_ticks >= [BLOCK.Start_tick]) ...
                & (D.Stim_ticks <= [BLOCK.End_tick]));
            ticks = D.Stim_ticks(STIMS)./2;
            ticks = ticks - (ticks(1)-1);
            sz = 1;
            if all(D.StimBlocks(idxs).([chan,'_ArtifactTemplate_t'])==0)
                %if isfield(D.StimBlocks(idxs-1),'AllArtTemps')
                si = size(D.StimBlocks(idxs-1).([chan,'_ArtifactTemplate_All']));
                if (si(1)<40) && (si(1) ~= 0) && (D.StimBlocks(idxs).Current_uA ~= 20)
                    sz = 1;
                end
                %end

                if sz
                    response = [];
                    noresponse = [];
                    d = [];
                    %             sbs = zeros(1,198);
                    vals = [];
                    maxls = [];
                    if ismember({chan},'NerveCon')
                        [pospks,lposocs,posw,posp] = findpeaks(ZAPPEDPlot,1:length(ZAPPEDPlot),'minpeakheight',max(trimAPTemp)*.5,'minpeakdistance',900);

                        if length(pospks)>180
                            testFP = 1;
                            response = ALLSTIM;
                        else
                            testFP = 0;
                        end
                    else
                        testFP = 0;
                    end
                    if ~testFP
                        for pos = 1:length(ticks)
                            noresponseFlag = 0;
                            responseFlag = 0;
                            if pos == length(ticks)
                                bound = round(ticks(pos)):length(ALLSTIMPlot);
                            else
                                bound = round(ticks(pos)):(round(ticks(pos+1))-1);
                            end
                            if length(bound) == 999
                                bound = [bound bound(end)+1];
                            elseif length(bound) > 1000
                                bound(1) = [];
                            end
                            newALLSTIMPlot = ZAPPEDPlot(bound);%ALLSTIMPlot(bound)-D.DAGAN2templ(:,idxs)';% 1.
                            newALLSTIMPlot = newALLSTIMPlot - mean(newALLSTIMPlot);

                            [maxv,maxl]=max(newALLSTIMPlot(40:end));
                            maxl = maxl +39;

                            %                 [maxv2,maxl2]=max(newALLSTIMPlot(40:end)-mean(newALLSTIMPlot));
                            %                 maxl2 = maxl2 +39;
                            if showFindFlg
                                figure
                                if ismember({chan},'NerveCon')
                                    findsignal(newALLSTIMPlot,trimAPTemp,'Normalization','zscore','NormalizationLength',10);
                                else
                                    findsignal(newALLSTIMPlot,trimAPTemp);
                                end

                                abc = gcf;
                                abc.Position = [1          41        1920         963];
                                hold on
                                plot(maxl,newALLSTIMPlot(maxl),'r*')
                                hold off
                            end

                            if ismember({chan},'NerveCon')
                                [istart,istop,dist] = findsignal(newALLSTIMPlot,trimAPTemp,'Normalization','zscore','NormalizationLength',10);
                                coLim = 195;
                            else
                                coLim = 150;
                                [istart,istop,dist] = findsignal(newALLSTIMPlot,trimAPTemp);
                            end
                            %                 [istart2,istop2,dist2] = findsignal(newALLSTIMPlot-mean(newALLSTIMPlot),trimAPTemp); %2.
                            d = [d dist];
                            dist;
                            val = istart-1;%ticks(pos)-istart
                            %val2 = istart2-1;
                            vals = [vals val];
                            t1 = newALLSTIMPlot;
                            %                     t2 = newALLSTIMPlot-mean(newALLSTIMPlot);
                            [t1mv, t1ml] = max(t1(istart:istop));
                            t1ml = t1ml+istart-1;
                            %                     [t2mv, t2ml] = max(t2(istart:istop));

                            % 1. If the max peak and the findsignal align
                            % 2. Check max within findsignal
                            if (maxl > istart) && (maxl < istop)
                                if abs(val) > coLim %%%%%%%% used to be 110
                                    noresponseFlag = 1;
                                    responseFlag = 0;
                                else
                                    if (maxl <= coLim) && (maxv/max(trimAPTemp) < 0.6)  %used to be 130
                                        %                             if (dist2 < 0.25) && (maxv2/max(trimAPTemp)>0.6) %0.5
                                        %                                 response = [response; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))]; % 3.
                                        %                                 sbs(pos) = 1;
                                        if (maxl <115) && (maxv/max(trimAPTemp) > 0.35) && (maxv/max(trimAPTemp)< 0.6) %&& (dist >0.25)
                                            noresponseFlag = 0;
                                            responseFlag = 1;
                                        else
                                            noresponseFlag = 1;
                                            responseFlag = 0;
                                        end
                                    elseif maxl > coLim %used to be 130
                                        noresponseFlag = 1;
                                        responseFlag = 0;
                                    else
                                        noresponseFlag = 0;
                                        responseFlag = 1;
                                    end
                                    %                     if dist > 0.5
                                    %                         noresponse = [noresponse; ALLSTIMPlot(bound)]; % 4.
                                    %                     else
                                    %                         response = [response; ALLSTIMPlot(bound)]; % 3.
                                    %                     end
                                end
                            else
                                %                     if (dist2 < 0.25) && (maxv2/max(trimAPTemp)>0.6) %0.5
                                %                                 response = [response; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))]; % 3.
                                %                                 sbs(pos) = 1;

                                if (t1ml <coLim) && (t1mv/max(trimAPTemp) > 0.6)
                                    noresponseFlag = 0;
                                    responseFlag = 1;
                                elseif (t1ml <115) && (t1mv/max(trimAPTemp) > 0.35) && (t1mv/max(trimAPTemp)< 0.6)
                                    noresponseFlag = 0;
                                    responseFlag = 1;
                                elseif ismember({chan},'NerveCon')
                                    [~,~,dd]=findsignal(newALLSTIMPlot,trimAPTemp);
                                    if D.StimBlocks(idxs).Current_uA > 100
                                        [~,~,dd]=findsignal(newALLSTIMPlot,trimAPTemp(1:100));
                                        thresh = 1;
                                    else
                                        thresh = 0.1;
                                    end

                                    if dd<thresh
                                        noresponseFlag = 0;
                                        responseFlag = 1;
                                    else
                                        noresponseFlag = 1;
                                        responseFlag = 0;
                                    end
                                else
                                    noresponseFlag = 1;
                                    responseFlag = 0;
                                end

                                %                     [maxl];
                                %                     [maxv/max(trimAPTemp)];
                                %                     [dist];
                                %                     if (maxl < 130) && (maxv/max(trimAPTemp) < 0.6)
                                %                         noresponse = [noresponse; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))]; % 4.
                                %                     elseif maxl > 130
                                %                         noresponse = [noresponse; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))];
                                %                     else
                                %                         response = [response; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))]; % 3.
                                %                         sbs(pos) = 1;
                                %                     end
                            end
                            if length(bound) ~= 1000
                                stpppp = 1;
                            end
                            if noresponseFlag
                                noresponse = [noresponse; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))]; % 4.
                            end
                            if responseFlag
                                response = [response; ALLSTIMPlot(bound)-mean(ALLSTIMPlot(bound))]; % 3.
                                sbs(pos) = 1;
                            end
                            if showFindFlg
                                close(abc)
                            end

                        end
                    end

                    if 1
                        if ~isempty(debugPlots)
                            sz = size(noresponse);
                            for i = 1:sz(1)
                                plot(debugPlots.NoResponsePlot,noresponse(i,:))
                                hold(debugPlots.NoResponsePlot, 'on')
                            end
                            plot(debugPlots.NoResponsePlot,mean(noresponse),'color','k','linewidth',3)
                            debugPlots.NoResponsePlot.Title.String = ['No Response ',num2str(sz(1))] ;
                            hold(debugPlots.NoResponsePlot, 'off')
                            drawnow

                            sz = size(response);
                            for i = 1:sz(1)
                                plot(debugPlots.ResponsePlot, response(i,:))
                                hold(debugPlots.ResponsePlot,'on')
                            end
                            plot(debugPlots.ResponsePlot,mean(response),'color','k','linewidth',3)
                            debugPlots.ResponsePlot.Title.String = ['Response ',num2str(sz(1))];
                            hold(debugPlots.ResponsePlot,'off')
                            drawnow
                        else
                            b=figure;
                            sz = size(noresponse);
                            ba = subplot(1,2,1);
                            for i = 1:sz(1)
                                plot(noresponse(i,:))
                                hold on
                            end
                            plot(mean(noresponse),'color','k','linewidth',3)
                            ba.Title.String = ['No Response ',num2str(sz(1))] ;
                            hold off

                            sz = size(response);
                            bb = subplot(1,2,2);
                            for i = 1:sz(1)
                                plot(response(i,:))
                                hold on
                            end
                            plot(mean(response),'color','k','linewidth',3)
                            bb.Title.String = ['Response ',num2str(sz(1))];
                            hold off
                            b.Position = [1          41        1920         963];
                        end
                        sz = size(noresponse);
                    end
                end
            else
                if 1
                    if ~isempty(debugPlots)
                        b=[];
                    else
                        b=figure;
                    end
                end
                response = ALLSTIM;
                noresponse = [];
                sz = 0;
            end
        end
		
        
        % Zap a block of data. Copy the original data from Values, and put
		% the zapped data at the same location in Values_zapped. FOR
		% CONVENIENCE ONLY, return ALLSTIM and ZAPPED, for easy plotting of
		% the "before" and "after" overlapped stim blocks.
		function [ALLSTIM, ZAPPED] = ZapBlock(obj, Block_idx, TEMPLATE, chan, Positivepks, cornerpt)
			DAT = obj.DAT; % Just for convenience.
			
			% If TEMPLATE is not provided, then we just "self zap" using this
			% block to create the zapping template.
			if nargin < 3
				TEMPLATE = obj.GetZapTemplate(Block_idx);
			end
			
			BLOCK = DAT.StimBlocks(Block_idx);
			STIMS = find((DAT.Stim_ticks >= [BLOCK.Start_tick]) ...
				& (DAT.Stim_ticks <= [BLOCK.End_tick]));
			
			NSAMPS = length(TEMPLATE);
			ALLSTIM = zeros(length(STIMS), NSAMPS);
			
            if ~isfield(DAT.StimBlocks(Block_idx),[chan,'_ALLSTIM'])
                % Subtract the TEMPLATE out of the source obj.Values data at each
                % stim pulse in this block, and store it in obj.Values_zapped.
                for i=1:length(STIMS)
                    %                 figure
                    %                 bc = gcf;
                    tick = DAT.Stim_ticks(STIMS(i));
                    [VALS, s_idx, e_idx] = obj.DAT.ValsAtTick(tick, NSAMPS, chan);
                    %                 plot(VALS)
                    %                 hold on
                    %                 plot(TEMPLATE)
                    if isempty(VALS); continue; end
                    %obj.Values_zapped(s_idx:e_idx) = VALS - TEMPLATE(:);
                    %                 plot(obj.Values_zapped(s_idx:e_idx))
                    %                 hold off
                    %                 close(bc)
                    ALLSTIM(i,:) = VALS;
                end
                ALLSTIM = ALLSTIM - mean(ALLSTIM(:,1));
            else
                ALLSTIM = DAT.StimBlocks(Block_idx).([chan,'_ALLSTIM']);
            end
            if iscolumn(TEMPLATE)
                TEMPLATE = TEMPLATE';
            end
            sz = size(ALLSTIM);
            ZAPPED = ALLSTIM - TEMPLATE;
            for ii = 1:sz(1)
                u = ALLSTIM(ii,:);
                [aa,ab] = max(u);
                [a2,ab2]=min(u);
                sameMin = find(u(1:150)-a2<0.005);
                sameMax = find(aa-u(1:150)<0.005);
                
                if length(sameMin)>4
                    ZAPPED(ii,sameMin) = ZAPPED(ii,sameMin(1)-1);
                end
                if length(sameMax)>4
                    if (sameMax(1)-1) < 1
                    else
                        sameMax = [sameMax(1)-1 sameMax];
                    end
                    ZAPPED(ii,sameMax) = ZAPPED(ii,sameMax(end)+1);
                end
                if ~any(isnan([ZAPPED(ii,1:55)]))
                    ZAPPED(ii,1:55) = filtfilt(ones(1,9)/9,1,ZAPPED(ii,1:55));
                end
            end
			
        end
		
        function ALLSTIM = GetAllStim(obj, Block_idx, chan)
			DAT = obj.DAT; % Just for convenience.
			
			BLOCK = DAT.StimBlocks(Block_idx);
			STIMS = find((DAT.Stim_ticks >= [BLOCK.Start_tick]) ...
				& (DAT.Stim_ticks <= [BLOCK.End_tick]));
			
			NSAMPS = 1000;
			ALLSTIM = zeros(length(STIMS), NSAMPS);
			

			for i=1:length(STIMS)
				tick = DAT.Stim_ticks(STIMS(i));
				[VALS, s_idx, e_idx] = obj.DAT.ValsAtTick(tick, NSAMPS, chan);

				if isempty(VALS); continue; end
				ALLSTIM(i,:) = VALS;
			end
			ALLSTIM = ALLSTIM - mean(ALLSTIM(:,1));
        end
        
		function [Start_idx, End_idx]=BlockSampleRange(obj, Block_idx)
			Start_idx = obj.DAT.StimBlocks(Block_idx).Start_tick / obj.DAT.TicksPerSample;
			End_idx = obj.DAT.StimBlocks(Block_idx).End_tick / obj.DAT.TicksPerSample;
		end
		
		function PlotBlock(obj, Block_idx, varargin)
			
			% How much data, in seconds, to show before and after block.
			PrePlot_sec = 0;
			PostPlot_sec = 0;
			
			for k=1:2:length(varargin)
				NAME = varargin{k};
				if k+1 > length(varargin); fprintf('Value not provided for argument "%s"\n', NAME); break; end
				VAL = varargin{k+1};
				switch(varargin{k})
					case {'PrePlot_sec'}; PrePlot_sec = double(VAL);
					case {'PostPlot_sec'}; PostPlot_sec = double(VAL);
					otherwise
						fprintf('BAD named parameter: %s\n', varargin{k});
				end
			end
			
			[istart, iend] = obj.BlockSampleRange(Block_idx);
			istart = max(1, istart-PrePlot_sec*obj.DAT.SampleRate);
			iend = min(length(obj.DAT.Values), iend+PostPlot_sec*obj.DAT.SampleRate);
			
			%plot([obj.DAT.Values(istart:iend) obj.Values_zapped(istart:iend)], '.-');
			plot([obj.DAT.Values(istart:iend) obj.Values_zapped(istart:iend)]);
			
			grid on;
			zoom on;
			legend('UnZapped', 'Zapped');
		end
		
		
	end
	
end
