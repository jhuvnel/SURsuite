classdef MotionStim_class <handle
    
    properties
        DAT
    end
    
    methods
        % Constructor takes a SpikeDATA object.
        function obj = MotionStim_class(SpikeDATA_object)
            obj.DAT = SpikeDATA_object;
        end
        
        function AlignSyncSignals(obj, debugFlg, handles)
            D = obj.DAT;
            FSYNC = D.GetChan('FSYNC');
            CEDSYNC = D.GetChan('CED-SYNC');
            
            tFSYNC = 0:1/FSYNC.SampleRate:(length(FSYNC.Data)-1)/FSYNC.SampleRate;
            yFSYNC = FSYNC.Data./max(FSYNC.Data);
            
            tCEDSYNC = CEDSYNC.Data/D.TicksPerSec;
            yCEDSYNC = zeros(length(tCEDSYNC),1);
            if CEDSYNC.FirstTimestampLevel
                yCEDSYNC(2:2:end) = 1;
            else
                yCEDSYNC(1:2:end) = 1;
            end
            
            tCEDSYNC = [0;tCEDSYNC];
            yCEDSYNC = [CEDSYNC.FirstTimestampLevel; yCEDSYNC];
            
            
            sharedT = [];
            newyCEDSYNC = [];
            newyFSYNC = [];
            for i = 1:length(tCEDSYNC)
                
                if i < length(tCEDSYNC)
                    tids = find(tFSYNC>=tCEDSYNC(i) & tFSYNC<tCEDSYNC(i+1));
                    if find(tFSYNC==tCEDSYNC(i+1))
                        sharedT = [sharedT tFSYNC(tids)];
                        if yCEDSYNC(i)
                            newyCEDSYNC = [newyCEDSYNC ones(1,length(tids))];
                        else
                            newyCEDSYNC = [newyCEDSYNC zeros(1,length(tids))];
                        end
                        newyFSYNC = [newyFSYNC yFSYNC(tids)'];
                    else
                        sharedT = [sharedT tFSYNC(tids) tCEDSYNC(i+1)];
                        if yCEDSYNC(i)
                            newyCEDSYNC = [newyCEDSYNC ones(1,length(tids)) 0];
                        else
                            newyCEDSYNC = [newyCEDSYNC zeros(1,length(tids)) 1];
                        end
                        newyFSYNC = [newyFSYNC yFSYNC(tids)' interp1(tFSYNC,yFSYNC,tCEDSYNC(i+1))];
                    end
                    handles.ChangeProgressBar(['Interpolating SYNC Data'], i, length(tCEDSYNC),0)
                end
            end
            handles.ChangeProgressBar('', 0, 0, 1);
            
            riseC = find(([0 diff(newyCEDSYNC)]>0));
            riseC = riseC(1:end-1);
            fallC = find(([0 diff(newyCEDSYNC)]<0));
            if fallC(1) < riseC(1)
                fallC(1) = [];
            end
            
            riseF = find(([0 diff(newyFSYNC)]>0.5));
            riseF = riseF(1:end-1);
            
            fallF = find(([0 diff(newyFSYNC)]<-0.5));
            if fallF(1) < riseF(1)
                fallF(1) = [];
            end

            
            if isempty(handles)
                a = figure;
                a.Position = [1 41 1920 963];
                aa = axes;
            else
                aa = axes('Parent',handles.AnalysisPlotsTab);
            end
             handles.ChangeProgressBar(['Check Paired SYNC Points, Plots will be adjusted at 150 points'], 0, 1,0)
            %aa.Title.String = 'Check Paired SYNC Points, Plots will be adjusted at 150 points';
            aa.YLim = [-0.35    1.15];
            cs = plot(aa,sharedT,newyCEDSYNC,'LineWidth',3);
            hold on
            fs = plot(aa,sharedT,newyFSYNC,'LineWidth',1.5);
            uukeep = plot(aa,0,0,'g*');
            uukeep.XData = [];
            uukeep.YData = [];
            uutemp = plot(aa,0,0,'k*');
            uutemp.XData = [];
            uutemp.YData = [];
            uutempF = plot(aa,0,0,'ko');
            uutempF.XData = [];
            uutempF.YData = [];
            hold off
            aa.YLim = [-0.35    1.15];
            
            rcV = 1;
            rfV = 1;
            keepRC = [];
            keepRF = [];
            typeDone = [];
            while (rcV < length(riseC)) && (rfV < length(riseF))
                titl = {'Check Paired SYNC Points, Plots will be adjusted at 150 points';...
                    [num2str(rcV),' of ',num2str(length(riseC)),' for CED and ',num2str(rfV),' of ',num2str(length(riseF)),' for FSYNC and ']};
                handles.ChangeProgressBar(titl, rcV, length(riseC) ,0)
%                 aa.Title.String = {'Check Paired SYNC Points, Plots will be adjusted at 150 points';...
%                     [num2str(rcV),' of ',num2str(length(riseC)),' for CED and ',num2str(rfV),' of ',num2str(length(riseF)),' for FSYNC and ']};
                uutempF.XData = [sharedT(riseF(rfV))];
                uutempF.YData = [newyFSYNC(riseF(rfV))];
%                 drawnow
                if length(keepRC) < 150
                    aa.XLim = [sharedT(riseC(rcV))-1.5    sharedT(riseC(rcV))+2];
                    uutemp.XData = [sharedT(riseC(rcV))];
                    uutemp.YData = [newyCEDSYNC(riseC(rcV))];
%                     drawnow
                elseif length(keepRC) == 150
                    tShift = sharedT(keepRF(1))-sharedT(keepRC(1));
                    tCED = sharedT+tShift;
                    x = sharedT(keepRF);
                    y = tCED(keepRC);
                    p = polyfit(x,y,1);
                    tCED = tCED./p(1);
                    cs.XData = tCED;
                    uukeep.XData(1:2:end) = tCED(keepRC);
                    aa.XLim = [tCED(riseC(rcV))-1.5    tCED(riseC(rcV))+2];
%                     drawnow
                    set(handles.figure1, 'currentaxes', aa)
                    refresh(handles.figure1)
                    uutemp.XData = [tCED(riseC(rcV))];
                    uutemp.YData = [newyCEDSYNC(riseC(rcV))];
%                     drawnow
                else
                    aa.XLim = [tCED(riseC(rcV))-1.5    tCED(riseC(rcV))+2];
                    uutemp.XData = [tCED(riseC(rcV))];
                    uutemp.YData = [newyCEDSYNC(riseC(rcV))];
%                     drawnow
                end
                handles.continueButton.BackgroundColor = 'y';
                handles.KeepCEDSYNCButton.Value = 1;
                handles.KeepFSYNCoButton.Value = 1;
                handles.continueButton.UserData = 0;
%                 set(handles.figure1, 'currentaxes', aa)
%                 refresh(handles.figure1)                
%                 fs.LineWidth = 1.5;
                handles.figure1.WindowKeyPressFcn = @(src,event)handles.KeyDownFcn;
                uiwait(handles.figure1);
                handles.KeyPressed = 0;
                handles.figure1.WindowKeyPressFcn = [];

                uutemp.XData = [];
                uutemp.YData = [];
                uutempF.XData = [];
                uutempF.YData = [];
%                 drawnow
                if handles.continueButton.UserData == -1
                    switch typeDone(end)
                        case 1
                            rcV = rcV - 1;
                            rfV = rfV - 1;
                            keepRC(keepRC==riseC(rcV)) = [];
                            keepRF(keepRF==riseF(rfV)) = [];
                            uukeep.XData(end-1:end) = [];
                            uukeep.YData(end-1:end) = [];
%                             drawnow
                        case 2
                            rcV = rcV - 1;
                            rfV = rfV - 1;
                        case 3
                            rcV = rcV - 1;
                        case 4
                            rfV = rfV - 1;
                    end
                    typeDone(end) = [];
                else
                    if handles.KeepCEDSYNCButton.Value && handles.KeepFSYNCoButton.Value
                        typeDone = [typeDone 1];
                        if length(keepRC) < 150
                            uukeep.XData = [uukeep.XData sharedT(riseC(rcV))];
                        else
                            uukeep.XData = [uukeep.XData tCED(riseC(rcV))];
                        end
                        uukeep.YData = [uukeep.YData newyCEDSYNC(riseC(rcV))];
                        keepRC = [keepRC riseC(rcV)];
%                         drawnow
                        rcV = rcV + 1;
                        uukeep.XData = [uukeep.XData sharedT(riseF(rfV))];
                        uukeep.YData = [uukeep.YData newyFSYNC(riseF(rfV))];
%                         drawnow
                        keepRF = [keepRF riseF(rfV)];
                        rfV = rfV + 1;
                    elseif  ~handles.KeepCEDSYNCButton.Value && ~handles.KeepFSYNCoButton.Value
                        typeDone = [typeDone 2];
                        rcV = rcV + 1;
                        rfV = rfV + 1;
                    elseif ~handles.KeepCEDSYNCButton.Value
                        typeDone = [typeDone 3];
                        rcV = rcV + 1;
                    elseif ~handles.KeepFSYNCoButton.Value
                        typeDone = [typeDone 4];
                        rfV = rfV + 1;
                    end
                end
                drawnow
            end
            handles.ChangeProgressBar('', 0, 0, 1);
            cla(aa)
            plot(aa,1:length(keepRC),sharedT(keepRC)-sharedT(keepRF))
            pause(0.3);
            hold(aa, 'on')
            tShift = sharedT(keepRC(1))-sharedT(keepRF(1));
            tF = sharedT+tShift;
            plot(aa,1:length(keepRC),sharedT(keepRC)-tF(keepRF))
            pause(0.3);
            x = tF(keepRF);
            y = sharedT(keepRC);
            p = polyfit(x,y,1);
            plot(aa,1:length(keepRC),sharedT(keepRC)-tF(keepRF)*p(1))
            pause(0.3);
            
            if debugFlg
                app = handles;
                delete(app.AnalysisPlotsTab.Children(isgraphics(app.AnalysisPlotsTab.Children,'axes')))
            end
            
            D.CHANS{FSYNC.Num}.Old_T = tFSYNC;
            D.CHANS{FSYNC.Num}.sharedT = sharedT;
            D.CHANS{FSYNC.Num}.sharedY = newyFSYNC;
            D.CHANS{FSYNC.Num}.RisePtLocs = keepRF; 
            D.CHANS{FSYNC.Num}.Adjusted_T = (tFSYNC+tShift)*p(1);
            D.CHANS{FSYNC.Num}.yFSYNC = yFSYNC;
            D.CHANS{FSYNC.Num}.tShift = tShift;
            D.CHANS{FSYNC.Num}.tMultiply = p(1);
            
            D.CHANS{CEDSYNC.Num}.tCEDSYNC = tCEDSYNC;
            D.CHANS{CEDSYNC.Num}.yCEDSYNC = yCEDSYNC;
            D.CHANS{CEDSYNC.Num}.sharedT = sharedT;
            D.CHANS{CEDSYNC.Num}.sharedY = newyCEDSYNC;
            D.CHANS{CEDSYNC.Num}.RisePtLocs = keepRC;
        end
        
        function StimStruct = FindMotionStim(obj, debugFlg, handles)
            if isempty(handles)
                
            else
                handles.ProgressBar.Title.String = ['Getting Motion Stim'];
                drawnow
            end
            APC = ActionPotential_class(obj.DAT);
            nineCycSinP5 = sin(linspace(0,18*pi,18000));
            nineCycSin1 = sin(linspace(0,18*pi,18000/2));
            nineCycSin2 = sin(linspace(0,18*pi,18000/4));

            sae = 4500;
            s = linspace(0,-3.05,sae);
            e = linspace(-3.05,0,sae);
            StaticTilt = [s linspace(-3.05,-3.05,29700-sae*2) e];

            StimStruct = struct();
            
            StimStruct.p5Hz.XVel.Data = [];
            StimStruct.p5Hz.XVel.start = [];
            StimStruct.p5Hz.XVel.stop = [];
            StimStruct.p5Hz.match = [];
            StimStruct.p5Hz.MotionTime = [];
            StimStruct.p5Hz.Direction = [];
            StimStruct.p5Hz.Magnitude = [];
            StimStruct.p5Hz.LARP = [];
            StimStruct.p5Hz.RALP = [];
            StimStruct.p5Hz.AP = [];
            StimStruct.p5Hz.APTime = [];
            StimStruct.p5Hz.APStart = [];
            StimStruct.p5Hz.APStop = [];
            StimStruct.p5Hz.APLocs = [];
            StimStruct.p5Hz.APSIsR = [];
            StimStruct.p5Hz.APSIsRt = [];
            StimStruct.p5Hz.StimID = [];
            
            StimStruct.p5Hz.YVel = StimStruct.p5Hz.XVel;
            StimStruct.p5Hz.ZVel = StimStruct.p5Hz.XVel;
            StimStruct.p5Hz.XAcc = StimStruct.p5Hz.XVel;
            StimStruct.p5Hz.YAcc = StimStruct.p5Hz.XVel;
            StimStruct.p5Hz.ZAcc = StimStruct.p5Hz.XVel;
            
            StimStruct.oneHz = StimStruct.p5Hz;
            StimStruct.twoHz = StimStruct.p5Hz;
            StimStruct.NoseDown = StimStruct.p5Hz;
            StimStruct.NoseUp = StimStruct.p5Hz;
            StimStruct.REarDown = StimStruct.p5Hz;
            StimStruct.LEarDown = StimStruct.p5Hz;

            D = obj.DAT;
            FSYNC = D.GetChan('FSYNC');
            CEDSYNC = D.GetChan('CED-SYNC');
            NerveCon = D.GetChan('NerveCon');
            NerveConT = 0:1/NerveCon.SampleRate:((length(NerveCon.Data)/NerveCon.SampleRate)-1/NerveCon.SampleRate);
%             chans = [{'X_VEL'} {'Y_VEL'} {'Z-VEL'} {'ACCEL_X'} {'ACCEL_Y'} {'ACCEL_Z'}];
            XVel = D.GetChan('X_VEL');
            YVel = D.GetChan('Y_VEL');
            ZVel = D.GetChan('Z-VEL');
            XAcc = D.GetChan('ACCEL_X');
            YAcc = D.GetChan('ACCEL_Y');
            ZAcc = D.GetChan('ACCEL_Z');
            fqTesting = 1;
            go = 1;
            tempSigXV = XVel.Data'-mean(XVel.Data);
            tempSigYV = YVel.Data'-mean(YVel.Data);
            tempSigZV = ZVel.Data'-mean(ZVel.Data);
            tempSigXA = (XAcc.Data')*9.8;
            tempSigXA = tempSigXA-mean(tempSigXA);
            tempSigYA = (YAcc.Data')*9.8;
            tempSigYA = tempSigYA-mean(tempSigYA);
            tempSigZA = (ZAcc.Data')*9.8;
            tempSigZA = tempSigZA-mean(tempSigZA);
            tempSigs = [tempSigXV; tempSigYV; tempSigZV; tempSigXA; tempSigYA; tempSigZA];
            if debugFlg && isempty(handles)
                ff = figure;
            elseif debugFlg && ~isempty(handles)
                ff = figure('Visible',0);
            else
                ff = [];
            end
            sig2Test = 1;
            StimID = 0;
            if isempty(handles)
                
            else
                handles.ProgressBar.Title.String = ['Getting Motion Stim: Testing X-Vel'];
                drawnow
            end
            while go
                switch fqTesting
                    case 1
                        if sig2Test > 3
                           dlim = 10000;
                        else
                            dlim = 150000;
                        end
                        [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, nineCycSinP5, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        fqTesting = 2;
                    case 2
                        if sig2Test > 3
                            dlim = 10000;
                        else
                            dlim = 150000;
                        end
                        [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, nineCycSin1, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        fqTesting = 3;
                    case 3
                        if sig2Test > 3
                            dlim = 10000;
                        else
                            dlim = 50000;
                        end
                        [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, nineCycSin2, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        fqTesting = 4;
                    case 4
                        if sig2Test == 3
                            dlim = 300000;
                        elseif sig2Test > 3
                            dlim = 10000;
                        else
                            dlim = 150000;
                        end
                        [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, -nineCycSinP5, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        fqTesting = 5;
                    case 5
                        if sig2Test == 3
                            dlim = 300000;
                        elseif sig2Test > 3
                            dlim = 10000;
                        else
                            dlim = 150000;
                        end
                        [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, -nineCycSin1, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        fqTesting = 6;
                    case 6
                        if sig2Test == 3
                            dlim = 100000;
                        elseif sig2Test > 3
                            dlim = 10000;
                        else
                            dlim = 50000;
                        end
                        [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, -nineCycSin2, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        fqTesting = 7;
                    case 7 %Nose Down
                        if sig2Test == 4
                            dlim = 10000;%5.0357e+03
                            [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, StaticTilt, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        end
                        fqTesting = 8;
                    case 8 %Nose Up
                        if sig2Test == 4
                            dlim = 10000;
                            [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, -StaticTilt, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        end
                        fqTesting = 9;
                    case 9 %Right Ear Down
                        if sig2Test ==5
                            dlim = 10000;
                            [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, StaticTilt, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        end
                        fqTesting = 10;
                    case 10 %Left Ear Down
                        if sig2Test == 5
                            dlim = 10000;
                            [StimStruct, tempSigs, StimID] = obj.FindandSaveStim(StimStruct, tempSigs, -StaticTilt, FSYNC.Adjusted_T, NerveCon.Data, NerveConT, dlim, fqTesting, sig2Test, APC, StimID, handles, debugFlg, ff);
                        end
                        sig2Test = sig2Test+1;
                        fqTesting = 1;
                        switch sig2Test
                            case 2
                                sigD = 'Y-Vel';
                            case 3
                                sigD = 'Z-Vel';
                            case 4
                                sigD = 'X-Acc';
                            case 5
                                sigD = 'Y-Acc';
                            case 6
                                sigD = 'Z-Acc';
                        end
                        if sig2Test > 6
                            go = 0;
                        end
                            if ~isempty(handles)
                                handles.PBarObj.Position(3) = (sig2Test-1)/6*1000;
                                handles.PBarTxt.String = [num2str(round((sig2Test-1)/6*100,1)),'%'];
                                handles.ProgressBar.Title.String = ['Getting Motion Stim: Testing ',sigD];
                                drawnow
                            end
                end
            end 
            if ~isempty(handles)
                handles.ProgressBar.Title.String = ['Finished'];
                handles.PBarObj.Position(3) = 0;
                handles.PBarTxt.String = 'Waiting';
                drawnow
            end
            if debugFlg && ~isempty(handles)
                app = handles;
                delete(app.AnalysisPlotsTab.Children(isgraphics(app.AnalysisPlotsTab.Children,'axes')))
                close(ff);
            end
        end
        
        function [StimStruct, tempSig, StimID] = FindandSaveStim(obj, StimStruct, tempSig, nineCyc, FSYNCTime, APSig, APTime, dlim, fqTest, sig2Test, APC, StimID, handles, debugFlg, ff)
            switch fqTest
                case {1, 4}
                    fq = 'p5Hz';
                    forList = '0.5 Hz';
                case {2, 5}
                    fq = 'oneHz';
                    forList = '1 Hz';
                case {3, 6}
                    fq = 'twoHz';
                    forList = '2 Hz';
                case 7
                    fq = 'NoseDown';
                    forList = 'Nose Down';
                case 8
                    fq = 'NoseUp';
                    forList = 'Nose Up';
                case 9
                    fq = 'REarDown';
                    forList = 'Right Ear Down';
                case 10
                    fq = 'LEarDown';
                    forList = 'Left Ear Down';
            end
            
            if sig2Test == 3
                stAmp = 12;
            elseif sig2Test > 3
                if fqTest < 7
                    stAmp = 2.5;
                else
                    stAmp = 1;
                end
            else
                stAmp = 6.5;
            end
            go2 = 1;
            while go2
                testSig = double(tempSig(sig2Test,:));
                if debugFlg && isempty(handles)
                    findsignal(testSig,nineCyc*stAmp)
                elseif debugFlg && ~isempty(handles)
                    set(0, 'CurrentFigure',ff)
                    findsignal(testSig,nineCyc*stAmp)
                    delete(handles.AnalysisPlotsTab.Children(isgraphics(handles.AnalysisPlotsTab.Children,'axes')))
                    set(findobj(ff,'type','axes'),'Parent',handles.AnalysisPlotsTab);
                end
                [istart,istop,dist]=findsignal(testSig,nineCyc*stAmp);
                if (sig2Test == 3) && (fqTest < 4)
                    [~,~,dist2]=findsignal(testSig,-nineCyc*stAmp);
                    if dist2<dist
                        dist = dlim*2;
                    end
                end
                if dist<dlim
                    StimID = StimID+1;
                    StimStruct.(fq).XVel.Data = [StimStruct.(fq).XVel.Data; tempSig(1,istart-50:istop+50)];
                    StimStruct.(fq).XVel.start = [StimStruct.(fq).XVel.start; istart-50];
                    StimStruct.(fq).XVel.stop = [StimStruct.(fq).XVel.stop; istop+50];
                    
                    StimStruct.(fq).YVel.Data = [StimStruct.(fq).YVel.Data; tempSig(2,istart-50:istop+50)];
                    StimStruct.(fq).YVel.start = [StimStruct.(fq).YVel.start; istart-50];
                    StimStruct.(fq).YVel.stop = [StimStruct.(fq).YVel.stop; istop+50];
                    
                    StimStruct.(fq).ZVel.Data = [StimStruct.(fq).ZVel.Data; tempSig(3,istart-50:istop+50)];
                    StimStruct.(fq).ZVel.start = [StimStruct.(fq).ZVel.start; istart-50];
                    StimStruct.(fq).ZVel.stop = [StimStruct.(fq).ZVel.stop; istop+50];
                    
                    StimStruct.(fq).XAcc.Data = [StimStruct.(fq).XAcc.Data; tempSig(4,istart-50:istop+50)];
                    StimStruct.(fq).XAcc.start = [StimStruct.(fq).XAcc.start; istart-50];
                    StimStruct.(fq).XAcc.stop = [StimStruct.(fq).XAcc.stop; istop+50];
                    
                    StimStruct.(fq).YAcc.Data = [StimStruct.(fq).YAcc.Data; tempSig(5,istart-50:istop+50)];
                    StimStruct.(fq).YAcc.start = [StimStruct.(fq).YAcc.start; istart-50];
                    StimStruct.(fq).YAcc.stop = [StimStruct.(fq).YAcc.stop; istop+50];
                    
                    StimStruct.(fq).ZAcc.Data = [StimStruct.(fq).ZAcc.Data; tempSig(6,istart-50:istop+50)];
                    StimStruct.(fq).ZAcc.start = [StimStruct.(fq).ZAcc.start; istart-50];
                    StimStruct.(fq).ZAcc.stop = [StimStruct.(fq).ZAcc.stop; istop+50];
                    
                    StimStruct.(fq).match = [StimStruct.(fq).match; fqTest];
                    StimStruct.(fq).MotionTime = [StimStruct.(fq).MotionTime; FSYNCTime(istart-50:istop+50)];
                    
                    rotMat = [cosd(-45) -sind(-45) 0;
                                      sind(-45) cosd(-45) 0;
                                      0 0 1];
                    vel = rotMat*[tempSig(1,istart-50:istop+50); tempSig(2,istart-50:istop+50); tempSig(3,istart-50:istop+50)];
                    
                    StimStruct.(fq).LARP = [StimStruct.(fq).LARP; vel(2,:)];
                    StimStruct.(fq).RALP = [StimStruct.(fq).RALP; vel(1,:)];
                    StimStruct.(fq).StimID = [StimStruct.(fq).StimID; StimID];
                    
                    lower = min(FSYNCTime(istart-50:istop+50));
                    upper = max(FSYNCTime(istart-50:istop+50));
                    lowI = find(APTime >= lower,1,'first');
                    upI = find(APTime <= upper,1,'last');
                    if ~isempty(StimStruct.(fq).AP)
                        sz = size(StimStruct.(fq).AP);
                        if length(APSig(lowI:upI)') > sz(2)
                            dd = length(APSig(lowI:upI)')-sz(2);
                            upI = upI - dd;
                        elseif length(APSig(lowI:upI)') < sz(2)
                            dd = sz(2)-length(APSig(lowI:upI)');
                            upI = upI + dd;
                        end
                    end
                    StimStruct.(fq).APStart = [StimStruct.(fq).APStart; lowI];
                    StimStruct.(fq).APStop = [StimStruct.(fq).APStop; upI];
                    StimStruct.(fq).AP = [StimStruct.(fq).AP; APSig(lowI:upI)'];
                    StimStruct.(fq).APTime = [StimStruct.(fq).APTime; APTime(lowI:upI)];
                    
                    [APT_unzeroed, APT] = APC.FindAPs(APSig(lowI:upI));
                    
                    apT = APTime(lowI:upI);
                    
                    StimStruct.(fq).APLocs = [StimStruct.(fq).APLocs; {APT.AP'}];
                    StimStruct.(fq).APSIsR = [StimStruct.(fq).APSIsR; {1./diff(apT(APT.AP))}];
                    StimStruct.(fq).APSIsRt = [StimStruct.(fq).APSIsRt; {apT(APT.AP(2:end))}];
                    
                    
                    [Xv,mXl]=max(smooth(tempSig(1,istart-50:istop+50),211,'sgolay'));
                    [Yv,mYl]=max(smooth(tempSig(2,istart-50:istop+50),211,'sgolay'));
                    [Zv,mZl]=max(smooth(tempSig(3,istart-50:istop+50),211,'sgolay'));
                    [XAv,mXAl]=max(smooth(tempSig(4,istart-50:istop+50),211,'sgolay'));
                    [YAv,mYAl]=max(smooth(tempSig(5,istart-50:istop+50),211,'sgolay'));
                    [ZAv,mZAl]=max(smooth(tempSig(6,istart-50:istop+50),211,'sgolay'));
                    
                    [~,maxDir] = max([Xv Yv Zv XAv YAv ZAv]);
                
                    maxLocs = [mXl mYl mZl mXAl mYAl mZAl];
                
                    Xv = tempSig(1,istart-49+maxLocs(maxDir));
                    Yv = tempSig(2,istart-49+maxLocs(maxDir));
                    Zv = tempSig(3,istart-49+maxLocs(maxDir));
                    XAv = tempSig(4,istart-49+maxLocs(maxDir));
                    YAv = tempSig(5,istart-49+maxLocs(maxDir));
                    ZAv = tempSig(6,istart-49+maxLocs(maxDir));
                    if fqTest > 6
                        maxDir = fqTest;
                    end
                switch maxDir
                    case 1
                        if Yv < 0
                            rotMat = [cosd(-45) -sind(-45) 0;
                                      sind(-45) cosd(-45) 0;
                                      0 0 1];
                            vel = rotMat*[Xv; Yv; Zv];
                            vel = abs(round(vel(2)));
                            %LARP
                            direct = 'LARP';
                        else
                            rotMat = [cosd(-45) -sind(-45) 0;
                                      sind(-45) cosd(-45) 0;
                                      0 0 1];
                            vel = rotMat*[Xv; Yv; Zv];
                            vel = abs(round(vel(1)));
                            %RALP
                            direct = 'RALP';
                        end
                        name = [direct,', ',forList,', ',num2str(vel),' m/s, ',num2str(StimID)];
                    case 2
                        if Xv < 0
                            rotMat = [cosd(-45) -sind(-45) 0;
                                      sind(-45) cosd(-45) 0;
                                      0 0 1];
                            vel = rotMat*[Xv; Yv; Zv];
                            %LARP
                            vel = abs(round(vel(2)));
                            direct = 'LARP';
                        else
                            rotMat = [cosd(-45) -sind(-45) 0;
                                      sind(-45) cosd(-45) 0;
                                      0 0 1];
                            vel = rotMat*[Xv; Yv; Zv];
                            %RALP
                            vel = abs(round(vel(1)));
                            direct = 'RALP';
                        end
                        name = [direct,', ',forList,', ',num2str(vel),' m/s, ',num2str(StimID)];
                    case 3
                        vel = abs(round(Zv));
                        direct = 'YAW';
                        name = [direct,', ',forList,', ',num2str(vel),' m/s, ',num2str(StimID)];
                        %Yaw
                    case 4
                        vel = abs(round(XAv));
                        direct = 'Surge';
                        %surge
                        name = [direct,', ',forList,', ',num2str(vel),' m/s2, ',num2str(StimID)];
                    case 5
                        vel = abs(round(YAv));
                        direct = 'Lateral';
                        %lateral
                        name = [direct,', ',forList,', ',num2str(vel),' m/s2, ',num2str(StimID)];
                    case 6
                        vel = abs(round(ZAv));
                        direct = 'Heave';
                        %heave
                        name = [direct,', ',forList,', ',num2str(vel),' m/s2, ',num2str(StimID)];
                    case 7
                        direct = forList;
                        vel = 20;
                        %Nose Down
                        name = [forList,', 20 Degrees, ',num2str(StimID)];
                    case 8
                        direct = forList;
                        vel = 20;
                        %Nose UP
                        name = [forList,', 20 Degrees, ',num2str(StimID)];
                    case 9
                        direct = forList;
                        vel = 20;
                        %Right Ear Down
                        name = [forList,', 20 Degrees, ',num2str(StimID)];
                    case 10
                        direct = forList;
                        vel = 20;
                        %Left Ear Down
                        name = [forList,', 20 Degrees, ',num2str(StimID)];
                end
                if ~isempty(handles)
                    
                    if any(ismember(handles.MotionStimList.Items, 'Listbox'))
                        handles.MotionStimList.Items = {name};
                    else
                        handles.MotionStimList.Items = [handles.MotionStimList.Items {name}];
                    end
                    drawnow
                end
                StimStruct.(fq).Direction = [StimStruct.(fq).Direction; {direct}];
                StimStruct.(fq).Magnitude = [StimStruct.(fq).Magnitude; vel];   
                    
                    tempSig(1,istart-50:istop+50) = mean(tempSig(1,:));
                    tempSig(2,istart-50:istop+50) = mean(tempSig(2,:));
                    tempSig(3,istart-50:istop+50) = mean(tempSig(3,:));
                    tempSig(4,istart-50:istop+50) = mean(tempSig(4,:));
                    tempSig(5,istart-50:istop+50) = mean(tempSig(5,:));
                    tempSig(6,istart-50:istop+50) = mean(tempSig(6,:));
                else
                    go2 = 0;
                end
            end
        end

    end
end