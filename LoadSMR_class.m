% This class handles loading of CED Spike2 "*.SMR" data files. It is
% "generic" in the sense that it is not specifically focused on particular
% channels that we use for stimulation. It is up to the caller to decide
% which channels to actually load from the file, using LoadChan().
%
classdef LoadSMR_class < handle
   %UNTITLED2 Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      F % File handle.
      Fname
      TickInterval_sec
      StartTime % Matlab datetime object with data start "wall clock" time.
      MaxTick
      AllChanList
   end
   
   methods
      function delete(obj)
         %fprintf('Closing SMR file %s!!\n', obj.Fname);
         CEDS64Close(obj.F);
         obj.F = 0;
      end
      
      % Open the SMR file, and read in all of the channel names and types,
      % and the basic file information. The constructor does not otherwise
      % load any of the actual channel data from the file.
      function obj=LoadSMR_class(fname, handles)
         if ~libisloaded('ceds64int')
            if isempty(which('CEDS64LoadLib'))
                % Find where THIS file (LoadSMR_class) is, and look there
                % for the CEDS64ML sub-directory.
                ThisPath = fileparts(which('LoadSMR_class'));
                addpath([ThisPath '/CEDS64ML']);
            end
            cedpath = fileparts(which('CEDS64LoadLib'));
            % Add the 'matson' path so CEDS64LoadLib.m is found
            %addpath( cedpath );
            % Load the library into the Matlab environment
            CEDS64LoadLib( cedpath ); % load ceds64int.dll
         end
         
         % Open file and get basic info, such as Tick Interval, etc.
         obj.Fname = fname;
         F = CEDS64Open( fname, 1 );  % Open Read-Only.
         obj.F = F;
         if (F <= 0); error(sprintf('Could not open file %s', fname)); end
         obj.TickInterval_sec = CEDS64TimeBase(F);
         obj.MaxTick = CEDS64MaxTime(F);
         [ok,dt]= CEDS64TimeDate(F);
         % 7 vals are: hundredths, secs, mins, hours, day, month, year.
         obj.StartTime = datetime(dt(7),dt(6),dt(5),dt(4),dt(3),dt(2),dt(1)*10);
         MAXchan = CEDS64MaxChan(F);
         
         obj.AllChanList=[];
         TypeStrings = {'Analog', 'EventFall', 'EventRise', 'EventBoth', ...
            'Marker', 'WaveMark', 'RealMark', 'TextMark', 'RealWave'};
         % Gather all channel information.
         for ch=1:MAXchan
            Type = CEDS64ChanType(F, ch);
            if Type == 0; continue; end
            CHAN.Num = ch;
            [ok, Title] = CEDS64ChanTitle(F, ch);
            CHAN.Name = Title;
            CHAN.Type = Type;
            CHAN.Type_str = TypeStrings{CHAN.Type};
            
            if isempty(obj.AllChanList); obj.AllChanList=CHAN;
            else obj.AllChanList(end+1) = CHAN; end
         end
      end % function LoadSMR()
      
      % Load several channels, by name, but don't throw an error if the
      % channel is not in the file. Pass in a Cell Array of channel names.
      % We return a Cell Array of Structures, one struct per channel. For
      % convenience, we also return a BOOL array of which channels were
      % present in the file (caller can ignore if not useful).
      function [CHANS,ChansPresent] = LoadChans(obj, ChanNames, handles)
         ChansPresent = ismember(ChanNames, {obj.AllChanList.Name});
         CHANS = {};
         if isempty(handles)
         ff = waitbar(0,'Please wait...');
         end
         for idx=1:length(ChanNames)
            Name = ChanNames(idx);
            
            if isempty(handles)
                waitbar((idx-1)/length(ChanNames),ff,['Loading Channel: ',Name{1}])
            else
                handles.ProgressBar.Title.String = ['Loading Channel: ',Name{1}];
                if any(ismember(handles.ChannelList.Items, 'Listbox'))
                    handles.ChannelList.Items = Name;
                else
                    handles.ChannelList.Items = [handles.ChannelList.Items Name];
                end
                drawnow
            end

            ChanIdx = find(strcmp(Name, {obj.AllChanList.Name}));
            if isempty(ChanIdx); continue; end
            CHAN = obj.LoadChan(Name);
            if ismember(CHAN.Name,{'NerveCon', 'ZapOut'})
                preD = CHAN.Data(1:5000000);
                preD = preD - mean(preD);
                SMOO = smooth(preD,21,'sgolay');
                [pospks,lposocs,posw,posp] = findpeaks(SMOO,1:length(SMOO),'minpeakheight',.05,'minpeakdistance',900);
                [negpks,neglocs,negw,negp] = findpeaks(-SMOO,1:length(SMOO),'minpeakheight',.05,'minpeakdistance',900);
                if mean(pospks) < mean(negpks)
                    CHAN.Data = -CHAN.Data;
                end

            end
            if isempty(CHAN); continue; end
            CHANS{end+1} = CHAN;
            
            if isempty(handles)
                waitbar((idx)/length(ChanNames),ff,['Loaded Channel: ',Name{1}])
            else
                handles.ProgressBar.Title.String = ['Loaded Channel: ',Name{1}];
                handles.PBarObj.Position(3) = (idx)/length(ChanNames)*1000;
                handles.PBarTxt.String = [num2str(round((idx)/length(ChanNames)*100)),'%'];
                drawnow
            end
         end
         if isempty(handles)
         waitbar((idx)/length(ChanNames),ff,['Finished'])
         close(ff)
         else
             handles.ProgressBar.Title.String = ['Finished'];
                handles.PBarObj.Position(3) = 0;
                handles.PBarTxt.String = 'Waiting';
                drawnow
         end
      end      
      
      % Load one channel, by name.
      function CHAN = LoadChan(obj, ChanName)
         if ~CEDS64IsOpen(obj.F); error('SMR file is closed!!'); end
         chidx = find(strcmp({obj.AllChanList.Name}, ChanName));
         if isempty(chidx); error('No such channel "%s"!!', ChanName); end
         
         CHAN = obj.AllChanList(chidx);
         
         if CHAN.Type == 1       % Analog ADC channel
            HandleAnalogChan();
            
         elseif CHAN.Type == 3   % EventRise Channel
            % Event channel, just returns the times of each event.
            MAX_EVENTS = 500000;
            [N,CHAN.Data]=CEDS64ReadEvents(obj.F,CHAN.Num,MAX_EVENTS,0);
            CHAN.Data = double(CHAN.Data); % We DO NOT want int64's.
            if N>=MAX_EVENTS
               warning(['Channel "%s" returned the MAX #Events of %d.\n' ...
                  '    (Consider increasing MAX_EVENTS in spikeload())\n'], ...
                  CHAN.Name, MAX_EVENTS);
            end
         elseif CHAN.Type == 2   % EventFall Channel
            % Event channel, just returns the times of each event.
            MAX_EVENTS = 500000;
            [N,CHAN.Data]=CEDS64ReadEvents(obj.F,CHAN.Num,MAX_EVENTS,0);
            CHAN.Data = double(CHAN.Data); % We DO NOT want int64's.
            if N>=MAX_EVENTS
               warning(['Channel "%s" returned the MAX #Events of %d.\n' ...
                  '    (Consider increasing MAX_EVENTS in spikeload())\n'], ...
                  CHAN.Name, MAX_EVENTS);
            end
         elseif CHAN.Type == 4   % EventRise/Fall Channel
            % Event channel, just returns the times of each event.
            MAX_EVENTS = 500000;
			% Returns times of ALL rising and falling edges, and
			% FirstTimestampLevel is either 0 or 1, giving the level after
			% the first timestamp.
            [N,CHAN.Data,CHAN.FirstTimestampLevel]=CEDS64ReadLevels(obj.F,CHAN.Num,MAX_EVENTS,0);
            CHAN.Data = double(CHAN.Data); % We DO NOT want int64's.
            if N>=MAX_EVENTS
               warning(['Channel "%s" returned the MAX #Events of %d.\n' ...
                  '    (Consider increasing MAX_EVENTS in spikeload())\n'], ...
                  CHAN.Name, MAX_EVENTS);
            end

         elseif CHAN.Type == 8   % TextMark channel
            % Text mark. For some reason, this is VERY slow.
            MAX_MARKS = 10000;
            [N,CHAN.Data] = CEDS64ReadExtMarks(obj.F,CHAN.Num,MAX_MARKS,0);
            % ??? FIX THIS!! Causes error. Dale 12/5/2019
            %CHAN.Data = double(CHAN.Data); % We DO NOT want int64's.
            if N>=MAX_MARKS
               warning(['WARNING: Channel "%s" returned the MAX #Marks of %d.\n' ...
                  '    (Consider increasing MAX_MARKS in spikeload())\n'], ...
                  CHAN.Name, MAX_MARKS);
            end
         else
            error('Invalid SMR channel type %d', CHAN.Type);
         end
         
         % The Analog Channel logic is confusing enough that we break it
         % out into a separate, nested, sub-function. AND THEN within that,
         % we have two sub-functions to handle the data SEGments.
         function HandleAnalogChan()
            CHAN.TicksPerSample = CEDS64ChanDiv(obj.F, CHAN.Num);
            CHAN.SampleRate = round(1/(CHAN.TicksPerSample*obj.TickInterval_sec));
            CHAN.TicksPerSec = CHAN.TicksPerSample*CHAN.SampleRate;
            CHAN.SampleInterval_sec = 1/CHAN.SampleRate;
            CHAN.TickInterval_sec = CHAN.SampleInterval_sec/CHAN.TicksPerSample;
            % Init, then do a "while" loop to load all of the "segments" of
            % the channel. Keep track of any "gaps/segments" in the data,
            % but otherwise just append all the data into a single Values[]
            % array.
            tick_start_read = 0; % Start first Read() from Ticks==0.
            segnum = 0;
            clear SEGS;
            SEGS=[];
            Values = [];

            % Read the data in small clumps. If this is too large, the
            % Read slows down.
            MAX_READ = 5e6;
            %MAX_READ = 548529; % Exact divisor for our test file's first segment.
            %MAX_READ = 548528; % One less.
            
            % Need to declare all variables shared with the nested functions.
            SEG=[]; SegVals=[];  % Initialize the SEG to empty.
            NRead=1; % Init to 1 to get into the "while" loop.
            while NRead > 0
               % Read each "clump" of MAX_READ samples.
               [ NRead, fVals, tick_start ] = CEDS64ReadWaveF(obj.F, CHAN.Num, MAX_READ, tick_start_read);
               UpdateSEG();
               tick_start_read = tick_start + CHAN.TicksPerSample*NRead;
            end
            CloseSEG(); % Be sure to close out last SEG, if it is not empty.

            if isempty(Values); CHAN=[]; return; end
            
            % Done reading! Save all of the analog data, and SEGment info
            % to the CHAN, to return to the caller.
            CHAN.Data = Values;
            CHAN.Segments = SEGS;

            % ============================================================
            % Sub-functions to handle SEGments (separated areas of data
            % that occur when "PAUSE" was used during data collection).
            function UpdateSEG()
               % Called after each Read() of a MAX_READ "clump" of samples,
               % to save the data, and check for gaps in the data.
               
               % If we get past this, we know we HAVE SOME DATA to deal with.
               if NRead <= 0; return; end % NRead of zero is End Of Channel.
               
               %fprintf('GAP %d, #Read %d, tick_start %d\n', tick_start-ticks_cur, NRead, tick_start);
               % Is this the first read of this segment? Save tick_start.
               if isempty(SEG)
                  SEG.tick_start = tick_start;
                  
               elseif (tick_start_read>0) && (tick_start>tick_start_read)
                  % If this was not the first Read of the segment, then was
                  % there a "gap" in the data? If so, the data we just read
                  % belongs in the NEXT segment!!
                  fprintf('GAP!! ticks_cur %d, tick_start %d\n', tick_start_read, tick_start);
                  CloseSEG();  % Close previous SEG.
                  SEG.tick_start = tick_start; %  Save tick_start to this new SEG.
               end
               
               % Now we simply save any data we have.
               SegVals = [SegVals; fVals];  % Append the data.
            end
            
            function CloseSEG()
               % Called when a break in the data is detected, or when we
               % reach the end of the file.
               if isempty(SEG); return; end  % Nothing to do!
               
               % Save the SEGment information to the SEGS array.
               SEG.NSamp = length(SegVals);
               segnum = segnum+1;
               if isempty(SEGS)
                   SEGS=SEG;
               else
                   SEGS(segnum) = SEG;
               end
               %fprintf('      SEGment #%d of %d samples is CLOSED!!\n', segnum, SEG.NSamp);
               
               % Between two segments? Put a marker pulse in the analog data.
               if segnum > 1
                  Values(end) = 6.0;
                  SegVals(1) = -6.0;
               end
               Values = [Values; SegVals]; % Append segment data to the overall array of data samples.
               SEG=[]; SegVals=[];  % Initialize the next SEG to empty.
            end
            
         end % Function HandleAnalogChan()
      end % LoadChan()
   end % Methods
   
   
   methods(Static)      
      % "Fix" event "tick" times to match up with an analog data channel.
      % Subtract out the data start times, and account for any gaps in the
      % data. We make use of the SEGment information that was saved by
      % HandleAnalogChan() above. Of course, the Ticks MUST be from the
      % same file as the analog CHAN.
      %
      % CHAN - The CHAN struct returned from LoadChan().
      % Ticks_in - Tick event times from any event channel in the same file.
      %
      % The returned Ticks_fixed will have their Tick times re-aligned to
      % match the saved analog data in CHAN. Any Ticks that occurred
      % outside of the Analog data range are deleted.
      %
      % NOTE that since Matlab array indexing STARTS AT ONE instead of
      % zero, and we want the Analog sample index to map directly to Ticks,
      % we ADD ONE SAMPLE TIME to all of the events in Ticks_fixed, so that
      % the lowest Tick number possible will be one sample time (i.e.,
      % TickDiv).
      %
      function Ticks_fixed = FixTicks(CHAN, Ticks_in)
         Ticks_fixed = [];
         if isempty(Ticks_in); return; end
         LastSEGEndTick = 0;
         
         % Make sure they are int64 type.
         Ticks_in = int64(round(Ticks_in));
         
         % Loop over the analog data segments, subtract out the start tick
         % times of each segment, and eliminate any Ticks that are outside
         % the Analog data range.
         for seg_i=1:length(CHAN.Segments)
            SEG = CHAN.Segments(seg_i);
            % Number of "ticks" for this number of samples.
            NSegTicks = SEG.NSamp * CHAN.TicksPerSample;
            % Find end "tick" of this segment.
            endtick = SEG.tick_start + NSegTicks;
            % Find all Tick_in events < this segment's end tick value, and
            % then subtract out the tick_start time for this Analog data
            % segment.
            Ticks_seg = Ticks_in(Ticks_in < endtick) - SEG.tick_start;
            Ticks_in(Ticks_in < endtick) = []; % Remove from input Ticks array.
            Ticks_seg(Ticks_seg<0) = []; % Remove negative Tick values.
            Ticks_fixed = [Ticks_fixed; (Ticks_seg+LastSEGEndTick)];
            LastSEGEndTick = LastSEGEndTick + NSegTicks;
         end
         
         % Add in one Analog sample time, to account for Matlab array
         % indexing starting from 1. That way, analog sample indexes will
         % directly map to ticks.
         Ticks_fixed = Ticks_fixed + CHAN.TicksPerSample;
         
      end % FixTicks()
   end % methods(Static)
end
