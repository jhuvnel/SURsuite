classdef StimTextMark < handle
   
   %StimTextMark - Class to hold array of CED TextMarks for
   %Keithley/Crosspoint stimulation.
   %
   %  Pass in the arrray of CEDTextMark objects.
   
   properties
      % The original array of TextMark objects, converted to a struct array
      % from an array of TextMark objects. This allows us to add new
      % fields.
      TM
      
      % Array, same size as TM, to hold additional information about each
      % entry in the array.
      
      % The CED Ticks per second for this data, so we can convert the
      % TextMark "times" (which are in Ticks) to seconds.
      TicksPerSec
      
      % Array of stim block Start times, in seconds. This is separated out
      % as a simple array of double numbers, so we can do a fast Binary
      % Search.
      StartTimes_sec
   end
   
   methods(Static)      
      % Helper function (should make more generally available) to extract a
      % particular channel from a CED data struct returned from
      % spikeload(). Pertinent struct fields are:
      %
      %   Chans: 1xN cell array of channel data arrays
      %   ChanNames: 1xN cell array of channel name strings
      %
      %   (where N is the number of channels in the CED data struct)      
      function CHAN = CEDGetChan(CED_DAT, ChanName)
         CHAN = [];
         chnum = find(strcmp(CED_DAT.ChanNames, ChanName));
         if isempty(chnum); return; end
         CHAN = CED_DAT.Chans{chnum};
      end
   end

   
   methods
      
      % Constructor.
      function obj=StimTextMark(TM_OR_CEDDAT, TicksPerSec, Ticks_fixed)

         % If the first argument is a data struct returned by spikeload(),
         % then extract the TextMark channel, and TicksPerSec.
         if isfield(TM_OR_CEDDAT, 'TicksPerSec')
            TM = StimTextMark.CEDGetChan(TM_OR_CEDDAT, 'TextMark');
            if isempty(TM)
               warning('No "TextMark" channel in CED data!!');
               return;
            end
            TicksPerSec = TM_OR_CEDDAT.TicksPerSec;
         elseif isa(TM_OR_CEDDAT, 'CEDTextMark')
            % We also accept an array of CEDTextMark objects.
            TM = TM_OR_CEDDAT;
         else
            error('Must pass CED data, or CEDTextMark array!!')
         end
         
         % First, set up TM as a struct (this takes just the first item in
         % the array of TextMark objects).
         obj.TM = struct(TM);  % Set up first object as a struct.
         
         % Now fill in all of the items.
         for i=2:length(TM)
            obj.TM(i) = struct(TM(i));
         end

         % If Ticks_fixed is passed in, then overwrite the ticks.
         if nargin >= 3
            % First, just for good measure, save the original ticks.
            [obj.TM.m_Time_UNFIXED] = obj.TM.m_Time;
            
            % "Trick" to distribute an array out to fields of struct array.
            ticks_cells = num2cell(Ticks_fixed);
            [obj.TM.m_Time] = ticks_cells{:};
         end
         
         obj.TicksPerSec = TicksPerSec;
         obj.ParseStrings();
      end
         
      % Parse each of the TextMark strings to extract information about the
      % stimulation parameters used. The text mark string looks like this
      % (must handle BOTH styles, eventually):
      %
      % Old style:
      %    'Stm1 Ref11 2.0uA 1x22x5mS 100/25/100uSec'
      %
      % New style:
      %    'Stm1 Ref10 20.0uA 1x(22x5mS+50mS_Gap)  100/25/100uSec [14:16:04.931]'
      %
      function ParseStrings(obj)
         for i=1:length(obj.TM)            
            parts = strsplit(obj.TM(i).m_Data);
            
            % Make sure we have a valid stim TextMark
            if (length(parts) < 5) || ~strncmp(parts{1}, 'Stm', 3)
               fprintf('Skipping TextMark at %.1fsec: "%s"\n', ...
			        double(obj.TM(i).m_Time)/obj.TicksPerSec, obj.TM(i).m_Data);
               continue;
            end
            
            Stim_E = str2double(parts{1}(4:end));
            Ref_E = str2double(parts{2}(4:end));
            Current_uA = str2double(parts{3}(1:end-2));
            

            % This should work for both old and new styles.
            P = strsplit(parts{4}, {'x','m','(', '+', 'S'});
            NGroups = str2double(P{1});
            NPulsePerGroup = str2double(P{2});
            PulseInterval_msec = str2double(P{3});
            
            % Old style does not have the GAP duration, so assume 50mSec.
            % Note that the old style can return a 0x0 empty char array as
            % the fourth part, so we check for that.
            Gap_msec = 50;
            if (length(P) >= 4) && ~isempty(P{4})
               Gap_msec = str2double(P{4});
            end

            % Total stim duration, in seconds, minus last gap time.
            Duration_sec = (NGroups*(NPulsePerGroup*PulseInterval_msec+Gap_msec) - Gap_msec)/1000.0;
            %fprintf('%3d: %d/%d %4.0fuA  Dur:%.3f\n', i, Stim_E, Ref_E, Current_uA, Duration_sec);
            
            obj.TM(i).Start_sec = double(obj.TM(i).m_Time) / obj.TicksPerSec;
            obj.TM(i).Duration_sec = Duration_sec;
            obj.TM(i).Stim_E = Stim_E;
            obj.TM(i).Ref_E = Ref_E;
            obj.TM(i).Current_uA = Current_uA;
         end
         
         % Separately save just the array of start times, so we can do a
         % very fast Binary Search in StimNearTime(). I use the transpose
         % operator at the end, just because I like column vectors, rather
         % than row vectors!
         obj.StartTimes_sec = [obj.TM.Start_sec]';
      end
      
      % Search the TM array to see if Time_sec is within a stimulation
      % block. Return empty if not.
      function [TM_near, idx] = StimNearTime(obj, Time_sec)
         TM_near = [];
         % BinarySearch_LE() will return the index of the time that is Less
         % Than or Equal To the passed Time_sec.
         idx = BinarySearch_LE(obj.StartTimes_sec, Time_sec);
         if isempty(idx); return; end
         TM = obj.TM(idx);
         
         % Make sure the time is less than the end time of the stim block.
         if Time_sec > (TM.Start_sec + TM.Duration_sec); return; end;
         TM_near = TM;
      end
      
      % This one returns a concise string summarizing the stimulation
      % block, instead of returning the raw TM structure.
      function [Stim_str, Block_idx] = StimNearTime_str(obj, Time_sec)
         Stim_str = '';
         [TM, Block_idx] = obj.StimNearTime(Time_sec);
         if isempty(TM); return; end
         Stim_str = sprintf('Stm_%d:Ref_%d %.0fuA', TM.Stim_E, TM.Ref_E, TM.Current_uA);
      end
      
   end  %methods

end  %class
