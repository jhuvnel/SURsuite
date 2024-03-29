classdef StimSet_class < handle
   % This class manages an "array" of itself. There is one StimSet_class
   % object per Data File.
   
   properties
      dir
      dirfull
      dirdate
      fname
      TicksPerSec
      STM  % StimTextMark objects, if the file has TextMarks. Otherwise, empty.
      
      % Stim_ticks from the SpikeData_class. Contains simply tick times for
      % every stimulation pulse event from the original Spike SMR data
      % file.
      Stim_ticks % Array of individual stim pulses in this file.
      StimBlocks % Array of all "blocks" of stimulation in this file.
      StimSets   % Array of all "sets" of stimulation in this file.
   end
   
   properties (Dependent)
      matnamefull  % Combines dirfull and fname.
   end
   
   % Does not completely prevent loss of static/constant/persistent array
   % when Class M-file is updated. Darn.
   %     properties(Constant)
   %        STATIC = StaticHelper_class
   %     end
   
   methods
      function value=get.matnamefull(obj)
         value = [obj.dirfull '\' obj.fname];
      end
      
      % Load the MAT data file associated with this object.
      function DATA = LoadMAT(obj)
         DATA = SpikeDATA_class.LoadMAT([obj.dirfull '\' obj.fname]);
      end
      
      % Return the BLOCKS array from the indexed StimSet, optionally with
      % Zero_tick subtracted out.
      function BLOCKS = GetBlocksZeroed(obj, SSidx, Zero_tick)
         BLOCKS = obj.StimSets(SSidx).BLOCKS;
         if nargin<3; return; end
         [BLOCKS.Start_tick] = dealvec([BLOCKS.Start_tick] - Zero_tick);
         [BLOCKS.End_tick] = dealvec([BLOCKS.End_tick] - Zero_tick);
         [BLOCKS.Start_sec] = dealvec([BLOCKS.Start_tick] ./ obj.TicksPerSec);
      end
      
      % Find Stim Sets with "glitches", where there are an unexpected
      % number of stimulations in a stim block. We can also take an array
      % of objects.
      function FindGlitches(obj)
         for OBJidx=1:length(obj)
            O = obj(OBJidx);
            for SSidx=1:length(O.StimSets);
               SS = O.StimSets(SSidx);
               mode_NStim = mode([SS.BLOCKS.NStim]); % Find most common value.
               NGLITCH = sum([SS.BLOCKS.NStim] ~= mode_NStim);
               if NGLITCH > 0
                  fprintf('%s %s  StimSet %d/%d has %d blocks with stim event glitches\n', ...
                     O.dirdate, strtrim(O.fname(1:min(6,length(O.fname)))), ...
                     SSidx, length(O.StimSets), NGLITCH);
                     
               end
            end
         end
      end
      
      
      % This is a simpler "Match" function which just works on all the
      % blocks in a single SS StimSet object. This one is "unaware" of any
      % StimSet "groups" of blocks.
      %
      % The info from the matching StimSet is stored directly into the
      % StimBlock as additional columns.
      function MatchStimBlocks(obj, varargin)
          SS = obj;
          MinBlockDuration_sec = 0.05; % Duration in seconds, else block is skipped.
          
          if ismember({class(varargin{1})},'SURsuite_App')
              handles = varargin{1};
          else
              handles = [];
              for k=1:2:length(varargin)
                  NAME = varargin{k};
                  if k+1 > length(varargin); fprintf('Value not provided for argument "%s"\n', NAME); break; end
                  VAL = varargin{k+1};
                  switch(lower(varargin{k}))
                      case {'MinBlockDuration_sec'}
                          MinBlockDuration_sec = double(VAL);
                      otherwise
                          fprintf('BAD named parameter: %s\n', varargin{k});
                  end
              end
          end

         
         
         
         % Make sure the new fields are added to the StimBlocks array.
         SS.StimBlocks(1).Stim_E=[];
         SS.StimBlocks(1).Ref_E=[];
         SS.StimBlocks(1).Current_uA=[];
         SS.StimBlocks(1).ToUse=[];

         % Loop through all blocks in the SS.StimBlocks array.
         if isempty(handles)
            ff = waitbar(0,'Parsing Stim Blocks...');
         else
                handles.ProgressBar.Title.String = ['Parsing Stim Blocks'];
                drawnow
         end
         combs = [];

         for k=1:length(SS.StimBlocks)
            NEW_BLOCK = FindTextMark(SS.StimBlocks(k));
            if isempty(NEW_BLOCK); continue; end;
            SS.StimBlocks(k) = NEW_BLOCK;

            if SS.StimBlocks(k).NStim > 100
                if ~isfield(SS.StimBlocks(k),'ToUse') || isempty(SS.StimBlocks(k).ToUse)
                    s = SS.StimBlocks(k).Stim_E;
                    r = SS.StimBlocks(k).Ref_E;
                    c = SS.StimBlocks(k).Current_uA;
                    if isempty(combs)
                        combs = [combs; [s r c]];
                        SS.StimBlocks(k).ToUse = 1;
                    else
                        if any(ismember(combs,[s r c],'row'))
                            SEL = logical(1:k);

                            % Now filter by each given parameter.
                            SEL = SEL & ismember([SS.StimBlocks(1:k).NStim], 198);
                            SEL = SEL & ismember([SS.StimBlocks(1:k).Stim_E], s);
                            SEL = SEL & ismember([SS.StimBlocks(1:k).Ref_E], r);
                            SEL = SEL & ismember([SS.StimBlocks(1:k).Current_uA], c);
                            idxsT = find(SEL);
                            if length(idxsT)>1
                                if length([SS.StimBlocks(idxsT).ToUse]) < length(idxsT)
                                    txt = ['Duplicates for Stim ', num2str(s),', Ref ', num2str(r), ', At ',num2str(c),'uA were found. Please choose one to use (Stim Block #)'];
                                    answer = nbuttondlg(txt,cellfun(@num2str,num2cell(idxsT),'UniformOutput',false),'PromptTextHeight',50);
                                    SS.StimBlocks(str2num(answer)).ToUse = 1;
                                    SS.StimBlocks(idxsT(idxsT~=str2num(answer))).ToUse = 0;
                                end
                            end
                        else
                            combs = [combs; [s r c]];
                            SS.StimBlocks(k).ToUse = 1;
                        end
                    end
                end
            else
                if isempty(SS.StimBlocks(k).ToUse)
                    SS.StimBlocks(k).ToUse = 1;
                end
            end



            if isempty(handles)
               waitbar(k/length(SS.StimBlocks),ff,'Parsing Stim Blocks...');
            else
                handles.PBarObj.Position(3) = (k)/length(SS.StimBlocks)*1000;
                handles.PBarTxt.String = [num2str(round((k)/length(SS.StimBlocks)*100)),'%'];
                if SS.StimBlocks(k).NStim > 100
                    if ismember(handles.StimList.Items,'Listbox')
                        handles.StimList.Items = {num2str(SS.StimBlocks(k).Stim_E)};
                    else
                        if ~ismember(handles.StimList.Items,num2str(SS.StimBlocks(k).Stim_E))
                            handles.StimList.Items = [handles.StimList.Items {num2str(SS.StimBlocks(k).Stim_E)}];
                        end
                    end
                    
                    if ismember(handles.RefList.Items,'Listbox')
                        handles.RefList.Items = {num2str(SS.StimBlocks(k).Ref_E)};
                    else
                        if ~ismember(handles.RefList.Items,num2str(SS.StimBlocks(k).Ref_E))
                            handles.RefList.Items = [handles.RefList.Items {num2str(SS.StimBlocks(k).Ref_E)}];
                        end
                    end
                    
                    if ismember(handles.CurrentList.Items,'Listbox')
                        handles.CurrentList.Items = {num2str(SS.StimBlocks(k).Current_uA)};
                    else
                        if ~ismember(handles.CurrentList.Items,num2str(SS.StimBlocks(k).Current_uA))
                            handles.CurrentList.Items = [handles.CurrentList.Items {num2str(SS.StimBlocks(k).Current_uA)}];
                        end
                    end
                end
                drawnow
            end
            
         end
         if isempty(handles)
             waitbar(k/length(SS.StimBlocks),ff,'Finished');
             close(ff)
         else
             handles.ProgressBar.Title.String = ['Finished'];
             handles.PBarObj.Position(3) = 0;
             handles.PBarTxt.String = 'Waiting';
             drawnow
         end
         

         % Sub-function to fine the Text Mark that matches this Stim Block.
         function BLOCK = FindTextMark(BLOCK)
            % Do it the "slow", non-binary-search way. Make sure start
            % times are w/in 0.3 sec. Make sure that the Text Mark is
            % BEFORE the end of the block. And that the block duration is
            % above the requested minimum.
            tm_idx = find(  (BLOCK.Start_sec-[SS.STM.TM.Start_sec] < 0.3) ...
                          & (BLOCK.Start_sec-[SS.STM.TM.Start_sec] > -0.1) ...
                          & ([SS.STM.TM.Start_sec] < BLOCK.End_sec) ...
                          & ([SS.STM.TM.Duration_sec] > MinBlockDuration_sec));
            if isempty(tm_idx)
               %BLOCK=[];
               BLOCK.Stim_E=nan; BLOCK.Ref_E=nan; BLOCK.Current_uA=nan;
            else
               % Take the "end" one. Hopefully there is only ONE match, but
               % if there is more than one, then the "end" TextMark is
               % closest to the start of the Stim Block.
               TM = SS.STM.TM(tm_idx(end));
               BLOCK.Stim_E = TM.Stim_E;
               BLOCK.Ref_E = TM.Ref_E;
               BLOCK.Current_uA = TM.Current_uA;
            end
         end
      end
      
      
      
   end % methods section
  
   
   methods(Static)
      
      % See if the given file is already present in the StimSet array.
      function found = ExistsInArray(ARY, dirdate, fname)
         found = false;
         for i=1:length(ARY)
            if strcmp(ARY(i).dirdate, dirdate) && strcmp(ARY(i).fname, fname)
               found=true; return;
            end
         end
      end
      
      % Save or get the persistent copy of the full array of StimSet_class
      % objects. Use a separate class for the persistence, so as we edit
      % this class file, we don't keep losing the persistent value. For
      % convenience, returns StimSet_class.empty if the array is empty or
      % has not yet been initialized.
      function ARY = GetSetARY(ARY)
         if nargin>0;
            PersistentClass.GetSet(ARY);
         else
            ARY = PersistentClass.GetSet;
         end
         if isempty(ARY); ARY=StimSet_class.empty; end
      end
      
      % Convenience function to convert this list of ARY objects, into a
      % struct array. This makes it easier to browse through using Matlab's
      % workspace variable viewer. Otherwise, "objects" are just displayed
      % as a single cell with the object type - you cannot see all the
      % "properties", whereas with a struct, it shows the fields as columns
      % in a table.
      function SS_struct = SS(ARY)
         if nargin < 1; ARY=StimSet_class.GetSetARY; end
         warning('off', 'MATLAB:structOnObject'); % hide annoying warning.
         SS_struct = arrayfun(@(x) struct(x), ARY);
         % Add extra column showing NumStimBlocks for each stim set.
         for i=1:length(SS_struct);
            SSS = SS_struct(i).StimSets;
            if isempty(SSS); continue; end;
            SS_struct(i).SetCounts = [SSS.NumStimBlocks];
         end
      end
      
      % The persistent variable keeps getting wiped out as we edit this
      % class file. Pain, during development.
      %       function ARY_out = GetSetARY(ARY_in)
      %          persistent ARY_persist;
      %          if nargin>0; ARY_persist=ARY_in; end
      %          if isempty(ARY_persist); ARY_persist=StimSet_class.empty; end
      %          ARY_out = ARY_persist;
      %       end

      % Match up TextMarks with StimSet blocks.
      function MatchStimBlockInfo()
         % Loop over all of the files.
         ARY=StimSet_class.GetSetARY;
         for i=1:length(ARY)
            FILE = ARY(i);
            fprintf('%d;', i);
            if ~mod(i,20); fprintf('\n'); end
            if isempty(FILE.STM);
               NoTextMarks(FILE);
               continue;
            end
            
            % This part applies to files with Text Marks.
            % Loop over StimSets in each file
            for j=1:length(FILE.StimSets)
               SS = FILE.StimSets(j);
               % Add new fields to hold electrode and stim info.
               FILE.StimSets(j).BLOCKS(1).Stim_E=[];
               FILE.StimSets(j).BLOCKS(1).Ref_E=[];
               FILE.StimSets(j).BLOCKS(1).Current_uA=[];
               % Loop over each block in the set.
               for k=1:length(SS.BLOCKS)
                  NEW_BLOCK = FindTextMark(SS.BLOCKS(k));
                  % Make sure to write back to the handle object.
                  % If we just update BLOCK, it will not get saved.
                  if isempty(NEW_BLOCK); continue; end;
                  FILE.StimSets(j).BLOCKS(k) = NEW_BLOCK;
               end
            end
         end
         
         function BLOCK = FindTextMark(BLOCK)
            % Do it the "slow", non-binary-search way.
            tm_idx = find( (abs([FILE.STM.TM.Start_sec]-BLOCK.Start_sec) < 0.3) ...
                          & ([FILE.STM.TM.Duration_sec] > 0.5));
            if isempty(tm_idx)
               BLOCK=[];
            else
               TM = FILE.STM.TM(tm_idx(1));
               BLOCK.Stim_E = TM.Stim_E;
               BLOCK.Ref_E = TM.Ref_E;
               BLOCK.Current_uA = TM.Current_uA;
            end
         end
         
         % This is called if the file has no TextMarks.
         % For now, we just recognize one particular stim set, if it is
         % exactly 108 blocks long. It is a unipolar set, where the 10 and
         % 11 reference alternates every other block. Order is:
         % StimE 1,4,7,3,6,9,2,5,8. Currents are 2,20,50,100,150,200.
         % 9*2*6==108.
         %
         function NoTextMarks(FILE)
            %Ref_array = [10 11];
            Ref_array = [11 10]; % Ref 11 is FIRST!! Not 10!!
            Cur_array = [2 20 50 100 150 200];
            Stim_array = [1 4 7 3 6 9 2 5 8];

            for j=1:length(FILE.StimSets)
               SS = FILE.StimSets(j);
               if SS.NumStimBlocks ~= 108; continue; end
               %fprintf('\n  === Found 108 block  %s ===\n', FILE.dirdate);
               % Init Ref, Current, and Stim indices.
               Ridx=1; Cidx=1; Sidx=1;
               
               % Add new fields to hold electrode and stim info.
               FILE.StimSets(j).BLOCKS(1).Stim_E=[];
               FILE.StimSets(j).BLOCKS(1).Ref_E=[];
               FILE.StimSets(j).BLOCKS(1).Current_uA=[];
               
               % Loop over each block in the set.
               for k=1:length(SS.BLOCKS)
                  NEW_BLOCK = GetBlockInfo(SS.BLOCKS(k));
                  % Make sure to write back to the handle object.
                  % If we just update BLOCK, it will not get saved.
                  if isempty(NEW_BLOCK); break; end;
                  FILE.StimSets(j).BLOCKS(k) = NEW_BLOCK;
               end
            end
            
            function BLOCK = GetBlockInfo(BLOCK)
               if Ridx > length(Ref_array); Ridx=1; Cidx=Cidx+1; end
               if Cidx > length(Cur_array); Cidx=1; Sidx=Sidx+1; end
               if Sidx > length(Stim_array); BLOCK=[]; return; end
               
               BLOCK.Ref_E = Ref_array(Ridx);
               BLOCK.Current_uA = Cur_array(Cidx);
               BLOCK.Stim_E = Stim_array(Sidx);
               
               Ridx = Ridx+1;
            end
            
         end
         
      end %MatchStimBlockInfo()
      
      
      % ======================================================================
      % This calls the 3 main functions below to 1) parse all stim events
      % into stim blocks, 2) parse stim blocks into stim sets, 3) combine
      % stim sets into a single array for easy viewing.
      function ALL_SETS = ParseCombineSetsBlocks(ARY)
         if nargin < 1; ARY=StimSet_class.GetSetARY; end
         StimSet_class.ParseAllStimBlocks();
         StimSet_class.ParseAllStimSets();
         ALL_SETS = StimSet_class.CombineAllSets();
      end
      
      % Combine all of the STIM_SETS from each file in the SBI_Array, into a
      % single output set.
      function ALL_SETS = CombineAllSets(ARY)
         if nargin < 1; ARY=StimSet_class.GetSetARY; end
         
         % Create empty struct array, based on the STIM_SETS struct.
         % Be sure to add the
         ALL_SETS = []; %EmptyStructArray(SBIA(1).STIM_SETS);
         
         for i=1:length(ARY)
            SS = ARY(i);
            SSETS = SS.StimSets;
            if isempty(SSETS); continue; end
            % Add a couple more fields, to keep track of which file this came
            % from. deal() distributes same value to each array item.
            [SSETS.dirdate] = deal(SS.dirdate);
            [SSETS.fname] = deal(SS.fname);
            
            % Append STIM_SETS from this file to the full output array.
            if isempty(ALL_SETS); ALL_SETS = SSETS;
            else ALL_SETS = [ALL_SETS SSETS]; end
         end
      end
      
      
      
      % ======================================================================
      % ======================================================================
      %
      % Step 1, parse the raw STIM events (literally one event for every single
      % stimulation pulse) into "blocks", where a BLOCK is a consecutive series
      % of stim pulses, usually about 1 second long at 200Hz (which, due to
      % timing, almost always comes out to 198 pulses), so spaced about 5mSec
      % apart.
      %
      function ParseAllStimBlocks()
         % This just loops over the StimSet array (there is one entry for
         % each data file) and calls ParseStimBlocks for each file.
         ARY = StimSet_class.GetSetARY;
         
         for i=1:length(ARY)
            SS = ARY(i);
            fprintf('%d;', i);
            if ~mod(i,20); fprintf('\n'); end
            SS.StimBlocks = StimSet_class.ParseStimBlocks(SS);
         end
         fprintf('\n');
      end
      
      function BLOCKS=ParseStimBlocks(SS, varargin)
          if ~isempty(varargin{3})
              handles = varargin{3};
              handles.ProgressBar.Title.String = 'Loading Stim Blocks';
              drawnow
          else
              handles = [];

          end
         % Parse stim pulses in each file in into an array of stim blocks.
         % This is an "intermediate/temporary" array that we will then
         % parse into Stim Sets.
         %
         % To keep the outer for() loop simple and neat, we break up the
         % functionality into separate sub-fuctions, NewBlock(), AddToBlock(),
         % and EndBlock(). Makes it much easier to understand what is
         % happening.
         
         BLOCKS=[];
         if isempty(SS.Stim_ticks); return; end
         
         % A gap of >= 10mSec ends the block. Typically, our stim interval is
         % 5mSec (200Hz).
         NTICKS_gap = 0.010 * SS.TicksPerSec; % Default 10mSec gap.

         % Min and max ticks allowed in a block, else it is rejected.
         MinNStim = 104;
         MaxNStim = 229;
         
         for k=1:2:length(varargin(1:2))
            NAME = varargin{k};
            if k+1 > length(varargin); fprintf('Value not provided for argument "%s"\n', NAME); break; end
            VAL = varargin{k+1};
            switch(lower(varargin{k}))
               case {'maxstimgap_ms'}
                  NTICKS_gap = double(VAL) * SS.TicksPerSec;
               case {'minnstim'}
                  MinNStim = double(VAL);
               case {'maxnstim'}
                  MaxNStim = double(VAL);
               otherwise
                  fprintf('BAD named parameter: %s\n', varargin{k});
            end
         end
         
         ticks = SS.Stim_ticks;
         %ampl = SBI.STIM_INFO.ampl;
         
         BLOCK=[];
         LastTick = ticks(1);
         % Loop across all tick events in this data file.
         for i=1:length(ticks)
            AddToBlock(i);
            LastTick = ticks(i);
            if ~isempty(handles)
                handles.PBarObj.Position(3) = i/length(ticks)*1000;
                handles.PBarTxt.String = [num2str(round(i/length(ticks)*100)),'%'];
                drawnow
            end
         end
         if ~isempty(handles)
            handles.ProgressBar.Title.String = ['Finished'];
            handles.PBarObj.Position(3) = 0;
            handles.PBarTxt.String = 'Waiting';
            drawnow
         end
         EndBlock();
         
         % SUB-FUNCTIONS to process each stim event into the blocks.
         function NewBlock(tick_idx)
            BLOCK=[];
            % In EndBlock() we will replace the idx with the tick number.
            BLOCK.Start_tick = tick_idx;
            BLOCK.End_tick = tick_idx; % This will keep getting updated.
            BLOCK.NStim = 0;
            BLOCK.AmplAvg = 0;
         end
         
         function EndBlock()
            % Too short? (likely a "learning" block)
            % Note that some files have a short "first block", so we've
            % lowered this from 170 down to 105.
            if BLOCK.NStim <= MinNStim || BLOCK.NStim > MaxNStim; return; end
            
            % Initially, these "_tick" items are actually indices, so that
            % we can index into the amplitudes array to find an average
            % artifact amplitude (using the middle 50% of the sorted
            % amplitudes, to eliminate outliers).
            %BLOCK.AmplAvg = trimmean(ampl(BLOCK.Start_tick:BLOCK.End_tick), 50);
            
            % Now we can convert the index values into ticks.
            BLOCK.Start_tick = ticks(BLOCK.Start_tick);
            BLOCK.End_tick = ticks(BLOCK.End_tick);
            
            % For convenience, compute the start time in seconds, so it will
            % be easier to later look for gaps between blocks.
            BLOCK.Start_sec = double(BLOCK.Start_tick) / SS.TicksPerSec;
            BLOCK.End_sec = double(BLOCK.End_tick) / SS.TicksPerSec;
            
            if isempty(BLOCKS); BLOCKS = BLOCK; else; BLOCKS(end+1) = BLOCK; end
         end
         
         function AddToBlock(tick_idx)
            if isempty(BLOCK); NewBlock(tick_idx); end
            
            if ticks(tick_idx)-LastTick > NTICKS_gap
               EndBlock();
               NewBlock(tick_idx);
            end
            
            BLOCK.End_tick = tick_idx;
            BLOCK.NStim = BLOCK.NStim+1;
         end
      end %ParseStimBlocks()
      
      
      % ======================================================================
      % ======================================================================
      %
      % Step 2, parse the STIM blocks into Stim Sets. Typically a "Set" is a
      % single run through a particular stimulation script (unipolar, bipolar,
      % whatever).
      %
      % Strategy and implementation is similar to parsing for Stim Blocks, except
      % we look for a gap of around 3.5 seconds, instead of 10mSec.
      
      function ParseAllStimSets()
         ARY = StimSet_class.GetSetARY;
         % Loop through the StimSet array (essentially looping the data
         % files), and for each of the StimBlocks, parse and save the
         % StimSets struct.
         for i=1:length(ARY)
            SS = ARY(i);
            fprintf('%d;', i);
            if ~mod(i,20); fprintf('\n'); end
            SS.StimSets = StimSet_class.ParseStimSets(SS);
         end
         fprintf('\n');
      end
      
      
      function SETS=ParseStimSets(SS)
         % Parse the STIM_BLOCKS into groups, "Sets", based on gaps in the
         % stimulation. Each SET struct will also contain the list of BLOCKS
         % that are in the set.
         SETS=[];
         if isempty(SS.StimBlocks); return; end
         
         % Parse each STIM_BLOCKS item (one per CED data file file) and look
         % for whole "sets" of stimulation blocks. A gap of >3.5 seconds
         % terminates a set.
         SetGap_sec = 3.5;  % Look for >3.5 second gap.
         LastTime = SS.StimBlocks(1).Start_sec;
         SET=[];
         
         for i=1:length(SS.StimBlocks)
            AddToSet(SS.StimBlocks(i));
            LastTime = SS.StimBlocks(i).Start_sec;
         end
         EndSet();
         
         % ========================================================
         % Sub-functions to handle the Set array.
         function NewSet(SBLK)
            SET=[];
            SET.BLOCKS=StimSet_class.EmptyStructArray(SBLK);
            SET.Start_sec = SBLK.Start_sec;
            SET.End_sec = [];
            SET.TotalDur_sec=0;
            SET.NumStimBlocks=0;
         end
         
         function EndSet()
            if SET.NumStimBlocks < 20; return; end;  % Don't save small sets.
            SET.TotalDur_sec = SET.End_sec - SET.Start_sec;
            if isempty(SETS); SETS = SET;
            else SETS(end+1) = SET; end
         end
         
         function AddToSet(SBLK)
            if isempty(SET); NewSet(SBLK); end
            if SBLK.Start_sec - LastTime > SetGap_sec
               EndSet();
               NewSet(SBLK);
            end
            
            SET.BLOCKS(end+1) = SBLK; % Add to the array of blocks.
            SET.End_sec = SBLK.Start_sec; % We "fudge" this and do not add in the block duration. Not critical.
            SET.NumStimBlocks = SET.NumStimBlocks+1;
         end
      end
      
      % Quick trick by Walter Robson using typically obscure Matlab array /
      % struct / cell jiu jitsu to create a zero-length struct array, based on
      % existing struct.
      %
      %   https://www.mathworks.com/matlabcentral/answers/122057-how-to-create-an-empty-struc-with-fields-of-a-given-struct
      %
      function ES = EmptyStructArray(InStruct)
         f = fieldnames(InStruct)';
         f{2,1}={};  % Add a second empty row, for the "value" of each field.
         ES = struct(f{:});
      end

      
      % Save the master array of StimSet objects to a file.
      function SaveMasterFile(ARY, fname)
         if nargin < 1; ARY = StimSet_class.GetSetARY; end
         if ~isa(ARY, 'StimSet_class') || (length(ARY) < 2)
            error('SaveMasterFile() requires a vector of StimSet_class objects');
         end
         if nargin < 2; fname='StimSet_Master'; end
         StimSet_array = ARY;
         save(fname, 'StimSet_array');
      end
      
      % Load the master file, and set the persistent ARY value.
      function ARY=LoadMasterFile(fname)
         if nargin < 1; fname='StimSet_Master'; end
         FLOAD = load(fname);
         ARY = FLOAD.StimSet_array;
         StimSet_class.GetSetARY(ARY);
      end
      
      
      % Take either a SpikeDATA MAT file name, or SpikeDATA object,
      % and return a new StimSet object.
      function ITEM = NewStimSetObject(fname_or_obj)
         
         if isstr(fname_or_obj)
            DATA = SpikeDATA_class.LoadMAT(LOAD_fname);
         else
            DATA = fname_or_obj;
         end
         
         ITEM = StimSet_class();  % Allocate new SBI_Array object.
         ITEM.TicksPerSec = DATA.TicksPerSec;
         ITEM.STM = DATA.STM;
         ITEM.Stim_ticks = DATA.Stim_ticks;
      end

      % Loop through all of the data folders, and load the
      % *.MAT files into one Master MAT file.
      function ARY=LoadAllFiles(ARY)
         MATroot='c:\prj\Ross7\SpikeMark\Data\MoMo\';
         if nargin < 1; ARY = StimSet_class.GetSetARY; end
         if isempty(ARY); ARY=StimSet_class.empty; end

         % Do a "wildcard" dir, looking for just the single-unit data folders, which
         % all end with the same words.
         dirs = dir([MATroot '*single unit recording']);
         
         DateStart = datetime(2000,1,1);
         DateEnd = datetime(2050, 1, 1); % Stop before we get to files with TextMark.
         
         % *** LOOP over Directories ***
         %for d_idx = 1:length(dirs)
         for d_idx = length(dirs):-1:1   % Reverse Loop, newest to oldest.
            dname = dirs(d_idx).name;
            fprintf('Looking at dir: %s\n', dname);
            
            if length(dname) < 7; continue; end
            YEAR = str2double(dname(1:4)); MONTH = str2double(dname(6:7)); DAY = str2double(dname(9:10));
            DIRDate = datetime(dname(1:10)); % Create datetime object from directory string.
            
            % Make sure directory date is in the user's desired range.
            if DIRDate < DateStart; continue; end
            if ~isempty(DateEnd) && DIRDate > DateEnd; continue; end
            dirdate_str = sprintf('%d-%02d-%02d', YEAR, MONTH, DAY);
            
            % Get list of ONLY the MAT files.
            files = dir([MATroot dname '\*.MAT']);
            
            % *** LOOP over Files ***
            % LOOP through all of the SBI MAT files in the folder, and add
            % their data to the Master SBI array.
            for f_idx = 1:length(files)
               if files(f_idx).bytes < 100; continue; end
               fname = files(f_idx).name;
               
               if StimSet_class.ExistsInArray(ARY, dirdate_str, fname)
                  fprintf('  File exists!! Skipping: %d-%02d-%02d: %s\n', YEAR, MONTH, DAY, fname);
                  continue;
               else
                  fprintf('  Processing file: %d-%02d-%02d: %s\n', YEAR, MONTH, DAY, fname);
               end
               
               LOAD_fname = [files(f_idx).folder '\' fname];
               ITEM = StimSet_class.NewStimSetObject(LOAD_fname);
               
               DATA = SpikeDATA_class.LoadMAT(LOAD_fname);

               % Create entry in ARY[] array.
               ITEM = StimSet_class();  % Allocate new SBI_Array object.
               ITEM.dir = dname;
               ITEM.dirfull = [MATroot dname];
               ITEM.dirdate = dirdate_str;
               ITEM.fname = fname;
               ITEM.TicksPerSec = DATA.TicksPerSec;
               ITEM.STM = DATA.STM;
               ITEM.Stim_ticks = DATA.Stim_ticks;
               %ITEM.SBI = SBI;
               %ITEM.STIM_INFO = FILE_LOAD.STIM_INFO;
               
               % Append to array.
               ARY(end+1) = ITEM;
               
            end % for-files
         end % for-directory
         
         % Save array persistently.
         StimSet_class.GetSetARY(ARY);
      end %LoadAllFiles
      
   end %methods(Static)
end

