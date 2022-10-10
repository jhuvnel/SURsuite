% Convenience function to load the Spike2 SMR data file, and also parse the
% stimulation blocks and add the 'StimBlocks' field to the DAT object. Also
% returns the Z=ZapQuick_class object.
%
% Also returns AZT, All Zap Templates. This is a matrix of all of the
% templates for each stimulation block in the file.
%
% The filename input can be a full path to the file, or a cell array of
% path partial components.
%
function [DAT,Z, AZT]=QuickLoadBrian(filename, handles)

warning('off','signal:findpeaks:largeMinPeakHeight');

if iscell(filename)
   BASE = '\\10.16.39.7\labdata\Monkey Single Unit Recording\';
   %BASE = 'c:\downloads\SURData\';
   SUBJECT = [filename{1} ' single unit recording\'];
   DATE = [filename{2} ' ' filename{1} ' Single Unit Recording\'];
   FILE = [filename{3} '.smr'];
   filename = [BASE SUBJECT DATE FILE];
end

%fprintf('Filename: "%s"\n', filename);
%return;

% Load CED/Spike2 *.SMR the data file.
%DAT = SpikeDATA_class.LoadSMRFile(filename);
%DAT = SpikeDATA_class.LoadSMRFile('..\Nancy\2019-12-04\Data2.smr');
if isempty(handles)
f = msgbox({'Loading SMR File';'Message Will Close When Complete'});
else
    handles.ProgressBar.Title.String = 'Loading SMR File';
    drawnow
end
DAT = SpikeDATA_class.LoadSMRFile(filename, handles);
if ~isprop(DAT,'APorientation')
    DAT.addprop('APorientation');
    if isempty(handles)
        DAT.APorientation = 0;
    else
        DAT.APorientation = handles.APorientation;
    end
    
end

if isempty(handles)
delete(f)
end

SS = StimSet_class.NewStimSetObject(DAT);

% Parse the individual stimulation pulses into "blocks", based on the time
% between pulses. Any gap >10mSec will mark the end of a block.
SS.StimBlocks = StimSet_class.ParseStimBlocks(SS, 'MinNStim', 20, handles);
% Match up the stimulation Text Marks with the Blocks that we just parsed.
% We add the Stim_E, Ref_E, and Current_uA to the columns of information in
% the SS.StimBlocks array.
SS.MatchStimBlocks(handles);

% Now copy the newly parsed "StimBlocks" back into the DAT structure, so we
% can ignore SS afterwards and just deal with the DAT object.
DAT.addprop('StimBlocks');
DAT.StimBlocks = SS.StimBlocks;
clear SS

if isempty(handles)
TC = Template_class(DAT);

[DAGANtempl_unzeroed, DAGANtempl] = TC.GetAllArtifactTemplate('NerveCon', []);


% [DAGAN2templ_unzeroed, DAGAN2templ] = TC.GetAllArtifactTemplate('ZapOut');



[DAGANAPtempl_unzeroed, DAGANAPtempl] = TC.GetAllAPTemplate('NerveCon', []);
end


% [DAGAN2APtempl_unzeroed, DAGAN2APtempl] = TC.GetAllAPTemplate('ZapOut');

% Make our ZapQuick object.

%Z = ZapQuick_class(DAT);
%ZapAllBlocks(Z);
Z=[];
AZT=[];
% Return numeric array of All Zap Templates.

%AZT=GetAllTemplates(Z);
%DATA=DAT;
%save([DAT.fname,'.SUR'], 'DATA');

% 
% answer = questdlg('Would you like to run Splie Zap on all data?',...
%     'Spline Zap',...
%     'Yes','No','Yes');
% switch answer
%     case 'Yes'
%         ZapAllBlocks_Spline(Z);
%     case 'No'
%         Z.Values_zapped_spline = [];
%         Z.APLoc_spline = [];
% end

end
