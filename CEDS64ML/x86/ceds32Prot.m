function [methodinfo,structs,enuminfo,ThunkLibName]=ceds32Prot
%CEDS32PROT Create structures to define interfaces found in 'ceds64int'.

%This function was generated by loadlibrary.m parser version 1.1.6.38 on Mon Nov 24 17:34:24 2014
%perl options:'ceds64int.i -outfile=ceds32Prot.m'
ival={cell(1,0)}; % change 0 to the actual number of functions to preallocate the data.
structs=[];enuminfo=[];fcnNum=1;
fcns=struct('name',ival,'calltype',ival,'LHS',ival,'RHS',ival,'alias',ival);
ThunkLibName=[];
%  int S64FileCount (); 
fcns.name{fcnNum}='S64FileCount'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}=[];fcnNum=fcnNum+1;
%  int S64CloseAll (); 
fcns.name{fcnNum}='S64CloseAll'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}=[];fcnNum=fcnNum+1;
%  int S64Create ( const char * FileName , const int nChans , const int nBig ); 
fcns.name{fcnNum}='S64Create'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'cstring', 'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64Open ( const char * FileName , const int nFlag ); 
fcns.name{fcnNum}='S64Open'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'cstring', 'int32'};fcnNum=fcnNum+1;
%  int S64IsOpen ( const int nFid ); 
fcns.name{fcnNum}='S64IsOpen'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64Close ( const int nFid ); 
fcns.name{fcnNum}='S64Close'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64Empty ( const int nFid ); 
fcns.name{fcnNum}='S64Empty'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64GetFileComment ( const int nFid , const int nInd , char * Comment , const int nGetSize ); 
fcns.name{fcnNum}='S64GetFileComment'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring', 'int32'};fcnNum=fcnNum+1;
%  int S64SetFileComment ( const int nFid , const int nInd , const char * Comment ); 
fcns.name{fcnNum}='S64SetFileComment'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring'};fcnNum=fcnNum+1;
%  int S64GetFreeChan ( const int nFid ); 
fcns.name{fcnNum}='S64GetFreeChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64MaxChans ( const int nFid ); 
fcns.name{fcnNum}='S64MaxChans'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  double S64GetTimeBase ( const int nFid ); 
fcns.name{fcnNum}='S64GetTimeBase'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='double'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64SetTimeBase ( const int nFid , const double dSecPerTick ); 
fcns.name{fcnNum}='S64SetTimeBase'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'double'};fcnNum=fcnNum+1;
%  long long S64SecsToTicks ( const int nFid , const double dSec ); 
fcns.name{fcnNum}='S64SecsToTicks'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'double'};fcnNum=fcnNum+1;
%  double S64TicksToSecs ( const int nFid , const long long tSec ); 
fcns.name{fcnNum}='S64TicksToSecs'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='double'; fcns.RHS{fcnNum}={'int32', 'int64'};fcnNum=fcnNum+1;
%  int S64GetVersion ( const int nFid ); 
fcns.name{fcnNum}='S64GetVersion'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  long long S64FileSize ( const int nFid ); 
fcns.name{fcnNum}='S64FileSize'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  long long S64MaxTime ( const int nFid ); 
fcns.name{fcnNum}='S64MaxTime'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64TimeDate ( const int nFid , long long * pTDGet , const long long * pTDSet , int iMode ); 
fcns.name{fcnNum}='S64TimeDate'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int64Ptr', 'int64Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64AppID ( const int nFid , int * pTDGet , const int * pTDSet , int iMode ); 
fcns.name{fcnNum}='S64AppID'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32Ptr', 'int32Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64GetExtraData ( const int nFid , void * pData , unsigned int nBytes , unsigned int nOffset ); 
fcns.name{fcnNum}='S64GetExtraData'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'voidPtr', 'uint32', 'uint32'};fcnNum=fcnNum+1;
%  int S64SetExtraData ( const int nFid , const void * pData , unsigned int nBytes , unsigned int nOffset ); 
fcns.name{fcnNum}='S64SetExtraData'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'voidPtr', 'uint32', 'uint32'};fcnNum=fcnNum+1;
%  int S64ChanType ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64ChanType'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  long long S64ChanDivide ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64ChanDivide'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  double S64GetIdealRate ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64GetIdealRate'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='double'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  double S64SetIdealRate ( const int nFid , const int nChan , const double dRate ); 
fcns.name{fcnNum}='S64SetIdealRate'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='double'; fcns.RHS{fcnNum}={'int32', 'int32', 'double'};fcnNum=fcnNum+1;
%  int S64GetChanComment ( const int nFid , const int nChan , char * Comment , const int nGetSize ); 
fcns.name{fcnNum}='S64GetChanComment'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring', 'int32'};fcnNum=fcnNum+1;
%  int S64SetChanComment ( const int nFid , const int nChan , const char * Comment ); 
fcns.name{fcnNum}='S64SetChanComment'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring'};fcnNum=fcnNum+1;
%  int S64GetChanTitle ( const int nFid , const int nChan , char * Title , const int nGetSize ); 
fcns.name{fcnNum}='S64GetChanTitle'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring', 'int32'};fcnNum=fcnNum+1;
%  int S64SetChanTitle ( const int nFid , const int nChan , const char * Title ); 
fcns.name{fcnNum}='S64SetChanTitle'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring'};fcnNum=fcnNum+1;
%  int S64GetChanScale ( const int nFid , const int nChan , double * dScale ); 
fcns.name{fcnNum}='S64GetChanScale'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'doublePtr'};fcnNum=fcnNum+1;
%  int S64SetChanScale ( const int nFid , const int nChan , const double dScale ); 
fcns.name{fcnNum}='S64SetChanScale'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double'};fcnNum=fcnNum+1;
%  int S64GetChanOffset ( const int nFid , const int nChan , double * dOffset ); 
fcns.name{fcnNum}='S64GetChanOffset'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'doublePtr'};fcnNum=fcnNum+1;
%  int S64SetChanOffset ( const int nFid , const int nChan , const double dOffset ); 
fcns.name{fcnNum}='S64SetChanOffset'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double'};fcnNum=fcnNum+1;
%  int S64GetChanUnits ( const int nFid , const int nChan , char * Units , const int nGetSize ); 
fcns.name{fcnNum}='S64GetChanUnits'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring', 'int32'};fcnNum=fcnNum+1;
%  int S64SetChanUnits ( const int nFid , const int nChan , const char * Units ); 
fcns.name{fcnNum}='S64SetChanUnits'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'cstring'};fcnNum=fcnNum+1;
%  long long S64ChanMaxTime ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64ChanMaxTime'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  long long S64PrevNTime ( const int nFid , const int nChan , long long tFrom , long long tTo , const int n , const int nMask , const int nAsWave ); 
fcns.name{fcnNum}='S64PrevNTime'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64', 'int64', 'int32', 'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64ChanDelete ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64ChanDelete'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64ChanUndelete ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64ChanUndelete'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64GetChanYRange ( const int nFid , const int nChan , double * dLow , double * dHigh ); 
fcns.name{fcnNum}='S64GetChanYRange'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'doublePtr', 'doublePtr'};fcnNum=fcnNum+1;
%  int S64SetChanYRange ( const int nFid , const int nChan , const double dLow , const double dHigh ); 
fcns.name{fcnNum}='S64SetChanYRange'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double', 'double'};fcnNum=fcnNum+1;
%  int S64ItemSize ( const int nFid , const int nChan ); 
fcns.name{fcnNum}='S64ItemSize'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64SetEventChan ( const int nFid , const int nChan , const double dRate , const int iType ); 
fcns.name{fcnNum}='S64SetEventChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double', 'int32'};fcnNum=fcnNum+1;
%  int S64WriteEvents ( const int nFid , const int nChan , const long long * pData , const int nCount ); 
fcns.name{fcnNum}='S64WriteEvents'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64ReadEvents ( const int nFid , const int nChan , long long * pData , int nMax , const long long tFrom , const long long tTo , const int nMask ); 
fcns.name{fcnNum}='S64ReadEvents'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64Ptr', 'int32', 'int64', 'int64', 'int32'};fcnNum=fcnNum+1;
%  int S64SetMarkerChan ( const int nFid , const int nChan , const double dRate , const int nkind ); 
fcns.name{fcnNum}='S64SetMarkerChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double', 'int32'};fcnNum=fcnNum+1;
%  int S64WriteMarkers ( const int nFid , const int nChan , const S64Marker * pData , const int count ); 
fcns.name{fcnNum}='S64WriteMarkers'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'int32'};fcnNum=fcnNum+1;
%  int S64ReadMarkers ( const int nFid , const int nChan , S64Marker * pData , const int nMax , const long long tFrom , const long long tUpto , const int nMask ); 
fcns.name{fcnNum}='S64ReadMarkers'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'int32', 'int64', 'int64', 'int32'};fcnNum=fcnNum+1;
%  int S64EditMarker ( const int nFid , const int nChan , long long t , const S64Marker * pM ); 
fcns.name{fcnNum}='S64EditMarker'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64', 'S64MarkerPtr'};fcnNum=fcnNum+1;
%  int S64GetMaskCodes ( const int nMask , int * iCode , int * nMode ); 
fcns.name{fcnNum}='S64GetMaskCodes'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32Ptr', 'int32Ptr'};fcnNum=fcnNum+1;
%  int S64SetMaskCodes ( const int nMask , const int * nCode ); 
fcns.name{fcnNum}='S64SetMaskCodes'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32Ptr'};fcnNum=fcnNum+1;
%  int S64SetMaskMode ( const int nMask , const int nMode ); 
fcns.name{fcnNum}='S64SetMaskMode'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64GetMaskMode ( const int nMask ); 
fcns.name{fcnNum}='S64GetMaskMode'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64ResetMask ( const int nMask ); 
fcns.name{fcnNum}='S64ResetMask'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32'};fcnNum=fcnNum+1;
%  int S64ResetAllMasks (); 
fcns.name{fcnNum}='S64ResetAllMasks'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}=[];fcnNum=fcnNum+1;
%  int S64SetLevelChan ( const int nFid , const int nChan , const double dRate ); 
fcns.name{fcnNum}='S64SetLevelChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double'};fcnNum=fcnNum+1;
%  int S64SetInitLevel ( const int nFid , const int nChan , const int nLevel ); 
fcns.name{fcnNum}='S64SetInitLevel'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int32'};fcnNum=fcnNum+1;
%  int S64WriteLevels ( const int nFid , const int nChan , const long long * pData , const int nCount ); 
fcns.name{fcnNum}='S64WriteLevels'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64ReadLevels ( const int nFid , const int nChan , long long * pData , int nMax , const long long tFrom , const long long tTo , int * nLevel ); 
fcns.name{fcnNum}='S64ReadLevels'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64Ptr', 'int32', 'int64', 'int64', 'int32Ptr'};fcnNum=fcnNum+1;
%  int S64SetTextMarkChan ( const int nFid , const int nChan , double dRate , const int nMax ); 
fcns.name{fcnNum}='S64SetTextMarkChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double', 'int32'};fcnNum=fcnNum+1;
%  int S64SetExtMarkChan ( const int nFid , const int nChan , double dRate , const int nType , const int nRows , const int nCols , const long long tDiv ); 
fcns.name{fcnNum}='S64SetExtMarkChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'double', 'int32', 'int32', 'int32', 'int64'};fcnNum=fcnNum+1;
%  int S64GetExtMarkInfo ( const int nFid , const int nChan , int * nRows , int * nCols ); 
fcns.name{fcnNum}='S64GetExtMarkInfo'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int32Ptr', 'int32Ptr'};fcnNum=fcnNum+1;
%  int S64Write1TextMark ( const int nFid , const int nChan , const S64Marker * pData , const char * text , const int nSize ); 
fcns.name{fcnNum}='S64Write1TextMark'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'cstring', 'int32'};fcnNum=fcnNum+1;
%  int S64Write1RealMark ( const int nFid , const int nChan , const S64Marker * pData , const float * real , const int nSize ); 
fcns.name{fcnNum}='S64Write1RealMark'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'singlePtr', 'int32'};fcnNum=fcnNum+1;
%  int S64Write1WaveMark ( const int nFid , const int nChan , const S64Marker * pData , const short * wave , const int nSize ); 
fcns.name{fcnNum}='S64Write1WaveMark'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'int16Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64Read1TextMark ( const int nFid , const int nChan , S64Marker * pData , char * text , const long long tFrom , const long long tUpto , const int nMask ); 
fcns.name{fcnNum}='S64Read1TextMark'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'cstring', 'int64', 'int64', 'int32'};fcnNum=fcnNum+1;
%  int S64Read1RealMark ( const int nFid , const int nChan , S64Marker * pData , float * real , const long long tFrom , const long long tUpto , const int nMask ); 
fcns.name{fcnNum}='S64Read1RealMark'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'singlePtr', 'int64', 'int64', 'int32'};fcnNum=fcnNum+1;
%  int S64Read1WaveMark ( const int nFid , const int nChan , S64Marker * pData , short * real , const long long tFrom , const long long tUpto , const int nMask ); 
fcns.name{fcnNum}='S64Read1WaveMark'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'S64MarkerPtr', 'int16Ptr', 'int64', 'int64', 'int32'};fcnNum=fcnNum+1;
%  int S64SetWaveChan ( const int nFid , const int nChan , const long long tDiv , const int nType , const double dRate ); 
fcns.name{fcnNum}='S64SetWaveChan'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int64', 'int32', 'double'};fcnNum=fcnNum+1;
%  long long S64WriteWaveS ( const int nFid , const int nChan , const short * pData , const int count , const long long tFrom ); 
fcns.name{fcnNum}='S64WriteWaveS'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'int32', 'int16Ptr', 'int32', 'int64'};fcnNum=fcnNum+1;
%  long long S64WriteWaveF ( const int nFid , const int nChan , const float * pData , const int count , const long long tFrom ); 
fcns.name{fcnNum}='S64WriteWaveF'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'int32', 'singlePtr', 'int32', 'int64'};fcnNum=fcnNum+1;
%  long long S64WriteWave64 ( const int nFid , const int nChan , const double * pData , const int count , const long long tFrom ); 
fcns.name{fcnNum}='S64WriteWave64'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int64'; fcns.RHS{fcnNum}={'int32', 'int32', 'doublePtr', 'int32', 'int64'};fcnNum=fcnNum+1;
%  int S64ReadWaveS ( const int nFid , const int nChan , short * pData , const int nMax , const long long tFrom , const long long tUpto , long long * tFirst , const int nMask ); 
fcns.name{fcnNum}='S64ReadWaveS'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'int16Ptr', 'int32', 'int64', 'int64', 'int64Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64ReadWaveF ( const int nFid , const int nChan , float * pData , const int nMax , const long long tFrom , const long long tUpto , long long * tFirst , const int nMask ); 
fcns.name{fcnNum}='S64ReadWaveF'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'singlePtr', 'int32', 'int64', 'int64', 'int64Ptr', 'int32'};fcnNum=fcnNum+1;
%  int S64ReadWave64 ( const int nFid , const int nChan , double * pData , const int nMax , const long long tFrom , const long long tUpto , long long * tFirst , const int nMask ); 
fcns.name{fcnNum}='S64ReadWave64'; fcns.calltype{fcnNum}='cdecl'; fcns.LHS{fcnNum}='int32'; fcns.RHS{fcnNum}={'int32', 'int32', 'doublePtr', 'int32', 'int64', 'int64', 'int64Ptr', 'int32'};fcnNum=fcnNum+1;
structs.S64Marker.members=struct('m_Time', 'int64', 'm_Code1', 'uint8', 'm_Code2', 'uint8', 'm_Code3', 'uint8', 'm_Code4', 'uint8');
methodinfo=fcns;