% File name: SetUpL7_4FlashPAwUST_ExternalTrig.m - Example of acquisition with ExtTrig
%
% Description:
%   Sequence programming file for L7-4 VVTnear array, using a plane wave
%   transmit and 2-1 synthetic aperture acquisition with an external trigger.
%   All 128 transmit and 64 receive channels are active for each synthetic
%   aperture acquisition. The external trigger-in, reconstruction/processing
%   are synchronoused with respect to acquisition. To run this script with
%   the VDAS hardware, one needs to connect an external trigger source to
%   the first trigger input.
%
% Last update:
% 04/25/2023
close all
clear all

P.path = 'C:\Users\Administrator\Documents\DataTemp\Ultrasound Data';
P.filePrefix = 'SaveIQData';

P.time = clock;
P.dateStr = strcat('_',num2str(P.time(2)), '-', num2str(P.time(3)), '-',...
    num2str(P.time(1)));

P.saveAcquisition = 1; %Default doesn't save
P.settingsNumber = 1; %Which version of the settings are you on?
P.settingsChanged = 1; %Whether the settings have changed since the last save.  Starts at 1 so that it doesn't automatically iterate to 2

P.runNumber = 1; %What run on the current setting?
P.itNumber = 1; %What iteration on the current run?


% Specify P structure array.
%% === Set commonly modified parameters ========================================

P(1).startDepth = 0;          % Acquisition start depth in wavelengths
P(1).endDepth = 192;          % Acquisition end depth
P(2).startDepth = 0;          % Acquisition start depth in wavelengths
P(2).endDepth = 192;          % Acquisition end depth may be different for PA
% Set 2D parameters
na = 7;      % Set na = number of flash angles for 2D.
    if na>1
        dtheta2D = (36*pi/180)/(na-1); % set dtheta2D to range over +/- 18 degrees.
    else
        dtheta2D = 0; % for a single angle, the angle will be zero.
    end

% Set PA parameters
oneway = 1;     % (logical) oneway=1 turns off the transmitters by setting TX.Apod to zero for all transmitters
flash2Qdelay = 188; % microseconds between trigger input and start of acquisition (which outputs a trigger pulse at time=0)
ne = 2;         % ne = number of acquisitions in PA ensemble for coherent addition in I/Q buffer. Capture the left and the right
PA_Angle = 0;   % angle of transmit plane wave for use in testing the PA mode in simulation
PA_PRF = 10;   % PA PRF in Hz. To be set in accordance with laser rep rate, when not using the input trigger mode.
                % When using the input trigger mode, remove the TTNA (SeqControl(4)) from the PA events to avoid
                % "missed TTNA" messages.

    if oneway==1
        disp(' *** PhotoAcoustic mode: Using one-way receive-only reconstruction ***')
    else
        disp(' *** Ultrasound Transmit mode: Using conventional T/R reconstruction ***')
    end
% ===============================================================================
% Specify system parameters
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.simulateMode = 2; %forces simulate mode, even if hardware is present.
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.verbose = 0;
Resource.Parameters.attenuation = -0.5;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.sizeApod = 64;

% Shows hardware memory allocation while running VSX
%Resource.VDAS.halDebugLevel = 0;
%Resource.VDAS.dmaTimeout = 20000; % in ms
%Resource.System.UTA = '260-S';
% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans.frequency = 5.208; % default for L7-4
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

%% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1,1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(1,3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
% - PA PData structure
PData(2).PDelta = [Trans.spacing, 0, 0.5];
PData(2).Size(1,1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3)); % PA window rows
PData(2).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1));  % PA window columns
PData(2).Size(1,3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(2).startDepth]; % x,y,z of upper lft crnr.

%% Specify Media object and point displacement function
pt1;
%Media.attenuation = Resource.Parameters.attenuation;
Media.function = 'movePoints';

%% Specify Resources.
% - RcvBuffer(1) is for both 2D and PA acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*(na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 10;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for PA reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for PA.
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';    % image buffer for 2D
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for PA image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for PA
Resource.ImageBuffer(2).numFrames = 1;          % 
% DisplayWindow is for 2D combined with PA
Resource.DisplayWindow(1).Title = mfilename;
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);

Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap; %gray(256) grayscaleCFImap
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
%% Specify Transmit waveforms structure
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - PA transmit waveform
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,0.67,6,1];

%% Specify Transmit beams structure
% 1st - PA left, na - US, na+1 = PA right
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [PA_Angle,0.0], ...
                   'Apod' , zeros(1,Trans.numelements,'int16'), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 1); % na TXs for 2D + 1 for PA

% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta2D;
else
    startAngle = -fix(na/2)*dtheta2D;
end

for n = 2:na+1   % na transmit events for 2D
    TX(n).waveform  = 2;
    TX(n).Apod  = ones(1,Trans.numelements);
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

%% -- only 02 TX struct needed for PA
TX(na+ne).waveform = 1;
TX(na+ne).Apod = zeros(1,Trans.numelements,'int16');     % THIS COMMAND TURNS OFF ALL TRANSMITTERS AND INVOKES THE RECEIVE-ONLY BEAMFORMER
TX(na+ne).Steer = [PA_Angle,0.0];            % only relevant when transmitters are active
TX(na+ne).Delay = computeTXDelays(TX(na+ne)); % only relevant when transmitters are active

%% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 3;

% This allows one to use different transmit profile for PA ... only relevant if transmitters are active
% --- currently TPC(2) is not used ---
TPC(2).name = 'PA';
TPC(2).maxHighVoltage = 1.7;

%% Analog front end gain settings.
RcvProfile(1).LnaGain = 12;     % 12, 18, or 24 dB  (18=default)
RcvProfile(1).condition = 'immediate';

RcvProfile(2).LnaGain = 24;
RcvProfile(2).condition = 'immediate';


%% Specify Receive structure arrays.
%   We need to acquire all the 2D and PA data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need 2*na Receives for a 2D frame and ne Receives for a PA frame.
maxAcqLngth2D = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxAcqLngthPA = ceil(sqrt(P(2).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wl4sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.

Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', maxAcqLngth2D, ...
                        'TGC', 1, ... 
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0,...
                        'inputFilter',[+0.00000 +0.00122 +0.00400 +0.00687 +0.00677 +0.00378 +0.00323 ...
                        +0.00931 +0.01569 +0.00861 -0.01431 -0.03320 -0.02631 -0.00491 ...
                        -0.01303 -0.07413 -0.14075 -0.12170 +0.02014 +0.20425 +0.28894]), 1, (na+ne)*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1; % move points before doing ensemble of different angle plane waves
    
    % PA acquisitions
    % -- 1st synthetic aperture acquisition for PA left.
    j=1;
    Receive(j+k).Apod(1:Resource.Parameters.numRcvChannels) = 1.0;
    Receive(j+k).framenum = i;
    Receive(j+k).acqNum = j;   % two acquisitions per frame
    Receive(j+k).startDepth = P(2).startDepth;
    Receive(j+k).endDepth = P(2).startDepth + wl4sPer128*ceil(maxAcqLngthPA/wl4sPer128);
    Receive(j+k).TGC = 2;           % TGC(1) is tied to the GUI sliders
    if j==1, Receive(j+k).callMediaFunc = 1; end % move points between 2D and PA to see difference in simulation
    
    % acquisitions for 2D
    for j = 2:na+1
        if mod(j,2) == 1 
            % -- 1st synthetic aperture acquisition for full frame.
            Receive(j+k).Apod(1:Resource.Parameters.numRcvChannels) = 1.0;
        else
            % -- 2nd synthetic aperture acquisition for full frame.
            Receive(j+k).Apod((Resource.Parameters.numRcvChannels+1):Trans.numelements) = 1.0;
        end
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    
    % PA acquisitions
    % -- 2nd synthetic aperture acquisition for PA right.
    j=na+ne;
    Receive(j+k).Apod((Resource.Parameters.numRcvChannels+1):Trans.numelements) = 1.0;

    Receive(j+k).framenum = i;
    Receive(j+k).acqNum = j;   % two acquisitions per frame
    Receive(j+k).startDepth = P(2).startDepth;
    Receive(j+k).endDepth = P(2).startDepth + wl4sPer128*ceil(maxAcqLngthPA/wl4sPer128);
    Receive(j+k).TGC = 2;           % TGC(1) is tied to the GUI sliders
end

%% Specify TGC Waveform structures.
% - 2D TGC
TGC(1).CntrlPts = [50,141,275,404,510,603,702,782];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - PA TGC
TGC(2).CntrlPts = [50,141,275,404,510,603,702,782]; % TO BE MODIFIED HERE AFTER COLLECTING PA DATA
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

%% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for PA. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,1)), 1, 2);

% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:na) = 2:na+1;  % na ReconInfos needed for na angles
% - Set Recon values for PA ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [2,-1];
Recon(2).RINums(1,1:ne) = [1,na+ne];   % 'ne' ReconInfos needed for PA ensemble.

%% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For PA, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % 4=accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).mode = 'replaceIQ';          % 3=replace IQ data (expect to use mode 5 on last acquisition)
for j = 2:na+1
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
end
if na>1
    ReconInfo(na+1).mode = 'accumIQ_replaceIntensity';     % 5=Reconstruct IQ data, add values to InterBuffer and compute magnitude, replacing data in ImageBuffer.
else
    ReconInfo(na+1).mode = 'replaceIntensity';     % 0=replace IQ data, detect, and replace Intensity data in ImageBuffer. (single acquisition)
end

%  - ReconInfos for PA ensemble.

for j = [1,na+ne]
    if j==1, ReconInfo(j).mode = 'replaceIQ'; end
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
end
if ne>1
    ReconInfo(na+ne).mode = 'accumIQ_replaceIntensity'; % 5=accum and detect
else
    ReconInfo(na+ne).mode = 'replaceIntensity'; % 0=replace IQ data, detect, and replace Intensity data;  1=Add the new reconstructed intensity data to the data in the ImageBuffer
end

%% Specify Process structure arrays.
cpt = 22;       % define here so we can use in UIControl below
cpers = 80;     % define here so we can use in UIControl below

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',10,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',50,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};

Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',10,...      % reject level
                         'persistMethod','dynamic',...
                         'persistLevel',cpers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

Process(3).classname = 'External';
Process(3).method = 'showSingleRF';
Process(3).Parameters = {'srcbuffer','receive',... % name of buffer to process. 
                         'srcbufnum',1,... 
                         'srcframenum',1,... 
                         'dstbuffer','none'}; 

Process(4).classname  = 'External';
Process(4).method = 'saveData';
Process(4).Parameters =  {'srcbuffer', ...
                            'inter',... 
                            'srcbufnum',1,...
                            'srcframenum',1,...
                            'dstbuffer','none'};

%% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;
% -- Change to Profile 2 (PA)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and PA ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 3000; % time in usec
% -- PRF for PA ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(PA_PRF*1e-5)); % (10 msecs for PA_PRF=100 Hz)
% -- Change to Profile 1 (2D)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;
% -- Time between PA and next 2D acquisition. Set to allow time for profile change.
SeqControl(6).command  = 'timeToNextAcq';
SeqControl(6).argument = 3000; % time in usec
% -- Jump back to start.
SeqControl(7).command  = 'jump';
SeqControl(7).argument = 1;
% set receive profile
SeqControl(8).command  = 'setRcvProfile';
SeqControl(8).argument = 1;
SeqControl(9).command  = 'setRcvProfile';
SeqControl(9).argument = 2;
% input trigger
SeqControl(10).command   = 'triggerIn';
SeqControl(10).condition = 'Trigger_1_Falling'; % Trigger input 1, enable with rising edge %% 'Trigger_1_Rising'
SeqControl(10).argument  = 1; % 500 msec timeout delay
    % (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
% noop delay between trigger in and start of acquisition
SeqControl(11).command  = 'noop';
SeqControl(11).argument = fix(flash2Qdelay)*1; % noop counts are in 2 microsec increments fix(flash2Qdelay)*5
% sync command
SeqControl(12).command  = 'sync';
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 13;  % next SeqControl number

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 2;

%% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire PA ensemble left.
    j=1;
    % Wait for input trigger from flash lamp firing
    Event(n).info = 'Wait for Trigger IN';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 10;
    n = n+1;

    % PA acquisition left
    Event(n).info = 'Acquire PA event';
    Event(n).tx = 1;
    Event(n).rcv = (na+ne)*(i-1)+j;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;

%     % Delay after PA left
%     Event(n).info = 'Delay';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = [11,12];
%     n = n+1;

    % Acquire 2D frame
    for j = 2:(na+1)
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j; 
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = [3,9];   % replace last 2D acquisition Event's seqControl (longer TTNA and new RCV profile)
    
    % Acquire PA ensemble right.
    % Wait for input trigger from flash lamp firing
    Event(n).info = 'Wait for Trigger IN';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 10;
    n = n+1;

    % PA acquisition right
    j = na+ne;
    Event(n).info = 'Acquire PA event';
    Event(n).tx = na+ne;
    Event(n).rcv = (na+ne)*(i-1)+j;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n-1).seqControl = [6,8]; % replace last PA acquisition Event's seqControl with longer TTNA and RCV profile change
  
    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer (each needs a different value of nsc)
    nsc = nsc+1;

    Event(n).info = 'Show single RF';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 12;
    n = n+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'PA image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab only every 3rd frame to prevent slowdown
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

%if P.saveAcquisition == 1
Event(n).info = 'Save image data';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 4;
Event(n).seqControl = 0;
%end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;


%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');

% - Color Priority Threshold Slider
UI(2).Control = {'UserB2','Style','VsSlider','Label','Color Priority','SliderMinMaxVal',[0,255,cpt],...
                 'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%-UI#2Callback');

% - Color Persistence Slider
UI(3).Control = {'UserB1','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,cpers],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%-UI#3Callback');

% - Save on/off button
UI(4).Control = {'UserC3','Style','VsToggleButton','Label','Save On/Off'};
UI(4).Callback = text2cell('%SaveToggle');

% - Show single RF
import vsv.seq.uicontrol.VsSliderControl
nr = Resource.Parameters.numRcvChannels; 
UI(5).Control = {'UserB3','Style','VsSlider',...
                          'Label','Plot Channel',...
                          'SliderMinMaxVal',[1,64,32],...
                          'SliderStep', [1/nr,8/nr],...
                          'ValueFormat', '%3.0f'};
UI(5).Callback = @setChannelNum; 

EF(1).Function = text2cell('%EF#1%');
EF(2).Function = vsv.seq.function.ExFunctionDef('showSingleRF',@showSingleRF); 

%% File end

%% Save all the structures to a .mat file, and run VSX automatically
filename = ('L7-4vFlashPA');   % define variable 'filename' to permit VSX to skip user query for matfile
save (filename)                 % save the structures to a matfile
VSX                             % invoke VSX automatically when running this Setup script

return


%% **** Callback routines to be encoded by text2cell function. ****
% ---------------------------------------------------------------------
%-UI#1Callback - Sensitivity cutoff change
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base','Control', Control);
    return
%-UI#1Callback


%-UI#2Callback - Color Threshold change
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'threshold'), Process(2).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.threshold.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',2,'threshold',UIValue};
    assignin('base','Control', Control);
%-UI#2Callback

%-UI#3Callback - Color Persistence change
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'persistLevel'), Process(2).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',2,'persistLevel',UIValue};
    assignin('base','Control', Control);
%-UI#3Callback

%SaveToggle
% Turns on/off saving.  Gets the base file name and destination folder if
% saving has been turned on.
P = evalin('base','P');
P.saveAcquisition = ~P.saveAcquisition;

if P.saveAcquisition
    
    %Dialogue box to make a new base file name
    prompt={'File Name:'};
    dlgTitle = 'Save Parameters';
    numLines = 2;
    defaultAns = {P.filePrefix};
    
    %NO CHECKS on the file name
    userInput = inputdlg(prompt,dlgTitle,numLines,defaultAns);
    P.filePrefix = userInput{1};
    
    %Assign a new directory
    P.path = strcat(uigetdir(P.path),'/');
else
    %Every toggle is a new run
    P.runNumber = P.runNumber +1;
    
    %And restarts the iterations
    P.itNumber = 1;
end

assignin('base','P',P);

%SaveToggle


%EF#1%
saveData(IQData)
    
    if evalin('base','P.saveAcquisition')
                %% File Naming and IQ data
        %Read in the misc variables struct 
        P = evalin('base','P');

        %Now we want to handle the settings file, to make sure we are saving
        %correctly

        %If the settings have changed since the last time, reset the boolean to
        %false so that new changes will propogate.  Also reset the run number
        %since it's the first run on the new settings
        if P.itNumber == 1
            %Check for previous settings
            while exist(strcat(P.path,P.filePrefix,P.dateStr,...
            '_Run-',int2str(P.runNumber),'_It-',int2str(P.itNumber),'_IQ.mat'),'file')
                P.runNumber = P.runNumber+1;
            end
            
           %TODO: Line to invoke save preSet here.
        end

        %Calculate the file name for any iteration specific file
        fileName = strcat(P.path,P.filePrefix,P.dateStr,...
            '_Run-',int2str(P.runNumber),'_It-',int2str(P.itNumber));
        
        %Save the IQ data for the run.
        save(strcat(fileName,'_IQ'),'IQData'); %Save the IQ data     
        
        %Modify the iteration number
        P.itNumber = P.itNumber+1;
        assignin('base','P',P);
    end
return
%EF#1%


function setChannelNum(~,~,UIValue) 
    assignin('base', 'myPlotChnl', round(UIValue));
end

function showSingleRF(RData)
    persistent myHandle
    % If ‘myPlotChnl’ exists, read it for the channel to plot.
    if evalin('base','exist(''myPlotChnl'',''var'')')
        channel = evalin('base','myPlotChnl');
    else
        channel = 32; % Channel no. to plot
    end
    % Create the figure if it doesn’t exist.
    if isempty(myHandle)||~ishandle(myHandle)
        figure;
        myHandle = axes('XLim',[0,1500],'YLim',[-16384 16384], ...
                        'NextPlot','replacechildren');
    end
    % Plot the RF data.
    plot(myHandle,RData(:,channel));
    drawnow
end