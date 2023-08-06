%Interleaved Photoacoustic and Ultrasound Imaging
clear all
%% constants :
% filename - used by VSX to automatically load.mat file
% P - local constants
filename ='L7-4vFlashPA';

% Photoacoustic
P.startDepth = 5;
P.endDepth = 150;
P.paReconMethod ='saft';
P.paOverlayToggled = 0;
P.overlayImgHandle = 0;
P.nAvg = 5;

% Define system parameters.
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 64;
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.verbose = 3;
Resource.Parameters.attenuation = -0.5;
% Shows hardware memory allocation while running VSX
%Resource.VDAS.halDebugLevel = 0;
%Resource.VDAS.dmaTimeout = 20000; % in ms
Resource.System.UTA = '260-S';
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.


%% Media : Specify media points for simulation -
% MP - [x, y, z reflectivity ]
pt1; % External code , which sets a range of media points for simulate mode
Media.attenuation = Resource.Parameters.attenuation;
Media.function ='movePoints';

%% Trans : Properties of the transducer -
% name - Name of the transducer
% connType - type of scan head interface
% frequency - center frequency of the transducer in MHz
% units - wavelengths or mm
% computeTrans - for known transducer names
Trans.name ='L7-4';
Trans.frequency = 5.208; % default for L7-4
Trans.units ='mm';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 10;  % set maximum high voltage limit for pulser supply.

%% PData : Define pixel grid for image reconstruction -
% PDelta - Spacing between pixels( wavelenghts )
% Origin - x,y,z of upper left corner.
% Size - number of rows , columns , z- planes(in case of 3d imaging )
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2, ...
                    0, ...
                    P.startDepth];
PData(1).Size = [ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)), ...
                 ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)),1];

PData(1).Region = repmat(struct('Shape', ...
                         struct('Name', 'Rectangle', ...
                         'Position', [0, 0, P.startDepth], ...
                         'width', Trans.numelements * Trans.spacing, ...
                         'height', P.endDepth-P.startDepth)), 1, 1);
PData(1).Region = computeRegions(PData(1));
PData(2) = PData(1);

%% Resource.RcvBuffer : Stores raw RF data -
% rowsPerFrame - number of rows in frame
% colsPerFrame - number of columns in frame
% numFrames - number of frames(1 or even number )
Resource.RcvBuffer = repmat(struct(...
        'datatype', 'int16', ...
        'rowsPerFrame', 4096, ...
        'colsPerFrame', Resource.Parameters.numTransmit, ...
        'numFrames', 20), 1, 2);
Resource.RcvBuffer(2).rowsPerFrame = 2*4096;

%% Resource.InterBuffer : Storing I,Q pixels -
% datatype -'complex double'( default ) or'single'
% numFrames - number of frames
Resource.InterBuffer = repmat(struct(...
        'datatype', 'complex', ...
        'numFrames', 1), 1, 2);

%% Resource.ImageBuffer : Storing intensity pixels -
% datatype - double( default )
% numFrames - number of frames
Resource.ImageBuffer(1).datatype ='double';
Resource.ImageBuffer(1).numFrames = Resource.RcvBuffer(1).numFrames;
Resource.ImageBuffer(2).datatype ='double';
Resoucre.ImageBuffer(2).numFrames = Resource.RcvBuffer(2).numFrames;

%% Resource.DisplayWindow -
% pdelta - interpolating pixel data to this size for display
% Position - [x, y, width , height ],(x,y): lower left corner
% ReferencePt - upper left corner relative to transducer(in wavelengths )
% numFrames - number of frames stored in the cineloop image buffer
% Colormap - standard Matlab colormap
ScrnSize = get(0, 'ScreenSize');
DwPdelta = 0.35;
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/DwPdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/DwPdelta);
DwPos{1} = [100., (ScrnSize(4) -(DwHeight + 150.))/2., DwWidth, DwHeight];
DwPos{2} = [100.+1.5*DwWidth, (ScrnSize(4) -(DwHeight+150.)) / 2., DwWidth, DwHeight];
DwRefPt{1} = [PData(1).Origin(1), PData(1).Origin(2), PData(1).Origin(3)];
DwRefPt{2} = [PData(2).Origin(1), PData(2).Origin(2), PData(2).Origin(3)];

%Photoacoustic window
Resource.DisplayWindow = repmat(struct(...
    'Title', 'Photoacoustic', ...
    'ScrnSize', ScrnSize, ...
    'pdelta', DwPdelta, ...
    'Position', DwPos{1}, ...
    'ReferencePt', DwRefPt{1}, ...
    'numFrames', Resource.ImageBuffer(1).numFrames, ...
    'AxesUnits', 'mm', ...
    'Colormap', gray(256)), 1, 2);

%Ultrasound window
Resource.DisplayWindow(2).Position = DwPos{2};
Resource.DisplayWindow(2).ReferencePt = DwRefPt{2};
Resource.DisplayWindow(2).Title ='Ultrasound';
Resource.DisplayWindow(2).numFrames = Resource.ImageBuffer(2).numFrames;
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).splitPalette = 1;

%% TW : Specify transmit waveform structure -
TW(1).type ='parametric'; % simple generic pulse shape
TW(1).Parameters = [Trans.frequency, 0.67, 2, 1];
% - PA transmit waveform
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency, 0.67, 2, 1];

%% TX : Setting waveform and timing for each transmit channel -
% waveform - number of transmit waveform(TW structure ).
% Origin - effective origin of transmit beam at transducer
% focus - focal distance in wavelengths(for computing delay )
% Steer - [theta , alpha ]: Beam steering angles in radians.
% Apod - Weight of the channels from -1.0 to 1.0(0 for off ).
% Delay - Delay times of transmitter in wavelengths.Use computeTXDelays
% utility to calculate from set attributes(focus , Steer ,...).
% to the second index of the Trans.HVMux.Aperture array.
TX = repmat(struct('waveform', 1, ...
                     'Origin', [0.0 , 0.0 , 0.0] , ...
                      'focus', 0.0, ...
                      'Steer', [0.0 , 0.0] , ...              
                       'Apod', zeros(1, Resource.Parameters.numTransmit), ...
                      'Delay', zeros(1, Resource.Parameters.numTransmit)), ... % plane wave
                            1, 2); 

% Switch on transmit waves for ultrasound mode
TX(2).Apod = ones(1, Resource.Parameters.numTransmit);


%% Receive : Specify Receive structure arrays -
% Apod - setting receiving transducer channels active or inactive
% startDepth - starting depth of acquisition in wavelengths
% endDepth - end depth of acquisition in wavelengths
% TGC - which TGC waveform to use(1 for TGC(1))
% mode - 0= replace data , 1= accumulate
% bufnum - number of receive buffer to use
% framenum - number of frame in the receive buffer to which to write data
% acqNum - number of acquisition in frame sequence ;
% multiple acquisitions per frame possible.
% sampleMode - method of sampling for receive data stored in RcvBuffer.
%'NS200BW','NS200BWI','BS100BW','BS67BW','BS50BW',' custom'
% LowPassCoef - [1 x12] double sym.coefs for 23 tap FIR following A/D
% InputFilter - [1 x21] double sym.coefs of 42 tap input filter
% --------------------------------------------------------------------------
% Acquisition needs to be longer than the image because of signal paths ,
% which are inclined towards the normal of the transducer plane.
% Maximum length is the path from bottom right point to outer left
% transducer channel or vice versa.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements -1) * Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
            'startDepth', P.startDepth, ...
            'endDepth', maxAcqLength, ...
            'TGC', 1, ...
            'mode', 0, ...
            'bufnum', 1, ...
            'framenum', 1, ...
            'acqNum', 1, ...
            'sampleMode', 'NS200BW', ... % Nyquist 200 bandwidth
            'LowPassCoef', [], ...
            'InputFilter', []), ...
            1, 2 * Resource.RcvBuffer(1).numFrames + ...
            2 * Resource.RcvBuffer(2).numFrames);

% PA Receives
for i = 1:Resource.RcvBuffer(1).numFrames
    % left side receive events
    Receive(2*i-1).framenum = i;
    Receive(2*i-1).acqNum = 1;

    % right side receive events
    Receive(2*i).framenum = i;
    Receive(2*i).acqNum = 2;
end

% US Receives
nPaRcvs = 2*Resource.RcvBuffer(1).numFrames;
for i = 1:Resource.RcvBuffer(2).numFrames
    % TX Receive left
    Receive(2 * i -1 + nPaRcvs).bufnum = 2;
    Receive(2 * i -1 + nPaRcvs).TGC = 2;
    Receive(2 * i -1 + nPaRcvs).framenum = i;
    Receive(2 * i -1 + nPaRcvs).acqNum = 1;

    % TX Receive right
    Receive(2 * i -0 + nPaRcvs).bufnum = 2;
    Receive(2 * i -0 + nPaRcvs).TGC = 2;
    Receive(2 * i -0 + nPaRcvs).framenum = i;
    Receive(2 * i -0 + nPaRcvs).acqNum = 2;
end

%% TGC : Time Gain Compensation curve for receiving.-
% CntrlPts - [1x8] Control points for specifying curve.
% rangeMax - Maximum range in wavelengths
% Waveform - Gain values for TGC curve

% PA
TGC(1).CntrlPts = [50 141 275 404 510 603 702 782];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

% US
TGC(2).CntrlPts = [50 141 275 404 510 603 702 782];
TGC(2).rangeMax = P.endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

%% Recon : Specify Recon structure arrays.-
% senscutoff - Sensitivity threshold below which to exclude element from
% reconstruction.Between 0.0 and 1.0.
% pdatanum - Number of PData structure to use
% rcvBufFrame - Overrides the frame no.in ReconInfo(if given ).
% newFrameTimeout - Time(in ms) to wait for new frame.
% IntBufDest - [ InterBuffer number , frame number ]
% ImgBufDest - [ ImageBuffer number , frame number ]
% RINums - Row vector of ReconInfo structure numbers.
Recon = repmat(struct('senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'newFrameTimeout', 20000, ...
    'IntBufDest', [2, 1], ...
    'ImgBufDest', [2, -1], ...
    'RINums', 1:2 * Resource.RcvBuffer(2).numFrames), 1, 1);

%% ReconInfo : Unique ReconInfo structures referenced by Recon.
% Holds information on how to process a Region of pixel data.-
% mode - Mode of reconstruction , possible options :
%'replaceIntensity','addIntensity','multiplyIntensity',
%'replaceIQ','accumIQ','accumIQ_replaceIntensity',
%'accumIQ_addIntensity','accumIQ_multiplyIntensity'
% note: only the last three options write to the ImageBuffer
% txnum - Number of transmit(TX) object used for acquisition
% rcvnum - Number of Receive object used for acquisition
% regionnum - number of region to be reconstructed by ReconInfo object
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'scaleFactor', 1.0, ...
    'regionnum', 1), 1, 2 * Resource.RcvBuffer(2).numFrames);
% - US
for i = 1:Resource.RcvBuffer(2).numFrames
    % TX right Receive left
    ReconInfo(2 * i -1).txnum = 1;
    ReconInfo(2 * i -1).rcvnum = 2 * i -1 + nPaRcvs;

    % TX right Receive right
    ReconInfo(2 * i).txnum = 1;
    ReconInfo(2 * i).rcvnum = 2 * i + nPaRcvs;
end

ReconInfo(1).mode ='replaceIQ';
ReconInfo(end).mode ='accumIQ_replaceIntensity';

%% Process : Type of processing for acquired or reconstructed data -
% classname -'Image','Doppler' or'External'
% method -'procFunc' name of external processing function to call
% Parameters - Properties defining attributes of the processing object.
% Are assigned to processing object at start of VSX.
% E.g.: VSX creates Image object and sets these parameters.
Process(1).classname ='External';
Process(1).method ='recon_worker';

% Process.Parameters :
% srcbuffer - buffer type to read from: receive , image , imageP , inter
% srcbufnum - buffer number to read from
% srcframenum - starting frame nunmber
% dstbuffer - buffer type to write to: receive , image , imageP , inter
% dstbufnum - buffer number to write to
% dstframenum - destination frame number in the buffer
Process(1).Parameters = {'srcbuffer', 'receive', ... % read from Rcv Buffer
                        'srcbufnum', 1, ...
                        'srcframenum', -1, ... % most recent frame
                        'dstbuffer', 'image', ... % write to ImageBuffer
                        'dstbufnum', 1, ...
                        'dstframenum', -2}; % last frame , increment on completion

Process(2).classname ='External';
Process(2).method ='updatePaDisplayWindow';
Process(2).Parameters = {'srcbuffer', 'image', ...
                        'srcbufnum', 1, ...
                        'srcframenum', -1, ...
                        'dstbuffer', 'imageP', ...
                        'dstbufnum', 1, ...
                        'dstframenum', -1};

% - US
Process(3).classname ='Image';
Process(3).method ='imageDisplay';
Process(3).Parameters = {'imgbufnum', 2, ... % number of buffer to process.
                        'framenum', -1, ... %(-1 => most recent )
                        'pdatanum', 2, ... % number of PData structure
                        'pgain', 1.0, ... % image processing gain
                        'reject', 3, ... % reject level
                        'persistMethod', 'simple', ...
                        'persistLevel', 20, ...
                        'interpMethod', '4pt', ...
                        'grainRemoval', 'none', ...
                        'processMethod', 'none', ...
                        'averageMethod', 'none', ...
                        'compressMethod', 'power', ...
                        'compressFactor', 40, ...
                        'mappingMethod', 'full', ...
                        'display', 1, ... % display image
                        'displayWindow', 2};

%% SeqControl : Sequence control event called from Event.seqControl
SeqControl(1).command   ='triggerIn';
SeqControl(1).condition ='Trigger_1_Falling'; % FIXED SYNC OUT
SeqControl(2).command   ='returnToMatlab';
SeqControl(3).command   ='jump';
SeqControl(3).argument  = 1;
SeqControl(4).command   ='sync';
SeqControl(4).argument  = 1e5; % in us
nsc = size(SeqControl, 2)+1;

%% Event : Specify sequence events.
% info - descriptive text( optional )
% tx - number of TX structure array to use
% rcv - number of Receive structure array to use
% recon - number of Recon structure array for reconstruction
% process - number of Process structure array to use
% seqControl - number of Sequence Control structure array to use
n = 1;

for i = 1:Resource.RcvBuffer(1).numFrames
    % Photoacoustic frame
    Event(n).info ='PA imaging for left half';
    Event(n).tx = 1;
    Event(n).rcv = 2 * i -1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 1; % Wait for trigger
    n = n + 1;

    Event(n).info ='PA imaging for right half';
    Event(n).tx = 1;
    Event(n).rcv = 2 * i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 1; % Wait for trigger
    n = n + 1;

    Event(n).info ='Transfer PA frame to host.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command ='transferToHost'; % transfer frame to host buffer
    nsc = nsc +1;
    n = n + 1;

    Event(n).info ='Reconstruct and process PA frame';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n + 1;

    Event(n).info ='Update Display Window';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n + 1;

    %% Ultrasound frame
    Event(n).info ='US imaging for TX right and Receive left';
    Event(n).tx = 2;
    Event(n).rcv = 2 * i -1 + nPaRcvs;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n + 1;

    Event(n).info ='US imaging for TX right and Receive right';
    Event(n).tx = 2;
    Event(n).rcv = 2 * i + nPaRcvs;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n + 1;

    Event(n).info ='Transfer US frame to host.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command ='transferToHost'; % transfer frame to host buffer
    nsc = nsc +1;
    n = n + 1;

    Event(n).info ='Reconstruct and process US frame';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n + 1;

    Event(n).info ='Synchronize hardware and software sequencers';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;

    if (mod(i, 5) == 0) && (i ~= Resource.RcvBuffer(1).numFrames)
        Event(n).seqControl = 2;
    end

    n = n + 1;
end

% Jump back to start Event after capturing on each frame of the RcvBuffer
Event(n).info ='Jump back to start event.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;

%% Define in -file external function callbacks

% Call to reconstruction method
EF(1).Function = text2cell('%recon_worker');
EF(2).Function = text2cell('%updatePaDisplayWindow');

% UI control for changing reconstruction method
UIButtonLabels{1} = {'SAFT', 'FSAFT', 'Fourier'};
UI(1).Control = {'UserB4', ...
                'Style', 'VsButtonGroup', ...
                'Title', 'Recon Method', ...
                'NumButtons', 3, ...
                'Labels', UIButtonLabels{1}};
UI(1).Callback = text2cell('%chooseRecon');

% UI control for PA overlay
UI(2).Control = {'UserC2', ...
                'Style', 'VsToggleButton', ...
                'Label', 'PAT Overlay'};
UI(2).Callback = text2cell('%togglePatOverlay');

% UI control for changing the number of acquisitions to average over
UI(3).Control = {'UserB2', ...
                'Style', 'VsSlider', ...
                'Label', '# Averages', ...
                'SliderMinMaxVal', [1, 10, 1], ...
                'SliderStep', [0.1, 1], ...
                'ValueFormat', '%2.0i'};
UI(3).Callback = text2cell('%numberOfAverages');

%% Save all the structures to a.mat file
save(filename);
return % Prevent Matlab from executing any code beyond this point

%% Functions to be read by text2cell

%chooseRecon
    switch UIState
        case 1
            reconMethod ='saft';
        case 2
            reconMethod ='fsaft';
        case 3
            reconMethod ='fourier';
        otherwise
            error('VSX:UI:reconMethod', 'Unrecognized reconstruction method !');
    end
    
    assignin('base', 'temp', reconMethod);
    evalin('base', 'P.paReconMethod=temp;');
        
    Control = evalin('base', 'Control');
    Control.Command ='update &Run';
    Control.Parameters = {'Event', 'SeqControl', 'Process'};
    assignin('base', 'Control', Control);
%chooseRecon

%numberOfAverages
    nAverages = round(UIValue);
    assignin('base', 'temp', nAverages);
    evalin('base', 'P.nAvg = temp;');
    Control = evalin('base', 'Control');
    Control.Command ='update&Run';
    Control.Parameters = {'Event', 'Process'};
    assignin('base', 'Control', Control);
%numberOfAverages

%togglePatOverlay
    global overlayImgHandle
    
    if UIState ~= 0
        
        if isempty(overlayImgHandle) || ~ishandle(overlayImgHandle)
            displayWindowHandle = evalin('base', 'Resource.DisplayWindow(1).imageHandle');
            axesHandle = findobj(displayWindowHandle, 'type', 'axes');
            h = figure;
            overlayImgHandle = copyobj(axesHandle, h);
        end
        
    else
        
        if ~isempty(overlayImgHandle) && ishandle(overlayImgHandle)
            close(overlayImgHandle);
        end
        
    end
        
    assignin('base', 'temp', UIState);
    evalin('base', 'P.paOverlayToggled = temp;');
    assignin('base', 'temp', overlayImgHandle);
    evalin('base', 'P.overlayImgHandle = temp;');
%togglePatOverlay

%updatePaDisplayWindow
    imgPData = updatePaDisplayWindow(ImgData)
    global paImgHandle usImgHandle overlayImgHandle
    
    if isempty(imageHandle) || ~ishandle(imageHandle)
        paImgHandle = evalin('base', ...
                             'Resource.DisplayWindow(1).imageHandle');
    end
    % Normalize to be between 0 and 255
    imgPData = scale(imgData, [0, 255]); % mex -file
        
    if ishandle(paImgHandle)% Check if window has been closed by the user
        set(paImgHandle, 'CData', imgPData);
    end
        
    % PA overlay
    if evalin('base', 'P.paOverlayToggled')
        overlayImgHandle = evalin('base', 'P.overlayImgHandle');
        if ~isempty(h) && ishandle(h)
            usImgHandle = evalin('base', ...
                                 'Resource.DisplayWindow(2).imageHandle');
        end
    end
%updatePaDisplayWindow

%recon_worker
    ImgData = recon_worker(RData)
    global loopCnt RData_buf prevRData
    
    Receive = evalin('base', 'Receive');
    Trans = evalin('base', 'Trans');
    PData = evalin('base', 'PData');
    nAvg = evalin('base', 'P.nAvg');
    Parameters = evalin('base', 'Resource.Parameters');
    reconMethod = evalin('base', 'P.paReconMethod');
    if (isempty(RData_buf)||isempty(prevRData)|| ...
           size(RData_buf,1) ~= nAvg || isempty(loopCnt))
       RData_buf = cell(nAvg, 1);
       prevRData = [];
       loopCnt = 1;
    end
    
    temp = [RData(Receive(1).startSample:Receive(1).endSample, :) ...
            RData(Receive(2).startSample:Receive(2).endSample, :) ];
    
    % Check if reprocessing previous frame
    if ~isequal(prevRData, temp)
        RData_buf{loopCnt} = temp;
        prevRData = temp;
        loopCnt = loopCnt + 1;
    end
    
    % Average RF data with previous acquisitions
    RData_avg = zeros(size(RData_buf{1}));
    avgCnt = 0;
    
    for i = 1:nAvg
        if ~isempty(RData_buf{i})
            RData_avg = RData_avg + RData_buf{i};
            avgCnt = avgCnt + 1;
        end
    
    end
    
    if avgCnt ~= 0
        RData_avg = RData_avg/avgCnt;
    end
    
    switch reconMethod
        case 'saft'
            ImgData = saft2d(RData_avg, Receive(1), Trans, PData(1), Parameters);
        case 'fsaft'
            ImgData = fsaft2d_pComp_mex(RData_avg, Receive(1), Trans, PData(1), Parameters);
        case 'fourier'
            ImgData = fdrecon2d(RData_avg, Receive(1), Trans, PData(1), Parameters);
        otherwise
            error('VSX:Process:reconMethod', [reconMethod ...
            'is not recognized as a reconstruction method !']);
    end
    
    % wrap around loop counter
    if (loopCnt > nAvg)
        loopCnt = 1;
    end
%recon_worker
