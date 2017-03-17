function catcherror = FallingFruit

catcherror = 0; % default to no error code, if error, this variable will contain details, and a stack trace
[subj, runnum, sex] = getSessionInfo; if isempty(subj); return; end

% ------------------------- user defined vars -----------------------------
% -------------------------------------------------------------------------

keys = {'escape', 'q'};
if isnumeric(subj); subj = num2str(subj); end;
if isnumeric(runnum); runnum = num2str(runnum); end
s.subj = subj;
s.runnum = runnum;
s.sex = sex;
rootdir = fileparts(mfilename('fullpath')); %set up directories
datadir = fullfile(rootdir,'data');
picdir = fullfile(rootdir,'cropped-big');
s.picdir = picdir;
s = loadImgsFromDisk(s);
s.stimList = readtable(fullfile(rootdir,'stims.csv'));
datestring = getDateAndTime; % get date string down to the minute
s.datafile = fullfile(datadir,sprintf('%s_%s_%s_%s.csv',subj, runnum,mfilename, datestring)); % make data file name for saving
veryPreciseTimingNeeded = false;
if ~exist(datadir,'dir') % if data folder doesn't exist, make it
    mkdir(datadir);
end
sca; % clear the screen, and other PTB backend things
maxNumItems = 8;
s.maxTime = 120; % seconds
s.imin = 5; % minimum rate of fall
s.imax = 8; % maximum rate of fall 
s.estTableSize = round(s.maxTime+((s.maxTime*0.25)*mean([s.imin s.imax]))); % preallocate for speed! Zoom zoom....



% ------ don't edit below this line unless you know what you're doing -----
% -------------------------------------------------------------------------
instruct = sprintf(['In this game, you will see different fruits and\n',...
                    'objects falling from the top of the screen.\n',...
                    'You are to use the mouse to move the basket\n',...
                    'around the screen and catch as many apples as possible.\n',...
                    'Do your best to avoid the other objects as you try to catch the apples.\n\n',...
                    '[Press the spacebar to begin]']);
s.T = table; % init empty table for saving
s.tbidx = 0;
s.timeRemaining = s.maxTime;
KbName('UnifyKeyNames'); % make all keyboards similar
Screen('Preference', 'VisualDebugLevel', 1);
if veryPreciseTimingNeeded
    Screen('Preference', 'SkipSyncTests', 0); %#ok
else
    Screen('Preference', 'SkipSyncTests',2);
end
PsychDefaultSetupPlus(3); % set up PTB with some standard values
params = PsychSetupParams(0,1); % set some standard parameters for this session
params.titleBarRect = [0 0 params.maxXpixels 80];
HideCursor;
params.pointerColor = [0 0 0];
SetMouse(params.Xc, params.Yc, params.win);
waitframes = 1;
numCornersInRect = 4;
s.maxNumItems = maxNumItems;
s.words = s.stimList.word(1:s.maxNumItems)';
s.wi = 1:maxNumItems;
s.colors = repmat([254 95 85]/255,s.maxNumItems,1);
s.isHit(1:maxNumItems) = 0;
s.cantHitYet(1:maxNumItems) = 0;
s.newPos(1:maxNumItems) = 1;
s.numHits = 0;
s.flipText(1:maxNumItems) = 0;
ratesVec = setRatesOfFall(s.imin, s.imax, maxNumItems); 
ratesVec
s.ratesVec = ratesVec;
[s.xGridLocs, s.yGridLocs] = createScreenGrid(params.maxXpixels, params.maxYpixels, maxNumItems+1);
s.xGridLocs = s.xGridLocs(1,2:end);
s.yGridLocs = 1:params.maxYpixels; % override screen grid output
s.words = s.stimList.word(1:s.maxNumItems)';
s.wi = 1:maxNumItems;
s = balanceXaxisLocations(s);
s.xlocs = s.xGridLocs(s.stimList.xlocs(1:s.maxNumItems)');
s.picSizeXDirection = s.xlocs(1)*0.3; % make pics take up 30% of the column they are in
s.basketSize = s.picSizeXDirection*2;
s = resizeImgsForThisScreen(s);
uipointerImg = imread(fullfile(s.picdir, 'basket.jpg'));
uipointerImgWidthToHeighRatio = size(uipointerImg,1)/size(uipointerImg,2);
uipointerImg = imresize(uipointerImg,[s.basketSize, s.basketSize*uipointerImgWidthToHeighRatio]);
uipointer = [0 0 size(uipointerImg,1) size(uipointerImg,2)];
uipointerOffset = uipointer(end)/2;
RestrictKeysForKbCheck(cellfun(@KbName, keys));
s.textBounds = zeros(maxNumItems, numCornersInRect);
ShowInstructions(params, s, instruct, {'space', 'escape'});
WaitSecs(0.5); % wait for 500 ms just to have a smooth transition from instruct to task
% Sync to get a time stamp
vbl = Screen('Flip', params.win);
s.vbl = vbl;
s.startTime = vbl;

while ~KbCheck & s.timeRemaining > 1 %#ok
    [mx, my, s.uipointer] = updateUIPointer(params, uipointer, uipointerImg);
    s = showStim(params, s, mx, my-uipointerOffset);
    s = showSessionStatsInTitleBar(params,s);
    vbl  = Screen('Flip', params.win, vbl + (waitframes - 0.5) * params.ifi); % Flip to the screen
    s.vbl = vbl;
end
filledCells = ~cellfun(@isempty,s.T.subj); % only use table rows that have subj values (gets rid of preallocated shit at the end)
onlyHits = find(s.T.hit(filledCells) > 0); % get the hit data using all the filled cells
onlyTargs = find(s.T.isTarget(filledCells) > 0); % get the target data using all the filled cells
% targsHit = find(onlyHits > 0 & onlyTargs > 0);
% alltargs = find(onlyTargs > 0);
% acc = length(targsHit)/length(onlyHits); % get average for the newly created vector with hit and target data
writetable(s.T,s.datafile);
CleanUp;
WaitSecs(1);
syncFilesToCloud({s.datafile}, subj, runnum);
ShowCursor;
% fprintf('\n\n***************\n\nAccuracy: %2.2f%%\nPoints: %d\n***************\n\n',acc*100,s.numHits);
% uiwait(msgbox(sprintf('Accuracy: %2.2f%%\nPoints: %d',acc*100,s.numHits)));
end























function PsychDefaultSetupPlus(featureLevel)
% PsychDefaultSetup(featureLevel) - Perform standard setup for Psychtoolbox.

% Default colormode to use: 0 = clamped, 0-255 range. 1 = unclamped 0-1 range.
global psych_default_colormode;
psych_default_colormode = 0;

% Reset KbName mappings:
clear KbName;

% Define maximum supported featureLevel for this Psychtoolbox installation:
maxFeatureLevel = 3;

% Sanity check featureLevel argument:
if nargin < 1 || isempty(featureLevel) || ~isscalar(featureLevel) || ~isnumeric(featureLevel) || featureLevel < 0
    error('Mandatory featureLevel argument missing or invalid (not a scalar number or negative).');
end

% Always AssertOpenGL:
AssertOpenGL;

% Level 1+ requested?
if featureLevel >= 1
    % Unify keycode to keyname mapping across operating systems:
    KbName('UnifyKeyNames');
end

% Level 2+ requested?
if featureLevel >= 2
    % Initial call to timing functions
    % Set global environment variable to ask PsychImaging() to enable
    % normalized color range for all drawing commands and Screen('MakeTexture'):
    psych_default_colormode = 1;
    GetSecs; WaitSecs(0.001);
end

% Level 2+ requested?
if featureLevel >= 3
    %suppress keypress to command window,
    %and hide the mouse pointer (usefull is most visual experiments)
    ListenChar(2);
    HideCursor;
end


if featureLevel > maxFeatureLevel
    error('This installation of Psychtoolbox can not execute scripts at the requested featureLevel of %i, but only up to level %i ! UpdatePsychtoolbox!', featureLevel, maxFeatureLevel);
end
return;
end


function params = PsychSetupParams(doAlphaBlending,doMultiSample)
%sets up some normal values used in experiments such as a gray background
%and Arial font, and a large text size, etc...
%saves all relevant screen info to the 'params' structure so that the
%entire structure can be passed in and out of functions, rather than
%zillions of variables. Also makes it expandable.
%
% History:
% 29-May-2015   th     made initial version of the function

global psych_default_colormode;
%make params structure
params = struct;
%set some defualt, common colors
params.colors.white = [1 1 1];
params.colors.black = [0 0 0];
params.colors.gray = [0.5 0.5 0.5];
params.colors.red = [1 0 0];
params.colors.green = [0 1 0];
params.headerColor = [184 216 216]/255;
%check if using normalized color values or not
if psych_default_colormode == 0
    params.colors.white = [255 255 255];
    params.colors.gray = [128 128 128];
end
%choose max screen number (will be the external monitor if connected)
params.screen = max(Screen('Screens'));
params.font = 'Arial'; %set the global font for PTB to use
params.tsize = 28; %set text size
params.TextColor = [79 99 103]/255; %set global text color
%set the background color of the screen (defaults to gray)
%params.background = [238 245 219]/255;
params.background = [1 1 1];
params.multiSample = [];
if doMultiSample
    params.multiSample = 4;%set to a value greater than 0 if you want super sampling
end
%open the PTB window
[params.win, params.rect] = PsychImaging('OpenWindow', params.screen, params.background,[],[],[],[],params.multiSample);
%get screen width and height
[params.maxXpixels, params.maxYpixels] = Screen('WindowSize', params.win);
if doAlphaBlending
    %Set blend function for alpha blending
    Screen('BlendFunction', params.win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
end
%find center of screen
[params.Xc,params.Yc] = RectCenter([0 0 params.maxXpixels params.maxYpixels]);
%now that the window pointer exists, set some values from earlier
Screen('TextSize', params.win, params.tsize);
Screen('TextFont',params.win, params.font);
Screen('TextSize',params.win, params.tsize);
Screen('TextStyle', params.win, 1);

%Maximum priority level
params.topPriorityLevel = MaxPriority(params.win);
Priority(params.topPriorityLevel);
%Query the frame duration
params.ifi = Screen('GetFlipInterval', params.win);
end


function ratesVec = setRatesOfFall(imin, imax, maxNumItems)
ratesVec = round(rand(1,maxNumItems)+randi([imin imax],1,maxNumItems));
end


function s = updateRateOfFall(s)
r = s.ratesVec;
rv = rand(1,1)+randi([s.imin s.imax],1,1);
r(s.i) = rv;
s.ratesVec = round(r);
end


function [xGridLocs, yGridLocs] = createScreenGrid(maxX, maxY, maxNumItems)
[xGridLocs,yGridLocs] = meshgrid(linspace(0,maxX,maxNumItems+1),linspace(0,maxY,maxY));
xGridLocs = xGridLocs(:,1:end-1); % remove the last column because it starts at the max x position (we wouldn't see it after it is drawn on the screen)
yGridLocs = yGridLocs(:,1:end-1);
end


function [mx, my, uipointer] = updateUIPointer(params, uipointer, uipointerImg, offset)
if nargin < 4
    offset = 0;
end
[mx, my] = GetMouse(params.win); % Get the current position of the mouse
%my = params.maxYpixels-offset;
my = my-offset;
%my = my-uipointer(end); % use top edge of uipointer for hitting targets instead of center
uipointer = CenterRectOnPointd(uipointer, mx, my); % Center the rectangle on the mouse
%Screen('FillRect', params.win, params.pointerColor, uipointer); % Draw the rect to the screen
uiTexture=Screen('MakeTexture', params.win, uipointerImg);
Screen('DrawTexture', params.win, uiTexture, [], uipointer)
Screen('Close', uiTexture);
end


function s = showStim(params, s, mx, my)
s.mx = mx;
s.my = my;
if ~isfield(s,'lastYpos')
    s.lastYpos(1:s.maxNumItems) = 1;
end
for i = 1:s.maxNumItems
    s.i = i;
    s = checkIfUIPointerIsInObject(mx,my,s);
    s = reverseDirectionForHitObjects(s);
    s = drawItemsToScreen(params, s);
    s = checkIfObjectWentOffScreen(params, s);
end
end


function datestring = getDateAndTime
d = fix(clock);
datestring=sprintf('Y%04d_M%02d_D%02d_H%02d_M%02d_S%02d',...
    d(1),...
    d(2),...
    d(3),...
    d(4),...
    d(5),...
    d(6));
end


function CleanUp(files)
if nargin < 1
    files = [];
end
ShowCursor;
ListenChar(0);
sca;
RestrictKeysForKbCheck([]);
end %CleanUp



function [subj, runnum, sex] = getSessionInfo
subj = [];
runnum = [];
sex = [];
prompt={'Participant: ','Session: ', 'Sex(M, F, NA): '};
   name='HAALT Sem';
   numlines=1;
   defaultanswer={'0','0', ''};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(answer); return; end
subj = answer{1};
runnum = answer{2};
sex = answer{3};
end


function s = showSessionStatsInTitleBar(params,s)
s.timeRemaining = s.maxTime-(s.vbl-s.startTime);
Screen('FillRect', params.win, params.headerColor, params.titleBarRect);
x = Screen('DrawText', params.win, sprintf('Time remaining: %3.0f ', s.timeRemaining),5,30, params.TextColor);
Screen('DrawText', params.win, sprintf('Points: %03d ',s.numHits),x+10,30, params.TextColor);
end



function s = checkIfObjectWentOffScreen(params,s)
% check if object went below last pixel in y direction or higher than first pixel (negative)
if (s.lastYpos(s.i) > params.maxYpixels) | s.lastYpos(s.i) <= 1 %#ok (for not using ||)
    %save data every time an object leaves the screen (includes hits and misses)
    if ~s.isHit(s.i)
        s.tbidx = s.tbidx+1;
        s = updateDataTable(s);
    end
    s.lastYpos(s.i) = 1;
    s = updateRateOfFall(s);
    s = updateWords(s);
    s.isHit(s.i) = 0;
    s.cantHitYet(s.i) = 0;
end
end


function s = checkIfUIPointerIsInObject(x,y,s)
pointerTopLeftX = s.uipointer(RectLeft);
pointerTopLeftY = s.uipointer(RectTop);
pointerTopRightX = s.uipointer(RectRight);
pointerTopRightY = s.uipointer(RectTop);
% if s.flipText(s.i)
%     rB = s.textBounds(s.i,RectBottom);
%     rT = s.textBounds(s.i,RectTop);
%     s.textBounds(s.i,RectBottom) = s.textBounds(s.i,RectBottom) + (rB-rT) + 15;
%     s.textBounds(s.i,RectTop) = s.textBounds(s.i,RectTop) + (rB-rT);
% end
if IsInRect(x,y,s.textBounds(s.i,:)) | IsInRect(pointerTopLeftX,pointerTopLeftY,s.textBounds(s.i,:)) | IsInRect(pointerTopRightX,pointerTopRightY,s.textBounds(s.i,:)) %#ok
    if ~s.cantHitYet(s.i)
        s.isHit(s.i) = 1;
        if s.stimList.isTarget(s.wi(s.i))
            s.numHits = s.numHits+1;
        end
        s.cantHitYet(s.i) = 1;
        s.tbidx = s.tbidx+1;
        s = updateDataTable(s);
    end
else
    %s.isHit(i) = 0;
end
end


function s = reverseDirectionForHitObjects(s)
if s.isHit(s.i) == 1
    hitAwayMultiplier = 10;
    %s.newPos(s.i) = round(s.yGridLocs(s.lastYpos(s.i))-s.ratesVec(s.i)-hitAwayMultiplier);
    s.newPos(s.i) = -1;
    if s.newPos(s.i) < 0;
        s.newPos(s.i) = 1;
        s.cantHitYet(s.i) = 1;
        %s.isHit(s.i) = 0;
        s = updateRateOfFall(s);
    end
else
    s.newPos(s.i) = round(s.yGridLocs(round(s.lastYpos(s.i)))+s.ratesVec(s.i));
end
end

function s = drawItemsToScreen(params,s)
% if ~isfield(s,'newPos')
%     s.newPos(s.i) = round(s.yGridLocs(s.lastYpos(s.i))+s.ratesVec(s.i));
% end
% if s.stimList.isTarget(s.wi(s.i));
%     flipText = 1;
%     revText = 0;
% else
%     flipText = 0;
%     revText = 0;
% end
%s.flipText(s.i) = flipText;
% [~,s.lastYpos(s.i),s.textBounds(s.i,:)] = DrawFormattedText(params.win, s.words{s.i}, s.xlocs(1,s.i),...
%     s.newPos(s.i), s.colors(s.i,:), [], revText, flipText);
%[nx, ny, textbounds] = DrawFormattedText(win, tstring [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical][, vSpacing][, righttoleft][, winRect])


%[wordRects, s.textBounds(s.i,:)] = makeRectsFromWord(s.words{s.i}, s.xlocs(1,s.i), s.newPos(s.i));
%s.lastYpos(s.i) = s.textBounds(s.i,RectTop);

img = s.(s.words{s.i});
%img = s.('apple');
imgRect = CenterRectOnPoint([0 0 size(img,1) size(img,2)], s.xlocs(1,s.i), s.newPos(s.i)+(size(img,2)/2));
s.textBounds(s.i,:) = imgRect;
s.lastYpos(s.i) = imgRect(RectTop);
thisTexture=Screen('MakeTexture', params.win, img);
Screen('DrawTexture', params.win, thisTexture, [], imgRect);
Screen('Close', thisTexture);
%Screen('FrameRect', params.win, s.colors(s.i,:), wordRects, penWidthPixels);
if s.lastYpos(s.i) <= 0; s.lastYpos(s.i) = 1; end;
end

function s = updateDataTable(s)
%vars = {'word', 'hit', 'numLetters', 'numVows', 'numCons', 'rateOfFall', 'numSyl', 'isTarget', 'scrXpos', 'scrYpos', 'dxFromMouse'};
subj = deblank(s.subj);
runnum = deblank(s.runnum);
sex = deblank(s.sex);
word = deblank(s.words{s.i});
hit = s.isHit(s.i);
rateOfFall = s.ratesVec(s.i);
scrXpos = s.xlocs(s.i);
scrYpos = s.lastYpos(s.i); % of word
dxFromMouse = getDistanceFromMouse(s);
mouseX = s.mx;
mouseY = s.my;
%isTarget = s.stimList.isTarget(s.wi(s.i));
if strcmpi(word,'apple')
    isTarget = 1;
else
    isTarget = 0;
end
T = table({subj}, {runnum}, {sex}, {word}, hit, rateOfFall, isTarget, scrXpos, scrYpos, mouseX, mouseY, dxFromMouse);
T.Properties.VariableNames = {'subj', 'runNum', 'sex', 'word', 'hit', 'rateOfFall', 'isTarget', 'scrXpos', 'scrYpos', 'mouseX', 'mouseY', 'dxFromMouse'};
if size(s.T,1) < 1
    s.T = T;
    s.T(s.estTableSize,:) = T(1,:);
    s.T(s.estTableSize,:) = [];
else
    s.T(s.tbidx,:) = T(1,:);
end
end


function n = numberOfVowelsInWord(w)
v = ('aeiou');
n = sum(ismember(w,v));
end


function n = numberOfConsInWord(w)
c = 'a':'z';
v = 'aeiou';
m = ismember(c,v);
c = c(~m);
n = sum(ismember(w,c));
end


function val = isTallLetter(letter)
tl = ('bdfhklt');
if ismember(letter,tl)
    val = 1;
else
    val = 0;
end
end


function val = isShortLetter(letter)
sl = ('aceimnorsuvwxz');
if ismember(letter,sl)
    val = 1;
else
    val = 0;
end
end


function val = isHangingLetter(letter)
hl = ('gjpqy');
if ismember(letter,hl)
    val = 1;
else
    val = 0;
end
end


function c = getDistanceFromMouse(s)
%when onbject enters viewing area
[x,~] = RectCenter(s.textBounds(s.i,:));
a = abs(x-s.mx);
b = abs(s.my-1);
c = sqrt(a^2 + b^2);
end


function s = updateWords(s)
s.wi(s.i) = s.wi(s.i)+s.maxNumItems;
if s.wi(s.i) > size(s.stimList,1);
    s.wi(s.i) = s.i;
end
s.words{1, s.i} = deblank(char(s.stimList.word(s.wi(s.i))));
s.xlocs(1,s.i) = s.xGridLocs(s.stimList.xlocs(s.wi(s.i)));
s.onsetTime(s.i) = s.vbl;
end


function syncFilesToCloud(files, subj, runnum)
if isempty(files); return; end;
cloudDir = ('~/Box Sync/MUSC_POLAR'); % check for box (MUSC computers)
if ~isdir(cloudDir)
    cloudDir = ('~/Dropbox (C-STAR)'); % if that doesn't exist check this
    if ~isdir(cloudDir) 
        cloudDir = ('~/Dropbox'); % check this last
        if ~isdir(cloudDir)
            warning('Data syncing not available. No Drop(box) folder detected at %s', cloudDir);
            return;
        end
    end
end;
taskFolder = fullfile(cloudDir,'PolarData', mfilename, subj, runnum);
if ~isdir(taskFolder)
    mkdir(taskFolder);
end
n = size(files,1);
h = waitbar(0,'Copying files for data syncing...');
steps = n;
for i = 1:steps
    [~,nm,ext] = fileparts(files{i});
    copyfile(files{i},fullfile(taskFolder,[nm ext]));
    waitbar(i / steps)
end
close(h) 
end


function startTime = ShowInstructions(params, s, instruct,keysToWaitFor)
if nargin < 3; keysToWaitFor = {'space', 'escape'}; end;
Screen('Flip',params.win); % for windows
DrawFormattedText(params.win, instruct, 'center', 'center',params.TextColor);
Screen('Flip',params.win);
RestrictKeysForKbCheck(cellfun(@KbName, keysToWaitFor));
deviceN = -1;
[startTime, keyCode] = KbWait(deviceN);
if strcmpi(KbName(keyCode), 'escape')
    CleanUp;
end
RestrictKeysForKbCheck([]);
Screen('Flip',params.win);
end


function s = balanceXaxisLocations(s)
sz = size(s.stimList,1); % query number of rows
nTargets = sum(s.stimList.isTarget); % get number of targets
nDistractors = sz - nTargets; % get number of distractors
rmT = rem(nTargets,s.maxNumItems); % see if nTargets is divisible by number of items on the screen
rmD = rem(nDistractors,s.maxNumItems); % see if nDistractors is divisible by number of items on the screen
nTargets = nTargets - rmT; % if it wasn't divisible, make it so
nDistractors = nDistractors - rmD; % if it wasn't divisible, make it so
tLocs = BalanceTrials(nTargets,1,1:s.maxNumItems); % balance target locations for number of items on screen
dLocs = BalanceTrials(nDistractors,1,1:s.maxNumItems); % balance distractor locations for number of items on screen
tTar = s.stimList;
tTar(~tTar.isTarget,:) = [];
tDis = s.stimList;
tDis(tDis.isTarget>0,:) = [];
tTar = tTar(1:nTargets,:);
tDis = tDis(1:nDistractors,:);
tTar.xlocs = tLocs; %locations produce overlapping stimuli, don't use for now
tDis.xlocs = dLocs;
s.tTar = tTar;
s.tDis = tDis;
s.stimList = [s.tTar; s.tDis];

%%% randomly permute the stimulus list (happens at the start of each
%%% session)
sz = size(s.stimList,1); % query number of rows
r = randperm(sz); % make a random permutation array based on number of rows
cols = s.stimList.Properties.VariableNames; % get the table's column names
ncols = numel(cols); % get number of columns
for i = 1:ncols
    s.stimList.(cols{i}) = s.stimList.(cols{i})(r); % give each column the same random permutation
end

%%% get indices of each stimuli's column it will fall from
idx1 = find(s.stimList.xlocs == 1);
idx2 = find(s.stimList.xlocs == 2);
idx3 = find(s.stimList.xlocs == 3);
idx4 = find(s.stimList.xlocs == 4);
idx5 = find(s.stimList.xlocs == 5);
idx6 = find(s.stimList.xlocs == 6);
idx7 = find(s.stimList.xlocs == 7);
idx8 = find(s.stimList.xlocs == 8);
newxlocs = [];
for i = 1:length(1:s.maxNumItems:sz)
    newxlocs = [newxlocs; idx1(i); idx2(i); idx3(i); idx4(i); idx5(i); idx6(i); idx7(i); idx8(i);];
end

%%% reorder stimulus list by screen location (1,2,3,4 -- 1,2,3,4 --
%%% 1,2,3,4) will be repeating like that.
cols = s.stimList.Properties.VariableNames; % get the table's column names
ncols = numel(cols); % get number of columns
for i = 1:ncols
    s.stimList.(cols{i}) = s.stimList.(cols{i})(newxlocs); % give each column the same random permutation
end
end


function [wordRects, textBounds] = makeRectsFromWord(word, beginX, beginY)
nLetters = numel(word);
wordRects = zeros(4,nLetters);
rectW = 12;
smallRectH = 12;
bigRectH = 24;
smallRect = [0 0 rectW smallRectH];
bigRect = [0 0 rectW bigRectH];
for i = 1:nLetters
    if isTallLetter(word(i))
        thisRect = bigRect;
    elseif isShortLetter(word(i))
        thisRect = smallRect;
    elseif isHangingLetter(word(i))
        thisRect = bigRect;
    else
        error('%s is not a valid character for this task',word(i));
    end
    if i > 1
        if isTallLetter(word(i-1)) & isShortLetter(word(i)) %#ok %tall letter then short letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-smallRectH);
        elseif isShortLetter(word(i-1)) & isTallLetter(word(i)) %#ok % short letter then tall letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-bigRectH);
        elseif isShortLetter(word(i-1)) & isHangingLetter(word(i)) %#ok short letter then hanging letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-smallRectH);
        elseif isTallLetter(word(i-1)) & isHangingLetter(word(i)) %#ok tall letter then hanging letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-smallRectH);
        elseif isHangingLetter(word(i-1)) & isShortLetter(word(i)) %#ok hanging letter then short letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-bigRectH);
        elseif isHangingLetter(word(i-1)) & isTallLetter(word(i)) %#ok hanging letter then tall letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-bigRectH-smallRectH);
        elseif isTallLetter(word(i-1)) & isTallLetter(word(i)) %#ok tall letter then tall letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-bigRectH);
        elseif isShortLetter(word(i-1)) & isShortLetter(word(i)) %#ok short letter then short letter
            thisRect = OffsetRect(thisRect,wordRects(3,i-1),wordRects(4,i-1)-smallRectH);
        end
    else
        thisRect = OffsetRect(thisRect, beginX, beginY);
    end
    wordRects(:,i) = thisRect;
end
textBounds = [min(wordRects(RectLeft,:)), min(wordRects(RectTop,:)), max(wordRects(RectRight,:)), min(wordRects(RectBottom,:))];
end


function s = loadImgsFromDisk(s)
s.apple = imread(fullfile(s.picdir, 'apple.jpg'));
s.orange = imread(fullfile(s.picdir, 'orange.jpg'));
s.peach = imread(fullfile(s.picdir, 'peach.png'));
s.basketball = imread(fullfile(s.picdir, 'basketball.png'));
s.button = imread(fullfile(s.picdir, 'button.jpg'));
s.wheel = imread(fullfile(s.picdir, 'wheel.jpg'));
end

function s = resizeImgsForThisScreen(s)
appleWtoHRatio = size(s.apple,1)/size(s.apple,2);
orangeWtoHRatio = size(s.orange,1)/size(s.orange,2);
peachWtoHRatio = size(s.peach,1)/size(s.peach,2);
basketballWtoHRatio = size(s.basketball,1)/size(s.basketball,2);
buttonWtoHRatio = size(s.button,1)/size(s.button,2);
wheelWtoHRatio = size(s.wheel,1)/size(s.wheel,2);

newAppleSize = [s.picSizeXDirection s.picSizeXDirection*appleWtoHRatio];
newOrangeSize = [s.picSizeXDirection s.picSizeXDirection*orangeWtoHRatio];
newPeachSize = [s.picSizeXDirection s.picSizeXDirection*peachWtoHRatio];
newBasketballSize = [s.picSizeXDirection s.picSizeXDirection*basketballWtoHRatio];
newButtonSize = [s.picSizeXDirection s.picSizeXDirection*buttonWtoHRatio];
newWheelSize = [s.picSizeXDirection s.picSizeXDirection*wheelWtoHRatio];

s.apple = imresize(s.apple, newAppleSize);
s.orange = imresize(s.orange, newOrangeSize);
s.peach = imresize(s.peach, newPeachSize);
s.basketball = imresize(s.basketball, newBasketballSize);
s.button = imresize(s.button, newButtonSize);
s.wheel = imresize(s.wheel, newWheelSize);
end









