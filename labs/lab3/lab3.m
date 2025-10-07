% File Name: lab2_arodriguez.m
% Author: Andrew Rodriguez
% Created: September 23rd 2025
% Description: Lab 3 Figure 1 / 5 Table 3

clear; clc;


%% Figure 1: ECG/BVP segments from comments, with matched durations
% Expected LabChart MAT fields: data, datastart, dataend, tickrate,
% samplerate, comtext, com, titles

matFile = 'Lab3_Ex1_Defne.mat';   % change if needed
S = load(matFile);

% Normalize strings
if ischar(S.comtext),   S.comtext = cellstr(S.comtext); end
if ischar(S.titles),    S.titles  = cellstr(S.titles);  end

data = S.data(:);
datastart = S.datastart(:);
dataend = S.dataend(:);     % 1-based (inclusive)
samplerate = S.samplerate(:);
tickrate = S.tickrate;
comtext = S.comtext(:);
com = S.com;

% --- choose channel: try channel 2 first, else fall back to ECG by title ---
reqChan = 2;  % your instruction: "second channel of data"
chanIdx = reqChan;

hasData = (chanIdx <= numel(datastart)) && datastart(chanIdx) > 0 && dataend(chanIdx) > 0;
hasFs   = chanIdx <= numel(samplerate) && samplerate(chanIdx) > 0;

if ~hasData || ~hasFs
    % Fall back to a channel titled ECG (common in this lab export)
    ecgIdx = find(contains(upper(S.titles),'ECG'), 1, 'first');
    if ~isempty(ecgIdx)
        chanIdx = ecgIdx;
    else
        % otherwise, first valid channel with data
        valid = find(datastart>0 & dataend>0, 1, 'first');
        if ~isempty(valid), chanIdx = valid; end
    end
end

% Extract signal and sampling rate for the chosen channel
i0 = datastart(chanIdx);         % 1-based
i1 = dataend(chanIdx);           % inclusive
sig = data(i0:i1);
fs  = samplerate(chanIdx); if fs<=0, fs = tickrate; end

% --- helper to get time (s) of a comment for this channel ---
getCommentTime = @(label) local_comment_time(label, comtext, com, chanIdx, tickrate);

t_hand = getCommentTime('Hand movement');
t_jump = getCommentTime('Jumping');

% --- build segments ---
% Segment A: 0 -> Hand movement
lenA = max(0, floor(t_hand * fs));           % samples
segA = sig(1 : min(lenA, numel(sig)));

% Segment B: Jumping -> end, trimmed to same duration as A
iJump = max(1, 1 + floor(t_jump * fs));      % 1-based index into 'sig'
availB = numel(sig) - (iJump - 1);
lenB   = min(lenA, availB);
segB   = sig(iJump : iJump + lenB - 1);

tA = (0:numel(segA)-1) / fs;
tB = (0:numel(segB)-1) / fs;                 % re-zero for overlay

% --- plot ---
figure('Color','w'); clf
plot(tA, segA, 'DisplayName','0 \rightarrow Hand movement'); hold on
plot(tB, segB, 'DisplayName','Jumping \rightarrow end (trimmed)');
grid on
xlabel('Time (s)');
ylabel('ECG/BVP (V)');
ttl = 'Figure 1: Segments from comments (matched duration)';
if exist('titles','var') && chanIdx<=numel(S.titles)
    ttl = sprintf('%s  |  Channel: %s', ttl, strtrim(S.titles{chanIdx}));
end
title(ttl);
legend('Location','best');

% Uncomment to save
% print(gcf, 'Figure1_comments_matched.png', '-dpng', '-r300');
% print(gcf, 'Figure1_comments_matched.pdf', '-dpdf');

%% -------- local function --------
function t = local_comment_time(label, comtext, com, chanIdx, tickrate)
    % label match is case/space-insensitive
    want = find(strcmpi(strtrim(comtext), strtrim(label)), 1, 'first');
    if isempty(want)
        error('Comment label "%s" not found in comtext.', label);
    end
    % com columns (LabChart): [channel, ?, tick_index, ?, text_index]
    rows = com(com(:,1)==chanIdx & com(:,5)==want, :);
    if isempty(rows)
        error('No "%s" comment found for channel %d.', label, chanIdx);
    end
    ticks = rows(1,3);    % use the first occurrence
    t = ticks / tickrate;
end