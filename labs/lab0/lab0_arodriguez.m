% File Name: lab0_arodriguez.m
% Author: Andrew Rodriguez
% Created: September 4th 2025
% Decription: Lab 0, learning how to use Matlab!

% Clear the output when initializing a new script
clear

%% I can make a new section like this
disp("hello world, andrews first matlab assignment")

%% Vector and Matrix Declaration and Operation

mat1 = [1 2 3 4 5];
mat2 = [11 12 13 14 15; 21 22 23 24 25];

disp(mat1 + mat2)
disp(mat1 - mat2)
disp(mat1 .* mat2)

%% Conditionals

% Take a number, n, and if n < 10 -> 0 elif n > 10 -> 1 else 2
n = 10;

if n > 10
    disp(1)
elseif n < 10
    disp(0)
else
    disp(2)
end

% Take numbers, A and B, print the relation (>, <, =) of the A and B
A = 10;
B = 10;

if A == B
    disp("A is equal to B")
elseif A > B
    disp("A is greater than B")
else
    disp("A is less than B")
end

%% Loops

% We have two scenarios summing all odd numbers from 1 to 3533
% and summing all even numbers until the sum is > 3533
% I think the first should be a for loop because we have a clearly defined
% range, and on the other hand the we have a while loop because we have
% a clearly defined conditional

sum_for_loop = 0;
for i = 1:3553
    if mod(i, 2) == 0
        continue
    end
    sum_for_loop = sum_for_loop + i;
end
disp(sum_for_loop)

sum_while_loop = 0;
incr = 0;
while sum_while_loop < 3553
    sum_while_loop = sum_while_loop + incr;
    incr = incr + 2;
end
disp(sum_while_loop)

%% Searching / Indexing

% Random Matrix
rand_mat = 100 * rand(100, 100);

in_range = rand_mat(rand_mat < 53 & rand_mat > 35);
% Not displayed for submission

% Indexing
rand_vec = 10 * rand(1, 100);
squared = rand_vec .^ 2;

indices = find(squared > 50);
smaller_than_50 = rand_vec(indices);

disp(size(smaller_than_50))

%% Functions

rand_vec = 10 * rand(1, 100); % Just redeclaring this to be safe
idx = 1:100;

[new_idx, new_vec] = greaterThan(idx, rand_vec, 5);
greater_than_5 = new_vec; % This is all numbers above 5

% Now for less than or equal to 5
[ids, new_vec_flipped] = greaterThan(idx, -rand_vec, -5);

% Now we just flip back 
leq_than_5 = -new_vec_flipped;

% Now to make a function called lessThan
function [newXs, newYs] = lessThan(xs, ys, threshold)
    lessThanInd = find(ys <= threshold);
    % Extracting the values where the y-values meet the condition
    newXs = xs(lessThanInd);
    newYs = ys(lessThanInd);
end

[ids50, below_50] = lessThan(idx, rand_vec, 50);
% So below 50 is our set below 50 which is all nums cause 0-10 < 50

%% Using Built-In MATLAB Functions and Documentation

fs = 4000; % 4khz 
T = 3; % 3 seconds
t = 0:1/fs:T-1/fs; % time vector
q = chirp(t-2,4,1/2,6,'quadratic',100,'convex').*exp(-4*(t-1).^2);
figure;
plot(t,q)

[up,lo] = envelope(q);
hold on
plot(t, up, 'r', 'LineWidth', 4)
plot(t, lo, 'b', 'LineWidth', 2)
legend('q','up','lo')
hold off


%% Snipping Data
load("SnippingDataExample.mat")
ch = 2;
bl = 2;

% So now we pull out the data for that block
pdata = data(datastart(ch,bl):dataend(ch,bl));
fs = samplerate(ch,bl);   

% So I want the comments in the block and channel I've picked
mask = (com(:,2) == bl) & (com(:,1) == ch | com(:,1) == -1);
idxs = find(mask);
% These are already sorted so I can just pull out the first (0 indexing)
idx_second = idxs(1);

t0 = com(idx_second,3) / tickrate(bl); % Seconds from start
start_samp = max(1, round(t0*fs) + 1);
N = round(10*fs); % Add 10 seconds
stop_samp  = min(start_samp + N - 1, numel(pdata)); % Then find the stop

y_snip = pdata(start_samp:stop_samp); % Snip the y data
t_snip = (0:numel(y_snip)-1) / fs; % Time relative to the comment

figure; 
plot(t_snip, y_snip, 'b-'); grid on
xlabel('Time since 2^{nd} comment (s)');
ylabel('Amplitude (v)')
title(sprintf('Channel 2, Block 2, 10 secs from second comment'));

%% Section to clear all output history
close all
clc