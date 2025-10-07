%%
% File Name: greaterThan.m
% Author: William Cho
% Created: Sept 6, 2017
% Description: A function that takes in x-values, y-values, and returns the
% x-values and y-values where the y-value exceeds a threshold.

function [newXs,newYs] = greaterThan(xs, ys, threshold)
    % Declaring boolean vector with indices where the y-values exceed the
    % threshold
    greaterThanInd = find(ys > threshold);
    % Extracting the values where the y-values exceed the threshold
    newXs = xs(greaterThanInd);
    newYs = ys(greaterThanInd);
end