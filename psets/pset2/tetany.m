% tetany.m - simulates the summation of muscle twitches to reach tetany
% For ES 53, designed by L. Moyer
% Last modified 9/16/22

% Students – please add comments to show that you understand 
% what is going on in each line of this code
 
clear all
close all
 
dt = 10; % time period between twitches in ms
frequency = 1/dt; % twitch frequency
numtwitches = 50; % number of twitches to initiate
x = [0:0.1:40]; % this is our time vector for a single twitch, so 400 time points
twamp = zeros(1,2001); % two twitch amplitude arrays, one for fast,
twamps = zeros(1,2001); % one for slow
 
% Generate time delayed twitches and sum
for i= 1:numtwitches;
y(i,:) = gampdf(x,3,1); % approximates a fast twitch response
start = round((i-1)*dt)+1;
twamp(i,start:(start+length(x)-1)) = y(i,:);
end
tetf = sum(twamp,1);
 
figure(5)
subplot(211)
plot(tetf,'r')
axis([0 500 0 1.1*max(max(tetf))])
xlabel('Time (ms)')
ylabel('Twitch Force (a.u.)')
legend('Fast Twitch','location','southeast')
 
 
% Now you generate time delayed slow twitches and sum them.
% Write your own code (or modify below) and plot your results in subplot(212)in blue
for j= 1:numtwitches;
y2(j,:) = gampdf(x,3,2); %approximates a slow twitch response
star = round((j-1)*dt)+1;
twamps(j,star:(star+length(x)-1)) = y2(j,:);
end
tets = sum(twamps,1);
  
subplot(212)
plot(tets,'b')
axis([0 500 0 1.1*max(max(tets))])
xlabel('Time (ms)')
ylabel('Twitch Force (a.u.)')
legend('Slow Twitch','location','southeast')
saveas(gcf, 'q2a.png');
 
