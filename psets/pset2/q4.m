%$ q4a
t = linspace(0,5,500); % over 5 seconds, 500 data points
y = exp(-t); 

figure('Color','w');
plot(t, y, 'LineWidth', 1.6); grid on; box on;
hold on;
plot(1, exp(-1), 'r*', 'MarkerSize', 10);                 % red asterisk
text(1, exp(-1), sprintf('  (t=1, y=e^{-1}\\approx %.3f)', exp(-1)), ...
     'VerticalAlignment','bottom');

xlabel('t (s)'); ylabel('y'); title('y = e^{-t}');
exportgraphics(gcf, 'q4a.png', 'Resolution', 300);

%% q4b
y2 = 1-y;
figure('Color','w');
plot(t, y2, 'LineWidth', 1.6, 'Color', 'red'); grid on; box on;
hold on;
plot(1, 1-exp(-1), 'r*', 'MarkerSize', 10);                 % red asterisk
text(1, 1-exp(-1), sprintf('  (t=1, y=1-e^{-1}\\approx %.3f)', 1-exp(-1)), ...
     'VerticalAlignment','bottom');

xlabel('t (s)'); ylabel('y'); title('y = 1 - e^{-t}');
exportgraphics(gcf, 'q4b.png', 'Resolution', 300);

%% q4c

A0  = 2;
tau = 3;   
t = linspace(0, 15, 600);
A = A0 * exp(-t/tau);

figure('Color','w');
plot(t, A, 'LineWidth', 1.6); grid on; box on; hold on;
plot(tau, A0*exp(-1), 'r*', 'MarkerSize', 10);      % marker at t = tau
xlabel('t (s)'); ylabel('A(t)'); title('A(t) = A_0 e^{-t/\tau},  A_0=2, \tau=3 s');
exportgraphics(gcf, 'q4c.png', 'Resolution', 300);

%% q4d
A0  = 2;            % from part (c)
tau = 3;            % seconds

t_half = tau*log(2);     % when A = A0/2
t_37   = tau;            % when A = A0/e ≈ 0.37*A0

fprintf('(d) t_half = tau*ln(2) = %.3f s\n', t_half);
fprintf('(d) t_(37%% of A0) = tau = %.3f s\n', t_37);
