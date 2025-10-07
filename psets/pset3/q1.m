%% q1a
a = 20;
b = 0.2; 
vmax = 1.0; 
T0 = a*vmax/b;
v = linspace(0, vmax, 400);
T = ((T0 + a)*b ./ (v + b)) - a;  

figure; hold on; box on; grid on;
plot(T, v, 'LineWidth', 2);
% scatter([0, vmax], [T0, 0], 60, 'filled'); 
xlabel('Velocity v (m/s)');
ylabel('Force T (N)');
title('Hill Force–Velocity Relationship');
% xlim([0 vmax]); ylim([0 T0*1.05]);
text(0, T0, '  (T_0 = 100)', 'VerticalAlignment','bottom');
text(vmax, 0, '  (v_{max}, 0)', 'VerticalAlignment','top');
saveas(gcf,"q1a.png")


%% q1b
T0 = 100;
a_vals = 5:5:30;
v = linspace(0, vmax, 800);

figure; hold on; box on; grid on;
for a = a_vals
    b = (a*vmax)/T0;           
    T = ((T0 + a)*b ./ (v + b)) - a;  
    plot(v, T, 'LineWidth', 2, ...
         'DisplayName', sprintf('a = %g N, b = %.3f m/s', a, b));
end
scatter([0, vmax], [T0, 0], 60, 'filled', 'MarkerFaceAlpha', 0.8);
text(0, T0, '  (0, T_0)', 'VerticalAlignment','bottom');
text(vmax, 0, '  (v_{max}, 0)', 'VerticalAlignment','top');

xlabel('Velocity v (m/s)'); ylabel('Force T (N)');
title('Hill Force–Velocity Curves for varying a (T0 = 100 N vmax =1 m/s)');
xlim([0 vmax]); ylim([0 T0*1.05]); legend('Location','best');
saveas(gcf,"q1b.png")