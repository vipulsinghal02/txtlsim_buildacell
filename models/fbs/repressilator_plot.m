% repressilator_plot.m - run the repressilator model
% RMM, 24 Jan 07

figure(1); clf; subplot(221);

sol = ode45(@repressilator, [0 20000], [1 0 0 200 0 0]);
xy = plot(...
  sol.x/60, sol.y(4,:), 'b-', ...
  sol.x/60, sol.y(5,:), 'r--', ...
  sol.x/60, sol.y(6,:), 'k-.');
set(xy, 'LineWidth', AM_data_linewidth);

axis([0, 300, 0, 5000]);
xlabel('Time {\itt} [min]');
ylabel('Proteins per cell');
lgh = legend(' cI', ' lacI', ' tetR');
legend(lgh, 'boxoff');

daspect([110 2500 1]);

% As a check, plot out the mRNA concentrations
figure(2); clf;
xy = plot(...
  sol.x/60, sol.y(1,:), 'b-', ...
  sol.x/60, sol.y(2,:), 'r--', ...
  sol.x/60, sol.y(3,:), 'k-.');
axis;
lgh = legend(' cI', ' lacI', ' tetR');
