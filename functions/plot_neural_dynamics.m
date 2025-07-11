function plot_neural_dynamics(projectedData, trialdef)

% Step 5: Visualize Neural Manifold
figure;
plot3(projectedData(:, 1), projectedData(:, 2), projectedData(:, 3));
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
title('Neural Manifold Visualization');
grid on;

end