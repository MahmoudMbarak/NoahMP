%Read Infiltration Rate
Infil = 10^-7;  % Can be read from file



noahmp.config.domain.NumSoilLayer = 30;
noahmp.config.domain.DepthSoilLayer = - (1:noahmp.config.domain.NumSoilLayer)/ noahmp.config.domain.NumSoilLayer; % Assumed Constant Spacing
noahmp.config.domain.ThicknessSnowSoilLayer = (1/noahmp.config.domain.NumSoilLayer) * ones(noahmp.config.domain.NumSoilLayer,1);

noahmp.water.param.SoilMoistureSat = 0.434*ones(noahmp.config.domain.NumSoilLayer,1);
noahmp.water.param.SoilWatConductivitySat = 3.47E-05*ones(noahmp.config.domain.NumSoilLayer,1); % Sand Saturated Conductivity from the Naohmp Table
noahmp.water.param.SoilExpCoeffB = 4.90*ones(noahmp.config.domain.NumSoilLayer); % Sand Coeffiecient from the NoahMP Table
noahmp.water.param.SoilWatDiffusivitySat = 8.05E-06*ones(noahmp.config.domain.NumSoilLayer); % Sand Saturated Diffusivity from the NoahMP Table
noahmp.water.param.SoilDrainSlope = 0.5;

noahmp.water.state.SoilMoisture = 0.1*ones(noahmp.config.domain.NumSoilLayer,1);


noahmp.water.state.SoilImpervFrac = zeros(noahmp.config.domain.NumSoilLayer,1); % For all soil layers
noahmp.water.state.SoilIce = zeros(noahmp.config.domain.NumSoilLayer,1);
noahmp.water.state.SoilLiqWater = noahmp.water.state.SoilMoisture;

noahmp.water.flux.InfilRateSfc =  0;
noahmp.water.flux.TranspWatLossSoilMean = zeros(noahmp.config.domain.NumSoilLayer,1);
noahmp.water.flux.TranspWatLossSoilMean(noahmp.config.domain.NumSoilLayer,1) = 0;
noahmp.water.flux.EvapSoilSfcLiqMean = 0; %.8899248E-10;

% Define the infiltration rate values
%infiltration_rates = [0.01, 0, 0, 0, 0, 0.005, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
[SoilWatDiffusivity, SoilWatConductivity] = SoilDiffusivityConductivity(noahmp, noahmp.water.state.SoilMoisture(1), 0, 1, 1);
infiltration_rates = Infil;


% Initialize SM to store soil moisture values
nts = 100 %number of timesteps 
SM = zeros(noahmp.config.domain.NumSoilLayer, nts);

% Store initial condition
SM(:, 1) = noahmp.water.state.SoilMoisture;

% Main loop
for i = 1:nts
    % Update infiltration rate
    noahmp.water.flux.InfilRateSfc = Infil;
    
    % Run your existing code
    [MatRight, MatLeft1, MatLeft2, MatLeft3, DrainSoilBot, SoilWatConductivity, SoilWatDiffusivity] = SoilWaterDiffusionRichards(noahmp);
    TimeStep = 3600;
    
    % Solve for soil moisture
    noahmp = SoilMoistureSolver(noahmp, TimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight);
    
    % Store the updated soil moisture
    SM(:, i+1) = noahmp.water.state.SoilMoisture;
end

% Create animation for theta vs depth
fig = figure('Position', [100, 100, 800, 600]);
ax = axes(fig);

% Define depth axis (convert to positive values for plotting)
depth = abs(noahmp.config.domain.DepthSoilLayer);

% Find overall min and max for consistent x-axis limits
theta_min = min(SM(:)) - 0.01;  % Add small buffer
theta_max = max(SM(:)) + 0.01;

% Color for the line
line_color = [0, 0.4470, 0.7410];  % MATLAB default blue

for i = 1:size(SM, 2)
    % Clear the axes for each frame
    cla(ax);
    
    % Plot theta vs depth
    plot(ax, SM(:,i), depth, '-o', 'LineWidth', 2.5, 'MarkerSize', 4, ...
         'Color', line_color, 'MarkerFaceColor', line_color);
    
    % Labels and title
    xlabel(ax, 'Soil Moisture (\theta)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Depth (normalized)', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax, sprintf('Soil Moisture Profile - Timestep %d', i-1), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Set consistent axis limits
    xlim(ax, [theta_min, theta_max]);
    ylim(ax, [0, 1]);  % Normalized depth from 0 to 1
    
    % Invert y-axis to show depth increasing downwards
    set(ax, 'YDir', 'reverse');
    
    % Add grid for better readability
    grid(ax, 'on');
    grid(ax, 'minor');
    
    % Improve appearance
    set(ax, 'FontSize', 10);
    set(ax, 'LineWidth', 1.2);
    
    % Add text box with current values
    textstr = sprintf('Surface \\theta: %.4f\nBottom \\theta: %.4f', SM(1,i), SM(end,i));
    annotation('textbox', [0.02, 0.02, 0.25, 0.15], 'String', textstr, ...
               'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Ensure figure is drawn
    drawnow;
    
    % Optional: pause to see the animation in real-time
    pause(0.1);
end

fprintf('Animation complete. Total frames: %d\n', size(SM, 2));

% Optional: Create a final static plot of all profiles overlaid
figure('Position', [200, 200, 800, 600]);
hold on;
colors = jet(size(SM, 2));  % Color map for different timesteps

for i = 1:size(SM, 2)
    plot(SM(:,i), depth, 'LineWidth', 1.5, 'Color', colors(i,:));
end

xlabel('Soil Moisture (\theta)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Depth (normalized)', 'FontSize', 12, 'FontWeight', 'bold');
title('All Soil Moisture Profiles Over Time', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse');
grid on;
colormap(jet);
c = colorbar;
c.Label.String = 'Timestep';
c.Label.FontSize = 12;
c.Label.FontWeight = 'bold';
caxis([0, size(SM, 2)-1]);
