%% Grasp_Force_Testing: GraspForceTestingAnalysis
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script is used to plot the differences between the modeled contact
%  force and the measureed contact force.
%  
%  Outputs:
%  N/A

clear all; close all; clc;

addpath '..\5B3D-Finger-Optimization'
addpath '..\5B3D-Finger-Optimization\Results'
% load('Data/Matched_Data_60x45.mat'); 
% load('Data/Matched_Data_60x60.mat'); 
% load('Data/Matched_Data_60x75.mat'); 
% load('Data/Matched_Data_75x75.mat'); 
% load('Data/Matched_Data_90x90.mat'); 


colorOrder = ['r', 'g', 'b', 'm', 'c'];
for ii = 1:5
    start = 1;
    jj = 0;
    if ii ==1
        t1_des{ii} = deg2rad(60);
        t2_des{ii} = deg2rad(45);
        load('Data/Matched_Data_60x45.mat'); 
        fig = figure(1);
        fig.Position = [100 100 550 400];
        hold on
        % plot([0; 30], [0; 30], '--k', 'LineWidth', 2)
        plot([0; 30], [0; 0], '--k', 'LineWidth', 2)
        fig2 = figure(2);
        fig2.Position = [100 600 550 400];
        hold on
        title('$\theta_1$ Error', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        xlabel('Desired $\theta_1$ (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        ylabel('Measured $\theta_1$ (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        plot([deg2rad(30), deg2rad(100)], [deg2rad(30), deg2rad(100)], ...
            'k', 'LineWidth', 2)
        fig3 = figure(3);
        fig3.Position = [650 600 550 400];
        hold on
        title('$\theta_2$ Error', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        xlabel('Desired $\theta_2$ (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        ylabel('Measured $\theta_2$ (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        plot([deg2rad(30), deg2rad(100)], [deg2rad(30), deg2rad(100)], ...
            'k', 'LineWidth', 2)
        fig5 = figure(5);
        fig5.Position = [1200 600 300 400];
        title('Joint Angle Error', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        ylabel('Angle Error (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
        hold on
    elseif ii == 2
        t1_des{ii} = deg2rad(60);
        t2_des{ii} = deg2rad(60);
        load('Data/Matched_Data_60x60.mat'); 
    elseif ii == 3
        t1_des{ii} = deg2rad(60);
        t2_des{ii} = deg2rad(75);
        load('Data/Matched_Data_60x75.mat'); 
    elseif ii == 4
        t1_des{ii} = deg2rad(75);
        t2_des{ii} = deg2rad(75);
        load('Data/Matched_Data_75x75.mat'); 
    else 
        t1_des{ii} = deg2rad(90);
        t2_des{ii} = deg2rad(90);
        load('Data/Matched_Data_90x90.mat'); 
    end
    for kk = 1:length(fin)
        try
            if start == 1
                t1{ii} = NaN*ones(length(fin{kk}.t1), 6);
                t2{ii} = NaN*ones(length(fin{kk}.t1), 6);
                t3{ii} = NaN*ones(length(fin{kk}.t1), 6);
                F{ii} = NaN*ones(length(fin{kk}.t1), 6);
                Fx{ii} = NaN*ones(length(fin{kk}.t1), 6);
                Fy{ii} = NaN*ones(length(fin{kk}.t1), 6);
                Fz{ii} = NaN*ones(length(fin{kk}.t1), 6);
                FT{ii} = NaN*ones(length(fin{kk}.t1), 6);
                Ref{ii} = NaN*ones(length(fin{kk}.t1), 6);
                Temp{ii} = NaN*ones(length(fin{kk}.t1), 6);
                trial{ii} = NaN*ones(length(fin{kk}.t1), 6);
                start = 0;
            end
            t1{ii}(:, kk) = fin{kk}.t1;
            t2{ii}(:, kk) = fin{kk}.t2;
            t3{ii}(:, kk) = fin{kk}.t3;
            F{ii}(:, kk) = movmean(fin{kk}.F, 10); % 1 second moving average
            Fx{ii}(:, kk) = movmean(rok{kk}.Fx, 10); % 1 second moving average
            Fy{ii}(:, kk) = movmean(rok{kk}.Fy, 10); % 1 second moving average
            Fz{ii}(:, kk) = movmean(rok{kk}.Fz, 10); % 1 second moving average
            FT{ii}(:, kk) = sqrt(Fx{ii}(:, kk).^2 + Fy{ii}(:, kk).^2);
            Ref{ii}(:, kk) = fin{kk}.Ref;
            Temp{ii}(:, kk) = rok{kk}.Temperature_C;
            trial{ii}(:, kk) = kk*ones(size(fin{kk}.t1));
        end
    end

    
    %% Angle Error Analysis
    figure(2)
    hold on
    temp = ~isnan(unique(t1{ii}));
    t1_unique = unique(t1{ii}); t1_unique = t1_unique(temp);
    scatter(t1_des{ii}*ones(size(t1_unique)), t1_unique, 30, colorOrder(ii), 'filled', 'o')
    figure(3)
    hold on
    temp = ~isnan(unique(t2{ii}));
    t2_unique = unique(t2{ii}); t2_unique = t2_unique(temp);
    scatter(t2_des{ii}*ones(size(t2_unique)), t2_unique, 30, colorOrder(ii), 'filled', 'o')
    t1_error{ii} = t1_unique-t1_des{ii}*ones(size(t1_unique));
    t2_error{ii} = t2_unique-t2_des{ii}*ones(size(t2_unique));
    %% Format plots
    if ii == 5
        figure(2)
        axis([deg2rad(55) deg2rad(95) deg2rad(55) deg2rad(95)])
        figure(3)
        axis([deg2rad(40) deg2rad(95) deg2rad(40) deg2rad(95)])
    end

    clear start jj
end
F_mean = mean(cell2mat(F), 2, 'omitnan');
idx = find(F_mean == max(F_mean), 1, 'first');
F_std = std(cell2mat(F), [], 2, 'omitnan');
FT_mean = mean(cell2mat(FT), 2, 'omitnan');
FT_std = std(cell2mat(FT), [], 2, 'omitnan');
F_max = (FT_mean+FT_std)./(F_mean-F_std).*F_mean;
F_min = (FT_mean-FT_std)./(F_mean+F_std).*F_mean;
x_patch = [F_mean(1:idx); flip(F_mean(1:idx))];
y_patch = [F_max(1:idx)-F_mean(1:idx); flip(F_min(1:idx)-F_mean(1:idx))];
figure(1)
plot(F_mean(1:idx), FT_mean(1:idx)-F_mean(1:idx), 'k', 'LineWidth', 2)
patch(x_patch, y_patch, 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none')

axis([6.5 25 -4 3])
xlabel('Estimated Force (N)', 'FontSize', 16, 'FontName', 'Times')
ylabel('Grasp Force Error (N)', 'FontSize', 16, 'FontName', 'Times')
lgd = legend('Ideal', 'Mean Error', '$\pm \ 1\sigma$', 'interpreter', 'latex', ...
    'FontSize', 16, 'FontName', 'Times', 'location', 'northwest');
lgd.EdgeColor = 'none'; lgd.Color = 'none';

%% Angle Error plot
t1Error = [t1_error{1}; t1_error{2}; t1_error{3}; t1_error{4}; t1_error{5}];
t2Error = [t2_error{1}; t2_error{2}; t2_error{3}; t2_error{4}; t2_error{5}];
combinedData = [t1Error; t2Error];
groups1 = ones(length(t1Error), 1);
groups2 = 2 * ones(length(t2Error), 1);
groupingVar = [groups1; groups2];

figure(5)
hold on
boxplot(combinedData, groupingVar, 'labels', {'t1', 't2'}, 'symbol', '')
xlim([0.75 2.25])
% ylim([0 0.08])
ylim([-0.06 -0.01])

%% Unpack the force plot
% figure(1)
% ax = gca; ax = ax.Children;
finger.x60x45 = F{1}(:);
rokubi.x60x45 = FT{1}(:);
finger.x60x60 = F{2}(:);
rokubi.x60x60 = FT{2}(:);
finger.x60x75 = F{3}(:);
rokubi.x60x75 = FT{3}(:);
finger.x75x75 = F{4}(:);
rokubi.x75x75 = FT{4}(:);
finger.x90x90 = F{5}(:);
rokubi.x90x90 = FT{5}(:);


%% Plot the boxplot of the grasp force error
Error = -[(finger.x60x45 - rokubi.x60x45); ...
        (finger.x60x60 - rokubi.x60x60); ...
        (finger.x60x75 - rokubi.x60x75); ...
        (finger.x75x75 - rokubi.x75x75); ...
        (finger.x90x90 - rokubi.x90x90)];
est = [finger.x60x45; finger.x60x60; finger.x60x75; ...
    finger.x75x75; finger.x90x90];
fig4 = figure(4);
fig4.Position = [650 100 250 400];
Error(isoutlier(Error, 'quartiles')) = NaN;
b = boxplot(Error, 'labels', 'N=1847', 'symbol', ''); %, 'FontSize', 12)
% title('Total Grasping Force Error', 'FontSize', 12)
ylabel('Grasp Force Error (N)', 'FontSize', 16, 'FontName', 'Times')
xlim([0.85 1.15])


%% Total Force Error Fitting with Mixed-Effect Model (low correlation with t1 and t2)
Fx = [Fx{1}; Fx{2}; Fx{3}; Fx{4}; Fx{5}];
Fy = [Fy{1}; Fy{2}; Fy{3}; Fy{4}; Fy{5}];
F = [F{1}; F{2}; F{3}; F{4}; F{5}];
Err = (sqrt(Fx.^2+Fy.^2)-F);
t1 = rad2deg([t1_des{1}*ones(size(t1{1})); t1_des{2}*ones(size(t1{2})); ...
    t1_des{3}*ones(size(t1{3})); t1_des{4}*ones(size(t1{4})); ...
    t1_des{5}*ones(size(t1{5}))]);
t2 = rad2deg([t2_des{1}*ones(size(t2{1})); t2_des{2}*ones(size(t2{2})); ...
    t2_des{3}*ones(size(t2{3})); t2_des{4}*ones(size(t2{4})); ...
    t2_des{5}*ones(size(t2{5}))]);
t3 = rad2deg([t3{1}; t3{2}; t3{3}; t3{4}; t3{5}]);
Temp = [Temp{1}; Temp{2}; Temp{3}; Temp{4}; Temp{5}];
trial = [trial{1}; trial{2}; trial{3}; trial{4}; trial{5}];
for ii = 1:numel(t1)
    if round(t1(ii)) == 60 && round(t2(ii)) == 45
        config(ii) = 1;
    elseif round(t1(ii)) == 60 && round(t2(ii)) == 60
        config(ii) = 2;
    elseif round(t1(ii)) == 60 && round(t2(ii)) == 75
        config(ii) = 3;
    elseif round(t1(ii)) == 75 && round(t2(ii)) == 75
        config(ii) = 4;
    elseif round(t1(ii)) == 90 && round(t2(ii)) == 90
        config(ii) = 5;        
    else
        config(ii) = NaN;
    end
end
tbl = array2table([Err(:), Fx(:), Fy(:), F(:), t1(:), t2(:), t3(:), Temp(:), trial(:), config(:)], ...
    'VariableNames', ["Error", "Fx", "Fy", "F", "t1", "t2", "t3", "Temp", "trial", "Config"]);
tbl = rmmissing(tbl);
lme = fitlme(tbl, 'Error ~ t1*t2')

fig6 = figure(6);
fig6.Position = [900 100 550 400];
title('')
xlabel('$\theta_1$ (deg)', 'interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times')
ylabel('$\theta_2$ (deg)', 'interpreter', 'latex', 'FontSize', 16, 'FontName', 'Times')
zlabel('Grasp Force Error (N)', 'FontSize', 16, 'FontName', 'Times')
boxPlot3D(tbl.Error(1:10:end),tbl.t1(1:10:end),tbl.t2(1:10:end)) % downsample
view([45 22.5])
ax = gca; ax = ax.Children;
jj=1;
for ii = 1:length(ax)
    if strcmp(ax(ii).Type, 'patch')
        if isscalar(unique(ax(ii).ZData)) == 1
            if (round(mean(ax(ii).XData, 'all')) == 60)
                if (round(mean(ax(ii).YData, 'all')) == 45)
                    l{1} = plot(NaN, NaN, 'r', 'LineWidth', 2);
                    ax(ii).FaceColor = 'red';
                elseif (round(mean(ax(ii).YData, 'all')) == 60)
                    l{2} = plot(NaN, NaN, 'g', 'LineWidth', 2);
                    ax(ii).FaceColor = 'green';
                else
                    l{3} = plot(NaN, NaN, 'b', 'LineWidth', 2);
                    ax(ii).FaceColor = 'blue';
                end
            elseif (round(mean(ax(ii).XData, 'all')) == 75)
                l{4} = plot(NaN, NaN, 'm', 'LineWidth', 2);
                ax(ii).FaceColor = 'magenta';
            else % 90
                l{5} = plot(NaN, NaN, 'c', 'LineWidth', 2);
                ax(ii).FaceColor = 'cyan';
            end
        end
    end
end
lgd = legend([l{1}, l{2}, l{3}, l{4} l{5}], ...
    {'(60$^\circ$, 45$^\circ$)', '(60$^\circ$, 60$^\circ$)', ...
    '(60$^\circ$, 75$^\circ$)', '(75$^\circ$, 75$^\circ$)', ...
    '(90$^\circ$, 90$^\circ$)'}, 'interpreter', 'latex', ...
    'FontSize', 14, 'FontName', 'Times');
lgd.EdgeColor = 'none'; lgd.Color = 'white'; lgd.Position = [0.4284 0.15 0.1782 0.290];
lgd.Box = 'off';

fig7 = figure(7);
fig7.Position = [1450 100 550 400];
boxplot(tbl.Error, tbl.Config, 'labels', {'(60°, 45°)', '(60°, 60°)', ...
    '(60°, 75°)', '(75°, 75°)', '(90°, 90°)'}, 'symbol', '');


exportgraphics(fig, 'Grasp Force Error.png', 'Resolution', 1000);
exportgraphics(fig4, 'Grasp Force Error-Box.png', 'Resolution', 1000);
exportgraphics(fig6, 'Grasp Force Error-Box-3D.png', 'Resolution', 1000);


anova(tbl, 'Error ~ Config')

% ppd = plotPartialDependence(lme, [1,2]);
% ppd.Children.EdgeColor = 'none'; ppd.Children.FaceColor = 'k'; ppd.Children.FaceAlpha = 0.5;
% hold on
% scatter3(tbl.t1, tbl.t2, tbl.Error, 20, '.k')
% title('')
% xlabel('$\theta_1$ (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
% ylabel('$\theta_2$ (rad)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
% zlabel('Grasp Force Error (N)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times')
