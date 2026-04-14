%% Grasp_Force_Testing: interpData
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script syncronizes the time series data from the finger controller
%  and the F/T sensor. The data stream from the finger has a much lower
%  sample rate than the F/T sensor, so the time series data are made the 
%  same length by interpolation.
%  
%  Outputs:
%  Matched_Data.mat -> Structure of trials with the finger data and the 
%  6-axis force and torque readings from the F/T sensor with time stamps.

clear all; close all; clc;
for jj = 1:3
    if jj == 1
        load('ANTERO_Data_60x45_5_Hz.mat'); % load the unpacked teensy data. produces a struct called 'finger'
    elseif jj == 2 
        clear ANTERO a0 ant f fa fieldNames id kk Name p_F p_Ref pw
        load('ANTERO_Data_60x60_5_Hz.mat'); % load the unpacked teensy data. produces a struct called 'finger'
    else
        clear ANTERO a0 ant f fa fieldNames id kk Name p_F p_Ref pw
        % load('ANTERO_Data_60x75_5_Hz.mat'); % load the unpacked teensy data. produces a struct called 'finger'
        load('ANTERO_Data_60x75_5_Hz.mat'); % load the unpacked teensy data. produces a struct called 'finger'
    end

    % load('Rokubi_Data_60x60_2_Hz.mat'); % load the unpacked rokubi. produces a struct called 'rokubi'
    
    fieldNames = fieldnames(ANTERO);
    kk = 1;
    fig = figure(1);
    fig.Position = [100 600 550 400];
    for ii = 1:numel(fieldNames)
        if jj ~= 2 || ii ~= 3 % omit this data
            tic;
            Name = fieldNames{ii};
            ant{ii} = ANTERO.(Name);
      
            a0 = find(ant{ii}.Ref < 12, 1, 'first'); 
            id{ii} = a0;
    
            ant{ii} = ant{ii}(end-28000:end, :); ant{ii}.Time = ant{ii}.Time - ant{ii}.Time(1);
            rf = 1./ant{ii}.Time(end); % minimum resolvable frequency
            rat = (10-min(ant{ii}.F))/2;
            pw_max = 20*log10(rat);
            
            %% FFT
            fa = round(1/mean(diff(ant{ii}.Time)));
            omega_comp = rf:rf:5.554;
            [p_F{ii}, f{ii}] = plomb(ant{ii}.F, ant{ii}.Time, omega_comp, 'power'); % max excitation frequency is around 5.5 Hz
            [p_Ref{ii}, ~] = plomb(ant{ii}.Ref, ant{ii}.Time, omega_comp, 'power');
            pw{ii} = p_F{ii}./p_Ref{ii};
            pw{ii} = 10*log10(pw{ii}); % convert power to magnitude for FRF
            kk = find(f{ii}>1, 1, 'first');
            
            pw{ii} = pw{ii} + (pw_max - max(pw{ii}));
    
            figure(1)
            % plot(f{ii}, pw{ii}, 'LineWidth', 1.0);
            hold on;
        end
    end
    temp = mean(cell2mat(f), 2);
    temp_min = min(cell2mat(pw), [], 2);
    temp_max = max(cell2mat(pw), [], 2);

    fs = 3;
    for ii = 1:floor(numel(temp)/fs) % downsample by factor of fs
        id1 = (ii-1)*fs + 1;
        id2 = ii*fs;
        freq{jj}(ii, 1) = median(temp(id1:id2));
        pw_temp = cell2mat(pw); % average trials together
        mag{jj}(ii, 1) = median(pw_temp(id1:id2, :), 'all');
        mag_std{jj}(ii, 1) = std(pw_temp(id1:id2, :), [], 'all');
        % mag_min{jj}(ii, 1) = median(temp_min(id1:id2), 'all');
        % mag_max{jj}(ii, 1) = median(temp_max(id1:id2), 'all');
    end

    % mag{jj} = movmean(mean(cell2mat(pw), 2), fs); % 1 mHz sampling -> about 
    % mag_std{jj} = movmean(std(cell2mat(pw), [], 2), 10*fs);
    x_patch{jj} = [freq{jj}; flip(freq{jj})];
    upper = mag{jj}+mag_std{jj};
    lower = mag{jj}-mag_std{jj};
    for kk = 1:numel(upper)
        if kk ~=1 && kk ~= numel(upper)
            temp_u(kk) = max([upper(kk-1), upper(kk), upper(kk+1), lower(kk-1), lower(kk), lower(kk+1)]);
            temp_l(kk) = min([upper(kk-1), upper(kk), upper(kk+1), lower(kk-1), lower(kk), lower(kk+1)]);
        else
            temp_u(kk) = upper(kk);
            temp_l(kk) = lower(kk);
        end
    end
    upper = temp_u'; lower = temp_l';
    y_patch{jj} = [upper; flip(lower)];

    if jj == 2
        temp = [ant{1}.Time, ant{2}.Time, ant{4}.Time, ant{5}.Time];
        Time{jj} = mean(temp, 2);
        temp = [ant{1}.Ref, ant{2}.Ref, ant{4}.Ref, ant{5}.Ref];
        Ref{jj} = mean(temp, 2);
        temp = [ant{1}.F, ant{2}.F, ant{4}.F, ant{5}.F];
        F{jj} = mean(temp, 2);
    else
        temp = [ant{1}.Time, ant{2}.Time, ant{3}.Time, ant{4}.Time, ant{5}.Time];
        Time{jj} = mean(temp, 2);
        temp = [ant{1}.Ref, ant{2}.Ref, ant{3}.Ref, ant{4}.Ref, ant{5}.Ref];
        Ref{jj} = mean(temp, 2);
        temp = [ant{1}.F, ant{2}.F, ant{3}.F, ant{4}.F, ant{5}.F];
        F{jj} = mean(temp, 2);
    end
    %% Plotting the average magnitude
    figure(1)
    plot(freq{jj}, mag{jj}, 'k', 'LineWidth', 1.5);
    grid on;
    title('Frequency Response of the Froce Control System');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([0.05 7])
end

fig2 = figure(2);
fig2.Position = [650 600 550 300];
%% 60x45

% semilogx([0.01; 6], [-3; -3], '--k', 'LineWidth', 2)
% hold on;
semilogx(freq{1}, movmean(mag{1}, 3), 'r', 'LineWidth', 2);
% s1 = semilogx(freq{1}, movmean(mag_max{1}, 1), '.r', 'LineWidth', 2);
hold on;
patch(x_patch{1}, y_patch{1}, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
%% 60x60
semilogx(freq{2}, movmean(mag{2}, 3), 'g', 'LineWidth', 2);
patch(x_patch{2}, y_patch{2}, 'g', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
%% 60x75
semilogx(freq{2}, movmean(mag{3}, 3), 'b', 'LineWidth', 2);
patch(x_patch{3}, y_patch{3}, 'b', 'FaceAlpha', 0.25, 'EdgeColor', 'none')

xlim([0.05 5.5])
ylim([-10 0])
grid on;
xlabel('Excitation Frequency (Hz)', 'FontName', 'Times', 'FontSize', 18)
ylabel('Magnitude (dB)', 'FontName', 'Times', 'FontSize', 18)
title("B)", 'FontName', 'times', 'FontSize', 20)
% lgd = legend('(60$^\circ$, 45$^\circ$)', '', '(60$^\circ$, 60$^\circ$)', '', ...
%     '(60$^\circ$, 75$^\circ$)', '', 'interpreter', 'latex', ...
%     'FontName', 'Times', 'FontSize', 18, 'location', 'southwest');
% lgd.EdgeColor = 'none'; lgd.Color = 'none'; lgd.Position = [0.2827 0.1829 0.2755 0.3047];

tf = 300;
rate = tf*tf*tf*tf/pi*0.35; % frequency over time
omega = Time{3}.^4/rate;

fig3 = figure(3);
fig3.Position = [650 100 550 400];
plot(Time{3}, Ref{3}, ':k', 'LineWidth', 2)
hold on
plot(Time{3}, F{3}, 'k', 'LineWidth', 2)
xlim([0 120])
ylim([7.5 12.5])
xlabel('Time (s)', 'FontName', 'Times', 'FontSize', 18)
ylabel('Total Grasp Force (N)', 'FontName', 'Times', 'FontSize', 18)
title("A)", 'FontName', 'times', 'FontSize', 20)
lgd2 = legend('Reference', 'Measured', 'interpreter', 'latex', ...
    'FontName', 'Times', 'FontSize', 18, 'location', 'southwest');
lgd2.EdgeColor = 'none'; lgd.Color = 'none';

fig4 = figure(4);
fig4.Position = [100 100 550 300];
clear upper lower temp_u temp_l
idx = 1;
for ii = 1:3
    time = Time{ii}(idx):0.0001:Time{3}(end);
    R = interp1(Time{ii}(idx:end), movmean(Ref{ii}(idx:end), 1), time, 'spline');
    f = interp1(Time{ii}(idx:end), movmean(F{ii}(idx:end), 1), time, 'spline');
    om = interp1(Time{ii}(idx:end), omega(idx:end), time, 'spline')';
    id1 = find((f(1:end-1)-10).*(f(2:end)-10) < 0); id1 = id1(:); id1f = movmean(id1(:), 1);
    id2 = find((R(1:end-1)-10).*(R(2:end)-10) < 0); id2 = id2(:); id2f = movmean(id2(:), 1);
    did = diff(id1); did = 2*movmean([did(1);did], 1); % two zero cross overs per period
    x_patch_2 = [om(id1); flip(om(id1))];
    upper = (id2f-id1f)./(did-100)*360; % uncertainty in the zero cross over point (0.01 s sampling)
    lower = (id2f-id1f)./(did+100)*360;
    for kk = 1:numel(upper)
        if kk ~=1 && kk ~= numel(upper)
            temp_u(kk) = max([upper(kk-1), upper(kk), upper(kk+1), lower(kk-1), lower(kk), lower(kk+1)]);
            temp_l(kk) = min([upper(kk-1), upper(kk), upper(kk+1), lower(kk-1), lower(kk), lower(kk+1)]);
        else
            temp_u(kk) = upper(kk);
            temp_l(kk) = lower(kk);
        end
    end
    upper = temp_u';
    lower = temp_l';
    y_patch_2 = [upper; flip(lower)];
    
    if ii == 1
        semilogx([1 5], [-180 -180])
        hold on
        semilogx(om(id1), movmean((id2f-id1f)./did*360, 15), 'r', 'LineWidth', 2)
        hold on;
        patch(x_patch_2, y_patch_2, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
    elseif ii == 2
        semilogx(om(id1), movmean((id2f-id1f)./did*360, 15), 'g', 'LineWidth', 2)
        patch(x_patch_2, y_patch_2, 'g', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
    else
        semilogx(om(id1), movmean((id2f-id1f)./did*360, 15), 'b', 'LineWidth', 2)
        patch(x_patch_2, y_patch_2, 'b', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
    end
    
end
figure(4)
xlim([0.05 5.5])
ylim([-90 0])
xlabel('Excitation Frequency (Hz)', 'FontName', 'Times', 'FontSize', 18)
ylabel('Phase Lag (deg)', 'FontName', 'Times', 'FontSize', 18)
title("C)", 'FontName', 'times', 'FontSize', 20)
grid on
lgd = legend('', '(60$^\circ$, 45$^\circ$)', '', '(60$^\circ$, 60$^\circ$)', ...
    '', '(60$^\circ$, 75$^\circ$)', '', 'interpreter', 'latex', ...
    'FontName', 'Times', 'FontSize', 18, 'location', 'southwest');
lgd.EdgeColor = 'none'; lgd.Color = 'none'; lgd.Position = [0.1664 0.2529 0.2755 0.3047];


%% Save the data for Grasp Force Analysis
% save('Data/Matched_Data_60x60_1_Hz.mat', 'fin', 'rok');

exportgraphics(fig2, 'Frequency Response vs Object Size.png', 'Resolution', 1000);
exportgraphics(fig3, 'Force Profile.png', 'Resolution', 1000);
exportgraphics(fig4, 'Phase Lag.png', 'Resolution', 1000);

