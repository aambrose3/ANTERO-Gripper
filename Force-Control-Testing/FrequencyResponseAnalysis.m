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

% Reiterate the Reference force signal
tf = 300;
rate = tf.^4/pi*0.35;
tt = 0:0.01:300; tt = tt(:);
omega = tt.^4/rate;
offset = 10;
R = 2*cos(omega.*tt);
% eliminate the first 20 seconds of data
tt = tt(end-26001:end); tt = tt - tt(1);
omega = omega(end-26001:end);
R = R(end-26001:end); dRdt = [diff(R)./diff(tt); NaN];


for jj = 1:3 % tested at three finger positions
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
    fig.Position = [100 500 550 400];
    for ii = 1:numel(fieldNames) % multiple trials per finger position
        % if jj ~= 2 || ii ~= 3 % omit this data -> outlier by vosual inspection
            Name = fieldNames{ii};
            ant{ii} = ANTERO.(Name);
            % eliminate useless data
            ant{ii} = ant{ii}(end-26001:end, :); ant{ii}.Time = ant{ii}.Time - ant{ii}.Time(1);
            rf = 1./floor(ant{ii}.Time(end)); % minimum resolvable frequency
            % Remove the 0 frequency offset
            ant{ii}.F = ant{ii}.F; ant{ii}.Ref = ant{ii}.Ref;

            % 1. Get Reference Indices
            [~, L_r_p] = findpeaks(ant{ii}.Ref);     % Peaks of Input
            [~, L_r_n] = findpeaks(-ant{ii}.Ref);    % Troughs of Input
            
            % 2. Pre-allocate Output indices
            L_f_p = zeros(size(L_r_p));
            L_f_n = zeros(size(L_r_n));
            
            % 3. Find Output Peaks/Troughs using Search Windows (No '==' or '-5')
            win = 10; % Search +/- 10 samples around the reference peak/trough
            for k = 1:numel(L_r_p)
                range = max(1, L_r_p(k)-win) : min(length(ant{ii}.F), L_r_p(k)+win);
                searchWin = 1./min(omega(range))/2; % search length for peak in time (up to a 180 deg phase lag)
                dt_est = (max(ant{ii}.Time(range)) - min(ant{ii}.Time(range)))/(numel(range)-1);
                searchWin = L_r_p(k):(L_r_p(k)+ceil(searchWin/dt_est));
                [~, relative_idx] = max(ant{ii}.F(searchWin));
                L_f_p(k) = L_r_p(k) + relative_idx - 1;
                %% Visualize Search
                % findpeaks(ant{ii}.Ref)
                % xlim([min(searchWin)-numel(searchWin) max(searchWin)+numel(searchWin)])
                % hold on
                % plot(searchWin, ant{ii}.F(searchWin))
                % scatter(L_f_p(k), ant{ii}.F(L_f_p(k)))
                % hold off
            end
            for k = 1:numel(L_r_n)
                range = max(1, L_r_n(k)-win) : min(length(ant{ii}.F), L_r_n(k)+win);
                searchWin = 1./min(omega(range))/2; % search length for peak in time (up to a 180 deg phase lag)
                dt_est = (max(ant{ii}.Time(range)) - min(ant{ii}.Time(range)))/(numel(range)-1);
                searchWin = L_r_n(k):(L_r_n(k)+ceil(searchWin/dt_est));
                if max(searchWin) > numel(ant{ii}.F)
                    searchWin = searchWin(1):numel(ant{ii}.F);
                end
                [~, relative_idx] = min(ant{ii}.F(searchWin));
                L_f_n(k) = L_r_n(k) + relative_idx - 1;
                %% Visualize Search
                % findpeaks(ant{ii}.Ref)
                % xlim([min(searchWin)-numel(searchWin) max(searchWin)+numel(searchWin)])
                % hold on
                % plot(searchWin, ant{ii}.F(searchWin))
                % scatter(L_f_p(k), ant{ii}.F(L_f_p(k)))
                % hold off
            end
            
            % 4. Pair them up correctly for "Rising" vs "Falling"
            % Let's assume the signal is Trough -> Peak -> Trough (Rising Segment)
            % We must ensure we have the same number of points for subtraction
            n_cycles = min([length(L_f_p), length(L_f_n)]);
            L_f_p = L_f_p(1:n_cycles);
            L_f_n = L_f_n(1:n_cycles);
            L_r_p = L_r_p(1:n_cycles);
            L_r_n = L_r_n(1:n_cycles);

            % 5. Calculate Dynamic Amplitudes
            % Amplitude = (Peak - Trough) / 2
            A_r{ii} = abs(R(L_r_p) - R(L_r_n)) / 2;
            A_f{ii} = abs(ant{ii}.F(L_f_p) - ant{ii}.F(L_f_n)) / 2;

            num_p = numel(L_r_p);
            freq_rising = zeros(num_p-1, 1);
            gain_rising = zeros(num_p-1, 1);
            
            % 6. Loop through segments to calculate frequency using Time
            for k = 1:num_p
                % Find the trough that sits between these two peaks
                % (We use the reference signal 't' to define the frequency)
                t_peak = ant{ii}.Time(L_r_p(k));
                t_trough = ant{ii}.Time(L_r_n(k)); 
                
                % Duration of this half-cycle (Rising)
                dt_rising = abs(t_peak - t_trough);
                
                % INSTANTANEOUS FREQUENCY (Hz)
                ff{ii}(k) = 1 / (2 * dt_rising);
                                % DYNAMIC GAIN
                % Use the values at the detected peak/trough indices
                amp_in = abs(R(L_r_p(k)) - R(L_r_n(k))) / 2;
                amp_out = abs(ant{ii}.F(L_f_p(k)) - ant{ii}.F(L_f_n(k))) / 2;
                gain_rising(k) = amp_out / amp_in;
            end 
            ff{ii} = ff{ii}(:);
            pw{ii} = 20*log10(A_f{ii}./A_r{ii});
            figure(1)
            semilogx([0.1 1], [-50 -50], 'w')
            hold on
            scatter(ff{ii}, pw{ii}, '.')
            ylim([-15 0])

            
        % end
    end   

    freq{jj} = median(cell2mat(ff), 2);
    mag{jj} = median(cell2mat(pw), 2);
    mag_std{jj} = std(cell2mat(pw), [], 2);
    mag_max{jj} = max(cell2mat(pw), [], 2);
    mag_min{jj} = min(cell2mat(pw), [], 2);
    %% Plot
    figure(1)
    plot(freq{jj}, mag{jj}, 'k', 'LineWidth', 1.5);
    grid on;
    title('Frequency Response of the Froce Control System');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([0.05 6])

    x_patch{jj} = [freq{jj}; flip(freq{jj})];
    upper = movmean(mag{jj}+mag_std{jj}, 10);
    lower = movmean(mag{jj}-mag_std{jj}, 10);
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

    temp = [ant{1}.Time, ant{2}.Time, ant{3}.Time, ant{4}.Time, ant{5}.Time];
    Time{jj} = mean(temp, 2);
    temp = [ant{1}.Ref, ant{2}.Ref, ant{3}.Ref, ant{4}.Ref, ant{5}.Ref];
    Ref{jj} = mean(temp, 2);
    temp = [ant{1}.F, ant{2}.F, ant{3}.F, ant{4}.F, ant{5}.F];
    F{jj} = mean(temp, 2);
    fig5 = figure(5);
    fig5.Position = [100 1000 1100 350]
    subplot(1, 3, jj)
    axis square
    xlim([7.9 12.1])
    C = linspace(0, 1, ceil(numel(ant{ii}.Time)/100))';
    if jj == 1
        C = [C, zeros(size(C)), zeros(size(C))]; % Red
    elseif jj == 2
        C = [zeros(size(C)), C, zeros(size(C))]; % Green
    else
        C = [zeros(size(C)), zeros(size(C)), C]; % Blue
    end
    for kk = 1:floor(numel(ant{ii}.Time)/100) 
        if kk < floor(numel(ant{ii}.Time)/100)-1
            plot(Ref{jj}(100*(kk-1)+1:100*kk), F{jj}(100*(kk-1)+1:100*kk), 'Color', C(kk, :), 'LineWidth', 2)
            hold on
        else
            plot(Ref{jj}(100*kk:end), F{jj}(100*kk:end), 'Color', C(kk, :), 'LineWidth', 2)
        end
    end
    if jj == 1
        ylabel('Measured Force (N)', 'FontName', 'Times', 'FontSize', 18)
        title('60$^\circ$x45$^\circ$', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 18)
        axis([7.5 12.5 7.5 12.5])
    elseif jj == 2
        title('60$^\circ$x60$^\circ$', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 18)
        axis([7.5 12.5 7.5 12.5])
    else
        title('60$^\circ$x75$^\circ$', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 18)
        axis([7.5 12.5 7.5 12.5])
    end
    xlabel('Reference Force (N)', 'FontName', 'Times', 'FontSize', 18)
    
    
        
    % title("A)", 'FontName', 'times', 'FontSize', 20)
end
fig2 = figure(2);
fig2.Position = [650 600 550 300];
%% 60x45

% semilogx([0.01; 6], [-3; -3], '--k', 'LineWidth', 2)
% hold on;
semilogx(freq{1}, movmean(mag{1}, 1), 'r', 'LineWidth', 2);
% s1 = semilogx(freq{1}, movmean(mag_max{1}, 1), '.r', 'LineWidth', 2);
hold on;
patch(x_patch{1}, y_patch{1}, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
%% 60x60
semilogx(freq{2}, movmean(mag{2}, 1), 'g', 'LineWidth', 2);
patch(x_patch{2}, y_patch{2}, 'g', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
%% 60x75
semilogx(freq{2}, movmean(mag{3}, 1), 'b', 'LineWidth', 2);
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
% plot(Time{3}, F{3}, 'k', 'LineWidth', 2)
C = linspace(0, 1, 121)';
C = [zeros(size(C)), zeros(size(C)), C]; % Blue
for kk = 1:121
    if kk < 120
        plot(Time{3}(100*(kk-1)+1:100*kk+2), F{jj}(100*(kk-1)+1:100*kk+2), 'Color', C(kk, :), 'LineWidth', 2)
        hold on
    else
        plot(Time{3}(100*kk:end), F{jj}(100*kk:end), 'Color', C(kk, :), 'LineWidth', 2)
    end
end

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
    time = Time{ii}(idx):0.0001:Time{3}(end); % interpolate to acheive better accuracy on zero cross over point
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
exportgraphics(fig5, 'Input-Output Portrait.png', 'Resolution', 1000);
% 
