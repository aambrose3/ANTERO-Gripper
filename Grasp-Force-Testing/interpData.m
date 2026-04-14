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
%  Matched_Data.mat -> Strucutre of trials with the finger data and the 
%  6-axis force and torque readings from the F/T sensor with time stamps.

clear all; close all; clc;

load('Data/Finger_Data_60x45.mat'); % load the unpacked teensy data. produces a struct called 'finger'

load('Data/Rokubi_Data_60x45.mat'); % load the unpacked rokubi. produces a struct called 'rokubi'

fieldNames = fieldnames(finger);
kk = 1;
for ii = 1:numel(fieldNames)
    if ii ~= 7 && ii ~= 8 && ii ~= 9 && ii ~= 10
        tic;
        Name = fieldNames{ii};
        fin{ii} = finger.(Name);
        rok{ii} = rokubi.(Name);
        mean(rok{ii}.temperature_C);
        std(rok{ii}.temperature_C);
        ff = 1/round(mean(diff(fin{ii}.time)), 3);
        fr = 1/round(mean(diff(rok{ii}.device_timestamp/1000000)), 3);
        num = ceil(ff/2); % 0.5 second moving average filter
        coeff = ones(1, num)/num;
        fin{ii}.F = filter(coeff, 1, fin{ii}.F);
        num = ceil(fr/2); % 0.5 second moving average filter
        coeff = ones(1, num)/num;
        rok{ii}.Fx = filter(coeff, 1, rok{ii}.Fx);
        rok{ii}.Fy = filter(coeff, 1, rok{ii}.Fy);
        rok{ii}.Fz = filter(coeff, 1, rok{ii}.Fz);
        rok{ii} = downsample(rok{ii}, fr/ff);
        % offset = 105;
        % stop = 415;
        % idf = find(fin{ii}.F(offset:stop) == max(fin{ii}.F(offset:stop)), 1, 'first')+offset;
        % fin{ii} = fin{ii}(idf:end, :); fin{ii}.time = fin{ii}.time - fin{ii}.time(1);
        offset = find(fin{ii}.F == max(fin{ii}.F), 1, 'last') - 5;
        temp = diff(fin{ii}.F(offset:offset+100));
        idf = find(temp < -0.1, 1, 'first')+offset-1; % first decrease from peak force
        fin{ii} = fin{ii}(201:900, :); fin{ii}.time = fin{ii}.time - fin{ii}.time(1);
    
        % stop = 543;
        % idr = find(rok{ii}.Fx(1:stop) == max(rok{ii}.Fx(1:stop)), 1, 'first');
        % rok{ii} = rok{ii}(idr:end, :); 
        offset = find(rok{ii}.Fy == max(rok{ii}.Fy), 1, 'last') - 5;
        temp = diff(rok{ii}.Fy(offset:offset+100));
        if find(temp == min(temp), 1, 'first') > 30
            offset = offset+find(temp == min(temp), 1, 'first')-1;
            temp = diff(rok{ii}.Fy(offset:offset+100));
            idr = find(temp < -0.06, 1, 'first')+offset-1; % first decrease from peak force
        else
            idr = find(temp < -0.06, 1, 'first')+offset-1; % first decrease from peak force
        end
        offset = idr-idf+200;
        temp = diff(rok{ii}.Fy);
        offset = find(temp >0.1, 1, 'first')+80;
        id1 = find(diff(rok{ii}.Fy(offset:offset+200)) > 0.1, 1, 'first')+offset-1;
        offset = offset+650;
        id2 = find(diff(rok{ii}.Fy(offset:offset+100)) < -0.1, 1, 'first')+offset-1;
        rok{ii} = rok{ii}(id1:id2, :); % still no precise
    
        rok{ii}.device_timestamp = (rok{ii}.device_timestamp - rok{ii}.device_timestamp(1))/1E6;
        rok{ii} = removevars(rok{ii}, ["host_time_s", "status"]);
        rok{ii}.Properties.VariableNames{1} = 'time';

        
        % dfx = diff(rok{ii}.Fx);
        % offset = 3600;
        % [npks, nlocs] = findpeaks(-dfx(offset:end), ...
              % 'MinPeakDistance', 300, 'MinPeakHeight', 0.01);
        % ix(ii) = nlocs(end)+offset;
        % idx(ii) = find(abs(dfx(1:ix(ii))) < 0.001, 1, 'last');
        % rok{ii} = rok{ii}(1:idx(ii), :);
        Temperature_C = NaN*ones(size(fin{ii}.time));
        tmp = fin{ii}; tmp = addvars(tmp, Temperature_C);
        for jj = 1:numel(tmp.Properties.VariableNames)
            xq = linspace(0, fin{ii}.time(end), size(fin{ii}, 1))';
            x = linspace(0, fin{ii}.time(end), size(rok{ii}, 1))';
            tmp{:, jj} = interp1(x, rok{ii}{:, jj}, xq);
        end
        tmp.Properties.VariableNames{1} = 'time';
        tmp.Properties.VariableNames{2} = 'a';
        tmp.Properties.VariableNames{3} = 'b';
        tmp.Properties.VariableNames{4} = 'cc';
        tmp.Properties.VariableNames{5} = 'd';
        tmp.Properties.VariableNames{6} = 'e';
        tmp.Properties.VariableNames{7} = 'f';
        tmp.Properties.VariableNames{2} = 'Fz'; tmp.Fz = abs(tmp.Fz);
        tmp.Properties.VariableNames{3} = 'Fx'; tmp.Fx = abs(tmp.Fx);
        tmp.Properties.VariableNames{4} = 'Fy'; tmp.Fy = abs(tmp.Fy);
        tmp.Properties.VariableNames{5} = 'Mz'; tmp.Mz = abs(tmp.Mz);
        tmp.Properties.VariableNames{6} = 'Mx'; tmp.Mx = abs(tmp.Mx);
        tmp.Properties.VariableNames{7} = 'My'; tmp.My = abs(tmp.My);
        tmp.Properties.VariableNames{8} = 'Temperature_C';
        rok{ii} = tmp;
        rok{ii}.time = fin{ii}.time;
        % for jj = 1:7 %% exclude the transient data
        %     if jj == 1
        %         finTemp = table2array(fin{ii}(30+(100*(jj-1)):90+(100*(jj-1)), :));
        %         rokTemp = table2array(rok{ii}(30+(100*(jj-1)):90+(100*(jj-1)), :));
        %     else
        %         finTemp = vertcat(finTemp, table2array(fin{ii}(30+(100*(jj-1)):90+(100*(jj-1)), :)));
        %         rokTemp = vertcat(rokTemp, table2array(rok{ii}(30+(100*(jj-1)):90+(100*(jj-1)), :)));
        %     end
        % end
        % fin{ii} = array2table(finTemp, 'VariableNames', ["time", "t1", "t2", "t3", "c", "F", "Ref"])
        % rok{ii} = array2table(rokTemp, 'VariableNames', ["time", "Fz", "Fx", "Fy", "Mz", "Mx", "My", "Temperature_C"])
        clear tmp;
        % plot(fin{ii}.time, fin{ii}.F, 'k', 'LineWidth', 2)
        % hold on
        % plot(rok{ii}.time, sqrt(rok{ii}.Fx.^2 + rok{ii}.Fy.^2 + rok{ii}.Fz.^2), 'r', 'LineWidth', 2)
        % xlabel('Time (s)', 'FontName', 'Times', 'FontSize', 12)
        % ylabel('Eucliden Norm of Grasp Force', 'FontName', 'Times', 'FontSize', 12)
        % legend('Estimated Force', 'Measured Force')
        % ylim([0 30])
        % hold off
        % fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
    end
end
%% Save the data for Grasp Force Analysis
% save('Data/Matched_Data_60x75.mat', 'fin', 'rok');

%% plotting the force profile
for ii = 1:length(fin)
    if ii == 1
        time = fin{ii}.time;
        f = fin{ii}.F;
       r = sqrt(rok{ii}.Fx.^2+rok{ii}.Fy.^2+rok{ii}.Fz.^2);
    else
        time = horzcat(time, fin{ii}.time);
        f = horzcat(f, fin{ii}.F);
        r = horzcat(r, sqrt(rok{ii}.Fx.^2+rok{ii}.Fy.^2+rok{ii}.Fz.^2));
    end
end
time = mean(time, 2, 'omitnan');
f_std = std(f, [], 2, 'omitnan');
f = mean(f, 2, 'omitnan');
r_std = std(r, [], 2, 'omitnan');
r = mean(r, 2, 'omitnan');
fig = figure(1)
p1 = plot(time, f, 'k', 'LineWidth', 2);
hold on
p2 = plot(time, r, '--k', 'LineWidth', 1);
xlabel('Time (s)', 'FontName', 'Times', 'FontSize', 14)
ylabel('Grasp Force Magnitude (N)', 'FontName', 'Times', 'FontSize', 14)
F_max = (r+r_std)./(f-f_std).*f;
F_min = (r-r_std)./(f+f_std).*f;
x_patch = [time; flip(time)];
y_patch = [F_max; flip(F_min)];
patch(x_patch, y_patch, 'k', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
legend('Estimated Force', 'Measured Force', '$\pm \ 1\sigma$', ...
    'interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12, 'location', ...
    'northwest', 'EdgeColor', 'none')
axis([0 70 4 26])

exportgraphics(fig, 'Grasp Force Profile.png', 'Resolution', 1000);