%% Grasp_Force_Testing: unpackRokubiData
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script unpacks the data streams from the 6-axis BOTA Systems Rokubi 
%  Serial Force-Torque (F/T) sensor. The sensor communicates with Python on a
%  host PC and this data is stored in a csv file with time stamps
%  
%  Outputs:
%  Rokubi_Data.mat -> Strucutre of trials with the 6-axis force and torque
%  readings from the F/T sensor with time stamps.

clear all; close all; clc;

%% Import Rokubi Data
Data.t1 = readtable('60x75\5_HZ\T1-FT.csv', 'FileType', 'text');
Data.t2 = readtable('60x75\5_HZ\T2-FT.csv', 'FileType', 'text');
Data.t3 = readtable('60x75\5_HZ\T3-FT.csv', 'FileType', 'text');
Data.t4 = readtable('60x75\5_HZ\T4-FT.csv', 'FileType', 'text');
Data.t5 = readtable('60x75\5_HZ\T5-FT.csv', 'FileType', 'text');
% Data.t6 = readtable('60x60\5_HZ\T6-FT.csv', 'FileType', 'text');


fieldNames = fieldnames(Data);
kk = 1;
% figure
for ii = 1:numel(fieldNames)
    tic;
    Name = fieldNames{ii};
    data = Data.(Name);
    data = removevars(data, {'host_time_s', 'status', 'uart_age_s', 'uart_host_time_s', 'uart_line'});
    data = renamevars(data, {'device_timestamp', 'temperature_C'}, {'Time', 'Temp'});
    F = NaN*ones(size(data.Fx));
    data = addvars(data, F, 'Before', 'Fx');
    % Zero the sensor out
    Zfx = mean(data.Fx(1:2000));
    Zfy = mean(data.Fy(1:2000));
    Zfz = mean(data.Fz(1:2000));
    Zmx = mean(data.Mx(1:2000));
    Zmy = mean(data.My(1:2000));
    Zmz = mean(data.Mz(1:2000));
    % assemble the structure
    data.Fx = abs(data.Fx - Zfx);
    data.Fy = abs(data.Fy - Zfy);
    data.Fz = abs(data.Fz - Zfz);
    data.Mx = data.Mx - Zmx;
    data.My = data.My - Zmy;
    data.Mz = data.Mz - Zmz;
    data.F = sqrt(data.Fx.^2+data.Fy.^2+data.Fz.^2);
    data.F = data.F - mean(data.F(1:2000));
    % Zero time
    data.Time = (data.Time - data.Time(1))/1E6;   
    % hold on
    % plot(data.Time, data.F)
    if ii == 1
        data.Time(9145) = (data.Time(9146) - data.Time(9144))/2+data.Time(9144);
    end
    rokubi.(Name) = data;
    clear x data
    fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
end
% save the assembled structure to a mat file
save('Rokubi_Data_60x75_5_Hz.mat', 'rokubi')