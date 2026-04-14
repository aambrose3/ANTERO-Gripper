%% Grasp_Force_Testing: unpackRokubiData
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

% Description:
%  This script unpacks the data streams from the 6-axis BOTA Systems Rokubi 
%  Serial Force-Torque (F/T) sensor. The sensor communicates with Python on a
%  host PC and this data is stored in a csv file with time stamps
%  
%  Outputs:
%  Rokubi_Data.mat -> Strucutre of trials with the 6-axis force and torque
%  readings from the F/T sensor with time stamps.

clear all; close all; clc;

% Import Rokubi Data
Data.t1 = readtable('60x75\5_HZ\T1.csv', 'FileType', 'text');
Data.t2 = readtable('60x75\5_HZ\T2.csv', 'FileType', 'text');
Data.t3 = readtable('60x75\5_HZ\T3.csv', 'FileType', 'text');
Data.t4 = readtable('60x75\5_HZ\T4.csv', 'FileType', 'text');
Data.t5 = readtable('60x75\5_HZ\T5.csv', 'FileType', 'text');
% Data.t6 = readtable('60x75\test.csv', 'FileType', 'text');



fieldNames = fieldnames(Data);
kk = 1;
for ii = 1:numel(fieldNames)
    tic;
    Name = fieldNames{ii};
    data = Data.(Name);
    data = removevars(data, {'host_time_rel_s', 'N1', 'N2', 't3'});
    data = renamevars(data, {'host_time_s', 'reportForce', 'reportFlag', 'E1'}, {'Time', 'F', 'Flag', 'E'});
    % zero time
    data.Time = data.Time - data.Time(1);
    % data.host_time_rel_s = data.host_time_rel_s - data.host_time_rel_s(1);
    ANTERO.(Name) = data;

    clear x data
    fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
end
% save the assembled structure to a mat file
save('ANTERO_Data_60x75_5_Hz.mat', 'ANTERO')