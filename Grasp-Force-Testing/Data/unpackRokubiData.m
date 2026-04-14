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
RData.t1 = readtable('Grasp_Force_Testing_Rokubi\60x75\T1.csv', 'FileType', 'text');
RData.t2 = readtable('Grasp_Force_Testing_Rokubi\60x75\T2.csv', 'FileType', 'text');
RData.t3 = readtable('Grasp_Force_Testing_Rokubi\60x75\T3.csv', 'FileType', 'text');
RData.t4 = readtable('Grasp_Force_Testing_Rokubi\60x75\T4.csv', 'FileType', 'text');
RData.t5 = readtable('Grasp_Force_Testing_Rokubi\60x75\T5.csv', 'FileType', 'text');
RData.t6 = readtable('Grasp_Force_Testing_Rokubi\60x75\T6.csv', 'FileType', 'text');


fieldNames = fieldnames(RData);
kk = 1;
for ii = 1:numel(fieldNames)
    tic;
    Name = fieldNames{ii};
    data = RData.(Name);
    if size(data, 1) > 25000
        data = data(1:25000, :)
    end
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
    % determine the average sample rate of the data
    fr = 1/round(mean(diff(data.device_timestamp/1000000)), 3);
    %% filter and crop?
    num = ceil(fr/4); % 0.26 second moving average filter
    df = diff(movmean(data.Fy, num)); % filter the data
    offset = 400; % find when the sensor makes contact with the finger
    idx = find(df > 0.01, 1, 'first');
    idx2 = find(df(idx-100:idx) < 0.001, 1, 'last');
    idx = idx-100+idx2-2; idx2 = idx+200*70;

    rokubi.(Name) = data;
    clear x data
    fprintf('-> Done Unpacking Trial %d in %.2f s\n', ii, toc);
end
% save the assembled structure to a mat file
save('Rokubi_Data_60x75.mat', 'rokubi')