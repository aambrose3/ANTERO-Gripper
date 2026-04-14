%% 5B3D-Finger-Model-Results: Optimization
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% This script helps analyze the results from the 5B3D-Model: Main

clear all; close all; clc;
addpath("../");
load('OptimizationFunctions.mat');

n1 = numel(t_1);
n2 = numel(t_2);
n3 = numel(t_3);

x = 46.018/1000;
y = 0.118342;
F1 = num2cell(nan(n1, n2, n3)); F2 = F1;
count = 0;
for ii = 1:n1
    for jj = 1:n2
        for kk = 1:n3
            % update parameters for this configuration
            localParam = updateParam(param, t_1(ii), t_2(jj), t_3(kk));
            % compute spring deflection & force
            [cmp, F_cmp, TA, flag] = getSpring(localParam, t_1(ii), t_2(jj), t_3(kk), true);
            solve(cmp - 25 < 0, [param.d, param.l1_r], 'ReturnConditions', true)
            cmp = matlabFunction(cmp);
            TA = matlabFunction(TA);
            if cmp(x, y) < 25 % limit spring compression
                if TA(x, y) < param.TA_max && TA(x, y) > 0 % limit actuator effort
                    if F1{ii, jj, kk} > 0 && F2{ii, jj, kk} > 0 % limit contact forces to be positive
                        F1{ii, jj, kk} = F1{ii, jj, kk}(x, y);
                        F2{ii, jj, kk} = F2{ii, jj, kk}(x, y);
                    end % else leave as NaN
                end % else leave as NaN
            end % else leave as NaN
            count=count+1;
            fprintf('%d\n',count)
        end
    end
end
