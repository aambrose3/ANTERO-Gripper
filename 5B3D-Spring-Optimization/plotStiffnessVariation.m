clear all; close all; clc;
load('OptimalSpring.mat');
f_s = linspace(-0.001, -150, 51);
figure
for kk = 1:2
    try
        clear param
    end
    load('OptimalSpring.mat');
    d = [param.d-1.2700e-04; param.d+1.2700e-04];
    if kk == 2 % variation from model uncertainty
        param.L_n = param.L_n - 2*0.0719*25.4/1000;
        param.ecc = param.ecc + 0.3041*25.4/1000;
        flag = 1;
    else 
        flag = 0;
    end
    for jj = 1:2 % variation from manufacturing
        param.d = d(jj); 
        param.d_min = d(jj);
        for ii = 1:length(f_s)
            dL(ii, jj, kk) = returnCompression(f_s(ii), param, flag);
        end
        p(:, jj, kk) = polyfit(dL(:, jj, kk), -f_s, 2);
    end
end
x = linspace(0, max(dL, [], 'all'), length(f_s))';
Fp = p(1, :, :).*x.^2 + p(2, :, :).*x;
F_avg = mean(Fp, [2, 3]);
Kp = Fp./x;
K_avg = mean(Kp, [2, 3]);
Fp_end = reshape(Fp(end, :, :), [2, 2]); % jj is row, kk is column
[jj, kk] = find(Fp_end == min(Fp_end, [], 'all'));
[JJ, KK] = find(Fp_end == max(Fp_end, [], 'all'));

figure
hold on
plot(x*1000, Fp(:, jj, kk), '--k', 'Linewidth', 1)
plot(x*1000, Fp(:, JJ, KK), '--k', 'Linewidth', 1)
plot(x*1000, F_avg, 'k', 'LineWidth', 2)
title('Beam Spring Force Uncertainty')
xlabel('Beam Spring Deflection (mm)')
ylabel('Beam Spring Compression Force (N)')

figure
hold on
plot(x*1000, Kp(:, jj, kk)/1000, '--k', 'Linewidth', 1)
plot(x*1000, Kp(:, JJ, KK)/1000, '--k', 'Linewidth', 1)
plot(x*1000, K_avg/1000, 'k', 'LineWidth', 2)
title('Beam Spring Stiffness Uncertainty')
xlabel('Beam Spring Deflection (mm)')
ylabel('Beam Spring Stiffness (N/mm)')
