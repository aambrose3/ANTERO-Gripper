%% 5B3D-Finger-Model-Results: genVideo
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
%  Generate an .mp4 video of the stable grasp configurations

%  Inputs:
%  param: basic stucture of fixed finger parameters
%  t_1 -> a vector of possible \theta_1 anlges (rad)
%  t_2 -> a vector of possible \theta_2 angles (rad)
%  t_3 -> a vector of possible crank/actaotr angles/\theta_3 angles (rad)
%  rate -> the desired fram rate of the video (fps) - determines length as well 
%  name -> the desired name (with file extensions) of the created video

function genVideo(param, N, t_1, t_2, t_3, rate, name)
    n1 = numel(t_1);
    d1 = 5*round(1/(rad2deg(t_1(2)) - rad2deg(t_1(1))));
    d1 = 1;
    n2 = numel(t_2);
    d2 = 5*round(1/(rad2deg(t_2(2)) - rad2deg(t_2(1))));
    d2 = 1;
    n3 = numel(t_3);
    d3 = 2*round(1/(rad2deg(t_3(2)) - rad2deg(t_3(1))));
    % d3 = 1;

    video_fileName = name;
    v = VideoWriter(video_fileName, 'MPEG-4');
    v.FrameRate = rate;
    open(v);
    figure(1)
    for ii = 1:d1:n1
        r12 = [param.a*cos(t_1(ii)); param.a*sin(t_1(ii))];
        for jj = 1:d2:n2
            r2_5 = r12 + [param.a*cos(t_1(ii)+t_2(jj)); ...
                param.a*sin(t_1(ii)+t_2(jj))];
            r23 = r12 + [param.b*cos(t_1(ii)+t_2(jj)-param.alpha); ...
                param.b*sin(t_1(ii)+t_2(jj)-param.alpha)];
            tmp = isnan(squeeze(N(ii, jj, :, 1)));
            if tmp(1) == 1
                start = find(tmp == 0, 1, 'first');
                stop = find(tmp == 0, 1, 'last');
                if isempty(start) == 1
                    start = n3;
                end
                if isempty(stop) == 1
                    stop = n3;
                end
            else
                start = 1;
                stop = n3;
            end
            for kk = start:d3:stop
    
                % if any(isnan(N(ii, jj, kk, :))) == 1 && kk > d3
                %     if any(isnan(N(ii, jj, kk-d3, :))) == 0
                %         break;
                %     end
                % end
                % Update Configuration
                lParam = updateParam(param, t_1(ii), t_2(jj), t_3(kk));
                if lParam.c < lParam.l1_o 
                    break;
                end
                r15 = param.e*[cos(-param.gamma);sin(-param.gamma)];
                r54 = r15 + param.d*[cos(t_3(kk));sin(t_3(kk))];
                rp1 = lParam.p1*[cos(t_1(ii)); sin(t_1(ii))] + ...
                    param.t*[cos(t_1(ii)+pi/2); sin(t_1(ii)+pi/2)];
                rp2 = r12 + lParam.p2*[cos(t_1(ii)+t_2(jj)); sin(t_1(ii)+t_2(jj))] + ...
                    param.t*[cos(t_1(ii)+t_2(jj)+pi/2); sin(t_1(ii)+t_2(jj)+pi/2)];
                %% Plot Finger
                % Plot Palm
                plot([0 0 param.O(1)], [0 param.O(2) param.O(2)], ...
                    'k', 'LineWidth', 4)
                hold on
                title('Finger Display')
                xlabel('Distance (mm)')
                ylabel('Distance (mm)')
                % Plot link a
                plot([0, r12(1)], [0, r12(2)], ...
                    'b', 'LineWidth', 4)
                % Plot link b
                plot([r12(1), r2_5(1)], [r12(2), r2_5(2)], ...
                    'g', 'LineWidth', 4)
                plot([r12(1), r23(1)], [r12(2), r23(2)], ...
                    'g', 'LineWidth', 4)
                % Plot link D
                plot([0, r15(1)], [0, r15(2)], ...
                    'r', 'LineWidth', 4)
                % Plot Link E
                plot([r15(1), r54(1)], [r15(2), r54(2)], ...
                    'y', 'LineWidth', 4)
                % Plot Link C
                plot([r54(1) r23(1)], [r54(2) r23(2)], ...
                    '--m', 'LineWidth', 4)
                % Plot Contact Points
                plot(rp1(1), rp1(2), '.k', 'MarkerSize', 10)
                plot(rp2(1), rp2(2), '.k', 'MarkerSize', 10)
                % Plot Contact Forces
                tmp = lParam.n1*N(ii, jj, kk, 1)/100;
                quiver(rp1(1)-tmp(1), rp1(2)-tmp(2), tmp(1), tmp(2), 'k', 'LineWidth', 2)
                tmp = lParam.n2*N(ii, jj, kk, 2)/100;
                quiver(rp2(1)-tmp(1), rp2(2)-tmp(2), tmp(1), tmp(2), 'k', 'LineWidth', 2)
                % Plot pin joints
                plot(0, 0, '.k', 'MarkerSize', 40)
                plot(r12(1), r12(2), '.k', 'MarkerSize', 40)
                plot(r23(1), r23(2), '.k', 'MarkerSize', 40)
                plot(r54(1), r54(2), '.k', 'MarkerSize', 40)
                plot(r15(1), r15(2), '.k', 'MarkerSize', 40)
                xlim([-0.1 0.25])
                ylim([-0.1 0.25])
                % axis equal
                drawnow
                pause(0.01);
                frame = getframe(gcf);
                writeVideo(v, frame);
                hold off
            end
        end
    end
    close(v)
    disp(['Video Saved to ' video_fileName]);
    close all

end