
poly_main = @(x, Root_pos_L, Root_pos_R, Scale, Asym, Scale_basic) ...
    Scale_basic*Scale*(1/4*x.^4 - 1/3*(Root_pos_L + Root_pos_R)*x.^3 + 1/2*Root_pos_L*Root_pos_R*x.^2) - Asym*x;

clc

fig = figure('position', [548 117 680 838]);
subplot('position', [0.1300 0.6635 0.7750 0.3031])
subplot('position', [0.1300 0.2243 0.7750 0.3574])

slider_pos_l = uicontrol('style','slider', 'InnerPosition', [20 100 300 20]);
uicontrol('style','text', 'string', 'Root_pos_L', 'InnerPosition', [20 120 300 20]);
slider_pos_l.Min = -30;
slider_pos_l.Max = -0.1;
slider_pos_l.Value = -20;

slider_pos_r = uicontrol('style','slider', 'InnerPosition', [340 100 300 20]);
uicontrol('style','text', 'string', 'Root_pos_R', 'InnerPosition', [340 120 300 20]);
slider_pos_r.Min = 0.1;
slider_pos_r.Max = 30;
slider_pos_r.Value = 20;

slider_sclae = uicontrol('style','slider', 'InnerPosition', [340 50 300 20]);
uicontrol('style','text', 'string', 'Scale', 'InnerPosition', [340 70 300 20]);
slider_sclae.Min = log10(2);
slider_sclae.Max = log10(30);
slider_sclae.Value = log10(5);

slider_asym = uicontrol('style','slider', 'InnerPosition', [20 50 300 20]);
uicontrol('style','text', 'string', 'Asym', 'InnerPosition', [20 70 300 20]);
slider_asym.Min = -0.8;
slider_asym.Max = 0.8;
slider_asym.Value = 0;

stop_button = uicontrol('style', 'togglebutton', 'String', 'Stop', 'InnerPosition', [20 20 60 20]);
grid_button = uicontrol('style', 'togglebutton', 'String', 'Full_grid', 'InnerPosition', [300 20 60 20]);


Scale_basic = 1e-4;
Far_away_right = 80;
alpha_array_full = [0:0.1:13 13:-0.1:-13 -13:0.1:0]; % kV/cm
alpha_array_lite = [ 0:0.7:2.1      2.6:0.5:7.1     7.75:1:13.75 ...
               13.0:-1:-1.75  -1.7:-0.5:-7.1   -7.75:-1:-13 ...
               -12:1:-7.75    -7.75:0.5:-2     -1.5:0.5:0 ];
filter_for_full = 0.11;
filter_for_lite = 0.4;

full_grid = 0;
stop = 0;
while ~stop
    Root_pos_L = get(slider_pos_l, 'value');
    Root_pos_R = get(slider_pos_r, 'value');
    Scale = 10^get(slider_sclae, 'value');
    Asym = get(slider_asym, 'value');
    stop = get(stop_button, 'value');
    full_grid = get(grid_button, 'value');

    fun = @(x) poly_main(x, Root_pos_L, Root_pos_R, Scale, Asym, Scale_basic);
    Well_value_R = fun(Root_pos_R);
    Well_value_L = fun(Root_pos_L);

    subplot('position', [0.1300 0.6635 0.7750 0.3031])
    cla
    x = -50*1.2:0.05:50*1.2;
    plot(x, fun(x))
    xlabel('P, uC/cm^2')
    ylabel('F, mJ/cm^3')
    xline(0)
    yline(0)
    if Well_value_L < 0 || Well_value_R < 0
        Ylim = abs(min([Well_value_L Well_value_R]));
        ylim([-Ylim*1.5 Ylim*.5])
    else
        ylim([-1 1])
    end
    drawnow

    if full_grid
        alpha_array = alpha_array_full;
        filter = filter_for_full;
    else
        alpha_array = alpha_array_lite;
        filter = filter_for_lite;
    end
    
    k = 0;
    x_min_out = [];
    y_min_out = [];
    for i = 1:numel(alpha_array)
        alpha = alpha_array(i);
        k = k + 1;

        fun = @(x) poly_main(x, Root_pos_L, Root_pos_R, Scale, Asym, Scale_basic) + alpha.*x;

        oprions = optimoptions('fminunc', 'Display', 'none');
        if i == 1
            [x_min, y_min] = fminunc(fun, Far_away_right, oprions);
            x_min_out(k) = x_min;
            y_min_out(k) = y_min;
        else
            [x_min, y_min] = fminunc(fun, x_min_out(k-1), oprions);
            x_min_out(k) = x_min*filter + x_min_out(k-1)*(1-filter);
            y_min_out(k) = y_min;
        end
    end


    subplot('position', [0.1300 0.2243 0.7750 0.3574])
    hold on
    cla
    plot(alpha_array, -x_min_out, '.-')
    plot(Einit_p, Pinit_p, '.r', 'linewidth', 1)
    plot(Einit_n, Pinit_n, '.r', 'linewidth', 1)
    xline(0)
    yline(0)
    xlabel('E, kV/cm')
    ylabel('P, uC/cm^2')
    ylim([min(Pinit_n)*1.2 max(Pinit_p)*1.2])
    title(['L: ' num2str(Root_pos_L) ' | R: ' num2str(Root_pos_R) ...
        ' | Scale: ' num2str(Scale) ' | Asym: ' num2str(Asym)])
%     xline(Ec_p, '--')
%     xline(Ec_n)
%     yline(Pr_p, '--')
%     yline(Pr_n)
%     ylim([-1.5 1.5])
%     xlim([-0.6 0.6])
    drawnow


    
end

close(fig)


disp([num2str(Root_pos_L) ' ' num2str(Root_pos_R) ' ' num2str(Scale) ' ' num2str(Asym)])


% Root_pos_L Root_pos_R Scale Asym   
output_coeffs = [-21.3053 20.2545 5.7016 0.31168; % 1
                 -20.4291 20.5633 5.4622 0.20178; % 2
                 -19.4743 20.4574 5.7315 0.10986; % 3
                 -18.2413 19.4903 5.935 0.070221; % 4
                 -16.9459 20.1765 5.7463 -0.17018; % 5
                 -16.4676 18.8617 5.9815 -0.097485; % 6
                 -16.7573 18.0863 6.7896 0.043822; % 7
                 -11.8182 13.6557 12.1856 0.10677; % 8
                 -10.5445 10.3912 15.264 0.15938; % 9
                 -8.0116 5.5295 17.5065 0.057253; % 10
                 -5.5364 0.1 20.6827 0.44101; % 11
                 ];





