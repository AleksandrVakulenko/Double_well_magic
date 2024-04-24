
% ! LOAD output_coeffs from model_03.m !
Temp_array = [22 35 45 55 65 75 85 100 115 125 130];

poly_main = @(x, Root_pos_L, Root_pos_R, Scale, Asym, Scale_basic) ...
    Scale_basic*Scale*(1/4*x.^4 - 1/3*(Root_pos_L + Root_pos_R)*x.^3 + 1/2*Root_pos_L*Root_pos_R*x.^2) - Asym*x;

clc

fig = figure('position', [548   196   543   758]);

Scale_basic = 1e-4;
Far_away_right = 80;
alpha_array_full = [0:0.1:13 13:-0.1:-13 -13:0.1:0]; % kV/cm
filter_for_full = 0.11;

shades = linspace(0.0, 1, 11);
colors = zeros(3, 11);
for i = 1:11
    colors(:, i) = [shades(i), 0, shades(11-i+1)];
end

for N = 1:numel(Temp_array)

    Root_pos_L = output_coeffs(N, 1);
    Root_pos_R = output_coeffs(N, 2);
    Scale = output_coeffs(N, 3);
    Asym = output_coeffs(N, 4);

    fun = @(x) poly_main(x, Root_pos_L, Root_pos_R, Scale, Asym, Scale_basic);
    Well_value_R = fun(Root_pos_R);
    Well_value_L = fun(Root_pos_L);
    
    subplot('Position', [0.13 0.58 0.84 0.37])
    hold on
%     cla
    x = -50*1.2:0.05:50*1.2;
    plot(-x, fun(x), 'Color', colors(:, N), 'LineWidth', 1.1)
    xlabel('P, uC/cm^2')
    ylabel('F, mJ/cm^3')
    xline(0)
    yline(0)
    xlim([-35 35])
    ylim([-35 20])
    title('Double-well energy landscape')
    set(gca, 'fontsize', 13)
    box('on')

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


    subplot('Position', [0.13 0.08 0.84 0.37])
    hold on
%     cla
    plot(alpha_array, -x_min_out, '-', 'Color', colors(:, N), 'LineWidth', 0.8)
%     plot(Einit_p, Pinit_p, '.r', 'linewidth', 1)
%     plot(Einit_n, Pinit_n, '.r', 'linewidth', 1)
    p = xline(0); set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    p = yline(0); set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    xlabel('E, kV/cm')
    ylabel('P, uC/cm^2')
    ylim([-40 40])
    title('FE hysteresis loops')
    set(gca, 'fontsize', 13)
    box('on')
    legend({'22 °C', '35 °C', '45 °C', '55 °C', '65 °C', ...
        '75 °C', '85 °C', '100 °C', '115 °C', '125 °C', '130 °C'}, 'Location', 'best')
    %     xline(Ec_p, '--')
    %     xline(Ec_n)
    %     yline(Pr_p, '--')
    %     yline(Pr_n)
    %     ylim([-1.5 1.5])
    %     xlim([-0.6 0.6])
    drawnow

end