
% addpath('./include')


% % ------ PMN-33PT fall 2023 100um ------
folder = 'Results 2023/Results_2023_10_11_PMN_33PT/';
Sample.h = 100e-6; % m
Sample.s = 0.003*0.003; % m^2

names = dir(folder);
names = {names.name};
names(1:2) = [];
names = string(names)';

% fig = figure('position', [416   272   785   739]);
figure

k = 0;
for i = 11%[1:11 13]%numel(names)%3:4:720 %[3:4:720]

load([folder char(names(i))])
k = k + 1;

temp = Loops.temp; % K

feloop = Loops.feloop;

% corrected = feloop_processing(feloop, Sample, fig); %FIXME: WRONG SAMPLE SIZE

% cla
hold on
LW_init = 0.8;
LW_ref = 1;

Einit = feloop.init.E.p;
Pinit = feloop.init.P.p;
Einit_p = Einit/1000/(Sample.h*100);
Pinit_p = Pinit*1e6/(Sample.s*100^2);

Einit = feloop.init.E.n;
Pinit = feloop.init.P.n;
Einit_n = Einit/1000/(Sample.h*100);
Pinit_n = Pinit*1e6/(Sample.s*100^2);

Shift_p = (Pinit_p(end)-Pinit_p(1))/2;
Pinit_p = Pinit_p - Shift_p;

Shift_n = (Pinit_n(end) - Pinit_n(1))/2;
Pinit_n = Pinit_n - Shift_n;

plot(Einit_p, Pinit_p, '-b', 'linewidth', LW_init)
plot(Einit_n, Pinit_n, '-b', 'linewidth', LW_init)
xline(0)
yline(0)

ylim([-40 40])

title([num2str(i) ' T = ' num2str(temp) ' K'])
% xlim([-50 50])
drawnow
pause(0.5)
end