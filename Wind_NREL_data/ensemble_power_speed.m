clear all
data_NREL = load_wind_para('model_data.mat',8);
wt = data_NREL.wt;
env = data_NREL.env;

n = wt.ctrl.gen.effeciency;
Pb = 5e6;
Pb_m = Pb/n;
w_LS = 12.1*2*pi/60; %rad/s from rpm
w_min = 6.9*2*pi/60; %rad/s from rpm

radius = wt.rotor.radius;
rho = env.rho;

path = '';

%% Draw power/speed curves
cp = wt.cp.table;
cp = cp(1,:);
lambda = wt.cp.tsr;
N = length(cp);

figureLatex;
figWidth = 2.5;%3.5;
set(hfig,'position',[pos(1:2),figWidth,figWidth/sqrt(2)])

eff = 0.944;
mult = eff*30/5;


v_ = [8];
for i = 1
    v = v_(i);
    v_3 = v^3*rho*radius^2*pi/2;
    w = lambda/radius*v;
    Pm = mult*v_3.*cp/1e6;
    [P_opt,idx] = max(Pm);
    w_opt = w(idx);
    plot(w,Pm), hold on
end
P_opt_8 = P_opt
w_opt_8 = w_opt
box off

N=30;
v_ = linspace(0,11.4,N);
P_opt = zeros(1,N);
w_opt = zeros(1,N);
for i = 1:N
    v = v_(i);
    v_3 = v^3*rho*radius^2*pi/2;
    w = lambda/radius*v;
    Pm = mult*v_3.*cp/1e6;
    [P_opt(i),idx] = max(Pm);
    w_opt(i) = w(idx);
end
% P_opt(N+1) = P_opt(N);
% w_opt(N+1) = 10;
plot(w_opt,P_opt,'k--')
ylim([0,13])
xlim([0,1.5])

ylabel('Power [MW]')
xlabel('Rotor speed [rad/s]')
%Legend
l = legend('$v = 8$ [m/s]',...
    'MPP');
set(l,'Interpreter','LaTeX');
set(l,'Box','off')
set(l,'location','best')
set(l,'FontSize',9)
set(l,'Orientation','vertical')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'wind_curve'),'epsc');
end


