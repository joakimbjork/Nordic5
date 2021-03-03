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

v_ = [11.4,10,8];
for i = 1:3
    v = v_(i);
    v_3 = v^3*rho*radius^2*pi/2;
    w = lambda/radius*v;
    Pm = v_3.*cp/1e6;
    [P_opt,idx] = max(Pm);
    w_opt = w(idx);
    plot(w,Pm), hold on
end
box off

N=30;
v_ = linspace(0,11.4,N);
P_opt = zeros(1,N);
w_opt = zeros(1,N);
for i = 1:N
    v = v_(i);
    v_3 = v^3*rho*radius^2*pi/2;
    w = lambda/radius*v;
    Pm = v_3.*cp/1e6;
    [P_opt(i),idx] = max(Pm);
    w_opt(i) = w(idx);
end
% P_opt(N+1) = P_opt(N);
% w_opt(N+1) = 10;
plot(w_opt,P_opt,'k--')
ylim([0,1.1*Pb_m/1e6])
xlim([0,2.0])

ylabel('Turbine power, $P_m$ [MW]')
xlabel('Turbine speed, $\Omega$ [rad/s]')
%Legend
l = legend('$v = 11.4$ [m/s]','$v = 10$ [m/s]','$v = 8$ [m/s]',...
    'MPP');
set(l,'Interpreter','LaTeX');
set(l,'Box','off')
set(l,'location','best')
set(l,'FontSize',9)
set(l,'Orientation','vertical')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'wind_curves'),'epsc');
end


%% Draw power coefficient derivative
cp = wt.cp.table;
cp = cp(1,:);
lambda = wt.cp.tsr;
N = length(lambda);
[c_opt,idx] = max(cp);
lambda_opt = lambda(idx);

figureLatex;
x = lambda/lambda_opt;
x_mod = x(2:N-1);

for i = 2:N-1
    i1 = i-1;
    i2 = i+1;
    dx = x(i2)-x(i1);
    d_cp = cp(i2)-cp(i1);
    d_cp_x(i-1) = d_cp/dx;
end

x_min = 0.8;
[~,idx] = min(abs(x_mod-x_min));
slope = d_cp_x(idx);
yline(slope,'r'), hold on
plot(x_mod,d_cp_x); hold all
plot(x_mod,d_cp_x./x_mod,'--','color',[0.7,0.7,0.7]);

xline(x_min,'r')



xline(1,'k')
yline(0,'k')

xticks(-0.8:0.2:2)
xlim([0.5,1.5])
ylim([-0.4,1.2])
% yticks([-0.5,0,round(slope,2),0.5,1])

l = legend([num2str(round(slope,2))],...
           '$\frac{\partial}{\partial x}c_p$',...
           '$\frac{1}{x}\frac{\partial}{\partial x}c_p$');

set(l,'Interpreter','LaTeX');
set(l,'Box','off')
set(l,'location','best')
set(l,'FontSize',9)
set(l,'Orientation','vertical')

ylabel('Power coefficient derivative')
xlabel('Normalized speed ratio, $x$ ')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'wind_Pm_w_ratio'),'epsc');
end

%% Draw normalized power/speed curve
figureLatex;
cp = wt.cp.table;
cp = cp(1,:);
[c_opt,idx] = max(cp);
lambda_opt = lambda(idx);
x = lambda/lambda_opt;

plot(x,cp,'k'), hold on
box off
xline(1,'k'); hold on
yline(c_opt,'k')


[~,idx] = min(abs(x-x_min));
cp_08 = cp(idx);
x2 = [0,x_min,2];
y2 = [-x_min*slope,0,1.2*slope]+cp_08;
plot(x2,y2,'r','LineWidth',0.5), hold on
xline(0.8,'r')


ylim([0.0,0.55])
xlim([0.0,2])

xticks([0,0.50,x_min,1,1.5,2])

l = legend('$c_p$',...
           ['$\lambda_\mathrm{opt} = $', num2str(round(lambda_opt,2))],...
           ['$c_\mathrm{opt} = $', num2str(round(c_opt,2))],...
           ['$\frac{\partial}{\partial x}c_p =$', num2str(round(slope,2))]);

set(l,'Interpreter','LaTeX');
set(l,'Box','off')
set(l,'location','best')
set(l,'FontSize',9)
set(l,'Orientation','vertical')


ylabel('Power coefficient, $c_p$ ')
xlabel('Normalized speed ratio, $x$ ')

saveFigure = true;
if saveFigure
    saveas(hfig,fullfile(path, 'wind_Pm_w_slope'),'epsc');
end
