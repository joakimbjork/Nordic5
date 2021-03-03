hfig = figure;
co = ([0;0.6]*[1,1,1]);
co = ([0,0,0;1,0,0;0,0.5,0.7]);

set(groot,'defaultAxesColorOrder',co)
set(0,'DefaultTextInterpreter','latex',...
    'DefaultTextFontSize',10,...
    'DefaultAxesFontName','Times',...
    'DefaultAxesFontSize',10,...
    'DefaultLineMarkerSize',6)
set(gca,'TickLabelInterpreter','latex')
set(hfig,'DefaultLineLineWidth',0.8)

set(hfig,'units','inches',...
    'NumberTitle','off','Name','Frequency');
if exist('figtitle','Var')
    set(hfig,'Name',figtitle);
end
pos = get(hfig,'position');

%% Figure size 3.5 (inches) good for two-column article
figWidth0 = 3.5;
figWidth = figWidth0;
if exist('doubleColumn','var')
    if doubleColumn == true
        figWidth = 5;
        doubleColumn = false;
    end
end

if exist('h','var')
set(hfig,'position',[pos(1:2),figWidth,figWidth0*h])
clear h
else
set(hfig,'position',[pos(1:2),figWidth,figWidth0/sqrt(3)])
end
%% Label
% ylabel('SM $f$s [Hz]')
% xlabel('Time [s]')
% set(gca,'YTick',[49.8, 49.9, 50])

%% Legends with Latex text interpreter
% l = legend('$\Delta f_1$', '', '$\Delta f_2$');
%         set(l,'Interpreter','LaTeX');
%         set(l,'Box','off')
%         set(l,'location','best')
%         set(l,'FontSize',9)
%         set(l,'Orientation','vertical')

%% Size of figure
% set(gca,'OuterPosition',[0.2 0.2 0.8 0.4]);


% set(hfig,'position',[pos(1:2),3.5,2.5])

