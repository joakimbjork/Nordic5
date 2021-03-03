% Copyright (C) 2021  Joakim Björk
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function hfig = plot_4_compass(w_all,e_all)
h=1;
doubleColumn =true;
figureLatex

for i = 1:4
    w = w_all(:,i);
    e = e_all(i);
    
    [~,idx] = max(abs(w));
    w_ = compass(w(idx),'w--'); hold all
    set(w_, 'Visible', 'on')
    f = abs(e)/(2*pi);
    d = -real(e)/abs(e);
    
    subplot(2,2,i)
    if length(w) == 5
        compass(w(1),'r');
        compass(w(2),'k'); 
        compass(w(3),'b'); 
        compass(w(4),'k--'); 
        compass(w(5),'b--'); 
        l = legend([num2str(round(f,2)), ' Hz',...
                   ' (' num2str(round(d,3)*100),' \% )'],...
                   '1','2','3','4','5') ;
        set(l,'Interpreter','LaTeX');
        set(l,'Box','on')
        set(l,'location','best')
        set(l,'FontSize',9)
        set(l,'Orientation','vertical')
    else
        compass(w);
        l = legend([num2str(round(f,2)), ' Hz',...
                   ' (' num2str(round(d,3)*100), ' \% )']) ;
        set(l,'Interpreter','LaTeX');
        set(l,'Box','on')
        set(l,'location','best')
        set(l,'FontSize',9)
        set(l,'color',[1,0,0])
        set(l,'Orientation','vertical')
    end
end

