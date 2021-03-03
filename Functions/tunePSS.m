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

function [F] = tunePSS(e,residu)
% G = ss(G);
% [v,e,w] = modal_analysis(G.A);
% 
% Observability and controllability
% v_obsv = G.C*v;
% w_ctrb = w*G.B;
% 
% %%%%%%%%%%%%%% Tune PSS1 %%%%%%%%%%%%%%%%%%%%
% Give the residue of the selected mode, see equations (6.52) and (6.81)
% residu=  v_obsv*w_ctrb;
% See Examples 6.4 to 6.6 in the textbook
% Real and imaginary parts of mode
sig=real(e);
wp=abs(imag(e));
% The angle of the residue
ang_R=angle(residu);
ang_phi=sign(ang_R)*pi-ang_R;
%%
% See equation (6.82)
if abs(ang_phi) > pi/3 && abs(ang_phi) <= 2*pi/3
    num_of_fiters=  2;
    sgnK= 1;
elseif abs(ang_phi) >= 2*pi/3 && abs(ang_phi) <= pi
    ang_phi=-ang_R;
    sgnK=  -1;
    num_of_fiters=  1;
else
    sgnK= 1;
    num_of_fiters=   1;
end
% Find alpha given in equation (6.83)
alph= (1+sin(ang_phi/num_of_fiters))/(1-sin(ang_phi/num_of_fiters)) ;
% Find T given in equation (6.83)
TT=  1/(wp*sqrt(alph)) ;
% Find T1-T4 given in equation (6.83)
T1= alph*TT ;
T2= TT ;
s = tf('s');
F = (T1*s + 1)/(T2*s +1);
if num_of_fiters==2
    F = F*F;
end
dep_ang = ang_R+angle(freqresp(F,wp));
disp(['Angle of departure ', num2str(dep_ang*180/pi),...
      ' = ', num2str(ang_R*180/pi) , '+' ,...
      num2str(angle(freqresp(F,wp))*180/pi) , ' deg '])

F = sgnK*F;
end

