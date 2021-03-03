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

function Gr = modelred_hsv(G,varargin)
%
% sys_r = modelred_hsv(sys)
%
% makes a balanced model reduction of the system G such that Gr is of 
% order n, where n has to be less than or equal to the order of G.

    
[HSV,BALDATA] = hsvd(G);
m = BALDATA.Split;
% hsv_unstable = HSV( 1:m(1) );
hsv_stable = HSV( m(1) + (1:m(2)) );

  
if nargin == 1
figure()
bar(hsv_stable)
if hsv_stable(1)-hsv_stable(2) > 1e3
    ylim([0,hsv_stable(2)]);
end
title('Hankel Singular Values (State Contributions of Stable Part)')
xlabel('State')
ylabel('State Energy')

disp(['The system has ', num2str(m(1)),...
      ' unstable poles. These cannot be reduced.'])
  
prompt = 'What order should the reduced stable part be? n_stab = ';
n = input(prompt);

else
    n = varargin{1};
    if n < 1 % Automatically select number of states
        eps = n;
        idx = find(hsv_stable<eps);
        for i = 1:length(idx)
            ii = max(idx(i)-1,1);
            ratio = hsv_stable(idx(i))/hsv_stable(ii);
            if ratio < 0.5            
                break
            end
        end
        n = idx(ii);
    end
end
disp( ['The total order is ', num2str( m(1)+n ) ]) 
Gr = balred(G,m(1)+n,BALDATA);
