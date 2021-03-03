% Copyright (C) 2021  Joakim Björk
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
function F = leadFilter(phase_add,wc)
% adds phase_add radians at frequency wc.
x = tan(phase_add);
p = -2-4*x^2;
q = 1;
beta = -p/2 - sqrt(p^2/4 -q);
tD = 1/(wc*sqrt(beta));

s= tf('s');
F = (tD*s+1)/(beta*tD*s+1);
end

