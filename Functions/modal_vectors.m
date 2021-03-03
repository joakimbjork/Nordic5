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

function [E,W,C,data] = modal_vectors(A,varargin)
% Sort eigenvalues end eigenvectors. Starting with oscilatory modes up to
% an eigenfrequency of w_osc_max rad/s. Modes are sorted according to
% damping ratio.

[~,E,W] = eig(A); 
% W'*A*inv(W') = E;  
E = diag(E);
n = length(E);

%  sort eigenvalues slower than w_osc_max
w_osc_min = 0;
if nargin == 1
    w_osc_max = 2*2*pi;
else
    w_osc_max = varargin{1};
    if length(w_osc_max) == 2
        temp = w_osc_max;
        w_osc_max = max(temp);
        w_osc_min = min(temp);
    end
end

idx1 = find(abs(imag(E))<w_osc_max &...
            abs(imag(E))>w_osc_min &...
            abs(imag(E))>1e-6);
E1 = E(idx1);
W1 = W(:,idx1); %eigenvectors

% sort after damping ratio
damping = -real(E1)./abs(E1); 
[damping1,idx_sort] = sort(damping);

% sort after frequency
% [~,idx_sort] = sort(abs(E1));

E1 = E1(idx_sort);
W1 = W1(:,idx_sort);

% sort remainding eigenvalues
idx2 = 1:n;
idx2(idx1) = [];
E2 = E(idx2);
W2 = W(:,idx2);
% sort after damping ratio
damping = -real(E2)./abs(E2);
[damping2,idx_sort] = sort(damping);
E2 = E2(idx_sort);
W2 = W2(:,idx_sort);
E = [E1;E2];
W = [W1,W2];

C = zeros(n,n);

i = 1;
while i <= n
    if imag(E(i))==0
        C(:,i) = real(W(:,i));
        i = i+1;
    else
        C(:,[0,1]+i) = [real(W(:,i)),imag(W(:,i+1))];
        i = i+2;
    end
end

C = C'; % modes = C*states

data.freq = abs(E)/(2*pi);
data.damping = -real(E)./abs(E);
end

