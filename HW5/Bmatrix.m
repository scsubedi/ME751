
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtains the Bmatrix from euler parameters. 
% Input: p(4x1) is the Euler Parameter
%        a(3x1) is the vector in local coordinate system
% Output: B(4 x 3) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = Bmatrix(p,a)

B = 2*[(p(1)*[1 0 0;0 1 0;0 0 1] + tilda(p(2:end)))*a...
    p(2:end)*a' - (p(1)*eye(3) + tilda(p(2:end)))*tilda(a)];

end