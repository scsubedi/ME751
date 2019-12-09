
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtains a skew symmetric matrix associated with algebraic vector a. 
% Input: a(3 x 1) is the Euler Parameter
% Output: A(3 x 3) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = tilda(a)
    A = [0 -a(3) a(2)
        a(3) 0 -a(1)
        -a(2) a(1) 0];
end
