function [R_skew] = skew(x)
% Skew symmetric matrix - cross product as matrix multipliation
% R_skew(X)Y = X x Y
R_skew = [0,-x(3),x(2);
          x(3),0,-x(1);
          -x(2),x(1),0];
end

