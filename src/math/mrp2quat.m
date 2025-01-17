function quat = mrp2quat(mrp)
    b = vecnorm(mrp,2,2);
    theta = 4*atan(b);
    n = zeros(size(mrp));

    n(b~=0,:) = mrp(b~=0,:)./b(b~=0,:);
    quat = [cos(theta./2),n.*sin(theta./2)];
end

