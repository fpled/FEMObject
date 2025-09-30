function [direction,angle] = calcrotation_direction_angle(u,v,n)
% function [direction,angle] = calcrotation_direction_angle(u,v,n)
% Calculate single axis direction and angle of rotation corresponding to
%   R = (tilt from z to n) * (twist about z by phi), where phi = atan2(vy,vx).
%
% Inputs
%   v = [vx, vy]     : in-plane tangent vector (twist direction)
%   n = [nx, ny, nz] : out-of-plane target normal vector (tilt direction)
%
% Outputs
%   direction : 1x3 unit axis vector
%   angle     : scalar angle in [0, 2*pi)

v = v(:).';
v = v / norm(v); % unit tangent vector
z  = [0, 0, 1];  % unit z axis

switch u.indim
    case 2
        % Rotate in-plane about z = [0, 0, 1] by angle theta = atan2(vy, vx)
        % using tangent vector v = [vx, vy] = [cos(theta), sin(theta)]
        direction = z;                  % unit z axis
        angle = atan2(v(2), v(1)); % in-plane twist angle theta

    case 3
        n = n(:).';
        n = n / norm(n); % unit normal vector
        
        % Use of quaternions [w, x, y, z]
        % Twist the XY-plane about z = [0, 0, 1] by angle phi = atan2(vy, vx)
        % using tangent vector v = [vx, vy] = [cos(phi), sin(phi)]
        phi = atan2(v(2), v(1));             % in-plane twist angle phi in the XY plane
        qz = [cos(phi/2), 0, 0, sin(phi/2)];
        
        % Tilt the XY-plane normal vector z = [0, 0, 1] into target
        % normal vector n = [nx, ny, nz] using unit tilt axis k = z x n / ||z x n||
        % (rotate about tilt axis k = z x n by tilt angle theta)
        k = cross(z, n); % tilt axis k = z x n that is orthogonal to both z and n
        st = norm(k);    % sin(theta)
        ct = dot(z, n);  % cos(theta)
        if st > eps
            k  = k / st;           % unit tilt axis
            theta = atan2(st, ct); % out-of-plane tilt angle theta
            qt = [cos(theta/2), k*sin(theta/2)];
        else                       % n = [nx, ny, nz] aligned or anti-aligned with z = [0, 0, 1]
            if ct < 0              % n = [nx, ny, nz] aligned with -z = [0, 0, -1]
                qt = [0, 1, 0, 0]; % rotate by pi = 180° about x = [1, 0, 0] (or any in-plane axis)
            else                   % n = [nx, ny, nz] aligned with z = [0, 0, 1]
                qt = [1, 0, 0, 0]; % identity
            end
        end

        % compose: rot = tilt * twist
        q = quatmul(qt, qz);
        q = q / norm(q);

        % convert quaternion to axis/angle
        w = q(1);     % cos(angle/2)
        vq = q(2:4);
        s = norm(vq); % sin(angle/2)
        if s > eps
            direction = vq / s;
            angle = 2*atan2(s, w);
            if angle < 0, angle = angle + 2*pi; end
        else
            direction = z; % arbitrary when angle = 0
            angle = 0; % zero angle
        end
end

end

function q = quatmul(q1, q2)
% quaternion product q = q1 ⨂ q2, with q = [w, x, y, z]
w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
q = [ w1*w2 - x1*x2 - y1*y2 - z1*z2, ...
      w1*x2 + x1*w2 + y1*z2 - z1*y2, ...
      w1*y2 - x1*z2 + y1*w2 + z1*x2, ...
      w1*z2 + x1*y2 - y1*x2 + z1*w2  ];
end
