function R = calcrotation(u,v,n)
% function R = calcrotation(u,v,n)

switch u.indim
    case 2
        % Rotate in-plane around z = [0, 0, 1] by angle of rotation theta =
        % atan2(vy, vx) using tangent vector v = [vx, vy]
        v = v / norm(v); % unit tangent vector v = [cos(theta), sin(theta)]
        R = [v(1), -v(2);
             v(2),  v(1)];
        % theta = atan2(v(2), v(1)); % in-plane twist angle theta
        % R = [cos(theta), -sin(theta);
        %      sin(theta),  cos(theta)];
    case 3
        %% Old version
        % Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
        % atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
        % rotation formula
        % n = n / norm(n);
        % Q = [    0, -n(3),  n(2);
        %       n(3),     0, -n(1);
        %      -n(2),  n(1),    0]; % skew for n = [nx, ny, nz]
        % v = [C.vx, C.vy];
        % v = v / norm(v); % unit tangent vector v = [cos(phi), sin(phi)]
        % % R = eye(3) + v(2)*Q + (1-v(1))*Q^2;
        % R = v(1)*eye(3) + v(2)*Q + (1-v(1))*(n'*n);
        % % phi = atan2(v(2), v(1)); % in-plane twist angle phi
        % % R = eye(3) + sin(phi)*Q + (1-cos(phi))*Q^2;
        % % R = cos(phi)*eye(3) + sin(phi)*Q + (1-cos(phi))*(n'*n);
        
        %% New version
        % Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
        % using tangent vector v = [vx, vy], then tilt from z axis to normal
        % vector n = [nx, ny, nz] so that the circle's normal is n regardless
        % of v = [vx, vy]
        
        % Twist XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
        % using tangent vector v = [vx, vy] = [cos(phi), sin(phi)]
        v = v / norm(v); % unit tangent vector v = [cos(phi), sin(phi)]
        z = [0, 0, 1];
        % Qz = [    0, -z(3),  z(2);
        %        z(3),     0, -z(1);
        %       -z(2),  z(1),    0]; % skew for z = [0, 0, 1]
        Qz = [0, -1, 0;
              1,  0, 0;
              0,  0, 0]; % skew for z = [0, 0, 1]
        % R_twist = eye(3) + v(2)*Qz + (1-v(1))*Qz^2;
        R_twist = v(1)*eye(3) + v(2)*Qz + (1-v(1))*(z'*z);
        % phi = atan2(v(2), v(1)); % in-plane twist angle phi in the XY plane
        % R_twist = eye(3) + sin(phi)*Qz + (1-cos(phi))*Qz^2;
        % R_twist = cos(phi)*eye(3) + sin(phi)*Qz + (1-cos(phi))*(z'*z);
        
        % Tilt the XY-plane normal vector z = [0, 0, 1] into target
        % normal vector n = [nx, ny, nz] using tilt axis k = z x n / ||z x n||
        n = n / norm(n); % unit normal vector
        k = cross(z, n); % tilt axis k = z x n that is orthogonal to both z and n
        st = norm(k);    % sin(theta)
        ct = dot(z, n);  % cos(theta)
        if st > eps
            k = k / st;  % unit tilt axis
            Qk = [    0, -k(3),  k(2);
                   k(3),     0, -k(1);
                  -k(2),  k(1),    0]; % skew for k = [kx, ky, kz]
            % R_tilt = eye(3) + st*Qk + (1-ct)*Qk^2;
            R_tilt = ct*eye(3) + st*Qk + (1-ct)*(k'*k);
            % theta = acos(ct); % out-of-plane angle theta
            % R_tilt = eye(3) + sin(theta)*Qk + (1-cos(theta))*Qk^2;
            % R_tilt = cos(theta)*eye(3) + sin(theta)*Qk + (1-cos(theta))*(k'*k);
        else % n = [nx, ny, nz] exactly aligned or anti-aligned with z = [0, 0, 1]
            if ct < 0 % n = [nx, ny, nz] aligned with -z = [0, 0, -1]
                R_tilt = diag([1,-1,-1]);  % rotate by pi about x = [1, 0, 0] (or any in-plane axis)
            else      % n = [nx, ny, nz] aligned with z = [0, 0, 1]
                R_tilt = eye(3);
            end
        end
        
        % Compose: first twist in the XY‐plane, then tilt into the plane with normal n
        R = R_tilt * R_twist; % constructed for column vectors
    otherwise
        error('Wrong space dimension');
end

R = R'; % constructed for row vectors
