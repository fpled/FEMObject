function D = calc_opmat(mat,elem,xnode,xgauss)
% function D = calc_opmat(mat,elem,xnode,xgauss)

if nargin<=2
    xnode = [];
    xgauss = [];
end

if israndom(mat)
    D = calc_opmatpc(mat,elem,xnode,xgauss);
    return
end

switch getdim(elem)
    case 1
        EL = evalparam(mat,'EL',elem,xnode,xgauss); % longitudinal Young modulus
        S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
        D = EL*S; % stiffness operator
    case 2
        e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
        EL = evalparam(mat,'EL',elem,xnode,xgauss); % longitudinal Young modulus
        ET = evalparam(mat,'ET',elem,xnode,xgauss); % transverse Young modulus
        nuL = evalparam(mat,'NUL',elem,xnode,xgauss); % longitudinal Poisson ratio
        GL = evalparam(mat,'GL',elem,xnode,xgauss); % longitudinal shear modulus
        n1 = evalparam(mat,'AXISL',elem,xnode,xgauss); % axis of rotational symmetry
        n1 = normalize(VECTOR(n1));
        switch getindim(n1)
            case 2
                if isparam(mat,'AXIST')
                    n2 = evalparam(mat,'AXIST',elem,xnode,xgauss);
                    n2 = normalize(VECTOR(n2));
                    if dot(n1,n2)~=0
                        lastwarn('Wrong axis : change axis T')
                        n2 = rot2D(n1,pi/2);
                    end
                else
                    n2 = rot2D(n1,pi/2);
                end
                syscoord = CARTESIAN2D(n1,n2);
            case 3
                if isparam(mat,'AXIST')
                    n2 = evalparam(mat,'AXIST',elem,xnode,xgauss);
                    n2 = normalize(VECTOR(n2));
                    n3 = cross(n1,n2);
                    n3 = normalize(n3);
                    n2 = cross(n3,n1);
                else
                    [n2,n3] = planortho(n1);
                end
                syscoord = CARTESIAN3D(n1,n2,n3);
            otherwise
                error('Wrong axis of rotational symmetry')
        end
        
        switch getoption(elem)
            case 'DEFO'
                nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % transverse Poisson ratio
                GT = ET/(1+nuT)/2; % transverse shear modulus
                kT = EL*ET/(2*(1-nuT)*EL-4*nuL^2*ET); % transverse bulk modulus
                switch getindim(n1)
                    case 2
                        % [11 22 12]
                        D = e*[EL+4*kT*nuL^2,2*kT*nuL,0;...
                            2*kT*nuL,kT+GT,0;...
                            0,0,GL]; % stiffness operator
                    case 3
                        % [11 22 33 23 13 12]
                        D = e*[EL+4*kT*nuL^2,2*kT*nuL,2*kT*nuL,0,0,0;...
                            2*kT*nuL,kT+GT,kT-GT,0,0,0;...
                            2*kT*nuL,kT-GT,kT+GT,0,0,0;...
                            0,0,0,GT,0,0;...
                            0,0,0,0,GL,0;...
                            0,0,0,0,0,GL]; % stiffness operator
                    otherwise
                        error('Wrong axis of rotational symmetry')
                end
                Pe = calc_Pmat(syscoord);
                % [xx yy xy] in 2D
                % [xx yy zz yz xz xy] in 3D
                D = Pe'*D*Pe;
                if getindim(n1)==3
                    % [xx yy xy]
                    D = [D(1,1),D(1,2),D(1,6);...
                        D(2,1),D(2,2),D(2,6);...
                        D(6,1),D(6,2),D(6,6)];
                end
                
            otherwise
                switch getindim(n1)
                    case 2
                        % [11 22 12]
                        S = 1/e*[1/EL,-nuL/EL,0;...
                            -nuL/EL,1/ET,0;...
                            0,0,1/GL]; % compliance operator
                    case 3
                        nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % transverse Poisson ratio
                        GT = ET/(1+nuT)/2; % transverse shear modulus
                        % [11 22 33 23 13 12]
                        S = 1/e*[1/EL,-nuL/EL,-nuL/EL,0,0,0;...
                            -nuL/EL,1/ET,-nuT/ET,0,0,0;...
                            -nuL/EL,-nuT/ET,1/ET,0,0,0;...
                            0,0,0,1/GT,0,0;...
                            0,0,0,0,1/GL,0;...
                            0,0,0,0,0,1/GL]; % compliance operator
                    otherwise
                        error('Wrong axis of rotational symmetry')
                end
                [~,Ps] = calc_Pmat(syscoord);
                % [xx yy xy] in 2D
                % [xx yy zz yz xz xy] in 3D
                S = Ps'*S*Ps;
                if getindim(n1)==3
                    % [xx yy xy]
                    S = [S(1,1),S(1,2),S(1,6);...
                        S(2,1),S(2,2),S(2,6);...
                        S(6,1),S(6,2),S(6,6)];
                end
                D = inv(S); % stiffness operator
        end
        
    case 3
        EL = evalparam(mat,'EL',elem,xnode,xgauss); % longitudinal Young modulus
        ET = evalparam(mat,'ET',elem,xnode,xgauss); % transverse Young modulus
        nuL = evalparam(mat,'NUL',elem,xnode,xgauss); % longitudinal Poisson ratio
        nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % transverse Poisson ratio
        GL = evalparam(mat,'GL',elem,xnode,xgauss); % longitudinal shear modulus
        GT = ET/(1+nuT)/2; % transverse shear modulus
        kT = EL*ET/(2*(1-nuT)*EL-4*nuL^2*ET); % transverse bulk modulus
        
        n1 = evalparam(mat,'AXISL',elem,xnode,xgauss); % axis of rotational symmetry
        n1 = normalize(VECTOR(n1));
        if isparam(mat,'AXIST')
            n2 = evalparam(mat,'AXIST',elem,xnode,xgauss);
            n2 = normalize(VECTOR(n2));
            n3 = cross(n1,n2);
            n3 = normalize(n3);
            n2 = cross(n3,n1);
        else
            [n2,n3] = planortho(n1);
        end
        syscoord = CARTESIAN3D(n1,n2,n3);
        Pe = calc_Pmat(syscoord);
        % [11 22 33 23 13 12]
        D = [EL+4*kT*nuL^2,2*kT*nuL,2*kT*nuL,0,0,0;...
            2*kT*nuL,kT+GT,kT-GT,0,0,0;...
            2*kT*nuL,kT-GT,kT+GT,0,0,0;...
            0,0,0,GT,0,0;...
            0,0,0,0,GL,0;...
            0,0,0,0,0,GL]; % stiffness operator
        % [xx yy zz yz xz xy]
        D = Pe'*D*Pe;
        
        P = [1,0,0,0,0,0;...
             0,1,0,0,0,0;...
             0,0,1,0,0,0;...
             0,0,0,0,0,1;...
             0,0,0,0,1,0;...
             0,0,0,1,0,0];
        % [xx yy zz xy xz yz]
        D = P'*D*P;
        
end


