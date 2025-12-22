%% Circular plate - Deterministic linear elasticity problem %%
%%----------------------------------------------------------%%
% Code_Aster v3.03.100.pdf
% SSLS100 - Plaque circulaire encastree soumise a une pression uniforme
% Code_Aster v3.03.101.pdf
% SSLS101 - Plaque circulaire posee soumise a une pression uniforme

% clc
clearvars
close all

%% Input data
solveProblem = true;
displaySolution = false;
displayCv = true;

% boundaries = {'SimplySupported'};
% boundaries = {'Clamped'};
boundaries = {'SimplySupported','Clamped'};
% loadings = {'Uniform'};
% loadings = {'Concentrated'};
loadings = {'Uniform','Concentrated'};
% elemtypes = {'DKT'};
% elemtypes = {'DKQ'};
% elemtypes = {'DST'};
% elemtypes = {'DSQ'};
% elemtypes = {'COQ4'};
% elemtypes = {'DKT','DKQ'}; % Kirchhoff-Love (classical) plate theory
% elemtypes = {'DST','DSQ','COQ4'}; % Reissner-Mindlin (first-order shear) plate theory
elemtypes = {'DKT','DKQ','DST','DSQ','COQ4'}; % Both plate theories
nbelems = 2.^(1:6);

fontsize = 16;
linewidth = 1;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

for ib=1:length(boundaries)
    boundary = boundaries{ib};
    
for il=1:length(loadings)
    loading = loadings{il};
    filename = ['plateCircDetLinElas' boundary loading];
    if displayCv
        hcvUz = figure('Name','Evolution of error indicator for Uz w.r.t number of elements');
        clf
        hcvRt = figure('Name','Evolution of error indicator for Rt w.r.t number of elements');
        clf
        hcvRx = figure('Name','Evolution of error indicator for Rx w.r.t number of elements');
        clf
        hcvRy = figure('Name','Evolution of error indicator for Ry w.r.t number of elements');
        clf
        htime = figure('Name','Evolution of CPU time w.r.t number of elements');
        clf
        leg = cell(1,length(elemtypes));
    end
    
for ie=1:length(elemtypes)
    elemtype = elemtypes{ie};
    pathname = fullfile(getfemobjectoptions('path'),'BASIC',...
        'EXAMPLES','SHELLS',filename,elemtype);
    if ~exist(pathname,'dir')
        mkdir(pathname);
    end
    leg{ie} = elemtype;
    
err_Uz = zeros(1,length(nbelems));
err_Rt = zeros(1,length(nbelems));
err_Rx = zeros(1,length(nbelems));
err_Ry = zeros(1,length(nbelems));
time = zeros(1,length(nbelems));
Nbelem = zeros(1,length(nbelems));
for i=1:length(nbelems)
    
%% Problem
if solveProblem
    %% Domains and meshes
    r = 1; % [m]
    C = CIRCLE(0.0,0.0,0.0,r);
    
    P_load = getcenter(C);
    x_load = double(P_load);
    
    cl = r./nbelems(i);
    switch lower(loading)
        case 'uniform'
            S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]));
            % S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]),'all'); % option 'all' for facets, ridges, peaks
        case 'concentrated'
            S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]),'points',x_load);
            % S = build_model(C,'cl',cl,'elemtype',elemtype,'filename',fullfile(pathname,['gmsh_plate_circ_' elemtype '_cl_' num2str(cl)]),'points',x_load,'all'); % option 'all' for facets, ridges, peaks
    end
    
    %% Materials
    % Gravitational acceleration
    g = 10; % [m/s2]
    % Young modulus
    E = 1; % [Pa]
    % Poisson ratio
    NU = 0.3;
    % Density
    RHO = 1; % [kg/m3]
    % Thickness
    h = 0.1; % [m]
    % Extensional stiffness (or Membrane rigidity)
    A_rig = E*h/(1-NU^2);
    % Bending stiffness (or Flexural rigidity)
    D_rig = E*h^3/(12*(1-NU^2));
    
    % Material
    mat = ELAS_SHELL('E',E,'NU',NU,'RHO',RHO,'DIM3',h,'k',5/6);
    mat = setnumber(mat,1);
    S = setmaterial(S,mat);
    
    %% Dirichlet boundary conditions
    x_support = cellfun(@(x) x',getvertices(C),'UniformOutput',false);
    x_support = [x_support{:}]';
    P_support = POINT(x_support);
    
    S = final(S);
    switch lower(boundary)
        case 'clamped'
            % No locking
            S = addcl(S,[]);
            % Shear locking for element COQ4
            % S = addcl(S,P_support);
            % Partial shear locking for element COQ4
            % S = addcl(S,getfacet(S,2));
            % S = addcl(S,P_support([1;4]));
        case 'simplysupported'
            % No locking
            S = addcl(S,[],'U');
            % Shear locking for element COQ4
            % S = addcl(S,P_support,'U');
            % Partial shear locking for element COQ4
            % S = addcl(S,getfacet(S,2),'U');
            % S = addcl(S,P_support([1;4]),'U');
    end
    % S = addcl(S,[],'R');
    
    %% Stiffness matrix and sollicitation vector
    switch lower(loading)
        case 'uniform' % Uniform transverse load per unit area applied on the plate surface
            p = RHO*g*h; % surface load (body load for plates) [N/m2]
        case 'concentrated' % Concentrated transverse load applied at point P_load
            Sec = pi*r^2;
            p = RHO*g*h*Sec; % point load [N]
    end
    % Moment per unit length applied on the plate boundary
    % (only for simply supported plate)
    c = 0; % [N.m/m]
    
    A = calc_rigi(S);
    
    switch lower(loading)
        case 'uniform'
            f = bodyload(S,[],'FZ',-p);
        case 'concentrated'
            f = nodalload(S,P_load,'FZ',-p);
            if isempty(ispointin(P_load,POINT(S.node)))
                error('Pointwise load must be applied to a node of the mesh')
            end
    end
    if strcmpi(boundary,'simplysupported')
        f = f + surfload(S,[],{'MX','MY'},-c*[1;1]);
    end
    
    %% Solution
    t = tic;
    u = A\f;
    time_i = toc(t);
    
    x = getcoord(S.node);
    t = cart2pol(x(:,1),x(:,2),x(:,3));
    funr = @(x,y,theta) dot([cos(theta),sin(theta)],[x,y],2);
    funt = @(x,y,theta) dot([-sin(theta),cos(theta)],[x,y],2);
    funx = @(r,t,theta) dot([cos(theta),-sin(theta)],[r,t],2);
    funy = @(r,t,theta) dot([sin(theta),cos(theta)],[r,t],2);
    
    u = unfreevector(S,u);
    
    % U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    % Ux = u(findddl(S,'UX'),:);
    % Uy = u(findddl(S,'UY'),:);
    Uz = u(findddl(S,'UZ'),:);
    % Ur = funr(Ux,Uy,t);
    % Ut = funt(Ux,Uy,t);
    
    % R = u(findddl(S,DDL(DDLVECT('R',S.syscoord,'ROTA'))),:);
    Rx = u(findddl(S,'RX'),:);
    Ry = u(findddl(S,'RY'),:);
    % Rz = u(findddl(S,'RZ'),:);
    % Rr = funr(Rx,Ry,t);
    Rt = funt(Rx,Ry,t);
    
    %% Reference solution
    switch elemtype
        case {'DKT','DKQ'} % Kirchhoff-Love
            phi = 0;
        case {'DST','DSQ','COQ4'} % Reissner-Mindlin (first-order shear) plate theory
            phi = 8/(3*getparam(mat,'k'))*(h/r)^2/(1-NU); % 16/5*(h/r)^2/(1-NU);
    end
    switch lower(loading)
        case 'uniform'
            switch lower(boundary)
                case 'clamped'
                    fun_Uz = @(x) -p/(64*D_rig) * (r^2 - (x(:,1).^2+x(:,2).^2)) .* (r^2*(1+phi) - (x(:,1).^2+x(:,2).^2));
                    fun_Rt = @(x) -p/(16*D_rig) * sqrt(x(:,1).^2+x(:,2).^2) .* (r^2 - (x(:,1).^2+x(:,2).^2));
                case 'simplysupported'
                    fun_Uz = @(x) -1/(2*D_rig) * (r^2 - (x(:,1).^2+x(:,2).^2)) .* (p/32*(((5+NU)/(1+NU)+phi)*r^2 - (x(:,1).^2+x(:,2).^2)) + c/(1+NU));
                    fun_Rt = @(x) -1/D_rig * sqrt(x(:,1).^2+x(:,2).^2) .* (p/16*(((3+NU)/(1+NU))*r^2 - (x(:,1).^2+x(:,2).^2)) + c/(1+NU));
            end
        case 'concentrated'
            switch lower(boundary)
                case 'clamped'
                    fun_Uz = @(x) -p/(16*pi*D_rig) * (r^2 - (x(:,1).^2+x(:,2).^2) + (2*(x(:,1).^2+x(:,2).^2) - 1/2*phi*r^2) .* log(sqrt(x(:,1).^2+x(:,2).^2))./r);
                    fun_Rt = @(x) p/(4*pi*D_rig) * sqrt(x(:,1).^2+x(:,2).^2) .* log(sqrt(x(:,1).^2+x(:,2).^2)./r);
                case 'simplysupported'
                    fun_Uz = @(x) -p/(16*pi*D_rig) * ((3+NU)/(1+NU)*(r^2 - (x(:,1).^2+x(:,2).^2)) + (2*(x(:,1).^2+x(:,2).^2) - 1/2*phi*r^2) .* log(sqrt(x(:,1).^2+x(:,2).^2)./r)) - c/(2*D_rig*(1+NU))*(r^2 - (x(:,1).^2+x(:,2).^2));
                    fun_Rt = @(x) p/(4*pi*D_rig) * sqrt(x(:,1).^2+x(:,2).^2) .* (log(sqrt(x(:,1).^2+x(:,2).^2)./r) - 1/(1+NU)) - c/(D_rig*(1+NU))*sqrt(x(:,1).^2+x(:,2).^2);
            end
    end
    fun_Uz = UserDefinedFunction(fun_Uz,3);
    fun_Rt = UserDefinedFunction(fun_Rt,3);
    fun_Uz.evaluationAtMultiplePoints = true;
    fun_Rt.evaluationAtMultiplePoints = true;
    
    Ux_ex = sparse(getnbnode(S),1);
    Uy_ex = sparse(getnbnode(S),1);
    Uz_ex = fun_Uz(x);
    Rt_ex = fun_Rt(x);
    if strcmpi(loading,'concentrated')
        switch lower(boundary)
            case 'clamped'
                if phi==0 % Kirchhoff-Love
                    Uz_ex(isnan(Uz_ex)) = -p/(16*pi*D_rig) * r^2;
                else % Reissner-Mindlin (first-order shear) plate theory
                    Uz_ex(isnan(Uz_ex)) = -Inf;
                end
                Rt_ex(isnan(Rt_ex)) = 0;
            case 'simplysupported'
                if phi==0 % Kirchhoff-Love
                    Uz_ex(isnan(Uz_ex)) = -p/(16*pi*D_rig) * (3+NU)/(1+NU)*r^2 - c/(2*D_rig*(1+NU))*r^2;
                else % Reissner-Mindlin (first-order shear) plate theory
                    Uz_ex(isnan(Uz_ex)) = -Inf;
                end
                Rt_ex(isnan(Rt_ex)) = 0;
        end
    end
    Rx_ex = funx(zeros(size(Rt_ex)),Rt_ex,t);
    Ry_ex = funy(zeros(size(Rt_ex)),Rt_ex,t);
    Rz_ex = sparse(getnbnode(S),1);
    u_ex = [Ux_ex Uy_ex Uz_ex Rx_ex Ry_ex Rz_ex]';
    u_ex = u_ex(:);
    
    ind_Uz = find(~isinf(Uz_ex));
    err_Uz_i = norm(Uz(ind_Uz)-Uz_ex(ind_Uz))/norm(Uz_ex(ind_Uz));
    err_Rt_i = norm(Rt-Rt_ex)/norm(Rt_ex);
    err_Rx_i = norm(Rx-Rx_ex)/norm(Rx_ex);
    err_Ry_i = norm(Ry-Ry_ex)/norm(Ry_ex);
    
    %% Test solution
    P = getcenter(C);
    xP = double(P);
    tP = cart2pol(xP(:,1),xP(:,2),xP(:,3));
    
    ux = eval_sol(S,u,P,'UX');
    uy = eval_sol(S,u,P,'UY');
    uz = eval_sol(S,u,P,'UZ');
    ur = funr(ux,uy,tP);
    ut = funt(ux,uy,tP);
    
    rx = eval_sol(S,u,P,'RX');
    ry = eval_sol(S,u,P,'RY');
    rz = eval_sol(S,u,P,'RZ');
    rr = funr(rx,ry,tP);
    rt = funt(rx,ry,tP);
    
    uz_ex = fun_Uz(xP);
    rt_ex = fun_Rt(xP);
    if eq(P,getcenter(C)) && strcmpi(loading,'concentrated')
        switch lower(boundary)
            case 'clamped'
                if phi==0 % Kirchhoff-Love
                    uz_ex = -p/(16*pi*D_rig) * r^2;
                else % Reissner-Mindlin (first-order shear) plate theory
                    uz_ex = -Inf;
                end
            case 'simplysupported'
                if phi==0 % Kirchhoff-Love
                    uz_ex = -p/(16*pi*D_rig) * (3+NU)/(1+NU)*r^2 - c/(2*D_rig*(1+NU))*r^2;
                else % Reissner-Mindlin (first-order shear) plate theory
                    uz_ex = -Inf;
                end
        end
        rt_ex = 0;
    end
    rx_ex = funx(0,rt_ex,tP);
    ry_ex = funy(0,rt_ex,tP);
    
    err_uz = norm(uz-uz_ex)/norm(uz_ex);
    err_rt = norm(rt-rt_ex)/norm(rt_ex);
    err_rx = norm(rx-rx_ex)/norm(rx_ex);
    err_ry = norm(ry-ry_ex)/norm(ry_ex);
    
    %% Save variables
    save(fullfile(pathname,['problem_' num2str(i) '.mat']),'S','C','r','h','f');
    save(fullfile(pathname,['solution_' num2str(i) '.mat']),'u','time_i',...
        'Uz','Rt','Rx','Ry');
    save(fullfile(pathname,['reference_solution_' num2str(i) '.mat']),'u_ex',...
        'Uz_ex','Rt_ex','Rx_ex','Ry_ex',...
        'err_Uz_i','err_Rt_i','err_Rx_i','err_Ry_i');
    save(fullfile(pathname,['test_solution_' num2str(i) '.mat']),'P',...
        'ux','uy','uz','ur','ut',...
        'rx','ry','rz','rr','rt',...
        'uz_ex','rt_ex','rx_ex','ry_ex',...
        'err_uz','err_rt','err_rx','err_ry');
else
    load(fullfile(pathname,['problem_' num2str(i) '.mat']),'S','C','r','h','f');
    load(fullfile(pathname,['solution_' num2str(i) '.mat']),'u','time_i',...
        'Uz','Rt','Rx','Ry');
    load(fullfile(pathname,['reference_solution_' num2str(i) '.mat']),'u_ex',...
        'Uz_ex','Rt_ex','Rx_ex','Ry_ex',...
        'err_Uz_i','err_Rt_i','err_Rx_i','err_Ry_i');
    load(fullfile(pathname,['test_solution_' num2str(i) '.mat']),'P',...
        'ux','uy','uz','ur','ut',...
        'rx','ry','rz','rr','rt',...
        'uz_ex','rt_ex','rx_ex','ry_ex',...
        'err_uz','err_rt','err_rx','err_ry');
end

%% Outputs
err_Uz(i) = err_Uz_i;
err_Rt(i) = err_Rt_i;
err_Rx(i) = err_Rx_i;
err_Ry(i) = err_Ry_i;
time(i) = time_i;
Nbelem(i) = getnbelem(S);

fprintf('Circular plate\n');
fprintf('\n');
fprintf(['boundary : ' boundary '\n']);
fprintf(['load     : ' loading '\n']);
fprintf(['mesh     : ' elemtype ' elements\n']);
fprintf('nb elements = %g\n',Nbelem(i));
fprintf('nb nodes    = %g\n',getnbnode(S));
fprintf('nb dofs     = %g\n',getnbddl(S));
fprintf('span-to-thickness ratio = %g\n',r/h);
fprintf('error = %.3e for Uz\n',err_Uz(i));
fprintf('      = %.3e for Rt\n',err_Rt(i));
fprintf('      = %.3e for Rx\n',err_Rx(i));
fprintf('      = %.3e for Ry\n',err_Ry(i));
fprintf('elapsed time = %f s\n',time(i));

fprintf('\n');
fprintf('Displacement u at point (%g,%g,%g)\n',double(P));
fprintf('ux    = %g\n',ux);
fprintf('uy    = %g\n',uy);
fprintf('uz    = %g\n',uz);
fprintf('uz_ex = %g, error = %.3e\n',uz_ex,err_uz);
fprintf('ur    = %g\n',ur);
fprintf('ut    = %g\n',ut);

fprintf('\n');
fprintf('Rotation r at point (%g,%g,%g)\n',double(P));
fprintf('rx    = %g\n',rx);
fprintf('rx_ex = %g, error = %.3e\n',rx_ex,err_rx);
fprintf('ry    = %g\n',ry);
fprintf('ry_ex = %g, error = %.3e\n',ry_ex,err_ry);
fprintf('rz    = %g\n',rz);
fprintf('rr    = %g\n',rr);
fprintf('rt    = %g\n',rt);
fprintf('rt_ex = %g, error = %.3e\n',rt_ex,err_rt);

%% Display
if displaySolution
    %% Display domains, boundary conditions and meshes
    plotDomain(C,'solid',true,'legend',false);
    mysaveas(pathname,['domain_' num2str(i)],formats,renderer);
    mymatlab2tikz(pathname,['domain_' num2str(i) '.tex']);
    
    [hD,legD] = plotBoundaryConditions(S,'legend',false);
    switch lower(loading)
        case 'uniform'
            ampl = 2;
        case 'concentrated'
            ampl = 0.2;
    end
    [hN,legN] = vectorplot(S,'F',f,ampl,'r','LineWidth',linewidth);
    % legend([hD,hN],[legD,legN],'Location','NorthEastOutside')
    mysaveas(pathname,['boundary_conditions_' num2str(i)],formats,renderer);
    
    plotModel(S,'Color','k','FaceColor','k','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,['mesh_' num2str(i)],formats,renderer);
    
    U = u(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
    ampl = getsize(S)/max(abs(U))/10;
    
    plotModelDeflection(S,u,'ampl',ampl,'Color','b','FaceColor','b','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,['mesh_deflected' num2str(i)],formats,renderer);
    
    plotModelDeflection(S,u_ex,'ampl',ampl,'Color','r','FaceColor','r','FaceAlpha',0.1,'legend',false);
    mysaveas(pathname,['mesh_deflected_ex' num2str(i)],formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*u,'Color','b','FaceColor','b','FaceAlpha',0.1);
    mysaveas(pathname,['meshes_deflected_' num2str(i)],formats,renderer);
    
    figure('Name','Meshes')
    clf
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1);
    plot(S+ampl*u_ex,'Color','r','FaceColor','r','FaceAlpha',0.1);
    mysaveas(pathname,['meshes_deflected_ex_' num2str(i)],formats,renderer);
    
    %% Display solution
    % ampl = 0;
    ampl = getsize(S)/max(abs(U))/10;
    options = {'solid',true};
    % options = {};
    
    plotSolution(S,u,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,['Uz_' num2str(i)],formats,renderer);
    
    plotSolution(S,u_ex,'displ',3,'ampl',ampl,options{:});
    mysaveas(pathname,['Uz_ex_' num2str(i)],formats,renderer);
    
%     plotSolution(S,norm(u-u_ex),'displ',3,'ampl',ampl,options{:});
%     mysaveas(pathname,['err_Uz_' num2str(i)],formats,renderer);
    
    plotSolution(S,u,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,['Rx_' num2str(i)],formats,renderer);
    
    plotSolution(S,u_ex,'rotation',1,'ampl',ampl,options{:});
    mysaveas(pathname,['Rx_ex_' num2str(i)],formats,renderer);
    
%     plotSolution(S,norm(u-u_ex),'rotation',1,'ampl',ampl,options{:});
%     mysaveas(pathname,['err_Rx_' num2str(i)],formats,renderer);
    
    plotSolution(S,u,'rotation',2,'ampl',ampl,options{:});
    mysaveas(pathname,['Ry_' num2str(i)],formats,renderer);
    
    plotSolution(S,u_ex,'rotation',2,'ampl',ampl,options{:});
    mysaveas(pathname,['Ry_ex_' num2str(i)],formats,renderer);
    
%     plotSolution(S,norm(u-u_ex),'rotation',2,'ampl',ampl,options{:});
%     mysaveas(pathname,['err_Ry_' num2str(i)],formats,renderer);
end

end

if displayCv
    figure(hcvUz)
    loglog(Nbelem,err_Uz,'LineStyle','-','Color',getfacecolor(ie+1),'LineWidth',linewidth);
    hold on
    
    figure(hcvRt)
    loglog(Nbelem,err_Rt,'LineStyle','-','Color',getfacecolor(ie+1),'LineWidth',linewidth);
    hold on
    
    figure(hcvRx)
    loglog(Nbelem,err_Rx,'LineStyle','-','Color',getfacecolor(ie+1),'LineWidth',linewidth);
    hold on
    
    figure(hcvRy)
    loglog(Nbelem,err_Ry,'LineStyle','-','Color',getfacecolor(ie+1),'LineWidth',linewidth);
    hold on
    
    figure(htime)
    loglog(Nbelem,time,'LineStyle','-','Color',getfacecolor(ie+1),'LineWidth',linewidth);
    hold on
end

end

if displayCv
    pathname = fullfile(getfemobjectoptions('path'),'BASIC',...
        'EXAMPLES','SHELLS',filename);
    
    figure(hcvUz)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of elements')
    ylabel('Error')
    %xlabel('Nombre d''éléments')
    %ylabel('Erreur')
    legend(leg{:})
    mysaveas(pathname,'error_Uz',formats,renderer);
    mymatlab2tikz(pathname,'error_Uz.tex');
    
    figure(hcvRt)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of elements')
    ylabel('Error')
    %xlabel('Nombre d''éléments')
    %ylabel('Erreur')
    legend(leg{:})
    mysaveas(pathname,'error_Rt',formats,renderer);
    mymatlab2tikz(pathname,'error_Rt.tex');
    
    figure(hcvRx)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of elements')
    ylabel('Error')
    %xlabel('Nombre d''éléments')
    %ylabel('Erreur')
    legend(leg{:})
    mysaveas(pathname,'error_Rx',formats,renderer);
    mymatlab2tikz(pathname,'error_Rx.tex');
    
    figure(hcvRy)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of elements')
    ylabel('Error')
    %xlabel('Nombre d''éléments')
    %ylabel('Erreur')
    legend(leg{:})
    mysaveas(pathname,'error_Ry',formats,renderer);
    mymatlab2tikz(pathname,'error_Ry.tex');
    
    figure(htime)
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Number of elements')
    ylabel('CPU time [s]')
    %xlabel('Nombre d''éléments')
    %ylabel('Temps CPU [s]')
    legend(leg{:})
    mysaveas(pathname,'cputime',formats,renderer);
    mymatlab2tikz(pathname,'cputime.tex');
end

end
end
