classdef THERMALHOMOGNBCPEN
    % THERMALHOMOGNBC class
    %   Homogenization of thermal conductivity using NATURAL boundary
    %   conditions, the rigid body motion is blocked via penalization
    properties
        size_domain;
        nb_elem;
        dim;
        domain;
        model;
        ls;
        A_unfree;
        A_free;
        u_unfree;
        u_free;
        f;
        C;
    end
    
    methods
        function obj=THERMALHOMOGNBCPEN(size_domain,nb_elem,ls,neg_side,pos_side)
            obj.dim=length(size_domain);
            obj.size_domain=size_domain;
            if (obj.dim ~= 2) && (obj.dim ~=3)
                error('Wrong number of dimensions in size_domain')
            end
            if (length(size_domain)~=length(nb_elem)) && (length(nb_elem)~=1)
                error('Error in the number of dimensions in nb_elem')
            end
            if length(nb_elem)==1
                obj.nb_elem=nb_elem*eyes(obj.dim,1);
            else
                obj.nb_elem=nb_elem;
            end
            
            obj=obj.buildmodels();
            obj=obj.definebcs();
            
            obj=obj.processls(ls,neg_side,pos_side);
            
            obj=obj.buildoperators();
            
            obj.u_unfree=cell([obj.dim 1]);
            obj.u_free=cell([obj.dim 1]);
            obj.C=zeros(obj.dim);
        end
        
        function obj=solve(obj,epsi)
            for i=1:obj.dim
                obj=obj.solve_mode(i,epsi);
            end
        end
        
        function obj=solve_mode(obj,dim,epsi)
            f0=calc_nonhomogeneous_vector(obj.model{dim},obj.A_unfree);
            %penalization
            lx = setfree( LINFORM(0,1),0);
            Lx = calc_vector(lx,obj.model{dim});
            obj.u_free{dim}=solvesingular(obj.A_free+epsi*freematrix(Lx*Lx',obj.model{dim}),obj.f{dim}-f0);
            obj.u_unfree{dim}=unfreevector(obj.model{dim},obj.u_free{dim});
        end
        
        function obj=findeffectiveconductivity(obj)
            for i=1:obj.dim
                for j=1:obj.dim
                    obj.C(i,j)=obj.u_unfree{i}'*obj.A_unfree*obj.u_unfree{j};
                end
            end
            obj.C=inv(obj.C/getvolume(obj.domain));
        end
        
        function []=plot_mode(obj,dim,varargin)
            dim_str=['x' 'y' 'z'];
            figure('Name',['Mode ' dim_str(dim)]);
            plot(obj.u_unfree{dim},obj.model{dim},varargin{:});
            colorbar;
        end
        
        function []=plot_ls(obj,varargin)
            figure('Name','Level-set');
            plot(obj.ls,obj.model{1},varargin{:});
            colorbar;
        end
    end
    
    methods (Access = private)
        function obj=buildmodels(obj)
            obj.domain=DOMAIN(obj.dim, zeros(1,obj.dim), obj.size_domain);
            obj.model=cell(obj.dim,1);
            if obj.dim == 2
                S=mesh(obj.domain,obj.nb_elem(1),obj.nb_elem(2));
            else
                S=mesh(obj.domain,obj.nb_elem(1),obj.nb_elem(2),obj.nb_elem(3));
            end
            for i=1:obj.dim
                obj.model{i}=S;
                obj.model{i}=createddlnode(obj.model{i},DDL('T'),DDL('QN'));
            end
        end
        
        function obj=definebcs(obj)
            point=obj.size_domain/2;
            %Neumann
            if obj.dim==2
                P=cell(4,1);
                L=cell(4,1);
                
                P{1}=POINT([0 0]);
                P{2}=POINT([obj.size_domain(1) 0]);
                P{3}=POINT([obj.size_domain(1) obj.size_domain(2)]);
                P{4}=POINT([0 obj.size_domain(2)]);
                
                L{1}=LINE(P{1},P{2});
                L{2}=LINE(P{2},P{3});
                L{3}=LINE(P{3},P{4});
                L{4}=LINE(P{4},P{1});
                
                obj.f{1}=surfload(obj.model{1},L{2},'QN',1);
                obj.f{1}=obj.f{1}+surfload(obj.model{1},L{4},'QN',-1);
                obj.f{2}=surfload(obj.model{2},L{1},'QN',-1);
                obj.f{2}=obj.f{2}+surfload(obj.model{2},L{3},'QN',1);
            elseif obj.dim==3
                P=cell(8,1);
                H=cell(6,1);
                
                x=obj.size_domain(1);
                y=obj.size_domain(2);
                z=obj.size_domain(3);
                
                P{1}=POINT([0 0 0]);
                P{2}=POINT([0 0 z]);
                P{3}=POINT([0 y 0]);
                P{4}=POINT([0 y z]);
                P{5}=POINT([x 0 0]);
                P{6}=POINT([x 0 z]);
                P{7}=POINT([x y 0]);
                P{8}=POINT([x y z]);
                
                H{1}=PLANE(P{1},P{2},P{3});
                H{4}=PLANE(P{5},P{6},P{7});
                
                H{2}=PLANE(P{1},P{2},P{5});
                H{5}=PLANE(P{3},P{4},P{7});
                
                H{3}=PLANE(P{1},P{3},P{5});
                H{6}=PLANE(P{2},P{4},P{6});
                
                for i=1:obj.dim
                    obj.f{i}=surfload(obj.model{i},H{i},'QN',-1)+surfload(obj.model{i},H{i+3},'QN',1);
                end
            end
        end
        
        function obj=processls(obj,ls,neg_side,pos_side)
            ls=lseval(ls,obj.model{1});
            ls=getvalue(ls);
            obj.ls=neg_side+(pos_side-neg_side)*double(ls>0);
        end
        
        function obj=buildoperators(obj)
            a=BILINFORM(1,1,obj.ls(:),0);
            a=setfree(a,0);
            obj.A_unfree=a{obj.model{1}}(:,:);
            obj.A_free=freematrix(obj.model{1},obj.A_unfree);
            
        end
        
        
    end
end

