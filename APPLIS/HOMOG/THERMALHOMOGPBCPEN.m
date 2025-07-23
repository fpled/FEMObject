classdef THERMALHOMOGPBCPEN
    % THERMALHOMOGPBC class
    %   Homogenization of thermal conductivity using PERIODIC boundary
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
        macro_mode_unfree;
        macro_mode_free;
        u_unfree_tot;
        u_unfree;
        u_free
        C;
    end
    
    methods
        function obj=THERMALHOMOGPBCPEN(size_domain,nb_elem,ls,neg_side,pos_side)
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
            
            obj=obj.buildmacromodes();
            
            obj=obj.buildoperators();
            
            obj.u_unfree_tot=cell([obj.dim 1]);
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
            f=obj.A_unfree*obj.macro_mode_unfree{dim};
            f=freevector(obj.model,f);
            f0=calc_nonhomogeneous_vector(obj.model,obj.A_unfree);
            
            %penalization
            lx = setfree( LINFORM(0,1),0);
            Lx = calc_vector(lx,obj.model);
            obj.u_free{dim}=solvesingular(obj.A_free+epsi*freematrix(Lx*Lx',obj.model),f-f0);
            obj.u_unfree{dim}=unfreevector(obj.model,obj.u_free{dim});
            obj.u_unfree_tot{dim}=obj.macro_mode_unfree{dim}-obj.u_unfree{dim};
        end
        
        function obj=findeffectiveconductivity(obj)
            for i=1:obj.dim
                for j=1:obj.dim
                    obj.C(i,j)=obj.u_unfree_tot{i}'*obj.A_unfree*obj.u_unfree_tot{j}/getvolume(obj.domain);
                end
            end
        end
        
        function []=plot_mode_micro(obj,dim,varargin)
            dim_str=['x' 'y' 'z'];
            figure('Name',['Mode micro ' dim_str(dim)]);
            plot(obj.u_unfree{dim},obj.model,varargin{:});
            colorbar;
        end
        
        function []=plot_mode_macro(obj,dim,varargin)
            dim_str=['x' 'y' 'z'];
            figure('Name',['Mode macro ' dim_str(dim)]);
            plot(obj.macro_mode_unfree{dim},obj.model,varargin{:});
            colorbar;
        end
        
        function []=plot_mode_tot(obj,dim,varargin)
            dim_str=['x' 'y' 'z'];
            figure('Name',['Mode total ' dim_str(dim)]);
            plot(obj.u_unfree_tot{dim},obj.model,varargin{:});
            colorbar;
        end
        
        function []=plot_ls(obj,varargin)
            figure('Name','Level-set');
            plot(obj.ls,obj.model,varargin{:});
            colorbar;
        end
        
    end
    methods (Access = private)
        function obj=buildmodels(obj)
            obj.domain=DOMAIN(obj.dim, zeros(1,obj.dim), obj.size_domain);
            if obj.dim == 2
                obj.model=mesh(obj.domain,obj.nb_elem(1),obj.nb_elem(2));
            else
                obj.model=mesh(obj.domain,obj.nb_elem(1),obj.nb_elem(2),obj.nb_elem(3));
            end
            obj.model=createddlnode(obj.model,DDL('T'));
        end
        
        function obj=definebcs(obj)
            if obj.dim==2
                obj.model=addcl(obj.model,POINT([0 0; obj.size_domain(1) 0;0 obj.size_domain(2);obj.size_domain(1) obj.size_domain(2)]),'T');
                obj.model=addclperiodic(obj.model,getedge(obj.domain,1),getedge(obj.domain,3),'T');
                obj.model=addclperiodic(obj.model,getedge(obj.domain,2),getedge(obj.domain,4),'T');
            elseif obj.dim==3
                c1=[0 obj.size_domain(1)];
                c2=[0 obj.size_domain(2)];
                c3=[0 obj.size_domain(3)];
                for i=1:2
                    for j=1:2
                        for k=1:2
                            obj.model=addcl(obj.model,POINT([c1(i) c2(j) c3(k)]),'T');
                        end
                    end
                end
                obj.model = addclperiodic(obj.model,getface(obj.domain,1),getface(obj.domain,3),'T');
                obj.model = addclperiodic(obj.model,getface(obj.domain,2),getface(obj.domain,4),'T');
                obj.model = addclperiodic(obj.model,getface(obj.domain,5),getface(obj.domain,6),'T');
            end
        end
        
        function obj=processls(obj,ls,neg_side,pos_side)
            ls=lseval(ls,obj.model);
            ls=getvalue(ls);
            obj.ls=neg_side+(pos_side-neg_side)*double(ls>0);
        end
        
        function obj=buildmacromodes(obj)
            obj.macro_mode_unfree=cell(obj.dim,1);
            obj.macro_mode_free=cell(obj.dim,1);
            if obj.dim==2
                ls=LSHYPERPLANE(2,0,0,1,0);
                ls=lseval(ls,obj.model);
                obj.macro_mode_unfree{1}=getvalue(ls);
                obj.macro_mode_free{1}=freevector(obj.model,obj.macro_mode_unfree{1});
                ls=LSHYPERPLANE(2,0,0,0,1);
                ls=lseval(ls,obj.model);
                obj.macro_mode_unfree{2}=getvalue(ls);
                obj.macro_mode_free{2}=freevector(obj.model,obj.macro_mode_unfree{2});
            else
                ls=LSHYPERPLANE(3,0,0,0,1,0,0);
                ls=lseval(ls,obj.model);
                obj.macro_mode_unfree{1}=getvalue(ls);
                obj.macro_mode_free{1}=freevector(obj.model,obj.macro_mode_unfree{1});
                ls=LSHYPERPLANE(3,0,0,0,0,1,0);
                ls=lseval(ls,obj.model);
                obj.macro_mode_unfree{2}=getvalue(ls);
                obj.macro_mode_free{2}=freevector(obj.model,obj.macro_mode_unfree{2});
                ls=LSHYPERPLANE(3,0,0,0,0,0,1);
                ls=lseval(ls,obj.model);
                obj.macro_mode_unfree{3}=getvalue(ls);
                obj.macro_mode_free{3}=freevector(obj.model,obj.macro_mode_unfree{3});
            end
        end
        
        function obj=buildoperators(obj)
            a=BILINFORM(1,1,obj.ls(:),0);
            a=setfree(a,0);
            obj.A_unfree=a{obj.model}(:,:);
            obj.A_free=freematrix(obj.model,obj.A_unfree);
        end
        
        
    end
end

