function Indices=PCbase_indices(M,p,typebase)
% function Indices=PCbase_indices(M,p,typebase)
% M variables aleatoires independantes
% p ordre du chaos polynomial

if nargin<3
typebase = 1 ;
end

NbRva = M ; OrderPol = p;
MM = NbRva -1 ;
n = 1;
switch typebase
case 1
P=round(prod(p+1:M+p)/factorial(M)-1);

Indices=zeros(P+1,M+1);
case 2
P=(p+1)^M-1;
Indices=zeros((p+1)^M,M+1);
end

% Particular case 
if (NbRva > 0)

switch NbRva
	
 case 1 % One-dimensional Polynomial Chaos <=> Hermite Polynomials
	for i = 1 : OrderPol
	  Indices(i+1,1) = i;
	  Indices(i+1,2) = i;
    end
	
 otherwise % Higher dimensional Polynomial Chaoses
	
	switch typebase
	case 1
	for CurrentOrder = 1 : OrderPol 
	  EndGenere = 0;
	  FirstThisOrder = 0; 
	  
	  while (EndGenere == 0)
		n = n +1 ;
		% First list t for order CurrentOrder
		if (FirstThisOrder == 0) 
		  for i=1 : MM
			t(i) = i;
          end
		  FirstThisOrder = 1; 
		else
		  % Regular incrementation
		  if (t(MM) < (MM + CurrentOrder))
			t(MM) = t(MM) + 1;  
		  else  % t(MM) = tmax = MM + CurrentOrder
			j = MM;
			while (t(j) == j + CurrentOrder)
			  j = j - 1;
            end
			t(j) = t(j) + 1 ;
			for k =(j + 1) :  MM 
			  t(k) = t(j) + k - j;
            end
          end
        end
		
		% Direct Translating t into PsiBasis{n}
		Indices(n,M) = t(1) -1;
		for i=2 : MM 
		  Indices(n,M+1-i) = t(i) - t(i-1) -1;
        end
		Indices(n,1) = NbRva + CurrentOrder - t(MM) -1;
		Indices(n,M+1) = CurrentOrder;
		% End of generation of order CurrentOrder
		if (t(1) == (CurrentOrder+1))
		  EndGenere = EndGenere + 1;
        end
      end
    end
     case 2
     
    for k=1:M
	for j=0:p
	for l=1:(p+1)^(k-1)
	Indices(l+j*(p+1)^(k-1):(p+1)^(k):end,k)=j;
	end
	end
	end
	Indices(:,end)=max(Indices(:,1:M)')';
	Indices=sortrows(Indices,[size(Indices,2):-1:1]);

     end

end
else
    Indices=zeros(0,1);
end
  
