function [ MM, Maor ] = MistuningMass( Ma, pm , dofm, N, sine, per, leftb, alt )
%Apply a mass mistuning in a certain degree of freedom

%  pm= percentage of change of M
%  dofM= degree of freedom to be changed on M
%  N= Number of dof
%  Maor= Original Ma matrix
%  leftb=1, apply a sinusoidal mistuning to the LEFT/UPPER blades (0: the right/down)
%  ones), if there is just one blade then that blade is LEFT

Maor=Ma;
ma=diag(Ma);
n=length(ma);

Mad=diag(Maor);




if alt==1
    
    for i=1:N
        if leftb==1
            Mad(N+i)=Mad(N+i)*(1+((-1)^(i+1))*pm/100);
        end
        if leftb==0
            Mad(2*N+i)=Mad(2*N+i)*(1+((-1)^(i+1))*pm/100);
        end
    end
end
if alt==0
    if dofm>n
        error(' DOF out of range.(Mass)')
    end
    if dofm<=N
        error(' Not a blade (Mass).')
    end
    
    Ma(dofm, dofm)=Ma(dofm, dofm)*(1+pm/100);
    Mad=diag(Ma);
    
    if sine==1
        x=linspace(0, 2*pi, N);
        y= (pm/100)*sin(per*x);%
        
        if leftb==1
            Mad(N+1:2*N)=Mad(N+1:2*N).*(1+y)';
        else
            Mad((2*N)+1:end)=Mad((2*N)+1:end).*(1+y)';
        end
        
        
    end
end
MM=diag(Mad);    
end

