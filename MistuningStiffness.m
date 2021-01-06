function [ kb, kbor ] = MistuningStiffness( kb, pk, dofk, N, sine, per, alt )
%Apply a stiffness mistuning in a certain degree of freedom

%  pk= percentage of change of K (if -ve, then reduce stiffness)
%  dofk= degree of freedom to be changed on K
%  N= Number of dof(between 1 to N blades
%  kbor= Original K matrix Tuned

%  for a sin wave mistuning
%  sin=1 then apply sinusoidal mistuning
%  pe= Period of the sin, or nodal diameter
%  leftb=1, apply a sinusoidal excitation to the LEFT blades (0: the right
%  ones)

kbor=kb;

if dofk>N
    error(' DOF out of range.(Stiffness: N blades)')
end



if alt==1
    kb=kbor;
    for i=1:N
        kb(i)=kb(i)*(1+((-1)^(i+1))*pk/100);
    end
else
    kb(dofk)= kbor(dofk)*(1+pk/100);
    
    if sine==1
        x=linspace(0, 2*pi, N);
        y= (pk/100)*sin(per*x);
        kb=kbor.*(1+y)';
        
    end
end
    

end

