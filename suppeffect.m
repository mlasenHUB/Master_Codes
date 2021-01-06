function [ Vsup, Vnsup, y2, y3 ] = suppeffect(Vs, Vns, DAMP, PHs, PHns )
%UNTITLED Summary of this function goes here
%   Compares the change from matrices Vs to Vns in columns(modes) and
%   nodes(rows)
% R is the ratio of change with and without suppression.
% y3 displays the merged Vs and Vns.
%R=[];
%R=Vs./Vns;
%Vs and Vns = Suppressed/UnSuppressed displacement at each resonant frequency
n=length(Vs(1,:));

if DAMP==1 %If there is DAMPING then consider the phase to see the direction of the displacement
  Vsup=Vs.*sign(sin(deg2rad(PHs)));
  Vnsup=Vns.*sign(sin(deg2rad(PHns)));
end

if DAMP==0
    Vsup=Vs;
    Vnsup=Vns;
end

y3=[];

for i=1:n
y2=[];%clean the vector of comparison y2 every iteration
for j=1:n
    y=[Vsup(j,i) Vnsup(j,i) ];
    y2= [y2 ; y];
end
y3=[y3, y2];

end

