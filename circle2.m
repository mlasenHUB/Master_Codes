function H=circle(radius,Vnd, Vnb)
NOP=length(Vnd);

if length(Vnd)~= length(Vnb)
    error('Different size of blades and disks');    
end

%NORMALIZE VECTOR
N1=max(abs(Vnd));
N2=max(abs(Vnb));
N=max(N1,N2);
NVnd=Vnd./N;
NVnb=Vnb./N;

THETA=linspace(0,2*pi,NOP+1);
RHOd=ones(1,NOP+1)*radius;%disks radious.
RHOb=ones(1,NOP+1)*2*radius;%blades radious.

[Xd,Yd] = pol2cart(THETA,RHOd);
[Xb,Yb] = pol2cart(THETA,RHOb);

Di=plot(Xd,Yd,'r');
hold on
Bl=plot(Xb,Yb,'b');
hold on

%For Plot numbers on the points
d=[1:NOP]';
b=[NOP+1:2*NOP]';
dn = num2str(d);
bn = num2str(b);
for i=1:NOP
if i<=NOP/2   
dc{i} = cellstr(dn(i));
else
bc{i} = cellstr(bn(i));
end
dx = 0.1; dy = 0.1;
end

%
for i=1:NOP
    
plot([Xb(i), Xb(i)+NVnb(i)*sin(THETA(i))],[Yb(i),Yb(i)-NVnb(i)*cos(THETA(i))], 'b-x');%Blades

text(Xb(i)+dx, Yb(i)+dy, cellstr(num2str(i+NOP)));
plot([Xd(i), Xd(i)+NVnd(i)*sin(THETA(i))],[Yd(i),Yd(i)-NVnd(i)*cos(THETA(i))], 'r-o');%Disks

text(Xd(i)+dx, Yd(i)+dy, cellstr(num2str(i))); 


end

hold off
r=2*radius*1.25;
axis([-r r -r r])
end
