syms x(t)
m=1;
c=50;
k=1;
F1=1;
F2=0;
f1=1;
f2=1.5;
eqn = m*diff(x,t,2)+c*diff(x,t,1)+k*x == F1*sin(f1*t)+F2*sin(f2*t);
%cond1 = x(0) == 1;
%cond2 = Dy(0) == 0;

%conds = [cond1 cond2];
%ySol(x) = dsolve(ode,conds);
%ySol = simplify(ySol)
%xSol(t) = dsolve(eqn);
V = odeToVectorField(eqn);
M = matlabFunction(V,'vars', {'t','Y'});

interval = [0 300];
x0 = [0 0];
xSol = ode45(M,interval,x0);

tValues = linspace(0,300,1003);
xValues = deval(xSol,tValues,1);

plot(tValues,xValues)

xlabel('Time (s)')
ylabel('Displacement (m)')