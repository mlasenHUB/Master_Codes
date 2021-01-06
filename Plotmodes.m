function [ ] = Plotmodes( c, d,V, Vnd, Vnb, t)
%Doing the plotting
%t= 1, use torsion modes else:0
N=length(Vnd(:,1));
Nb=length(Vnb(:,1));

if t==0 || t==1
%Lines
figure('units','normalized','outerposition',[0 0 1 1])
suptitle([num2str(Nb+N),' Mode Shapes', ' Disks: -o ','Blades: -x'  ])
bb=(N+Nb)/c;

for i=1:Nb+N

subplot(bb,c,i)       
stem(Vnd(:,i), 'r')
hold on
stem(Vnb(1:N,i),'r-x')%Left/up blades
stem(Vnb(N+1:Nb,i), 'b-x')%Right/down blades

yli=-(real(max(V(:))));
yls=(real(max(V(:))));
ylim([yli yls])
label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
title(label)
set(get(gca,'title'),'Position',[N/2 yls 0.0])
hold off
set(gca,'xtick',0:2*N)
end
legend('Disk','Base-blade', 'Tip-blade')

end

%Squares
%{
figure('units','normalized','outerposition',[0 0 1 1])
suptitle([num2str(2*N),' Mode Shapes', ' Disks: -o ','Blades: -x','. Positive clockwise, normalise to the max. amplitude'  ])
bb=2*N/c;

for i=1:2*N

subplot(bb,c,i) 

circle(2,Vnd(:,i),Vnb(:,i));
label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
title(label)
hold on
end
hold off
%}


%second mass-blade
if t==1
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle([num2str(Nb+N),' Mode Shapes', ' Disks: -o ','Blades: -x'  ])
    bb=(N+Nb)/c;
    %circles
    THETA=2*pi/N:(2*pi/N):(2*pi);
    theta=linspace(0, 2*pi,100);
    for i=1:Nb+N
    
    subplot(bb,c,i)

    x=sin(THETA);
    y=cos(THETA);
    stem3(x,y,[Vnd(:,i)],'r-o')
    hold on
    stem3(x,y,[Vnb(1:N,i)], 'r-x')%first mass-blade
    plot(sin(theta),cos(theta), 'r')
    r=1.2;%External radious
    stem3(r*x,r*y,[Vnb(N+1:end,i)],'b-x')
    %plot3(x,y,[Vnd(:,i)])
    plot(r*sin(theta),r*cos(theta), 'b')
    hold off
    end
    legend('Disk','Base-blade', 'Tip-blade')
end

if t==0 %can be model 1 or model 2
    
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle([num2str(Nb+N),' Mode Shapes', ' Disks: -o ','Blades: -x'  ])
    bb=(N+Nb)/c;
    %circles
    THETA=2*pi/N:(2*pi/N):(2*pi);
    theta=linspace(0, 2*pi,100);
    for i=1:Nb+N
        subplot(bb,c,i)

        x=sin(THETA);
        y=cos(THETA);
        stem3(x,y,[Vnd(:,i)],'r-o')
        hold on
        stem3(x,y,[Vnb(1:N,i)], 'r-x')%first mass-blade
        plot(sin(theta),cos(theta), 'r')

        if Nb > N
            stem3(x,y,[Vnb(N+1:end,i)],'b-x')
        else
            stem3(x,y,[Vnb(1:end,i)],'b-x')
        end
    %plot3(x,y,[Vnd(:,i)])
        hold off
    end
    legend('Disk','Base-blade', 'Tip-blade')
end


if t==2 % linear +torsional modes
  
    %Lines
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle([num2str(2*Nb+N),' Mode Shapes', ' Disks: -o ',' Blades Linear: -x ', ' Blades Rotation: -+'  ])
    bb=(N+Nb)/c;

    for i=1:Nb+N

        subplot(bb,c,i)       
        stem(Vnd(:,i), 'r')
        hold on
        stem(Vnb(1:N,i),'r-x')%Left/up blades linear
        stem(Vnb(N+1:2*N,i), 'b-x')%Right/down blades linear
        stem(Vnb(2*N+1:3*N,i),'r-s')%Left/up blades torsion
        stem(Vnb(3*N+1:end,i), 'b-s')%Right/down blades torsion
        

        yli=-(real(max(V(:))));
        yls=(real(max(V(:))));
        ylim([yli yls])
        label=(['F', num2str(d(i)^0.5), ' P', num2str(i)]);
        title(label)
        set(get(gca,'title'),'Position',[N/2 yls 0.0])
        hold off
        set(gca,'xtick',0:2*N)
    end
    legend('Disk','Base-blade', 'Tip-blade')
    
    
    %circles
    figure('units','normalized','outerposition',[0 0 1 1])
    suptitle([num2str(Nb+N),' Mode Shapes', ' Disks: -o ','Blades: -x'  ])
    bb=(N+Nb)/c;
    %circles
    THETA=2*pi/N:(2*pi/N):(2*pi);
    theta=linspace(0, 2*pi,100);
    for i=1:Nb+N
    
        subplot(bb,c,i)

        x=sin(THETA);
        y=cos(THETA);
        stem3(x,y,[Vnd(:,i)],'r-o')
        hold on
        stem3(x,y,[Vnb(1:N,i)], 'r-x')%first mass-blade linear
        stem3(x,y,[Vnb(N+1:2*N,i)],'b-x')%2nd mass-blade linear
        plot(sin(theta),cos(theta), 'r')
        r=1.2;%External radious
        
        stem3(r*x,r*y,[Vnb(2*N+1:3*N,i)],'r-s')%first mass-blade torsion
        stem3(r*x,r*y,[Vnb(3*N+1:end,i)],'b-s')%2nd mass-blade torsion
        
        plot(r*sin(theta),r*cos(theta), 'b')
        hold off
    end
    legend('Disk','Base-blade', 'Tip-blade')
    
    
end



end





