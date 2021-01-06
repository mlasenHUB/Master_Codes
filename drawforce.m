function [ ] = drawforce(fs)
%DRAW THE FORCE IN A SCHEMATIC SYSTEM OF 2 DOFS PER BLADE PLUS DISKS
%   Detailed explanation goes here
n=length(fs);
%basic representation

line([-3 3], [0 0], 'Color','black')
xlim([-3.5 3.5])
ylim([-3.5 3.5])
hold on
line([0 0], [-3 3], 'Color','black')

line([-0.5 0.5], [1 1], 'Color','black')
line([-0.5 0.5], [2 2], 'Color','black')
line([-0.5 0.5], [3 3], 'Color','black')
line([-0.5 0.5], [-1 -1], 'Color','black')
line([-0.5 0.5], [-2 -2], 'Color','black')
line([-0.5 0.5], [-3 -3], 'Color','black')

line([1 1], [-0.5 0.5], 'Color','black')
line([2 2],[-0.5, 0.5], 'Color','black')
line([3 3],[-0.5, 0.5], 'Color','black')
line([-1 -1], [-0.5 0.5], 'Color','black')
line([-2 -2],[-0.5, 0.5], 'Color','black')
line([-3 -3],[-0.5, 0.5], 'Color','black')

for i=1:n
    if i==1 && fs(i)>0
        rectangle('Position', [-0.5 0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==1 && fs(i)<0
        rectangle('Position', [-0.5 0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==2 && fs(i)>0
        rectangle('Position', [0.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==2 && fs(i)<0
        rectangle('Position', [0.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==3 && fs(i)>0
        rectangle('Position', [-0.5 -1.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==3 && fs(i)<0
        rectangle('Position', [-0.5 -1.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==4 && fs(i)>0
        rectangle('Position', [-1.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==4 && fs(i)<0
        rectangle('Position', [-1.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    
    
    %%%
    if i==5 && fs(i)>0
        rectangle('Position', [-0.5 1.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==5 && fs(i)<0
        rectangle('Position', [-0.5 1.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==6 && fs(i)>0
        rectangle('Position', [1.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==6 && fs(i)<0
        rectangle('Position', [1.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==7 && fs(i)>0
        rectangle('Position', [-0.5 -2.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==7 && fs(i)<0
        rectangle('Position', [-0.5 -2.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==8 && fs(i)>0
        rectangle('Position', [-2.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==8 && fs(i)<0
        rectangle('Position', [-2.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    
    %%%%
    
    
    if i==9 && fs(i)>0
        rectangle('Position', [-0.5 2.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==9 && fs(i)<0
        rectangle('Position', [-0.5 2.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==10 && fs(i)>0
        rectangle('Position', [2.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==10 && fs(i)<0
        rectangle('Position', [2.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==11 && fs(i)>0
        rectangle('Position', [-0.5 -3.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==11 && fs(i)<0
        rectangle('Position', [-0.5 -3.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end
    if i==12 && fs(i)>0
        rectangle('Position', [-3.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [0 0 0])
    end
    if i==12 && fs(i)<0
        rectangle('Position', [-3.5 -0.5 1 1], 'Curvature', 1, 'FaceColor', [1 1 1])
    end 
        
        
end
    

