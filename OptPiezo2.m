function [ Mm, Mm2, Pt ,P] = OptPiezo2( Mmin )
%Finds the optimal position for 2 piezos suppressing various modes
%   Mm: deliver the minimum force value, pos of 1st and 2nd piezo
%   Mm2: eliminated the redundants positions

Ma=abs(Mmin);
[r, c, deep]=size(Ma);%number of row: piezo1 pos, c=piezo2 pos, deep=Modes
for i = 1:deep %modes
    aux=Ma(:,:,i);
    aux(logical(eye(r)))=nan;
    Ma(:,:,i) = aux;  %eliminate the trivial cases of 2 piezos: one over the other one
end
[min_val,idx]=min(Ma(:,:,:));



Mm=zeros(r,3,deep);


for k=1:deep %number of modes[
    
    m =min(min_val(:,:,k));
    count=1;
    for i=1:r %position of the first piezo (rows)
        for j=1:r% position of the second piezo (columns)
        
        if Ma(i,j,k) < m+1e-5 && Ma(i,j,k)> m-1e-5 
            Mm(count,1,k)= Ma(i,j,k);
            Mm(count,2,k)= i;
            Mm(count,3,k)=j;
            
            count=count+1;
        end
        end
    end
end

ll=length(Mm(:,1,1));
Mm2=Mm;
for k=1:deep
    for i=1:ll
        v1=Mm2(i,2,k);
        v2=Mm2(i,3,k);
        
        for j=1:ll
            if i~=j && Mm2(j,3,k)==v1 && Mm2(j,2,k)==v2
                Mm2(j,:,k)=0;
            end
        end
    end
end

Pp1=Mm2(:, 2,:);%just the coordinates for piezo 1
Pp1=Pp1(:);%in one row
Pp1=Pp1(Pp1~=0);%get rid of 0 values

Pp2=Mm2(:, 3,:);%just the coordinates for piezo 2
Pp2=Pp2(:);%in one row
Pp2=Pp2(Pp2~=0);%get rid of 0 values


Pt=zeros(deep,3);%Matriz of both piezo positions and force
Pp2aux=Pp2;

for n=1:r
    Pp2=Pp2aux;
    for m=1:r
        for k=1:deep
            
            %count=1;
            [pmp2, frepp2]=mode(Pp2);
            [pmp1, frepp1]=mode(Pp1);

            for i =1:ll
                if pmp1==Mm2(i,2,k) && pmp2==Mm2(i,3,k) && Pt(k,1)==0
                    Pt(k,1)=pmp1;
                    Pt(k,2)=pmp2;
                    Pt(k,3)=Mm2(i,1,k);
                end
            end
            %count=count+2;
        end
        Pp2=Pp2(Pp2~=pmp2);

        if isempty(Pp2)
           break
        end
    end
    Pp1=Pp1(Pp1~=pmp1);

    if isempty(Pp1)
       break
    end
end
P=Pt(:,1:2,:);
P=P(:);
P=unique(P);%just the position of all piezos.
     

end

            
        
        
    




