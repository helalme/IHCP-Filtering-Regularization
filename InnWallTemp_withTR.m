clear;
%stacksize('max');
temp_ref = 100.0; %reference temperature of the simulation run with ANSYS
temp_ref_inv = 1.0/100.0;

OutWall_temp = importdata('/data/OuterWallTempHistory_Perturbed.txt'); %fscanfMat*() returns a matrix after reading
OutT_of_UT = importdata('/data/Ref_Out_Temp.txt');

[Time,x] = size(OutWall_temp);        % 'time' is the total time for the simulation, size() returns nos of rows or columns from a matrix
[responseTimeOfUT,x] = size(OutT_of_UT);
[x,Layer]= size(OutWall_temp);
Layer=Layer-1;

initialTemp(Layer)=0;
for i = 1 : Layer
    initialTemp(i)=OutWall_temp(1,i+1);
end
UT_influenced_outwall_T(Time,Layer*Layer)= 0;  %3D matrix for storing calculated valued per timestep 
UT_influenced_outwall_T_sum(Time,Layer*Layer)= 0;
UT_influenced_outwall_T(:,:,:)= 0; %initializing with 0
DelInnWall_temp(Time-1,Layer)=0;
DelInnWall_temp(:,:)=0; %2D matrix, for storing final results

A(Layer,Layer)=0;
for i = 1 : Layer
    for j = 1 : Layer
        A(i,j) = OutT_of_UT(2,j+i*Layer-Layer+1);
    end
end

b(Layer)=0;
alf_list(Time-1,1)=0;       
for i = 1 : Time-1
    for j = 1 : Layer %to calculate temperature of a layer at time i
        temporary=0.0;
        if i>1 
            for l = 1 : Layer
                temporary = temporary + UT_influenced_outwall_T_sum(i+1,l*Layer-Layer+j);
            end
        end
                
        b(j)=(OutWall_temp(i+1,j+1)-initialTemp(j)- temporary)*temp_ref; %-OutWall_temp(i,j+1) is it necessary, try to use it, it could reduce computations
    end
    
   %T = A\b';%linsolve returns a warning and does not return a solution, pcg also returns some NaN values so A|b is better
  [T,alf]=TikhonovReg(A,b',i);

  alf_list(i)= alf;

    for m = 1:Layer
        DelInnWall_temp(i,m)= T(m);
    end
       
    UT_influenced_outwall_T(:,:)=0;
    for j = 1 : Layer*Layer
        
        for k = i+1 : i+responseTimeOfUT-1
            if k <= Time 
                if Layer>1 
                    UT_influenced_outwall_T(k,j)= OutT_of_UT(k-i+1,j+1) *DelInnWall_temp(i,floor((j-1)/Layer+1))* temp_ref_inv ;
                end
            end
        end
    end
    
    for it = 1 : Time
        for jt = 1 : Layer*Layer
            UT_influenced_outwall_T_sum(it,jt)= UT_influenced_outwall_T_sum(it,jt)+ UT_influenced_outwall_T(it,jt);
        end
    end
    
end

InnWall_Temp(Time,Layer)=0;
for it = 1 : Time
    for jt = 1 : Layer
        InnWall_Temp(it,jt)=initialTemp(jt);
    end
end

for m = 2 : Time
    for n =1 : Layer
        InnWall_Temp(m,n) = DelInnWall_temp(m-1,n);
        InnWall_Temp(m,n)= InnWall_Temp(m,n)+ InnWall_Temp(m-1,n);
    end
end

save('/data/InnWallTemp_withTikhonovReg.txt','InnWall_Temp','-ASCII');
%save('/data/optimum_alpha.txt','alf_list','-ASCII');
