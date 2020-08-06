% Algorithm for IHCP, implemented in Matlab. Computed innerWall temperature history is saved into a file

clear;
temp_ref = 100.0; %reference temperature of the simulation run with ANSYS
temp_ref_inv = 1.0/100.0;

OutWall_temp = importdata('/data/OuterWallTempHistory_Perturbed.txt'); % returns a matrix after reading from file
OutT_of_UT = importdata('/data/Ref_Out_Temp.txt');

[Time,x] = size(OutWall_temp);        % 'time' is the total time for the simulation, size() returns nos of rows or columns from a matrix
[responseTimeOfUT,x] = size(OutT_of_UT);
[x,Layer]= size(OutWall_temp);
Layer=Layer-1;

%setting initial condition
initialTemp(Layer)=0;
for i = 1 : Layer
    initialTemp(i)=OutWall_temp(1,i+1);
end

%required matrices memory allocation
UT_influenced_outwall_T(Time,Layer*Layer,Time)= 0;  %3D matrix for storing calculated valued per timestep 
UT_influenced_outwall_T(:,:,:)= 0; %initializing with 0
DelInnWall_temp(Time-1,Layer)=0;
DelInnWall_temp(:,:)=0; %2D matrix, for storing final results

%constructing matrix A from Ax=b
A(Layer,Layer)=0;
for i = 1 : Layer
    for j = 1 : Layer
        A(i,j) = OutT_of_UT(2,j+i*Layer-Layer+1);
    end
end

%initializing b with 0
b(Layer)=0;

%computations starts        
for i = 1 : Time-1 %3rd dimension of the 3D matrix
    for j = 1 : Layer %Constructing b starts
        temporary=0.0;
        if i>1 
            for k = 1 : i-1
                for l = 1 : Layer
                    temporary = temporary + UT_influenced_outwall_T(i+1,l*Layer-Layer+j,k);
                end
            end
        end
                
        b(j)=(OutWall_temp(i+1,j+1)-initialTemp(j)- temporary)*temp_ref; 
    end
    
    T = A\b';	%solver will automatically be selected, solving Ax=b
    
    %storing solutions
    for m = 1:Layer
        
        DelInnWall_temp(i,m)= T(m);
    end
       
    %forward algorithm starts
    for j = 1 : Layer*Layer
        
        for k = i+1 : i+responseTimeOfUT-1
            if k <= Time 
                if Layer>1 
                    UT_influenced_outwall_T(k,j,i)= OutT_of_UT(k-i+1,j+1) *DelInnWall_temp(i,floor((j-1)/Layer+1))* temp_ref_inv ;
                end
            end
        end
    end
    
end

%setting initial condition to the result matrix
InnWall_Temp(Time,Layer)=0;
for it = 1 : Time
    for jt = 1 : Layer
        InnWall_Temp(it,jt)=initialTemp(jt);
    end
end

%completing resultant matrix with each second changed Temperature
for m = 2 : Time
    for n =1 : Layer
        InnWall_Temp(m,n) = DelInnWall_temp(m-1,n);
        InnWall_Temp(m,n)= InnWall_Temp(m,n)+ InnWall_Temp(m-1,n);
    end
end

%Writing into a text file
save('/data/InnWallTemp_AlgorithmComputed.txt','InnWall_Temp','-ASCII');
