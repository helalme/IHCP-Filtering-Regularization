% This Matlab script adds perturbation to the original simulation results.
% i.e., perturbation to the outerWall temperature.
% Because sensor-measured temperature history (outerWall temperature) contains noises and perturbation   
clear;

InnWall_temp = importdata('/data/InnerWallTempHistory.txt');
OutWall_temp = importdata('/data/OuterWallTempHistory.txt');

[Time,Layer] = size(OutWall_temp);        

InnWall_temp_2s((Time+1)/2,Layer)=0;
OutWall_temp_2s((Time+1)/2,Layer)=0;

for i = 1 : Time
    for j = 1 : Layer 
        if (rem(i+1,2)==0)
            InnWall_temp_2s(floor((i+1)/2),j)=InnWall_temp(i,j);
            OutWall_temp_2s(floor((i+1)/2),j)=OutWall_temp(i,j);
        end
    end
end

OutWall_temp_2s(:,1)=[];
Prtrb = randint((Time+1)/2,Layer-1,[-5,5]);
Prtrb1=0.01*Prtrb;
Prtb_OutWall_temp_2s = OutWall_temp_2s + Prtrb1;

%save('ANSYS_InnT_2s.txt','InnWall_temp_2s','-ASCII');
save('/data/OuterWallTempHistory_Perturbed.txt','Prtb_OutWall_temp_2s','-ASCII');