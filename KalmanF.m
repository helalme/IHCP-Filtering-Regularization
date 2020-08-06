% This script takes a noisy data, e.g. computed solution of IHCP, and applies Kalman filter to get smooth solution.

clear;

NoisyProfile = importdata('/data/InnWallTemp_AlgorithmComputed.txt'); % returns a matrix after reading from file
[Time,Layer] = size(NoisyProfile);        % 'time' is the total time for the simulation, size() returns nos of rows or columns from a matrix
Layer=Layer-1;

%required matrices memory allocation
OriginalProfile= NoisyProfile;  

%computation of MA
for i = 1 :Layer %will start from 2nd column
    for j = 1 : Time-1 % row
        if j>1 
            mean = (NoisyProfile(j-1,i+1)+NoisyProfile(j,i+1)+NoisyProfile(j+1,i+1))/3.0;
            NoisyProfile(j,i+1)= mean;
        end
        if j==Time-1
            NoisyProfile(j+1,i+1)=NoisyProfile(j,i+1);
        end
    end
end

for i = 1 :Layer %will start from 2nd column
    for j = 1 : Time % row
        if j>1 
            varianceCurr=(NoisyProfile(j,i+1)-OriginalProfile(j,i+1))* (NoisyProfile(j,i+1)-OriginalProfile(j,i+1));
            variancePrev=(NoisyProfile(j-1,i+1)-OriginalProfile(j-1,i+1))* (NoisyProfile(j-1,i+1)-OriginalProfile(j-1,i+1));
            Kgain=varianceCurr/(varianceCurr+variancePrev);
            NoisyProfile(j,i+1)= NoisyProfile(j,i+1)+ Kgain*(OriginalProfile(j,i+1)-NoisyProfile(j,i+1));
           
        end
    end
end


%Writing into a text file
save('/data/InnWallTemp_KalmanFiltered.txt','NoisyProfile','-ASCII');

