
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves the results from computation for different geometric constraints. 
% Input: Results (array of all the computed results)
%        Constraint (name of the goemetric constraint
%        groundcondition (describes the status of body J)
% Output: File with all the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function file = saveFile(Results, constraint,groundCondition,fileName)

 fileID = fopen([fileName '.txt'],'w');

% str1 = Results;
% fileID = fopen('Results.txt', 'w');

fprintf(fileID,'- - - - - - - - - - - - \n');

fprintf(fileID,'The bodies I and J have a %s constraint. \n', constraint);
fprintf(fileID,'The body J is %s. \n', groundCondition);

if size(Results,2) ==7

constraint = Results{1,1};
Velocity = Results{1,2};
Acceleration = Results{1,3};
PartialPhi_ri= Results{1,4};
PartialPhi_rj= Results{1,5};
PartialPhi_pi= Results{1,6};
PartialPhi_pj= Results{1,7};

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
% fprintf(fileID,'The value of the expression of the constraint is %f \n', constraint);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The right-hand side of the velocity equation(mu) is %.5f \n', Velocity);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The right-hand side of the acceleration equation(gamma) is %.5f \n', Acceleration);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The expression for partial derivatives phi_R:\n');
fprintf(fileID, 'phi_Ri = %d.\n', PartialPhi_ri);
fprintf(fileID, 'phi_Rj = %d.\n', PartialPhi_rj);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The expression for partial derivatives phi_P:\n');
fprintf(fileID, 'phi_Ri = %d.\n', PartialPhi_pi);
fprintf(fileID, 'phi_Rj = %d.\n', PartialPhi_pj);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
else
constraint = Results{1,1};
Velocity = Results{1,2};
Acceleration = Results{1,3};
PartialPhi_ri= Results{1,4};
PartialPhi_pi= Results{1,5};

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The value of the expression of the constraint is %f \n', constraint);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The right-hand side of the velocity equation(mu) is %.5f \n', Velocity);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The right-hand side of the acceleration equation(gamma) is %.5f \n', Acceleration);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');

fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The expression for partial derivatives phi_R:\n');
fprintf(fileID, 'phi_Ri = %d.\n', PartialPhi_ri);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf(fileID,'The expression for partial derivatives phi_P:\n');
fprintf(fileID, 'phi_Ri = %d.\n', PartialPhi_pi);
fprintf(fileID,'- - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
end 
fclose(fileID);true
end

