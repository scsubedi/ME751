%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position, Velocity and Acceleration analysis of different joints and
% Intermediate level Gcons on two bodies. One of the bodies is grounded. 
% Input: input all the attributes in the 'input.txt' file
% Output: Position, velocity and acceleration of ungrounded body based on the inputs. 
% For Revolute joint, there's an option to save a video and plots for position, velocity and acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the input file from input.txt and saves them as variables in
% workspace

readInput;

%%
if strcmp(joint,'Revolute')
    revoluteJoint;
elseif strcmp(joint,'Spherical')
    sphericalJoint;
elseif strcmp(joint,'Universal')
    universalJoint;
elseif strcmp(joint,'Cylindrical')
    cylindricalJoint;
elseif strcmp(joint,'Perpendicular1')
    perpendicular1;
elseif strcmp(joint,'Perpendicular2')
    perpendicular2;
elseif strcmp(joint, 'TSDA')
    TSDA;
elseif strcmp(joint,'RSDA')
    RSDA;
else
    disp('Not a well-defined joint')
end


if strcmp(type,'None')
    disp(' ')
elseif strcmp(type, 'TSDA')
    TSDA;
elseif strcmp(type,'RSDA')
    RSDA;
else
    disp('Unknown Forces or Torque')
end

%%%% Results and Post-Processing
if strcmp(joint,'Revolute')
    if strcmp((val{11,2}),'Video')
        disp('Video of the motion requested.')
        makeVideo(output)
        disp('Video has been saved.')
    elseif strcmp((val{11,2}),'Plots')
        disp('Plot has been requested')
        saveImage(L, output)
        disp('Plot has been saved')
    elseif strcmp((val{11,2}),'Both')
        disp('Video and Plot have been requested.')
        makeVideo(output)
        saveImage(L, output)
        disp('Video and plot have been saved.')
    else
        disp('Unrecognized Option, please check the input file')
    end
    
end

