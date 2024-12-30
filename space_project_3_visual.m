clear all
clc
UI=0;%THESE ARE REFRENCE PARAMETERS I USED DURING THE CODE
UO=0;
S=0;
while S==0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   HERE I ASK FOR THE USER INPUT 

fprintf('------------------------------------------\n')
fprintf('Enter 1 for euler angles as an input:\n')
fprintf('Enter 2 for Direct cosine matrix as an input:\n')
fprintf('Enter 3 for quatirions as an input:\n')
fprintf('Enter 4 for principle axis as an input:\n')
fprintf('Enter 5 to stop:\n')
fprintf('------------------------------------------\n')
 while UI==0 
user_input=input('select INPUT atitude representation:');
if user_input==5
    fprintf('------------------------------------------\n')
    disp('Exited')
    S=1;
    break
    
elseif user_input==1||user_input==2|| user_input==3|| user_input==4
UI=user_input;
elseif size(user_input)~=size([1 1])
    fprintf('invalid try again\n')
else
    fprintf('invalid try again\n')
 end
 end
if S==1
    break 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%HERE I ASK THE USER FOR THE REQUIRED OUTPUT

fprintf('------------------------------------------\n')
fprintf('Enter 1 for euler angles as an output:\n')
fprintf('Enter 2 for Direct cosine matrix as an output:\n')
fprintf('Enter 3 for quatirions as an output:\n')
fprintf('Enter 4 for principle axis as an output:\n')
fprintf('Enter 5 to stop:\n')
fprintf('------------------------------------------\n')

 while UO==0 
user_output=input('select OUTPUT atitude  representation:');
if user_output==5
    fprintf('------------------------------------------\n')
    disp('Exited')
    S=1;
    break
    
elseif user_output==user_input
    fprintf('------------------------------------------\n')
    fprintf('The chosen output is similar to the input\n')
    fprintf('------------------------------------------\n')
elseif user_output==1||user_output==2|| user_output==3||user_output==4
UO=user_output;

else
    fprintf('invalid try again\n')
 end
 end
 fprintf('------------------------------------------\n')
if S==1
    break 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%HERE WE ASK THE USER FOR THE INPUT BASED ON PREVIOUS INPUT AND OUTPUT

switch UI
    case 1 %HERE WE ASK FOR THE EULER ANGLES INPUT
S=1;
        X=input('enter the euler angles [1/pitch/x ,2/roll/y ,3/yaw/z] in degree:');
         
        while size(X,1)~=1  || size(X,2)~=3
            fprintf('------------------------------------------\n')
           X=input('Re-enter the euler angles [1/pitch/x ,2/roll/y ,3/yaw/z] in degree:');
    
           
        end
case 2 %HERE WE ASK FOR THE DCM INPUT
    S=1;
      seq=input('enter the DCM matrix [3-by-3] matrix:');
     
        while size(seq,1)~=3 || size(seq,2)~=3
            fprintf('------------------------------------------\n')
           seq=input('Re-enter the DCM matrix [3-by-3] matrix:') ;
        end
case 3 %HERE WE ASK FOR THE QUATERNIONS INPUT
    S=1;
                 quat=input('enter the quaternions elements[1-by-4]matrix [q1,q2,q3,q4]:');
                 
        while size(quat,1)~=1 ||size(quat,2)~=4
            fprintf('------------------------------------------\n')
           quat=input('Re-enter the quaternions elements [1-by-4] matrix [q1,q2,q3,q4]:'); 
        
        
        end

        case 4 %HERE WE ASK FOR THE PRINCIPLE AXIS AND PRINCIPLE ANGLE INPUT\
            S=1;
               prinp=input('enter the princple matrix [4-by-1] matrix:');
        while size(prinp,1)~=1 || size(prinp,2)~=4
            fprintf('------------------------------------------\n')
            S=1;
           prinp=input('Re-enter the principle matrix [4-by-1] matrix:') ;
        end
        end
        



%HERE WE CONVERT BASED ON USER INPUTS
if UI==1 && UO==2 %CONVERT BETWEEN EULER ANGLES TO DCM

rotation=askrotation();
seq=euler2dcm(rotation,X);
fprintf('The DCM MATRIX IS:\n')
disp(seq)
elseif  UI==1 && UO==3 %CONVERT BETWEEN EULER ANGLES TO QUATERNIONS
rotation=askrotation();
seq=euler2dcm(rotation,X);
[quat,q1,q2,q3,q4]=dcm2quat(seq);
fprintf('The QUATERNION MATRIX WHERE ELEMENT 1 TO 3 ARE VECTORS AND 4 IS A SCALER IS:\n')
disp(quat)
elseif UI==1 && UO==4 %CONVERT BETWEEN EULER ANGLES TO PRINCIPLE AXIS
 rotation=askrotation();
seq=euler2dcm(rotation,X);
[quat,q1,q2,q3,q4]=dcm2quat(seq);
[principal_line] = quaternion_to_principal(quat);
fprintf('The PRINCILE AXIS AND ANGLE MATRIX IS:\n')
disp(principal_line);
elseif UI==2 && UO==3 %CONVERT BETWEEN DCM  TO QUATERNIONS
[quat,q1,q2,q3,q4]=dcm2quat(seq);
fprintf('The QUATERNION MATRIX WHERE ELEMENT 1 TO 3 ARE VECTORS AND 4 IS A SCALER IS:\n')
disp(quat)
elseif UI==2 && UO==4 %CONVERT BETWEEN DCM  TO PRINCIPLE AXIS
    quat=dcm2quat(seq); 
[principal_line] = quaternion_to_principal(quat);
fprintf('The PRINCILE AXIS AND ANGLE MATRIX IS:\n')
disp(principal_line)
    elseif UI==3 && UO==4 %CONVERT BETWEEN QUATERNIONS TO PRINCIPLE AXIS
[principal_line] = quaternion_to_principal(quat);
fprintf('The PRINCILE AXIS AND ANGLE MATRIX IS:\n')
disp(principal_line)

elseif UO==1 && UI==2 %HERE WE CONVERT FROM DCM TO DCM TO EULER ANGLES 
rotation=askrotation();
X=dcm2euler(rotation,seq);
fprintf('The euler angles are [xyz]\n')
disp(X)

elseif UO==1 && UI==3 %HERE WE CONVERT FROM QUATERNIONS TO EULER ANGLES 
 
 seq=quat2dcm(quat);
rotation=askrotation();
X=dcm2euler(rotation,seq);

fprintf('The euler angles are [xyz]\n')
disp(X)
elseif UO==1 && UI==4 %HERE WE CONVERT FROM PRINCIPLE TO EULER ANGLES 
 [quat] = prinp_to_quat(prinp);
 seq=quat2dcm(quat);
rotation=askrotation();
X=dcm2euler(rotation,seq);
fprintf('The euler angles are [xyz]\n')
disp(X)



elseif UO==2 && UI==3 %HERE WE CONVERT FROM QUATERNIONS TO DCM
    seq=quat2dcm(quat);
rotation=askrotation()
fprintf('The DCM MATRIX IS:\n')
disp(seq);
elseif UO==2 && UI==4 %HERE WE CONVERT FROM PRINCIPLE AXIS TO DCM
[quat] = prinp_to_quat(prinp);
 seq=quat2dcm(quat);
 fprintf('The DCM MATRIX IS:\n')
disp(seq)
elseif UO==3 && UI==4 %HERE WE CONVERT FROM PRINCIPLE AXIS TO QUATERNIONS
[quat] = prinp_to_quat(prinp);

fprintf('The QUATERNION MATRIX WHERE ELEMENT 1 TO 3 ARE VECTORS AND 4 IS A SCALER IS:\n')
disp(quat)


end
end

function seq=euler2dcm(rotation,X) 
%HERE WE DO CONVERSION BETWEEN EULER TO DCM 
   psi=X(3);
           theta=X(2);
           phi=X(1);
R1=[[1 0 0];[0 cosd(phi) sind(phi)];[0 -sind(phi) cosd(phi)]];
R2=[[cosd(theta) 0 -sind(theta)];[0 1 0];[sind(theta) 0 cosd(theta)]];
R3=[cosd(psi) sind(psi) 0;-sind(psi) cosd(psi) 0;0 0 1];
switch rotation
    case 1
        seq=R3*R2*R1;
    case 2
      seq=R2*R3*R1 ; 
    case 3
        seq=R3*R1*R2;
    case 4
        seq=R1*R3*R2;
    case 5
        seq=R1*R2*R3;
    case 6
        seq=R2*R1*R3;
end
seq=seq';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPECAIL CASE FOR WHEN SWITCHING BETWEEN EULER AND OTHER TYPES AND VISE
%VERSA
%HERE WE ASK THE USER WHICH ROTATION SEQUENCE IS NEEDED
function rotation=askrotation()
fprintf('------------------------------------------\n')
    fprintf('1 For R123\n')
    fprintf('2 For R132\n')
    fprintf('3 For R213\n')
    fprintf('4 For R231\n')
    fprintf('5 For R321\n')
    fprintf('6 For R312\n')   
fprintf('note the 3=yaw=z 2=roll=y 1=pitch=x\n')
fprintf('------------------------------------------\n')
rotation=input('Select the Rotation Sequence of the DCM:');
while rotation>6 || rotation<1

    rotation=input('Re-enter the Rotation Sequence of the DCM:');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HERE WE DO THE CONVERSION BETWEEN DCM AND QUATERNIONS
function [quat,q1,q2,q3,q4]=dcm2quat(seq);
    a11=seq(1,1);
    a22=seq(2,2);
    a33=seq(3,3);
    a12=seq(1,2);
    a21=seq(2,1);
    a31=seq(3,1);
    a13=seq(1,3);
    a23=seq(2,3);
    a32=seq(3,2);
q4=0.5*((1+a11+a22+a33)^(1/2));
q3=(a12-a21)*(1/(q4*4));
q2=(a31-a13)*(1/(4*q4));
q1=(a23-a32)*(1/(4*q4));

quat=[-q1 -q2 -q3 q4];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HERE WE CHANGE FROM QUATERMIONS TO PRINCIPLE AXIS
function [principal_line] = quaternion_to_principal(quat)
    quat = quat ;
    angle = 2 * acosd(quat(4));
    axis = quat(1:3) / sind(angle / 2);
    principal_line = [axis, angle];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HERE WE CHANGE FROM DCM TO EULER ANGLES
function X=dcm2euler(rotation,seq)
syms phi1 theta1 psi1 %USING SYMBOLIC MATH TOOLBOX IS HELPFUL TO MULTIPLY THE MATRIX DEPEINDING ON THE ROTATION MATRIX CHOSEN
R1=[[1 0 0];[0 cosd(phi1) sind(phi1)];[0 -sind(phi1) cosd(phi1)]];
R2=[[cosd(theta1) 0 -sind(theta1)];[0 1 0];[sind(theta1) 0 cosd(theta1)]];
R3=[cosd(psi1) sind(psi1) 0;-sind(psi1) cosd(psi1) 0;0 0 1];
switch rotation
    case 1
        seq1=R3*R2*R1;
        theta1=asind(seq(1,3));
        phi1=-atand(seq(2,3)/seq(3,3));
        psi1=-atand(seq(1,2)/seq(1,1));
       
    case 2
      seq1=R2*R3*R1;
 psi1=-asind(seq(1,2));
phi1=atand(seq(3,2)/seq(2,2));
theta1=atand((seq(1,3)/seq(1,1)));


    case 3
        seq1=R3*R1*R2;
        seq1=seq1';
  
        phi1=-asind(seq(2,3));
        psi1=atand(seq(2,2)/seq(2,1));
        theta1=atand(seq(1,3)/seq(3,3));
    case 4
        seq1=R1*R3*R2;
        seq=seq'
        seq1=seq1';

         psi1=asind(seq(1,2));
        phi1=-atand(seq(3,2)/seq(2,2));
        theta1=-atand((seq(1,3)/seq(1,1)));
    case 5
        seq1=R1*R2*R3;
            seq=seq'
        seq1=seq1';
        theta1=-asind(seq(1,3));
        phi1=atand(seq(2,3)/seq(3,3));
        psi1=atand(seq(1,2)/seq(1,1));
    case 6
        seq1=R2*R1*R3;

      
        phi1=asind(seq(3,2));
        psi1=-atand(seq(1,2)/seq(2,2));
        theta1=-atand(seq(3,1)/seq(3,3));
end

X=[phi1 theta1 psi1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HERE WE CHANGE FROM QUATERNIONS TO DCM
function seq=quat2dcm(quat)
q1=quat(1);
q2=quat(2);
q3=quat(3);
q4=quat(4);
% I CHANGED THE NEGATIVE SGINS HERE TO TEST THE OUTPUT a13 a21 a32
seq=[q1^2-q2^2-q3^2+q4^2,2*(q1*q2-q4*q3),2*(q1*q3+(q4*q2));2*(q2*q1+(+q4*q3)),-q1^2+q2^2-q3^2+q4^2,2*(q2*q3-q4*q1);2*(q3*q1-q4*q2),2*(q3*q2+q4*q1),-q1^2-q2^2+q3^2+q4^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function[ quat] = prinp_to_quat(prinp);
%HERE WE CHANGE FROM QUATERNIONS TO PRINCIPLE AXIS

p_axis = prinp(1:3); 
p_angle = prinp(4); 

p_axis = p_axis / norm(p_axis);

q4 = cosd(p_angle/2);
q1 = p_axis(1) * sind(p_angle/2);
q2 = p_axis(2) * sind(p_angle/2);
q3 = p_axis(3) * sind(p_angle/2);


quat = [q1 q2 q3 q4];
end