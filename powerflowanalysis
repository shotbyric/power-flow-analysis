% Title:    Matlab Program to Solve for the state variables in a load flow problem
% Method:   Newton-Raphson's Method
% Author:   Ricardo Castillo
% Date:     April 2020

clear

fprintf("\nInitiating Code to Solve for State Variables in a Load Flow Problem Applying Newton-Raphson's Method for...\n");


%           Bus#    |V|         Angle   P           Q           Bus Type
linedata=[  1,      1.060,      0,      0,          0,          1;
            2,      1.040,      0,      1.50,       0,          3; 
            3,      1.030,      0,      1.00,       0,          3;
            4,      0,          0,      -1.00,      -0.70,      2; 
            5,      0,          0,      -0.90,      -0.30,      2; 
            6,      0,          0,      -1.60,      -1.10,      2];

        % Last column shows Bus Type: 1.Slack Bus 2.PQ Bus 3.PV Bus

%load the admittance matrix from the ybusmatrix.m file
ybusmatrix

%create a copy of the input matrix to preserve the original values
linedatacopy = linedata;
linedata;


%get the number of rows and columns of the matrix
[rownumlinedata,colnumlinedata]=size(linedata);

w=0; % total no of PV busses
u=0; % total no of PQ busses

for i=1:rownumlinedata % total no of PV busses
    if linedata(i,6)== 3
        w=w+1;
    end
end

for i=2:rownumlinedata % total no of PQ busses
    if linedata(i,6)== 2
        u=u+1;
    end
end
fprintf("\nOverview:\n-------------------------\n");
fprintf("Total Number of buses = %d\n", rownumlinedata);
fprintf("Total Number of PV buses = %d\n", w);
fprintf("Total Number of PQ buses = %d\n", u);

%------------------------------------------------------------------------------------------------------------
%
%
%------------------------------------------------------------------------------------------------------------
%make it so all |V| that are unknown get a default value of 1 for the first iteration
%also determine size of Jacobian based on the number of unknown |V| and del values

unknown_voltages = 0;
unknown_angles = -1; %ASUMING SLACK BUS REFERENCE ANGLE IS 0 SO DO NOT COUNT THIS ANGLE AS AN UNKNOWN

for i=1:rownumlinedata
    if linedata(i,2) == 0
        linedata(i,2) = 1;
        unknown_voltages = unknown_voltages + 1;
    end
    if linedata(i,3) == 0
        unknown_angles = unknown_angles + 1;
    end
end
size_of_Jacobian = unknown_voltages + unknown_angles;

%------------------------------------------------------------------------------------------------------------
%
%---------------------------Calculating Psch and Qsch for the first time-------------------------------------
%------------------------------------------------------------------------------------------------------------

%creating new arrays that will contain the values of Psch and Qsch
Psch = zeros(unknown_angles, 1);
Qsch = zeros(unknown_voltages, 1);

% adding the values of P and Q into the new arrays
for i=1:rownumlinedata
    Psch(i,1) = linedata(i,4);
    Qsch(i,1) = linedata(i,5);
end

%create a S scheduled matrix combined from the data extracted by Psch and Qsch (to be subtracted later by S calculated)
Ssch = zeros(size_of_Jacobian,1);
SPY = 1;
for i=1:rownumlinedata
    if Psch(i,1) ~= 0
        Ssch(SPY,1) = Psch(i,1);
        SPY = SPY + 1;
    end
end

for i=1:rownumlinedata
    if Qsch(i,1) ~= 0
        Ssch(SPY, 1) = Qsch(i,1);
        SPY = SPY + 1;
    end
end

%create a matrix that will hold the final value of the state variables
IXIC = zeros(size_of_Jacobian, 1);

%solution must converge with a tolerance of 0.0001
tol = 0.0001;

%flag used to identify if tolerance level has been achieved
flag = false;
hola = 1;
%------------------------------------------------------------------------------------------------------------
%
%----------------------------------Start of Iteration Loop:--------------------------------------------------
%------------------------------------------------------------------------------------------------------------

% while the tolerance condition is not met keep iterating
while flag == false;
    %creating the mismatch vector
    mismatch_func = zeros(unknown_angles + unknown_voltages,1);
    %initialize a temp vector that will store the dels temporally
    mismatch_func_angles = zeros(unknown_angles,1);
    %initialize a temp vector that will store the voltages temporally
    mismatch_func_voltages = zeros(unknown_voltages,1);
    %{

if hola ~= 1
for i = 1: size_of_Jacobian
if i <= unknown_angles
mismatch_func_angles(i,1) = IXIC (i,1);
end
if i > unknown_angles
mismatch_func_voltages(i-unknown_angles,1) = IXIC (i,1);
end
end
end
    %}
    %------------------------------------------------------------------------------------------------------------
    %
    %--------------------------------Calculating F(P) -> # of dels-----------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    for i=2:rownumlinedata
        RealP = 0;
        for j = 1:rownumlinedata
            RealP = RealP +linedata(i,2)*linedata(j,2)*abs(ybus(i,j))*cosd(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
            %mismatch_func_angles(i-1,1) = mismatch_func_angles(i-1,1) + linedata(i,2)*linedata(j,2)*abs(ybus(i,j))*cosd(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
        end
        mismatch_func_angles(i-1) = Ssch(i-1) - RealP;
    end
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Calculating F(Q) -> # of Vs-----------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    for i=4:rownumlinedata
        ImgQ = 0;
        for j=1:rownumlinedata
            ImgQ = ImgQ - abs(linedata(i,2))*abs(linedata(j,2))*abs(ybus(i,j))*sind(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
        end
        mismatch_func_voltages(i-3) = Ssch(i+2) - ImgQ;
    end
    %load the values of the Pcalc and Qcalc onto the mismatch matrix
    mismatch_func = [mismatch_func_angles;
        mismatch_func_voltages];
    mismatch_func;
    %if the tolerance condition has been meet, print the state variables,
    %mismatched function values and the iterations needed to complete
    %convergence
    if max(abs(mismatch_func)) < tol
        flag == true;
        fprintf("\n Solution:\n The system converges within %d iterations given an error tolerance of %d\n\n %d\n--------------------------\n", hola, tol);
        fprintf("Value of state variables after %d iterations:\n", hola);
        fprintf("[ Voltage Magnitudes ]\t\t[ Angles (degrees) ]\n", i,linedata(i,2),i,linedata(i,3));
        for i = 1 : rownumlinedata
            
            fprintf(" | V_%d | = %.7f \t\t del_%d = %.7f \n", i,linedata(i,2),i,linedata(i,3));
        end
        fprintf(" \nMismatch Function Values:  \n");
        for i = 1 : rownumlinedata
            fprintf(" %.7f \n", mismatch_func(i));
        end
        
        break
    end
    hola = hola +1;
    
    %create submatrices to derive the jacobian when put together
    J11=zeros(rownumlinedata-1,rownumlinedata-1);
    J12=zeros(rownumlinedata-1,u);
    J21=zeros(u,rownumlinedata-1);
    J22=zeros(u,u);
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Calculating J11-----------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    for i=2:rownumlinedata
        for j=2:rownumlinedata
            if i~=j
                J11(i-1,j-1)= -abs(linedata(i,2))*abs(linedata(j,2))*abs(ybus(i,j))*sind(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
            end
            if i==j
                for n=1:rownumlinedata
                    if n~=i
                        J11(i-1,i-1)=J11(i-1,i-1)+( -abs(linedata(i,2))*abs(linedata(n,2))*abs(ybus(n,j)))*sind(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
                    end
                end
            end
        end
    end
    J11;
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Calculating J12-----------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    for i=w:rownumlinedata
        for j=4:rownumlinedata
            %if linedata(j,6)== 2
            if i~=j
                J12(i-1,j-3)= abs(linedata(i,2))*abs(ybus(i,j))*cosd(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
            end
            yv12=0;
            if i==j
                for n=1:rownumlinedata
                    if n~=i
                        yv12=yv12 + abs(ybus(i,n))*linedata(n,2)*cosd(rad2deg(angle(ybus(i,n))) + linedata(n,3) - linedata(i,3));
                        %abs(ybus(i,j))
                    end
                end
                %J12(i-1,j-1) = abs(linedata(i,2))*(2*abs(linedata(i,2))*real(ybus(i,i))+yv12);
                J12(i-1,j-3) = yv12+ linedata(i,2)*2*abs(ybus(i,i))*cos(angle(ybus(i,i)));
            end
            %end
        end
    end
    J12;
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Calculating J21-----------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    for i=4:rownumlinedata
        %if linedata(i,6)== 2
        for j=w:rownumlinedata
            if i~=j
                J21(i-3,j-1)= -linedata(i,2)*linedata(j,2)*abs(ybus(i,j))*cosd(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
            end
            if i==j
                for n=1:rownumlinedata
                    if n~=i
                        J21(i-3,j-1)= J21(i-3,j-1)+ (( linedata(i,2)*linedata(n,2)*abs(ybus(i,n)))*cosd(rad2deg(angle(ybus(i,n))) + linedata(n,3) - linedata(i,3)));
                    end
                end
            end
        end
        %end
    end
    J21;
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Calculating J22-----------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    for i=4:rownumlinedata
        for j=4:rownumlinedata
            %if linedata(i,6)== 2 && linedata(j,6)== 2
            if i~=j
                J22(i-3,j-3)= -(linedata(i,2)*abs(ybus(i,j)))*sind(rad2deg(angle(ybus(i,j))) + linedata(j,3) - linedata(i,3));
            end
            yv22=0;
            if i==j
                for n=1:rownumlinedata
                    if n~=i
                        yv22= yv22 - linedata(n,2)*abs(ybus(i,n))*sind(rad2deg(angle(ybus(i,n))) + linedata(n,3) - linedata(i,3));
                        
                    end
                end
                J22(i-3,j-3)=  yv22 - 2*linedata(i,2)*abs(ybus(i,i))*sind(rad2deg(angle(ybus(i,i))));
            end
            
            %end
        end
    end
    J22;
    
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Computing the actual Jacobian---------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    
    Jacobian_Matrix = [J11, J12;
        J21, J22];
    Jacobian_Matrix;
    
    
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Computing the inverse Jacobian--------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    
    BATS = inv(Jacobian_Matrix);
    
    %------------------------------------------------------------------------------------------------------------
    %
    %----------------------------------Solving for the State variables-------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    ARCA = BATS*mismatch_func;
    
    
    %write the new values into the table
    
    for i= 1: size_of_Jacobian
        if i <= unknown_angles
            linedata(i+1,3) = linedata(i+1,3) + rad2deg(ARCA(i,1));
        end
        if i > unknown_angles
            
            linedata(i-2 ,2) = linedata(i-2,2) + ARCA(i,1);
        end
        
    end
    
    linedata;
    
end
