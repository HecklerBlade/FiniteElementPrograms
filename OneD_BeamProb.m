%=======================================================================%
% Problem - 1D Beam                                                     %   
%=======================================================================%
%--------------------------------START----------------------------------%
tic               % For computing computational time 
clc               % Command Window clear
clear             % Workspace variable clear
%-----------------------------USER INPUTS-------------------------------%
% Any consistant unit system can be used.

%----------------------------Material def.------------------------------%
% E = [Modulus]
E = 210*10^9;

%---------------------------Beam properties-----------------------------%
% Total length
len = 5;
% Crossection type
% Use 'nan' to suppress
    % Square a = [side length]
    a = nan;
    % Circle r = [radius]
    r = nan;
    % rectangle R = [length, height]
    R = [nan,nan];
    % Manual MOI value for custom crossection
    I = 4*10^-4;

%---------------------------------Mesh----------------------------------%
% Number of elements
e = 19;
% Total number of nodes = e + 1

%-------------------------Boundary condition----------------------------%
% Bc = [Node number, Displacement val. [mm], Rotation val. [Deg]]
% Use '0' to suppress DOF
% For Fixed [Node number, 0, 0]
% For Pinned [Node number, 0, nan]

Bc = [1, 0, 0
    20, 0, 0];

%--------------------------------Load-----------------------------------%
% L = [Node number, Force, Moment]
% Use 'nan' to calculate reaction Force & Moment
% Use '0' to suppress Force/Moment
% For at Fixed [Node number, nan, nan]
% For Pinned [Node number, nan, 0]

L = [1, nan, nan
    10, -5000, 0
    20, nan, nan];

%-------------------------END OF USER INPUTS----------------------------%

%----------------------Assembly and calculation-------------------------%

%--------------------Moment of inertia calculation----------------------%
if (isnan(a) == 0)         % Reading Moment of inertia setting 
    sh = 1;                % Square
elseif (isnan(r) == 0)
    sh = 2;                % Circle
elseif (isnan(R(:,:)) == 0)
    sh = 3;                % Rectangle
else
    sh = 4;                % Manual
end

switch sh                  % Calculating Moment of inertia
    case 1
        I = (1/12) * a^4;
    case 2
        I = (1/4) * pi * r^4;
    case 3
        I = (1/12) * R(2) * R(1)^3;
    case 4                 
end

%----------------------Beam Stiffnes Matrix ref.------------------------%
el = len/e;
KK = [12 6*el -12 6*el
    6*el 4*el^2 -6*el 2*el^2
    -12 -6*el 12 -6*el
    6*el 2*el^2 -6*el 4*el^2];

%------------------------------LM array---------------------------------%
% DOF per node is 2 & 4 per element

LM = zeros(4,e);
n = 1;
for i=1:e
    LM(:,i)=LM(:,i)+(n:n+3)';
    n = n+2;
end

%--------------------Global Stiffness matrix assembly-------------------%
K = zeros((e+1)*2);                 %Pre-assigning zeros
for i = 1:((e+1)*2)
    for j = 1:((e+1)*2)
        for r = 1:e
            LMi = 0;
            LMj = 0;
            for t = 1:4
                if LM(t,r) == i     %Searching i index of K in LM array
                    LMi = t;
                end
                if LM(t,r) == j     %Searching j index of K in LM array
                    LMj = t;
                end
            end
            if (LMi ~= 0 && LMj ~= 0)
                K(i,j) = K(i,j) + KK(LMi,LMj); % Assigning values from 
                                               % KK stiffness ref.Mat.
            end
        end
    end
end

%---------------------Global force vector assembly----------------------%
f = zeros(e+1,2);        % For User load input reading
F = zeros(2*(e+1),1);    % For global force vector
Ls = size(L);

for i = 1:Ls(1)                 % Reading User input and storing
    for j = 1:e+1
        if L(i,1) == j
            f(j,:) = [L(i,2),L(i,3)];
        end
    end
end

fl = 1;                         % Flag variable
for i = 1:e+1                   % Assigning values in global F vector
    F(fl) = f(i,1);
    F(fl+1) = f(i,2);
    fl = fl + 2;
end 

%----------------------Global DOF vector assembly-----------------------%
d = nan .* zeros(e+1,2);        % For User Bc input reading
D = nan .* zeros(2*(e+1),1);    % For global d vector
Bcs = size(Bc);

for i = 1:Bcs(1)                % Reading User input and storing
    for j = 1:e+1
        if Bc(i,1) == j
            d(j,:) = [Bc(i,2),Bc(i,3)];
        end
    end
end

fl = 1;                         % Flag variable
for i = 1:e+1                   % Assigning values in global d vector
    D(fl) = d(i,1);
    D(fl+1) = d(i,2);
    fl = fl + 2;
end 

%------------------Finding unknown degrees of freedom-------------------%
%------------------------Row-Column elimination-------------------------%
Kd = K;
fl = 0;
for i = 1:(e+1)*2              % Making rows and column Zero
    if isnan(F(i)) == 1
        Kd(i,:) = Kd(i,:)*0;
        Kd(:,i) = Kd(:,i)*0;
        fl = fl + 1;
    end
end

i = fl;
j = 1;
while (i>0)                    % Removing Zero rows and columns
    if Kd(j,:) == 0
        Kd(j,:) = [];
        Kd(:,j) = [];
        i = i - 1;
    else
        j = j+1;
    end
end

%-----------------Force vector trimming for calculation-----------------%
Fd = zeros(((e+1)*2 - fl),2);
j = 1;
for i = 1:(e+1)*2
    if isnan(F(i)) == 0
        Fd(j,1) = i;
        Fd(j,2) = F(i);
        j = j + 1;
    end
end

%-----------------------Unknown DOF calculation-------------------------%
dd = Fd;
dd(:,2) = inv(((E*I)/(len/e)^3) .* Kd) * Fd(:,2);

%-----------------Finding reaction forces and moments-------------------%
%-------------------------Assembling D vector---------------------------%
for i = 1:((e+1)*2 - fl)
    D(dd(i,1)) = dd(i,2);
end

%-------------------------Calculating forces----------------------------%
RF = (((E*I)/(len/e)^3) .* K)*D;
toc

%--------------------------------POST-----------------------------------%
P1 = (1:2:(e+1)*2)';    % [1,3,5, ..... (e+1)*2] For Disp. and Force.
P2 = (2:2:(e+1)*2)';    % [2,4,6, ..... (e+1)*2] For Rot. and Moment.
figure('Name','Displacement');
plot(1:e+1,D(P1));
axis([1 e+1 min(D(P1))-10e-6 max(D(P1))+10e-6 ]) % Comment to suppress
figure('Name','Rotation');
plot(1:e+1,D(P2));
axis([1 e+1 min(D(P2))-10e-6 max(D(P2))+10e-6 ]) % Comment to suppress
figure('Name','Force');
plot(1:e+1,RF(P1));
axis([1 e+1 min(RF(P1))-100 max(RF(P1))+100])    % Comment to suppress
figure('Name','Moment');
plot(1:e+1,RF(P2));
axis([1 e+1 min(RF(P2))-100 max(RF(P2))+100 ])   % Comment to suppress

%---------------------------END OF PROGRAM------------------------------%
