% ================================================================ %
% Program for two dimensional steady state diffusion
% ================================================================ %

clear           % clears variables
clc             % clears command window
cl = 0.1/8;     % sq.cell side length in meters
nl = cl;        % node to node distance. Same as side length
ny = 80;        % no of cells in y dir.
nx = 80;        % no of cells in x dir.
k = 1000;       % k value W/mK
a = k*(cl/nl);  % calculating 'a' value


% Boundary Conditions
% BC's syntax : [temp., flux, Temp. Surrounding] use - or + for choosing direction
% nan - not a number
h = -1000;              % Convective coefficient
N = [nan,0,nan];        % north wall
E = [nan,500,nan];      % east wall
W = [200,nan,nan];      % west wall
S = [nan,nan,50];       % south wall
D = [N;E;W;S];          % assigning given bc's in row wice in a matrix for calculation
d = zeros(4,1);         % 4*[source term , a value]

% Seperating boundary conditions & computing source terms
for i = 1:4                     % four boundaries - one to four
    if isnan(D(i,1)) == 0       % checking for temp boundary
        e = D(i,1)*2*a;         % calculating source term
        d(i,1) = e;
        d(i,2) = 2*a;           % additional ap value for temperature boundary
    elseif isnan(D(i,2)) == 0   % checking for flux boundary
        e = D(i,2)*cl;          % flux source term
        d(i,1) = e; 
        d(i,2) = 0;             % additional ap valur for flux boundary
    else                        % else convective boundary
        e = -D(i,3)*h*cl;       % calculating source term
        d(i,1) = e;
        d(i,2) = -h*cl;         % additional ap value due to convective boundary
    end
end

% Generating node grid
grid = zeros(ny+2,nx+2);        % Actual node grid
num = 1;
for i = 2:ny+1                  % Assigning node numbers in the grid with zeros as boundary
    for j = 2:nx+1
        grid(i,j) = num;
        num = num + 1;
    end
end

% Generating coefficient matrix
A = zeros(ny*nx,nx*ny);         % Coefficent grid
c = zeros(ny*nx,1);             % Source term collecting matrix --- RHS in equation
g = zeros(1,5);                 % Variable for collecting node numbers 
                                % 5 - current node number 1,2,3,4 - surrounding node numbers
for i = 2:ny+1
    for j = 2:nx+1
        g(1) = grid(i-1,j);     % collecting north node number
        g(2) = grid(i,j+1);     % collecting east node number
        g(3) = grid(i,j-1);     % collecting west node number
        g(4) = grid(i+1,j);     % collecting south node number
        g(5) = grid(i,j);       % collecting center node number (current position)
        ac = zeros(1,ny*nx);    % variable for collecting 'a' values --- equation LHS coefficients of particular node
        for p = 1:4             % iterating four boundaries 
            if g(p) == 0        % if the node number is zero ---> boundary
                ac(g(5)) = ac(g(5)) - d(p,2);   % adding ap term due to boundary (ap terms are negative)
                c(g(5)) = c(g(5)) - d(p,1);     % adding source term due to boundary (RHS of node equation)
            else                % else ---> surrounding node
                ac(g(p)) = ac(g(p)) + a;        % adding 'a' value of surrounding node
                ac(g(5)) = ac(g(5)) - a;        % adding 'a' value to ap term due to surrounding node
            end
        end
        A(grid(i,j),:) = ac;    % transfering nodel equation coefficients to the global 'A' matrix
    end
end

% Calculating solution
AI = inv(A);        % calculating inverse
t = AI * c;         % calculating temperature
num = 1;
T = zeros (ny,nx);  % placing temperature values in matrix
for i = 1:ny 
    for j = 1:nx
        T(i,j) = t(num);
        num = num + 1;
    end
end
pcolor(T)
shading interp;
colorbar;
colormap jet;
%contourf(T)
heatmap(T)

