global e invDirections Lx Ly b dx dt isObstacle isObstacle_full omega cs cs2inv cs4inv X Y XX YY;
Lx = 1;
Ly = 1;
% Nt = 10000;
endT = 1;
Nx = 100;
Ny = 100;
dt = 0.0001;
dx = Lx/Nx; %Attalums starp shunam
dy = Ly/Ny;
cs = sqrt(dx^2/dt^2/3);
cs2inv = 1/cs^2;
cs4inv = 1/cs^4;

X = linspace(-dx/2, Lx+dx/2, Nx+2);
Y = linspace(-dy/2, Ly+dy/2, Ny+2);
[XX, YY] = meshgrid(X, Y);

inner = meshgrid(2:Nx-1, 2:Ny-1);

isObstacle = padarray(zeros(Nx, Ny), [1, 1], 1);
D = 1; % SV koeficients
tau = thermalRelaxationTime(D);

T = sin(pi*XX'/Lx).*sin(pi*YY'/Ly);
e = [0, 0; 1, 0; 0, 1; -1, 0; 0, -1; 1, 1; -1, 1; -1, -1; 1, -1];
b = length(e); %2DQ9
isObstacle_full = repmat(isObstacle, [1, 1, b]);

invDirectionIndex = @(ie) find(ismember(e, -e(ie, :), 'rows'), 1, 'first');
%precalculate inverse directions for each ie
invDirections = zeros(1, b);
for ie=1:b
    invDirections(ie) = invDirectionIndex(ie); 
end

omega = [16, 4, 4, 4, 4, 1, 1, 1, 1]/36;

colormap jet;

% Sakuma atrodas lidzsvara stavokli
f = equilibrium(T);
% atrodam stabilu sakuma stavokli
% while 1
%     
%     if %konverge
%         break;
% end
initOutput(T, T);

TT = {};
% katra laika solii
ErrAbs = [];
ErrRel = [];
ErrNorm = [];
ErrNormRel = [];
tt = 0:dt:endT;
for t=tt
    tFrame = tic;
        tic
            % Calculate RHO and VX, VY
            T = calculate_macro(f);
        TT.macro = toc;

        tic
            sol = analyticSolution(t);    
            visualize(T, sol);
            
            ErrAbs(end+1) = max(abs((T(2:end-1, 2:end-1) - sol(2:end-1,2:end-1))), [], 'all');
            ErrRel(end+1) = max(abs((T(2:end-1, 2:end-1) - sol(2:end-1, 2:end-1))./sol(2:end-1, 2:end-1)), [], 'all');
            ErrNorm(end+1) = norm(T(2:end-1, 2:end-1)- sol(2:end-1, 2:end-1), 2);
            ErrNormRel(end+1) = norm(T(2:end-1, 2:end-1)- sol(2:end-1, 2:end-1), 2)./norm(sol(2:end-1, 2:end-1), 2);
        TT.visualize = toc;

        tic
            % calculate equilibrium
            fEq = equilibrium(T);
        TT.equilibirum = toc;

        tic
            % collision step
            ftemp = collision(f, fEq, tau);
        TT.collision = toc;

        tic
            %streaming (propagation) step
            f = streaming(ftemp);
        TT.streaming = toc;
    TT.total = toc(tFrame)
end

function f_res = stream(f, e)
    global b;
    % streamo katru no virzieniem
    f_res = zeros(size(f));
    for k = 1:b
        f_res(:, :, k) = circshift(f(:,:,k),e(k, :));
    end
end

function f = streaming(ftemp)
    global invDirections isObstacle e b;
    
    ftemp(isnan(ftemp)) = 0;
    f = stream(ftemp, e);
    
    % Atspogulojam atpakal no malam un skersliem
    for ie=1:b
        % vertibas, kuras iegaja ieksa "sienaas"
        toAdd = isObstacle.*f(:, :, ie);
        
        % Dirihle, anti-bounceback
        % nobidam atpakal
        toAdd = -circshift(toAdd, -e(ie, :));
        f(:, :, invDirections(ie)) = f(:, :, invDirections(ie)) + toAdd;
    end
end

function T = calculate_macro(f)
    T = sum(f, 3);
end

function fEq = equilibrium(T)
    global b omega;
    fEq = zeros([size(T), b]);
    for ie = 1:b
       fEq(:, :, ie) = omega(ie)*T;
    end
end

function ftemp = collision(f, fEq, tau)
    global b isObstacle_full dt;
    
    %collision (relaxation) step
    ftemp = zeros(size(f));
    for ie = 1:b
        ftemp(2:end-1, 2:end-1, ie) = f(2:end-1, 2:end-1, ie)*(1-dt/tau) + fEq(2:end-1, 2:end-1, ie)*dt/tau;
    end
    ftemp(isObstacle_full == 1) = NaN;
end

function initOutput(T1, T2)
    global isObstacle I1 I2 X Y;
    figure(1);
    I1 = imagesc(X, Y, T1', [0, 1]);
    set(I1, 'AlphaData', ~isObstacle');
    colormap(jet);
    colorbar;
    title('Temperature');
    figure(2);
    I2 = imagesc(X, Y, T2', [0, 1]);
    set(I2, 'AlphaData', ~isObstacle');
    colormap(jet);
    colorbar;
    title('Analytic');
end

function visualize(T1, T2)
    global I1 I2;
    set(I1, 'CData', T1');
    set(I2, 'CData', T2');
    drawnow limitrate;  
end

function tauG = thermalRelaxationTime(D)
    global cs dt;
    syms tau;
    tauG = double(solve(D == cs^2*(tau - dt/2)));
end

function T = analyticSolution(t)
    global Lx Ly XX YY;
    T = sin(pi*XX'/Lx).*sin(pi*YY'/Ly).*exp(-pi^2*(1/Lx^2 + 1/Ly^2)*t);
end

function [dx, dy] = initialGradient()
    global Lx Ly XX YY;
    dx = pi/Lx*cos(pi*XX'/Lx).*sin(pi*YY'/Ly);
    dy = pi/Ly*sin(pi*XX'/Lx).*cos(pi*YY'/Ly) ;
end