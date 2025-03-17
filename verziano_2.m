% Calculation of the gate opening height needed to pass a discharge of 0.8m3/s

% Let’s focus on a single entering discharge equal to 2 mc/s and
% let’s start from the sill((weir?) configuration you designed in the
% previous exercise. Let’s assume that the waste water treatment
% plant should reduce to 0.8 the discharge that can be treated.
% Determine the opening of the sluice gate which allows reducing
% the discharge to this value and the resulting water surface
% profile
Qin = 2;
Qout = 0.8;
c = 1.335;
hds = 0.5;
Sb = 0.001; % bed slope
ks = 81;    % Strickler's coefficient
mu = 0.41;  % discharge coefficient of lateral weir
B = 2;      % main channel width
b = 1.3;    % contracted channel width
L = 10.95;
Cc = 0.7;
a = 0.25;

%% Physical constants
g = 9.806;

%% Section properties
[Aw, WPw, HRw, HDw, DCw] = rectangular_geometry(B);     % properties of the wider section
[An, WPn, HRn, HDn, Dcvar] = rectangular_geometry(b);   % properties of the narrower section

%% Computation of the sluice gate height
dx = 0.1;
dh = 0.0001;
hj = a*Cc;  % the jet depth
% Define arrays with a fixed size "steps" and initialize them 
% with temporary values to prevent array growth on every 
% loop iteration and thus improve performance
steps = ceil(L/dx);
xP = zeros(1,steps);
zP = zeros(1,steps);
hP = zeros(1,steps);
Q = zeros(1,steps);
Q(1) = Qout;
xP(1) = 0;
zP(1) = 0;
overflow = Qin - Qout;
% Vary the height until the right amount of overflow
% corresponding to Qout is found
found = false;  % loop control boolean variable
while ~found
    % take some initial value of jet depth downstream and
    % calculate the constant specific energy
    Econst = hds + Qout^2/(2*g*An(hj)^2);
    % calculate depth just downstream of the weir
    fhdw =@(h) Econst - h - Qout^2/(2*g*B^2*h^2);
    hdw = newton(fhdw,Econst);
    % integrate upstream and find the discharge upstream the weir
    k = 1;
    hP(1) = hdw;
    while xP(k) < L
        xP(k+1) = xP(k) + dx;
        zP(k+1) = zP(k) + Sb*dx;
        hw = hP(k)-c;  % water depth above the weir
        if hw > 0
            dq = mu*sqrt(2*g)*(hw)^(3/2)*dx;
        else
            break
        end
        Q(k+1) = Q(k) + dq;    % discharge at subsequent section
        fh =@(h) Econst - h - Q(k+1)^2/(2*g*Aw(h)^2);
        hP(k+1) = newton(fh,hP(k));
        k = k + 1;
    end

    tentative_overflow = Q(end)-Qout;
    disp(hj)
    if abs(overflow-tentative_overflow) < 0.01
        found = true;
    elseif tentative_overflow < 1.2
        hj = hj - dh;
    else
        hj = hj + dh;
    end
end

figure(1)
grid on
hold on

plot(xP,hP,'-r',LineWidth=2);
plot(xP,zP+.4,'k','LineWidth',0.5)
crest = linspace(c,c,length(xP));
plot(xP,crest,'k','LineWidth',0.50);   % plot the weir

a = hj/Cc;
fprintf("To get an output discharge of 0.8m3/s with an input discharge of 2m3/s, the gate opening should be %.2f m\n", a);
