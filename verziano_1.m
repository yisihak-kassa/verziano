close all
clear all
%% Channel properties
Sb = 0.001; % bed slope
ks = 81;    % Strickler's coefficient
mu = 0.41;  % discharge coefficient of lateral weir
B = 2;      % main channel width
b = 1.3;    % contracted channel width
Lv1 = 1.5;  % extent of converging section
Lv2 = 1.65; % extent of narrowest section
L = 10.95;	% extent of overflow-weir
Cc = 0.7;   % sluice gate coefficient of contraction
a = 0.25;   % sluice gate opening height

zbed = 106.2;       % channel bed height
zcrest = 107.75;    % lateral weir crest height
initial_crest_height = zcrest - zbed;  % height of the weir

hds = 0.5;  % water level downstream sluice gate
Q_out = 1;   % discharge to be delivered to the plant

%% Physical constants
g = 9.806;  % acceleration due to gravity
rho = 998;  % density

%% Section properties
[Aw, WPw, HRw, HDw, DCw] = rectangular_geometry(B);     % properties of the wider section
[An, WPn, HRn, HDn, Dcvar] = rectangular_geometry(b);   % properties of the narrower section

%% Accept a single or a series of input discharges from the user.
while true
    user_input = input('Enter a series of Q_in (space-separated): ', 's');    
    % Split the input string into an array of strings
    input_values = strsplit(user_input);    
    % Initialize a flag to check if all values are valid
    all_valid = true;    
    % Convert each input value to a number and check if it's valid
    Q_in = [];
    for i = 1:length(input_values)
        value = str2double(input_values{i});
        if isnan(value) || value < 1.5 || value > 2.0
            all_valid = false;
            break;
        else
            Q_in(i) = value;
        end
    end
    
    if all_valid
        break;
    else
        disp('Invalid input. Some Q_in are not valid or between 1.5 and 2.0.');
        % Ask if the user wants to continue
        choice = input('Do you want to continue? (yes/no): ', 's');
        if strcmpi(choice, 'no')
            return; % Exit the loop if the user does not want to continue
        end
    end
end

%% Verify the channel is mild sloped
for i = 1 : length(Q_in)
    [Q,V] = Q_chezy(Aw,HRw,ks,Sb);
    fn =@(h) Q(h) - Q_in(i);
    x0 = (Q_in(i)/ks/Sb^0.5/B)^(3/5);   % trial depth, ideal cross section
    h03 = newton(fn,x0);                % normal depth channel 3
    
    [Fr] = froude(Aw,HDw,Q_in(i));
    fc =@(h) Fr(h)-1;
    hc3 = newton(fc,0.5);               % critical depth, channel 3
    
    if h03 > hc3
        disp("The channel is mild")
    else
        disp("The channel is steep")
    end
end
%% Compute water depth immediately upstream the sluice gate assuming drowned outflow
[PH,VH,E] = energy(An, Q_out);
% Here, the piezometric depth is the total submerged depth "hds" and not the
% jet depth "a*Cc", while the kienetic head is computed using the jet
% depth. % For more details refer to "Open Channel Flow by 
% F.M Henderson, 1966, pp. 208)
Eds = hds + Q_out^2/(2*g*An(a*Cc)^2);
fh =@(h) E(h) - Eds;
hus = newton(fh,Eds);  % water depth upstream the sluice gate

%% Compute water depth just downstream of the weir
Econst1 = hus + Q_out^2/(2*g*b^2*hus^2);
fhdw =@(h) Econst1 - h - Q_out^2/(2*g*B^2*h^2);
hdw = newton(fhdw,Econst1);

%% Determine the crest height for a specific incoming discharge; and compute and plot the water surface profile above the weir
% table datastructure which stores the output
nrow = 4*length(Q_in);
ncol = 4;
varnames = {'Q','Number of lowered gates (cumecs)','Crest height (m)','Lowered by (m)'};
vartypes = {'double','int8','double','double'};
output = table('size',[nrow,ncol],'VariableTypes',vartypes,'VariableNames',varnames);

% boundary conditions
xP(1) = Lv2/2 + Lv1;
dc = 0.001;
dx = 0.1;
Lx = L/4*[1,2,3,4]; % array of lengths of the four weir configurations
for i=1:length(Q_in)    % for every incoming discharge
    overflow = Q_in(i) - Q_out;
    for m=1:length(Lx)  % for every configuration of the gate starting from the initial one
        clear xP zP hP Q2;
        % Define arrays with a fixed size "steps" and initialize them 
        % with temporary values to prevent array growth on every 
        % loop iteration and thus improve performance
        steps = ceil(Lx(m)/dx);
        xP = zeros(1,steps);
        zP = zeros(1,steps);
        hP = zeros(1,steps);
        Q2 = zeros(1,steps);
        
        xP(1) = Lv2/2 + Lv1;
        zP(1) = (Lv2/2 + Lv1)*Sb;
        hP(1) = hdw;
        Q2(1) = Q_out;
        c_var = initial_crest_height; % the varying depth of the crest
        found = false;
        while ~found
            j = 1;
            % determine discharge at the upstream corresponding to the current crest height
            % assuming constant specific energy along the weir
            while xP(j) < xP(1) + Lx(m)
                xP(j+1) = xP(j) + dx;
                zP(j+1) = zP(j) + Sb*dx;
                hw = hP(j)-c_var;  % water depth above the weir
                if hw > 0
                    dq = mu*sqrt(2*g)*(hw)^(3/2)*dx;
                else
                    break
                end
                Q2(j+1) = Q2(j) + dq;    % discharge at subsequent section
                fhP =@(h) E(hus) - h - Q2(j+1)^2/(2*g*Aw(h)^2);
                hP(j+1) = newton(fhP,hP(j));
                j = j + 1;
            end
            % if the discharges do not match adjust the discharge
            % by a small increment
            tentative_overflow = Q2(end)-Q_out;
            if abs(overflow-tentative_overflow) < 0.01
                found = true;
            elseif tentative_overflow < 1
                c_var = c_var - dc;
            else
                c_var = c_var + dc;
            end
        end
        currentrow = 4*(i-1) + m;
        output(currentrow, :) = {Q_in(i),m,c_var,(initial_crest_height-c_var)*100};
        
       % plot the profiles above the weir
        figure(i)
        titlei = sprintf("With %i gates lowered", m);
        subplot(2,2,m)
        hold on
        grid on;
        title(titlei);        
        yticks(0:0.1:1.5)
        plot(xP,zP + hP,'-r',LineWidth=1);
        plot(xP,mean(hP),'-r',LineWidth=1);
        plot(xP,zP,'k','LineWidth',0.5);
        crest = linspace(c_var,c_var,length(xP));
        plot(xP,crest,'k','LineWidth',0.50);   % plot the lateral weir
        xlabel('Station (m)','FontName','Arial','FontSize',10); % add the label for the x axis
        ylabel('Water depth (m)','FontName','Arial','FontSize',10); % add the label for the y axis
    end
end
disp(output)

