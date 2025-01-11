% ENME 570 Assignment 2 - Jared Crebo 30085839

clc
clear
close all

%% Read Airfoil Data
% Symmetric Airfoils
naca0012 = readmatrix("NACA0012.csv");
naca0025 = readmatrix("NACA0025.csv");

% Cambered Airfoil
naca2412 = readmatrix("naca2412.csv");

%% Read Literature CL vs AoA Data
naca0012_CLvAoA_exp = readmatrix("naca0012_exp.csv");
naca0025_CLvAoA_exp = readmatrix("naca0025_exp.csv");
naca2412_CLvAoA_exp = readmatrix("NACA2412_exp.csv");

%% Parameters
rho = 1.225;
c = 1;
nu = 1.5e-5;
panels = 100; % Must be even number of panels
excluded_panel = panels/2 - 1; % Panel omitted to apply Kutta condition

% NACA0012 Experimental Re = 200,000
Re = 200000;
V_naca0012 = Re * nu / rho / c; %m/s

% NACA0025 Experimental Re = 5,000,000
Re = 3200000;
V_naca0025 =  Re * nu / rho / c; %m/s

% NACA2412 Experimental Re = 3,110,000
Re = 3.11e6;
V_naca2412 = Re * nu / rho / c;

%% Execute VPM
alpha = linspace(0,16,17);

% Panel Exclusion Sensitivity Analysis
num_panels = 36; % Number of panels for sensitivity analysis only
counter = 0;
for i = 13:18 % Change this to the indexes of panels to remove
    counter = counter + 1;
    for AoA = 0:16
        Cl_naca2412(AoA+1,counter) = vortexPanelMethod("NACA 2412", naca2412, AoA, num_panels, V_naca2412, i);
    end
end
figure('Position',[500,100, 1700, 1100]);
sgtitle("Sensitivity Analysis of Panel Omission with NACA 2412 and " + num2str(num_panels) + " Panels");
for j = 1:6
    subplot(2,3,j);
    hold on;
    grid on;
    plot(alpha, Cl_naca2412(:,j),'o-m');
    plot(naca2412_CLvAoA_exp(:,1),naca2412_CLvAoA_exp(:,2),'x-b');
    xlabel("Angle of Attack (deg)");
    ylabel("Coefficient of Lift");
    title("NACA 2412 Cl vs AoA (Re = 3110k, Panels = "+num2str(num_panels)+", Panel #"+num2str(i-6+j)+" excluded");
    xlim([0 16]);
    hold off;
end



%CL vs Angle of Attack Plots
for AoA = 0:16
    Cl_naca0012(AoA+1,1) = vortexPanelMethod("NACA 0012", naca0012, AoA, panels, V_naca0012, excluded_panel);
    Cl_naca0025(AoA+1,1) = vortexPanelMethod("NACA 0025", naca0025, AoA, panels, V_naca0025, excluded_panel);
    Cl_naca2412(AoA+1,1) = vortexPanelMethod("NACA 2412", naca2412, AoA, panels, V_naca2412, excluded_panel);
end

% Plot Cl v alpha for all three airfoils
figure('Position',[500, 400, 1700,500]);
subplot(1,3,1);
hold on;
grid on;
plot(alpha, Cl_naca0012(:,1),'o-m');
plot(naca0012_CLvAoA_exp(:,1),naca0012_CLvAoA_exp(:,2),'x-b');
plot(alpha, 2*pi*deg2rad(alpha), '-r');
xlabel("Angle of Attack (deg)");
ylabel("Coefficient of Lift");
title("NACA 0012 CL vs AoA (Re = 200k, Panels = " + num2str(panels) + ")");
legend('NACA0012 VPM','NACA0012 Lit','slope = 2\pi','Location','northwest');
xlim([0 16]);

subplot(1,3,2);
hold on;
grid on;
plot(alpha, Cl_naca0025(:,1),'o-m');
plot(naca0025_CLvAoA_exp(:,1),naca0025_CLvAoA_exp(:,2),'x-b');
plot(alpha, 2*pi*deg2rad(alpha), '-r');
xlabel("Angle of Attack (deg)");
ylabel("Coefficient of Lift");
title("NACA 0025 CL vs AoA (Re = 5000k, Panels = " + num2str(panels) + ")");
legend('NACA0025 VPM','NACA0025 Lit','slope = 2\pi','Location','northwest');
xlim([0 16]);

subplot(1,3,3);
hold on;
grid on;
plot(alpha, Cl_naca2412(:,1),'o-m');
plot(naca2412_CLvAoA_exp(:,1),naca2412_CLvAoA_exp(:,2),'x-b');
xlabel("Angle of Attack (deg)");
ylabel("Coefficient of Lift");
title("NACA 2412 CL vs AoA (Re = 3110k, Panels = " + num2str(panels) + ")");
legend('NACA2412 VPM','NACA2412 Lit','Location','northwest');
xlim([0 16]);

%% Vortex Panel Method (Airfoil Data, Angle of Attack (deg), Num of Panels, Freestream Velocity)
function Cl = vortexPanelMethod(name, airfoil, alpha, panels, V_inf, excluded_panel)
    % Change degrees to radians
    alpha = deg2rad(alpha);
    
    % Convert total number of panels to panel nodes per top/bottom surface
    total_panels = panels;
    panels_per_surface = panels / 2 + 1;
    
    % Separate airfoil data into xy coords of top and bottom surfaces
    x_top = airfoil(1:length(airfoil)/2+1, 1);
    y_top = airfoil(1:length(airfoil)/2+1, 2);
    x_bot = airfoil(length(airfoil)/2:end, 1);
    y_bot = airfoil(length(airfoil)/2:end, 2);

    % Interpolate b/w data points
    top_spline = spline(x_top, y_top);
    bot_spline = spline(x_bot, y_bot);

    % Discretize points along splines at specified resolution (# panels)
    fitted_top_spline_x = linspace(min(x_top), max(x_top), panels_per_surface);
    fitted_top_spline_y = ppval(top_spline, fitted_top_spline_x);
    fitted_bot_spline_x = linspace(min(x_bot), max(x_bot), panels_per_surface);
    fitted_bot_spline_y = ppval(bot_spline, fitted_bot_spline_x);

    % Discretize airfoil into panels    
    top_panel_centers = [];
    bot_panel_centers = [];
    top_panel_nodes = [];
    bot_panel_nodes = [];
    for i = 1:panels_per_surface
        % Record coords of edges of panels from discretized splines
        top_panel_nodes = [top_panel_nodes; [fitted_top_spline_x(i), fitted_top_spline_y(i)]];
        bot_panel_nodes = [bot_panel_nodes; [fitted_bot_spline_x(i), fitted_bot_spline_y(i)]];

        if i < panels_per_surface % There is one less center than nodes of panels
            % Calculate (x,y) coords of panel centers at midpoint b/w nodes
            top_panel_centers = [top_panel_centers; [(fitted_top_spline_x(i)+fitted_top_spline_x(i+1))/2, (fitted_top_spline_y(i)+fitted_top_spline_y(i+1))/2]];
            bot_panel_centers = [bot_panel_centers; [(fitted_bot_spline_x(i)+fitted_bot_spline_x(i+1))/2, (fitted_bot_spline_y(i)+fitted_bot_spline_y(i+1))/2]];
        end
    end

    % Correct the trailing edge to meet at (1,0)
    top_panel_nodes(end,2) = 0;
    bot_panel_nodes(end,2) = 0;
    top_panel_centers = [top_panel_centers; [(top_panel_nodes(end-1,1) + top_panel_nodes(end,1))/2, (top_panel_nodes(end-1,2) + top_panel_nodes(end,2))/2]];
    bot_panel_centers = [bot_panel_centers; [(bot_panel_nodes(end-1,1) + bot_panel_nodes(end,1))/2, (bot_panel_nodes(end-1,2) + bot_panel_nodes(end,2))/2]];

    % Plot airfoil shape
    if alpha == 0
        figure('Position',[700, 700, 1000, 400]);
        subplot(1,2,1);
        hold on;
        plot(x_top,y_top,'Color','red');
        plot(x_bot,y_bot,'Color','red');
        axis equal;
        hold off;
    % Plot to validate discretization method
        subplot(1,2,2);
        hold on;
        axis equal;
        plot(top_panel_nodes(:,1),top_panel_nodes(:,2),'b-o','DisplayName','Nodes');
        plot(bot_panel_nodes(:,1),bot_panel_nodes(:,2),'b-o');
        plot(top_panel_centers(:,1),top_panel_centers(:,2),'r*','DisplayName','Control Points');
        plot(bot_panel_centers(:,1),bot_panel_centers(:,2),'r*');
        legend('Nodes','','Control Points','');
        hold off;
        sgtitle("Discretization of " + name);
    end
    
    % Recombine calculated data of top and bottom surfaces
    % Order of panels from upper TE to LE to lower TE
    % Reorganize into (X,Y) and (x,y) coords for easier implementation to formulas
    panel_nodes = [top_panel_nodes(end:-1:2,:); bot_panel_nodes];
    X = panel_nodes(:,1)';
    Y = panel_nodes(:,2)';
    panel_centers = [top_panel_centers(end-1:-1:1,:); bot_panel_centers(1:end-1,:)];
    x = panel_centers(:,1)';
    y = panel_centers(:,2)';
    
    % Define panel angles (Phi) and panel lengths (S)
    for i = 1:length(x)
        % Panel angle relative to x-axis
        Phi(i) = mod(atan2(Y(i+1) - Y(i), X(i+1) - X(i)),2*pi);        
        % Panel lengths
        S(i) = sqrt((Y(i+1) - Y(i))^2 + (X(i+1) - X(i))^2);
    end

    for i = 1:length(x)
        for j = 1:length(x)
            % Coefficients from lecture example
            A(i,j) = -(x(i) - X(j))*cos(Phi(j)) - (y(i) - Y(j))*sin(Phi(j));
            B(i,j) = (x(i) - X(j))^2 + (y(i) - Y(j))^2;
            C(i,j) = -cos(Phi(i) - Phi(j));
            D(i,j) = (x(i) - X(j))*cos(Phi(i)) + (y(i) - Y(j))*sin(Phi(i));
            E(i,j) = sqrt(B(i,j) - A(i,j)^2);
            if ~isreal(E(i,j))
                E(i,j) = real(E(i,j));
            end
            t1 = atan2(S(j) + A(i,j) , E(i,j));
            t2 = atan2(A(i,j) , E(i,j));
            if i ~= j
                % Influence of panel j on panel i
                J(i,j) = (1/(2*pi)) * C(i,j)/2 * log((S(j)^2+2*A(i,j)*S(j)+B(i,j))/B(i,j)) + (D(i,j)-A(i,j)*C(i,j))/E(i,j) * (t1 - t2);
            else
                % Influence of panel on itself
                J(i,j) = 0;
            end
        end
        b(i,1) = - V_inf * sin(Phi(i) - alpha);
    end

    % Kutta Condition (upper TE + lower TE = 0)
    % First and last indexes are the TE surfaces
    % Exclude specified panel
    b(excluded_panel,1) = 0;
    J(excluded_panel,:) = 0;
    J(excluded_panel,1) = 1;
    J(excluded_panel,end) = 1;
    
    % Solve for vortex strengths
    gamma = J \ b;  

    for i = 1:length(x)
        influence = 0;
        for j = 1:2*panels_per_surface - 2
            influence = influence + (J(i,j) * gamma(j,1));
        end
        
        Vc(i,1) = V_inf * cos(alpha - Phi(i)) + influence;

        % Calculate pressure coefficient at control point
        Cp(i) = 1 - (Vc(i,1) / V_inf)^2;
    end

    % Circulation
    Gamma = gamma .* S';

    % Kutta-Jowkouski Theorem
    Cl_i = 2*Gamma / (V_inf * 1);
    Cl = sum(Cl_i);
    
    % To add a breakpoint at specific AoA and airfoils for debugging
    if alpha == deg2rad(3) && name == "NACA 2412"
        disp("Breakpoint");
    end
    
    % Plot Cp for specific AoA
    % if alpha == 0 || alpha == deg2rad(4) || alpha == deg2rad(8)
    %     figure;
    %     plot(x(1:length(x)/2),Cp(1:length(x)/2),'-r','DisplayName','Cp vs x/c');
    %     axis([0 1 min(Cp) max(Cp)]);
    %     xlim([0 1]);
    %     yline(0,'--');
    %     xlabel("x/c");
    %     ylabel("Coefficient of Pressure");
    %     title("Cp vs x/c of " + name + " (AoA = " + num2str(rad2deg(alpha)) + "deg, " + num2str((panels-1)*2) + " panels)");
    % end
end