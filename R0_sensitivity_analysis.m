%==========================================================================
% File:        R0_sensitivity_analysis.m
%
% Title:       Sensitivity Analysis of Basic Reproduction Number (R0)
%              for Environmental Melioidosis Model with Recurrence
%
% Author:      Dr. Habtamu Ayalew Engida
%              Assistant Professor, Department of Mathematics
%              Debre Markos University, Ethiopia
%              Email:  hayalew21@gmail.com
%              ORCID:  0000-0003-2291-6207
%
% Publication: Engida, H.A. (2024). Modelling environmental-born
%              melioidosis dynamics with recurrence: An application
%              of optimal control.
%              Results in Control and Optimization.
%              DOI: https://doi.org/10.1016/j.rico.2024.100476
%
% Description: This script produces two figures:
%              Figure 1 - 3D surface plot of R0 as a function of
%                         transmission rate (alpha) and bacterial
%                         death rate (mu_b)
%              Figure 2 - 2D filled contour plot of R0 over the
%                         same parameter space
%
% Software:    MATLAB R2020 or later
% Date:        2024
%==========================================================================

clc;    % Clear command window
clear;  % Clear all workspace variables

%% ========================================================================
%  SECTION 1: Initial Conditions
%% ========================================================================

S0  = 450;   % Initial susceptible human population
E0  = 50;    % Initial exposed human population
I0  = 20;    % Initial infected human population
R10 = 10;    % Initial first-stage recovered population
R20 = 5;     % Initial second-stage recovered population
Bm0 = 200;   % Initial bacterial concentration in environment

% Total initial human population
N0 = S0 + E0 + I0 + R10 + R20;

%% ========================================================================
%  SECTION 2: Model Parameters
%% ========================================================================

% Demographic parameters
m1 = 1/(65 * 365);  % Natural human death rate (per day)
                    % Assumes average lifespan of 65 years
p  = m1 * N0;       % Human recruitment (birth) rate

% Disease parameters
a  = 0.0185;  % Environmental transmission rate (alpha)
t  = 1/9;     % Progression rate from exposed (E) to infected (I)
s  = 0.13;    % Bacterial shedding/environmental transmission factor
g  = 0.0037;  % Human recovery rate
m2 = 0.02;    % Bacterial death rate in environment (mu_b)
d  = 0.005;   % Disease-induced (melioidosis) death rate
c  = 5000;    % Environmental bacterial carrying capacity

% Recurrence parameters
e1 = 0.069;   % Rate of first recurrence/relapse
e2 = 0.035;   % Rate of second recurrence/relapse
tu = 0.26;    % Treatment/control rate

%% ========================================================================
%  SECTION 3: Figure 1 — 3D Surface Plot of R0
%% ========================================================================

figure(1)

% Define R0 as anonymous function of alpha (a) and mu_b (m2)
% R0 derived using next-generation matrix approach
f1 = @(a, m2) (p .* a .* t .* s .* (e1 + m1)) ./ ...
              (c .* m1 .* m2 .* (t + m1) .* ...
              ((e1 + m1) .* (d + m1) + m1 .* g));

% Generate surface over biologically realistic parameter ranges
% alpha (a)  : [0.007, 0.09]
% mu_b  (m2) : [0.001, 0.09]
fsurf(f1, [0.007 0.09 0.001 0.09])

% Figure formatting
grid on;
box on;
xlabel('$\alpha$', ...
       'FontSize', 20, ...
       'Interpreter', 'Latex');
ylabel('$\mu_b$', ...
       'FontSize', 20, ...
       'Interpreter', 'Latex');
zlabel('$\mathcal{R}_0$', ...
       'FontSize', 20, ...
       'Interpreter', 'Latex');
title('3D Surface: Sensitivity of $\mathcal{R}_0$', ...
      'FontSize', 14, ...
      'Interpreter', 'Latex');

% Save figure as EPS file
print -depsc Fig_R0_surface_3D

%% ========================================================================
%  SECTION 4: Figure 2 — 2D Filled Contour Plot of R0
%% ========================================================================

figure(2)

% Define parameter ranges for contour evaluation
a_range  = 0.0008 : 0.01 : 0.03;   % Range of alpha values
m2_range = 0.008  : 0.01 : 0.03;   % Range of mu_b values

% Create 2D grid of parameter combinations
[a_grid, m2_grid] = meshgrid(a_range, m2_range);

% Compute R0 at each grid point
R0 = (p .* a_grid .* t .* s .* (e1 + m1)) ./ ...
     (c .* m1 .* m2_grid .* (t + m1) .* ...
     ((e1 + m1) .* (d + m1) + m1 .* g));

% Plot filled contour map with 15 contour levels
contourf(a_grid, m2_grid, R0, 15, 'ShowText', 'on')

% Figure formatting
grid on;
box on;
colorbar;
xlabel('$\alpha$', ...
       'FontSize', 20, ...
       'Interpreter', 'Latex');
ylabel('$\mu_b$', ...
       'FontSize', 20, ...
       'Interpreter', 'Latex');
title('Contour Plot: Sensitivity of $\mathcal{R}_0$', ...
      'FontSize', 14, ...
      'Interpreter', 'Latex');

% Save figure as EPS file
print -depsc Fig_R0_contour_2D

%% ========================================================================
%  END OF SCRIPT
%  For questions or collaboration: hayalew21@gmail.com
%==========================================================================