%% Assignment 4: Additive & Shunting Networks

clc; clear all; close all;

% Constants
A = .1;
B = 1;
x0 = 0;
I = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
neurons = 10;
dt = 0.01; 

[x, t] = additive(A, B, x0, I, neurons, dt); % Additive network
[x1, t1] = shunting(A, B, x0, I, neurons, dt); % Shunting network

x_norm = x;
x1_norm = x;

sums = sum(x_norm, 1);

% Normalizing additive network cell activities
for i = 1:length(t)
    x_norm(:, i) = x_norm(:, i) / sums(i);
end

%% Plotting additive newtork 

figure(1)

for i = 1:neurons
    
    plot(t, x(i, :))
    hold on
    
end

xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Additive STM Network Responses', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')

figure(2)

for i = 1:neurons
    
    plot(t, x_norm(i, :))
    hold on
    
end

xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Additive STM Network Responses (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')

%% Plotting shunting network

figure(3)

for i = 1:neurons
    
    plot(t1, x1(i, :))
    hold on
    
end

xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Shunting STM Network Responses', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')

figure(4)

sums = sum(x1_norm, 1);

% Normalizing shunting network cell activities
for i = 1:length(t)
    x1_norm(:, i) = x1_norm(:, i) / sums(i);
end

for i = 1:neurons
    
    plot(t, x1_norm(i, :))
    hold on
    
end

xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Shunting STM Network Responses (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')

%% Distance-Dependent Network

% First, the Constants

k_i = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

% k_i = [3, 2, 1, 0, 1, 2, 3, 4, 5; ...
%     2, 1, 0, 1, 2, 3, 4, 5, 6; ...
%     1, 0, 1, 2, 3, 4, 5, 6, 7; ...
%     0, 1, 2, 3, 4, 5, 6, 7, 8; ...
%     1, 2, 3, 4, 5, 6, 7, 8, 9; ...
%     2, 3, 4, 5, 6, 7, 8, 9, 10; ...
%     3, 4, 5, 6, 7, 8, 9, 10, 11; ...
%     4, 5, 6, 7, 8, 9, 10, 11, 12; ...
%     5, 6, 7, 8, 9, 10, 11, 12, 13; ...
%     6, 7, 8, 9, 10, 11, 12, 13, 14];

C_plt = C_plot(k_i);
E_plt = E_plot(k_i);

figure(5)

subplot(1, 2, 1)

plot(k_i, C_plt)
xlabel('Distance', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Value of Excitation Coefficient (C)', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Distance-Dependent Excitation in an STM Network', 'FontName', 'Times New Roman', 'FontSize', 14);

subplot(1, 2, 2)

plot(k_i, E_plt)
xlabel('Distance', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Value of Inhibition Coefficient (E)', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Distance-Dependent Inhibition in an STM Network', 'FontName', 'Times New Roman', 'FontSize', 14);
%% Now the Network

A_inp = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.8, 0.8, 0.8, 0.8, 0.8 ];
B_inp = [ 0.1, 0.1, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1, 0.1, 0.1 ];
C_inp = [ 0.1, 0.1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.1, 0.1 ];
D_inp = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ];

[V_dist1, t_dist1] = distanceDep(A, B, x0, A_inp, neurons, dt);
[V_dist2, t_dist2] = distanceDep(A, B, x0, B_inp, neurons, dt);
[V_dist3, t_dist3] = distanceDep(A, B, x0, C_inp, neurons, dt);
[V_dist4, t_dist4] = distanceDep(A, B, x0, D_inp, neurons, dt);

[V_shunt1, t_shunt1] = shunting(A, B, x0, A_inp, neurons, dt);
[V_shunt2, t_shunt2] = shunting(A, B, x0, B_inp, neurons, dt);
[V_shunt3, t_shunt3] = shunting(A, B, x0, C_inp, neurons, dt);
[V_shunt4, t_shunt4] = shunting(A, B, x0, D_inp, neurons, dt);

figure(6)

subplot(1, 2, 1)

plot(t_dist1, V_dist1)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Distance-Dependent STM Network Responses for input A (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.7])

subplot(1, 2, 2)

plot(t_shunt1, V_shunt1)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Shunting STM Network Responses for input A (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.2])

figure(7)

subplot(1, 2, 1)

plot(t_dist2, V_dist2)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Distnace-Dependent STM Network Responses for input B (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.7])

subplot(1, 2, 2)
plot(t_shunt2, V_shunt2)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Shunting STM Network Responses for input B (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.35])

figure(8)

subplot(1, 2, 1)

plot(t_dist3, V_dist3)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Distnace-Dependent STM Network Responses for input C (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.7])

subplot(1, 2, 2)

plot(t_shunt3, V_shunt3)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Shunting STM Network Responses for input C (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.2])

figure(9)

subplot(1, 2, 1)

plot(t_dist4, V_dist4)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Distnace-Dependent STM Network Responses for input D (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'Southwest')
axis([800, 1000, 0, 0.7])

subplot(1, 2, 2)

plot(t_shunt4, V_shunt4)
xlabel('Time (t)', 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel('Equilibrium STM Activity of neuron', 'FontName', 'Times New Roman', 'FontSize', 14)
title('Shunting STM Network Responses for input D (EQ)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', 'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10', 'Times New Roman', 'FontSize', 14, 'Location', 'SouthWest')
axis([800, 1000, 0, 0.2])

%% modAdditive

B_B = 100;
A_A = 0.1;

[x_madd, t_madd] = modAdditive(A_A, B_B, x0, I, neurons, dt);

figure(10)

subplot(1, 2, 1)

plot(t_madd, x_madd)

[x_add, t_Add] = additive(A_A, B_B, x0, I, neurons, dt); 

subplot(1, 2, 2)

plot(t_Add, x_add)

%% functions

% Additive neuron model
function [x, t] = additive(A, B, x0, I, neurons, dt)

    t = 0:6000;
    x = zeros(neurons, length(t));
    x(:, 1) = x0;
    
    I_inh = zeros(1, length(I)) + sum(I);
    I_inh = I_inh - I;
    
    for i = 1:neurons
        for j = 1:length(t) - 1
            
            x(i, j+1) = x(i, j) + dt*(-A*x(i, j) + B*I(i) - I_inh(i));
         
        end
    end
    
end

% Shunting neuron model
function [x, t] = shunting(A, B, x0, I, neurons, dt)

    t = 1:1000;
    x = zeros(neurons, length(t));
    x(:, 1) = x0;
    
    I_inh = zeros(1, length(I)) + sum(I);
    I_inh = I_inh - I;
    
    for i = 1:neurons
        for j = 1:length(t) - 1
            
            x(i, j+1) = x(i, j) + dt*(-A*x(i, j) + (B - x(i, j))*I(i) - x(i, j)*I_inh(i));
         
        end
    end

end

% Modified additive neuron model (activity-dependent input)
function [x, t] = modAdditive(A, B, x0, I, neurons, dt)

    t = 1:1000;
    x = zeros(neurons, length(t));
    x(:, 1) = x0;
    
    I_inh = zeros(1, length(I)) + sum(I);
    I_inh = I_inh - I;
    
    for i = 1:neurons
        for j = 1:length(t) - 1
            
            x(i, j+1) = x(i, j) + dt*(-A*x(i, j) + (B - x(i, j))*I(i) - I_inh(i));
         
        end
    end

end

% Inhibition Constant
function [C] = C_plot(k_i)
    C = ones(size(k_i));
    
    for i = 1:length(k_i)
        
            C(i) = exp(-(k_i(i).^2) ./ 4);
    end
end

% Excitation Constant
function [E] = E_plot(k_i)
    E = ones(size(k_i));
    
    for i = 1: length(k_i)

            E(i) = .5*exp(-(k_i(i).^2) ./16);

    end
end

function [C_sum, i] = C(i, I)

    C_sum = 0;

    for k = (i-4:i+4)
        if i - 4 <= 0
            ind = 1;
            
        elseif i + 4 > 10
            ind = 10;
        else
            ind = k;
        end
        
        C = exp(-(ind - i).^2 ./ 4)*I(ind);
        C_sum = C_sum + C;
        
    end

end

function [E_sum, i] = E(i, I)

    E_sum = 0;

    for k = (i-4:i+4)
        
        if i - 4 <= 0
            ind = 1;
        elseif i + 4 > 10
            ind = 10;
        else
            ind = k;
        end
        
        E = 0.5.*exp(-(ind - i).^2 ./16)*I(ind);
        E_sum = E_sum + E;
        
    end
    
end

% Lateral, 'distance-dependent' inhibition of neurons 
function [x, t] = distanceDep(A, B, x0, I, neurons, dt)

    t = 1:1000;
    x = zeros(neurons, length(t));
    x(:, 1) = x0;
    
    for i = 1:neurons
        
        for j = 1:length(t) - 1
            
            x(i, j+1) = x(i, j) + dt*(-A*x(i, j) + (B - x(i, j))*C(i,I) - x(i, j)*E(i, I));
         
        end
    end

end
    