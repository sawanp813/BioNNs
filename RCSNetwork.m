%% Recurrent, Competetive Shunting Network - Sawan Patel

clc; clear; close all;

% Constants
I = [ 0.2, 0.6, 0.9, 0.6, 0.2, 0.1, 0.4, 0.8, 0.4, 0.1 ];
A = 1;
B = 3;
T = 10;
x0 = 0;
neurons = 10;
dt = 0.001;

[x1, t1] = network(A,B,T,dt,I,x0,neurons,2);

for i = 1:length(t1)
    x1(:, i) = x1(:, i) / sum(x1(:, i));
end

figure(1)

for i = 1:neurons

    hold on
    plot(t1, x1(i, :))

end

xlabel('Time (t)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('STM Activity (x_i)', 'FontSize', 12, 'FontName', 'Times New Roman')
title('Recurrent Competetive Network Activity for Input #1 (SF #4)', 'FontSize', 12, 'FontName', 'Times New Roman')
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', ...
    'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10')

figure(2)

% for i = 1:neurons
%     
%     hold on
%     plot(t1(1901:length(t1)-1), x1(i, 1901:length(t1)-1))
%     
% end

plot(1:10, x1(:, 2001))

xlabel('Neuron Index', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Final STM Activity (x_i)', 'FontSize', 12, 'FontName', 'Times New Roman')
title('Recurrent Competetive Network Activity (at EQ) for Input #1 (SF #4)', 'FontSize', 12, 'FontName', 'Times New Roman')
axis([1 10 0 0.2])

I = [ 0.7, 0.6, 0.8, 0.9, 0.5, 0.3, 0.5, 0.7, 0.8, 0.4 ];

[x2, t2] = network(A,B,T,dt,I,x0,neurons,2);

for i = 1:length(t1)
    x2(:, i) = x2(:, i) / sum(x2(:, i));
end

figure(3)

for i = 1:neurons

    hold on
    plot(t2, x2(i, :))

end

xlabel('Time (t)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('STM Activity (x_i)', 'FontSize', 12, 'FontName', 'Times New Roman')
title('Recurrent Competetive Network Activity for Input #2 (SF #4)', 'FontSize', 12, 'FontName', 'Times New Roman')
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Neuron 4', 'Neuron 5', 'Neuron 6', ...
    'Neuron 7', 'Neuron 8', 'Neuron 9', 'Neuron 10')

figure(4)

% for i = 1:neurons
%     
%     hold on
%     plot(t2(1901:length(t2)-1), x2(i, 1901:length(t2)-1))
%     
% end

plot(1:10, x2(:, 2001))

xlabel('Time (t)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('STM Activity (x_i)', 'FontSize', 12, 'FontName', 'Times New Roman')
title('Recurrent Competetive Network Activity (at EQ) for Input #2 (SF #4)', 'FontSize', 12, 'FontName', 'Times New Roman')
axis([1 10 0 0.2])






%% Functions

% Network simulation
function [x, t] = network(A, B, T, dt, I, x0, neurons, Signal)

   t = 0:dt:T;
   x = zeros(neurons, length(t));
   x(:, 1) = x0;
   
   otherAct = zeros(1, 10);
   
   for j = 1:length(t) - 1 % Iterating over time
       
       if t(j) > 1
           
           I = zeros(1, neurons);
           
       end
       
       for i = 1:neurons % Iterating over neurons
           
           otherAct(i) = 0; % Reset inhibitory inputs to 0
           
           for k = 1:neurons % Needed to calculate inhibitory inputs for each neuron
               
               if (i ~= k) % Sums activities for all k != i
                   otherAct(i) = otherAct(i) + signal(x(i, j), Signal);  
               end
               
           end
           
           x(i, j+1) = x(i, j) + dt*(-A*x(i,j) + (B - x(i,j))*(signal(x(i,j), Signal) + ...
               I(i)) - x(i,j)*otherAct(i)); % Euler's method
           
       end
       
   end
                             
end

% Activation functions
function [f] = signal(x, i)

    if i == 1
        f = x;
    end
    
    if i == 2
        f = x^2;
    end
    
    if i == 3
        f = x/(0.25 + x);
    end
    
    if i == 4
        f = (x^2)/(0.25 + x^2);
    end
    
end
        