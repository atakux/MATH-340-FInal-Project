%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Runge-Kutta 4th Order Methods - Lotka-volterra Model           %
%                                                                         %
%       Names:                   Emails:                                  %
%       Kris Melgar Morales      cmelgarmorales@csu.fullerton.edu         %
%       Armanul Ambia            arman714@csu.fullerton.edu               %
%       Angela DeLeo             atakux707@csu.fullerton.edu              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default Parameters
a = 1; %Growth rate per capita of prey
b = .01; %Predator population effects on growth rate of prey
c = 1; %Effect of prey population on preadtor growth rate
d = .02; %Death rate of preadtor per capita

% Initial conditions
R0 = 20;
F0 = 20;

% Time interval
tspan = [0 15];
h = 1/300;
N = (tspan(2)-tspan(1))/h; % Total Steps

% System of ODEs
dYdt = @(t,Y) [a*Y(1) - b*Y(1)*Y(2); -c*Y(2) + d*Y(1)*Y(2)];

% 4th order Runge-Kutta method
[tW,W] = ode45(dYdt, tspan, [R0; F0], h);
[tY,Y] = RK4(dYdt, tspan, [R0; F0], N);
length(tY)

% Plot results of Density
figure('units','normalized','outerposition',[0 0 1 1])
plot(tY,Y,tW,W)
legend('Prey RK4','Predator RK4','Prey ODE45','Predator ODE45','FontSize',48,'interpreter','latex','Location','north')
xlabel('Time','FontSize',48,'interpreter','latex')
ylabel('Population Density','FontSize',48,'interpreter','latex')
print -depsc2 Density.eps

% Plot results of Phase
figure('units','normalized','outerposition',[0 0 1 1])
plot(Y(1,:),Y(2,:),'-',W(:,1),W(:,2),'-')
legend('RK4','ODE45','FontSize',48)
xlabel('Prey Population Density','FontSize',48,'interpreter','latex')
ylabel('Predator Population Density','FontSize',48,'interpreter','latex')
print -depsc2 Phase.eps

% Plot results on subplot
% Create labels for sliders
%hTextA = uicontrol('Style', 'text', 'String', strcat('a: ',num2str(a)), 'Position', [230 20 40 20]);
%hTextB = uicontrol('Style', 'text', 'String', strcat('b: ',num2str(b)), 'Position', [230 50 40 20]);
%hTextC = uicontrol('Style', 'text', 'String', strcat('c: ',num2str(c)), 'Position', [230 80 40 20]);
%hTextD = uicontrol('Style', 'text', 'String', strcat('d: ',num2str(d)), 'Position', [230 110 40 20]);

%Init sliders
%hSliderA = uicontrol('Style', 'slider', 'Min', 0, 'Max', 100, 'Value', a, 'SliderStep', [1/100 1/100], 'Position', [20 20 200 20]);
%hSliderB = uicontrol('Style', 'slider', 'Min', 0, 'Max', 100, 'Value', b, 'SliderStep', [1/100 1/100], 'Position', [20 50 200 20]);
%hSliderC = uicontrol('Style', 'slider', 'Min', 0, 'Max', 100, 'Value', c, 'SliderStep', [1/100 1/100], 'Position', [20 80 200 20]);
%hSliderD = uicontrol('Style', 'slider', 'Min', 0, 'Max', 100, 'Value', d, 'SliderStep', [1/100 1/100], 'Position', [20 110 200 20]);

%Init callback for sliders
%set(hSliderA,'Callback',{@update,tspan,R0,F0,N,h,hSliderA,hSliderB,hSliderC,hSliderD,hTextA,1});
%set(hSliderB,'Callback',{@update,tspan,R0,F0,N,h,hSliderA,hSliderB,hSliderC,hSliderD,hTextB,2});
%set(hSliderC,'Callback',{@update,tspan,R0,F0,N,h,hSliderA,hSliderB,hSliderC,hSliderD,hTextC,3});
%set(hSliderD,'Callback',{@update,tspan,R0,F0,N,h,hSliderA,hSliderB,hSliderC,hSliderD,hTextD,4});

%Slider update function
function update(hObj,event,tspan,R0,F0,N,h,hSliderA,hSliderB,hSliderC,hSliderD,hText,id)
    st = ["a:","b:","c:","d:"];
    a = get(hSliderA,'Value');
    b = get(hSliderB,'Value');
    c = get(hSliderC,'Value');
    d = get(hSliderD,'Value');
    params = [a,b,c,d];
    set(hText,'String',strcat(st(id),num2str(params(id))));
    updatePlot(tspan,R0,F0,N,h,a,b,c,d);
end
    function updatePlot(tspan,R0,F0,N,h,a,b,c,d)
        % Recalculate solution with new parameter values
        dYdt = @(t,Y) [a*Y(1) - b*Y(1)*Y(2); -c*Y(2) + d*Y(1)*Y(2)];
        
        %RK4 methods
        [tW,W] = ode45(dYdt, tspan, [R0; F0], h);
        [tY,Y] = RK4(dYdt, tspan, [R0; F0], N);
        
        % Update plot with new solution
        plot(tY,Y,tW,W)
        legend('Rabbits','Foxes')
        xlabel('Time')
        ylabel('Population')
    end
    
%Runga Kutta 4th Order Function
%input: dYdt as vector of anonymous function, tspan as an array of length 2
%containing an start point and end point,
%Y0 as an array matching the length of the dYdt, containing the intital conditions
%N as the amount of places to evaluate the interval at
%output: t as an array containg the different timesteps evaluated at,
%Y as the matrix where each row contains the solution of the system of
%equations with respect to t array.
function [t,y] = RK4(f, tspan, y0, N)
    h = (tspan(2)-tspan(1))/N;
    t = (tspan(1):h:tspan(2))';
    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    for i = 1:N
        k1 = h * f(t(i), y(:,i));
        k2 = h * f(t(i) + h/2, y(:,i) + k1/2);
        k3 = h * f(t(i) + h/2, y(:,i) + k2/2);
        k4 = h * f(t(i) + h, y(:,i) + k3);
        y(:,i+1) = y(:,i) + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end
