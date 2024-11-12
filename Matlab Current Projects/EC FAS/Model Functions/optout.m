%% optout
function stop = optout(x, optimValues, state)
% Output for running optimization

stop = false;

% Show current guess
disp(x)

% % Play sound after each iteration
%load handel
%sound(y,Fs)