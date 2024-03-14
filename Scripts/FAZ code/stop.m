function F = stop(string)

% Messagebox:
M = msgbox(string,'Stop');
pos1 = [50, 50,150,80];
set(M,'OuterPosition',pos1) 

% create the two anonymous functions
F.Stop = @() stopfun(M) ; % false if message box still exists
F.Clear = @() clearfun(M) ; % delete message box

function r = stopfun(M)
drawnow ;          % ensure that button presses are recorded
r = ~ishandle(M) ; % false if message box still exists

function clearfun(M)
% clear the message box if it still exists
if ishandle(M),
    delete(M) ;
end