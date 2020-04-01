function [th_s,th_n] = MCMT(sl,st,fn)
% Minimum controled maximum thresholds
% wl: smoothing length (in frame)
% st: smoothing increment (in frame)
% fn: frame number 

%% Erlang Distribution
 intval = 0.01;
 cd = zeros(length(0.00:intval:3),1);
 it = 0;
 for x = sl*(0.00:intval:3)
     it = it+1;
     value = 1;
     itr = 1;
     for n = 1:sl-1
         value = value*x/n;
         itr = itr+value;
     end
     cd(it) = 1-exp(-x)*itr;
 end
    
 G = sl*(0.00:intval:3);
 GC = cd';
 
%% Minimum controled maximum thresholds 

 bl = fn*st*(1+log(sl/st))/sl;
 pl = log(G)*(sl-1)-G-log(factorial(sl-1))+log(1-GC)*(bl-1);
 pu = log(G)*(sl-1)-G-log(factorial(sl-1))+log(GC)*(bl-1);       
 el = sum(exp(pl).*G)/sum(exp(pl));
%  eu = sum(exp(pu).*G)/sum(exp(pu));    
 Fu = cumsum(exp(pu))/sum(exp(pu));
 th_s = G(find(Fu>0.9,1))/el;
 th_n = G(find(Fu>0.5,1))/el;
 