function [X,Y] = Observing_filter(public,CASE,S_temp,mu)

X = feval(public.function.HX,  CASE.C, S_temp, mu);
Y = feval(public.function.HV,  CASE.C, CASE.A,  S_temp);