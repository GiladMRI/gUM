%	function [r1,r2] = qdf(a,b,c)
%
%	Outputs quadratic roots of ax^2+bx+c = 0.
%


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: qdf.m,v $
%	Revision 1.1  2002/03/28 01:27:46  bah
%	Added to CVS
%	
%
% ===========================================================


function [roots] = qdf(a,b,c)

d = b^2 - 4*a*c;

roots(1) = (-b + sqrt(d))/(2*a);
roots(2) = (-b - sqrt(d))/(2*a);




