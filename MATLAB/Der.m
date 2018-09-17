function [dy,xc] = Der(F,x,varargin)
% D calculates the derivative of the function F(x) using one of three
% methods: forward difference, central difference, and expolated difference
% D(F,x,method) returns the x and y coordinate of the derivative function.
% Specify which method in the optional third argument as 'fd', 'cd', or
% 'ed'. If F is a numeric vector, it is treated as data and a simple
% derivative is calculated

if isempty(varargin)
   % Default method
    varargin{1} = 'cd'; 
end
if ~isnumeric(x)
    error('Second input must be numeric vector')
end
if ~isa(F,'function_handle')
    if isnumeric(F)
      xc = chop(x,0); % method for data depends on how you assign dy to xc
      dy = diff(F)./diff(x); 
      varargin{1} = 'data';
    else
      error('First input must be a function handle or a numeric array')
    end
end

n=length(x);
dx = diff(x);

switch varargin{1}

    case 'df'     % 1st derivative by forward differentiation
        dy = (F(x(2:n)) - F(x(1:n-1)))./dx;
        xc = x(1:n-1);
    case 'db'     % 1st derivative by backward differentiation
        dy = (F(x(2:n)) - F(x(1:n-1)))./dx;
        xc = x(2:n);
    case 'dc'     % 1st derivative by central differentation
        dy = (F(x(1:n-1)+0.5*dx) - F(x(1:n-1)-0.5*dx))./dx;
        xc = x(1:n-1);
    case 'de'     
        half = (F(x(1:n-1)+0.25*dx) - F(x(1:n-1)-0.25*dx))./(0.5*dx); 
        full = (F(x(1:n-1)+0.5*dx) - F(x(1:n-1)-0.5*dx))./dx;
        dy = (4/3).*half - (1/3).*full;
        xc = chop(x);
    case '2dc'    % 2nd derivative by central differentiation
        dy = (0.5*F(x(1:n-1)+ dx) - 0.5*F(x(1:n-1)-dx))./dx;
        xc = x(1:n-1);
    case '3dc'    % 3rd derivative by central differentiation
        dy = (27*F(x(1:n-1)+0.5*dx) - 27*F(x(1:n-1)-0.5*dx) + F(x(1:n-1)-1.5*dx) - F(x(1:n-1)+1.5*dx))./(24*dx);
        xc = x(1:n-1);
    case '4dc'    % 4th derivative by central differentiation
        dy = (F(x(1:n-1)-2*dx) -8*F(x(1:n-1)-dx) + 8*F(x(1:n-1)+dx) - F(x(1:n-1)+2*dx) )./(12.*dx);
        xc = x(1:n-1);
    case '5dc'    % 5th derivative by central differentiation
        dy = ((3*F(x(1:n-1)+2.5*dx)./640) -(25*F(x(1:n-1)+1.5*dx)./384)+(75*F(x(1:n-1)+0.5*dx)/64)-(3*F(x(1:n-1)-2.5*dx)./640) +(25*F(x(1:n-1)-1.5*dx)./384)-(75*F(x(1:n-1)-0.5*dx)/64))./dx;
        xc = x(1:n-1);
    case '2df'    % 2nd derivative by forward differentiation
        dy = (F(x(1:n-1)+dx)+F(x(1:n-1)-dx)-2 .*F(x(1:n-1))) ./ dx .^2;
        xc = x(1:n-1);
    case '2db'    % 2nd derivative by backward differentiation
        dy = (F(x(1:n-1)+dx)+F(x(1:n-1)-dx)-2 .*F(x(1:n-1))) ./ dx .^2;
        xc = x(2:n);
    case 'data'
        return;
    otherwise
        error('method not recognized');
        
end

end