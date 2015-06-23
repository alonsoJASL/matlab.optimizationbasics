function [xstar, it, kfail] = BFGS(fname, x0, maxITER, varargin)
%           BROYDEN-FLETCHER-GOLDFARB-SHANNO (BFGS)
%
% Genera la actualizacion de Broyden-Fletcher-Goldfarb-Shanno simetrico
% positivo definido. Los sistemas lineales se resuelven utilizando el
% metodo de Gradiente conjugado. Usa SPARSE MATRICES, y toma en cuenta los
% argumentos extra que pueda necesitar la funcion fname.
%
% usage:     [xout, it, kfail] = BFGS(fname, xin, maxITER)
%           donde:
%                SALIDA
%                  xstar := aproximacion al minimo local.
%                     it := numero de iteraciones tomadas por el metodo
%                  kfail := numero de veces en que el metodo falla en
%                           generar una direccion de descenso.
%                ENTRADA
%                  fname := nombre de la funcion a evaluar (string)
%                     x0 := punto de inicio del metodo.
%                maxiter := numero maximo de iteraciones 
%                           (DEFAULT maxiter = 100)
%
% 

%
% Jose Alonso Solis Lemus (2012. ITAM) 
% for license and more code check: 
% https://github.com/alonsoJASL/matlab.optimizationbasics.git
%

%---- Algunas variables necesarias:-----
if nargin < 3
    maxITER = 100;
end
tol = 10e-6; % tolerancia
n = length(x0); 

%---- Para las condiciones de Wolfe ----:
maxW = 30;
c = [10e-4;0.9];

%---- Inicializacion --------------------
xk = x0;

B = speye(length(x0)); %Hessiana(fname, x0);
Bk = B;
k = 0;
kfail = 0;
termine = false;

while ~termine
    
    [fxk, gk] = feval(fname, xk, varargin{:});
    
    sk = GradienteConjugado(Bk,-gk); 
    
    %--------- Condiciones de Wolfe ------------
    [alpha] = condicionesWOLFE(fname, xk, fxk, sk, gk, c, maxW, ...
        varargin{:});
    xk = xk + alpha*sk;
    
    [~, yaux] = feval(fname, xk, varargin{:});
    yk = yaux - gk;
    
    % Actualizacion de BFGS
    if (cond(B)>100)
        Bk = speye(n);
        kfail = kfail+1;
    else
        Bk = Bk - ((Bk*(sk*sk')*Bk)/(sk'*Bk*sk))...
            + ((yk*yk')/(yk'*sk));
    end
    k = k+1;
    termine = (k >= maxITER || norm(gk) <= tol);
    disp(norm(gk));
end

xstar = xk;
it = k;

