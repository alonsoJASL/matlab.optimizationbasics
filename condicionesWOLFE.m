function [alpha] = condicionesWOLFE(fname, xk, fxk, sk, gk, c, maxW, ...
    varargin)
%             CONDICIONES DE WOLFE
%
% Metodo auxiliar a BFGS para encontrar el valor alpha necesario para
% que se tenga un paso apropiado en la busqueda de linea. Dado que es 
% un metodo auxiliar, toma en cuenta valores como sk, necesario para 
% obtener la derivada direccional. 
%
% usage:       [alpha] = condicionesWOLFE(xk, sk, gk, maxW)
%          donde:
%               SALIDA
%                   alpha := valor entre 0 y 1 que define el paso
%                            a dar, dadas las condiciones de Wolfe
%               ENTRADA
%                   fname := Nombre de la funcion (string)
%                     fxk := Valor de f en el punto xk
%                      xk := punto de inicio 
%                      sk := vector de direccion de descenso
%                      gk := gradiente de f en el punto xk
%                       c := vector que contiene las constantes 
%                            c1 y c2 de las condiciones de Wolfe
%                    maxW := maximo numero de iteraciones
%                   
%

%
% Jose Alonso Solis Lemus (2012. ITAM) 
% for license and more code check: 
% https://github.com/alonsoJASL/matlab.optimizationbasics.git
%


alpha = 1; 
xp = xk+ alpha*sk;       % Punto de prueba
[fxp, gradXp] = feval(fname, xp, varargin{:});  % Funcion en el punto de prueba
s = sk'*gk; % Derivada direccional
iterW = 0;

termine = (fxp > fxk+(alpha*c(1)*s)) || ...
        (gradXp'*sk < c(2)*gk'*sk); % condicion de paro;


    
while ~termine
    alpha = alpha/2;
    xp = xk + alpha*sk;      % Punto de prueba
    fxp = feval(fname, xp, varargin{:});
    iterW = iterW+1;       
    
    termine = ~(fxp > fxk+(alpha*c(1)*s)) || ...
        (gradXp'*sk < c(2)*gk'*sk)&&...
        (iterW<maxW); % condicion de paro
end