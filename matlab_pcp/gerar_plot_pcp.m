%%
% This is an example of how to create a 3D bar chart in MATLAB&#174;.
% 
% Read about the <http://www.mathworks.com/help/matlab/ref/bar3.html |bar3|> function in the MATLAB documentation.
%
% For more examples, go to <http://www.mathworks.com/discovery/gallery.html MATLAB Plot Gallery>
%
% Copyright 2012-2014 The MathWorks, Inc.

M = csvread('tempos.csv')
matrix = reshape (M, [32,64])
figure

b = bar3( matrix)

axis tight

% Add title and axis labels
title('{Rela\c{c}\~ao entre Processos MPI e Threads OpenMP para Algoritmo H\''ibrido :: Matriz 2048*2048}','interpreter','latex')
%title('Rela??o entre Processos MPI e Threads OpenMP para Algoritmo H?brido :: Matriz 2048*2048')
xlabel('# Processos MPI');
ylabel('# Threads OpenMP');
zlabel('Tempo Total (ms)');



for k = 1:length(b)
   zdata = b(k).ZData;
  b(k).CData = zdata;
 b(k).FaceColor = 'interp';
end

colorbar;

