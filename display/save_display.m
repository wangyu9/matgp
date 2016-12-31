function [im] = save_display(filename)
% save what is currently showed to image files
% example" save_display('vector-field.gif');


  im = myaa('raw');
  [imind,cm] = rgb2ind(im,256);
  
  imwrite(imind,cm,filename);
  
  %if w == 1;
  %imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
  %else
    % imwrite(imind,cm,filename,'gif','WriteMode','append');
  %end

