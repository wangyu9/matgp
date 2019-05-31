function writeECS(fnameECS, fnameOBJ, argu)


f = fopen( fnameECS, 'w' );

fprintf( f, '-i %s\n', fnameOBJ);
fprintf( f, '\n');
fprintf( f, sprintf('-size %d %d \n',[argu.image_size(1),argu.image_size(2)]));

fprintf( f, '#-vp 0 1 -5 -vi 0.0 0.0 0 -vu 0 1 -0.2 -fov 11\n');
%fprintf( f, '#-vp 0 5 0 -vi 0.0 0.0 0 -vu 1 0 0 -fov 11\n');
fprintf( f, '#-ambientlight 1 1 1\n');
fprintf( f, '#-directionallight 0 -1 0 .2 .2 .2\n');

fclose(f);
