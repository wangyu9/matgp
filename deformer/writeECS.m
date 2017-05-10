function writeECS(fnameECS, fnameOBJ)


f = fopen( fnameECS, 'w' );

fprintf( f, '-i %s\n', fnameOBJ);
fprintf( f, '\n');
fprintf( f, '-vp 1.2 0.1 -0.6 -vi 0.0 0.0 0 -vu 0 1 0 -fov 37\n');
fprintf( f, '-size 1600 1200\n');
fprintf( f, '-ambientlight 1 1 1\n');

fclose(f);
