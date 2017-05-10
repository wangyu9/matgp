function writeMTL(fnameMTL)

f = fopen( fnameMTL, 'w' );

fprintf( f, 'newmtl lightred\n');
fprintf( f, 'Ka 0.1 0.0 0.0\n');
fprintf( f, 'Kd 0.4 0.0 0.0\n');
fprintf( f, 'Ks 0.0 0.0 0.0\n');
fprintf( f, '\n');

fprintf( f, 'newmtl lambert4SGtexture\n');
fprintf( f, 'Kd 0.4 0.0 0.0\n');
fprintf( f, 'Ka 0.0 0.0 0.0\n');
fprintf( f, 'Ks 0.0 0.0 0.0\n');
%fprintf( f, 'Tf 1 1 1\n');
fprintf( f, 'map_Kd checkerboard.png\n');
fprintf( f, 'Ni 1.00\n');
fprintf( f, '\n');

fprintf( f, 'newmtl white\n');
fprintf( f, 'Ka 0 0 0\n');
fprintf( f, 'Kd 1 1 1\n');
fprintf( f, 'Ks 0 0 0\n');
fprintf( f, '\n');

fclose(f);
