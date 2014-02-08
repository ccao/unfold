set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette defined (0  1  1  1, \
                     1  0.875  0.875  1.000, \
                     2  0.750  0.750  1.000, \
                     3  0.625  0.625  1.000, \
                     4  0.500  0.500  1.000, \
                     5  0.375  0.375  1.000, \
                     6  0.250  0.250  1.000, \
                     7  0.125  0.125  1.000, \
                     8  0.000  0.000  1.000 )
#set palette defined ( 0 '#000090',\
#                      1 '#000fff',\
#                      2 '#0090ff',\
#                      3 '#0fffee',\
#                      4 '#90ff70',\
#                      5 '#ffee00',\
#                      6 '#ff7000',\
#                      7 '#ee0000',\
#                      8 '#7f0000')
plot 'spec.dat' u 1:($2)/1000:3 matrix with image t ''
