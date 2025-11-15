% buildit
% mex -setup C++

arch = computer('arch');

if strcmp(arch,'maci64')
    mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka1D.cpp
    mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka2D.cpp
    mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka3D.cpp
    mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand1D.cpp
    mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand2D.cpp
    mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand3D.cpp
elseif strcmp(arch,'maca64')
    mex -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -L/opt/homebrew/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' LDFLAGS='$LDFLAGS -ld_classic' shinozuka1D.cpp
    mex -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -L/opt/homebrew/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' LDFLAGS='$LDFLAGS -ld_classic' shinozuka2D.cpp
    mex -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -L/opt/homebrew/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' LDFLAGS='$LDFLAGS -ld_classic' shinozuka3D.cpp
    mex -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -L/opt/homebrew/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' LDFLAGS='$LDFLAGS -ld_classic' shinozukaRand1D.cpp
    mex -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -L/opt/homebrew/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' LDFLAGS='$LDFLAGS -ld_classic' shinozukaRand2D.cpp
    mex -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 -L/opt/homebrew/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' LDFLAGS='$LDFLAGS -ld_classic' shinozukaRand3D.cpp
end