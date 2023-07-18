# surrogate_models

Steps to run SDS_main.jl

1. Install Julia and ITensor:
https://itensor.github.io/ITensors.jl/dev/getting_started/Installing.html

2. To enable multi-threading in Julia use in a terminal e.g.

julia -t 8

or, if having previously compiled the ITensor library

julia -t 8 --sysimage C:\Users\your_username\.julia\sysimages\sys_itensors.so

See also https://itensor.github.io/ITensors.jl/dev/getting_started/RunningCodes.html

3. In the Julia session run

julia> include("SDS_main.jl")

Note: the code has been written using VS Code https://code.visualstudio.com/ 
and uses unicode math symbols that may not be displayed correctly in other editors.

