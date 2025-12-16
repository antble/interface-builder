# Water contact angle simulation 

<!-- ![passivated_silica](../_static/wca_system.png) -->

```{image} ../_static/wca_system.png
:alt: passivated_silica
:width: 500px
:align: center
```
Figure above shows the final system obtained using the following script for water contact angle simulation. 

```python 
from interface_builder import Silica
import os

current_directory = os.getcwd()
print("Current working directory:", current_directory)


sim_params = {
    'lx' : 3, 
    'ly' : 3, 
    'lz' : 3, 
    'vacuum' : 20,
    'filename' : 'silica_amorphous.data',
    'output_folder' : f'{current_directory}',
    'input_folder' : None, 
    'sio2_potential' : f'{current_directory}/SiO2.vashishta',
    'sio2_h2o_potential' : f"{current_directory}/potential_hoga_5_0.mod",  
    'h2o_potential' : f'{current_directory}/potential_shrink_11_4.mod'
}

# initilize the silica system
silica = Silica(**sim_params)

# Replicate [lx,ly,lz], then do initial annealing
sio2_file = silica.build_amorphous(16, "lmp")
sio2_surface = silica.add_vacuum(sio2_file, 20, 
                                "silica_surface+vacuum.data")
sio2_relax = silica.thermalize(sio2_surface, time=1000, init_temp=300, final_temp=300,
                                    mpirun_n=16, lmp_exec="lmp")

sio2_replicated = silica.replicate_v2((2,3,1),sio2_relax, "silica_replicated-2x3x1")
sio2_vacuum = silica.resize_z(sio2_replicated, 80, "silica+vacuum.data")
sio2_vacuum_shifted = silica.shift_z_origin(sio2_vacuum, z_shift=4)

h2o_sio2_system = silica.add_water(sio2_vacuum_shifted, dimensions=[None,20,20], 
                                    num_mol=800,
                                    thermalized_h2o=0, 
                                    geometry="sphere",
                                    mpirun_n=16)
```