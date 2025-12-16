# Passivated Amorphous Silica Surface


```{image} ../_static/passivated_asilica0002.png
:alt: passivated_silica
:width: 200px
:align: center
```


```python
from interface_builder import Silica

# 3 --> 21.36x21.36x21.36
silica_params = {
    'lx' : 3, 
    'ly' : 3, 
    'lz' : 3, 
    'filename' : 'silica_amorphous.data',
    'output_folder' : './passivation_test/',
    'input_folder' : None, 
    'sio2_potential' : './data/SiO2.vashishta',
    'sio2_h2o_potential' : "./passivation_test/potential_hoga_5_0.mod",
    'h2o_potential' : './data/potential_initial_0_0-A.mod' 
}

# initilize the silica system
silica = Silica(**silica_params)

silica_file = silica.build_amorphous(16, "lmp", run=False)
silica_surface = silica.add_vacuum(silica_file, 20, 
                                "silica_surface+vacuum.data", 
                                which="both")

silica_relax = silica.thermalize(silica_surface, 
                                time=2000, 
                                init_temp=300, final_temp=300, 
                                output_filename="amorphous_surface_thermalized-300ps.data",
                                mpirun_n=16, 
                                lmp_exec="lmp")


passivated_silica = silica.passivate_new(silica_relax, 
                                'passivated_silica.data')
```
