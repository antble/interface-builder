# Adding Vacuum 

```{image} ../_static/vacuum_asilica0000.png
:alt: passivated_silica
:width: 200px
:align: center
```


```python
from interface_builder import Silica

silica_params = {
    'lx' : 3, 
    'ly' : 3, 
    'lz' : 3, 
    'filename' : 'silica_amorphous.data',
    'output_folder' : './passivation_test/',
    'input_folder' : None, 
    'sio2_potential' : './data/SiO2.vashishta',
}

# initilize the silica system
silica = Silica(**silica_params)

# Replicate [lx,ly,lz], then do initial annealing
silica_file = silica.build_amorphous(16, "lmp")
silica_surface = silica.add_vacuum(silica_file, 20, 
                            "silica_surface+vacuum.data", 
                            which="both")
```