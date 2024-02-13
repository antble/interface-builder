from interface_builder import Silica


# prepare the surface 
silica_params = {
    'vacuum' : 0,
    'filename' : 'Q4_siloxane.data',
    'output_folder' : 'Q4',
    'input_folder' : 'test',
    'sio2_potential' : 'SiO2.vashishta',
    'sio2_h2o_potential' : "SiOH2O_199_16_adjusted.vashishta",
    'h2o_potential' : 'potential_nm_0_0.mod'
}

silica = Silica(**silica_params)
silica_filename = silica() 
print(silica_filename)

h2o_silica_system = silica.add_water(silica_filename, [None,20,20], 
                                    center=(33.3675/2, 104.526/2, 20+21),
                                    num_mol=512, 
                                    thermalized_h2o=50, 
                                    mpirun_n=16)

