from interface_builder import Silica

# Generate 10 sample of silica surface! 
for indx in range(2):
    # prepare the surface 
    silica_params = {
        'lx' : 4, 
        'ly' : 4, 
        'lz' : 4, 
        'vacuum' : 30,
        'filename' : 'silica_amorphous.data',
        'output_folder' : f'amorphous_{indx}',
        'input_folder' : 'test', 
        'sio2_potential' : './data/SiO2.vashishta',
        'sio2_h2o_potential' : "./data/SiOH2O_199_16_adjusted.vashishta",
        'h2o_potential' : './data/potential_nm_0_0.mod'
    }
    #================================= PASSIVATION OF THE SILICA SYSTEM ==================== 
    # initilize the silica system 
    silica = Silica(**silica_params)

    silica_file = silica.build_amorphous(16, "lmp", run=True)
    
    # silica_surface = silica.resize_z(silica_file, 90, "silica_surface.data")
    # silica_annealed = silica.anneal(silica_surface, mpirun_n=16, lmp_exec="lmp", run=True)


