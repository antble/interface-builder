from interface_builder import Silica


# prepare the surface 
silica_params = {
    # 'polymorph' : 'quartz1',
    'lx' : 30, 
    'ly' : 100, 
    'lz' : 30, 
    'vacuum' : 60,
    'filename' : 'silica_quartz.data',
    'output_folder' : 'quartz1',
    'input_folder' : 'test',
    'sio2_potential' : 'SiO2.vashishta',
    'sio2_h2o_potential' : "SiOH2O_199_16_adjusted.vashishta",
    'h2o_potential' : 'potential_nm_0_0.mod'
}

silica = Silica(**silica_params)
silica_file = silica.build_quartz()


# annealing procedure  
potential = 'SiO2.vashishta'
T = 5 # ps
silica_thermalize = silica.thermalize(silica_file, time=T, temp=300, output_filename= "silica_thermalize.data",
                                    mpirun_n=16, lmp_exec="lmp", run=False)


### Iteratively passivate the system in short burst of simulation to avoid breakage of siloxane bond
T = 1 # ps 
passivated_silica = silica_thermalize 
for _ in range(3):
    passivated_silica = silica.passivate(passivated_silica , T, 'passivated_quartz.data',
                                    water_thickness=5, mpirun_n=16, lmp_exec='lmp_usc', run=False)

silica_passivated_thermalize = silica.thermalize_passivated(passivated_silica, time=50, temp=300, 
                                                            output_filename= "silica_passivated_thermalize.data",
                                                            mpirun_n=16, lmp_exec="lmp_usc", run=True)


# silica_surface = silica.resize_z(passivated_silica,  "passivated_quartz.data")

print("Silanol concentration:", silica.get_silanol(silica_passivated_thermalize))
        
# set the desired silanol number 
# silanol_silica = silica.set_silanol(0.2, passivated_silica, "passivated_silanol-2.5.data", 
#                                     lmp_exec='lmp_usc', delete_oxygen="yes")


h2o_silica_system = silica.add_water(silica_passivated_thermalize, [None,20,20], num_mol=512, thermalized_h2o=50, mpirun_n=16)

