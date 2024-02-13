from interface_builder import Silica


# prepare the surface 
silica_params = {
    'lx' : 15, 
    'ly' : 15, 
    'lz' : 15, 
    'vacuum' : 60,
    'filename' : 'silica_quartz.data',
    'output_folder' : 'quartz2',
    'input_folder' : 'test',
    'sio2_potential' : 'SiO2.vashishta',
    'sio2_h2o_potential' : "SiOH2O_199_16_adjusted.vashishta",
    'h2o_potential' : 'potential_nm_0_0.mod'
}

silica = Silica(**silica_params)
silica_file = silica.build_quartz()


# annealing procedure  
potential = 'SiO2.vashishta'
T = 100 # ps
silica_thermalize = silica.thermalize(silica_file, time=T, temp=300, output_filename= "silica_thermalize.data",
                                    mpirun_n=16, lmp_exec="lmp", run=True)


### Iteratively passivate the system in short burst of simulation to avoid breakage of siloxane bond
T = 1 # ps 
passivated_silica = silica_thermalize 
for i in range(2):
    passivated_silica = silica.passivate(passivated_silica , T, 'passivated_quartz.data',
                                    water_thickness=5, mpirun_n=1, lmp_exec='lmp_usc', run=True)
    print(i, passivated_silica)
    
    
silica_passivated_thermalize = silica.thermalize_passivated(passivated_silica, time=50, temp=300, 
                                                            output_filename= "silica_passivated_thermalize.data",
                                                            mpirun_n=16, lmp_exec="lmp_usc", run=True)


# silica_surface = silica.resize_z(passivated_silica,  "passivated_quartz.data")

       
# # set the desired silanol number 
# passivated_silica = 'passivated_quartz.data'
# silanol_silica = silica.set_silanol(0.2, passivated_silica, "passivated_quartz-2.5.data", 
#                                     lmp_exec='/home/users/anthonca/.bin/lmp_usc', delete_oxygen="yes")
# print("Silanol concentration:", silica.get_silanol(silanol_silica))
 

h2o_silica_system = silica.add_water(silica_passivated_thermalize, [None, 10, 10], num_mol=50, 
                                     thermalized_h2o=10, output_filename="silica_water_passivated.data", mpirun_n=16)

