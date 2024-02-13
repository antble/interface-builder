from interface_builder import Silica

# 3 --> 21.36x21.36x21.36
# 4 --> 28.48x28.48x28.48 



# Generate 10 sample of silica surface! 
# for indx in range(10):

# prepare the surface 
silica_params = {
    'lx' : 3, 
    'ly' : 3, 
    'lz' : 3, 
    'vacuum' : 30,
    'filename' : 'silica_amorphous.data',
    'output_folder' : 'amorphous_passivated',
    'input_folder' : 'test', 
    'sio2_potential' : 'SiO2.vashishta',
    'sio2_h2o_potential' : "SiOH2O_199_16_adjusted.vashishta",
    'h2o_potential' : 'potential_nm_0_0.mod'
}
#================================= PASSIVATION OF THE SILICA SYSTEM ==================== 
# initilize the silica system
silica = Silica(**silica_params)

# replicate [lx,ly,lz], then do initial annealing
silica_file = silica.build_amorphous(16, "lmp", run=False)
silica_surface = silica.add_vacuum(silica_file, "silica_surface+vacuum.data")

silica_relax = silica.thermalize(silica_surface, time=50, temp=300,
                                output_filename="amorphous_surface_thermalized-300.data",
                                mpirun_n=16, lmp_exec="lmp", run=False)

minimized_silica = silica.minimize_sio2(silica_relax, mpirun_n=16, lmp_exec="lmp", run=False)
# time=[ps], themp=[K] | relax the system for 250 ps
# silica_surface = silica.resize_z(silica_file, 90, "silica_surface.data")

silica_annealed = silica.anneal(minimized_silica, mpirun_n=16, lmp_exec="lmp", run=False)


# # #================================= PASSIVATION OF THE SILICA SYSTEM ==================== 
# # passivate a surface
# '''WARNING: there maybe some excess atoms 
# manually remove the excess atoms in the passivated.data 
#     - still require human intervention! :)
# '''

passivated_silica = silica_annealed
# for indx in range(1):
passivated_silica = silica.passivate(passivated_silica, 1, 'passivated_silica.data',
                                    water_thickness=5, mpirun_n=1, lmp_exec='lmp_usc', run=False)


# # # set the desired silanol number 
silanol_silica = silica.set_silanol(1.0, passivated_silica, "passivated_silanol-2.7.data", delete_oxygen="yes")


# #================================= WATER + SILICA SYSTEM =================================
# h2o_silica_system = silica.add_water(passivated_silica, [10,None,10])




