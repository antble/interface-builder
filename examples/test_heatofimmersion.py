from interface_builder import Silica
'''Methodology: 
An amorphous silica was used in the interface simulation. We used N number of water molecule to calculate the heat of immersion.
The water was thermalized for 20ps. Then a pure silica with a vacuum of X angstrom was thermalized for 100 ps. 
'''
# prepare the surface 
silica_params = {
    'lx' : 3, 
    'ly' : 3, 
    'lz' : 3, 
    'filename' : 'silica_amorphous.data',
    'output_folder' : f'heatofimmersion5',
    'sio2_potential' : 'SiO2.vashishta',
    'sio2_h2o_potential' : "potential_hoga_0_0.mod",#"SiOH2O_nm_0_0.vashishta", #"SiOH2O_199_16_adjusted.vashishta",#
    'h2o_potential' : "potential_nm_0_0.mod", #"hoga_199_16.mod"
}
#================================= CALCULATION OF THE WATER/SILICA SYSTEM HEAT OF IMMERSION==================== 
# initilize the silica system 
silica = Silica(**silica_params)

# replicate [lx,ly,lz], then do initial annealing
silica_file = silica.build_amorphous(16, "lmp", run=False)


N = int(586/2)   # number of water 
dH_imm = silica.get_heat_of_immersion(silica_file, 
                                      num_h2o=N, 
                                      h2o_therm_time=20,
                                      sio2_therm_time=100,
                                      output_filename="water-silica_system.data",
                                      mpirun_n=16,
                                      run=True
                                      )
print("Calculated heat of immersion:", dH_imm, "J/m^2")
