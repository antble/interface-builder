from molecular_builder.core import carve_geometry, create_bulk_crystal
from molecular_builder.geometry import SphereGeometry, CylinderGeometry, BoxGeometry
from molecular_builder import pack_water, write 

import ase.io as io
import numpy as np 
import os 



script_dir = os.path.dirname(os.path.abspath(__file__))



class Silica:
    def __init__(self, lx, ly, lz, output_folder=None, input_folder=None, 
                vacuum=None, sio2_potential=None, sio2_h2o_potential=None, 
                h2o_potential=None,  filename="silica_quartz.data"):
        self.lx = lx 
        self.ly = ly 
        self.lz = lz 
        self.vacuum = vacuum
        self.filename = filename
        self.output_folder = output_folder
        self.input_folder = input_folder 
        self.h2o_potential = h2o_potential
        self.sio2_potential = sio2_potential
        self.sio2_h2o_potential = sio2_h2o_potential
        
        # create the output directory
        os.system(f'mkdir {output_folder}')
        
        
    def build_quartz(self):
        '''Create a quartz silica slab. 
        Process the quartz crystal structure and write it to a file in LAMMPS data format.
        see: https://henriasv.github.io/molecular-builder for more info. 
        
        :return: The filename of the output LAMMPS data file.
        :rtype: str
        '''
        L = [self.lx, self.ly, self.lz]
        quartz_structure = create_bulk_crystal("alpha_quartz", L) 
        cell = quartz_structure.get_cell()
        # resetting off-diagonal elements 
        cell[0,1] = 0 
        cell[1,0] = 0 
        cell[0,2] = 0 
        cell[2,0] = 0 
        cell[1,2] = 0
        cell[2,1] = 0
        quartz_structure.set_cell(cell) 
        quartz_structure.wrap()
        # adjusting cell size
        quartz_structure.cell[2,2] = self.lz + self.vacuum
        # write to a file
        quartz_structure.write(f'{self.output_folder}/{self.filename}', format="lammps-data")
        return self.filename
    
    
    
    def thermalize(self, datafile, time, temp, output_filename="silica-thermalize.data", mpirun_n=4, lmp_exec='lmp', run=True):
        '''Run a thermalization procedure on the system.

        :param datafile: The filename of the input data file.
        :type datafile: str
        :param time: The time duration for thermalization.
        :type time: float
        :param temp: The temperature for thermalization.
        :type temp: float
        :param output_filename: (Optional) The filename for the output data file.Default is silica_thermaliz.data.
        :type output_filename: str or None
        :param mpirun_n: (Optional) The number of MPI processes to use for the thermalization simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is 'lmp'.
        :type lmp_exec: str
        :param run: (Optional) Whether to run the thermalization simulation. Default is True.
        :type run: bool
        
        Returns
        :return: The filename of the output data file, if provided.
        :rtype: str or None
        '''
        input_filepath = os.path.join(self.output_folder, datafile)
        output_filepath = os.path.join(self.output_folder, output_filename) if output_filename else None

        print(f'{self.output_folder}/{datafile}')
        lmp_args = {
            'input_filename' : input_filepath, 
            'potential_filename' : self.sio2_potential, 
            'output_filename' : output_filepath,
            'therm_time' : time,
            'therm_temp' : temp, 
        }
        if run:
            script_path = os.path.join(script_dir, 'script', 'in.thermalize') 
            self.execute_lammps(lmp_args, mpirun_n, lmp_exec, script=script_path)
        return output_filename
    
    
    
    def thermalize_passivated(self, datafile, time, temp, output_filename="silica_passivated-thermalize.data", mpirun_n=4, lmp_exec='lmp', run=True):
        """Run a thermalization procedure on the system.

        :param datafile: The filename of the input data file.
        :type datafile: str
        :param time: The time duration for thermalization.
        :type time: float
        :param temp: The temperature for thermalization.
        :type temp: float
        :param output_filename: (Optional) The filename for the output data file. Default is silica_passivated-thermalize.data. 
        :type output_filename: str or None
        :param mpirun_n: (Optional) The number of MPI processes to use for the thermalization simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is 'lmp'.
        :type lmp_exec: str
        :param run: (Optional) Whether to run the thermalization simulation. Default is True.
        :type run: bool

        :return: The filename of the output data file, if provided.
        :rtype: str or None
        """
        input_filepath = os.path.join(self.output_folder, datafile)
        output_filepath = os.path.join(self.output_folder, output_filename) if output_filename else None

        lmp_args = {
            'input_filename' : input_filepath, 
            'potential_filename' : self.sio2_h2o_potential, 
            'output_filename' : output_filepath,
            'therm_time' : time,
            'therm_temp' : temp, 
        }
        if run:
            script_path = os.path.join(script_dir, 'script', 'in.thermalize_passivated') 
            self.execute_lammps(lmp_args, mpirun_n, lmp_exec, script=script_path)
        return output_filename
    
    
    
    def passivate(self, datafile, passivation_time, output_file, water_thickness=10, run=True, mpirun_n=4, lmp_exec='lmp_mpi'):
        '''Run a passivation simulation.

        :param datafile: The filename of the input data file.
        :type datafile: str
        :param passivation_time: The duration of the passivation simulation.
        :type passivation_time: float
        :param output_file: The filename for the output data file.
        :type output_file: str
        :param thickness: (Optional) The thickness of the water slab to be immersed in the silica surface. Default is 10.
        :type thickness: int
        :param run: (Optional) Whether to run the passivation simulation. Default is True.
        :type run: bool
        :param mpirun_n: (Optional) The number of MPI processes to use for the passivation simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is 'lmp_mpi'.
        :type lmp_exec: str

        :return: The filename of the output data file.
        :rtype: str
        '''
        input_filepath = os.path.join(self.output_folder, datafile)
        # 0. extract dimension of the silica slab
        silica_atoms = io.read(input_filepath, format="lammps-data", style="atomic")
        dimensions = silica_atoms.cell.cellpar()
        Lx, Ly, Lz = dimensions[:3] 
        
        # 1. add water slab on top of the silica system 
        # geometry = BoxGeometry((Lx/2, Ly/2, self.lz + Wz/2 + 6), (Lx, Ly, Wz))
        geometry = BoxGeometry((Lx/2, Ly/2, (Lz-self.vacuum) + water_thickness/2 + 6), (Lx, Ly, water_thickness))
        # water = pack_water(nummol=512, geometry=geometry)
        water_slab_filepath = os.path.join(self.output_folder,"water_slab.data")
        water_slab = pack_water(volume=Lx*Ly*water_thickness, geometry=geometry)
        water_slab.cell[0,0] = Lx
        water_slab.cell[1,1] = Ly
        water_slab.cell[2,2] = Lz
        water_slab.write(water_slab_filepath, format="lammps-data")
        
        # 2. add the water slab to the silica surface 
        silica_slab_filepath = os.path.join(self.output_folder,datafile)
        print("XXX",silica_slab_filepath )
        silica_atoms = io.read(silica_slab_filepath, format="lammps-data", style="atomic")
        atom_types = set(atom.symbol for atom in silica_atoms)  # Collect unique atom symbols
        num_atom_types = len(atom_types)  # Count the number of unique atom types
        ## why do i need to read it again, why not use water_slab??? 
        water_atoms = io.read(water_slab_filepath, format="lammps-data", style="atomic")
        if num_atom_types == 2: # to avoid conflict in atomic number when merging the two system
            silica_atoms.set_chemical_symbols(np.where(silica_atoms.get_atomic_numbers()==1, "O", "Si"))
            water_atoms.set_chemical_symbols(np.where(water_atoms.get_atomic_numbers()==1, "H", "O"))
            
        system = silica_atoms + water_atoms 
        water_silica_filepath = os.path.join(self.output_folder, "water_silica_system.data")
        system.write(water_silica_filepath, format="lammps-data")
        
        # 3. execute the passivation run
        output_filepath = os.path.join(self.output_folder, output_file)
        lmps_args = {
            'passivate_T' : passivation_time, 
            'input_filename' : water_silica_filepath,
            'potential_filename' : self.sio2_h2o_potential,
            'output_filename' : output_filepath,
            'water_zInitial' : (Lz-self.vacuum) + water_thickness/2 + 5,
            'water_zFinal' :  Lz, 
            'silica_zInitial' : 0, 
            'silica_zFinal' : Lz - self.vacuum,
            'si_slabDepth' : (Lz - self.vacuum)-1 
        }
        if run:
            script_path = os.path.join(script_dir, 'script', 'in.passivate') 
            self.execute_lammps(lmps_args, mpirun_n, lmp_exec, script=script_path)
        return output_file
    
    
    
    def anneal(self, datafile, anneal_time=None, mpirun_n=4, lmp_exec='lmp', run=True):
        '''Run a basic annealing procedure on the system.

        :param datafile: The filename of the input data file containing the system configuration.
        :type datafile: str
        :param anneal_time: The duration of the annealing procedure.
        :type anneal_time: float
        :param run: (Optional) Whether to run the annealing procedure. Default is True.
        :type run: bool
        :param mpirun_n: (Optional) The number of MPI processes to use for the annealing simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is 'lmp'.
        :type lmp_exec: str

        :return: The filename of the output data file containing the annealed system configuration.
        :rtype: str
        '''
        input_filepath = os.path.join(self.output_folder, datafile)
        output_filepath = os.path.join(self.output_folder, "annealed_silica.data")
        
        lmp_args = {
            'input_filename' : input_filepath, 
            'potential_filename' : self.sio2_potential, 
            'output_filename' : output_filepath,
            'anneal_T' : anneal_time,
            'seed' : 123 
        }
        if run:
            script_path = os.path.join(script_dir, 'script', 'in.anneal') 
            self.execute_lammps(lmp_args, mpirun_n, lmp_exec, script=script_path)
        return 'annealed_silica.data'
    
    
    def minimize_sio2(self, datafile, output_filename=None, mpirun_n=4, lmp_exec='lmp', run=True):
        '''Run a minimization procedure to optimize the structure of silica.

        :param datafile: The filename of the input data file containing the initial structure of silica.
        :type datafile: str
        :param mpirun_n: (Optional) The number of MPI processes to use for the minimization simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is 'lmp'.
        :type lmp_exec: str
        :param run: (Optional) Whether to run the minimization simulation. Default is True.
        :type run: bool

        :return: The filename of the output data file containing the minimized silica structure.
        :rtype: str
        '''
        if output_filename is None:
            output_filename = "minimized_silica.data"
        input_filepath = os.path.join(self.output_folder, datafile)
        output_filepath = os.path.join(self.output_folder, output_filename)
       
        lmp_args = {
            'input_filename' : input_filepath, 
            'potential_filename' : self.sio2_potential, 
            'output_filename' : output_filepath,
        }
        if run:
            script_path = os.path.join(script_dir, 'script', 'in.minimize_sio2') 
            self.execute_lammps(lmp_args, mpirun_n, lmp_exec, script=script_path)
        return output_filename
    

    def set_silanol(self, silanol_concentration, input_filename, output_filename=None, lmp_exec="lmp", delete_oxygen=True, run=True):    
        '''Set the concentration of silanol groups in the silica structure.

        :param silanol_concentration: The concentration of silanol groups to set in the structure.
        :type silanol_concentration: float
        :param input_filename: The filename of the input structure file.
        :type input_filename: str
        :param output_filename: The filename for the output structure file with modified silanol concentration.
        :type output_filename: str
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is "lmp".
        :type lmp_exec: str
        :param delete_oxygen: (Optional) Whether to delete oxygen atoms. Default is True.
        :type delete_oxygen: bool
        :param run: (Optional) Whether to run the simulation. Default is True.
        :type run: bool

        :return: The filename of the output structure file with modified silanol concentration.
        :rtype: str
        '''
        input_filepath = os.path.join(self.output_folder, datafile)
        output_filepath = os.path.join(self.output_folder, output_filename)
        if output_filepath is None:
            output_filepath = os.path.join(self.output_folder, f"silica_silanol={silanol_concentration}.data")
        
        lmps_args = {
            'silanol' : silanol_concentration, 
            'input_filename' : input_filepath, 
            'output_filename' : output_filepath,
            'potential_filename' : self.sio2_h2o_potential, 
        }
        if run:
            script_path = os.path.join(script_dir, 'script', 'in.set_silanol') 
            self.execute_lammps(lmps_args, lmp_exec=lmp_exec, mpirun_n=1, script=script_path)
        return output_filename
    
    
    
    
    def build_amorphous(self, mpirun_n=None, lmp_exec="lmp", cristobalite=None, seed=None, cluster=None, run=True, partition=None):
        '''Build an amorphous system using molecular dynamics simulation.
        This method constructs an amorphous system by annealing the beta cristobalite structure using molecular dynamics simulation.

        :param mpirun_n: (Optional) The number of MPI processes to use for the simulation.
        :type mpirun_n: int or None
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is "lmp".
        :type lmp_exec: str
        :param cristobalite: (Optional) The filename of the beta cristobalite structure data file. If not provided, it will be generated.
        :type cristobalite: str or None
        :param cluster: (Optional) Cluster information for parallel execution. Default is None.
        :type cluster: str or None
        :param run: (Optional) Whether to run the simulation. Default is True.
        :type run: bool
        :param partition: (Optional) The partition information for parallel execution. Default is None.
        :type partition: str or None

        :return: The filename of the output amorphous system data file.
        :rtype: str
        '''
        import random
        # Generate a random number between 100 and 9999 (inclusive)
        if seed is None:
            seed = random.randint(100, 9999)
        # Create the beta cristobalite
        if cristobalite is None:
            cristobalite_filepath = os.path.join(script_dir, 'data', 'beta_cristobalite_atomic.data') 
        else:
            cristobalite_filepath = os.path.join(self.output_folder, cristobalite)
            
        output_filepath = os.path.join(self.output_folder,"amorphous.data")
        lmps_args = {
            'lx' : self.lx, 
            'ly' : self.ly, 
            'lz' : self.lz, 
            'cristobalite_data' : cristobalite_filepath, 
            'potential_filename' : self.sio2_potential,
            'output_filename' : output_filepath,
            'seed' : seed
        }
        
        if run:
            # Save the seed for future reference
            seed_filepath = os.path.join(self.output_folder,"seed.txt")
            with open(seed_filepath,"a") as f: 
                f.write(str(seed)+"\n")
                
            script_path = os.path.join(script_dir, 'script', 'in.temp_quench') 
            self.execute_lammps(lmps_args, mpirun_n, lmp_exec,  script=script_path) 
        return f"amorphous.data"
    
    
        
    def execute_lammps(self, lmps_args, mpirun_n=4, lmp_exec='lmp', script=None, slurm=False):
        '''Execute a LAMMPS simulation with the given arguments.
        This method executes a LAMMPS simulation using the specified LAMMPS executable and input script, along with additional arguments provided in the lmps_args dictionary.

        :param lmps_args: A dictionary containing the LAMMPS variables and their corresponding values.
        :type lmps_args: dict
        :param mpirun_n: (Optional) The number of MPI processes to use for the simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The path to the LAMMPS executable. Default is 'lmp'.
        :type lmp_exec: str
        :param script: (Optional) The path to the input script for LAMMPS. If provided, this takes precedence over the slurm parameter. Default is None.
        :type script: str or None
        :param slurm: (Optional) If True, run the simulation using Slurm on a HPC cluster. Default is False.
        :type slurm: bool

        :return: None
        '''
        lmps_args_list = ' '.join(f'-var {key} {val}' for key, val in lmps_args.items())
        if slurm:
            computer = SlurmCPU(1, lmp_exec=lmp_exec, 
                                slurm_args={'job-name': f"SilicaSim",  "partition": partition})
            sim = Simulator(directory=data_dir, overwrite=True)
            sim.set_input_script(f"{data_dir}/in.lammps", copy=False, **lmps_args)
            sim.run(computer=computer)
        else:
            os.system(f"mpirun -n {mpirun_n} {lmp_exec} -in {script} {lmps_args_list}")
            
    
    def add_water(self, datafile,  dimensions, num_mol=None, output_filename="water-silica_system-new.data", thermalized_h2o=None, geometry="rectangle", mpirun_n=4, lmp_exec='lmp', run=True):
        '''Add a water slab on top of a silica system.
        This method adds a water slab on top of a silica system and runs a thermalization simulation if specified.

        :param datafile: The filename of the input data file containing the silica system.
        :type datafile: str
        :param dimensions: The dimensions of the silica system in the format (lx, ly, lz).
        :type dimensions: tuple
        :param num_molecules: (Optional) The number of water molecules to add. Default is None.
        :type num_molecules: int or None
        :param thermalized_water: (Optional) The temperature for thermalization of water molecules. Default is None.
        :type thermalized_water: float or None
        :param geometry: (Optional) The geometry of the water slab. Default is "rectangle".
        :type geometry: str
        :param mpirun_n: (Optional) The number of MPI processes to use for the simulation. Default is 4.
        :type mpirun_n: int
        :param lmp_exec: (Optional) The LAMMPS executable path or command. Default is "lmp".
        :type lmp_exec: str
        :param run: (Optional) Whether to run the simulation. Default is True.
        :type run: bool

        :return: The filename of the output data file containing the silica-water system.
        :rtype: str
        '''
        
        # Extract dimension of the silica slab
        lx, ly, lz = dimensions[:3] 
        input_filepath = os.path.join(self.output_folder, datafile)
        silica_atoms = io.read(input_filepath, format="lammps-data", style="atomic")
        si_dimensions = silica_atoms.get_cell_lengths_and_angles()
        si_Lx, si_Ly, si_Lz = si_dimensions[:3] 
        
        # Set default dimensions if not provided:
        # To create a rectangle that is periodic on one side, 1.0 is a buffer to avoid hydronium formation at the start
        lx = lx if lx is not None else si_Lx - 1.0
        ly = ly if ly is not None else si_Ly - 1.0
        lz = lz if lz is not None else si_Lz - 1.0
        
        # 1. add water slab on top of the silica system 
        center = (si_Lx/2, si_Ly/2, (si_Lz-self.vacuum)  + lz/2 + 9)
        dimension = (lx, ly, lz)
        geometry = BoxGeometry(center, dimension)
        
        if num_mol is None:
            water = pack_water(volume=lx*ly*lz, geometry=geometry)
        else:
            water = pack_water(nummol=num_mol, geometry=geometry)
        water.cell[0,0] = si_Lx
        water.cell[1,1] = si_Ly
        water.cell[2,2] = si_Lz
        
        water_slab_filepath = os.path.join(self.output_folder, "water_slab.data")
        silica_slab_filepath = os.path.join(self.output_folder, datafile)
        water.write(water_slab_filepath , format="lammps-data")
        # 2. add the water slab to the silica surface 
        silica_atoms = io.read(silica_slab_filepath, format="lammps-data", style="atomic")
        water_atoms = io.read(water_slab_filepath , format="lammps-data", style="atomic")
        
        # 2.1 run a thermalization on the water first before adding it
        if thermalized_h2o is not None:
            water_slab_thermalized_filepath = os.path.join(self.output_folder, "water_slab_thermalized.data")
            lmp_args = {
                'input_filename' : water_slab_filepath, 
                'potential_filename' : self.h2o_potential, 
                'output_filename' : water_slab_thermalized_filepath,
                'therm_T' : thermalized_h2o,
                'therm_temp' : 300, 
            }
            script_path = os.path.join(script_dir, 'script', 'in.thermalize_h2o') 
            self.execute_lammps(lmp_args, mpirun_n, lmp_exec, script=script_path)
            # read the thermalized water structure 
            water_atoms = io.read(water_slab_thermalized_filepath, format="lammps-data", style="atomic")
        system = silica_atoms + water_atoms 
        output_filepath = os.path.join(self.output_folder, output_filename)
        system.write(output_filepath, format="lammps-data")
        return output_filename
        
        
    
    def resize_z(self, datafile, z_value, output_filename=None, overwrite=False, run=True):
        '''Resize the z-dimension of the system in the input data file.
        This method resizes the z-dimension of the system in the input data file to the specified value.

        :param datafile: The filename of the input data file containing the system.
        :type datafile: str
        :param z_value: The new value for the z-dimension of the system.
        :type z_value: float
        :param output_filename: (Optional) The filename for the output data file with the resized system. If not provided, the input data file is overwritten. Default is None.
        :type output_filename: str or None
        :param overwrite: (Optional) Whether to overwrite the input data file with the resized system. If False, a new output file with the resized system is created. Default is False.
        :type overwrite: bool
        :param run: (Optional) Whether to perform the resizing operation. Default is True.
        :type run: bool

        :return: The filename of the data file containing the resized system.
        :rtype: str
        '''
        input_filepath = os.path.join(self.output_folder,datafile)
        output_filepath = os.path.join(self.output_folder,output_filename)
        silica_atoms = io.read(input_filepath, format="lammps-data", style="atomic")
        silica_atoms.cell[2,2] = z_value
        if run:
            if overwrite:
                silica_atoms.write(input_filepath, format="lammps-data")
                return datafile
            else:
                silica_atoms.write(output_filepath, format="lammps-data")
                return output_filename
            
    
    def add_vacuum(self, datafile, vacuum, output_filename=None, overwrite=False, run=True):
        '''Add vacuum space to the z-dimension of the system in the input data file.
        This method adds vacuum space to the z-dimension of the system in the input data file, effectively increasing its size along the z-axis.

        :param datafile: The filename of the input data file containing the system.
        :type datafile: str
        :param vacuum: The amount of vacuum space to add to the z-dimension of the system.
        :type vacuum: float
        :param output_filename: (Optional) The filename for the output data file with the added vacuum space. If not provided, the input data file is overwritten. Default is None.
        :type output_filename: str or None
        :param overwrite: (Optional) Whether to overwrite the input data file with the system containing added vacuum space. If False, a new output file with the added vacuum space is created. Default is False.
        :type overwrite: bool
        :param run: (Optional) Whether to perform the operation of adding vacuum space. Default is True.
        :type run: bool

        :return: The filename of the data file containing the system with added vacuum space.
        :rtype: str
        '''
        input_filepath = os.path.join(self.output_folder,datafile)
        output_filepath = os.path.join(self.output_folder,output_filename)
        silica_atoms = io.read(input_filepath, format="lammps-data", style="atomic")
        dimensions = silica_atoms.cell.cellpar() # get dimensions of the system
        _, _, Lz = dimensions[:3] 
        self.lz_new = Lz # used in passivation | it assumes that add_vacuum is going to be used!
        silica_atoms.cell[2,2] = Lz + vacuum
        if run:
            if overwrite:
                silica_atoms.write(input_filepath, format="lammps-data")
                return datafile
            else:
                silica_atoms.write(output_filepath, format="lammps-data")
                return output_filename
    
    
    def get_silanol(self, datafile):
        '''Calculate the density of silanol groups in the silica system.

        This method calculates the density of silanol groups in the silica system described by the provided datafile. 
        Silanol groups are hydroxyl (-OH) functional groups attached to the silica surface.

        :param datafile: The filename of the input data file containing the silica system.
        :type datafile: str

        :return: The density of silanol groups in the silica system, measured in silanol groups per square nanometer (nm^-2).
        :rtype: float
        '''
        input_filepath = os.path.join(self.output_folder,datafile)
        silica_atoms = io.read(input_filepath, format="lammps-data", style="atomic")
        dimensions = silica_atoms.cell.cellpar()
        # Get surface area of the system in nm
        Lx, Ly, _ = dimensions[:3]
        # Count the number of hydrogen atoms in the system
        num_hydrogen = sum(1 for atom in silica_atoms if atom.symbol == 'H')
        return (num_hydrogen)/(Lx*Ly/100)