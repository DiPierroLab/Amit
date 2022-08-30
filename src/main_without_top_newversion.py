import hoomd, hoomd.md
import mdtraj  as md
import mbuild  as mb
import numpy   as np


def find_nearest(array, value):
    '''Util function to bin values '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class TFparticle(mb.Compound):
    
    '''Simple particle class to add to chromsome configuration via packing/solvation'''
    
    def __init__(self):
        
        super(TFparticle, self).__init__()
        
        self.add( mb.Particle(name='TF', pos=[0, 0, 0]) )
        

class ChromoChip(object):

    '''A simple polymer class to load chromosome structure and assign bead types based on Chip-seq dictionary
    INPUT: traj and topology input for mdtraj (gro and dcd/h5) AND snapshot index AND chip_seq affinity list of energy multiples of kT, e.g [1.2, 1.4, 1.5, 100, ...]
        '''

    def __init__(self, traj, snap_index, num_bins, chip_seq):
        
        # read Chi-seq patterns
        self.snap_index   = snap_index
        self.chip_seq     = chip_seq
        self.chrom_beadtypesN      = num_bins
        #self.chrom_beadsN = md.load(top).n_atoms
                    
        # Load mtraj and pass it to mbuild for adding TF particles
        
        newtraj = md.load(traj);
        self.chromosome = mb.load ( newtraj[snap_index] )
        
        count, bins = np.histogram(self.chip_seq, num_bins)   
        
        bin_centers = (bins[:-1] + bins[1:])/2
        
        for i in range(self.chromosome.n_particles): 
           
                # Name chromosome loci as X + unqiue bin_value
                self.chromosome[i].name = 'X'+ str( find_nearest(bin_centers, self.chip_seq[i]) )
            
        self.coarse_chip_seq = bin_centers            

        # Must provide Chip-seq binding affinities for all beads
        #assert len(self.chip_seq) == self.chromosome.n_particles

        # Create a TF particle
#         self.tf_bead  = TFparticle()
        
    def __repr__(self):
        
        return f'Nucleus with {self.chromosome.n_particles} number of beads extracted from snapshot {self.snap_index}.\n Chip-seq average binding affinity is {np.mean(self.chip_seq)}.\n Now you may use sprinkle_tf() method to add TFs to your taste.'
        

    def visualize(self):
        
        '''Use mbuild visualize functionality'''
        
#         for i in range(self.chromosome.n_particles):            
#                 self.chromosome.visualize({: 'green'} )
                
        self.chromosome.visualize(color_scheme = {'TF': 'red'})
        
        
    def sprinkle_tf(self, nTF, boxL,  overlap=0.2, edge=0.1, seed=1235, inpfile = 'input_no_top.gsd'):
        
        '''Use mbuild packmol functionality to sprinkle tf around chromsome and save input file for simulation'''
        
        
        self.chromosome_tf = mb.packing.solvate(self.chromosome, self.tf_bead, nTF, 
                             box=(boxL,boxL,boxL),  
                             overlap=0.2, seed=12345, 
                             edge=0.1, fix_orientation=False, 
                             temp_file=None, update_port_locations=False)
        
       
        self.chromosome_tf_pbc = mb.packing.fill_box(compound = self.chromosome_tf, n_compounds = 1, 
                             box=(-boxL,-boxL,-boxL,boxL,boxL,boxL), density=None, 
                             overlap=0.1, seed=12345, 
                             edge=0.2, compound_ratio=None, 
                             aspect_ratio=None, fix_orientation=False, 
                             temp_file=None, update_port_locations=False)
        self.chromosome_tf_pbc.save(inpfile)
        
    def make_nuclear_shell(self, nTF, chr_shell_radius, boxL, origin = [0,0,0], overlap=0.2, edge=0.1, seed=1235, inpfile = 'input_no_top.gsd'):
        
#         '''Use mbuild packmol functionality to sprinkle tf outside a spherical shell around the chromsome and save input file for simulation'''
                
        tf_bead  = TFparticle()

        class tf_shell(mb.Compound):
            def __init__(self):
                super(tf_shell, self).__init__()
        #         lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

                n = 0
                while n < nTF:
                    random_point = mb.Random3DPattern(1)
                    random_point.scale(chr_shell_radius)
                    pos = random_point.points[0]
                    pos = [2*x - chr_shell_radius for x in pos]
                    dd = np.linalg.norm(pos)
                    if dd > chr_shell_radius:
                        lj_particle = mb.clone(tf_bead)
                        pos[0]-=1.0
                        lj_particle.translate(pos)
                        self.add(lj_particle)
                        n += 1

        self.chrom_tf = tf_shell()
        self.chrom_tf.translate_to(origin)

        chr_mb = self.chromosome
        chr_mb.translate_to(origin)

        self.chrom_tf.add(chr_mb)
        self.chromosome_tf_pbc = mb.packing.fill_box(compound = self.chrom_tf, n_compounds = 1, 
                             box=(-boxL,-boxL,-boxL,boxL,boxL,boxL), density=None, 
                             overlap=0.1, seed=12345, 
                             edge=0.2, compound_ratio=None, 
                             aspect_ratio=None, fix_orientation=False, 
                             temp_file=None, update_port_locations=False)

        
        self.chromosome_tf_pbc.save(inpfile)
        
    
#     def visualize_shell(self):
        
#         try: 
#             self.chromosome_shell_tf.visualize(color_scheme = {'TF': 'green'})
#         except: 
#             print("Oops! Need to add TF before visualizing")
        
#         self.chromosome_shell_tf.visualize(color_scheme = {'TF': 'green'})
        
#     def visualize_tf(self):
        
#         try: 
#             self.chromosome_tf.visualize(color_scheme = {'TF': 'blue'})
#         except: 
#             print("Oops! Need to add TF before visualizing")
        
#         self.chromosome_tf.visualize(color_scheme = {'TF': 'blue'})
        
        
    def sim_param(self, sigmaX=1, sigmaTF=1, epsTF = 1, rcut_tf = 2.5, rcut_chtf = 7.5, inpfile = 'input_no_top.gsd'):
        '''Provide simulation parameters'''
        
        hoomd.context.initialize("")
        system = hoomd.init.read_gsd(inpfile)
        
        nl    = hoomd.md.nlist.cell()
        lj    = hoomd.md.pair.lj(rcut_tf, nlist=nl)
        morse = hoomd.md.pair.morse(rcut_chtf, nlist=nl)
        
        print('done')
        
        for i in range(self.chrom_beadtypesN): 
            
            # TF-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'TF', epsilon= epsTF, sigma = 0.5 * (sigmaTF + sigmaX), r_cut = rcut_chtf )
            morse.pair_coeff.set('X'+str(i), 'TF', D0 = self.coarse_chip_seq[i], alpha = 3.0, 
                              r0 = 0.5 * (sigmaTF + sigmaX), r_cut = rcut_chtf )

            for j in range(self.chrom_beadtypesN):

                lj.pair_coeff.set('X'+str(i), 'X'+str(j), epsilon= 0, sigma = 0, r_cut = 0 )
                morse.pair_coeff.set('X'+str(i), 'X'+str(j), D0 = 0, alpha = 0.0, r0 = 0, r_cut = 0 )
            
        # TF-TF interactions?
        lj.pair_coeff.set('TF', 'TF', epsilon=epsTF, sigma = sigmaTF, r_cut = rcut_tf )        
        morse.pair_coeff.set('TF', 'TF', D0 = 0, alpha = 0.0, r0 = 0, r_cut = 0 )
           
        
    def simulate(self, gam=0.01, kT=1, dt = 0.01, tsim=10e4, outfile='sim_out', freqout=1000, seed =12373):
        
        tf_group = hoomd.group.type('TF')
        hoomd.md.integrate.mode_standard(dt=dt)
        
        integrator = hoomd.md.integrate.langevin(group = tf_group, kT = kT, seed=seed)
        integrator.set_gamma('TF', gamma=gam)
        
        
        #Logging sim data
        hoomd.analyze.log(filename = outfile + '.log', quantities= ['pair_lj_energy','pair_morse_energy','potential_energy','kinetic_energy','temperature'], 
                          period=freqout, overwrite=True, header_prefix='#')
                                                         
        # Dump gsd traj
        hoomd.dump.gsd(outfile + '.gsd',  period=freqout, group=tf_group, overwrite=True)
        
  


        #Run MD
        hoomd.run(tsim) 
