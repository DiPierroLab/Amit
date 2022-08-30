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

    def __init__(self, traj, top, snap_index, chip_seq):
        
        # read Chi-seq patterns
        self.snap_index   = snap_index
        self.chip_seq     = chip_seq
        self.chrom_beadsN = md.load(top).n_atoms
            
        # Must provide Chip-seq binding affinities for all beads
        assert len(self.chip_seq) == self.chrom_beadsN
        
        # Load mtraj and pass it to mbuild for adding TF particles
        self.chromosome = mb.load ( md.load_frame(traj, snap_index, top) )
        
        count, bins = np.histogram(self.chip_seq, bins=100)   
        
        bin_centers = (bins[:-1] + bins[1:])/2
        
        for i in range(self.chromosome.n_particles): 
           
                # Name chromosome loci as X + unqiue bin_value
                self.chromosome[i].name = 'X'+ str( find_nearest(bin_centers, self.chip_seq[i]) )
            

        # Create a TF particle
        self.tf_bead  = TFparticle()
        
    def __repr__(self):
        
        return f'Nucleus with {self.chrom_beadsN} number of beads extracted from snapshot {self.snap_index}.\n Chip-seq average binding affinity is {np.mean(self.chip_seq)}.\n Now you may use sprinkle_tf() method to add TFs to your taste.'
        

    def visualize(self):
        
        '''Use mbuild visualize functionality'''
        
        self.chromosome.visualize({'TF': 'blue'})
        
        
    def sprinkle_tf(self, nTF, boxL,  overlap=0.9, edge=0.1, seed=1235, inpfile = 'input.gsd'):
        
        '''Use mbuild packmol functionality to sprinkle tf around chromsome and save input file for simulation'''
        
        
        self.chromosome_tf =        mb.packing.fill_box(compound = [self.chromosome, self.tf_bead], n_compounds = [1,  nTF], 
                             box=(boxL,boxL,boxL), density=None, 
                             overlap=0.9, seed=12345, 
                             edge=0.1, compound_ratio=None, 
                             aspect_ratio=None, fix_orientation=False, 
                             temp_file=None, update_port_locations=False)
        
        self.chromosome_tf.save(inpfile)
        
        
    def visualize_tf(self):
        
        try: 
            self.chromosome_tf.visualize({'TF': 'blue'})
        except: 
            print("Oops! Need to add TF before visualizing")
        
        self.chromosome_tf.visualize({'TF': 'blue'})
        
        
    def sim_param(self, sigmaX=1, sigmaTF=1, epsTF = 1, r_c = 3.5, inpfile = 'input.gsd'):
        '''Provide simulation parameters'''
        
        hoomd.context.initialize("")
        system = hoomd.init.read_gsd(inpfile)
        
        nl    = hoomd.md.nlist.cell()
        lj    = hoomd.md.pair.lj(r_c, nlist=nl)
        
        print('done')
        
        for i in range(100): 
            
            # TF-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'TF', epsilon= epsTF * self.chip_seq[i], 
                              sigma = 0.5 * (sigmaTF + sigmaX), r_cut = r_c )

            # chromosome-chromsome interactions
            #lj.pair_coeff.set('X'+str(i), 'X'+str(i), epsilon=0, sigma = 0, r_cut = 0 )

            for j in range(100):

                lj.pair_coeff.set('X'+str(i), 'X'+str(j), epsilon= 0, sigma = 0, r_cut = 0 )
            
        # TF-TF interactions?
        lj.pair_coeff.set('TF', 'TF', epsilon=epsTF, sigma = sigmaTF, r_cut = r_c )
            
        
    def simulate(self, gam=0.01, kT=1, dt = 0.01, tsim=10e4, outfile='sim_out', freqout=1000, seed =12373):
        
        tf_group = hoomd.group.type('TF')
        hoomd.md.integrate.mode_standard(dt=dt)
        
        integrator = hoomd.md.integrate.langevin(group = tf_group, kT = kT, seed=seed)
        integrator.set_gamma('TF', gamma=gam)
        
        
        #Logging sim data
        hoomd.analyze.log(filename = outfile + '.log', quantities= ['pair_lj_energy','pair_yukawa_energy','potential_energy','kinetic_energy','temperature'], 
                          period=freqout, overwrite=True, header_prefix='#')
                                                         
        # Dump gsd traj
        hoomd.dump.gsd(outfile + '.gsd',  period=freqout, group=tf_group, overwrite=True)
        
  


        #Run MD
        hoomd.run(tsim) 
