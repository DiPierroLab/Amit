import numpy as np
import hoomd
import mdtraj as md
import mbuild  as mb



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
        
        # Load mtraj and pass it to mbuild
        self.chromosome = mb.load ( md.load_frame(traj, snap_index, top) )
        
        # Name chromosome loci as X0, X1, X2, ... etc
        for i in range(self.chromosome.n_particles): 
            self.chromosome.name = 'X'+str(i) 

        # TF particle
        self.tf_bead  = TFparticle()
        

    def visualize(self):
        
        '''Use mbuild visualize functionality'''
        
        self.chromosome.visualize()
        
        
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
        
        try self.packed.visualize():
        except ValueError:
            print("Oops! Need to add TF before visualizing")
        
        
    def sim_param(self, sigmaX=1, sigmaTF=1, epsTF = 1, r_c = 3.5):
        
        '''Provide simulation parameters'''
        
        nlist = hoomd.md.nlist.cell()
        lj    = hoomd.md.pair.lj(r_c, nlist)
        
        for i in range(len(chip_seq)): 
            
            # TF-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'TF', epsilon= epsTF * chip_seq[i], sigma = 0.5 * (sigmaTF + sigmaX), r_cut = r_c )
            
            # TF-TF interactions
            lj.pair_coeff.set('TF', 'TF', epsilon=epsTF, sigma = sigmaTF, r_cut = r_c )

            # chromosome-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'X'+str(i), epsilon=0, sigma = 0, r_cut = 0 )
            
        
    def simulate(mass = 10, gam = 0.01, kT=1, dt = 0.01, tsim=10e4, outfile='sim_out', freqout=1000):
        
        tf_group = groupA = group.type('TF')
        
        integrator.set_gamma('TF', gamma = 0.01*mass)
        
        
        #Logging sim data
        hoomd.analyze.log(filename = outfile + '.log', quantities=['pair_lj_energy','pair_yukawa_energy','potential_energy','kinetic_energy','temperature'], 
                          period=freqout, overwrite=True, header_prefix='#')
                                                         
        # Dump gsd traj
        hoomd.dump.gsd(outfile + '.gsd',  period=freqout, group=tf_group, overwrite=True)


        #Run MD
        hoomd.run(tsim)

        
        