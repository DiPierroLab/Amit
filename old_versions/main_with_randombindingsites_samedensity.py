import hoomd, hoomd.md
import mdtraj  as md
import mbuild  as mb
import numpy   as np

def find_nearest(array, value):
    '''Util function to bin values '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def pad_zero_to_fit(a):
    max_len = max(len(line) for line in a)
    my_array = np.zeros([len(a), max_len])
    for i, line in enumerate(a):
        my_array[i, :len(line)] = line

    return my_array 

def pair_sigmoid(r,rmin,rmax,epsilon,a,s,sigma):
    import numpy as np
    V =  -0.5*epsilon*(1.0-np.tanh((r/sigma-a)/s))
    F = 0.5*epsilon/s*(1.0-np.tanh((r/sigma-a)/s)**2.0)
    return (V, F)

class TFparticle(mb.Compound):
    
    '''Simple particle class to add to chromsome configuration via packing/solvation'''
    
    def __init__(self):
        
        super(TFparticle, self).__init__()
        
        self.add( mb.Particle(name='TF', pos=[0, 0, 0]) )
                
class ChromoChip5(object):

    '''A simple polymer class to load chromosome structure and assign bead types based on Chip-seq dictionary
    INPUT: traj and topology input for mdtraj (gro and dcd/h5) AND snapshot index AND chip_seq affinity list of energy multiples of kT, e.g [1.2, 1.4, 1.5, 100, ...]. This particular class also organizes the binding sites randomly in space keeping the packing density of the polymer same.
        '''

    def __init__(self, traj, snap_index, energy_bin_width, chipfile, chromatin_bead_size):
        
        # process chip_seq data
        datafile = open(chipfile,'r')
        lines = datafile.readlines()
        lines2 = [line.strip().split(',') for line in lines]
        chip_seq = pad_zero_to_fit(lines2)        
        
        # read Chromosome configuration: load mtraj and pass it to mbuild for adding TF particles                

        self.snap_index = snap_index
        newtraj = md.load(traj);
        self.chromosome = mb.load ( newtraj[snap_index] )

        # Chip-seq patterns                

        self.all_chip_seq = chip_seq
                                    
        # process chip_seq input to generate a coarse-grained list of all binding energies
        
        all_be = self.all_chip_seq[:,1:]  # collecting all binding energies
        input_bins = np.arange(0, all_be.max()+energy_bin_width, energy_bin_width)
        count, bins = np.histogram(all_be.flat, bins = input_bins)           
        bin_centers = (bins[:-1] + bins[1:])/2        

        self.chrom_beadtypesN = len(bin_centers)  # number of particle types = # of coarse-grainining bins for binding energies           
        self.coarse_chip_seq = bin_centers 
        
        # generate chromosome configuration with the binding sites for each bead
        
        self.chr_binding_sites = mb.Compound() 
        CBS_proto = mb.Particle(name='CBS', pos=[0, 0, 0]) # CBS - chromatin binding site
        print(chromatin_bead_size)
        ll=0
        
        xl = self.chromosome.xyz.max(0)-self.chromosome.xyz.min(0)
        dim = np.min(xl)
        print(xl,dim)
        for i in range(self.chromosome.n_particles):
            xpos = mb.Random3DPattern(1,i)     #self.chromosome.xyz[i,:]
            xpos.scale(dim)
            cpos = xpos.points[0,:]
#             print(cpos)
            xa = all_be[i,:]
            ii = np.where(xa > 0)
            nn = len(ii[0])
            chr_bead_sites = mb.Compound()
            if nn > 1:
                ll += 1    
                pattern_sphere = mb.SpherePattern(nn)
                pattern_sphere.scale(chromatin_bead_size)
                m = 0
                for pos in pattern_sphere:
                    CBS_particle = mb.clone(CBS_proto)
                    CBS_particle.name = 'X'+ str( find_nearest(self.coarse_chip_seq, xa[m]) )
                    CBS_particle.translate(pos)        
                    chr_bead_sites.add(CBS_particle)
                    m += 1
                chr_bead_sites.translate(cpos)
            else:
                CBS_particle = mb.clone(CBS_proto)
                CBS_particle.name = 'X'+ str( find_nearest(self.coarse_chip_seq, xa[0]) )
                CBS_particle.translate(cpos)        
                chr_bead_sites.add(CBS_particle)

            self.chr_binding_sites.add(chr_bead_sites)
#         chr_binding_sites.translate_to(origin)  


    def __repr__(self):
        
        return f'Nucleus with {self.chromosome.n_particles} number of beads extracted from snapshot {self.snap_index}.\n Chip-seq average binding affinity is {np.mean(self.all_chip_seq[:,1:].flat)}.\n Now you may use sprinkle_tf() method to add TFs to your taste.'
        

    def visualize(self):
        
        '''Use mbuild visualize functionality'''
        
#         for i in range(self.chromosome.n_particles):            
#                 self.chromosome.visualize({: 'green'} )
                
        self.chromosome.visualize(color_scheme = {'TF': 'red'})                     
        
    def make_nuclear_shell_with_binding_sites(self, gridsize, TF_diameter, chr_shell_radius, box_length, origin = [0,0,0], seed=1235, inpfile = 'input_binding_sites.gsd'):
        
#         '''Use mbuild packmol functionality to sprinkle tf outside a spherical shell around the chromsome and save input file for simulation'''
                
        # to convert all positions from MiChroM (in microns) in units of TF_diameter 
        
        chr_mb = self.chr_binding_sites
        chr_mb.xyz = chr_mb.xyz/TF_diameter                                                          
        r_nuc = chr_shell_radius/TF_diameter
        boxl = box_length/TF_diameter
        
        # for i in range(chr_binding_sites.n_particles):
        #     print(chr_binding_sites[i].name)

        tf_bead  = TFparticle()
        self.Chrom_TF = mb.Compound()
        pattern = mb.Random3DPattern(gridsize,seed)  # A random arrangement of pieces inside a 3D region.
        pattern.scale(boxl)
        n = 0
        for pos in pattern:
            pos = pos - boxl/2
            dd = np.linalg.norm(pos)
            if dd > r_nuc and abs(pos[0]) < boxl/2-0.5 and abs(pos[1]) < boxl/2-0.5 and abs(pos[2]) < boxl/2-0.5:
                TF_particle = mb.clone(tf_bead)
                TF_particle.translate(pos)        
                self.Chrom_TF.add(TF_particle)
                n += 1
#         print(n)
        
        chr_mb.translate_to(origin)  # place chromatin at the center
        self.Chrom_TF.add(chr_mb)
        self.Chrom_TF.periodicity = [boxl,boxl,boxl]
        self.Chrom_TF.save(inpfile)   

    def make_uniform_tf_with_binding_sites(self, nTF, TF_diameter, box_length, seed=1235, inpfile = 'input_binding_sites.gsd'):
         
        chr_mb = self.chr_binding_sites
        chr_mb.xyz = chr_mb.xyz/TF_diameter                                                          
        boxl2 = 0.5*box_length/TF_diameter
        
        tf_bead  = TFparticle()        
        self.Chrom_TF_uni = mb.packing.fill_box(compound = [chr_mb, tf_bead], n_compounds = [1,nTF], 
                             box=(-boxl2,-boxl2,-boxl2,boxl2,boxl2,boxl2), 
                             overlap=0.2, seed=seed, sidemax=boxl2-TF_diameter*0.1,
                             edge=0.1, compound_ratio=None, 
                             aspect_ratio=None, fix_orientation=True, 
                             temp_file=None, update_port_locations=False)
        self.Chrom_TF_uni.save(inpfile)
        
    def sim_param(self, sigmaX=1, sigmaTF=1, epsTF = 1, rcut_tf = 2.5, rcut_chtf = 2.5, ch_tf_alpha = 3, inpfile = 'input_no_top.gsd'):
        '''Provide simulation parameters'''
        
        hoomd.context.initialize("")
        system = hoomd.init.read_gsd(inpfile)

        snapshot = system.take_snapshot()
        x = snapshot.particles.position
        nn = x.shape[0]
        n1 = self.chr_binding_sites.n_particles
        xx = x[0:n1-1,:]
        print(snapshot.box.Lx,xx.max(0)-xx.min(0))
        yy = x[n1:,:]
        print(snapshot.box.Ly,yy.max(0)-yy.min(0))
        
        nl    = hoomd.md.nlist.cell()
        lj    = hoomd.md.pair.lj(rcut_tf, nlist=nl)
        morse = hoomd.md.pair.morse(rcut_chtf, nlist=nl)
        
        print('done')
        
        for i in range(self.chrom_beadtypesN): 
            
            # TF-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'TF', epsilon= epsTF, sigma = 0.5 * (sigmaTF + sigmaX), r_cut = rcut_chtf )
            morse.pair_coeff.set('X'+str(i), 'TF', D0 = self.coarse_chip_seq[i], alpha = ch_tf_alpha, 
                              r0 = 0.5*(sigmaTF + sigmaX), r_cut = rcut_chtf )

            for j in range(self.chrom_beadtypesN):

                lj.pair_coeff.set('X'+str(i), 'X'+str(j), epsilon= 0, sigma = 0, r_cut = 0 )
                morse.pair_coeff.set('X'+str(i), 'X'+str(j), D0 = 0, alpha = 0.0, r0 = 0, r_cut = 0 )
            
        # TF-TF interactions?
        lj.pair_coeff.set('TF', 'TF', epsilon=epsTF, sigma = sigmaTF, r_cut = rcut_tf )        
        morse.pair_coeff.set('TF', 'TF', D0 = 0, alpha = 0.0, r0 = 0, r_cut = 0 )
           
        
    def simulate(self, gam=0.01, kT=1, dt = 0.01, tsim=1e5, outfile='sim_out', freqout=1000, seed =12345):
        
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
        hoomd.run_upto(tsim) 

        
    def sim_param_full(self, sigmaX=1, sigmaTF=1, epsTF = 1, rcut_tf = 2.5, rcut_chtf = 2.5, ch_tf_alpha = 3, inpfile = 'input_no_top.gsd', restartfile='restart_no_top.gsd'):
        '''Provide simulation parameters'''
        
        hoomd.context.initialize("")
        system = hoomd.init.read_gsd(filename=inpfile, restart=restartfile)
       
        nl    = hoomd.md.nlist.cell()
        lj    = hoomd.md.pair.lj(rcut_tf, nlist=nl)
        morse = hoomd.md.pair.morse(rcut_chtf, nlist=nl)
        
        print('done')
        
        for i in range(self.chrom_beadtypesN): 
            
            # TF-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'TF', epsilon= epsTF, sigma = 0.5 * (sigmaTF + sigmaX), r_cut = rcut_chtf )
            morse.pair_coeff.set('X'+str(i), 'TF', D0 = self.coarse_chip_seq[i], alpha = ch_tf_alpha, 
                              r0 = 0.5*(sigmaTF + sigmaX), r_cut = rcut_chtf )

            for j in range(self.chrom_beadtypesN):

                lj.pair_coeff.set('X'+str(i), 'X'+str(j), epsilon= 0, sigma = 0, r_cut = 0 )
                morse.pair_coeff.set('X'+str(i), 'X'+str(j), D0 = 0, alpha = 0.0, r0 = 0, r_cut = 0 )
            
        # TF-TF interactions?
        lj.pair_coeff.set('TF', 'TF', epsilon=epsTF, sigma = sigmaTF, r_cut = rcut_tf )        
        morse.pair_coeff.set('TF', 'TF', D0 = 0, alpha = 0.0, r0 = 0, r_cut = 0 )
        
        
    def simulate_full(self, gam=0.01, kT=1, dt = 0.01, teq=1e5, tsim=1e5, outfile='sim_out', freqout=1000, restart_freq=10000, seed =12373):
        
        tf_group = hoomd.group.type('TF')
        hoomd.md.integrate.mode_standard(dt=dt)
        
        integrator = hoomd.md.integrate.langevin(group = tf_group, kT = kT, seed=seed)
        integrator.set_gamma('TF', gamma=gam)
        
        
        #Logging sim data
        hoomd.analyze.log(filename = outfile + '.log', quantities= ['pair_lj_energy','pair_morse_energy','potential_energy','kinetic_energy','temperature'], 
                          period=freqout, header_prefix='#', phase=0)

        #Run MD for equilibration, before dumping
        hoomd.run(teq)
        
        # Dump gsd traj
        hoomd.dump.gsd(outfile + '.gsd',  period=freqout, group=tf_group, phase=0)         
        gsd_restart = hoomd.dump.gsd(filename=outfile + '_restart.gsd', group=tf_group, truncate=True, period=restart_freq, phase=0)

        #Run MD
#         hoomd.run_upto(teq+tsim) 
        
        try:
            hoomd.run_upto(teq+tsim, limit_multiple=10000)
        except WalltimeLimitReached:
            pass

        gsd_restart.write_restart()
        
    def sim_param_pair_sigmoid(self, sigmaX=1, sigmaTF=1, epsTF = 1, rcut_tf = 2.5, rcut_chtf = 2.5, ch_tf_range = 3, s = 1, inpfile = 'input_no_top.gsd'):
        '''Provide simulation parameters'''
        
        hoomd.context.initialize("")
        system = hoomd.init.read_gsd(inpfile)

#         snapshot = system.take_snapshot()
#         x = snapshot.particles.position
#         nn = x.shape[0]
#         n1 = self.chr_binding_sites.n_particles
#         xx = x[0:n1-1,:]
#         print(snapshot.box.Lx,xx.max(0)-xx.min(0))
#         yy = x[n1:,:]
#         print(snapshot.box.Ly,yy.max(0)-yy.min(0))
        
        nl    = hoomd.md.nlist.cell()
        lj    = hoomd.md.pair.lj(rcut_tf, nlist=nl)
        table = hoomd.md.pair.table(width=1000, nlist=nl)
        
        print('done')
        
        for i in range(self.chrom_beadtypesN): 
            
            # TF-chromsome interactions
            lj.pair_coeff.set('X'+str(i), 'TF', epsilon= epsTF, sigma = 0.5 * (sigmaTF + sigmaX), r_cut = rcut_chtf )
            table.pair_coeff.set('X'+str(i), 'TF', func=pair_sigmoid, rmin=0.8, rmax=rcut_chtf,
                    coeff=dict(epsilon=self.coarse_chip_seq[i], a=ch_tf_range, s=1, sigma=0.5 * (sigmaTF + sigmaX)))

            for j in range(self.chrom_beadtypesN):

                lj.pair_coeff.set('X'+str(i), 'X'+str(j), epsilon= 0, sigma = 0, r_cut = 0 )
                table.pair_coeff.set('X'+str(i), 'X'+str(j), func=pair_sigmoid, rmin=0.8, rmax=5.0,
                                     coeff=dict(epsilon=0, a=3, s=1, sigma=1.0))
            
        # TF-TF interactions?
        lj.pair_coeff.set('TF', 'TF', epsilon=epsTF, sigma = sigmaTF, r_cut = rcut_tf )        
        table.pair_coeff.set('TF', 'TF', func=pair_sigmoid, rmin=0.8, rmax=5.0,
                             coeff=dict(epsilon=0, a=3, s=1, sigma=sigmaTF))

           
        
    def simulate_pair_sigmoid(self, gam=0.01, kT=1, dt = 0.01, tsim=1e5, outfile='sim_out', freqout=1000, seed =12345):
        
        tf_group = hoomd.group.type('TF')
        hoomd.md.integrate.mode_standard(dt=dt)
        
        integrator = hoomd.md.integrate.langevin(group = tf_group, kT = kT, seed=seed)
        integrator.set_gamma('TF', gamma=gam)
        
        
        #Logging sim data
        hoomd.analyze.log(filename = outfile + '.log', quantities= ['pair_lj_energy','pair_table_energy','potential_energy','kinetic_energy','temperature'], 
                          period=freqout, overwrite=True, header_prefix='#')
        
        # Dump gsd traj
        hoomd.dump.gsd(outfile + '.gsd',  period=freqout, group=tf_group, overwrite=True)         

        #Run MD
        hoomd.run_upto(tsim) 