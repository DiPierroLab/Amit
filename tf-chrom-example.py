#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import hoomd, hoomd.md
import mdtraj as md
import mbuild as mb
import numpy  as np

import warnings
warnings.filterwarnings('ignore')


# In[64]:


get_ipython().system(' rm *gsd')


# In[60]:


from src.main_without_top import ChromoChip
#help(ChromoChip)


# In[67]:


nucleus = ChromoChip('./data/GM12878-OliveiraJunior2020_chr21.pdb.pdb', 5, 100, chip_seq = np.random.rand(964)*10)


# In[5]:


nucleus


# In[6]:


nucleus.chromosome.visualize()


# In[68]:


nucleus.sprinkle_tf(nTF=1000, boxL=2, seed=1235 )


# In[19]:


nucleus.chromosome_tf


# In[75]:


nucleus.chromosome_tf.visualize(color_scheme={'TF':'green'})


# In[70]:


nucleus.sim_param(sigmaX=1, sigmaTF=1, epsTF = 1, r_c = 3.5)


# In[71]:


nucleus.simulate(gam=1, kT=1, dt = 0.01, tsim=10e4, outfile='sim_out', freqout=1000, seed =12383)


# In[81]:


xx = mb.load(md.load('sim_out.gsd'))


# In[93]:


nucleus.chromosome_tf.visualize()


# In[92]:


xx.visualize(color_scheme={'TF':'green'})


# In[39]:


yy.visualize(color_scheme={'TF':'green'})

