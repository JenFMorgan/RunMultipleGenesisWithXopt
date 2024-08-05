from xopt import Xopt
import os
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import json
from impact import Impact
from distgen import Generator
from matplotlib import pyplot as plt
from genesis.version4 import Genesis4
from genesis import tools
import shutil

def archive_infolder(self, h5=None, archive_file=None):
        """
        Dump inputs and outputs into HDF5 file.

        Parameters
        ----------
        h5 : str or h5py.File
            The filename or handle to HDF5 file in which to write the information.
            If not in informed, a new file is generated.

        Returns
        -------
        h5 : h5py.File
            Handle to the HDF5 file.
        """
        if not h5:
            h5 = "genesis4_" + self.fingerprint() + ".h5"

        if archive_file:
            h5 = archive_file



        if isinstance(h5, str):
            if "outfile" in self.output:
                shutil.copy(self.output["outfile"], h5)
                self.vprint(f"Archiving to file {h5}")

            # fname = os.path.expandvars(h5)
            # g = h5py.File(fname, 'w')
            # self.vprint(f'Archiving to file {fname}')
        else:
            g = h5

        return h5


def runGenesis(input_dict, run_file_location, run_file, workdir, nproc=2):


    
    file=run_file_location + run_file
    G = Genesis4(run_file_location + run_file, workdir=workdir, verbose = True)
    if 'alphax' in input_dict:
        G.input['main'][5]['alphax'] = input_dict["alphax"] #alpha_x
        alpha_x=input_dict["alphax"] 
    else: 
        alpha_x=G.input['main'][5]['alphax']

    if "alphay" in input_dict:
        G.input['main'][5]['alphay'] = input_dict["alphay"]
        alpha_y=input_dict["alphay"]
    else:
        alpha_y=G.input['main'][5]['alphay']
        
    if "betax" in input_dict:
        G.input['main'][5]['betax'] = input_dict["betax"]
        beta_x=input_dict["betax"]
    else:
        beta_x= G.input['main'][5]['betax']

    if "betay" in input_dict:
        G.input['main'][5]['betay'] = input_dict["betay"]
        beta_y=input_dict["betay"]
    else:
        beta_y= G.input['main'][5]['betay']

    
    if 'taper' in input_dict:

        taper=input_dict['taper']

        if 'ustart' in input_dict:

            ustart=input_dict['ustart']
        else:
            ustart = 12

        Kstart=1.7017
        nwig=130
        xlamdu=0.026
        ustop = 33
        l_und = 0.026
        
        G.apply_taper(write_filename = 'lattice' + '.lat', Kstart = Kstart, dKbyK = taper, ustart = ustart, ustop = ustop, nwig =nwig, uperiod = xlamdu, order = 1)
    
   # G.input['main'][4]['zstop'] = 3.5
    G.write_input(path=G.path, input_filename=run_file)
    print(nproc)
    G.nproc = nproc

    runscript = "salloc --partition milano --account ad:beamphysics --mem-per-cpu=4g  -N 1 -n {nproc} mpirun {command_mpi} {inputfile}".format(nproc = G.nproc, command_mpi = G.command_mpi, inputfile = G.path + '/' + run_file)

    log = []
    for path in tools.execute(runscript.split(), cwd=G.path):
        G.vprint(path, end="")
        log.append(path)
    G.vprint("Finished.")
         
    G.log = log
    G.load_output()

    twiss=[alpha_x, alpha_y, beta_x, beta_y]
 
    return G, twiss

def merrit_genesis(G):
    # Check for error

    m={}

    #if G.output['run_info']['error']:
    #    return {'error':True}
    #else:
    #    m= {'error':False}
        
    m['energy']=G.stat('field_energy')[-1]

    return m

def full_path(path):
    """
    Helper function to expand enviromental variables and return the absolute path
    """
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def evaluate_genesis(input_dict, run_file_location, run_file, workdir, archive_path, nproc=2):
    """
     If an archive_path is given, the complete evaluated Genesiss and Generator objects will be archived
     to a file named using a fingerprint from both objects.
    """
    
    d={}
    print('input dictonary', input_dict)
    G , twiss = runGenesis(input_dict=input_dict, run_file_location=run_file_location, run_file=run_file, workdir=workdir, nproc=nproc)
    
    # Evaluate merit of Genesis output
    output=merrit_genesis(G)  

    d['inputs'] = [input_dict]

#    d['twiss']=twiss

 #   if 'error' in output and output['error']:

  #       raise ValueError('run_genesis returned error in output')
   
        
    #Recreate Generator object for fingerprint, proper archiving
    # TODO: make this cleaner
    #Gen = Generator(G.input)  ## don't need this unless distgen file is used. 

    fingerprint = G.fingerprint()
    output['fingerprint'] = fingerprint

    if archive_path:
        path = full_path(archive_path)
        assert os.path.exists(path), f'archive path does not exist: {path}'
        archive_file = os.path.join(path, fingerprint+'.h5')
        output['archive'] = archive_file

       # Call the composite archive method
        archive_infolder(G, archive_file=archive_file)

    d['outputs'] = output
    return output
