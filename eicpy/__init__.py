'''
eicpy package-level code.
Used for initialising the ROOT environment with eic-smear code.
'''

import os, subprocess, sys

def initialise(logon=False):
    '''
    Initialise the ROOT environment and load eic-smear code.
    Optionally load the ROOT logon script.
    '''
    try:
        set_rootlib_path()
        import ROOT
        # Set this before we do anthing ROOT-related to ROOT
        # doesn't grab the command line arguments instead
        # of the user's own programmes.
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        # Load ROOT logon script if requested.
        if logon:
            root_logon()
        if ROOT.gSystem.Load('libeicsmear') < 0:
            raise IOError('libeicsmear cannot be located')
    except Exception as e:
        print 'initialise Error:', e

def set_rootlib_path():
    '''
    Add the ROOT library path to the Python module search path
    if it is not already in it.
    Raises an exception if the path cannot be determined or
    ROOT cannot be imported.
    '''
    try:
        # Test if the path is already set by trying to import ROOT
        import ROOT
    except:
        # That failed, so try to determine the ROOT path and add
        # it to the Python module search path.
        rootlib = find_root_lib_dir()
        if rootlib not in sys.path:
            sys.path.append(rootlib)
        # Test that it worked by trying to import ROOT again.
        # If not, this will raise an ImportError to be handled by the caller.
        import ROOT

def find_root_lib_dir():
    '''
    Attempt to locate the ROOT library directory.
    Return the path to it, or raise an exception if it cannot be found.
    '''
    try:
        rootlib = ''
        # First check with the ROOTSYS environment variable.
        rootsys = os.path.expandvars('$ROOTSYS')
        if os.path.exists(rootsys):
            # This will give the ROOT top-level directory, while
            # we want the lib subdirectory.
            rootlib = '/'.join([rootsys, 'lib'])
        else:
            # If that's empty, attempt to locate via the location
            # of the ROOT binary.
            # If 'which root' fails a CalledProcessError is raised.
            rootbin = subprocess.check_output(['which', 'root'])
            # rootbin should be of the form '$ROOTSYS/bin/root\n'.
            # We want the corresponding lib directory.
            rootlib = rootbin.replace('bin/root\n', 'lib')
        # Check that the ROOT library directory we found does
        # indeed exist. If not, raise an exception.
        if not os.path.exists(rootlib):
            raise IOError(
                ' '.join(['ROOT lib directory', rootlib, 'does not exist']))
        # Check that the lib directory contains the Python library.
        libpyroot = '/'.join([rootlib, 'libPyROOT.so'])
        if not os.path.exists(libpyroot):
            raise IOError('PyROOT library not found at ' + libpyroot)
        return rootlib
    # If 'which root' fails, catch the exception and re-raise it,
    # printing something more informative.
    except subprocess.CalledProcessError as e:
        print 'Failed to locate ROOT; is it installed?'
        raise

def root_logon():
    '''
    Load the user's ROOT logon script, if it can be located.
    Raise an IOError if it cannot.
    '''
    import ROOT
    logon = ROOT.gEnv.GetValue('Rint.Logon', '')
    if os.path.exists(logon):
        ROOT.gROOT.Macro(logon)
    else:
        raise IOError('Unable to load user\'s ROOT logon script')
