#!/usr/bin/env python

"""Script for building ROOT trees.

Thomas Burton, BNL, 5th September 2013, tpb@bnl.gov

For usage run
    build.py --help

Prerequisites:

-- ROOT installed with Python support (active by default for recent versions)
-- ROOT.py should be accessible via $PYTHONPATH
-- libeicsmear should be accessible vi $LD_LIBRARY_PATH

If in doubt, adding the ROOT library directory (e.g. $ROOTSYS/lib)
to both LD_LIBRARY_PATH and PYTHONPATH should suffice.

"""

import os
import Queue

# Names of command to zip and unzip each type of zipped file
ZIP =   {'.gz': 'gzip', '.bz2': 'bzip2'}
UNZIP = {'.gz': 'gunzip', '.bz2': 'bunzip2'}

class File:
    """Processes an input file into a ROOT tree file."""
    
    # Use Queue to manage threading.
    # This Queue stores all File objects.
    queue = Queue.Queue()

    # List of allowed input file extensions
    supported_extensions = {'.txt', '.dat'}

    def __init__(self, filename, outdir, nevents, rezip):
        """Constructor.
        
        Initialise with the input file and output information.
        Determines whether the input file is zipped.
        
        """
        name, extension = os.path.splitext(filename)
        # self.zipext stores the zipped extension if the file was zipped
        # initially, or None if it was not zipped.
        if extension in ZIP:
            self.zipext = extension
            self.rezip = rezip
            # Get the name and extension now that the zip extension is gone
            name, extension = os.path.splitext(name)
        else:
            self.zipext = None
            self.rezip = False
        self.name = name # Without extension
        self.ext = extension # File extension, not zipped extension
        self.outdir = outdir
        self.nevents = nevents

    def process(self):
        """Build the tree for the input the File was initialised with.
        
        If the file is zipped, unzip it first and rezip it after making
        the ROOT tree (unless requested not to rezip).
        
        """
        import subprocess
        if self.ext not in File.supported_extensions:
            return
        fullname = ''.join([self.name, self.ext])
        # Unzip a zipped file
        if self.zipext:
            zipped_name = ''.join([fullname, self.zipext])
            # subprocess.call should return 0 if unzipping is successful
            unzipped = not subprocess.call(
                                [UNZIP[self.zipext], '-v', zipped_name])
        else:
            unzipped = False
        # Catch any errors from tree building to we always rezip if asked to
        try:
            # Here's what we actually came to do!
            ROOT.BuildTree(fullname, self.outdir, self.nevents)
        except:
            print 'Error encountered building tree'
        if unzipped and self.rezip and self.zipext:
            print 'Rezipping', fullname
            subprocess.call([ZIP[self.zipext], '-v', fullname])


def processor():
    """Function run by each thread.
    
    Pops Files from the queue and processes them until there are no files
    remaining.
    
    """
    while True:
        file = File.queue.get()     # Get the next file from the Queue
        file.process()              # Process the file
        File.queue.task_done()      # Inform the Queue that file is done

def build_list(filelist, outdir='.', nevents=-1, nthreads = 1, rezip=True):
    """Build ROOT tree files from all the listed files.
    
    Arguments:
    filelist -- a list of file names
    outdir -- the directory to which to write ROOT files [current directory]
    nevents -- the maximum number of events per file to process [all]
    nthreads -- the maximum number of threads permitted to run [1]
    
    Zipped files (either gzip or bzip2) are unzipped then rezipped
    after the ROOT file has been made.
    
    """
    import threading
    # Generate our threads, each running the processing function
    for i in range(nthreads):
        t = threading.Thread(target=processor)
        t.daemon = True
        t.start()
    # Populate the Queue with Files
    for i in filelist:
        File.queue.put(File(i, outdir, nevents, rezip))
    # Wait for all the Queue elements to finish processing
    File.queue.join()

def parse():
    """Process command line arguments and return argparse Namespace."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
        description='build trees from files in a list or directory')
    parser.add_argument('source', help='name of file list or directory')
    parser.add_argument('-o', '--outdir', default='.', help='output directory')
    parser.add_argument('-e', '--events', type=int, default=-1,
                        help='number of events per file')
    parser.add_argument('-p', '--processes', type=int, default=1,
                        help='number of parallel processes')
    parser.add_argument('--norezip', action='store_true',
                        help='do not rezip files')
    return parser.parse_args()

# Execute buildlist on all the files in a directory, the current
# directory by default
if __name__ == "__main__":
    """Executes buildlist on all the files in a list or directory."""
    # Parse command line options
    options = parse()
    # Import ROOT after parsing arguments so build.py --help will still
    # work even if ROOT cannot be located.
    try:
        import ROOT
        # Try to tell ROOT to ignore command line arguments.
        # Otherwise ROOT will intercept command line arguments and
        # build.py --help won't work. This isn't supported in older versions
        # of ROOT, so don't complain if it can't be done.
        try:
            ROOT.PyConfig.IgnoreCommandLineOptions = True
        except AttributeError:
            pass
        # Load libeicsmear, raise and error if it can't be found
        if ROOT.gSystem.Load('libeicsmear') < 0:
            raise IOError('libeicsmear could not be located')
    # If importing ROOT and libeicsmear failed print an error and quit
    except Exception as error:
        print 'Error:', error
        quit()
    # Try to get list of file names from an input file
    if os.path.isfile(options.source):
        with open(options.source) as file:
            files = file.read().splitlines()
    # Try to get list of file names from all files in a directory
    elif os.path.isdir(options.source):
        files = [os.path.join(options.source, file)
                 for file in os.listdir(options.source)]
    # Got some bad input, complain
    else:
        print options.source, 'not recognized... quitting'
        quit()
    # Build everything
    build_list(files, options.outdir, options.events, options.processes,
               not options.norezip)