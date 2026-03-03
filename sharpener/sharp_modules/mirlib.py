# Monday, 28 September 2015
import logging
from configparser import ConfigParser
import subprocess
import os
import sys
from astropy import units as u
from astropy.io import fits, ascii
import numpy as np
import imp
deg2rad = np.pi/180.
# Its rather messy to reload the logging library, but is necessary if the logger is going to work.
imp.reload(logging)

print('Setup logger with lib.setup_logger()')

#THIS EXCEPTIONER DOES NOT WORK
def exceptioner(outputs,err):
    search_term='Fatal'
    for e in outputs:
        print(e)
        stdout_content = e.decode('utf-8').lower()
        print(stdout_content)
    print(search_term)
    if search_term in stdout_content:
            raise Exception(f"Found '{search_term}' - exiting process.")
    


def masher(task=None, **kwargs):
    '''
    masher - Miriad Task Runner
    Usage: masher(task='sometask', arg1=val1, arg2=val2)
    Example: masher(task='invert', vis='/home/frank/test.uv/', options='mfs,double', ...)
    Each argument is passed to the task through the use of the keywords.
    '''
    logger = logging.getLogger('masher')
    if task != None:
        argstr = " "
        for k in kwargs.keys():
            if str(kwargs[k]).upper() != 'NONE':
                if k == 'in_':
                    argstr += 'in=' + str(kwargs[k]) + ' '
                else:
                    k = k
                    argstr += k + '=' + str(kwargs[k]) + ' '
        cmd = task + argstr
        logger.debug(cmd)
        if ("-k" in cmd) is True:
            out = basher(cmd, showasinfo=True)
        else:
            out = basher(cmd, showasinfo=False)
        return out
    else:
        logger.critical(
            "Usage = masher(task='sometask', arg1=val1, arg2=val2...)")
        sys.exit(0)


def basher(cmd, showasinfo=False):
    '''
    basher: shell run - helper function to run commands on the shell.
    '''
    logger = logging.getLogger('basher')
    logger.debug(cmd)
    # Replacing brackets so that bash won't complain.
    #cmd = cmd.replace('""','"')
    if 'window' or 'uvrange' or 'percentage' in cmd:
        # pos = cmd.find('window')
        # cmdlist = list(cmd)
        # cmdlist.insert(pos, '')
        # cmd = ''.join(cmdlist)
        # br = cmd.find(')')
        # cmdlist = list(cmd)
        # cmdlist.insert(br+1, '')
        # cmd = ''.join(cmdlist)
        # print(cmd)
        pass
    else:
        cmd = cmd.replace("(", "\(")
        cmd = cmd.replace(")", "\)")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell=True)
    out, err = proc.communicate()
    print(out,err)
    if len(out) > 0:
        if showasinfo:
            logger.debug("Command = {} \n {}".format(cmd,out))
        else:
            logger.debug("Command = {} \n {}".format(cmd,out))
    if len(err) > 0:
        logger.debug(err)
    # NOTE: Returns the STD output.
    #iexceptioner(out, err)
    logger.debug("Returning output.")
    # Standard output error are returned in a more convenient way
    return out.splitlines()[0:-1]


class miriad:
    def __init__(self, task, **kwargs):
        '''
        DO NOT DEFINE ANY OTHER VARIABLES HERE
        '''
        self.__dict__.update(kwargs)
        self.task = task

    def __getitem__(self, key):
        return getattr(self, key)

    def keywords(self):
        masher(task=self.task+" -kw")
        masher(task=self.task + " -w")

    def help(self):
        masher(task=self.task+" -k")

    def rmfiles(self):
        logger = logging.getLogger('miriad '+self.task)
        logger.debug("Cleanup - files will be DELETED.")
        if self.task == 'invert':
            if os.path.exists(self.map):
                basher("rm -r "+self.map)
            if os.path.exists(self.beam):
                basher("rm -r "+self.beam)
        elif self.task == 'clean':
            if os.path.exists(self.out):
                basher('rm -r '+self.out)
        elif self.task == 'restor':
            if os.path.exists(self.out):
                basher('rm -r '+self.out)
        elif self.task == 'maths':
            if os.path.exists(self.out):
                basher('rm -r '+self.out)
        elif self.task == 'uvlin':
            if os.path.exists(self.out):
                basher('rm -r '+self.out)
        elif self.task == 'uvcat':
            if os.path.exists(self.out):
                basher('rm -r '+self.out)
        else:
            if os.path.exists(self.out):
                basher('rm -r '+self.out)

    def inp(self):
        logger = logging.getLogger('miriad '+self.task)
        attrs = vars(self)
        logger.info(', '.join("%s=%s" %
                              item for item in attrs.items() if item[0] is not 'logger'))

    def go(self, rmfiles=False):
        logger = logging.getLogger('miriad '+self.task)
        if rmfiles:
            self.rmfiles()
        output = masher(**self.__dict__)
#        logger.info('Completed.')
        return output
