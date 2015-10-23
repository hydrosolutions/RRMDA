#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Test matlab_oda_batcher_offline.py
#
# details
#
# author    Beatrice Marti, hydrosolutions ltd. (marti@hydrosolutions.ch)
#
# copyright hydrosolutions ltd, 2015
#
# license   see LICENSE
#

import unittest
import sys
import time        # Time references.
import os          # OS commands.
import shlex
from datetime import datetime  # To print current time.
import shutil      # Hihg level file operations, e.g. copy file from source to destination.
import re
# Append custom python modules to the path.
sys.path.append(os.path.realpath(os.path.join(os.path.dirname(__file__),'..','..','src','python_modules')))
import jdutil


def observationsAvailable(homeDirectory,modelName):
  '''
    Check if observations are available for the current time step.
    
    - Figure out todays date in matlab datenum format to get the filename
    - If there is a file with todays measurements read it in and figure out if
    there is discharge data available.
    
    If there are discharge observations available return true, else return false.
    '''
  date = jdutil.date_to_jd(int(datetime.now().strftime('%Y')),
                           int(datetime.now().strftime('%m')),
                           float(datetime.now().strftime('%d')))
  date = jdutil.jd_to_mjd(date)
  # Convert to matlab datenum
  date = date + 678942
  datenum = "%d" % date
  obsDir = os.path.join(homeDirectory,'app',modelName,'data','raw','imomo')
  filename = os.path.join(obsDir,(datenum+'DLoad.csv'))
  if os.path.exists(filename):
    # Scan line by line until the second entry (variable ID) of a line matches 25
    # (discharge).
    for line in open(filename):
      if re.search('^\d+\.*\d*,25,.+',line):
        return True
    return False  # No match while reading the lines.
  else:  # Return False if there is no new data.
    return False



class TestMatlabOdaBatcherMethods(unittest.TestCase):
  
  def runTest(self):
    
    self.modelName = 'ThemiTest'
    self.homeDir = os.path.realpath(os.path.join(os.path.dirname(__file__),'testEnvironment'))
    
    # Rename first observation file to have todays date.
    imomoDir = os.path.join(self.homeDir,'app',self.modelName,'data','raw','imomo')
    date = jdutil.date_to_jd(int(datetime.now().strftime('%Y')),
                             int(datetime.now().strftime('%m')),
                             float(datetime.now().strftime('%d')))
    date = jdutil.jd_to_mjd(date)
    # Convert to matlab datenum
    date = date + 678942
    datenum = "%d" % date
    filename = os.path.join(imomoDir,(datenum+'DLoad.csv'))
    shutil.copy(os.path.join(imomoDir,'736239DLoad.csv'),os.path.join(imomoDir,filename))
    
    # Call to observationsAvailable.
    self.bol = observationsAvailable(self.homeDir, self.modelName)

    # AssertTrue
    self.assertEqual(self.bol,True)

    # Test false case.
    shutil.copy(os.path.join(imomoDir,'736238DLoad.csv'),os.path.join(imomoDir,filename))
  
    # Call to observationsAvailable.
    self.bol = observationsAvailable(self.homeDir, self.modelName)
    
    # AssertTrue
    self.assertEqual(self.bol,False)


# Main.
if __name__=='__main__':
  sys.path.append(os.path.realpath(os.path.join(os.path.dirname(__file__),'..','..','src','setup')))
  import matlab_oda_batcher_offline

  unittest.main()
  



