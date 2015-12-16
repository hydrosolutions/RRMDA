#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# matlab_oda_batcher.py is a wrapper for matlab functions and openDA that can be
# started with a launch daemon on mac machines.
#
# details
#
# author    Beatrice Marti, hydrosolutions ltd. (marti@hydrosolutions.ch)
#
# copyright hydrosolutions ltd, 2015
#
# license   see LICENSE
#

import logging     # Allows in-line logging with different levels.
import sys
import smtplib     # Sending e-mails via smtp.
from email.mime.text import MIMEText  # e-mail modules for message preparation.
import subprocess  # Allows system commands.
import time        # Time references.
import os          # OS commands.
import glob        # File name pattern matching.
import tempfile    # Temporary files.
import shlex
from datetime import datetime  # To print current time.
from datetime import timedelta
import shutil      # Hihg level file operations, e.g. copy file from source to destination.
import xml.dom.minidom  # Parsing of entire xml files in memory.
import re
import scipy.io    # Reading of .mat files.
# Append custom python modules to the path.
sys.path.append(os.path.realpath(os.path.join(os.path.dirname(__file__),'..','..','..','src','python_modules')))
import jdutil

# Global fields, yahoo settings.
SMTP_SERVER = "smtp.mail.yahoo.com"
SMTP_PORT = 587
SMTP_SENDER = str('imomowb@yahoo.com')
SMTP_USERNAME = str('imomowb@yahoo.com')
SMTP_PASSWORD = str('iMoMo567')

logger = None


# Local methods.
def preprocessObservations(homeDirectory,modelName):
  # Preprocess observation data. Write it to noos format for each observation
  # location.
  # Here an example of the noos format:
  #
  # #------------------------------------------------------
  # # Timeseries retrieved from the MATROOS series database
  # # Created at Mon Mar 17 10:37:58CET 2008
  # #------------------------------------------------------
  # # Location :den helder
  # # Position : (4.745356,52.966001)
  # # Source : observed
  # # Unit : waterlevel_astro
  # # Analyse time: most recent
  # # Timezone : GMT
  # #------------------------------------------------------
  # 200801010000 0.7300
  # 200801010010 0.7300
  # 200801010020 0.7200E
  #

  # Read in the file.
  date = jdutil.date_to_jd(int(datetime.now().strftime('%Y')),
                           int(datetime.now().strftime('%m')),
                           float(datetime.now().strftime('%d')))
  date = jdutil.jd_to_mjd(date)
  # Convert to matlab datenum
  date = date + 678942
  datenum = "%d" % date
  inputDir = os.path.join(homeDirectory,'app',modelName,'src','oda_model','input')
  obsDir = os.path.join(inputDir,'app',modelName,'data','raw','imomo')
  logger.info('Preprocessing observation data.')
  filename = os.path.join(obsDir,(datenum+'DLoad.csv'))
  logger.debug('filename = %s',filename)
  dischargeData = []
  numberOfSubCatchments = 8
  if os.path.exists(filename):
    logger.debug('File exists.')
    # Scan line by line until the second entry (variable ID) of a line matches 25
    # (discharge).
    for line in open(filename):
      #print line
      if re.search('.+,25,.+',line):
        dischargeData.append(line)
        #print 'yeah'

  
    # Get maximum value -> reverse order.
    xmllist = []
    index = '9'
    for item in reversed(dischargeData):
      columns = re.split(',',item)
      #if index == columns[5]:
      #  continue  # Skip double entries. ??
      #else:
      index = columns[5]
      time = columns[2]  # Time in format: '2015-10-01 07:05:25,3'
      time = time[0:4]+time[5:7]+time[8:10]+'0000'  # reformat time to yyyymmddHHMM
      value = columns[0]
      #logger.debug('index = %s, value = %s',index,value)
      newfilename = 'imomoQ'+str(index)+'.txt'
      #logger.debug('newfilename = %s',newfilename)
      newfile = os.path.join(homeDirectory,'app',modelName,'src','oda_observer',newfilename)
      file = open(newfile,'w')
      writestringcom = '#--------------'+os.linesep
      file.write(writestringcom)
      writestring = '# Location : Q'+index+os.linesep
      file.write(writestring)
      file.write(writestringcom)
      writestring = time+' '+value+os.linesep
      file.write(writestring)
      file.close()

  else:
    logger.error('File not found: %s.',filename)


def generateObservations(homeDirectory, modelName,numberOfSubCatchments):
  '''
  Use last computed discharge and add some noise to it (currently + 20).
  '''
  # Get simulation results of previous time step from E.mat.
  resultdir = os.path.join(homeDirectory,'app',modelName,'src','oda_model','input','app',modelName,'results','nowcast')
  logger.debug('resultdir = %s.',resultdir)
  filename = max(glob.iglob(os.path.join(resultdir,'*.mat')), key=os.path.getctime)
  Emat = scipy.io.loadmat(os.path.join(resultdir,filename))
  E = Emat.get('E')
  for i in range(0,numberOfSubCatchments):
    E[0 + i*numberOfSubCatchments][0] = E[0 + i*numberOfSubCatchments][0] + 20  # Here "noise" is added.

  # Emat['E'] = E
  # scipy.io.savemat(filename,Emat)
 
  # Get todays observation file.
  date = jdutil.date_to_jd(int(datetime.now().strftime('%Y')),
                           int(datetime.now().strftime('%m')),
                           float(datetime.now().strftime('%d')))
  date = jdutil.jd_to_mjd(date)
  # Convert to matlab datenum
  date = date + 678942
  datenum = "%d" % date
  obsDir = os.path.join(homeDirectory,'app',modelName,'src','oda_model','input','app',modelName,'data','raw','imomo')
  logger.debug('obsDir = %s',obsDir)
  filename = os.path.join(obsDir,(datenum+'DLoad.csv'))
  if os.path.exists(filename):
    # Scan line by line until the second entry (variable ID) of a line matches 25
    # (discharge).
    with open(filename,'a') as obsfile:
      for i in range(0,numberOfSubCatchments):
        val = str(E[0+i*numberOfSubCatchments][0])
        line = val+',25,'+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'.0,-3.51185,36.783066,'+str(i+1)+',400,offline_location,offline_user'+os.linesep
        obsfile.write(line)
  return True


def checkTimeMat(homeDirectory,modelName):
  ''''
  Check contents of time.mat and adapt if the dates are older than 1 day or the 
  model doesn't run.
  '''
  logger.info('Checking time.mat . . . ')
  restartDir = os.path.join(homeDirectory,'app',modelName,'src','oda_model','input')
  filename = os.path.join(restartDir,'time.mat')
  endtime = []
  if os.path.exists(filename):
    timemat = scipy.io.loadmat(filename)
    endtime = timemat.get('endTime')
  year = float(endtime[0][0:4])
  month = float(endtime[0][4:6])
  day = float(endtime[0][6:8])
  endtime = jdutil.date_to_jd(year,month,day)
  endtime = jdutil.jd_to_mjd(endtime)
  logger.debug('endTime = %d',endtime)

  date = jdutil.date_to_jd(int(datetime.now().strftime('%Y')),
                           int(datetime.now().strftime('%m')),
                           float(datetime.now().strftime('%d')))
  date = jdutil.jd_to_mjd(date)
  if (endtime == date):
    logger.debug('Yeah, dates match.')
  else:
    logger.debug('Adapting times . . . ')
    outstarttime = date - 1
    outendtime = jdutil.mjd_to_jd(date)
    outendtimes = jdutil.jd_to_date(outendtime)
    outendtime = "%d%d%d0000"%(outendtimes[0],outendtimes[1],round(outendtimes[2]))
    outstarttime = jdutil.mjd_to_jd(outstarttime)
    outstarttimes = jdutil.jd_to_date(outstarttime)
    outstarttime = "%d%d%d0000"%(outstarttimes[0],outstarttimes[1],round(outstarttimes[2]))
    timemat['startTime'] = outstarttime
    timemat['endTime'] = outendtime
    step = float(timemat.get('step'))
    timemat['step'] = step
    scipy.io.savemat(filename,timemat)


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
  obsDir = os.path.join(homeDirectory,'app',modelName,'src','oda_model','input','app',modelName,'data','raw','imomo')
  filename = os.path.join(obsDir,(datenum+'DLoad.csv'))
  if os.path.exists(filename):
    # Scan line by line until the second entry (variable ID) of a line matches 25
    # (discharge).
    logger.debug('Reading todays observation file %s.',filename)
    for line in open(filename):
      if re.search('^\d+\.*\d*,25,.+',line):
        return True
    logger.debug('No discharge observation today.')
    return False  # No match while reading the lines.

  else:  # Return False if there is no new data.
    return False
    logger.debug('No observation file for today.')


def setupOpenDaRun(homeDirectory,modelName):
  '''
  Setup call to openDA.
  
  Cleaning up input directory
  Copy input files to model/input/
  Copy observation file to model/observation/
  Write RRMDA.oda
  Write observation.xml
  
  returns 0 upon success, 1 uppon failure.
  '''
  
  logger.info('Setting up openDA run.')
  inputDir = os.path.join(homeDirectory,'app',modelName,'src','oda_model','input')
  
  ## Cleaning up: Removing all (old) files and directories in inputDir.
  files = [f for f in os.listdir(inputDir)]
  for f in files:
    filePath = os.path.join(inputDir,f)
    if os.path.isfile(filePath):
      try:
        os.remove(filePath)
        logger.info('Removing file %s.',filePath)
      except Exception, e:
        logger.error('Problem removing file %s e: %s',f,e)
        return 1
    elif os.path.isdir(filePath):
      try:
        shutil.rmtree(filePath)
        logger.info('Removing directory %s.',filePath)
      except Exception, e:
        logger.error('Problem removing directory %s e: %s',filePath,e)
        return 1

  
  ## Copy input files from resources/restart
  restartPath = os.path.join(homeDirectory,'app',modelName,'resources','restart')
  srcPath = os.path.join(homeDirectory,'src')
  appPath = os.path.join(homeDirectory,'app',modelName)
  try:
    shutil.copy(os.path.join(restartPath,'E_oda.mat'),inputDir)
    #shutil.copy(os.path.join(restartPath,'S0G0.mat'),inputDir)
    shutil.copy(os.path.join(restartPath,'time.mat'),inputDir)
    shutil.copytree(os.path.join(srcPath),os.path.join(inputDir,'src'))
    shutil.copytree(os.path.join(appPath,'data'),os.path.join(inputDir,'app',modelName,'data'))
    shutil.copytree(os.path.join(appPath,'prm'),os.path.join(inputDir,'app',modelName,'prm'))
    shutil.copytree(os.path.join(appPath,'resources'),os.path.join(inputDir,'app',modelName,'resources'))
    shutil.copytree(os.path.join(appPath,'results'),os.path.join(inputDir,'app',modelName,'results'))
    runmodelname = 'runModelOpenDA_'+modelName+'.m'
    appsrcpath = os.path.join(inputDir,'app',modelName,'src')
    if not os.path.exists(appsrcpath):
      os.makedirs(appsrcpath)  # Create directory.
    shutil.copy(os.path.join(appPath,'src',runmodelname),os.path.join(inputDir,'app',modelName,'src'))
    shutil.copy(os.path.join(appPath,'src','setup.mat'),os.path.join(inputDir,'app',modelName,'src'))
    shutil.copy(os.path.join(appPath,'src','oda_adaptSetupMat.m'),os.path.join(inputDir,'app',modelName,'src'))
    logger.info('Copying restart files.')
  except Exception, e:
    logger.error('Problem copying restart files. Error: %s.',e)
    return 1

  ''' Not needed.
  ## Copy most recent observation file to model/observation/
  # Get the list of files in the source directory.
  observation_path = os.path.join(inputDir,'app',modelName,'data','raw','imomo')
  # In case there are characters in the path that need to be escaped put quotes around the path.
  observation_path = observation_path + '/'
  logger.debug('setupOpenDArun : observation_path = %s' , observation_path)
  try:
    files = [f for f in os.listdir(observation_path)]
    # Get the extension of the last (= most recent) file in the list.
    index = -1
    filename, file_extension = os.path.splitext(files[index])
    while not (file_extension == '.csv'):
      index = index - 1
      filename, file_extension = os.path.splitext(files[index])
    # Copy this file to homeDirectory/app/modelName/src/model/input/
    try:
      shutil.copy(os.path.join(observation_path,(filename+file_extension)),inputDir)
      logger.info('Copying file %s to directory %s.', os.path.join(observation_path,(filename+file_extension)),inputDir)
    except IOError as ioe:
      logger.error('IOError %d: %s.', ioe.errno, ioe.strerror)
      return 1
  except OSError as ose:
    logger.error('OSError %d: %s.', ose.errno, ose.strerror)
    return 1
  '''

  ## Write RRMDA.oda
  try:
    DOMTree = xml.dom.minidom.parse("RRMDA_Themi.oda")
    openDaApplication = DOMTree.documentElement
    restartInFile = openDaApplication.getElementsByTagName("restartInFile")
    restartInFileName = restartInFile[0].childNodes[0].nodeValue
    restartOutFile = openDaApplication.getElementsByTagName("restartOutFilePrefix")
    restartOutFilePrefix = restartOutFile[0].childNodes[0].nodeValue
    logger.info('Parsing RRMDA_Themi.oda')
    logger.debug('restartOutFilePrefix = %s',restartOutFilePrefix)
    # Change fileName
    date = jdutil.date_to_jd(int(datetime.now().strftime('%Y')),
                             int(datetime.now().strftime('%m')),
                             float(datetime.now().strftime('%d')))
    date = jdutil.jd_to_mjd(date)
    # Convert to matlab datenum
    date = date + 678942
    datenum = "%d" % date
    logger.debug('datenum of today = %d', date)
    restartInFileName = restartOutFilePrefix + datetime.now().strftime('%Y%m%d') + '.zip'
    logger.debug('new restartInFileName = %s',restartInFileName)
    restartInFile[0].childNodes[0].replaceWholeText(restartInFileName)
    # Adapt result file name.
    resultWriter = openDaApplication.getElementsByTagName("resultWriter")[0]
    configFile = resultWriter.getElementsByTagName("configFile")[0]
    resultFileName = configFile.childNodes[0].nodeValue
    resultFileName = 'RRMDA_EnKF_' + datetime.now().strftime('%Y%m%d') + '.m'
    configFile.childNodes[0].replaceWholeText(resultFileName)
    logger.debug('resultFileName: %s',resultFileName)
    # Write.
    file = open('RRMDA_Themi.oda','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing RRMDA.oda. Error: %s.',e)
    return 1

  ## Check if restart file is present. Create empty restart file if not.
  fn = os.path.join(homeDirectory,'app',modelName,'src',restartInFileName)
  if not os.path.exists(fn):
    try:
      f = open(fn,'a')
      f.close()
    except Exception, e:
      logger.error('Problem creating file %s.',fn)

  ## Write current simulation time span to noise model configuration files.
  now_date = datetime.now() - timedelta(days=1)
  next_date = now_date + timedelta(days=1)
  end_date = now_date + timedelta(days = 4)
  simTimespan = (now_date.strftime('%Y')+
                 now_date.strftime('%m')+
                 now_date.strftime('%d')+
                 '0000,'+
                 next_date.strftime('%Y')+
                 next_date.strftime('%m')+
                 next_date.strftime('%d')+
                 '0000,...,'+
                 end_date.strftime('%Y')+
                 end_date.strftime('%m')+
                 end_date.strftime('%d')+
                 '0000')
  
  try:
    logger.info('Parsing oda_model/NoiseModel_alpha1.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_alpha1.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_alpha1.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_alpha1.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_alpha2.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_alpha2.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_alpha2.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_alpha2.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_d.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_d.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_d.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_d.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_Smax.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_Smax.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_Smax.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_Smax.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_discharge.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_discharge.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_discharge.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_discharge.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_ETact.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_ETact.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_ETact.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_ETact.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_GroundwaterStorage.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_GroundwaterStorage.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_GroundwaterStorage.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_GroundwaterStorage.xml')
    return 1

  try:
    logger.info('Parsing oda_model/NoiseModel_soilMoisture.xml')
    DOMTree = xml.dom.minidom.parse("oda_model/NoiseModel_soilMoisture.xml")
    nm = DOMTree.documentElement
    simTimespanObj = nm.getElementsByTagName("simulationTimespan")
    simTimespanObj[0].childNodes[0].replaceWholeText(simTimespan)
    file = open('oda_model/NoiseModel_soilMoisture.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing NoiseModel_soilMoisture.xml')
    return 1


  ''' Not needed.
  # Write Observations.xml
  try:
    DOMTree = xml.dom.minidom.parse("observer/Observations.xml")
    observer = DOMTree.documentElement
    timeSeries = observer.getElementsByTagName("timeSeries")
    logger.info('Parsing Observations.xml')
    for node in timeSeries:
      node.childNodes[0].replaceWholeText(datenum + 'DLoad.csv')
    file = open('observer/Observations.xml','w')
    DOMTree.writexml(file)
    file.close()
  except Exception, e:
    logger.error('Problem parsing and writing Observations.xml. Error: %s.',e)
    return 1
  '''
  
  return 0

def systemCall(args, numberOfTrials, output_dir, out_buffer, out_path, recipients, command):
  '''
  Call to system
  
  @arg commandLine String with system command line.
  @arg numberOfTrials Integer specifying how many times the command is executed.
  @arg output_dir Path where the log is stored.
  @arg out_buffer Temporary output buffer.
  @arg out_path Path where temporary output is stored.
  @arg recipients List of strings with mail recipients.
  @arg command String with short name for command.
  
  @return returnValue Boolean 0 for successfull call to command and 1 for failure.
  '''
  logger.debug('Arguments: %s', args)
  
  for i in range(0,numberOfTrials):  # Try to execute args numberOfTrials times.
    try:
      process = subprocess.check_call(args, stdout=out_buffer, stderr=subprocess.STDOUT,shell=True)
      logger.debug('process = %s',process)
      logger.debug('Command %s successful, continuing.', command)
      returnValue = process
      return returnValue
    except (KeyboardInterrupt, SystemExit):  # Allows keyboard interrupt in shell.
      sys.exit("Keyboard interrupt while calling matlab. Exiting.")
    except OSError as e:
      logger.error('OSError %d: %s.', e.errno, e.strerror)
      returnValue = 1
      return returnValue
    except subprocess.CalledProcessError:  # check_call returned non-zero returncode.
      logger.warning('Command %s failed for the %d(st,nd,rd,th) time.', command, i+1)
      returnValue = 2
      if command == "RRMDA_Themi.oda":  # For some reason oda_run does not return 0 for success.
        return 0
      with open(os.path.join(output_dir, 'error.txt'), 'w') as ferr:
        ferr.write("Error calling system with %s%s"%(args,os.linesep))

    if (returnValue != 0):
      if i < numberOfTrials-1:
        logger.info('Waiting for 30 minutes before trying again.')
        try:
          for i in range(1800):
            time.sleep(1)  # Wait for 30 minutes and try again. In for-loop to allow fast exiting after keyboard interrupt.
        except KeyboardInterrupt:
          sys.exit("Keyboard interrupt. Exiting.")
  
  if (returnValue != 0) :
    logger.error('Command %s failed after %d times. Sending error mail.' ,command, i+1)
    message = """Failed in call to %s in %d trials.""" % (command,i+1)
    subject = """imomowb ERROR"""
    sendEmail(recipients,subject,message)

  return returnValue


def cleanUp():
  '''
  Removes files that are no longer needed in order to avoid filling of disk.
  
  '''
  
  '''
  for filename in glob.glob('data/raw/imomo/*') :
    os.remove( filename )
  for filename in glob.glob('data/raw/ndvi/*') :
    os.remove( filename )
  for filename in glob.glob('data/raw/*.mat') :
    os.remove( filename )
  '''



def testType(variable,expected_type):
  """
  Test type of variable against expected_type. Prints an error message 
  in case the types do not match and exits the script.
  
  @param variable to be type-tested.
  @param expected_type Python type.
  """
  if (type(variable) == expected_type) == False:
    logger.error('Error: Type missmatch. Found type(%s)=%s but expected %s' ,eval(variable),type(variable),expected_type)
    exit()


def sendEmail(recipients,subject,message):
  """ 
  send_email(recipients,subject,message)
  
  @param recipients List of strings containing e-mail
                    addresses of recipients.
  @param message    String with the body of the e-mail.
  @throws SMTPException
  
  Connect to server and send e-mail. Two trials in case internet is down temporarily.
  """
  
  # Test arguments, number and type.
  testType(recipients,list)
  testType(recipients[0],str)
  testType(subject,str)
  testType(message,str)
  
  msg = MIMEText(message)
  msg['Subject'] = subject
  msg['From'] = SMTP_SENDER
  # If there are multiple recipients join the addresses.
  if len(recipients) > 1:
    msg['To'] = ", ".join(recipients)
  else:
    msg['To'] = recipients

  server = smtplib.SMTP(SMTP_SERVER,SMTP_PORT)
  try:
    server.starttls()  # Encript connection.
    server.login(SMTP_USERNAME,SMTP_PASSWORD)
    server.sendmail(SMTP_SENDER, recipients, msg.as_string())
    logger.info("Successfully sent e-mail.")
  except:
    # Sleep and try again.
    try:
      logger.info("Failed to send e-mail the first time. Sleep for 10 min and try again.")
      for i in range(600):
        time.sleep(1)  # Sleep for 10 minutes.
      server.starttls()
      server.login(SMTP_USERNAME,SMTP_PASSWORD)
      server.sendmail(SMTP_SENDER,recipients,msg.as_string())
      logger.info("Successfully sent e-mail.")
    except KeyboardInterrupt:
      logger.info("Keyboard interrupt. Exiting.")
      sys.exit("Keyboard interrupt. Exiting.")
    except:
      logger.error("Error: unable to send e-mail")
  finally:
    server.quit

def main():
  """
  main()
  
  Specification of e-mail-recipients and model name.
  Setting up of logger and matlab commands.
  Procedure of main with system calls (sc):
    - sc matlab getRaw_Themi
    - sc matlab processRaw_Themi
    - sc openDA runModel_Themi
   (- sc matlab sendtoDB_Themi)
   (- sc matlab matlabMail_Themi)
   (- sc matlab controlData_Themi)
    - wrapping up
  System calls that are commented out in the offline version are in brackets.
  """
  global logger
  
  ## Specify e-mail recipients and model name.
  # Set up e-mail.
  recipients = [str('marti@hydrosolutions.ch'),str('martibeatrice@gmail.com')] # Comma-separate multiple recipients.
  subject = """imomowb SUCCESS"""
  
  # Model name.
  modelName = "Themi"

  #---------------------
  # Do not edit below!
  
  
  ## Setting up of logger.
  # Error messages are piped to the directory where this file is stored.
  output_dir = os.path.realpath(os.path.join(os.path.dirname(__file__),os.pardir,'src'))
  
  logger = logging.getLogger('matlab_oda_batcher')
  logger.setLevel(logging.DEBUG)
  consoleHandler = logging.StreamHandler()
  consoleHandler.setLevel(logging.DEBUG)
  consoleHandler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
  logger.addHandler(consoleHandler)

  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

  logger.info('-------------------------------------------------')
  logger.info('%s Starting off-line model . . . ',currentDateTime)
  logger.info('-------------------------------------------------')

  logger.debug('Error message and log directory: %s', str(output_dir))
  
  (out_buffer, out_path) = tempfile.mkstemp()  # Generates temporary output file in os-dependent path.
  logger.debug('tempfile out_buffer: %s',out_buffer)
  logger.debug('tempfile out_path: %s',str(out_path))
  
  
  ## Setting up of matlab commands.
  # The matlab user path needs to be adapted.
  homeDir = os.path.realpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
  logger.debug('homeDir: %s',homeDir)
  
  # Test if modelName describes a directory name.
  if not os.path.isdir(os.path.join(homeDir,'app',modelName)):
    logger.error('Error: %s is not a directory.',os.path.join(homeDir,'app',modelName))
    print('Error: %s is not a directory.'%os.path.join(homeDir,'app',modelName))
    sys.exit()
  
  addMatlabPath = os.path.join(homeDir,'app',modelName,'src') + ":" + \
                  os.path.join(homeDir,'src') + ":" + \
                  os.path.join(homeDir,'src','data') + ":" + \
                  os.path.join(homeDir,'src','dbase') + ":" + \
                  os.path.join(homeDir,'src','enkf') + ":" + \
                  os.path.join(homeDir,'src','enkf','prm') + ":" + \
                  os.path.join(homeDir,'src','rrm') + ":" + \
                  os.path.join(homeDir,'src','rrm','budyko') + ":" + \
                  os.path.join(homeDir,'src','setup') + ":"
  # In case there are characters in the path that need to be escaped put quotes around the path.
  addMatlabPath = "\"" + addMatlabPath + "\""
  matlabPathCommand = "export MATLABPATH="+addMatlabPath+"$MATLABPATH"
  logger.debug('MatlabPath: %s',matlabPathCommand)

  ret = 1  # Initialize ret with value for "failure".

  ## Call getRaw.m
  # Set up command.
  numberOfTrials = 2
  commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"getRaw_Themi\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
  logger.debug('commandLine: %s',commandLine)
  command = "getRaw_Themi"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s . . . ',command, currentDateTime)
  # Calling matlab.
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")
  finally:
    os.remove(out_path)  # Remove temporary directory.
    logger.debug('removing %s',str(out_path))
  
  if (ret==1):
    #sys.exit("Call to getRaw_Themi.m failed. Exiting.")
    logger.warning('Call to failed. Continue anyway.',command)


  ## Call processRaw.m
  numberOfTrials = 1
  commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"processRaw_Themi\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
  command = "processRaw_Themi"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s. . . ',command, currentDateTime)
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")

  if (ret==1):
    #sys.exit("Call to processRaw_Themi.m fialed. Exiting.")
    logger.warning('Call to %s failed. Continue anyway.',command)



  ## Setup call to openDA.
  # Copy input files to model/input/
  # Copy observation file to observer/
  # Write RRMDA.oda
  # Write observation.xml

  ret = setupOpenDaRun(homeDir,modelName)
 
  logger.debug('Checking if observations are available, generating if not.')
  ret = 0;
  if not observationsAvailable(homeDir,modelName):  # If True assimilate data. Simple forecast if False.
  
    ## Generate observations.
    ret = 0;
    logger.info('Generating Observations . . . ')
    numberOfSubCatchments = 8
    ret = generateObservations(homeDir,modelName,numberOfSubCatchments)
    if (ret==0 or ret==True):
      logger.info('. . . done.')
    else:
      logger.warning('Failed.')

  # Checking time.mat
  checkTimeMat(homeDir, modelName)

  # Copy observations to nooa-data format.
  preprocessObservations(homeDir,modelName)



  ## Call openDA.
  # openDA needs to be propperly installed! OpenDA calls matlab with runModelOpenDA_Themi.m
  numberOfTrials = 1
  commandLine = "oda_run.sh RRMDA_Themi.oda"
  command = "RRMDA_Themi.oda"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling openDA with %s at %s. . . ',command,currentDateTime)
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")

  if (ret != 0):
    logger.warning('Call to %s failed. Continues anyway.', command)

  
  
  ## Clean up openDA run.
  # Call matlab function that concatenates the model/output/workX/ files to
  # re-create a complete input file.
  numberOfTrials = 1
  commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"wrappingUpOpenDARun\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
  command = "wrappingUpOpenDARun"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s . . . ',command, currentDateTime)
  ret = 2
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")
  if (ret==1):
    logger.warning('Call to %s failed. Continue anyway.',command)
  if (ret==2):
    logger.debug('ret=2')
  
  
  '''
  else:
  
    ## Call runModel.m
    numberOfTrials = 1
    commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"runModel_Themi\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
    command = "runModel_Themi"
    currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    logger.info('Calling matlab with %s at %s . . . ',command, currentDateTime)
    ret = 2
    try:
      ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
      if (ret == 0):
        logger.info('. . . done.')
    except (KeyboardInterrupt, SystemExit):
      sys.exit("Keyboard interrupt. Exiting.")

    if (ret==1):
      #sys.exit("Call to runModel_Themi.m fialed. Exiting.")
      logger.warning('Call to %s failed. Continue anyway.',command)

    if (ret==2):
      logger.debug('ret=2')
  '''


  ## Offline mode: no sending of data to data base, no sending of matlab mails,
  ## and not controling of data.
  '''
  ## Call sendtoDB_Themi.m
  numberOfTrials = 1
  commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"sendtoDB_Themi\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
  command = "sendtoDB_Themi"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s . . . ',command, currentDateTime)
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")
        
  if (ret==1):
    #sys.exit("Call to sendtoDB_Themi.m fialed. Exiting.")
    logger.warning('Call to %s failed. Done.', command)


  ## Call MatlabMail_Themi.m
  numberOfTrials = 1
  commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"MatlabMail_Themi\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
  command = "MatlabMail_Themi"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s . . . ',command, currentDateTime)
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")
  
  if (ret==1):
    #sys.exit("Call to MatlabMail_Themi.m fialed. Exiting.")
    logger.warning('Call to %s failed. Continue anyway.',command)


  ## Call controlData_Themi.m
  numberOfTrials = 1
  commandLine = matlabPathCommand+" && "+"/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nosplash -r \"controlData_Themi\"".format(os.path.realpath(os.path.dirname(__file__)),output_dir)
  command = "controlData_Themi"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s . . . ',command, currentDateTime)
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")
  '''

  ## Wrapping up.
  if (ret == 0):
    currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    message = 'Successfully ran off-line model RRMDA_Themi on imac.'
    sendEmail(recipients,subject,message)
    logger.info('Done at %s. ',currentDateTime)
  

  '''
  elif (ret==1):
    #sys.exit("Call to controlData_Themi.m fialed. Exiting.")
    logger.warning('Call to %s failed. Done.',command)
  '''


# Main.
if __name__=='__main__':
  main()
  



