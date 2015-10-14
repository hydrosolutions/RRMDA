#!/usr/bin/python
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
# matlab.engine has trouble starting when user inactive. Use system commands instead.
# import matlab.engine  # Calling matlab from python.
import subprocess  # Allows system commands.
import time        # Time references.
import os          # OS commands.
import glob        # File name pattern matching.
import tempfile    # Temporary files.
import shlex
from datetime import datetime  # To print current time.
import shutil


# Global fields, yahoo settings.
SMTP_SERVER = "smtp.mail.yahoo.com"
SMTP_PORT = 587
SMTP_SENDER = str('imomowb@yahoo.com')
SMTP_USERNAME = str('imomowb@yahoo.com')
SMTP_PASSWORD = str('iMoMo567')


# Local methods.
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
      with open(os.path.join(output_dir, 'error.txt'), 'w') as ferr:
        ferr.write("Error calling system with %s\n" % args)

    if (returnValue != 0):
      if i < numberOfTrials-1:
        logger.info('Waiting for 30 minutes before trying again.')
        try:
          for i in range(1800):
            time.sleep(1)  # Wait for 30 minutes and try again. In for-loop to allow fast exiting after keyboard interrupt.
        except KeyboardInterrupt:
          sys.exit("Keyboard interrupt. Exiting.")
  
  if (returnValue != 0) :
    logger.error('Command %s failed after %d times. Sending error mail. Aborting run.' ,command, i+1)
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
    logger.error("Error: Type missmatch. Found type(%s)=%s but expected %s" eval(variable),type(variable),expected_type)
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
  ## Specify e-mail recipients and model name.
  # Set up e-mail.
  recipients = [str('marti@hydrosolutions.ch')] # Comma-separate multiple recipients.
  subject = """imomowb SUCCESS"""
  
  # Model name.
  modelName = "Themi"

  #---------------------
  # Do not edit below!
  
  
  ## Setting up of logger.
  # Error messages are piped to the directory where this file is stored.
  output_dir = os.path.realpath(os.path.join(os.path.dirname(__file__),os.pardir))
  
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

  ## Call openDA.
  # openDA needs to be propperly installed!
  numberOfTrials = 1
  commandLine = "oda_run.sh RRMDA_Themi.oda"
  command = "RRMDA_Themi.oda"
  currentDateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  logger.info('Calling matlab with %s at %s. . . ',command,currentDateTime)
  try:
    ret = systemCall(commandLine, numberOfTrials, output_dir, out_buffer, out_path, recipients, command)
    if (ret == 0):
      logger.info('. . . done.')
  except (KeyboardInterrupt, SystemExit):
    sys.exit("Keyboard interrupt. Exiting.")

  if (ret == 1):
    logger.warning('Call to %s failed. Continues anyway.', command)


  '''
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
    logger.info('Done at %s. Sending success mail.',currentDateTime)
    message = 'Successfully ran on-line model.'
    sendEmail(recipients,subject,message)

  '''
  elif (ret==1):
    #sys.exit("Call to controlData_Themi.m fialed. Exiting.")
    logger.warning('Call to %s failed. Done.',command)
  '''


# Main.
if __name__=='__main__':
  sys.exit(main())  # Exit python upon execution of main().
  
  




