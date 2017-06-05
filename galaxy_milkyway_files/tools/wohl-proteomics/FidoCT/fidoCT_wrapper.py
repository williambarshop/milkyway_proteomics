#from os.path import basename
import subprocess
import tempfile
import os, time
import sys, shutil
argus=sys.argv[1:]
output = argus.pop(0)
roc_output = argus.pop(0)
std_output = argus.pop(0)
err_output = argus.pop(0)
cmd=argus
working_directory = os.getcwd()
from subprocess import call


#tmp_stderr_name = tempfile.NamedTemporaryFile(dir = working_directory, suffix = '.stderr').name
#tmp_stdout_name = tempfile.NamedTemporaryFile(dir = working_directory, suffix = '.stdout').name




#with open(tmp_stderr_name, 'wb') as tmp_stderr:
#    with open(tmp_stdout_name, 'wb') as tmp_stdout:
subprocess.call(cmd, stderr=sys.stdout.fileno())

time.sleep(30)

#shutil.copy(os.path.join(working_directory,tmp_stdout_name),std_output)
#shutil.copy(os.path.join(working_directory,tmp_stderr_name),err_output)
shutil.copy(os.path.join(working_directory,"posteriors.txt"),output)
shutil.copy(os.path.join(working_directory,"rocCurve.txt"),roc_output)
