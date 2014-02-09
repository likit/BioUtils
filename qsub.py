'''The code is a modified version of cluster_utils.py

from MISO package.

'''

import time
import subprocess

def check_job(jobid):
    '''Returns True is a job is finished,

    otherwise False.
    '''

    output = subprocess.Popen('qstat %i' %(jobid),
            shell = True,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE).communicate()
    if 'Unknown Job Id Error' in output[1]:
        return True
    else:
        return False

def wait_on_job(jobids, delay=60):
    '''Wait until a job is finished.'''

    completed = []
    while(len(completed) < len(jobids)):
        for jid in jobids:
            if jid in completed:
                continue
            status = check_job(jid)
            if status:  # job is finished
                print >> sys.stderr, 'Job #%d is finished [%d/%d]' \
                                % (jid, len(completed), len(jobids))
                completed.append(jid)
            else:
                print >> sys.stderr, 'Waiting for job #%d' % (jid)
        time.sleep(delay)

def launch_job(commands):
    '''Launches a job by executing a command.'''
    jobids = []

    for cmd in commands:
        proc = subprocess.Popen(cmd, shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE)
        # Read the job ID if it's a known cluster
        # submission system
        output = proc.communicate()
        if "." in output[0][:-1] and ">" not in output[0]:
            jobids.append(int(output[0].split(".")[0]))

    wait_on_job(jobids)

if __name__=='__main__':
    import sys
    input = sys.argv[1]
    commands = ['qsub -v input=%s blat_job.sh' % (input)]
    print commands
    launch_job(commands)
