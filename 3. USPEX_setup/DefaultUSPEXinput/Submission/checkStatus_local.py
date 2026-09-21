from __future__ import absolute_import
import argparse
import glob
import os

from subprocess import check_output, CalledProcessError, STDOUT

_author_ = 'etikhonov'


def checkStatus_local(jobID):
    u"""
    This function is to check if the submitted job is done or not
    One needs to do a little edit based on your own case.
    1   : whichCluster (0: no-job-script, 1: local submission, 2: remote submission)
    Step1: the command to check job by ID. 
    Step2: to find the keywords from screen message to determine if the job is done
    Below is just a sample:
    -------------------------------------------------------------------------------
    Job id                    Name             User            Time Use S Queue
    ------------------------- ---------------- --------------- -------- - -----
    2455453.nano              USPEX            qzhu            02:28:42 R cfn_gen04 
    -------------------------------------------------------------------------------
    If the job is still running, it will show as above.
    
    If there is no key words like 'R/Q Cfn_gen04', it indicates the job is done.
    :param jobID: 
    :return: doneOr
    """

    # Step 1
    # SLURM only keeps a finished job visible to squeue for a short while
    # (MinJobAge, 300 s by default). After that - and for any id it does not
    # recognise - squeue exits non-zero with "Invalid job id specified". That
    # is the ordinary way this function learns the job is done, so it must be
    # treated as completion, not as a failure.
    #
    # Any other failure (squeue not on PATH, slurmctld unreachable) means we
    # genuinely cannot tell. The original code returned True there, which
    # advances USPEX past results that were never produced, so we raise
    # instead of guessing.
    try:
        output = str(check_output('squeue -j {}'.format(jobID),
                                  shell=True, stderr=STDOUT))
    except CalledProcessError as exc:
        message = str(exc.output).lower()
        if 'invalid job id' in message or 'invalid job specified' in message:
            output = ''  # job has left the queue: finished
        else:
            raise RuntimeError(
                'could not determine job status for {0}: {1} (squeue said: {2})'
                .format(jobID, exc, str(exc.output).strip()))
    except Exception as exc:
        raise RuntimeError(
            'could not determine job status for {0}: {1}'.format(jobID, exc))
    # Step 2
    doneOr = True
    if ' R ' in output or ' Q ' in output or ' PD ' in output or ' CG ' in output:
        doneOr = False
    if doneOr:
        for file in glob.glob('USPEX*'):
            os.remove(file)  # to remove the log file
    return doneOr

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', dest='jobID', type=int)
    args = parser.parse_args()

    isDone = checkStatus_local(jobID=args.jobID)
    print('<CALLRESULT>')
    print(int(isDone))
