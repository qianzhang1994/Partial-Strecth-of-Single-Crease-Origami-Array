"""
Run Abaqus using inp file.
When one job has been completed, other procedure are carried out.
The Cpus Information cpu = 2
"""
import os
class RunABAQUS(object):
    def __init__(self, jobname):
        self.jobname = jobname

    def Simulate(self):
        cmd = "abaqus interactive job=" + self.jobname + " cpus=2"
        # cmd = "abaqus job=" + self.jobname + " cpus=4"
        print(cmd)
        os.system(cmd)
        name = self.jobname + ".odb"
        if os.path.exists(name):
            pass
        else:
            raise ValueError("try again for the" + str(name))


