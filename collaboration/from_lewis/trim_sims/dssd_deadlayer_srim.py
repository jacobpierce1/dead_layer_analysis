import subprocess
from contextlib import contextmanager
import os
import shutil

# directory where SRIM folder is located
direc = "C:/Users/LV/Desktop/uchicago/anl/SRIM-2013/"

# define a way to change the current working directory
@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(prevdir)

def main():
	#energies = list(range(5722,5901))
	energies=[3182.69, 5770.0] # approximate means of the distribution?
	thicknesses = list(range(22,51))
	for energy in energies:
		for thickness in thicknesses:

			Z = 2 # charge of alpha
			mass = 4.0015 # mass of alpha (u)
			angle = 0.0 # angle from normal of surface
			numions = 10000 # number of ions to run
			correctionfactor = 1.0 # 1.0 for pure elements, most non-organics
			autosavenum = 100000 # TRIM will autosave after this many ions

			filename = 'C:/Users/LV/Desktop/uchicago/anl/source_deadlayer/dssd_deadlayer/%ikev_%iA_Si_transmit.txt' %(energy, thickness*thickness)
			params3 = [Z, mass, energy, angle, numions, correctionfactor, autosavenum]

			params11 = [thickness*thickness]

			params16 = [thickness*thickness]

			# edit the TRIMAUTO file to run TRIM in batch mode
			editTRIMAUTO()

			#edit the TRIM.IN file
			editTRIMIN(params3, params11, params16)

			#call TRIM.exe
			runTRIM()

			#move and rename the TRANSMIT.txt files
			savedata(filename)
			print(filename)

def editTRIMAUTO():
	with cd(direc):
		fin = open("defaultTRIMAUTO",'r')
		fout = open("TRIMAUTO",'w')

		fin.readline()

		# run in batch mode
		print("1",file=fout)

		for line in fin:
			print(line,file=fout)

		fin.close()
		fout.close()

def editTRIMIN(params3, params11, params16):
	with cd(direc):
		fin = open("defaultTRIM.IN",'r')
		fout = open("TRIM.IN",'w')

		for i in range(1,29):
			if i == 3:
				fin.readline()
				line = '%i %f %f %f %i %f %i' % tuple(params3)
				print(line,file=fout)
			elif i == 11:
				fin.readline()
				line = '       5                         0             %i' % tuple(params11)
				print(line,file=fout)
			elif i == 16:
				fin.readline()
				line = ' 1      "Layer 1"           %i  2.3212       1' % tuple(params16)
				print(line,file=fout)
			else:
				line = fin.readline().strip()
				print(line,file=fout)

		fin.close()
		fout.close()

def runTRIM():
	with cd(direc):
		p1 = subprocess.Popen(["TRIM.exe"],shell=True)
		p1.wait()

def savedata(filename):
	shutil.move(direc+"SRIMOU~1/TRANSMIT.txt",filename)

main()