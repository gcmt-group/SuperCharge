#!/usr/bin/env python

#!/usr/bin/env python
import os, sys, getopt
import numpy as np
from decimal import Decimal

def extend(parameter, outname = 'POSCAR_Extend'):
	outputpath = str(os.getcwd())+"/"+str(outname)
	Na = int(parameter[0])
	Nb = int(parameter[1])
	Nc = int(parameter[2])
	index = []
	x = []
	y = []
	z = []
	gridx = []
	gridy = []
	gridz = []
	poscar = open(str(os.getcwd())+"/POSCAR").readlines()
	rows = poscar[6].split()

	atoms = len(rows)
	for i in range(atoms):
		index.append(int(rows[i]))
	totoalatoms =  sum(index)
	newindex = [str(i*Na*Nb*Nc) for i in index]
	line7 = "     " + "     ".join(newindex)+'\n'
	for i in range(totoalatoms):
		x.append(Decimal(poscar[i+8].split()[0]))
		y.append(Decimal(poscar[i+8].split()[1]))
		z.append(Decimal(poscar[i+8].split()[2]))
	newx = []
	newy = []
	newz = []
	for i in range(len(x)):
		x[i] = Decimal(x[i])/Decimal(Na)
		y[i] = Decimal(y[i])/Decimal(Nb)
		z[i] = Decimal(z[i])/Decimal(Nc)
		for ia in range(Na):
			for ib in range(Nb):
				for ic in range(Nc):
					# print ("X:"+str(x[i])+"Y:"+str(y[i])+"Z:"+str(z[i])+'\n')
					newx.append(x[i] + Decimal(ia/Na))
					newy.append(y[i] + Decimal(ib/Nb))
					newz.append(z[i] + Decimal(ic/Nc))
					# print ("NX:"+str(newx[i])+"NY: "+str(newy[i])+"NZ:"+str(newz[i])+'\n')
	print (newindex)
	
	fout = open(outputpath,'w')
	for i in range (2):
		fout.writelines(poscar[i])
	for i in range (3):
		gridx.append(poscar[2].split()[i])
		gridy.append(poscar[3].split()[i])
		gridz.append(poscar[4].split()[i])
	line3 = ' '
	for strn in [Decimal(gridx[0])*Na,Decimal(gridx[1])*Na,Decimal(gridx[2])*Na]:
		line3 += '{0:>22.16f}'.format(strn)
	line3 += '\n '
	for strn in [Decimal(gridy[0])*Nb,Decimal(gridy[1])*Nb,Decimal(gridy[2])*Nb]:
		line3 += '{0:>22.16f}'.format(strn)
	line3 += '\n '
	for strn in [Decimal(gridz[0])*Nc,Decimal(gridz[1])*Nc,Decimal(gridz[2])*Nc]:
		line3 += '{0:>22.16f}'.format(strn)
	line3 += '\n '
	# for i in range(3):
	# 	for strn in [Decimal(gridx[i])*Na,Decimal(gridy[i])*Nb,Decimal(gridz[i])*Nc]:
	# 		line3 += '{0:>22.16f}'.format(strn)
	# 	line3 += '\n '
	fout.writelines(line3[:-1])
	fout.writelines(poscar[5])
	fout.writelines(line7)
	fout.writelines('Direct\n')
	for i in range(len(newx)):
		fout.writelines(' '+"{:.16f}".format(newx[i])+' '+"{:.16f}".format(newy[i])+' '+"{:.16f}".format(newz[i])+'\n')
	# fout.writelines('     '.join(str(int(rows)*Na*Nb*Nc)))
	fout.close()

def main(argv):
	outputfile_name = ''
	parameter = []
	try:
		opts, args = getopt.getopt(argv, "ci:o:")
	except getopt.GetoptError:
		print ("Error parameter")
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-c':
			print ('Copyright CaoZheng @ 2018')
			sys.exit()
		elif opt in ("-i","--inputfile"):
			pass
		elif opt in ("-o","--outputfile"):
			outputfile_name = arg
	try:
		parametercontent = open(str(os.getcwd())+"/parameter.txt").readline()
		parameter = parametercontent.split()
		if len(parameter) != 3 :
			print ("Parameter Error")
			sys.exit(1)
	except:
		print ("Cannot find Parameter")
	try:
		open(str(os.getcwd())+"/POSCAR").readline()
	except:
		print ("Cannot find POSCAR")

	if outputfile_name != '':
		extend(parameter,outputfile_name)
	else:
		extend(parameter)



if __name__ == "__main__":
	main(sys.argv[1:])