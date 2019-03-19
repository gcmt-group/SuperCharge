#!/usr/bin/env python
# Version 1.0.2 Created by Z.Cao 2/Feb/2019
import os, sys, getopt
import numpy as np
def cut(line):
	newline = []
	for i in line.split():
		if i[0] == '-':
			newstr = i[0]+i[2]+i[1]
			if i[-3:] == '-01':
				newstr = i[0]+i[2]+i[1]+i[3:13]+'E'+'+00'
			elif i[-3] == '+':
				newstr = i[0]+i[2]+i[1]+i[3:13]+'E+'+str(int(i[-2:])+1).zfill(2)
			elif i[-3] == '-':
				newstr = i[0]+i[2]+i[1]+i[3:13]+'E-'+str(int(i[-2:])-1).zfill(2)
			newline.append(newstr)
		else:
			newline.append(i)
	final = ' ' + ' '.join(newline)+'\n'
	return final

def makeuplist():
	templist = []
	file_path = str(os.getcwd())+ '/atom_selection.txt'
	CHGCAR_path = str(os.getcwd())+ '/CHGCAR'
	file_content = open(file_path, 'r') 
	sourcelines = file_content.readlines()
	temp = open(CHGCAR_path).readlines()[6].split()
	tempnum = [int(i) for i in temp]
	totallines = sum(tempnum)

	if sourcelines[0] == 'list\n' or sourcelines[0] == 'List\n' or sourcelines[0] == 'LIST\n':
		for index in range(len(sourcelines)-1):
			templist.append(sourcelines[index+1].strip('\n').zfill(4))
	elif sourcelines[0] == 'range\n' or sourcelines[0] == 'Range\n'or sourcelines[0] == 'RANGE\n':
		print ("Finding Selected Atom...")
		rangemin = float(sourcelines[2].strip('\n'))
		rangemax = float(sourcelines[3].strip('\n'))
		if sourcelines[1] == 'x\n' or sourcelines[1] == 'X\n':
			for j in range(totallines):
				if float(open(CHGCAR_path).readlines()[8+j].split()[0])>=rangemin and float(open(CHGCAR_path).readlines()[8+j].split()[0])<rangemax:
					templist.append(str(j+1).zfill(4))
		elif sourcelines[1] == 'y\n' or sourcelines[1] == 'Y\n':
			for j in range(totallines):
				if float(open(CHGCAR_path).readlines()[8+j].split()[1])>=rangemin and float(open(CHGCAR_path).readlines()[8+j].split()[1])<rangemax:
					templist.append(str(j+1).zfill(4))
		elif sourcelines[1] == 'z\n' or sourcelines[1] == 'Z\n':
			for j in range(totallines):
				if float(open(CHGCAR_path).readlines()[8+j].split()[2])>=rangemin and float(open(CHGCAR_path).readlines()[8+j].split()[2])<rangemax:
					templist.append(str(j+1).zfill(4))
	else:
		print ("Atom_selection.txt Error: MUST LABEL range or list in first line")
		
	return templist


def addlines(linelist, lineindex, needlineindex):
	output = ''
	count = 1
	for i in needlineindex:
		startline = lineindex[int(i)-1]
		endline = lineindex[int(i)]-1
		flag = startline

		templine = linelist[flag].split()
		if count>9:
			templinecontent = templine[0] + ' '+ templine[1] + ' ' + str(count) + ' ' +templine[3] + '\n'
		else:
			templinecontent = templine[0] + ' '+ templine[1] + '  ' + str(count) + ' ' +templine[3] + '\n'
		count += 1
		output += templinecontent

		for j in range(endline-startline):
			output += linelist[flag+1]
			flag += 1
	return output

def mergy(filename, outname = 'CHGCAR_OUTPUT'):
	headercontent = ''
	footercontent = ''
	p = len(filename)
	CHGCAR_path = str(os.getcwd())+ '/CHGCAR'
	#Add the first 6 lines in head
	for i in range(6):
		headercontent += open(CHGCAR_path).readlines()[i]
	line_7 = open(CHGCAR_path).readlines()[6]
	rows = 0
	for s in line_7.split():
		rows += int(s)
	factor = int(rows/p)
	#Add line #7
	linenumber = [int(int(i)/factor+0.5) for i in line_7.split()]
	print (linenumber)
	# headercontent += '     ' + str(int(int(line_7.split()[0])/factor)) + '     '+ str(int(int(line_7.split()[1])/factor)) + '\n'
	line7 = [str(i) for i in linenumber]
	headercontent += '     '+'     '.join(line7) + '\n'
	# open(CHGCAR_path).readlines()[0:5]
	# Add line #8
	headercontent += open(CHGCAR_path).readlines()[7]

	for i in range(p):
		nextrow = int(filename[i]) + 7
		headercontent += open(CHGCAR_path).readlines()[nextrow]
	# Add final line
	headercontent += open(CHGCAR_path).readlines()[8+rows]
	finalrow = open(CHGCAR_path).readlines()[9+rows].splitlines()
	headercontent += finalrow[0]
	# Caculate the main martix
	flag = 0


	for i in range(len(filename)):
		path = str(os.getcwd())+"/BvAt"+filename[i]+".dat"
		try:
			x = open(path).readline()
			print("Processing Atom #"+str(filename[i]))
			a = np.loadtxt(path , skiprows = rows + 10)
			if flag == 0:
				res = a
				flag = 1
			else:
				res += a
		except:
			print("Error: Problem File "+"/BvAt"+filename[i]+".dat")

	# Add footer
	print ('Shape is '+str(res.shape[0]))
	footerlist = open(CHGCAR_path).readlines()[res.shape[0]+10+rows : len(open(CHGCAR_path).readlines())]
	footerindex = []
	for i, val in enumerate(footerlist):
		if val.find('augmentation') == 0:
				footerindex.append(i)
	footerindex.append(len(footerlist))

	footercontent = addlines(footerlist,footerindex,filename)
	# print (footercontent)


	# np.savetxt('chgcar.temp',res,fmt=' %.11E %.11E %.11E %.11E %.11E',header=headercontent, footer=footercontent[:-1], comments='')
	np.savetxt('chgcar.temp',res,fmt=' %.11E %.11E %.11E %.11E %.11E',comments='')
	print ("Final Check Negative Value")

	temppath = open(str(os.getcwd())+"/chgcar.temp").readlines()

	fout = open(str(os.getcwd())+"/"+outname,'w')
	fout.writelines(headercontent)
	fout.writelines('\n')
	for i in temppath:
		fout.writelines(cut(i))
	fout.writelines(footercontent[:-1])	
	fout.close()
	os.remove(str(os.getcwd())+'/chgcar.temp')

def check(filename):
	for i in range(len(filename)):
		path = str(os.getcwd())+"/BvAt"+filename[i]+".dat"
		content = open(path).readlines()

		lastrowlist = content[-1].split()
		# print (lastrowlist)

		#Stupid Method haha
		if len(lastrowlist)==5:
			print ("Atom #" + str(filename[i]) + " Matrix Check Pass")
		else:
			if len(lastrowlist)==4:
				lastrowlist.append(lastrowlist[-1])
			elif len(lastrowlist)==3:
				lastrowlist.append(lastrowlist[-1])
				lastrowlist.append(lastrowlist[-1])
			elif len(lastrowlist)==2:
				lastrowlist.append(lastrowlist[-1])
				lastrowlist.append(lastrowlist[-1])
				lastrowlist.append(lastrowlist[-1])
			elif len(lastrowlist)==1:
				lastrowlist.append(lastrowlist[-1])
				lastrowlist.append(lastrowlist[-1])
				lastrowlist.append(lastrowlist[-1])
				lastrowlist.append(lastrowlist[-1])
			# print (lastrowlist)
			tempxx = " ".join(lastrowlist)
			strlastrow = " " + tempxx + '\n'
			# print (strlastrow)
			finalcontent = content[:-1]
			finalcontent.append (strlastrow)
			fout = open(path,'w')
			fout.writelines(finalcontent)
			fout.close()
			print ("Atom #" + str(filename[i]) + " Matrix Makeup")

def main(argv):
	Atom_selection = []
	outputfile_name = ''

	try:
		opts, args = getopt.getopt(argv, "ci:o:",["inputfile=","outputfile="])

	except getopt.GetoptError:
		print ("Error parameter: BaderMerge.py -i <AtomSelection_1> -i <AtomSelection_2> -o <Outputfile_name>")
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-c':
			print ('Copyright CaoZheng @ 2018')
			sys.exit()

		elif opt in ("-i","--inputfile"): #Atom-selection number
			# Atom_selection.append(arg)
			pass

		elif opt in ("-o","--outputfile"):
			outputfile_name = arg

	Atom_selection = makeuplist()
	#Atom_selection should be list ['0001','0002',etc.]
	if Atom_selection == []:
		print ("Please add -i atom-selection!")
	else:
		print ("Atom selection:" + str(Atom_selection))
		if outputfile_name != '':
			check(Atom_selection)
			mergy(Atom_selection,outputfile_name)
		else:
			check(Atom_selection)
			mergy(Atom_selection)


if __name__ == "__main__":
	main(sys.argv[1:])