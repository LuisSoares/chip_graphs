'''
Data file template:
	mean1	IRRELEVANT	stderror1
	mean2	IRRELEVANT	stderror2
	mean3	IRRELEVANT	stderror3
	mean4	IRRELEVANT	stderror4
	mean5	IRRELEVANT	stderror5
	#
	mean1	IRRELEVANT	stderror1
	mean2	IRRELEVANT	stderror2
	mean3	IRRELEVANT	stderror3
	mean4	IRRELEVANT	stderror4
	mean4	IRRELEVANT	stderror5

Currently requires 5 primer pairs in normal orientation look down if different orientation
'''
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from coord import * #coord contains the coordinates of primers of used genes
from scipy.interpolate import spline

#global variables definitions 
primers=('5\'', '.','..','...','3\'')
color_scheme=["#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84",]
line_color=['Blue','Yellow','Red','Orange','#00CC66','Green','Orange','Yellow','Yellow']
os.chdir('/Users/Luis/Desktop')

'''global definitions from command line
	argument 1- should be file name, will be split at dot to extract GENE NAMe
	argument 2- should be a string with IP used
	argument 3- string with different conditions separated by white space
	'''
file_name=sys.argv[1]
Condition=sys.argv[2]
Conditions=sys.argv[3].split()
Gene=file_name.split('.')[0]

coordinates=eval(Gene) #uses the coord defined variables and the parsed Gene name

#Data retrieval
def load_data(ChIP):
	'''Function to Parse the Data
		right now it is not very organized but gets the job done
	'''
	with open(ChIP) as file:
		data=file.read()
	try:
		assert len(data.split('#'))==len(Conditions)
	except:
		print ('Number of conditions in Command line and Data file different')
		sys.exit()
		
	data2=[item.split() for item in data.split('#')]
	means=[] #list of lists (each internal list contains the means for each condition)
	std=[] #list of lists (each internal list contains the stderror for each condition)
	for item in data2:
		means.append([float(item[0]),float(item[3]),float(item[6]),float(item[9]),float(item[12])])
	for item in data2:
		std.append([float(item[2]),float(item[5]),float(item[8]),float(item[11]),float(item[14])])
	return (means,std)
	
	#TODO: get this part as different function
def parse_data(means,std):
	means_per_primer=[] #list of lists (each internal list contains the means for each condition)
	std_per_primer=[]
	n=0
	while n<=4:
		temp=[]
		for item in means:
			temp.append(item[n])
		n+=1
		means_per_primer.append(temp)
	
	n=0
	while n<=4:
		temp=[]
		for item in std:
			temp.append(item[n])
		n+=1
		std_per_primer.append(temp)
	#in case primer order is inverted
	#means_per_primer=means_per_primer[::-1]
	#std_per_primer=std_per_primer[::-1]
	return(means_per_primer,std_per_primer)
	
means,std=load_data(file_name)
means_per_primer,std_per_primer=parse_data(means,std)
	
#To bar graphs	
ind=np.arange(len(Conditions))
width=0.15
ind=ind+0.1
fig = plt.figure(1,figsize=(len(Conditions)*1.5,5))
ax = fig.add_subplot(111, axisbg='0.92')
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth('2')
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_color('1')
ax.spines['bottom'].set_zorder(4)
ax.grid(axis='y',linestyle="-",color='#D8D8D8',linewidth=1)
n=0
for item in means_per_primer:
	ax.bar(ind+width*n,item,width,color=color_scheme[n],zorder=3,yerr=std_per_primer[n],linewidth=1.5,edgecolor='0.92',label=primers[n],error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1,zorder=4))
	n+=1

ax.set_ylim(0,max([max(item) for item in means_per_primer])*1.2)
ax.set_title(Gene,fontsize=20,fontstyle='italic',fontweight="medium",color='0.35')
ax.set_ylabel('IP',fontsize=16,fontweight='medium',color='0.35')
ax.set_xlabel(Condition,fontsize=16,fontweight="medium",color='0.35')
ax.set_xticks(ind+2.5*width)
ax.set_xticklabels( Conditions,fontsize=16,fontstyle='italic',fontweight='light',color='0.35')
ax.xaxis.set_tick_params(pad=8)
ax.yaxis.set_tick_params(labelsize=12, labelcolor='0.25',color='w')
ax.yaxis.tick_left()
plt.tick_params(\
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off' )
plt.legend(('5\'', '.','..','...','3\'') ,bbox_to_anchor=(1.2, 1),frameon=False)
plt.show()

#interpolation
fig = plt.figure(1,figsize=(6,4))
ax = fig.add_subplot(111, axisbg='0.92')
ax.set_title(Gene+' '+Condition,fontsize=18,fontstyle='italic',fontweight="medium",color='0.35')
ax.yaxis.set_tick_params(labelsize=0)
ax.xaxis.set_tick_params(labelsize=14)
plt.subplots_adjust(bottom=0.13)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth('2')
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_color('1')
ax.spines['bottom'].set_zorder(4)
ax.grid(axis='y',linestyle="-",color='#D8D8D8',linewidth=1)
xdata=coordinates
xnew=np.linspace(min(xdata),max(xdata),300)
n=0
ax.set_ylim(0,25)
ax.set_xlim(min(xdata)-300,max(xdata)+300)
for item in means:
	plt.plot(xnew,spline(xdata,item,xnew),color=line_color[n],zorder=4,linewidth=3,alpha=0.80)
	n+=1
plt.tick_params(\
    #axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',
	top='off')
plt.legend(Conditions,frameon=False, bbox_to_anchor=(1.35, 1))
plt.show()
