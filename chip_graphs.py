import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
import re
import sys
from coord import *
from scipy.interpolate import spline

primers=('5\'', '.','..','...','3\'')
color_scheme=["#c7e9b4",
"#7fcdbb",
"#41b6c4",
"#1d91c0",
"#225ea8",
"#0c2c84",]
os.chdir('/Users/Luis/Desktop')
file_name=sys.argv[1]
Condition=sys.argv[2]
Conditions=sys.argv[3].split()
Gene=file_name.split('.')[0]
coordinates=eval(Gene)
#print(coordinates)
file=open(file_name)
data=file.read()
file.close()
N=data.count('#')
data_split=data.split('#')
data2=[item.split() for item in data_split[:-1]]
#print (data2,len(data2),N)
means=[]
std=[]
for item in data2:
	means.append([float(item[0]),float(item[3]),float(item[6]),float(item[9]),float(item[12])])
for item in data2:
	std.append([float(item[2]),float(item[5]),float(item[8]),float(item[11]),float(item[14])])
#print (means,std)
means_per_primer=[]
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
ind=np.arange(N)
width=0.15
ind=ind+0.1
fig = plt.figure(1,figsize=(N*1.5,5))
ax = fig.add_subplot(111, axisbg='0.92')
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth('2')
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_color('1')
ax.spines['bottom'].set_zorder(4)
ax.grid(axis='y',linestyle="-",color='#D8D8D8',linewidth=2)
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
plt.legend(('5\'', '.','..','...','3\'') ,bbox_to_anchor=(1.3, 1),frameon=False)
plt.show()

#interpolation
fig,ax=plt.subplots(figsize=(6,4))
ax.set_title(Gene,fontsize=24,fontstyle='italic',fontweight="medium")
ax.yaxis.set_tick_params(labelsize=0)
ax.xaxis.set_tick_params(labelsize=14)
plt.subplots_adjust(bottom=0.13)
ax.spines['bottom'].set_linewidth('2')
ax.spines['left'].set_linewidth('2')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
line_color=['#e41a1c',
'#377eb8',
'#4daf4a',
'#984ea3',
'#ff7f00',
'#ffff33',
'#a65628',
'#f781bf',
'#999999']
xdata=coordinates
xnew=np.linspace(min(xdata),max(xdata),300)
n=0
ax.set_xlim(min(xdata)-300,max(xdata)+300)
for item in means:
	plt.plot(xnew,spline(xdata,item,xnew),color=line_color[n],linewidth=3,alpha=0.80)
	n+=1
plt.tick_params(\
    #axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',
	top='off')
plt.legend(Conditions,frameon=False, bbox_to_anchor=(1.3, 1))
plt.show()
