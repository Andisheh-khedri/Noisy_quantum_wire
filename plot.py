import matplotlib.pyplot as plt
import scipy as sy
import pylab as plb 
import matplotlib as mpl
from pylab import *
from scipy import special
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

import scipy.special as sc

matplotlib.rc('font', family='serif')

matplotlib.rc('text.latex', 
              preamble=[r'\usepackage[T1]{fontenc}',
                        r'\usepackage{amsmath}',
                        r'\usepackage{txfonts}',
                        r'\usepackage{textcomp}'])

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True 
#-----------------------------------------------------

#---------------------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------
#---------------------------------------------------------------
# PROPAGATOR-input data:
#---------------------------------------------------------------

data = plb.loadtxt('prop_s=1_k=0.2.txt')

ws1k0p2= data[:,0]
Gs1k0p2= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1_k=0.6.txt')

ws1k0p6= data[:,0]
Gs1k0p6= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1_k=1.0.txt')

ws1k1= data[:,0]
Gs1k1= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1_k=1.4.txt')

ws1k1p4= data[:,0]
Gs1k1p4= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1_k=1.8.txt')

ws1k1p8= data[:,0]
Gs1k1p8= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=0.8_k=0.2.txt')

ws0p8k0p2= data[:,0]
Gs0p8k0p2= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=0.8_k=0.6.txt')

ws0p8k0p6= data[:,0]
Gs0p8k0p6= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=0.8_k=1.0.txt')

ws0p8k1= data[:,0]
Gs0p8k1= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=0.8_k=1.4.txt')

ws0p8k1p4= data[:,0]
Gs0p8k1p4= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=0.8_k=1.8.txt')

ws0p8k1p8= data[:,0]
Gs0p8k1p8= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1.2_k=0.2.txt')

ws1p2k0p2= data[:,0]
Gs1p2k0p2= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1.2_k=0.6.txt')

ws1p2k0p6= data[:,0]
Gs1p2k0p6= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1.2_k=1.0.txt')

ws1p2k1= data[:,0]
Gs1p2k1= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1.2_k=1.4.txt')

ws1p2k1p4= data[:,0]
Gs1p2k1p4= data[:,1]
#-------------------------------------------
data = plb.loadtxt('prop_s=1.2_k=1.8.txt')

ws1p2k1p8= data[:,0]
Gs1p2k1p8= data[:,1]
#---------------------------------------------------------------
#---------------------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, gridspec_kw={'hspace':0}, figsize=(30, 19))

rcParams['axes.titlepad'] = 20

gs = gridspec.GridSpec(5, 3)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Green's functions:
#---------------------------------------------------------------------------------
ax4= plt.subplot(gs[0:2,0:3])

#-----------------------------------------------------
wc=100.0
#-----------------------------------------------------

s=1 


plt.plot(ws1k1,Gs1k1*(ws1k1*wc), color='red', ms=8, lw=6)


plt.fill_between(ws1k0p6,Gs1k0p6*(ws1k0p6*wc),Gs1k1*(ws1k0p6*wc), color='red', alpha=0.2) 
plt.fill_between(ws1k0p6,Gs1k1*(ws1k0p6*wc),Gs1k1p4*(ws1k0p6*wc), color='red', alpha=0.4) 


plt.text(0.00004, 2,r'Ohmic',fontsize=55,color='red',alpha=0.8)


s=0.8 


plt.plot(ws0p8k1,Gs0p8k1*(ws0p8k1*wc), color='blue', ms=8, lw=6)


plt.fill_between(ws0p8k0p6,Gs0p8k0p6*(ws0p8k0p6*wc),Gs0p8k1*(ws0p8k0p6*wc), color='blue', alpha=0.2) 
plt.fill_between(ws0p8k0p6,Gs0p8k1*(ws0p8k0p6*wc),Gs0p8k1p4*(ws0p8k0p6*wc), color='blue', alpha=0.4) 

plt.text(0.01, 10.0,r'sub-Ohmic',fontsize=55,color='blue',alpha=0.8)




s=1.2 

plt.plot(ws1p2k1,Gs1p2k1*(ws1p2k1*wc), color='green', ms=8, lw=6)


plt.fill_between(ws1p2k0p6,Gs1p2k0p6*(ws1p2k0p6*wc),Gs1p2k1*(ws1p2k0p6*wc), color='green', alpha=0.2) 
plt.fill_between(ws1p2k0p6,Gs1p2k1*(ws1p2k0p6*wc),Gs1p2k1p4*(ws1p2k0p6*wc), color='green', alpha=0.4) 


plt.text(0.01, 0.1,r'super-Ohmic',fontsize=55,color='green',alpha=0.8)

plt.text(0.6, 10.0,r'$(a)$',fontsize=75)


plt.xscale('log')
plt.yscale('log')

ylabel(r'$K(\omega)$',fontsize=75, labelpad=30)
xlabel(r'$\omega/\omega_c$',fontsize=75, labelpad=-40)


plt.xticks(fontsize = 75)
plt.yticks(fontsize = 75)

gca().set_xlim(0.00001,1.2)

gca().set_ylim(0.04,40.0)

plt.tick_params(direction='in', which='major', pad=20, size=30, width=2.5)
plt.tick_params(direction='in', which='minor', pad=20, size=15, width=2.5)

for axis in ['top','bottom','left','right']:
    ax4.spines[axis].set_linewidth(2.5)

#-------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------


ax1= plt.subplot(gs[2:5, 0])



number = 9
cmap = plt.get_cmap('RdYlBu')
colors = [cmap(i) for i in np.linspace(0, 1, number)]

p=np.arange(-3.0,5.0,0.1)
w = (10.0**(p))


style = "Simple, tail_width=3.0, head_width=20, head_length=20"
#-----------------------------------------------------
data = plb.loadtxt('flow_s=1.txt')


plt.plot([1.0,1.0], [0.05,50.0], 'k--', lw=2.0)


n=8600
for i in range(25):
    print i
    lams1= data[(i)*n:(i+1)*n,1]
    Vres1= data[(i)*n:(i+1)*n,2]
    Kres1= data[(i)*n:(i+1)*n,11]
    plt.plot(Kres1,Vres1, color='red', ms=8, lw=8)
    for j in range (1):
        m=(800)*(j+1)
        if np.absolute(Kres1[500]-1.0)>0.1:
           a = patches.FancyArrowPatch((Kres1[m],Vres1[m]), (Kres1[m+1],Vres1[m+1]), **dict(arrowstyle=style, color='red'))
           plt.gca().add_patch(a)
#---------------------------------------------------------------
plt.text(1.5, 20.0,r'$(b)$',fontsize=75)
#-----------------------------------------------------
k=np.arange(1.0,2.0,0.001)
plt.fill_between(k,0.05*k/k,50.0*k/k, color='orange', alpha=0.1) 
k=np.arange(0.001,1.0,0.001)
plt.fill_between(k,0.05*k/k,50.0*k/k, color='darkslategray', alpha=0.1) 
#-----------------------------------------------------

plt.yscale('log')

for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2.5)


ylabel(r'$V(\Lambda)$',fontsize=75, labelpad=30)
xlabel(r'$K(\Lambda)$',fontsize=75, labelpad=30)


plt.xticks([0.5,1.0,1.5],[r'$0.5$',r'$1.0$',r'$1.5$'])


plt.xticks(fontsize = 75)
plt.yticks(fontsize = 75)

gca().set_xlim(0.01,2.0)

gca().set_ylim(0.05,50.0)

plt.tick_params(direction='in', which='major', pad=20, size=30, width=2.5)
plt.tick_params(direction='in', which='minor', pad=20, size=15, width=2.5)


for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2.5)


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
ax2= plt.subplot(gs[2:5, 1])


number = 9
cmap = plt.get_cmap('RdYlBu')
colors = [cmap(i) for i in np.linspace(0, 1, number)]

p=np.arange(-3.0,5.0,0.1)
w = (10.0**(p))


style = "Simple, tail_width=3.0, head_width=20, head_length=20"
#-----------------------------------------------------
data = plb.loadtxt('flow_s=0.8.txt')

n=8600
for i in range(8):
    print i
    lams1= data[(i+1)*n:(i+2)*n,1]
    Vres1= data[(i+1)*n:(i+2)*n,2]
    Kres1= data[(i+1)*n:(i+2)*n,11]
    plt.plot(Kres1,Vres1, color='blue', ms=8, lw=8)
    for j in range (8):
        m=(900)*(j+1)
        a = patches.FancyArrowPatch((Kres1[m],Vres1[m]), (Kres1[m+1],Vres1[m+1]), **dict(arrowstyle=style, color='blue'))
        plt.gca().add_patch(a)



k=np.arange(1.0,2.0,0.001)
plt.fill_between(k,0.05*k/k,50.0*k/k, color='orange', alpha=0.1) 
k=np.arange(0.001,1.0,0.001)
plt.fill_between(k,0.05*k/k,50.0*k/k, color='darkslategray', alpha=0.1) 

#---------------------------------------------------------------
plt.text(1.5, 20.0,r'$(c)$',fontsize=75)
#-----------------------------------------------------
#-----------------------------------------------------

plt.setp(ax2.get_yticklabels(), visible=False)

plt.yscale('log')

for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(2.5)

xlabel(r'$K(\Lambda)$',fontsize=75, labelpad=30)


plt.xticks([0.5,1.0,1.5],[r'$0.5$',r'$1.0$',r'$1.5$'])

plt.xticks(fontsize = 75)
plt.yticks(fontsize = 75)

gca().set_xlim(0.01,2.0)

gca().set_ylim(0.05,50.0)

plt.tick_params(direction='in', which='major', pad=20, size=30, width=2.5)
plt.tick_params(direction='in', which='minor', pad=20, size=15, width=2.5)

for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(2.5)


#--------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
ax3= plt.subplot(gs[2:5, 2])


number = 9
cmap = plt.get_cmap('RdYlBu')
colors = [cmap(i) for i in np.linspace(0, 1, number)]

p=np.arange(-3.0,5.0,0.1)
w = (10.0**(p))


style = "Simple, tail_width=3.0, head_width=20, head_length=20"
#-----------------------------------------------------
data = plb.loadtxt('flow_s=1.2.txt')

n=8600
for i in range(8):
    print i
    lams1= data[(i+1)*n:(i+2)*n,1]
    Vres1= data[(i+1)*n:(i+2)*n,2]
    Kres1= data[(i+1)*n:(i+2)*n,11]
    plt.plot(Kres1,Vres1, color='green', ms=8, lw=8)
    for j in range (8):
        m=(600)*(j+1)
        a = patches.FancyArrowPatch((Kres1[m],Vres1[m]), (Kres1[m+1],Vres1[m+1]), **dict(arrowstyle=style, color='green'))
        plt.gca().add_patch(a)
#---------------------------------------------------------------
plt.text(1.6, 20.0,r'$(d)$',fontsize=75)
#-----------------------------------------------------
k=np.arange(1.0,2.0,0.001)
plt.fill_between(k,0.05*k/k,50.0*k/k, color='orange', alpha=0.1) 
k=np.arange(0.001,1.0,0.001)
plt.fill_between(k,0.05*k/k,50.0*k/k, color='darkslategray', alpha=0.1) 
#-----------------------------------------------------

plt.setp(ax3.get_yticklabels(), visible=False)


plt.yscale('log')

for axis in ['top','bottom','left','right']:
    ax3.spines[axis].set_linewidth(2.5)

xlabel(r'$K(\Lambda)$',fontsize=75, labelpad=30)


plt.xticks([0.5,1.0,1.5],[r'$0.5$',r'$1.0$',r'$1.5$'])



plt.xticks(fontsize = 75)
plt.yticks(fontsize = 75)

gca().set_xlim(0.01,2.0)

gca().set_ylim(0.05,50.0)

plt.tick_params(direction='in', which='major', pad=20, size=30, width=2.5)
plt.tick_params(direction='in', which='minor', pad=20, size=15, width=2.5)

for axis in ['top','bottom','left','right']:
    ax3.spines[axis].set_linewidth(2.5)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

plt.subplots_adjust(left=0.15, right=0.99, top=0.95, bottom=0.16)


f.tight_layout(h_pad=0.1, w_pad=-0.5)
f.subplots_adjust(wspace=0,hspace=0.7)

savefig('flow.pdf')
show()
