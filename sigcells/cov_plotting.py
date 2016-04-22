
import matplotlib.pyplot as plt
import numpy as np

ions = ['HI', 'MgII', 'CIV', 'OVI']

def plotfunc(prop,m):

    markList = ['$1$', '$2$', '$3$', '$4$', '$5$', '$6$', '$7$', '$8$', '$9$', '$10$'] 
    markList = ['o', 's', 'x', 'v', '^', '<', '>', '+', '*', 'd'] 
    colorList = ['blue', 'green', 'red', 'black', 'brown', 'cyan', 'crimson', 'azure', 
                    'coral', 'darkgreen']

#    prop = filename.split('_')[2]
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2,2, figsize=(10,10))
    axes = (ax1,ax2,ax3,ax4)
    for ion, ax in zip(ions,axes):
        filename = '{0:s}_clump_{1:s}.out'.format(ion, prop)
        k, cov = [], []
        with open(filename, 'r') as f:
            for line in f:
                l = line.split()
                for i in range(1,len(l)):
                    k.append(int(l[0]))
                    cov.append(float(l[i]))

        filename = '{0:s}_clump_member_count.out'.format(ion)
        number = []
        marks = []
        cols = []
        with open(filename, 'r') as f:
            for line in f:
                memberSum = 0
                l = line.split()
                for i in range(1,len(l)):
                    #number.append(np.log10(float(l[i])))
                    number.append(i)
                    memberSum += int(l[i])
                    marks.append(markList[i-1])
                    cols.append(colorList[i-1])
        for a, b, c, d in zip(k, cov, cols, marks):
            ax.scatter(a,b,marker=d,c=c,s=10)

#        ax.scatter(k,cov,c=number,marker=marks,s=50,cmap='viridis')

        ax.set_xlabel('k')
        if 'cov' in prop:
            p1 = prop.split('_')[0]
            p2 = prop.split('_')[1]
            ax.set_ylabel('{0:s} {1:s} for group'.format(p1,p2))
        else:
            ax.set_ylabel('{0:s} for group'.format(prop))
        ax.set_xlim([0,max(k)+1])
        ax.set_ylim([0,m])
        ax.set_title(ion)
    fig.tight_layout()
    s = 'cluster_{0:s}.png'.format(prop)
    fig.savefig(s,bbox_inches='tight')


props = ['temperature_cov', 'density_cov', 'snII_cov', 'radius']
maxes = [5,25,2.5, 200]
for p,m in zip(props,maxes):
    plotfunc(p,m)

