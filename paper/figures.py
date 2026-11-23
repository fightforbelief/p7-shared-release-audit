"""Render the returned Windows result tables for the manuscript."""
import csv
from pathlib import Path
from statistics import median
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
with (HERE.parent/"results/profiles.csv").open() as f:
    rows = list(csv.DictReader(f))
fig, axes = plt.subplots(1,3,figsize=(7.3,2.65),sharex=True)
for ax,size in zip(axes,("small","medium","large")):
    for method,label,style,color in (("independent_","Separate tables","--","#ac5928"),
                                     ("centered_","Centered bound",":","#655397"),
                                     ("","Shared table","-","#006875")):
        xs = (0,.01,.1,.25)
        vals = [median(100*(float(r[method+"hi"])-float(r[method+"lo"])) for r in rows
                       if r["size"]==size and int(r["step"])==10 and float(r["rho"])==x) for x in xs]
        ax.plot(xs,vals,style,marker='o',markersize=3,lw=1.7,color=color,label=label)
    ax.set_title(size.capitalize()+' institutions',fontsize=10)
    ax.set_xticks([0,.1,.25],['0','10','25'])
    ax.set_ylim(bottom=0)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.2)
axes[0].set_ylabel("Median width (pp)",fontsize=9)
fig.supxlabel("Replacement share (%)",y=.14,fontsize=10)
fig.legend(*axes[0].get_legend_handles_labels(),loc="lower center",ncol=3,frameon=False,fontsize=9)
fig.subplots_adjust(left=.095,right=.985,bottom=.33,top=.88,wspace=.36)
fig.savefig(HERE/"profiles.pdf",metadata={"CreationDate":None})
